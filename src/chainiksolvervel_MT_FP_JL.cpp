#include <simple_chain_ik/solvers/chainiksolvervel_MT_FP_JL.hpp>
#include <iostream>
#include <limits>

#include "kdl/utilities/svd_eigen_HH.hpp"

// ATTENTION: these have to be removed as soon as KDL released version will be >=1.4.0
#define E_SIZE_MISMATCH -4
#define E_SVD_FAILED -8

#define DEBUG 2

using namespace KDL;

ChainIkSolverVel_MT_FP_JL::ChainIkSolverVel_MT_FP_JL(const Chain& chain, double eps, int maxiter) : 
    chain(chain),
    eps(eps),
    maxiter(maxiter),
    nj(JS_dim),
    ts_dim(TS_dim),
    jnt2jac(chain),
    weightTS(Eigen::Matrix<double,TS_dim,TS_dim>::Identity()),
    weightJS(Eigen::Matrix<double,JS_dim,JS_dim>::Identity()),
    lambda(0.0),
    svdResult(0),
    q_dot_lb(Eigen::Matrix<double,JS_dim,1>::Ones()*std::numeric_limits<double>::min()),
    q_dot_ub(Eigen::Matrix<double,JS_dim,1>::Ones()*std::numeric_limits<double>::max()),
    is_jac_weighted(false),
    q_lb(Eigen::Matrix<double,JS_dim,1>::Ones()*std::numeric_limits<double>::min()),
    q_ub(Eigen::Matrix<double,JS_dim,1>::Ones()*std::numeric_limits<double>::min()), 
    svdU(Eigen::Matrix<double,TS_dim,TS_dim>::Zero()),
    svdV(Eigen::Matrix<double,JS_dim,JS_dim>::Zero()),
    svdS(Eigen::Matrix<double,TS_dim,1>::Zero())
{
    assert(nj == chain.getNrOfJoints());
}

void ChainIkSolverVel_MT_FP_JL::setJointLimits(const Eigen::ArrayXd& lower_bound, const Eigen::ArrayXd& upper_bound)
{
    // either asser here or do everything in header file
    assert(lower_bound.rows() == nj);
    assert(upper_bound.rows() == nj);
    
    q_lb = lower_bound;
    q_ub = upper_bound;
}

int ChainIkSolverVel_MT_FP_JL::setWeightJS(const Eigen::MatrixXd& Mq)
{
    // either asser here or do everything in header file
    assert(Mq.rows() == nj);
    assert(Mq.cols() == nj);
    
    if( Mq == Eigen::Matrix<double,JS_dim,JS_dim>::Identity() )
    {
        is_jac_weighted = false;
        return 0;
    }
    
    weightJS = Mq;
    is_jac_weighted = true;
    return 0;
}

int ChainIkSolverVel_MT_FP_JL::setWeightTS(const Eigen::MatrixXd& weights)
{
    // either asser here or do everything in header file
    assert(weights.rows() == ts_dim);
    assert(weights.cols() == ts_dim);
    
    double w(1.0),last_w(0.0);
    uint task_nr(0);
    Eigen::Matrix<double,TS_dim,1> task_list(Eigen::Matrix<double,TS_dim,1>::Zero());
    Eigen::Matrix<double,TS_dim,1> Mx(weights.diagonal());
    
    // find current most-important task weight
    w = Mx.maxCoeff();
#if DEBUG>1
    std::cout << "Mx = " << Mx.transpose() << std::endl;
    std::cout << "w = " << w << std::endl;
#endif
    
    // base case: no other task is present
    while(std::abs(w) > 1e-10)
    {
        // task_list will have 1's in 1st task elements, 2's in 2nd task elements...
        task_nr++;
        task_list = (Mx.array() == w).select(Eigen::Matrix<double,TS_dim,1>::Ones()*task_nr,task_list);
        
        // adjust weights for a new cycle
        Mx = (Mx.array() == w).select(Eigen::Matrix<double,TS_dim,1>::Zero(),Mx);
        
        // find current most-important task weight
        w = Mx.maxCoeff();
        
        // debug output
#if DEBUG>1
        std::cout << "task_list = " << task_list.transpose() << std::endl;
        std::cout << "Mx = " << Mx.transpose() << std::endl;
        std::cout << "w = " << w << std::endl;
#endif
    }
    
#if DEBUG>1
    std::cout << "task_nr = " << task_nr << std::endl;
    #endif
    task_nr_ = task_nr;
    task_list_ = task_list;
    return 0;
}

const char* ChainIkSolverVel_MT_FP_JL::strError(const int error) const
{
    if (E_CONVERGE_PINV_SINGULAR == error) return "Converged put pseudo inverse of jacobian is singular.";
    else return SolverI::strError(error);
}

int ChainIkSolverVel_MT_FP_JL::CartToJnt(const JntArray& q_in, const Twist& v_in, JntArray& qdot_out)
{
    SetToZero(qdot_out);
    return 0;
    
    if(nj != q_in.rows() || nj != qdot_out.rows())
        return (error = E_SIZE_MISMATCH);
    error = jnt2jac.JntToJac(q_in,jac_kdl);
    if ( error < E_NOERROR)
        return error;
    
    // initializations
    // total solution
    Eigen::Vector<double,JS_dim,1> S_0(Eigen::Vector<double,JS_dim,1>::Zero());
    // solution at step k-1
    Eigen::Vector<double,JS_dim,1> S_k_old(Eigen::Vector<double,JS_dim,1>::Zero());
    // solution at step k
    Eigen::Vector<double,JS_dim,1> S_k(Eigen::Vector<double,JS_dim,1>::Zero());
    // null-space projector at step k
    Eigen::Matrix<double,JS_dim,JS_dim> N_k(Eigen::Matrix<double,JS_dim,JS_dim>::Identity());
    // complete task specification
    Eigen::Matrix<double,TS_dim,1> xi(Eigen::Matrix<double,TS_dim,1>::Zero());
    
    //S_0.setZero();
    //S_k_old.setZero();
    //S_k.setZero();
    //N_k.setIdentity();
    jac = jac_kdl.data;
    for(int i=0; i<TS_dim; ++i)
        xi(i) = v_in(i);
    
    // for each task
    for(uint k=0; k<task_nr_; ++k)
    {
        // jacobian at k-th step
        Eigen::MatrixXd J_k;
        // task at k-th step
        Eigen::VectorXd xi_k;
        // null-projected jacobian at k-th step
        Eigen::MatrixXd NJ_k;
        // pseudo-inverted null-projected jacobian at k-th step
        Eigen::MatrixXd NJ_k_pinv;
        
        // get local jacobian (dimension may vary for each task...)
        selectMatrixRows(task_list_,k+1,jac,J_k);
        // get local task (dimension may vary for each task...)
        selectMatrixRows(task_list_,k+1,xi,xi_k);
        
        // comupute null-projected jacobian and its Pseudo-Inverse
        NJ_k.noalias() = J_k * N_k;
        pinvDLS(NJ_k,NJ_k_pinv);
        
        // compute q_dot at this task
        S_k_old = S_k;
        S_k = NJ_k_pinv*(xi_k - J_k * S_k_old);
        
        // check limits based on q_in
        int worst_joint = -1;
        worst_joint = checkVelocityLimits(q_in,S_0 + S_k);
        
        // SNS - De Luca : dispense pag. 75+
        // enforce limits - if failed, don't go on with tasks
        if(worst_joint >= 0)
        {
            if(!enforceLimits())
            {
                double s;
                s = computeMaxScaling();
                S_0 += NJ_k_pinv*(s*xi_k - J_k * S_k_old);
                // return some inaccuracy error, but do not fail
                return E_SNS_NEEDED;
            }
        }
        
        // if it's not last task, compute iterative step
        S_0 += S_k;
        N_k -= NJ_k_pinv*NJ_k;
    }
    
    // return q_dot
    qdot_out.data = S_0;
    return E_NOERROR;
    
    ///////////////// OLD SOLVER
    
    Eigen::MatrixXd tmp_jac_weight1;
    Eigen::MatrixXd tmp_jac_weight2;
    Eigen::MatrixXd tmp_ts;
    Eigen::MatrixXd tmp_js;
    double lambda_scaled;
    unsigned int nrZeroSigmas ;
    double sigmaMin;
    Eigen::VectorXd tmp(VectorXd::Zero(nj));
    
    if(nj != q_in.rows() || nj != qdot_out.rows())
        return (error = E_SIZE_MISMATCH);
    error = jnt2jac.JntToJac(q_in,jac_kdl);
    if ( error < E_NOERROR) return error;
    
    double sum;
    unsigned int i,j;
    
    // Initialize (internal) return values
    nrZeroSigmas = 0 ;
    sigmaMin = 0.;
    lambda_scaled = 0.;
    
    /*
     *        for (i=0;i<jac_kdl.rows();i++) {
     *            for (j=0;j<jac_kdl.columns();j++)
     *                tmp_jac(i,j) = jac(i,j);
}
*/
    
    // Create the Weighted jacobian
    tmp_jac_weight1.noalias() = jac_kdl.data * weightJS;
    tmp_jac_weight2.noalias() = weightTS.asDiagonal() * tmp_jac_weight1;
    
    // Compute the SVD of the weighted jacobian
    Eigen::MatrixXd svdU(MatrixXd::Zero(6,nj));
    Eigen::VectorXd svdS(VectorXd::Zero(nj));
    Eigen::MatrixXd svdV(MatrixXd::Zero(nj,nj));
    svdResult = svd_eigen_HH(tmp_jac_weight2,svdU,svdS,svdV,tmp,maxiter);
    if (0 != svdResult)
    {
        qdot_out.data.setZero() ;
        return (error = E_SVD_FAILED);
    }
    
    //Pre-multiply U and V by the task space and joint space weighting matrix respectively
    tmp_ts.noalias() = weightTS.asDiagonal() * svdU.topLeftCorner(6,6);
    tmp_js = weightJS.lazyProduct(svdV);
    
    // Minimum of six largest singular values of J is S(5) if number of joints >=6 and 0 for <6
    if ( jac_kdl.columns() >= 6 ) {
        sigmaMin = svdS(5);
    }
    else {
        sigmaMin = 0.;
    }
    
    // tmp = (Si*U'*Ly*y),
    for (i=0;i<jac_kdl.columns();i++) {
        sum = 0.0;
        for (j=0;j<jac_kdl.rows();j++) {
            if(i<6)
                sum+= tmp_ts(j,i)*v_in(j);
            else
                sum+=0.0;
        }
        // If sigmaMin > eps, then wdls is not active and lambda_scaled = 0 (default value)
        // If sigmaMin < eps, then wdls is active and lambda_scaled is scaled from 0 to lambda
        // Note:  singular values are always positive so sigmaMin >=0
        if ( sigmaMin < eps )
        {
            lambda_scaled = sqrt(1.0-(sigmaMin/eps)*(sigmaMin/eps))*lambda ;
        }
        if(fabs(svdS(i))<eps) {
            if (i<6) {
                // Scale lambda to size of singular value sigmaMin
                tmp(i) = sum*((svdS(i)/(svdS(i)*svdS(i)+lambda_scaled*lambda_scaled)));
            }
            else {
                tmp(i)=0.0;  // S(i)=0 for i>=6 due to cols>rows
            }
            //  Count number of singular values near zero
            ++nrZeroSigmas ;
        }
        else {
            tmp(i) = sum/svdS(i);
        }
    }
    
    /*
     *        // x = Lx^-1*V*tmp + x
     *        for (i=0;i<jac_kdl.columns();i++) {
     *            sum = 0.0;
     *            for (j=0;j<jac_kdl.columns();j++) {
     *                sum+=tmp_js(i,j)*tmp(j);
}
qdot_out(i)=sum;
}
*/
    qdot_out.data=tmp_js.lazyProduct(tmp);
    
    // If number of near zero singular values is greater than the full rank
    // of jac, then wdls is active
    if ( nrZeroSigmas > (jac_kdl.columns()-jac_kdl.rows()) )    {
        return (error = E_CONVERGE_PINV_SINGULAR);  // converged but pinv singular
    } else {
        return (error = E_NOERROR);                 // have converged
    }
}

int ChainIkSolverVel_MT_FP_JL::pinvDLS(const Matrix< double, -1, 7 >& NJ_k, Matrix< double, 7, -1 >& NJ_k_pinv)
{
    Eigen::MatrixXd weigthedJ;
    if(is_jac_weighted)
        weigthedJ.noalias() = NJ_k * weightJS;
    else
        weigthedJ = NJ_k;
    
    // Compute the SVD of the weighted jacobian
    // TODO: understand why this doesn't work properly with regular, non-dynamic matrixes
    Eigen::MatrixXd svdU(MatrixXd::Zero(6,nj));
    Eigen::VectorXd svdS(VectorXd::Zero(nj));
    Eigen::MatrixXd svdV(MatrixXd::Zero(nj,nj));
    Eigen::VectorXd tmp(VectorXd::Zero(nj));
    svdResult = svd_eigen_HH(weigthedJ,svdU,svdS,svdV,tmp,maxiter);
    if (0 != svdResult)
        return (error = E_SVD_FAILED);
    
    // TODO finish computing SVD... or we probably want to compute qdot here?
}

