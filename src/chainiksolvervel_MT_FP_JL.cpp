#include <simple_chain_ik/solvers/chainiksolvervel_MT_FP_JL.hpp>
#include <iostream>
#include <limits>
// #include <eigen3/Eigen/src/Core/Matrix.h>
// #include <eigen3/Eigen/src/SVD/JacobiSVD.h>
// #include <eigen3/Eigen/Jacobi>

#include "kdl/utilities/svd_eigen_HH.hpp"

// ATTENTION: these have to be removed as soon as KDL released version will be >=1.4.0
#define E_SIZE_MISMATCH -4
#define E_SVD_FAILED -8

#define DEBUG 3
#define CLASS_NAMESPACE "ChainIkSolverVel_MT_FP_JL::"

#define QDOT_ZERO 0e-6

using namespace KDL;

ChainIkSolverVel_MT_FP_JL::ChainIkSolverVel_MT_FP_JL(const Chain& chain, double eps, int maxiter) : 
    chain(chain),
    eps(eps),
    maxiter(maxiter),
    nj(JS_dim),
    jac_kdl(nj),
    ts_dim(TS_dim),
    jnt2jac(chain),
    weightTS(MatrixT::Identity()),
    weightW(MatrixJ::Identity()),
    weightJS(MatrixJ::Identity()),
    lambda(0.0),
    q_dot_lb(VectorJ::Ones()*(-1.0*std::numeric_limits<double>::max())),
    q_dot_ub(VectorJ::Ones()*std::numeric_limits<double>::max()),
    is_jac_weighted(false),
    q_lb(VectorJ::Ones()*(-1.0*std::numeric_limits<double>::max())),
    q_ub(VectorJ::Ones()*std::numeric_limits<double>::max()),
    S_k_old(VectorJ::Zero()),
    S_k(VectorJ::Zero()),
    N_k(MatrixJ::Identity()),
    xi(VectorT::Zero())
{
    assert(nj == chain.getNrOfJoints());
}

void ChainIkSolverVel_MT_FP_JL::setJointLimits(const VectorJ& lower_bound, const VectorJ& upper_bound)
{
//     // either asser here or do everything in header file
//     assert(lower_bound.rows() == nj);
//     assert(upper_bound.rows() == nj);
    q_lb = lower_bound;
    q_ub = upper_bound;
}

int ChainIkSolverVel_MT_FP_JL::setWeightJS(const MatrixJ& Mq)
{
//     // either asser here or do everything in header file
//     assert(Mq.rows() == nj);
//     assert(Mq.cols() == nj);
    
    if( Mq == MatrixJ::Identity() )
    {
        is_jac_weighted = false;
        return 0;
    }
    
    weightJS = Mq;
    is_jac_weighted = true;
    return 0;
}

int ChainIkSolverVel_MT_FP_JL::setWeightTS(const MatrixT& weights)
{
//     // either asser here or do everything in header file
//     assert(weights.rows() == ts_dim);
//     assert(weights.cols() == ts_dim);
    
    double w(1.0),last_w(0.0);
    uint task_nr(0);
    task_list_.setZero();
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
        task_list_ = (Mx.array() == w).select(VectorTi::Ones()*task_nr,task_list_);
        
        // adjust weights for a new cycle
        Mx = (Mx.array() == w).select(Eigen::Matrix<double,TS_dim,1>::Zero(),Mx);
        
        // find current most-important task weight
        w = Mx.maxCoeff();
        
        // debug output
#if DEBUG>1
        std::cout << "task_list = " << task_list_.transpose() << std::endl;
        std::cout << "Mx = " << Mx.transpose() << std::endl;
        std::cout << "w = " << w << std::endl;
#endif
    }
    
#if DEBUG>1
    std::cout << "task_nr = " << task_nr << std::endl;
#endif
    task_nr_ = task_nr;
    return 0;
}

const char* ChainIkSolverVel_MT_FP_JL::strError(const int error) const
{
    if (E_CONVERGE_PINV_SINGULAR == error) return "Converged put pseudo inverse of jacobian is singular.";
    else return SolverI::strError(error);
}

int ChainIkSolverVel_MT_FP_JL::CartToJnt(const JntArray& q_in, const Twist& v_in, JntArray& qdot_out)
{
    if(nj != q_in.rows() || nj != qdot_out.rows())
        return (error = E_SIZE_MISMATCH);
    error = jnt2jac.JntToJac(q_in,jac_kdl);
    if ( error < E_NOERROR)
    {
        std::cout << CLASS_NAMESPACE << __func__ << "@" << __LINE__ << " : error=" << error << std::endl;
        return error;
    }
    
    // initializations
    S_k_old.setZero();
    S_k.setZero();
    N_k.setIdentity();
    jac = jac_kdl.data;
    for(int i=0; i<TS_dim; ++i)
        xi(i) = v_in(i);
    
    // for each task
    for(uint k=0; k<task_nr_; ++k)
    {
        // jacobian at k-th step
        MatrixXJ J_k;
        // task at k-th step
        VectorX xi_k;
        // null-projected jacobian at k-th step
        MatrixXJ NJ_k;
        // pseudo-inverted null-projected jacobian at k-th step
        MatrixJX NJ_k_pinv;
        
        // get local jacobian (dimension may vary for each task...)
        selectMatrixRows(task_list_,k+1,jac,J_k);
        // get local task (dimension may vary for each task...)
        selectMatrixRows(task_list_,k+1,xi,xi_k);
        
        // comupute null-projected jacobian and its Pseudo-Inverse
        NJ_k.resize(J_k.rows(),NJ_k.ColsAtCompileTime);
        NJ_k.noalias() = J_k * N_k;
        pinvDLS(NJ_k,NJ_k_pinv);
        
        // compute q_dot at this task
        S_k_old = S_k;
        S_k.noalias() = NJ_k_pinv*(xi_k - J_k * S_k_old);
        S_k += S_k_old;
        
        // check limits based on q_in
        bool respecting_limits = false;
        respecting_limits = checkVelocityLimits(q_in.data,S_k);
        
        // Saturation in the Null-Space (SNS) - De Luca : class notes pag. 75+
        // initialize variables to perform SNS
        // weight matrix
        weightW.setIdentity();
        // current task scale factor
        double s = 1.0;
        // largest scale factor found
        double sStar = 0.0;
        // velocity vector to saturate, used in internal loop
        VectorJ qN = S_k_old;
        // saturated velocities in the best task scaling scenario
        VectorJ qNstar;
        // weight matrix in the best task scaling scenario
        MatrixJ weightWstar;
        // pseudo-inverse of jacobian projected in the null-space with saturations (initially equal to non-saturated pseudo-inverse)
        MatrixJX JNW_k_pinv = NJ_k_pinv;
        
        // enforce limits - if failed, don't go on with tasks
        while(!respecting_limits)
        {
            // ATTENTION: not sure about this > is it 
            // a = pinv(J_k*N_k*W)*xi_k OR a = S_k_old + pinv(J_k*N_k*W)*(xi - Jk*S_k_old) ? I will opt for the second, it makes more sense
            VectorJ a;
            pinvDLS(NJ_k * weightW, JNW_k_pinv);
            a.noalias() = JNW_k_pinv*(xi_k - J_k*S_k_old);
            a += S_k_old;
            // compute b - S_k is the variable which gets updated also in the internal cycle
            VectorJ b = S_k - a; 
            
            s = computeMaxScaling(q_in.data,a,b,&worst_joint);
            if(s > sStar)
            {
                sStar = s;
                weightWstar = weightW;
                qNstar = qN;
            }
            
            weightW(worst_joint,worst_joint) = 0.0;
            // apply in qN the limit resulting from S_k
            if(S_k(worst_joint) >= q_dot_ub(worst_joint))
                qN(worst_joint) = q_dot_ub(worst_joint);
            else if(S_k(worst_joint) <= q_dot_lb(worst_joint))
                qN(worst_joint) = q_dot_lb(worst_joint);
            else
            {
                std::cout << CLASS_NAMESPACE << __func__ << " : This is really strange: shouldn't happen!" << std::endl;
                assert(false);
            }
            
            // the task went out of scope of the jacobian (deficient rank): scale it and exit
            
            if( (NJ_k * weightW).colPivHouseholderQr().rank() < xi_k.rows() )
            {
                pinvDLS(NJ_k * weightWstar, JNW_k_pinv);
                S_k = qNstar + JNW_k_pinv*(s*xi_k - J_k * qNstar);
#if DEBUG>1
                std::cout << CLASS_NAMESPACE << __func__ << " : exiting after enforcing joint limits at task #" << k+1 << " out of " << task_nr_ << std::endl;
#endif
                // return some inaccuracy error, but do not fail
                qdot_out.data = S_k;
                return E_SNS_NEEDED;
            }
            
            // compute an extra step and check again for limits
            // if I didn't return, this is the same as above >> pinvDLS(NJ_k * weightW, JNW_k_pinv);
            S_k = qN + JNW_k_pinv*(xi_k - J_k * qN);
            
            respecting_limits = checkVelocityLimits(q_in.data,S_k);
        }
        
        // if it's not last task, compute iterative step
        N_k -= NJ_k_pinv*NJ_k;
    }
    
    // return q_dot
    qdot_out.data = S_k;
    return E_NOERROR;
}

int ChainIkSolverVel_MT_FP_JL::pinvDLS(const MatrixXJ& NJ_k, MatrixJX& NJ_k_pinv)
{
    MatrixXJ weigthedJ;
    weigthedJ.resize(NJ_k.rows(),weigthedJ.ColsAtCompileTime);
    if(is_jac_weighted)
        weigthedJ.noalias() = NJ_k * weightJS;
    else
        weigthedJ = NJ_k;
    
    int zero_sv_counter = 0;
    
    // Compute the SVD of the weighted jacobian
    JacobiSVD<MatrixXd> svd(weigthedJ, ComputeThinU | ComputeThinV);
    Eigen::VectorXd sigma = svd.singularValues();
    
    double lambda_scaled = 0.0;
    double sigmaMin = sigma.minCoeff();
    if ( sigmaMin < eps )
        lambda_scaled = sqrt(1.0-(sigmaMin/eps)*(sigmaMin/eps))*lambda ;
    for(int i=0; i<sigma.rows(); ++i)
    {
        // Damp the singular value
        if(fabs(sigma(i))<eps)
        {
            ++zero_sv_counter;
            sigma(i) = (sigma(i)/(sigma(i)*sigma(i)+lambda_scaled*lambda_scaled));
        }
        else
            sigma(i) = 1.0/sigma(i);
    }
    
    // resize the output matrix and compute the Pinv
    NJ_k_pinv.resize(JS_dim,sigma.rows());
    NJ_k_pinv.noalias() = svd.matrixV() * sigma.asDiagonal() * svd.matrixU().transpose();
    
    return zero_sv_counter;
}

bool ChainIkSolverVel_MT_FP_JL::checkVelocityLimits(const VectorJ& q_in, const VectorJ& q_dot_k)
{
    // condition to satisfy:
    // lb <= q+qd <= ub
    
    int ru,cu,rl,cl;
    double maxu,maxl;
    maxu = ((q_in + q_dot_k) - q_ub).maxCoeff(&ru,&cu);
    maxl = (q_lb - (q_in + q_dot_k)).maxCoeff(&rl,&cl);
    
    // if either one is positive (greater than a small number treated as zero), return the index of the worst one
    if ((maxu > QDOT_ZERO) || (maxl > QDOT_ZERO))
        return false;
    
    // none is positive, return -1
    return true;
}

double ChainIkSolverVel_MT_FP_JL::computeMaxScaling(const VectorJ& q_in, const VectorJ& a, const VectorJ& b, int* r)
{
    updateVelocityLimits(q_in);
    
    VectorJ sMin, sMax;
    // the use of a for cycle is maybe the best thing here...
    for(int i = 0; i<JS_dim; ++i)
    {
        sMin(i) = (q_dot_lb(i) - b(i)) / a(i);
        sMax(i) = (q_dot_ub(i) - b(i)) / a(i);
        if(sMin(i) > sMax(i))
            std::swap(sMin(i),sMax(i));
    }
    
    double smin, smax;
    int c;
    smax = sMax.minCoeff(r,&c);
    smin = sMin.maxCoeff();
    
    if((smin > smax) || (smax < 0.0) || (smin > 1.0))
        return 0.0;
    else
        return smax;
}

void ChainIkSolverVel_MT_FP_JL::updateVelocityLimits(const VectorJ& q_in)
{
    // this is not the most generic way of doing this...
    // would it make sense to consider also velocity and acceleration limits of the real robot?
    q_dot_lb = q_lb - q_in;
    q_dot_ub = q_ub - q_in;
}

void ChainIkSolverVel_MT_FP_JL::selectMatrixRows(const VectorTi& task_list_, uint k, const MatrixTJ& jac, MatrixXJ& jac_k) const
{
    VectorTi indexer = (task_list_.array() == k).select(VectorTi::Ones(),VectorTi::Zero());
    uint task_dim = indexer.sum();
#if DEBUG > 3
    std::cout << CLASS_NAMESPACE << __func__ << " : task #" << k << " | indexer = [" << indexer.transpose() << "] | task_dim = " << task_dim << std::endl;
#endif
    jac_k.resize(task_dim,JS_dim);
    uint jac_row = 0;
    
    for(int i=0; i<TS_dim; ++i)
    {
        if(indexer(i))
            jac_k.row(jac_row++) = jac.row(i);
    }
    
#if DEBUG > 3
    std::cout << "jacobian : " << std::endl << jac << std::endl;
    std::cout << "jacobian([" << task_list_.transpose() << "] == " << k << ") : " << std::endl << jac_k << std::endl;
#endif
}

void ChainIkSolverVel_MT_FP_JL::selectMatrixRows(const VectorTi& task_list_, uint k, const VectorT& xi, VectorX& xi_k) const
{
    VectorTi indexer = (task_list_.array() == k).select(VectorTi::Ones(),VectorTi::Zero());
    uint task_dim = indexer.sum();
#if DEBUG > 2
    std::cout << CLASS_NAMESPACE << __func__ << " : task #" << k << " | indexer = [" << indexer.transpose() << "] | task_dim = " << task_dim << std::endl;
#endif
    xi_k.resize(task_dim,1);
    uint xi_row = 0;
    
    for(int i=0; i<TS_dim; ++i)
    {
        if(indexer(i))
            xi_k(xi_row++) = xi(i);
    }
    
#if DEBUG > 3
    std::cout << "xi = [" << xi.transpose() << "]" << std::endl;
    std::cout << "xi([" << task_list_.transpose() << "] == " << k << ") : " << xi_k.transpose() << std::endl;
#endif
}
