#include <simple_chain_ik/solvers/chainiksolvervel_MT_FP_JL.hpp>
#include <iostream>
#include <limits>
// #include <eigen3/Eigen/src/Core/Matrix.h>
// #include <eigen3/Eigen/src/SVD/JacobiSVD.h>
// #include <eigen3/Eigen/Jacobi>
#include <ros/ros.h>

#include "kdl/utilities/svd_eigen_HH.hpp"

// ATTENTION: these have to be removed as soon as KDL released version will be >=1.4.0
#define E_SIZE_MISMATCH -4
#define E_SVD_FAILED -8

#define DEBUG 2
#define CLASS_NAMESPACE "ChainIkSolverVel_MT_FP_JL::"

// value used when checking for a value to be zero: the smaller the lambda, the higher this should be
// (found to be working empirically when QDOT_ZERO * lambda = 1e-10)
#define QDOT_ZERO 1e-5
// task weight margin to avoid oscillations in velocity commands
#define S_margin 0.0
// choose whether to scale only the task or also the contribution to the last task given by the previous (k-1) tasks
#define SCALE_PREVIOUS_TASK_CONTRIBUTION 0

using namespace KDL;

ChainIkSolverVel_MT_FP_JL::ChainIkSolverVel_MT_FP_JL(const Chain& chain, double eps, int maxiter) : 
    chain(chain),
    eps(eps),
    maxiter(maxiter),
    nj(JS_dim),
    jac_kdl(nj),
    ts_dim(TS_dim),
    fksolver(chain),
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
    xi(VectorT::Zero()),
    model_tolerance_(0.0),
    use_ee_task_(false)
{
    assert(nj == chain.getNrOfJoints());
}

void ChainIkSolverVel_MT_FP_JL::setJointLimits(const VectorJ& lower_bound, const VectorJ& upper_bound)
{
//     // either assert here or do everything in header file
//     assert(lower_bound.rows() == nj);
//     assert(upper_bound.rows() == nj);
    q_lb = lower_bound;
    q_ub = upper_bound;
}

int ChainIkSolverVel_MT_FP_JL::setWeightJS(const MatrixJ& Mq)
{
//     // either assert here or do everything in header file
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
//     // either assert here or do everything in header file
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
    
    KDL::Twist v_ref(v_in);
    // check whether a change of frame is needed, and in such case perform it
    if(use_ee_task_)
    {
        KDL::Frame F;
        fksolver.JntToCart(q_in,F);
        v_ref = F.M.Inverse() * v_in;
        jac_kdl.changeBase(F.M.Inverse());
    }
    
    // initializations
    S_k_old.setZero();
    S_k.setZero();
    N_k.setIdentity();
    jac = jac_kdl.data;
    for(int i=0; i<TS_dim; ++i)
        xi(i) = v_ref(i);
    updateVelocityLimits(q_in.data);
    
    // count the number of iterations we are using here
    int counter = 0;
    
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
        S_k = S_k_old + NJ_k_pinv*(xi_k - J_k * S_k_old);
        // S_k.noalias() = NJ_k_pinv*(xi_k - J_k * S_k_old);
        // S_k += S_k_old;
        
        // Saturation in the Null-Space (SNS) - De Luca : class notes pag. 75+
        // initialize variables to perform SNS
        // weight matrix
        weightW.setIdentity();
        // current task scale factor
        double s = 1.0;
        // largest scale factor found
        double sStar = 0.0;
        // Null-projected jacobian considering also saturations
        MatrixJ Nbar_k = N_k;
        // Null-projected jacobian for projectin saturations in the previous (k-1) task null-spaces
        MatrixJ Ntilde_k;
        // velocity vector to saturate, used in internal loop
        VectorJ qN; // = S_k_old;
        qN.setZero();
        // saturated velocities in the best task scaling scenario
        VectorJ qNstar = qN;
        // weight matrix in the best task scaling scenario
        MatrixJ weightWstar = weightW;
        // null-projected jacobian considering also saturations in the best task scaling scenario
        MatrixJ Nbar_star = Nbar_k;
        // pseudo-inverse of jacobian projected in the null-space with saturations (initially equal to non-saturated pseudo-inverse)
        MatrixXJ JNbar_k = NJ_k;
        // pseudo-inverse of jacobian projected in the null-space with saturations (initially equal to non-saturated pseudo-inverse)
        MatrixJX JNbar_k_pinv = NJ_k_pinv;
        // pseudo-inverse of the auxiliary null-projector, which cosiders also joint saturation
        MatrixJX IWN_k_pinv;
        
        // check limits based on q_in
        bool respecting_limits = false;
        respecting_limits = checkVelocityLimits(q_in.data,S_k);
        
        // enforce limits - if failed, don't go on with tasks
        while(!respecting_limits)
        {
            VectorJ a;
#if (SCALE_PREVIOUS_TASK_CONTRIBUTION > 0)
            {
                pinvDLS((MatrixJ::Identity() - weightW)*N_k,IWN_k_pinv);
                Nbar_k = N_k - IWN_k_pinv*N_k;
                JNbar_k = J_k*Nbar_k;
                pinvDLS(JNbar_k,JNbar_k_pinv);
                VectorJ q_tilde_k = S_k_old + IWN_k_pinv*qN;
                a = JNbar_k_pinv*(xi_k - J_k*q_tilde_k);
            }
#else
            a.noalias() = JNbar_k_pinv*xi_k;
#endif
            // ensure delta_q_dot is zero where it should be multiplying by weightW (needed for numerical reasons)
            a = weightW*a;
            // compute b -- S_k is the variable which gets updated also in the internal cycle
            VectorJ b = S_k - a; 
            
            int worst_joint = -1;
            s = computeMaxScaling(a,b,&worst_joint);
            if(s > sStar)
            {
                sStar = s;
                weightWstar = weightW;
                qNstar = qN;
                Nbar_star = Nbar_k;
            }
            
            weightW(worst_joint,worst_joint) = 0.0;
            // apply in qN the limit resulting from S_k (and a contribution from a if S_k was not out of limits!)
            if(S_k(worst_joint) >= q_dot_ub(worst_joint))
            {
#if DEBUG>2
                std::cout << "saturating ub..." << std::endl;
#endif
                qN(worst_joint) = q_dot_ub(worst_joint) - S_k_old(worst_joint);
            }
            else if(S_k(worst_joint) <= q_dot_lb(worst_joint))
            {
#if DEBUG>2
                std::cout << "saturating lb..." << std::endl;
#endif
                qN(worst_joint) = q_dot_lb(worst_joint) - S_k_old(worst_joint);
            }
            else
            {
#if DEBUG>2
                std::cout << "NOT saturating anything..." << std::endl;
#endif
            }

#if DEBUG>1
            std::cout << "qN = [" << qN.transpose() << "]" << std::endl;
#endif
            
            pinvDLS((MatrixJ::Identity() - weightW)*N_k,IWN_k_pinv);
            Nbar_k = N_k - IWN_k_pinv*N_k;
            JNbar_k = J_k*Nbar_k;
            
            // the task went out of scope of the jacobian (deficient rank): scale it and exit the k-th task cycle
            if( (JNbar_k).colPivHouseholderQr().rank() < xi_k.rows() )
            {
                pinvDLS((MatrixJ::Identity() - weightWstar)*N_k,IWN_k_pinv);
                VectorJ q_tilde_k = S_k_old + IWN_k_pinv*qNstar;
                // enforce limits on q_tilde_k
                enforceWLimits(q_tilde_k);
                Nbar_k = N_k - IWN_k_pinv*N_k;
                JNbar_k = J_k*Nbar_k;
                pinvDLS(JNbar_k, JNbar_k_pinv);
                if(k>0 && sStar <= S_margin)
                {
                    S_k = S_k_old;
                    // as we break, and potentially update the task-space projector, avoid including the last task
                    NJ_k.setZero();
                    NJ_k_pinv.setZero();
                }
                else
                {
#if SCALE_PREVIOUS_TASK_CONTRIBUTION>0
                    // test the idea of scaling also effect of k-1 tasks onto k-th task
                    S_k = q_tilde_k + weightWstar*(JNbar_k_pinv*((sStar-S_margin)*(xi_k - J_k * q_tilde_k)));
#else
                    S_k = q_tilde_k + weightWstar*(JNbar_k_pinv*((sStar-S_margin)*xi_k - J_k * q_tilde_k));
#endif
                }
#if DEBUG>1
                ROS_WARN_STREAM(CLASS_NAMESPACE << __func__ << " : enforcing joint limits at task #" << k+1 << " out of " << task_nr_ << " after performing " << 100*sStar << "% of last task - continuing with next one...");
#endif
                break;
            }
            else
            {
#if DEBUG>1
                std::cout << CLASS_NAMESPACE << __func__ << " : NOT YET exiting (cycle #" << counter++ << "), rank=" << (JNbar_k).colPivHouseholderQr().rank() << " | task #" << k+1 << " out of " << task_nr_ << " | dimension=" << xi_k.rows() << std::endl;
#endif
            }

            // compute an extra step and check again for limits
            pinvDLS(JNbar_k,JNbar_k_pinv);
            VectorJ q_tilde_k = S_k_old + IWN_k_pinv*qN;
            // enforce limits on q_tilde_k
            enforceWLimits(q_tilde_k);
            
            // ensure delta_q_dot is zero where it should be multiplying by weightW (needed for numerical reasons)
            S_k = q_tilde_k + weightW*(JNbar_k_pinv*(xi_k - J_k * q_tilde_k));
            
            respecting_limits = checkVelocityLimits(q_in.data,S_k);
        }
        
        // if it's not last task, compute iterative step
        N_k -= NJ_k_pinv*NJ_k;
    }
    
    // return q_dot
    qdot_out.data = S_k;
    
    // compute max scaling for the full q_dot vector using line-search
    if(model_tolerance_ > 0.0)
    {
        double alpha = computeBestAlphaLineSearch(q_in,v_in,qdot_out,jac_kdl);
        qdot_out.data = alpha*S_k;
    }
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
    
#if DEBUG>1
    std::cout << CLASS_NAMESPACE << __func__ << std::endl;
    std::cout << " : sigma = [" << sigma.transpose() << "]" << std::endl;
    std::cout << " : zero_sv_counter = " << zero_sv_counter << std::endl;
#endif
    
    return zero_sv_counter;
}

bool ChainIkSolverVel_MT_FP_JL::checkVelocityLimits(const VectorJ& q_in, const VectorJ& q_dot_k)
{
    // condition to satisfy:
    // qd_lb <= qd <= qd_ub
    
    int ru,cu,rl,cl;
    double maxu,maxl;
    maxu = (q_dot_k - q_dot_ub).maxCoeff(&ru,&cu);
    maxl = (q_dot_lb - q_dot_k).maxCoeff(&rl,&cl);
    
    to_be_checked_for_limits_.setZero();
    for(int i=0; i<JS_dim; ++i)
    {
        if((q_dot_k(i) - q_dot_ub(i) > QDOT_ZERO) || (q_dot_lb(i) - q_dot_k(i) > QDOT_ZERO))
            to_be_checked_for_limits_(i) = 1.0;
    }
    
#if DEBUG>1
    std::cout << CLASS_NAMESPACE << __func__ << std::endl;
    std::cout << " : to_be_checked_for_limits_=[" << to_be_checked_for_limits_.transpose() << "]" << std::endl;
#endif
#if DEBUG>2
    std::cout << " : q_in=[" << q_in.transpose() << "]" << std::endl;
    std::cout << " : q_dot_k=[" << q_dot_k.transpose() << "]" << std::endl;
    std::cout << " : q_lb=[" << q_lb.transpose() << "]" << std::endl;
    std::cout << " : q_ub=[" << q_ub.transpose() << "]" << std::endl;
    std::cout << " : maxu=" << maxu << " | maxl=" << maxl << " | ru=" << ru << " | rl=" << rl << std::endl;
#endif
    
    // if either one is positive (greater than a small number treated as zero), return false
    if ((maxu > QDOT_ZERO) || (maxl > QDOT_ZERO))
        return false;
    
    // none is positive, return true
    return true;
}

double ChainIkSolverVel_MT_FP_JL::computeMaxScaling(const VectorJ& a, const VectorJ& b, int* r)
{
#if DEBUG>1
    std::cout << CLASS_NAMESPACE << __func__ << std::endl;
    std::cout << " : a        = [" << a.transpose() << "]" << std::endl;
    std::cout << " : b        = [" << b.transpose() << "]" << std::endl;
    std::cout << " : weightW  = [" << weightW.diagonal().transpose() << "]" << std::endl;
#endif
#if DEBUG>2
    std::cout << " : q_dot_lb = [" << q_dot_lb.transpose() << "]" << std::endl;
    std::cout << " : q_dot_ub = [" << q_dot_ub.transpose() << "]" << std::endl;
#endif
    
    VectorJ sMin, sMax;
    // the use of a for cycle is maybe the best thing here...
    for(int i = 0; i<JS_dim; ++i)
    {
        // if(weightW.diagonal()(i) == 0.0 || to_be_checked_for_limits_(i) == 0.0)
        if(weightW.diagonal()(i) == 0.0)
        {
            sMin(i) = -1.0*std::numeric_limits<double>::max();
            sMax(i) = std::numeric_limits<double>::max();
#if DEBUG>2
            if(!((b(i) <= (q_dot_ub(i) + QDOT_ZERO)) && (b(i) >= (q_dot_lb(i) - QDOT_ZERO))))
            {
                std::cout << " - looking at component #" << i << std::endl;
                std::cout << " - (q_dot_ub(i) - b(i)) = " << (q_dot_ub(i) - b(i)) << std::endl;
                std::cout << " - (b(i) - q_dot_lb(i)) = " << (b(i) - q_dot_lb(i)) << std::endl;
            }
#endif
        }
        else
        {
            sMin(i) = (q_dot_lb(i) - b(i)) / a(i);
            sMax(i) = (q_dot_ub(i) - b(i)) / a(i);
        }
        
        bool swapped = false;
        if(sMin(i) > sMax(i))
        {
            swapped = true;
            std::swap(sMin(i),sMax(i));
        }
#if DEBUG>2
        std::cout << " > \t" << sMin(i) << "\t" << sMax(i) << (swapped?"\t>>\tswapping!":"") << std::endl;
#endif
    }
    
    double smin, smax;
    int c,r2,c2;
    smax = sMax.minCoeff(r,&c);
    smin = sMin.maxCoeff(&r2,&c2);
    
#if DEBUG>1
    std::cout << " : sMin     = [" << sMin.transpose() << "]" << std::endl;
    std::cout << " : sMax     = [" << sMax.transpose() << "]" << std::endl;
    std::cout << " : worst_joint=" << *r << std::endl;
#endif
    
    if((smin > smax) || (smax < 0.0) || (smin > 1.0 + S_margin))
    {
#if DEBUG>1
        std::cout << " : returning 0" << std::endl;
#endif
        return 0.0;
    }
    else
    {
        assert(to_be_checked_for_limits_(*r) != 0.0);
#if DEBUG>1
        std::cout << " : returning smax=" << smax << std::endl;
#endif
        return smax;
    }
}

void ChainIkSolverVel_MT_FP_JL::updateVelocityLimits(const VectorJ& q_in)
{
    // this is not the most generic way of doing this...
    // would it make sense to consider also velocity and acceleration limits of the real robot?
    q_dot_lb = q_lb - q_in;
    q_dot_ub = q_ub - q_in;
    
//     // TODO: make this come from the outside: set velocity limits (between iterations - avoid too high jumps)
//     for(int i=0; i<JS_dim; ++i)
//     {
//         if(q_dot_lb(i) < -0.2)
//             q_dot_lb(i) = -0.2;
//         if(q_dot_ub(i) > 0.2)
//             q_dot_ub(i) = 0.2;
//     }
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

bool ChainIkSolverVel_MT_FP_JL::enforceWLimits(VectorJ& q_dot)
{
    // enforce limits due to W in S_k
    bool respecting_limits = true;
    for(int i=0; i<JS_dim; ++i)
    {
        if(weightW(i,i) == 0.0)
        {
            if(q_dot(i) < q_dot_lb(i))
            {
                respecting_limits = false;
                q_dot(i) = q_dot_lb(i);
            }
            else if(q_dot(i) > q_dot_ub(i))
            {
                respecting_limits = false;
                q_dot(i) = q_dot_ub(i);
            }
        }
    }
    return respecting_limits;
}

double ChainIkSolverVel_MT_FP_JL::computeBestAlphaLineSearch(const KDL::JntArray& q, const KDL::Twist& xi, const KDL::JntArray& q_dot, const KDL::Jacobian& jac)
{
    // compute originally desired pose and beginning pose
    KDL::Frame pStar,p;
    fksolver.JntToCart(q,p);
    pStar = KDL::addDelta(p,xi,1.0);
    // compute twist due to q_dot
    VectorT vec_p = jac.data*q_dot.data;
    KDL::Twist delta_p;
    for(int j=0; j<TS_dim; ++j)
        delta_p[j] = vec_p(j);
    if(use_ee_task_)
        delta_p = p.M*delta_p;
    
    double alpha = 1.0;
    double alpha_prev_low = 0.0;
    double alpha_prev_high = 1.0;
    KDL::JntArray q_a(JS_dim);
    KDL::Frame p_alpha,p_alpha_lin;
    KDL::Twist diff_t;
    
    while(alpha_prev_high - alpha_prev_low > 0.125) // i.e. either 1 or 4 steps
    {
        // compute REAL pose at alpha distance from q, in the direction of q_dot
        q_a.data = q.data + (q_dot.data * alpha);
        fksolver.JntToCart(q_a,p_alpha);
        
        // compute LINEAR pose at alpha distance from q, in the direction of q_dot
        p_alpha_lin = KDL::addDelta(p,delta_p,alpha);
        
        // take the difference, compute the norm, check for tolerance reached
        diff_t = KDL::diff(p_alpha,p_alpha_lin);
        
        // consider if we are weighting in end-effector frame
        if(use_ee_task_)
            diff_t = p.M.Inverse()*diff_t;
        
        for(int j=0; j<TS_dim; ++j)
            diff_t[j] *= weightTS(j,j);
        if(Equal(diff_t,Twist::Zero(),model_tolerance_))
        {
            alpha_prev_low = alpha;
        }
        else
        {
            alpha_prev_high = alpha;
        }
        alpha = (alpha_prev_high + alpha_prev_low)/2.0;
    }
    
#if DEBUG>1
    std::cout << CLASS_NAMESPACE << __func__ << " : scaling the task (if alpha is 1.0 > 1 iteration, else always max iterations) > alpha=" << alpha << std::endl;
#endif
#if DEBUG>2
    char y;
    std::cin >> y;
#endif
    
    return alpha;
}
