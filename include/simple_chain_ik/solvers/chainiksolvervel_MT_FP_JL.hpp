// Copyright (C)  2016 Hamal Marino <hamal dot marino at gmail dot com>
// license: BSD

#ifndef KDL_CHAIN_IKSOLVERVEL_MT_FP_JL_HPP
#define KDL_CHAIN_IKSOLVERVEL_MT_FP_JL_HPP

#include "kdl/chainiksolver.hpp"
#include "kdl/chainjnttojacsolver.hpp"
#include <Eigen/Dense>

namespace KDL
{

/**
    * Implementation of an IK velocity solver algorithm.
    * The (quite cumbersome) name comes from:
    * Multi-Task, Fixed-Priority, with Joint Limits
    * 
    * It uses an iterative approach to solve the tasks, and internally uses
    * a Selectively Damped Least Square (SDLS) algorithm to compute the various
    * components of qDot.
    *
    * For more details, see
    * 1) Siciliano, Slotine (ICAR - 1991)
    * 2) Buss, Kim (2004)
    */

static const int TS_dim = 6;
static const int JS_dim = 7;

typedef Eigen::Matrix<double,TS_dim,TS_dim> MatrixT;
typedef Eigen::Matrix<double,JS_dim,JS_dim> MatrixJ;
typedef Eigen::Matrix<double,-1,-1> MatrixX;
typedef Eigen::Matrix<double,TS_dim,JS_dim> MatrixTJ;
typedef Eigen::Matrix<double,JS_dim,TS_dim> MatrixJT;
typedef Eigen::Matrix<double,TS_dim,-1> MatrixTX;
typedef Eigen::Matrix<double,JS_dim,-1> MatrixJX;
typedef Eigen::Matrix<double,-1,JS_dim> MatrixXJ;
typedef Eigen::Matrix<double,-1,TS_dim> MatrixXT;
typedef Eigen::Matrix<double,TS_dim,1> VectorT;
typedef Eigen::Matrix<int,TS_dim,1> VectorTi;
typedef Eigen::Matrix<double,JS_dim,1> VectorJ;
typedef Eigen::Matrix<double,-1,1> VectorX;

// template<int TS_dim, int JS_dim>
class ChainIkSolverVel_MT_FP_JL : public ChainIkSolverVel
{
public:
    /// solution converged but (pseudo)inverse is singular
    static const int E_CONVERGE_PINV_SINGULAR = +100;
    /// found only a partial solution due to joint limits
    static const int E_SNS_NEEDED = +101;
    
    /**
    * Constructor
    *
    * @param chain the chain to calculate the inverse velocity
    * kinematics for
    * @param eps if a singular value is below this value, its
    * inverse is set to zero, default: 0.00001
    * @param maxiter maximum iterations for the SVD calculation,
    * default: 150
    */
    explicit ChainIkSolverVel_MT_FP_JL(const Chain& chain,double eps=0.00001,int maxiter=150);
    
    /// Destructor
    ~ChainIkSolverVel_MT_FP_JL() {};
    
    /**
    * Find an output joint velocity \a qdot_out, given a starting joint pose
    * \a q_init and a desired cartesian velocity \a v_in
    *
    * @return
    *  E_NOERROR=svd solution converged in maxiter
    *  E_SVD_FAILED=svd solution failed
    *  E_CONVERGE_PINV_SINGULAR=svd solution converged but (pseudo)inverse singular
    *
    * @note if E_CONVERGE_PINV_SINGULAR returned then converged and can
    * continue motion, but have degraded solution
    *
    * @note If E_SVD_FAILED returned, then getSvdResult() returns the error
    * code from the SVD algorithm.
    */
    virtual int CartToJnt(const JntArray& q_in, const Twist& v_in, JntArray& qdot_out);
    
    /**
        * not (yet) implemented.
        */
    virtual int CartToJnt(const JntArray& q_init, const FrameVel& v_in, JntArrayVel& q_out){return -1;};
    
    /**
    * @brief Set the joint space weighting matrix
    *
    * @param weight_js joint space weighting matrix, must be symmetric and positive definite
    *
    * @return succes/error code
    */
    int setWeightJS(const MatrixJ& Mq);
    
    /**
    * @brief Set the task space weighting
    *
    * The multi-task priority comes from the diagonal of this matrix: equal weights correspond to equal priority tasks,
    * higher weights correspond to higher priority.
    * 
    * @param weight_ts task space weighting matrix (only diagonal vector considered), all weights have to be positive.
    *
    * @return succes/error code
    */
    int setWeightTS(const MatrixT& weights);
    
    /**
    * @brief Set lambda parameter used in damping the eigenvalues of the SVD
    */
    void setLambda(const double lambda_in) {lambda = lambda_in;};

    /**
    * @brief Set eps, the tolerance for damping the SVD
    */
    void setEps(const double eps_in) {eps = eps_in;};

    /**
    * @brief Set the maximum number of iterations for computing the SVD
    */
    void setMaxIter(const int maxiter_in) {maxiter = maxiter_in;};
    
    /**
     * @brief Set bounds for the q's of this chain
     * 
     * @param lower_bound Minimum for q
     * @param upper_bound Maximum for q
     */
    void setJointLimits(const VectorJ& lower_bound, const VectorJ& upper_bound);
    
    /**
    * @brief Return the value of eps
    */
    double getEps()const {return eps;};
    
    /**
    * @brief Return the value of lambda
    */
    double getLambda()const {return lambda;};
    
    /// @copydoc KDL::SolverI::strError()
    virtual const char* strError(const int error) const;
    
private:
    /// internal copy of the chain we are solving for
    const Chain chain;
    /// number of joints
    unsigned int nj;
    /// task space dimension
    unsigned int ts_dim;
    
    // Jacobian parameters
    /// Joint to Jacobian solver
    ChainJntToJacSolver jnt2jac;
    /// KDL jacobian
    Jacobian jac_kdl;
    /// Eigen jacobian
    MatrixTJ jac;
    
    // SVD parameters
    /// tolerance for damping singular values
    double eps;
    /// maximum number of iterations for computing SVD
    int maxiter;
    /// damping parameter for the SVD
    double lambda;
    
    // variables for computing IK
    /// solution at step k-1
    VectorJ S_k_old;
    /// solution at step k
    VectorJ S_k;
    /// null-space projector at step k
    MatrixJ N_k;
    /// complete task specification
    VectorT xi;
    
    // weighting-related parameters
    /// flag to say whether it is needed to multiply the Jacobian for the Joint Weight Matrix
    bool is_jac_weighted;
    /// weigth matrix in Task-space (only diagonal elements are considered)
    MatrixT weightTS;
    /// weight matrix used for limiting joint in null-space (Saturation in the Null-Space - SNS)
    MatrixJ weightW;
    /// weight matrix in Joint-space
    MatrixJ weightJS;
    /// number of tasks to perform (based on task-space weighting)
    uint task_nr_;
    /// list of indexes of tasks to perform
    VectorTi task_list_;
    /// list of joints which are out of limits
    VectorJ to_be_checked_for_limits_;
    
    // joint limits
    /// joint lower bounds
    VectorJ q_lb;
    /// joint upper bounds
    VectorJ q_ub;
    /// joint velocity lower bounds
    VectorJ q_dot_lb;
    /// joint velocity upper bounds
    VectorJ q_dot_ub;
    
private:
    /**
     * @brief Select the submatrix associated to the rows of @p jac which correspond to rows where @p task_list_ == @p k
     * 
     * @return jac_k
     */
    void selectMatrixRows(const VectorTi& task_list_, uint k, const MatrixTJ& jac, MatrixXJ& jac_k) const;
    void selectMatrixRows(const VectorTi& task_list_, uint k, const VectorT& xi, VectorX& xi_k) const;
    
    /**
     * @brief Check joint limits based on the newly received joint position and the currently computed velocities
     * 
     * @return true if all joints are inside the limits, false otherwise
     */
    bool checkVelocityLimits(const KDL::VectorJ& q_in, const KDL::VectorJ& q_dot_k);
    
    /**
     * @brief Compute pseudo-inverse with damped-least-square
     * 
     * @return success/error code
     */
    int pinvDLS(const MatrixXJ& NJ_k, MatrixJX& NJ_k_pinv);
    
    /**
     * @brief Compute maximum scaling parameter, given limits and currently computed solution
     * 
     * @param a velocity contribution of the current task, with limited joints
     * @param b difference between total task velocity and @p a
     * @param r index of the most critical joint, which makes sense only if the return is non-zero
     * 
     * @return maximum scaling factor
     */
    double computeMaxScaling(const VectorJ& a, const VectorJ& b, int* r);
    
    /**
     * @brief Update the internal values of joint velocity limits, based on position threshold
     */
    void updateVelocityLimits(const VectorJ& q_in);
    
    /**
     * @brief Enforce the limits on the vector @p q_dot dictated by weightW: if an element of the diagonal of weightW is zero, the corresponding value must be saturated
     * 
     * @param q_dot the joint velocity vector to check for limits
     * 
     * @return true if all joints were already inside limits, false otherwise
     */
    bool enforceWLimits(KDL::VectorJ& q_dot);
    
};

/// required forward declaration of template class for it to be instanciated in the library
// template <> class ChainIkSolverVel_MT_FP_JL<6,7>{};

}
#endif // KDL_CHAIN_IKSOLVERVEL_MT_FP_JL_HPP
