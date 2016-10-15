// Copyright (C)  2016 Hamal Marino <hamal dot marino at gmail dot com>
// license: BSD

#ifndef KDL_CHAIN_IKSOLVERVEL_MT_FP_JL_HPP
#define KDL_CHAIN_IKSOLVERVEL_MT_FP_JL_HPP

#include "kdl/chainiksolver.hpp"
#include "kdl/chainjnttojacsolver.hpp"
#include <Eigen/Core>

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
    int setWeightJS(const Eigen::MatrixXd& Mq);
    
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
    int setWeightTS(const Eigen::MatrixXd& weights);
    
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
    void setJointLimits(const Eigen::ArrayXd& lower_bound, const Eigen::ArrayXd& upper_bound);
    
    /**
    * @brief Return the value of eps
    */
    double getEps()const {return eps;};
    
    /**
    * @brief Return the value of lambda
    */
    double getLambda()const {return lambda;};
    
    /**
    * @brief Return the latest return code from the SVD algorithm.
    * 
    * @return 0 if CartToJnt() not yet called, otherwise latest SVD result code.
    */
    int getSVDResult()const {return svdResult;};
    
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
    Eigen::Matrix<double,TS_dim,JS_dim> jac;
    
    // SVD parameters
    /// matrixes to use for the Jacobian SVD: they will change size at each iteration
    Eigen::MatrixXd svdU;
    Eigen::VectorXd svdS;
    Eigen::Matrix<double,JS_dim,JS_dim> svdV;
    /// tolerance for damping singular values
    double eps;
    /// maximum number of iterations for computing SVD
    int maxiter;
    /// damping parameter for the SVD
    double lambda;
    /// latest return value while computing the SVD
    int svdResult;
    
    // weighting-related parameters
    /// flag to say whether it is needed to multiply the Jacobian for the Joint Weight Matrix
    bool is_jac_weighted;
    /// weigth matrix in Task-space (only diagonal elements are considered)
    Eigen::Matrix<double,TS_dim,TS_dim> weightTS;
    /// weight matrix in Joint-space
    Eigen::Matrix<double,JS_dim,JS_dim> weightJS;
    /// number of tasks to perform (based on task-space weighting)
    uint task_nr_;
    /// list of indexes of tasks to perform
    Eigen::Matrix<double,TS_dim,1> task_list_;
    
    // joint limits
    /// joint lower bounds
    Eigen::Matrix<double,JS_dim,1> q_lb;
    /// joint upper bounds
    Eigen::Matrix<double,JS_dim,1> q_ub;
    /// joint velocity lower bounds
    Eigen::Matrix<double,JS_dim,1> q_dot_lb;
    /// joint velocity upper bounds
    Eigen::Matrix<double,JS_dim,1> q_dot_ub;
    
private:
    /**
     * @brief Select the submatrix associated to the rows of @p jac which correspond to rows where @p task_list_ == @p k
     * 
     * @return jac_k
     */
    void selectMatrixRows(const Eigen::Matrix<double,TS_dim,1>& task_list_, uint k, const Eigen::Matrix<double,TS_dim,JS_dim>& jac, Eigen::Matrix<double,-1,JS_dim>& jac_k) const;
    
    /**
     * @brief Check joint limits based on the newly received joint position and the currently computed velocities
     * 
     * @return -1 if all joints are inside the limits, otherwise the index of the "worst" joint (further apart from the limits)
     */
    int checkVelocityLimits(const KDL::JntArray& q_in, const Eigen::Matrix<double,TS_dim,1>& q_dot_k);
    
    /**
     * @brief Compute pseudo-inverse with damped-least-square
     * 
     * @return success/error code
     */
    int pinvDLS(const Eigen::Matrix<double,-1,JS_dim>& NJ_k, Eigen::Matrix<double,JS_dim,-1>& NJ_k_pinv);
    
    
};

/// required forward declaration of template class for it to be instanciated in the library
// template <> class ChainIkSolverVel_MT_FP_JL<6,7>{};

}
#endif // KDL_CHAIN_IKSOLVERVEL_MT_FP_JL_HPP
