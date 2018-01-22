// Copyright  (C)  2016 Hamal Marino <hamal dot marino at gmail dot com>
// Copyright  (C)  2007-2008  Ruben Smits <ruben dot smits at mech dot kuleuven dot be>
// Copyright  (C)  2008  Mikael Mayer
// Copyright  (C)  2008  Julia Jesse

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

// This class is based on KDL::ChainIkSolverPos_NR_JL, extended with task weighting

#ifndef CHAINIKSOLVERPOS_RELAXED_HPP
#define CHAINIKSOLVERPOS_RELAXED_HPP

#include <kdl/chainiksolver.hpp>
#include <kdl/chainfksolver.hpp>
#include <kdl/solveri.hpp>
#include <Eigen/Core>

inline void skew_symmetric(KDL::Vector &v_, Eigen::Matrix<double,3,3> &skew_mat_)
{
    skew_mat_ = Eigen::Matrix<double,3,3>::Zero();
    
    skew_mat_(0,1) = -v_(2);
    skew_mat_(0,2) =  v_(1);
    skew_mat_(1,0) =  v_(2);
    skew_mat_(1,2) = -v_(0);
    skew_mat_(2,0) = -v_(1);
    skew_mat_(2,1) =  v_(0);
}

namespace KDL
{

/**
    * @brief Implementation of a position IK solver based on Newton-Raphson and using relaxation coefficients to weight the task space
    */
class ChainIkSolverPos_relaxed : public ChainIkSolverPos
{
public:
    
    static const int E_IKSOLVERVEL_FAILED = -100; //! Child IK solver vel failed
    static const int E_FKSOLVERPOS_FAILED = -101; //! Child FK solver failed
    
    // ATTENTION remove these when KDL version >= 1.4.0 will be available
    static const int E_SIZE_MISMATCH = -4;
    static const int E_MAX_ITERATIONS_EXCEEDED = -5;
    
    /**
     * @brief Constructor; needs the chain, a fwd position and velocity IK solvers.
     *
     * @param chain the chain to calculate the inverse position for
     * @param q_min the minimum joint positions
     * @param q_max the maximum joint positions
     * @param fksolver a forward position kinematics solver
     * @param iksolver an inverse velocity kinematics solver
     * @param maxiter the maximum Newton-Raphson iterations,
     * default: 100
     * @param eps the tolearnce to use for the position, used to end the
     * iterations (default: 1e-6)
     */
    ChainIkSolverPos_relaxed(const Chain& chain, const JntArray& q_min, const JntArray& q_max, ChainFkSolverPos& fksolver, ChainIkSolverVel& iksolver, unsigned int maxiter=100, double eps=1e-6);
    
    /**
     * @brief Constructor without joint limits: uses numerical limits to initialize them
     */
    ChainIkSolverPos_relaxed(const Chain& chain, ChainFkSolverPos& fksolver, ChainIkSolverVel& iksolver, unsigned int maxiter=100, double eps=1e-6);
    
    ~ChainIkSolverPos_relaxed();
    
    /**
     * @brief Computes angular velocity from the rotation matrices for two frames
     * 
     * @param frame_1 first frame
     * @param frame_2 second frame
     * @param ang_vel computed angular velocity
     * 
     * @return false if something goes wrong during the procedure, else return true
     */
    void getAngularVelFromPoses(const Frame& frame_1, const Frame& frame_2, Vector& ang_vel);
    
    /**
     * @brief Compute joint values for the input pose, starting from the initial guess.
     * 
     * @param q_init joints initial guess
     * @param p_in pose for the end-effector (chain-tip)
     * @param q_out resulting joints
     * 
     * @return E_MAX_ITERATIONS_EXCEEDED if the maximum number of iterations was exceeded before a result was found
     *         E_NOT_UP_TO_DATE if the internal data is not up to date with the chain
     *         E_SIZE_MISMATCH if the size of the input/output data does not match the chain.
     */
    virtual int CartToJnt(const JntArray& q_init, const Frame& p_in, JntArray& q_out);
    
    /**
     * @brief Set the joint limits
     * 
     * @param q_min lower limits
     * @param q_max upper limits
     * 
     * @return E_SIZE_MISMATCH if input sizes do not match the chain
     */
    int setJointLimits(const JntArray& q_min, const JntArray& q_max);
    
    /**
     * @brief Set the task weight, a 6x1 vector which will be used as square roots of actual weighting coefficients
     * 
     * @param W the weight vector, must be 6x1 (will multiply the twist when performing IK)
     */
    void setTaskWeight(const Eigen::Matrix<double,6,1>& W_);
    
    /// @copydoc KDL::SolverI::strError()
    const char* strError(const int error) const;
    
    /**
     * @brief Decide whether to use the task weight as expressed in world-frame (default) or End-effector frame.
     * 
     * @param s false to use tolerances in world frame, true to use tolerances in end-effector frame
     * 
     * @return The old value for this flag.
     */
    bool useWeigthEndEffector(bool s = false) {bool tmp(use_ee_task_); use_ee_task_ = s; return tmp;}

    
private:
    const Chain chain;
    /// number of joints
    unsigned int nj;
    /// joint lower limits
    JntArray q_min;
    /// joint upper limits
    JntArray q_max;
    /// internal IK velocity solver
    ChainIkSolverVel& iksolver;
    /// internal FK position solver
    ChainFkSolverPos& fksolver;
    JntArray delta_q;
    unsigned int maxiter;
    double eps;
    /// the task weighting vector
    Eigen::Matrix<double,6,1> W;
    
    /// variables used in internal computation
    Frame f;
    Twist delta_twist;
    
    // flag to say whether using end-effector or world-frame task - for this to work properly, the fksolver needs to be instanciated.
    bool use_ee_task_;
};

}

#endif // CHAINIKSOLVERPOS_RELAXED_HPP
