#ifndef CHAIN_AND_SOLVERS_H_
#define CHAIN_AND_SOLVERS_H_

#include <kdl/frames.hpp>
#include <kdl/chain.hpp>
#include <kdl/jntarray.hpp>
// #include <kdl/chainiksolverpos_nr_jl.hpp>
#include <simple_chain_ik/solvers/chainiksolverpos_relaxed.hpp>
// #include <kdl/chainiksolvervel_wdls.hpp>
#include <my_chainiksolvervel_wdls.hpp>
#include <kdl/chainfksolverpos_recursive.hpp>
#include <kdl/chainidsolver_recursive_newton_euler.hpp>
#include <kdl/tree.hpp>
#include <kdl/treefksolver.hpp>
#include <memory>
#include <eigen3/Eigen/Dense>

namespace ChainAndSolversTypes
{
    typedef KDL::ChainIkSolverPos_relaxed ChainIkSolverPos;
    typedef KDL::ChainFkSolverPos_recursive ChainFkSolverPos;
    typedef KDL::my_ChainIkSolverVel_wdls ChainIkSolverVel;
    // typedef KDL::ChainIkSolverVel_wdls ChainIkSolverVel;
    typedef KDL::ChainIdSolver_RNE ChainIdSolver;
}

class ChainAndSolvers
{
public:
    /**
     * @brief Full constructor - all other will construct members these parameters internally
     * 
     * @param tree_ A shared pointer to a KDL tree which will be used to get the chain
     * @param tree_fk_ A shared pointer to a KDL tree FK position solver
     * @param chain_root_ The root name of the chain to be used
     * @param chain_tip_ The tip of the chain to be used
     * @param tree_root_gravity_ The gravity expressed in the root of the tree: will be later converted in the root of the chain
     */
    ChainAndSolvers(const std::shared_ptr< KDL::Tree >& tree_, const std::shared_ptr< KDL::TreeFkSolverPos >& tree_fk_, const std::string& chain_root_, const std::string& chain_tip_, const KDL::Vector& tree_root_gravity_);
    
    /// @brief Constructor
    ChainAndSolvers(const std::shared_ptr< KDL::Tree >& tree_, const std::string& chain_root_, const std::string& chain_tip_, const KDL::Vector& tree_root_gravity_);
    
    /// @brief Constructor
    ChainAndSolvers(const KDL::Chain& chain_, const KDL::Vector& tree_root_gravity_);
    
    /**
     * @brief Set all parameters needed by the solvers
     * 
     * @param q_min_ vector of min values for the joints
     * @param q_max_ vector of max values for the joints
     * @param pos_IK_max_iter_ maximum number of iterations to be used in the position IK solver
     * @param pos_IK_eps_ tolerance to be used in the position IK solver
     * @param vel_IK_max_iter_ maximum number of iterations to be used in the velocity IK solver
     * @param vel_IK_eps_ tolerance to be used in the velocity IK solver (for computing the pseudo-inverse)
     * @param vel_IK_lambda_ damping parameter to be used in the velocity IK solver (for computing the pseudo-inverse)
     */
    bool setSolverParameters(const KDL::JntArray& q_min_, const KDL::JntArray& q_max_, uint pos_IK_max_iter_, double pos_IK_eps_, uint vel_IK_max_iter_, double vel_IK_eps_, double vel_IK_lambda_);
    
    /**
     * @brief Initialize the solvers with the currently set parameters
     * 
     * @return false if something went wrong
     */
    bool initSolvers();
    
    /**
     * @brief Get the names of the joints this class is managing at the moment
     */
    const std::vector<std::string>& jointNames();
    
    /**
     * @brief Function to get a vector of joint torques due to gravity; such vector is rectified using tau_multiplier's in order to have inverted torques if the chain has been obtained reverting (even part of) the tree from tip to root
     * 
     * @param j joint values for which we are computing gravity
     * @param w_ext map of external wrenches applied to some links
     * @param tau computed joint torques
     * 
     * @return false if something went wrong
     */
    bool getGravity(const KDL::JntArray& j, const std::map< std::string, KDL::Wrench >& w_ext, KDL::JntArray& tau);
    
    /// Direct access to internal members; can be used for calling member function(s), but it advised not to keep a reference, as the solver may be reinstanciated at run-time.
    std::unique_ptr< ChainAndSolversTypes::ChainIkSolverPos >& getIKSolver();
    std::unique_ptr< ChainAndSolversTypes::ChainFkSolverPos >& getFKSolver();
    std::unique_ptr< ChainAndSolversTypes::ChainIkSolverVel >& getIKVelSolver();
    
    /**
     * @brief Get a random valid joint array (between minimum and maximum)
     * 
     * @return A joint array, granted to have values between minimum and maximum allowed, as passed in @p setSolverParameters
     */
    KDL::JntArray getValidRandomJoints();
    
    /**
     * @brief Changes the constant frame at the tip of the chain used for computing IK and FK. When this function is called, all solvers get erased and created as new, so any property which was set on the solver directly, and not through this class interface, needs to be set again.
     * 
     * @param ee_tip_ a frame for the new tip, to be set w.r.t. the chain's end-effector
     */
    void changeTip(const KDL::Frame& ee_tip_);
    
    /**
     * @brief Change the weight of the task for computing inverse kinematics, i.e. how each direction is weighted in terms of the error: default to Identity. A weight of zero on a direction means it is not considered (the error can be very large), while a higher value means the task is weighted more along that direction.
     * 
     * @param Wx task space weights: will be used as diagonal of a symetric matrix, (semi-)positive definite
     * 
     * @return true on success
     * @see KDL::ChainIkSolverVel_wdls::setWeightTS for further information
     * @see KDL::ChainIkSolverPos_relaxed::setTaskWeight for further information
     */
    bool changeIkTaskWeigth(const Eigen::Matrix<double,6,1>& Wx_);
    
private:
    
    /**
     * @brief Actions needed when the Chain gets changed
     * NOTE: this can support chaining root/tip of the chain
     * 
     * @return false if something went wrong
     */
    bool onChainChanged();
    
    /**
     * @brief compute the gravity in a frame local to the chain root, starting from gravity in a frame local to the tree root
     * 
     * @param gravity The gravity expressed in the base frame of the tree
     * 
     * @return false if something went wrong
     */
    bool computeLocalGravity(const KDL::Vector& tree_root_gravity_);
    
    /**
     * @brief Compute the vector of tau multipliers used for multiplying the torques got out of the inverse dynamics solver, needed when getting a chain out of a tree which goes (even partly) from tip to root of the original tree
     * 
     * @return false if something went wrong
     */
    bool computeTauMultipliers();
    
private:
    /// tree and FK solver got from the constructor
    std::shared_ptr<KDL::Tree> tree;
    std::shared_ptr<KDL::TreeFkSolverPos> tree_fk;
    /// original chain
    KDL::Chain chain_orig;
    /// chain which may contain an extra segment, needed for referencing a constant frame which is not the end-effector
    KDL::Chain chain;
    
    std::unique_ptr<ChainAndSolversTypes::ChainFkSolverPos> fksolver; /// FK position solver
    std::unique_ptr<ChainAndSolversTypes::ChainIkSolverPos> iksolver; /// IK position solver
    std::unique_ptr<ChainAndSolversTypes::ChainIkSolverVel> ikvelsolver; /// IK velocity solver
    std::unique_ptr<ChainAndSolversTypes::ChainIdSolver> idsolver; /// ID solver
    std::vector<std::string> joint_names;
    /// joint limits
    KDL::JntArray q_min, q_max;
    /// vector used for multiplying the torques got out of the inverse dynamics solver, needed when getting a chain out of a tree which goes (even partly) from tip to root of the original tree
    std::vector<double> tau_multiplier;
    
    /// stores whether the function @p initSolvers has been called already; will make other calls fail otherwise
    bool initialized;
    
    /// maximum number of iteration to be used in computing IK (used by IK position solver)
    uint pos_IK_max_iter;
    /// tolerance for computing IK (used by IK position solver)
    double pos_IK_eps;
    /// maximum number of iteration to be used in computing IK (used by IK velocity solver)
    uint vel_IK_max_iter;
    /// tolerance for computing pseudo-inverse (used by IK velocity solver)
    double vel_IK_eps;
    /// damping parameter in pseudo-inverse (used by IK velocity solver)
    double vel_IK_lambda;
    /// names of root and tip of the chain to be used
    std::string chain_root, chain_tip;
    /// local gravity expressed in the root of the chain
    KDL::Vector gravity;
    /// gravity expressed in the root of the tree
    const KDL::Vector tree_root_gravity;
    /// frame containing the tip to use for IK/FK computation, w.r.t. the end-effector
    KDL::Frame ee_tip;
    /// matrix to be used for weighting the task for IKvel computation; has to be simmetric and (semi-)positive definite
    Eigen::Matrix<double,6,6> Mx;
    /// vector to be used for weighting the task for IKpos computation
    Eigen::Matrix<double,6,1> Wx;
};

#endif // CHAIN_AND_SOLVERS_H_
