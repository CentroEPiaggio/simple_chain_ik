#ifndef QB_LEGS_IK_H
#define QB_LEGS_IK_H
#include <kdl/frames.hpp>
#include <kdl/chain.hpp>
#include <kdl/jntarray.hpp>
#include <kdl/chainiksolvervel_pinv.hpp>
#include <kdl/chainiksolverpos_nr_jl.hpp>
#include <kdl/chainfksolverpos_recursive.hpp>
#include <kdl/tree.hpp>
#include <XmlRpcValue.h>
#include <geometry_msgs/Pose.h>
#include <urdf/model.h>
#include <kdl/chainidsolver_recursive_newton_euler.hpp>

#define BIG_NR 10000

class chain_and_solvers
{
public:
    KDL::Chain chain;
    KDL::ChainFkSolverPos_recursive* fksolver=0;
    KDL::ChainIkSolverPos_NR_JL* iksolver=0;
    KDL::ChainIkSolverVel_pinv* ikvelsolver=0;
    std::vector<std::string> joint_names;
    int index;
    KDL::JntArray q_min, q_max;
};


class qb_legs_ik
{
public:
    qb_legs_ik();
    bool normalizePose(geometry_msgs::Pose& pose);
private:
    bool check_ik(std::string ee_name, KDL::Frame World_EE) const;
    void initialize_solvers(chain_and_solvers* container) const;
//     void parseParameters(XmlRpc::XmlRpcValue& params);
    bool publishConfig(const std::vector<std::string>& joint_names, const KDL::JntArray& q) const;
    bool normalizePoses(std::vector< geometry_msgs::Pose >& poses);
private:
    std::map<int,KDL::Frame> fine_tuning;
    std::vector<KDL::Rotation> sphere_sampling;
//     mutable dual_manipulation::ik_control::ikCheckCapability *ik_check_capability;
    mutable chain_and_solvers double_arm_solver;
    std::string robot_urdf;
    urdf::Model urdf_model;
    KDL::Tree robot_kdl;
    mutable std::default_random_engine generator;
    mutable std::uniform_real_distribution<double> distribution;
    // managing external parameters
    XmlRpc::XmlRpcValue ik_control_params;
    std::vector<std::string> chain_names_list;
    std::map<std::string,KDL::Chain> chains;
    std::map<std::string,KDL::Chain> chains_reverse;
    std::map<std::string,KDL::ChainIdSolver_RNE> chain_id_solvers;
};

#endif // QB_LEGS_IK_H
