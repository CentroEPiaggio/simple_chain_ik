#ifndef QB_LEGS_IK_H
#define QB_LEGS_IK_H
#include <ros/ros.h>
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
    KDL::ChainIdSolver_RNE* idsolver=0;
    std::vector<std::string> joint_names;
    int index;
    KDL::JntArray q_min, q_max;
};


class qb_legs_ik
{
public:
    qb_legs_ik();
    bool get_ik(const std::string& chain, const KDL::Frame& ee, const KDL::JntArray& q_init, KDL::JntArray& q_out, bool publish = false);
    bool get_fk(const std::string& chain, const KDL::JntArray& j, KDL::Frame& ee, bool publish = false);
    bool get_gravity(const std::string& chain, const KDL::JntArray& j, KDL::JntArray& tau, bool publish = false);
private:
    void initialize_solvers(chain_and_solvers* container) const;
    void parseParameters(XmlRpc::XmlRpcValue& params);
    bool publishConfig(const std::vector< std::string >& joint_names, const KDL::JntArray& q);
    bool normalizePoses(std::vector< geometry_msgs::Pose >& poses);
    bool normalizePose(geometry_msgs::Pose& pose);
private:
    ros::NodeHandle nh;
    std::string robot_urdf;
    urdf::Model urdf_model;
    KDL::Tree robot_kdl;
    std::vector<std::string> chain_names_list;
    std::vector<std::string> chain_roots_list;
    std::vector<std::string> chain_ees_list;
    std::map<std::string,KDL::Chain> chains;
    std::map<std::string,KDL::Chain> chains_reverse;
    std::map<std::string,chain_and_solvers> solvers;
    ros::Publisher joint_state_pub;
};

#endif // QB_LEGS_IK_H
