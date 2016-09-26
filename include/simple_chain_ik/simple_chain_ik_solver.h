#ifndef SIMPLE_CHAIN_IK_SOLVER_H
#define SIMPLE_CHAIN_IK_SOLVER_H
#include <ros/ros.h>
#include <XmlRpcValue.h>
#include <geometry_msgs/Pose.h>
#include <urdf/model.h>
#include <kdl/treefksolverpos_recursive.hpp>

#include "simple_chain_ik/chain_and_solvers.h"

class simple_chain_ik_solver
{
public:
    simple_chain_ik_solver();
    bool get_ik(const std::string& chain, const KDL::Frame& ee, const KDL::JntArray& q_init, KDL::JntArray& q_out, bool publish = false);
    bool get_fk(const std::string& chain, const KDL::JntArray& j, KDL::Frame& ee, bool publish = false);
    bool get_gravity(const std::string& chain, const KDL::JntArray& j, KDL::JntArray& tau, bool publish = false);
    bool get_gravity(const std::string& chain, const KDL::JntArray& j, const std::map< std::string, KDL::Wrench >& w_ext, KDL::JntArray& tau, bool publish = false);
private:
    void initialize_solvers(ChainAndSolvers& container) const; //, const KDL::Tree& robot_kdl, int chain_index) const;
    void parseParameters(XmlRpc::XmlRpcValue& params);
    bool publishConfig(const std::vector< std::string >& joint_names, const KDL::JntArray& q);
    bool normalizePoses(std::vector< geometry_msgs::Pose >& poses);
    bool normalizePose(geometry_msgs::Pose& pose);
private:
    ros::NodeHandle nh;
    std::string robot_urdf;
    urdf::Model urdf_model;
    std::vector<std::string> chain_names_list;
    std::vector<std::string> chain_roots_list;
    std::vector<std::string> chain_ees_list;
    std::map<std::string,KDL::Chain> chains;
    std::map<std::string,std::unique_ptr<ChainAndSolvers>> solvers;
    ros::Publisher joint_state_pub;
    KDL::Vector gravity;
    std::shared_ptr<KDL::TreeFkSolverPos_recursive> tree_fk;
    std::shared_ptr<KDL::Tree> tree;
};

#endif // SIMPLE_CHAIN_IK_SOLVER_H
