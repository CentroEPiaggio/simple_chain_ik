#include <simple_chain_ik/chain_and_solvers.h>
#include <ros/ros.h>
#include <urdf/model.h>
#include <kdl_parser/kdl_parser.hpp>
#include <kdl/treefksolverpos_recursive.hpp>

#define CLASS_NAMESPACE "chain_and_solvers_test::"

int main(int argc, char** argv)
{
    while(!ros::isInitialized())
    {
        ros::init(argc,argv,"simple_chain_ik_node");
    }
    ros::NodeHandle nh;
    ros::AsyncSpinner aspin(1);
    aspin.start();
    
    urdf::Model urdf_model;
    if (!urdf_model.initParam("robot_description"))
    {
        ROS_ERROR_STREAM(CLASS_NAMESPACE << " : cannot load robot_description");
        return -1;
    }
    if(urdf_model.getName() != "vito")
    {
        ROS_ERROR_STREAM(CLASS_NAMESPACE << " : I was expecting 'vito' robot, got instead '" << urdf_model.getName() << "'...");
        return -2;
    }
    
    KDL::Tree robot_kdl;
    if (!kdl_parser::treeFromUrdfModel(urdf_model, robot_kdl))
    {
        ROS_ERROR_STREAM(CLASS_NAMESPACE << " : Failed to construct kdl tree");
        return -3;
    }
    std::shared_ptr<KDL::Tree> tree = std::make_shared<KDL::Tree>(robot_kdl);
    
    std::shared_ptr<KDL::TreeFkSolverPos_recursive> tree_fk;
    tree_fk.reset(new KDL::TreeFkSolverPos_recursive(*tree));
    
    ChainAndSolvers rh_solver(tree,tree_fk,"vito_anchor","right_hand_palm_link",KDL::Vector(0.0,0.0,-9.81));
    
    std::cout << "chain_and_solvers_test running!!!" << std::endl;
    
    ros::Rate rate(0.1);
    while(ros::ok())
    {
        rate.sleep();
    }
    return 0;
}
