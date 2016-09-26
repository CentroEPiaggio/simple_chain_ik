#include "simple_chain_ik/simple_chain_ik_solver.h"
#include <kdl_conversions/kdl_msg.h>
#include <math.h>
#include <algorithm>    // std::min_element, std::max_element
#include <kdl_parser/kdl_parser.hpp>
#include <kdl/frames_io.hpp>
#include <kdl/kinfam_io.hpp>
#include <dual_manipulation_shared/parsing_utils.h>
#include <sensor_msgs/JointState.h>

#define CLASS_NAMESPACE "simple_chain_ik_solver::"
#define DEBUG 1 // if 1, 2, ... print some more information
#define MAX_ITER 100
#define EPS 5e-4

double eps = EPS;

simple_chain_ik_solver::simple_chain_ik_solver() : gravity(0.0,0.0,-9.81)
{
    std::string global_name, relative_name, default_param;
    if (nh.getParam("/robot_description", global_name)) //The checks here are redundant, but they may provide more debug info
    {
        robot_urdf=global_name;
    }
    else if (nh.getParam("robot_description", relative_name))
    {
        robot_urdf=relative_name;
    }
    else 
    {
        ROS_ERROR_STREAM(CLASS_NAMESPACE << "constructor : cannot find robot_description");
        abort();
    }
    if (!urdf_model.initParam("robot_description"))
    {
        ROS_ERROR_STREAM(CLASS_NAMESPACE << "constructor : cannot load robot_description");
        abort();
    }
    
    KDL::Tree robot_kdl;
    if (!kdl_parser::treeFromUrdfModel(urdf_model, robot_kdl))
    {
        ROS_ERROR_STREAM("Failed to construct kdl tree");
        abort();
    }
    tree = std::make_shared<KDL::Tree>(robot_kdl);
    
    tree_fk.reset(new KDL::TreeFkSolverPos_recursive(*tree));
    
    // managing external parameters
    XmlRpc::XmlRpcValue ext_params;
    if (nh.getParam("simple_chain_ik_params", ext_params))
        parseParameters(ext_params);
    else
    {
        ROS_ERROR_STREAM(CLASS_NAMESPACE << "constructor : cannot find external parameters, did you load the YAML file?");
        abort();
    }
    
    joint_state_pub = nh.advertise<sensor_msgs::JointState>("simple_chain_ik/joint_states",10);
    
    KDL::Chain temp;
    #if DEBUG
    std::string robot_root = tree->getRootSegment()->first;
    std::cout << "root: " << robot_root << std::endl;
    #endif
    
    for (int i=0; i<chain_ees_list.size(); i++)
    {
        std::string end_effector = chain_ees_list.at(i);
        std::string root = chain_roots_list.at(i);
        tree->getChain(root,end_effector,temp);
        
        chains[chain_names_list.at(i)] = temp;
        
        #if DEBUG>0
        std::cout << "chains[" << chain_names_list.at(i) << "].getNrOfJoints(): " << chains[chain_names_list.at(i)].getNrOfJoints() << std::endl;
        #endif
        #if DEBUG>1
        std::cout << "chains[" << chain_names_list.at(i) << "].segments: | ";
        for(auto segs:chains[chain_names_list.at(i)].segments)
            std::cout << segs.getName() << " | ";
        std::cout << std::endl;
        std::cout << "chains[" << chain_names_list.at(i) << "].segments masses: | ";
        for(auto segs:chains[chain_names_list.at(i)].segments)
            std::cout << segs.getInertia().getMass() << " | ";
        std::cout << std::endl;
        #endif
        
        solvers[chain_names_list.at(i)].reset(new ChainAndSolvers(tree,tree_fk,chain_roots_list.at(i),chain_ees_list.at(i)));
        
        initialize_solvers(*solvers[chain_names_list.at(i)]);
    }
    
    // check which robot I am using
    std::string robot_name = urdf_model.getName();
    std::cout << CLASS_NAMESPACE << __func__ << " : I am using \"" << robot_name << "\" robot!" << std::endl;
}

void simple_chain_ik_solver::parseParameters(XmlRpc::XmlRpcValue& params)
{
    ROS_ASSERT(params.getType() == XmlRpc::XmlRpcValue::TypeStruct);
    parseSingleParameter(params,chain_names_list,"chain_names",1);
    parseSingleParameter(params,chain_roots_list,"chain_roots",1);
    parseSingleParameter(params,chain_ees_list,"chain_ees",1);
    parseSingleParameter(params,eps,"eps");
    std::vector<double> base_gravity({gravity(0),gravity(1),gravity(2)});
    parseSingleParameter(params,base_gravity,"base_gravity",3);
    for(int i=0; i<3; i++) gravity(i) = base_gravity[i];
}

void simple_chain_ik_solver::initialize_solvers(ChainAndSolvers& container) const
{
    KDL::JntArray q_min, q_max;
    q_min.resize(container.jointNames().size());
    q_max.resize(container.jointNames().size());
    int j=0;
    for (auto& joint_name:container.jointNames())
    {
        if(urdf_model.joints_.at(joint_name)->safety)
        {
            q_max(j)=urdf_model.joints_.at(joint_name)->safety->soft_upper_limit;
            q_min(j)=urdf_model.joints_.at(joint_name)->safety->soft_lower_limit;
        }
        else
        {
            q_max(j)=urdf_model.joints_.at(joint_name)->limits->upper;
            q_min(j)=urdf_model.joints_.at(joint_name)->limits->lower;
        }
        j++;
    }
    
    if(!container.setSolverParameters(q_min,q_max,gravity,MAX_ITER,eps,150,1e-5,1e-5) || !container.initSolvers())
    {
        std::cout << CLASS_NAMESPACE << __func__ << " : unable to initialize the solvers! Returning..." << std::endl;
        abort();
    }
}

bool simple_chain_ik_solver::publishConfig(const std::vector< std::string >& joint_names, const KDL::JntArray& q)
{
    sensor_msgs::JointState js_msg;
    js_msg.name = joint_names;
    js_msg.header.stamp = ros::Time::now();
    js_msg.position.clear();
    for(int i=0; i<js_msg.name.size(); i++)
    {
        js_msg.position.push_back(q(i));
    }
    joint_state_pub.publish(js_msg);
    
    return true;
}

bool simple_chain_ik_solver::normalizePose(geometry_msgs::Pose& pose)
{
    geometry_msgs::Quaternion& q(pose.orientation);
    double q_norm = std::sqrt(q.x*q.x+q.y*q.y+q.z*q.z+q.w*q.w);
    q.x = q.x/q_norm;
    q.y = q.y/q_norm;
    q.z = q.z/q_norm;
    q.w = q.w/q_norm;
    
    bool ok = (q_norm < 1.01 && q_norm > 0.99);
    
    if(!ok)
        std::cout << "Pose not properly normalized, quaternion norm was " << q_norm << std::endl;
    
    return ok;
}

bool simple_chain_ik_solver::normalizePoses(std::vector< geometry_msgs::Pose >& poses)
{
    bool ok = true;
    for (auto& p:poses)
        ok = ok & normalizePose(p);
    
    return ok;
}

bool simple_chain_ik_solver::get_ik(const std::string& chain, const KDL::Frame& ee, const KDL::JntArray& q_init, KDL::JntArray& q_out, bool publish)
{
    int res = solvers.at(chain)->getIKSolver()->CartToJnt(q_init,ee,q_out);
    
    if(publish)
    {
        publishConfig(solvers.at(chain)->jointNames(),q_out);
    }
    
    return (res>=0);
}

bool simple_chain_ik_solver::get_fk(const std::string& chain, const KDL::JntArray& j, KDL::Frame& ee, bool publish)
{
    int res = solvers.at(chain)->getFKSolver()->JntToCart(j,ee);
    
    if(publish)
    {
        publishConfig(solvers.at(chain)->jointNames(),j);
    }
    
    return (res>=0);
}

bool simple_chain_ik_solver::get_gravity(const std::string& chain, const KDL::JntArray& j, KDL::JntArray& tau, bool publish)
{
    std::map<std::string,KDL::Wrench> w_ext;
    
    return get_gravity(chain,j,w_ext,tau,publish);
}

bool simple_chain_ik_solver::get_gravity(const std::string& chain, const KDL::JntArray& j, const std::map<std::string,KDL::Wrench>& w_ext, KDL::JntArray& tau, bool publish)
{
    if (!solvers.at(chain)->getGravity(j,w_ext,tau))
        return false;
    
    if(publish)
    {
        publishConfig(solvers.at(chain)->jointNames(),j);
    }
    
    return true;
}
