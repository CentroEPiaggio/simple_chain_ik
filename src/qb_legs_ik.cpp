#include "qb_legs/qb_legs_ik.h"
#include <vector>
#include <kdl_conversions/kdl_msg.h>
#include <math.h>
#include <algorithm>    // std::min_element, std::max_element
#include <std_msgs/String.h>
#include <kdl_parser/kdl_parser.hpp>
#include <random>
#include <kdl/frames_io.hpp>
#include <moveit/robot_model/joint_model_group.h>
#include <kdl/kinfam_io.hpp>
#include <dual_manipulation_shared/parsing_utils.h>
#include <sensor_msgs/JointState.h>

#define CLASS_NAMESPACE "qb_legs_ik::"
#define DEBUG 1 // if 1, 2, ... print some more information
#define MAX_ITER 100
#define EPS 5e-4

double eps = EPS;

qb_legs_ik::qb_legs_ik() : gravity(0.0,0.0,-9.81)
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
    
    if (!kdl_parser::treeFromUrdfModel(urdf_model, robot_kdl))
    {
        ROS_ERROR_STREAM("Failed to construct kdl tree");
        abort();
    }
    
    // managing external parameters
    XmlRpc::XmlRpcValue ext_params;
    if (nh.getParam("qb_legs_params", ext_params))
        parseParameters(ext_params);
    else
    {
        ROS_ERROR_STREAM(CLASS_NAMESPACE << "constructor : cannot find external parameters, did you load the YAML file?");
        abort();
    }
    
    joint_state_pub = nh.advertise<sensor_msgs::JointState>("qb_legs/joint_states",10);
    
    KDL::Chain temp;
    #if DEBUG
    std::string robot_root = robot_kdl.getRootSegment()->first; //ik_check_capability->get_robot_state().getRobotModel()->getRootLinkName();
    std::cout << "root: " << robot_root << std::endl;
    #endif
    
    for (int i=0; i<chain_ees_list.size(); i++)
    {
        std::string end_effector = chain_ees_list.at(i);
        std::string root = chain_roots_list.at(i);
        robot_kdl.getChain(root,end_effector,temp);
        
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
        
        solvers[chain_names_list.at(i)].chain = chains[chain_names_list.at(i)];
        initialize_solvers(&(solvers[chain_names_list.at(i)]),robot_kdl,i);
    }
    
    // check which robot I am using
    std::string robot_name = urdf_model.getName();
    std::cout << CLASS_NAMESPACE << " : I am using \"" << robot_name << "\" robot!" << std::endl;
}

void qb_legs_ik::parseParameters(XmlRpc::XmlRpcValue& params)
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

void qb_legs_ik::initialize_solvers(chain_and_solvers* container, const KDL::Tree& robot_kdl, int chain_index) const
{
    delete container->fksolver;
    delete container->iksolver;
    delete container->ikvelsolver;
    delete container->idsolver;
    container->joint_names.clear();
    container->tau_multiplier.clear();
    container->tau_multiplier.resize(container->chain.getNrOfJoints(),1.0);
    for (KDL::Segment& segment: container->chain.segments)
    {
        if (segment.getJoint().getType()==KDL::Joint::None) continue;
        #if DEBUG>2
        std::cout<<segment.getJoint().getName()<<std::endl;
        #endif
        container->joint_names.push_back(segment.getJoint().getName());
    }
    assert(container->joint_names.size()==container->chain.getNrOfJoints());
    container->q_max.resize(container->chain.getNrOfJoints());
    container->q_min.resize(container->chain.getNrOfJoints());
    container->fksolver=new KDL::ChainFkSolverPos_recursive(container->chain);
    container->ikvelsolver = new KDL::ChainIkSolverVel_pinv(container->chain);
    KDL::Vector gravity(0,0,9.81);
    container->idsolver = new KDL::ChainIdSolver_RNE(container->chain,gravity);
    int j=0;
    for (auto joint_name:container->joint_names)
    {
        if(urdf_model.joints_.at(joint_name)->safety)
        {
            container->q_max(j)=urdf_model.joints_.at(joint_name)->safety->soft_upper_limit;
            container->q_min(j)=urdf_model.joints_.at(joint_name)->safety->soft_lower_limit;
        }
        else
        {
            container->q_max(j)=urdf_model.joints_.at(joint_name)->limits->upper;
            container->q_min(j)=urdf_model.joints_.at(joint_name)->limits->lower;
        }
        j++;
    }
    uint max_iter = MAX_ITER;
    container->iksolver= new KDL::ChainIkSolverPos_NR_JL(container->chain,container->q_min,container->q_max,*container->fksolver,*container->ikvelsolver,max_iter,eps);
    
    std::string chain_root = container->chain.getSegment(0).getName();
    std::string robot_root = robot_kdl.getRootSegment()->first;
    // code from KDL::tree.cpp
    // walk down from chain_root to the root of the tree
    std::vector<KDL::SegmentMap::key_type> parents_chain_root;
    for (KDL::SegmentMap::const_iterator s=robot_kdl.getSegment(chain_root); s!=robot_kdl.getSegments().end(); s = s->second.parent){
        parents_chain_root.push_back(s->first);
        if (s->first == robot_root) break;
    }
    if (parents_chain_root.empty() || parents_chain_root.back() != robot_root)
    {
        ROS_ERROR_STREAM(CLASS_NAMESPACE << __func__ << " : there has been an error while looking for robot_root in chain_root ancestors...");
        abort();
    }
    
    #if DEBUG>1
    std::cout << "parents_chain_root: [";
    #endif
    int j_count = 0;
    for(int i=0; i<parents_chain_root.size() && i<container->chain.getNrOfSegments(); i++)
    {
        #if DEBUG>1
        std::cout << parents_chain_root[i] << " ";
        #endif
        
        const KDL::Segment& seg=container->chain.getSegment(i);
        if( seg.getName() != parents_chain_root[i] )
        {
            #if DEBUG>1
            continue;
            #else
            break;
            #endif
        }
        else if( seg.getJoint().getType() != KDL::Joint::None )
            container->tau_multiplier.at(j_count++) = -1;
    }
    #if DEBUG>1
    std::cout << "]" << std::endl;
    #endif
    #if DEBUG>0
    std::cout << "Changed " << j_count << " joint torques directions!" << std::endl;
    #endif
}

bool qb_legs_ik::publishConfig(const std::vector< std::string >& joint_names, const KDL::JntArray& q)
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

bool qb_legs_ik::normalizePose(geometry_msgs::Pose& pose)
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

bool qb_legs_ik::normalizePoses(std::vector< geometry_msgs::Pose >& poses)
{
    bool ok = true;
    for (auto& p:poses)
        ok = ok & normalizePose(p);
    
    return ok;
}

bool qb_legs_ik::get_ik(const std::string& chain, const KDL::Frame& ee, const KDL::JntArray& q_init, KDL::JntArray& q_out, bool publish)
{
    int res = solvers.at(chain).iksolver->CartToJnt(q_init,ee,q_out);
    
    if(publish)
    {
        publishConfig(solvers.at(chain).joint_names,q_out);
    }
    
    return (res>=0);
}

bool qb_legs_ik::get_fk(const std::string& chain, const KDL::JntArray& j, KDL::Frame& ee, bool publish)
{
    int res = solvers.at(chain).fksolver->JntToCart(j,ee);
    
    if(publish)
    {
        publishConfig(solvers.at(chain).joint_names,j);
    }
    
    return (res>=0);
}

bool qb_legs_ik::get_gravity(const std::string& chain, const KDL::JntArray& j, KDL::JntArray& tau, bool publish)
{
    int nj = solvers.at(chain).chain.getNrOfJoints();
    KDL::JntArray qzero(nj);
    KDL::Wrenches f_ext(solvers.at(chain).chain.getNrOfSegments(),KDL::Wrench(KDL::Vector::Zero(),KDL::Vector::Zero()));

    int res = solvers.at(chain).idsolver->CartToJnt(j,qzero,qzero,f_ext,tau);
    if(res<0)
    {
        ROS_ERROR_STREAM(CLASS_NAMESPACE << __func__ << " : unable to get the right ID, did you use the right dimensions?");
        return false;
    }
    
    // rectify the sign of torques for switched joint axes
    for(int i=0; i<tau.rows(); i++)
    {
        tau(i) *= solvers[chain].tau_multiplier[i];
    }
    
    if(publish)
    {
        publishConfig(solvers.at(chain).joint_names,j);
    }
    
    return true;
}
