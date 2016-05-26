#include "qb_legs/qb_legs_ik.h"
#include <ros/ros.h>
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
// #include <dual_manipulation_shared/parsing_utils.h>
// #include "dual_manipulation_shared/serialization_utils.h"

#define HIGH 0.35
#define LOW 0.06
#define BOX_CONSTR_SIDE 0.2
#define ANGLE_STEPS 4.0 // 6.0
#define BIMANUAL_IK_ATTEMPTS 3
#define BIMANUAL_IK_TIMEOUT 0.005
#define OBJ_GRASP_FACTOR 1000
#define EXTRA_TIME_TO_GO_HOME 0.0

#define DEBUG 1 // if 1, print some more information
#define SHOW_IK 2
#define MAX_ITER 100
#define EPS 5e-3
#define MULTI_OBJECT_PLANNING 1 // HOME commands should be blocking with multi-object planner
#define SINGLE_MOVEMENT_DURATION 10 // consider 10 seconds for each movement... maybe less?

bool am_I_Vito = false;
double eps = EPS;

qb_legs_ik::qb_legs_ik()
{
    std::cout << "building qb_legs_ik!!" << std::endl;
    //   this->database=database;
    ros::NodeHandle nh;
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
        ROS_ERROR_STREAM("qb_legs_ik::constructor : cannot find robot_description");
        abort();
    }
    if (!urdf_model.initParam("robot_description"))
    {
        ROS_ERROR_STREAM("qb_legs_ik::constructor : cannot load robot_description");
        abort();
    }
    
    if (!kdl_parser::treeFromUrdfModel(urdf_model, robot_kdl))
    {
        ROS_ERROR_STREAM("Failed to construct kdl tree");
        abort();
    }
    
//     if (nh.getParam("ik_control_parameters", ik_control_params))
//         parseParameters(ik_control_params);
    
    KDL::Chain temp;
    std::vector<std::string> link_names({"qb_legs_torso_box","qb_legs_left_mid_shaft","qb_legs_right_foot"});
    
    for(auto& seg:robot_kdl.getSegments())
        std::cout << seg.first << " ";
    std::cout << std::endl;
    
//     const moveit::core::JointModelGroup* jmg = ik_check_capability->get_robot_state().getRobotModel()->getJointModelGroup("full_robot");
//     jmg->getEndEffectorTips(link_names);
    std::vector<std::string> ee_names({"torso","left_foot","right_foot"}); // = jmg->getAttachedEndEffectorNames();
    std::string root = robot_kdl.getRootSegment()->first; //ik_check_capability->get_robot_state().getRobotModel()->getRootLinkName();
    
    #if DEBUG
    std::cout << "root: " << root << std::endl;
    #endif
    
    for (int i=0; i<link_names.size(); i++)
    {
        std::string end_effector = link_names.at(i);
        robot_kdl.getChain(root,end_effector,temp);
        std::string seg_fake_name;
        bool fake_next=true;
        for (auto s: temp.segments)
        {
            KDL::Joint j = s.getJoint();
            
            #if DEBUG>1
            std::cout << "s.getName(): " << s.getName() << std::endl;
            std::cout << "s.getJoint().getName(): " << s.getJoint().getName() << std::endl;
            std::cout << "s.getJoint().JointOrigin(): " << s.getJoint().JointOrigin().x() << " " << s.getJoint().JointOrigin().y() << " " << s.getJoint().JointOrigin().z() << std::endl;
            std::cout << "s.getJoint().JointAxis(): " << s.getJoint().JointAxis().x() << " "<< s.getJoint().JointAxis().y() << " "<< s.getJoint().JointAxis().z() << std::endl;
            std::cout << "s.getJoint().getType(): " << s.getJoint().getType() << std::endl;
            KDL::Frame f = s.getFrameToTip();
            std::cout << f.p.data[0] << " "<< f.p.data[1] << " "<< f.p.data[2] << std::endl;
            double r,p,y; f.M.GetRPY(r,p,y);
            std::cout << r << " " << p << " " << y << std::endl;
            #endif
            
            if (fake_next)
            {
                if (s.getJoint().getType()==KDL::Joint::None)
                {
                    seg_fake_name=end_effector;
                    fake_next=true;
                    j=KDL::Joint(s.getJoint().getName()+end_effector);
                }
                else if (s.getJoint().getType()!=KDL::Joint::None && fake_next)
                {
                    seg_fake_name="";
                    fake_next=false;
                }
            }
            else
            {
                seg_fake_name="";
                fake_next=false;
            }
            KDL::Segment b(s.getName()+seg_fake_name,j,s.getFrameToTip(),s.getInertia());
            chains[ee_names.at(i)].addSegment(b);
        }
        KDL::Tree t("fake_root");
        bool done = t.addChain(chains[ee_names.at(i)],"fake_root");
        assert(done);
        done = t.getChain(end_effector,"fake_root",chains_reverse[ee_names.at(i)]);
        
        #if DEBUG>1
        temp = chains_reverse.at(ee_names.at(i));
        for (auto s: temp.segments)
        {
            
            std::cout << "s.getName(): " << s.getName() << std::endl;
            KDL::Joint j = s.getJoint();
            
            std::cout << "s.getJoint().getName(): " << s.getJoint().getName() << std::endl;
            KDL::Frame f = s.getFrameToTip();
            std::cout << f.p.data[0] << " "<< f.p.data[1] << " "<< f.p.data[2] << std::endl;
            double r,p,y; f.M.GetRPY(r,p,y);
            std::cout << r << " " << p << " " << y << std::endl;
        }
        #endif
        #if DEBUG
        // std::cout << "t.getNrOfJoints(): " << t.getNrOfJoints() << std::endl;
        std::cout << "chains[ee_names.at(i)].getNrOfJoints(): " << chains[ee_names.at(i)].getNrOfJoints() << std::endl;
        std::cout << "chains_reverse[ee_names.at(i)].getNrOfJoints(): " << chains_reverse[ee_names.at(i)].getNrOfJoints() << std::endl;
        std::cout << "chains[" << ee_names.at(i) << "].segments: | ";
        for(auto segs:chains[ee_names.at(i)].segments)
            std::cout << segs.getName() << " | ";
        std::cout << std::endl;
        std::cout << "chains_reverse[" << ee_names.at(i) << "].segments: | ";
        for(auto segs:chains_reverse[ee_names.at(i)].segments)
            std::cout << segs.getName() << " | ";
        std::cout << std::endl;
//         std::cout << "t.segments: | ";
//         for(auto segs:t.getSegments())
//             std::cout << segs.first << " | ";
//         std::cout << std::endl;
        std::cout << "chains[" << ee_names.at(i) << "].segments masses: | ";
        for(auto segs:chains[ee_names.at(i)].segments)
            std::cout << segs.getInertia().getMass() << " | ";
        std::cout << std::endl;
        std::cout << "chains_reverse[" << ee_names.at(i) << "].segments masses: | ";
        for(auto segs:chains_reverse[ee_names.at(i)].segments)
            std::cout << segs.getInertia().getMass() << " | ";
        std::cout << std::endl;
        if(!done)
        {
            std::cout << "qb_legs_ik : unable to construct chains_reverse[" << ee_names.at(i) << "] - aborting" << std::endl;
            abort();
        }
        #endif
        
        if(ee_names.at(i) == "torso")
            continue;
        // do some tests with gravity
        KDL::Vector gravity(0,0,-9.81);
        KDL::ChainIdSolver_RNE cids(chains[ee_names.at(i)],gravity);
        KDL::ChainIdSolver_RNE cids2(chains_reverse[ee_names.at(i)],gravity);
        
        int nj = chains[ee_names.at(i)].getNrOfJoints();
        KDL::JntArray qs(nj),qzero(nj),tau(nj);
        qs.data[0] = 0.785*2.0;
//         qs.data[1] = 0.75;
//         qs.data[2] = 0.75;
        
//         qs.rows();
        
        int res;
        KDL::Chain tmp; tmp.getNrOfSegments();
        res = cids.CartToJnt(qs,qzero,qzero,KDL::Wrenches(chains[ee_names.at(i)].getNrOfSegments(),KDL::Wrench(KDL::Vector::Zero(),KDL::Vector::Zero())),tau);
        std::cout << "Tau straight (res=" << res << "): " << tau << std::endl;
        
        
//         qs.data[2]=qs.data[0];
//         qs.data[0]=0.0;
        res = cids2.CartToJnt(qs,qzero,qzero,KDL::Wrenches(chains[ee_names.at(i)].getNrOfSegments(),KDL::Wrench(KDL::Vector::Zero(),KDL::Vector::Zero())),tau);
        std::cout << "Tau reverse (res=" << res << "): " << tau << std::endl;
        
    }
    
    // check whether I am using Vito
    std::string robot_name = urdf_model.getName();
    std::cout << "qb_legs_ik : I am using \"" << robot_name << "\" robot!" << std::endl;
}

// void qb_legs_ik::parseParameters(XmlRpc::XmlRpcValue& params)
// {
//     double goal_position_tolerance = EPS, goal_orientation_tolerance = EPS;
//     ROS_ASSERT(params.getType() == XmlRpc::XmlRpcValue::TypeStruct);
//     parseSingleParameter(params,chain_names_list,"chain_group_names",1);
//     parseSingleParameter(params,goal_position_tolerance,"goal_position_tolerance");
//     parseSingleParameter(params,goal_orientation_tolerance,"goal_orientation_tolerance");
//     
//     eps = std::min(goal_position_tolerance,goal_orientation_tolerance);
// }

void qb_legs_ik::initialize_solvers(chain_and_solvers* container) const
{
    delete container->fksolver;
    delete container->iksolver;
    delete container->ikvelsolver;
    container->joint_names.clear();
    for (KDL::Segment& segment: container->chain.segments)
    {
        if (segment.getJoint().getType()==KDL::Joint::None) continue;
        #if DEBUG
        std::cout<<segment.getJoint().getName()<<std::endl;
        #endif
        container->joint_names.push_back(segment.getJoint().getName());
    }
    assert(container->joint_names.size()==container->chain.getNrOfJoints());
    container->q_max.resize(container->chain.getNrOfJoints());
    container->q_min.resize(container->chain.getNrOfJoints());
    container->fksolver=new KDL::ChainFkSolverPos_recursive(container->chain);
    container->ikvelsolver = new KDL::ChainIkSolverVel_pinv(container->chain);
    int j=0;
    for (auto joint_name:container->joint_names)
    {
        #if IGNORE_JOINT_LIMITS
        container->q_max(j)=M_PI/3.0;
        container->q_min(j)=-M_PI/3.0;
        #else
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
        #endif
        j++;
    }
    if (am_I_Vito) //Particular case of Vito robot
    {
        //TODO impose some new limits on shoulder joints of both arms
        // double LSh0_mean = 0.5, LSh1_mean = 1.0;
        // double RSh0_mean = -0.5, RSh1_mean = 1.0;
        int start_ind, inc;
        if (container->joint_names.at(0).find("left") == 0)
        {
            start_ind = 10; inc = 1;
        }
        else if (container->joint_names.at(0).find("right") == 0)
        {
            start_ind = 3; inc = -1;
        }
        else
        {
            std::cout << "Vito should start with something else! (container->joint_names.at(0) = " << container->joint_names.at(0) << ")" << std::endl;
            abort();
        }
        double allowed_range = 0.25;
        container->q_min(start_ind) = 1.4 - allowed_range;
        container->q_max(start_ind) = 1.4 + allowed_range;
        container->q_min(start_ind+inc) = 1.8 - allowed_range;
        container->q_max(start_ind+inc) = 1.8 + allowed_range;
        container->q_min(start_ind+2*inc) = 1.0 - allowed_range/2.0;
        container->q_max(start_ind+2*inc) = 1.0 + allowed_range/2.0;
        container->q_min(start_ind+3*inc) = 0.5 - allowed_range/2.0;
        container->q_max(start_ind+3*inc) = 0.5 + allowed_range/2.0;
    }
    uint max_iter = MAX_ITER;
    container->iksolver= new KDL::ChainIkSolverPos_NR_JL(container->chain,container->q_min,container->q_max,*container->fksolver,*container->ikvelsolver,max_iter,eps);
}

bool qb_legs_ik::publishConfig(const std::vector< std::string >& joint_names, const KDL::JntArray& q) const
{
//     static ros::Publisher joint_state_pub_;
//     static bool pub_initialized(false);
//     static ros::NodeHandle node;
//     if (!pub_initialized)
//     {
//         joint_state_pub_ = node.advertise<sensor_msgs::JointState>("sem2cart/joint_states",10);
//         pub_initialized = true;
//     }
//     sensor_msgs::JointState js_msg;
//     js_msg.name = joint_names;
//     js_msg.header.stamp = ros::Time::now();
//     js_msg.position.clear();
//     for(int i=0; i<js_msg.name.size(); i++)
//     {
//         js_msg.position.push_back(q(i));
//     }
//     joint_state_pub_.publish(js_msg);
//     
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
