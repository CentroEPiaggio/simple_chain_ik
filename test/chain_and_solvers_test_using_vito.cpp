#include <simple_chain_ik/chain_and_solvers.h>
#include <ros/ros.h>
#include <urdf/model.h>
#include <kdl_parser/kdl_parser.hpp>
#include <kdl/treefksolverpos_recursive.hpp>
#include <sensor_msgs/JointState.h>
#include <kdl/frames_io.hpp>

#define CLASS_NAMESPACE "chain_and_solvers_test::"

KDL::Frame ee_tip = KDL::Frame(KDL::Vector(0.0,0.0,0.2));

void initialize_solvers(ChainAndSolvers& container, const urdf::Model& urdf_model)
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
    
    if(!container.setSolverParameters(q_min,q_max,100,5e-4,150,1e-5,1e-5) || !container.initSolvers())
    {
        std::cout << CLASS_NAMESPACE << __func__ << " : unable to initialize the solvers! Returning..." << std::endl;
        abort();
    }
}

bool publishConfig(const std::vector< std::string >& joint_names, const KDL::JntArray& q, const ros::Publisher& joint_state_pub)
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

void computeAndDisplayDifference(const KDL::Frame& target, const KDL::JntArray& q_actual, ChainAndSolvers& solvers)
{
    KDL::Frame actual;
    solvers.getFKSolver()->JntToCart(q_actual,actual);
    KDL::Twist xi;
    xi = KDL::diff(target,actual);
    std::cout << "Difference: " << xi << std::endl;
    std::cout << "Difference norm [rot+lin]: " << xi.rot.Norm() << "+" << xi.vel.Norm() << std::endl;
}

int main(int argc, char** argv)
{
    while(!ros::isInitialized())
    {
        ros::init(argc,argv,"simple_chain_ik_node");
    }
    ros::NodeHandle nh;
    ros::AsyncSpinner aspin(1);
    aspin.start();
    ros::Publisher joint_state_pub = nh.advertise<sensor_msgs::JointState>("simple_chain_ik/joint_states",10,true);
    
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
    
    ChainAndSolvers lh_solver(tree,tree_fk,"vito_anchor","left_hand_palm_link",KDL::Vector(0.0,0.0,-9.81));
    
    std::cout << "chain_and_solvers_test running!!!" << std::endl;
    
    lh_solver.changeTip(ee_tip);
    initialize_solvers(lh_solver,urdf_model);
    
    KDL::Frame target;
    double init_y(-0.1), end_y(0.5);
    target.p = KDL::Vector(-0.8 , init_y , 0.07);
    target.M = KDL::Rotation::RotY(M_PI); // x forward, z downward
    target.M.DoRotZ(M_PI/2.0);
    
    int nj = lh_solver.jointNames().size();
    KDL::JntArray q_out(nj),q_init(nj);
    q_init = lh_solver.getValidRandomJoints();
    int solver_error = lh_solver.getIKSolver()->CartToJnt(q_init,target,q_out);
    if(solver_error < 0)
    {
        ROS_ERROR_STREAM(CLASS_NAMESPACE << " : unable to get first IK! Error " << solver_error << " > " << lh_solver.getIKSolver()->strError(solver_error));
        return -10;
    }
    publishConfig(lh_solver.jointNames(),q_out,joint_state_pub);
    computeAndDisplayDifference(target,q_out,lh_solver);
    
    ros::Duration tsleep(0.1);
    tsleep.sleep();
    
    // iterate with some frames
    int counter=0;
    while(target.p.y() < end_y)
    {
        ROS_INFO_STREAM(CLASS_NAMESPACE << " : computing IK #" << ++counter);
        target.p.y( target.p.y() + 0.025 );
        // use old q_out as new seed
        solver_error = lh_solver.getIKSolver()->CartToJnt(q_out,target,q_out);
        if(solver_error < 0)
        {
            ROS_ERROR_STREAM(CLASS_NAMESPACE << " : unable to get IK! Error " << solver_error << " > " << lh_solver.getIKSolver()->strError(solver_error));
            ROS_WARN_STREAM(CLASS_NAMESPACE << " : was able to move until y=" << target.p.y());
            break;
        }
        publishConfig(lh_solver.jointNames(),q_out,joint_state_pub);
        computeAndDisplayDifference(target,q_out,lh_solver);
        
        tsleep.sleep();
    }
    
    ROS_INFO_STREAM(CLASS_NAMESPACE << " : trying again adding some more tolerance for rotation...");
    Eigen::Matrix<double,6,1> Mx;
    Mx.setOnes(Mx.RowsAtCompileTime,Mx.ColsAtCompileTime);
//     Mx(1,0) = 0.0;
//     Mx(2,0) = 0.0;
    Mx(3,0) = 0.0; // care less about x-rotation
    Mx(4,0) = 0.0; // care less about y-rotation
//     Mx(5,0) = 0.0; // care less about z-rotation
    
    target.p.y( init_y );
    
    // iterate with some frames
    counter=0;
    solver_error = lh_solver.getIKSolver()->CartToJnt(q_init,target,q_out);
    if(solver_error < 0)
    {
        ROS_ERROR_STREAM(CLASS_NAMESPACE << " : unable to get first IK! Error " << solver_error << " > " << lh_solver.getIKSolver()->strError(solver_error));
        return -10;
    }
    publishConfig(lh_solver.jointNames(),q_out,joint_state_pub);
    computeAndDisplayDifference(target,q_out,lh_solver);
    tsleep.sleep();
    
    if(!lh_solver.changeIkTaskWeigth(Mx))
        ROS_WARN_STREAM(CLASS_NAMESPACE << " : could not change TS weight as the matrix is not positive (semi-)definite!");
    
    while(target.p.y() < end_y)
    {
        ROS_INFO_STREAM(CLASS_NAMESPACE << " : computing relaxed IK #" << ++counter);
        target.p.y( target.p.y() + 0.025 );
        // use old q_out as new seed
        solver_error = lh_solver.getIKSolver()->CartToJnt(q_out,target,q_out);
        if(solver_error < 0)
        {
            ROS_ERROR_STREAM(CLASS_NAMESPACE << " : unable to get IK! Error " << solver_error << " > " << lh_solver.getIKSolver()->strError(solver_error));
            ROS_WARN_STREAM(CLASS_NAMESPACE << " : was able to move until y=" << target.p.y());
            break;
        }
        publishConfig(lh_solver.jointNames(),q_out,joint_state_pub);
        computeAndDisplayDifference(target,q_out,lh_solver);
        
        tsleep.sleep();
    }
    
    return 0;
}
