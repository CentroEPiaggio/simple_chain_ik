#include <ros/ros.h>
#include <qb_legs/qb_legs_ik.h>
#include <kdl/jntarray.hpp>
#include <kdl/kinfam_io.hpp>
#include <kdl/frames_io.hpp>
#include <kdl_conversions/kdl_msg.h>

qb_legs_ik qbik;

void arrayToJntArray(const std::vector<double>& in, KDL::JntArray& out)
{
    out.resize(in.size());
    for(int i=0; i<in.size(); i++)
        out(i) = in[i];
}

void jntArrayToArray(const KDL::JntArray& in, std::vector<double>& out)
{
    out.resize(in.rows());
    for(int i=0; i<in.rows(); i++)
        out[i] = in(i);
}

void qb_srv_handler()
{
    
}

int main(int argc, char** argv)
{
    ros::init(argc,argv,"qb_legs_node");
    ros::NodeHandle nh;
    ros::AsyncSpinner aspin(1);
    aspin.start();
    
    KDL::JntArray q_in(4),q_out(4),tau(4);
    q_in(0) = M_PI/2.0;
    q_in(1) = 0.0*(-M_PI/2.0);
    
    KDL::Frame ee;
    bool publish = true;
    qbik.get_fk("lfoot_torso",q_in,ee,publish);
    
    qbik.get_gravity("lfoot_torso",q_in,tau);
    
    std::cout << "q_in: " << q_in << std::endl;
    std::cout << "ee_pose: " << ee << std::endl;
    std::cout << "tau: " << tau << std::endl;
    
    ros::Rate rate(0.1);
    while(ros::ok())
    {
        std::cout << "I'm running in the main...!" << std::endl;
        rate.sleep();
    }
    return 0;
}