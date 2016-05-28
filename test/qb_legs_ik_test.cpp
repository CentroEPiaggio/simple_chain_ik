#include <ros/ros.h>
#include <qb_legs/qb_legs_ik.h>
#include <qb_legs/qb_legs_srv.h>
#include <kdl/jntarray.hpp>
#include <kdl/kinfam_io.hpp>
#include <kdl/frames_io.hpp>
#include <kdl_conversions/kdl_msg.h>

qb_legs_ik* qbik;

int main(int argc, char** argv)
{
    while(!ros::isInitialized())
    {
        ros::init(argc,argv,"qb_legs_ik_test");
    }
    ros::NodeHandle nh;
    ros::AsyncSpinner aspin(1);
    aspin.start();
    
    sleep(2);
    
    qbik = new qb_legs_ik;
    
    std::cout << "qb_legs_ik_test running!!!" << std::endl;
    
    // test 1
    std::string chain("lfoot_rfoot");
    KDL::JntArray j(6);
    j(3) = 1.57;
    std::map<std::string,KDL::Wrench> w_ext({{"qb_legs_right_foot",KDL::Wrench(KDL::Vector(1.0,0.0,0.0),KDL::Vector::Zero())}});
    KDL::JntArray tau(j.rows());
    bool publish(true);
    qbik->get_gravity(chain,j,w_ext,tau,publish);
    std::cout << "Test 1: chain \'" << chain << "\' with external wrench on " << w_ext.begin()->first << " of " << w_ext.begin()->second << std::endl;
    std::cout << "tau: " << tau << std::endl;
    
    // test 2
    chain = "lfoot_torso";
    j.resize(4);
    KDL::SetToZero(j);
    w_ext.clear();
    w_ext["qb_legs_torso_link"] = KDL::Wrench(KDL::Vector(1.0,0.0,0.0),KDL::Vector::Zero());
    tau.resize(j.rows());
    qbik->get_gravity(chain,j,w_ext,tau,publish);
    std::cout << "Test 2: chain \'" << chain << "\' with external wrench on " << w_ext.begin()->first << " of " << w_ext.begin()->second << std::endl;
    std::cout << "tau: " << tau << std::endl;
    
    ros::Rate rate(0.1);
    while(ros::ok())
    {
        rate.sleep();
    }
    return 0;
}