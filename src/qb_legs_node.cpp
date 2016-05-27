#include <ros/ros.h>
#include <qb_legs/qb_legs_ik.h>

int main(int argc, char** argv)
{
    ros::init(argc,argv,"qb_legs_node");
    ros::NodeHandle nh;
    
    qb_legs_ik qbik;

    ros::Rate rate(0.1);
    while(ros::ok())
    {
        std::cout << "I'm running in the main...!" << std::endl;
        rate.sleep();
        ros::spinOnce();
    }
    return 0;
}