#include <ros/ros.h>
#include <qb_legs/qb_legs_ik.h>

int main(int argc, char** argv)
{
    ros::init(argc,argv,"qb_legs_node");
    ros::NodeHandle nh;
    
    qb_legs_ik qbik;
    geometry_msgs::Pose p;
    p.position.x = 0;
    p.position.y = 0;
    p.position.z = 0;
    p.orientation.x = 10;
    p.orientation.y = 10;
    p.orientation.z = 10;
    p.orientation.w = 10;
    qbik.normalizePose(p);

    ros::Rate rate(0.1);
    while(ros::ok())
    {
        std::cout << "I'm running in the main...!" << std::endl;
        rate.sleep();
        ros::spinOnce();
    }
    return 0;
}