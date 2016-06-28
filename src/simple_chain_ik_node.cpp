#include <ros/ros.h>
#include <simple_chain_ik/simple_chain_ik_solver.h>
#include <kdl/jntarray.hpp>
#include <kdl/kinfam_io.hpp>
#include <kdl/frames_io.hpp>
#include <kdl_conversions/kdl_msg.h>
#include <simple_chain_ik/simple_chain_ik_srv.h>

simple_chain_ik_solver* scik;

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

bool scik_srv_handler(simple_chain_ik::simple_chain_ik_srv::Request& req, simple_chain_ik::simple_chain_ik_srv::Response& res)
{
    KDL::Frame f(KDL::Frame::Identity());
    int nj = req.q_in.size();
    KDL::JntArray j_in(nj),j_out(nj),jzero(nj);
    arrayToJntArray(req.q_in,j_in);
    
    if(req.command == "ik")
    {
        tf::poseMsgToKDL(req.pose_in,f);
        res.ok = scik->get_ik(req.chain_name,f,j_in,j_out,req.publish);
        jntArrayToArray(j_out,res.q_out);
    }
    else if(req.command == "fk")
    {
        res.ok = scik->get_fk(req.chain_name,j_in,f,req.publish);
        tf::poseKDLToMsg(f,res.pose_out);
    }
    else if(req.command == "grav")
    {
        std::map<std::string,KDL::Wrench> w_ext;
        KDL::Wrench w_tmp;
        for(const geometry_msgs::WrenchStamped& w:req.w_ext)
        {
            tf::wrenchMsgToKDL(w.wrench,w_tmp);
            w_ext[w.header.frame_id] = w_tmp;
        }
        
        res.ok = scik->get_gravity(req.chain_name,j_in,w_ext,j_out,req.publish);
        jntArrayToArray(j_out,res.tau);
    }
    else
    {
        ROS_ERROR_STREAM(__func__ << " : unrecognized command \'" << req.command << "\'...");
        tf::poseKDLToMsg(f,res.pose_out);
        res.q_out.clear();
        res.tau.clear();
        res.ok = false;
    }
    if(!res.ok)
        ROS_ERROR_STREAM(__func__ << " : There was an error!");
    
    return res.ok;
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
    
    scik = new simple_chain_ik_solver;
    
    ros::ServiceServer service_srv;
    service_srv = nh.advertiseService("simple_chain_ik_srv", &scik_srv_handler);
    
    std::cout << "simple_chain_ik_node running!!!" << std::endl;
    
    ros::Rate rate(0.1);
    while(ros::ok())
    {
        rate.sleep();
    }
    return 0;
}