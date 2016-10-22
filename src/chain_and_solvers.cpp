#include "simple_chain_ik/chain_and_solvers.h"
#include <iostream>
#include <kdl/frames_io.hpp>
#include <kdl/kinfam_io.hpp>
#include <kdl/treefksolverpos_recursive.hpp>

#define CLASS_NAMESPACE "ChainAndSolvers::"
#define DEBUG 0
#define DEBUG_STRING {std::cout << CLASS_NAMESPACE << __func__ << "@" << __LINE__ << std::endl;}
#define USING_CUSTOM_SOLVERS 1

ChainAndSolvers::ChainAndSolvers(const std::shared_ptr< KDL::Tree >& tree_, const std::shared_ptr< KDL::TreeFkSolverPos >& tree_fk_, const std::string& chain_root_, const std::string& chain_tip_, const KDL::Vector& tree_root_gravity_) : tree(tree_), tree_fk(tree_fk_), chain_root(chain_root_), chain_tip(chain_tip_), tree_root_gravity(tree_root_gravity_)
{
    if(!onChainChanged()) abort();
}

ChainAndSolvers::ChainAndSolvers(const std::shared_ptr< KDL::Tree >& tree_, const std::string& chain_root_, const std::string& chain_tip_, const KDL::Vector& tree_root_gravity_) : tree(tree_), chain_root(chain_root_), chain_tip(chain_tip_), tree_root_gravity(tree_root_gravity_)
{
    tree_fk.reset(new KDL::TreeFkSolverPos_recursive(*tree));
    
    if(!onChainChanged()) abort();
}

ChainAndSolvers::ChainAndSolvers(const KDL::Chain& chain_, const KDL::Vector& tree_root_gravity_) : tree_root_gravity(tree_root_gravity_)
{
    /// a fake_root is needed in the tree in order to get a chain which is equal to the one passed as input (as it will not contain the first segment)
    tree.reset(new KDL::Tree("fake_root"));
    tree->addChain(chain_,"fake_root");
    tree_fk.reset(new KDL::TreeFkSolverPos_recursive(*tree));
    chain_root = "fake_root";
    chain_tip = chain_.segments.back().getName();
    
    if(!onChainChanged()) abort();
}

bool ChainAndSolvers::onChainChanged()
{
    initialized = false;
    
    /// reset variables connected to the extra frame at the end-effector
    ee_tip = KDL::Frame::Identity();
    Mx.setIdentity(Mx.RowsAtCompileTime,Mx.ColsAtCompileTime);
    Wx.setOnes(Wx.RowsAtCompileTime,Wx.ColsAtCompileTime);
    
    if (!tree->getChain(chain_root,chain_tip,chain_orig))
    {
        std::cout << CLASS_NAMESPACE << __func__ << " : unable to get Chain from Tree - aborting!" << std::endl;
        return false;
    }
#if DEBUG>2
    std::cout << "chain_root: " << chain_root << std::endl;
    std::cout << "chain_tip: " << chain_tip << std::endl;
    std::cout << "chain segment names: ";
    for(auto& seg:chain_orig.segments)
        std::cout << seg.getName() << " ";
    std::cout << std::endl;
#endif
    
    joint_names.clear();
    for (const KDL::Segment& segment: chain_orig.segments)
    {
        if (segment.getJoint().getType()==KDL::Joint::None)
            continue;
        
        #if DEBUG>2
        std::cout << segment.getJoint().getName() << " " << std::endl;
        #endif
        joint_names.push_back(segment.getJoint().getName());
    }
    
    assert(joint_names.size() == chain_orig.getNrOfJoints());
    
    if(!computeTauMultipliers() || !computeLocalGravity(tree_root_gravity))
        return false;
    
    return true;
}

bool ChainAndSolvers::setSolverParameters(const KDL::JntArray& q_min_, const KDL::JntArray& q_max_, uint pos_IK_max_iter_, double pos_IK_eps_, uint vel_IK_max_iter_, double vel_IK_eps_, double vel_IK_lambda_)
{
    initialized = false;
    
    pos_IK_max_iter = pos_IK_max_iter_;
    pos_IK_eps = pos_IK_eps_;
    vel_IK_max_iter = vel_IK_max_iter_;
    vel_IK_eps = vel_IK_eps_;
    vel_IK_lambda = vel_IK_lambda_;
    q_max = q_max_;
    q_min = q_min_;
    
    return true;
}

bool ChainAndSolvers::computeLocalGravity(const KDL::Vector& tree_root_gravity_)
{
    KDL::JntArray tree_j(tree->getNrOfJoints());
    KDL::Frame robotroot_chainroot;
    if(tree_fk->JntToCart(tree_j,robotroot_chainroot,chain_root) < 0)
    {
        std::cout << CLASS_NAMESPACE << __func__ << " : unable to compute chain_root position in robot_root frame! Returning..." << std::endl;
        return false;
    }
    
    gravity = robotroot_chainroot.M.Inverse(tree_root_gravity_);
    #if DEBUG>2
    std::cout << "robotroot_chainroot: " << robotroot_chainroot << std::endl;
    #endif
    #if DEBUG>1
    std::cout << "local gravity in the chain root: " << gravity << std::endl;
    #endif
    
    return true;
}

bool ChainAndSolvers::computeTauMultipliers()
{
    tau_multiplier.clear();
    tau_multiplier.resize(chain_orig.getNrOfJoints(),1.0);
    
    std::string chain_root = chain_orig.getSegment(0).getName();
    std::string robot_root = tree->getRootSegment()->first;
    // code from KDL::tree.cpp
    // walk down from chain_root to the root of the tree
    std::vector<KDL::SegmentMap::key_type> parents_chain_root;
    for (KDL::SegmentMap::const_iterator s = tree->getSegment(chain_root); s != tree->getSegments().end(); s = s->second.parent){
        parents_chain_root.push_back(s->first);
        if (s->first == robot_root) break;
    }
    if (parents_chain_root.empty() || parents_chain_root.back() != robot_root)
    {
        std::cout << CLASS_NAMESPACE << __func__ << " : there has been an error while looking for robot_root in chain_root ancestors..." << std::endl;
        return false;
    }
    
    #if DEBUG>1
    std::cout << "parents_chain_root: [";
    #endif
    int j_count = 0;
    for(int i=0; i<parents_chain_root.size() && i<chain_orig.getNrOfSegments(); i++)
    {
        #if DEBUG>1
        std::cout << parents_chain_root[i] << " ";
        #endif
        
        const KDL::Segment& seg=chain_orig.getSegment(i);
        if( seg.getName() != parents_chain_root[i] )
        {
            #if DEBUG>1
            continue;
            #else
            break;
            #endif
        }
        else if( seg.getJoint().getType() != KDL::Joint::None )
            tau_multiplier.at(j_count++) = -1;
    }
    #if DEBUG>1
    std::cout << "]" << std::endl;
    #endif
    #if DEBUG>0
    std::cout << "Changed " << j_count << " joint torques directions!" << std::endl;
    #endif
    
    return true;
}

bool ChainAndSolvers::initSolvers()
{
    if(q_max.rows() != q_min.rows() || q_max.rows() != chain_orig.getNrOfJoints())
    {
        std::cout << CLASS_NAMESPACE << __func__ << " : q_max (" << q_max.rows() << ") and q_min (" << q_min.rows() << ") must have the same size, equal to the number of joints in the chain (" << chain_orig.getNrOfJoints() << ")!" << std::endl;
        return false;
    }
    
    chain = chain_orig;
    chain.addSegment(KDL::Segment("ee_tip",KDL::Joint(KDL::Joint::None),ee_tip,KDL::RigidBodyInertia::Zero()));
    
    fksolver.reset(new ChainAndSolversTypes::ChainFkSolverPos(chain));
    ikvelsolver.reset(new ChainAndSolversTypes::ChainIkSolverVel(chain,vel_IK_eps,vel_IK_max_iter));
    ikvelsolver->setLambda(vel_IK_lambda);
    ikvelsolver->setWeightTS(Mx);
    ikvelsolver->setJointLimits(q_min.data,q_max.data);
    iksolver.reset(new ChainAndSolversTypes::ChainIkSolverPos(chain,q_min,q_max,*fksolver,*ikvelsolver,pos_IK_max_iter,pos_IK_eps));
#if USING_CUSTOM_SOLVERS>0
    iksolver->setTaskWeight(Wx);
#endif
    idsolver.reset(new ChainAndSolversTypes::ChainIdSolver(chain,gravity));
    
    initialized = true;
    
    return true;
}

const std::vector< std::string >& ChainAndSolvers::jointNames()
{
    return joint_names;
}

bool ChainAndSolvers::getGravity(const KDL::JntArray& j, const std::map< std::string, KDL::Wrench >& w_ext, KDL::JntArray& tau)
{
    if(!initialized)
        return false;
    
    int nj = chain.getNrOfJoints();
    KDL::JntArray qzero(nj);
    KDL::Wrenches f_ext(chain.getNrOfSegments(),KDL::Wrench(KDL::Vector::Zero(),KDL::Vector::Zero()));
    
    // add the external wrench in the appropriate position
    int counter = 0, w_counter = w_ext.size();
    for(KDL::Segment seg:chain.segments)
    {
        if(w_counter <= 0)
            break;
        if(w_ext.count(seg.getName()))
        {
            f_ext.at(counter) = w_ext.at(seg.getName());
            w_counter--;
        }
        counter++;
    }
    
    int res = idsolver->CartToJnt(j,qzero,qzero,f_ext,tau);
    if(res<0)
    {
        std::cout << CLASS_NAMESPACE << __func__ << " : unable to get the right ID, did you use the right dimensions?" << std::endl;
        return false;
    }
    
    // rectify the sign of torques for switched joint axes
    for(int i=0; i<tau.rows(); i++)
    {
        tau(i) *= tau_multiplier[i];
    }
    
    return true;
}

std::unique_ptr< ChainAndSolversTypes::ChainFkSolverPos >& ChainAndSolvers::getFKSolver()
{
    if(initialized)
        return fksolver;
    
    std::cout << CLASS_NAMESPACE << __func__ << " : getting non-initialized solver not allowed!" << std::endl;
    abort();
}

std::unique_ptr< ChainAndSolversTypes::ChainIkSolverPos >& ChainAndSolvers::getIKSolver()
{
    if(initialized)
        return iksolver;
    
    std::cout << CLASS_NAMESPACE << __func__ << " : getting non-initialized solver not allowed!" << std::endl;
    abort();
}

std::unique_ptr< ChainAndSolversTypes::ChainIkSolverVel >& ChainAndSolvers::getIKVelSolver()
{
    if(initialized)
        return ikvelsolver;
    
    std::cout << CLASS_NAMESPACE << __func__ << " : getting non-initialized solver not allowed!" << std::endl;
    abort();
}

KDL::JntArray ChainAndSolvers::getValidRandomJoints()
{
    if(!initialized)
        return KDL::JntArray();
    
    KDL::JntArray ret(q_max.rows());
    ret.data.setRandom(ret.data.rows(),ret.data.cols());
    ret.data.array() += 1;
    ret.data /= 2.0;
    ret.data = q_min.data + ret.data.cwiseProduct(q_max.data - q_min.data);
    return ret;
}

void ChainAndSolvers::changeTip(const KDL::Frame& ee_tip_)
{
    initialized = false;
    ee_tip = ee_tip_;
}

bool ChainAndSolvers::changeIkTaskWeigth(const Eigen::Matrix<double,6,1>& Wx_)
{
    Wx = Wx_;
    Mx = Wx_.asDiagonal();
    if(initialized)
    {
#if USING_CUSTOM_SOLVERS>0
        iksolver->setTaskWeight(Wx);
#endif
        return (ikvelsolver->setWeightTS(Mx) == KDL::SolverI::E_NOERROR);
    }
    
    return true;
}
