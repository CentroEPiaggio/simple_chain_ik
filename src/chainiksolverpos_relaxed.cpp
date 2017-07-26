// Copyright (c) 2017, Centro di Ricerca "E. Piaggio", University of Pisa
// All rights reserved.
// Copyright (c) 2017, George Jose Pollayil <gpollayil at gmail dot com>
// Copyright (c) 2017, Mathew Jose Pollayil <mathewjosepollayil at gmail dot com>
// Copyright  (C)  2016 Hamal Marino <hamal dot marino at gmail dot com>
// Copyright  (C)  2007-2008  Ruben Smits <ruben dot smits at mech dot kuleuven dot be>
// Copyright  (C)  2008  Mikael Mayer
// Copyright  (C)  2008  Julia Jesse

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

// This class is based on KDL::ChainIkSolverPos_NR_JL, extended with task weighting

#include <simple_chain_ik/solvers/chainiksolverpos_relaxed.hpp>

#include <limits>
#include <iostream>
#include <kdl/frames_io.hpp>
#include <math.h>
#include <Eigen/Core>

#define DEBUG 0

using namespace KDL;

ChainIkSolverPos_relaxed::ChainIkSolverPos_relaxed(const Chain& chain, const JntArray& q_min, const JntArray& q_max, ChainFkSolverPos& fksolver, ChainIkSolverVel& iksolver, unsigned int maxiter, double eps) : chain(chain), nj(chain.getNrOfJoints()), q_min(q_min), q_max(q_max), iksolver(iksolver), fksolver(fksolver), delta_q(chain.getNrOfJoints()), maxiter(maxiter), eps(eps), W(Eigen::MatrixXd::Ones(W.RowsAtCompileTime,W.ColsAtCompileTime)), use_ee_task_(false)
{
}

ChainIkSolverPos_relaxed::ChainIkSolverPos_relaxed(const Chain& chain, ChainFkSolverPos& fksolver, ChainIkSolverVel& iksolver, unsigned int maxiter, double eps) : chain(chain), nj(chain.getNrOfJoints()), q_min(nj), q_max(nj), iksolver(iksolver), fksolver(fksolver), delta_q(chain.getNrOfJoints()), maxiter(maxiter), eps(eps), W(Eigen::MatrixXd::Ones(W.RowsAtCompileTime,W.ColsAtCompileTime)), use_ee_task_(false)
{
    q_min.data.setConstant(std::numeric_limits<double>::min());
    q_max.data.setConstant(std::numeric_limits<double>::max());
}

void ChainIkSolverPos_relaxed::getAngularVelFromPoses(const Frame& frame_1, const Frame& frame_2, Vector& ang_vel){
    // Getting rotation matrices from frames
    Rotation r_1 = frame_1.M;
    Rotation r_2 = frame_2.M;

    // Computing R(1)R(0)'
    Rotation R_Rtranspose = r_2*r_1.Inverse();

    // Getting axis and angle of R(1)R(0)'
    Vector axis = R_Rtranspose.GetRot();        // rotation axis
    double angle = R_Rtranspose.GetRotAngle(axis);     // redundant: both angle and rotation axis

    // Computing skew matrix of axis
    Eigen::Matrix<double,3,3> skew_axis;
    skew_symmetric(axis, skew_axis);

    // Computing skew angular velocity
    Eigen::Matrix<double,3,3> skew_ang = skew_axis*angle;

    // Extraction angular velocity from previous skew matrix
    Vector ang_vel_tmp;
    ang_vel_tmp(0) = skew_ang(2, 1);
    ang_vel_tmp(1) = skew_ang(0, 2);
    ang_vel_tmp(2) = skew_ang(1, 0);

    ang_vel = ang_vel_tmp;
}

int ChainIkSolverPos_relaxed::CartToJnt(const JntArray& q_init, const Frame& p_in, JntArray& q_out)
{
    if(nj != q_init.rows() || nj != q_out.rows())
        return (error = E_SIZE_MISMATCH);
    
    q_out = q_init;
    
    unsigned int i;
    Twist delta_twist_chk;
    for(i=0;i<maxiter;i++)
    {
        if ( fksolver.JntToCart(q_out,f) < 0)
            return (error = E_FKSOLVERPOS_FAILED);
        
        delta_twist = diff(f,p_in);

        // WARNING: diff does not give a twist -> check KDL Reference for more info.
        if(f.M != p_in.M){
            #if DEBUG>1
            std::cout << "I ENTERED ANGULAR VELOCITY COMPUTATION AND CHANGE!" << std::endl;
            #endif
            Vector angular_vel;
            getAngularVelFromPoses(f, p_in, angular_vel);
            delta_twist.rot = angular_vel;
        }

        
        /// apply weighting to the task for checking the completion
        if(use_ee_task_)
        {
            delta_twist_chk = p_in.M.Inverse()*delta_twist;
        }
        else
        {
            delta_twist_chk = delta_twist;
        }
        for(int j=0; j<6; ++j)
            delta_twist_chk[j] *= W(j);
        
#if DEBUG>1
        std::cout << "delta_twist = " << delta_twist << std::endl;
        std::cout << "|delta_twist| = " << sqrt(delta_twist.rot.Norm()*delta_twist.rot.Norm() + delta_twist.vel.Norm()*delta_twist.vel.Norm()) << std::endl;
        std::cout << "delta_twist_chk = " << delta_twist_chk << std::endl;
        std::cout << "|delta_twist_chk| = " << sqrt(delta_twist_chk.rot.Norm()*delta_twist_chk.rot.Norm() + delta_twist_chk.vel.Norm()*delta_twist_chk.vel.Norm()) << std::endl;
#endif
        
        if(Equal(delta_twist_chk,Twist::Zero(),eps))
            break;
        
        int err;
        KDL::SetToZero(delta_q);
        if ( (err=iksolver.CartToJnt(q_out,delta_twist,delta_q)) < 0)
        {
#if DEBUG>1
            std::cout << "Error #" << err << std::endl;
#endif
            return (error = E_IKSOLVERVEL_FAILED);
        }
        JntArray tmp();
#if DEBUG>1
        std::cout << "IKpos Solver delta_q = [" << delta_q.data.transpose() << "]" << std::endl;
#endif
//         assert(!delta_q.data.isNaN());
        for(int i=0; i<delta_q.rows(); ++i)
            assert(!std::isnan(delta_q.data(i)));
        if(Equal(delta_q,JntArray(7),1e-15))
        {
#if DEBUG>1
            std::cout << __func__ << " : exiting in cycle #" << i << " (out of " << maxiter << ") because |q_dot|<1e-15 - error on pose = " << delta_twist << std::endl;
#endif
            return (error = E_NO_CONVERGE);
        }
        
        Add(q_out,delta_q,q_out);
        
        /// set limits
        for(unsigned int j=0; j<q_min.rows(); j++) {
            if(q_out(j) < q_min(j))
                q_out(j) = q_min(j);
        }
        for(unsigned int j=0; j<q_max.rows(); j++) {
            if(q_out(j) > q_max(j))
                q_out(j) = q_max(j);
        }
    }
    
    if(i!=maxiter)
        return (error = E_NOERROR);
    else
        return (error = E_MAX_ITERATIONS_EXCEEDED);
}

int ChainIkSolverPos_relaxed::setJointLimits(const JntArray& q_min_in, const JntArray& q_max_in) {
    if (q_min_in.rows() != nj || q_max_in.rows() != nj)
        return (error = E_SIZE_MISMATCH);
    q_min = q_min_in;
    q_max = q_max_in;
    return (error = E_NOERROR);
}

ChainIkSolverPos_relaxed::~ChainIkSolverPos_relaxed()
{
}

const char* ChainIkSolverPos_relaxed::strError(const int error) const
{
    if (E_FKSOLVERPOS_FAILED == error) return "Internal forward position solver failed.";
    else if (E_IKSOLVERVEL_FAILED == error) return "Internal inverse velocity solver failed.";
// ATTENTION remove these two when KDL version >= 1.4
    else if (E_MAX_ITERATIONS_EXCEEDED == error) return "Maximum number of iterations exceeded";
    else if (E_SIZE_MISMATCH == error) return "Size mismatch";
    else return SolverI::strError(error);
}

void ChainIkSolverPos_relaxed::setTaskWeight(const Eigen::Matrix< double, 6, 1 >& W_)
{
    W = W_;
    W.array().square();
}
