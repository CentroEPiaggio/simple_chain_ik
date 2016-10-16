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

using namespace KDL;

ChainIkSolverPos_relaxed::ChainIkSolverPos_relaxed(const Chain& chain, const JntArray& q_min, const JntArray& q_max, ChainFkSolverPos& fksolver, ChainIkSolverVel& iksolver, unsigned int maxiter, double eps) : chain(chain), nj(chain.getNrOfJoints()), q_min(q_min), q_max(q_max), iksolver(iksolver), fksolver(fksolver), delta_q(chain.getNrOfJoints()), maxiter(maxiter), eps(eps), W(Eigen::MatrixXd::Ones(W.RowsAtCompileTime,W.ColsAtCompileTime))
{
}

ChainIkSolverPos_relaxed::ChainIkSolverPos_relaxed(const Chain& chain, ChainFkSolverPos& fksolver, ChainIkSolverVel& iksolver, unsigned int maxiter, double eps) : chain(chain), nj(chain.getNrOfJoints()), q_min(nj), q_max(nj), iksolver(iksolver), fksolver(fksolver), delta_q(chain.getNrOfJoints()), maxiter(maxiter), eps(eps), W(Eigen::MatrixXd::Ones(W.RowsAtCompileTime,W.ColsAtCompileTime))
{
    q_min.data.setConstant(std::numeric_limits<double>::min());
    q_max.data.setConstant(std::numeric_limits<double>::max());
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
        
        /// apply weighting to the task for checking the completion
        delta_twist_chk = delta_twist;
        for(int j=0; j<6; ++j)
            delta_twist_chk[j] *= W(j);
        if(Equal(delta_twist_chk,Twist::Zero(),eps))
            break;
        
        std::cout << "BLA1 : try" << std::endl;
        int err;
        if ( (err=iksolver.CartToJnt(q_out,delta_twist,delta_q)) < 0)
            std::cout << "Error #" << err << std::endl;
            return (error = E_IKSOLVERVEL_FAILED);
        std::cout << "BLA2 : catch" << std::endl;
        
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
