#include <simple_chain_ik/solvers/my_eigen_indexing.h>
#include <iostream>

using namespace Eigen;

int main()
{
    Eigen::MatrixXi A = Eigen::MatrixXi::Random(4,4);
    std::cout << "A =" << std::endl;
    std::cout << A << std::endl << std::endl;
    Array3i ri(1,2,1);
    ArrayXi ci(6); ci << 3,2,1,0,0,2;
    
    indexing_functor<decltype(A),Array3i,ArrayXi> indf(A,ri,ci);
    
    std::cout << "indf(2,3) [aka A(1,0)] = " << indf(2,3) << std::endl;
    
//     Eigen::MatrixXi B = indexing<decltype(A),Array3i,ArrayXi>(A, ri, ci);
//     std::cout << "A([" << ri.transpose() << "], [" << ci.transpose() << "]) =" << std::endl;
//     std::cout << B << std::endl;
//     
//     B =  indexing(A, ri+1, ci);
//     std::cout << "A(ri+1,ci) =" << std::endl;
//     std::cout << B << std::endl << std::endl;
// #if __cplusplus >= 201103L
//     B =  indexing(A, ArrayXi::LinSpaced(13,0,12).unaryExpr([](int x){return x%4;}), ArrayXi::LinSpaced(4,0,3));
//     std::cout << "A(ArrayXi::LinSpaced(13,0,12).unaryExpr([](int x){return x%4;}), ArrayXi::LinSpaced(4,0,3)) =" << std::endl;
//     std::cout << B << std::endl << std::endl;
// #endif
    
    return 0;
}

