#include <Eigen/Eigen>

// This should go into Eigen 3.4, but there is already an example at
// http://eigen.tuxfamily.org/dox-devel/TopicCustomizing_NullaryExpr.html#title1

namespace Eigen
{

typedef EIGEN_DEFAULT_DENSE_INDEX_TYPE Index;

template<class ArgType, class RowIndexType, class ColIndexType>
class indexing_functor {
    const ArgType &m_arg;
    const RowIndexType &m_rowIndices;
    const ColIndexType &m_colIndices;
public:
//     typedef Matrix<typename ArgType::Scalar,
//     RowIndexType::SizeAtCompileTime,
//     ColIndexType::SizeAtCompileTime,
//     ArgType::Flags&RowMajorBit?RowMajor:ColMajor,
//     RowIndexType::MaxSizeAtCompileTime,
//     ColIndexType::MaxSizeAtCompileTime> MatrixType;
    typedef MatrixXi MatrixType;
    
    indexing_functor(const ArgType& arg, const RowIndexType& row_indices, const ColIndexType& col_indices)
    : m_arg(arg), m_rowIndices(row_indices), m_colIndices(col_indices)
    {};
    
    const typename ArgType::Scalar& operator() (Index row, Index col) const {
        return m_arg(m_rowIndices[row], m_colIndices[col]);
    };
};

template <class ArgType, class RowIndexType, class ColIndexType>
CwiseNullaryOp<indexing_functor<ArgType,RowIndexType,ColIndexType>, typename indexing_functor<ArgType,RowIndexType,ColIndexType>::MatrixType>
indexing(const Eigen::MatrixBase<ArgType>& arg, const RowIndexType& row_indices, const ColIndexType& col_indices)
{
    typedef indexing_functor<ArgType,RowIndexType,ColIndexType> Func;
    typedef typename Func::MatrixType MatrixType;
    return MatrixType::NullaryExpr(row_indices.size(), col_indices.size(), Func(arg.derived(), row_indices, col_indices));
};

}