/*
 * Copyright (c) 2025-2026 Simulisoft
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *     http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

// this file inc
#include <culite/proxies/cusparse_proxy.hpp>

// system

// 3rd

// culite
#include <culite/error/cuda.hpp>
#include <culite/error/exceptions.hpp>

/*-------------------------------------------------*/
namespace culite {
namespace cusparse {
/*-------------------------------------------------*/
op_t cusparseOp2cla3pOp(cusparseOperation_t op)
{
    switch(op) {
        case CUSPARSE_OPERATION_NON_TRANSPOSE: return op_t::N;
        case CUSPARSE_OPERATION_TRANSPOSE:     return op_t::T;
        case CUSPARSE_OPERATION_CONJUGATE_TRANSPOSE: return op_t::C;
        default:
            throw err::CudaException("Invalid cusparseOperation_t value.");
    }
}
/*-------------------------------------------------*/
cusparseOperation_t cla3pOp2cusparseOp(op_t op)
{
    switch(op) {
        case op_t::N: return CUSPARSE_OPERATION_NON_TRANSPOSE;
        case op_t::T: return CUSPARSE_OPERATION_TRANSPOSE;
        case op_t::C: return CUSPARSE_OPERATION_CONJUGATE_TRANSPOSE;
        default:
            throw err::CudaException("Invalid (op_t) value.");
    }
}
/*-------------------------------------------------*/
uplo_t cusparseUplo2cla3pUplo(cusparseFillMode_t uplo)
{
    switch(uplo) {
        case CUSPARSE_FILL_MODE_LOWER: return uplo_t::Lower;
        case CUSPARSE_FILL_MODE_UPPER: return uplo_t::Upper;
        default:
            throw err::CudaException("Invalid cusparseFillMode_t value.");
    }
}
/*-------------------------------------------------*/
cusparseFillMode_t cla3pUplo2cusparseUplo(uplo_t uplo)
{
    switch(uplo) {
        case uplo_t::Lower: return CUSPARSE_FILL_MODE_LOWER;
        case uplo_t::Upper: return CUSPARSE_FILL_MODE_UPPER;
        default:
            throw err::CudaException("Invalid (uplo_t) value.");
    }
}
/*-------------------------------------------------*/
diag_t cusparseDiag2cla3pDiag(cusparseDiagType_t diag)
{
    switch(diag) {
        case CUSPARSE_DIAG_TYPE_NON_UNIT: return diag_t::NonUnit;
        case CUSPARSE_DIAG_TYPE_UNIT:     return diag_t::Unit;
        default:
            throw err::CudaException("Invalid cusparseDiagType_t value.");
    }
}
/*-------------------------------------------------*/
cusparseDiagType_t cla3pDiag2cusparseDiag(diag_t diag)
{
    switch(diag) {
        case diag_t::NonUnit: return CUSPARSE_DIAG_TYPE_NON_UNIT;
        case diag_t::Unit:    return CUSPARSE_DIAG_TYPE_UNIT;
        default:
            throw err::CudaException("Invalid (diag_t) value.");
    }
}
/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/
template <typename T_Scalar>
DnVec<T_Scalar>::DnVec(cuSparseInt size, T_Scalar *vals)
{
    cusparseStatus_t cusparseStatus = 
    cusparseCreateDnVec(&m_descr, 
                        size, 
                        vals, 
                        TypeTraits<T_Scalar>::cuda_type());
    err::check_cusparse(cusparseStatus);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
DnVec<T_Scalar>::DnVec(cuSparseInt size, const T_Scalar *vals)
{
    cusparseStatus_t cusparseStatus = 
    cusparseCreateConstDnVec(const_cast<cusparseConstDnVecDescr_t*>(&m_descr), 
                             size, 
                             vals, 
                             TypeTraits<T_Scalar>::cuda_type());
    err::check_cusparse(cusparseStatus);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
DnVec<T_Scalar>::~DnVec()
{
    cusparseStatus_t cusparseStatus = cusparseDestroyDnVec(m_descr);
    err::check_cusparse(cusparseStatus);
}
/*-------------------------------------------------*/
template class DnVec<real_t>;
template class DnVec<real4_t>;
template class DnVec<complex_t>;
template class DnVec<complex8_t>;
/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/
template <typename T_Scalar>
DnMat<T_Scalar>::DnMat(cuSparseInt rows, cuSparseInt cols, T_Scalar *vals, cuSparseInt ld, 
                       cusparseOrder_t order)
{
    cusparseStatus_t cusparseStatus = 
    cusparseCreateDnMat(&m_descr, 
                        rows, 
                        cols, 
                        ld, 
                        vals, 
                        TypeTraits<T_Scalar>::cuda_type(), 
                        order);
    err::check_cusparse(cusparseStatus);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
DnMat<T_Scalar>::DnMat(cuSparseInt rows, cuSparseInt cols, const T_Scalar *vals, cuSparseInt ld, 
                       cusparseOrder_t order)
{
    cusparseStatus_t cusparseStatus = 
    cusparseCreateConstDnMat(const_cast<cusparseConstDnMatDescr_t*>(&m_descr), 
                             rows, 
                             cols, 
                             ld, 
                             vals, 
                             TypeTraits<T_Scalar>::cuda_type(), 
                             order);
    err::check_cusparse(cusparseStatus);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
DnMat<T_Scalar>::~DnMat()
{
    cusparseStatus_t cusparseStatus = cusparseDestroyDnMat(m_descr);
    err::check_cusparse(cusparseStatus);
}
/*-------------------------------------------------*/
template class DnMat<real_t>;
template class DnMat<real4_t>;
template class DnMat<complex_t>;
template class DnMat<complex8_t>;
/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/
SpMatBase::~SpMatBase()
{
    if(m_descr != nullptr) {
        cusparseStatus_t cusparseStatus = cusparseDestroySpMat(m_descr);
        err::check_cusparse(cusparseStatus);
    }
}
/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/
template <typename T_Scalar>
SpMatCsr<T_Scalar>::SpMatCsr(cuSparseInt         rows,
                              cuSparseInt         cols,
                              cuSparseInt         nnz,
                              int_t*              rowptr,
                              int_t*              colidx,
                              T_Scalar*           values,
                              cusparseIndexBase_t idxBase)
{
    cusparseStatus_t cusparseStatus = 
    cusparseCreateCsr(&m_descr,
                      rows,
                      cols,
                      nnz,
                      rowptr,
                      colidx,
                      values,
                      indexType(),
                      indexType(),
                      idxBase,
                      TypeTraits<T_Scalar>::cuda_type());
    err::check_cusparse(cusparseStatus);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
SpMatCsr<T_Scalar>::SpMatCsr(cuSparseInt         rows,
                              cuSparseInt         cols,
                              cuSparseInt         nnz,
                              const int_t*        rowptr,
                              const int_t*        colidx,
                              const T_Scalar*     values,
                              cusparseIndexBase_t idxBase)
{
    cusparseStatus_t cusparseStatus = 
    cusparseCreateConstCsr(const_cast<cusparseConstSpMatDescr_t*>(&m_descr),
                           rows,
                           cols,
                           nnz,
                           rowptr,
                           colidx,
                           values,
                           indexType(),
                           indexType(),
                           idxBase,
                           TypeTraits<T_Scalar>::cuda_type());
    err::check_cusparse(cusparseStatus);
}
/*-------------------------------------------------*/
template class SpMatCsr<real_t>;
template class SpMatCsr<real4_t>;
template class SpMatCsr<complex_t>;
template class SpMatCsr<complex8_t>;
/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/
template <typename T_Scalar>
SpMatCsc<T_Scalar>::SpMatCsc(cuSparseInt         rows,
                              cuSparseInt         cols,
                              cuSparseInt         nnz,
                              int_t*              colptr,
                              int_t*              rowidx,
                              T_Scalar*           values,
                              cusparseIndexBase_t idxBase)
{
    cusparseStatus_t cusparseStatus = 
    cusparseCreateCsc(&m_descr,
                      rows,
                      cols,
                      nnz,
                      colptr,
                      rowidx,
                      values,
                      indexType(),
                      indexType(),
                      idxBase,
                      TypeTraits<T_Scalar>::cuda_type());
    err::check_cusparse(cusparseStatus);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
SpMatCsc<T_Scalar>::SpMatCsc(cuSparseInt         rows,
                              cuSparseInt         cols,
                              cuSparseInt         nnz,
                              const int_t*        colptr,
                              const int_t*        rowidx,
                              const T_Scalar*     values,
                              cusparseIndexBase_t idxBase)
{
    cusparseStatus_t cusparseStatus = 
    cusparseCreateConstCsc(const_cast<cusparseConstSpMatDescr_t*>(&m_descr),
                           rows,
                           cols,
                           nnz,
                           colptr,
                           rowidx,
                           values,
                           indexType(),
                           indexType(),
                           idxBase,
                           TypeTraits<T_Scalar>::cuda_type());
    err::check_cusparse(cusparseStatus);
}
/*-------------------------------------------------*/
template class SpMatCsc<real_t>;
template class SpMatCsc<real4_t>;
template class SpMatCsc<complex_t>;
template class SpMatCsc<complex8_t>;
/*-------------------------------------------------*/
} // namespace cusparse
} // namespace culite
/*-------------------------------------------------*/
