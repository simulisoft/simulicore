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

#ifndef CULITE_CUSPARSE_PROXY_HPP_
#define CULITE_CUSPARSE_PROXY_HPP_

/**
 * @file
 */

 #include <cusparse.h>
 #include <cla3p/types/enums.hpp>

 #include "culite/types/integer.hpp"
 #include "culite/types/scalar.hpp"
 #include "culite/error/cuda.hpp"

/*-------------------------------------------------*/
namespace culite {
namespace cusparse {
/*-------------------------------------------------*/

::cla3p::op_t cusparseOp2cla3pOp(cusparseOperation_t op);
cusparseOperation_t cla3pOp2cusparseOp(::cla3p::op_t op);

::cla3p::uplo_t cusparseUplo2cla3pUplo(cusparseFillMode_t uplo);
cusparseFillMode_t cla3pUplo2cusparseUplo(::cla3p::uplo_t uplo);

::cla3p::diag_t cusparseDiag2cla3pDiag(cusparseDiagType_t diag);
cusparseDiagType_t cla3pDiag2cusparseDiag(::cla3p::diag_t diag);

/*-------------------------------------------------*/

inline constexpr cusparseIndexType_t indexType()
{
#if defined (CULITE_I64)
    return cusparseIndexType_t::CUSPARSE_INDEX_64I; 
#else
    return cusparseIndexType_t::CUSPARSE_INDEX_32I;
#endif
}

/*-------------------------------------------------*/

template <typename T_Scalar>
class DnVec {

    public:
        DnVec(int64_t size, T_Scalar *vals) 
        {
            cusparseStatus_t cusparseStatus = cusparseCreateDnVec(m_descr, size, vals, TypeTraits<T_Scalar>::cuda_type());
            err::check_cusparse(cusparseStatus);
        }

        DnVec(int64_t size, const T_Scalar *vals) 
        {
            cusparseStatus_t cusparseStatus = cusparseCreateConstDnVec(m_descr, size, vals, TypeTraits<T_Scalar>::cuda_type());
            err::check_cusparse(cusparseStatus);
        }

        ~DnVec()
        {
            cusparseStatus_t cusparseStatus = cusparseDestroyDnVec(m_descr);
            err::check_cusparse(cusparseStatus);
        }

        cusparseDnVecDescr_t descr() { return m_descr; }
        cusparseConstDnVecDescr_t descr() const { return m_descr; }

    private:
        cusparseDnVecDescr_t m_descr{nullptr};
};

/*-------------------------------------------------*/

template <typename T_Scalar>
class DnMat {

    public:
        DnMat(int64_t rows, int64_t cols, T_Scalar *vals, int64_t ld, 
              cusparseOrder_t order = cusparseOrder_t::CUSPARSE_ORDER_COL)
        {
            cusparseStatus_t cusparseStatus = cusparseCreateDnMat(&m_descr, rows, cols, ld, vals, TypeTraits<T_Scalar>::cuda_type(), order);
            err::check_cusparse(cusparseStatus);
        }

        DnMat(int64_t rows, int64_t cols, const T_Scalar *vals, int64_t ld, 
              cusparseOrder_t order = cusparseOrder_t::CUSPARSE_ORDER_COL)
        {
            cusparseStatus_t cusparseStatus = cusparseCreateConstDnMat(&m_descr, rows, cols, ld, vals, TypeTraits<T_Scalar>::cuda_type(), order);
            err::check_cusparse(cusparseStatus);
        }

        ~DnMat()
        {
            cusparseStatus_t cusparseStatus = cusparseDestroyDnMat(m_descr);
            err::check_cusparse(cusparseStatus);
        }

        cusparseDnMatDescr_t descr() { return m_descr; }
        cusparseConstDnMatDescr_t descr() const { return m_descr; }

    private:
        cusparseDnMatDescr_t m_descr{nullptr};
};

/*-------------------------------------------------*/

class SpMatBase {

    public:
        SpMatBase() = default;
        
        ~SpMatBase()
        {
            if(m_descr != nullptr) {
                cusparseStatus_t cusparseStatus = cusparseDestroySpMat(m_descr);
                err::check_cusparse(cusparseStatus);
            }
        }

        cusparseSpMatDescr_t descr() { return m_descr; };
        cusparseConstSpMatDescr_t descr() const { return m_descr; };

    protected:
        cusparseSpMatDescr_t m_descr{nullptr};
};

/*-------------------------------------------------*/

template <typename T_Scalar>
class SpMatCsr : public SpMatBase {

    public:
        SpMatCsr(int64_t             rows,
                 int64_t             cols,
                 int64_t             nnz,
                 int_t*              rowptr,
                 int_t*              colidx,
                 T_Scalar*           values,
                 cusparseIndexBase_t idxBase = cusparseIndexBase_t::CUSPARSE_INDEX_BASE_ZERO)
        {
            cusparseStatus_t cusparseStatus = cusparseCreateCsr(&m_descr,
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

        SpMatCsr(int64_t             rows,
                 int64_t             cols,
                 int64_t             nnz,
                 const int_t*        rowptr,
                 const int_t*        colidx,
                 const T_Scalar*     values,
                 cusparseIndexBase_t idxBase = cusparseIndexBase_t::CUSPARSE_INDEX_BASE_ZERO)
        {
            cusparseStatus_t cusparseStatus = cusparseCreateConstCsr(&m_descr,
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

        ~SpMatCsr() = default;
};

/*-------------------------------------------------*/

template <typename T_Scalar>
class SpMatCsc : public SpMatBase {

    public:
        SpMatCsc(int64_t             rows,
                 int64_t             cols,
                 int64_t             nnz,
                 int_t*              colptr,
                 int_t*              rowidx,
                 T_Scalar*           values,
                 cusparseIndexBase_t idxBase = cusparseIndexBase_t::CUSPARSE_INDEX_BASE_ZERO)
        {
            cusparseStatus_t cusparseStatus = cusparseCreateCsc(&m_descr,
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

        SpMatCsc(int64_t             rows,
                 int64_t             cols,
                 int64_t             nnz,
                 const int_t*        colptr,
                 const int_t*        rowidx,
                 const T_Scalar*     values,
                 cusparseIndexBase_t idxBase = cusparseIndexBase_t::CUSPARSE_INDEX_BASE_ZERO)
        {
            cusparseStatus_t cusparseStatus = cusparseCreateConstCsc(&m_descr,
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

        ~SpMatCsc() = default;
};

/*-------------------------------------------------*/
} // namespace cusparse
} // namespace culite
/*-------------------------------------------------*/

#endif // CULITE_CUSPARSE_PROXY_HPP_
