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

#include <culite/types/integer.hpp>
#include <culite/types/scalar.hpp>
#include <culite/types/enums.hpp>

/*-------------------------------------------------*/
namespace culite {
namespace cusparse {
/*-------------------------------------------------*/

using cuSparseInt = int64_t;

/*-------------------------------------------------*/

op_t cusparseOp2cla3pOp(cusparseOperation_t op);
cusparseOperation_t cla3pOp2cusparseOp(op_t op);

uplo_t cusparseUplo2cla3pUplo(cusparseFillMode_t uplo);
cusparseFillMode_t cla3pUplo2cusparseUplo(uplo_t uplo);

diag_t cusparseDiag2cla3pDiag(cusparseDiagType_t diag);
cusparseDiagType_t cla3pDiag2cusparseDiag(diag_t diag);

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

        // disable copying
        DnVec(const DnVec&) = delete;
        DnVec& operator=(const DnVec&) = delete;

        DnVec(cuSparseInt size, T_Scalar *vals);
        DnVec(cuSparseInt size, const T_Scalar *vals);
        ~DnVec();

        cusparseDnVecDescr_t descr() { return m_descr; }
        cusparseConstDnVecDescr_t descr() const { return m_descr; }

    private:
        cusparseDnVecDescr_t m_descr{nullptr};
};

/*-------------------------------------------------*/

template <typename T_Scalar>
class DnMat {

    public:

        // disable copying
        DnMat(const DnMat&) = delete;
        DnMat& operator=(const DnMat&) = delete;

        DnMat(cuSparseInt     rows, 
              cuSparseInt     cols, 
              T_Scalar*       vals, 
              cuSparseInt     ld, 
              cusparseOrder_t order = cusparseOrder_t::CUSPARSE_ORDER_COL);
        
        DnMat(cuSparseInt     rows, 
              cuSparseInt     cols, 
              const T_Scalar* vals, 
              cuSparseInt     ld, 
              cusparseOrder_t order = cusparseOrder_t::CUSPARSE_ORDER_COL);
        
        ~DnMat();

        cusparseDnMatDescr_t descr() { return m_descr; }
        cusparseConstDnMatDescr_t descr() const { return m_descr; }

    private:
        cusparseDnMatDescr_t m_descr{nullptr};
};

/*-------------------------------------------------*/

class SpMatBase {

    public:

        // disable copying
        SpMatBase(const SpMatBase&) = delete;
        SpMatBase& operator=(const SpMatBase&) = delete;

        SpMatBase() = default;
        ~SpMatBase();

        cusparseSpMatDescr_t descr() { return m_descr; };
        cusparseConstSpMatDescr_t descr() const { return m_descr; };

    protected:
        cusparseSpMatDescr_t m_descr{nullptr};
};

/*-------------------------------------------------*/

template <typename T_Scalar>
class SpMatCsr : public SpMatBase {

    public:

        // disable copying
        SpMatCsr(const SpMatCsr&) = delete;
        SpMatCsr& operator=(const SpMatCsr&) = delete;

        SpMatCsr(cuSparseInt         rows,
                 cuSparseInt         cols,
                 cuSparseInt         nnz,
                 int_t*              rowptr,
                 int_t*              colidx,
                 T_Scalar*           values,
                 cusparseIndexBase_t idxBase = cusparseIndexBase_t::CUSPARSE_INDEX_BASE_ZERO);

        SpMatCsr(cuSparseInt         rows,
                 cuSparseInt         cols,
                 cuSparseInt         nnz,
                 const int_t*        rowptr,
                 const int_t*        colidx,
                 const T_Scalar*     values,
                 cusparseIndexBase_t idxBase = cusparseIndexBase_t::CUSPARSE_INDEX_BASE_ZERO);

        ~SpMatCsr() = default;
};

/*-------------------------------------------------*/

template <typename T_Scalar>
class SpMatCsc : public SpMatBase {

    public:

        // disable copying
        SpMatCsc(const SpMatCsc&) = delete;
        SpMatCsc& operator=(const SpMatCsc&) = delete;

        SpMatCsc(cuSparseInt         rows,
                 cuSparseInt         cols,
                 cuSparseInt         nnz,
                 int_t*              colptr,
                 int_t*              rowidx,
                 T_Scalar*           values,
                 cusparseIndexBase_t idxBase = cusparseIndexBase_t::CUSPARSE_INDEX_BASE_ZERO);

        SpMatCsc(cuSparseInt         rows,
                 cuSparseInt         cols,
                 cuSparseInt         nnz,
                 const int_t*        colptr,
                 const int_t*        rowidx,
                 const T_Scalar*     values,
                 cusparseIndexBase_t idxBase = cusparseIndexBase_t::CUSPARSE_INDEX_BASE_ZERO);

        ~SpMatCsc() = default;
};

/*-------------------------------------------------*/
} // namespace cusparse
} // namespace culite
/*-------------------------------------------------*/

#endif // CULITE_CUSPARSE_PROXY_HPP_
