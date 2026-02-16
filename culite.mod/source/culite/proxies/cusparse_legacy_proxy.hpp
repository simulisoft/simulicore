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

#ifndef CULITE_CUSPARSE_LEGACY_PROXY_HPP_
#define CULITE_CUSPARSE_LEGACY_PROXY_HPP_

/**
 * @file
 */

#include <cusparse.h>

#include "culite/types/integer.hpp"
#include "culite/types/scalar.hpp"
#include "culite/error/cuda.hpp"

/*-------------------------------------------------*/
namespace culite {
namespace cusparse {
/*-------------------------------------------------*/

#if !defined (CULITE_I64)

class LegacyMatDescr {

    public:

        // disable copying
        LegacyMatDescr(const LegacyMatDescr&) = delete;
        LegacyMatDescr& operator=(const LegacyMatDescr&) = delete;

        LegacyMatDescr()
        {
            cusparseStatus_t cusparseStatus = cusparseCreateMatDescr(&m_descr);
            err::check_cusparse(cusparseStatus);
        }

        ~LegacyMatDescr()
        {
            cusparseStatus_t cusparseStatus = cusparseDestroyMatDescr(m_descr);
            err::check_cusparse(cusparseStatus);
        }

        cusparseMatDescr_t descr() const { return m_descr; }

    private:
        cusparseMatDescr_t m_descr{nullptr};
};

//
// At the moment only 32bit integers are supported by legacy cuSparse
//
#define csr_geam2_reserve_macro(typein) \
void csr_geam2_reserve(cusparseHandle_t         handle, \
                       int_t                    m, \
                       int_t                    n, \
                       const typein*            alpha, \
                       const cusparseMatDescr_t descrA, \
                       int_t                    nnzA, \
                       const typein*            csrSortedValA, \
                       const int_t*             csrSortedRowPtrA, \
                       const int_t*             csrSortedColIndA, \
                       const typein*            beta, \
                       const cusparseMatDescr_t descrB, \
                       int_t                    nnzB, \
                       const typein*            csrSortedValB, \
                       const int_t*             csrSortedRowPtrB, \
                       const int_t*             csrSortedColIndB, \
                       const cusparseMatDescr_t descrC, \
                       const typein*            csrSortedValC, \
                       const int_t*             csrSortedRowPtrC, \
                       const int_t*             csrSortedColIndC, \
                       size_t*                  pBufferSizeInBytes)
csr_geam2_reserve_macro(real_t);
csr_geam2_reserve_macro(real4_t);
csr_geam2_reserve_macro(complex_t);
csr_geam2_reserve_macro(complex8_t);
#undef csr_geam2_reserve_macro

void csr_geam2_nnz(cusparseHandle_t         handle,
                   int_t                    m,
                   int_t                    n,
                   const cusparseMatDescr_t descrA,
                   int_t                    nnzA,
                   const int_t*             csrSortedRowPtrA,
                   const int_t*             csrSortedColIndA,
                   const cusparseMatDescr_t descrB,
                   int_t                    nnzB,
                   const int_t*             csrSortedRowPtrB,
                   const int_t*             csrSortedColIndB,
                   const cusparseMatDescr_t descrC,
                   int_t*                   csrSortedRowPtrC,
                   int_t*                   nnzTotalDevHostPtr,
                   void*                    workspace);

#define csr_geam2_macro(typein) \
void csr_geam2(cusparseHandle_t            handle, \
                  int_t                    m, \
                  int_t                    n, \
                  const typein*            alpha, \
                  const cusparseMatDescr_t descrA, \
                  int_t                    nnzA, \
                  const typein*            csrSortedValA, \
                  const int_t*             csrSortedRowPtrA, \
                  const int_t*             csrSortedColIndA, \
                  const typein*            beta, \
                  const cusparseMatDescr_t descrB, \
                  int_t                    nnzB, \
                  const typein*            csrSortedValB, \
                  const int_t*             csrSortedRowPtrB, \
                  const int_t*             csrSortedColIndB, \
                  const cusparseMatDescr_t descrC, \
                  typein*                  csrSortedValC, \
                  int_t*                   csrSortedRowPtrC, \
                  int_t*                   csrSortedColIndC, \
                  void*                    pBuffer)
csr_geam2_macro(real_t);
csr_geam2_macro(real4_t);
csr_geam2_macro(complex_t);
csr_geam2_macro(complex8_t);
#undef csr_geam2_macro

#endif // CULITE_I64

/*-------------------------------------------------*/
} // namespace cusparse
} // namespace culite
/*-------------------------------------------------*/

#endif // CULITE_CUSPARSE_LEGACY_PROXY_HPP_
