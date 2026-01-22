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
#include "culite/proxies/cusparse_proxy.hpp"

// system

// 3rd

// culite
#include "culite/error/exceptions.hpp"

/*-------------------------------------------------*/
namespace culite {
namespace cusparse {
/*-------------------------------------------------*/
::cla3p::op_t cusparseOp2cla3pOp(cusparseOperation_t op)
{
    switch(op) {
        case CUSPARSE_OPERATION_NON_TRANSPOSE: return ::cla3p::op_t::N;
        case CUSPARSE_OPERATION_TRANSPOSE:     return ::cla3p::op_t::T;
        case CUSPARSE_OPERATION_CONJUGATE_TRANSPOSE: return ::cla3p::op_t::C;
        default:
            throw err::CudaException("Invalid cusparseOperation_t value.");
    }
}
/*-------------------------------------------------*/
cusparseOperation_t cla3pOp2cusparseOp(::cla3p::op_t op)
{
    switch(op) {
        case ::cla3p::op_t::N: return CUSPARSE_OPERATION_NON_TRANSPOSE;
        case ::cla3p::op_t::T: return CUSPARSE_OPERATION_TRANSPOSE;
        case ::cla3p::op_t::C: return CUSPARSE_OPERATION_CONJUGATE_TRANSPOSE;
        default:
            throw err::CudaException("Invalid (::cla3p::op_t) value.");
    }
}
/*-------------------------------------------------*/
::cla3p::uplo_t cusparseUplo2cla3pUplo(cusparseFillMode_t uplo)
{
    switch(uplo) {
        case CUSPARSE_FILL_MODE_LOWER: return ::cla3p::uplo_t::Lower;
        case CUSPARSE_FILL_MODE_UPPER: return ::cla3p::uplo_t::Upper;
        default:
            throw err::CudaException("Invalid cusparseFillMode_t value.");
    }
}
/*-------------------------------------------------*/
cusparseFillMode_t cla3pUplo2cusparseUplo(::cla3p::uplo_t uplo)
{
    switch(uplo) {
        case ::cla3p::uplo_t::Lower: return CUSPARSE_FILL_MODE_LOWER;
        case ::cla3p::uplo_t::Upper: return CUSPARSE_FILL_MODE_UPPER;
        default:
            throw err::CudaException("Invalid (::cla3p::uplo_t) value.");
    }
}
/*-------------------------------------------------*/
::cla3p::diag_t cusparseDiag2cla3pDiag(cusparseDiagType_t diag)
{
    switch(diag) {
        case CUSPARSE_DIAG_TYPE_NON_UNIT: return ::cla3p::diag_t::NonUnit;
        case CUSPARSE_DIAG_TYPE_UNIT:     return ::cla3p::diag_t::Unit;
        default:
            throw err::CudaException("Invalid cusparseDiagType_t value.");
    }
}
/*-------------------------------------------------*/
cusparseDiagType_t cla3pDiag2cusparseDiag(::cla3p::diag_t diag)
{
    switch(diag) {
        case ::cla3p::diag_t::NonUnit: return CUSPARSE_DIAG_TYPE_NON_UNIT;
        case ::cla3p::diag_t::Unit:    return CUSPARSE_DIAG_TYPE_UNIT;
        default:
            throw err::CudaException("Invalid (::cla3p::diag_t) value.");
    }
}
/*-------------------------------------------------*/
} // namespace cusparse
} // namespace culite
/*-------------------------------------------------*/
