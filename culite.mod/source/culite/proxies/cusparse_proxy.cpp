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
} // namespace cusparse
} // namespace culite
/*-------------------------------------------------*/
