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
#include <culite/algebra/functional_add.hpp>

// system

// 3rd

// culite
#include <culite/error/exceptions.hpp>

// forwards
#include <culite/types/cla3p_forwards.hpp>
#include <culite/checks/cla3p_forwards.hpp>

/*-------------------------------------------------*/
namespace culite {
namespace ops {
/*-------------------------------------------------*/
template <typename T_Scalar>
void add(T_Scalar alpha, const dns::XxVector<T_Scalar>& x,
         T_Scalar beta , const dns::XxVector<T_Scalar>& y,
         dns::XxVector<T_Scalar>& z,
         CuBlasHandler& cublasHandler)
{
    ::cla3p::similarity_dim_check(x.size(), y.size());
    ::cla3p::similarity_dim_check(x.size(), z.size());
    ::cla3p::similarity_dim_check(y.size(), z.size());

    cublasHandler.geam(op_t::N,
                       op_t::N,
                       z.size(), 1,
                       &alpha, x.values(), x.size(),
                       &beta,  y.values(), y.size(),
                       z.values(), z.size());
}
/*-------------------------------------------------*/
#define instantiate_add(T_Scl) \
template void add( \
        T_Scl, const dns::XxVector<T_Scl>&, \
        T_Scl, const dns::XxVector<T_Scl>&, \
        dns::XxVector<T_Scl>&, CuBlasHandler&)
instantiate_add(real_t);
instantiate_add(real4_t);
instantiate_add(complex_t);
instantiate_add(complex8_t);
#undef instantiate_add
/*-------------------------------------------------*/
template <typename T_Scalar>
void add(op_t opA, T_Scalar alpha, const dns::XxMatrix<T_Scalar>& A,
         op_t opB, T_Scalar beta , const dns::XxMatrix<T_Scalar>& B,
         dns::XxMatrix<T_Scalar>& C,
         CuBlasHandler& cublasHandler)
{
    op_t opC = op_t::N;

    ::cla3p::similarity_check(opA, A, opB, B);
    ::cla3p::similarity_check(opA, A, opC, C);
    ::cla3p::similarity_check(opB, B, opC, C);

    if(A.prop().isGeneral() && B.prop().isGeneral() && C.prop().isGeneral()) {

        cublasHandler.geam(opA,
                           opB,
                           C.nrows(), C.ncols(),
                           &alpha, A.values(), A.ld(),
                           &beta,  B.values(), B.ld(),
                           C.values(), C.ld());

    } else {

        throw err::CudaException("Unsupported matrix property for addition. Only General (dense) matrices are currently supported.");

    } // prop
}
/*-------------------------------------------------*/
#define instantiate_add(T_Scl) \
template void add( \
        op_t, T_Scl, const dns::XxMatrix<T_Scl>&, \
        op_t, T_Scl, const dns::XxMatrix<T_Scl>&, \
        dns::XxMatrix<T_Scl>&, CuBlasHandler&)
instantiate_add(real_t);
instantiate_add(real4_t);
instantiate_add(complex_t);
instantiate_add(complex8_t);
#undef instantiate_add
/*-------------------------------------------------*/
} // namespace ops
} // namespace culite
/*-------------------------------------------------*/
