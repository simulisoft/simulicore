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
#include "culite/algebra/functional_outer.hpp"

// system

// 3rd

// culite
#include "culite/types/scalar.hpp"
#include "culite/error/exceptions.hpp"
#include "culite/dense/dns_xxvector.hpp"
#include "culite/dense/dns_xxmatrix.hpp"

// forwards
#include "culite/checks/cla3p_forwards.hpp"

/*-------------------------------------------------*/
namespace culite {
namespace ops {
/*-------------------------------------------------*/
template <typename T_Scalar>
static void outerx(bool conjop, T_Scalar alpha, 
                   const dns::XxVector<T_Scalar>& x, 
                   const dns::XxVector<T_Scalar>& y, 
                   dns::XxMatrix<T_Scalar>& A,
                   CuBlasHandler& cublasHandler)
{
    conjop = (TypeTraits<T_Scalar>::is_real() ? false : conjop);

    ::cla3p::outer_product_consistency_check(conjop, A, x, y);
    hermitian_coeff_check2<T_Scalar>(A.prop(), alpha);

    if(A.prop().isGeneral()) {

        if(conjop) {
            cublasHandler.gerc(x.size(), y.size(), &alpha, x.values(), 1, y.values(), 1, A.values(), A.ld());
        } else {
            cublasHandler.ger(x.size(), y.size(), &alpha, x.values(), 1, y.values(), 1, A.values(), A.ld());
        }

    } else if (A.prop().isSymmetric()) {

        if(x.values() == y.values()) {

            cublasHandler.syr(A.prop().uplo(),
                              A.ncols(),
                              &alpha,
                              x.values(), 1,
                              A.values(), A.ld());
        } else {

            T_Scalar beta = makeScalar<T_Scalar>(1);
            cublasHandler.syrkx(A.prop().uplo(),
                                op_t::N,
                                A.ncols(),
                                1,
                                &alpha,
                                x.values(), x.size(),
                                y.values(), y.size(),
                                &beta,
                                A.values(), A.ld());

        }

    } else if (A.prop().isHermitian()) {

        using T_RScalar = typename TypeTraits<T_Scalar>::real_type;
        T_RScalar realAlpha = arith::getRe(alpha);
        T_RScalar realBeta = makeScalar<T_RScalar>(1);

        if(x.values() == y.values()) {

            cublasHandler.her(A.prop().uplo(),
                              A.ncols(),
                              &realAlpha,
                              x.values(), 1,
                              A.values(), A.ld());
        } else {

            cublasHandler.herkx(A.prop().uplo(),
                                op_t::N,
                                A.ncols(),
                                1,
                                &alpha,
                                x.values(), x.size(),
                                y.values(), y.size(),
                                &realBeta,
                                A.values(), A.ld());
        }

    } else {

        throw err::CudaException("Unsupported matrix property for outer product update (" + A.prop().name() + ").");

    } // valid props
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void outer(T_Scalar alpha, 
           const dns::XxVector<T_Scalar>& x, 
           const dns::XxVector<T_Scalar>& y, 
           dns::XxMatrix<T_Scalar>& A,
           CuBlasHandler& cublasHandler)
{
    outerx(false, alpha, x, y, A, cublasHandler);
}
/*-------------------------------------------------*/
#define instantiate_outer(T_Scl) \
template void outer(T_Scl, const dns::XxVector<T_Scl>&, const dns::XxVector<T_Scl>&, dns::XxMatrix<T_Scl>&, CuBlasHandler&)
instantiate_outer(real_t);
instantiate_outer(real4_t);
instantiate_outer(complex_t);
instantiate_outer(complex8_t);
#undef instantiate_outer
/*-------------------------------------------------*/
template <typename T_Scalar>
void outerc(T_Scalar alpha, 
            const dns::XxVector<T_Scalar>& x, 
            const dns::XxVector<T_Scalar>& y, 
            dns::XxMatrix<T_Scalar>& A,
            CuBlasHandler& cublasHandler)
{
    outerx(true, alpha, x, y, A, cublasHandler);
}
/*-------------------------------------------------*/
#define instantiate_outerc(T_Scl) \
template void outerc(T_Scl, const dns::XxVector<T_Scl>&, const dns::XxVector<T_Scl>&, dns::XxMatrix<T_Scl>&, CuBlasHandler&)
instantiate_outerc(real_t);
instantiate_outerc(real4_t);
instantiate_outerc(complex_t);
instantiate_outerc(complex8_t);
#undef instantiate_outerc
/*-------------------------------------------------*/
} // namespace ops
} // namespace culite
/*-------------------------------------------------*/
