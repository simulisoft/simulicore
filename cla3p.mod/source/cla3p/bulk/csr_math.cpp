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
#include "cla3p/bulk/csr_math.hpp"

// system

// 3rd

// cla3p
#include "cla3p/error/exceptions.hpp"
#include "cla3p/error/literals.hpp"
#include "cla3p/proxies/blas_proxy.hpp"
#if defined(CLA3P_INTEL_MKL)
#include "cla3p/proxies/mkl_sparse_proxy.hpp"
#elif defined(CLA3P_ARMPL)
#include "cla3p/proxies/armpl_sparse_proxy.hpp"
#endif

/*-------------------------------------------------*/
namespace cla3p {
namespace blk {
namespace csr {
/*-------------------------------------------------*/
template <typename T_Scalar>
void add(int_t m, int_t n,
        T_Scalar alpha, const int_t *rowptrA, const int_t *colidxA, const T_Scalar *valuesA,
        T_Scalar beta, const int_t *rowptrB, const int_t *colidxB, const T_Scalar *valuesB,
        int_t **rowptrC, int_t **colidxC, T_Scalar **valuesC)
{
#if defined(CLA3P_INTEL_MKL)
    // Multiply everything with (1/beta) to calculate
    // (1/beta)*C = (alpha/beta) * A + B
    // then scale C to get proper result
    if(beta != T_Scalar(1)) { alpha = alpha / beta; }
    mkl::csr_add(m, n, alpha, op_t::N,
            rowptrA, colidxA, valuesA,
            rowptrB, colidxB, valuesB,
            rowptrC, colidxC, valuesC);
    if(beta != T_Scalar(1)) { blas::scal((*rowptrC)[n], beta, *valuesC, 1); }
#elif defined(CLA3P_ARMPL)
    armpl::csr_add(m, n,
            alpha, op_t::N, rowptrA, colidxA, valuesA,
            beta, op_t::N, rowptrB, colidxB, valuesB,
            rowptrC, colidxC, valuesC);
#else
    throw err::InvalidOp(msg::MissingThirdParty());
#endif
}
/*-------------------------------------------------*/
#define instantiate_add(T_Scl) \
template void add(int_t, int_t, \
        T_Scl, const int_t*, const int_t*, const T_Scl*, \
        T_Scl, const int_t*, const int_t*, const T_Scl*, \
        int_t**, int_t**, T_Scl**)
instantiate_add(real_t);
instantiate_add(real4_t);
instantiate_add(complex_t);
instantiate_add(complex8_t);
#undef instantiate_add
/*-------------------------------------------------*/
template <typename T_Scalar>
void gem_x_vec(op_t opA, int_t m, int_t n, T_Scalar alpha,
        const int_t *rowptr, const int_t *colidx, const T_Scalar *values,
        const T_Scalar *x, T_Scalar beta, T_Scalar *y)
{
#if defined(CLA3P_INTEL_MKL)
    Property pr = Property::General();
    mkl::csr_mv(pr.type(), pr.uplo(), m, n, alpha, opA, rowptr, colidx, values, x, beta, y);
#elif defined(CLA3P_ARMPL)
    Property pr = Property::General();
    armpl::csr_mv(pr.type(), pr.uplo(), m, n, alpha, opA, rowptr, colidx, values, x, beta, y);
#else
    throw err::InvalidOp(msg::MissingThirdParty());
#endif
}
/*-------------------------------------------------*/
#define instantiate_gem_x_vec(T_Scl) \
template void gem_x_vec(op_t, int_t, int_t, T_Scl, \
        const int_t*, const int_t*, const T_Scl*, \
        const T_Scl*, T_Scl, T_Scl*)
instantiate_gem_x_vec(real_t);
instantiate_gem_x_vec(real4_t);
instantiate_gem_x_vec(complex_t);
instantiate_gem_x_vec(complex8_t);
#undef instantiate_gem_x_vec
/*-------------------------------------------------*/
template <typename T_Scalar>
void sym_x_vec(uplo_t uplo, int_t n, T_Scalar alpha,
        const int_t *rowptr, const int_t *colidx, const T_Scalar *values,
        const T_Scalar *x, T_Scalar beta, T_Scalar *y)
{
#if defined(CLA3P_INTEL_MKL)
    Property pr = Property(prop_t::Symmetric, uplo);
    mkl::csr_mv(pr.type(), pr.uplo(), n, n, alpha, op_t::N, rowptr, colidx, values, x, beta, y);
#elif defined(CLA3P_ARMPL)
    Property pr = Property(prop_t::Symmetric, uplo);
    armpl::csr_mv(pr.type(), pr.uplo(), n, n, alpha, op_t::N, rowptr, colidx, values, x, beta, y);
#else
    throw err::InvalidOp(msg::MissingThirdParty());
#endif
}
/*-------------------------------------------------*/
#define instantiate_sym_x_vec(T_Scl) \
template void sym_x_vec(uplo_t, int_t, T_Scl, \
        const int_t*, const int_t*, const T_Scl*, \
        const T_Scl*, T_Scl, T_Scl*)
instantiate_sym_x_vec(real_t);
instantiate_sym_x_vec(real4_t);
instantiate_sym_x_vec(complex_t);
instantiate_sym_x_vec(complex8_t);
#undef instantiate_sym_x_vec
/*-------------------------------------------------*/
template <typename T_Scalar>
void hem_x_vec(uplo_t uplo, int_t n, T_Scalar alpha,
        const int_t *rowptr, const int_t *colidx, const T_Scalar *values,
        const T_Scalar *x, T_Scalar beta, T_Scalar *y)
{
#if defined(CLA3P_INTEL_MKL)
    Property pr = sanitizeProperty<T_Scalar>(Property(prop_t::Hermitian, uplo));
    mkl::csr_mv(pr.type(), pr.uplo(), n, n, alpha, op_t::N, rowptr, colidx, values, x, beta, y);
#elif defined(CLA3P_ARMPL)
    Property pr = sanitizeProperty<T_Scalar>(Property(prop_t::Hermitian, uplo));
    armpl::csr_mv(pr.type(), pr.uplo(), n, n, alpha, op_t::N, rowptr, colidx, values, x, beta, y);
#else
    throw err::InvalidOp(msg::MissingThirdParty());
#endif
}
/*-------------------------------------------------*/
#define instantiate_hem_x_vec(T_Scl) \
template void hem_x_vec(uplo_t, int_t, T_Scl, \
        const int_t*, const int_t*, const T_Scl*, \
        const T_Scl*, T_Scl, T_Scl*)
instantiate_hem_x_vec(real_t);
instantiate_hem_x_vec(real4_t);
instantiate_hem_x_vec(complex_t);
instantiate_hem_x_vec(complex8_t);
#undef instantiate_hem_x_vec
/*-------------------------------------------------*/
template <typename T_Scalar>
void gem_x_gem(op_t opA, int_t m, int_t n, int_t k, T_Scalar alpha,
        const int_t *rowptr, const int_t *colidx, const T_Scalar *values,
        const T_Scalar *b, int_t ldb, T_Scalar beta, T_Scalar *c, int_t ldc)
{
#if defined(CLA3P_INTEL_MKL)
    int_t mA = (opA == op_t::N ? m : k);
    int_t nA = (opA == op_t::N ? k : m);
    Property pr = Property::General();
    mkl::csr_mm(pr.type(), pr.uplo(), mA, nA, alpha, opA, rowptr, colidx, values, n, b, ldb, beta, c, ldc);
#elif defined(CLA3P_ARMPL)
    int_t mA = (opA == op_t::N ? m : k);
    int_t nA = (opA == op_t::N ? k : m);
    Property pr = Property::General();
    armpl::csr_mm(pr.type(), pr.uplo(), mA, nA, alpha, opA, rowptr, colidx, values, n, b, ldb, beta, c, ldc);
#else
    throw err::InvalidOp(msg::MissingThirdParty());
#endif
}
/*-------------------------------------------------*/
#define instantiate_gem_x_gem(T_Scl) \
template void gem_x_gem(op_t, int_t, int_t, int_t, T_Scl, \
        const int_t*, const int_t*, const T_Scl*, \
        const T_Scl*, int_t, T_Scl, T_Scl*, int_t )
instantiate_gem_x_gem(real_t);
instantiate_gem_x_gem(real4_t);
instantiate_gem_x_gem(complex_t);
instantiate_gem_x_gem(complex8_t);
#undef instantiate_gem_x_gem
/*-------------------------------------------------*/
template <typename T_Scalar>
void sym_x_gem(uplo_t uplo, int_t m, int_t n, T_Scalar alpha,
        const int_t *rowptr, const int_t *colidx, const T_Scalar *values,
        const T_Scalar *b, int_t ldb, T_Scalar beta, T_Scalar *c, int_t ldc)
{
#if defined(CLA3P_INTEL_MKL)
    Property pr = Property(prop_t::Symmetric, uplo);
    mkl::csr_mm(pr.type(), pr.uplo(), m, m, alpha, op_t::N, rowptr, colidx, values, n, b, ldb, beta, c, ldc);
#elif defined(CLA3P_ARMPL)
    Property pr = Property(prop_t::Symmetric, uplo);
    armpl::csr_mm(pr.type(), pr.uplo(), m, m, alpha, op_t::N, rowptr, colidx, values, n, b, ldb, beta, c, ldc);
#else
    throw err::InvalidOp(msg::MissingThirdParty());
#endif
}
/*-------------------------------------------------*/
#define instantiate_sym_x_gem(T_Scl) \
    template void sym_x_gem(uplo_t, int_t, int_t, T_Scl, \
            const int_t*, const int_t*, const T_Scl*, \
            const T_Scl*, int_t, T_Scl, T_Scl*, int_t)
instantiate_sym_x_gem(real_t);
instantiate_sym_x_gem(real4_t);
instantiate_sym_x_gem(complex_t);
instantiate_sym_x_gem(complex8_t);
#undef instantiate_sym_x_gem
/*-------------------------------------------------*/
template <typename T_Scalar>
void hem_x_gem(uplo_t uplo, int_t m, int_t n, T_Scalar alpha,
        const int_t *rowptr, const int_t *colidx, const T_Scalar *values,
        const T_Scalar *b, int_t ldb, T_Scalar beta, T_Scalar *c, int_t ldc)
{
#if defined(CLA3P_INTEL_MKL)
    Property pr = sanitizeProperty<T_Scalar>(Property(prop_t::Hermitian, uplo));
    mkl::csr_mm(pr.type(), pr.uplo(), m, m, alpha, op_t::N, rowptr, colidx, values, n, b, ldb, beta, c, ldc);
#elif defined(CLA3P_ARMPL)
    Property pr = sanitizeProperty<T_Scalar>(Property(prop_t::Hermitian, uplo));
    armpl::csr_mm(pr.type(), pr.uplo(), m, m, alpha, op_t::N, rowptr, colidx, values, n, b, ldb, beta, c, ldc);
#else
    throw err::InvalidOp(msg::MissingThirdParty());
#endif
}
/*-------------------------------------------------*/
#define instantiate_hem_x_gem(T_Scl) \
template void hem_x_gem(uplo_t, int_t, int_t, T_Scl, \
        const int_t*, const int_t*, const T_Scl*, \
        const T_Scl*, int_t, T_Scl, T_Scl*, int_t)
instantiate_hem_x_gem(real_t);
instantiate_hem_x_gem(real4_t);
instantiate_hem_x_gem(complex_t);
instantiate_hem_x_gem(complex8_t);
#undef instantiate_hem_x_gem
/*-------------------------------------------------*/
template <typename T_Scalar>
void gem_x_gem(int_t m, int_t n, int_t k, T_Scalar alpha,
        op_t opA, const int_t *rowptrA, const int_t *colidxA, const T_Scalar *valuesA,
        op_t opB, const int_t *rowptrB, const int_t *colidxB, const T_Scalar *valuesB,
        T_Scalar beta, T_Scalar *c, int_t ldc)
{
    int_t mA = (opA == op_t::N ? m : k);
    int_t nA = (opA == op_t::N ? k : m);
    int_t mB = (opB == op_t::N ? k : n);
    int_t nB = (opB == op_t::N ? n : k);

#if defined(CLA3P_INTEL_MKL)
    mkl::csr_sp2md(alpha,
            opA, mA, nA, rowptrA, colidxA, valuesA,
            opB, mB, nB, rowptrB, colidxB, valuesB,
            beta, c, ldc);
#elif defined(CLA3P_ARMPL)
    armpl::csr_spmm(alpha,
            opA, mA, nA, rowptrA, colidxA, valuesA,
            opB, mB, nB, rowptrB, colidxB, valuesB,
            beta, c, ldc);
#else
    throw err::InvalidOp(msg::MissingIntelMKL());
#endif
}
/*-------------------------------------------------*/
#define instantiate_gem_x_gem(T_Scl) \
template void gem_x_gem(int_t, int_t, int_t, T_Scl, \
        op_t, const int_t*, const int_t*, const T_Scl*, \
        op_t, const int_t*, const int_t*, const T_Scl*, \
        T_Scl, T_Scl*, int_t)
instantiate_gem_x_gem(real_t);
instantiate_gem_x_gem(real4_t);
instantiate_gem_x_gem(complex_t);
instantiate_gem_x_gem(complex8_t);
#undef instantiate_gem_x_gem
/*-------------------------------------------------*/
template <typename T_Scalar>
void gem_x_gem(int_t m, int_t n, int_t k,
        op_t opA, const int_t *rowptrA, const int_t *colidxA, const T_Scalar *valuesA,
        op_t opB, const int_t *rowptrB, const int_t *colidxB, const T_Scalar *valuesB,
        int_t **rowptrC, int_t **colidxC, T_Scalar **valuesC)
{
    int_t mA = (opA == op_t::N ? m : k);
    int_t nA = (opA == op_t::N ? k : m);
    int_t mB = (opB == op_t::N ? k : n);
    int_t nB = (opB == op_t::N ? n : k);

#if defined(CLA3P_INTEL_MKL)
    mkl::csr_sp2m(
            opA, mA, nA, rowptrA, colidxA, valuesA,
            opB, mB, nB, rowptrB, colidxB, valuesB,
            rowptrC, colidxC, valuesC);
#elif defined(CLA3P_ARMPL)
    armpl::csr_spmm(T_Scalar(1),
            opA, mA, nA, rowptrA, colidxA, valuesA,
            opB, mB, nB, rowptrB, colidxB, valuesB,
            rowptrC, colidxC, valuesC);
#else
    throw err::InvalidOp(msg::MissingThirdParty());
#endif
}
/*-------------------------------------------------*/
#define instantiate_gem_x_gem(T_Scl) \
template void gem_x_gem(int_t, int_t, int_t, \
        op_t, const int_t*, const int_t*, const T_Scl*, \
        op_t, const int_t*, const int_t*, const T_Scl*, \
        int_t**, int_t**, T_Scl**)
instantiate_gem_x_gem(real_t);
instantiate_gem_x_gem(real4_t);
instantiate_gem_x_gem(complex_t);
instantiate_gem_x_gem(complex8_t);
#undef instantiate_gem_x_gem
/*-------------------------------------------------*/
} // namespace csr
} // namespace blk
} // namespace cla3p
/*-------------------------------------------------*/
