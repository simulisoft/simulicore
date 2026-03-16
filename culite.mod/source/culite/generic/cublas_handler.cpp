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
#include <culite/generic/cublas_handler.hpp>

// system

// 3rd

// culite
#include <culite/error/cuda.hpp>

/*-------------------------------------------------*/
namespace culite {
/*-------------------------------------------------*/
CuBlasSpm::CuBlasSpm(CuBlasHandler& cuBlasHandler, cublasPointerMode_t mode)
    : m_handler(cuBlasHandler)
{
    m_oldMode = m_handler.setPointerMode(mode);
}
/*-------------------------------------------------*/
CuBlasSpm::~CuBlasSpm()
{
    m_handler.setPointerMode(m_oldMode);
}
/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/
CuBlasHandler::CuBlasHandler()
{
    cublasStatus_t cublasStatus = cublasCreate(&m_handle);
    err::check_cublas(cublasStatus);
}
/*-------------------------------------------------*/
CuBlasHandler::~CuBlasHandler()
{
    cublasStatus_t cublasStatus = cublasDestroy(m_handle);
    err::check_cublas(cublasStatus);
}
/*-------------------------------------------------*/
cublasPointerMode_t CuBlasHandler::setPointerMode(cublasPointerMode_t mode)
{
    cublasPointerMode_t ret = pointerMode();
    if (ret != mode) {
        cublasStatus_t status = cublasSetPointerMode(m_handle, mode);
        err::check_cublas(status);
    }
    return ret;
}
/*-------------------------------------------------*/
cublasPointerMode_t CuBlasHandler::pointerMode()
{
    cublasPointerMode_t ret;
    cublasStatus_t status = cublasGetPointerMode(m_handle, &ret);
    err::check_cublas(status);
    return ret;
}
/*-------------------------------------------------*/
cublasHandle_t CuBlasHandler::handle()
{
    return m_handle;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void CuBlasHandler::iamax(int_t n, const T_Scalar *x, int_t incx, int_t *result)
{
    cublas::iamax(handle(), n, x, incx, result);
}
/*-------------------------------------------------*/
template void CuBlasHandler::iamax<real_t>(int_t, const real_t*, int_t, int_t*);
template void CuBlasHandler::iamax<real4_t>(int_t, const real4_t*, int_t, int_t*);
template void CuBlasHandler::iamax<complex_t>(int_t, const complex_t*, int_t, int_t*);
template void CuBlasHandler::iamax<complex8_t>(int_t, const complex8_t*, int_t, int_t*);
/*-------------------------------------------------*/
template <typename T_Scalar>
void CuBlasHandler::asum(int_t n, const T_Scalar *x, int_t incx, typename TypeTraits<T_Scalar>::real_type *result)
{
    cublas::asum(handle(), n, x, incx, result);
}
/*-------------------------------------------------*/
template void CuBlasHandler::asum<real_t>(int_t, const real_t*, int_t, typename TypeTraits<real_t>::real_type*);
template void CuBlasHandler::asum<real4_t>(int_t, const real4_t*, int_t, typename TypeTraits<real4_t>::real_type*);
template void CuBlasHandler::asum<complex_t>(int_t, const complex_t*, int_t, typename TypeTraits<complex_t>::real_type*);
template void CuBlasHandler::asum<complex8_t>(int_t, const complex8_t*, int_t, typename TypeTraits<complex8_t>::real_type*);
/*-------------------------------------------------*/
template <typename T_Scalar>
void CuBlasHandler::nrm2(int_t n, const T_Scalar *x, int_t incx, typename TypeTraits<T_Scalar>::real_type *result)
{
    cublas::nrm2(handle(), n, x, incx, result);
}
/*-------------------------------------------------*/
template void CuBlasHandler::nrm2<real_t>(int_t, const real_t*, int_t, typename TypeTraits<real_t>::real_type*);
template void CuBlasHandler::nrm2<real4_t>(int_t, const real4_t*, int_t, typename TypeTraits<real4_t>::real_type*);
template void CuBlasHandler::nrm2<complex_t>(int_t, const complex_t*, int_t, typename TypeTraits<complex_t>::real_type*);
template void CuBlasHandler::nrm2<complex8_t>(int_t, const complex8_t*, int_t, typename TypeTraits<complex8_t>::real_type*);
/*-------------------------------------------------*/
template <typename T_Scalar>
void CuBlasHandler::copy(int_t n, const T_Scalar *x, int_t incx, T_Scalar *y, int_t incy)
{
    cublas::copy(handle(), n, x, incx, y, incy);
}
/*-------------------------------------------------*/
template void CuBlasHandler::copy<real_t>(int_t, const real_t*, int_t, real_t*, int_t);
template void CuBlasHandler::copy<real4_t>(int_t, const real4_t*, int_t, real4_t*, int_t);
template void CuBlasHandler::copy<complex_t>(int_t, const complex_t*, int_t, complex_t*, int_t);
template void CuBlasHandler::copy<complex8_t>(int_t, const complex8_t*, int_t, complex8_t*, int_t);
/*-------------------------------------------------*/
template <typename T_Scalar>
void CuBlasHandler::scal(int_t n, const T_Scalar *alpha, T_Scalar *x, int_t incx)
{
    cublas::scal(handle(), n, alpha, x, incx);
}
/*-------------------------------------------------*/
template void CuBlasHandler::scal<real_t>(int_t, const real_t*, real_t*, int_t);
template void CuBlasHandler::scal<real4_t>(int_t, const real4_t*, real4_t*, int_t);
template void CuBlasHandler::scal<complex_t>(int_t, const complex_t*, complex_t*, int_t);
template void CuBlasHandler::scal<complex8_t>(int_t, const complex8_t*, complex8_t*, int_t);
/*-------------------------------------------------*/
template <typename T_Scalar>
void CuBlasHandler::dot(int_t n, const T_Scalar *x, int_t incx, const T_Scalar *y, int_t incy, T_Scalar *result)
{
    cublas::dot(handle(), n, x, incx, y, incy, result);
}
/*-------------------------------------------------*/
template void CuBlasHandler::dot<real_t>(int_t, const real_t*, int_t, const real_t*, int_t, real_t*);
template void CuBlasHandler::dot<real4_t>(int_t, const real4_t*, int_t, const real4_t*, int_t, real4_t*);
template void CuBlasHandler::dot<complex_t>(int_t, const complex_t*, int_t, const complex_t*, int_t, complex_t*);
template void CuBlasHandler::dot<complex8_t>(int_t, const complex8_t*, int_t, const complex8_t*, int_t, complex8_t*);
/*-------------------------------------------------*/
template <typename T_Scalar>
void CuBlasHandler::dotc(int_t n, const T_Scalar *x, int_t incx, const T_Scalar *y, int_t incy, T_Scalar *result)
{
    cublas::dot(handle(), n, x, incx, y, incy, result);
}
/*-------------------------------------------------*/
template void CuBlasHandler::dotc<real_t>(int_t, const real_t*, int_t, const real_t*, int_t, real_t*);
template void CuBlasHandler::dotc<real4_t>(int_t, const real4_t*, int_t, const real4_t*, int_t, real4_t*);
template void CuBlasHandler::dotc<complex_t>(int_t, const complex_t*, int_t, const complex_t*, int_t, complex_t*);
template void CuBlasHandler::dotc<complex8_t>(int_t, const complex8_t*, int_t, const complex8_t*, int_t, complex8_t*);
/*-------------------------------------------------*/
template <typename T_Scalar>
void CuBlasHandler::axpy(int_t n, const T_Scalar *alpha, const T_Scalar *x, int_t incx, T_Scalar *y, int_t incy)
{
    cublas::axpy(handle(), n, alpha, x, incx, y, incy);
}
/*-------------------------------------------------*/
template void CuBlasHandler::axpy<real_t>(int_t, const real_t*, const real_t*, int_t, real_t*, int_t);
template void CuBlasHandler::axpy<real4_t>(int_t, const real4_t*, const real4_t*, int_t, real4_t*, int_t);
template void CuBlasHandler::axpy<complex_t>(int_t, const complex_t*, const complex_t*, int_t, complex_t*, int_t);
template void CuBlasHandler::axpy<complex8_t>(int_t, const complex8_t*, const complex8_t*, int_t, complex8_t*, int_t);
/*-------------------------------------------------*/
template <typename T_Scalar>
void CuBlasHandler::geam(op_t opA,
                          op_t opB,
                          int_t m, int_t n,
                          const T_Scalar *alpha, const T_Scalar *a, int_t lda,
                          const T_Scalar *beta,  const T_Scalar *b, int_t ldb,
                          T_Scalar *c, int_t ldc)
{
    cublas::geam(handle(),
                 cublas::cla3pOp2cublasOp(opA),
                 cublas::cla3pOp2cublasOp(opB),
                 m, n,
                 alpha, a, lda,
                 beta, b, ldb,
                 c, ldc);
}
/*-------------------------------------------------*/
template void CuBlasHandler::geam<real_t>(op_t, op_t, int_t, int_t, const real_t*, const real_t*, int_t, const real_t*, const real_t*, int_t, real_t*, int_t);
template void CuBlasHandler::geam<real4_t>(op_t, op_t, int_t, int_t, const real4_t*, const real4_t*, int_t, const real4_t*, const real4_t*, int_t, real4_t*, int_t);
template void CuBlasHandler::geam<complex_t>(op_t, op_t, int_t, int_t, const complex_t*, const complex_t*, int_t, const complex_t*, const complex_t*, int_t, complex_t*, int_t);
template void CuBlasHandler::geam<complex8_t>(op_t, op_t, int_t, int_t, const complex8_t*, const complex8_t*, int_t, const complex8_t*, const complex8_t*, int_t, complex8_t*, int_t);
/*-------------------------------------------------*/
template <typename T_Scalar>
void CuBlasHandler::dgmm(side_t side,
                          int_t m, int_t n,
                          const T_Scalar *a, int_t lda,
                          const T_Scalar *x, int_t incx,
                          T_Scalar *c, int_t ldc)
{
    cublas::dgmm(handle(),
                 cublas::cla3pSide2cublasSide(side),
                 m, n,
                 a, lda,
                 x, incx,
                 c, ldc);
}
/*-------------------------------------------------*/
template void CuBlasHandler::dgmm<real_t>(side_t, int_t, int_t, const real_t*, int_t, const real_t*, int_t, real_t*, int_t);
template void CuBlasHandler::dgmm<real4_t>(side_t, int_t, int_t, const real4_t*, int_t, const real4_t*, int_t, real4_t*, int_t);
template void CuBlasHandler::dgmm<complex_t>(side_t, int_t, int_t, const complex_t*, int_t, const complex_t*, int_t, complex_t*, int_t);
template void CuBlasHandler::dgmm<complex8_t>(side_t, int_t, int_t, const complex8_t*, int_t, const complex8_t*, int_t, complex8_t*, int_t);
/*-------------------------------------------------*/
template <typename T_Scalar>
void CuBlasHandler::ger(int_t m, int_t n,
                         const T_Scalar *alpha,
                         const T_Scalar *x, int_t incx,
                         const T_Scalar *y, int_t incy,
                         T_Scalar *a, int_t lda)
{
    cublas::ger(handle(), m, n, alpha, x, incx, y, incy, a, lda);
}
/*-------------------------------------------------*/
template void CuBlasHandler::ger<real_t>(int_t, int_t, const real_t*, const real_t*, int_t, const real_t*, int_t, real_t*, int_t);
template void CuBlasHandler::ger<real4_t>(int_t, int_t, const real4_t*, const real4_t*, int_t, const real4_t*, int_t, real4_t*, int_t);
template void CuBlasHandler::ger<complex_t>(int_t, int_t, const complex_t*, const complex_t*, int_t, const complex_t*, int_t, complex_t*, int_t);
template void CuBlasHandler::ger<complex8_t>(int_t, int_t, const complex8_t*, const complex8_t*, int_t, const complex8_t*, int_t, complex8_t*, int_t);
/*-------------------------------------------------*/
template <typename T_Scalar>
void CuBlasHandler::gerc(int_t m, int_t n,
                          const T_Scalar *alpha,
                          const T_Scalar *x, int_t incx,
                          const T_Scalar *y, int_t incy,
                          T_Scalar *a, int_t lda)
{
    cublas::gerc(handle(), m, n, alpha, x, incx, y, incy, a, lda);
}
/*-------------------------------------------------*/
template void CuBlasHandler::gerc<real_t>(int_t, int_t, const real_t*, const real_t*, int_t, const real_t*, int_t, real_t*, int_t);
template void CuBlasHandler::gerc<real4_t>(int_t, int_t, const real4_t*, const real4_t*, int_t, const real4_t*, int_t, real4_t*, int_t);
template void CuBlasHandler::gerc<complex_t>(int_t, int_t, const complex_t*, const complex_t*, int_t, const complex_t*, int_t, complex_t*, int_t);
template void CuBlasHandler::gerc<complex8_t>(int_t, int_t, const complex8_t*, const complex8_t*, int_t, const complex8_t*, int_t, complex8_t*, int_t);
/*-------------------------------------------------*/
template <typename T_Scalar>
void CuBlasHandler::syr(uplo_t uplo,
                         int_t n,
                         const T_Scalar *alpha,
                         const T_Scalar *x, int_t incx,
                         T_Scalar *a, int_t lda)
{
    cublas::syr(handle(),
                cublas::cla3pUplo2cublasUplo(uplo),
                n, alpha, x, incx, a, lda);
}
/*-------------------------------------------------*/
template void CuBlasHandler::syr<real_t>(uplo_t, int_t, const real_t*, const real_t*, int_t, real_t*, int_t);
template void CuBlasHandler::syr<real4_t>(uplo_t, int_t, const real4_t*, const real4_t*, int_t, real4_t*, int_t);
template void CuBlasHandler::syr<complex_t>(uplo_t, int_t, const complex_t*, const complex_t*, int_t, complex_t*, int_t);
template void CuBlasHandler::syr<complex8_t>(uplo_t, int_t, const complex8_t*, const complex8_t*, int_t, complex8_t*, int_t);
/*-------------------------------------------------*/
template <typename T_Scalar>
void CuBlasHandler::her(uplo_t uplo,
                         int_t n,
                         const typename TypeTraits<T_Scalar>::real_type *alpha,
                         const T_Scalar *x, int_t incx,
                         T_Scalar *a, int_t lda)
{
    cublas::her(handle(),
                cublas::cla3pUplo2cublasUplo(uplo),
                n, alpha, x, incx, a, lda);
}
/*-------------------------------------------------*/
template void CuBlasHandler::her<real_t>(uplo_t, int_t, const typename TypeTraits<real_t>::real_type*, const real_t*, int_t, real_t*, int_t);
template void CuBlasHandler::her<real4_t>(uplo_t, int_t, const typename TypeTraits<real4_t>::real_type*, const real4_t*, int_t, real4_t*, int_t);
template void CuBlasHandler::her<complex_t>(uplo_t, int_t, const typename TypeTraits<complex_t>::real_type*, const complex_t*, int_t, complex_t*, int_t);
template void CuBlasHandler::her<complex8_t>(uplo_t, int_t, const typename TypeTraits<complex8_t>::real_type*, const complex8_t*, int_t, complex8_t*, int_t);
/*-------------------------------------------------*/
template <typename T_Scalar>
void CuBlasHandler::gemv(op_t op,
                          int_t m, int_t n, const T_Scalar* alpha,
                          const T_Scalar *a, int_t lda,
                          const T_Scalar *x, int_t incx,
                          const T_Scalar* beta, T_Scalar *y, int_t incy)
{
    cublas::gemv(handle(),
                 cublas::cla3pOp2cublasOp(op),
                 m, n, alpha,
                 a, lda,
                 x, incx,
                 beta, y, incy);
}
/*-------------------------------------------------*/
template void CuBlasHandler::gemv<real_t>(op_t, int_t, int_t, const real_t*, const real_t*, int_t, const real_t*, int_t, const real_t*, real_t*, int_t);
template void CuBlasHandler::gemv<real4_t>(op_t, int_t, int_t, const real4_t*, const real4_t*, int_t, const real4_t*, int_t, const real4_t*, real4_t*, int_t);
template void CuBlasHandler::gemv<complex_t>(op_t, int_t, int_t, const complex_t*, const complex_t*, int_t, const complex_t*, int_t, const complex_t*, complex_t*, int_t);
template void CuBlasHandler::gemv<complex8_t>(op_t, int_t, int_t, const complex8_t*, const complex8_t*, int_t, const complex8_t*, int_t, const complex8_t*, complex8_t*, int_t);
/*-------------------------------------------------*/
template <typename T_Scalar>
void CuBlasHandler::symv(uplo_t uplo,
                          int_t n, const T_Scalar* alpha,
                          const T_Scalar *a, int_t lda,
                          const T_Scalar *x, int_t incx,
                          const T_Scalar* beta, T_Scalar *y, int_t incy)
{
    cublas::symv(handle(),
                 cublas::cla3pUplo2cublasUplo(uplo),
                 n, alpha,
                 a, lda,
                 x, incx,
                 beta, y, incy);
}
/*-------------------------------------------------*/
template void CuBlasHandler::symv<real_t>(uplo_t, int_t, const real_t*, const real_t*, int_t, const real_t*, int_t, const real_t*, real_t*, int_t);
template void CuBlasHandler::symv<real4_t>(uplo_t, int_t, const real4_t*, const real4_t*, int_t, const real4_t*, int_t, const real4_t*, real4_t*, int_t);
template void CuBlasHandler::symv<complex_t>(uplo_t, int_t, const complex_t*, const complex_t*, int_t, const complex_t*, int_t, const complex_t*, complex_t*, int_t);
template void CuBlasHandler::symv<complex8_t>(uplo_t, int_t, const complex8_t*, const complex8_t*, int_t, const complex8_t*, int_t, const complex8_t*, complex8_t*, int_t);
/*-------------------------------------------------*/
template <typename T_Scalar>
void CuBlasHandler::hemv(uplo_t uplo,
                          int_t n, const T_Scalar* alpha,
                          const T_Scalar *a, int_t lda,
                          const T_Scalar *x, int_t incx,
                          const T_Scalar* beta, T_Scalar *y, int_t incy)
{
    cublas::hemv(handle(),
                 cublas::cla3pUplo2cublasUplo(uplo),
                 n, alpha,
                 a, lda,
                 x, incx,
                 beta, y, incy);
}
/*-------------------------------------------------*/
template void CuBlasHandler::hemv<real_t>(uplo_t, int_t, const real_t*, const real_t*, int_t, const real_t*, int_t, const real_t*, real_t*, int_t);
template void CuBlasHandler::hemv<real4_t>(uplo_t, int_t, const real4_t*, const real4_t*, int_t, const real4_t*, int_t, const real4_t*, real4_t*, int_t);
template void CuBlasHandler::hemv<complex_t>(uplo_t, int_t, const complex_t*, const complex_t*, int_t, const complex_t*, int_t, const complex_t*, complex_t*, int_t);
template void CuBlasHandler::hemv<complex8_t>(uplo_t, int_t, const complex8_t*, const complex8_t*, int_t, const complex8_t*, int_t, const complex8_t*, complex8_t*, int_t);
/*-------------------------------------------------*/
template <typename T_Scalar>
void CuBlasHandler::trmv(uplo_t uplo,
                          op_t trans,
                          diag_t diag,
                          int_t n, const T_Scalar* a, int_t lda,
                          T_Scalar *x, int_t incx)
{
    cublas::trmv(handle(),
                 cublas::cla3pUplo2cublasUplo(uplo),
                 cublas::cla3pOp2cublasOp(trans),
                 cublas::cla3pDiag2cublasDiag(diag),
                 n, a, lda,
                 x, incx);
}
/*-------------------------------------------------*/
template void CuBlasHandler::trmv<real_t>(uplo_t, op_t, diag_t, int_t, const real_t*, int_t, real_t*, int_t);
template void CuBlasHandler::trmv<real4_t>(uplo_t, op_t, diag_t, int_t, const real4_t*, int_t, real4_t*, int_t);
template void CuBlasHandler::trmv<complex_t>(uplo_t, op_t, diag_t, int_t, const complex_t*, int_t, complex_t*, int_t);
template void CuBlasHandler::trmv<complex8_t>(uplo_t, op_t, diag_t, int_t, const complex8_t*, int_t, complex8_t*, int_t);
/*-------------------------------------------------*/
template <typename T_Scalar>
void CuBlasHandler::trsv(uplo_t uplo,
                          op_t trans,
                          diag_t diag,
                          int_t n, const T_Scalar* a, int_t lda,
                          T_Scalar *x, int_t incx)
{
    cublas::trsv(handle(),
                 cublas::cla3pUplo2cublasUplo(uplo),
                 cublas::cla3pOp2cublasOp(trans),
                 cublas::cla3pDiag2cublasDiag(diag),
                 n, a, lda,
                 x, incx);
}
/*-------------------------------------------------*/
template void CuBlasHandler::trsv<real_t>(uplo_t, op_t, diag_t, int_t, const real_t*, int_t, real_t*, int_t);
template void CuBlasHandler::trsv<real4_t>(uplo_t, op_t, diag_t, int_t, const real4_t*, int_t, real4_t*, int_t);
template void CuBlasHandler::trsv<complex_t>(uplo_t, op_t, diag_t, int_t, const complex_t*, int_t, complex_t*, int_t);
template void CuBlasHandler::trsv<complex8_t>(uplo_t, op_t, diag_t, int_t, const complex8_t*, int_t, complex8_t*, int_t);
/*-------------------------------------------------*/
template <typename T_Scalar>
void CuBlasHandler::gemm(op_t transa,
                          op_t transb,
                          int_t m, int_t n, int_t k,
                          const T_Scalar* alpha,
                          const T_Scalar* a, int_t lda,
                          const T_Scalar* b, int_t ldb,
                          const T_Scalar* beta,
                          T_Scalar* c, int_t ldc)
{
    cublas::gemm(handle(),
                 cublas::cla3pOp2cublasOp(transa),
                 cublas::cla3pOp2cublasOp(transb),
                 m, n, k,
                 alpha,
                 a, lda,
                 b, ldb,
                 beta,
                 c, ldc);
}
/*-------------------------------------------------*/
template void CuBlasHandler::gemm<real_t>(op_t, op_t, int_t, int_t, int_t, const real_t*, const real_t*, int_t, const real_t*, int_t, const real_t*, real_t*, int_t);
template void CuBlasHandler::gemm<real4_t>(op_t, op_t, int_t, int_t, int_t, const real4_t*, const real4_t*, int_t, const real4_t*, int_t, const real4_t*, real4_t*, int_t);
template void CuBlasHandler::gemm<complex_t>(op_t, op_t, int_t, int_t, int_t, const complex_t*, const complex_t*, int_t, const complex_t*, int_t, const complex_t*, complex_t*, int_t);
template void CuBlasHandler::gemm<complex8_t>(op_t, op_t, int_t, int_t, int_t, const complex8_t*, const complex8_t*, int_t, const complex8_t*, int_t, const complex8_t*, complex8_t*, int_t);
/*-------------------------------------------------*/
template <typename T_Scalar>
void CuBlasHandler::symm(side_t side,
                          uplo_t uplo,
                          int_t m, int_t n,
                          const T_Scalar* alpha,
                          const T_Scalar* a, int_t lda,
                          const T_Scalar* b, int_t ldb,
                          const T_Scalar* beta,
                          T_Scalar* c, int_t ldc)
{
    cublas::symm(handle(),
                 cublas::cla3pSide2cublasSide(side),
                 cublas::cla3pUplo2cublasUplo(uplo),
                 m, n,
                 alpha,
                 a, lda,
                 b, ldb,
                 beta,
                 c, ldc);
}
/*-------------------------------------------------*/
template void CuBlasHandler::symm<real_t>(side_t, uplo_t, int_t, int_t, const real_t*, const real_t*, int_t, const real_t*, int_t, const real_t*, real_t*, int_t);
template void CuBlasHandler::symm<real4_t>(side_t, uplo_t, int_t, int_t, const real4_t*, const real4_t*, int_t, const real4_t*, int_t, const real4_t*, real4_t*, int_t);
template void CuBlasHandler::symm<complex_t>(side_t, uplo_t, int_t, int_t, const complex_t*, const complex_t*, int_t, const complex_t*, int_t, const complex_t*, complex_t*, int_t);
template void CuBlasHandler::symm<complex8_t>(side_t, uplo_t, int_t, int_t, const complex8_t*, const complex8_t*, int_t, const complex8_t*, int_t, const complex8_t*, complex8_t*, int_t);
/*-------------------------------------------------*/
template <typename T_Scalar>
void CuBlasHandler::hemm(side_t side,
                          uplo_t uplo,
                          int_t m, int_t n,
                          const T_Scalar* alpha,
                          const T_Scalar* a, int_t lda,
                          const T_Scalar* b, int_t ldb,
                          const T_Scalar* beta,
                          T_Scalar* c, int_t ldc)
{
    cublas::hemm(handle(),
                 cublas::cla3pSide2cublasSide(side),
                 cublas::cla3pUplo2cublasUplo(uplo),
                 m, n,
                 alpha,
                 a, lda,
                 b, ldb,
                 beta,
                 c, ldc);
}
/*-------------------------------------------------*/
template void CuBlasHandler::hemm<real_t>(side_t, uplo_t, int_t, int_t, const real_t*, const real_t*, int_t, const real_t*, int_t, const real_t*, real_t*, int_t);
template void CuBlasHandler::hemm<real4_t>(side_t, uplo_t, int_t, int_t, const real4_t*, const real4_t*, int_t, const real4_t*, int_t, const real4_t*, real4_t*, int_t);
template void CuBlasHandler::hemm<complex_t>(side_t, uplo_t, int_t, int_t, const complex_t*, const complex_t*, int_t, const complex_t*, int_t, const complex_t*, complex_t*, int_t);
template void CuBlasHandler::hemm<complex8_t>(side_t, uplo_t, int_t, int_t, const complex8_t*, const complex8_t*, int_t, const complex8_t*, int_t, const complex8_t*, complex8_t*, int_t);
/*-------------------------------------------------*/
template <typename T_Scalar>
void CuBlasHandler::trmm(side_t side,
                          uplo_t uplo,
                          op_t trans,
                          diag_t diag,
                          int_t m, int_t n,
                          const T_Scalar* alpha,
                          const T_Scalar* a, int_t lda,
                          const T_Scalar* b, int_t ldb,
                          T_Scalar* c, int_t ldc)
{
    cublas::trmm(handle(),
                 cublas::cla3pSide2cublasSide(side),
                 cublas::cla3pUplo2cublasUplo(uplo),
                 cublas::cla3pOp2cublasOp(trans),
                 cublas::cla3pDiag2cublasDiag(diag),
                 m, n,
                 alpha,
                 a, lda,
                 b, ldb,
                 c, ldc);
}
/*-------------------------------------------------*/
template void CuBlasHandler::trmm<real_t>(side_t, uplo_t, op_t, diag_t, int_t, int_t, const real_t*, const real_t*, int_t, const real_t*, int_t, real_t*, int_t);
template void CuBlasHandler::trmm<real4_t>(side_t, uplo_t, op_t, diag_t, int_t, int_t, const real4_t*, const real4_t*, int_t, const real4_t*, int_t, real4_t*, int_t);
template void CuBlasHandler::trmm<complex_t>(side_t, uplo_t, op_t, diag_t, int_t, int_t, const complex_t*, const complex_t*, int_t, const complex_t*, int_t, complex_t*, int_t);
template void CuBlasHandler::trmm<complex8_t>(side_t, uplo_t, op_t, diag_t, int_t, int_t, const complex8_t*, const complex8_t*, int_t, const complex8_t*, int_t, complex8_t*, int_t);
/*-------------------------------------------------*/
template <typename T_Scalar>
void CuBlasHandler::trsm(side_t side,
                          uplo_t uplo,
                          op_t trans,
                          diag_t diag,
                          int_t m, int_t n,
                          const T_Scalar* alpha,
                          const T_Scalar* a, int_t lda,
                          T_Scalar* b, int_t ldb)
{
    cublas::trsm(handle(),
                 cublas::cla3pSide2cublasSide(side),
                 cublas::cla3pUplo2cublasUplo(uplo),
                 cublas::cla3pOp2cublasOp(trans),
                 cublas::cla3pDiag2cublasDiag(diag),
                 m, n,
                 alpha,
                 a, lda,
                 b, ldb);
}
/*-------------------------------------------------*/
template void CuBlasHandler::trsm<real_t>(side_t, uplo_t, op_t, diag_t, int_t, int_t, const real_t*, const real_t*, int_t, real_t*, int_t);
template void CuBlasHandler::trsm<real4_t>(side_t, uplo_t, op_t, diag_t, int_t, int_t, const real4_t*, const real4_t*, int_t, real4_t*, int_t);
template void CuBlasHandler::trsm<complex_t>(side_t, uplo_t, op_t, diag_t, int_t, int_t, const complex_t*, const complex_t*, int_t, complex_t*, int_t);
template void CuBlasHandler::trsm<complex8_t>(side_t, uplo_t, op_t, diag_t, int_t, int_t, const complex8_t*, const complex8_t*, int_t, complex8_t*, int_t);
/*-------------------------------------------------*/
template <typename T_Scalar>
void CuBlasHandler::syrk(uplo_t uplo,
                          op_t trans,
                          int_t n, int_t k,
                          const T_Scalar* alpha,
                          const T_Scalar* a, int_t lda,
                          const T_Scalar* beta,
                          T_Scalar* c, int_t ldc)
{
    cublas::syrk(handle(),
                 cublas::cla3pUplo2cublasUplo(uplo),
                 cublas::cla3pOp2cublasOp(trans),
                 n, k,
                 alpha,
                 a, lda,
                 beta,
                 c, ldc);
}
/*-------------------------------------------------*/
template void CuBlasHandler::syrk<real_t>(uplo_t, op_t, int_t, int_t, const real_t*, const real_t*, int_t, const real_t*, real_t*, int_t);
template void CuBlasHandler::syrk<real4_t>(uplo_t, op_t, int_t, int_t, const real4_t*, const real4_t*, int_t, const real4_t*, real4_t*, int_t);
template void CuBlasHandler::syrk<complex_t>(uplo_t, op_t, int_t, int_t, const complex_t*, const complex_t*, int_t, const complex_t*, complex_t*, int_t);
template void CuBlasHandler::syrk<complex8_t>(uplo_t, op_t, int_t, int_t, const complex8_t*, const complex8_t*, int_t, const complex8_t*, complex8_t*, int_t);
/*-------------------------------------------------*/
template <typename T_Scalar>
void CuBlasHandler::herk(uplo_t uplo,
                          op_t trans,
                          int_t n, int_t k,
                          const typename TypeTraits<T_Scalar>::real_type* alpha,
                          const T_Scalar* a, int_t lda,
                          const typename TypeTraits<T_Scalar>::real_type* beta,
                          T_Scalar* c, int_t ldc)
{
    cublas::herk(handle(),
                 cublas::cla3pUplo2cublasUplo(uplo),
                 cublas::cla3pOp2cublasOp(trans),
                 n, k,
                 alpha,
                 a, lda,
                 beta,
                 c, ldc);
}
/*-------------------------------------------------*/
template void CuBlasHandler::herk<real_t>(uplo_t, op_t, int_t, int_t, const typename TypeTraits<real_t>::real_type*, const real_t*, int_t, const typename TypeTraits<real_t>::real_type*, real_t*, int_t);
template void CuBlasHandler::herk<real4_t>(uplo_t, op_t, int_t, int_t, const typename TypeTraits<real4_t>::real_type*, const real4_t*, int_t, const typename TypeTraits<real4_t>::real_type*, real4_t*, int_t);
template void CuBlasHandler::herk<complex_t>(uplo_t, op_t, int_t, int_t, const typename TypeTraits<complex_t>::real_type*, const complex_t*, int_t, const typename TypeTraits<complex_t>::real_type*, complex_t*, int_t);
template void CuBlasHandler::herk<complex8_t>(uplo_t, op_t, int_t, int_t, const typename TypeTraits<complex8_t>::real_type*, const complex8_t*, int_t, const typename TypeTraits<complex8_t>::real_type*, complex8_t*, int_t);
/*-------------------------------------------------*/
template <typename T_Scalar>
void CuBlasHandler::syrkx(uplo_t uplo,
                           op_t trans,
                           int_t n, int_t k,
                           const T_Scalar* alpha,
                           const T_Scalar* a, int_t lda,
                           const T_Scalar* b, int_t ldb,
                           const T_Scalar* beta,
                           T_Scalar* c, int_t ldc)
{
    cublas::syrkx(handle(),
                  cublas::cla3pUplo2cublasUplo(uplo),
                  cublas::cla3pOp2cublasOp(trans),
                  n, k,
                  alpha,
                  a, lda,
                  b, ldb,
                  beta,
                  c, ldc);
}
/*-------------------------------------------------*/
template void CuBlasHandler::syrkx<real_t>(uplo_t, op_t, int_t, int_t, const real_t*, const real_t*, int_t, const real_t*, int_t, const real_t*, real_t*, int_t);
template void CuBlasHandler::syrkx<real4_t>(uplo_t, op_t, int_t, int_t, const real4_t*, const real4_t*, int_t, const real4_t*, int_t, const real4_t*, real4_t*, int_t);
template void CuBlasHandler::syrkx<complex_t>(uplo_t, op_t, int_t, int_t, const complex_t*, const complex_t*, int_t, const complex_t*, int_t, const complex_t*, complex_t*, int_t);
template void CuBlasHandler::syrkx<complex8_t>(uplo_t, op_t, int_t, int_t, const complex8_t*, const complex8_t*, int_t, const complex8_t*, int_t, const complex8_t*, complex8_t*, int_t);
/*-------------------------------------------------*/
template <typename T_Scalar>
void CuBlasHandler::herkx(uplo_t uplo,
                           op_t trans,
                           int_t n, int_t k,
                           const T_Scalar* alpha,
                           const T_Scalar* a, int_t lda,
                           const T_Scalar* b, int_t ldb,
                           const typename TypeTraits<T_Scalar>::real_type* beta,
                           T_Scalar* c, int_t ldc)
{
    cublas::herkx(handle(),
                  cublas::cla3pUplo2cublasUplo(uplo),
                  cublas::cla3pOp2cublasOp(trans),
                  n, k,
                  alpha,
                  a, lda,
                  b, ldb,
                  beta,
                  c, ldc);
}
/*-------------------------------------------------*/
template void CuBlasHandler::herkx<real_t>(uplo_t, op_t, int_t, int_t, const real_t*, const real_t*, int_t, const real_t*, int_t, const typename TypeTraits<real_t>::real_type*, real_t*, int_t);
template void CuBlasHandler::herkx<real4_t>(uplo_t, op_t, int_t, int_t, const real4_t*, const real4_t*, int_t, const real4_t*, int_t, const typename TypeTraits<real4_t>::real_type*, real4_t*, int_t);
template void CuBlasHandler::herkx<complex_t>(uplo_t, op_t, int_t, int_t, const complex_t*, const complex_t*, int_t, const complex_t*, int_t, const typename TypeTraits<complex_t>::real_type*, complex_t*, int_t);
template void CuBlasHandler::herkx<complex8_t>(uplo_t, op_t, int_t, int_t, const complex8_t*, const complex8_t*, int_t, const complex8_t*, int_t, const typename TypeTraits<complex8_t>::real_type*, complex8_t*, int_t);
/*-------------------------------------------------*/
} // namespace culite
/*-------------------------------------------------*/
