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

#ifndef CULITE_CUBLAS_PROXY_HPP_
#define CULITE_CUBLAS_PROXY_HPP_

/**
 * @file
 */

#include <cublas_v2.h>

#include <cla3p/types/enums.hpp>

#include "culite/types/integer.hpp"
#include "culite/types/scalar.hpp"

/*-------------------------------------------------*/
namespace culite {
namespace cublas {
/*-------------------------------------------------*/

::cla3p::op_t cublasOp2cla3pOp(cublasOperation_t op);
cublasOperation_t cla3pOp2cublasOp(::cla3p::op_t op);

::cla3p::side_t cublasSide2cla3pSide(cublasSideMode_t side);
cublasSideMode_t cla3pSide2cublasSide(::cla3p::side_t side);

::cla3p::uplo_t cublasUplo2cla3pUplo(cublasFillMode_t uplo);
cublasFillMode_t cla3pUplo2cublasUplo(::cla3p::uplo_t uplo);

::cla3p::diag_t cublasDiag2cla3pDiag(cublasDiagType_t diag);
cublasDiagType_t cla3pDiag2cublasDiag(::cla3p::diag_t diag);

/*-------------------------------------------------*/

template <typename T_Scalar>
void VectorH2D(int_t n, const T_Scalar *src, T_Scalar *dest);

template <typename T_Scalar>
void VectorD2H(int_t n, const T_Scalar *src, T_Scalar *dest);

/*-------------------------------------------------*/

template <typename T_Scalar>
void MatrixH2D(int_t nr, int_t nc, 
               const T_Scalar *src, int_t lds, 
               T_Scalar *dest, int_t ldd);

template <typename T_Scalar>
void MatrixD2H(int_t nr, int_t nc, 
               const T_Scalar *src, int_t lds, 
               T_Scalar *dest, int_t ldd);

/*------------------ Level 1 ----------------------*/

#define iamax_macro(typein) \
void iamax(cublasHandle_t handle, int_t n, const typein *x, int_t incx, int_t *result)
iamax_macro(real_t);
iamax_macro(real4_t);
iamax_macro(complex_t);
iamax_macro(complex8_t);
#undef iamax_macro

template <typename T_Scalar>
int_t iamax(cublasHandle_t handle, int_t n, const T_Scalar *x, int_t incx)
{
	int_t result = 0;
	iamax(handle, n, x, incx, &result);
	return result;
}

/*-------------------------------------------------*/

#define asum_macro(typein) \
void asum(cublasHandle_t handle, int_t n, const typein *x, int_t incx, TypeTraits<typein>::real_type *result)
asum_macro(real_t);
asum_macro(real4_t);
asum_macro(complex_t);
asum_macro(complex8_t);
#undef asum_macro

template <typename T_Scalar>
typename TypeTraits<T_Scalar>::real_type asum(cublasHandle_t handle, int_t n, const T_Scalar *x, int_t incx)
{
	typename TypeTraits<T_Scalar>::real_type result = 0;
	asum(handle, n, x, incx, &result);
	return result;
}

/*-------------------------------------------------*/

#define copy_macro(typein) \
void copy(cublasHandle_t handle, int_t n, const typein *x, int_t incx, typein *y, int_t incy)
copy_macro(real_t);
copy_macro(real4_t);
copy_macro(complex_t);
copy_macro(complex8_t);
#undef copy_macro

/*-------------------------------------------------*/

#define swap_macro(typein) \
void swap(cublasHandle_t handle, int_t n, typein *x, int_t incx, typein *y, int_t incy)
swap_macro(real_t);
swap_macro(real4_t);
swap_macro(complex_t);
swap_macro(complex8_t);
#undef swap_macro

/*-------------------------------------------------*/

#define scal_macro(typein) \
void scal(cublasHandle_t handle, int_t n, const typein *alpha, typein *x, int_t incx)
scal_macro(real_t);
scal_macro(real4_t);
scal_macro(complex_t);
scal_macro(complex8_t);
#undef scal_macro

template <typename T_Scalar>
void scal(cublasHandle_t handle, int_t n, T_Scalar alpha, T_Scalar *x, int_t incx)
{
	scal(handle, n, &alpha, x, incx);
}

#define scal_macro(typein) \
void scal(cublasHandle_t handle, int_t n, const TypeTraits<typein>::real_type *alpha, typein *x, int_t incx)
scal_macro(complex_t);
scal_macro(complex8_t);
#undef scal_macro

template <typename T_Scalar>
void scal(cublasHandle_t handle, int_t n, typename TypeTraits<T_Scalar>::real_type alpha, T_Scalar *x, int_t incx)
{
	scal(handle, n, &alpha, x, incx);
}

/*-------------------------------------------------*/

#define nrm2_macro(typein) \
void nrm2(cublasHandle_t handle, int_t n, const typein *x, int_t incx, TypeTraits<typein>::real_type *result)
nrm2_macro(real_t);
nrm2_macro(real4_t);
nrm2_macro(complex_t);
nrm2_macro(complex8_t);
#undef nrm2_macro

template <typename T_Scalar>
typename TypeTraits<T_Scalar>::real_type nrm2(cublasHandle_t handle, int_t n, const T_Scalar *x, int_t incx)
{
	typename TypeTraits<T_Scalar>::real_type result = 0;
	nrm2(handle, n, x, incx, &result);
	return result;
}

/*-------------------------------------------------*/

#define dot_macro(typeio) \
void dot(cublasHandle_t handle, int_t n, const typeio *x, int_t incx, const typeio *y, int_t incy, typeio *result)
dot_macro(real_t);
dot_macro(real4_t);
dot_macro(complex_t);
dot_macro(complex8_t);
#undef dot_macro

template <typename T_Scalar>
T_Scalar dot(cublasHandle_t handle, int_t n, const T_Scalar *x, int_t incx, const T_Scalar *y, int_t incy)
{
	T_Scalar result = 0;
	dot(handle, n, x, incx, y, incy, &result);
	return result;
}

#define dotc_macro(typeio) \
void dotc(cublasHandle_t handle, int_t n, const typeio *x, int_t incx, const typeio *y, int_t incy, typeio *result)
dotc_macro(real_t);
dotc_macro(real4_t);
dotc_macro(complex_t);
dotc_macro(complex8_t);
#undef dotc_macro

template <typename T_Scalar>
T_Scalar dotc(cublasHandle_t handle, int_t n, const T_Scalar *x, int_t incx, const T_Scalar *y, int_t incy)
{
	T_Scalar result = 0;
	dotc(handle, n, x, incx, y, incy, &result);
	return result;
}

/*-------------------------------------------------*/

#define axpy_macro(typein) \
void axpy(cublasHandle_t handle, int_t n, const typein *alpha, const typein *x, int_t incx, typein *y, int_t incy)
axpy_macro(real_t);
axpy_macro(real4_t);
axpy_macro(complex_t);
axpy_macro(complex8_t);
#undef axpy_macro

template <typename T_Scalar>
void axpy(cublasHandle_t handle, int_t n, T_Scalar alpha, const T_Scalar *x, int_t incx, T_Scalar *y, int_t incy)
{
	axpy(handle, n, &alpha, x, incx, y, incy);
}

/*------------------ Level 2 ----------------------*/

#define geam_macro(typein) \
void geam(cublasHandle_t handle, \
	      cublasOperation_t transa, \
		  cublasOperation_t transb, \
		  int_t m, int_t n, \
		  const typein *alpha, const typein *a, int_t lda, \
		  const typein *beta, const typein *b, int_t ldb, \
		  typein *c, int_t ldc)
geam_macro(real_t);
geam_macro(real4_t);
geam_macro(complex_t);
geam_macro(complex8_t);
#undef geam_macro

/*-------------------------------------------------*/

#define dgmm_macro(typein) \
void dgmm(cublasHandle_t handle, \
	      cublasSideMode_t mode, \
          int_t m, int_t n, \
          const typein *a, int_t lda, \
          const typein *x, int_t incx, \
          typein *c, int_t ldc)
dgmm_macro(real_t);
dgmm_macro(real4_t);
dgmm_macro(complex_t);
dgmm_macro(complex8_t);
#undef dgmm_macro

/*-------------------------------------------------*/

#define gemv_macro(typein) \
void gemv(cublasHandle_t handle, \
	      cublasOperation_t trans, \
		  int_t m, int_t n, const typein* alpha, \
		  const typein *a, int_t lda, \
		  const typein *x, int_t incx, \
		  const typein* beta, typein *y, int_t incy)
gemv_macro(real_t);
gemv_macro(real4_t);
gemv_macro(complex_t);
gemv_macro(complex8_t);
#undef gemv_macro

/*-------------------------------------------------*/

#define symv_macro(typein) \
void symv(cublasHandle_t handle, \
	      cublasFillMode_t uplo, \
		  int_t n, const typein* alpha, \
		  const typein *a, int_t lda, \
		  const typein *x, int_t incx, \
		  const typein* beta, typein *y, int_t incy)
symv_macro(real_t);
symv_macro(real4_t);
symv_macro(complex_t);
symv_macro(complex8_t);
#undef symv_macro

/*-------------------------------------------------*/

#define hemv_macro(typein) \
void hemv(cublasHandle_t handle, \
	      cublasFillMode_t uplo, \
		  int_t n, const typein* alpha, \
		  const typein *a, int_t lda, \
		  const typein *x, int_t incx, \
		  const typein* beta, typein *y, int_t incy)
hemv_macro(real_t);
hemv_macro(real4_t);
hemv_macro(complex_t);
hemv_macro(complex8_t);
#undef hemv_macro

/*-------------------------------------------------*/

#define trmv_macro(typein) \
void trmv(cublasHandle_t handle, \
	      cublasFillMode_t uplo, \
		  cublasOperation_t trans, \
		  cublasDiagType_t diag, \
		  int_t n, const typein *a, int_t lda, \
		  typein *x, int_t incx)
trmv_macro(real_t);
trmv_macro(real4_t);
trmv_macro(complex_t);
trmv_macro(complex8_t);
#undef trmv_macro

/*-------------------------------------------------*/

#define trsv_macro(typein) \
void trsv(cublasHandle_t handle, \
	      cublasFillMode_t uplo, \
		  cublasOperation_t trans, \
		  cublasDiagType_t diag, \
		  int_t n, const typein *a, int_t lda, \
		  typein *x, int_t incx)
trsv_macro(real_t);
trsv_macro(real4_t);
trsv_macro(complex_t);
trsv_macro(complex8_t);
#undef trsv_macro

/*------------------ Level 3 ----------------------*/

#define gemm_macro(typein) \
void gemm(cublasHandle_t handle, \
          cublasOperation_t transa, \
		  cublasOperation_t transb, \
          int_t m, int_t n, int_t k, \
          const typein* alpha, \
          const typein* a, int_t lda, \
          const typein* b, int_t ldb, \
          const typein* beta, \
          typein* c, int_t ldc)
gemm_macro(real_t);
gemm_macro(real4_t);
gemm_macro(complex_t);
gemm_macro(complex8_t);
#undef gemm_macro

/*-------------------------------------------------*/

#define symm_macro(typein) \
void symm(cublasHandle_t handle, \
          cublasSideMode_t side, \
		  cublasFillMode_t uplo, \
          int_t m, int_t n, \
          const typein* alpha, \
          const typein* a, int_t lda, \
          const typein* b, int_t ldb, \
          const typein* beta, \
          typein* c, int_t ldc)
symm_macro(real_t);
symm_macro(real4_t);
symm_macro(complex_t);
symm_macro(complex8_t);
#undef symm_macro

/*-------------------------------------------------*/

#define hemm_macro(typein) \
void hemm(cublasHandle_t handle, \
          cublasSideMode_t side, \
		  cublasFillMode_t uplo, \
          int_t m, int_t n, \
          const typein* alpha, \
          const typein* a, int_t lda, \
          const typein* b, int_t ldb, \
          const typein* beta, \
          typein* c, int_t ldc)
hemm_macro(real_t);
hemm_macro(real4_t);
hemm_macro(complex_t);
hemm_macro(complex8_t);
#undef hemm_macro

/*-------------------------------------------------*/

#define trmm_macro(typein) \
void trmm(cublasHandle_t handle, \
          cublasSideMode_t side, \
		  cublasFillMode_t uplo, \
		  cublasOperation_t trans, \
		  cublasDiagType_t diag , \
          int_t m, int_t n, \
          const typein* alpha, \
          const typein* a, int_t lda, \
          const typein* b, int_t ldb, \
          typein* c, int_t ldc)
trmm_macro(real_t);
trmm_macro(real4_t);
trmm_macro(complex_t);
trmm_macro(complex8_t);
#undef trmm_macro

/*-------------------------------------------------*/

#define trsm_macro(typein) \
void trsm(cublasHandle_t handle, \
          cublasSideMode_t side, \
		  cublasFillMode_t uplo, \
		  cublasOperation_t trans, \
		  cublasDiagType_t diag , \
          int_t m, int_t n, \
          const typein* alpha, \
          const typein* a, int_t lda, \
          typein* b, int_t ldb)
trsm_macro(real_t);
trsm_macro(real4_t);
trsm_macro(complex_t);
trsm_macro(complex8_t);
#undef trsm_macro

/*-------------------------------------------------*/

#define syrk_macro(typein) \
void syrk(cublasHandle_t handle, \
          cublasFillMode_t uplo, \
		  cublasOperation_t trans, \
          int_t n, int_t k, \
          const typein* alpha, \
          const typein* a, int_t lda, \
          const typein* beta, \
          typein* c, int_t ldc)
syrk_macro(real_t);
syrk_macro(real4_t);
syrk_macro(complex_t);
syrk_macro(complex8_t);
#undef syrk_macro

/*-------------------------------------------------*/

#define herk_macro(typein) \
void herk(cublasHandle_t handle, \
          cublasFillMode_t uplo, \
		  cublasOperation_t trans, \
          int_t n, int_t k, \
          const typename TypeTraits<typein>::real_type* alpha, \
          const typein* a, int_t lda, \
          const typename TypeTraits<typein>::real_type* beta, \
          typein* c, int_t ldc)
herk_macro(real_t);
herk_macro(real4_t);
herk_macro(complex_t);
herk_macro(complex8_t);
#undef herk_macro

/*-------------------------------------------------*/

#define syrkx_macro(typein) \
void syrkx(cublasHandle_t handle, \
          cublasFillMode_t uplo, \
		  cublasOperation_t trans, \
          int_t n, int_t k, \
          const typein* alpha, \
          const typein* a, int_t lda, \
          const typein* b, int_t ldb, \
          const typein* beta, \
          typein* c, int_t ldc)
syrkx_macro(real_t);
syrkx_macro(real4_t);
syrkx_macro(complex_t);
syrkx_macro(complex8_t);
#undef syrkx_macro

/*-------------------------------------------------*/

#define herkx_macro(typein) \
void herkx(cublasHandle_t handle, \
          cublasFillMode_t uplo, \
		  cublasOperation_t trans, \
          int_t n, int_t k, \
          const typein* alpha, \
          const typein* a, int_t lda, \
          const typein* b, int_t ldb, \
          const typename TypeTraits<typein>::real_type* beta, \
          typein* c, int_t ldc)
herkx_macro(real_t);
herkx_macro(real4_t);
herkx_macro(complex_t);
herkx_macro(complex8_t);
#undef herkx_macro

/*-------------------------------------------------*/
} // namespace cublas
} // namespace culite
/*-------------------------------------------------*/

#endif // CULITE_CUBLAS_PROXY_HPP_
