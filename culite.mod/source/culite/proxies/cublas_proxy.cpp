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
#include "culite/proxies/cublas_proxy.hpp"

// system

// 3rd

// culite
#include "culite/error/cuda.hpp"
#include "culite/error/exceptions.hpp"

/*-------------------------------------------------*/
#if defined(CULITE_I64)
#define cublas_func_name(name) cublas##name##_64
#else
#define cublas_func_name(name) cublas##name
#endif
/*-------------------------------------------------*/
namespace culite {
namespace cublas {
/*-------------------------------------------------*/
::cla3p::op_t cublasOp2cla3pOp(cublasOperation_t op)
{
    switch(op) {
        case CUBLAS_OP_N: return ::cla3p::op_t::N;
        case CUBLAS_OP_T: return ::cla3p::op_t::T;
        case CUBLAS_OP_C: return ::cla3p::op_t::C;
        default:
            throw err::CudaException("Invalid cublasOperation_t value.");
    }
}
/*-------------------------------------------------*/
cublasOperation_t cla3pOp2cublasOp(::cla3p::op_t op)
{
    switch(op) {
        case ::cla3p::op_t::N: return CUBLAS_OP_N;
        case ::cla3p::op_t::T: return CUBLAS_OP_T;
        case ::cla3p::op_t::C: return CUBLAS_OP_C;
        default:
            throw err::CudaException("Invalid (::cla3p::op_t) value.");
    }
}
/*-------------------------------------------------*/
::cla3p::side_t cublasSide2cla3pSide(cublasSideMode_t side)
{
	switch(side) {
		case CUBLAS_SIDE_LEFT:  return ::cla3p::side_t::Left;
		case CUBLAS_SIDE_RIGHT: return ::cla3p::side_t::Right;
		default:
			throw err::CudaException("Invalid cublasSideMode_t value.");
	}
}
/*-------------------------------------------------*/
cublasSideMode_t cla3pSide2cublasSide(::cla3p::side_t side)
{
	switch(side) {
		case ::cla3p::side_t::Left: return CUBLAS_SIDE_LEFT;
		case ::cla3p::side_t::Right: return CUBLAS_SIDE_RIGHT;
		default:
			throw err::CudaException("Invalid (::cla3p::side_t) value.");
	}
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void VectorH2D(int_t n, const typename TypeTraits<T_Scalar>::host_type *src, T_Scalar *dest)
{
    cublasStatus_t cublasStatus = cublas_func_name(SetVector)(n, sizeof(T_Scalar), src, 1, dest, 1);
    err::check_cublas(cublasStatus);
}
/*-------------------------------------------------*/
#define vector_h2d_macro(typein) \
template void VectorH2D<typein>(int_t, const typename TypeTraits<typein>::host_type*, typein*);
vector_h2d_macro(int_t);
vector_h2d_macro(real_t);
vector_h2d_macro(real4_t);
vector_h2d_macro(complex_t);
vector_h2d_macro(complex8_t);
#undef vector_h2d_macro
/*-------------------------------------------------*/
template <typename T_Scalar>
void VectorD2H(int_t n, const T_Scalar *src, typename TypeTraits<T_Scalar>::host_type *dest)
{
    cublasStatus_t cublasStatus = cublas_func_name(GetVector)(n, sizeof(T_Scalar), src, 1, dest, 1);
    err::check_cublas(cublasStatus);
}
/*-------------------------------------------------*/
#define vector_d2h_macro(typein) \
template void VectorD2H<typein>(int_t, const typein*, typename TypeTraits<typein>::host_type*);
vector_d2h_macro(int_t);
vector_d2h_macro(real_t);
vector_d2h_macro(real4_t);
vector_d2h_macro(complex_t);
vector_d2h_macro(complex8_t);
#undef vector_d2h_macro
/*-------------------------------------------------*/
template <typename T_Scalar>
void MatrixH2D(int_t nr, int_t nc, 
               const typename TypeTraits<T_Scalar>::host_type *src, int_t lds, 
               T_Scalar *dest, int_t ldd)
{
    cublasStatus_t cublasStatus = cublas_func_name(SetMatrix)(nr, nc, sizeof(T_Scalar), src, lds, dest, ldd);
    err::check_cublas(cublasStatus);
}
/*-------------------------------------------------*/
#define matrix_h2d_macro(typein) \
template void MatrixH2D<typein>(int_t, int_t, const typename TypeTraits<typein>::host_type*, int_t, typein*, int_t);
matrix_h2d_macro(real_t);
matrix_h2d_macro(real4_t);
matrix_h2d_macro(complex_t);
matrix_h2d_macro(complex8_t);
#undef matrix_h2d_macro
/*-------------------------------------------------*/
template <typename T_Scalar>
void MatrixD2H(int_t nr, int_t nc, 
               const T_Scalar *src, int_t lds, 
               typename TypeTraits<T_Scalar>::host_type *dest, int_t ldd)
{
    cublasStatus_t cublasStatus = cublas_func_name(GetMatrix)(nr, nc, sizeof(T_Scalar), src, lds, dest, ldd);
    err::check_cublas(cublasStatus);
}
/*-------------------------------------------------*/
#define matrix_d2h_macro(typein) \
template void MatrixD2H<typein>(int_t, int_t, const typein*, int_t, typename TypeTraits<typein>::host_type*, int_t);
matrix_d2h_macro(real_t);
matrix_d2h_macro(real4_t);
matrix_d2h_macro(complex_t);
matrix_d2h_macro(complex8_t);
#undef matrix_d2h_macro
/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/
#define iamax_macro(typein, prefix) \
void iamax(cublasHandle_t handle, int_t n, const typein *x, int_t incx, int_t *result) \
{ \
	cublasStatus_t cublasStatus = cublas_func_name(prefix##amax)(handle, n, x, incx, result); \
	err::check_cublas(cublasStatus); \
}
iamax_macro(real_t    , Id);
iamax_macro(real4_t   , Is);
iamax_macro(complex_t , Iz);
iamax_macro(complex8_t, Ic);
#undef iamax_macro
/*-------------------------------------------------*/
#define asum_macro(typein, prefix) \
void asum(cublasHandle_t handle, int_t n, const typein *x, int_t incx, TypeTraits<typein>::real_type *result) \
{ \
	cublasStatus_t cublasStatus = cublas_func_name(prefix##asum)(handle, n, x, incx, result); \
	err::check_cublas(cublasStatus); \
}
asum_macro(real_t    ,  D);
asum_macro(real4_t   ,  S);
asum_macro(complex_t , Dz);
asum_macro(complex8_t, Sc);
#undef asum_macro
/*-------------------------------------------------*/
#define copy_macro(typein, prefix) \
void copy(cublasHandle_t handle, int_t n, const typein *x, int_t incx, typein *y, int_t incy) \
{ \
	cublasStatus_t cublasStatus = cublas_func_name(prefix##copy)(handle, n, x, incx, y, incy); \
	err::check_cublas(cublasStatus); \
}
copy_macro(real_t    , D)
copy_macro(real4_t   , S)
copy_macro(complex_t , Z)
copy_macro(complex8_t, C)
#undef copy_macro
/*-------------------------------------------------*/
#define swap_macro(typein, prefix) \
void swap(cublasHandle_t handle, int_t n, typein *x, int_t incx, typein *y, int_t incy) \
{ \
	cublasStatus_t cublasStatus = cublas_func_name(prefix##swap)(handle, n, x, incx, y, incy); \
	err::check_cublas(cublasStatus); \
}
swap_macro(real_t    , D)
swap_macro(real4_t   , S)
swap_macro(complex_t , Z)
swap_macro(complex8_t, C)
#undef swap_macro
/*-------------------------------------------------*/
#define scal_macro(typein, prefix) \
void scal(cublasHandle_t handle, int_t n, const typein *alpha, typein *x, int_t incx) \
{ \
	cublasStatus_t cublasStatus = cublas_func_name(prefix##scal)(handle, n, alpha, x, incx); \
	err::check_cublas(cublasStatus); \
}
scal_macro(real_t    , D)
scal_macro(real4_t   , S);
scal_macro(complex_t , Z);
scal_macro(complex8_t, C);
#undef scal_macro
/*-------------------------------------------------*/
#define scal_macro(typein, prefix) \
void scal(cublasHandle_t handle, int_t n, const TypeTraits<typein>::real_type *alpha, typein *x, int_t incx) \
{ \
	cublasStatus_t cublasStatus = cublas_func_name(prefix##scal)(handle, n, alpha, x, incx); \
	err::check_cublas(cublasStatus); \
}
scal_macro(complex_t , Zd);
scal_macro(complex8_t, Cs);
#undef scal_macro
/*-------------------------------------------------*/
#define nrm2_macro(typein, prefix) \
void nrm2(cublasHandle_t handle, int_t n, const typein *x, int_t incx, TypeTraits<typein>::real_type *result) \
{ \
	cublasStatus_t cublasStatus = cublas_func_name(prefix##nrm2)(handle, n, x, incx, result); \
	err::check_cublas(cublasStatus); \
}
nrm2_macro(real_t    ,  D)
nrm2_macro(real4_t   ,  S)
nrm2_macro(complex_t , Dz)
nrm2_macro(complex8_t, Sc)
#undef nrm2_macro
/*-------------------------------------------------*/
#define dot_macro(typeio, prefix, suffix) \
void dot(cublasHandle_t handle, int_t n, const typeio *x, int_t incx, const typeio *y, int_t incy, typeio *result) \
{ \
	cublasStatus_t cublasStatus = cublas_func_name(prefix##dot##suffix)(handle, n, x, incx, y, incy, result); \
	err::check_cublas(cublasStatus); \
}
dot_macro(real_t    , D,  )
dot_macro(real4_t   , S,  )
dot_macro(complex_t , Z, u)
dot_macro(complex8_t, C, u)
#undef dot_macro
/*-------------------------------------------------*/
#define dot_macro(typeio, prefix, suffix) \
void dotc(cublasHandle_t handle, int_t n, const typeio *x, int_t incx, const typeio *y, int_t incy, typeio *result) \
{ \
	cublasStatus_t cublasStatus = cublas_func_name(prefix##dot##suffix)(handle, n, x, incx, y, incy, result); \
	err::check_cublas(cublasStatus); \
}
dot_macro(real_t    , D,  )
dot_macro(real4_t   , S,  )
dot_macro(complex_t , Z, c)
dot_macro(complex8_t, C, c)
#undef dot_macro
/*-------------------------------------------------*/
#define axpy_macro(typein, prefix) \
void axpy(cublasHandle_t handle, int_t n, const typein *alpha, const typein *x, int_t incx, typein *y, int_t incy) \
{ \
	cublasStatus_t cublasStatus = cublas_func_name(prefix##axpy)(handle, n, alpha, x, incx, y, incy); \
	err::check_cublas(cublasStatus); \
}
axpy_macro(real_t    , D)
axpy_macro(real4_t   , S)
axpy_macro(complex_t , Z)
axpy_macro(complex8_t, C)
#undef axpy_macro
/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/
#define geam_macro(typein, prefix) \
void geam(cublasHandle_t handle, \
	      cublasOperation_t transa, \
		  cublasOperation_t transb, \
		  int_t m, int_t n, \
		  const typein *alpha, const typein *a, int_t lda, \
		  const typein *beta, const typein *b, int_t ldb, \
		  typein *c, int_t ldc) \
{ \
	cublasStatus_t cublasStatus = cublas_func_name(prefix##geam)(handle, transa, transb, m, n, alpha, a, lda, beta, b, ldb, c, ldc); \
	err::check_cublas(cublasStatus); \
}
geam_macro(real_t    , D);
geam_macro(real4_t   , S);
geam_macro(complex_t , Z);
geam_macro(complex8_t, C);
#undef geam_macro
/*-------------------------------------------------*/
#define dgmm_macro(typein, prefix) \
void dgmm(cublasHandle_t handle, \
	      cublasSideMode_t mode, \
          int_t m, int_t n, \
          const typein *a, int_t lda, \
          const typein *x, int_t incx, \
          typein *c, int_t ldc) \
{ \
	cublasStatus_t cublasStatus = cublas_func_name(prefix##dgmm)(handle, mode, m, n, a, lda, x, incx, c, ldc); \
	err::check_cublas(cublasStatus); \
}
dgmm_macro(real_t    , D);
dgmm_macro(real4_t   , S);
dgmm_macro(complex_t , Z);
dgmm_macro(complex8_t, C);
#undef dgmm_macro
/*-------------------------------------------------*/




#if 0

/*-------------------------------------------------*/
/*-------------------------------------------------*/

/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/
#define real_ger_macro(typein, prefix, suffix) \
void ger##suffix(int_t m, int_t n, typein alpha, \
		const typein *x, int_t incx, \
		const typein *y, int_t incy, \
		typein *a, int_t lda) \
{ \
	blas_func_name(prefix##ger)(&m, &n, &alpha, x, &incx, y, &incy, a, &lda); \
}
real_ger_macro(real_t , d,  )
real_ger_macro(real4_t, s,  )
real_ger_macro(real_t , d, c)
real_ger_macro(real4_t, s, c)
#undef real_ger_macro
/*-------------------------------------------------*/
#define complex_ger_macro(typein, prefix, suffix, blas_suffix) \
void ger##suffix(int_t m, int_t n, typein alpha, \
    const typein *x, int_t incx, \
    const typein *y, int_t incy, \
    typein *a, int_t lda) \
{ \
	blas_func_name(prefix##ger##blas_suffix)(&m, &n, &alpha, x, &incx, y, &incy, a, &lda); \
}
complex_ger_macro(complex_t , z,  , u)
complex_ger_macro(complex8_t, c,  , u)
complex_ger_macro(complex_t , z, c, c)
complex_ger_macro(complex8_t, c, c, c)
#undef complex_ger_macro
/*-------------------------------------------------*/
#define syr_macro(typein, prefix) \
void syr(char uplo, int_t n, typein alpha, \
		const typein *x, int_t incx, \
		typein *a, int_t lda) \
{ \
	blas_func_name(prefix##syr)(&uplo, &n, &alpha, x, &incx, a, &lda); \
}
syr_macro(real_t    , d)
syr_macro(real4_t   , s)
syr_macro(complex_t , z)
syr_macro(complex8_t, c)
#undef syr_macro
/*-------------------------------------------------*/
#define her_macro(typein, prefix) \
void her(char uplo, int_t n, TypeTraits<typein>::real_type alpha, \
		const typein *x, int_t incx, \
		typein *a, int_t lda) \
{ \
	blas_func_name(prefix##syr)(&uplo, &n, &alpha, x, &incx, a, &lda); \
}
her_macro(real_t , d)
her_macro(real4_t, s)
#undef her_macro
/*-------------------------------------------------*/
#define her_macro(typein, prefix) \
void her(char uplo, int_t n, TypeTraits<typein>::real_type alpha, \
		const typein *x, int_t incx, \
		typein *a, int_t lda) \
{ \
	blas_func_name(prefix##her)(&uplo, &n, &alpha, x, &incx, a, &lda); \
}
her_macro(complex_t , z)
her_macro(complex8_t, c)
#undef her_macro
/*-------------------------------------------------*/
#define gemv_macro(typein, prefix) \
void gemv(char trans, int_t m, int_t n, typein alpha, \
		const typein *a, int_t lda, const typein *x, int_t incx, \
		typein beta, typein *y, int_t incy) \
{ \
	blas_func_name(prefix##gemv)(&trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy); \
}
gemv_macro(real_t    , d)
gemv_macro(real4_t   , s)
gemv_macro(complex_t , z)
gemv_macro(complex8_t, c)
#undef gemv_macro
/*-------------------------------------------------*/
#define symv_macro(typein, prefix) \
void symv(char uplo, int_t n, typein alpha, const typein *a, int_t lda, \
		const typein *x, int_t incx, typein beta, typein *y, int_t incy) \
{ \
	blas_func_name(prefix##symv)(&uplo, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy); \
}
symv_macro(real_t    , d)
symv_macro(real4_t   , s)
symv_macro(complex_t , z)
symv_macro(complex8_t, c)
#undef symv_macro
/*-------------------------------------------------*/
#define hemv_macro(typein, prefix) \
void hemv(char uplo, int_t n, typein alpha, \
		const typein *a, int_t lda, const typein *x, int_t incx, \
		typein beta, typein *y, int_t incy) \
{ \
	blas_func_name(prefix##hemv)(&uplo, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy); \
}
hemv_macro(complex_t , z)
hemv_macro(complex8_t, c)
#undef hemv_macro
/*-------------------------------------------------*/
#define trmv_macro(typein, prefix) \
void trmv(char uplo, char transa, char diag, int_t n, \
		const typein *a, int_t lda, typein *b, int_t incx) \
{ \
	blas_func_name(prefix##trmv)(&uplo, &transa, &diag, &n, a, &lda, b, &incx); \
}
trmv_macro(real_t    , d)
trmv_macro(real4_t   , s)
trmv_macro(complex_t , z)
trmv_macro(complex8_t, c)
#undef trmv_macro
/*-------------------------------------------------*/
#define trsv_macro(typein, prefix) \
void trsv(char uplo, char transa, char diag, int_t n, \
		const typein *a, int_t lda, typein *b, int_t incx) \
{ \
	blas_func_name(prefix##trsv)(&uplo, &transa, &diag, &n, a, &lda, b, &incx); \
}
trsv_macro(real_t    , d)
trsv_macro(real4_t   , s)
trsv_macro(complex_t , z)
trsv_macro(complex8_t, c)
#undef trsv_macro
/*-------------------------------------------------*/
#define gemm_macro(typein, prefix) \
void gemm(char transa, char transb, int_t m, int_t n, int_t k, \
		typein alpha, const typein *a, int_t lda, const typein *b, int_t ldb, \
		typein beta, typein *c, int_t ldc) \
{ \
	blas_func_name(prefix##gemm)(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc); \
}
gemm_macro(real_t    , d)
gemm_macro(real4_t   , s)
gemm_macro(complex_t , z)
gemm_macro(complex8_t, c)
#undef gemm_macro
/*-------------------------------------------------*/
template <typename T_Scalar>
static void gemmt_recursive(char uplo, char transa, char transb, int_t n, int_t k,
		T_Scalar alpha, const T_Scalar *a, int_t lda, const T_Scalar *b, int_t ldb,
		T_Scalar beta, T_Scalar *c, int_t ldc)
{
	constexpr int_t gemmtRecLimit = 64;

	if(n < gemmtRecLimit) {

		for(int_t j = 0; j < n; j++) {

			int_t incb = (transb == 'N' ? 1 : ldb);
			const T_Scalar *bj = (transb == 'N' ? blk::dns::ptrmv(ldb,b,0,j) : blk::dns::ptrmv(ldb,b,j,0));

			if(uplo == 'L') {

				const T_Scalar *al = (transa == 'N' ? blk::dns::ptrmv(lda,a,j,0) : blk::dns::ptrmv(lda,a,0,j));

				int_t ma = (transa == 'N' ? n-j : k  );
				int_t na = (transa == 'N' ? k   : n-j);
				int_t incc = 1;
				T_Scalar *cj = blk::dns::ptrmv(ldc,c,j,j);
				gemv(transa, ma, na, alpha, al, lda, bj, incb, beta, cj, incc);

			} else if(uplo == 'U') {

				const T_Scalar *au = a;

				int_t ma = (transa == 'N' ? j+1 : k  );
				int_t na = (transa == 'N' ? k   : j+1);
				int_t incc = 1;
				T_Scalar *cj = blk::dns::ptrmv(ldc,c,0,j);
				gemv(transa, ma, na, alpha, au, lda, bj, incb, beta, cj, incc);

			} // U/L

		} // j

	} else {

		int_t n1 = n/2;
		int_t n2 = n - n1;

		const T_Scalar *a1 = a;
		const T_Scalar *a2 = (transa == 'N' ? blk::dns::ptrmv(lda,a,n1,0) : blk::dns::ptrmv(lda,a,0,n1));

		const T_Scalar *b1 = b;
		const T_Scalar *b2 = (transb == 'N' ? blk::dns::ptrmv(ldb,b,0,n1) : blk::dns::ptrmv(ldb,b,n1,0));

		T_Scalar *c11 = c;
		T_Scalar *c22 = blk::dns::ptrmv(ldc,c,n1,n1);

		gemmt_recursive(uplo, transa, transb, n1, k, alpha, a1, lda, b1, ldb, beta, c11, ldc);
		gemmt_recursive(uplo, transa, transb, n2, k, alpha, a2, lda, b2, ldb, beta, c22, ldc);

		if(uplo == 'L') {

			T_Scalar *c21 = blk::dns::ptrmv(ldc,c,n1,0);
			gemm(transa, transb, n2, n1, k, alpha, a2, lda, b1, ldb, beta, c21, ldc);

		} else if(uplo == 'U') {

			T_Scalar *c12 = blk::dns::ptrmv(ldc,c,0,n1);
			gemm(transa, transb, n1, n2, k, alpha, a1, lda, b2, ldb, beta, c12, ldc);

		} // U/L

	} // n
}
/*-------------------------------------------------*/
#if defined(CLA3P_INTEL_MKL) || defined(CLA3P_ARMPL)
#define gemmt_macro(typein, prefix) \
void gemmt(char uplo, char transa, char transb, int_t n, int_t k, \
		typein alpha, const typein *a, int_t lda, const typein *b, int_t ldb, \
		typein beta, typein *c, int_t ldc) \
{ \
	blas_func_name(prefix##gemmt)(&uplo, &transa, &transb, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc); \
}
#else
#define gemmt_macro(typein, prefix) \
void gemmt(char uplo, char transa, char transb, int_t n, int_t k, \
		typein alpha, const typein *a, int_t lda, const typein *b, int_t ldb, \
		typein beta, typein *c, int_t ldc) \
{ \
	gemmt_recursive<typein>(uplo, transa, transb, n, k, alpha, a, lda, b, ldb, beta, c, ldc); \
}
#endif
gemmt_macro(real_t    , d)
gemmt_macro(real4_t   , s)
gemmt_macro(complex_t , z)
gemmt_macro(complex8_t, c)
#undef gemmt_macro
/*-------------------------------------------------*/
#define symm_macro(typein, prefix) \
void symm(char side, char uplo, int_t m, int_t n, \
		typein alpha, const typein *a, int_t lda, const typein *b, int_t ldb, \
		typein beta, typein *c, int_t ldc) \
{ \
	blas_func_name(prefix##symm)(&side, &uplo, &m, &n, &alpha, a, &lda, b, &ldb, &beta, c, &ldc); \
}
symm_macro(real_t    , d)
symm_macro(real4_t   , s)
symm_macro(complex_t , z)
symm_macro(complex8_t, c)
#undef symm_macro
/*-------------------------------------------------*/
#define hemm_macro(typein, prefix) \
void hemm(char side, char uplo, int_t m, int_t n, \
		typein alpha, const typein *a, int_t lda, const typein *b, int_t ldb, \
		typein beta, typein *c, int_t ldc) \
{ \
	blas_func_name(prefix##hemm)(&side, &uplo, &m, &n, &alpha, a, &lda, b, &ldb, &beta, c, &ldc); \
}
hemm_macro(complex_t , z)
hemm_macro(complex8_t, c)
#undef hemm_macro
/*-------------------------------------------------*/
#define trmm_macro(typein, prefix) \
void trmm(char side, char uplo, char transa, char diag, \
		int_t m, int_t n, typein alpha, const typein *a, int_t lda, \
		typein *b, int_t ldb) \
{ \
	blas_func_name(prefix##trmm)(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb); \
}
trmm_macro(real_t    , d)
trmm_macro(real4_t   , s)
trmm_macro(complex_t , z)
trmm_macro(complex8_t, c)
#undef trmm_macro
/*-------------------------------------------------*/
#define trsm_macro(typein, prefix) \
void trsm(char side, char uplo, char transa, char diag, \
		int_t m, int_t n, typein alpha, const typein *a, int_t lda, \
		typein *b, int_t ldb) \
{ \
	blas_func_name(prefix##trsm)(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb); \
}
trsm_macro(real_t    , d)
trsm_macro(real4_t   , s)
trsm_macro(complex_t , z)
trsm_macro(complex8_t, c)
#undef trsm_macro


#endif // 0

/*-------------------------------------------------*/
} // namespace cublas
} // namespace culite
/*-------------------------------------------------*/
#undef cublas_func_name
/*-------------------------------------------------*/
