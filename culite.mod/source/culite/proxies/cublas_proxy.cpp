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
::cla3p::uplo_t cublasUplo2cla3pUplo(cublasFillMode_t uplo)
{
    switch(uplo) {
        case CUBLAS_FILL_MODE_UPPER: return ::cla3p::uplo_t::Upper;
        case CUBLAS_FILL_MODE_LOWER: return ::cla3p::uplo_t::Lower;
        case CUBLAS_FILL_MODE_FULL: return ::cla3p::uplo_t::Full;
        default:
            throw err::CudaException("Invalid cublasFillMode_t value.");
    }
}
/*-------------------------------------------------*/
cublasFillMode_t cla3pUplo2cublasUplo(::cla3p::uplo_t uplo)
{
    switch(uplo) {
        case ::cla3p::uplo_t::Upper: return CUBLAS_FILL_MODE_UPPER;
        case ::cla3p::uplo_t::Lower: return CUBLAS_FILL_MODE_LOWER;
        case ::cla3p::uplo_t::Full: return CUBLAS_FILL_MODE_FULL;
        default:
            throw err::CudaException("Invalid (::cla3p::uplo_t) value.");
    }
}
/*-------------------------------------------------*/
::cla3p::diag_t cublasDiag2cla3pDiag(cublasDiagType_t diag)
{
    switch(diag) {
        case CUBLAS_DIAG_UNIT: return ::cla3p::diag_t::Unit;
        case CUBLAS_DIAG_NON_UNIT: return ::cla3p::diag_t::NonUnit;
        default:
            throw err::CudaException("Invalid cublasDiagType_t value.");
    }
}
/*-------------------------------------------------*/
cublasDiagType_t cla3pDiag2cublasDiag(::cla3p::diag_t diag)
{
    switch(diag) {
        case ::cla3p::diag_t::Unit: return CUBLAS_DIAG_UNIT;
        case ::cla3p::diag_t::NonUnit: return CUBLAS_DIAG_NON_UNIT;
        default:
            throw err::CudaException("Invalid (::cla3p::diag_t) value.");
    }
}
/*-------------------------------------------------*/

template <typename T_Scalar>
void VectorH2D(int_t n, const T_Scalar *src, T_Scalar *dest)
{
    cublasStatus_t cublasStatus = cublas_func_name(SetVector)(n, sizeof(T_Scalar), src, 1, dest, 1);
    err::check_cublas(cublasStatus);
}
/*-------------------------------------------------*/
#define vector_h2d_macro(typein) \
template void VectorH2D<typein>(int_t, const typein*, typein*);
vector_h2d_macro(int_t);
vector_h2d_macro(real_t);
vector_h2d_macro(real4_t);
vector_h2d_macro(complex_t);
vector_h2d_macro(complex8_t);
#undef vector_h2d_macro
/*-------------------------------------------------*/
template <typename T_Scalar>
void VectorD2H(int_t n, const T_Scalar *src, T_Scalar *dest)
{
    cublasStatus_t cublasStatus = cublas_func_name(GetVector)(n, sizeof(T_Scalar), src, 1, dest, 1);
    err::check_cublas(cublasStatus);
}
/*-------------------------------------------------*/
#define vector_d2h_macro(typein) \
template void VectorD2H<typein>(int_t, const typein*, typein*);
vector_d2h_macro(int_t);
vector_d2h_macro(real_t);
vector_d2h_macro(real4_t);
vector_d2h_macro(complex_t);
vector_d2h_macro(complex8_t);
#undef vector_d2h_macro
/*-------------------------------------------------*/
template <typename T_Scalar>
void MatrixH2D(int_t nr, int_t nc, 
               const T_Scalar *src, int_t lds, 
               T_Scalar *dest, int_t ldd)
{
    cublasStatus_t cublasStatus = cublas_func_name(SetMatrix)(nr, nc, sizeof(T_Scalar), src, lds, dest, ldd);
    err::check_cublas(cublasStatus);
}
/*-------------------------------------------------*/
#define matrix_h2d_macro(typein) \
template void MatrixH2D<typein>(int_t, int_t, const typein*, int_t, typein*, int_t);
matrix_h2d_macro(real_t);
matrix_h2d_macro(real4_t);
matrix_h2d_macro(complex_t);
matrix_h2d_macro(complex8_t);
#undef matrix_h2d_macro
/*-------------------------------------------------*/
template <typename T_Scalar>
void MatrixD2H(int_t nr, int_t nc, 
               const T_Scalar *src, int_t lds, 
               T_Scalar *dest, int_t ldd)
{
    cublasStatus_t cublasStatus = cublas_func_name(GetMatrix)(nr, nc, sizeof(T_Scalar), src, lds, dest, ldd);
    err::check_cublas(cublasStatus);
}
/*-------------------------------------------------*/
#define matrix_d2h_macro(typein) \
template void MatrixD2H<typein>(int_t, int_t, const typein*, int_t, typein*, int_t);
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
#define ger_macro(typein, prefix, suffix) \
void ger(cublasHandle_t handle, \
         int_t m, int_t n, \
         const typein *alpha, \
         const typein *x, int_t incx, \
         const typein *y, int_t incy, \
         typein *a, int_t lda) \
{ \
    cublasStatus_t cublasStatus = cublas_func_name(prefix##ger##suffix)(handle, m, n, alpha, x, incx, y, incy, a, lda); \
    err::check_cublas(cublasStatus); \
}
ger_macro(real_t    , D,  );
ger_macro(real4_t   , S,  );
ger_macro(complex_t , Z, u);
ger_macro(complex8_t, C, u);
#undef ger_macro
/*-------------------------------------------------*/
#define gerc_macro(typein, prefix, suffix) \
void gerc(cublasHandle_t handle, \
         int_t m, int_t n, \
         const typein *alpha, \
         const typein *x, int_t incx, \
         const typein *y, int_t incy, \
         typein *a, int_t lda) \
{ \
    cublasStatus_t cublasStatus = cublas_func_name(prefix##ger##suffix)(handle, m, n, alpha, x, incx, y, incy, a, lda); \
    err::check_cublas(cublasStatus); \
}
gerc_macro(real_t    , D,  );
gerc_macro(real4_t   , S,  );
gerc_macro(complex_t , Z, c);
gerc_macro(complex8_t, C, c);
#undef gerc_macro
/*-------------------------------------------------*/
#define syr_macro(typein, prefix) \
void syr(cublasHandle_t handle, \
         cublasFillMode_t uplo, \
         int_t n, \
         const typein *alpha, \
         const typein *x, int_t incx, \
         typein *a, int_t lda) \
{ \
    cublasStatus_t cublasStatus = cublas_func_name(prefix##syr)(handle, uplo, n, alpha, x, incx, a, lda); \
    err::check_cublas(cublasStatus); \
}
syr_macro(real_t    , D);
syr_macro(real4_t   , S);
syr_macro(complex_t , Z);
syr_macro(complex8_t, C);
#undef syr_macro
/*-------------------------------------------------*/
#define her_macro(typein, prefix) \
void her(cublasHandle_t handle, \
         cublasFillMode_t uplo, \
         int_t n, \
         const typename TypeTraits<typein>::real_type *alpha, \
         const typein *x, int_t incx, \
         typein *a, int_t lda) \
{ \
    cublasStatus_t cublasStatus = cublas_func_name(prefix##r)(handle, uplo, n, alpha, x, incx, a, lda); \
    err::check_cublas(cublasStatus); \
}
her_macro(real_t    , Dsy);
her_macro(real4_t   , Ssy);
her_macro(complex_t , Zhe);
her_macro(complex8_t, Che);
#undef her_macro
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
#define gemv_macro(typein, prefix) \
void gemv(cublasHandle_t handle, \
          cublasOperation_t  trans, \
          int_t m, int_t n, const typein* alpha, \
          const typein *a, int_t lda, \
          const typein *x, int_t incx, \
          const typein* beta, typein *y, int_t incy) \
{ \
    cublasStatus_t cublasStatus = cublas_func_name(prefix##gemv)(handle, trans, m, n, alpha, a, lda, x, incx, beta, y, incy); \
    err::check_cublas(cublasStatus); \
}
gemv_macro(real_t    , D);
gemv_macro(real4_t   , S);
gemv_macro(complex_t , Z);
gemv_macro(complex8_t, C);
#undef gemv_macro
/*-------------------------------------------------*/
#define symv_macro(typein, prefix) \
void symv(cublasHandle_t handle, \
          cublasFillMode_t uplo, \
          int_t n, const typein* alpha, \
          const typein *a, int_t lda, \
          const typein *x, int_t incx, \
          const typein* beta, typein *y, int_t incy) \
{ \
    cublasStatus_t cublasStatus = cublas_func_name(prefix##symv)(handle, uplo, n, alpha, a, lda, x, incx, beta, y, incy); \
    err::check_cublas(cublasStatus); \
}
symv_macro(real_t    , D);
symv_macro(real4_t   , S);
symv_macro(complex_t , Z);
symv_macro(complex8_t, C);
#undef symv_macro
/*-------------------------------------------------*/
#define hemv_macro(typein, prefix) \
void hemv(cublasHandle_t handle, \
          cublasFillMode_t uplo, \
          int_t n, const typein* alpha, \
          const typein *a, int_t lda, \
          const typein *x, int_t incx, \
          const typein* beta, typein *y, int_t incy) \
{ \
    cublasStatus_t cublasStatus = cublas_func_name(prefix##mv)(handle, uplo, n, alpha, a, lda, x, incx, beta, y, incy); \
    err::check_cublas(cublasStatus); \
}
hemv_macro(real_t    , Dsy);
hemv_macro(real4_t   , Ssy);
hemv_macro(complex_t , Zhe);
hemv_macro(complex8_t, Che);
#undef hemv_macro
/*-------------------------------------------------*/
#define trmv_macro(typein, prefix) \
void trmv(cublasHandle_t handle, \
          cublasFillMode_t uplo, \
          cublasOperation_t trans, \
          cublasDiagType_t diag, \
          int_t n, const typein* a, int_t lda, \
          typein *x, int_t incx) \
{ \
    cublasStatus_t cublasStatus = cublas_func_name(prefix##trmv)(handle, uplo, trans, diag, n, a, lda, x, incx); \
    err::check_cublas(cublasStatus); \
}
trmv_macro(real_t    , D);
trmv_macro(real4_t   , S);
trmv_macro(complex_t , Z);
trmv_macro(complex8_t, C);
#undef trmv_macro
/*-------------------------------------------------*/
#define trsv_macro(typein, prefix) \
void trsv(cublasHandle_t handle, \
          cublasFillMode_t uplo, \
          cublasOperation_t trans, \
          cublasDiagType_t diag, \
          int_t n, const typein* a, int_t lda, \
          typein *x, int_t incx) \
{ \
    cublasStatus_t cublasStatus = cublas_func_name(prefix##trsv)(handle, uplo, trans, diag, n, a, lda, x, incx); \
    err::check_cublas(cublasStatus); \
}
trsv_macro(real_t    , D);
trsv_macro(real4_t   , S);
trsv_macro(complex_t , Z);
trsv_macro(complex8_t, C);
#undef trsv_macro
/*-------------------------------------------------*/
#define gemm_macro(typein, prefix) \
void gemm(cublasHandle_t handle, \
          cublasOperation_t transa, \
          cublasOperation_t transb, \
          int_t m, int_t n, int_t k, \
          const typein* alpha, \
          const typein* a, int_t lda, \
          const typein* b, int_t ldb, \
          const typein* beta, \
          typein* c, int_t ldc) \
{ \
    cublasStatus_t cublasStatus = cublas_func_name(prefix##gemm)(handle, transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc); \
    err::check_cublas(cublasStatus); \
}
gemm_macro(real_t    , D)
gemm_macro(real4_t   , S)
gemm_macro(complex_t , Z)
gemm_macro(complex8_t, C)
#undef gemm_macro
/*-------------------------------------------------*/
#define symm_macro(typein, prefix) \
void symm(cublasHandle_t handle, \
          cublasSideMode_t side, \
          cublasFillMode_t uplo, \
          int_t m, int_t n, \
          const typein* alpha, \
          const typein* a, int_t lda, \
          const typein* b, int_t ldb, \
          const typein* beta, \
          typein* c, int_t ldc) \
{ \
    cublasStatus_t cublasStatus = cublas_func_name(prefix##symm)(handle, side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc); \
    err::check_cublas(cublasStatus); \
}
symm_macro(real_t    , D)
symm_macro(real4_t   , S)
symm_macro(complex_t , Z)
symm_macro(complex8_t, C)
#undef symm_macro
/*-------------------------------------------------*/
#define hemm_macro(typein, prefix) \
void hemm(cublasHandle_t handle, \
          cublasSideMode_t side, \
          cublasFillMode_t uplo, \
          int_t m, int_t n, \
          const typein* alpha, \
          const typein* a, int_t lda, \
          const typein* b, int_t ldb, \
          const typein* beta, \
          typein* c, int_t ldc) \
{ \
    cublasStatus_t cublasStatus = cublas_func_name(prefix##mm)(handle, side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc); \
    err::check_cublas(cublasStatus); \
}
hemm_macro(real_t    , Dsy)
hemm_macro(real4_t   , Ssy)
hemm_macro(complex_t , Zhe)
hemm_macro(complex8_t, Che)
#undef hemm_macro
/*-------------------------------------------------*/
#define trmm_macro(typein, prefix) \
void trmm(cublasHandle_t handle, \
          cublasSideMode_t side, \
          cublasFillMode_t uplo, \
          cublasOperation_t trans, \
          cublasDiagType_t diag , \
          int_t m, int_t n, \
          const typein* alpha, \
          const typein* a, int_t lda, \
          const typein* b, int_t ldb, \
          typein* c, int_t ldc) \
{ \
    cublasStatus_t cublasStatus = cublas_func_name(prefix##trmm)(handle, side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb, c, ldc); \
    err::check_cublas(cublasStatus); \
}
trmm_macro(real_t    , D)
trmm_macro(real4_t   , S)
trmm_macro(complex_t , Z)
trmm_macro(complex8_t, C)
#undef trmm_macro
/*-------------------------------------------------*/
#define trsm_macro(typein, prefix) \
void trsm(cublasHandle_t handle, \
          cublasSideMode_t side, \
          cublasFillMode_t uplo, \
          cublasOperation_t trans, \
          cublasDiagType_t diag , \
          int_t m, int_t n, \
          const typein* alpha, \
          const typein* a, int_t lda, \
          typein* b, int_t ldb) \
{ \
    cublasStatus_t cublasStatus = cublas_func_name(prefix##trsm)(handle, side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb); \
    err::check_cublas(cublasStatus); \
}
trsm_macro(real_t    , D)
trsm_macro(real4_t   , S)
trsm_macro(complex_t , Z)
trsm_macro(complex8_t, C)
#undef trsm_macro
/*-------------------------------------------------*/
#define syrk_macro(typein, prefix) \
void syrk(cublasHandle_t handle, \
          cublasFillMode_t uplo, \
          cublasOperation_t trans, \
          int_t n, int_t k, \
          const typein* alpha, \
          const typein* a, int_t lda, \
          const typein* beta, \
          typein* c, int_t ldc) \
{ \
    cublasStatus_t cublasStatus = cublas_func_name(prefix##syrk)(handle, uplo, trans, n, k, alpha, a, lda, beta, c, ldc); \
    err::check_cublas(cublasStatus); \
}
syrk_macro(real_t    , D)
syrk_macro(real4_t   , S)
syrk_macro(complex_t , Z)
syrk_macro(complex8_t, C)
#undef syrk_macro
/*-------------------------------------------------*/
#define herk_macro(typein, prefix) \
void herk(cublasHandle_t handle, \
          cublasFillMode_t uplo, \
          cublasOperation_t trans, \
          int_t n, int_t k, \
          const typename TypeTraits<typein>::real_type* alpha, \
          const typein* a, int_t lda, \
          const typename TypeTraits<typein>::real_type* beta, \
          typein* c, int_t ldc) \
{ \
    cublasStatus_t cublasStatus = cublas_func_name(prefix##rk)(handle, uplo, trans, n, k, alpha, a, lda, beta, c, ldc); \
    err::check_cublas(cublasStatus); \
}
herk_macro(real_t    , Dsy)
herk_macro(real4_t   , Ssy)
herk_macro(complex_t , Zhe)
herk_macro(complex8_t, Che)
#undef herk_macro
/*-------------------------------------------------*/
#define syrkx_macro(typein, prefix) \
void syrkx(cublasHandle_t handle, \
          cublasFillMode_t uplo, \
          cublasOperation_t trans, \
          int_t n, int_t k, \
          const typein* alpha, \
          const typein* a, int_t lda, \
          const typein* b, int_t ldb, \
          const typein* beta, \
          typein* c, int_t ldc) \
{ \
    cublasStatus_t cublasStatus = cublas_func_name(prefix##syrkx)(handle, uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc); \
    err::check_cublas(cublasStatus); \
}
syrkx_macro(real_t    , D)
syrkx_macro(real4_t   , S)
syrkx_macro(complex_t , Z)
syrkx_macro(complex8_t, C)
#undef syrkx_macro
/*-------------------------------------------------*/
#define herkx_macro(typein, prefix) \
void herkx(cublasHandle_t handle, \
          cublasFillMode_t uplo, \
          cublasOperation_t trans, \
          int_t n, int_t k, \
          const typein* alpha, \
          const typein* a, int_t lda, \
          const typein* b, int_t ldb, \
          const typename TypeTraits<typein>::real_type* beta, \
          typein* c, int_t ldc) \
{ \
    cublasStatus_t cublasStatus = cublas_func_name(prefix##rkx)(handle, uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc); \
    err::check_cublas(cublasStatus); \
}
herkx_macro(real_t    , Dsy)
herkx_macro(real4_t   , Ssy)
herkx_macro(complex_t , Zhe)
herkx_macro(complex8_t, Che)
#undef herkx_macro
/*-------------------------------------------------*/
} // namespace cublas
} // namespace culite
/*-------------------------------------------------*/
#undef cublas_func_name
/*-------------------------------------------------*/
