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
#include "culite/bulk/dns2d_impl.hpp"

// system
#include <cmath>

// 3rd
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

// culite
#include "culite/types/scalar.hpp"
#include "culite/support/utils.hpp"
#include "culite/support/grid.hpp"
#include "culite/support/imalloc.hpp"

/*-------------------------------------------------*/
namespace culite {
namespace blk {
namespace dns {
/*-------------------------------------------------*/
template <typename T_Scalar>
__global__ void copy_2d_kernel(::cla3p::uplo_t uplo, int_t m, int_t n, 
                               const T_Scalar* a, int_t lda, 
                               T_Scalar* b, int_t ldb)
{
    int_t i = blockIdx.x * blockDim.x + threadIdx.x;
    int_t j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i < m && j < n) {
        if(uplo == ::cla3p::uplo_t::Upper && i > j) return;
        if(uplo == ::cla3p::uplo_t::Lower && i < j) return;
        b[j * ldb + i] = a[j * lda + i];
    }
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void launch_copy_kernel(::cla3p::uplo_t uplo, int_t m, int_t n, 
                        const T_Scalar *a, int_t lda, 
                        T_Scalar *b, int_t ldb)
{
    if(m <= 0 || n <= 0) return;

    Grid2D grid(m, n);
    copy_2d_kernel<T_Scalar><<<grid.numBlocks(), grid.threadsPerBlock()>>>(uplo, m, n, a, lda, b, ldb);

    syncDevice();
}
/*-------------------------------------------------*/
template void launch_copy_kernel(::cla3p::uplo_t, int_t, int_t, const real_t*, int_t, real_t*, int_t);
template void launch_copy_kernel(::cla3p::uplo_t, int_t, int_t, const real4_t*, int_t, real4_t*, int_t);
template void launch_copy_kernel(::cla3p::uplo_t, int_t, int_t, const complex_t*, int_t, complex_t*, int_t);
template void launch_copy_kernel(::cla3p::uplo_t, int_t, int_t, const complex8_t*, int_t, complex8_t*, int_t);
/*-------------------------------------------------*/
template <typename T_Scalar>
__global__ void fill_2d_kernel(::cla3p::uplo_t uplo, int_t m, int_t n, T_Scalar *a, int_t lda, T_Scalar val)
{
    int_t i = blockIdx.x * blockDim.x + threadIdx.x;
    int_t j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i < m && j < n) {
        if(uplo == ::cla3p::uplo_t::Upper && i > j) return;
        if(uplo == ::cla3p::uplo_t::Lower && i < j) return;
        a[j * lda + i] = val;
    }
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void launch_fill_kernel(::cla3p::uplo_t uplo, int_t m, int_t n, T_Scalar *a, int_t lda, T_Scalar val)
{
    if(m <= 0 || n <= 0) return;

    // TODO: consider using memSetZero() for zero values --- IGNORE ---

    Grid2D grid(m, n);
    fill_2d_kernel<T_Scalar><<<grid.numBlocks(), grid.threadsPerBlock()>>>(uplo, m, n, a, lda, val);
    syncDevice();
}
/*-------------------------------------------------*/
template void launch_fill_kernel(::cla3p::uplo_t, int_t, int_t, real_t*, int_t, real_t);
template void launch_fill_kernel(::cla3p::uplo_t, int_t, int_t, real4_t*, int_t, real4_t);
template void launch_fill_kernel(::cla3p::uplo_t, int_t, int_t, complex_t*, int_t, complex_t);
template void launch_fill_kernel(::cla3p::uplo_t, int_t, int_t, complex8_t*, int_t, complex8_t);
/*-------------------------------------------------*/
template <typename T_Scalar>
__global__ void scale_2d_kernel(::cla3p::uplo_t uplo, int_t m, int_t n, T_Scalar* a, int_t lda, T_Scalar alpha)
{
    int_t i = blockIdx.x * blockDim.x + threadIdx.x;
    int_t j = blockIdx.y * blockDim.y + threadIdx.y;
    
    if (i < m && j < n) {
        if(uplo == ::cla3p::uplo_t::Upper && i > j) return;
        if(uplo == ::cla3p::uplo_t::Lower && i < j) return;
        a[j * lda + i] = alpha * a[j * lda + i];
    }
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void launch_scale_kernel(::cla3p::uplo_t uplo, int_t m, int_t n, T_Scalar* a, int_t lda, T_Scalar alpha)
{
    if(m <= 0 || n <= 0) return;

    Grid2D grid(m, n);
    scale_2d_kernel<T_Scalar><<<grid.numBlocks(), grid.threadsPerBlock()>>>(uplo, m, n, a, lda, alpha);
    syncDevice();
}
/*-------------------------------------------------*/
template void launch_scale_kernel<real_t>(::cla3p::uplo_t, int_t, int_t, real_t*, int_t, real_t);
template void launch_scale_kernel<real4_t>(::cla3p::uplo_t, int_t, int_t, real4_t*, int_t, real4_t);
template void launch_scale_kernel<complex_t>(::cla3p::uplo_t, int_t, int_t, complex_t*, int_t, complex_t);
template void launch_scale_kernel<complex8_t>(::cla3p::uplo_t, int_t, int_t, complex8_t*, int_t, complex8_t);
/*-------------------------------------------------*/
__device__ static double atomicMax(double* address, double val) {
    unsigned long long int* address_as_ull = (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(fmax(val, __longlong_as_double(assumed))));
    } while (assumed != old);
    return __longlong_as_double(old);
}
/*-------------------------------------------------*/
__device__ static float atomicMax(float* address, float val) {
    int* address_as_int = (int*)address;
    int old = *address_as_int, assumed;
    do {
        assumed = old;
        old = atomicCAS(address_as_int, assumed,
                        __float_as_int(fmaxf(val, __int_as_float(assumed))));
    } while (assumed != old);
    return __int_as_float(old);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
__global__ void matrix_one_norm_kernel(int_t m, int_t n, const T_Scalar* a, int_t lda, typename TypeTraits<T_Scalar>::real_type* ret)
{
    using T_RScalar = typename TypeTraits<T_Scalar>::real_type;

    int_t j = blockIdx.x * blockDim.x + threadIdx.x;

    if (j < n) {
        T_RScalar colSum = 0;
        for (int_t i = 0; i < m; ++i) {
            colSum += arith::abs(a[j * lda + i]);
        } // i
        atomicMax(ret, colSum);
    } // j
}
/*-------------------------------------------------*/
template <typename T_Scalar>
typename TypeTraits<T_Scalar>::real_type
launch_matrix_one_norm_kernel(int_t m, int_t n, const T_Scalar* a, int_t lda)
{
    using T_RScalar = typename TypeTraits<T_Scalar>::real_type;
    
    T_RScalar retHost = 0;
    T_RScalar* retDevice = device_alloc_t<T_RScalar>(1);
    cudaMemset(retDevice, 0, sizeof(T_RScalar));

    Grid1D grid(n);
    matrix_one_norm_kernel<T_Scalar><<<grid.numBlocks(), grid.threadsPerBlock()>>>(m, n, a, lda, retDevice);

    memCopyD2H(1, retDevice, &retHost);
    device_free(retDevice);
    syncDevice(); // double sync affects performance?

    return retHost;
}
/*-------------------------------------------------*/
template real_t launch_matrix_one_norm_kernel<real_t>(int_t, int_t, const real_t*, int_t);
template real4_t launch_matrix_one_norm_kernel<real4_t>(int_t, int_t, const real4_t*, int_t);
template real_t launch_matrix_one_norm_kernel<complex_t>(int_t, int_t, const complex_t*, int_t);
template real4_t launch_matrix_one_norm_kernel<complex8_t>(int_t, int_t, const complex8_t*, int_t);
/*-------------------------------------------------*/
template <typename T_Scalar>
__global__ void fill_row_sums(int_t m, int_t n, const T_Scalar* a, int_t lda, typename TypeTraits<T_Scalar>::real_type* rowSums) 
{
    int_t i = blockIdx.x * blockDim.x + threadIdx.x;
    int_t j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i < m && j < n) {
        atomicAdd(&rowSums[i], arith::abs(a[j * lda + i]));
    }
}
/*-------------------------------------------------*/
template <typename T_RScalar>
__global__ void find_max_sum(int_t m, T_RScalar* rowSums, T_RScalar* globalMax) 
{
    int_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < m) {
        atomicMax(globalMax, rowSums[i]);
    }
}
/*-------------------------------------------------*/
template <typename T_Scalar>
typename TypeTraits<T_Scalar>::real_type
launch_matrix_inf_norm_kernel(int_t m, int_t n, const T_Scalar* a, int_t lda) 
{
    using T_RScalar = typename TypeTraits<T_Scalar>::real_type;

    if (m <= 0 || n <= 0) return 0;
    
    T_RScalar retHost = 0;
    T_RScalar *rowSums = device_alloc_t<T_RScalar>(m);
    T_RScalar *retDevice = device_alloc_t<T_RScalar>(1);

    cudaMemset(rowSums, 0, m * sizeof(T_RScalar));
    cudaMemset(retDevice, 0, sizeof(T_RScalar));

    Grid2D grid2d(m, n);
    fill_row_sums<T_Scalar><<<grid2d.numBlocks(), grid2d.threadsPerBlock()>>>(m, n, a, lda, rowSums);

    Grid1D grid1d(m);
    find_max_sum<T_RScalar><<<grid1d.numBlocks(), grid1d.threadsPerBlock()>>>(m, rowSums, retDevice);

    memCopyD2H(1, retDevice, &retHost);
    device_free(rowSums);
    device_free(retDevice);
    syncDevice(); // double sync affects performance?

    return retHost;
}
/*-------------------------------------------------*/
template real_t launch_matrix_inf_norm_kernel<real_t>(int_t, int_t, const real_t*, int_t);
template real4_t  launch_matrix_inf_norm_kernel<real4_t>(int_t, int_t, const real4_t*, int_t);
template real_t launch_matrix_inf_norm_kernel<complex_t>(int_t, int_t, const complex_t*, int_t);
template real4_t  launch_matrix_inf_norm_kernel<complex8_t>(int_t, int_t, const complex8_t*, int_t);
/*-------------------------------------------------*/
template <typename T_Scalar>
__global__ void matrix_max_norm_kernel(int_t m, int_t n, const T_Scalar* a, int_t lda, typename TypeTraits<T_Scalar>::real_type* globalMax) 
{
    using T_RScalar = typename TypeTraits<T_Scalar>::real_type;

    int_t j = blockIdx.x * blockDim.x + threadIdx.x;

    if (j < n) {
        T_RScalar localMax = 0;
        for (int_t i = 0; i < m; ++i) {
            localMax = fmax(localMax, arith::abs(a[j * lda + i]));
        }
        atomicMax(globalMax, localMax);
    }
}
/*-------------------------------------------------*/
template <typename T_Scalar>
typename TypeTraits<T_Scalar>::real_type
launch_matrix_max_norm_kernel(int_t m, int_t n, const T_Scalar* a, int_t lda) 
{
    using T_RScalar = typename TypeTraits<T_Scalar>::real_type;

    if (m <= 0 || n <= 0) return 0;

    T_RScalar retHost = 0;
    T_RScalar* retDevice = device_alloc_t<T_RScalar>(1);
    cudaMemset(retDevice, 0, sizeof(T_RScalar));

    Grid1D grid(n);
    matrix_max_norm_kernel<T_Scalar><<<grid.numBlocks(), grid.threadsPerBlock()>>>(m, n, a, lda, retDevice);

    memCopyD2H(1, retDevice, &retHost);
    device_free(retDevice);
    syncDevice(); // double sync affects performance?

    return retHost;
}
/*-------------------------------------------------*/
template real_t launch_matrix_max_norm_kernel<real_t>(int_t, int_t, const real_t*, int_t);
template real4_t  launch_matrix_max_norm_kernel<real4_t>(int_t, int_t, const real4_t*, int_t);
template real_t launch_matrix_max_norm_kernel<complex_t>(int_t, int_t, const complex_t*, int_t);
template real4_t  launch_matrix_max_norm_kernel<complex8_t>(int_t, int_t, const complex8_t*, int_t);
/*-------------------------------------------------*/
template <typename T_Scalar>
__global__ void matrix_fro_norm_kernel(int_t m, int_t n, const T_Scalar* a, int_t lda, typename TypeTraits<T_Scalar>::real_type* globalSum)
{
    using T_RScalar = typename TypeTraits<T_Scalar>::real_type;

    int_t j = blockIdx.x * blockDim.x + threadIdx.x;

    if (j < n) {
        T_RScalar colSum = 0;
        for (int_t i = 0; i < m; ++i) {
            colSum += arith::abs2(a[j * lda + i]);
        }
        atomicAdd(globalSum, colSum);
    }
}
/*-------------------------------------------------*/
template <typename T_Scalar>
typename TypeTraits<T_Scalar>::real_type
launch_matrix_fro_norm_kernel(int_t m, int_t n, const T_Scalar* a, int_t lda) 
{
    using T_RScalar = typename TypeTraits<T_Scalar>::real_type;

    if (m <= 0 || n <= 0) return 0;

    T_RScalar retHost = 0;
    T_RScalar* retDevice = device_alloc_t<T_RScalar>(1);
    cudaMemset(retDevice, 0, sizeof(T_RScalar));

    Grid1D grid(n);
    matrix_fro_norm_kernel<T_Scalar><<<grid.numBlocks(), grid.threadsPerBlock()>>>(m, n, a, lda, retDevice);

    memCopyD2H(1, retDevice, &retHost);
    device_free(retDevice);
    syncDevice(); // double sync affects performance?
    
    return std::sqrt(retHost);
}
/*-------------------------------------------------*/
template real_t launch_matrix_fro_norm_kernel<real_t>(int_t, int_t, const real_t*, int_t);
template real4_t  launch_matrix_fro_norm_kernel<real4_t>(int_t, int_t, const real4_t*, int_t);
template real_t launch_matrix_fro_norm_kernel<complex_t>(int_t, int_t, const complex_t*, int_t);
template real4_t  launch_matrix_fro_norm_kernel<complex8_t>(int_t, int_t, const complex8_t*, int_t);
/*-------------------------------------------------*/
template <typename T_Scalar>
__global__ void get_real_2d_kernel(::cla3p::uplo_t uplo, int_t m, int_t n, 
                                   const T_Scalar* a, int_t lda, 
                                   typename TypeTraits<T_Scalar>::real_type* b, int_t ldb)
{
    int_t i = blockIdx.x * blockDim.x + threadIdx.x;
    int_t j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i < m && j < n) {
        if(uplo == ::cla3p::uplo_t::Upper && i > j) return;
        if(uplo == ::cla3p::uplo_t::Lower && i < j) return;
        b[j * ldb + i] = arith::getRe(a[j * lda + i]);
    }
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void launch_get_real_kernel(::cla3p::uplo_t uplo, int_t m, int_t n, 
                            const T_Scalar* a, int_t lda, 
                            typename TypeTraits<T_Scalar>::real_type* b, int_t ldb)
{
    if(m <= 0 || n <= 0) return;

    Grid2D grid(m, n);
    get_real_2d_kernel<<<grid.numBlocks(), grid.threadsPerBlock()>>>(uplo, m, n, a, lda, b, ldb);

    syncDevice();
}
/*-------------------------------------------------*/
template void launch_get_real_kernel<complex_t>(::cla3p::uplo_t, int_t, int_t, const complex_t*, int_t, real_t*, int_t);
template void launch_get_real_kernel<complex8_t>(::cla3p::uplo_t, int_t, int_t, const complex8_t*, int_t, real4_t*, int_t);
/*-------------------------------------------------*/
template <typename T_Scalar>
__global__ void get_imag_2d_kernel(::cla3p::uplo_t uplo, int_t m, int_t n, 
                                   const T_Scalar* a, int_t lda, 
                                   typename TypeTraits<T_Scalar>::real_type* b, int_t ldb)
{
    int_t i = blockIdx.x * blockDim.x + threadIdx.x;
    int_t j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i < m && j < n) {
        if(uplo == ::cla3p::uplo_t::Upper && i > j) return;
        if(uplo == ::cla3p::uplo_t::Lower && i < j) return;
        b[j * ldb + i] = arith::getIm(a[j * lda + i]);
    }
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void launch_get_imag_kernel(::cla3p::uplo_t uplo, int_t m, int_t n, 
                            const T_Scalar* a, int_t lda, 
                            typename TypeTraits<T_Scalar>::real_type* b, int_t ldb)
{
    if(m <= 0 || n <= 0) return;

    Grid2D grid(m, n);
    get_imag_2d_kernel<<<grid.numBlocks(), grid.threadsPerBlock()>>>(uplo, m, n, a, lda, b, ldb);

    syncDevice();
}
/*-------------------------------------------------*/
template void launch_get_imag_kernel<complex_t>(::cla3p::uplo_t, int_t, int_t, const complex_t*, int_t, real_t*, int_t);
template void launch_get_imag_kernel<complex8_t>(::cla3p::uplo_t, int_t, int_t, const complex8_t*, int_t, real4_t*, int_t);
/*-------------------------------------------------*/
template <typename T_Scalar>
__global__ void conjugate_2d_kernel(::cla3p::uplo_t uplo, int_t m, int_t n, T_Scalar* a, int_t lda)
{
    int_t i = blockIdx.x * blockDim.x + threadIdx.x;
    int_t j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i < m && j < n) {
        if(uplo == ::cla3p::uplo_t::Upper && i > j) return;
        if(uplo == ::cla3p::uplo_t::Lower && i < j) return;
        a[j * lda + i] = arith::conj(a[j * lda + i]);
    }
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void launch_conjugate_kernel(::cla3p::uplo_t uplo, int_t m, int_t n, T_Scalar* a, int_t lda)
{
    if(m <= 0 || n <= 0) return;

    Grid2D grid(m, n);
    conjugate_2d_kernel<<<grid.numBlocks(), grid.threadsPerBlock()>>>(uplo, m, n, a, lda);

    syncDevice();
}
/*-------------------------------------------------*/
template void launch_conjugate_kernel<real_t>(::cla3p::uplo_t, int_t, int_t, real_t*, int_t);
template void launch_conjugate_kernel<real4_t>(::cla3p::uplo_t, int_t, int_t, real4_t*, int_t);
template void launch_conjugate_kernel<complex_t>(::cla3p::uplo_t, int_t, int_t, complex_t*, int_t);
template void launch_conjugate_kernel<complex8_t>(::cla3p::uplo_t, int_t, int_t, complex8_t*, int_t);
/*-------------------------------------------------*/
template <typename T_Scalar>
__global__ void geev_calculate_complex_eigenvectors_kernel(int_t n, 
                                                           const typename TypeTraits<T_Scalar>::real_type* w, 
                                                           const typename TypeTraits<T_Scalar>::real_type* vr, int_t ldvr,
                                                           T_Scalar *vc, int_t ldvc)
{
    using T_RScalar = typename TypeTraits<T_Scalar>::real_type;

    const T_RScalar *wr = w;
    const T_RScalar *wi = wr + n;
 
    int_t i = blockIdx.x * blockDim.x + threadIdx.x;
    int_t j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i < n && j < n) {

        T_RScalar im = wi[j];

        if (im == 0) {

            vc[j * ldvc + i] = makeScalar<T_Scalar>(vr[j * ldvr + i]);

        } else if (im > 0) {

            T_RScalar u = vr[ j      * ldvr + i];
            T_RScalar v = vr[(j + 1) * ldvr + i];
            vc[j * ldvc + i] = makeComplex(u, v);

        } else {

            T_RScalar u = vr[(j - 1) * ldvr + i];
            T_RScalar v = vr[ j      * ldvr + i];
            vc[j * ldvc + i] = makeComplex(u, -v);

        } // im cases

    } // dim check
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void launch_geev_calculate_complex_eigenvectors_kernel(int_t n, 
                                                       const typename TypeTraits<T_Scalar>::real_type* w, 
                                                       const typename TypeTraits<T_Scalar>::real_type* vr, int_t ldvr,
                                                       T_Scalar *vc, int_t ldvc)
{
    if(n <= 0) return;

    Grid2D grid(n, n);
    geev_calculate_complex_eigenvectors_kernel<T_Scalar><<<grid.numBlocks(), grid.threadsPerBlock()>>>(n, w, vr, ldvr, vc, ldvc);

    syncDevice();
}
/*-------------------------------------------------*/
template void launch_geev_calculate_complex_eigenvectors_kernel<complex_t>(int_t, const real_t*, const real_t*, int_t, complex_t*, int_t);
template void launch_geev_calculate_complex_eigenvectors_kernel<complex8_t>(int_t, const real4_t*, const real4_t*, int_t, complex8_t*, int_t);
/*-------------------------------------------------*/
} // namespace dns
} // namespace blk
} // namespace culite
/*-------------------------------------------------*/
