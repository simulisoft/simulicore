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
#include "culite/bulk/dns1D_impl.hpp"

// system

// 3rd
#include <cuda_runtime.h>

// culite
#include "culite/types/scalar.hpp"
#include "culite/support/utils.hpp"

/*-------------------------------------------------*/
namespace culite {
namespace blk {
namespace dns {
/*-------------------------------------------------*/
template <typename T_Scalar>
__global__ void fill_1d_kernel(int_t n, T_Scalar* x, T_Scalar val)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) x[idx] = val;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void launch_fill_kernel(int_t n, T_Scalar* x, T_Scalar val)
{ 
    if(n <= 0) return;

	int threads = 256;
	int blocks = (n + threads - 1) / threads;
	fill_1d_kernel<T_Scalar><<<blocks, threads>>>(n, x, val);

    syncDevice();
}
/*-------------------------------------------------*/
template void launch_fill_kernel(int_t, real4_t*, real4_t);
template void launch_fill_kernel(int_t, real_t*, real_t);
template void launch_fill_kernel(int_t, complex_t*, complex_t);
template void launch_fill_kernel(int_t, complex8_t*, complex8_t);
/*-------------------------------------------------*/
template <typename T_Scalar>
__global__ void conjugate_1d_kernel(int_t n, T_Scalar* x)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) x[idx] = arith::conj(x[idx]);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void launch_conjugate_kernel(int_t n, T_Scalar* x)
{ 
    if(n <= 0 || TypeTraits<T_Scalar>::is_real()) return;

	int threads = 256;
	int blocks = (n + threads - 1) / threads;
	conjugate_1d_kernel<T_Scalar><<<blocks, threads>>>(n, x);

    syncDevice();
}
/*-------------------------------------------------*/
template void launch_conjugate_kernel(int_t, real4_t*);
template void launch_conjugate_kernel(int_t, real_t*);
template void launch_conjugate_kernel(int_t, complex_t*);
template void launch_conjugate_kernel(int_t, complex8_t*);
/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/
template <typename T_Scalar>
__global__ void get_real_1d_kernel(int_t n, const T_Scalar* x, 
                                   typename TypeTraits<T_Scalar>::real_type *y)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) y[idx] = arith::getRe(x[idx]);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void launch_get_real_kernel(int_t n, const T_Scalar* x, 
                            typename TypeTraits<T_Scalar>::real_type* y)
{
    if(n <= 0) return;

    int threads = 256;
    int blocks = (n + threads - 1) / threads;
    get_real_1d_kernel<T_Scalar><<<blocks, threads>>>(n, x, y);
    syncDevice();
}
/*-------------------------------------------------*/
template void launch_get_real_kernel<complex_t>(int_t, const complex_t*, real_t*);
template void launch_get_real_kernel<complex8_t>(int_t, const complex8_t*, real4_t*);
/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/
template <typename T_Scalar>
__global__ void get_imag_1d_kernel(int_t n, const T_Scalar* x, 
                                   typename TypeTraits<T_Scalar>::real_type *y)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) y[idx] = arith::getIm(x[idx]);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void launch_get_imag_kernel(int_t n, const T_Scalar* x, 
                            typename TypeTraits<T_Scalar>::real_type* y)
{
    if(n <= 0) return;

    int threads = 256;
    int blocks = (n + threads - 1) / threads;
    get_imag_1d_kernel<T_Scalar><<<blocks, threads>>>(n, x, y);
    syncDevice();
}
/*-------------------------------------------------*/
template void launch_get_imag_kernel<complex_t>(int_t, const complex_t*, real_t*);
template void launch_get_imag_kernel<complex8_t>(int_t, const complex8_t*, real4_t*);
/*-------------------------------------------------*/
template <typename T_Scalar>
__global__ void geev_calculate_complex_eigenvalues_kernel(int_t n, 
                                                          const typename TypeTraits<T_Scalar>::real_type* wri, 
                                                          T_Scalar *w)
{
    using T_RScalar = typename TypeTraits<T_Scalar>::real_type;

    const T_RScalar *wr = wri;
    const T_RScalar *wi = wr + n;
 
    int_t i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < n) {
        w[i] = makeComplex(wr[i], wi[i]);
    } // dim check
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void launch_geev_calculate_complex_eigenvalues_kernel(int_t n, 
                                                      const typename TypeTraits<T_Scalar>::real_type* wri, 
                                                      T_Scalar *w)
{
    if(n <= 0) return;
    
    int threads = 256;
    int blocks = (n + threads - 1) / threads;
    geev_calculate_complex_eigenvalues_kernel<T_Scalar><<<blocks, threads>>>(n, wri, w);

    syncDevice();
}
/*-------------------------------------------------*/
template void launch_geev_calculate_complex_eigenvalues_kernel<complex_t>(int_t, const real_t*, complex_t*);
template void launch_geev_calculate_complex_eigenvalues_kernel<complex8_t>(int_t, const real4_t*, complex8_t*);
/*-------------------------------------------------*/
} // namespace dns
} // namespace blk
} // namespace culite
/*-------------------------------------------------*/
