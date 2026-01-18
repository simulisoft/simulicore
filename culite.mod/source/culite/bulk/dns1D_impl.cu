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
__global__ void conjugate_kernel_1d(std::size_t n, real_t* d) { }
/*-------------------------------------------------*/
__global__ void conjugate_kernel_1d(std::size_t n, real4_t* s) { }
/*-------------------------------------------------*/
__global__ void conjugate_kernel_1d(std::size_t n, complex_t* z) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) z[idx] = cuConj(z[idx]);
}
/*-------------------------------------------------*/
__global__ void conjugate_kernel_1d(std::size_t n, complex8_t* c) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) c[idx] = cuConjf(c[idx]);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void launch_conjugate_kernel_1d(std::size_t n, T_Scalar* z)
{ 
    if(n <= 0) return;
    if(TypeTraits<T_Scalar>::is_real()) return;

	int threads = 256;
	int blocks = (n + threads - 1) / threads;
	conjugate_kernel_1d<<<blocks, threads>>>(n, z);

    syncDevice();
}
/*-------------------------------------------------*/
template void launch_conjugate_kernel_1d(std::size_t, real4_t*);
template void launch_conjugate_kernel_1d(std::size_t, real_t*);
template void launch_conjugate_kernel_1d(std::size_t, complex_t*);
template void launch_conjugate_kernel_1d(std::size_t, complex8_t*);
/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/
__global__ void get_real_kernel_1d(std::size_t n, const complex_t* z, real_t *d) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) d[idx] = cuCreal(z[idx]);
}
/*-------------------------------------------------*/
__global__ void get_real_kernel_1d(std::size_t n, const complex8_t* c, real4_t *s) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) s[idx] = cuCrealf(c[idx]);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void launch_get_real_kernel_1d(std::size_t n, const T_Scalar* x, typename TypeTraits<T_Scalar>::real_type* y)
{
    int threads = 256;
    int blocks = (n + threads - 1) / threads;
    get_real_kernel_1d<<<blocks, threads>>>(n, x, y);
    syncDevice();
}
/*-------------------------------------------------*/
template void launch_get_real_kernel_1d<complex_t>(std::size_t, const complex_t*, real_t*);
template void launch_get_real_kernel_1d<complex8_t>(std::size_t, const complex8_t*, real4_t*);
/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/
__global__ void get_imag_kernel_1d(std::size_t n, const complex_t* z, real_t *d) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) d[idx] = cuCimag(z[idx]);
}
/*-------------------------------------------------*/
__global__ void get_imag_kernel_1d(std::size_t n, const complex8_t* c, real4_t *s) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) s[idx] = cuCimagf(c[idx]);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void launch_get_imag_kernel_1d(std::size_t n, const T_Scalar* x, typename TypeTraits<T_Scalar>::real_type* y)
{
    int threads = 256;
    int blocks = (n + threads - 1) / threads;
    get_imag_kernel_1d<<<blocks, threads>>>(n, x, y);
    syncDevice();
}
/*-------------------------------------------------*/
template void launch_get_imag_kernel_1d<complex_t>(std::size_t, const complex_t*, real_t*);
template void launch_get_imag_kernel_1d<complex8_t>(std::size_t, const complex8_t*, real4_t*);
/*-------------------------------------------------*/
} // namespace dns
} // namespace blk
} // namespace culite
/*-------------------------------------------------*/
