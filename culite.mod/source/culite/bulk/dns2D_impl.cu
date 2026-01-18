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
#include "culite/bulk/dns2D_impl.hpp"

// system

// 3rd

// culite
#include "culite/types/scalar.hpp"
#include "culite/support/utils.hpp"

/*-------------------------------------------------*/
namespace culite {
namespace blk {
namespace dns {
/*-------------------------------------------------*/
template <typename T_Scalar>
__global__ void scale_matrix_kernel(int_t m, int_t n, T_Scalar* v, int_t ldv, T_Scalar alpha)
{
    int_t i = blockIdx.x * blockDim.x + threadIdx.x;
    int_t j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i < m && j < n) {
        v[j * ldv + i] = alpha * v[j * ldv + i];
    }
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void launch_scale_matrix_kernel(int_t m, int_t n, T_Scalar* v, int_t ldv, T_Scalar alpha)
{
    if(m <= 0 || n <= 0) return;

    dim3 threads(16, 16);
    dim3 blocks((m + threads.x - 1) / threads.x, (n + threads.y - 1) / threads.y);
    scale_matrix_kernel<T_Scalar><<<blocks, threads>>>(m, n, v, ldv, alpha);
    syncDevice();
}
/*-------------------------------------------------*/
template void launch_scale_matrix_kernel<real_t>(int_t m, int_t n, real_t* v, int_t ldv, real_t alpha);
template void launch_scale_matrix_kernel<real4_t>(int_t m, int_t n, real4_t* v, int_t ldv, real4_t alpha);
template void launch_scale_matrix_kernel<complex_t>(int_t m, int_t n, complex_t* v, int_t ldv, complex_t alpha);
template void launch_scale_matrix_kernel<complex8_t>(int_t m, int_t n, complex8_t* v, int_t ldv, complex8_t alpha);
/*-------------------------------------------------*/





/*-------------------------------------------------*/
} // namespace dns
} // namespace blk
} // namespace culite
/*-------------------------------------------------*/
