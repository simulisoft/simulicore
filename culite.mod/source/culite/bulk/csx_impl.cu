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
#include "culite/bulk/csx_impl.hpp"

// system

// 3rd
#include <cuda_runtime.h>

// culite
#include "culite/types/integer.hpp"
#include "culite/types/scalar.hpp"
#include "culite/support/utils.hpp"

/*-------------------------------------------------*/
namespace culite {
namespace blk {
namespace csx {
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
__global__ void csx_diag_times_vec_kernel(const T_Scalar* alpha, T_Int np, 
                                          const T_Int* xxxptr, 
                                          const T_Int* xxxidx, 
                                          const T_Scalar* values, 
                                          const T_Scalar*x, T_Int incx,
                                          T_Scalar *y, T_Int incy)
{
    T_Int k = blockIdx.x * blockDim.x + threadIdx.x;

    if (k < np) {
        for (T_Int xidx = xxxptr[k]; xidx < xxxptr[k+1]; xidx++) {
            if (xxxidx[xidx] == k) {
                y[k] = y[k] + (*alpha) * values[xidx] * x[k];
                break; // assume only one diagonal element per row
            }
        } // xidx
    } // k
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
void launch_csx_diag_times_vec_kernel(const T_Scalar* alpha, T_Int np, 
                                      const T_Int* xxxptr, 
                                      const T_Int* xxxidx, 
                                      const T_Scalar* values, 
                                      const T_Scalar*x, T_Int incx,
                                      T_Scalar *y, T_Int incy)
{ 
    int threadsPerBlock = 256;
    int blocksPerGrid = (np + threadsPerBlock - 1) / threadsPerBlock;

    csx_diag_times_vec_kernel<<<blocksPerGrid, threadsPerBlock>>>(alpha, np, xxxptr, xxxidx, values, x, incx, y, incy);

    syncDevice();
}
/*-------------------------------------------------*/
template void launch_csx_diag_times_vec_kernel(const real_t    *, int_t, const int_t*, const int_t*, const real_t    *, const real_t    *, int_t, real_t    *, int_t);
template void launch_csx_diag_times_vec_kernel(const real4_t   *, int_t, const int_t*, const int_t*, const real4_t   *, const real4_t   *, int_t, real4_t   *, int_t);
template void launch_csx_diag_times_vec_kernel(const complex_t *, int_t, const int_t*, const int_t*, const complex_t *, const complex_t *, int_t, complex_t *, int_t);
template void launch_csx_diag_times_vec_kernel(const complex8_t*, int_t, const int_t*, const int_t*, const complex8_t*, const complex8_t*, int_t, complex8_t*, int_t);
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
__global__ void csx_diag_times_mat_kernel(const T_Scalar* alpha, T_Int np, 
                                          const T_Int* __restrict__ xxxptr, 
                                          const T_Int* __restrict__ xxxidx, 
                                          const T_Scalar* __restrict__ values, 
                                          T_Int nc,
                                          const T_Scalar*__restrict__ b, T_Int ldb,
                                          T_Scalar* __restrict__ c, T_Int ldc)
{
    T_Int k = blockIdx.x * blockDim.x + threadIdx.x;
    T_Int l = blockIdx.y * blockDim.y + threadIdx.y;

    if (k < np && l < nc) {
        for (T_Int xidx = xxxptr[k]; xidx < xxxptr[k+1]; xidx++) {
            if (xxxidx[xidx] == k) {
                c[ldc * l + k] = c[ldc * l + k] + (*alpha) * values[xidx] * b[ldb * l + k];
                break;
            } // diag
        } // xidx
    } // dim check
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
void launch_csx_diag_times_mat_kernel(const T_Scalar* alpha, T_Int np, 
                                      const T_Int* xxxptr, 
                                      const T_Int* xxxidx, 
                                      const T_Scalar* values, 
                                      T_Int nc,
                                      const T_Scalar*b, T_Int ldb,
                                      T_Scalar *c, T_Int ldc)
{ 
    // TODO: check if nc is large enough to warrant a 2D grid, otherwise use a 1D grid
    dim3 blockSize(16, 16);
    dim3 gridSize((np + blockSize.x - 1) / blockSize.x, 
                  (nc  + blockSize.y - 1) / blockSize.y);

    csx_diag_times_mat_kernel<<<gridSize, blockSize>>>(alpha, np, xxxptr, xxxidx, values, nc, b, ldb, c, ldc);

    syncDevice();
}
/*-------------------------------------------------*/
template void launch_csx_diag_times_mat_kernel(const real_t    *, int_t, const int_t*, const int_t*, const real_t    *, int_t, const real_t    *, int_t, real_t    *, int_t);
template void launch_csx_diag_times_mat_kernel(const real4_t   *, int_t, const int_t*, const int_t*, const real4_t   *, int_t, const real4_t   *, int_t, real4_t   *, int_t);
template void launch_csx_diag_times_mat_kernel(const complex_t *, int_t, const int_t*, const int_t*, const complex_t *, int_t, const complex_t *, int_t, complex_t *, int_t);
template void launch_csx_diag_times_mat_kernel(const complex8_t*, int_t, const int_t*, const int_t*, const complex8_t*, int_t, const complex8_t*, int_t, complex8_t*, int_t);
/*-------------------------------------------------*/

} // namespace csx
} // namespace blk
} // namespace culite
/*-------------------------------------------------*/
