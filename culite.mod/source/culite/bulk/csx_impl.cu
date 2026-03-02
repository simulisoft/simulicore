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
#include "culite/support/grid.hpp"

/*-------------------------------------------------*/
namespace culite {
namespace blk {
namespace csx {
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
__global__ void diag_times_vec_kernel(const T_Scalar* alpha, T_Int np, 
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
void launch_diag_times_vec_kernel(const T_Scalar* alpha, T_Int np, 
                                  const T_Int* xxxptr, 
                                  const T_Int* xxxidx, 
                                  const T_Scalar* values, 
                                  const T_Scalar*x, T_Int incx,
                                  T_Scalar *y, T_Int incy)
{
    if(np <= 0) return;

    Grid1D grid(np);
    diag_times_vec_kernel<<<grid.numBlocks(), 
                            grid.threadsPerBlock()>>>(alpha, np, xxxptr, xxxidx, values, x, incx, y, incy);

    syncDevice();
}
/*-------------------------------------------------*/
template void launch_diag_times_vec_kernel(const real_t    *, int_t, const int_t*, const int_t*, const real_t    *, const real_t    *, int_t, real_t    *, int_t);
template void launch_diag_times_vec_kernel(const real4_t   *, int_t, const int_t*, const int_t*, const real4_t   *, const real4_t   *, int_t, real4_t   *, int_t);
template void launch_diag_times_vec_kernel(const complex_t *, int_t, const int_t*, const int_t*, const complex_t *, const complex_t *, int_t, complex_t *, int_t);
template void launch_diag_times_vec_kernel(const complex8_t*, int_t, const int_t*, const int_t*, const complex8_t*, const complex8_t*, int_t, complex8_t*, int_t);
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
__global__ void diag_times_mat_kernel(const T_Scalar* alpha, T_Int np, 
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
void launch_diag_times_mat_kernel(const T_Scalar* alpha, T_Int np, 
                                  const T_Int* xxxptr, 
                                  const T_Int* xxxidx, 
                                  const T_Scalar* values, 
                                  T_Int nc,
                                  const T_Scalar*b, T_Int ldb,
                                  T_Scalar *c, T_Int ldc)
{ 
    if(np <= 0 || nc <= 0) return;

    // TODO: check if nc is large enough to warrant a 2D grid, otherwise use a 1D grid
    Grid2D grid(np, nc);
    diag_times_mat_kernel<<<grid.numBlocks(), 
                            grid.threadsPerBlock()>>>(alpha, np, xxxptr, xxxidx, values, nc, b, ldb, c, ldc);

    syncDevice();
}
/*-------------------------------------------------*/
template void launch_diag_times_mat_kernel(const real_t    *, int_t, const int_t*, const int_t*, const real_t    *, int_t, const real_t    *, int_t, real_t    *, int_t);
template void launch_diag_times_mat_kernel(const real4_t   *, int_t, const int_t*, const int_t*, const real4_t   *, int_t, const real4_t   *, int_t, real4_t   *, int_t);
template void launch_diag_times_mat_kernel(const complex_t *, int_t, const int_t*, const int_t*, const complex_t *, int_t, const complex_t *, int_t, complex_t *, int_t);
template void launch_diag_times_mat_kernel(const complex8_t*, int_t, const int_t*, const int_t*, const complex8_t*, int_t, const complex8_t*, int_t, complex8_t*, int_t);
/*-------------------------------------------------*/
template <typename T_Int>
__global__ void roll_xxptr_kernel(T_Int np, T_Int* xxxptr)
{
    T_Int tid = blockIdx.x * blockDim.x + threadIdx.x;

    if (tid == 0) {
        for(T_Int k = 0; k < np; k++) {
            xxxptr[k+1] = xxxptr[k+1] + xxxptr[k];
        } // k
    }
}
/*-------------------------------------------------*/
template <typename T_Int>
__global__ void unroll_xxptr_kernel(T_Int np, T_Int* xxxptr)
{
    T_Int tid = blockIdx.x * blockDim.x + threadIdx.x;

    if (tid == 0) {
        for(T_Int k = np-1; k > 0; k--) {
            xxxptr[k] = xxxptr[k-1];
        } // k
        xxxptr[0] = 0;
    }
}
/*-------------------------------------------------*/
template <typename T_Int>
__global__ void add_xxptr_kernel(T_Int np, 
                                 const T_Int* xxxptrA, const T_Int *xxxidxA,
                                 const T_Int* xxxptrB, const T_Int *xxxidxB,
                                 T_Int* xxxptrC)
{
    T_Int k = blockIdx.x * blockDim.x + threadIdx.x;

    if (k < np) {

        T_Int xbgnA = xxxptrA[k];
        T_Int xendA = xxxptrA[k+1];
        T_Int xbgnB = xxxptrB[k];
        T_Int xendB = xxxptrB[k+1];

        while (xbgnA < xendA && xbgnB < xendB) {

            T_Int xidxA = xxxidxA[xbgnA];
            T_Int xidxB = xxxidxB[xbgnB];

            if (xidxA < xidxB) {
                xxxptrC[k+1]++;
                xbgnA++;
            } else if (xidxA > xidxB) {
                xxxptrC[k+1]++;
                xbgnB++;
            } else {
                xxxptrC[k+1]++;
                xbgnA++;
                xbgnB++;
            }

        } // merge the two sorted lists of column indices

        while (xbgnA < xendA) {
            xxxptrC[k+1]++;
            xbgnA++;
        } // append rest of A

        while (xbgnB < xendB) {
            xxxptrC[k+1]++;
            xbgnB++;
        } // append rest of B

    } // dim check
}
/*-------------------------------------------------*/
template <typename T_Int>
void launch_add_xxptr_kernel(T_Int np, 
                             const T_Int* xxxptrA, const T_Int *xxxidxA,
                             const T_Int* xxxptrB, const T_Int *xxxidxB,
                             T_Int* nnzC, T_Int* xxxptrC)
{
    if(np <= 0) return;

    memSetZero(np + 1, xxxptrC);

    Grid1D grid(np);
    add_xxptr_kernel<<<grid.numBlocks(), 
                       grid.threadsPerBlock()>>>(np, 
                                                 xxxptrA, xxxidxA, 
                                                 xxxptrB, xxxidxB, 
                                                 xxxptrC);

    roll_xxptr_kernel<<<1, 1>>>(np, xxxptrC);

    memCopyD2H(1, xxxptrC + np, nnzC);
    syncDevice();
}
/*-------------------------------------------------*/
template void launch_add_xxptr_kernel(int_t, const int_t*, const int_t*, const int_t*, const int_t*, int_t*, int_t*);
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
__global__ void add_kernel(T_Int np, T_Scalar alpha, T_Scalar beta,
                           const T_Int* xxxptrA, const T_Int *xxxidxA, const T_Scalar* valuesA,
                           const T_Int* xxxptrB, const T_Int *xxxidxB, const T_Scalar* valuesB,
                           T_Int* xxxptrC, T_Int *xxxidxC, T_Scalar* valuesC)
{
    T_Int k = blockIdx.x * blockDim.x + threadIdx.x;

    if (k < np) {

        T_Int xbgnA = xxxptrA[k];
        T_Int xendA = xxxptrA[k+1];
        T_Int xbgnB = xxxptrB[k];
        T_Int xendB = xxxptrB[k+1];

        while (xbgnA < xendA && xbgnB < xendB) {

            T_Int xidxA = xxxidxA[xbgnA];
            T_Int xidxB = xxxidxB[xbgnB];

            if (xidxA < xidxB) {
                xxxidxC[xxxptrC[k]] = xidxA;
                valuesC[xxxptrC[k]] = alpha * valuesA[xbgnA];
                xxxptrC[k]++;
                xbgnA++;
            } else if (xidxA > xidxB) {
                xxxidxC[xxxptrC[k]] = xidxB;
                valuesC[xxxptrC[k]] = beta * valuesB[xbgnB];
                xxxptrC[k]++;
                xbgnB++;
            } else {
                xxxidxC[xxxptrC[k]] = xidxA; // or xidxB, they are equal
                valuesC[xxxptrC[k]] = alpha * valuesA[xbgnA] + beta * valuesB[xbgnB];
                xxxptrC[k]++;
                xbgnA++;
                xbgnB++;
            }

        } // merge the two sorted lists of column indices

        while (xbgnA < xendA) {
            xxxidxC[xxxptrC[k]] = xxxidxA[xbgnA];
            valuesC[xxxptrC[k]] = alpha * valuesA[xbgnA];
            xxxptrC[k]++;
            xbgnA++;
        } // append rest of A

        while (xbgnB < xendB) {
            xxxidxC[xxxptrC[k]] = xxxidxB[xbgnB];
            valuesC[xxxptrC[k]] = beta * valuesB[xbgnB];
            xxxptrC[k]++;
            xbgnB++;
        } // append rest of B

    } // dim check
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
void launch_add_kernel(T_Int np, T_Scalar alpha, T_Scalar beta,
                       const T_Int* xxxptrA, const T_Int *xxxidxA, const T_Scalar* valuesA, 
                       const T_Int* xxxptrB, const T_Int *xxxidxB, const T_Scalar* valuesB,
                       T_Int* xxxptrC, T_Int *xxxidxC, T_Scalar* valuesC)
{
    if(np <= 0) return;
    
    Grid1D grid(np);
    add_kernel<<<grid.numBlocks(), 
                 grid.threadsPerBlock()>>>(np, alpha, beta, 
                                           xxxptrA, xxxidxA, valuesA, 
                                           xxxptrB, xxxidxB, valuesB, 
                                           xxxptrC, xxxidxC, valuesC);

    unroll_xxptr_kernel<<<1, 1>>>(np, xxxptrC);

    syncDevice();
}
/*-------------------------------------------------*/
template void launch_add_kernel(int_t, real_t    , real_t    , const int_t*, const int_t*, const real_t    *, const int_t*, const int_t*, const real_t    *, int_t*, int_t*, real_t    *);
template void launch_add_kernel(int_t, real4_t   , real4_t   , const int_t*, const int_t*, const real4_t   *, const int_t*, const int_t*, const real4_t   *, int_t*, int_t*, real4_t   *);
template void launch_add_kernel(int_t, complex_t , complex_t , const int_t*, const int_t*, const complex_t *, const int_t*, const int_t*, const complex_t *, int_t*, int_t*, complex_t *);
template void launch_add_kernel(int_t, complex8_t, complex8_t, const int_t*, const int_t*, const complex8_t*, const int_t*, const int_t*, const complex8_t*, int_t*, int_t*, complex8_t*);
/*-------------------------------------------------*/
} // namespace csx
} // namespace blk
} // namespace culite
/*-------------------------------------------------*/
