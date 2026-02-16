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

#ifndef CULITE_BULK_CSX_HPP_
#define CULITE_BULK_CSX_HPP_

/**
 * @file
 */

#include "culite/bulk/csx_impl.hpp"

/*-------------------------------------------------*/
namespace culite {
namespace blk {
namespace csx {
/*-------------------------------------------------*/

// 
// A(np x np)
// x(np)
// y(np)
// y = y + alpha * diag(A) * x
//
template <typename T_Int, typename T_Scalar>
void diag_times_vec(const T_Scalar* alpha, int_t np, 
                        const T_Int* xxxptr, 
                        const T_Int* xxxidx, 
                        const T_Scalar* values, 
                        const T_Scalar*x, T_Int incx,
                        T_Scalar *y, T_Int incy)
{
	if(np > 0) {
		launch_diag_times_vec_kernel(alpha, np, xxxptr, xxxidx, values, x, incx, y, incy);
	}
}

/*-------------------------------------------------*/

// 
// A(np x np)
// B(np x nc)
// C(np x nc)
// C = C + alpha * diag(A) * B
//
template <typename T_Int, typename T_Scalar>
void diag_times_mat(const T_Scalar* alpha, int_t np, 
                    const T_Int* xxxptr, 
                    const T_Int* xxxidx, 
                    const T_Scalar* values, 
                    T_Int nc, 
                    const T_Scalar*b, T_Int ldb, 
                    T_Scalar *c, T_Int ldc)
{
	if(nc > 0 && np > 0) {
		launch_diag_times_mat_kernel(alpha, np, xxxptr, xxxidx, values, nc, b, ldb, c, ldc);
	}
}

/*-------------------------------------------------*/

template <typename T_Int>
void add_xxptr(T_Int np, 
               const T_Int* xxxptrA, const T_Int *xxxidxA,
               const T_Int* xxxptrB, const T_Int *xxxidxB,
               T_Int* nnzC, T_Int* xxxptrC)
{
    if(np > 0) {
        launch_add_xxptr_kernel(np, xxxptrA, xxxidxA, xxxptrB, xxxidxB, nnzC, xxxptrC);
    }
}

template <typename T_Int, typename T_Scalar>
void add(T_Int np, T_Scalar alpha, T_Scalar beta,
         const T_Int* xxxptrA, const T_Int *xxxidxA, const T_Scalar* valuesA, 
         const T_Int* xxxptrB, const T_Int *xxxidxB, const T_Scalar* valuesB,
         T_Int* xxxptrC, T_Int *xxxidxC, T_Scalar* valuesC)
{
    if(np > 0) {
        launch_add_kernel(np, alpha, beta,
                          xxxptrA, xxxidxA, valuesA, 
                          xxxptrB, xxxidxB, valuesB, 
                          xxxptrC, xxxidxC, valuesC);
    }
}

/*-------------------------------------------------*/
} // namespace csx
} // namespace blk
} // namespace culite
/*-------------------------------------------------*/

#endif // CULITE_BULK_CSX_HPP_
