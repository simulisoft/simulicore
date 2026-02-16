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

#ifndef CULITE_BULK_CSX_IMPL_HPP_
#define CULITE_BULK_CSX_IMPL_HPP_

/**
 * @file
 */

/*-------------------------------------------------*/
namespace culite {
namespace blk {
namespace csx {
/*-------------------------------------------------*/

template <typename T_Int, typename T_Scalar>
void launch_diag_times_vec_kernel(const T_Scalar* alpha, T_Int np, 
                                  const T_Int* xxxptr, 
                                  const T_Int* xxxidx, 
                                  const T_Scalar* values, 
                                  const T_Scalar*x, T_Int incx,
                                  T_Scalar *y, T_Int incy);

template <typename T_Int, typename T_Scalar>
void launch_diag_times_mat_kernel(const T_Scalar* alpha, T_Int np, 
                                  const T_Int* xxxptr, 
                                  const T_Int* xxxidx, 
                                  const T_Scalar* values, 
                                  T_Int nc,
                                  const T_Scalar*b, T_Int ldb,
                                  T_Scalar *c, T_Int ldc);

template <typename T_Int>
void launch_add_xxptr_kernel(T_Int np, 
                             const T_Int* xxxptrA, const T_Int *xxxidxA,
                             const T_Int* xxxptrB, const T_Int *xxxidxB,
                             T_Int* nnzC, T_Int* xxxptrC);

template <typename T_Int, typename T_Scalar>
void launch_add_kernel(T_Int np, T_Scalar alpha, T_Scalar beta,
                       const T_Int* xxxptrA, const T_Int *xxxidxA, const T_Scalar* valuesA, 
                       const T_Int* xxxptrB, const T_Int *xxxidxB, const T_Scalar* valuesB,
                       T_Int* xxxptrC, T_Int *xxxidxC, T_Scalar* valuesC);


/*-------------------------------------------------*/
} // namespace csx
} // namespace blk
} // namespace culite
/*-------------------------------------------------*/

#endif // CULITE_BULK_CSX_IMPL_HPP_
