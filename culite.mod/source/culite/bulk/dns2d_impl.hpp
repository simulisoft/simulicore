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

#ifndef CULITE_BULK_DNS2D_IMPL_HPP_
#define CULITE_BULK_DNS2D_IMPL_HPP_

/**
 * @file
 */

#include "culite/types/traits.hpp"
#include "culite/types/integer.hpp"
#include "culite/types/enums.hpp"

// forwards
#include "culite/types/cla3p_forwards.hpp"

/*-------------------------------------------------*/
namespace culite {
namespace blk {
namespace dns {
/*-------------------------------------------------*/

template <typename T_Scalar>
void launch_copy_kernel(uplo_t uplo, int_t m, int_t n, 
                        const T_Scalar *a, int_t lda, 
                        T_Scalar *b, int_t ldb);

template <typename T_Scalar>
void launch_fill_kernel(uplo_t uplo, int_t m, int_t n, T_Scalar *a, int_t lda, T_Scalar val);

template <typename T_Scalar>
void launch_scale_kernel(uplo_t uplo, int_t m, int_t n, T_Scalar* a, int_t lda, T_Scalar alpha);

template <typename T_Scalar>
typename TypeTraits<T_Scalar>::real_type
launch_matrix_one_norm_kernel(int_t m, int_t n, const T_Scalar* a, int_t lda);

template <typename T_Scalar>
typename TypeTraits<T_Scalar>::real_type
launch_matrix_inf_norm_kernel(int_t m, int_t n, const T_Scalar* a, int_t lda);

template <typename T_Scalar>
typename TypeTraits<T_Scalar>::real_type
launch_matrix_max_norm_kernel(int_t m, int_t n, const T_Scalar* a, int_t lda);

template <typename T_Scalar>
typename TypeTraits<T_Scalar>::real_type
launch_matrix_fro_norm_kernel(int_t m, int_t n, const T_Scalar* a, int_t lda);

template <typename T_Scalar>
void launch_get_real_kernel(uplo_t uplo, int_t m, int_t n, 
                            const T_Scalar* a, int_t lda, 
                            typename TypeTraits<T_Scalar>::real_type* b, int_t ldb);

template <typename T_Scalar>
void launch_get_imag_kernel(uplo_t uplo, int_t m, int_t n, 
                            const T_Scalar* a, int_t lda, 
                            typename TypeTraits<T_Scalar>::real_type* b, int_t ldb);

template <typename T_Scalar>
void launch_conjugate_kernel(uplo_t uplo, int_t m, int_t n, T_Scalar* a, int_t lda);

template <typename T_Scalar>
void launch_geev_calculate_complex_eigenvectors_kernel(int_t n, 
                                                       const typename TypeTraits<T_Scalar>::real_type* w, 
                                                       const typename TypeTraits<T_Scalar>::real_type* vr, int_t ldvr,
                                                       T_Scalar *vc, int_t ldvc);

/*-------------------------------------------------*/
} // namespace dns
} // namespace blk
} // namespace culite
/*-------------------------------------------------*/

#endif // CULITE_BULK_DNS2D_IMPL_HPP_
