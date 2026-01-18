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
#include "culite/types/integer.hpp"

/*-------------------------------------------------*/
namespace culite {
namespace blk {
namespace dns {
/*-------------------------------------------------*/

template <typename T_Scalar>
void launch_scale_matrix_kernel(int_t m, int_t n, T_Scalar* v, int_t ldv, T_Scalar alpha);

/*-------------------------------------------------*/
} // namespace dns
} // namespace blk
} // namespace culite
/*-------------------------------------------------*/

#endif // CULITE_BULK_DNS2D_IMPL_HPP_
