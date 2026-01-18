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

#ifndef CULITE_BULK_DNS1D_IMPL_CUH_
#define CULITE_BULK_DNS1D_IMPL_CUH_

/**
 * @file
 */
#include <cuda_runtime.h>

#include "culite/types/scalar.hpp"

/*-------------------------------------------------*/
namespace culite {
namespace blk {
namespace dns {
/*-------------------------------------------------*/

template <typename T_Scalar>
void launch_conjugate_kernel_1d(std::size_t n, T_Scalar* x);

/*-------------------------------------------------*/
} // namespace dns
} // namespace blk
} // namespace culite
/*-------------------------------------------------*/

#endif // CULITE_BULK_DNS1D_IMPL_CUH_
