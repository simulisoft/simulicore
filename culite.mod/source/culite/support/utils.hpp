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

#ifndef CULITE_UTILS_HPP_
#define CULITE_UTILS_HPP_

/** 
 * @file
 */

#include <cstddef>
#include <cuda_runtime.h>

#include "culite/types/traits.hpp"

/*-------------------------------------------------*/
namespace culite {
/*-------------------------------------------------*/

void memCopyX2X(std::size_t count, const void *src, void *dest, cudaMemcpyKind kind);

template <typename T>
inline void memCopyD2D(std::size_t n, const T *src, T *dest)
{
	memCopyX2X(n * sizeof(T), src, dest, cudaMemcpyDeviceToDevice);
}

template <typename T>
inline void memCopyH2D(std::size_t n, const typename TypeTraits<T>::host_type *src, T *dest)
{
	memCopyX2X(n * sizeof(T), src, dest, cudaMemcpyHostToDevice);
}

template <typename T>
inline void memCopyD2H(std::size_t n, const T *src, typename TypeTraits<T>::host_type *dest)
{
	memCopyX2X(n * sizeof(T), src, dest, cudaMemcpyDeviceToHost);
}

template <typename T>
inline void memCopyH2H(std::size_t n, const T *src, T *dest)
{
	memCopyX2X(n * sizeof(T), src, dest, cudaMemcpyHostToHost);
}

/*-------------------------------------------------*/
} // namespace culite
/*-------------------------------------------------*/

#endif // CULITE_UTILS_HPP_
