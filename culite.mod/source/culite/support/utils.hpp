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

/**
 * @brief Synchronizes the CUDA device.
 * @details Blocks until the device has completed all preceding requested tasks.
 */
void syncDevice();

/**
 * @ingroup culite_module_index_copying
 * @brief Copies memory between different memory spaces.
 * @param count Number of bytes to copy.
 * @param src Source memory pointer.
 * @param dest Destination memory pointer.
 * @param kind Type of memory copy operation (device-to-device, host-to-device, etc.).
 */
void memCopyX2X(std::size_t count, const void *src, void *dest, cudaMemcpyKind kind);

/**
 * @ingroup culite_module_index_copying
 * @brief Copies a vector from device memory to device memory.
 * @tparam T_Scalar The scalar type of the vector elements.
 * @param n Number of elements to copy.
 * @param src Source device memory pointer.
 * @param dest Destination device memory pointer.
 */
template <typename T_Scalar>
inline void memCopyD2D(std::size_t n, const T_Scalar *src, T_Scalar *dest)
{
	memCopyX2X(n * sizeof(T_Scalar), src, dest, cudaMemcpyDeviceToDevice);
}

/**
 * @ingroup culite_module_index_copying
 * @brief Copies a vector from host memory to device memory.
 * @tparam T_Scalar The scalar type of the vector elements.
 * @param n Number of elements to copy.
 * @param src Source host memory pointer.
 * @param dest Destination device memory pointer.
 */
template <typename T_Scalar>
inline void memCopyH2D(std::size_t n, const T_Scalar *src, T_Scalar *dest)
{
	memCopyX2X(n * sizeof(T_Scalar), src, dest, cudaMemcpyHostToDevice);
}

/**
 * @ingroup culite_module_index_copying
 * @brief Copies a vector from device memory to host memory.
 * @tparam T_Scalar The scalar type of the vector elements.
 * @param n Number of elements to copy.
 * @param src Source device memory pointer.
 * @param dest Destination host memory pointer.
 */
template <typename T_Scalar>
inline void memCopyD2H(std::size_t n, const T_Scalar *src, T_Scalar *dest)
{
	memCopyX2X(n * sizeof(T_Scalar), src, dest, cudaMemcpyDeviceToHost);
}

/**
 * @ingroup culite_module_index_copying
 * @brief Copies a vector from host memory to host memory.
 * @tparam T_Scalar The scalar type of the vector elements.
 * @param n Number of elements to copy.
 * @param src Source host memory pointer.
 * @param dest Destination host memory pointer.
 */
template <typename T_Scalar>
inline void memCopyH2H(std::size_t n, const T_Scalar *src, T_Scalar *dest)
{
	memCopyX2X(n * sizeof(T_Scalar), src, dest, cudaMemcpyHostToHost);
}

/*-------------------------------------------------*/

/**
 * @ingroup culite_module_index_copying
 * @brief Copies a vector from host memory (cla3p type) to device memory.
 * @tparam T_Scalar The scalar type of the vector elements.
 * @param n Number of elements to copy.
 * @param src Source host memory pointer with cla3p type.
 * @param dest Destination device memory pointer.
 */
template <typename T_Scalar>
inline void memCopyH2D(std::size_t n, const typename TypeTraits<T_Scalar>::cla3p_type *src, T_Scalar *dest)
{
	memCopyX2X(n * sizeof(T_Scalar), src, dest, cudaMemcpyHostToDevice);
}

/**
 * @ingroup culite_module_index_copying
 * @brief Copies a vector from device memory to host memory (cla3p type).
 * @tparam T_Scalar The scalar type of the vector elements.
 * @param n Number of elements to copy.
 * @param src Source device memory pointer.
 * @param dest Destination host memory pointer with cla3p type.
 */
template <typename T_Scalar>
inline void memCopyD2H(std::size_t n, const T_Scalar *src, typename TypeTraits<T_Scalar>::cla3p_type *dest)
{
	memCopyX2X(n * sizeof(T_Scalar), src, dest, cudaMemcpyDeviceToHost);
}

/*-------------------------------------------------*/

/**
 * @ingroup culite_module_index_copying
 * @brief Copies a 2D matrix between different memory spaces.
 * @param elemSize Size of each element in bytes.
 * @param m Number of rows to copy.
 * @param n Number of columns to copy.
 * @param src Source memory pointer.
 * @param lds Leading dimension of the source matrix.
 * @param dest Destination memory pointer.
 * @param ldd Leading dimension of the destination matrix.
 * @param kind Type of memory copy operation (device-to-device, host-to-device, etc.).
 */
void memCopyX2X(std::size_t elemSize, std::size_t m, std::size_t n, 
	            const void *src, std::size_t lds, 
				void *dest, std::size_t ldd, cudaMemcpyKind kind);

/**
 * @ingroup culite_module_index_copying
 * @brief Copies a 2D matrix from device memory to device memory.
 * @tparam T_Scalar The scalar type of the matrix elements.
 * @param m Number of rows to copy.
 * @param n Number of columns to copy.
 * @param src Source device memory pointer.
 * @param lds Leading dimension of the source matrix.
 * @param dest Destination device memory pointer.
 * @param ldd Leading dimension of the destination matrix.
 */
template <typename T_Scalar>
inline void memCopyD2D(std::size_t m, std::size_t n, 
	                   const T_Scalar *src, std::size_t lds, 
					   T_Scalar *dest, std::size_t ldd)
{
	memCopyX2X(sizeof(T_Scalar), m, n, src, lds, dest, ldd, cudaMemcpyDeviceToDevice);
}

/**
 * @ingroup culite_module_index_copying
 * @brief Copies a 2D matrix from host memory to device memory.
 * @tparam T_Scalar The scalar type of the matrix elements.
 * @param m Number of rows to copy.
 * @param n Number of columns to copy.
 * @param src Source host memory pointer.
 * @param lds Leading dimension of the source matrix.
 * @param dest Destination device memory pointer.
 * @param ldd Leading dimension of the destination matrix.
 */
template <typename T_Scalar>
inline void memCopyH2D(std::size_t m, std::size_t n, 
	                   const T_Scalar *src, std::size_t lds, 
					   T_Scalar *dest, std::size_t ldd)
{
	memCopyX2X(sizeof(T_Scalar), m, n, src, lds, dest, ldd, cudaMemcpyHostToDevice);
}

/**
 * @ingroup culite_module_index_copying
 * @brief Copies a 2D matrix from device memory to host memory.
 * @tparam T_Scalar The scalar type of the matrix elements.
 * @param m Number of rows to copy.
 * @param n Number of columns to copy.
 * @param src Source device memory pointer.
 * @param lds Leading dimension of the source matrix.
 * @param dest Destination host memory pointer.
 * @param ldd Leading dimension of the destination matrix.
 */
template <typename T_Scalar>
inline void memCopyD2H(std::size_t m, std::size_t n, 
	                   const T_Scalar *src, std::size_t lds, 
					   T_Scalar *dest, std::size_t ldd)
{
	memCopyX2X(sizeof(T_Scalar), m, n, src, lds, dest, ldd, cudaMemcpyDeviceToHost);
}

/**
 * @ingroup culite_module_index_copying
 * @brief Copies a 2D matrix from host memory to host memory.
 * @tparam T_Scalar The scalar type of the matrix elements.
 * @param m Number of rows to copy.
 * @param n Number of columns to copy.
 * @param src Source host memory pointer.
 * @param lds Leading dimension of the source matrix.
 * @param dest Destination host memory pointer.
 * @param ldd Leading dimension of the destination matrix.
 */
template <typename T_Scalar>
inline void memCopyH2H(std::size_t m, std::size_t n, 
	                   const T_Scalar *src, std::size_t lds, 
					   T_Scalar *dest, std::size_t ldd)
{
	memCopyX2X(sizeof(T_Scalar), m, n, src, lds, dest, ldd, cudaMemcpyHostToHost);
}

/*-------------------------------------------------*/

/**
 * @ingroup culite_module_index_copying
 * @brief Copies a 2D matrix from host memory (cla3p type) to device memory.
 * @tparam T_Scalar The scalar type of the matrix elements.
 * @param m Number of rows to copy.
 * @param n Number of columns to copy.
 * @param src Source host memory pointer with cla3p type.
 * @param lds Leading dimension of the source matrix.
 * @param dest Destination device memory pointer.
 * @param ldd Leading dimension of the destination matrix.
 */
template <typename T_Scalar>
inline void memCopyH2D(std::size_t m, std::size_t n, 
	                   const typename TypeTraits<T_Scalar>::cla3p_type *src, std::size_t lds, 
					   T_Scalar *dest, std::size_t ldd)
{
	memCopyX2X(sizeof(T_Scalar), m, n, src, lds, dest, ldd, cudaMemcpyHostToDevice);
}

/**
 * @ingroup culite_module_index_copying
 * @brief Copies a 2D matrix from device memory to host memory (cla3p type).
 * @tparam T_Scalar The scalar type of the matrix elements.
 * @param m Number of rows to copy.
 * @param n Number of columns to copy.
 * @param src Source device memory pointer.
 * @param lds Leading dimension of the source matrix.
 * @param dest Destination host memory pointer with cla3p type.
 * @param ldd Leading dimension of the destination matrix.
 */
template <typename T_Scalar>
inline void memCopyD2H(std::size_t m, std::size_t n, 
	                   const T_Scalar *src, std::size_t lds, 
					   typename TypeTraits<T_Scalar>::cla3p_type *dest, std::size_t ldd)
{
	memCopyX2X(sizeof(T_Scalar), m, n, src, lds, dest, ldd, cudaMemcpyDeviceToHost);
}

/*-------------------------------------------------*/

/**
 * @ingroup culite_module_index_copying
 * @brief Sets a vector in device memory to zero.
 * @param n Number of elements to set to zero.
 * @param x Pointer to device memory vector.
 * @param elemSize Size of each element in bytes.
 */
void memSetZeroX(std::size_t n, void *x, std::size_t elemSize);

/**
 * @ingroup culite_module_index_copying
 * @brief Sets a vector in device memory to zero.
 * @tparam T_Scalar The scalar type of the vector elements.
 * @param n Number of elements to set to zero.
 * @param x Pointer to device memory vector.
 */
template <typename T_Scalar>
void memSetZero(std::size_t n, T_Scalar *x)
{
    memSetZeroX(n, x, sizeof(T_Scalar));
}

/**
 * @ingroup culite_module_index_copying
 * @brief Sets a 2D matrix in device memory to zero.
 * @param m Number of rows to set to zero.
 * @param n Number of columns to set to zero.
 * @param a Pointer to device memory matrix.
 * @param lda Leading dimension of the matrix.
 * @param elemSize Size of each element in bytes.
 */
void memSetZeroX(std::size_t m, std::size_t n, void *a, std::size_t lda, std::size_t elemSize);

/**
 * @ingroup culite_module_index_copying
 * @brief Sets a 2D matrix in device memory to zero.
 * @tparam T_Scalar The scalar type of the matrix elements.
 * @param m Number of rows to set to zero.
 * @param n Number of columns to set to zero.
 * @param a Pointer to device memory matrix.
 * @param lda Leading dimension of the matrix.
 */
template <typename T_Scalar>
void memSetZero(std::size_t m, std::size_t n, T_Scalar *a, std::size_t lda)
{
    memSetZeroX(m, n, a, lda, sizeof(T_Scalar));
}

/*-------------------------------------------------*/
} // namespace culite
/*-------------------------------------------------*/

#endif // CULITE_UTILS_HPP_
