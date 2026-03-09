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

#ifndef CLA3P_IMALLOC_HPP_
#define CLA3P_IMALLOC_HPP_

/** 
 * @file
 */

#include "cla3p/types.hpp"

/*-------------------------------------------------*/
namespace cla3p {
/*-------------------------------------------------*/

/**
 * @ingroup cla3p_module_index_allocators
 * @brief Default allocator for raw byte storage.
 * @details Allocates @p size bytes of uninitialized storage.
 * @param[in] size The requested size in bytes.
 * @return Pointer to the allocated storage on success, otherwise @c nullptr.
 */
void* i2malloc(std::size_t size);

/**
 * @ingroup cla3p_module_index_allocators
 * @brief Default zero-initializing allocator for raw byte storage.
 * @details Allocates @p nmemb objects of @p size bytes each and zero-initializes the memory.
 * @param[in] nmemb The number of elements to allocate.
 * @param[in] size The size of each element in bytes.
 * @return Pointer to the allocated storage on success, otherwise @c nullptr.
 */
void* i2calloc(std::size_t nmemb, std::size_t size);

/**
 * @ingroup cla3p_module_index_allocators
 * @brief Resizes a previously allocated memory block.
 * @details Reallocates the memory referenced by @p ptr to @p size bytes, preserving existing data up to the new size.
 * @param[in] ptr Pointer to the memory block to resize (may be @c nullptr).
 * @param[in] size The new size of the allocation in bytes.
 * @return Pointer to the resized memory on success, otherwise @c nullptr.
 */
void* i2realloc(void *ptr, std::size_t size);

/**
 * @ingroup cla3p_module_index_allocators
 * @brief Type-safe allocator for uninitialized storage.
 * @details Allocates @p nmemb elements of type @p T_Elem without initialization.
 * @param[in] nmemb The number of elements to allocate.
 * @return Pointer to the allocated storage on success, otherwise @c nullptr.
 */
template <typename T_Elem>
T_Elem* i_malloc_t(std::size_t nmemb)
{
    return static_cast<T_Elem*>(i2malloc(nmemb * sizeof(T_Elem)));
}

/**
 * @ingroup cla3p_module_index_allocators
 * @brief Type-safe allocator for zero-initialized storage.
 * @details Allocates @p nmemb elements of type @p T_Elem and initializes all bits to zero.
 * @param[in] nmemb The number of elements to allocate.
 * @return Pointer to the allocated storage on success, otherwise @c nullptr.
 */
template <typename T_Elem>
T_Elem* i_calloc_t(std::size_t nmemb)
{
    return static_cast<T_Elem*>(i2calloc(nmemb, sizeof(T_Elem)));
}

/**
 * @ingroup cla3p_module_index_allocators
 * @brief Type-safe reallocator for typed storage.
 * @details Resizes the memory block referenced by @p ptr to store @p nmemb elements of type @p T_Elem,
 *          preserving existing data up to the new size.
 * @param[in] ptr Pointer to the memory block to resize (may be @c nullptr).
 * @param[in] nmemb The new number of elements requested.
 * @return Pointer to the resized storage on success, otherwise @c nullptr.
 */
template <typename T_Elem>
T_Elem* i_realloc_t(T_Elem* ptr, std::size_t nmemb)
{
    return static_cast<T_Elem*>(i2realloc(ptr, nmemb * sizeof(T_Elem)));
}

/**
 * @ingroup cla3p_module_index_allocators
 * @brief Default deallocator for cla3p allocations.
 * @details Releases memory previously allocated by @ref i_malloc_t(), @ref i_calloc_t(), or @ref i_realloc_t().
 * @param[in,out] ptr Pointer to the memory block to deallocate (may be @c nullptr).
 */
void i_free(void *ptr);

/*-------------------------------------------------*/
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_IMALLOC_HPP_
