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

#ifndef CULITE_IMALLOC_HPP_
#define CULITE_IMALLOC_HPP_

/** 
 * @file
 */

#include <cstddef>

/*-------------------------------------------------*/
namespace culite {
/*-------------------------------------------------*/

/**
 * @ingroup culite_module_index_allocators
 * @brief Allocates memory on the device.
 * @details This function reserves a block of memory of the specified size on the 
 *          underlying device.
 * @param[in] size The number of bytes to allocate.
 * @return A pointer to the allocated memory on success; @c nullptr if the allocation fails.
 * @note Memory allocated by this function must be freed using the corresponding device free function.
 * @see device_free()
 */
void* device_alloc(std::size_t size);

/**
 * @ingroup culite_module_index_allocators
 * @brief Frees a block of memory on the device.
 * @details This function releases the memory previously allocated by @ref device_alloc.
 *          If @p ptr is @c nullptr, no action is performed.
 * @param[in] ptr Pointer to the device memory to be freed.
 * @note This function is marked @c noexcept and is guaranteed not to throw 
 *       exceptions during the deallocation process.
 * @see device_alloc()
 */
void device_free(void *ptr) noexcept;

/**
 * @ingroup culite_module_index_allocators
 * @brief Allocates typed memory on the device.
 * @details This function reserves a block of memory for @p n elements of type @p T 
 *          on the underlying device.
 * @tparam T The type of elements to allocate.
 * @param[in] n The number of elements to allocate.
 * @return A pointer to the allocated memory on success; @c nullptr if the allocation fails.
 * @note Memory allocated by this function must be freed using the corresponding device free function.
 * @see device_free()
 */
template <typename T>
T* device_alloc_t(std::size_t n) {
	return static_cast<T*>(device_alloc(n * sizeof(T)));
}

/**
 * @ingroup culite_module_index_allocators
 * @brief Allocates page-locked (pinned) host memory.
 * @details Reserves a block of host memory that is page-locked and accessible to the device. 
 * @param[in] size The number of bytes to allocate.
 * @return A pointer to the pinned host memory on success; @c nullptr if the allocation fails.
 * @note Must be freed using the corresponding pinned memory deallocation function.
 * @see pinned_free()
 */

void* pinned_alloc(std::size_t size);
/**
 * @ingroup culite_module_index_allocators
 * @brief Frees a block of pinned host memory.
 * @details Releases the memory previously allocated by @ref pinned_alloc. 
 *          If @p ptr is @c nullptr, no action is performed.
 * @param[in] ptr Pointer to the pinned host memory to be freed.
 * @note This function is marked @c noexcept and is guaranteed not to throw 
 *       exceptions during the deallocation process.
 * @see pinned_alloc()
 */
void pinned_free(void *ptr) noexcept;

/**
 * @ingroup culite_module_index_allocators
 * @brief Allocates typed pinned host memory.
 * @details This function reserves a block of host memory for @p n elements of type @p T 
 *          that is page-locked and accessible to the device.
 * @tparam T The type of elements to allocate.
 * @param[in] n The number of elements to allocate.
 * @return A pointer to the allocated pinned host memory on success; @c nullptr if the allocation fails.
 * @note Must be freed using the corresponding pinned memory deallocation function.
 * @see pinned_free()
 */
template <typename T>
T* pinned_alloc_t(std::size_t n) {
	return static_cast<T*>(pinned_alloc(n * sizeof(T)));
}

/*-------------------------------------------------*/

/* // TODO: add doc ?
 * @ingroup culite_module_index_allocators
 * @brief The device memory allocator.
 */
template <typename T>
class DeviceAllocator {

	public:
    	using value_type = T;

    	DeviceAllocator() noexcept = default;

    	template <typename U>
    	DeviceAllocator(const DeviceAllocator<U>&) noexcept {}

    	T* allocate(std::size_t n) { return static_cast<T*>(device_alloc(n * sizeof(T))); }
    	void deallocate(T* p, std::size_t n) noexcept { device_free(static_cast<void*>(p)); }
};

template <typename T, typename U>
bool operator==(const DeviceAllocator<T>&, const DeviceAllocator<U>&) { return true; }

template <typename T, typename U>
bool operator!=(const DeviceAllocator<T>&, const DeviceAllocator<U>&) { return false; }

/*-------------------------------------------------*/

/* // TODO: add doc ?
 * @ingroup culite_module_index_allocators
 * @brief The pinned memory allocator.
 */
template <typename T>
class PinnedAllocator {

	public:
    	using value_type = T;

    	PinnedAllocator() noexcept = default;

    	template <typename U>
    	PinnedAllocator(const PinnedAllocator<U>&) noexcept {}

    	T* allocate(std::size_t n) { return static_cast<T*>(pinned_alloc(n * sizeof(T))); }
    	void deallocate(T* p, std::size_t n) noexcept { pinned_free(p); }
};

template <typename T, typename U>
bool operator==(const PinnedAllocator<T>&, const PinnedAllocator<U>&) { return true; }

template <typename T, typename U>
bool operator!=(const PinnedAllocator<T>&, const PinnedAllocator<U>&) { return false; }

/*-------------------------------------------------*/
} // namespace culite
/*-------------------------------------------------*/

#endif // CULITE_IMALLOC_HPP_
