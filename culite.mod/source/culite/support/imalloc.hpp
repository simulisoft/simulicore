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

void* device_alloc(std::size_t size);
void device_free(void *ptr) noexcept;

template <typename T>
T* device_alloc_t(std::size_t n) {
	return static_cast<T*>(device_alloc(n * sizeof(T)));
}

void* pinned_alloc(std::size_t size);
void pinned_free(void *ptr) noexcept;

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
