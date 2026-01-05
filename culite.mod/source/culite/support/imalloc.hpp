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

#include "culite/types/enums.hpp"

/*-------------------------------------------------*/
namespace culite {
/*-------------------------------------------------*/

/**
 * @ingroup culite_module_index_allocators
 * @brief Detects the allocation type of a given memory pointer.
 * @details This function queries the CUDA runtime to determine the allocation type of the
 *          provided pointer. It identifies whether the memory is allocated on the device,
 *          pinned on the host, managed (unified memory), or unregistered (unknown origin).
 * @param[in] ptr Pointer to the memory whose allocation type is to be detected.
 * @return The @ref culite::alloc_t enumeration value indicating the allocation type:
 *         - @ref culite::alloc_t::Device for device memory
 *         - @ref culite::alloc_t::Pinned for pinned host memory
 *         - @ref culite::alloc_t::Managed for managed (unified) memory
 *         - @ref culite::alloc_t::Unregistered for unknown or unregistered memory
 * @note This function incurs significant performance overhead due to CUDA runtime queries
 *       and should be avoided in performance-critical loops. Cache the result when possible.
 * @see alloc_t, auto_free()
 */
alloc_t detect_allocation_type(const void* ptr);

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

/**
 * @ingroup culite_module_index_allocators
 * @brief Automatically detects and frees memory based on allocation type.
 * @details This function detects the allocation type of the given pointer (device memory,
 *          pinned host memory, or managed memory) and automatically calls the appropriate
 *          deallocation function. If @p ptr is @c nullptr or if the allocation type is
 *          unregistered (unknown origin), no action is performed.
 * @param[in] ptr Pointer to the memory to be freed.
 * @note This function incurs performance overhead due to allocation type detection and
 *       should be avoided in performance-critical loops. Prefer using explicit
 *       @ref device_free() or @ref pinned_free() when the allocation type is known.
 * @see device_free(), pinned_free(), detect_allocation_type()
 */
void auto_free(void *ptr);

/**
 * @ingroup culite_module_index_allocators
 * @brief A RAII wrapper for untyped device memory allocation.
 * @details This class manages a dynamically allocated buffer of raw bytes on the device.
 *          It automatically handles memory allocation and deallocation following RAII principles,
 *          ensuring that device memory is properly freed when the buffer goes out of scope.
 *          Unlike @ref DeviceBuffer, this class manages untyped memory (void*).
 */
class DeviceBufferVoid {

	public:
		/**
		 * @brief Default constructor.
		 * @details Constructs an empty device buffer with no allocated memory.
		 */
		DeviceBufferVoid() 
		{ 
			defaults(); 
		}

		/**
		 * @brief Constructs a device buffer with specified capacity.
		 * @details Constructs a device buffer and reserves memory for @p n bytes on the device.
		 * @param[in] n The number of bytes to allocate.
		 */
		explicit DeviceBufferVoid(std::size_t n) 
			: DeviceBufferVoid()
		{
			reserve(n);
		}

		/**
		 * @brief Destructor.
		 * @details Automatically frees the allocated device memory.
		 */
		~DeviceBufferVoid()
		{
			clear();
		}
		
		/**
		 * @brief Reserves device memory for at least the specified number of bytes.
		 * @details If @p n is greater than the current size, the existing memory is freed and 
		 *          new memory is allocated. If the allocation fails, the buffer is left empty.
		 * @param[in] n The number of bytes to reserve.
		 */
		void reserve(std::size_t n)
		{
			if(n > size()) {
				clear();
				m_data = device_alloc(n);
				m_size = (data() != nullptr) ? n : 0;
			}
		}
		/**
		 * @brief Frees the allocated device memory.
		 * @details Releases all device memory and resets the buffer to an empty state.
		 * @note This function is marked @c noexcept and is guaranteed not to throw exceptions.
		 */
		void clear() noexcept
		{
			device_free(data());
			defaults();
		}
		/**
		 * @brief Returns a pointer to the device memory buffer.
		 * @details Provides direct access to the underlying untyped device memory.
		 * @return A pointer to the device buffer.
		 */
		void* data() noexcept { return m_data; }

		/**
		 * @brief Returns a const pointer to the device memory buffer.
		 * @details Provides direct read-only access to the underlying untyped device memory.
		 * @return A const pointer to the device buffer.
		 */
		const void* data() const noexcept { return m_data; }

		/**
		 * @brief Returns the number of bytes in the buffer.
		 * @details Returns the capacity of the buffer (number of bytes allocated).
		 * @return The number of bytes the buffer can hold.
		 */
		std::size_t size() const noexcept { return m_size; }

	private:
		void *m_data;
		std::size_t m_size;
		
		void defaults()
		{
			m_data = nullptr;
			m_size = 0;
		}
};

/**
 * @ingroup culite_module_index_allocators
 * @brief A RAII wrapper for device memory allocation.
 * @details This class manages a dynamically allocated buffer on the device.
 *          It automatically handles memory allocation and deallocation following RAII principles,
 *          ensuring that device memory is properly freed when the buffer goes out of scope.
 * @tparam T The type of elements stored in the device buffer.
 */
template <typename T>
class DeviceBuffer {

	public:
		/**
		 * @brief Default constructor.
		 * @details Constructs an empty device buffer with no allocated memory.
		 */
		DeviceBuffer() {}

		/**
		 * @brief Constructs a device buffer with specified capacity.
		 * @details Constructs a device buffer and reserves memory for @p n elements on the device.
		 * @param[in] n The number of elements to allocate.
		 */
		explicit DeviceBuffer(std::size_t n) 
			: DeviceBuffer()
		{
			reserve(n);
		}

		/**
		 * @brief Destructor.
		 * @details Automatically frees the allocated device memory.
		 */
		~DeviceBuffer()
		{
			clear();
		}

		/**
		 * @brief Reserves device memory for at least the specified number of elements.
		 * @details If @p n is greater than the current size, the existing memory is freed and 
		 *          new memory is allocated. If the allocation fails, the buffer is left empty.
		 * @param[in] n The number of elements to reserve.
		 */
		void reserve(std::size_t n)
		{
			m_voidBuffer.reserve(n * sizeof(T));
		}

		/**
		 * @brief Frees the allocated device memory.
		 * @details Releases all device memory and resets the buffer to an empty state.
		 * @note This function is marked @c noexcept and is guaranteed not to throw exceptions.
		 */
		void clear() noexcept
		{
			m_voidBuffer.clear();
		}

		/**
		 * @brief Returns a pointer to the device memory buffer.
		 * @details Provides direct access to the underlying device memory.
		 * @return A pointer to the device buffer.
		 */
		T* data() noexcept { return static_cast<T*>(m_voidBuffer.data()); }

		/**
		 * @brief Returns a const pointer to the device memory buffer.
		 * @details Provides direct read-only access to the underlying device memory.
		 * @return A const pointer to the device buffer.
		 */
		const T* data() const noexcept { return static_cast<const T*>(m_voidBuffer.data()); }

		/**
		 * @brief Returns the number of elements in the buffer.
		 * @details Returns the capacity of the buffer (number of elements allocated).
		 * @return The number of elements the buffer can hold.
		 */
		std::size_t size() const noexcept { return m_voidBuffer.size(); }

	private:
		DeviceBufferVoid m_voidBuffer;
};

/*-------------------------------------------------*/
} // namespace culite
/*-------------------------------------------------*/

#endif // CULITE_IMALLOC_HPP_
