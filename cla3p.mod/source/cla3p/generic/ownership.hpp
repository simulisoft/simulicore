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

#ifndef CLA3P_OWNERSHIP_HPP_
#define CLA3P_OWNERSHIP_HPP_

/**
 * @file
 */

/*-------------------------------------------------*/
namespace cla3p { 
/*-------------------------------------------------*/

/**
 * @nosubgrouping
 * @brief Memory ownership metadata.
 * @details Tracks whether an object is responsible for deallocating its
 *          associated memory. When the ownership flag is true, the object
 *          will deallocate its contents upon destruction.
 */
class Ownership {

	public:
		/**
		 * @brief Default constructor.
		 * @details Initializes ownership state to default.
		 */
		Ownership();
		
		/**
		 * @brief Ownership constructor.
		 * @details Initializes ownership state with the specified flag.
		 * @param[in] owner If @c true, the object takes ownership of its contents.
		 */
		Ownership(bool owner);
		
		/**
		 * @brief Destructor.
		 * @details Destroys the ownership metadata object.
		 */
		~Ownership();

		/**
		 * @brief Tests memory ownership.
		 * @details Returns @c true if the object is responsible for deallocating
		 *          its contents. When true, memory is automatically freed upon
		 *          object destruction.
		 * @return @c true if the object owns its contents, @c false otherwise.
		 */
		bool owner() const;

		/**
		 * @brief Relinquishes memory ownership.
		 * @details Makes the object no longer responsible for content deallocation.
		 *          Use with caution as this may lead to memory leaks if the
		 *          contents are not deallocated manually. Auto-allocated data
		 *          should be freed using @ref i_free().
		 */
		void unbind();

	protected:
		void clear();

	private:
		bool m_owner;

		void setOwner(bool owner);

		void defaults();
};

/*-------------------------------------------------*/
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_OWNERSHIP_HPP_
