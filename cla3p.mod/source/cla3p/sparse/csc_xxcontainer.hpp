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

#ifndef CLA3P_CSC_XXCONTAINER_HPP_
#define CLA3P_CSC_XXCONTAINER_HPP_

/**
 * @file
 */

#include <cstddef>

#include "cla3p/types.hpp"
#include "cla3p/generic/ownership.hpp"
#include "cla3p/generic/guard.hpp"

/*-------------------------------------------------*/
namespace cla3p { 
namespace csc {
/*-------------------------------------------------*/

/**
 * @nosubgrouping 
 * @brief The sparse container class (compressed sparse column).
 */
template <typename T_Int, typename T_Scalar>
class XxContainer : public Ownership {

	public:
		using index_type = T_Int;
		using value_type = T_Scalar;

	public:
		XxContainer();
		explicit XxContainer(std::size_t nc, std::size_t nz);
		explicit XxContainer(T_Int *cptr, T_Int *ridx, T_Scalar *vals, bool bind);
		~XxContainer();

		/**
		 * @copydoc standard_docs::colptr()
		 */
		T_Int* colptr();

		/**
		 * @copydoc standard_docs::colptr()
		 */
		const T_Int* colptr() const;

		/**
		 * @copydoc standard_docs::rowidx()
		 */
		T_Int* rowidx();

		/**
		 * @copydoc standard_docs::rowidx()
		 */
		const T_Int* rowidx() const;

		/**
		 * @copydoc standard_docs::values()
		 */
		T_Scalar* values();

		/**
		 * @copydoc standard_docs::values()
		 */
		const T_Scalar* values() const;

	protected:
		void clear();

	private:

		T_Int*    m_colptr;
		T_Int*    m_rowidx;
		T_Scalar* m_values;

		void setColptr(T_Int*);
		void setRowidx(T_Int*);
		void setValues(T_Scalar*);

		void defaults();
};

/*-------------------------------------------------*/
} // namespace csc
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_CSC_XXCONTAINER_HPP_
