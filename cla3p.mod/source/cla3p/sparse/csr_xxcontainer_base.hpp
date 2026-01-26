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

#ifndef CLA3P_CSR_XXCONTAINER_BASE_HPP_
#define CLA3P_CSR_XXCONTAINER_BASE_HPP_

/**
 * @file
 */

 #include "cla3p/sparse/csx_xxcontainer_base.hpp"

/*-------------------------------------------------*/
namespace cla3p { 
namespace csr {
/*-------------------------------------------------*/

/**
 * @nosubgrouping 
 * @brief The sparse container base class (compressed sparse row).
 */
template <typename T_Int, typename T_Scalar>
class XxContainerBase : public csx::XxContainerBase<T_Int,T_Scalar> {

	public:
		XxContainerBase() {}

		explicit XxContainerBase(T_Int *rptr, T_Int *cidx, T_Scalar *vals)
        : csx::XxContainerBase<T_Int,T_Scalar>(rptr, cidx, vals) {}

		XxContainerBase(XxContainerBase<T_Int,T_Scalar>&) = delete;
		XxContainerBase<T_Int,T_Scalar>& operator=(XxContainerBase<T_Int,T_Scalar>&) = delete;

		XxContainerBase(XxContainerBase<T_Int,T_Scalar>&& other)
        : csx::XxContainerBase<T_Int,T_Scalar>(std::move(other)) {}
		XxContainerBase<T_Int,T_Scalar>& operator=(XxContainerBase<T_Int,T_Scalar>&& other)
		{
            csx::XxContainerBase<T_Int,T_Scalar>::operator=(std::move(other));
            return *this;
        }

		~XxContainerBase() { clear(); }

		/**
		 * @brief Access the row pointer array.
		 * @details Returns a pointer to the row pointer array (CSR format).
		 * @return Pointer to the row pointer array.
		 */
		T_Int* rowptr() { return this->xxxptr(); }

		/**
		 * @brief Access the row pointer array.
		 * @details Returns a pointer to the row pointer array (CSR format).
		 * @return Pointer to the row pointer array.
		 */
		const T_Int* rowptr() const { return this->xxxptr(); }

		/**
		 * @brief Access the column index array.
		 * @details Returns a pointer to the column index array (CSR format).
		 * @return Pointer to the column index array.
		 */
		T_Int* colidx() { return this->xxxidx(); }

		/**
		 * @brief Access the column index array.
		 * @details Returns a pointer to the column index array (CSR format).
		 * @return Pointer to the column index array.
		 */
		const T_Int* colidx() const { return this->xxxidx(); }

	protected:
		void clear()
		{
			csx::XxContainerBase<T_Int,T_Scalar>::clear();
		}
};

/*-------------------------------------------------*/
} // namespace csr
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_CSR_XXCONTAINER_BASE_HPP_
