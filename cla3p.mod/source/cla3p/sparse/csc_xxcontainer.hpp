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

#include "cla3p/sparse/csc_xxcontainer_base.hpp"
#include "cla3p/generic/ownership.hpp"
#include "cla3p/support/imalloc.hpp"

/*-------------------------------------------------*/
namespace cla3p { 
namespace csc {
/*-------------------------------------------------*/

/**
 * @nosubgrouping 
 * @brief The sparse container class (compressed sparse column).
 */
template <typename T_Int, typename T_Scalar>
class XxContainer : public XxContainerBase<T_Int, T_Scalar>, public Ownership {

	public:
		XxContainer() {}

		explicit XxContainer(std::size_t nc, std::size_t nz)
			: XxContainerBase<T_Int,T_Scalar>(nc ? i_malloc<T_Int>(nc+1) : nullptr,
			                                  nz ? i_malloc<T_Int>(nz) : nullptr,
			                                  nz ? i_malloc<T_Scalar>(nz) : nullptr),
			  Ownership(nc ? true : false)
		{
			this->colptr()[nc] = static_cast<T_Int>(nz);
		}

		explicit XxContainer(T_Int *cptr, T_Int *ridx, T_Scalar *vals, bool bind)
			: XxContainerBase<T_Int,T_Scalar>(cptr ? cptr : nullptr, 
											  cptr ? ridx : nullptr, 
											  cptr ? vals : nullptr),
			  Ownership(cptr ? bind : false) {}

		XxContainer(XxContainer<T_Int,T_Scalar>&) = delete;
		XxContainer<T_Int,T_Scalar>& operator=(XxContainer<T_Int,T_Scalar>&) = delete;

		XxContainer(XxContainer<T_Int,T_Scalar>&& other)
		{
			moveFrom(other);
		}

		XxContainer<T_Int,T_Scalar>& operator=(XxContainer<T_Int,T_Scalar>&& other)
		{
			moveFrom(other);
			return *this;
		}

		~XxContainer()
		{
			clear();
		}

	protected:
		void clear()
		{
			if(owner()) {
				i_free(this->colptr());
				i_free(this->rowidx());
				i_free(this->values());
			} // owner
			XxContainerBase<T_Int,T_Scalar>::clear();
			Ownership::clear();
		}

	private:
		void moveFrom(XxContainer<T_Int,T_Scalar>& other)
		{
			if(this != &other) {
				clear();
				XxContainerBase<T_Int,T_Scalar>::operator=(std::move(other));
				Ownership::operator=(std::move(other));
				other.unbind();
				other.clear();
			} // do not apply on self
		}
};

/*-------------------------------------------------*/
} // namespace csc
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_CSC_XXCONTAINER_HPP_
