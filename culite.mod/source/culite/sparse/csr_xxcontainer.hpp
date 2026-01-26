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

#ifndef CULITE_CSR_XXCONTAINER_HPP_
#define CULITE_CSR_XXCONTAINER_HPP_

/**
 * @file
 */

#include <cstddef>

#include <cla3p/sparse/csr_xxcontainer_base.hpp>
#include <cla3p/generic/ownership.hpp>

#include "culite/support/imalloc.hpp"
#include "culite/support/utils.hpp"

/*-------------------------------------------------*/
namespace culite { 
namespace csr {
/*-------------------------------------------------*/

/**
 * @nosubgrouping 
 * @brief The sparse container class (compressed sparse row).
 */
template <typename T_Int, typename T_Scalar>
class XxContainer : public ::cla3p::csr::XxContainerBase<T_Int, T_Scalar>, public ::cla3p::Ownership {

	public:
		XxContainer() {}

		explicit XxContainer(std::size_t nr, std::size_t nz)
			: ::cla3p::csr::XxContainerBase<T_Int,T_Scalar>(nr ? device_alloc_t<T_Int>(nr+1) : nullptr,
			                                                nr ? device_alloc_t<T_Int>(nz) : nullptr,
			                                                nr ? device_alloc_t<T_Scalar>(nz) : nullptr),
			  Ownership(nr ? true : false)
		{
            // TODO: rethink about copying
            T_Int intNz = static_cast<T_Int>(nz);
            memCopyH2D(1, &intNz, this->rowptr() + nr);
		}

		explicit XxContainer(T_Int *rptr, T_Int *cidx, T_Scalar *vals, bool bind)
			: ::cla3p::csr::XxContainerBase<T_Int,T_Scalar>(rptr ? rptr : nullptr, 
											                rptr ? cidx : nullptr, 
											                rptr ? vals : nullptr),
			  Ownership(rptr ? bind : false) {}

		XxContainer(XxContainer<T_Int,T_Scalar>&) = delete;
		XxContainer<T_Int,T_Scalar>& operator=(XxContainer<T_Int,T_Scalar>&) = delete;

		XxContainer(XxContainer<T_Int,T_Scalar>&& other) { moveFrom(other); }
		XxContainer<T_Int,T_Scalar>& operator=(XxContainer<T_Int,T_Scalar>&& other) { return moveFrom(other); }

		~XxContainer() { clear(); }

	protected:
		void clear()
		{
			if(owner()) {
				device_free(this->rowptr());
				device_free(this->colidx());
				device_free(this->values());
			} // owner
			::cla3p::csr::XxContainerBase<T_Int,T_Scalar>::clear();
			Ownership::clear();
		}

	private:
		XxContainer<T_Int,T_Scalar>& moveFrom(XxContainer<T_Int,T_Scalar>& other)
		{
			if(this != &other) {
				clear();
				::cla3p::csr::XxContainerBase<T_Int,T_Scalar>::operator=(std::move(other));
				Ownership::operator=(std::move(other));
				other.unbind();
				other.clear();
			} // do not apply on self
            return *this;
		}
};

/*-------------------------------------------------*/
} // namespace csr
} // namespace culite
/*-------------------------------------------------*/

#endif // CULITE_CSR_XXCONTAINER_HPP_