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

#ifndef CLA3P_CSX_XXCONTAINER_HPP_
#define CLA3P_CSX_XXCONTAINER_HPP_

/**
 * @file
 */

#include <cstddef>

#include "cla3p/sparse/csx_xxcontainer_base.hpp"
#include "cla3p/generic/ownership.hpp"
#include "cla3p/support/imalloc.hpp"

/*-------------------------------------------------*/
namespace cla3p { 
namespace csx {
/*-------------------------------------------------*/

/**
 * @nosubgrouping 
 * @brief The sparse container class (compressed sparse generic).
 */
template <typename T_Int, typename T_Scalar>
class XxContainer : public XxContainerBase<T_Int, T_Scalar>, public Ownership {

    public:
        XxContainer() {}

        explicit XxContainer(std::size_t np, std::size_t nz)
            : XxContainerBase<T_Int,T_Scalar>(np ? i_malloc_t<T_Int>(np+1) : nullptr,
                                              np ? i_malloc_t<T_Int>(nz) : nullptr,
                                              np ? i_malloc_t<T_Scalar>(nz) : nullptr),
              Ownership(np ? true : false)
        {
            this->xxxptr()[np] = static_cast<T_Int>(nz);
        }

        explicit XxContainer(T_Int *xptr, T_Int *xidx, T_Scalar *vals, bool bind)
            : XxContainerBase<T_Int,T_Scalar>(xptr ? xptr : nullptr, 
                                              xptr ? xidx : nullptr, 
                                              xptr ? vals : nullptr),
              Ownership(xptr ? bind : false) {}

        XxContainer(XxContainer<T_Int,T_Scalar>&) = delete;
        XxContainer<T_Int,T_Scalar>& operator=(XxContainer<T_Int,T_Scalar>&) = delete;

        XxContainer(XxContainer<T_Int,T_Scalar>&& other) { moveFrom(other); }
        XxContainer<T_Int,T_Scalar>& operator=(XxContainer<T_Int,T_Scalar>&& other) { return moveFrom(other); }

        ~XxContainer() { clear(); }

    protected:
        void clear()
        {
            if(owner()) {
                i_free(this->xxxptr());
                i_free(this->xxxidx());
                i_free(this->values());
            } // owner
            XxContainerBase<T_Int,T_Scalar>::clear();
            Ownership::clear();
        }

    private:
        XxContainer<T_Int,T_Scalar>& moveFrom(XxContainer<T_Int,T_Scalar>& other)
        {
            if(this != &other) {
                clear();
                XxContainerBase<T_Int,T_Scalar>::operator=(std::move(other));
                Ownership::operator=(std::move(other));
                other.unbind();
                other.clear();
            } // do not apply on self
            return *this;
        }
};

/*-------------------------------------------------*/
} // namespace csx
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_CSX_XXCONTAINER_HPP_
