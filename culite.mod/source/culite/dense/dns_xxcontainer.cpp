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

// this file inc
#include "culite/dense/dns_xxcontainer.hpp"

// system

// 3rd

// cla3p

// culite
#include "culite/types/integer.hpp"
#include "culite/types/scalar.hpp"
#include "culite/support/imalloc.hpp"

/*-------------------------------------------------*/
namespace culite {
namespace dns {
/*-------------------------------------------------*/
template <typename T_Scalar>
XxContainer<T_Scalar>::XxContainer()
{
  defaults();
}
/*-------------------------------------------------*/
template <typename T_Scalar>
XxContainer<T_Scalar>::XxContainer(std::size_t numElements)
{
	if(numElements) {
		T_Scalar *vals = device_alloc_t<T_Scalar>(numElements);
		::cla3p::Ownership::operator=(::cla3p::Ownership(true));
		setValues(vals);
	}
}
/*-------------------------------------------------*/
template <typename T_Scalar>
XxContainer<T_Scalar>::XxContainer(T_Scalar *vals, bool bind)
{
	if(vals) {
		Ownership::operator=(Ownership(bind));
		setValues(vals);
	}
}
/*-------------------------------------------------*/
template <typename T_Scalar>
XxContainer<T_Scalar>::XxContainer(XxContainer<T_Scalar>&& other)
{
	moveFrom(other);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
XxContainer<T_Scalar>& XxContainer<T_Scalar>::operator=(XxContainer<T_Scalar>&& other)
{
	moveFrom(other);
	return *this;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
XxContainer<T_Scalar>::~XxContainer()
{
  clear();
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void XxContainer<T_Scalar>::defaults()
{
	setValues(nullptr);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void XxContainer<T_Scalar>::setValues(T_Scalar *vals)
{
	m_values = vals;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
T_Scalar* XxContainer<T_Scalar>::values()
{
	return m_values;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
const T_Scalar* XxContainer<T_Scalar>::values() const
{
	return m_values;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void XxContainer<T_Scalar>::clear()
{
	if(owner()) {
		device_free(static_cast<void*>(values()));
	} // owner

	::cla3p::Ownership::clear();

	defaults();
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void XxContainer<T_Scalar>::moveFrom(XxContainer<T_Scalar>& other)
{
	if(this != &other) {
		clear();
		::cla3p::Ownership::operator=(std::move(other));
		setValues(other.values());
		other.unbind();
		other.clear();
	} // do not apply on self
}
/*-------------------------------------------------*/
template class XxContainer<int_t>;
template class XxContainer<real_t>;
template class XxContainer<real4_t>;
template class XxContainer<complex_t>;
template class XxContainer<complex8_t>;
/*-------------------------------------------------*/
} // namespace dns
} // namespace culite
/*-------------------------------------------------*/
