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
#include "cla3p/dense/dns_xivector.hpp"

// system
#include <algorithm>

// 3rd

// cla3p
#include "cla3p/bulk/dns_io.hpp"

#include "cla3p/error/exceptions.hpp"
#include "cla3p/error/literals.hpp"

#include "cla3p/support/utils.hpp"
#include "cla3p/checks/basic_checks.hpp"
#include "cla3p/checks/dns_checks.hpp"
#include "cla3p/checks/block_ops_checks.hpp"

/*-------------------------------------------------*/
namespace cla3p {
namespace dns {
/*-------------------------------------------------*/
template <typename T_Scalar>
XiVector<T_Scalar>::XiVector()
{
}
/*-------------------------------------------------*/
template <typename T_Scalar>
XiVector<T_Scalar>::XiVector(int_t n)
	: Meta1D<int_t>(n), XxContainer<T_Scalar>(n)
{
	if(n > 0) {
		checker();
	} else {
		clear();
	}
}
/*-------------------------------------------------*/
template <typename T_Scalar>
XiVector<T_Scalar>::XiVector(int_t n, T_Scalar *vals, bool bind)
	: Meta1D<int_t>(n), XxContainer<T_Scalar>(vals, bind)
{
	if(n > 0) {
		checker();
	} else {
		clear();
	}
}
/*-------------------------------------------------*/
template <typename T_Scalar>
XiVector<T_Scalar>::~XiVector()
{
	clear();
}
/*-------------------------------------------------*/
template <typename T_Scalar>
XiVector<T_Scalar>::XiVector(const XiVector<T_Scalar>& other)
	: XiVector(other.size())
{
	copyFromExisting(other);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
XiVector<T_Scalar>& XiVector<T_Scalar>::operator=(const XiVector<T_Scalar>& other)
{
	if(!(*this)) {
		*this = XiVector<T_Scalar>(other.size());
	}
	copyFromExisting(other);
	return *this;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
XiVector<T_Scalar>::XiVector(XiVector<T_Scalar>&& other)
{
	moveFrom(other);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
XiVector<T_Scalar>& XiVector<T_Scalar>::operator=(XiVector<T_Scalar>&& other)
{
	moveFrom(other);
	return *this;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void XiVector<T_Scalar>::clear()
{
	Meta1D<int_t>::clear();
	XxContainer<T_Scalar>::clear();
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void XiVector<T_Scalar>::fill(T_Scalar val)
{
	std::fill(this->values(), this->values() + size(), val);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
T_Scalar& XiVector<T_Scalar>::operator()(int_t i)
{
	if(i >= size()) {
		throw err::OutOfBounds(msg::IndexOutOfBounds(size(),i));
	} // out-of-bounds

	return (this->values())[i];
}
/*-------------------------------------------------*/
template <typename T_Scalar>
const T_Scalar& XiVector<T_Scalar>::operator()(int_t i) const
{
	if(i >= size()) {
		throw err::OutOfBounds(msg::IndexOutOfBounds(size(),i));
	} // out-of-bounds

	return (this->values())[i];
}
/*-------------------------------------------------*/
template <typename T_Scalar>
XiVector<T_Scalar> XiVector<T_Scalar>::copy() const
{
	XiVector<T_Scalar> ret(size());
	ret.copyFromExisting(*this);
	return ret;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
XiVector<T_Scalar> XiVector<T_Scalar>::rcopy()
{
	return XiVector<T_Scalar>(size(), this->values(), false);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
Guard<XiVector<T_Scalar>> XiVector<T_Scalar>::rcopy() const
{
	return view(size(), this->values());
}
/*-------------------------------------------------*/
template <typename T_Scalar>
XiVector<T_Scalar> XiVector<T_Scalar>::move()
{
	XiVector<T_Scalar> ret;
	ret.moveFrom(*this);
	return ret;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
std::string XiVector<T_Scalar>::info(const std::string& header) const
{ 
	std::string top;
	std::string bottom;
	fill_info_margins(header, top, bottom);

	std::ostringstream ss;

	ss << top << "\n";

	ss << "  Datatype............. " << TypeTraits<T_Scalar>::type_name() << "\n";
	ss << "  Precision............ " << TypeTraits<T_Scalar>::prec_name() << "\n";
	ss << "  Size................. " << size() << "\n";
	ss << "  Values............... " << this->values() << "\n";
	ss << "  Owner................ " << boolToYesNo(this->owner()) << "\n";

	ss << bottom << "\n";

	return ss.str();
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void XiVector<T_Scalar>::print(std::streamsize prec) const
{
	toStream(std::cout, prec);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void XiVector<T_Scalar>::toStream(std::ostream& os, std::streamsize prec) const
{
	blk::dns::print_to_stream(os, uplo_t::Full, size(), 1, this->values(), size(), prec);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void XiVector<T_Scalar>::checker() const
{
	dns_consistency_check(size(), 1, this->values(), size());
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void XiVector<T_Scalar>::moveFrom(XiVector<T_Scalar>& other)
{
	if(this != &other) {

		if(*this) {
			*this = other;
		} else {
			Meta1D<int_t>::operator=(std::move(other));
			XxContainer<T_Scalar>::operator=(std::move(other));
			other.unbind();
		} // similar

		other.clear();

	} // do not apply on self
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void XiVector<T_Scalar>::copyFromExisting(const XiVector<T_Scalar>& other)
{
	if(this != &other) {
		similarity_dim_check(size(), other.size());
		std::copy(other.values(), other.values() + size(), this->values());
	} // do not apply on self
}
/*-------------------------------------------------*/
template <typename T_Scalar>
Guard<XiVector<T_Scalar>> XiVector<T_Scalar>::view(int_t n, const T_Scalar *vals)
{
	XiVector<T_Scalar> tmp(n, const_cast<T_Scalar*>(vals), false);
	Guard<XiVector<T_Scalar>> ret(tmp);
	return ret;
}
/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/
template class XiVector<int_t>;
template class XiVector<real_t>;
template class XiVector<real4_t>;
template class XiVector<complex_t>;
template class XiVector<complex8_t>;
/*-------------------------------------------------*/
} // namespace dns
} // namespace cla3p
/*-------------------------------------------------*/
