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
#include "culite/dense/dns_xxvector.hpp"

// system

// 3rd
#include <cla3p/types/property.hpp>
#include <cla3p/checks/block_ops_checks.hpp>

// culite
#include "culite/types/scalar.hpp"
#include "culite/bulk/dns1D.hpp"

/*-------------------------------------------------*/
namespace culite {
namespace dns {
/*-------------------------------------------------*/
template <typename T_Scalar>
XxVector<T_Scalar>::XxVector()
{
}
/*-------------------------------------------------*/
template <typename T_Scalar>
XxVector<T_Scalar>::XxVector(int_t n)
	: XiVector<T_Scalar>(n)
{
}
/*-------------------------------------------------*/
template <typename T_Scalar>
XxVector<T_Scalar>::XxVector(int_t n, T_Scalar *vals, bool bind)
	: XiVector<T_Scalar>(n, vals, bind)
{
}
/*-------------------------------------------------*/
template <typename T_Scalar>
XxVector<T_Scalar>::~XxVector()
{
}
/*-------------------------------------------------*/
template <typename T_Scalar>
XxVector<T_Scalar>::XxVector(const XiVector<T_Scalar>& other)
	: XiVector<T_Scalar>(other)
{
}
/*-------------------------------------------------*/
template <typename T_Scalar>
XxVector<T_Scalar>& XxVector<T_Scalar>::operator=(const XiVector<T_Scalar>& other)
{
	XiVector<T_Scalar>::operator=(other);
	return *this;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
XxVector<T_Scalar>::XxVector(XiVector<T_Scalar>&& other)
	: XiVector<T_Scalar>(std::move(other))
{
}
/*-------------------------------------------------*/
template <typename T_Scalar>
XxVector<T_Scalar>& XxVector<T_Scalar>::operator=(XiVector<T_Scalar>&& other)
{
	XiVector<T_Scalar>::operator=(std::move(other));
	return *this;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
alias::VirtualScal_vec<T_Scalar> XxVector<T_Scalar>::operator-() const
{
    return alias::VirtualScal_vec<T_Scalar>(virtualize(), makeScalar<T_Scalar>(-1));
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void XxVector<T_Scalar>::iscale(const T_Scalar& val)
{
	blk::dns::scale1D(this->size(), val, this->values());
}
/*-------------------------------------------------*/
template <typename T_Scalar>
VirtualRowvec<T_Scalar> XxVector<T_Scalar>::transpose() const
{
	return VirtualRowvec<T_Scalar>(this->size(), this->values(), 1, false);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
VirtualRowvec<T_Scalar> XxVector<T_Scalar>::ctranspose() const
{
	return VirtualRowvec<T_Scalar>(this->size(), this->values(), 1, true);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
alias::VirtualConj_vec<T_Scalar> XxVector<T_Scalar>::conjugate() const
{
    return alias::VirtualConj_vec<T_Scalar>(*this);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void XxVector<T_Scalar>::iconjugate()
{
	blk::dns::conjugate1D(this->size(), this->values());
}
/*-------------------------------------------------*/
template <typename T_Scalar>
typename XxVector<T_Scalar>::T_RScalar XxVector<T_Scalar>::normOne() const
{
	return blk::dns::normOne1D(this->size(), this->values());
}
/*-------------------------------------------------*/
template <typename T_Scalar>
typename XxVector<T_Scalar>::T_RScalar XxVector<T_Scalar>::normInf() const
{
	return blk::dns::normInf1D(this->size(), this->values());
}
/*-------------------------------------------------*/
template <typename T_Scalar>
typename XxVector<T_Scalar>::T_RScalar XxVector<T_Scalar>::normEuc() const
{ 
	return blk::dns::normEuc1D(this->size(), this->values());
}
/*-------------------------------------------------*/
template <typename T_Scalar>
XxVector<T_Scalar> XxVector<T_Scalar>::block(int_t ibgn, int_t ni) const
{
	XiVector<T_Scalar> tmp = rblock(ibgn,ni).get().copy();
	XxVector<T_Scalar> ret(tmp);
	return ret;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
XxVector<T_Scalar> XxVector<T_Scalar>::rblock(int_t ibgn, int_t ni)
{
	::cla3p::block_op_consistency_check(::cla3p::Property::General(), this->size(), 1, ibgn, 0, ni, 1);

	return XxVector<T_Scalar>(ni, this->values() + ibgn, false);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
::cla3p::Guard<XxVector<T_Scalar>> XxVector<T_Scalar>::rblock(int_t ibgn, int_t ni) const
{
	::cla3p::block_op_consistency_check(::cla3p::Property::General(), this->size(), 1, ibgn, 0, ni, 1);

	::cla3p::Guard<XxVector<T_Scalar>> ret = this->view(ni, this->values() + ibgn);
	return ret;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void XxVector<T_Scalar>::setBlock(int_t ibgn, const XxVector<T_Scalar>& src)
{
	rblock(ibgn, src.size()) = src;
}
/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/
template class XxVector<real_t>;
template class XxVector<real4_t>;
template class XxVector<complex_t>;
template class XxVector<complex8_t>;
/*-------------------------------------------------*/
} // namespace dns
} // namespace culite
/*-------------------------------------------------*/
