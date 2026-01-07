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
#include "cla3p/dense/dns_xxvector.hpp"

// system
#include <algorithm>

// 3rd

// cla3p
#include "cla3p/perms.hpp"
#include "cla3p/bulk/dns.hpp"
#include "cla3p/checks/basic_checks.hpp"
#include "cla3p/checks/perm_checks.hpp"
#include "cla3p/checks/block_ops_checks.hpp"
#include "cla3p/algebra/operators_scale.hpp"

/*-------------------------------------------------*/
namespace cla3p {
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
VirtualScale<XxVector<T_Scalar>,VirtualObject<XxVector<T_Scalar>>> XxVector<T_Scalar>::operator-() const
{	
	return (T_Scalar(-1) * (*this));
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void XxVector<T_Scalar>::operator=(T_Scalar val)
{
	this->fill(val);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void XxVector<T_Scalar>::iscale(T_Scalar val)
{
	blk::dns::scale(uplo_t::Full, this->size(), 1, this->values(), this->size(), val);
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
VirtualConjugate<XxVector<T_Scalar>> XxVector<T_Scalar>::conjugate() const
{
	return VirtualConjugate<XxVector<T_Scalar>>(*this);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void XxVector<T_Scalar>::iconjugate()
{
	blk::dns::conjugate(uplo_t::Full, this->size(), 1, this->values(), this->size());
}
/*-------------------------------------------------*/
template <typename T_Scalar>
typename XxVector<T_Scalar>::T_RScalar XxVector<T_Scalar>::normOne() const
{
	return blk::dns::norm_one(
			prop_t::General,
			uplo_t::Full,
			this->size(),
			1,
			this->values(),
			this->size());
}
/*-------------------------------------------------*/
template <typename T_Scalar>
typename XxVector<T_Scalar>::T_RScalar XxVector<T_Scalar>::normInf() const
{
	return blk::dns::norm_inf(
			prop_t::General,
			uplo_t::Full,
			this->size(),
			1,
			this->values(),
			this->size());
}
/*-------------------------------------------------*/
template <typename T_Scalar>
typename XxVector<T_Scalar>::T_RScalar XxVector<T_Scalar>::normEuc() const
{ 
	return blk::dns::norm_euc(this->size(), this->values());
}
/*-------------------------------------------------*/
template <typename T_Scalar>
XxVector<T_Scalar> XxVector<T_Scalar>::permuteLeft(const prm::PiMatrix& P) const
{
	XxVector<T_Scalar> ret(this->size());
	permuteLeft(P, ret);
	return ret;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void XxVector<T_Scalar>::permuteLeft(const prm::PiMatrix& P, XxVector<T_Scalar>& dest) const
{
	perm_op_consistency_check(this->size(), 1, P.size(), 1);
	similarity_dim_check(this->size(), dest.size());
	blk::dns::permute(prop_t::General, uplo_t::Full, this->size(), 1, this->values(), this->size(), dest.values(), dest.size(), P.values(), nullptr);
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
	block_op_consistency_check(Property::General(), this->size(), 1, ibgn, 0, ni, 1);

	return XxVector<T_Scalar>(ni, this->values() + ibgn, false);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
Guard<XxVector<T_Scalar>> XxVector<T_Scalar>::rblock(int_t ibgn, int_t ni) const
{
	block_op_consistency_check(Property::General(), this->size(), 1, ibgn, 0, ni, 1);

	Guard<XxVector<T_Scalar>> ret = this->view(ni, this->values() + ibgn);
	return ret;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void XxVector<T_Scalar>::setBlock(int_t ibgn, const XxVector<T_Scalar>& src)
{
	rblock(ibgn, src.size()) = src;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
XxVector<T_Scalar> XxVector<T_Scalar>::random(int_t n, T_RScalar lo, T_RScalar hi)
{
	XxVector<T_Scalar> ret(n);
	blk::dns::rand(uplo_t::Full, ret.size(), 1, ret.values(), ret.size(), lo, hi);
	return ret;
}
/*-------------------------------------------------*/
template class XxVector<real_t>;
template class XxVector<real4_t>;
template class XxVector<complex_t>;
template class XxVector<complex8_t>;
/*-------------------------------------------------*/
} // namespace dns
} // namespace cla3p
/*-------------------------------------------------*/
