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
XxVector<T_Scalar> XxVector<T_Scalar>::operator-() const
{	
	XxVector<T_Scalar> ret = *this;
	ret.iscale(makeScalar<T_Scalar>(-1));
	return ret;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void XxVector<T_Scalar>::iscale(T_Scalar /*val*/)
{
	// TODO: implement
	// blk::dns::scale(uplo_t::Full, this->size(), 1, this->values(), this->size(), val);
}
/*-------------------------------------------------*/
#if 0
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
#endif // 0
/*-------------------------------------------------*/
template <typename T_Scalar>
void XxVector<T_Scalar>::iconjugate()
{
	// TODO: implement
	// blk::dns::conjugate(uplo_t::Full, this->size(), 1, this->values(), this->size());
}
/*-------------------------------------------------*/
template <typename T_Scalar>
typename XxVector<T_Scalar>::T_RScalar XxVector<T_Scalar>::normOne() const
{
	// TODO: implement
	// return blk::dns::norm_one(
	// 		prop_t::General,
	// 		uplo_t::Full,
	// 		this->size(),
	// 		1,
	// 		this->values(),
	// 		this->size());
	return 0;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
typename XxVector<T_Scalar>::T_RScalar XxVector<T_Scalar>::normInf() const
{
	// TODO: implement
	// return blk::dns::norm_inf(
	// 		prop_t::General,
	// 		uplo_t::Full,
	// 		this->size(),
	// 		1,
	// 		this->values(),
	// 		this->size());
	return 0;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
typename XxVector<T_Scalar>::T_RScalar XxVector<T_Scalar>::normEuc() const
{ 
	// TODO: implement
	// return blk::dns::norm_euc(this->size(), this->values());
	return 0;
}
/*-------------------------------------------------*/
#if 0
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
#endif // 0
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
Guard<XxVector<T_Scalar>> XxVector<T_Scalar>::rblock(int_t ibgn, int_t ni) const
{
	::cla3p::block_op_consistency_check(::cla3p::Property::General(), this->size(), 1, ibgn, 0, ni, 1);

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
XxVector<T_Scalar> XxVector<T_Scalar>::random(int_t n, T_RScalar /*lo*/, T_RScalar /*hi*/)
{
	XxVector<T_Scalar> ret(n);
	// TODO: implement
	// blk::dns::rand(uplo_t::Full, ret.size(), 1, ret.values(), ret.size(), lo, hi);
	return ret;
}
/*-------------------------------------------------*/
template class XxVector<real_t>;
template class XxVector<real4_t>;
template class XxVector<complex_t>;
template class XxVector<complex8_t>;
/*-------------------------------------------------*/
} // namespace dns
} // namespace culite
/*-------------------------------------------------*/
