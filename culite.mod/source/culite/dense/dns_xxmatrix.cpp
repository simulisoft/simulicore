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
#include "culite/dense/dns_xxmatrix.hpp"

// system

// 3rd
#include <cla3p/checks/basic_checks.hpp>
#include <cla3p/checks/hermitian_coeff_checks.hpp>
#include <cla3p/checks/block_ops_checks.hpp>
#include <cla3p/checks/dns_checks.hpp>
#include <cla3p/support/utils.hpp>

// culite
#include "culite/types/scalar.hpp"

/*-------------------------------------------------*/
namespace culite {
namespace dns {
/*-------------------------------------------------*/
template <typename T_Scalar>
XxMatrix<T_Scalar>::XxMatrix()
{
	defaults();
}
/*-------------------------------------------------*/
template <typename T_Scalar>
XxMatrix<T_Scalar>::XxMatrix(int_t nr, int_t nc, const ::cla3p::Property& pr)
	: MatrixMeta<int_t>(nr, nc, ::cla3p::sanitizeProperty<T_ScalarHost>(pr)), XxContainer<T_Scalar>(nr * nc)
{
	if(nr > 0 && nc > 0) {
		setLd(nr);
		checker();
	} else {
		clear();
	}
}
/*-------------------------------------------------*/
template <typename T_Scalar>
XxMatrix<T_Scalar>::XxMatrix(int_t nr, int_t nc, T_Scalar *vals, int_t ldv, bool bind, const ::cla3p::Property& pr)
	: MatrixMeta<int_t>(nr, nc, ::cla3p::sanitizeProperty<T_ScalarHost>(pr)), XxContainer<T_Scalar>(vals, bind)
{
	if(nr > 0 && nc > 0) {
		setLd(ldv);
		checker();
	} else {
		clear();
	}
}
/*-------------------------------------------------*/
template <typename T_Scalar>
XxMatrix<T_Scalar>::~XxMatrix()
{
	clear();
}
/*-------------------------------------------------*/
template <typename T_Scalar>
XxMatrix<T_Scalar>::XxMatrix(const XxMatrix<T_Scalar>& other)
	: XxMatrix(other.nrows(), other.ncols(), other.prop())
{
	copyFromExisting(other);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
XxMatrix<T_Scalar>& XxMatrix<T_Scalar>::operator=(const XxMatrix<T_Scalar>& other)
{
	if(!(*this)) {
		*this = XxMatrix<T_Scalar>(other.nrows(), other.ncols(), other.prop());
	}
	copyFromExisting(other);
	return *this;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
XxMatrix<T_Scalar>::XxMatrix(XxMatrix<T_Scalar>&& other)
{
	moveFrom(other);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
XxMatrix<T_Scalar>& XxMatrix<T_Scalar>::operator=(XxMatrix<T_Scalar>&& other)
{
	moveFrom(other);
	return *this;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void XxMatrix<T_Scalar>::defaults()
{
	setLd(0);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void XxMatrix<T_Scalar>::setLd(int_t ld)
{
	m_ld = ld;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
int_t XxMatrix<T_Scalar>::ld() const
{
	return m_ld;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void XxMatrix<T_Scalar>::clear()
{
	::cla3p::MatrixMeta<int_t>::clear();
	XxContainer<T_Scalar>::clear();
	defaults();
}
/*-------------------------------------------------*/
template <typename T_Scalar>
XxMatrix<T_Scalar> XxMatrix<T_Scalar>::operator-() const
{
	XxMatrix<T_Scalar> ret = *this;
	ret.iscale(makeScalar<T_Scalar>(-1));
	return ret;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
XxMatrix<T_Scalar> XxMatrix<T_Scalar>::copy() const
{
	XxMatrix<T_Scalar> ret(nrows(), ncols(), prop());
	ret.copyFromExisting(*this);
	return ret;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
XxMatrix<T_Scalar> XxMatrix<T_Scalar>::rcopy()
{
	return XxMatrix<T_Scalar>(nrows(), ncols(), this->values(), ld(), false, prop());
}
/*-------------------------------------------------*/
template <typename T_Scalar>
::cla3p::Guard<XxMatrix<T_Scalar>> XxMatrix<T_Scalar>::rcopy() const
{
	return view(nrows(), ncols(), this->values(), ld(), prop());
}
/*-------------------------------------------------*/
template <typename T_Scalar>
XxMatrix<T_Scalar> XxMatrix<T_Scalar>::move()
{
	XxMatrix<T_Scalar> ret;
	ret.moveFrom(*this);
	return ret;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void XxMatrix<T_Scalar>::iscale(T_Scalar val)
{
	T_ScalarHost valHost = TypeTraits<T_Scalar>::toHostType(val);
	::cla3p::hermitian_coeff_check<T_ScalarHost>(prop(), valHost);
	// TODO: implement
	// blk::dns::scale(
	//		prop().uplo(),
	//		nrows(),
	//		ncols(),
	//		this->values(),
	//		ld(), val);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
std::string XxMatrix<T_Scalar>::info(const std::string& header) const
{ 
	std::string top;
	std::string bottom;
	::cla3p::fill_info_margins(header, top, bottom);

	std::ostringstream ss;

	ss << top << "\n";

	ss << "  Datatype............. " << TypeTraits<T_Scalar>::type_name() << "\n";
	ss << "  Precision............ " << TypeTraits<T_Scalar>::prec_name() << "\n";
	ss << "  Number of rows....... " << nrows() << "\n";
	ss << "  Number of columns.... " << ncols() << "\n";
	ss << "  Leading dimension.... " << ld() << "\n";
	ss << "  Values............... " << this->values() << "\n";
	ss << "  Property............. " << prop() << "\n";
	ss << "  Owner................ " << ::cla3p::boolToYesNo(this->owner()) << "\n";

	ss << bottom << "\n";

	return ss.str();
}
/*-------------------------------------------------*/
#if 0
template <typename T_Scalar>
VirtualTranspose<XxMatrix<T_Scalar>> XxMatrix<T_Scalar>::transpose() const
{
	return VirtualTranspose<XxMatrix<T_Scalar>>(*this, false);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
VirtualTranspose<XxMatrix<T_Scalar>> XxMatrix<T_Scalar>::ctranspose() const
{
	return VirtualTranspose<XxMatrix<T_Scalar>>(*this, true);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
VirtualConjugate<XxMatrix<T_Scalar>> XxMatrix<T_Scalar>::conjugate() const
{
	return VirtualConjugate<XxMatrix<T_Scalar>>(*this);
}
#endif // 0
/*-------------------------------------------------*/
template <typename T_Scalar>
void XxMatrix<T_Scalar>::iconjugate()
{
	// TODO: implement
	// blk::dns::conjugate(prop().uplo(), nrows(), ncols(), this->values(), ld());
}
/*-------------------------------------------------*/
template <typename T_Scalar>
typename XxMatrix<T_Scalar>::T_RScalar XxMatrix<T_Scalar>::normOne() const
{
	// TODO: implement
	//return blk::dns::norm_one(
	//		prop().type(),
	//		prop().uplo(),
	//		nrows(),
	//		ncols(),
	//		this->values(),
	//		ld());
	return 0;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
typename XxMatrix<T_Scalar>::T_RScalar XxMatrix<T_Scalar>::normInf() const
{
	// TODO: implement
	//return blk::dns::norm_inf(
	//		prop().type(),
	//		prop().uplo(),
	//		nrows(),
	//		ncols(),
	//		this->values(),
	//		ld());
	return 0;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
typename XxMatrix<T_Scalar>::T_RScalar XxMatrix<T_Scalar>::normMax() const
{
	// TODO: implement
	//return blk::dns::norm_max(
	//		prop().type(),
	//		prop().uplo(),
	//		nrows(),
	//		ncols(),
	//		this->values(),
	//		ld());
	return 0;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
typename XxMatrix<T_Scalar>::T_RScalar XxMatrix<T_Scalar>::normFro() const
{
	// TODO: implement
	//return blk::dns::norm_fro(
	//		prop().type(),
	//		prop().uplo(),
	//		nrows(),
	//		ncols(),
	//		this->values(),
	//		ld());
	return 0;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
XxMatrix<T_Scalar> XxMatrix<T_Scalar>::general() const
{
	XxMatrix<T_Scalar> ret = *this;
	ret.igeneral();
	return ret;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void XxMatrix<T_Scalar>::igeneral()
{
	// TODO: implement
	//if(prop().isGeneral()) {
//
	//	return;
//
	//} else if(prop().isSymmetric()) {
//
	//	blk::dns::sy2ge(prop().uplo(), ncols(), this->values(), ld());
//
	//} else if(prop().isHermitian()) {
//
	//	blk::dns::he2ge(prop().uplo(), ncols(), this->values(), ld());
//
	//} else if(prop().isTriangular()) {
//
	//	blk::dns::tr2ge(prop().uplo(), nrows(), ncols(), this->values(), ld());
//
	//} else if(prop().isSkew()) {
//
	//	blk::dns::sk2ge(prop().uplo(), ncols(), this->values(), ld());
//
	//} else {
//
	//	throw err::Exception();
//
	//} // property 
//
	//setProp(Property::General());
}
/*-------------------------------------------------*/
#if 0
template <typename T_Scalar>
XxMatrix<T_Scalar> XxMatrix<T_Scalar>::permuteLeftRight(const prm::PiMatrix& P, const prm::PiMatrix& Q) const
{
	XxMatrix<T_Scalar> ret(nrows(), ncols(), prop());
	permuteLeftRight(P, Q, ret);
	return ret;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void XxMatrix<T_Scalar>::permuteLeftRight(const prm::PiMatrix& P, const prm::PiMatrix& Q, XxMatrix<T_Scalar>& dest) const
{
	perm_ge_op_consistency_check(prop().type(), nrows(), ncols(), P.size(), Q.size());
	similarity_check(prop(), nrows(), ncols(), dest.prop(), dest.nrows(), dest.ncols());
	blk::dns::permute(prop().type(), prop().uplo(), nrows(), ncols(), this->values(), ld(), dest.values(), dest.ld(), P.values(), Q.values());
}
/*-------------------------------------------------*/
template <typename T_Scalar>
XxMatrix<T_Scalar> XxMatrix<T_Scalar>::permuteLeft(const prm::PiMatrix& P) const
{
	XxMatrix<T_Scalar> ret(nrows(), ncols(), prop());
	permuteLeft(P, ret);
	return ret;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void XxMatrix<T_Scalar>::permuteLeft(const prm::PiMatrix& P, XxMatrix<T_Scalar>& dest) const
{
	perm_ge_op_consistency_check(prop().type(), nrows(), ncols(), P.size(), ncols());
	similarity_check(prop(), nrows(), ncols(), dest.prop(), dest.nrows(), dest.ncols());
	blk::dns::permute(prop().type(), prop().uplo(), nrows(), ncols(), this->values(), ld(), dest.values(), dest.ld(), P.values(), nullptr);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
XxMatrix<T_Scalar> XxMatrix<T_Scalar>::permuteRight(const prm::PiMatrix& Q) const
{
	XxMatrix<T_Scalar> ret(nrows(), ncols(), prop());
	permuteRight(Q, ret);
	return ret;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void XxMatrix<T_Scalar>::permuteRight(const prm::PiMatrix& Q, XxMatrix<T_Scalar>& dest) const
{
	perm_ge_op_consistency_check(prop().type(), nrows(), ncols(), nrows(), Q.size());
	similarity_check(prop(), nrows(), ncols(), dest.prop(), dest.nrows(), dest.ncols());
	blk::dns::permute(prop().type(), prop().uplo(), nrows(), ncols(), this->values(), ld(), dest.values(), dest.ld(), nullptr, Q.values());
}
/*-------------------------------------------------*/
template <typename T_Scalar>
XxMatrix<T_Scalar> XxMatrix<T_Scalar>::permuteMirror(const prm::PiMatrix& P) const
{
	XxMatrix<T_Scalar> ret(nrows(), ncols(), prop());
	permuteMirror(P, ret);
	return ret;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void XxMatrix<T_Scalar>::permuteMirror(const prm::PiMatrix& P, XxMatrix<T_Scalar>& dest) const
{
	perm_op_consistency_check(nrows(), ncols(), P.size(), P.size());
	similarity_check(prop(), nrows(), ncols(), dest.prop(), dest.nrows(), dest.ncols());

	prm::PiMatrix iP;
	if(prop().isGeneral())
		iP = P.inverse();

	blk::dns::permute(prop().type(), prop().uplo(), nrows(), ncols(), this->values(), ld(), dest.values(), dest.ld(), P.values(), iP.values());
}
#endif // 0
/*-------------------------------------------------*/
template <typename T_Scalar>
XxMatrix<T_Scalar> XxMatrix<T_Scalar>::block(int_t ibgn, int_t jbgn, int_t ni, int_t nj) const
{
	return rblock(ibgn,jbgn,ni,nj).get().copy();
}
/*-------------------------------------------------*/
template <typename T_Scalar>
XxMatrix<T_Scalar> XxMatrix<T_Scalar>::rblock(int_t ibgn, int_t jbgn, int_t ni, int_t nj)
{
	::cla3p::Property pr = ::cla3p::block_op_consistency_check(
			prop(), nrows(), ncols(),
			ibgn, jbgn, ni, nj);
	
	// TODO: implement bulks
	T_Scalar *p_vij = nullptr; // blk::dns::ptrmv(ld(), this->values(), ibgn, jbgn);
	
	return XxMatrix<T_Scalar>(ni, nj, p_vij, ld(), false, pr);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
::cla3p::Guard<XxMatrix<T_Scalar>> XxMatrix<T_Scalar>::rblock(int_t ibgn, int_t jbgn, int_t ni, int_t nj) const
{
	::cla3p::Property pr = ::cla3p::block_op_consistency_check(
			prop(), nrows(), ncols(),
			ibgn, jbgn, ni, nj);
	
	// TODO: implement bulks
	const T_Scalar *p_vij = nullptr; //blk::dns::ptrmv(ld(), this->values(), ibgn, jbgn);
	
	return view(ni, nj, p_vij, ld(), pr);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void XxMatrix<T_Scalar>::setBlock(int_t ibgn, int_t jbgn, const XxMatrix<T_Scalar>& src)
{
	rblock(ibgn, jbgn, src.nrows(), src.ncols()) = src;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
XxVector<T_Scalar> XxMatrix<T_Scalar>::column(int_t j) const
{
	XiVector<T_Scalar> tmp = rcolumn(j).get().copy();
	XxVector<T_Scalar> ret(tmp);
	return ret;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
XxVector<T_Scalar> XxMatrix<T_Scalar>::rcolumn(int_t j)
{
	XxMatrix<T_Scalar> tmp = rblock(0, j, nrows(), 1);
	return XxVector<T_Scalar>(tmp.nrows(), tmp.values(), false);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
::cla3p::Guard<XxVector<T_Scalar>> XxMatrix<T_Scalar>::rcolumn(int_t j) const 
{ 
	::cla3p::Guard<XxMatrix<T_Scalar>> tmpMat = rblock(0, j, nrows(), 1);
	::cla3p::Guard<XxVector<T_Scalar>> ret = XiVector<T_Scalar>::view(tmpMat.get().nrows(), tmpMat.get().values());
	return ret;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
XxMatrix<T_Scalar> XxMatrix<T_Scalar>::row(int_t i) const
{ 
	return rrow(i).get().copy();
}
/*-------------------------------------------------*/
template <typename T_Scalar>
XxMatrix<T_Scalar> XxMatrix<T_Scalar>::rrow(int_t i)
{
	return rblock(i, 0, 1, ncols());
}
/*-------------------------------------------------*/
template <typename T_Scalar>
::cla3p::Guard<XxMatrix<T_Scalar>> XxMatrix<T_Scalar>::rrow(int_t i) const 
{ 
	return rblock(i, 0, 1, ncols());
}
/*-------------------------------------------------*/
// TODO: use virtuals
//template <typename T_Scalar>
//VirtualRowvec<T_Scalar> XxMatrix<T_Scalar>::rrowvec(int_t i) const
//{
//	::cla3p::Guard<XxMatrix<T_Scalar>> tmp = rrow(i);
//	return VirtualRowvec<T_Scalar>(tmp.get().ncols(), tmp.get().values(), tmp.get().ld(), false);
//}
/*-------------------------------------------------*/
template <typename T_Scalar>
XxMatrix<T_Scalar> XxMatrix<T_Scalar>::random(int_t nr, int_t nc, const ::cla3p::Property& pr, T_RScalar /*lo*/, T_RScalar /*hi*/)
{
	XxMatrix<T_Scalar> ret(nr, nc, pr);
	// TODO: implement
	//blk::dns::rand(ret.prop().uplo(), ret.nrows(), ret.ncols(), ret.values(), ret.ld(), lo, hi);
	//blk::dns::set_diag_zeros(ret.prop().type(), ret.ncols(), ret.values(), ret.ld());
	return ret;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
::cla3p::Guard<XxMatrix<T_Scalar>> XxMatrix<T_Scalar>::view(int_t nr, int_t nc, const T_Scalar *vals, int_t ldv, const ::cla3p::Property& pr)
{
	XxMatrix<T_Scalar> tmp(nr, nc, const_cast<T_Scalar*>(vals), ldv, false, pr);
	::cla3p::Guard<XxMatrix<T_Scalar>> ret(tmp);
	return ret;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void XxMatrix<T_Scalar>::checker() const
{
	::cla3p::dns_consistency_check(prop(), nrows(), ncols(), this->values(), ld());
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void XxMatrix<T_Scalar>::moveFrom(XxMatrix<T_Scalar>& other)
{
	if(this != &other) {

		if(*this) {
			*this = other;
		} else {
			::cla3p::MatrixMeta<int_t>::operator=(std::move(other));
			XxContainer<T_Scalar>::operator=(std::move(other));
			setLd(other.ld());
			other.unbind();
		} // similar

		other.clear();

	} // do not apply on self
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void XxMatrix<T_Scalar>::copyFromExisting(const XxMatrix<T_Scalar>& other)
{
	if(this != &other) {
		::cla3p::similarity_check(prop(), nrows(), ncols(), other.prop(), other.nrows(), other.ncols());
		// TODO: implement bulks
		//blk::dns::copy(other.prop().uplo(), other.nrows(), other.ncols(), other.values(), other.ld(), this->values(), ld());
	} // do not apply on self
}
/*-------------------------------------------------*/
template class XxMatrix<real_t>;
template class XxMatrix<real4_t>;
template class XxMatrix<complex_t>;
template class XxMatrix<complex8_t>;
/*-------------------------------------------------*/
} // namespace dns
} // namespace culite
/*-------------------------------------------------*/
