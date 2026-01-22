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
#include "culite/bulk/dns2D.hpp"
#include "culite/support/utils.hpp"

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
	: MatrixMeta<int_t>(nr, nc, ::cla3p::sanitizeProperty<T_Cla3pScalar>(pr)), XxContainer<T_Scalar>(nr * nc)
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
	: MatrixMeta<int_t>(nr, nc, ::cla3p::sanitizeProperty<T_Cla3pScalar>(pr)), XxContainer<T_Scalar>(vals, bind)
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
	T_Cla3pScalar valHost = TypeTraits<T_Scalar>::toHostType(val);
	::cla3p::hermitian_coeff_check<T_Cla3pScalar>(prop(), valHost);
	blk::dns::scale2D(
		prop().uplo(),
		nrows(),
		ncols(),
		this->values(),
		ld(), val);
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
template <typename T_Scalar>
XxMatrix<T_Scalar> XxMatrix<T_Scalar>::transpose() const
{
	XxMatrix<T_Scalar> ret(ncols(), nrows(), prop().transpose());
	blk::dns::transpose2D(
		nrows(),
		ncols(),
		this->values(),
		ld(),
		ret.values(),
		ret.ld());
	return ret;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
XxMatrix<T_Scalar> XxMatrix<T_Scalar>::ctranspose() const
{
	XxMatrix<T_Scalar> ret(ncols(), nrows(), prop().transpose());
	blk::dns::ctranspose2D(
		nrows(),
		ncols(),
		this->values(),
		ld(),
		ret.values(),
		ret.ld());
	return ret;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
XxMatrix<T_Scalar> XxMatrix<T_Scalar>::conjugate() const
{
	XxMatrix<T_Scalar> ret(nrows(), ncols(), prop().transpose());
	ret.iconjugate();
	return ret;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void XxMatrix<T_Scalar>::iconjugate()
{
	blk::dns::conjugate2D(
		prop().uplo(),
		nrows(),
		ncols(),
		this->values(),
		ld());
}
/*-------------------------------------------------*/
template <typename T_Scalar>
typename XxMatrix<T_Scalar>::T_RScalar XxMatrix<T_Scalar>::normOne() const
{
	return blk::dns::normOne2D(
		prop().type(),
		prop().uplo(),
		nrows(),
		ncols(),
		this->values(),
		ld());
}
/*-------------------------------------------------*/
template <typename T_Scalar>
typename XxMatrix<T_Scalar>::T_RScalar XxMatrix<T_Scalar>::normInf() const
{
	return blk::dns::normInf2D(
		prop().type(),
		prop().uplo(),
		nrows(),
		ncols(),
		this->values(),
		ld());
}
/*-------------------------------------------------*/
template <typename T_Scalar>
typename XxMatrix<T_Scalar>::T_RScalar XxMatrix<T_Scalar>::normMax() const
{
	return blk::dns::normMax2D(
		prop().type(),
		prop().uplo(),
		nrows(),
		ncols(),
		this->values(),
		ld());
}
/*-------------------------------------------------*/
template <typename T_Scalar>
typename XxMatrix<T_Scalar>::T_RScalar XxMatrix<T_Scalar>::normFro() const
{
	return blk::dns::normFro2D(
		prop().type(),
		prop().uplo(),
		nrows(),
		ncols(),
		this->values(),
		ld());
}
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
	
	T_Scalar *p_vij = blk::dns::ptrmv(ld(), this->values(), ibgn, jbgn);
	
	return XxMatrix<T_Scalar>(ni, nj, p_vij, ld(), false, pr);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
::cla3p::Guard<XxMatrix<T_Scalar>> XxMatrix<T_Scalar>::rblock(int_t ibgn, int_t jbgn, int_t ni, int_t nj) const
{
	::cla3p::Property pr = ::cla3p::block_op_consistency_check(
			prop(), nrows(), ncols(),
			ibgn, jbgn, ni, nj);
	
	const T_Scalar *p_vij = blk::dns::ptrmv(ld(), this->values(), ibgn, jbgn);
	
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
	if(!prop().isGeneral()) {
		throw err::CudaException("Only the 'General' matrix property is supported currently.");
	}
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

		::cla3p::similarity_check(prop(), 
		                          nrows(), 
								  ncols(), 
								  other.prop(), 
								  static_cast<int_t>(other.nrows()), 
								  static_cast<int_t>(other.ncols()));

		blk::dns::copy2D(other.prop().uplo(), 
				         other.nrows(), 
						 other.ncols(), 
						 other.values(), 
						 other.ld(), 
						 this->values(), 
						 ld());

	} // do not apply on self
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void XxMatrix<T_Scalar>::copyToHost(::cla3p::dns::XxMatrix<T_Cla3pScalar>& dest) const
{
	if(!dest) {
		dest = ::cla3p::dns::XxMatrix<T_Cla3pScalar>(nrows(), ncols(), prop());
	}
	::cla3p::similarity_check(prop(), 
	                          nrows(), 
							  ncols(), 
							  dest.prop(), 
							  static_cast<int_t>(dest.nrows()), 
							  static_cast<int_t>(dest.ncols()));

    memCopyD2H(nrows(), ncols(),
               this->values(), ld(),
               dest.values(), dest.ld());
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void XxMatrix<T_Scalar>::copyFromHost(const ::cla3p::dns::XxMatrix<T_Cla3pScalar>& src)
{
	if(!(*this)) {
		*this = XxMatrix<T_Scalar>(src.nrows(), src.ncols(), src.prop());
	}
	::cla3p::similarity_check(prop(), 
	                          nrows(), 
							  ncols(), 
							  src.prop(), 
							  static_cast<int_t>(src.nrows()), 
							  static_cast<int_t>(src.ncols()));

    memCopyH2D<T_Scalar>(src.nrows(), src.ncols(), 
	                     src.values(), src.ld(), 
					     this->values(), ld());

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
