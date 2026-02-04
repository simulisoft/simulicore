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
#include "culite/sparse/csr_xxmatrix.hpp"

// system

// 3rd
#include <cla3p/support/utils.hpp>
#include <cla3p/checks/basic_checks.hpp>
#include <cla3p/checks/hermitian_coeff_checks.hpp>
#include <cla3p/checks/csx_checks.hpp>

// culite
#include "culite/types/integer.hpp"
#include "culite/types/scalar.hpp"
#include "culite/error/exceptions.hpp"
#include "culite/bulk/dns1D.hpp"

/*-------------------------------------------------*/
namespace culite {
namespace csr {
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
XxMatrix<T_Int,T_Scalar>::XxMatrix()
{
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
XxMatrix<T_Int,T_Scalar>::XxMatrix(T_Int nr, T_Int nc, T_Int nz, const ::cla3p::Property& pr)
	: ::cla3p::MatrixMeta<T_Int>(nr, nc, ::cla3p::sanitizeProperty<T_Cla3pScalar>(pr)), 
      csx::XxContainer<T_Int,T_Scalar>(nr, nz)
{
	if(nr > 0 && nc > 0) {
		checker();
	} else {
		clear();
	}
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
XxMatrix<T_Int,T_Scalar>::XxMatrix(T_Int nr, T_Int nc, T_Int *rptr, T_Int *cidx, T_Scalar *vals, bool bind, const ::cla3p::Property& pr)
	: ::cla3p::MatrixMeta<T_Int>(nr, nc, ::cla3p::sanitizeProperty<T_Cla3pScalar>(pr)), 
      csx::XxContainer<T_Int,T_Scalar>(rptr, cidx, vals, bind)
{
	if(nr > 0 && nc > 0) {
		checker();
	} else {
		clear();
	}
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
XxMatrix<T_Int,T_Scalar>::~XxMatrix()
{
	clear();
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
XxMatrix<T_Int,T_Scalar>::XxMatrix(const XxMatrix<T_Int,T_Scalar>& other)
	: XxMatrix(other.nrows(), other.ncols(), other.nnz(), other.prop())
{
	copyFromExisting(other);
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
XxMatrix<T_Int,T_Scalar>& XxMatrix<T_Int,T_Scalar>::operator=(const XxMatrix<T_Int,T_Scalar>& other)
{
	if(!(*this)) {
		*this = XxMatrix<T_Int,T_Scalar>(other.nrows(), other.ncols(), other.nnz(), other.prop());
	}
	return copyFromExisting(other);
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
XxMatrix<T_Int,T_Scalar>::XxMatrix(XxMatrix<T_Int,T_Scalar>&& other)
{
	moveFrom(other);
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
XxMatrix<T_Int,T_Scalar>& XxMatrix<T_Int,T_Scalar>::operator=(XxMatrix<T_Int,T_Scalar>&& other)
{
	return moveFrom(other);
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
void XxMatrix<T_Int,T_Scalar>::clear()
{
	::cla3p::MatrixMeta<T_Int>::clear();
	csx::XxContainer<T_Int,T_Scalar>::clear();
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
XxMatrix<T_Int,T_Scalar> XxMatrix<T_Int,T_Scalar>::operator-() const
{
    T_Scalar coeff = makeScalar<T_Scalar>(-1);
    XxMatrix<T_Int,T_Scalar> ret = *this;
    ret.iscale(coeff);
    return ret;
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
T_Int* XxMatrix<T_Int,T_Scalar>::rowptr() { return this->xxxptr(); }
template <typename T_Int, typename T_Scalar>
T_Int* XxMatrix<T_Int,T_Scalar>::colidx() { return this->xxxidx(); }
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
const T_Int* XxMatrix<T_Int,T_Scalar>::rowptr() const { return this->xxxptr(); }
template <typename T_Int, typename T_Scalar>
const T_Int* XxMatrix<T_Int,T_Scalar>::colidx() const { return this->xxxidx(); }
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
T_Int XxMatrix<T_Int,T_Scalar>::nnz() const
{
    T_Int nz = 0;
	if(!this->empty()) {
        // TODO: rethink about copying
        memCopyD2H(1, this->rowptr() + this->nrows(), &nz);
	}
	return nz;
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
std::string XxMatrix<T_Int,T_Scalar>::info(const std::string& header) const
{ 
	std::string top;
	std::string bottom;
	::cla3p::fill_info_margins(header, top, bottom);

	std::ostringstream ss;

	ss << top << "\n";

	ss << "  Datatype............. " << TypeTraits<T_Scalar>::type_name() << "\n";
	ss << "  Precision............ " << TypeTraits<T_Scalar>::prec_name() << "\n";
	ss << "  Number of rows....... " << this->nrows() << "\n";
	ss << "  Number of columns.... " << this->ncols() << "\n";
	ss << "  Number of non zeros.. " << nnz() << "\n";
	ss << "  Rowptr............... " << this->rowptr() << "\n";
	ss << "  Colidx............... " << this->colidx() << "\n";
	ss << "  Values............... " << this->values() << "\n";
	ss << "  Property............. " << this->prop() << "\n";
	ss << "  Owner................ " << ::cla3p::boolToYesNo(this->owner()) << "\n";

	ss << bottom << "\n";

	return ss.str();
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
XxMatrix<T_Int,T_Scalar>& XxMatrix<T_Int,T_Scalar>::copyFromExisting(const XxMatrix<T_Int,T_Scalar>& other)
{
	if(this != &other) {

		::cla3p::similarity_check(*this, other);
		::cla3p::similarity_dim_check(nnz(), other.nnz());

		T_Int nr = other.nrows() + 1;
		T_Int nz = other.nnz();

        memCopyD2D(nr, other.rowptr(), this->rowptr());
        memCopyD2D(nz, other.colidx(), this->colidx());
        memCopyD2D(nz, other.values(), this->values());

	} // do not apply on self

    return *this;
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
XxMatrix<T_Int,T_Scalar>& XxMatrix<T_Int,T_Scalar>::moveFrom(XxMatrix<T_Int,T_Scalar>& other)
{
	if(this != &other) {

		if(*this) {
			*this = other;
		} else {
			::cla3p::MatrixMeta<T_Int>::operator=(std::move(other));
			csx::XxContainer<T_Int,T_Scalar>::operator=(std::move(other));
			other.unbind();
		} // similar

		other.clear();

	} // do not apply on self

    return *this;
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
XxMatrix<T_Int,T_Scalar> XxMatrix<T_Int,T_Scalar>::copy() const
{
	XxMatrix<T_Int,T_Scalar> ret(this->nrows(), this->ncols(), nnz(), this->prop());
	ret.copyFromExisting(*this);
	return ret;
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
XxMatrix<T_Int,T_Scalar> XxMatrix<T_Int,T_Scalar>::rcopy()
{
	return XxMatrix<T_Int,T_Scalar>(this->nrows(), this->ncols(), this->rowptr(), this->colidx(), this->values(), false, this->prop());
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
::cla3p::Guard<XxMatrix<T_Int,T_Scalar>> XxMatrix<T_Int,T_Scalar>::rcopy() const
{
	return view(this->nrows(), this->ncols(), this->rowptr(), this->colidx(), this->values(), this->prop());
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
XxMatrix<T_Int,T_Scalar> XxMatrix<T_Int,T_Scalar>::move()
{
	XxMatrix<T_Int,T_Scalar> ret;
	ret.moveFrom(*this);
	return ret;
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
void XxMatrix<T_Int,T_Scalar>::iscale(T_Scalar val)
{
    T_Cla3pScalar cla3pVal = TypeTraits<T_Scalar>::toCla3pType(val);
    ::cla3p::hermitian_coeff_check<T_Cla3pScalar>(this->prop(), cla3pVal);
	blk::dns::scale1D(nnz(), val, this->values());
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
XxMatrix<T_Int,T_Scalar> XxMatrix<T_Int,T_Scalar>::conjugate() const
{
    XxMatrix<T_Int,T_Scalar> ret = *this;
    ret.iconjugate();
    return ret;
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
void XxMatrix<T_Int,T_Scalar>::iconjugate()
{
    blk::dns::conjugate1D(nnz(), this->values());
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
void XxMatrix<T_Int,T_Scalar>::checker() const
{
    if(!this->prop().isGeneral()) {
		throw err::CudaException("Only the 'General' matrix property is supported currently.");
	}

	::cla3p::csx_consistency_check(this->nrows(), this->ncols(), nnz(), 
                                   this->rowptr(), 
                                   this->colidx(), 
                                   this->values(), 
                                   this->prop());
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
::cla3p::Guard<XxMatrix<T_Int,T_Scalar>> XxMatrix<T_Int,T_Scalar>::view(T_Int nr, T_Int nc, 
                                                                        const T_Int *cptr, 
                                                                        const T_Int *ridx, 
                                                                        const T_Scalar *vals, 
                                                                        const ::cla3p::Property& pr)
{
	XxMatrix<T_Int,T_Scalar> tmp(nr, nc,
			const_cast<T_Int   *>(cptr),
			const_cast<T_Int   *>(ridx),
			const_cast<T_Scalar*>(vals), false, pr);
	::cla3p::Guard<XxMatrix<T_Int,T_Scalar>> ret(tmp);
	return ret;
}
/*-------------------------------------------------*/
#if 0
template <typename T_Int, typename T_Scalar>
void XxMatrix<T_Int,T_Scalar>::copyToHost(::cla3p::csr::XxMatrix<T_Cla3pInt,T_Cla3pScalar>& dest) const
{
	if(!dest) {
		dest = ::cla3p::csr::XxMatrix<T_Cla3pInt,T_Cla3pScalar>(this->nrows(), this->ncols(), nnz(), this->prop());
	}

	::cla3p::similarity_check<T_Int>(this->prop(), this->nrows(), this->ncols(), 
                                     dest.prop(), dest.nrows(), dest.ncols());
    ::cla3p::similarity_dim_check<T_Int>(nnz(), dest.nnz());

    T_Int nr = nrows() + 1;
    T_Int nz = nnz();

    memCopyD2H(nr, this->rowptr(), dest.rowptr());
    memCopyD2H(nz, this->colidx(), dest.colidx());
    memCopyD2H(nz, this->values(), dest.values());
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
void XxMatrix<T_Int,T_Scalar>::copyFromHost(const ::cla3p::csr::XxMatrix<T_Cla3pInt,T_Cla3pScalar>& src)
{
	if(!(*this)) {
		*this = XxMatrix<T_Int,T_Scalar>(src.nrows(), src.ncols(), src.nnz(), src.prop());
	}

	::cla3p::similarity_check<T_Int>(this->prop(), this->nrows(), this->ncols(), 
                                     src.prop(), src.nrows(), src.ncols());
    ::cla3p::similarity_dim_check<T_Int>(nnz(), src.nnz());

    T_Int nr = src.nrows() + 1;
    T_Int nz = src.nnz();

    memCopyH2D(nr, src.rowptr(), this->rowptr());
    memCopyH2D(nz, src.colidx(), this->colidx());
    memCopyH2D(nz, src.values(), this->values());
}
#endif // 0
/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/
template class XxMatrix<int_t,real_t>;
template class XxMatrix<int_t,real4_t>;
template class XxMatrix<int_t,complex_t>;
template class XxMatrix<int_t,complex8_t>;
/*-------------------------------------------------*/
} // namespace csr
} // namespace culite
/*-------------------------------------------------*/
