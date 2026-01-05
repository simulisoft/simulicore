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
#include "cla3p/sparse/csc_xxmatrix.hpp"

// system

// 3rd

// cla3p
#include "cla3p/dense.hpp"
#include "cla3p/sparse.hpp"
#include "cla3p/perms.hpp"
#include "cla3p/bulk/dns.hpp"
#include "cla3p/bulk/csc.hpp"
#include "cla3p/bulk/csc_math.hpp"
#include "cla3p/support/imalloc.hpp"
#include "cla3p/support/rand.hpp"
#include "cla3p/support/utils.hpp"

#include "cla3p/checks/basic_checks.hpp"
#include "cla3p/checks/csc_checks.hpp"
#include "cla3p/checks/block_ops_checks.hpp"
#include "cla3p/checks/transp_checks.hpp"
#include "cla3p/checks/perm_checks.hpp"
#include "cla3p/checks/hermitian_coeff_checks.hpp"
#include "cla3p/algebra/operators_scale.hpp"

/*-------------------------------------------------*/
namespace cla3p {
namespace csc {
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
XxMatrix<T_Int,T_Scalar>::XxMatrix()
{
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
XxMatrix<T_Int,T_Scalar>::XxMatrix(int_t nr, int_t nc, int_t nz, const Property& pr)
	: MatrixMeta<T_Int>(nr, nc, sanitizeProperty<T_Scalar>(pr)), XxContainer<T_Int,T_Scalar>(nc, nz)
{
	if(nr > 0 && nc > 0) {
		checker();
	} else {
		clear();
	}
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
XxMatrix<T_Int,T_Scalar>::XxMatrix(int_t nr, int_t nc, T_Int *cptr, T_Int *ridx, T_Scalar *vals, bool bind, const Property& pr)
	: MatrixMeta<T_Int>(nr, nc, sanitizeProperty<T_Scalar>(pr)), XxContainer<T_Int,T_Scalar>(cptr, ridx, vals, bind)
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
	copyFromExisting(other);
	return *this;
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
	moveFrom(other);
	return (*this);
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
void XxMatrix<T_Int,T_Scalar>::clear()
{
	MatrixMeta<T_Int>::clear();
	XxContainer<T_Int,T_Scalar>::clear();
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
VirtualScale<XxMatrix<T_Int,T_Scalar>,VirtualObject<XxMatrix<T_Int,T_Scalar>>> XxMatrix<T_Int,T_Scalar>::operator-() const
{
	return (T_Scalar(-1) * (*this));
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
int_t XxMatrix<T_Int,T_Scalar>::nnz() const
{
	if(!this->empty()) {
		return (this->colptr())[this->ncols()];
	}
	return 0;
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
std::string XxMatrix<T_Int,T_Scalar>::info(const std::string& header) const
{ 
	std::string top;
	std::string bottom;
	fill_info_margins(header, top, bottom);

	std::ostringstream ss;

	ss << top << "\n";

	ss << "  Datatype............. " << TypeTraits<T_Scalar>::type_name() << "\n";
	ss << "  Precision............ " << TypeTraits<T_Scalar>::prec_name() << "\n";
	ss << "  Number of rows....... " << this->nrows() << "\n";
	ss << "  Number of columns.... " << this->ncols() << "\n";
	ss << "  Number of non zeros.. " << nnz() << "\n";
	ss << "  Colptr............... " << this->colptr() << "\n";
	ss << "  Rowidx............... " << this->rowidx() << "\n";
	ss << "  Values............... " << this->values() << "\n";
	ss << "  Property............. " << this->prop() << "\n";
	ss << "  Owner................ " << boolToYesNo(this->owner()) << "\n";

	ss << bottom << "\n";

	return ss.str();
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
void XxMatrix<T_Int,T_Scalar>::copyFromExisting(const XxMatrix<T_Int,T_Scalar>& other)
{
	if(this != &other) {

		similarity_check(this->prop(), this->nrows(), this->ncols(), other.prop(), other.nrows(), other.ncols());
		similarity_dim_check(nnz(), other.nnz());

		int_t nc = other.ncols() + 1;
		int_t nz = other.nnz();

		// 
		// TODO: perhaps use a copy for 1D arrays
		//
		blk::dns::copy(uplo_t::Full, nc, 1, other.colptr(), nc, this->colptr(), nc);
		blk::dns::copy(uplo_t::Full, nz, 1, other.rowidx(), nz, this->rowidx(), nz);
		blk::dns::copy(uplo_t::Full, nz, 1, other.values(), nz, this->values(), nz);

	} // do not apply on self
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
void XxMatrix<T_Int,T_Scalar>::moveFrom(XxMatrix<T_Int,T_Scalar>& other)
{
	if(this != &other) {

		if(*this) {
			*this = other;
		} else {
			MatrixMeta<T_Int>::operator=(std::move(other));
			XxContainer<T_Int,T_Scalar>::operator=(std::move(other));
			other.unbind();
		} // similar

		other.clear();

	} // do not apply on self
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
void XxMatrix<T_Int,T_Scalar>::toStream(std::ostream& os, std::streamsize prec) const
{
	blk::csc::print_to_stream(os, this->ncols(), this->colptr(), this->rowidx(), this->values(), prec);
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
	return XxMatrix<T_Int,T_Scalar>(this->nrows(), this->ncols(), this->colptr(), this->rowidx(), this->values(), false, this->prop());
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
Guard<XxMatrix<T_Int,T_Scalar>> XxMatrix<T_Int,T_Scalar>::rcopy() const
{
	return view(this->nrows(), this->ncols(), this->colptr(), this->rowidx(), this->values(), this->prop());
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
	hermitian_coeff_check(this->prop(), val);
	blk::dns::scale(uplo_t::Full, nnz(), 1, this->values(), nnz(), val);
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
VirtualTranspose<XxMatrix<T_Int,T_Scalar>> XxMatrix<T_Int,T_Scalar>::transpose() const
{
	return VirtualTranspose<XxMatrix<T_Int,T_Scalar>>(*this, false);
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
VirtualTranspose<XxMatrix<T_Int,T_Scalar>> XxMatrix<T_Int,T_Scalar>::ctranspose() const
{
	return VirtualTranspose<XxMatrix<T_Int,T_Scalar>>(*this, true);
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
VirtualConjugate<XxMatrix<T_Int,T_Scalar>> XxMatrix<T_Int,T_Scalar>::conjugate() const
{
	return VirtualConjugate<XxMatrix<T_Int,T_Scalar>>(*this);
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
void XxMatrix<T_Int,T_Scalar>::iconjugate()
{
	//
	// TODO: perhaps use a conjugate for 1D arrays
	//
	blk::dns::conjugate(uplo_t::Full, nnz(), 1, this->values(), nnz());
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
typename XxMatrix<T_Int,T_Scalar>::T_RScalar XxMatrix<T_Int,T_Scalar>::normOne() const
{
	return blk::csc::norm_one(
			this->prop().type(),
			this->ncols(), 
			this->colptr(), 
			this->rowidx(), 
			this->values());
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
typename XxMatrix<T_Int,T_Scalar>::T_RScalar XxMatrix<T_Int,T_Scalar>::normInf() const
{
	return blk::csc::norm_inf(
			this->prop().type(),
			this->nrows(),
			this->ncols(),
			this->colptr(), 
			this->rowidx(), 
			this->values());
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
typename XxMatrix<T_Int,T_Scalar>::T_RScalar XxMatrix<T_Int,T_Scalar>::normMax() const
{
	return blk::csc::norm_max(
			this->ncols(),
			this->colptr(), 
			this->values());
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
typename XxMatrix<T_Int,T_Scalar>::T_RScalar XxMatrix<T_Int,T_Scalar>::normFro() const
{
	return blk::csc::norm_fro(
			this->prop().type(),
			this->ncols(),
			this->colptr(), 
			this->rowidx(), 
			this->values());
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
XxMatrix<T_Int,T_Scalar> XxMatrix<T_Int,T_Scalar>::general() const
{
	XxMatrix<T_Int,T_Scalar> ret;

	T_Int    *colptr_ge = nullptr;
	T_Int    *rowidx_ge = nullptr;
	T_Scalar *values_ge = nullptr;

	if(this->prop().isGeneral()) {

		ret = copy();

	} else if(this->prop().isSymmetric()) {

		colptr_ge = i_malloc<T_Int>(this->ncols() + 1);
		blk::csc::uplo2ge_colptr(this->prop().uplo(), this->ncols(), this->colptr(), this->rowidx(), colptr_ge);
		T_Int nz = colptr_ge[this->ncols()];
		rowidx_ge = i_malloc<T_Int>(nz);
		values_ge = i_malloc<T_Scalar>(nz);
		blk::csc::sy2ge(this->prop().uplo(), this->ncols(), this->colptr(), this->rowidx(), this->values(), colptr_ge, rowidx_ge, values_ge);
		ret = XxMatrix<T_Int,T_Scalar>(this->nrows(), this->ncols(), colptr_ge, rowidx_ge, values_ge, true);

	} else if(this->prop().isHermitian()) {

		colptr_ge = i_malloc<T_Int>(this->ncols() + 1);
		blk::csc::uplo2ge_colptr(this->prop().uplo(), this->ncols(), this->colptr(), this->rowidx(), colptr_ge);
		T_Int nz = colptr_ge[this->ncols()];
		rowidx_ge = i_malloc<T_Int>(nz);
		values_ge = i_malloc<T_Scalar>(nz);
		blk::csc::he2ge(this->prop().uplo(), this->ncols(), this->colptr(), this->rowidx(), this->values(), colptr_ge, rowidx_ge, values_ge);
		ret = XxMatrix<T_Int,T_Scalar>(this->nrows(), this->ncols(), colptr_ge, rowidx_ge, values_ge, true);

	} else if(this->prop().isTriangular()) {

		ret = copy();
		ret.setProp(Property::General());

	} else {

		throw err::Exception();

	} // property 

	return ret;
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
dns::XxMatrix<T_Scalar> XxMatrix<T_Int,T_Scalar>::toDns() const
{
	dns::XxMatrix<T_Scalar> ret(this->nrows(), this->ncols(), this->prop());
	ret = 0;
	for(int_t j = 0; j < this->ncols(); j++) {
		for(T_Int irow = (this->colptr())[j]; irow < (this->colptr())[j+1]; irow++) {
			ret((this->rowidx())[irow],j) = (this->values())[irow];
		} // irow
	} // j

	return ret;
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
XxMatrix<T_Int,T_Scalar> XxMatrix<T_Int,T_Scalar>::permuteLeftRight(const prm::PiMatrix& P, const prm::PiMatrix& Q) const
{
	perm_ge_op_consistency_check(this->prop().type(), this->nrows(), this->ncols(), P.size(), Q.size());

	XxMatrix<T_Int,T_Scalar> ret(this->nrows(), this->ncols(), nnz(), this->prop());
	blk::csc::permute(this->prop().type(), this->prop().uplo(), this->nrows(), this->ncols(), 
			this->colptr(), this->rowidx(), this->values(), 
			ret.colptr(), ret.rowidx(), ret.values(), 
			P.values(), Q.values());

	return ret;
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
XxMatrix<T_Int,T_Scalar> XxMatrix<T_Int,T_Scalar>::permuteLeft(const prm::PiMatrix& P) const
{
	perm_ge_op_consistency_check(this->prop().type(), this->nrows(), this->ncols(), P.size(), this->ncols());

	XxMatrix<T_Int,T_Scalar> ret(this->nrows(), this->ncols(), nnz(), this->prop());
	blk::csc::permute(this->prop().type(), this->prop().uplo(), this->nrows(), this->ncols(), 
			this->colptr(), this->rowidx(), this->values(), 
			ret.colptr(), ret.rowidx(), ret.values(), 
			P.values(), nullptr);

	return ret;
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
XxMatrix<T_Int,T_Scalar> XxMatrix<T_Int,T_Scalar>::permuteRight(const prm::PiMatrix& Q) const
{
	perm_ge_op_consistency_check(this->prop().type(), this->nrows(), this->ncols(), this->nrows(), Q.size());

	XxMatrix<T_Int,T_Scalar> ret(this->nrows(), this->ncols(), nnz(), this->prop());
	blk::csc::permute(this->prop().type(), this->prop().uplo(), this->nrows(), this->ncols(), 
			this->colptr(), this->rowidx(), this->values(), 
			ret.colptr(), ret.rowidx(), ret.values(), 
			nullptr, Q.values());

	return ret;
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
XxMatrix<T_Int,T_Scalar> XxMatrix<T_Int,T_Scalar>::permuteMirror(const prm::PiMatrix& P) const
{
	perm_op_consistency_check(this->nrows(), this->ncols(), P.size(), P.size());

	XxMatrix<T_Int,T_Scalar> ret(this->nrows(), this->ncols(), nnz(), this->prop());
	blk::csc::permute(this->prop().type(), this->prop().uplo(), this->nrows(), this->ncols(), 
			this->colptr(), this->rowidx(), this->values(), 
			ret.colptr(), ret.rowidx(), ret.values(), 
			P.values(), nullptr);

	return ret;
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
XxMatrix<T_Int,T_Scalar> XxMatrix<T_Int,T_Scalar>::block(int_t ibgn, int_t jbgn, int_t ni, int_t nj) const
{
	Property pr = block_op_consistency_check(this->prop(), this->nrows(), this->ncols(), ibgn, jbgn, ni, nj);

	if(!ni || !nj) return XxMatrix<T_Int,T_Scalar>();

	T_Int *cptr = i_calloc<T_Int>(nj + 1);

	int_t iend = ibgn + ni;
	int_t jend = jbgn + nj;

	for(int_t j = jbgn; j < jend; j++) {
		int_t jlocal = j - jbgn;
		for(T_Int irow = (this->colptr())[j]; irow < (this->colptr())[j+1]; irow++) {
			T_Int i = (this->rowidx())[irow];
			if(static_cast<T_Int>(ibgn) <= i && i < static_cast<T_Int>(iend)) {
				cptr[jlocal+1]++;
			} // i in range
		} // irow
	} // j

	blk::csc::roll(nj, cptr);

	T_Int     nnz  = cptr[nj];
	T_Int    *ridx = nullptr;
	T_Scalar *vals = nullptr;

	if(nnz) {

		ridx = i_malloc<T_Int>(nnz);
		vals = i_malloc<T_Scalar>(nnz);

		for(int_t j = jbgn; j < jend; j++) {
			int_t jlocal = j - jbgn;
			for(T_Int irow = (this->colptr())[j]; irow < (this->colptr())[j+1]; irow++) {
				T_Int i = (this->rowidx())[irow];
				T_Int ilocal = i - ibgn;
				T_Scalar v = (this->values())[irow];
				if(static_cast<T_Int>(ibgn) <= i && i < static_cast<T_Int>(iend)) {
					ridx[cptr[jlocal]] = ilocal;
					vals[cptr[jlocal]] = v;
					cptr[jlocal]++;
				} // i in range
			} // irow
		} // j

		blk::csc::unroll(nj, cptr);

	} // nnz

	XxMatrix<T_Int,T_Scalar> ret(ni, nj, cptr, ridx, vals, true, pr);

	return ret;
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
void XxMatrix<T_Int,T_Scalar>::checker() const
{
	csc_consistency_check(this->prop(), this->nrows(), this->ncols(), nnz(), this->colptr(), this->rowidx(), this->values());
}
/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
XxMatrix<T_Int,T_Scalar> XxMatrix<T_Int,T_Scalar>::random(int_t nr, int_t nc, int_t nz, const Property& pr, T_RScalar lo, T_RScalar hi)
{
	if(!nr || !nc)
		return XxMatrix<T_Int,T_Scalar>();

	coo::XxMatrix<T_Int,T_Scalar> Acoo(nr, nc, pr);

	int_t offDiagNnz = nz;

	/*
	 * Fill diagonal if needed
	 */
	if(pr.isSymmetric() || pr.isHermitian() || pr.isTriangular() || (pr.isGeneral() && nr == nc)) {

		int_t diagNnz = std::min(std::min(nr,nc),nz);
		offDiagNnz = nz - diagNnz;

		for(int_t j = 0; j < diagNnz; j++) {

			T_Scalar Ajj = rand<T_Scalar>(lo,hi);

			if(pr.isHermitian())
				arith::setIm(Ajj,0);

			Acoo.insert(j,j,Ajj);

		} // j

	} // sy/he

	/*
	 * Fill off-diagonal
	 * Do not use while loops
	 * Do not treat cases where i == j
	 * Trivial cases like 1x1 Skew are insignificant
	 */
	for(int_t k = 0; k < offDiagNnz; k++) {

		int_t iend = nr - 1;
		int_t jend = nc - 1;

		if(pr.isTriangular() && pr.isUpper() && nr > nc) iend = jend; 
		if(pr.isTriangular() && pr.isLower() && nr < nc) jend = iend; 

		T_Int i = rand<T_Int>(0,iend);
		T_Int j = rand<T_Int>(0,jend);

		if((pr.isUpper() && i > j) || (pr.isLower() && i < j))
			std::swap(i,j);

		T_Scalar Aij = rand<T_Scalar>(lo,hi);

		if(i == j && pr.isHermitian())
			arith::setIm(Aij,0);

		if(i == j && pr.isSkew())
			continue;

		Acoo.insert(i,j,Aij);

	} // off diag

	return Acoo.toCsc();
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
Guard<XxMatrix<T_Int,T_Scalar>> XxMatrix<T_Int,T_Scalar>::view(int_t nr, int_t nc, const T_Int *cptr, const T_Int *ridx, const T_Scalar *vals, const Property& pr)
{
	XxMatrix<T_Int,T_Scalar> tmp(nr, nc,
			const_cast<T_Int   *>(cptr),
			const_cast<T_Int   *>(ridx),
			const_cast<T_Scalar*>(vals), false, pr);
	Guard<XxMatrix<T_Int,T_Scalar>> ret(tmp);
	return ret;
}
/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/
template class XxMatrix<int_t,real_t>;
template class XxMatrix<int_t,real4_t>;
template class XxMatrix<int_t,complex_t>;
template class XxMatrix<int_t,complex8_t>;
/*-------------------------------------------------*/
} // namespace csc
} // namespace cla3p
/*-------------------------------------------------*/
