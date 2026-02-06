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
#include "cla3p/sparse/csr_xxmatrix.hpp"

// system

// 3rd

// cla3p
#include "cla3p/sparse/coo_xxmatrix.hpp"
#include "cla3p/perms/pxmatrix.hpp"

#include "cla3p/bulk/dns.hpp"
#include "cla3p/bulk/csr.hpp"
#include "cla3p/support/imalloc.hpp"
#include "cla3p/support/rand.hpp"
#include "cla3p/support/utils.hpp"

#include "cla3p/checks/basic_checks.hpp"
#include "cla3p/checks/csx_checks.hpp"
#include "cla3p/checks/block_ops_checks.hpp"
#include "cla3p/checks/perm_checks.hpp"
#include "cla3p/checks/hermitian_coeff_checks.hpp"

/*-------------------------------------------------*/
namespace cla3p {
namespace csr {
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
XxMatrix<T_Int,T_Scalar>::XxMatrix()
{
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
XxMatrix<T_Int,T_Scalar>::XxMatrix(T_Int nr, T_Int nc, T_Int nz, const Property& pr)
	: MatrixMeta<T_Int>(nr, nc, sanitizeProperty<T_Scalar>(pr)), csx::XxContainer<T_Int,T_Scalar>(nr, nz)
{
	if(nr > 0 && nc > 0) {
		checker();
	} else {
		clear();
	}
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
XxMatrix<T_Int,T_Scalar>::XxMatrix(T_Int nr, T_Int nc, T_Int *rptr, T_Int *cidx, T_Scalar *vals, bool bind, const Property& pr)
	: MatrixMeta<T_Int>(nr, nc, sanitizeProperty<T_Scalar>(pr)), csx::XxContainer<T_Int,T_Scalar>(rptr, cidx, vals, bind)
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
	MatrixMeta<T_Int>::clear();
	csx::XxContainer<T_Int,T_Scalar>::clear();
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
alias::VirtualScal_csr<T_Int,T_Scalar> XxMatrix<T_Int,T_Scalar>::operator-() const
{
	return alias::VirtualScal_csr<T_Int,T_Scalar>(virtualize(), T_Scalar(-1));
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
	if(!this->empty()) {
		return rowptr()[this->nrows()];
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
    ss << "  Index Precision...... " << TypeTraits<T_Int>::prec_name() << "\n";
	ss << "  Number of rows....... " << this->nrows() << "\n";
	ss << "  Number of columns.... " << this->ncols() << "\n";
	ss << "  Number of non zeros.. " << nnz() << "\n";
	ss << "  Rowptr............... " << rowptr() << "\n";
	ss << "  Colidx............... " << colidx() << "\n";
	ss << "  Values............... " << this->values() << "\n";
	ss << "  Property............. " << this->prop() << "\n";
	ss << "  Owner................ " << boolToYesNo(this->owner()) << "\n";

	ss << bottom << "\n";

	return ss.str();
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
XxMatrix<T_Int,T_Scalar>& XxMatrix<T_Int,T_Scalar>::copyFromExisting(const XxMatrix<T_Int,T_Scalar>& other)
{
	if(this != &other) {

        similarity_check(*this, other);
		similarity_dim_check(nnz(), other.nnz());

		T_Int nr = other.nrows() + 1;
		T_Int nz = other.nnz();

		// 
		// TODO: perhaps use a copy for 1D arrays
		//
		blk::dns::copy(uplo_t::Full, nr, 1, other.rowptr(), nr, rowptr(), nr);
		blk::dns::copy(uplo_t::Full, nz, 1, other.colidx(), nz, colidx(), nz);
		blk::dns::copy(uplo_t::Full, nz, 1, other.values(), nz, this->values(), nz);

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
			MatrixMeta<T_Int>::operator=(std::move(other));
			csx::XxContainer<T_Int,T_Scalar>::operator=(std::move(other));
			other.unbind();
		} // similar

		other.clear();

	} // do not apply on self

    return *this;
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
void XxMatrix<T_Int,T_Scalar>::toStream(std::ostream& os, std::streamsize prec) const
{
	blk::csr::print_to_stream(os, this->nrows(), rowptr(), colidx(), this->values(), prec);
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
	return XxMatrix<T_Int,T_Scalar>(this->nrows(), 
                                    this->ncols(), 
                                    rowptr(), 
                                    colidx(), 
                                    this->values(), 
                                    false, 
                                    this->prop());
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
Guard<XxMatrix<T_Int,T_Scalar>> XxMatrix<T_Int,T_Scalar>::rcopy() const
{
	return view(this->nrows(), 
                this->ncols(), 
                rowptr(), 
                colidx(), 
                this->values(), 
                this->prop());
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
alias::VirtualTrans_csr<T_Int,T_Scalar> XxMatrix<T_Int,T_Scalar>::transpose() const
{
	return alias::VirtualTrans_csr<T_Int,T_Scalar>(*this, false);
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
alias::VirtualTrans_csr<T_Int,T_Scalar> XxMatrix<T_Int,T_Scalar>::ctranspose() const
{
	return alias::VirtualTrans_csr<T_Int,T_Scalar>(*this, true);
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
alias::VirtualConj_csr<T_Int,T_Scalar> XxMatrix<T_Int,T_Scalar>::conjugate() const
{
	return alias::VirtualConj_csr<T_Int,T_Scalar>(*this);
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
	return blk::csr::norm_one(
			this->prop().type(),
			this->nrows(),
			this->ncols(),
			rowptr(), 
			colidx(), 
			this->values());
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
typename XxMatrix<T_Int,T_Scalar>::T_RScalar XxMatrix<T_Int,T_Scalar>::normInf() const
{
	return blk::csr::norm_inf(
			this->prop().type(),
			this->nrows(),
			rowptr(), 
			colidx(), 
			this->values());
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
typename XxMatrix<T_Int,T_Scalar>::T_RScalar XxMatrix<T_Int,T_Scalar>::normMax() const
{
	return blk::csr::norm_max(
			this->nrows(),
			rowptr(), 
			this->values());
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
typename XxMatrix<T_Int,T_Scalar>::T_RScalar XxMatrix<T_Int,T_Scalar>::normFro() const
{
	return blk::csr::norm_fro(
			this->prop().type(),
			this->nrows(),
			rowptr(), 
			colidx(), 
			this->values());
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
XxMatrix<T_Int,T_Scalar> XxMatrix<T_Int,T_Scalar>::general() const
{
	XxMatrix<T_Int,T_Scalar> ret;

	T_Int    *rowptr_ge = nullptr;
	T_Int    *colidx_ge = nullptr;
	T_Scalar *values_ge = nullptr;

	if(this->prop().isGeneral()) {

		ret = copy();

	} else if(this->prop().isSymmetric()) {

		rowptr_ge = i_malloc_t<T_Int>(this->nrows() + 1);
		blk::csr::uplo2ge_rowptr(this->prop().uplo(), this->nrows(), rowptr(), colidx(), rowptr_ge);
		T_Int nz = rowptr_ge[this->nrows()];
		colidx_ge = i_malloc_t<T_Int>(nz);
		values_ge = i_malloc_t<T_Scalar>(nz);
		blk::csr::sy2ge(this->prop().uplo(), this->nrows(), rowptr(), colidx(), this->values(), rowptr_ge, colidx_ge, values_ge);
		ret = XxMatrix<T_Int,T_Scalar>(this->nrows(), this->ncols(), rowptr_ge, colidx_ge, values_ge, true);

	} else if(this->prop().isHermitian()) {

		rowptr_ge = i_malloc_t<T_Int>(this->nrows() + 1);
		blk::csr::uplo2ge_rowptr(this->prop().uplo(), this->nrows(), rowptr(), colidx(), rowptr_ge);
		T_Int nz = rowptr_ge[this->nrows()];
		colidx_ge = i_malloc_t<T_Int>(nz);
		values_ge = i_malloc_t<T_Scalar>(nz);
		blk::csr::he2ge(this->prop().uplo(), this->nrows(), rowptr(), colidx(), this->values(), rowptr_ge, colidx_ge, values_ge);
		ret = XxMatrix<T_Int,T_Scalar>(this->nrows(), this->ncols(), rowptr_ge, colidx_ge, values_ge, true);

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
	for(T_Int i = 0; i < this->nrows(); i++) {
		for(T_Int jcol = rowptr()[i]; jcol < rowptr()[i+1]; jcol++) {
			ret(i, colidx()[jcol]) = (this->values())[jcol];
		} // jcol
	} // i
	return ret;
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
XxMatrix<T_Int,T_Scalar> XxMatrix<T_Int,T_Scalar>::permuteLeftRight(const prm::PxMatrix<T_Int>& P, const prm::PxMatrix<T_Int>& Q) const
{
	perm_ge_op_consistency_check(this->prop().type(), this->nrows(), this->ncols(), P.size(), Q.size());

	XxMatrix<T_Int,T_Scalar> ret(this->nrows(), this->ncols(), nnz(), this->prop());

	blk::csr::permute(this->prop().type(), 
                      this->prop().uplo(), 
                      this->nrows(), 
                      this->ncols(), 
			          rowptr(), 
                      colidx(), 
                      this->values(), 
			          ret.rowptr(), 
                      ret.colidx(), 
                      ret.values(), 
			          P.values(), 
                      Q.values());

	return ret;
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
XxMatrix<T_Int,T_Scalar> XxMatrix<T_Int,T_Scalar>::permuteLeft(const prm::PxMatrix<T_Int>& P) const
{
	perm_ge_op_consistency_check(this->prop().type(), this->nrows(), this->ncols(), P.size(), this->ncols());

	XxMatrix<T_Int,T_Scalar> ret(this->nrows(), this->ncols(), nnz(), this->prop());

	blk::csr::permute(this->prop().type(), 
                      this->prop().uplo(), 
                      this->nrows(), 
                      this->ncols(), 
			          rowptr(), 
                      colidx(), 
                      this->values(), 
			          ret.rowptr(), 
                      ret.colidx(), 
                      ret.values(), 
			          P.values(), 
                      nullptr);

	return ret;
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
XxMatrix<T_Int,T_Scalar> XxMatrix<T_Int,T_Scalar>::permuteRight(const prm::PxMatrix<T_Int>& Q) const
{
	perm_ge_op_consistency_check(this->prop().type(), this->nrows(), this->ncols(), this->nrows(), Q.size());

	XxMatrix<T_Int,T_Scalar> ret(this->nrows(), this->ncols(), nnz(), this->prop());

	blk::csr::permute(this->prop().type(), 
                      this->prop().uplo(), 
                      this->nrows(), 
                      this->ncols(), 
			          rowptr(), 
                      colidx(), 
                      this->values(), 
			          ret.rowptr(), 
                      ret.colidx(), 
                      ret.values(), 
			          nullptr, 
                      Q.values());

	return ret;
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
XxMatrix<T_Int,T_Scalar> XxMatrix<T_Int,T_Scalar>::permuteMirror(const prm::PxMatrix<T_Int>& P) const
{
	perm_op_consistency_check(this->nrows(), this->ncols(), P.size(), P.size());

	XxMatrix<T_Int,T_Scalar> ret(this->nrows(), this->ncols(), nnz(), this->prop());

	blk::csr::permute(this->prop().type(), 
                      this->prop().uplo(),
                      this->nrows(), 
                      this->ncols(), 
			          rowptr(), 
                      colidx(), 
                      this->values(), 
			          ret.rowptr(), 
                      ret.colidx(), 
                      ret.values(), 
			          P.values(), 
                      nullptr);

	return ret;
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
XxMatrix<T_Int,T_Scalar> XxMatrix<T_Int,T_Scalar>::block(T_Int ibgn, T_Int jbgn, T_Int ni, T_Int nj) const
{
	Property pr = block_op_consistency_check(this->prop(), this->nrows(), this->ncols(), ibgn, jbgn, ni, nj);

	if(ni <= 0 || nj <= 0) return XxMatrix<T_Int,T_Scalar>();

	T_Int *rptr = i_calloc_t<T_Int>(ni + 1);

	T_Int iend = ibgn + ni;
	T_Int jend = jbgn + nj;

	for(T_Int i = ibgn; i < iend; i++) {
		T_Int ilocal = i - ibgn;
		for(T_Int jcol = rowptr()[i]; jcol < rowptr()[i+1]; jcol++) {
			T_Int j = colidx()[jcol];
			if(jbgn <= j && j < jend) {
				rptr[ilocal+1]++;
			} // j in range
		} // jcol
	} // i

	blk::csx::roll(ni, rptr);

	T_Int     nnz  = rptr[ni];
	T_Int    *cidx = nullptr;
	T_Scalar *vals = nullptr;

	if(nnz) {

		cidx = i_malloc_t<T_Int>(nnz);
		vals = i_malloc_t<T_Scalar>(nnz);

		for(T_Int i = ibgn; i < iend; i++) {
			T_Int ilocal = i - ibgn;
			for(T_Int jcol = rowptr()[i]; jcol < rowptr()[i+1]; jcol++) {
				T_Int j = colidx()[jcol];
				T_Int jlocal = j - jbgn;
				T_Scalar v = this->values()[jcol];
				if(jbgn <= j && j < jend) {
                    cidx[rptr[ilocal]] = jlocal;
					vals[rptr[ilocal]] = v;
					rptr[ilocal]++;
				} // j in range
			} // jcol
		} // i

		blk::csx::unroll(ni, rptr);

	} // nnz

	XxMatrix<T_Int,T_Scalar> ret(ni, nj, rptr, cidx, vals, true, pr);
	return ret;
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
void XxMatrix<T_Int,T_Scalar>::checker() const
{
	csx_consistency_check(this->nrows(), 
                          this->ncols(), 
                          nnz(), 
                          rowptr(), 
                          colidx(), 
                          this->values(), 
                          this->prop());
}
/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
XxMatrix<T_Int,T_Scalar> XxMatrix<T_Int,T_Scalar>::random(T_Int nr, T_Int nc, T_Int nz, const Property& pr, T_RScalar lo, T_RScalar hi)
{
	if(nr <= 0 || nc <= 0 || nz < 0)
		return XxMatrix<T_Int,T_Scalar>();

	coo::XxMatrix<T_Int,T_Scalar> Acoo(nr, nc, pr);

	T_Int offDiagNnz = nz;

	/*
	 * Fill diagonal if needed
	 */
	if(pr.isSymmetric() || pr.isHermitian() || pr.isTriangular() || (pr.isGeneral() && nr == nc)) {

		T_Int diagNnz = std::min(std::min(nr,nc),nz);
		offDiagNnz = nz - diagNnz;

		for(T_Int j = 0; j < diagNnz; j++) {

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
	for(T_Int k = 0; k < offDiagNnz; k++) {

		T_Int iend = nr - 1;
		T_Int jend = nc - 1;

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

	return Acoo.toCsr();
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
Guard<XxMatrix<T_Int,T_Scalar>> XxMatrix<T_Int,T_Scalar>::view(T_Int nr, T_Int nc, const T_Int *rptr, const T_Int *cidx, const T_Scalar *vals, const Property& pr)
{
	XxMatrix<T_Int,T_Scalar> tmp(nr, nc,
			const_cast<T_Int   *>(rptr),
			const_cast<T_Int   *>(cidx),
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
} // namespace csr
} // namespace cla3p
/*-------------------------------------------------*/
