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
#include "cla3p/dense/dns_cxmatrix.hpp"

// system

// 3rd

// cla3p
#include "cla3p/bulk/dns.hpp"

/*-------------------------------------------------*/
namespace cla3p {
namespace dns {
/*-------------------------------------------------*/
template <typename T_Scalar>
CxMatrix<T_Scalar>::CxMatrix(const XxMatrix<T_Scalar>& other)
	: XxMatrix<T_Scalar>(other)
{
}
/*-------------------------------------------------*/
template <typename T_Scalar>
CxMatrix<T_Scalar>& CxMatrix<T_Scalar>::operator=(const XxMatrix<T_Scalar>& other)
{
	XxMatrix<T_Scalar>::operator=(other);
	return *this;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
CxMatrix<T_Scalar>::CxMatrix(XxMatrix<T_Scalar>&& other)
	: XxMatrix<T_Scalar>(std::move(other))
{
}
/*-------------------------------------------------*/
template <typename T_Scalar>
CxMatrix<T_Scalar>& CxMatrix<T_Scalar>::operator=(XxMatrix<T_Scalar>&& other)
{
	XxMatrix<T_Scalar>::operator=(std::move(other));
	return *this;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
CxMatrix<T_Scalar>::CxMatrix()
{
}
/*-------------------------------------------------*/
template <typename T_Scalar>
CxMatrix<T_Scalar>::CxMatrix(int_t nr, int_t nc, const Property& pr)
	: CxMatrix<T_Scalar>::XxMatrix(nr, nc, pr)
{
}
/*-------------------------------------------------*/
template <typename T_Scalar>
CxMatrix<T_Scalar>::CxMatrix(int_t nr, int_t nc, T_Scalar *vals, int_t ldv, bool bind, const Property& pr)
	: CxMatrix<T_Scalar>::XxMatrix(nr, nc, vals, ldv, bind, pr)
{
}
/*-------------------------------------------------*/
template <typename T_Scalar>
CxMatrix<T_Scalar>::~CxMatrix()
{
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void CxMatrix<T_Scalar>::operator=(T_Scalar val)
{
	this->fill(val);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
XxMatrix<typename TypeTraits<T_Scalar>::real_type> CxMatrix<T_Scalar>::real() const
{
	Property ret_prop = (this->prop().isHermitian() ? Property(prop_t::Symmetric, this->prop().uplo()) : this->prop());

	XxMatrix<T_RScalar> ret(this->nrows(), this->ncols(), ret_prop);

	blk::dns::get_real(
			this->prop().uplo(), 
			this->nrows(), 
			this->ncols(), 
			this->values(), 
			this->ld(), 
			ret.values(), ret.ld());

	return ret;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
XxMatrix<typename TypeTraits<T_Scalar>::real_type> CxMatrix<T_Scalar>::imag() const
{
	Property ret_prop = (this->prop().isHermitian() ? Property(prop_t::Skew, this->prop().uplo()) : this->prop());

	XxMatrix<T_RScalar> ret(this->nrows(), this->ncols(), ret_prop);

	blk::dns::get_imag(
			this->prop().uplo(), 
			this->nrows(), 
			this->ncols(), 
			this->values(), 
			this->ld(), 
			ret.values(), ret.ld());

	return ret;
}
/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/
template class CxMatrix<complex_t>;
template class CxMatrix<complex8_t>;
/*-------------------------------------------------*/
} // namespace dns
} // namespace cla3p
/*-------------------------------------------------*/

