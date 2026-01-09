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
#include "culite/dense/dns_cxvector.hpp"

// system

// 3rd

// culite
#include "culite/types/scalar.hpp"

/*-------------------------------------------------*/
namespace culite {
namespace dns {
/*-------------------------------------------------*/
template <typename T_Scalar>
CxVector<T_Scalar>::CxVector(const XiVector<T_Scalar>& other)
	: XxVector<T_Scalar>(other)
{
}
/*-------------------------------------------------*/
template <typename T_Scalar>
CxVector<T_Scalar>& CxVector<T_Scalar>::operator=(const XiVector<T_Scalar>& other)
{
	XxVector<T_Scalar>::operator=(other);
	return *this;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
CxVector<T_Scalar>::CxVector(XiVector<T_Scalar>&& other)
	: XxVector<T_Scalar>(std::move(other))
{
}
/*-------------------------------------------------*/
template <typename T_Scalar>
CxVector<T_Scalar>& CxVector<T_Scalar>::operator=(XiVector<T_Scalar>&& other)
{
	XxVector<T_Scalar>::operator=(std::move(other));
	return *this;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
CxVector<T_Scalar>::CxVector()
{
}
/*-------------------------------------------------*/
template <typename T_Scalar>
CxVector<T_Scalar>::CxVector(int_t n)
	: CxVector<T_Scalar>::XxVector(n)
{
}
/*-------------------------------------------------*/
template <typename T_Scalar>
CxVector<T_Scalar>::CxVector(int_t n, T_Scalar *vals, bool bind)
	: CxVector<T_Scalar>::XxVector(n, vals, bind)
{
}
/*-------------------------------------------------*/
template <typename T_Scalar>
CxVector<T_Scalar>::~CxVector()
{
}
/*-------------------------------------------------*/
#if 0
template <typename T_Scalar>
XxVector<typename TypeTraits<T_Scalar>::real_type> CxVector<T_Scalar>::real() const
{
	XxVector<T_RScalar> ret(this->size());
	blk::dns::get_real(uplo_t::Full, this->size(), 1, this->values(), this->size(), ret.values(), ret.size());
	return ret;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
XxVector<typename TypeTraits<T_Scalar>::real_type> CxVector<T_Scalar>::imag() const
{
	XxVector<T_RScalar> ret(this->size());
	blk::dns::get_imag(uplo_t::Full, this->size(), 1, this->values(), this->size(), ret.values(), ret.size());
	return ret;
}
#endif // 0
/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/
template class CxVector<complex_t>;
template class CxVector<complex8_t>;
/*-------------------------------------------------*/
} // namespace dns
} // namespace culite
/*-------------------------------------------------*/

