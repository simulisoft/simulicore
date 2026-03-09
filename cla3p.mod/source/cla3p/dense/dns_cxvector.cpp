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
#include <cla3p/dense/dns_cxvector.hpp>

// system

// 3rd

// cla3p
#include <cla3p/bulk/dns.hpp>

/*-------------------------------------------------*/
namespace cla3p {
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
template <typename T_Scalar>
void CxVector<T_Scalar>::operator=(T_Scalar val)
{
    this->fill(val);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
alias::VirtualStrided_vec<typename TypeTraits<T_Scalar>::real_type> 
CxVector<T_Scalar>::real()
{
    T_RScalar *realValues = reinterpret_cast<T_RScalar*>(this->values());
    return alias::VirtualStrided_vec<T_RScalar>(this->size(), realValues, 2);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
alias::VirtualStrided_vec<typename TypeTraits<T_Scalar>::real_type> 
CxVector<T_Scalar>::imag()
{
    T_RScalar *imagValues = reinterpret_cast<T_RScalar*>(this->values()) + 1;
    return alias::VirtualStrided_vec<T_RScalar>(this->size(), imagValues, 2);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
alias::GuardedStrided_vec<typename TypeTraits<T_Scalar>::real_type> 
CxVector<T_Scalar>::real() const
{
    const T_RScalar *realValues = reinterpret_cast<const T_RScalar*>(this->values());
    return alias::VirtualStrided_vec<T_RScalar>::view(this->size(), realValues, 2);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
alias::GuardedStrided_vec<typename TypeTraits<T_Scalar>::real_type> 
CxVector<T_Scalar>::imag() const
{
    const T_RScalar *imagValues = reinterpret_cast<const T_RScalar*>(this->values()) + 1;
    return alias::VirtualStrided_vec<T_RScalar>::view(this->size(), imagValues, 2);
}
/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/
template class CxVector<complex_t>;
template class CxVector<complex8_t>;
/*-------------------------------------------------*/
} // namespace dns
} // namespace cla3p
/*-------------------------------------------------*/

