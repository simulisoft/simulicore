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
#include "culite/virtuals/virtual_rowvec.hpp"

// system

// 3rd
#include <cla3p/checks/basic_checks.hpp>
#include <cla3p/checks/outer_checks.hpp>
#include <cla3p/checks/hermitian_coeff_checks.hpp>

// culite
#include "culite/support/utils.hpp"
#include "culite/generic/cublas_handler.hpp"
#include "culite/dense/dns_xxvector.hpp"
#include "culite/dense/dns_xxmatrix.hpp"

/*-------------------------------------------------*/
namespace culite {
/*-------------------------------------------------*/
template <typename T_Scalar>
VirtualRowvec<T_Scalar>::VirtualRowvec(int_t n, const T_Scalar *vals, int_t incv, bool conj)
    : VirtualRowvecBase<int_t,T_Scalar>(n, vals, incv, conj)
{
}
/*-------------------------------------------------*/
template <typename T_Scalar>
VirtualRowvec<T_Scalar>::~VirtualRowvec()
{
}
/*-------------------------------------------------*/
template <typename T_Scalar>
T_Scalar VirtualRowvec<T_Scalar>::evaluateInner(const dns::XxVector<T_Scalar>& vec) const
{
  ::cla3p::similarity_dim_check(this->size(), vec.size());

  T_Scalar ret = makeScalar<T_Scalar>(0);
  if(this->isConj()) {
    globalCuBlasHandler().dotc(this->size(), this->values(), this->incv(), vec.values(), 1, &ret);
  } else {
    globalCuBlasHandler().dot(this->size(), this->values(), this->incv(), vec.values(), 1, &ret);
  }

  return ret;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void VirtualRowvec<T_Scalar>::evaluateOuterOnNew(T_Scalar coeff, const dns::XxVector<T_Scalar>& vec, dns::XxMatrix<T_Scalar>& dest) const
{
    dest.clear();
    dest = dns::XxMatrix<T_Scalar>(this->size(), vec.size());
    evaluateOuterOnExisting(coeff, vec, dest);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void VirtualRowvec<T_Scalar>::evaluateOuterOnExisting(T_Scalar coeff, const dns::XxVector<T_Scalar>& vec, dns::XxMatrix<T_Scalar>& dest) const
{
  //memSetZero(dest.nrows(), dest.ncols(), dest.values(), dest.ld());
  dest = makeScalar<T_Scalar>(0);
  accumulateOuterOnExisting(coeff, vec, dest);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void VirtualRowvec<T_Scalar>::accumulateOuterOnExisting(T_Scalar coeff, const dns::XxVector<T_Scalar>& vec, dns::XxMatrix<T_Scalar>& dest) const
{
  ::cla3p::outer_product_consistency_check(this->isConj(), dest.nrows(), dest.ncols(), dest.prop(), vec.size(), this->size());
  ::cla3p::hermitian_coeff_check(dest.prop(), TypeTraits<T_Scalar>::toCla3pType(coeff));

  if(dest.prop().isGeneral()) {

    if(this->isConj()) {
      globalCuBlasHandler().gerc(dest.nrows(), dest.ncols(), &coeff, vec.values(), 1, this->values(), this->incv(), dest.values(), dest.ld());
    } else {
      globalCuBlasHandler().ger(dest.nrows(), dest.ncols(), &coeff, vec.values(), 1, this->values(), this->incv(), dest.values(), dest.ld());
    }

  } else if(dest.prop().isSymmetric()) {

    T_Scalar beta = makeScalar<T_Scalar>(1);
    globalCuBlasHandler().syrkx(dest.prop().uplo(), ::cla3p::op_t::N, dest.nrows(), 1, 
                                &coeff, 
                                this->values(), this->incv(), 
                                vec.values(), 1, 
                                &beta, dest.values(), dest.ld());

  } else if(dest.prop().isHermitian()) { 

    using T_RScalar = typename TypeTraits<T_Scalar>::real_type;
    T_RScalar beta = makeScalar<T_RScalar>(1);
    globalCuBlasHandler().herkx(dest.prop().uplo(), ::cla3p::op_t::N, dest.nrows(), 1, 
                                &coeff, 
                                this->values(), this->incv(), 
                                vec.values(), 1, 
                                &beta, dest.values(), dest.ld());

  } // prop
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void VirtualRowvec<T_Scalar>::evaluateOnNew(dns::XxMatrix<T_Scalar>& dest) const
{
  dest = dns::XxMatrix<T_Scalar>(1, this->size());
  evaluateOnExisting(dest);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void VirtualRowvec<T_Scalar>::evaluateOnExisting(dns::XxMatrix<T_Scalar>& dest) const
{
    ::cla3p::similarity_dim_check<int_t>(this->size(), dest.ncols());
    ::cla3p::similarity_dim_check<int_t>(1, dest.nrows());
    globalCuBlasHandler().copy(this->size(), this->values(), this->incv(), dest.values(), dest.ld());
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void VirtualRowvec<T_Scalar>::accumulateOnExisting(T_Scalar coeff, dns::XxMatrix<T_Scalar>& dest) const
{
    ::cla3p::similarity_dim_check<int_t>(this->size(), dest.ncols());
    ::cla3p::similarity_dim_check<int_t>(1, dest.nrows());
    globalCuBlasHandler().axpy(this->size(), &coeff, this->values(), this->incv(), dest.values(), dest.ld());
}
/*-------------------------------------------------*/
template class VirtualRowvec<real_t>;
template class VirtualRowvec<real4_t>;
template class VirtualRowvec<complex_t>;
template class VirtualRowvec<complex8_t>;
/*-------------------------------------------------*/
} // namespace culite
/*-------------------------------------------------*/
