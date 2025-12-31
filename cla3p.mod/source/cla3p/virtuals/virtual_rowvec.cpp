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
#include "cla3p/virtuals/virtual_rowvec.hpp"

// system

// 3rd

// cla3p
#include "cla3p/checks/basic_checks.hpp"
#include "cla3p/checks/outer_checks.hpp"
#include "cla3p/checks/hermitian_coeff_checks.hpp"
#include "cla3p/proxies/blas_proxy.hpp"
#include "cla3p/dense/dns_xxvector.hpp"
#include "cla3p/dense/dns_xxmatrix.hpp"

/*-------------------------------------------------*/
namespace cla3p {
/*-------------------------------------------------*/
template <typename T_Scalar>
VirtualRowvec<T_Scalar>::VirtualRowvec(int_t n, const T_Scalar *vals, int_t incv, bool conj)
    : m_size(n), m_vals(vals), m_incv(incv), m_conj(conj)
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
  similarity_dim_check(m_size, vec.size());

  T_Scalar ret = 0;
  if(m_conj) {
    ret = blas::dotc(m_size, m_vals, m_incv, vec.values(), 1);
  } else {
    ret = blas::dot(m_size, m_vals, m_incv, vec.values(), 1);
  }

  return ret;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void VirtualRowvec<T_Scalar>::evaluateOuterOnNew(T_Scalar coeff, const dns::XxVector<T_Scalar>& vec, dns::XxMatrix<T_Scalar>& dest) const
{
	dest.clear();
	dest = dns::XxMatrix<T_Scalar>(m_size, vec.size());
	evaluateOuterOnExisting(coeff, vec, dest);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void VirtualRowvec<T_Scalar>::evaluateOuterOnExisting(T_Scalar coeff, const dns::XxVector<T_Scalar>& vec, dns::XxMatrix<T_Scalar>& dest) const
{
  dest = 0;
  accumulateOuterOnExisting(coeff, vec, dest);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void VirtualRowvec<T_Scalar>::accumulateOuterOnExisting(T_Scalar coeff, const dns::XxVector<T_Scalar>& vec, dns::XxMatrix<T_Scalar>& dest) const
{
  outer_product_consistency_check(m_conj, dest.nrows(), dest.ncols(), dest.prop(), vec.size(), m_size);
  hermitian_coeff_check(dest.prop(), coeff);

  T_Scalar beta = T_Scalar(1);

  if(dest.prop().isGeneral()) {

    if(m_conj) {
      blas::gerc(dest.nrows(), dest.ncols(), coeff, vec.values(), 1, m_vals, m_incv, dest.values(), dest.ld());
    } else {
      blas::ger(dest.nrows(), dest.ncols(), coeff, vec.values(), 1, m_vals, m_incv, dest.values(), dest.ld());
    }

  } else if(dest.prop().isSymmetric()) {

    blas::gemmt(dest.prop().cuplo(), 'N', 'T', dest.nrows(), 1, coeff, m_vals, m_incv, vec.values(), 1, beta, dest.values(), dest.ld());

  } else if(dest.prop().isHermitian()) { 

    blas::gemmt(dest.prop().cuplo(), 'N', 'C', dest.nrows(), 1, coeff, m_vals, m_incv, vec.values(), 1, beta, dest.values(), dest.ld());

  } // prop
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void VirtualRowvec<T_Scalar>::evaluateOnNew(dns::XxMatrix<T_Scalar>& dest) const
{
  dest = dns::XxMatrix<T_Scalar>(1, m_size);
  evaluateOnExisting(dest);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void VirtualRowvec<T_Scalar>::evaluateOnExisting(dns::XxMatrix<T_Scalar>& dest) const
{
  blas::copy(m_size, m_vals, m_incv, dest.values(), dest.ld());
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void VirtualRowvec<T_Scalar>::accumulateOnExisting(T_Scalar coeff, dns::XxMatrix<T_Scalar>& dest) const
{
  blas::axpy(m_size, coeff, m_vals, m_incv, dest.values(), dest.ld());
}
/*-------------------------------------------------*/
template class VirtualRowvec<real_t>;
template class VirtualRowvec<real4_t>;
template class VirtualRowvec<complex_t>;
template class VirtualRowvec<complex8_t>;
/*-------------------------------------------------*/
} // namespace cla3p
/*-------------------------------------------------*/
