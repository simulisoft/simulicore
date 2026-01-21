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
#include "cla3p/virtuals/virtual_strided_block.hpp"

// system

// 3rd

// cla3p
#include "cla3p/types/scalar.hpp"
#include "cla3p/bulk/dns.hpp"
#include "cla3p/checks/basic_checks.hpp"
#include "cla3p/dense/dns_xxvector.hpp"
#include "cla3p/dense/dns_xxmatrix.hpp"
#include "cla3p/sparse/csc_xxmatrix.hpp"

/*-------------------------------------------------*/
namespace cla3p {
/*-------------------------------------------------*/
template <typename T_Scalar>
VirtualStridedDnsVector<T_Scalar>::VirtualStridedDnsVector(int_t n, const T_Scalar *vals, int_t incv)
    : Meta1D<int_t>(n), m_values(vals), m_incv(incv)
{
}
/*-------------------------------------------------*/
template <typename T_Scalar>
VirtualStridedDnsVector<T_Scalar>::~VirtualStridedDnsVector()
{
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void VirtualStridedDnsVector<T_Scalar>::evaluateOnNew(dns::XxVector<T_Scalar>& dest) const
{
    dest = dns::XxVector<T_Scalar>(this->size());
    evaluateOnExisting(dest);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void VirtualStridedDnsVector<T_Scalar>::evaluateOnExisting(dns::XxVector<T_Scalar>& dest) const
{
    similarity_dim_check(this->size(), dest.size());
    blk::dns::copy_stride(uplo_t::Full, this->size(), 1, 
                          m_values, 1, m_incv,
                          dest.values(), 1, 1);
}
/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/
template <typename T_Scalar>
VirtualStridedDnsMatrix<T_Scalar>::VirtualStridedDnsMatrix(int_t m, int_t n, const T_Scalar *vals, int_t ldv, 
                                                           const Property& pr, int_t incv)
    : MatrixMeta<int_t>(m,n,pr), m_values(vals), m_incv(incv), m_ldv(ldv)
{
}
/*-------------------------------------------------*/
template <typename T_Scalar>
VirtualStridedDnsMatrix<T_Scalar>::~VirtualStridedDnsMatrix()
{
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void VirtualStridedDnsMatrix<T_Scalar>::evaluateOnNew(dns::XxMatrix<T_Scalar>& dest) const
{
    dest = dns::XxMatrix<T_Scalar>(this->nrows(), this->ncols(), this->prop());
    evaluateOnExisting(dest);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void VirtualStridedDnsMatrix<T_Scalar>::evaluateOnExisting(dns::XxMatrix<T_Scalar>& dest) const
{
    similarity_check(*this, dest);
    blk::dns::copy_stride(this->prop().uplo(), this->nrows(), this->ncols(), 
                          m_values, m_ldv, m_incv,
                          dest.values(), dest.ld(), 1);
}
/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
VirtualStridedCscMatrix<T_Int,T_Scalar>::VirtualStridedCscMatrix(T_Int m, T_Int n, 
                                                                 const T_Int *cptr, 
                                                                 const T_Int *ridx,
                                                                 const T_Scalar *vals, 
                                                                 const Property& pr, T_Int incv)
    : MatrixMeta<T_Int>(m,n,pr), m_colptr(cptr), m_rowidx(ridx), m_values(vals), m_incv(incv)
{
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
VirtualStridedCscMatrix<T_Int,T_Scalar>::~VirtualStridedCscMatrix()
{
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
void VirtualStridedCscMatrix<T_Int,T_Scalar>::evaluateOnNew(csc::XxMatrix<T_Int,T_Scalar>& dest) const
{
    dest = csc::XxMatrix<T_Int,T_Scalar>(this->nrows(), this->ncols(), nnz(), this->prop());
    evaluateOnExisting(dest);
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
void VirtualStridedCscMatrix<T_Int,T_Scalar>::evaluateOnExisting(csc::XxMatrix<T_Int,T_Scalar>& dest) const
{
    similarity_check(*this, dest);
    similarity_dim_check(nnz(), dest.nnz());

    blk::dns::copy_stride(uplo_t::Full, nnz(), 1, 
                          m_values, nnz(), m_incv,
                          dest.values(), dest.nnz(), 1);
}
/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/
template class VirtualStridedDnsVector<real_t>;
template class VirtualStridedDnsVector<real4_t>;
template class VirtualStridedDnsVector<complex_t>;
template class VirtualStridedDnsVector<complex8_t>;
/*-------------------------------------------------*/
template class VirtualStridedDnsMatrix<real_t>;
template class VirtualStridedDnsMatrix<real4_t>;
template class VirtualStridedDnsMatrix<complex_t>;
template class VirtualStridedDnsMatrix<complex8_t>;
/*-------------------------------------------------*/
template class VirtualStridedCscMatrix<int_t, real_t>;
template class VirtualStridedCscMatrix<int_t, real4_t>;
template class VirtualStridedCscMatrix<int_t, complex_t>;
template class VirtualStridedCscMatrix<int_t, complex8_t>;
/*-------------------------------------------------*/
} // namespace cla3p
/*-------------------------------------------------*/
