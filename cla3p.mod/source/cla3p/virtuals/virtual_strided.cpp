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
#include "cla3p/virtuals/virtual_strided.hpp"

// system

// 3rd

// cla3p
#include "cla3p/types/scalar.hpp"
#include "cla3p/bulk/dns1D.hpp"
#include "cla3p/bulk/dns.hpp"
#include "cla3p/checks/basic_checks.hpp"
#include "cla3p/dense/dns_xxvector.hpp"
#include "cla3p/dense/dns_xxmatrix.hpp"
#include "cla3p/sparse/csr_xxmatrix.hpp"
#include "cla3p/sparse/csc_xxmatrix.hpp"

/*-------------------------------------------------*/
namespace cla3p {
/*-------------------------------------------------*/
template <typename T_Scalar>
VirtualStridedDnsVector<T_Scalar>::VirtualStridedDnsVector()
{
    defaults();
}
/*-------------------------------------------------*/
template <typename T_Scalar>
VirtualStridedDnsVector<T_Scalar>::VirtualStridedDnsVector(int_t n, T_Scalar *vals, int_t incv)
    : Meta1D<int_t>(n), m_values(vals), m_incv(incv)
{
}
/*-------------------------------------------------*/
template <typename T_Scalar>
Guard<VirtualStridedDnsVector<T_Scalar>> 
VirtualStridedDnsVector<T_Scalar>::view(int_t n, const T_Scalar *vals, int_t incv)
{
    VirtualStridedDnsVector<T_Scalar> tmp(n, const_cast<T_Scalar*>(vals), incv);
    Guard<VirtualStridedDnsVector<T_Scalar>> ret(tmp);
    return ret;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
VirtualStridedDnsVector<T_Scalar>::~VirtualStridedDnsVector()
{
    clear();
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void VirtualStridedDnsVector<T_Scalar>::defaults() 
{
    m_values = nullptr;
    m_incv = 0;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void VirtualStridedDnsVector<T_Scalar>::clear() 
{ 
    Meta1D<int_t>::clear();
    defaults(); 
}
/*-------------------------------------------------*/
template <typename T_Scalar>
VirtualStridedDnsVector<T_Scalar> VirtualStridedDnsVector<T_Scalar>::rcopy()
{
    return VirtualStridedDnsVector<T_Scalar>(this->size(), m_values, m_incv);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
VirtualStridedDnsVector<T_Scalar>& 
VirtualStridedDnsVector<T_Scalar>::operator=(const VirtualStridedDnsVector& other)
{
    if(this != &other) {
        similarity_dim_check(this->size(), other.size());
        blk::dns::copy(this->size(), other.m_values, m_values, other.m_incv, m_incv);
    }
    return *this;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void VirtualStridedDnsVector<T_Scalar>::operator=(const dns::XxVector<T_Scalar>& vec)
{
    *this = view(vec.size(), vec.values(), 1).get();
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void VirtualStridedDnsVector<T_Scalar>::operator=(T_Scalar val)
{
    blk::dns::fill(this->size(), m_values, val, m_incv);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void VirtualStridedDnsVector<T_Scalar>::evaluateOnNew(dns::XxVector<T_Scalar>& dest) const
{
    dest.clear();
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
VirtualStridedDnsMatrix<T_Scalar>::VirtualStridedDnsMatrix()
{
    defaults();
}
/*-------------------------------------------------*/
template <typename T_Scalar>
VirtualStridedDnsMatrix<T_Scalar>::VirtualStridedDnsMatrix(int_t m, int_t n, 
                                                           T_Scalar *vals, int_t ldv, 
                                                           const Property& pr, int_t incv)
    : MatrixMeta<int_t>(m,n,pr), m_values(vals), m_ldv(ldv), m_incv(incv)
{
}
/*-------------------------------------------------*/
template <typename T_Scalar>
Guard<VirtualStridedDnsMatrix<T_Scalar>> 
VirtualStridedDnsMatrix<T_Scalar>::view(int_t m, int_t n, 
                                        const T_Scalar *vals, int_t ldv, 
                                        const Property& pr, int_t incv)
{
    VirtualStridedDnsMatrix<T_Scalar> tmp(m, n, const_cast<T_Scalar*>(vals), ldv, pr, incv);
    Guard<VirtualStridedDnsMatrix<T_Scalar>> ret(tmp);
    return ret;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
VirtualStridedDnsMatrix<T_Scalar>::~VirtualStridedDnsMatrix()
{
    clear();
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void VirtualStridedDnsMatrix<T_Scalar>::defaults()
{
    m_values = nullptr;
    m_ldv = 0;
    m_incv = 0;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void VirtualStridedDnsMatrix<T_Scalar>::clear() 
{ 
    MatrixMeta<int_t>::clear();
    defaults(); 
}
/*-------------------------------------------------*/
template <typename T_Scalar>
VirtualStridedDnsMatrix<T_Scalar>
VirtualStridedDnsMatrix<T_Scalar>::rcopy()
{
    return VirtualStridedDnsMatrix<T_Scalar>(this->nrows(), this->ncols(), 
                                             m_values, m_ldv, 
                                             this->prop(), m_incv);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
VirtualStridedDnsMatrix<T_Scalar>& 
VirtualStridedDnsMatrix<T_Scalar>::operator=(const VirtualStridedDnsMatrix& other)
{
    if(this != &other) {
        similarity_check(*this, other);
        blk::dns::copy_stride(other.prop().uplo(), other.nrows(), other.ncols(), 
                              other.m_values, other.m_ldv, other.m_incv,
                              m_values, m_ldv, m_incv);
    }
    return *this;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void VirtualStridedDnsMatrix<T_Scalar>::operator=(const dns::XxMatrix<T_Scalar>& mat)
{
    *this = view(mat.nrows(), mat.ncols(), mat.values(), mat.ld(), mat.prop(), 1).get();
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void VirtualStridedDnsMatrix<T_Scalar>::operator=(T_Scalar val)
{
    blk::dns::fill_stride(this->prop().uplo(), this->nrows(), this->ncols(), 
                          m_values, m_ldv, m_incv, val);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void VirtualStridedDnsMatrix<T_Scalar>::evaluateOnNew(dns::XxMatrix<T_Scalar>& dest) const
{
    dest.clear();
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
VirtualStridedCsrMatrix<T_Int,T_Scalar>::VirtualStridedCsrMatrix()
{
    defaults();
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
VirtualStridedCsrMatrix<T_Int,T_Scalar>::VirtualStridedCsrMatrix(T_Int m, T_Int n, 
                                                                 T_Int *rptr, 
                                                                 T_Int *cidx,
                                                                 T_Scalar *vals, 
                                                                 const Property& pr, T_Int incv)
    : MatrixMeta<T_Int>(m,n,pr), m_rowptr(rptr), m_colidx(cidx), m_values(vals), m_incv(incv)
{
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
VirtualStridedCsrMatrix<T_Int,T_Scalar>::~VirtualStridedCsrMatrix()
{
    clear();
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
void VirtualStridedCsrMatrix<T_Int,T_Scalar>::defaults()
{
    m_rowptr = nullptr;
    m_colidx = nullptr;
    m_values = nullptr;
    m_incv = 0;
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
void VirtualStridedCsrMatrix<T_Int,T_Scalar>::clear()
{
    MatrixMeta<T_Int>::clear();
    defaults();
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
VirtualStridedCsrMatrix<T_Int,T_Scalar>
VirtualStridedCsrMatrix<T_Int,T_Scalar>::rcopy()
{
    return VirtualStridedCsrMatrix<T_Int,T_Scalar>(this->nrows(), this->ncols(), 
                                                     m_rowptr, m_colidx, m_values, 
                                                     this->prop(), m_incv);
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
Guard<VirtualStridedCsrMatrix<T_Int,T_Scalar>> 
VirtualStridedCsrMatrix<T_Int,T_Scalar>::view(T_Int m, T_Int n, 
                                              const T_Int *rptr, 
                                              const T_Int *cidx, 
                                              const T_Scalar *vals, 
                                              const Property& pr, T_Int incv)
{
    VirtualStridedCsrMatrix<T_Int,T_Scalar> tmp(m, n, 
                                                const_cast<T_Int*>(rptr), 
                                                const_cast<T_Int*>(cidx), 
                                                const_cast<T_Scalar*>(vals), 
                                                pr, incv);
    Guard<VirtualStridedCsrMatrix<T_Int,T_Scalar>> ret(tmp);
    return ret;
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
VirtualStridedCsrMatrix<T_Int,T_Scalar>& 
VirtualStridedCsrMatrix<T_Int,T_Scalar>::operator=(const VirtualStridedCsrMatrix& other)
{
    if(this != &other) {
        similarity_check(*this, other);
        similarity_dim_check(nnz(), other.nnz());
        blk::dns::copy(other.nrows() + 1, other.m_rowptr, m_rowptr);
        blk::dns::copy(other.nnz(), other.m_colidx, m_colidx);
        blk::dns::copy(other.nnz(), other.m_values, m_values, other.m_incv, m_incv);
    }
    return *this;
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
void VirtualStridedCsrMatrix<T_Int,T_Scalar>::operator=(const csr::XxMatrix<T_Int,T_Scalar>& mat)
{
    *this = view(mat.nrows(), 
                 mat.ncols(), 
                 mat.rowptr(), 
                 mat.colidx(), 
                 mat.values(),
                 mat.prop(), 1).get();
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
void VirtualStridedCsrMatrix<T_Int,T_Scalar>::operator=(T_Scalar val)
{
    blk::dns::fill(nnz(), m_values, val, m_incv);
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
void VirtualStridedCsrMatrix<T_Int,T_Scalar>::evaluateOnNew(csr::XxMatrix<T_Int,T_Scalar>& dest) const
{
    dest.clear();
    dest = csr::XxMatrix<T_Int,T_Scalar>(this->nrows(), this->ncols(), nnz(), this->prop());
    evaluateOnExisting(dest);
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
void VirtualStridedCsrMatrix<T_Int,T_Scalar>::evaluateOnExisting(csr::XxMatrix<T_Int,T_Scalar>& dest) const
{
    similarity_check(*this, dest);
    similarity_dim_check(nnz(), dest.nnz());
    blk::dns::copy(this->nrows() + 1, m_rowptr, dest.rowptr());
    blk::dns::copy(nnz(), m_colidx, dest.colidx());
    blk::dns::copy(nnz(), m_values, dest.values(), m_incv, 1);
}
/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
VirtualStridedCscMatrix<T_Int,T_Scalar>::VirtualStridedCscMatrix()
{
    defaults();
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
VirtualStridedCscMatrix<T_Int,T_Scalar>::VirtualStridedCscMatrix(T_Int m, T_Int n, 
                                                                 T_Int *cptr, 
                                                                 T_Int *ridx,
                                                                 T_Scalar *vals, 
                                                                 const Property& pr, T_Int incv)
    : MatrixMeta<T_Int>(m,n,pr), m_colptr(cptr), m_rowidx(ridx), m_values(vals), m_incv(incv)
{
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
VirtualStridedCscMatrix<T_Int,T_Scalar>::~VirtualStridedCscMatrix()
{
    clear();
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
void VirtualStridedCscMatrix<T_Int,T_Scalar>::defaults()
{
    m_colptr = nullptr;
    m_rowidx = nullptr;
    m_values = nullptr;
    m_incv = 0;
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
void VirtualStridedCscMatrix<T_Int,T_Scalar>::clear()
{
    MatrixMeta<T_Int>::clear();
    defaults();
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
VirtualStridedCscMatrix<T_Int,T_Scalar> 
VirtualStridedCscMatrix<T_Int,T_Scalar>::rcopy()
{
    return VirtualStridedCscMatrix<T_Int,T_Scalar>(this->nrows(), this->ncols(), 
                                                   m_colptr, m_rowidx, m_values, 
                                                   this->prop(), m_incv);
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
Guard<VirtualStridedCscMatrix<T_Int,T_Scalar>> 
VirtualStridedCscMatrix<T_Int,T_Scalar>::view(T_Int m, T_Int n, 
                                              const T_Int *cptr, 
                                              const T_Int *ridx, 
                                              const T_Scalar *vals, 
                                              const Property& pr, T_Int incv)
{
    VirtualStridedCscMatrix<T_Int,T_Scalar> tmp(m, n, 
                                                const_cast<T_Int*>(cptr), 
                                                const_cast<T_Int*>(ridx), 
                                                const_cast<T_Scalar*>(vals), 
                                                pr, incv);
    Guard<VirtualStridedCscMatrix<T_Int,T_Scalar>> ret(tmp);
    return ret;
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
VirtualStridedCscMatrix<T_Int,T_Scalar>& 
VirtualStridedCscMatrix<T_Int,T_Scalar>::operator=(const VirtualStridedCscMatrix<T_Int,T_Scalar>& other)
{
    if(this != &other) {
        similarity_check(*this, other);
        similarity_dim_check(nnz(), other.nnz());
        blk::dns::copy(other.ncols() + 1, other.m_colptr, m_colptr);
        blk::dns::copy(other.nnz(), other.m_rowidx, m_rowidx);
        blk::dns::copy(other.nnz(), other.m_values, m_values, other.m_incv, m_incv);
    }
    return *this;
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
void VirtualStridedCscMatrix<T_Int,T_Scalar>::operator=(const csc::XxMatrix<T_Int,T_Scalar>& mat)
{
    *this = view(mat.nrows(), 
                 mat.ncols(), 
                 mat.colptr(), 
                 mat.rowidx(), 
                 mat.values(),
                 mat.prop(), 1).get();
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
void VirtualStridedCscMatrix<T_Int,T_Scalar>::operator=(T_Scalar val)
{
    blk::dns::fill(nnz(), m_values, val, m_incv);
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
void VirtualStridedCscMatrix<T_Int,T_Scalar>::evaluateOnNew(csc::XxMatrix<T_Int,T_Scalar>& dest) const
{
    dest.clear();
    dest = csc::XxMatrix<T_Int,T_Scalar>(this->nrows(), this->ncols(), nnz(), this->prop());
    evaluateOnExisting(dest);
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
void VirtualStridedCscMatrix<T_Int,T_Scalar>::evaluateOnExisting(csc::XxMatrix<T_Int,T_Scalar>& dest) const
{
    similarity_check(*this, dest);
    similarity_dim_check(nnz(), dest.nnz());
    blk::dns::copy(this->ncols() + 1, m_colptr, dest.colptr());
    blk::dns::copy(nnz(), m_rowidx, dest.rowidx());
    blk::dns::copy(nnz(), m_values, dest.values(), m_incv, 1);
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
template class VirtualStridedCsrMatrix<int_t, real_t>;
template class VirtualStridedCsrMatrix<int_t, real4_t>;
template class VirtualStridedCsrMatrix<int_t, complex_t>;
template class VirtualStridedCsrMatrix<int_t, complex8_t>;
/*-------------------------------------------------*/
template class VirtualStridedCscMatrix<int_t, real_t>;
template class VirtualStridedCscMatrix<int_t, real4_t>;
template class VirtualStridedCscMatrix<int_t, complex_t>;
template class VirtualStridedCscMatrix<int_t, complex8_t>;
/*-------------------------------------------------*/
} // namespace cla3p
/*-------------------------------------------------*/
