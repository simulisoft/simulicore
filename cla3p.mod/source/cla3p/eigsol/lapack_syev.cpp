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
#include <cla3p/eigsol/lapack_syev.hpp>

// system

// 3rd

// cla3p
#include <cla3p/types/scalar.hpp>
#include <cla3p/proxies/lapack_proxy.hpp>
#include <cla3p/checks/basic_checks.hpp>
#include <cla3p/checks/lapack_checks.hpp>
#include <cla3p/error/exceptions.hpp>
#include <cla3p/bulk/dns.hpp>
#include <cla3p/dense.hpp>

/*-------------------------------------------------*/
namespace cla3p {
/*-------------------------------------------------*/
template <typename T_Matrix>
LapackSyev<T_Matrix>::LapackSyev()
{
    defaults();
}
/*-------------------------------------------------*/
template <typename T_Matrix>
LapackSyev<T_Matrix>::LapackSyev(bool calcEigenvectors)
    : LapackSyev<T_Matrix>()
{
    setCalcEigenvectors(calcEigenvectors);
}
/*-------------------------------------------------*/
template <typename T_Matrix>
LapackSyev<T_Matrix>::~LapackSyev()
{
    clear();
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void LapackSyev<T_Matrix>::defaults()
{
    m_calcEigenvectors = true;
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void LapackSyev<T_Matrix>::clearInternalWrappers()
{
    m_eigenvalues.clear();
    m_eigenvectors.clear();
    m_matrixBackup.clear();

    m_realBuffers.rewind();
    m_scalarBuffers.rewind();
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void LapackSyev<T_Matrix>::clearInternalBuffers()
{
    m_realBuffers.clear();
    m_scalarBuffers.clear();
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void LapackSyev<T_Matrix>::clear()
{
    clearInternalWrappers();
    clearInternalBuffers();

    defaults();
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void LapackSyev<T_Matrix>::reserve(int_t n)
{
    m_realBuffers.resize(n);
    m_scalarBuffers.resize(n * n);
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void LapackSyev<T_Matrix>::resizeInternalObjects(int_t n)
{
    reserve(n);
    clearInternalWrappers();

    int_t n2 = n * n;
    
    m_matrixBackup = T_Matrix(n, n, m_scalarBuffers.request(n2), n, false);
    m_eigenvalues = T_RVector(n, m_realBuffers.request(n), false);
    if(m_calcEigenvectors) m_eigenvectors = m_matrixBackup.rcopy();
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void LapackSyev<T_Matrix>::setCalcEigenvectors(bool calcEigenvectors)
{
    m_calcEigenvectors = calcEigenvectors;
}
/*-------------------------------------------------*/
template <typename T_Matrix>
bool LapackSyev<T_Matrix>::getCalcEigenvectors() const
{
    return m_calcEigenvectors;
}
/*-------------------------------------------------*/
template <typename T_Matrix>
const typename LapackSyev<T_Matrix>::T_RVector& 
LapackSyev<T_Matrix>::eigenvalues() const
{
    return m_eigenvalues;
}
/*-------------------------------------------------*/
template <typename T_Matrix>
const T_Matrix& LapackSyev<T_Matrix>::eigenvectors() const
{
    return m_eigenvectors;
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void LapackSyev<T_Matrix>::decompose(const T_Matrix& mat)
{
    if(TypeTraits<T_Scalar>::is_real() && !mat.prop().isSymmetric()) {
        throw err::NoConsistency("LapackSyev: Real input matrix must have 'Symmetric' property.");
    }
    if(TypeTraits<T_Scalar>::is_complex() && !mat.prop().isHermitian()) {
        throw err::NoConsistency("LapackSyev: Complex input matrix must have 'Hermitian' property.");
    }
    square_check(mat);

    int_t n = mat.nrows();

    resizeInternalObjects(n);

    // Use bulk copying to avoid property mismatch error trigger
    blk::dns::copy(mat.prop().uplo(), 
                   mat.nrows(), mat.ncols(), mat.values(), mat.ld(),
                   m_matrixBackup.values(), m_matrixBackup.ld());

    char jobz = m_calcEigenvectors  ? 'V' : 'N';

    int_t info = lapack::heev(jobz, mat.prop().cuplo(), n, 
                              m_matrixBackup.values(), m_matrixBackup.ld(), 
                              m_eigenvalues.values());

    lapack_info_check(info);
}
/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/
template class LapackSyev<dns::RdMatrix>;
template class LapackSyev<dns::RfMatrix>;
template class LapackSyev<dns::CdMatrix>;
template class LapackSyev<dns::CfMatrix>;
/*-------------------------------------------------*/
} // namespace cla3p
/*-------------------------------------------------*/
