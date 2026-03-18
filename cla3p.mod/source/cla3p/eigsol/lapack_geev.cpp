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
#include <cla3p/eigsol/lapack_geev.hpp>

// system

// 3rd

// cla3p
#include <cla3p/types/scalar.hpp>
#include <cla3p/proxies/lapack_proxy.hpp>
#include <cla3p/proxies/blas_proxy.hpp>
#include <cla3p/checks/basic_checks.hpp>
#include <cla3p/checks/lapack_checks.hpp>
#include <cla3p/error/exceptions.hpp>
#include <cla3p/bulk/dns1d.hpp>
#include <cla3p/bulk/dns.hpp>
#include <cla3p/dense.hpp>

/*-------------------------------------------------*/
namespace cla3p {
/*-------------------------------------------------*/
template <typename T_Matrix>
LapackGeev<T_Matrix>::LapackGeev()
{
    defaults();
}
/*-------------------------------------------------*/
template <typename T_Matrix>
LapackGeev<T_Matrix>::LapackGeev(bool calcLeft, bool calcRight)
    : LapackGeev<T_Matrix>()
{
    setCalcLeft(calcLeft);
    setCalcRight(calcRight);
}
/*-------------------------------------------------*/
template <typename T_Matrix>
LapackGeev<T_Matrix>::~LapackGeev()
{
    clear();
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void LapackGeev<T_Matrix>::defaults()
{
    m_calcLeft = false;
    m_calcRight = true;
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void LapackGeev<T_Matrix>::clearInternalWrappers()
{
    m_eigenvalues.clear();
    m_leftEigenvectors.clear();
    m_rightEigenvectors.clear();
    m_matrixBackup.clear();

    m_scalarBuffers.rewind();
    m_complexBuffers.rewind();
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void LapackGeev<T_Matrix>::clearInternalBuffers()
{
    m_scalarBuffers.clear();
    m_complexBuffers.clear();
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void LapackGeev<T_Matrix>::clear()
{
    clearInternalWrappers();
    clearInternalBuffers();

    defaults();
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void LapackGeev<T_Matrix>::reserve(int_t n)
{
    int_t n2 = n * n;

    // scalar buffer: matrix backup (n*n) 
    //              + geev w array (2*n for real, n for complex) 
    //              + vl (n*n if calcLeft) 
    //              + vr (n*n if calcRight)
    std::size_t aSize = n2;
    std::size_t wSize = TypeTraits<T_Scalar>::is_real() ? (2 * n) : n;
    std::size_t vlSize = m_calcLeft ? n2 : 0;
    std::size_t vrSize = m_calcRight ? n2 : 0;
    m_scalarBuffers.resize(aSize + wSize + vlSize + vrSize);

    // complex buffer: eigenvalues (n) 
    //               + left eigenvectors (n*n if calcLeft) 
    //               + right eigenvectors (n*n if calcRight)
    std::size_t cwSize = n;
    std::size_t cvlSize = m_calcLeft ? n2 : 0;
    std::size_t cvrSize = m_calcRight ? n2 : 0;
    m_complexBuffers.resize(cwSize + cvlSize + cvrSize);
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void LapackGeev<T_Matrix>::resizeInternalObjects(int_t n)
{
    reserve(n);
    clearInternalWrappers();

    int_t n2 = n * n;

    m_matrixBackup = T_Matrix(n, n, m_scalarBuffers.request(n2), n, false);

    // complex output wrappers
    m_eigenvalues = T_CVector(n, m_complexBuffers.request(n), false);
    if(m_calcLeft)  m_leftEigenvectors  = T_CMatrix(n, n, m_complexBuffers.request(n2), n, false);
    if(m_calcRight) m_rightEigenvectors = T_CMatrix(n, n, m_complexBuffers.request(n2), n, false);
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void LapackGeev<T_Matrix>::setCalcLeft(bool calcLeft)
{
    m_calcLeft = calcLeft;
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void LapackGeev<T_Matrix>::setCalcRight(bool calcRight)
{
    m_calcRight = calcRight;
}
/*-------------------------------------------------*/
template <typename T_Matrix>
bool LapackGeev<T_Matrix>::getCalcLeft() const
{
    return m_calcLeft;
}
/*-------------------------------------------------*/
template <typename T_Matrix>
bool LapackGeev<T_Matrix>::getCalcRight() const
{
    return m_calcRight;
}
/*-------------------------------------------------*/
template <typename T_Matrix>
const typename LapackGeev<T_Matrix>::T_CVector& 
LapackGeev<T_Matrix>::eigenvalues() const
{
    return m_eigenvalues;
}
/*-------------------------------------------------*/
template <typename T_Matrix>
const typename LapackGeev<T_Matrix>::T_CMatrix& 
LapackGeev<T_Matrix>::leftEigenvectors() const
{
    return m_leftEigenvectors;
}
/*-------------------------------------------------*/
template <typename T_Matrix>
const typename LapackGeev<T_Matrix>::T_CMatrix& 
LapackGeev<T_Matrix>::rightEigenvectors() const
{
    return m_rightEigenvectors;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
static void copyEigenvaluesReal(int_t n, const T_Scalar *w, 
                                typename TypeTraits<T_Scalar>::complex_type *eigs)
{
    using T_CScalar = typename TypeTraits<T_Scalar>::complex_type;

    const T_Scalar *wr = w;
    const T_Scalar *wi = w + n;

    for(int_t i = 0; i < n; i++) {
        eigs[i] = T_CScalar(wr[i], wi[i]);
    }
}
/*-------------------------------------------------*/
template <typename T_Scalar>
static void copyEigenvaluesComplex(int_t n, const T_Scalar *w, T_Scalar *eigs)
{
    blk::dns::copy(n, w, eigs);
}
/*-------------------------------------------------*/
static void copyEigenvalues(int_t n, const real_t *w, complex_t *eigs) { copyEigenvaluesReal(n, w, eigs); }
static void copyEigenvalues(int_t n, const real4_t *w, complex8_t *eigs) { copyEigenvaluesReal(n, w, eigs); }
static void copyEigenvalues(int_t n, const complex_t *w, complex_t *eigs) { copyEigenvaluesComplex(n, w, eigs); }
static void copyEigenvalues(int_t n, const complex8_t *w, complex8_t *eigs) { copyEigenvaluesComplex(n, w, eigs); }
/*-------------------------------------------------*/
template <typename T_Scalar>
static void copyEigenvectorsReal(int_t n, const T_Scalar *w, const T_Scalar *v, int_t ldv,
                                 typename TypeTraits<T_Scalar>::complex_type *cv, int_t ldcv)
{
    blk::dns::zero(uplo_t::Full, n, n, cv, ldcv);

    const T_Scalar *wi = w + n;

    // dgeev/sgeev ordering: wi[j] > 0 => first of conjugate pair (j, j+1)
    //                       wi[j] < 0 => second of conjugate pair, skip
    //                       wi[j] = 0 => real eigenvalue
    #pragma omp parallel for schedule(static)
    for(int_t j = 0; j < n; j++) {

        if(wi[j] < T_Scalar(0)) continue; // second of conjugate pair, handled at j-1

        const T_Scalar *vj   = blk::dns::ptrmv(ldv, v, 0, j);
        T_Scalar       *cvRe = reinterpret_cast<T_Scalar*>(blk::dns::ptrmv(ldcv, cv, 0, j));
        T_Scalar       *cvIm = cvRe + 1;

        if(wi[j] == T_Scalar(0)) {

            // real eigenvalue: copy real column into real parts
            blas::copy(n, vj, 1, cvRe, 2);

        } else {

            // wi[j] > 0: first of conjugate pair, also fill column j+1
            const T_Scalar *vj1   = blk::dns::ptrmv(ldv, v, 0, j + 1);
            T_Scalar       *cvRe1 = reinterpret_cast<T_Scalar*>(blk::dns::ptrmv(ldcv, cv, 0, j + 1));
            T_Scalar       *cvIm1 = cvRe1 + 1;

            // column j:   re = v[j], im = +v[j+1]
            blas::copy(n, vj , 1, cvRe, 2);
            blas::copy(n, vj1, 1, cvIm, 2);

            // column j+1: re = v[j], im = -v[j+1]
            blas::copy(n, vj, 1, cvRe1, 2);
            blas::axpy(n, T_Scalar(-1), vj1, 1, cvIm1, 2);

        } // wi

    } // j
}
/*-------------------------------------------------*/
template <typename T_Scalar>
static void copyEigenvectorsComplex(int_t n, const T_Scalar *v, int_t ldv,
                                    T_Scalar *cv, int_t ldcv)
{
    blk::dns::copy(uplo_t::Full, n, n, v, ldv, cv, ldcv);
}
/*-------------------------------------------------*/
static void copyEigenvectors(int_t n, const real_t     *w, const real_t     *v, int_t ldv, complex_t  *cv, int_t ldcv) { copyEigenvectorsReal(n, w, v, ldv, cv, ldcv); }
static void copyEigenvectors(int_t n, const real4_t    *w, const real4_t    *v, int_t ldv, complex8_t *cv, int_t ldcv) { copyEigenvectorsReal(n, w, v, ldv, cv, ldcv); }
static void copyEigenvectors(int_t n, const complex_t  * , const complex_t  *v, int_t ldv, complex_t  *cv, int_t ldcv) { copyEigenvectorsComplex(n, v, ldv, cv, ldcv); }
static void copyEigenvectors(int_t n, const complex8_t * , const complex8_t *v, int_t ldv, complex8_t *cv, int_t ldcv) { copyEigenvectorsComplex(n, v, ldv, cv, ldcv); }
/*-------------------------------------------------*/
template <typename T_Matrix>
void LapackGeev<T_Matrix>::decompose(const T_Matrix& mat)
{
    if(!mat.prop().isGeneral()) {
        throw err::NoConsistency("LapackGeev: Input matrix must have 'General' property.");
    }
    square_check(mat);

    int_t n = mat.nrows();

    resizeInternalObjects(n);

    m_matrixBackup = mat;

    std::size_t wSize = TypeTraits<T_Scalar>::is_real() ? (2 * n) : n;
    T_Scalar *w  = m_scalarBuffers.request(wSize);
    T_Scalar *vl = m_calcLeft  ? m_scalarBuffers.request(n * n) : nullptr;
    T_Scalar *vr = m_calcRight ? m_scalarBuffers.request(n * n) : nullptr;

    char jobvl = m_calcLeft  ? 'V' : 'N';
    char jobvr = m_calcRight ? 'V' : 'N';

    int_t info = lapack::geev(jobvl, jobvr, n,
                              m_matrixBackup.values(), m_matrixBackup.ld(),
                              w,
                              vl, n,
                              vr, n);

    lapack_info_check(info);

    copyEigenvalues(n, w, m_eigenvalues.values());
    if(m_calcLeft)  copyEigenvectors(n, w, vl, n, m_leftEigenvectors.values(), m_leftEigenvectors.ld());
    if(m_calcRight) copyEigenvectors(n, w, vr, n, m_rightEigenvectors.values(), m_rightEigenvectors.ld());
}
/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/
template class LapackGeev<dns::RdMatrix>;
template class LapackGeev<dns::RfMatrix>;
template class LapackGeev<dns::CdMatrix>;
template class LapackGeev<dns::CfMatrix>;
/*-------------------------------------------------*/
} // namespace cla3p
/*-------------------------------------------------*/
