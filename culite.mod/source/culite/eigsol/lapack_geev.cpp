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
#include <culite/eigsol/lapack_geev.hpp>

// system

// 3rd

// culite
#include <culite/dense.hpp>

/*-------------------------------------------------*/
namespace culite {
/*-------------------------------------------------*/
template <typename T_Matrix>
LapackGeev<T_Matrix>::LapackGeev(CuSolverHandler& cusolver, bool calcLeft, bool calcRight)
    : m_cusolver(cusolver), m_calcLeft(calcLeft), m_calcRight(calcRight)
{
}
/*-------------------------------------------------*/
template <typename T_Matrix>
LapackGeev<T_Matrix>::~LapackGeev()
{
    clear();
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void LapackGeev<T_Matrix>::reserve(const T_Matrix& mat)
{
    m_cusolver.reserveGeev(mat, m_calcLeft, m_calcRight);

    clearOutput();

    std::size_t n = mat.nrows();
    std::size_t sizeN2 = n * n;

    std::size_t sizeE = n;
    std::size_t sizeVL = m_calcLeft ? sizeN2 : 0;
    std::size_t sizeVR = m_calcRight ? sizeN2 : 0;

    m_deviceBuffer.reserve(sizeE + sizeVL + sizeVR);

    T_CScalar* currentBuff = m_deviceBuffer.data();

    m_eigenvalues = dns::CxVector<T_CScalar>(n, currentBuff, false); currentBuff += sizeE;

    if(m_calcLeft ) { m_leftEigenvectors  = dns::CxMatrix<T_CScalar>(n, n, currentBuff, n, false); currentBuff += sizeVL; }
    if(m_calcRight) { m_rightEigenvectors = dns::CxMatrix<T_CScalar>(n, n, currentBuff, n, false); currentBuff += sizeVR; }
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void LapackGeev<T_Matrix>::clearOutput()
{
    m_eigenvalues.clear();
    m_leftEigenvectors.clear();
    m_rightEigenvectors.clear();
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void LapackGeev<T_Matrix>::clear()
{
    clearOutput();
    m_deviceBuffer.clear();
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void LapackGeev<T_Matrix>::decompose(const T_Matrix& mat)
{
    reserve(mat);
    
    m_cusolver.executeGeev(mat, m_calcLeft, m_calcRight);
    m_cusolver.geevGetEigenvalues(m_eigenvalues);
    m_cusolver.geevGetEigenvectors(m_calcLeft, m_calcRight, m_leftEigenvectors, m_rightEigenvectors);
}
/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/
template class LapackGeev<dns::RdMatrix>;
template class LapackGeev<dns::RfMatrix>;
template class LapackGeev<dns::CdMatrix>;
template class LapackGeev<dns::CfMatrix>;
/*-------------------------------------------------*/
} // namespace culite
/*-------------------------------------------------*/
