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
#include <culite/svd/default_svd.hpp>

// system

// 3rd

// culite
#include <culite/dense.hpp>

/*-------------------------------------------------*/
namespace culite {
/*-------------------------------------------------*/
template <typename T_Matrix>
DefaultSVD<T_Matrix>::DefaultSVD(CuSolverHandler& cusolver)
    : m_cusolver(cusolver)
{
    defaults();
}
/*-------------------------------------------------*/
template <typename T_Matrix>
DefaultSVD<T_Matrix>::DefaultSVD(svdPolicy_t leftPolicy, svdPolicy_t rightPolicy, CuSolverHandler& cusolver)
    : DefaultSVD<T_Matrix>::DefaultSVD(cusolver)
{
    setLeftPolicy(leftPolicy);
    setRightPolicy(rightPolicy);
}
/*-------------------------------------------------*/
template <typename T_Matrix>
DefaultSVD<T_Matrix>::~DefaultSVD()
{
    clear();
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void DefaultSVD<T_Matrix>::defaults()
{
    m_LeftPolicy = svdPolicy_t::Economy;
    m_RightPolicy = svdPolicy_t::Economy;
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void DefaultSVD<T_Matrix>::reserve(const T_Matrix& mat)
{
    m_cusolver.reserveGesvd(mat, m_LeftPolicy, m_RightPolicy);

    clearOutput();

    int_t m = mat.nrows();
    int_t n = mat.ncols();
    int_t k = std::min(m, n);

    int_t sizeS = k;
    int_t sizeLeft = m_cusolver.svdVectorSize(m_LeftPolicy, m, k);
    int_t sizeRight = m_cusolver.svdVectorSize(m_RightPolicy, n, k);

    m_deviceValuesBuffer.reserve(sizeS);
    m_deviceVectorBuffer.reserve(sizeLeft + sizeRight);

    m_singularValues = T_RVector(k, m_deviceValuesBuffer.data(), false);

    T_Scalar *currentVectorBuff = m_deviceVectorBuffer.data();

    if(m_LeftPolicy != svdPolicy_t::NoCalculation) { 
        int_t nrowsU = m;
        int_t ncolsU = (m_LeftPolicy == svdPolicy_t::Full) ? m : k;
        m_leftSingularVectors = T_Matrix(nrowsU, ncolsU, currentVectorBuff, nrowsU, false); 
        currentVectorBuff += sizeLeft; 
    }

    if(m_RightPolicy != svdPolicy_t::NoCalculation) { 
        int_t nrowsV = (m_RightPolicy == svdPolicy_t::Full) ? n : k;
        int_t ncolsV = n;
        m_rightSingularVectors = T_Matrix(nrowsV, ncolsV, currentVectorBuff, nrowsV, false); 
        currentVectorBuff += sizeRight; 
    }
}
///*-------------------------------------------------*/
template <typename T_Matrix>
void DefaultSVD<T_Matrix>::clearOutput()
{
    m_singularValues.clear();
    m_leftSingularVectors.clear();
    m_rightSingularVectors.clear();
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void DefaultSVD<T_Matrix>::clear()
{
    clearOutput();
    m_deviceValuesBuffer.clear();
    m_deviceVectorBuffer.clear();

    defaults();
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void DefaultSVD<T_Matrix>::decompose(const T_Matrix& mat)
{
    reserve(mat);
    
    m_cusolver.executeGesvd(mat, m_LeftPolicy, m_RightPolicy);
    m_cusolver.gesvdGetSingularValues(m_singularValues);
    m_cusolver.gesvdGetSingularVectors(m_LeftPolicy, m_RightPolicy, m_leftSingularVectors, m_rightSingularVectors);
}
/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/
template class DefaultSVD<dns::RdMatrix>;
template class DefaultSVD<dns::RfMatrix>;
template class DefaultSVD<dns::CdMatrix>;
template class DefaultSVD<dns::CfMatrix>;
/*-------------------------------------------------*/
} // namespace culite
/*-------------------------------------------------*/
