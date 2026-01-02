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
#include "cla3p/qr/partial_qr.hpp"

// system
#include <cstring>
#include <algorithm>

// 3rd

// cla3p
#include <cla3p/dense.hpp>
#include <cla3p/proxies/lapack_proxy.hpp>

/*-------------------------------------------------*/
namespace cla3p {
/*-------------------------------------------------*/
template <typename T_Matrix>
PartialQR<T_Matrix>::PartialQR()
{
	defaults();
}
/*-------------------------------------------------*/
template <typename T_Matrix>
PartialQR<T_Matrix>::PartialQR(int_t m, int_t n)
	:
		PartialQR<T_Matrix>::PartialQR()
{
	reserve(m,n);
}
/*-------------------------------------------------*/
template <typename T_Matrix>
PartialQR<T_Matrix>::~PartialQR()
{
	clear();
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void PartialQR<T_Matrix>::clearInternalWrappers()
{
	m_matrixP.clear();
	m_matrixA.clear();
	m_matrixF.clear();
	m_matrixR.clear();
	m_tauVector.clear();
	m_auxVector.clear();
	m_vn1Vector.clear();
	m_vn2Vector.clear();

	m_scalarBuffer.rewind();
	m_realBuffer.rewind();
	m_intBuffer.rewind();
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void PartialQR<T_Matrix>::clearInternalBuffers()
{
	m_scalarBuffer.clear();
	m_realBuffer.clear();
	m_intBuffer.clear();
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void PartialQR<T_Matrix>::clear()
{
	clearInternalWrappers();
	clearInternalBuffers();

	defaults();
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void PartialQR<T_Matrix>::defaults()
{
	m_cutTolerance = 1.e-6;
	m_blockSize = 0;
	m_transFlag = false;
	m_offset = 0;
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void PartialQR<T_Matrix>::reserve(int_t m, int_t n)
{
	if(m < n)
		std::swap(m,n);

	setOptimumBlockSize(m, n);

	std::size_t scalarBufferSize = 
		m * n +           // matrixA
		n * m_blockSize + // matrixF
		n * n +           // matrixR
		n +               // tauVector
		m_blockSize;      // auxVector

	m_scalarBuffer.resize(scalarBufferSize);
	m_realBuffer.resize(2 * n);
	m_intBuffer.resize(n);
}
/*-------------------------------------------------*/
template <typename T_Matrix>
int_t PartialQR<T_Matrix>::totalColumnsProcessed() const
{
	return m_offset;
}
/*-------------------------------------------------*/
template <typename T_Matrix>
bool PartialQR<T_Matrix>::allColumnsProcessed() const
{
	return (totalColumnsProcessed() >= m_matrixA.ncols());
}
/*-------------------------------------------------*/
template <typename T_Matrix>
real_t PartialQR<T_Matrix>::cutTolerance() const
{
	return m_cutTolerance;
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void PartialQR<T_Matrix>::setCutTolerance(real_t tol)
{
	m_cutTolerance = tol;
}
/*-------------------------------------------------*/
template <typename T_Matrix>
const T_Matrix& PartialQR<T_Matrix>::R() const
{
	return m_matrixR;
}
/*-------------------------------------------------*/
template <typename T_Matrix>
const prm::PiMatrix& PartialQR<T_Matrix>::P() const
{
	return m_matrixP;
}
/*-------------------------------------------------*/
template <typename T_Matrix>
const T_Matrix& PartialQR<T_Matrix>::elementaryReflectors() const
{
	return m_matrixA;
}
/*-------------------------------------------------*/
template <typename T_Matrix>
const typename PartialQR<T_Matrix>::T_Vector& PartialQR<T_Matrix>::tau() const
{
	return m_tauVector;
}
/*-------------------------------------------------*/
template <typename T_Matrix>
bool PartialQR<T_Matrix>::transFlag() const
{
	return m_transFlag;
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void PartialQR<T_Matrix>::initializeDecompositionProcess(const T_Matrix& mat)
{
	resizeInternalObjects(mat);

	if(mat.nrows() >= mat.ncols()) {
		m_matrixA = mat;
		m_transFlag = false;
	} else {
		m_matrixA = mat.ctranspose();
		m_transFlag = true;
	} // m?n

	m_offset = 0;
	m_tauVector = 0;

	for(int_t j = 0; j < m_matrixA.ncols(); j++) {
		m_vn1Vector(j) = m_matrixA.rcolumn(j).normEuc();
	} // j

	m_vn2Vector = m_vn1Vector;

	for(int_t j = 0; j < m_matrixP.size(); j++) {
		m_matrixP(j) = j;
	} // j
}
/*-------------------------------------------------*/
template <typename T_Matrix>
int_t PartialQR<T_Matrix>::performDecompositionStep()
{
	int_t m = m_matrixA.nrows();
	int_t n = m_matrixA.ncols();
	int_t nb = m_blockSize;

	pivotSwitchBase1();

	int_t noff = totalColumnsProcessed();
	int_t rem = n - noff;
	int_t kb = 0;

	if(rem >= nb) {

		lapack::laqps(m, n - noff, noff, nb, &kb, 
				m_matrixA.rcolumn(noff).values(), m_matrixA.ld(),
				m_matrixP.values() + noff,
				m_tauVector.values() + noff,
				m_vn1Vector.values() + noff,
				m_vn2Vector.values() + noff,
				m_auxVector.values(), 
				m_matrixF.values(), n - noff);

	} else if(rem) {

		lapack::laqp2(m, n - noff, noff, 
				m_matrixA.rcolumn(noff).values(), m_matrixA.ld(),
				m_matrixP.values() + noff,
				m_tauVector.values() + noff,
				m_vn1Vector.values() + noff,
				m_vn2Vector.values() + noff,
				m_auxVector.values()); 

		kb = rem;

	} // rem?nb

	pivotSwitchBase0();

	m_offset += kb;

	return kb;
}
/*-------------------------------------------------*/
template <typename T_Matrix>
int_t PartialQR<T_Matrix>::toleranceCriterionReached(int_t numColsProcessedOnStep)
{
	int_t ret = 0;

	T_RScalar absQ00 = std::abs(m_matrixA(0,0));

	for(int_t i = totalColumnsProcessed() - numColsProcessedOnStep; i < totalColumnsProcessed(); i++) {

		T_Scalar Qii = m_matrixA(i,i);

		if(std::abs(Qii) < absQ00 * T_RScalar(m_cutTolerance)) {
			ret = i;
			break;
		} // tol

	} // i

	if(!ret && allColumnsProcessed()) {
		ret = m_matrixA.ncols();
	} // return if entire matrix is factored

	return ret;
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void PartialQR<T_Matrix>::fillMatrixR(int_t numRanks)
{
	numRanks = std::min(numRanks, totalColumnsProcessed());

	if(!numRanks)
		numRanks = totalColumnsProcessed();

	int_t n = m_matrixA.ncols();

	if(numRanks && n) {

		m_matrixR = T_Matrix(numRanks, n, m_matrixR.values(), numRanks, false);
		m_matrixR = m_matrixA.rblock(0,0,numRanks,n);

		for(int_t j = 0; j < numRanks; j++) {
			for(int_t i = j+1; i < numRanks; i++) {
				m_matrixR(i,j) = T_Scalar(0);
			} // i
		} // j

	} // valid dims
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void PartialQR<T_Matrix>::resume(bool ignoreTol)
{
	int_t ranksFound = 0;

	do {

		int_t numColsProcessedOnStep = performDecompositionStep();

		if(ignoreTol) {
			ranksFound = (allColumnsProcessed() ? m_matrixA.ncols() : 0);
		} else {
			ranksFound = toleranceCriterionReached(numColsProcessedOnStep);
		} // ignoreTol

	} while(!ranksFound);

	fillMatrixR(ranksFound);
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void PartialQR<T_Matrix>::decompose(const T_Matrix& mat)
{
	initializeDecompositionProcess(mat);
	resume();
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void PartialQR<T_Matrix>::setOptimumBlockSize(int_t m, int_t n)
{
	int_t n1 = std::max(m,n);
	int_t n2 = std::min(m,n);

	char funName[] = "XGEQRF";
	funName[0] = TypeTraits<T_Scalar>::netlibChar();
	funName[0] = std::toupper(funName[0]);
	
	m_blockSize = lapack::laenv(1, funName, " ", n1, n2, -1, -1);

	if(m_blockSize <= 0)
		m_blockSize = std::min(int_t(24), n2);
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void PartialQR<T_Matrix>::resizeInternalObjects(const T_Matrix& mat)
{
	int_t m = mat.nrows();
	int_t n = mat.ncols();

	if(m < n)
		std::swap(m,n);

	reserve(m,n);

	clearInternalWrappers();

	m_matrixP = prm::PiMatrix(n, m_intBuffer.request(n), false);

	m_matrixA = T_Matrix(m, n, m_scalarBuffer.request(m * n), m, false);
	m_matrixF = T_Matrix(n, m_blockSize, m_scalarBuffer.request(n * m_blockSize), n, false);
	m_matrixR = T_Matrix(n, n, m_scalarBuffer.request(n * n), n, false);
	m_tauVector = T_Vector(n, m_scalarBuffer.request(n), false);
	m_auxVector = T_Vector(m_blockSize, m_scalarBuffer.request(m_blockSize), false);

	m_vn1Vector = T_RVector(n, m_realBuffer.request(n), false);
	m_vn2Vector = T_RVector(n, m_realBuffer.request(n), false);
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void PartialQR<T_Matrix>::pivotSwitchBase0()
{
	std::for_each(m_matrixP.values() + totalColumnsProcessed(), m_matrixP.values() + m_matrixP.size(), [](int_t& j) { j--; });
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void PartialQR<T_Matrix>::pivotSwitchBase1()
{
	std::for_each(m_matrixP.values() + totalColumnsProcessed(), m_matrixP.values() + m_matrixP.size(), [](int_t& j) { j++; });
}
/*-------------------------------------------------*/
template class PartialQR<dns::RdMatrix>;
template class PartialQR<dns::RfMatrix>;
template class PartialQR<dns::CdMatrix>;
template class PartialQR<dns::CfMatrix>;
/*-------------------------------------------------*/
} // namespace cla3p
/*-------------------------------------------------*/

