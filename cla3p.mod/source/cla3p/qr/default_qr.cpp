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
#include "cla3p/qr/default_qr.hpp"

// system

// 3rd

// cla3p
#include "cla3p/dense.hpp"
#include "cla3p/proxies/lapack_proxy.hpp"
#include "cla3p/checks/lapack_checks.hpp"

/*-------------------------------------------------*/
namespace cla3p {
/*-------------------------------------------------*/
template <typename T_Matrix>
DefaultQR<T_Matrix>::DefaultQR()
{
	defaults();
}
/*-------------------------------------------------*/
template <typename T_Matrix>
DefaultQR<T_Matrix>::DefaultQR(int_t m, int_t n, qrPolicy_t qrPolicy)
	: 
		DefaultQR<T_Matrix>::DefaultQR()
{
	setQrPolicy(qrPolicy);
	reserve(m, n);
}
/*-------------------------------------------------*/
template <typename T_Matrix>
DefaultQR<T_Matrix>::~DefaultQR()
{
	clear();
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void DefaultQR<T_Matrix>::defaults()
{
	setQrPolicy(qrPolicy_t::Reflection);
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void DefaultQR<T_Matrix>::clearInternalWrappers()
{
	m_matrixA.clear();
	m_tauVector.clear();
	m_matrixR.clear();
	m_matrixQ.clear();

	m_scalarBuffers.rewind();
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void DefaultQR<T_Matrix>::clearInternalBuffers()
{
	m_scalarBuffers.clear();
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void DefaultQR<T_Matrix>::clear()
{
	clearInternalWrappers();
	clearInternalBuffers();

	defaults();
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void DefaultQR<T_Matrix>::reserve(int_t m, int_t n)
{
	int_t k = std::min(m,n);

	std::size_t bufferSize = 
		m * n + // matrixBackup
		k +     // tau
		k * n;  // matrixR

	if(getQrPolicy() == qrPolicy_t::Full)
		bufferSize += m * k; // matrixQ

	m_scalarBuffers.resize(bufferSize);
}
/*-------------------------------------------------*/
template <typename T_Matrix>
qrPolicy_t DefaultQR<T_Matrix>::getQrPolicy() const
{
	return m_qrPolicy;
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void DefaultQR<T_Matrix>::setQrPolicy(qrPolicy_t qrPolicy)
{
	m_qrPolicy = qrPolicy;
}
/*-------------------------------------------------*/
template <typename T_Matrix>
const typename DefaultQR<T_Matrix>::T_Vector& DefaultQR<T_Matrix>::tau() const
{
	return m_tauVector; 
}
/*-------------------------------------------------*/
template <typename T_Matrix>
const T_Matrix& DefaultQR<T_Matrix>::R() const
{ 
	return m_matrixR; 
}
/*-------------------------------------------------*/
template <typename T_Matrix>
const T_Matrix& DefaultQR<T_Matrix>::Q() const
{ 
	return m_matrixQ; 
}
/*-------------------------------------------------*/
template <typename T_Matrix>
const T_Matrix& DefaultQR<T_Matrix>::elementaryReflectors() const
{ 
	return m_matrixA; 
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void DefaultQR<T_Matrix>::resizeInternalObjects(const T_Matrix& mat)
{
	int_t m = mat.nrows();
	int_t n = mat.ncols();

	reserve(m, n);

	int_t k = std::min(m,n);

	clearInternalWrappers();

	{
		m_matrixA = T_Matrix(m, n, m_scalarBuffers.request(m * n), m, false, mat.prop());
		m_tauVector = T_Vector(k, m_scalarBuffers.request(k), false);
		m_matrixR = T_Matrix(k, n, m_scalarBuffers.request(k * n), k, false);

		if(getQrPolicy() == qrPolicy_t::Full) {
			m_matrixQ = T_Matrix(m, k, m_scalarBuffers.request(m * k), m, false);
		}
	}
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void DefaultQR<T_Matrix>::decompose(const T_Matrix& mat)
{ 
	resizeInternalObjects(mat);

	m_matrixA = mat;
	m_matrixA.igeneral();

	int_t info = lapack::geqrf(
			m_matrixA.nrows(), 
			m_matrixA.ncols(), 
			m_matrixA.values(),
			m_matrixA.ld(),
			m_tauVector.values());

	lapack_info_check(info);

	fillMatrixR();
	fillMatrixQ();
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void DefaultQR<T_Matrix>::fillMatrixR()
{
	int_t m = m_matrixR.nrows();
	int_t n = m_matrixR.ncols();
	int_t k = std::min(m,n);

	m_matrixR = m_matrixA.rblock(0,0,m,n);

	for(int_t j = 0; j < k; j++) {
		for(int_t i = j+1; i < k; i++) {
			m_matrixR(i,j) = T_Scalar(0);
		} // i
	} // j
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void DefaultQR<T_Matrix>::fillMatrixQ()
{
	if(m_matrixQ) {

		int_t m = m_matrixQ.nrows();
		int_t n = m_matrixQ.ncols();

		m_matrixQ = m_matrixA.rblock(0,0,m,n);

		int_t info = lapack::xxgqr(m, n, n, 
				m_matrixQ.values(), 
				m_matrixQ.ld(), 
				m_tauVector.values());

		lapack_info_check(info);

	} // explicit calculation of Q
}
/*-------------------------------------------------*/
template class DefaultQR<dns::RdMatrix>;
template class DefaultQR<dns::RfMatrix>;
template class DefaultQR<dns::CdMatrix>;
template class DefaultQR<dns::CfMatrix>;
/*-------------------------------------------------*/
} // namespace cla3p
/*-------------------------------------------------*/
