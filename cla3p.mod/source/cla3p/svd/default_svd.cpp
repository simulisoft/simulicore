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
#include "cla3p/svd/default_svd.hpp"

// system

// 3rd

// cla3p
#include "cla3p/dense.hpp"
#include "cla3p/error/exceptions.hpp"
#include "cla3p/proxies/lapack_proxy.hpp"
#include "cla3p/checks/lapack_checks.hpp"

/*-------------------------------------------------*/
namespace cla3p {
/*-------------------------------------------------*/
template <typename T_Matrix>
DefaultSVD<T_Matrix>::DefaultSVD()
{
	defaults();
}
/*-------------------------------------------------*/
template <typename T_Matrix>
DefaultSVD<T_Matrix>::DefaultSVD(int_t m, int_t n, svdPolicy_t leftPolicy, svdPolicy_t rightPolicy)
	: 
		DefaultSVD<T_Matrix>::DefaultSVD()
{
	setLeftPolicy(leftPolicy);
	setRightPolicy(rightPolicy);
	reserve(m, n);
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
	setLeftPolicy(svdPolicy_t::Limited);
	setRightPolicy(svdPolicy_t::Limited);
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void DefaultSVD<T_Matrix>::clearInternalWrappers()
{
	m_singularValues.clear();
	m_superbVector.clear();
	m_matrixBackup.clear();
	m_leftSingularVectors.clear();
	m_rightSingularVectors.clear();

	m_realBuffers.rewind();
	m_scalarBuffers.rewind();
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void DefaultSVD<T_Matrix>::clearInternalBuffers()
{
	m_realBuffers.clear();
	m_scalarBuffers.clear();
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void DefaultSVD<T_Matrix>::clear()
{
	clearInternalWrappers();
	clearInternalBuffers();

	defaults();
}
/*-------------------------------------------------*/
static int_t policyToNumCols(svdPolicy_t policy, int_t n, int_t k)
{
	int_t ret = 0;

	if(policy == svdPolicy_t::Limited)
		ret = k;
	else if(policy == svdPolicy_t::Full)
		ret = n;

	return ret;
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void DefaultSVD<T_Matrix>::reserve(int_t m, int_t n)
{
	int_t k = std::min(m,n);

	int_t leftCols = policyToNumCols(getLeftPolicy(), m, k);
	int_t rightCols = policyToNumCols(getRightPolicy(), n, k);

	m_realBuffers.resize(2 * k - 1);
	m_scalarBuffers.resize(m * n + m * leftCols + n * rightCols);
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void DefaultSVD<T_Matrix>::setLeftPolicy(svdPolicy_t leftPolicy)
{
	m_leftPolicy = leftPolicy; 
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void DefaultSVD<T_Matrix>::setRightPolicy(svdPolicy_t rightPolicy)
{
	m_rightPolicy = rightPolicy; 
}
/*-------------------------------------------------*/
template <typename T_Matrix>
svdPolicy_t DefaultSVD<T_Matrix>::getLeftPolicy() const
{
	return m_leftPolicy; 
}
/*-------------------------------------------------*/
template <typename T_Matrix>
svdPolicy_t DefaultSVD<T_Matrix>::getRightPolicy() const
{ 
	return m_rightPolicy; 
}
/*-------------------------------------------------*/
template <typename T_Matrix>
const typename DefaultSVD<T_Matrix>::T_RVector& DefaultSVD<T_Matrix>::singularValues() const
{
	return m_singularValues; 
}
/*-------------------------------------------------*/
template <typename T_Matrix>
const T_Matrix& DefaultSVD<T_Matrix>::leftSingularVectors() const
{ 
	return m_leftSingularVectors; 
}
/*-------------------------------------------------*/
template <typename T_Matrix>
const T_Matrix& DefaultSVD<T_Matrix>::rightSingularVectors() const
{ 
	return m_rightSingularVectors; 
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void DefaultSVD<T_Matrix>::resizeInternalObjects(const T_Matrix& mat)
{
	int_t m = mat.nrows();
	int_t n = mat.ncols();

	reserve(m, n);

	int_t k = std::min(m,n);

	int_t leftCols = policyToNumCols(getLeftPolicy(), m, k);
	int_t rightCols = policyToNumCols(getRightPolicy(), n, k);

	clearInternalWrappers();

	{
		m_singularValues = T_RVector(k, m_realBuffers.request(k), false);
		m_superbVector = T_RVector(k - 1, m_realBuffers.request(k-1), false);
	}

	{
		m_matrixBackup = T_Matrix(m, n, m_scalarBuffers.request(m * n), m, false, mat.prop());

		if(leftCols) {
			m_leftSingularVectors = T_Matrix(m, leftCols, m_scalarBuffers.request(m * leftCols), m, false);
		}

		if(rightCols) {
			m_rightSingularVectors = T_Matrix(rightCols, n, m_scalarBuffers.request(n * rightCols), rightCols, false);
		}
	}
}
/*-------------------------------------------------*/
template <typename T_Matrix>
static void recursive_transpose_right_singular_vectors(T_Matrix& Vt, T_Matrix& A)
{
	using T_Scalar = typename T_Matrix::value_type;

	int_t m = A.nrows();
	int_t n = A.ncols();
	int_t k = Vt.nrows();

	if(Vt.ncols() != n)
		throw err::NoConsistency("Invalid dimensions for right singular vectors matrix");

	if(k <= m) {

		T_Matrix Akn = A.rblock(0,0,k,n);
		Akn = Vt;

		T_Scalar *values = Vt.values();
		Vt.clear();
		Vt = T_Matrix(n, k, values, n, false);
		Vt = Akn.ctranspose();

	} else {

		// 
		// in this scope: k == n && n > m
		//

		T_Matrix Amm = A.rblock(0,0,m,m);
		T_Matrix Vmm = Vt.rblock(0,0,m,m);
		Amm = Vmm;
		Vmm = Amm.ctranspose();

		int_t l = n - m;

		if(l) {

			T_Matrix Aml = A.rblock(0,m,m,l);
			T_Matrix Vml = Vt.rblock(0,m,m,l);
			T_Matrix Vlm = Vt.rblock(m,0,l,m);
			T_Matrix Vll = Vt.rblock(m,m,l,l);
			Aml = Vml;
			Vml = Vlm.ctranspose();
			Vlm = Aml.ctranspose();

			recursive_transpose_right_singular_vectors(Vll, Aml);

		} // l

	} // m?n
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void DefaultSVD<T_Matrix>::decompose(const T_Matrix& mat)
{ 
	resizeInternalObjects(mat);

	m_matrixBackup = mat;
	m_matrixBackup.igeneral();

	int_t minld = 1;

	int_t info = lapack::gesvd(
			static_cast<char>(getLeftPolicy()), 
			static_cast<char>(getRightPolicy()), 
			m_matrixBackup.nrows(), 
			m_matrixBackup.ncols(), 
			m_matrixBackup.values(),
			m_matrixBackup.ld(),
			m_singularValues.values(), 
			m_leftSingularVectors.values(), 
			std::max(minld, m_leftSingularVectors.ld()), // ld == 0 leads to lapack error
			m_rightSingularVectors.values(), 
			std::max(minld, m_rightSingularVectors.ld()), // ld == 0 leads to lapack error
			m_superbVector.values());

	lapack_info_check(info);

	if(m_rightSingularVectors) 
		recursive_transpose_right_singular_vectors(m_rightSingularVectors, m_matrixBackup);
}
/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/
template class DefaultSVD<dns::RdMatrix>;
template class DefaultSVD<dns::RfMatrix>;
template class DefaultSVD<dns::CdMatrix>;
template class DefaultSVD<dns::CfMatrix>;
/*-------------------------------------------------*/
} // namespace cla3p
/*-------------------------------------------------*/
