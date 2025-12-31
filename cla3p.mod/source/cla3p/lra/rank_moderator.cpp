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
#include "cla3p/lra/rank_moderator.hpp"

// system

// 3rd

// cla3p
#include "cla3p/dense.hpp"
#include "cla3p/svd.hpp"
#include "cla3p/proxies/lapack_proxy.hpp"
#include "cla3p/checks/lapack_checks.hpp"
#include "cla3p/qr/partial_qr.hpp"

/*-------------------------------------------------*/
namespace cla3p {
namespace lra {
/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/
//template <typename T_Matrix>
//RankModerator<T_Matrix>& RankModerator<T_Matrix>::Default()
//{
//	return __defaultRankModerator;
//}
/*-------------------------------------------------*/
template <typename T_Matrix>
RankModerator<T_Matrix>::RankModerator()
{
	defaults();
}
/*-------------------------------------------------*/
template <typename T_Matrix>
RankModerator<T_Matrix>::RankModerator(T_RScalar tol, lraMethod_t method, T_RScalar ref)
	:
		m_cutTolerance(tol),
		m_reductionMethod(method),
		m_maxReference(ref)
{
}
/*-------------------------------------------------*/
template <typename T_Matrix>
RankModerator<T_Matrix>::~RankModerator()
{
	clear();
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void RankModerator<T_Matrix>::defaults()
{
#if 0
	m_cutTolerance = Default().cutTolerance();
	m_reductionMethod = Default().reductionMethod();
	m_maxReference = Default().maxReference();
#else
	m_cutTolerance = T_RScalar(1e-6);
	m_reductionMethod = lraMethod_t::StandardSVD;
	m_maxReference = autoReference();
#endif
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void RankModerator<T_Matrix>::clear()
{
	defaults();
}
/*-------------------------------------------------*/
template <typename T_Matrix>
typename RankModerator<T_Matrix>::T_RScalar RankModerator<T_Matrix>::cutTolerance() const
{
	return m_cutTolerance;
}
/*-------------------------------------------------*/
template <typename T_Matrix>
lraMethod_t RankModerator<T_Matrix>::reductionMethod() const
{
	return m_reductionMethod;
}
/*-------------------------------------------------*/
template <typename T_Matrix>
typename RankModerator<T_Matrix>::T_RScalar RankModerator<T_Matrix>::maxReference() const
{
	return m_maxReference;
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void RankModerator<T_Matrix>::setCutTolerance(T_RScalar tol)
{
	m_cutTolerance = tol;
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void RankModerator<T_Matrix>::setReductionMethod(lraMethod_t method)
{
	m_reductionMethod = method;
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void RankModerator<T_Matrix>::setMaxReference(T_RScalar ref)
{
	m_maxReference = ref;
}
/*-------------------------------------------------*/
template <typename T_Matrix>
int_t RankModerator<T_Matrix>::findRank(const dns::XxVector<T_RScalar>& S) const
{
	if(!S)
		return 0;

	T_RScalar tol = cutTolerance();
	T_RScalar ref = maxReference();
	int_t lbgn = 0;

	if(ref < 0) {
		ref = S(0);
		lbgn = 1;
	} // auto/fixed

	if(!ref)
		return 0;

	T_RScalar absRef = T_RScalar(ref * tol);

	for(int_t l = lbgn; l < S.size(); l++) {
		if(S(l) < absRef) {
			return l;
		}
	} // l

	return S.size();
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void RankModerator<T_Matrix>::accumulateSigma(const dns::XxVector<T_RScalar>& S, T_Matrix& /*A*/, T_Matrix& B) const
{
	for(int_t j = 0; j < B.ncols(); j++) {
		T_Scalar sj = T_Scalar(S(j));
		B.rcolumn(j).iscale(sj);
	} // j
}
/*-------------------------------------------------*/
template <typename T_Matrix>
XxMatrix<T_Matrix> RankModerator<T_Matrix>::reduceUsingDefaultSVD(const T_Matrix& mat, DefaultSVD<T_Matrix>& svd) const
{
	svd.decompose(mat);

	int_t k = findRank(svd.singularValues());

	XxMatrix<T_Matrix> ret(mat.nrows(), mat.ncols(), k);

	const T_Matrix& U = svd.leftSingularVectors();
	const T_Matrix& V = svd.rightSingularVectors();

	ret.A() = U.rblock(0, 0, U.nrows(), k).get();
	ret.B() = V.rblock(0, 0, V.nrows(), k).get();

	accumulateSigma(svd.singularValues(), ret.A(), ret.B());

	return ret;
}
/*-------------------------------------------------*/
template <typename T_Matrix>
static void calculateMatrixA(T_Matrix& A, const T_Matrix& Q, const T_Matrix& Uk, const dns::XxVector<typename T_Matrix::value_type>& tau)
{
	int_t m = Q.nrows();
	int_t l = Uk.nrows();
	int_t k = Uk.ncols();

	A = 0;
	A.rblock(0,0,l,k) = Uk;

	int_t info = lapack::xxmqr('L', 'N', m, k, l,
			Q.values(), Q.ld(),
			tau.values(),
			A.values(), A.ld());

	lapack_info_check(info);
}
/*-------------------------------------------------*/
template <typename T_Matrix>
static void calculateMatrixB(T_Matrix& B, const T_Matrix& Vk, const prm::PiMatrix& P)
{
	Vk.permuteLeft(P, B);
}
/*-------------------------------------------------*/
template <typename T_Matrix>
int_t RankModerator<T_Matrix>::reduceUntilAllValuesFound(const T_Matrix& mat, PartialQR<T_Matrix>& pqr, DefaultSVD<T_Matrix>& svd) const
{
	//
	// 1) Perform standard decomposition using tol
	// 2) Perform an extra step if needed
	// 3) Decompose entire matrix and return
	//
	
	int_t k = 0;

	pqr.decompose(mat);
	svd.decompose(pqr.R());
	k = findRank(svd.singularValues());

	if(k < pqr.R().nrows())
		return k;

	pqr.performDecompositionStep();
	pqr.fillMatrixR();
	svd.decompose(pqr.R());
	k = findRank(svd.singularValues());

	if(k < pqr.R().nrows())
		return k;

	pqr.resume(true);
	svd.decompose(pqr.R());
	k = findRank(svd.singularValues());

	return k;
}
/*-------------------------------------------------*/
template <typename T_Matrix>
XxMatrix<T_Matrix> RankModerator<T_Matrix>::reduceUsingFastQR(const T_Matrix& mat, PartialQR<T_Matrix>& pqr, DefaultSVD<T_Matrix>& svd) const
{
	int_t k = reduceUntilAllValuesFound(mat, pqr, svd);
	
	XxMatrix<T_Matrix> ret;

	if(k) {

		std::cout << pqr.tau();

		const T_Matrix& U = svd.leftSingularVectors();
		const T_Matrix& V = svd.rightSingularVectors();

		T_Matrix Uk = U.rblock(0,0,U.nrows(),k).get().copy();
		T_Matrix Vk = V.rblock(0,0,V.nrows(),k).get().copy();

		int_t m = pqr.elementaryReflectors().nrows();
		int_t n = pqr.elementaryReflectors().ncols();

		if(pqr.transFlag()) {
			ret = XxMatrix<T_Matrix>(n, m, k);
			calculateMatrixA(ret.B(), pqr.elementaryReflectors(), Uk, pqr.tau());
			calculateMatrixB(ret.A(), Vk, pqr.P());
			accumulateSigma(svd.singularValues(), ret.B(), ret.A());
		} else {
			ret = XxMatrix<T_Matrix>(m, n, k);
			calculateMatrixA(ret.A(), pqr.elementaryReflectors(), Uk, pqr.tau());
			calculateMatrixB(ret.B(), Vk, pqr.P());
			accumulateSigma(svd.singularValues(), ret.A(), ret.B());
		} // worked with transpose

	} // k
	
	return ret;
}
/*-------------------------------------------------*/
template <typename T_Matrix>
XxMatrix<T_Matrix> RankModerator<T_Matrix>::reduce(const T_Matrix& mat) const
{
	if(reductionMethod() == lraMethod_t::StandardSVD) {
		int_t m = mat.nrows();
		int_t n = mat.ncols();
		DefaultSVD<T_Matrix> svd(m,n);
		return reduceUsingDefaultSVD(mat, svd);
	}

	if(reductionMethod() == lraMethod_t::FastQR) {
		int_t m = std::max(mat.nrows(),mat.ncols());
		int_t n = std::min(mat.nrows(),mat.ncols());
		int_t k = std::max(int_t(1),m/2);
		DefaultSVD<T_Matrix> svd(k,n);
		PartialQR<T_Matrix> pqr(m,n);
		real_t pqrTolerance = cutTolerance() * std::sqrt(cutTolerance());
		pqr.setCutTolerance(pqrTolerance);
		return reduceUsingFastQR(mat, pqr, svd);
	}

	return XxMatrix<T_Matrix>();
}
/*-------------------------------------------------*/
template class RankModerator<dns::RdMatrix>;
template class RankModerator<dns::RfMatrix>;
template class RankModerator<dns::CdMatrix>;
template class RankModerator<dns::CfMatrix>;
/*-------------------------------------------------*/
} // namespace lra
} // namespace cla3p
/*-------------------------------------------------*/
template <typename T_Matrix>
std::ostream& operator<<(std::ostream& os, const cla3p::lra::RankModerator<T_Matrix>& rmod)
{
	using T_Scalar = typename T_Matrix::value_type;
	using T_RScalar = typename cla3p::TypeTraits<T_Scalar>::real_type;

	os << "====== Rank Reducer ======\n";
	os << "  tolerance......... " << rmod.cutTolerance() << "\n";
	os << "  reference......... " << (rmod.maxReference() < T_RScalar(0) ? "Auto" : std::to_string(rmod.maxReference())) << "\n";
	os << "  reduction method.. " << rmod.reductionMethod() << std::endl;
	return os;
}
/*-------------------------------------------------*/

