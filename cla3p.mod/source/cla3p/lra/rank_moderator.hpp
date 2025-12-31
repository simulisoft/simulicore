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

#ifndef CLA3P_RANK_MODERATOR_HPP_
#define CLA3P_RANK_MODERATOR_HPP_

/**
 * @file
 */

#include <string>
#include <iostream>

#include "cla3p/types/integer.hpp"
#include "cla3p/types/scalar.hpp"
#include "cla3p/types/enums.hpp"
#include "cla3p/generic/type_traits.hpp"
#include "cla3p/lra/lra_xxmatrix.hpp"

/*-------------------------------------------------*/

namespace cla3p { template <typename T_Matrix> class DefaultSVD; }
namespace cla3p { template <typename T_Matrix> class PartialQR; }

namespace cla3p { namespace dns { template <typename T_Scalar> class XxVector; } }

/*-------------------------------------------------*/
namespace cla3p { 
namespace lra {
/*-------------------------------------------------*/

/**
 * @nosubgrouping 
 * @brief The rank moderator class.
 */
template <typename T_Matrix>
class RankModerator {

	private:
		using T_Scalar = typename T_Matrix::value_type;
		using T_RScalar = typename TypeTraits<T_Scalar>::real_type;

	public:

		/**
		 * @name Constructors
		 * @{
		 */

		/**
		 * @copydoc standard_docs::constructor
		 */
		explicit RankModerator();

		/**
		 * @brief The parameterized constructor.
		 * @details Constructs a rank moderator with user defined parameters.
		 * @param[in] tol The cutoff tolerance. Ranks kept must satisfy &sigma;<sub>j</sub> > `ref` * `tol`
		 * @param[in] method The reduction method
		 * @param[in] ref Use an auto/fixed value as reference. @n
		 *                If negative, the maximum singular value will be used instead.
		 */
		explicit RankModerator(T_RScalar tol, lraMethod_t method, T_RScalar ref = autoReference());

		/**
		 * @copydoc standard_docs::destructor
		 */
		~RankModerator();

		/** @} */

		/**
		 * @name Operators
		 * @{
		 */

		/** @} */

		/**
		 * @name Arguments
		 * @{
		 */

		/**
		 * @brief The cutoff tolerance.
		 * @details Controls the level of detail in the low-rank matrix.
		 *          Ranks kept must satisfy &sigma;<sub>j</sub> > maxReference() * cutTolerance().
		 */
		T_RScalar cutTolerance() const;

		/**
		 * @brief The reduction method.
		 * @details Controls the full-rank to low-rank conversion method.
		 */
		lraMethod_t reductionMethod() const;

		/**
		 * @brief The reference value.
		 * @details Controls the level of detail in the low-rank matrix.
		 *          Ranks kept must satisfy &sigma;<sub>j</sub> > maxReference() * cutTolerance().
		 *          If negative, the maximum singular value will be used instead. 
		 */
		T_RScalar maxReference() const;

		/**
		 * @brief Sets the cutoff tolerance.
		 * @details Controls the level of detail in the low-rank matrix.
		 *          Ranks kept must satisfy &sigma;<sub>j</sub> > maxReference() * cutTolerance().
		 * @param[in] tol The desired cutoff tolerance.
		 */
		void setCutTolerance(T_RScalar tol);

		/**
		 * @brief Sets the reduction method.
		 * @details Controls the full-rank to low-rank conversion method.
		 * @param[in] method The desired reduction method.
		 */
		void setReductionMethod(lraMethod_t method);

		/**
		 * @brief The reference value.
		 * @details Controls the level of detail in the low-rank matrix.
		 *          Ranks kept must satisfy &sigma;<sub>j</sub> > maxReference() * cutTolerance().
		 *          If negative, the maximum singular value will be used instead. 
		 * @param[in] ref The desired auto/fixed reference value.
		 */
		void setMaxReference(T_RScalar ref);

		/** @} */

		/**
		 * @name Public Member Functions
		 * @{
		 */

		/**
		 * @brief Clears the rank moderator.
		 * @details Resets all members to match those of the default rank moderator.
		 */
		void clear();

		/**
		 * @brief Full-rank to low-rank convertor.
		 * @details Converts a full rank dense matrix to low-rank using internal parameterization.
		 * @param[in] mat The dense matrix to be reduced.
		 * @return The low-rank approximation of `mat`.
		 */
		XxMatrix<T_Matrix> reduce(const T_Matrix& mat) const;

		/** @} */

		/**
		 * @name Advanced Public Member Functions
		 * @{
		 */

		/**
		 * @brief Full-rank to low-rank convertor.
		 * @details Converts a full rank dense matrix to low-rank using internal the standard SVD method.
		 *          A user parameterized svd object is supplied in order to avoid repeated (de-)allocations.
		 * @param[in] mat The dense matrix to be reduced.
		 * @param[in] svd The Singular Value Decomposition object to be used.
		 * @return The low-rank approximation of `mat`.
		 */
		XxMatrix<T_Matrix> reduceUsingDefaultSVD(const T_Matrix& mat, DefaultSVD<T_Matrix>& svd) const;

		/**
		 * @brief Full-rank to low-rank convertor.
		 * @details Converts a full rank dense matrix to low-rank using internal the hybrid QR/SVD method.
		 *          User parameterized pqr/svd objects are supplied in order to avoid repeated (de-)allocations.
		 * @param[in] mat The dense matrix to be reduced.
		 * @param[in] pqr The Partial QR Decomposition object to be used.
		 * @param[in] svd The Singular Value Decomposition object to be used.
		 * @return The low-rank approximation of `mat`.
		 */
		XxMatrix<T_Matrix> reduceUsingFastQR(const T_Matrix& mat, PartialQR<T_Matrix>& pqr, DefaultSVD<T_Matrix>& svd) const;

		/** @} */

		/**
		 * @name Global Options
		 * @{
		 */

		/**
		 * @brief Reference value will be set automaticaly.
		 * @details The default value for automatic deduction of the reference value (&sigma;<sub>max</sub>).
		 */
		static constexpr T_RScalar autoReference() { return T_RScalar(-1); }

		/* TODO: maybe enable this later
		 * @brief The default rank moderator.
		 * @details Gets the global rank moderator.
		 *          The global moderator is fully customizable.
		 *          Every rank moderator will inherit the global moderator settings on instantiation.
		 */
		//static RankModerator<T_Matrix>& Default();

		/** @} */

	private:

		T_RScalar m_cutTolerance;
		lraMethod_t m_reductionMethod;
		T_RScalar m_maxReference;

		void defaults();

		int_t findRank(const dns::XxVector<T_RScalar>& S) const;
		void accumulateSigma(const dns::XxVector<T_RScalar>& S, T_Matrix& A, T_Matrix& B) const;
		int_t reduceUntilAllValuesFound(const T_Matrix& mat, PartialQR<T_Matrix>& pqr, DefaultSVD<T_Matrix>& svd) const;
};

/*-------------------------------------------------*/
} // namespace lra
} // namespace cla3p
/*-------------------------------------------------*/

/*
 @ingroup module_index_stream_operators
 @brief Writes to os the contents of rmod.
*/
template <typename T_Matrix>
std::ostream& operator<<(std::ostream& os, const cla3p::lra::RankModerator<T_Matrix>& rmod);

#endif // CLA3P_RANK_MODERATOR_HPP_
