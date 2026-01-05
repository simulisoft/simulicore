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

#ifndef CLA3P_PXMATRIX_HPP_
#define CLA3P_PXMATRIX_HPP_

/** 
 * @file
 */

#include <string>

#include "cla3p/types/literals.hpp"
#include "cla3p/generic/type_traits.hpp"
#include "cla3p/dense/dns_xivector.hpp"

/*-------------------------------------------------*/
namespace cla3p {
namespace prm {
/*-------------------------------------------------*/

/**
 * @nosubgrouping
 * @brief Permutation matrix utility.
 * @tparam T_Int Integer type for permutation indices.
 * @details Permutation matrices are stored as a 1D array of (unsigned) integers
 *          and reorder indices based on a predefined mapping. For a size @p n
 *          matrix, each value lies in the range [0, @p n-1] using 0-based indexing.
 */
template <typename T_Int>
class PxMatrix : public dns::XiVector<T_Int> {

	public:

		/** 
		 * @name Constructors
		 * @{
		 */

		/**
		 * @brief Default constructor.
		 * @details Constructs an empty permutation matrix.
		 */
		explicit PxMatrix();

		/**
		 * @brief Size constructor.
		 * @details Constructs a permutation matrix of size @p n with uninitialized values.
		 * @param[in] n The permutation matrix size.
		 */
		explicit PxMatrix(int_t n);

		/**
		 * @brief Buffer binding constructor.
		 * @details Creates a permutation matrix of size @p n bound to external
		 *          storage or copies from provided values.
		 * @param[in] n The permutation matrix size.
		 * @param[in] vals Pointer to @p n entries holding the permutation values.
		 * @param[in] bind If @c true, binds to @p vals without ownership;
		 *                 otherwise copies the values.
		 */
		explicit PxMatrix(int_t n, T_Int *vals, bool bind);

		/**
		 * @brief Copy constructor.
		 * @details Constructs a permutation matrix with a copy of @p other;
		 *          @p other remains unchanged.
		 */
		PxMatrix(const PxMatrix<T_Int>& other) = default;

		/**
		 * @brief Move constructor.
		 * @details Transfers the contents of @p other; @p other is left empty.
		 */
		PxMatrix(PxMatrix<T_Int>&& other) = default;

		/**
		 * @brief Destructor.
		 * @details Destroys the permutation matrix.
		 */
		~PxMatrix();

		/** @} */

		/** 
		 * @name Operators
		 * @{
		 */

		/**
		 * @brief Copy assignment operator.
		 * @details If @p (*this) is empty, copies the contents of @p other. If @p (*this)
		 *          is not empty, performs a deep copy of @p other into @p (*this). The two
		 *          matrices must have equal size.
		 */
		PxMatrix& operator=(const PxMatrix<T_Int>& other) = default;

		/**
		 * @brief Move assignment operator.
		 * @details Replaces the contents with those of @p other; @p other is left empty.
		 */
		PxMatrix<T_Int>& operator=(PxMatrix<T_Int>&& other) = default;

		/**
		 * @brief Value setter.
		 * @details Sets all entries of @p (*this) to @p val.
		 * @param[in] val The value to assign to all entries.
		 */
		void operator=(T_Int val);

		/** @} */

		/** 
		 * @name Public Member Functions
		 * @{
		 */

		/**
		 * @brief Inverse permutation.
		 * @details Computes and returns the inverse (transpose) of the permutation matrix.
		 * @return The inverse permutation matrix.
		 */
		PxMatrix<T_Int> inverse() const;

		/**
		 * @brief Left permutation composition.
		 * @details Computes and returns the composed permutation @f$ P \cdot (*this) @f$.
		 * @param[in] P The left-side permutation matrix.
		 * @return The composed permutation @f$ P \cdot (*this) @f$.
		 */
		PxMatrix<T_Int> permuteLeft(const PxMatrix<T_Int>& P) const;

		/**
		 * @brief In-place left permutation composition.
		 * @details Stores the composed permutation @f$ P \cdot (*this) @f$ into @p trg.
		 * @param[in] P The left-side permutation matrix.
		 * @param[out] trg Output permutation matrix receiving @f$ P \cdot (*this) @f$;
		 *          must match the size of @p (*this).
		 */
		void permuteLeft(const PxMatrix<T_Int>& P, PxMatrix<T_Int>& trg) const;

		/** @} */

		/** 
		 * @name Creators/Generators
		 * @{
		 */

		/**
		 * @brief Identity permutation generator.
		 * @details Creates an @p n-sized permutation matrix with @f$ P(i) = i @f$.
		 * @param[in] n The permutation matrix size.
		 * @return The identity permutation matrix.
		 */
		static PxMatrix<T_Int> identity(int_t n);

		/**
		 * @brief Random permutation generator.
		 * @details Creates an @p n-sized permutation matrix with randomly
		 *          rearranged indices.
		 * @param[in] n The permutation matrix size.
		 * @return The random permutation matrix.
		 */
		static PxMatrix<T_Int> random(int_t n);

		/** @} */
};

/*-------------------------------------------------*/
} // namespace prm
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_PXMATRIX_HPP_
