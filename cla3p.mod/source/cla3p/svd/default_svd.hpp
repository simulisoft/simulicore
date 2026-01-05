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

#ifndef CLA3P_DEFAULT_SVD_HPP_
#define CLA3P_DEFAULT_SVD_HPP_

/**
 * @file
 */

#include "cla3p/types/integer.hpp"
#include "cla3p/types/enums.hpp"
#include "cla3p/support/heap_buffer.hpp"

/*-------------------------------------------------*/
namespace cla3p {
/*-------------------------------------------------*/

namespace dns{ template <typename T_Scalar> class XxVector; }

/**
 * @nosubgrouping
 * @brief Singular Value Decomposition (SVD) implementation for dense matrices.
 * @tparam T_Matrix The dense matrix type to decompose.
 * @details This class provides a complete Singular Value Decomposition interface,
 *          decomposing a matrix @f$ A @f$ into @f$ A = U \cdot \Sigma \cdot V^H @f$,
 *          where @f$ U @f$ and @f$ V @f$ are unitary matrices and @f$ \Sigma @f$ is diagonal.
 */
template <typename T_Matrix>
class DefaultSVD {

	private:
		using T_Scalar = typename T_Matrix::value_type;
		using T_RScalar = typename TypeTraits<T_Scalar>::real_type;
		using T_RVector = dns::XxVector<T_RScalar>;

	public:

		/**
		 * @brief Constructs an uninitialized SVD decomposition object.
		 * @details Initializes an empty SVD object with no allocated memory for decomposition results.
		 */
		DefaultSVD();

		/**
		 * @brief Constructs an SVD object with pre-allocated buffers.
		 * @details Initializes the SVD decomposition object and allocates internal buffers
		 *          to accommodate matrices up to the specified dimensions.
		 * @param[in] m The maximum number of matrix rows.
		 * @param[in] n The maximum number of matrix columns.
		 * @param[in] leftPolicy The singular vector computation policy for @f$ U @f$.
		 * @param[in] rightPolicy The singular vector computation policy for @f$ V @f$.
		 */
		DefaultSVD(int_t m, int_t n, svdPolicy_t leftPolicy = svdPolicy_t::Limited, svdPolicy_t rightPolicy = svdPolicy_t::Limited);

		/**
		 * @brief Destroys the SVD decomposition object.
		 * @details Releases all allocated memory and clears internal data structures.
		 */
		~DefaultSVD();

		/**
		 * @brief Clears all decomposition results and resets the object.
		 * @details Deallocates decomposition results and returns the object to its default state.
		 */
		void clear();

		/**
		 * @brief Pre-allocates buffers for decomposition.
		 * @details Allocates internal buffers to accommodate matrices with up to @p m rows and @p n columns.
		 *          This avoids memory reallocation during subsequent decompositions.
		 * @param[in] m The maximum number of matrix rows to support.
		 * @param[in] n The maximum number of matrix columns to support.
		 */
		void reserve(int_t m, int_t n);

		/**
		 * @brief Retrieves the left singular vector computation policy.
		 * @details Returns the policy that determines whether the left singular vectors (@f$ U @f$ matrix)
		 *          will be computed in full, limited to @f$ \min(m,n) @f$ columns, or not at all.
		 * @return The current left singular vector computation policy.
		 */
		svdPolicy_t getLeftPolicy() const;

		/**
		 * @brief Retrieves the right singular vector computation policy.
		 * @details Returns the policy that determines whether the right singular vectors (@f$ V^H @f$ matrix)
		 *          will be computed in full, limited to @f$ \min(m,n) @f$ rows, or not at all.
		 * @return The current right singular vector computation policy.
		 */
		svdPolicy_t getRightPolicy() const;

		/**
		 * @brief Configures the left singular vector computation policy.
		 * @details Specifies whether the left singular vectors (@f$ U @f$ matrix) should be computed
		 *          in full, limited to @f$ \min(m,n) @f$ columns, or not computed at all.
		 * @param[in] leftPolicy The left singular vector computation policy.
		 */
		void setLeftPolicy(svdPolicy_t leftPolicy);

		/**
		 * @brief Configures the right singular vector computation policy.
		 * @details Specifies whether the right singular vectors (@f$ V^H @f$ matrix) should be computed
		 *          in full, limited to @f$ \min(m,n) @f$ rows, or not computed at all.
		 * @param[in] rightPolicy The right singular vector computation policy.
		 */
		void setRightPolicy(svdPolicy_t rightPolicy);

		/**
		 * @brief Performs Singular Value Decomposition on the input matrix.
		 * @details Computes the decomposition @f$ A = U \cdot \Sigma \cdot V^H @f$ of the input matrix.
		 *          The computation respects the singular vector policies set via @ref setLeftPolicy() and @ref setRightPolicy().
		 * @param[in] mat The matrix to decompose.
		 */
		void decompose(const T_Matrix& mat);

		/**
		 * @brief Retrieves the computed singular values.
		 * @details Returns the diagonal elements of @f$ \Sigma @f$ in descending order.
		 * @return A constant reference to a vector of @f$ \min(m,n) @f$ singular values.
		 */
		const T_RVector& singularValues() const;

		/**
		 * @brief Retrieves the left singular vectors matrix.
		 * @details Returns the matrix @f$ U @f$ of size @f$ m \times k @f$, where:
		 *          - @f$ k = m @f$ if left policy is Full
		 *          - @f$ k = \min(m,n) @f$ if left policy is Limited
		 *          - The matrix is empty if left policy is NoCalculation
		 * @return A constant reference to the left singular vectors matrix.
		 */
		const T_Matrix& leftSingularVectors() const;

		/**
		 * @brief Retrieves the right singular vectors matrix.
		 * @details Returns the matrix @f$ V^H @f$ of size @f$ k \times n @f$, where:
		 *          - @f$ k = n @f$ if right policy is Full
		 *          - @f$ k = \min(m,n) @f$ if right policy is Limited
		 *          - The matrix is empty if right policy is NoCalculation
		 * @return A constant reference to the right singular vectors matrix.
		 */
		const T_Matrix& rightSingularVectors() const;

	private:
		svdPolicy_t m_leftPolicy;
		svdPolicy_t m_rightPolicy;

		T_RVector m_singularValues;
		T_RVector m_superbVector;

		T_Matrix m_matrixBackup;
		T_Matrix m_leftSingularVectors;
		T_Matrix m_rightSingularVectors;

		HeapBuffer<T_RScalar> m_realBuffers;
		HeapBuffer<T_Scalar> m_scalarBuffers;

		void defaults();

		void clearInternalWrappers();
		void clearInternalBuffers();
		void resizeInternalObjects(const T_Matrix& mat);
};

/*-------------------------------------------------*/
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_DEFAULT_SVD_HPP_
