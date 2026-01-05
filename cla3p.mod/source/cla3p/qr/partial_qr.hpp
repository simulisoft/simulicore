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

#ifndef CLA3P_PARTIAL_QR_HPP_
#define CLA3P_PARTIAL_QR_HPP_

/**
 * @file
 */

#include "cla3p/types/integer.hpp"
#include "cla3p/types/scalar.hpp"
#include "cla3p/types/enums.hpp"
#include "cla3p/support/heap_buffer.hpp"
#include "cla3p/dense/dns_xxvector.hpp"
#include "cla3p/perms.hpp"

/*-------------------------------------------------*/
namespace cla3p { 
/*-------------------------------------------------*/

/**
 * @nosubgrouping
 * @brief Partial Pivoted QR Decomposition for rank reduction algorithms.
 * @tparam T_Matrix The dense matrix type to decompose.
 * @details This class implements a partial (stopped) QR decomposition with column pivoting,
 *          useful for low-rank approximation and rank reduction processes.
 *          The algorithm processes the matrix incrementally, stopping when singular values
 *          fall below a specified tolerance threshold. The input matrix is automatically
 *          transposed as needed to ensure tall-matrix orientation (@f$ m \geq n @f$).
 */
template <typename T_Matrix>
class PartialQR {

	private:
		using T_Scalar = typename T_Matrix::value_type;
		using T_RScalar = typename TypeTraits<T_Scalar>::real_type;
		using T_Vector = dns::XxVector<T_Scalar>;
		using T_RVector = dns::XxVector<T_RScalar>;

	public:

		/**
		 * @brief Constructs an uninitialized partial QR object.
		 * @details Initializes an empty partial QR object with no allocated memory for decomposition results.
		 */
		PartialQR();

		/**
		 * @brief Constructs a partial QR object with pre-allocated buffers.
		 * @details Initializes the partial QR decomposition object and allocates internal buffers
		 *          to accommodate matrices up to the specified dimensions.
		 * @param[in] m The maximum number of matrix rows.
		 * @param[in] n The maximum number of matrix columns.
		 */
		PartialQR(int_t m, int_t n);

		/**
		 * @brief Destroys the partial QR decomposition object.
		 * @details Releases all allocated memory and clears internal data structures.
		 */
		~PartialQR();

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
		 * @brief Retrieves the singular value tolerance threshold.
		 * @details Returns the relative tolerance parameter used to determine which columns are kept
		 *          in the partial QR decomposition. Diagonal elements @f$ R_{i,i} @f$ are retained
		 *          if @f$ |R_{i,i}| > |R_{0,0}| \cdot \text{cutTolerance}() @f$.
		 * @return The current cutoff tolerance value.
		 */
		real_t cutTolerance() const;

		/**
		 * @brief Configures the singular value tolerance threshold.
		 * @details Sets the relative tolerance parameter that determines stopping condition.
		 *          Diagonal elements @f$ R_{i,i} @f$ are retained if @f$ |R_{i,i}| > |R_{0,0}| \cdot \text{tol} @f$.
		 * @param[in] tol The desired cutoff tolerance value.
		 */
		void setCutTolerance(real_t tol);

		/**
		 * @brief Initiates the partial QR decomposition.
		 * @details Begins column-wise QR decomposition with column pivoting, stopping when the tolerance criterion
		 *          @f$ |R_{i,i}| > |R_{0,0}| \cdot \text{cutTolerance}() @f$ is violated.
		 * @param[in] mat The matrix to decompose.
		 */
		void decompose(const T_Matrix& mat);

		/**
		 * @brief Resumes the partial QR decomposition.
		 * @details Continues the decomposition from where @ref decompose() or the previous @ref resume() stopped.
		 *          The decomposition proceeds until the tolerance criterion @f$ |R_{i,i}| > |R_{0,0}| \cdot \text{cutTolerance}() @f$ is violated,
		 *          or all columns are processed if @p ignoreTol is `true`.
		 * @param[in] ignoreTol If `true`, ignores the tolerance criterion and decomposes all remaining columns.
		 */
		void resume(bool ignoreTol = false);

		/**
		 * @brief Performs a single decomposition step.
		 * @details Executes one decomposition step with an automatically determined block size.
		 * @return The number of columns processed in this decomposition step.
		 */
		int_t performDecompositionStep();

		/**
		 * @brief Retrieves the count of processed columns.
		 * @details Returns the total number of columns processed so far throughout the partial QR decomposition.
		 * @return The total number of columns processed.
		 */
		int_t totalColumnsProcessed() const;

		/**
		 * @brief Retrieves the upper triangular matrix.
		 * @details Returns the matrix @f$ R @f$ calculated so far in the partial QR decomposition process.
		 * @return A constant reference to the upper triangular matrix @f$ R @f$ in general form.
		 */
		const T_Matrix& R() const;

		/**
		 * @brief Retrieves the Householder reflection coefficients.
		 * @details Returns the vector @f$ \tau @f$ of size @f$ \min(m,n) @f$ containing
		 *          the scalar factors for elementary Householder reflectors.
		 * @return A constant reference to the Householder scalar coefficient vector.
		 */
		const T_Vector& tau() const;

		/**
		 * @brief Retrieves the column pivoting permutation matrix.
		 * @details Returns the permutation matrix @f$ P @f$ that encodes the column pivots applied during
		 *          the partial QR decomposition process.
		 * @return A constant reference to the column pivoting permutation matrix.
		 */
		const prm::PiMatrix& P() const;

		/**
		 * @brief Retrieves the elementary Householder reflector matrix.
		 * @details Returns a matrix containing the elementary Householder vectors used to implicitly represent
		 *          the orthogonal transformation matrix from the QR decomposition.
		 * @return A constant reference to the elementary reflector matrix.
		 */
		const T_Matrix& elementaryReflectors() const;

		/**
		 * @brief Determines if input matrix was transposed.
		 * @details Indicates whether the input matrix was automatically conjugate-transposed
		 *          during initialization to achieve tall-matrix orientation.
		 * @return `true` if the matrix was transposed, `false` otherwise.
		 */
		bool transFlag() const;

		/**
		 * @brief Extracts the upper triangular portion to form matrix @f$ R @f$.
		 * @details Constructs matrix @f$ R @f$ from the elementary reflector matrix,
		 *          optionally keeping only the first @p numRanks columns.
		 * @param[in] numRanks The number of columns to retain (0 means all columns).
		 */
		void fillMatrixR(int_t numRanks = 0);

	private:
		real_t m_cutTolerance;
		int_t m_blockSize;
		bool m_transFlag;
		int_t m_offset;

		prm::PiMatrix m_matrixP;
		T_Matrix m_matrixA;
		T_Matrix m_matrixF;
		T_Matrix m_matrixR;
		T_Vector m_tauVector;
		T_Vector m_auxVector;
		T_RVector m_vn1Vector;
		T_RVector m_vn2Vector;

		HeapBuffer<T_Scalar> m_scalarBuffer;
		HeapBuffer<T_RScalar> m_realBuffer;
		HeapBuffer<int_t> m_intBuffer;

		void defaults();

		void setOptimumBlockSize(int_t, int_t);

		void clearInternalWrappers();
		void clearInternalBuffers();
		void resizeInternalObjects(const T_Matrix&);

		void initializeDecompositionProcess(const T_Matrix&);
		int_t toleranceCriterionReached(int_t);

		void pivotSwitchBase0();
		void pivotSwitchBase1();

		bool allColumnsProcessed() const;
};

/*-------------------------------------------------*/
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_PARTIAL_QR_HPP_
