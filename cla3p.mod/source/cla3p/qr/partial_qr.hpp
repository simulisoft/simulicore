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
 * @brief The partial QR decomposition class.
 * @details Utility class for the rank reduction process.
 *          Uses a pivoted QR decomposition to process a chunk of the input matrix
 *          until the stopping criterion is reached. Can be resumed on request, 
 *          or perform a single step. The input matrix is automatically transposed
 *          if needed, in order for the decomposition to be applied on tall matrices
 *          (num rows >= num columns).
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
		 * @brief The default constructor.
		 * @details Constructs an empty partial qr object.
		 */
		PartialQR();

		/**
		 * @brief The dimensional constructor.
		 * @details Constructs a partial qr object with internal buffer pre-allocation.
		 * @param[in] m The maximum number of rows the buffer can fit.
		 * @param[in] n The maximum number of columns the buffer can fit.
		 */
		PartialQR(int_t m, int_t n);

		/**
		 * @brief Destroys the partial qr object.
		 * @details Clears all internal data and destroys the partial qr object.
		 */
		~PartialQR();

		/**
		 * @brief Clears the object internal data.
		 * @details Clears the object internal data and resets all settings.
		 */
		void clear();

		/**
		 * @brief Allocates internal buffers.
		 * @param[in] m The maximum number of rows the buffer can fit.
		 * @param[in] n The maximum number of columns the buffer can fit.
		 */
		void reserve(int_t m, int_t n);

		/**
		 * @brief The cutoff tolerance.
		 * @details Controls the number of rows in R to be kept.
		 *          Rows kept must satisfy |R<sub>i,i</sub>| > |R<sub>0,0</sub>| * cutTolerance().
		 */
		real_t cutTolerance() const;

		/**
		 * @brief Sets the cutoff tolerance.
		 * @details Controls the number of rows in R to be kept.
		 *          Rows kept must satisfy |R<sub>i,i</sub>| > |R<sub>0,0</sub>| * cutTolerance().
		 * @param[in] tol The desired cutoff tolerance.
		 */
		void setCutTolerance(real_t tol);

		/**
		 * @brief Performs matrix partial qr decomposition.
		 * @details Performs QR decomposition until the tolerance criterion |R<sub>i,i</sub>| > |R<sub>0,0</sub>| * cutTolerance() is reached.
		 * @param[in] mat The matrix to be decomposed.
		 */
		void decompose(const T_Matrix& mat);

		/**
		 * @brief Resumes the partial qr decomposition.
		 * @details In case decompose() returned earlier than expected, set a new tolerance value and call 
		 *          resume() that will run until the new tolerance criterion |R<sub>i,i</sub>| > |R<sub>0,0</sub>| * cutTolerance() is reached.@n
		 *          If the ignoreTol option is true, the decomposition will proceed until every column is processed.
		 * @param[in] ignoreTol The decomposition will ignore the tolerance criterion and proceed to full decomposition.
		 */
		void resume(bool ignoreTol = false);

		/**
		 * @brief Performs a single decomposition step.
		 * @details Performs a single decomposition step using a auto-generated block size.
		 * @return The number of columns processed in this step.
		 */
		int_t performDecompositionStep();

		/**
		 * @brief Informs about the number of columns processed.
		 * @details Returns the total number of columns processed so far in the partial QR decomposition process.
		 * @return The number of total columns processed.
		 */
		int_t totalColumnsProcessed() const;

		/**
		 * @brief The matrix R.
		 * @details Returns the matrix R caclulated so far in the partial QR decomposition process.
		 * @return The matrix R in general form.
		 */
		const T_Matrix& R() const;

		/**
		 * @brief The vector tau.
		 * @details The vector tau of size min(m,n).
		 */

		const T_Vector& tau() const;

		/**
		 * @brief The pivoting permutation matrix.
		 * @details Returns the permutation matrix that holds the column pivots of the partial QR decomposition process.
		 * @return The pivoting permutation matrix.
		 */
		const prm::PiMatrix& P() const;

		const T_Matrix& elementaryReflectors() const;
		bool transFlag() const;
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
