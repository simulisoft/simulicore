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

#ifndef CLA3P_DEFAULT_QR_HPP_
#define CLA3P_DEFAULT_QR_HPP_

/**
 * @file
 */

#include "cla3p/types/integer.hpp"
#include "cla3p/types/enums.hpp"
#include "cla3p/support/heap_buffer.hpp"
#include "cla3p/dense/dns_xxvector.hpp"

/*-------------------------------------------------*/
namespace cla3p {
/*-------------------------------------------------*/

/**
 * @nosubgrouping
 * @brief The default singular value decomposition class for dense matrices.
 */
template <typename T_Matrix>
class DefaultQR {

	private:
		using T_Scalar = typename T_Matrix::value_type;
		using T_Vector = dns::XxVector<T_Scalar>;

	public:

		/**
		 * @brief The default constructor.
		 * @details Constructs an empty qr object.
		 */
		DefaultQR();

		/**
		 * @brief The dimensional constructor.
		 * @details Constructs a qr object with internal buffer pre-allocation.
		 * @param[in] m The maximum number of rows the buffer can fit.
		 * @param[in] n The maximum number of columns the buffer can fit.
		 * @param[in] qrPolicy The calculation policy for matrix Q.
		 */
		DefaultQR(int_t m, int_t n, qrPolicy_t qrPolicy = qrPolicy_t::Reflection);

		/**
		 * @brief Destroys the qr object.
		 * @details Clears all internal data and destroys the qr object.
		 */
		~DefaultQR();

		/**
		 * @brief Clears the object internal data.
		 * @details Clears the object internal data and resets all settings.
		 */
		void clear();

		/**
		 * @brief Gets the policy for the Q matrix calculation mode.
		 * @details Informs about the Q matrix calculation options.@n
		 *          Meaning if matrix Q will explicitly be calculated or will be storead as an elementary reflector product.
		 */
		qrPolicy_t getQrPolicy() const;

		/**
		 * @brief Sets the policy for the Q matrix calculation.
		 * @details Informs the object about the Q matrix calculation options.@n
		 *          Meaning if matrix Q will explicitly be calculated or will be storead as an elementary reflector product.
		 */
		void setQrPolicy(qrPolicy_t qrPolicy);

		/**
		 * @brief Allocates internal buffers.
		 * @param[in] m The maximum number of rows the buffer can fit.
		 * @param[in] n The maximum number of columns the buffer can fit.
		 */
		void reserve(int_t m, int_t n);

		/**
		 * @brief Performs matrix QR decomposition.
		 * @param[in] mat The matrix to be decomposed.
		 */
		void decompose(const T_Matrix& mat);

		/**
		 * @brief The matrix R.
		 * @details The (min(m,n) x n) matrix R.
		 */
		const T_Matrix& R() const;

		/**
		 * @brief The matrix Q.
		 * @details The (m x min(m,n)) matrix Q.@n
		 *          Holds for policy qrPolicy_t::Full, else returns a reference to an empty object.
		 */
		const T_Matrix& Q() const;

		/**
		 * @brief The vector tau.
		 * @details The vector tau of size min(m,n).
		 */
		const T_Vector& tau() const;

		const T_Matrix& elementaryReflectors() const;

	private:
		qrPolicy_t m_qrPolicy;

		T_Matrix m_matrixA;
		T_Vector m_tauVector;
		T_Matrix m_matrixR;
		T_Matrix m_matrixQ;

		HeapBuffer<T_Scalar> m_scalarBuffers;

		void defaults();

		void clearInternalWrappers();
		void clearInternalBuffers();
		void resizeInternalObjects(const T_Matrix& mat);

		void fillMatrixR();
		void fillMatrixQ();
};

/*-------------------------------------------------*/
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_DEFAULT_QR_HPP_
