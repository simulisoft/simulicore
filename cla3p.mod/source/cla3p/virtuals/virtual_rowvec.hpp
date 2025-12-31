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

#ifndef CLA3P_VIRTUAL_ROWVEC_HPP_
#define CLA3P_VIRTUAL_ROWVEC_HPP_

/**
 * @file
 */

#include "cla3p/types/integer.hpp"

/*-------------------------------------------------*/
namespace cla3p { 
/*-------------------------------------------------*/

namespace dns { template <typename T_Scalar> class XxVector; }
namespace dns { template <typename T_Scalar> class XxMatrix; }

/**
 * @nosubgrouping
 * @brief The virtual row-vector class.
 */
template <typename T_Scalar>
class VirtualRowvec {

	public:
		VirtualRowvec(int_t n, const T_Scalar *vals, int_t incv, bool conj); 
		~VirtualRowvec();

		/* TODO: Consider enabling this later
		 * @brief Calculates the inner-product.
		 * @details Calculates the inner-product of `this` virtual row-vector and a given vector.
		 * @param[in] vec The vector to compute the inner-product with.
		 * @returns The inner-product <*this, vec>.
		 */
		T_Scalar evaluateInner(const dns::XxVector<T_Scalar>& vec) const;

		// TODO: doxy comments if needed
		void evaluateOuterOnNew(T_Scalar coeff, const dns::XxVector<T_Scalar>& vec, dns::XxMatrix<T_Scalar>& dest) const;
		void evaluateOuterOnExisting(T_Scalar coeff, const dns::XxVector<T_Scalar>& vec, dns::XxMatrix<T_Scalar>& dest) const;
		void accumulateOuterOnExisting(T_Scalar coeff, const dns::XxVector<T_Scalar>& vec, dns::XxMatrix<T_Scalar>& dest) const;

		void evaluateOnNew(dns::XxMatrix<T_Scalar>& dest) const;
		void evaluateOnExisting(dns::XxMatrix<T_Scalar>& dest) const;
		void accumulateOnExisting(T_Scalar coeff, dns::XxMatrix<T_Scalar>& dest) const;

	private:
		int_t m_size;
		const T_Scalar *m_vals;
		int_t m_incv;
		bool m_conj;
};

/*-------------------------------------------------*/
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_VIRTUAL_ROWVEC_HPP_
