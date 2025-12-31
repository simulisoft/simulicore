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

#ifndef CLA3P_LAPACK_LU_HPP_
#define CLA3P_LAPACK_LU_HPP_

/**
 * @file
 */

#include "cla3p/linsol/lapack_base.hpp"

/*-------------------------------------------------*/
namespace cla3p { 
/*-------------------------------------------------*/

namespace dns { template <typename T_Scalar> class XxVector; }

/**
 * @nosubgrouping
 * @brief The partial pivoting LU linear solver for dense matrices.
 */
template <typename T_Matrix>
class LapackLU : public LapackBase<T_Matrix> {

	using T_Scalar = typename T_Matrix::value_type;
	using T_Vector = dns::XxVector<T_Scalar>;

	public:

		// no copy
		LapackLU(const LapackLU&) = delete;
		LapackLU& operator=(const LapackLU&) = delete;

		/**
		 * @copydoc cla3p::LapackLLt::LapackLLt()
		 */
		LapackLU() : LapackBase<T_Matrix>(decomp_t::LU) {}

		/**
		 * @copydoc cla3p::LapackLLt::LapackLLt(int_t n)
		 */
		LapackLU(int_t n) : LapackBase<T_Matrix>(decomp_t::LU, n) {}

		/**
		 * @copydoc cla3p::LapackLLt::~LapackLLt()
		 */
		~LapackLU() = default;
};

/*-------------------------------------------------*/
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_LAPACK_LU_HPP_
