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

#ifndef CLA3P_LAPACK_AUTO_HPP_
#define CLA3P_LAPACK_AUTO_HPP_

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
 * @brief The linear solver for dense matrices with automatic method detection.
 */
template <typename T_Matrix>
class LapackAuto : public LapackBase<T_Matrix> {

	using T_Scalar = typename T_Matrix::value_type;
	using T_Vector = dns::XxVector<T_Scalar>;

	public:

		// no copy
		LapackAuto(const LapackAuto&) = delete;
		LapackAuto& operator=(const LapackAuto&) = delete;

		/**
		 * @copydoc cla3p::LapackLLt::LapackLLt()
		 */
		LapackAuto() : LapackBase<T_Matrix>(decomp_t::Auto) {}

		/**
		 * @copydoc cla3p::LapackLLt::LapackLLt(int_t n)
		 */
		LapackAuto(int_t n) : LapackBase<T_Matrix>(decomp_t::Auto, n) {}

		/**
		 * @copydoc cla3p::LapackLLt::~LapackLLt()
		 */
		~LapackAuto() = default;
};

/*-------------------------------------------------*/
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_LAPACK_AUTO_HPP_
