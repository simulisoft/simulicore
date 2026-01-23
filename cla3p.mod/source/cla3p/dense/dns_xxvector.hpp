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

#ifndef CLA3P_DNS_XXVECTOR_HPP_
#define CLA3P_DNS_XXVECTOR_HPP_

/**
 * @file
 */

#include <string>
#include <ostream>

#include "cla3p/generic/guard.hpp"
#include "cla3p/dense/dns_xivector.hpp"

#include "cla3p/virtuals/virtual_rowvec.hpp"
#include "cla3p/virtuals/virtual_object.hpp"
#include "cla3p/virtuals/virtual_transpose.hpp"
#include "cla3p/virtuals/virtual_conjugate.hpp"
#include "cla3p/virtuals/virtual_scale.hpp"

/*-------------------------------------------------*/

namespace cla3p { namespace prm { template <typename T_Int> class PxMatrix; } }

/*-------------------------------------------------*/
namespace cla3p { 
namespace dns {
/*-------------------------------------------------*/

/**
 * @nosubgrouping 
 * @brief The floating point dense vector class.
 */
template <typename T_Scalar>
class XxVector : public XiVector<T_Scalar> {

	private:
		using T_RScalar = typename TypeTraits<T_Scalar>::real_type;

	public:

		//
		// Convertors
		// Move convertors intentionally left as non-explicit
		//
		explicit XxVector(const XiVector<T_Scalar>& other);
		XxVector(XiVector<T_Scalar>&& other);

		XxVector<T_Scalar>& operator=(XiVector<T_Scalar>&& other);
		XxVector<T_Scalar>& operator=(const XiVector<T_Scalar>& other);

		template <typename T_Virtual>
		explicit XxVector(const alias::VirtualExpr_vec<T_Scalar,T_Virtual>& v) { operator=(v); }
		template <typename T_Virtual>
		XxVector(alias::VirtualExpr_vec<T_Scalar,T_Virtual>&& v) { operator=(std::move(v)); }

		template <typename T_Virtual>
		XxVector<T_Scalar>& operator=(const alias::VirtualExpr_vec<T_Scalar,T_Virtual>& v) { evaluateFrom(v); return *this; }
		template <typename T_Virtual>
		XxVector<T_Scalar>& operator=(alias::VirtualExpr_vec<T_Scalar,T_Virtual>&& v) { evaluateFrom(v); return *this; }

		alias::VirtualObj_vec<T_Scalar> virtualize() const { return alias::VirtualObj_vec<T_Scalar>(*this); }

		/**
		 * @name Constructors
		 * @{
		 */

		/**
		 * @brief Default constructor.
		 * @details Creates an empty vector with no allocated memory.
		 */
		XxVector();

		/**
		 * @brief Dimension constructor.
		 * @details Creates a vector of the specified size and allocates memory.
		 * @param[in] n The number of elements.
		 */
		explicit XxVector(int_t n);

		/**
		 * @brief Auxiliary constructor.
		 * @details Creates a vector using existing memory.
		 * @param[in] n The number of elements.
		 * @param[in] vals Pointer to existing memory.
		 * @param[in] bind If true, the vector does not take ownership of the memory.
		 */
		explicit XxVector(int_t n, T_Scalar *vals, bool bind);

		/**
		 * @brief Copy constructor.
		 * @details Creates a new vector by copying another vector.
		 * @param[in] other The vector to copy.
		 */
		XxVector(const XxVector<T_Scalar>& other) = default;

		/**
		 * @brief Move constructor.
		 * @details Creates a new vector by moving resources from another vector.
		 * @param[in] other The vector to move from.
		 */
		XxVector(XxVector<T_Scalar>&& other) = default;

		/**
		 * @brief Destructor.
		 * @details Destroys the vector and releases allocated memory.
		 */
		~XxVector();

		/** @} */

		/** 
		 * @name Operators
		 * @{
		 */

		/**
		 * @brief Copy assignment operator.
		 * @details Copies the contents of another vector to this vector.
		 * @param[in] other The vector to copy.
		 * @return Reference to this vector.
		 */
		XxVector<T_Scalar>& operator=(const XxVector<T_Scalar>& other) = default;

		/**
		 * @brief Move assignment operator.
		 * @details Moves resources from another vector to this vector.
		 * @param[in] other The vector to move from.
		 * @return Reference to this vector.
		 */
		XxVector<T_Scalar>& operator=(XxVector<T_Scalar>&& other) = default;

		/**
		 * @brief Fill operator.
		 * @details Fills all elements with the specified value.
		 * @param[in] val The scalar value to fill with.
		 */
		void operator=(T_Scalar val);

		/**
		 * @brief Unary negation operator.
		 * @details Returns a negated copy of the vector.
		 * @return A virtual expression containing the negated elements.
		 */
		alias::VirtualScal_vec<T_Scalar> operator-() const;

		/** @} */

		/** 
		 * @name Public Member Functions
		 * @{
		 */

		/**
		 * @brief Scale the vector in-place.
		 * @details Multiplies all elements of the vector by a scalar value.
		 * @param[in] val The scalar value to multiply by.
		 */
		void iscale(T_Scalar val);

		/**
		 * @brief Transpose the vector.
		 * @details Returns a row vector view of this column vector.
		 * @return A virtual row vector expression.
		 */
		VirtualRowvec<T_Scalar> transpose() const;

		/**
		 * @brief Conjugate transpose the vector.
		 * @details Returns a conjugate transposed row vector view of this column vector.
		 * @return A virtual row vector expression.
		 */
		VirtualRowvec<T_Scalar> ctranspose() const;

		/**
		 * @brief Compute the complex conjugate.
		 * @details Returns a virtual expression containing the complex conjugate of each element.
		 * @return A virtual expression with conjugated elements.
		 */
		alias::VirtualConj_vec<T_Scalar> conjugate() const;

		/**
		 * @brief Conjugate the vector in-place.
		 * @details Replaces all elements with their complex conjugates.
		 */
		void iconjugate();

		/**
		 * @brief Compute the 1-norm.
		 * @details Computes the sum of absolute values of all elements.
		 * @return The 1-norm of the vector.
		 */
		T_RScalar normOne() const;

		/**
		 * @brief Compute the infinity norm.
		 * @details Computes the maximum absolute value of all elements.
		 * @return The infinity norm of the vector.
		 */
		T_RScalar normInf() const;

		/**
		 * @brief Compute the Euclidean norm.
		 * @details Computes the square root of the sum of squared absolute values of all elements.
		 * @return The Euclidean norm of the vector.
		 */
		T_RScalar normEuc() const;

		/**
		 * @brief Permute the vector.
		 * @details Applies a left permutation to the vector.
		 * @param[in] P The permutation matrix.
		 * @return A new vector containing the permuted elements.
		 */
		XxVector<T_Scalar> permuteLeft(const prm::PxMatrix<int_t>& P) const;

		/**
		 * @brief Permute the vector into a destination.
		 * @details Applies a left permutation to the vector and stores the result in the destination vector.
		 * @param[in] P The permutation matrix.
		 * @param[out] dest The destination vector for the permuted elements.
		 */
		void permuteLeft(const prm::PxMatrix<int_t>& P, XxVector<T_Scalar>& dest) const;

		/**
		 * @brief Extract a block as a new vector.
		 * @details Creates a new vector containing a deep copy of a contiguous block of elements.
		 * @param[in] ibgn The starting index of the block.
		 * @param[in] ni The number of elements in the block.
		 * @return A new vector containing the block's data.
		 */
		XxVector<T_Scalar> block(int_t ibgn, int_t ni) const;

		/**
		 * @brief Extract a reference block.
		 * @details Creates a vector that references a contiguous block of this vector's memory.
		 * @param[in] ibgn The starting index of the block.
		 * @param[in] ni The number of elements in the block.
		 * @return A vector that shares memory with this vector.
		 */
		XxVector<T_Scalar> rblock(int_t ibgn, int_t ni);

		/**
		 * @brief Extract a guarded reference block.
		 * @details Creates a guarded vector that references a contiguous block of this vector's memory.
		 * @param[in] ibgn The starting index of the block.
		 * @param[in] ni The number of elements in the block.
		 * @return A guarded vector that shares memory with this vector.
		 */
		Guard<XxVector<T_Scalar>> rblock(int_t ibgn, int_t ni) const;

		/**
		 * @brief Set a block of elements.
		 * @details Copies elements from the source vector to a contiguous block of this vector.
		 * @param[in] ibgn The starting index where the block will be written.
		 * @param[in] src The source vector to copy from.
		 */
		void setBlock(int_t ibgn, const XxVector<T_Scalar>& src);

		/** @} */

		/** 
		 * @name Creators/Generators
		 * @{
		 */

		/**
		 * @brief Create a random vector.
		 * @details Creates a vector with random values uniformly distributed in the specified range.
		 * @param[in] n The number of elements.
		 * @param[in] lo The lower bound of the random values (default: 0).
		 * @param[in] hi The upper bound of the random values (default: 1).
		 * @return A vector containing random values.
		 */
		static XxVector<T_Scalar> random(int_t n, T_RScalar lo = T_RScalar(0), T_RScalar hi = T_RScalar(1));

		/** @} */

	protected:
		template <typename T_Virtual>
		void evaluateFrom(const alias::VirtualExpr_vec<T_Scalar,T_Virtual>& v)
		{
			if(*this) {
				v.evaluateOnExisting(*this);
			} else {
				v.evaluateOnNew(*this);
			}
		}

};

/*-------------------------------------------------*/
} // namespace dns
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_DNS_XXVECTOR_HPP_
