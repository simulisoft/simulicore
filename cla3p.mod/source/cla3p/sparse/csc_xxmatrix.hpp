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

#ifndef CLA3P_CSC_XXMATRIX_HPP_
#define CLA3P_CSC_XXMATRIX_HPP_

/**
 * @file
 */

#include <ostream>
#include <string>

#include "cla3p/generic/matrix_meta.hpp"
#include "cla3p/dense/dns_xxmatrix.hpp"
#include "cla3p/sparse/csc_xxcontainer.hpp"
#include "cla3p/generic/guard.hpp"

/*-------------------------------------------------*/
namespace cla3p { 
/*-------------------------------------------------*/

namespace prm { template <typename T_Int> class PxMatrix; }

/*-------------------------------------------------*/
namespace csc {
/*-------------------------------------------------*/

/**
 * @nosubgrouping 
 * @brief The sparse matrix class (compressed sparse column format).
 */
template <typename T_Int, typename T_Scalar>
class XxMatrix : public MatrixMeta<T_Int>, public XxContainer<T_Int,T_Scalar> {

	private:
		using T_RScalar = typename TypeTraits<T_Scalar>::real_type;

	public:

		//
		// Convertors
		// Move convertors intentionally left as non-explicit
		//
		template <typename T_Virtual>
		explicit XxMatrix(const VirtualExpression<XxMatrix<T_Int,T_Scalar>,T_Virtual>& v) { operator=(v); }
		template <typename T_Virtual>
		XxMatrix(VirtualExpression<XxMatrix<T_Int,T_Scalar>,T_Virtual>&& v) { operator=(std::move(v)); }

		template <typename T_Virtual>
		XxMatrix<T_Int,T_Scalar>& operator=(const VirtualExpression<XxMatrix<T_Int,T_Scalar>,T_Virtual>& v) { evaluateFrom(v); return *this; }
		template <typename T_Virtual>
		XxMatrix<T_Int,T_Scalar>& operator=(VirtualExpression<XxMatrix<T_Int,T_Scalar>,T_Virtual>&& v) { evaluateFrom(v); return *this; }

		VirtualObject<XxMatrix<T_Int,T_Scalar>> virtualize() const { return VirtualObject<XxMatrix<T_Int,T_Scalar>>(*this); }

		/**
		 * @name Constructors
		 * @{
		 */

		/**
		 * @copydoc standard_matrix_docs::constructor()
		 */
		XxMatrix();

		/**
		 * @copydoc standard_csx_docs::dim_constructor()
		 */
		explicit XxMatrix(int_t nr, int_t nc, int_t nz, const Property& pr = Property::General());

		/**
		 * @copydoc standard_csc_docs::aux_constructor()
		 */
		explicit XxMatrix(int_t nr, int_t nc, T_Int *cptr, T_Int *ridx, T_Scalar *vals, bool bind, const Property& pr = Property::General());

		/**
		 * @copydoc standard_docs::copy_constructor()
		 */
		XxMatrix(const XxMatrix<T_Int,T_Scalar>& other);

		/**
		 * @copydoc standard_docs::move_constructor()
		 */
		XxMatrix(XxMatrix<T_Int,T_Scalar>&& other);

		/**
		 * @copydoc standard_matrix_docs::destructor()
		 */
		~XxMatrix();

		/** @} */

		/**
		 * @name Operators
		 * @{
		 */

		/**
		 * @copydoc standard_docs::copy_assignment()
		 */
		XxMatrix<T_Int,T_Scalar>& operator=(const XxMatrix<T_Int,T_Scalar>& other);

		/**
		 * @copydoc standard_docs::move_assignment()
		 */
		XxMatrix<T_Int,T_Scalar>& operator=(XxMatrix<T_Int,T_Scalar>&& other);

		/**
		 * @copydoc standard_docs::virtual_negate_operator()
		 */
		VirtualScale<XxMatrix<T_Int,T_Scalar>,VirtualObject<XxMatrix<T_Int,T_Scalar>>> operator-() const;

		/** @} */

		/** 
		 * @name Arguments
		 * @{
		 */

		/**
		 * @copydoc standard_docs::nnz()
		 */
		int_t nnz() const;

		/** @} */

		/** 
		 * @name Public Member Functions
		 * @{
		 */

		/**
		 * @copydoc standard_docs::clear()
		 */
		void clear();

		/**
		 * @copydoc standard_matrix_docs::info()
		 */
		std::string info(const std::string& header = "") const;

		/**
		 * @copydoc standard_docs::copy()
		 */
		XxMatrix<T_Int,T_Scalar> copy() const;

		/**
		 * @copydoc standard_docs::rcopy()
		 */
		XxMatrix<T_Int,T_Scalar> rcopy();

		/**
		 * @copydoc standard_docs::rcopy_const()
		 */
		Guard<XxMatrix<T_Int,T_Scalar>> rcopy() const;

		/**
		 * @copydoc standard_docs::move()
		 */
		XxMatrix<T_Int,T_Scalar> move();

		/**
		 * @copydoc standard_docs::toStream()
		 */
		void toStream(std::ostream& os, std::streamsize prec = 0) const;

		/**
		 * @copydoc standard_docs::iscale()
		 */
		void iscale(T_Scalar val);

		/**
		 * @copydoc standard_matrix_docs::virtual_transpose()
		 */
		VirtualTranspose<XxMatrix<T_Int,T_Scalar>> transpose() const;

		/**
		 * @copydoc standard_matrix_docs::virtual_ctranspose()
		 */
		VirtualTranspose<XxMatrix<T_Int,T_Scalar>> ctranspose() const;

		/**
		 * @copydoc standard_matrix_docs::virtual_conjugate()
		 */
		VirtualConjugate<XxMatrix<T_Int,T_Scalar>> conjugate() const;

		/**
		 * @copydoc standard_matrix_docs::iconjugate()
		 */
		void iconjugate();

		/**
		 * @copydoc standard_docs::normOne()
		 */
		T_RScalar normOne() const;

		/**
		 * @copydoc standard_docs::normInf()
		 */
		T_RScalar normInf() const;

		/**
		 * @copydoc standard_docs::normMax()
		 */
		T_RScalar normMax() const;

		/**
		 * @copydoc standard_docs::normFro()
		 */
		T_RScalar normFro() const;

		/**
		 * @copydoc standard_matrix_docs::general()
		 */
		XxMatrix<T_Int,T_Scalar> general() const;

		/**
		 * @copydoc standard_matrix_docs::to_dns()
		 */
		dns::XxMatrix<T_Scalar> toDns() const;

		/**
		 * @copydoc standard_matrix_docs::permute_leftright()
		 */
		XxMatrix<T_Int,T_Scalar> permuteLeftRight(const prm::PxMatrix<int_t>& P, const prm::PxMatrix<int_t>& Q) const;

		/**
		 * @copydoc standard_matrix_docs::permute_left()
		 */
		XxMatrix<T_Int,T_Scalar> permuteLeft(const prm::PxMatrix<int_t>& P) const;

		/**
		 * @copydoc standard_matrix_docs::permute_right()
		 */
		XxMatrix<T_Int,T_Scalar> permuteRight(const prm::PxMatrix<int_t>& Q) const;

		/**
		 * @copydoc standard_matrix_docs::permute_mirror()
		 */
		XxMatrix<T_Int,T_Scalar> permuteMirror(const prm::PxMatrix<int_t>& P) const;

		/**
		 * @copydoc standard_matrix_docs::block()
		 */
		XxMatrix<T_Int,T_Scalar> block(int_t ibgn, int_t jbgn, int_t ni, int_t nj) const;

		/** @} */

		/** 
		 * @name Creators/Generators
		 * @{
		 */

		/**
		 * @copydoc standard_csx_docs::random()
		 */
		static XxMatrix<T_Int,T_Scalar> random(int_t nr, int_t nc, int_t nz, const Property& pr = Property::General(),
				T_RScalar lo = T_RScalar(0), T_RScalar hi = T_RScalar(1));

		/**
		 * @copydoc standard_csc_docs::view()
		 */
		static Guard<XxMatrix<T_Int,T_Scalar>> view(int_t nr, int_t nc, const T_Int *cptr, const T_Int *ridx, const T_Scalar *vals, const Property& pr = Property::General());

		/** @} */

	private:
		void copyFromExisting(const XxMatrix<T_Int,T_Scalar>& other);
		void moveFrom(XxMatrix<T_Int,T_Scalar>& other);
		void checker() const;

	protected:
		template <typename T_Virtual>
		void evaluateFrom(const VirtualExpression<XxMatrix<T_Int,T_Scalar>,T_Virtual>& v)
		{
			if(*this) {
				v.evaluateOnExisting(*this);
			} else {
				v.evaluateOnNew(*this);
			}
		}

};

/*-------------------------------------------------*/
} // namespace csc
} // namespace cla3p
/*-------------------------------------------------*/

/**
 * @ingroup cla3p_module_index_stream_operators
 * @brief Writes to os the contents of mat.
 */
template <typename T_Int, typename T_Scalar>
std::ostream& operator<<(std::ostream& os, const cla3p::csc::XxMatrix<T_Int,T_Scalar>& mat)
{
	mat.toStream(os);
	return os;
}

#endif // CLA3P_CSC_XXMATRIX_HPP_
