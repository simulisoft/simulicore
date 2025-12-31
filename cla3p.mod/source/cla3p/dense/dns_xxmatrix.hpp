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

#ifndef CLA3P_DNS_XXMATRIX_HPP_
#define CLA3P_DNS_XXMATRIX_HPP_

/**
 * @file
 */

#include <string>

#include "cla3p/generic/guard.hpp"
#include "cla3p/generic/matrix_meta.hpp"
#include "cla3p/dense/dns_xxcontainer.hpp"
#include "cla3p/dense/dns_xxvector.hpp"

#include "cla3p/virtuals/virtual_object.hpp"
#include "cla3p/virtuals/virtual_rowvec.hpp"
#include "cla3p/virtuals/virtual_transpose.hpp"
#include "cla3p/virtuals/virtual_conjugate.hpp"
#include "cla3p/virtuals/virtual_scale.hpp"

/*-------------------------------------------------*/
namespace cla3p { 
/*-------------------------------------------------*/

namespace prm { template <typename T_Int> class PxMatrix; }

/*-------------------------------------------------*/
namespace dns {
/*-------------------------------------------------*/

/**
 * @nosubgrouping 
 * @brief The dense matrix class.
 */
template <typename T_Scalar>
class XxMatrix : public MatrixMeta, public XxContainer<T_Scalar> {

	private:
		using T_RScalar = typename TypeTraits<T_Scalar>::real_type;

	public:

		//
		// Convertors
		// Move convertors intentionally left as non-explicit
		//
		template <typename T_Virtual>
		explicit XxMatrix(const VirtualExpression<XxMatrix<T_Scalar>,T_Virtual>& v) { operator=(v); }
		template <typename T_Virtual>
		XxMatrix(VirtualExpression<XxMatrix<T_Scalar>,T_Virtual>&& v) { operator=(std::move(v)); }

		template <typename T_Virtual>
		XxMatrix<T_Scalar>& operator=(const VirtualExpression<XxMatrix<T_Scalar>,T_Virtual>& v) { evaluateFrom(v); return *this; }
		template <typename T_Virtual>
		XxMatrix<T_Scalar>& operator=(VirtualExpression<XxMatrix<T_Scalar>,T_Virtual>&& v) { evaluateFrom(v); return *this; }

		explicit XxMatrix(const VirtualRowvec<T_Scalar>& rv) { operator=(rv); }
		XxMatrix(VirtualRowvec<T_Scalar>&& rv) { operator=(std::move(rv)); }

		XxMatrix<T_Scalar>& operator=(const VirtualRowvec<T_Scalar>& rv) { evaluateFrom(rv); return *this; }
		XxMatrix<T_Scalar>& operator=(VirtualRowvec<T_Scalar>&& rv) { evaluateFrom(rv); return *this;  }

		VirtualObject<XxMatrix<T_Scalar>> virtualize() const { return VirtualObject<XxMatrix<T_Scalar>>(*this); }
	
		/**
		 * @name Constructors
		 * @{
		 */

		/**
		 * @copydoc standard_matrix_docs::constructor()
		 */
		XxMatrix();

		/**
		 * @copydoc standard_matrix_docs::dim_constructor()
		 */
		explicit XxMatrix(int_t nr, int_t nc, const Property& pr = Property::General());

		/**
		 * @copydoc standard_matrix_docs::aux_constructor()
		 */
		explicit XxMatrix(int_t nr, int_t nc, T_Scalar *vals, int_t ldv, bool bind, const Property& pr = Property::General());

		/**
		 * @copydoc standard_docs::copy_constructor()
		 */
		XxMatrix(const XxMatrix<T_Scalar>& other);

		/**
		 * @copydoc standard_docs::move_constructor()
		 */
		XxMatrix(XxMatrix<T_Scalar>&& other);

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
		XxMatrix<T_Scalar>& operator=(const XxMatrix<T_Scalar>& other);

		/**
		 * @copydoc standard_docs::move_assignment()
		 */
		XxMatrix<T_Scalar>& operator=(XxMatrix<T_Scalar>&& other);

		/**
		 * @copydoc standard_docs::index_operator_2d()
		 */
		T_Scalar& operator()(int_t i, int_t j);

		/**
		 * @copydoc standard_docs::index_operator_2d()
		 */
		const T_Scalar& operator()(int_t i, int_t j) const;

		/**
		 * @copydoc standard_docs::virtual_negate_operator()
		 */
		VirtualScale<XxMatrix<T_Scalar>,VirtualObject<XxMatrix<T_Scalar>>> operator-() const;

		/**
		 * @copydoc standard_matrix_docs::fill()
		 */
		void operator=(T_Scalar val);

		/** @} */

		/** 
		 * @name Arguments
		 * @{
		 */

		/**
		 * @copydoc standard_matrix_docs::ld()
		 */
		int_t ld() const;

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
		 * @copydoc standard_matrix_docs::fill()
		 */
		void fill(T_Scalar val);

		/**
		 * @copydoc standard_docs::copy()
		 */
		XxMatrix<T_Scalar> copy() const;

		/**
		 * @copydoc standard_docs::rcopy()
		 */
		XxMatrix<T_Scalar> rcopy();

		/**
		 * @copydoc standard_docs::rcopy_const()
		 */
		Guard<XxMatrix<T_Scalar>> rcopy() const;

		/**
		 * @copydoc standard_docs::move()
		 */
		XxMatrix<T_Scalar> move();

		/**
		 * @copydoc standard_matrix_docs::info()
		 */
		std::string info(const std::string& header = "") const;

		/**
		 * @copydoc standard_docs::print()
		 */
		void print(std::streamsize prec = 0) const;

		/**
		 * @copydoc standard_docs::toStream()
		 */
		void toStream(std::ostream& os, std::streamsize prec = 0) const;

		/**
		 * @copydoc standard_docs::iscale()
		 */
		void iscale(T_Scalar val);

		/**
		 *  @copydoc standard_matrix_docs::virtual_transpose()
		 */
		VirtualTranspose<XxMatrix<T_Scalar>> transpose() const;

		/**
		 * @copydoc standard_matrix_docs::virtual_ctranspose()
		 */
		VirtualTranspose<XxMatrix<T_Scalar>> ctranspose() const;

		/**
		 * @copydoc standard_matrix_docs::virtual_conjugate()
		 */
		VirtualConjugate<XxMatrix<T_Scalar>> conjugate() const;

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
		XxMatrix<T_Scalar> general() const;

		/**
		 * @copydoc standard_matrix_docs::igeneral()
		 */
		void igeneral();

		/**
		 * @copydoc standard_matrix_docs::permute_leftright()
		 */
		XxMatrix<T_Scalar> permuteLeftRight(const prm::PxMatrix<int_t>& P, const prm::PxMatrix<int_t>& Q) const;

		/**
		 * @copydoc standard_matrix_docs::permute_leftright_dst()
		 */
		void permuteLeftRight(const prm::PxMatrix<int_t>& P, const prm::PxMatrix<int_t>& Q, XxMatrix<T_Scalar>& dest) const;

		/**
		 * @copydoc standard_matrix_docs::permute_left()
		 */
		XxMatrix<T_Scalar> permuteLeft(const prm::PxMatrix<int_t>& P) const;

		/**
		 * @copydoc standard_matrix_docs::permute_left_dst()
		 */
		void permuteLeft(const prm::PxMatrix<int_t>& P, XxMatrix<T_Scalar>& dest) const;

		/**
		 * @copydoc standard_matrix_docs::permute_right()
		 */
		XxMatrix<T_Scalar> permuteRight(const prm::PxMatrix<int_t>& Q) const;

		/**
		 * @copydoc standard_matrix_docs::permute_right_dst()
		 */
		void permuteRight(const prm::PxMatrix<int_t>& Q, XxMatrix<T_Scalar>& dest) const;

		/**
		 * @copydoc standard_matrix_docs::permute_mirror()
		 */
		XxMatrix<T_Scalar> permuteMirror(const prm::PxMatrix<int_t>& P) const;

		/**
		 * @copydoc standard_matrix_docs::permute_mirror_dst()
		 */
		void permuteMirror(const prm::PxMatrix<int_t>& P, XxMatrix<T_Scalar>& dest) const;

		/**
		 * @copydoc standard_matrix_docs::block()
		 */
		XxMatrix<T_Scalar> block(int_t ibgn, int_t jbgn, int_t ni, int_t nj) const;

		/**
		 * @copydoc standard_matrix_docs::rblock()
		 */
		XxMatrix<T_Scalar> rblock(int_t ibgn, int_t jbgn, int_t ni, int_t nj);

		/**
		 * @copydoc standard_matrix_docs::rblock()
		 */
		Guard<XxMatrix<T_Scalar>> rblock(int_t ibgn, int_t jbgn, int_t ni, int_t nj) const;

		/**
		 * @copydoc standard_matrix_docs::set_block()
		 */
		void setBlock(int_t ibgn, int_t jbgn, const XxMatrix<T_Scalar>& src);

		/**
		 * @copydoc standard_matrix_docs::column()
		 */
		XxVector<T_Scalar> column(int_t j) const;

		/**
		 * @copydoc standard_matrix_docs::rcolumn()
		 */
		XxVector<T_Scalar> rcolumn(int_t j);

		/**
		 * @copydoc standard_matrix_docs::rcolumn()
		 */
		Guard<XxVector<T_Scalar>> rcolumn(int_t j) const;

		/**
		 * @copydoc standard_matrix_docs::row()
		 */
		XxMatrix<T_Scalar> row(int_t i) const;

		/**
		 * @copydoc standard_matrix_docs::rrow()
		 */
		XxMatrix<T_Scalar> rrow(int_t i);

		/**
		 * @copydoc standard_matrix_docs::rrow()
		 */
		Guard<XxMatrix<T_Scalar>> rrow(int_t i) const;

		/**
		 * @copydoc standard_matrix_docs::rrowvec()
		 */
		VirtualRowvec<T_Scalar> rrowvec(int_t i) const;

		/** @} */

		/** 
		 * @name Creators/Generators
		 * @{
		 */

		/**
		 * @copydoc standard_matrix_docs::random()
		 */
		static XxMatrix<T_Scalar> random(int_t nr, int_t nc, const Property& pr = Property::General(), 
				T_RScalar lo = T_RScalar(0), 
				T_RScalar hi = T_RScalar(1));

		/**
		 * @copydoc standard_matrix_docs::view()
		 */
		static Guard<XxMatrix<T_Scalar>> view(int_t nr, int_t nc, const T_Scalar *vals, int_t ldv, const Property& pr = Property::General());

		/** @} */

	private:
		int_t m_ld;

		void defaults();

		void setLd(int_t ld);

		void moveFrom(XxMatrix<T_Scalar>& other);
		void copyFromExisting(const XxMatrix<T_Scalar>& other);
		void checker() const;

	protected:
		template <typename T_Virtual>
		void evaluateFrom(const VirtualExpression<XxMatrix<T_Scalar>,T_Virtual>& v)
		{
			if(*this) {
				v.evaluateOnExisting(*this);
			} else {
				v.evaluateOnNew(*this);
			}
		}

		void evaluateFrom(const VirtualRowvec<T_Scalar>& rv)
		{
			if(*this) {
				rv.evaluateOnExisting(*this);
			} else {
				rv.evaluateOnNew(*this);
			}
		}
};

/*-------------------------------------------------*/
} // namespace dns
} // namespace cla3p
/*-------------------------------------------------*/

/**
 * @ingroup cla3p_module_index_stream_operators
 * @brief Writes to os the contents of mat.
 */
template <typename T_Scalar>
std::ostream& operator<<(std::ostream& os, const cla3p::dns::XxMatrix<T_Scalar>& mat)
{
	mat.toStream(os);
	return os;
}

#endif // CLA3P_DNS_XXMATRIX_HPP_
