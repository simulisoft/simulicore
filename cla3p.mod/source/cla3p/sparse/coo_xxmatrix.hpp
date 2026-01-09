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

#ifndef CLA3P_COO_XXMATRIX_HPP_
#define CLA3P_COO_XXMATRIX_HPP_

/**
 * @file
 */

#include <ostream>
#include <string>
#include <vector>

#include "cla3p/types.hpp"
#include "cla3p/generic/matrix_meta.hpp"
#include "cla3p/generic/tuple.hpp"
#include "cla3p/sparse/csc_xxmatrix.hpp"

/*-------------------------------------------------*/
namespace cla3p { 
namespace coo {
/*-------------------------------------------------*/

/**
 * @nosubgrouping 
 * @brief The sparse matrix class (coordinate format).
 */
template <typename T_Int, typename T_Scalar>
class XxMatrix : public MatrixMeta<T_Int> {

	private:
		using TupleVec = std::vector<Tuple<T_Int,T_Scalar>>;

	public:
		using index_type = T_Int;
		using value_type = T_Scalar;

	public:

		/**
		 * @name Constructors
		 * @{
		 */

		/**
		 * @copydoc standard_matrix_docs::constructor()
		 */
		XxMatrix();

		/**
		 * @copydoc standard_coo_docs::dim_constructor()
		 */
		explicit XxMatrix(int_t nr, int_t nc, const Property& pr = Property::General());

		/**
		 * @copydoc standard_docs::copy_constructor()
		 */
		XxMatrix(const XxMatrix<T_Int,T_Scalar>& other) = default;

		/**
		 * @copydoc standard_docs::move_constructor()
		 */
		XxMatrix(XxMatrix<T_Int,T_Scalar>&& other) = default;

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
		XxMatrix<T_Int,T_Scalar>& operator=(const XxMatrix<T_Int,T_Scalar>& other) = default;

		/**
		 * @copydoc standard_docs::move_assignment()
		 */
		XxMatrix<T_Int,T_Scalar>& operator=(XxMatrix<T_Int,T_Scalar>&& other) = default;

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
		 * @copydoc standard_coo_docs::reserve()
		 */
		void reserve(int_t nz);

		/**
		 * @copydoc standard_coo_docs::insert_tup()
		 */
		void insert(const Tuple<T_Int,T_Scalar>& tuple);

		/**
		 * @copydoc standard_coo_docs::insert_ijv()
		 */
		void insert(T_Int i, T_Int j, T_Scalar v);

		/**
		 * @copydoc standard_matrix_docs::info()
		 */
		std::string info(const std::string& header = "") const;

		/**
		 * @copydoc standard_docs::toStream()
		 */
		void toStream(std::ostream& os, std::streamsize prec = 0) const;

		/**
		 * @copydoc standard_coo_docs::to_csc()
		 */
		csc::XxMatrix<T_Int,T_Scalar> toCsc(dup_t duplicatePolicy = dup_t::Sum) const;

		/** @} */

	private:
		TupleVec m_tuples;

		TupleVec& tupleVec();
		const TupleVec& tupleVec() const;

		void checker() const;
};

/*-------------------------------------------------*/
} // namespace coo
} // namespace cla3p
/*-------------------------------------------------*/

/**
 * @ingroup cla3p_module_index_stream_operators
 * @brief Writes to os the contents of mat.
 */
template <typename T_Int, typename T_Scalar>
std::ostream& operator<<(std::ostream& os, const cla3p::coo::XxMatrix<T_Int,T_Scalar>& mat)
{
	mat.toStream(os);
	return os;
}

#endif // CLA3P_COO_XXMATRIX_HPP_
