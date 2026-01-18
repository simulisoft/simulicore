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

#ifndef CULITE_DNS_XXMATRIX_HPP_
#define CULITE_DNS_XXMATRIX_HPP_

/**
 * @file
 */

//#include <string>

#include <cla3p/generic/guard.hpp>
#include <cla3p/generic/matrix_meta.hpp>
#include <cla3p/dense/dns_xxmatrix.hpp>

#include "culite/dense/dns_xxcontainer.hpp"
#include "culite/dense/dns_xxvector.hpp"

/*-------------------------------------------------*/
namespace culite { 
namespace dns {
/*-------------------------------------------------*/

/**
 * @nosubgrouping 
 * @brief The device dense matrix class.
 * @details Represents a dense matrix stored on the GPU device.
 * @tparam T_Scalar The scalar type (e.g., float, double, complex).
 */
template <typename T_Scalar>
class XxMatrix : public ::cla3p::MatrixMeta<int_t>, public XxContainer<T_Scalar> {

	private:
		using T_RScalar = typename TypeTraits<T_Scalar>::real_type;
		using T_Cla3pScalar = typename TypeTraits<T_Scalar>::cla3p_type;

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
		 * @copydoc standard_matrix_docs::dim_constructor()
		 */
		explicit XxMatrix(int_t nr, int_t nc, const ::cla3p::Property& pr = ::cla3p::Property::General());

		/**
		 * @copydoc standard_matrix_docs::aux_constructor()
		 */
		explicit XxMatrix(int_t nr, int_t nc, T_Scalar *vals, int_t ldv, bool bind, const ::cla3p::Property& pr = ::cla3p::Property::General());

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
		 * @copydoc standard_docs::virtual_negate_operator()
		 */
		XxMatrix<T_Scalar> operator-() const; // TODO: use virtuals

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
		::cla3p::Guard<XxMatrix<T_Scalar>> rcopy() const;

		/**
		 * @copydoc standard_docs::move()
		 */
		XxMatrix<T_Scalar> move();

		/**
		 * @copydoc standard_matrix_docs::info()
		 */
		std::string info(const std::string& header = "") const;

		/**
		 * @copydoc standard_docs::iscale()
		 */
		void iscale(T_Scalar val);

		/**
		 *  @copydoc standard_matrix_docs::virtual_transpose()
		 */
		XxMatrix<T_Scalar> transpose() const; // TODO: use virtuals

		/**
		 * @copydoc standard_matrix_docs::virtual_ctranspose()
		 */
		XxMatrix<T_Scalar> ctranspose() const; // TODO: use virtuals

		/**
		 * @copydoc standard_matrix_docs::virtual_conjugate()
		 */
		XxMatrix<T_Scalar> conjugate() const; // TODO: use virtuals

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
		::cla3p::Guard<XxMatrix<T_Scalar>> rblock(int_t ibgn, int_t jbgn, int_t ni, int_t nj) const;

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
		::cla3p::Guard<XxVector<T_Scalar>> rcolumn(int_t j) const;

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
		::cla3p::Guard<XxMatrix<T_Scalar>> rrow(int_t i) const;

		/* TODO: use virtuals
		 * @copydoc standard_matrix_docs::rrowvec()
		 */
		//VirtualRowvec<T_Scalar> rrowvec(int_t i) const;

		/**
		 * @brief Copies the device matrix to a host matrix.
		 * @details Performs a device-to-host copy operation.
		 *          If @p dest is empty, it is resized to match the dimensions and properties of the device matrix.
		 *          Otherwise, @p dest must be compatible with the device matrix.
		 * @param[out] dest The destination host matrix.
		 */
		void copyToHost(::cla3p::dns::XxMatrix<T_Cla3pScalar>& dest) const;

		/**
		 * @brief Copies a host matrix to the device matrix.
		 * @details Performs a host-to-device copy operation.
		 *          If the device matrix is empty, it is resized to match the dimensions and properties of @p src.
		 *          Otherwise, the device matrix must be compatible with @p src.
		 * @param[in] src The source host matrix.
		 */
		void copyFromHost(const ::cla3p::dns::XxMatrix<T_Cla3pScalar>& src);

		/** @} */

		/** 
		 * @name Creators/Generators
		 * @{
		 */

		/**
		 * @copydoc standard_matrix_docs::view()
		 */
		static ::cla3p::Guard<XxMatrix<T_Scalar>> view(int_t nr, int_t nc, const T_Scalar *vals, int_t ldv, const ::cla3p::Property& pr = ::cla3p::Property::General());

		/** @} */

	private:
		int_t m_ld;

		void defaults();

		void setLd(int_t ld);

		void moveFrom(XxMatrix<T_Scalar>& other);
		void copyFromExisting(const XxMatrix<T_Scalar>& other);
		void checker() const;
};

/*-------------------------------------------------*/
} // namespace dns
} // namespace culite
/*-------------------------------------------------*/

/**
 * @ingroup culite_module_index_device2host_operators
 * @brief Copies the contents of src to dest.
 * @details Performs a device-to-host copy from @p src to @p dest.
 *          If @p dest is empty, it is resized to match @p src, otherwise it must be compatible with @p src.
 * @tparam T_Scalar The scalar type (e.g., float, double, complex).
 * @param[in] src The source device matrix.
 * @param[out] dest The destination host matrix.
 */
template <typename T_Scalar>
void operator>>(const culite::dns::XxMatrix<T_Scalar>& src,
	            ::cla3p::dns::XxMatrix<typename culite::TypeTraits<T_Scalar>::cla3p_type>& dest)
{
	src.copyToHost(dest);
}

/*-------------------------------------------------*/

/**
 * @ingroup culite_module_index_host2device_operators
 * @brief Copies the contents of src to dest.
 * @details Performs a host-to-device copy from @p src to @p dest.
 *          If @p dest is empty, it is resized to match @p src, otherwise it must be compatible with @p src.
 * @tparam T_Scalar The scalar type (e.g., float, double, complex).
 * @param[in] src The source host matrix.
 * @param[out] dest The destination device matrix.
 */
template <typename T_Scalar>
void operator>>(const ::cla3p::dns::XxMatrix<typename culite::TypeTraits<T_Scalar>::cla3p_type>& src,
	            culite::dns::XxMatrix<T_Scalar>& dest)
{
	dest.copyFromHost(src);
}

/*-------------------------------------------------*/

/**
 * @ingroup culite_module_index_stream_operators
 * @brief Writes to os the contents of mat.
 * @details Copies the device matrix to host and writes it to the output stream.
 * @tparam T_Scalar The scalar type (e.g., float, double, complex).
 * @param[in,out] os The output stream.
 * @param[in] mat The device matrix to output.
 * @return The output stream.
 */
template <typename T_Scalar>
std::ostream& operator<<(std::ostream& os, const culite::dns::XxMatrix<T_Scalar>& mat)
{
	::cla3p::dns::XxMatrix<typename culite::TypeTraits<T_Scalar>::cla3p_type> hostMat;
	mat >> hostMat;
	os << hostMat;
	return os;
}

/*-------------------------------------------------*/

#endif // CULITE_DNS_XXMATRIX_HPP_
