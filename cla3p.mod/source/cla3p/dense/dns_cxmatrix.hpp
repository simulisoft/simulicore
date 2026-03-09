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

#ifndef CLA3P_DNS_CXMATRIX_HPP_
#define CLA3P_DNS_CXMATRIX_HPP_

/**
 * @file
 */

#include "cla3p/dense/dns_xxmatrix.hpp"
#include "cla3p/virtuals/virtual_strided.hpp"

/*-------------------------------------------------*/
namespace cla3p { 
namespace dns {
/*-------------------------------------------------*/

/**
 * @nosubgrouping 
 * @brief The dense complex matrix class.
 * @details Represents a dense complex-valued matrix stored in column-major format.
 * @tparam T_Scalar The complex scalar type (e.g., std::complex<float>, std::complex<double>).
 */
template <typename T_Scalar>
class CxMatrix : public XxMatrix<T_Scalar> {

    private:
        using T_RScalar = typename TypeTraits<T_Scalar>::real_type;

    public:

        //
        // Convertors
        // Move convertors intentionally left as non-explicit
        //
        explicit CxMatrix(const XxMatrix<T_Scalar>& other);
        CxMatrix<T_Scalar>& operator=(const XxMatrix<T_Scalar>& other);

        CxMatrix(XxMatrix<T_Scalar>&& other);
        CxMatrix<T_Scalar>& operator=(XxMatrix<T_Scalar>&& other);

        template <typename T_Virtual>
        CxMatrix(const alias::VirtualExpr_dns<T_Scalar,T_Virtual>& v) : XxMatrix<T_Scalar>(v) {}
        template <typename T_Virtual>
        CxMatrix<T_Scalar>& operator=(const alias::VirtualExpr_dns<T_Scalar,T_Virtual>& v) { XxMatrix<T_Scalar>::operator=(v); return *this; }

        CxMatrix(const VirtualRowvec<T_Scalar>& rv) : XxMatrix<T_Scalar>(rv) {}
        CxMatrix<T_Scalar>& operator=(const VirtualRowvec<T_Scalar>& rv) { XxMatrix<T_Scalar>::operator=(rv); return *this; }

        /**
         * @name Constructors
         * @{
         */

        /**
         * @brief Default constructor.
         * @details Creates an empty complex matrix with no allocated memory.
         */
        CxMatrix();

        /**
         * @brief Dimension constructor.
         * @details Creates a complex matrix of the specified dimensions and allocates memory.
         * @param[in] nr The number of rows.
         * @param[in] nc The number of columns.
         * @param[in] pr The matrix property (default: General).
         */
        explicit CxMatrix(int_t nr, int_t nc, const Property& pr = Property::General());

        /**
         * @brief Auxiliary constructor.
         * @details Creates a complex matrix using existing memory.
         * @param[in] nr The number of rows.
         * @param[in] nc The number of columns.
         * @param[in] vals Pointer to existing memory.
         * @param[in] ldv The leading dimension.
         * @param[in] bind If true, the matrix takes ownership of the memory.
         * @param[in] pr The matrix property (default: General).
         */
        explicit CxMatrix(int_t nr, int_t nc, T_Scalar *vals, int_t ldv, bool bind, const Property& pr = Property::General());

        /**
         * @brief Copy constructor.
         * @details Creates a new complex matrix by copying another complex matrix.
         * @param[in] other The complex matrix to copy.
         */
        CxMatrix(const CxMatrix<T_Scalar>& other) = default;

        /**
         * @brief Move constructor.
         * @details Creates a new complex matrix by moving resources from another complex matrix.
         * @param[in] other The complex matrix to move from.
         */
        CxMatrix(CxMatrix<T_Scalar>&& other) = default;

        /**
         * @brief Destructor.
         * @details Destroys the complex matrix and releases allocated memory.
         */
        ~CxMatrix();

        /** @} */

        /**
         * @name Operators
         * @{
         */

        /**
         * @brief Copy assignment operator.
         * @details Copies the contents of another complex matrix to this complex matrix.
         * @param[in] other The complex matrix to copy.
         * @return Reference to this complex matrix.
         */
        CxMatrix<T_Scalar>& operator=(const CxMatrix<T_Scalar>& other) = default;

        /**
         * @brief Move assignment operator.
         * @details Moves resources from another complex matrix to this complex matrix.
         * @param[in] other The complex matrix to move from.
         * @return Reference to this complex matrix.
         */
        CxMatrix<T_Scalar>& operator=(CxMatrix<T_Scalar>&& other) = default;

        /**
         * @brief Fill operator.
         * @details Fills all elements with the specified complex value.
         * @param[in] val The complex scalar value to fill with.
         */
        void operator=(T_Scalar val);

        /** @} */

        /** 
         * @name Public Member Functions
         * @{
         */

        /**
         * @brief Accesses the real part of the complex matrix.
         * @details Returns a strided view of the real components that can be used for both reading and writing.
         * @return A virtual strided view of the real part.
         */
        alias::VirtualStrided_dns<T_RScalar> real();

        /**
         * @brief Accesses the imaginary part of the complex matrix.
         * @details Returns a strided view of the imaginary components that can be used for both reading and writing.
         * @return A virtual strided view of the imaginary part.
         */
        alias::VirtualStrided_dns<T_RScalar> imag();

        /**
         * @brief Accesses the real part of the complex matrix (const).
         * @details Returns a read-only strided view of the real components.
         * @return A guarded strided view of the real part.
         */
        alias::GuardedStrided_dns<T_RScalar> real() const;

        /**
         * @brief Accesses the imaginary part of the complex matrix (const).
         * @details Returns a read-only strided view of the imaginary components.
         * @return A guarded strided view of the imaginary part.
         */
        alias::GuardedStrided_dns<T_RScalar> imag() const;

        /** @} */
};

/*-------------------------------------------------*/
} // namespace dns
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_DNS_CXMATRIX_HPP_
