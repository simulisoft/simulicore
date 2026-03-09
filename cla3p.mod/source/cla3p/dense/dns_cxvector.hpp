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

#ifndef CLA3P_DNS_CXVECTOR_HPP_
#define CLA3P_DNS_CXVECTOR_HPP_

/**
 * @file
 */

#include <cla3p/dense/dns_xxvector.hpp>
#include <cla3p/virtuals/virtual_strided.hpp>

/*-------------------------------------------------*/
namespace cla3p { 
namespace dns {
/*-------------------------------------------------*/

/**
 * @nosubgrouping 
 * @brief The dense complex vector class.
 * @details Represents a dense complex-valued vector stored in contiguous memory.
 * @tparam T_Scalar The complex scalar type (e.g., std::complex<float>, std::complex<double>).
 */
template <typename T_Scalar>
class CxVector : public XxVector<T_Scalar> {

    private:
        using T_RScalar = typename TypeTraits<T_Scalar>::real_type;

    public:

        //
        // Convertors
        // Move convertors intentionally left as non-explicit
        //
        explicit CxVector(const XiVector<T_Scalar>& other);
        CxVector<T_Scalar>& operator=(const XiVector<T_Scalar>& other);

        CxVector(XiVector<T_Scalar>&& other);
        CxVector<T_Scalar>& operator=(XiVector<T_Scalar>&& other);
        
        template <typename T_Virtual>
        CxVector(const alias::VirtualExpr_vec<T_Scalar,T_Virtual>& v) : XxVector<T_Scalar>(v) {}
        template <typename T_Virtual>
        CxVector<T_Scalar>& operator=(const alias::VirtualExpr_vec<T_Scalar,T_Virtual>& v) { XxVector<T_Scalar>::operator=(v); return *this; }

        /**
         * @name Constructors
         * @{
         */

        /**
         * @brief Default constructor.
         * @details Creates an empty complex vector with no allocated memory.
         */
        CxVector();

        /**
         * @brief Dimension constructor.
         * @details Creates a complex vector of the specified size and allocates memory.
         * @param[in] n The number of elements.
         */
        explicit CxVector(int_t n);

        /**
         * @brief Auxiliary constructor.
         * @details Creates a complex vector using existing memory.
         * @param[in] n The number of elements.
         * @param[in] vals Pointer to existing memory.
         * @param[in] bind If true, the vector takes ownership of the memory.
         */
        explicit CxVector(int_t n, T_Scalar *vals, bool bind);

        /**
         * @brief Copy constructor.
         * @details Creates a new complex vector by copying another complex vector.
         * @param[in] other The complex vector to copy.
         */
        CxVector(const CxVector<T_Scalar>& other) = default;

        /**
         * @brief Move constructor.
         * @details Creates a new complex vector by moving resources from another complex vector.
         * @param[in] other The complex vector to move from.
         */
        CxVector(CxVector<T_Scalar>&& other) = default;

        /**
         * @brief Destructor.
         * @details Destroys the complex vector and releases allocated memory.
         */
        ~CxVector();

        /** @} */

        /**
         * @name Operators
         * @{
         */

        /**
         * @brief Copy assignment operator.
         * @details Copies the contents of another complex vector to this complex vector.
         * @param[in] other The complex vector to copy.
         * @return Reference to this complex vector.
         */
        CxVector<T_Scalar>& operator=(const CxVector<T_Scalar>& other) = default;

        /**
         * @brief Move assignment operator.
         * @details Moves resources from another complex vector to this complex vector.
         * @param[in] other The complex vector to move from.
         * @return Reference to this complex vector.
         */
        CxVector<T_Scalar>& operator=(CxVector<T_Scalar>&& other) = default;

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
         * @brief Extract the real part (mutable).
         * @details Returns a strided view of the real parts of all complex elements.
         *          Modifications to the returned view will affect this vector.
         * @return A virtual strided vector view of the real components.
         */
        alias::VirtualStrided_vec<T_RScalar> real();

        /**
         * @brief Extract the imaginary part (mutable).
         * @details Returns a strided view of the imaginary parts of all complex elements.
         *          Modifications to the returned view will affect this vector.
         * @return A virtual strided vector view of the imaginary components.
         */
        alias::VirtualStrided_vec<T_RScalar> imag();

        /**
         * @brief Extract the real part (const).
         * @details Returns a guarded strided view of the real parts of all complex elements.
         * @return A guarded strided vector view of the real components.
         */
        alias::GuardedStrided_vec<T_RScalar> real() const;

        /**
         * @brief Extract the imaginary part (const).
         * @details Returns a guarded strided view of the imaginary parts of all complex elements.
         * @return A guarded strided vector view of the imaginary components.
         */
        alias::GuardedStrided_vec<T_RScalar> imag() const;

        /** @} */

};

/*-------------------------------------------------*/
} // namespace dns
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_DNS_CXVECTOR_HPP_
