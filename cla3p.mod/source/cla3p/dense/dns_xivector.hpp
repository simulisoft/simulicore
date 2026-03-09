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

#ifndef CLA3P_DNS_XIVECTOR_HPP_
#define CLA3P_DNS_XIVECTOR_HPP_

/**
 * @file
 */

#include <string>
#include <ostream>

#include <cla3p/generic/guard.hpp>
#include <cla3p/generic/meta1d.hpp>
#include <cla3p/dense/dns_xxcontainer.hpp>

/*-------------------------------------------------*/

namespace cla3p { namespace prm { template <typename T_Int> class PxMatrix; } }

/*-------------------------------------------------*/
namespace cla3p { 
namespace dns {
/*-------------------------------------------------*/

/**
 * @nosubgrouping 
 * @brief The general purpose dense vector class.
 * @details Represents a dense vector stored in contiguous memory.
 * @tparam T_Scalar The scalar type (e.g., float, double, complex).
 */
template <typename T_Scalar>
class XiVector : public Meta1D<int_t>, public XxContainer<T_Scalar> {

    public:
        XiVector();
        explicit XiVector(int_t n);
        explicit XiVector(int_t n, T_Scalar *vals, bool bind);
        ~XiVector();

        XiVector(const XiVector<T_Scalar>& other);
        XiVector<T_Scalar>& operator=(const XiVector<T_Scalar>& other);

        XiVector(XiVector<T_Scalar>&& other);
        XiVector<T_Scalar>& operator=(XiVector<T_Scalar>&& other);

        /** 
         * @name Operators
         * @{
         */

        /**
         * @brief Element access operator.
         * @details Accesses the element at the specified index.
         * @param[in] i The index of the element (0-based).
         * @return Reference to the element at index i.
         */
        T_Scalar& operator()(int_t i);

        /**
         * @brief Element access operator.
         * @details Accesses the element at the specified index.
         * @param[in] i The index of the element (0-based).
         * @return Const reference to the element at index i.
         */
        const T_Scalar& operator()(int_t i) const;

        /** @} */

        /** 
         * @name Public Member Functions
         * @{
         */

        /**
         * @brief Clear the vector.
         * @details Releases all memory and resets the vector to an empty state.
         */
        void clear();

        /**
         * @brief Fill all elements with a value.
         * @details Sets all elements in the vector to the specified scalar value.
         * @param[in] val The scalar value to fill with.
         */
        void fill(T_Scalar val);

        /**
         * @brief Create a deep copy of the vector.
         * @details Creates a new vector with its own memory allocation
         *          and copies all elements from this vector to the new vector.
         * @return A new vector containing a copy of this vector's data.
         */
        XiVector<T_Scalar> copy() const;

        /**
         * @brief Create a reference copy (shallow copy) of the vector.
         * @details Creates a new vector object that references the same memory
         *          as this vector. Changes to either vector will affect both.
         * @return A vector that shares memory with this vector.
         */
        XiVector<T_Scalar> rcopy();

        /**
         * @brief Create a guarded reference copy (shallow copy) of the vector.
         * @details Creates a guarded vector object that references the same memory
         *          as this vector. The guard ensures the reference is read-only.
         * @return A guarded vector that shares memory with this vector.
         */
        Guard<XiVector<T_Scalar>> rcopy() const;

        /**
         * @brief Move the vector's resources.
         * @details Transfers ownership of the memory to a new vector object,
         *          leaving this vector in an empty state.
         * @return A new vector containing this vector's memory.
         */
        XiVector<T_Scalar> move();

        /**
         * @brief Get information about the vector.
         * @details Returns a string containing information about the vector's dimensions.
         * @param[in] header Optional header string to prepend to the information.
         * @return A string containing vector information.
         */
        std::string info(const std::string& header = "") const;

        /**
         * @brief Print the vector to standard output.
         * @details Prints the vector elements to standard output with specified precision.
         * @param[in] prec The output precision (default: 0 for default precision).
         */
        void print(std::streamsize prec = 0) const;

        /**
         * @brief Output the vector to a stream.
         * @details Writes the vector elements to the specified output stream with specified precision.
         * @param[in,out] os The output stream to write to.
         * @param[in] prec The output precision (default: 0 for default precision).
         */
        void toStream(std::ostream& os, std::streamsize prec = 0) const;

        /** @} */

        /** 
         * @name Creators/Generators
         * @{
         */

        /**
         * @brief Create a view of existing memory.
         * @details Creates a guarded vector that references existing memory
         *          without taking ownership. The memory must remain valid for the lifetime
         *          of the returned view.
         * @param[in] n The number of elements.
         * @param[in] vals Pointer to the memory.
         * @return A guarded vector that views the specified memory.
         */
        static Guard<XiVector<T_Scalar>> view(int_t n, const T_Scalar *vals);

        /** @} */

    private:
        XiVector<T_Scalar>& moveFrom(XiVector<T_Scalar>& other);
        XiVector<T_Scalar>& copyFromExisting(const XiVector<T_Scalar>& other);
        void checker() const;
};

/*-------------------------------------------------*/
} // namespace dns
} // namespace cla3p
/*-------------------------------------------------*/

/**
 * @ingroup cla3p_module_index_stream_operators
 * @brief Writes to os the contents of vec.
 * @details Outputs the vector to the stream in a human-readable format.
 * @tparam T_Scalar The scalar type (e.g., float, double, complex).
 * @param[in,out] os The output stream.
 * @param[in] vec The vector to output.
 * @return The output stream.
 */
template <typename T_Scalar>
std::ostream& operator<<(std::ostream& os, const cla3p::dns::XiVector<T_Scalar>& vec)
{
    vec.toStream(os);
    return os;
}

/*-------------------------------------------------*/

#endif // CLA3P_DNS_XIVECTOR_HPP_
