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

#ifndef CLA3P_PROPERTY_HPP_
#define CLA3P_PROPERTY_HPP_

/** 
 * @file
 */

#include <ostream>
#include <string>
#include <vector>

#include "cla3p/types/scalar.hpp"
#include "cla3p/types/enums.hpp"

/*-------------------------------------------------*/
namespace cla3p {
/*-------------------------------------------------*/

/**
 * @brief The property class.
 * @details Used to define matrix properties.@n
 *          Constructing a matrix with property can accelerate certain operations, 
 *          but at the same time imposes equivalent constraints on others.
 * @include ex00_property_create.cpp
 */
class Property {

	public:

		/**
		 * @brief Default constructor.
		 * @details Constructs an empty property with default values.
		 */
		Property();

		/**
		 * @brief Constructs a property with specified type and fill pattern.
		 * @details Constructs a property object with the specified property type and upper/lower fill type.
		 * @param[in] ptype The property type (general, symmetric, hermitian, triangular, or skew).
		 * @param[in] ftype The matrix fill pattern (upper, lower, or full).
		 */
		explicit Property(prop_t ptype, uplo_t ftype);

		/**
		 * @brief Destructor.
		 * @details Destroys the property object.
		 */
		~Property();

		/**
		 * @brief Copy constructor.
		 * @details Creates a new property object by copying another property object.
		 * @param[in] other The property object to copy.
		 */
		Property(const Property& other);

		/**
		 * @brief Copy assignment operator.
		 * @details Copies the contents of another property object to this property object.
		 * @param[in] other The property object to copy.
		 * @return Reference to this property object.
		 */
		Property& operator=(const Property& other);

		/**
		 * @brief Equality comparison operator.
		 * @details Compares two property objects for equality.
		 * @param[in] other The property object to compare against.
		 * @return @c true if both properties have identical type and fill pattern, @c false otherwise.
		 */
		bool operator==(const Property& other) const;

		/**
		 * @brief Inequality comparison operator.
		 * @details Compares two property objects for inequality.
		 * @param[in] other The property object to compare against.
		 * @return @c true if properties differ in type or fill pattern, @c false otherwise.
		 */
		bool operator!=(const Property& other) const;

		/**
		 * @brief Resets the property to its default state.
		 * @details Clears all property settings and returns the object to its initial default configuration.
		 */
		void clear();

		/**
		 * @brief Retrieves the property type.
		 * @return The property type enumeration value.
		 */
		prop_t type() const;

		/**
		 * @brief Retrieves the matrix fill pattern.
		 * @return The upper/lower fill type enumeration value.
		 */
		uplo_t uplo() const;

		/**
		 * @brief Retrieves the matrix fill pattern as a character.
		 * @return The upper/lower fill type represented as a single character.
		 */
		char cuplo() const;

		/**
		 * @brief Toggles the matrix fill pattern.
		 * @details Switches between upper and lower fill patterns (upper becomes lower and vice versa).
		 *          The full fill pattern remains unchanged.
		 */
		void switchUplo();

		/**
		 * @brief Retrieves the property name.
		 * @return A string representation of the property type.
		 */
		std::string name() const;

		/**
		 * @brief Validates the property type.
		 * @return @c true if the property has a valid type, @c false otherwise.
		 */
		bool isValid() const;

		/**
		 * @brief Determines if the property requires a square matrix.
		 * @return @c true if the property necessitates a square matrix structure, @c false otherwise.
		 */
		bool isSquare() const;

		/**
		 * @brief Determines if the property type is general.
		 * @return @c true if the property represents a general matrix, @c false otherwise.
		 */
		bool isGeneral() const;

		/**
		 * @brief Determines if the property type is symmetric.
		 * @return @c true if the property represents a symmetric matrix, @c false otherwise.
		 */
		bool isSymmetric() const;

		/**
		 * @brief Determines if the property type is Hermitian.
		 * @return @c true if the property represents a Hermitian matrix, @c false otherwise.
		 */
		bool isHermitian() const;

		/**
		 * @brief Determines if the property type is triangular or trapezoidal.
		 * @return @c true if the property represents a triangular or trapezoidal matrix, @c false otherwise.
		 */
		bool isTriangular() const;

		/**
		 * @brief Determines if the property type is skew-symmetric.
		 * @return @c true if the property represents a skew-symmetric matrix, @c false otherwise.
		 */
		bool isSkew() const;

		/**
		 * @brief Determines if the entire matrix is utilized.
		 * @return @c true if both upper and lower matrix parts are used, @c false otherwise.
		 */
		bool isFull() const;

		/**
		 * @brief Determines if the upper triangular portion is utilized.
		 * @return @c true if the upper part of the matrix is used, @c false otherwise.
		 */
		bool isUpper() const;

		/**
		 * @brief Determines if the lower triangular portion is utilized.
		 * @return @c true if the lower part of the matrix is used, @c false otherwise.
		 */
		bool isLower() const;

		/**
		 * @brief Computes the transposed property.
		 * @details Determines the property of the transposed (or conjugate-transposed) matrix.@n
		 *          This operation is applicable to general and triangular matrix types.
		 * @return The property corresponding to the transposed matrix.
		 */
		Property transpose() const;

	public:

		/**
		 * @brief Factory method for general property.
		 * @details Creates a property object representing a general matrix with no special structure.
		 * @return A Property object configured for general matrices.
		 */
		static Property General();

		/**
		 * @brief Factory method for upper-triangular symmetric property.
		 * @details Creates a property object representing a symmetric matrix where the upper triangular part is stored and utilized.
		 * @return A Property object configured for symmetric matrices with upper storage.
		 */
		static Property SymmetricUpper();

		/**
		 * @brief Factory method for lower-triangular symmetric property.
		 * @details Creates a property object representing a symmetric matrix where the lower triangular part is stored and utilized.
		 * @return A Property object configured for symmetric matrices with lower storage.
		 */
		static Property SymmetricLower();

		/**
		 * @brief Factory method for upper-triangular Hermitian property.
		 * @details Creates a property object representing a Hermitian matrix where the upper triangular part is stored and utilized.
		 * @return A Property object configured for Hermitian matrices with upper storage.
		 */
		static Property HermitianUpper();

		/**
		 * @brief Factory method for lower-triangular Hermitian property.
		 * @details Creates a property object representing a Hermitian matrix where the lower triangular part is stored and utilized.
		 * @return A Property object configured for Hermitian matrices with lower storage.
		 */
		static Property HermitianLower();

		/**
		 * @brief Factory method for upper-triangular property.
		 * @details Creates a property object representing an upper triangular or trapezoidal matrix.
		 * @return A Property object configured for upper triangular matrices.
		 */
		static Property TriangularUpper();

		/**
		 * @brief Factory method for lower-triangular property.
		 * @details Creates a property object representing a lower triangular or trapezoidal matrix.
		 * @return A Property object configured for lower triangular matrices.
		 */
		static Property TriangularLower();

		/**
		 * @brief Factory method for upper-triangular skew-symmetric property.
		 * @details Creates a property object representing a skew-symmetric matrix where the upper triangular part is stored and utilized.
		 * @return A Property object configured for skew-symmetric matrices with upper storage.
		 */
		static Property SkewUpper();

		/**
		 * @brief Factory method for lower-triangular skew-symmetric property.
		 * @details Creates a property object representing a skew-symmetric matrix where the lower triangular part is stored and utilized.
		 * @return A Property object configured for skew-symmetric matrices with lower storage.
		 */
		static Property SkewLower();

	private:
		prop_t m_type;
		uplo_t m_uplo;

		void defaults();
		void check() const;
};

/*-------------------------------------------------*/

/**
 * @brief Sanitizes property for scalar type compatibility.
 * @tparam T_Scalar The scalar type of the matrix elements.
 * @param[in] pr The property to sanitize.
 * @return The sanitized property. For real scalar types, Hermitian properties are converted to symmetric properties.
 * @details Ensures property compatibility with the scalar type by converting Hermitian properties to
 *          symmetric properties when applied to real-valued matrices, as Hermitian and symmetric
 *          properties are equivalent for real matrices.
 */
template <typename T_Scalar>
inline Property sanitizeProperty(const Property& pr)
{
	if(pr.isHermitian() && TypeTraits<T_Scalar>::is_real()) {
		return Property(prop_t::Symmetric, pr.uplo());
	} // convert to symmetric
	return pr;
}

/**
 * @brief Retrieves all available property types.
 * @return A vector containing all predefined Property objects.
 */
std::vector<Property> allProperties();

/*-------------------------------------------------*/
} // namespace cla3p
/*-------------------------------------------------*/

/**
 * @ingroup cla3p_module_index_stream_operators
 * @brief Outputs the property information to an output stream.
 * @param[in,out] so The output stream to write to.
 * @param[in] pr The property object to output.
 * @return A reference to the output stream.
 */
std::ostream& operator<<(std::ostream& so, const cla3p::Property& pr);

/*-------------------------------------------------*/

#endif // CLA3P_PROPERTY_HPP_
