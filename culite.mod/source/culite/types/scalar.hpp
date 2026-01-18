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

#ifndef CULITE_SCALAR_HPP_
#define CULITE_SCALAR_HPP_

/** 
 * @file
 */

#include <cuComplex.h>
#include <library_types.h> // cuda
#include <cla3p/types/scalar.hpp>

#include "culite/types/cuda_macros.hpp"
#include "culite/types/traits.hpp"

/*-------------------------------------------------*/
namespace culite {
/*-------------------------------------------------*/

/**
 * @ingroup culite_module_index_datatypes
 * @brief Single precision real.
 */
using real4_t = float;

/**
 * @ingroup culite_module_index_datatypes
 * @brief Double precision real.
 */
using real_t = double;

/**
 * @ingroup culite_module_index_datatypes
 * @brief Single precision complex.
 */
using complex8_t = cuFloatComplex;

/**
 * @ingroup culite_module_index_datatypes
 * @brief Double precision complex.
 */
using complex_t = cuDoubleComplex;

/*-------------------------------------------------*/

template<> class TypeTraits<real_t> {
	public:
		using real_type = real_t;
		using complex_type = complex_t;
		using host_type = ::cla3p::real_t;
		static std::string type_name();
		static std::string prec_name();
		CULITE_HOST CULITE_DEVICE static host_type toHostType(const real_t v);
		CULITE_HOST CULITE_DEVICE static constexpr cudaDataType cuda_type() { return cudaDataType::CUDA_R_64F; }
		CULITE_HOST CULITE_DEVICE static constexpr char netlibChar() { return 'd'; }
		CULITE_HOST CULITE_DEVICE static constexpr bool is_real() { return true; }
		CULITE_HOST CULITE_DEVICE static constexpr bool is_complex() { return false; }
		CULITE_HOST CULITE_DEVICE static constexpr bool is_single_precision() { return false; }
		CULITE_HOST CULITE_DEVICE static constexpr bool is_double_precision() { return true; }
};

template<> class TypeTraits<real4_t> {
	public:
		using real_type = real4_t;
		using complex_type = complex8_t;
		using host_type = ::cla3p::real4_t;
		static std::string type_name();
		static std::string prec_name();
		CULITE_HOST CULITE_DEVICE static host_type toHostType(const real4_t v);
		CULITE_HOST CULITE_DEVICE static constexpr cudaDataType cuda_type() { return cudaDataType::CUDA_R_32F; }
		CULITE_HOST CULITE_DEVICE static constexpr char netlibChar() { return 's'; }
		CULITE_HOST CULITE_DEVICE static constexpr bool is_real() { return true; }
		CULITE_HOST CULITE_DEVICE static constexpr bool is_complex() { return false; }
		CULITE_HOST CULITE_DEVICE static constexpr bool is_single_precision() { return true; }
		CULITE_HOST CULITE_DEVICE static constexpr bool is_double_precision() { return false; }
};

template<> class TypeTraits<complex_t> {
	public:
		using real_type = real_t;
		using complex_type = complex_t;
		using host_type = ::cla3p::complex_t;
		static std::string type_name();
		static std::string prec_name();
		CULITE_HOST CULITE_DEVICE static host_type toHostType(const complex_t v);
		CULITE_HOST CULITE_DEVICE static constexpr cudaDataType cuda_type() { return cudaDataType::CUDA_C_64F; }
		CULITE_HOST CULITE_DEVICE static constexpr char netlibChar() { return 'z'; }
		CULITE_HOST CULITE_DEVICE static constexpr bool is_real() { return false; }
		CULITE_HOST CULITE_DEVICE static constexpr bool is_complex() { return true; }
		CULITE_HOST CULITE_DEVICE static constexpr bool is_single_precision() { return false; }
		CULITE_HOST CULITE_DEVICE static constexpr bool is_double_precision() { return true; }
};

template<> class TypeTraits<complex8_t> {
	public:
		using real_type = real4_t;
		using complex_type = complex8_t;
		using host_type = ::cla3p::complex8_t;
		static std::string type_name();
		static std::string prec_name();
		CULITE_HOST CULITE_DEVICE static host_type toHostType(const complex8_t v);
		CULITE_HOST CULITE_DEVICE static constexpr cudaDataType cuda_type() { return cudaDataType::CUDA_C_32F; }
		CULITE_HOST CULITE_DEVICE static constexpr char netlibChar() { return 'c'; }
		CULITE_HOST CULITE_DEVICE static constexpr bool is_real() { return false; }
		CULITE_HOST CULITE_DEVICE static constexpr bool is_complex() { return true; }
		CULITE_HOST CULITE_DEVICE static constexpr bool is_single_precision() { return true; }
		CULITE_HOST CULITE_DEVICE static constexpr bool is_double_precision() { return false; }
};

/*-------------------------------------------------*/

/**
 * @ingroup culite_module_index_datatypes
 * @brief Constructs a complex number from the associated real type.
 * @details Constructs a complex number from the associated real type.
 * @param re The real part of the complex number.
 * @param im The imaginary part of the complex number (defaults to 0).
 * @return The complex number initialized with the provided values
 */
CULITE_HOST CULITE_DEVICE inline complex_t
makeComplex(const real_t re, const real_t im = 0) { return make_cuDoubleComplex(re, im); }

/**
 * @ingroup culite_module_index_datatypes
 * @copydoc makeComplex(const real_t, const real_t)
 */
CULITE_HOST CULITE_DEVICE inline complex8_t
makeComplex(const real4_t re, const real4_t im = 0) { return make_cuFloatComplex(re, im); }

/**
 * @ingroup culite_module_index_datatypes
 * @brief Converts a real value to a scalar type.
 * @details Creates a scalar value (real or complex) from the associated real type.
 * @param v The input value of the underlying real type.
 * @return The value converted to T_Scalar.
 */
template <typename T_Scalar>
CULITE_HOST CULITE_DEVICE inline T_Scalar
makeScalar(const typename TypeTraits<T_Scalar>::real_type v) { return v; }

template <>
CULITE_HOST CULITE_DEVICE inline complex_t
makeScalar<complex_t>(const real_t re) { return makeComplex(re); }

template <>
CULITE_HOST CULITE_DEVICE inline complex8_t
makeScalar<complex8_t>(const real4_t re) { return makeComplex(re); }

/*-------------------------------------------------*/

namespace arith {
CULITE_HOST CULITE_DEVICE inline real_t     conj(const real_t d) { return d; }
CULITE_HOST CULITE_DEVICE inline real4_t    conj(const real4_t s) { return s; }
CULITE_HOST CULITE_DEVICE inline complex_t  conj(const complex_t z) { return cuConj(z); }
CULITE_HOST CULITE_DEVICE inline complex8_t conj(const complex8_t c) { return cuConjf(c); }

CULITE_HOST CULITE_DEVICE inline real_t  abs(const real_t d) { return d; }
CULITE_HOST CULITE_DEVICE inline real4_t abs(const real4_t s) { return s; }
CULITE_HOST CULITE_DEVICE inline real_t  abs(const complex_t z) { return cuCabs(z); }
CULITE_HOST CULITE_DEVICE inline real4_t abs(const complex8_t c) { return cuCabsf(c); }

CULITE_HOST CULITE_DEVICE inline real_t  abs2(const real_t d) { return d*d; }
CULITE_HOST CULITE_DEVICE inline real4_t abs2(const real4_t s) { return s*s; }
CULITE_HOST CULITE_DEVICE inline real_t  abs2(const complex_t z) { return cuCreal(z)*cuCreal(z) + cuCimag(z)*cuCimag(z); }
CULITE_HOST CULITE_DEVICE inline real4_t abs2(const complex8_t c) { return cuCrealf(c)*cuCrealf(c) + cuCimagf(c)*cuCimagf(c); }

CULITE_HOST CULITE_DEVICE inline real_t  getRe(const real_t d) { return d; }
CULITE_HOST CULITE_DEVICE inline real4_t getRe(const real4_t f) { return f; }
CULITE_HOST CULITE_DEVICE inline real_t  getRe(const complex_t z) { return cuCreal(z); }
CULITE_HOST CULITE_DEVICE inline real4_t getRe(const complex8_t c) { return cuCrealf(c); }

CULITE_HOST CULITE_DEVICE inline real_t  getIm(const real_t) { return 0; }
CULITE_HOST CULITE_DEVICE inline real4_t getIm(const real4_t) { return 0; }
CULITE_HOST CULITE_DEVICE inline real_t  getIm(const complex_t z) { return cuCimag(z); }
CULITE_HOST CULITE_DEVICE inline real4_t getIm(const complex8_t c) { return cuCimagf(c); }

CULITE_HOST CULITE_DEVICE inline void setIm(real_t*, const real_t) { }
CULITE_HOST CULITE_DEVICE inline void setIm(real4_t*, const real4_t) { }
CULITE_HOST CULITE_DEVICE inline void setIm(complex_t* z, const real_t d) { z->y = d; }
CULITE_HOST CULITE_DEVICE inline void setIm(complex8_t* c, const real4_t s) { c->y = s; }
} // namespace arith

/*-------------------------------------------------*/
} // namespace culite
/*-------------------------------------------------*/

CULITE_HOST CULITE_DEVICE inline culite::complex_t operator*(const culite::complex_t z1, const culite::complex_t z2) { return cuCmul(z1, z2); }
CULITE_HOST CULITE_DEVICE inline culite::complex8_t operator*(const culite::complex8_t c1, const culite::complex8_t c2) { return cuCmulf(c1, c2); }
CULITE_HOST CULITE_DEVICE inline culite::complex_t operator/(const culite::complex_t z1, const culite::complex_t z2) { return cuCdiv(z1, z2); }
CULITE_HOST CULITE_DEVICE inline culite::complex8_t operator/(const culite::complex8_t c1, const culite::complex8_t c2) { return cuCdivf(c1, c2); }

/*-------------------------------------------------*/

#endif // CULITE_SCALAR_HPP_
