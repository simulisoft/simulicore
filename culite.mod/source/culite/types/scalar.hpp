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
#include <cla3p/types/scalar.hpp>

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
		static host_type toHostType(const real_t& v);
		static std::string type_name();
		static std::string prec_name();
		constexpr static char netlibChar() { return 'd'; }
		constexpr static bool is_real() { return true; }
		constexpr static bool is_complex() { return false; }
		constexpr static bool is_single_precision() { return false; }
		constexpr static bool is_double_precision() { return true; }
};

template<> class TypeTraits<real4_t> {
	public:
		using real_type = real4_t;
		using complex_type = complex8_t;
		using host_type = ::cla3p::real4_t;
		static host_type toHostType(const real4_t& v);
		static std::string type_name();
		static std::string prec_name();
		constexpr static char netlibChar() { return 's'; }
		constexpr static bool is_real() { return true; }
		constexpr static bool is_complex() { return false; }
		constexpr static bool is_single_precision() { return true; }
		constexpr static bool is_double_precision() { return false; }
};

template<> class TypeTraits<complex_t> {
	public:
		using real_type = real_t;
		using complex_type = complex_t;
		using host_type = ::cla3p::complex_t;
		static host_type toHostType(const complex_t& v);
		static std::string type_name();
		static std::string prec_name();
		constexpr static char netlibChar() { return 'z'; }
		constexpr static bool is_real() { return false; }
		constexpr static bool is_complex() { return true; }
		constexpr static bool is_single_precision() { return false; }
		constexpr static bool is_double_precision() { return true; }
};

template<> class TypeTraits<complex8_t> {
	public:
		using real_type = real4_t;
		using complex_type = complex8_t;
		using host_type = ::cla3p::complex8_t;
		static host_type toHostType(const complex8_t& v);
		static std::string type_name();
		static std::string prec_name();
		constexpr static char netlibChar() { return 'c'; }
		constexpr static bool is_real() { return false; }
		constexpr static bool is_complex() { return true; }
		constexpr static bool is_single_precision() { return true; }
		constexpr static bool is_double_precision() { return false; }
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
inline complex_t makeComplex(const real_t& re, const real_t& im = 0) { return make_cuDoubleComplex(re, im); }

/**
 * @ingroup culite_module_index_datatypes
 * @copydoc makeComplex(const real_t&, const real_t&)
 */
inline complex8_t makeComplex(const real4_t& re, const real4_t& im = 0) { return make_cuFloatComplex(re, im); }

/**
 * @ingroup culite_module_index_datatypes
 * @brief Converts a real value to a scalar type.
 * @details Creates a scalar value (real or complex) from the associated real type.
 * @param v The input value of the underlying real type.
 * @return The value converted to T_Scalar.
 */
template <typename T_Scalar>
inline T_Scalar makeScalar(const typename TypeTraits<T_Scalar>::real_type& v) { return v; }

template <>
inline complex_t makeScalar<complex_t>(const real_t& re) { return makeComplex(re); }

template <>
inline complex8_t makeScalar<complex8_t>(const real4_t& re) { return makeComplex(re); }

/*-------------------------------------------------*/

namespace arith {
inline real_t conj(const real_t& d) { return d; }
inline real4_t conj(const real4_t& s) { return s; }
inline complex_t conj(const complex_t& z) { return cuConj(z); }
inline complex8_t conj(const complex8_t& c) { return cuConjf(c); }

inline real_t  getRe(const real_t& d) { return d; }
inline real4_t getRe(const real4_t& f) { return f; }
inline real_t  getRe(const complex_t& z) { return cuCreal(z); }
inline real4_t getRe(const complex8_t& c) { return cuCrealf(c); }

inline real_t  getIm(const real_t&) { return 0; }
inline real4_t getIm(const real4_t&) { return 0; }
inline real_t  getIm(const complex_t& z) { return cuCimag(z); }
inline real4_t getIm(const complex8_t& c) { return cuCimagf(c); }

inline void setIm(real_t&, real_t) { }
inline void setIm(real4_t&, real4_t) { }
inline void setIm(complex_t& z, real_t d) { z.y = d; }
inline void setIm(complex8_t& c, real4_t s) { c.y = s; }
} // namespace arith

/*-------------------------------------------------*/
} // namespace culite
/*-------------------------------------------------*/

#endif // CULITE_SCALAR_HPP_
