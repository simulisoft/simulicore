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

/*-------------------------------------------------*/
namespace culite {
/*-------------------------------------------------*/



	int test(); // TODO: DELETEME



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

#if 0
template<> class TypeTraits<real_t> {
	public:
		using real_type = real_t;
		using complex_type = complex_t;
		static std::string type_name();
		static std::string prec_name();
		constexpr static char netlibChar() { return 'd'; }
		constexpr static bool is_real() { return true; }
		constexpr static bool is_complex() { return false; }
		constexpr static bool is_single_precision() { return false; }
		constexpr static bool is_double_precision() { return true; }
		static real_type epsilon();
};

template<> class TypeTraits<real4_t> {
	public:
		using real_type = real4_t;
		using complex_type = complex8_t;
		static std::string type_name();
		static std::string prec_name();
		constexpr static char netlibChar() { return 's'; }
		constexpr static bool is_real() { return true; }
		constexpr static bool is_complex() { return false; }
		constexpr static bool is_single_precision() { return true; }
		constexpr static bool is_double_precision() { return false; }
		static real_type epsilon();
};

template<> class TypeTraits<complex_t> {
	public:
		using real_type = complex_t::value_type;
		using complex_type = complex_t;
		static std::string type_name();
		static std::string prec_name();
		constexpr static char netlibChar() { return 'z'; }
		constexpr static bool is_real() { return false; }
		constexpr static bool is_complex() { return true; }
		constexpr static bool is_single_precision() { return false; }
		constexpr static bool is_double_precision() { return true; }
		static real_type epsilon();
};

template<> class TypeTraits<complex8_t> {
	public:
		using real_type = complex8_t::value_type;
		using complex_type = complex8_t;
		static std::string type_name();
		static std::string prec_name();
		constexpr static char netlibChar() { return 'c'; }
		constexpr static bool is_real() { return false; }
		constexpr static bool is_complex() { return true; }
		constexpr static bool is_single_precision() { return true; }
		constexpr static bool is_double_precision() { return false; }
		static real_type epsilon();
};

/*-------------------------------------------------*/

namespace arith {
inline real_t conj(const real_t& d) { return d; }
inline real4_t conj(const real4_t& s) { return s; }
inline complex_t conj(const complex_t& z) { return std::conj(z); }
inline complex8_t conj(const complex8_t& c) { return std::conj(c); }

inline real_t  getRe(const real_t& d) { return d; }
inline real4_t getRe(const real4_t& f) { return f; }
inline real_t  getRe(const complex_t& z) { return z.real(); }
inline real4_t getRe(const complex8_t& c) { return c.real(); }

inline real_t  getIm(const real_t&) { return 0; }
inline real4_t getIm(const real4_t&) { return 0; }
inline real_t  getIm(const complex_t& z) { return z.imag(); }
inline real4_t getIm(const complex8_t& c) { return c.imag(); }

inline void setIm(real_t&, real_t) { }
inline void setIm(real4_t&, real4_t) { }
inline void setIm(complex_t& z, real_t d) { z.imag(d); }
inline void setIm(complex8_t& c, real4_t s) { c.imag(s); }

template <typename T_Scalar>
inline T_Scalar inv(T_Scalar val)
{
	return (T_Scalar(1) / val); // skip zero check
}

} // namespace arith
#endif // 0

/*-------------------------------------------------*/
} // namespace culite
/*-------------------------------------------------*/

#endif // CULITE_SCALAR_HPP_
