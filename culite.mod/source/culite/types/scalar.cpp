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

// this file inc
#include "culite/types/scalar.hpp"

// system

// 3rd
#include <cla3p/types/literals.hpp>

// culite

/*-------------------------------------------------*/
namespace culite {
/*-------------------------------------------------*/
TypeTraits<real_t>::host_type TypeTraits<real_t>::toHostType(const real_t& v) { return host_type(v); }
std::string TypeTraits<real_t>::type_name() { return ::cla3p::msg::Real(); }
std::string TypeTraits<real_t>::prec_name() { return ::cla3p::msg::Double(); }
/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/
TypeTraits<real4_t>::host_type TypeTraits<real4_t>::toHostType(const real4_t& v) { return host_type(v); }
std::string TypeTraits<real4_t>::type_name() { return ::cla3p::msg::Real(); }
std::string TypeTraits<real4_t>::prec_name() { return ::cla3p::msg::Single(); }
/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/
TypeTraits<complex_t>::host_type TypeTraits<complex_t>::toHostType(const complex_t& v)
{ 
    return host_type(arith::getRe(v), arith::getIm(v));
}
std::string TypeTraits<complex_t>::type_name() { return ::cla3p::msg::Complex(); }
std::string TypeTraits<complex_t>::prec_name() { return ::cla3p::msg::Double(); }
/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/
TypeTraits<complex8_t>::host_type TypeTraits<complex8_t>::toHostType(const complex8_t& v)
{ 
    return host_type(arith::getRe(v), arith::getIm(v));
}
std::string TypeTraits<complex8_t>::type_name() { return ::cla3p::msg::Complex(); }
std::string TypeTraits<complex8_t>::prec_name() { return ::cla3p::msg::Single(); }
/*-------------------------------------------------*/
} // namespace culite
/*-------------------------------------------------*/

