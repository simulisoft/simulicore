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

#ifndef CULITE_INTEGER_HPP_
#define CULITE_INTEGER_HPP_

/** 
 * @file
 */

#include <cstdint>
#include <library_types.h> // cuda
#include <cla3p/types/integer.hpp>

#include <culite/types/traits.hpp>

/*-------------------------------------------------*/
namespace culite {
/*-------------------------------------------------*/

#if defined (CULITE_I64)
using int_t = int64_t;
using uint_t = uint64_t;
const cudaDataType __cuda_itype = cudaDataType::CUDA_R_64I;
const cudaDataType __cuda_utype = cudaDataType::CUDA_R_64U;
#else
using int_t = int;
using uint_t = unsigned int;
const cudaDataType __cuda_itype = cudaDataType::CUDA_R_32I;
const cudaDataType __cuda_utype = cudaDataType::CUDA_R_32U;
#endif

template<> class TypeTraits<int_t> {
    public:
        using cla3p_type = ::cla3p::int_t;
        static std::string type_name();
        static std::string prec_name();
        static constexpr cudaDataType cuda_type() { return __cuda_itype; }
};

template<> class TypeTraits<uint_t> {
    public:
        using cla3p_type = ::cla3p::uint_t;
        static std::string type_name();
        static std::string prec_name();
        static constexpr cudaDataType cuda_utype() { return __cuda_utype; }
};

/*-------------------------------------------------*/
} // namespace culite
/*-------------------------------------------------*/

#endif // CULITE_INTEGER_HPP_
