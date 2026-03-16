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

#ifndef CULITE_ENUMS_HPP_
#define CULITE_ENUMS_HPP_

/** 
 * @file
 */

#include <iostream>

#include <cla3p/types/enums.hpp>

/*-------------------------------------------------*/
namespace culite {
/*-------------------------------------------------*/

/**
 * @ingroup culite_module_index_forwards
 * @brief Matrix algebraic property enumeration.
 * @details Alias of @ref cla3p::prop_t into the `culite` namespace.
 *          Defines the algebraic properties of matrices used to optimize computational operations.
 *          Each property type imposes specific structural constraints and enables specialized algorithms.
 */
using prop_t = ::cla3p::prop_t;

/**
 * @ingroup culite_module_index_forwards
 * @brief Matrix fill type enumeration.
 * @details Alias of @ref cla3p::uplo_t into the `culite` namespace.
 *          Specifies which triangular portion of a matrix is stored and utilized.
 *          This enumeration is used to define storage efficiency for symmetric, Hermitian,
 *          and triangular matrices where only one triangular portion needs to be stored.
 */
using uplo_t = ::cla3p::uplo_t;

/**
 * @ingroup culite_module_index_forwards
 * @brief Virtual matrix operation enumeration.
 * @details Alias of @ref cla3p::op_t into the `culite` namespace.
 *          Specifies virtual operations to be applied to matrices in computational routines.
 *          These operations are not explicitly performed; instead, the operation flag is passed
 *          to optimized algorithms that implicitly handle the transformation.
 */
using op_t = ::cla3p::op_t;

/**
 * @ingroup culite_module_index_forwards
 * @brief Diagonal element type enumeration.
 * @details Alias of @ref cla3p::diag_t into the `culite` namespace.
 *          Specifies whether the diagonal of a triangular matrix consists of unit elements
 *          or arbitrary values. This affects storage requirements and computational algorithms.
 */
using diag_t = ::cla3p::diag_t;

/**
 * @ingroup culite_module_index_forwards
 * @brief Matrix operand positioning enumeration.
 * @details Alias of @ref cla3p::side_t into the `culite` namespace.
 *          Specifies the position (left or right) of a particular matrix operand
 *          in binary matrix operations such as triangular solves and matrix multiplication.
 */
using side_t = ::cla3p::side_t;

/**
 * @ingroup culite_module_index_forwards
 * @brief Matrix decomposition method enumeration.
 * @details Alias of @ref cla3p::decomp_t into the `culite` namespace.
 *          Specifies the factorization algorithm to be used for matrix decomposition.
 *          The choice of decomposition method depends on matrix properties and computational requirements.
 */
using decomp_t = ::cla3p::decomp_t;

/**
 * @ingroup culite_module_index_forwards
 * @brief Singular vector computation policy enumeration.
 * @details Alias of @ref cla3p::svdPolicy_t into the `culite` namespace.
 *          Specifies the extent of singular vector computation in Singular Value Decomposition (SVD).
 *          Different policies trade off computational cost against the completeness of the decomposition.
 */
using svdPolicy_t = ::cla3p::svdPolicy_t;

/**
 * @ingroup culite_module_index_forwards
 * @brief Eigenvalue computation range enumeration.
 * @details Alias of @ref cla3p::eigRange_t into the `culite` namespace.
 *          Specifies which subset of eigenvalues to compute in eigenvalue decomposition routines.
 *          This allows selective computation to improve performance when only specific eigenvalues
 *          are needed, rather than computing the full spectrum.
 */
using eigRange_t = ::cla3p::eigRange_t;

/*-------------------------------------------------*/

/**
 * @ingroup culite_module_index_datatypes
 * @enum alloc_t
 * @brief Specifies the memory allocation strategy and physical residency.
 * @details This enum is used to categorize pointers by their origin, allowing the
 *          runtime to resolve the correct deallocation function (e.g., @c cudaFree vs
 *          @c cudaFreeHost) and optimize data transfer directions.
 */
enum class alloc_t {
    Device       = 0, /**< Device memory */
    Pinned          , /**< Pinned memory */
    Managed         , /**< Managed memory */
    Unregistered      /**< Unknown origin memory */
};

} // namespace culite

/*-------------------------------------------------*/

/**
 * @ingroup culite_module_index_stream_operators
 * @brief Stream insertion operator for @ref culite::alloc_t enumeration.
 * @details Writes a string representation of the memory allocation type to the output stream.
 * @param[in,out] os The output stream to write to.
 * @param[in] alloc_type The allocation type to be written.
 * @return Reference to the output stream @p os.
 */
std::ostream& operator<<(std::ostream& os, const culite::alloc_t& alloc_type);

/*-------------------------------------------------*/

#endif // CULITE_ENUMS_HPP_
