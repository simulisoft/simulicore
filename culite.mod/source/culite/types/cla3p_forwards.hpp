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

#ifndef CULITE_TYPES_CLA3P_FORWARDS_HPP_
#define CULITE_TYPES_CLA3P_FORWARDS_HPP_

/**
 * @file
 */

#include <cla3p/types/enums.hpp>
#include <cla3p/types/property.hpp>

/*-------------------------------------------------*/
namespace culite {
/*-------------------------------------------------*/

/**
 * @ingroup culite_module_index_datatypes
 * @brief The property class.
 * @details Alias of @ref cla3p::Property into the `culite` namespace.
 *          Used to define matrix properties.
 *          Constructing a matrix with property can accelerate certain operations,
 *          but at the same time imposes equivalent constraints on others.
 */
using Property = ::cla3p::Property;

/**
 * @ingroup culite_module_index_datatypes
 * @brief Matrix algebraic property enumeration.
 * @details Alias of @ref cla3p::prop_t into the `culite` namespace.
 *          Defines the algebraic properties of matrices used to optimize computational operations.
 *          Each property type imposes specific structural constraints and enables specialized algorithms.
 */
using prop_t = ::cla3p::prop_t;

/**
 * @ingroup culite_module_index_datatypes
 * @brief Matrix fill type enumeration.
 * @details Alias of @ref cla3p::uplo_t into the `culite` namespace.
 *          Specifies which triangular portion of a matrix is stored and utilized.
 *          This enumeration is used to define storage efficiency for symmetric, Hermitian,
 *          and triangular matrices where only one triangular portion needs to be stored.
 */
using uplo_t = ::cla3p::uplo_t;

/**
 * @ingroup culite_module_index_datatypes
 * @brief Virtual matrix operation enumeration.
 * @details Alias of @ref cla3p::op_t into the `culite` namespace.
 *          Specifies virtual operations to be applied to matrices in computational routines.
 *          These operations are not explicitly performed; instead, the operation flag is passed
 *          to optimized algorithms that implicitly handle the transformation.
 */
using op_t = ::cla3p::op_t;

/**
 * @ingroup culite_module_index_datatypes
 * @brief Diagonal element type enumeration.
 * @details Alias of @ref cla3p::diag_t into the `culite` namespace.
 *          Specifies whether the diagonal of a triangular matrix consists of unit elements
 *          or arbitrary values. This affects storage requirements and computational algorithms.
 */
using diag_t = ::cla3p::diag_t;

/**
 * @ingroup culite_module_index_datatypes
 * @brief Matrix operand positioning enumeration.
 * @details Alias of @ref cla3p::side_t into the `culite` namespace.
 *          Specifies the position (left or right) of a particular matrix operand
 *          in binary matrix operations such as triangular solves and matrix multiplication.
 */
using side_t = ::cla3p::side_t;

/**
 * @ingroup culite_module_index_datatypes
 * @brief Matrix decomposition method enumeration.
 * @details Alias of @ref cla3p::decomp_t into the `culite` namespace.
 *          Specifies the factorization algorithm to be used for matrix decomposition.
 *          The choice of decomposition method depends on matrix properties and computational requirements.
 */
using decomp_t = ::cla3p::decomp_t;

/**
 * @ingroup culite_module_index_datatypes
 * @brief Singular vector computation policy enumeration.
 * @details Alias of @ref cla3p::svdPolicy_t into the `culite` namespace.
 *          Specifies the extent of singular vector computation in Singular Value Decomposition (SVD).
 *          Different policies trade off computational cost against the completeness of the decomposition.
 */
using svdPolicy_t = ::cla3p::svdPolicy_t;

/**
 * @ingroup culite_module_index_datatypes
 * @brief Eigenvalue computation range enumeration.
 * @details Alias of @ref cla3p::eigRange_t into the `culite` namespace.
 *          Specifies which subset of eigenvalues to compute in eigenvalue decomposition routines.
 *          This allows selective computation to improve performance when only specific eigenvalues
 *          are needed, rather than computing the full spectrum.
 */
using eigRange_t = ::cla3p::eigRange_t;

/*-------------------------------------------------*/
} // namespace culite
/*-------------------------------------------------*/

#endif // CULITE_TYPES_CLA3P_FORWARDS_HPP_
