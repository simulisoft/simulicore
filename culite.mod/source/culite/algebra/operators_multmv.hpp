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

#ifndef CULITE_OPERATORS_MULTMV_HPP_
#define CULITE_OPERATORS_MULTMV_HPP_

/**
 * @file
 */

#include "culite/dense/dns_xxvector.hpp"
#include "culite/dense/dns_xxmatrix.hpp"
#include "culite/sparse/csr_xxmatrix.hpp"
#include "culite/sparse/csc_xxmatrix.hpp"

#include "culite/virtuals/virtual_product.hpp"

/*-------------------------------------------------*/

/*
 * Virtual x XxVector
 */
template <typename T_Left>
culite::VirtualProduct<
    culite::dns::XxVector<typename T_Left::result_type::value_type>,
    T_Left,
    culite::alias::VirtualObj_vec<typename T_Left::result_type::value_type>> 
operator*(
    const culite::VirtualExpression<typename T_Left::result_type, T_Left>& left, 
    const culite::dns::XxVector<typename T_Left::result_type::value_type>& x) 
{ 
    return culite::VirtualProduct<
        culite::dns::XxVector<typename T_Left::result_type::value_type>,
        T_Left,
        culite::alias::VirtualObj_vec<typename T_Left::result_type::value_type>>(left.self(), x.virtualize());
}

/*
 * XxMatrix x XxVector
 */

/**
 * @ingroup culite_module_index_math_operators_mult
 * @brief Matrix-vector multiplication for dense matrices.
 * @details Performs the operation @f$ A \cdot x @f$.
 * @param[in] A The input dense matrix.
 * @param[in] x The input vector.
 * @return The virtual product representing the multiplication.
 */
template <typename T_Scalar>
culite::alias::VirtualProd_dnsmv<T_Scalar>
operator*(
    const culite::dns::XxMatrix<T_Scalar>& A, 
    const culite::dns::XxVector<T_Scalar>& x) 
{ 
    return culite::alias::VirtualProd_dnsmv<T_Scalar>(A.virtualize(), x.virtualize());
}

/**
 * @ingroup culite_module_index_math_operators_mult
 * @brief Matrix-vector multiplication for sparse matrices.
 * @details Performs the operation @f$ A \cdot x @f$.
 * @param[in] A The input sparse matrix.
 * @param[in] x The input vector.
 * @return The virtual product representing the multiplication.
 */
template <typename T_Int, typename T_Scalar>
culite::alias::VirtualProd_csrmv<T_Int,T_Scalar>
operator*(
    const culite::csr::XxMatrix<T_Int,T_Scalar>& A, 
    const culite::dns::XxVector<T_Scalar>& x) 
{ 
    return culite::alias::VirtualProd_csrmv<T_Int,T_Scalar>(A.virtualize(), x.virtualize());
}

/**
 * @ingroup culite_module_index_math_operators_mult
 * @brief Matrix-vector multiplication for sparse matrices.
 * @details Performs the operation @f$ A \cdot x @f$.
 * @param[in] A The input sparse matrix.
 * @param[in] x The input vector.
 * @return The virtual product representing the multiplication.
 */
template <typename T_Int, typename T_Scalar>
culite::alias::VirtualProd_cscmv<T_Int,T_Scalar>
operator*(
    const culite::csc::XxMatrix<T_Int,T_Scalar>& A, 
    const culite::dns::XxVector<T_Scalar>& x) 
{ 
    return culite::alias::VirtualProd_cscmv<T_Int,T_Scalar>(A.virtualize(), x.virtualize());
}

/*-------------------------------------------------*/

#endif // CULITE_OPERATORS_MULTMV_HPP_
