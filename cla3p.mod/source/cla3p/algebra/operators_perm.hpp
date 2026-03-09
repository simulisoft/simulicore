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

#ifndef CLA3P_OPERATORS_PERM_HPP_
#define CLA3P_OPERATORS_PERM_HPP_

/**
 * @file
 */

#include <cla3p/virtuals/virtual_expression.hpp>

#include <cla3p/dense/dns_xxvector.hpp>
#include <cla3p/dense/dns_xxmatrix.hpp>
#include <cla3p/sparse/csr_xxmatrix.hpp>
#include <cla3p/sparse/csc_xxmatrix.hpp>
#include <cla3p/perms/pxmatrix.hpp>

/*-------------------------------------------------*/
//namespace cla3p {
//namespace prm { template <typename T_Int> class PxMatrix; }
//namespace dns { template <typename T_Scalar> class XxVector; }
//namespace dns { template <typename T_Scalar> class XxMatrix; }
//namespace csr { template <typename T_Int, typename T_Scalar> class XxMatrix; }
//namespace csc { template <typename T_Int, typename T_Scalar> class XxMatrix; }
//} // namespace cla3p
/*-------------------------------------------------*/

/**
 * @ingroup cla3p_module_index_math_operators_perm
 * @brief Permutation-vector multiplication.
 * @details Performs the operation @f$ P \cdot x @f$, reordering vector entries.
 * @param[in] P The permutation matrix.
 * @param[in] x The input vector.
 * @return The permuted vector.
 */
template <typename T_Int, typename T_Scalar>
cla3p::dns::XxVector<T_Scalar> operator*(
    const cla3p::prm::PxMatrix<T_Int>& P, 
    const cla3p::dns::XxVector<T_Scalar>& x);

/**
 * @ingroup cla3p_module_index_math_operators_perm
 * @brief Left permutation-matrix multiplication for dense matrices.
 * @details Performs the operation @f$ P \cdot A @f$, reordering matrix rows.
 * @param[in] P The permutation matrix.
 * @param[in] A The input dense matrix.
 * @return The permuted matrix.
 */
template <typename T_Int, typename T_Scalar>
cla3p::dns::XxMatrix<T_Scalar> operator*(
    const cla3p::prm::PxMatrix<T_Int>& P, 
    const cla3p::dns::XxMatrix<T_Scalar>& A);

/**
 * @ingroup cla3p_module_index_math_operators_perm
 * @brief Right permutation-matrix multiplication for dense matrices.
 * @details Performs the operation @f$ A \cdot P @f$, reordering matrix columns.
 * @param[in] A The input dense matrix.
 * @param[in] P The permutation matrix.
 * @return The permuted matrix.
 */

template <typename T_Int, typename T_Scalar>
cla3p::dns::XxMatrix<T_Scalar> operator*(
    const cla3p::dns::XxMatrix<T_Scalar>& A, 
    const cla3p::prm::PxMatrix<T_Int>& P);

/**
 * @ingroup cla3p_module_index_math_operators_perm
 * @brief Left permutation-matrix multiplication for CSR sparse matrices.
 * @details Performs the operation @f$ P \cdot A @f$, reordering matrix rows.
 * @param[in] P The permutation matrix.
 * @param[in] A The input CSR sparse matrix.
 * @return The permuted matrix.
 */
template <typename T_Int, typename T_Scalar>
cla3p::csr::XxMatrix<T_Int,T_Scalar> operator*(
    const cla3p::prm::PxMatrix<T_Int>& P, 
    const cla3p::csr::XxMatrix<T_Int,T_Scalar>& A);

/**
 * @ingroup cla3p_module_index_math_operators_perm
 * @brief Right permutation-matrix multiplication for CSR sparse matrices.
 * @details Performs the operation @f$ A \cdot P @f$, reordering matrix columns.
 * @param[in] A The input CSR sparse matrix.
 * @param[in] P The permutation matrix.
 * @return The permuted matrix.
 */
template <typename T_Int, typename T_Scalar>
cla3p::csr::XxMatrix<T_Int,T_Scalar> operator*(
    const cla3p::csr::XxMatrix<T_Int,T_Scalar>& A, 
    const cla3p::prm::PxMatrix<T_Int>& P);

/**
 * @ingroup cla3p_module_index_math_operators_perm
 * @brief Left permutation-matrix multiplication for CSC sparse matrices.
 * @details Performs the operation @f$ P \cdot A @f$, reordering matrix rows.
 * @param[in] P The permutation matrix.
 * @param[in] A The input CSC sparse matrix.
 * @return The permuted matrix.
 */
template <typename T_Int, typename T_Scalar>
cla3p::csc::XxMatrix<T_Int,T_Scalar> operator*(
    const cla3p::prm::PxMatrix<T_Int>& P, 
    const cla3p::csc::XxMatrix<T_Int,T_Scalar>& A);

/**
 * @ingroup cla3p_module_index_math_operators_perm
 * @brief Right permutation-matrix multiplication for CSC sparse matrices.
 * @details Performs the operation @f$ A \cdot P @f$, reordering matrix columns.
 * @param[in] A The input CSC sparse matrix.
 * @param[in] P The permutation matrix.
 * @return The permuted matrix.
 */
template <typename T_Int, typename T_Scalar>
cla3p::csc::XxMatrix<T_Int,T_Scalar> operator*(
    const cla3p::csc::XxMatrix<T_Int,T_Scalar>& A, 
    const cla3p::prm::PxMatrix<T_Int>& P);

/**
 * @ingroup cla3p_module_index_math_operators_perm
 * @brief Permutation composition.
 * @details Performs the operation @f$ P \cdot Q @f$, composing two permutations.
 * @param[in] P The left permutation matrix.
 * @param[in] Q The right permutation matrix.
 * @return The composed permutation matrix.
 */
template <typename T_Int>
cla3p::prm::PxMatrix<T_Int> operator*(
    const cla3p::prm::PxMatrix<T_Int>& P,
    const cla3p::prm::PxMatrix<T_Int>& Q);

/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/

template <typename T_Int, typename T_Virtual>
typename T_Virtual::result_type operator*(
    const cla3p::prm::PxMatrix<T_Int>& P, 
    const cla3p::VirtualExpression<typename T_Virtual::result_type,T_Virtual>& v)
{
    return (P * v.evaluate());
}

template <typename T_Int, typename T_Virtual>
typename T_Virtual::result_type operator*(
    const cla3p::VirtualExpression<typename T_Virtual::result_type,T_Virtual>& v,
    const cla3p::prm::PxMatrix<T_Int>& P) 
{
    return (v.evaluate() * P);
}

/*-------------------------------------------------*/

#endif // CLA3P_OPERATORS_PERM_HPP_
