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

#ifndef CULITE_OPERATORS_SCALE_HPP_
#define CULITE_OPERATORS_SCALE_HPP_

/**
 * @file
 */

#include <culite/dense/dns_xxvector.hpp>
#include <culite/dense/dns_xxmatrix.hpp>
#include <culite/sparse/csr_xxmatrix.hpp>
#include <culite/sparse/csc_xxmatrix.hpp>

#include <culite/virtuals/virtual_object.hpp>
#include <culite/virtuals/virtual_scale.hpp>

/*-------------------------------------------------*/

/*
 * Generic scale operator
 */
template <typename T_Virtual> 
culite::VirtualScale<typename T_Virtual::result_type, T_Virtual>
operator*(
    typename T_Virtual::result_type::value_type val, 
    const culite::VirtualExpression<typename T_Virtual::result_type, T_Virtual>& src)
{ 
    return culite::VirtualScale<typename T_Virtual::result_type, T_Virtual>(src.self(), val); 
}

/**
 * @ingroup culite_module_index_math_operators_scal
 * @brief Scalar-vector multiplication.
 * @details Performs the operation @f$ \alpha \cdot x @f$.
 * @param[in] val The scalar coefficient.
 * @param[in] x The input vector.
 * @return The virtually scaled vector.
 */
template <typename T_Scalar>
culite::alias::VirtualScal_vec<T_Scalar>
operator*(T_Scalar val, const culite::dns::XxVector<T_Scalar>& x) 
{ 
    return culite::alias::VirtualScal_vec<T_Scalar>(x.virtualize(), val);
}

/**
 * @ingroup culite_module_index_math_operators_scal
 * @brief Scalar-matrix multiplication for dense matrices.
 * @details Performs the operation @f$ \alpha \cdot A @f$.
 * @param[in] val The scalar coefficient.
 * @param[in] A The input dense matrix.
 * @return The virtually scaled matrix.
 */
template <typename T_Scalar>
culite::alias::VirtualScal_dns<T_Scalar>
operator*(T_Scalar val, const culite::dns::XxMatrix<T_Scalar>& A) 
{ 
    return culite::alias::VirtualScal_dns<T_Scalar>(A.virtualize(), val);
}

/**
 * @ingroup culite_module_index_math_operators_scal
 * @brief Scalar-matrix multiplication for sparse matrices.
 * @details Performs the operation @f$ \alpha \cdot A @f$.
 * @param[in] val The scalar coefficient.
 * @param[in] A The input sparse matrix.
 * @return The virtually scaled matrix.
 */
template <typename T_Int, typename T_Scalar>
culite::alias::VirtualScal_csr<T_Int,T_Scalar>
operator*(T_Scalar val, const culite::csr::XxMatrix<T_Int,T_Scalar>& A) 
{ 
    return culite::alias::VirtualScal_csr<T_Int,T_Scalar>(A.virtualize(), val);
}

/**
 * @ingroup culite_module_index_math_operators_scal
 * @brief Scalar-matrix multiplication for sparse matrices.
 * @details Performs the operation @f$ \alpha \cdot A @f$.
 * @param[in] val The scalar coefficient.
 * @param[in] A The input sparse matrix.
 * @return The virtually scaled matrix.
 */
template <typename T_Int, typename T_Scalar>
culite::alias::VirtualScal_csc<T_Int,T_Scalar>
operator*(T_Scalar val, const culite::csc::XxMatrix<T_Int,T_Scalar>& A) 
{ 
    return culite::alias::VirtualScal_csc<T_Int,T_Scalar>(A.virtualize(), val);
}

/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/

/**
 * @ingroup culite_module_index_math_operators_scal
 * @brief In-place vector scaling.
 * @details Scales @p x by @p val, performing @f$ x = \alpha \cdot x @f$.
 * @param[in,out] x The dense vector to scale.
 * @param[in] val The scaling coefficient.
 */
template <typename T_Scalar>
void operator*=(
        culite::dns::XxVector<T_Scalar>& x, 
        T_Scalar val)
{
    x.iscale(val);
}

/**
 * @ingroup culite_module_index_math_operators_scal
 * @brief In-place dense matrix scaling.
 * @details Scales @p A by @p val, performing @f$ A = \alpha \cdot A @f$.
 * @param[in,out] A The dense matrix to scale.
 * @param[in] val The scaling coefficient.
 */
template <typename T_Scalar>
void operator*=(
        culite::dns::XxMatrix<T_Scalar>& A, 
        T_Scalar val)
{
    A.iscale(val);
}

/**
 * @ingroup culite_module_index_math_operators_scal
 * @brief In-place sparse matrix scaling.
 * @details Scales @p A by @p val, performing @f$ A = \alpha \cdot A @f$.
 * @param[in,out] A The sparse matrix to scale.
 * @param[in] val The scaling coefficient.
 */
template <typename T_Int, typename T_Scalar>
void operator*=(
        culite::csr::XxMatrix<T_Int,T_Scalar>& A, 
        T_Scalar val)
{
    A.iscale(val);
}

/**
 * @ingroup culite_module_index_math_operators_scal
 * @brief In-place sparse matrix scaling.
 * @details Scales @p A by @p val, performing @f$ A = \alpha \cdot A @f$.
 * @param[in,out] A The sparse matrix to scale.
 * @param[in] val The scaling coefficient.
 */
template <typename T_Int, typename T_Scalar>
void operator*=(
        culite::csc::XxMatrix<T_Int,T_Scalar>& A, 
        T_Scalar val)
{
    A.iscale(val);
}

/*-------------------------------------------------*/

#endif // CULITE_OPERATORS_SCALE_HPP_
