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
#include <culite/linsol/lapack_base.hpp>

// system

// 3rd

// culite
#include <culite/dense.hpp>
#include <culite/error/exceptions.hpp>

// forwards
#include <culite/checks/cla3p_forwards.hpp>

/*-------------------------------------------------*/
namespace culite {
/*-------------------------------------------------*/
template <typename T_Matrix>
LapackBase<T_Matrix>::LapackBase(CuSolverHandler& cusolver, decomp_t decompType)
    : m_cusolver(cusolver)
    , m_decompType(decompType)
{
}
/*-------------------------------------------------*/
template <typename T_Matrix>
LapackBase<T_Matrix>::~LapackBase()
{
    clear();
}
/*-------------------------------------------------*/
template <typename T_Matrix>
decomp_t LapackBase<T_Matrix>::decompType() const
{
    return m_decompType;
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void LapackBase<T_Matrix>::clear()
{
    m_factorMeta2D = Meta2D<int_t>();
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void LapackBase<T_Matrix>::reserve(const T_Matrix& mat)
{
    switch(decompType()) {

    case decomp_t::LU:
        m_cusolver.reserveLU(mat);
        break;

    default:
        throw err::CudaException("Unsupported decomposition type in Culite LapackBase reserve");
    }
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void LapackBase<T_Matrix>::prepareForDecomposition(const T_Matrix& mat)
{
    ::cla3p::decomp_generic_check(mat);

    if(decompType() == decomp_t::LU && !mat.prop().isGeneral()) {
        throw err::CudaException("Only general matrices are supported for LU decomposition");
    }

    reserve(mat);
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void LapackBase<T_Matrix>::decompose(const T_Matrix& mat)
{
    prepareForDecomposition(mat);

    switch(decompType()) {

        case decomp_t::LU:
            m_cusolver.decomposeLU(mat);
            m_factorMeta2D = static_cast<Meta2D<int_t>>(mat);
            break;

        default:
            throw err::CudaException("Only LU decomposition is supported in Culite LapackBase");
    }
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void LapackBase<T_Matrix>::prepareForSolution(T_Matrix& rhs) const
{
    ::cla3p::default_solve_input_check(m_factorMeta2D.ncols(), rhs);
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void LapackBase<T_Matrix>::solve(T_Matrix& rhs) const
{
    prepareForSolution(rhs);

    switch(decompType()) {

        case decomp_t::LU:
            m_cusolver.solveLU(rhs);
            break;

        default:
            throw err::CudaException("Only LU back substitution is supported in Culite LapackBase");
    }
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void LapackBase<T_Matrix>::solve(T_Vector& rhs) const
{
    T_Matrix tmp(rhs.size(), 1, rhs.values(), rhs.size(), false);
    solve(tmp);
}
/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/
template class LapackBase<dns::RdMatrix>;
template class LapackBase<dns::RfMatrix>;
template class LapackBase<dns::CdMatrix>;
template class LapackBase<dns::CfMatrix>;
/*-------------------------------------------------*/
} // namespace culite
/*-------------------------------------------------*/
