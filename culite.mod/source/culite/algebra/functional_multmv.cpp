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
#include "culite/algebra/functional_multmv.hpp"

// system
#include <string>
#include <sstream>

// 3rd
#include <cla3p/checks/matrix_math_checks.hpp>

// culite
#include "culite/bulk/csx.hpp"
#include "culite/error/exceptions.hpp"
#include "culite/dense/dns_xxvector.hpp"
#include "culite/dense/dns_xxmatrix.hpp"
#include "culite/sparse/csr_xxmatrix.hpp"
#include "culite/sparse/csc_xxmatrix.hpp"

/*-------------------------------------------------*/
namespace culite {
namespace ops {
/*-------------------------------------------------*/
template <typename T_Scalar>
void mult(T_Scalar alpha, ::cla3p::op_t opA,
          const dns::XxMatrix<T_Scalar>& A,
          const dns::XxVector<T_Scalar>& x,
          T_Scalar beta,
          dns::XxVector<T_Scalar>& y,
          CuBlasHandler& cuBlasHandler)
{
    ::cla3p::mult_dim_check(opA, A, x, y);

    if(A.prop().isGeneral()) {

        cuBlasHandler.gemv(opA, 
                           A.nrows(), 
                           A.ncols(), 
                           &alpha, 
                           A.values(), A.ld(), 
                           x.values(), 1, 
                           &beta, 
                           y.values(), 1);

    } else if(A.prop().isSymmetric()) {

        cuBlasHandler.symv(A.prop().uplo(), 
                           A.ncols(), 
                           &alpha, 
                           A.values(), A.ld(), 
                           x.values(), 1, 
                           &beta, 
                           y.values(), 1);

    } else if(A.prop().isHermitian()) {

        cuBlasHandler.hemv(A.prop().uplo(), 
                           A.ncols(), 
                           &alpha, 
                           A.values(), A.ld(), 
                           x.values(), 1, 
                           &beta, 
                           y.values(), 1);
    } else {

        std::stringstream ss;
        ss << "Invalid matrix property for matrix-vector multiplication: " << A.prop();
        throw err::CudaException(ss.str());

    } // property 
}
/*-------------------------------------------------*/
#define instantiate_mult(T_Scl) \
template void mult(T_Scl, ::cla3p::op_t, \
                   const dns::XxMatrix<T_Scl>&, \
                   const dns::XxVector<T_Scl>&, \
                   T_Scl, \
                   dns::XxVector<T_Scl>&, \
                   CuBlasHandler&)
instantiate_mult(real_t);
instantiate_mult(real4_t);
instantiate_mult(complex_t);
instantiate_mult(complex8_t);
#undef instantiate_mult
/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
void mult(T_Scalar alpha, ::cla3p::op_t opA,
          const csr::XxMatrix<T_Int,T_Scalar>& A,
          const dns::XxVector<T_Scalar>& x,
          T_Scalar beta,
          dns::XxVector<T_Scalar>& y,
          CuSparseHandler& cuSparseHandler)
{
     // ignore opA for symmetric/hermitian matrices
    if(A.prop().isSymmetric() || A.prop().isHermitian()) opA = ::cla3p::op_t::N;

    ::cla3p::mult_dim_check(opA, A, x, y);

    cusparse::SpMatCsr<T_Scalar> csrA(A.nrows(), A.ncols(), A.nnz(), A.rowptr(), A.colidx(), A.values());
    cusparse::DnVec<T_Scalar> vecX(x.size(), x.values());
    cusparse::DnVec<T_Scalar> vecY(y.size(), y.values());

    if(A.prop().isGeneral()) {

        cuSparseHandler.reserveSpmv(opA, &alpha, csrA, vecX, &beta, vecY);
        cuSparseHandler.preprocessSpmv(opA, &alpha, csrA, vecX, &beta, vecY);
        cuSparseHandler.performSpmv(opA, &alpha, csrA, vecX, &beta, vecY);

    } else if(A.prop().isSymmetric() || A.prop().isHermitian()) {

        {
            // y = beta * y + alpha * A{uplo} * x
            ::cla3p::op_t op = ::cla3p::op_t::N;
            cuSparseHandler.reserveSpmv(op, &alpha, csrA, vecX, &beta, vecY);
            cuSparseHandler.preprocessSpmv(op, &alpha, csrA, vecX, &beta, vecY);
            cuSparseHandler.performSpmv(op, &alpha, csrA, vecX, &beta, vecY);
        }
        {
            // y = y + alpha * A{uplo}.transpose() * x
            T_Scalar betaOne = makeScalar<T_Scalar>(1);
            ::cla3p::op_t op = (A.prop().isSymmetric() ? ::cla3p::op_t::T : ::cla3p::op_t::C);
            cuSparseHandler.reserveSpmv(op, &alpha, csrA, vecX, &betaOne, vecY);
            cuSparseHandler.preprocessSpmv(op, &alpha, csrA, vecX, &betaOne, vecY);
            cuSparseHandler.performSpmv(op, &alpha, csrA, vecX, &betaOne, vecY);
        }
        {
            // y = y - alpha * A{diag} * x
            T_Scalar alphaMinus = -alpha;
            blk::csx::diag_times_vec<T_Int, T_Scalar>(&alphaMinus, 
                                                      A.nrows(), 
                                                      A.rowptr(), 
                                                      A.colidx(), 
                                                      A.values(), 
                                                      x.values(), 1, 
                                                      y.values(), 1);
        }

    } else {

        std::stringstream ss;
        ss << "Invalid matrix property for matrix-vector multiplication: " << A.prop();
        throw err::CudaException(ss.str());

    } // property 
}
/*-------------------------------------------------*/
#define instantiate_mult(T_Int, T_Scl) \
template void mult(T_Scl, ::cla3p::op_t, \
                   const csr::XxMatrix<T_Int, T_Scl>&, \
                   const dns::XxVector<T_Scl>&, \
                   T_Scl, \
                   dns::XxVector<T_Scl>&, \
                   CuSparseHandler&)
instantiate_mult(int_t, real_t);
instantiate_mult(int_t, real4_t);
instantiate_mult(int_t, complex_t);
instantiate_mult(int_t, complex8_t);
#undef instantiate_mult
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
void mult(T_Scalar alpha, ::cla3p::op_t opA,
          const csc::XxMatrix<T_Int,T_Scalar>& A,
          const dns::XxVector<T_Scalar>& x,
          T_Scalar beta, 
          dns::XxVector<T_Scalar>& y,
          CuSparseHandler& cuSparseHandler)
{
    // ignore opA for symmetric/hermitian matrices
    if(A.prop().isSymmetric() || A.prop().isHermitian()) opA = ::cla3p::op_t::N;

    ::cla3p::mult_dim_check(opA, A, x, y);

    cusparse::SpMatCsc<T_Scalar> cscA(A.nrows(), A.ncols(), A.nnz(), A.colptr(), A.rowidx(), A.values());
    cusparse::DnVec<T_Scalar> vecX(x.size(), x.values());
    cusparse::DnVec<T_Scalar> vecY(y.size(), y.values());


    if(A.prop().isGeneral()) {

        cuSparseHandler.reserveSpmv(opA, &alpha, cscA, vecX, &beta, vecY);
        cuSparseHandler.preprocessSpmv(opA, &alpha, cscA, vecX, &beta, vecY);
        cuSparseHandler.performSpmv(opA, &alpha, cscA, vecX, &beta, vecY);

    } else if(A.prop().isSymmetric() || A.prop().isHermitian()) {

        {
            // y = beta * y + alpha * A{uplo} * x
            ::cla3p::op_t op = ::cla3p::op_t::N;
            cuSparseHandler.reserveSpmv(op, &alpha, cscA, vecX, &beta, vecY);
            cuSparseHandler.preprocessSpmv(op, &alpha, cscA, vecX, &beta, vecY);
            cuSparseHandler.performSpmv(op, &alpha, cscA, vecX, &beta, vecY);
        }
        {
            // y = y + alpha * A{uplo}.transpose() * x
            T_Scalar betaOne = makeScalar<T_Scalar>(1);
            ::cla3p::op_t op = (A.prop().isSymmetric() ? ::cla3p::op_t::T : ::cla3p::op_t::C);
            cuSparseHandler.reserveSpmv(op, &alpha, cscA, vecX, &betaOne, vecY);
            cuSparseHandler.preprocessSpmv(op, &alpha, cscA, vecX, &betaOne, vecY);
            cuSparseHandler.performSpmv(op, &alpha, cscA, vecX, &betaOne, vecY);
        }
        {
            // y = y - alpha * A{diag} * x
            T_Scalar alphaMinus = -alpha;
            blk::csx::diag_times_vec<T_Int, T_Scalar>(&alphaMinus, 
                                                      A.ncols(), 
                                                      A.colptr(), 
                                                      A.rowidx(), 
                                                      A.values(), 
                                                      x.values(), 1, 
                                                      y.values(), 1);
        }

    } else {

        std::stringstream ss;
        ss << "Invalid matrix property for matrix-vector multiplication: " << A.prop();
        throw err::CudaException(ss.str());

    } // property 
}
/*-------------------------------------------------*/
#define instantiate_mult(T_Int, T_Scl) \
template void mult(T_Scl, ::cla3p::op_t, \
                   const csc::XxMatrix<T_Int, T_Scl>&, \
                   const dns::XxVector<T_Scl>&, \
                   T_Scl, \
                   dns::XxVector<T_Scl>&, \
                   CuSparseHandler&)
instantiate_mult(int_t, real_t);
instantiate_mult(int_t, real4_t);
instantiate_mult(int_t, complex_t);
instantiate_mult(int_t, complex8_t);
#undef instantiate_mult
/*-------------------------------------------------*/
} // namespace ops
} // namespace culite
/*-------------------------------------------------*/
