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
#include <culite/algebra/functional_multmm.hpp>

// system
#include <string>

// 3rd
#include <cla3p/error/exceptions.hpp>

// culite
#include <culite/bulk/csx.hpp>
#include <culite/error/exceptions.hpp>
#include <culite/dense/dns_xxmatrix.hpp>

// forwards
#include <culite/checks/cla3p_forwards.hpp>

/*-------------------------------------------------*/
namespace culite {
namespace ops {
/*-------------------------------------------------*/
template <typename T_Int>
static void throw_prop_compatibility_error(const MatrixMeta<T_Int>& A, 
                                           const MatrixMeta<T_Int>& B, 
                                           const MatrixMeta<T_Int>& C)
{
    std::string message = "Bad matrix property combo:\n";
    message = message + "A: " + A.prop().name() + "\n";
    message = message + "B: " + B.prop().name() + "\n";
    message = message + "C: " + C.prop().name();

    throw ::cla3p::err::InvalidOp(message);
}
/*-------------------------------------------------*/
template <typename T_Int>
static void throw_prop_compatibility_error(const MatrixMeta<T_Int>& A, 
                                           const MatrixMeta<T_Int>& B) 
{
    std::string message = "Bad matrix property combo:\n";
    message = message + "A: " + A.prop().name() + "\n";
    message = message + "B: " + B.prop().name();

    throw ::cla3p::err::InvalidOp(message);
}
/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/
template <typename T_Scalar>
void mult(T_Scalar alpha,
          op_t opA, const dns::XxMatrix<T_Scalar>& A,
          op_t opB, const dns::XxMatrix<T_Scalar>& B,
          T_Scalar beta, dns::XxMatrix<T_Scalar>& C,
          CuBlasHandler& cuBlasHandler)
{
    if(A.prop().isSymmetric() || A.prop().isHermitian()) opA = op_t::N;
    if(B.prop().isSymmetric() || B.prop().isHermitian()) opB = op_t::N;

    opA = (TypeTraits<T_Scalar>::is_real() && opA == op_t::C ? op_t::T : opA);
    opB = (TypeTraits<T_Scalar>::is_real() && opB == op_t::C ? op_t::T : opB);

    ::cla3p::mult_dim_check(opA, A, opB, B, C);
    hermitian_coeff_check2<T_Scalar>(C.prop(), alpha);
    hermitian_coeff_check2<T_Scalar>(C.prop(), beta);

    if(A.prop().isGeneral() && B.prop().isGeneral() && C.prop().isGeneral()) {

        int_t k = (opA == op_t::N ? A.ncols() : A.nrows());

        cuBlasHandler.gemm(opA,
                           opB,
                           C.nrows(), C.ncols(), k,
                           &alpha,
                           A.values(), A.ld(),
                           B.values(), B.ld(),
                           &beta,
                           C.values(), C.ld());

    } else if(A.prop().isSymmetric() && B.prop().isGeneral() && C.prop().isGeneral()) {

        cuBlasHandler.symm(side_t::Left,
                           A.prop().uplo(),
                           C.nrows(), C.ncols(),
                           &alpha,
                           A.values(), A.ld(),
                           B.values(), B.ld(),
                           &beta,
                           C.values(), C.ld());

    } else if(A.prop().isHermitian() && B.prop().isGeneral() && C.prop().isGeneral()) {

        cuBlasHandler.hemm(side_t::Left,
                           A.prop().uplo(),
                           C.nrows(), C.ncols(),
                           &alpha,
                           A.values(), A.ld(),
                           B.values(), B.ld(),
                           &beta,
                           C.values(), C.ld());

    } else if(A.prop().isGeneral() && B.prop().isGeneral() && C.prop().isSymmetric()) {

        op_t op = op_t::N;

        if(opA == op_t::N && opB == op_t::T) {
            op = op_t::N;
        } else if(opA == op_t::T && opB == op_t::N) {
            op = op_t::T;
        } else {
            std::stringstream ss;
            ss << "Invalid opA/opB combo for symmetric rank-k update: opA=" << opA << ", opB=" << opB;
            throw err::CudaException(ss.str());
        }

        int_t k = (opA == op_t::N ? A.ncols() : A.nrows());

        cuBlasHandler.syrkx(C.prop().uplo(),
                            op,
                            C.ncols(), k,
                            &alpha,
                            A.values(), A.ld(),
                            B.values(), B.ld(),
                            &beta,
                            C.values(), C.ld());

    } else if(A.prop().isGeneral() && B.prop().isGeneral() && C.prop().isHermitian()) {

        op_t op = op_t::N;

        if(opA == op_t::N && opB == op_t::C) {
            op = op_t::N;
        } else if(opA == op_t::C && opB == op_t::N) {
            op = op_t::C;
        } else {
            std::stringstream ss;
            ss << "Invalid opA/opB combo for Hermitian rank-k update: opA=" << opA << ", opB=" << opB;
            throw err::CudaException(ss.str());
        }

        int_t k = (opA == op_t::N ? A.ncols() : A.nrows());
        typename TypeTraits<T_Scalar>::real_type realBeta = arith::getRe(beta);

        cuBlasHandler.herkx(C.prop().uplo(),
                            op,
                            C.ncols(), k,
                            &alpha,
                            A.values(), A.ld(),
                            B.values(), B.ld(),
                            &realBeta,
                            C.values(), C.ld());

    } else {

        throw_prop_compatibility_error(A, B, C);

    } // property combos
}
/*-------------------------------------------------*/
#define instantiate_mult(T_Scl) \
template void mult(T_Scl, \
                   op_t, const dns::XxMatrix<T_Scl>&, \
                   op_t, const dns::XxMatrix<T_Scl>&, \
                   T_Scl, dns::XxMatrix<T_Scl>&, \
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
void mult(T_Scalar alpha, op_t opA,
          const csr::XxMatrix<T_Int,T_Scalar>& A,
          const dns::XxMatrix<T_Scalar>& B,
          T_Scalar beta, dns::XxMatrix<T_Scalar>& C,
          CuSparseHandler& cuSparseHandler)
{
    if(A.prop().isSymmetric() || A.prop().isHermitian()) opA = op_t::N;

    opA = (TypeTraits<T_Scalar>::is_real() && opA == op_t::C ? op_t::T : opA);
    op_t opB = op_t::N;

    ::cla3p::mult_dim_check(opA, A, opB, B, C);

    cusparse::SpMatCsr<T_Scalar> csrA(A.nrows(), A.ncols(), A.nnz(), A.rowptr(), A.colidx(), A.values());
    cusparse::DnMat<T_Scalar> dnsB(B.nrows(), B.ncols(), B.values(), B.ld());
    cusparse::DnMat<T_Scalar> dnsC(C.nrows(), C.ncols(), C.values(), C.ld());

    if(A.prop().isGeneral() && B.prop().isGeneral() && C.prop().isGeneral()) {

        cuSparseHandler.reserveSpmm(opA, opB, &alpha, csrA, dnsB, &beta, dnsC);
        cuSparseHandler.preprocessSpmm(opA, opB, &alpha, csrA, dnsB, &beta, dnsC);
        cuSparseHandler.performSpmm(opA, opB, &alpha, csrA, dnsB, &beta, dnsC);

    } else if((A.prop().isSymmetric() || A.prop().isHermitian()) && B.prop().isGeneral() && C.prop().isGeneral()) {

        {
            // C = beta * C + alpha * A{uplo} * B
            op_t op = op_t::N;
            cuSparseHandler.reserveSpmm(op, opB, &alpha, csrA, dnsB, &beta, dnsC);
            cuSparseHandler.preprocessSpmm(op, opB, &alpha, csrA, dnsB, &beta, dnsC);
            cuSparseHandler.performSpmm(op, opB, &alpha, csrA, dnsB, &beta, dnsC);
        }
        {
            // C = C + alpha * A{uplo}.transpose() * B
            T_Scalar betaOne = makeScalar<T_Scalar>(1);
            op_t op = (A.prop().isSymmetric() ? op_t::T : op_t::C);
            cuSparseHandler.reserveSpmm(op, opB, &alpha, csrA, dnsB, &betaOne, dnsC);
            cuSparseHandler.preprocessSpmm(op, opB, &alpha, csrA, dnsB, &betaOne, dnsC);
            cuSparseHandler.performSpmm(op, opB, &alpha, csrA, dnsB, &betaOne, dnsC);
        }
        {
            // y = y - alpha * A{diag} * x
            T_Scalar alphaMinus = -alpha;
            blk::csx::diag_times_mat<T_Int, T_Scalar>(&alphaMinus, 
                                                      A.nrows(), 
                                                      A.rowptr(), 
                                                      A.colidx(), 
                                                      A.values(), 
                                                      C.ncols(),
                                                      B.values(), B.ld(), 
                                                      C.values(), C.ld());
        }

    } else {

        throw_prop_compatibility_error(A, B, C);

    } // property combos
}
/*-------------------------------------------------*/
#define instantiate_mult(T_Int, T_Scl) \
template void mult(T_Scl, op_t, \
              const csr::XxMatrix<T_Int,T_Scl>&, \
              const dns::XxMatrix<T_Scl>&, \
              T_Scl, dns::XxMatrix<T_Scl>&, \
              CuSparseHandler&)
instantiate_mult(int_t, real_t);
instantiate_mult(int_t, real4_t);
instantiate_mult(int_t, complex_t);
instantiate_mult(int_t, complex8_t);
#undef instantiate_mult
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
void mult(T_Scalar alpha, op_t opA,
          const csc::XxMatrix<T_Int,T_Scalar>& A,
          const dns::XxMatrix<T_Scalar>& B,
          T_Scalar beta, dns::XxMatrix<T_Scalar>& C,
          CuSparseHandler& cuSparseHandler)
{
    if(A.prop().isSymmetric() || A.prop().isHermitian()) opA = op_t::N;

    opA = (TypeTraits<T_Scalar>::is_real() && opA == op_t::C ? op_t::T : opA);
    op_t opB = op_t::N;

    ::cla3p::mult_dim_check(opA, A, opB, B, C);

    cusparse::SpMatCsc<T_Scalar> cscA(A.nrows(), A.ncols(), A.nnz(), A.colptr(), A.rowidx(), A.values());
    cusparse::DnMat<T_Scalar> dnsB(B.nrows(), B.ncols(), B.values(), B.ld());
    cusparse::DnMat<T_Scalar> dnsC(C.nrows(), C.ncols(), C.values(), C.ld());

    if(A.prop().isGeneral() && B.prop().isGeneral() && C.prop().isGeneral()) {

        cuSparseHandler.reserveSpmm(opA, opB, &alpha, cscA, dnsB, &beta, dnsC);
        cuSparseHandler.preprocessSpmm(opA, opB, &alpha, cscA, dnsB, &beta, dnsC);
        cuSparseHandler.performSpmm(opA, opB, &alpha, cscA, dnsB, &beta, dnsC);

    } else if((A.prop().isSymmetric() || A.prop().isHermitian()) && B.prop().isGeneral() && C.prop().isGeneral()) {

        {
            // C = beta * C + alpha * A{uplo} * B
            op_t op = op_t::N;
            cuSparseHandler.reserveSpmm(op, opB, &alpha, cscA, dnsB, &beta, dnsC);
            cuSparseHandler.preprocessSpmm(op, opB, &alpha, cscA, dnsB, &beta, dnsC);
            cuSparseHandler.performSpmm(op, opB, &alpha, cscA, dnsB, &beta, dnsC);
        }
        {
            // C = C + alpha * A{uplo}.transpose() * B
            T_Scalar betaOne = makeScalar<T_Scalar>(1);
            op_t op = (A.prop().isSymmetric() ? op_t::T : op_t::C);
            cuSparseHandler.reserveSpmm(op, opB, &alpha, cscA, dnsB, &betaOne, dnsC);
            cuSparseHandler.preprocessSpmm(op, opB, &alpha, cscA, dnsB, &betaOne, dnsC);
            cuSparseHandler.performSpmm(op, opB, &alpha, cscA, dnsB, &betaOne, dnsC);
        }
        {
            // y = y - alpha * A{diag} * x
            T_Scalar alphaMinus = -alpha;
            blk::csx::diag_times_mat<T_Int, T_Scalar>(&alphaMinus, 
                                                      A.ncols(), 
                                                      A.colptr(), 
                                                      A.rowidx(), 
                                                      A.values(), 
                                                      C.ncols(),
                                                      B.values(), B.ld(), 
                                                      C.values(), C.ld());
        }

    } else {

        throw_prop_compatibility_error(A, B, C);

    } // property combos
}
/*-------------------------------------------------*/
#define instantiate_mult(T_Int, T_Scl) \
template void mult(T_Scl, op_t, \
                   const csc::XxMatrix<T_Int,T_Scl>&, \
                   const dns::XxMatrix<T_Scl>&, \
                   T_Scl, dns::XxMatrix<T_Scl>&, \
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
