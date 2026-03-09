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
#include <cla3p/algebra/functional_multmv.hpp>

// system

// 3rd

// cla3p
#include <cla3p/error.hpp>
#include <cla3p/checks/matrix_math_checks.hpp>
#include <cla3p/proxies/blas_proxy.hpp>
#include <cla3p/bulk/dns_math.hpp>
#include <cla3p/bulk/csr_math.hpp>
#include <cla3p/bulk/csc_math.hpp>
#include <cla3p/dense/dns_xxvector.hpp>
#include <cla3p/dense/dns_xxmatrix.hpp>
#include <cla3p/sparse/csr_xxmatrix.hpp>
#include <cla3p/sparse/csc_xxmatrix.hpp>
#include <cla3p/algebra/functional_update.hpp>

/*-------------------------------------------------*/
namespace cla3p {
namespace ops {
/*-------------------------------------------------*/
template <typename T_Scalar>
void mult(T_Scalar alpha, op_t opA,
    const dns::XxMatrix<T_Scalar>& A,
    const dns::XxVector<T_Scalar>& x,
        T_Scalar beta,
    dns::XxVector<T_Scalar>& y)
{
    Operation _opA(opA);
    mat_x_vec_mult_check(_opA, A.prop(), A.nrows(), A.ncols(), x.size(), y.size());

    if(A.prop().isGeneral()) {

        blk::dns::gem_x_vec(opA, A.nrows(), A.ncols(), alpha, A.values(), A.ld(), x.values(), beta, y.values());

    } else if(A.prop().isSymmetric()) {

        blk::dns::sym_x_vec(A.prop().uplo(), A.ncols(), alpha, A.values(), A.ld(), x.values(), beta, y.values());

    } else if(A.prop().isHermitian()) {

        blk::dns::hem_x_vec(A.prop().uplo(), A.ncols(), alpha, A.values(), A.ld(), x.values(), beta, y.values());

    } else if(A.prop().isTriangular()) {

        if(beta == T_Scalar(0)) {

            blk::dns::trm_x_vec(A.prop().uplo(), opA, A.nrows(), A.ncols(), alpha, A.values(), A.ld(), x.values(), y.values());

        } else {

            dns::XxVector<T_Scalar> tmp(y.size());
            blk::dns::trm_x_vec(A.prop().uplo(), opA, A.nrows(), A.ncols(), alpha, A.values(), A.ld(), x.values(), tmp.values());
            y.iscale(beta);
            ops::update(T_Scalar(1), tmp, y);

        } // beta

    } else {

        throw err::Exception();

    } // property 
}
/*-------------------------------------------------*/
#define instantiate_mult(T_Scl) \
template void mult(T_Scl, op_t, \
    const dns::XxMatrix<T_Scl>&, \
    const dns::XxVector<T_Scl>&, \
        T_Scl, \
    dns::XxVector<T_Scl>&)
instantiate_mult(real_t);
instantiate_mult(real4_t);
instantiate_mult(complex_t);
instantiate_mult(complex8_t);
#undef instantiate_mult
/*-------------------------------------------------*/
template <typename T_Scalar>
void trimult(op_t opA,
        const dns::XxMatrix<T_Scalar>& A,
        dns::XxVector<T_Scalar>& x)
{
    Operation _opA(opA);
    trivec_mult_replace_check(A.prop(), A.nrows(), A.ncols(), _opA, x.size());

    blas::trmv(A.prop().cuplo(), _opA.ctype(), 'N', A.ncols(), A.values(), A.ld(), x.values(), 1);
}
/*-------------------------------------------------*/
#define instantiate_trimult(T_Scl) \
template void trimult(op_t, \
    const dns::XxMatrix<T_Scl>&, \
    dns::XxVector<T_Scl>&)
instantiate_trimult(real_t);
instantiate_trimult(real4_t);
instantiate_trimult(complex_t);
instantiate_trimult(complex8_t);
#undef instantiate_trimult
/*-------------------------------------------------*/
template <typename T_Scalar>
void trisol(op_t opA,
    const dns::XxMatrix<T_Scalar>& A,
    dns::XxVector<T_Scalar>& b)
{
    Operation _opA(opA);
    trivec_mult_replace_check(A.prop(), A.nrows(), A.ncols(), _opA, b.size());

    blas::trsv(A.prop().cuplo(), _opA.ctype(), 'N', A.ncols(), A.values(), A.ld(), b.values(), 1);
}
/*-------------------------------------------------*/
#define instantiate_trisol(T_Scl) \
template void trisol(op_t, \
    const dns::XxMatrix<T_Scl>&, \
    dns::XxVector<T_Scl>&)
instantiate_trisol(real_t);
instantiate_trisol(real4_t);
instantiate_trisol(complex_t);
instantiate_trisol(complex8_t);
#undef instantiate_trisol
/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
void mult(T_Scalar alpha, op_t opA,
          const csr::XxMatrix<T_Int,T_Scalar>& A,
          const dns::XxVector<T_Scalar>& x,
          T_Scalar beta,
          dns::XxVector<T_Scalar>& y)
{
    Operation _opA(opA);
    mat_x_vec_mult_check(_opA, A.prop(), A.nrows(), A.ncols(), x.size(), y.size());

    if(A.prop().isGeneral() || A.prop().isTriangular()) {

        blk::csr::gem_x_vec(opA, A.nrows(), A.ncols(), alpha, 
                            A.rowptr(), A.colidx(), A.values(), 
                            x.values(), beta, y.values());

    } else if(A.prop().isSymmetric()) {

        blk::csr::sym_x_vec(A.prop().uplo(), A.ncols(), alpha, 
                            A.rowptr(), A.colidx(), A.values(), 
                            x.values(), beta, y.values());

    } else if(A.prop().isHermitian()) {

        blk::csr::hem_x_vec(A.prop().uplo(), A.ncols(), alpha, 
                            A.rowptr(), A.colidx(), A.values(), 
                            x.values(), beta, y.values());

    } else {

        throw err::Exception();

    } // property 
}
/*-------------------------------------------------*/
#define instantiate_mult(T_Int, T_Scl) \
template void mult(T_Scl, op_t, \
                   const csr::XxMatrix<T_Int, T_Scl>&, \
                   const dns::XxVector<T_Scl>&, \
                   T_Scl, \
                   dns::XxVector<T_Scl>&)
instantiate_mult(int_t, real_t);
instantiate_mult(int_t, real4_t);
instantiate_mult(int_t, complex_t);
instantiate_mult(int_t, complex8_t);
#undef instantiate_mult
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
void mult(T_Scalar alpha, op_t opA,
          const csc::XxMatrix<T_Int,T_Scalar>& A,
          const dns::XxVector<T_Scalar>& x,
          T_Scalar beta, 
          dns::XxVector<T_Scalar>& y)
{
    Operation _opA(opA);
    mat_x_vec_mult_check(_opA, A.prop(), A.nrows(), A.ncols(), x.size(), y.size());

    if(A.prop().isGeneral() || A.prop().isTriangular()) {

        blk::csc::gem_x_vec(opA, A.nrows(), A.ncols(), alpha, 
                            A.colptr(), A.rowidx(), A.values(), 
                            x.values(), beta, y.values());

    } else if(A.prop().isSymmetric()) {

        blk::csc::sym_x_vec(A.prop().uplo(), A.ncols(), alpha, 
                            A.colptr(), A.rowidx(), A.values(), 
                            x.values(), beta, y.values());

    } else if(A.prop().isHermitian()) {

        blk::csc::hem_x_vec(A.prop().uplo(), A.ncols(), alpha, 
                            A.colptr(), A.rowidx(), A.values(), 
                            x.values(), beta, y.values());

    } else {

        throw err::Exception();

    } // property 
}
/*-------------------------------------------------*/
#define instantiate_mult(T_Int, T_Scl) \
template void mult(T_Scl, op_t, \
                   const csc::XxMatrix<T_Int, T_Scl>&, \
                   const dns::XxVector<T_Scl>&, \
                   T_Scl, \
                   dns::XxVector<T_Scl>&)
instantiate_mult(int_t, real_t);
instantiate_mult(int_t, real4_t);
instantiate_mult(int_t, complex_t);
instantiate_mult(int_t, complex8_t);
#undef instantiate_mult
/*-------------------------------------------------*/
} // namespace ops
} // namespace cla3p
/*-------------------------------------------------*/
