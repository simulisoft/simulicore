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
#include <culite/generic/cusparse_handler.hpp>

// system

// 3rd

// culite
#include <culite/error/cuda.hpp>
#include <culite/types/scalar.hpp>

/*-------------------------------------------------*/
namespace culite {
/*-------------------------------------------------*/
CuSparseSpm::CuSparseSpm(CuSparseHandler& cuSparseHandler, cusparsePointerMode_t mode)
    : m_handler(cuSparseHandler)
{
    m_oldMode = m_handler.setPointerMode(mode);
}
/*-------------------------------------------------*/
CuSparseSpm::~CuSparseSpm()
{
    m_handler.setPointerMode(m_oldMode);
}
/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/
CuSparseHandler::CuSparseHandler()
{
    defaults();
    cusparseStatus_t cusparseStatus = cusparseCreate(&m_handle);
    err::check_cusparse(cusparseStatus);
}
/*-------------------------------------------------*/
CuSparseHandler::~CuSparseHandler()
{
    clear();
    cusparseStatus_t cusparseStatus = cusparseDestroy(m_handle);
    err::check_cusparse(cusparseStatus);
}
/*-------------------------------------------------*/
void CuSparseHandler::defaults()
{
    m_workspaceInBytes = 0;
}
/*-------------------------------------------------*/
void CuSparseHandler::clear()
{
    deviceWork().clear();

    defaults();
}
/*-------------------------------------------------*/
cusparsePointerMode_t CuSparseHandler::setPointerMode(cusparsePointerMode_t mode)
{
    cusparsePointerMode_t ret = pointerMode();
    if (ret != mode) {
        cusparseStatus_t status = cusparseSetPointerMode(m_handle, mode);
        err::check_cusparse(status);
    }
    return ret;
}
/*-------------------------------------------------*/
cusparsePointerMode_t CuSparseHandler::pointerMode()
{
    cusparsePointerMode_t ret;
    cusparseStatus_t status = cusparseGetPointerMode(m_handle, &ret);
    err::check_cusparse(status);
    return ret;
}
/*-------------------------------------------------*/
cusparseHandle_t CuSparseHandler::handle()
{
    return m_handle;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void CuSparseHandler::reserveSpmv(op_t                             opA,
                                  const T_Scalar*                  alpha,
                                  const cusparse::SpMatBase&       matA,
                                  const cusparse::DnVec<T_Scalar>& vecX,
                                  const T_Scalar*                  beta,
                                  cusparse::DnVec<T_Scalar>&       vecY,
                                  cusparseSpMVAlg_t                alg)
{
    cusparseStatus_t status = cusparseSpMV_bufferSize(handle(),
                                                      cusparse::cla3pOp2cusparseOp(opA),
                                                      alpha,
                                                      matA.descr(),
                                                      vecX.descr(),
                                                      beta,
                                                      vecY.descr(),
                                                      TypeTraits<T_Scalar>::cuda_type(),
                                                      alg,
                                                      &m_workspaceInBytes);
    err::check_cusparse(status);

    deviceWork().reserve(m_workspaceInBytes);
}
//-------------------------------------------------
template void CuSparseHandler::reserveSpmv<real_t>(op_t, const real_t*, const cusparse::SpMatBase&, const cusparse::DnVec<real_t>&, const real_t*, cusparse::DnVec<real_t>&, cusparseSpMVAlg_t);
template void CuSparseHandler::reserveSpmv<real4_t>(op_t, const real4_t*, const cusparse::SpMatBase&, const cusparse::DnVec<real4_t>&, const real4_t*, cusparse::DnVec<real4_t>&, cusparseSpMVAlg_t);
template void CuSparseHandler::reserveSpmv<complex_t>(op_t, const complex_t*, const cusparse::SpMatBase&, const cusparse::DnVec<complex_t>&, const complex_t*, cusparse::DnVec<complex_t>&, cusparseSpMVAlg_t);
template void CuSparseHandler::reserveSpmv<complex8_t>(op_t, const complex8_t*, const cusparse::SpMatBase&, const cusparse::DnVec<complex8_t>&, const complex8_t*, cusparse::DnVec<complex8_t>&, cusparseSpMVAlg_t);
/*-------------------------------------------------*/
template <typename T_Scalar>
void CuSparseHandler::preprocessSpmv(op_t                             opA,
                                     const T_Scalar*                  alpha,
                                     const cusparse::SpMatBase&       matA,
                                     const cusparse::DnVec<T_Scalar>& vecX,
                                     const T_Scalar*                  beta,
                                     cusparse::DnVec<T_Scalar>&       vecY,
                                     cusparseSpMVAlg_t                alg)
{
    cusparseStatus_t status = cusparseSpMV_preprocess(handle(),
                                                      cusparse::cla3pOp2cusparseOp(opA),
                                                      alpha,
                                                      matA.descr(),
                                                      vecX.descr(),
                                                      beta,
                                                      vecY.descr(),
                                                      TypeTraits<T_Scalar>::cuda_type(),
                                                      alg,
                                                      deviceWork().data());
    err::check_cusparse(status);
}
//-------------------------------------------------
template void CuSparseHandler::preprocessSpmv<real_t>(op_t, const real_t*, const cusparse::SpMatBase&, const cusparse::DnVec<real_t>&, const real_t*, cusparse::DnVec<real_t>&, cusparseSpMVAlg_t);
template void CuSparseHandler::preprocessSpmv<real4_t>(op_t, const real4_t*, const cusparse::SpMatBase&, const cusparse::DnVec<real4_t>&, const real4_t*, cusparse::DnVec<real4_t>&, cusparseSpMVAlg_t);
template void CuSparseHandler::preprocessSpmv<complex_t>(op_t, const complex_t*, const cusparse::SpMatBase&, const cusparse::DnVec<complex_t>&, const complex_t*, cusparse::DnVec<complex_t>&, cusparseSpMVAlg_t);
template void CuSparseHandler::preprocessSpmv<complex8_t>(op_t, const complex8_t*, const cusparse::SpMatBase&, const cusparse::DnVec<complex8_t>&, const complex8_t*, cusparse::DnVec<complex8_t>&, cusparseSpMVAlg_t);
/*-------------------------------------------------*/
template <typename T_Scalar>
void CuSparseHandler::performSpmv(op_t                             opA,
                                  const T_Scalar*                  alpha,
                                  const cusparse::SpMatBase&       matA,
                                  const cusparse::DnVec<T_Scalar>& vecX,
                                  const T_Scalar*                  beta,
                                  cusparse::DnVec<T_Scalar>&       vecY,
                                  cusparseSpMVAlg_t                alg)
{
    cusparseStatus_t status = cusparseSpMV(handle(),
                                           cusparse::cla3pOp2cusparseOp(opA),
                                           alpha,
                                           matA.descr(),
                                           vecX.descr(),
                                           beta,
                                           vecY.descr(),
                                           TypeTraits<T_Scalar>::cuda_type(),
                                           alg,
                                           deviceWork().data());
    err::check_cusparse(status);
}
//-------------------------------------------------
template void CuSparseHandler::performSpmv<real_t>(op_t, const real_t*, const cusparse::SpMatBase&, const cusparse::DnVec<real_t>&, const real_t*, cusparse::DnVec<real_t>&, cusparseSpMVAlg_t);
template void CuSparseHandler::performSpmv<real4_t>(op_t, const real4_t*, const cusparse::SpMatBase&, const cusparse::DnVec<real4_t>&, const real4_t*, cusparse::DnVec<real4_t>&, cusparseSpMVAlg_t);
template void CuSparseHandler::performSpmv<complex_t>(op_t, const complex_t*, const cusparse::SpMatBase&, const cusparse::DnVec<complex_t>&, const complex_t*, cusparse::DnVec<complex_t>&, cusparseSpMVAlg_t);
template void CuSparseHandler::performSpmv<complex8_t>(op_t, const complex8_t*, const cusparse::SpMatBase&, const cusparse::DnVec<complex8_t>&, const complex8_t*, cusparse::DnVec<complex8_t>&, cusparseSpMVAlg_t);
/*-------------------------------------------------*/
template <typename T_Scalar>
void CuSparseHandler::reserveSpmm(op_t                             opA,
                                  op_t                             opB,
                                  const T_Scalar*                  alpha,
                                  const cusparse::SpMatBase&       matA,
                                  const cusparse::DnMat<T_Scalar>& matB,
                                  const T_Scalar*                  beta,
                                  cusparse::DnMat<T_Scalar>&       matC,
                                  cusparseSpMMAlg_t                alg)
{
    cusparseStatus_t status = cusparseSpMM_bufferSize(handle(),
                                                      cusparse::cla3pOp2cusparseOp(opA),
                                                      cusparse::cla3pOp2cusparseOp(opB),
                                                      alpha,
                                                      matA.descr(),
                                                      matB.descr(),
                                                      beta,
                                                      matC.descr(),
                                                      TypeTraits<T_Scalar>::cuda_type(),
                                                      alg,
                                                      &m_workspaceInBytes);
    err::check_cusparse(status);

    deviceWork().reserve(m_workspaceInBytes);
}
//-------------------------------------------------
template void CuSparseHandler::reserveSpmm<real_t>(op_t, op_t, const real_t*, const cusparse::SpMatBase&, const cusparse::DnMat<real_t>&, const real_t*, cusparse::DnMat<real_t>&, cusparseSpMMAlg_t);
template void CuSparseHandler::reserveSpmm<real4_t>(op_t, op_t, const real4_t*, const cusparse::SpMatBase&, const cusparse::DnMat<real4_t>&, const real4_t*, cusparse::DnMat<real4_t>&, cusparseSpMMAlg_t);
template void CuSparseHandler::reserveSpmm<complex_t>(op_t, op_t, const complex_t*, const cusparse::SpMatBase&, const cusparse::DnMat<complex_t>&, const complex_t*, cusparse::DnMat<complex_t>&, cusparseSpMMAlg_t);
template void CuSparseHandler::reserveSpmm<complex8_t>(op_t, op_t, const complex8_t*, const cusparse::SpMatBase&, const cusparse::DnMat<complex8_t>&, const complex8_t*, cusparse::DnMat<complex8_t>&, cusparseSpMMAlg_t);
/*-------------------------------------------------*/
template <typename T_Scalar>
void CuSparseHandler::preprocessSpmm(op_t                             opA,
                                     op_t                             opB,
                                     const T_Scalar*                  alpha,
                                     const cusparse::SpMatBase&       matA,
                                     const cusparse::DnMat<T_Scalar>& matB,
                                     const T_Scalar*                  beta,
                                     cusparse::DnMat<T_Scalar>&       matC,
                                     cusparseSpMMAlg_t                alg)
{
    cusparseStatus_t status = cusparseSpMM_preprocess(handle(),
                                                      cusparse::cla3pOp2cusparseOp(opA),
                                                      cusparse::cla3pOp2cusparseOp(opB),
                                                      alpha,
                                                      matA.descr(),
                                                      matB.descr(),
                                                      beta,
                                                      matC.descr(),
                                                      TypeTraits<T_Scalar>::cuda_type(),
                                                      alg,
                                                      deviceWork().data());
    err::check_cusparse(status);
}
//-------------------------------------------------
template void CuSparseHandler::preprocessSpmm<real_t>(op_t, op_t, const real_t*, const cusparse::SpMatBase&, const cusparse::DnMat<real_t>&, const real_t*, cusparse::DnMat<real_t>&, cusparseSpMMAlg_t);
template void CuSparseHandler::preprocessSpmm<real4_t>(op_t, op_t, const real4_t*, const cusparse::SpMatBase&, const cusparse::DnMat<real4_t>&, const real4_t*, cusparse::DnMat<real4_t>&, cusparseSpMMAlg_t);
template void CuSparseHandler::preprocessSpmm<complex_t>(op_t, op_t, const complex_t*, const cusparse::SpMatBase&, const cusparse::DnMat<complex_t>&, const complex_t*, cusparse::DnMat<complex_t>&, cusparseSpMMAlg_t);
template void CuSparseHandler::preprocessSpmm<complex8_t>(op_t, op_t, const complex8_t*, const cusparse::SpMatBase&, const cusparse::DnMat<complex8_t>&, const complex8_t*, cusparse::DnMat<complex8_t>&, cusparseSpMMAlg_t);
/*-------------------------------------------------*/
template <typename T_Scalar>
void CuSparseHandler::performSpmm(op_t                             opA,
                                  op_t                             opB,
                                  const T_Scalar*                  alpha,
                                  const cusparse::SpMatBase&       matA,
                                  const cusparse::DnMat<T_Scalar>& matB,
                                  const T_Scalar*                  beta,
                                  cusparse::DnMat<T_Scalar>&       matC,
                                  cusparseSpMMAlg_t                alg)
{
    cusparseStatus_t status = cusparseSpMM(handle(),
                                           cusparse::cla3pOp2cusparseOp(opA),
                                           cusparse::cla3pOp2cusparseOp(opB),
                                           alpha,
                                           matA.descr(),
                                           matB.descr(),
                                           beta,
                                           matC.descr(),
                                           TypeTraits<T_Scalar>::cuda_type(),
                                           alg,
                                           deviceWork().data());
    err::check_cusparse(status);
}
//-------------------------------------------------
template void CuSparseHandler::performSpmm<real_t>(op_t, op_t, const real_t*, const cusparse::SpMatBase&, const cusparse::DnMat<real_t>&, const real_t*, cusparse::DnMat<real_t>&, cusparseSpMMAlg_t);
template void CuSparseHandler::performSpmm<real4_t>(op_t, op_t, const real4_t*, const cusparse::SpMatBase&, const cusparse::DnMat<real4_t>&, const real4_t*, cusparse::DnMat<real4_t>&, cusparseSpMMAlg_t);
template void CuSparseHandler::performSpmm<complex_t>(op_t, op_t, const complex_t*, const cusparse::SpMatBase&, const cusparse::DnMat<complex_t>&, const complex_t*, cusparse::DnMat<complex_t>&, cusparseSpMMAlg_t);
template void CuSparseHandler::performSpmm<complex8_t>(op_t, op_t, const complex8_t*, const cusparse::SpMatBase&, const cusparse::DnMat<complex8_t>&, const complex8_t*, cusparse::DnMat<complex8_t>&, cusparseSpMMAlg_t);
/*-------------------------------------------------*/
} // namespace culite
/*-------------------------------------------------*/
