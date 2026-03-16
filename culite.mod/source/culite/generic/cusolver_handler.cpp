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
#include <culite/generic/cusolver_handler.hpp>

// system

// 3rd

// culite
#include <culite/types/scalar.hpp>
#include <culite/error/cuda.hpp>
#include <culite/error/exceptions.hpp>
#include <culite/support/utils.hpp>
#include <culite/bulk/dns1d.hpp>
#include <culite/bulk/dns2d.hpp>
#include <culite/proxies/cusolver_proxy.hpp>
#include <culite/proxies/cublas_proxy.hpp>
#include <culite/dense.hpp>

// forwards
#include <culite/checks/cla3p_forwards.hpp>

/*-------------------------------------------------*/
namespace culite {
/*-------------------------------------------------*/
CuSolverHandler::InputMeta::InputMeta()
{
    defaults();
}
/*-------------------------------------------------*/
CuSolverHandler::InputMeta::InputMeta(cuSolverInt nrows, cuSolverInt ncols, cudaDataType_t cuda_type)
    : m_nrows(nrows), m_ncols(ncols), m_cuda_type(cuda_type)
{
}
/*-------------------------------------------------*/
CuSolverHandler::cuSolverInt CuSolverHandler::InputMeta::nrows() const
{
    return m_nrows;
}
/*-------------------------------------------------*/
CuSolverHandler::cuSolverInt CuSolverHandler::InputMeta::ncols() const
{
    return m_ncols;
}
/*-------------------------------------------------*/
cudaDataType_t CuSolverHandler::InputMeta::cuda_type() const
{
    return m_cuda_type;
}
/*-------------------------------------------------*/
void CuSolverHandler::InputMeta::clear()
{
    defaults();
}
/*-------------------------------------------------*/
void CuSolverHandler::InputMeta::defaults()
{
    m_nrows = 0;
    m_ncols = 0;
    m_cuda_type = CUDA_R_32F; // default to single precision real
}
/*-------------------------------------------------*/
CuSolverHandler::CuSolverHandler()
{
    defaults();
    cusolverStatus_t cusolverStatus = cusolverDnCreate(&m_handle);
    err::check_cusolver(cusolverStatus);
}
/*-------------------------------------------------*/
CuSolverHandler::~CuSolverHandler()
{
    clear();
    cusolverStatus_t cusolverStatus = cusolverDnDestroy(m_handle);
    err::check_cusolver(cusolverStatus);
}
/*-------------------------------------------------*/
void CuSolverHandler::defaults()
{
    m_workspaceInBytesOnDevice = 0;
    m_workspaceInBytesOnHost = 0;
}
/*-------------------------------------------------*/
CuSolverHandler::InputMeta& CuSolverHandler::inputMeta() { return m_inputMeta; }
DeviceBuffer<CuSolverHandler::cuSolverInt>& CuSolverHandler::ipiv() { return m_ipiv; }
DeviceBuffer<CuSolverHandler::infoInt>& CuSolverHandler::info() { return m_info; }
DeviceBufferVoid& CuSolverHandler::customWork() { return m_customBuffer; }
DeviceBufferVoid& CuSolverHandler::deviceWork() { return m_deviceBuffer; }
PinnedBufferVoid& CuSolverHandler::hostWork() { return m_hostBuffer; }
/*-------------------------------------------------*/
const CuSolverHandler::InputMeta& CuSolverHandler::inputMeta() const { return m_inputMeta; }
const DeviceBuffer<CuSolverHandler::cuSolverInt>& CuSolverHandler::ipiv() const { return m_ipiv; }
const DeviceBuffer<CuSolverHandler::infoInt>& CuSolverHandler::info() const { return m_info; }
const DeviceBufferVoid& CuSolverHandler::customWork() const { return m_customBuffer; }
const DeviceBufferVoid& CuSolverHandler::deviceWork() const { return m_deviceBuffer; }
const PinnedBufferVoid& CuSolverHandler::hostWork() const { return m_hostBuffer; }
/*-------------------------------------------------*/
void CuSolverHandler::clear()
{
    inputMeta().clear();
    ipiv().clear();
    info().clear();
    customWork().clear();
    deviceWork().clear();
    hostWork().clear();

    defaults();
}
/*-------------------------------------------------*/
cusolverDnHandle_t CuSolverHandler::handle()
{
    return m_handle;
}
/*-------------------------------------------------*/
cusolverDnParams_t CuSolverHandler::params() const
{
    return m_params.get();
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void CuSolverHandler::reserveLU(const T_Matrix& A)
{
    using T_Scalar = typename T_Matrix::value_type;

    cusolverStatus_t cusolverStatus =
    cusolverDnXgetrf_bufferSize(handle(),
                                params(),
                                A.nrows(), A.ncols(),
                                TypeTraits<T_Scalar>::cuda_type(), A.values(), A.ld(),
                                TypeTraits<T_Scalar>::cuda_type(),
                                &m_workspaceInBytesOnDevice,
                                &m_workspaceInBytesOnHost);

    err::check_cusolver(cusolverStatus);

    ipiv().reserve(std::min(A.nrows(), A.ncols()));
    info().reserve(1);
    customWork().reserve(A.nrows() * A.ncols() * sizeof(T_Scalar));
    deviceWork().reserve(m_workspaceInBytesOnDevice);
    hostWork().reserve(m_workspaceInBytesOnHost);
}
/*-------------------------------------------------*/
template void CuSolverHandler::reserveLU<dns::RdMatrix>(const dns::RdMatrix&);
template void CuSolverHandler::reserveLU<dns::RfMatrix>(const dns::RfMatrix&);
template void CuSolverHandler::reserveLU<dns::CdMatrix>(const dns::CdMatrix&);
template void CuSolverHandler::reserveLU<dns::CfMatrix>(const dns::CfMatrix&);
/*-------------------------------------------------*/
template <typename T_Matrix>
void CuSolverHandler::decomposeLU(const T_Matrix& A)
{
    using T_Scalar = typename T_Matrix::value_type;

    cuSolverInt n = A.nrows();

    memCopyD2D<T_Scalar>(A.nrows(), A.ncols(),
                         A.values(), A.ld(),
                         static_cast<T_Scalar*>(customWork().data()), n);

    cusolverStatus_t cusolverStatus =
    cusolverDnXgetrf(handle(),
                     params(),
                     A.nrows(), A.ncols(),
                     TypeTraits<T_Scalar>::cuda_type(), customWork().data(), n,
                     ipiv().data(),
                     TypeTraits<T_Scalar>::cuda_type(),
                     deviceWork().data(), m_workspaceInBytesOnDevice,
                     hostWork().data(), m_workspaceInBytesOnHost,
                     info().data());

    err::check_cusolver(cusolverStatus);

    inputMeta() = InputMeta(A.nrows(), A.ncols(), TypeTraits<T_Scalar>::cuda_type());
}
/*-------------------------------------------------*/
template void CuSolverHandler::decomposeLU<dns::RdMatrix>(const dns::RdMatrix&);
template void CuSolverHandler::decomposeLU<dns::RfMatrix>(const dns::RfMatrix&);
template void CuSolverHandler::decomposeLU<dns::CdMatrix>(const dns::CdMatrix&);
template void CuSolverHandler::decomposeLU<dns::CfMatrix>(const dns::CfMatrix&);
/*-------------------------------------------------*/
template <typename T_Matrix>
void CuSolverHandler::solveLU(T_Matrix& B)
{
    using T_Scalar = typename T_Matrix::value_type;

    cuSolverInt n = B.nrows();
    ::cla3p::similarity_dim_check(inputMeta().nrows(), n);

    cusolverStatus_t cusolverStatus =
    cusolverDnXgetrs(handle(),
                     params(),
                     CUBLAS_OP_N,
                     n, B.ncols(),
                     inputMeta().cuda_type(), customWork().data(), n,
                     ipiv().data(),
                     TypeTraits<T_Scalar>::cuda_type(), B.values(), B.ld(),
                     info().data());

    err::check_cusolver(cusolverStatus);
}
/*-------------------------------------------------*/
template void CuSolverHandler::solveLU<dns::RdMatrix>(dns::RdMatrix&);
template void CuSolverHandler::solveLU<dns::RfMatrix>(dns::RfMatrix&);
template void CuSolverHandler::solveLU<dns::CdMatrix>(dns::CdMatrix&);
template void CuSolverHandler::solveLU<dns::CfMatrix>(dns::CfMatrix&);
/*-------------------------------------------------*/
template <typename T_Matrix>
void CuSolverHandler::reserveGeev(const T_Matrix& A, bool calcLeft, bool calcRight)
{
    using T_Scalar = typename T_Matrix::value_type;
    bool realCase = TypeTraits<T_Scalar>::is_real();

    cuSolverInt n = A.nrows();
    cuSolverInt sizeW = (realCase ? 2 * n : n) * sizeof(T_Scalar);
    cuSolverInt sizeN2 = n * n * sizeof(T_Scalar);
    cuSolverInt sizevA = sizeN2;
    cuSolverInt sizeVL = calcLeft ? sizeN2 : 0;
    cuSolverInt sizeVR = calcRight ? sizeN2 : 0;

    customWork().reserve(sizeW + sizeVL + sizeVR + sizevA);

    T_Scalar *W  = nullptr;
    T_Scalar *VL = nullptr;
    T_Scalar *VR = nullptr;
    geevAssignInternalPointers<T_Scalar>(calcLeft, calcRight, n, &W, nullptr, &VL, &VR);

    cusolverStatus_t cusolverStatus =
    cusolverDnXgeev_bufferSize(handle(),
                               params(),
                               cusolver::bool2cusolverEigMode(calcLeft),
                               cusolver::bool2cusolverEigMode(calcRight),
                               n,
                               TypeTraits<T_Scalar>::cuda_type(), A.values(), A.ld(),
                               TypeTraits<T_Scalar>::cuda_type(), W,
                               TypeTraits<T_Scalar>::cuda_type(), VL, n,
                               TypeTraits<T_Scalar>::cuda_type(), VR, n,
                               TypeTraits<T_Scalar>::cuda_type(),
                               &m_workspaceInBytesOnDevice,
                               &m_workspaceInBytesOnHost);

    err::check_cusolver(cusolverStatus);

    info().reserve(1);
    deviceWork().reserve(m_workspaceInBytesOnDevice);
    hostWork().reserve(m_workspaceInBytesOnHost);
}
/*-------------------------------------------------*/
template void CuSolverHandler::reserveGeev<dns::RdMatrix>(const dns::RdMatrix&, bool, bool);
template void CuSolverHandler::reserveGeev<dns::RfMatrix>(const dns::RfMatrix&, bool, bool);
template void CuSolverHandler::reserveGeev<dns::CdMatrix>(const dns::CdMatrix&, bool, bool);
template void CuSolverHandler::reserveGeev<dns::CfMatrix>(const dns::CfMatrix&, bool, bool);
/*-------------------------------------------------*/
template <typename T_Matrix>
void CuSolverHandler::executeGeev(const T_Matrix& A, bool calcLeft, bool calcRight)
{
    using T_Scalar = typename T_Matrix::value_type;

    cuSolverInt n = A.nrows();

    T_Scalar *W  = nullptr;
    T_Scalar *vA = nullptr;
    T_Scalar *VL = nullptr;
    T_Scalar *VR = nullptr;
    geevAssignInternalPointers<T_Scalar>(calcLeft, calcRight, n, &W, &vA, &VL, &VR);

    memCopyD2D<T_Scalar>(n, n, A.values(), A.ld(), vA, n);

    cusolverStatus_t cusolverStatus =
    cusolverDnXgeev(handle(),
                    params(),
                    cusolver::bool2cusolverEigMode(calcLeft),
                    cusolver::bool2cusolverEigMode(calcRight),
                    n,
                    TypeTraits<T_Scalar>::cuda_type(), vA, n, // vA is overwritten
                    TypeTraits<T_Scalar>::cuda_type(), W,
                    TypeTraits<T_Scalar>::cuda_type(), VL, n,
                    TypeTraits<T_Scalar>::cuda_type(), VR, n,
                    TypeTraits<T_Scalar>::cuda_type(),
                    deviceWork().data(), m_workspaceInBytesOnDevice,
                    hostWork().data(), m_workspaceInBytesOnHost,
                    info().data());

    err::check_cusolver(cusolverStatus);

    inputMeta() = InputMeta(A.nrows(), A.ncols(), TypeTraits<T_Scalar>::cuda_type());
}
/*-------------------------------------------------*/
template void CuSolverHandler::executeGeev<dns::RdMatrix>(const dns::RdMatrix&, bool, bool);
template void CuSolverHandler::executeGeev<dns::RfMatrix>(const dns::RfMatrix&, bool, bool);
template void CuSolverHandler::executeGeev<dns::CdMatrix>(const dns::CdMatrix&, bool, bool);
template void CuSolverHandler::executeGeev<dns::CfMatrix>(const dns::CfMatrix&, bool, bool);
/*-------------------------------------------------*/
template <typename T_Vector>
void CuSolverHandler::geevGetEigenvalues(T_Vector& eigs)
{
    using T_Scalar = typename T_Vector::value_type;

    if(TypeTraits<T_Scalar>::is_real()) {
        throw err::CudaException("Geev only returns complex eigenpairs.");
    }

    int_t n = inputMeta().nrows();

    if(!eigs) eigs = T_Vector(n);
    ::cla3p::similarity_dim_check(eigs.size(), n);

    if(inputMeta().cuda_type() == cudaDataType::CUDA_R_32F || inputMeta().cuda_type() == cudaDataType::CUDA_R_64F) {

        using T_RScalar = typename TypeTraits<T_Scalar>::real_type;

        T_RScalar *W = nullptr;
        geevAssignInternalPointers<T_RScalar>(false, false, n, &W, nullptr, nullptr, nullptr);

        blk::dns::geevCalculateComplexEigenvalues(n, W, eigs.values());

    } else if(inputMeta().cuda_type() == cudaDataType::CUDA_C_32F || inputMeta().cuda_type() == cudaDataType::CUDA_C_64F) {

        T_Scalar *W = nullptr;
        geevAssignInternalPointers<T_Scalar>(false, false, n, &W, nullptr, nullptr, nullptr);

        memCopyD2D<T_Scalar>(n, W, eigs.values());

    } else {

        throw err::CudaException("Unsupported data type for eigenvalue retrieval.");

    } // problem cuda type
}
/*-------------------------------------------------*/
template void CuSolverHandler::geevGetEigenvalues<dns::CdVector>(dns::CdVector&);
template void CuSolverHandler::geevGetEigenvalues<dns::CfVector>(dns::CfVector&);
/*-------------------------------------------------*/
template <typename T_Matrix>
void CuSolverHandler::geevGetEigenvectors(bool calcLeft, bool calcRight, T_Matrix& leftEigenvectors, T_Matrix& rightEigenvectors)
{
    using T_Scalar = typename T_Matrix::value_type;

    if(TypeTraits<T_Scalar>::is_real()) {
        throw err::CudaException("Geev only returns complex eigenpairs.");
    }

    int_t n = inputMeta().nrows();

    if(calcLeft) {
        if(!leftEigenvectors) leftEigenvectors = T_Matrix(n, n);
        ::cla3p::similarity_dim_check(leftEigenvectors.nrows(), n);
        ::cla3p::similarity_dim_check(leftEigenvectors.ncols(), n);
    }

    if(calcRight) {
        if(!rightEigenvectors) rightEigenvectors = T_Matrix(n, n);
        ::cla3p::similarity_dim_check(rightEigenvectors.nrows(), n);
        ::cla3p::similarity_dim_check(rightEigenvectors.ncols(), n);
    }

    if(inputMeta().cuda_type() == cudaDataType::CUDA_R_32F || inputMeta().cuda_type() == cudaDataType::CUDA_R_64F) {

        using T_RScalar = typename TypeTraits<T_Scalar>::real_type;

        T_RScalar *W  = nullptr;
        T_RScalar *VL = nullptr;
        T_RScalar *VR = nullptr;
        geevAssignInternalPointers<T_RScalar>(calcLeft, calcRight, n, &W, nullptr, &VL, &VR);

        if(calcLeft) {
            blk::dns::geevCalculateComplexEigenvectors(n, W, VL, n,
                                                       leftEigenvectors.values(),
                                                       leftEigenvectors.ld());
        } // calcLeft

        if(calcRight) {
            blk::dns::geevCalculateComplexEigenvectors(n, W, VR, n,
                                                       rightEigenvectors.values(),
                                                       rightEigenvectors.ld());
        } // calcRight

    } else if(inputMeta().cuda_type() == cudaDataType::CUDA_C_32F || inputMeta().cuda_type() == cudaDataType::CUDA_C_64F) {

        T_Scalar *VL = nullptr;
        T_Scalar *VR = nullptr;
        geevAssignInternalPointers<T_Scalar>(calcLeft, calcRight, n, nullptr, nullptr, &VL, &VR);

        if(calcLeft)  memCopyD2D<T_Scalar>(n, n, VL, n, leftEigenvectors.values(), leftEigenvectors.ld());
        if(calcRight) memCopyD2D<T_Scalar>(n, n, VR, n, rightEigenvectors.values(), rightEigenvectors.ld());

    } else {

        throw err::CudaException("Unsupported data type for eigenvector retrieval.");

    } // problem cuda type
}
/*-------------------------------------------------*/
template void CuSolverHandler::geevGetEigenvectors<dns::CdMatrix>(bool, bool, dns::CdMatrix&, dns::CdMatrix&);
template void CuSolverHandler::geevGetEigenvectors<dns::CfMatrix>(bool, bool, dns::CfMatrix&, dns::CfMatrix&);
/*-------------------------------------------------*/
template <typename T_Matrix>
void CuSolverHandler::reserveSyevd(const T_Matrix& A, bool calcVectors)
{
    using T_Scalar = typename T_Matrix::value_type;
    using T_RScalar = typename TypeTraits<T_Scalar>::real_type;

    cuSolverInt n = A.nrows();
    cuSolverInt sizeW = n * sizeof(T_RScalar);
    cuSolverInt sizeV = n * n * sizeof(T_Scalar);

    customWork().reserve(sizeW + sizeV);

    T_RScalar *W = nullptr;
    T_Scalar *V = nullptr;
    syevdAssignInternalPointers<T_Scalar>(n, &W, &V);

    cusolverStatus_t cusolverStatus =
    cusolverDnXsyevd_bufferSize(handle(),
                                params(),
                                cusolver::bool2cusolverEigMode(calcVectors),
                                cublas::cla3pUplo2cublasUplo(A.prop().uplo()),
                                n,
                                TypeTraits<T_Scalar>::cuda_type(), A.values(), A.ld(),
                                TypeTraits<T_RScalar>::cuda_type(), W,
                                TypeTraits<T_Scalar>::cuda_type(),
                                &m_workspaceInBytesOnDevice,
                                &m_workspaceInBytesOnHost);

    err::check_cusolver(cusolverStatus);

    info().reserve(1);
    deviceWork().reserve(m_workspaceInBytesOnDevice);
    hostWork().reserve(m_workspaceInBytesOnHost);
}
/*-------------------------------------------------*/
template void CuSolverHandler::reserveSyevd<dns::RdMatrix>(const dns::RdMatrix&, bool);
template void CuSolverHandler::reserveSyevd<dns::RfMatrix>(const dns::RfMatrix&, bool);
template void CuSolverHandler::reserveSyevd<dns::CdMatrix>(const dns::CdMatrix&, bool);
template void CuSolverHandler::reserveSyevd<dns::CfMatrix>(const dns::CfMatrix&, bool);
/*-------------------------------------------------*/
template <typename T_Matrix>
void CuSolverHandler::executeSyevd(const T_Matrix& A, bool calcVectors)
{
    using T_Scalar = typename T_Matrix::value_type;
    using T_RScalar = typename TypeTraits<T_Scalar>::real_type;

    cuSolverInt n = A.nrows();

    T_RScalar *W = nullptr;
    T_Scalar *V = nullptr;
    syevdAssignInternalPointers<T_Scalar>(n, &W, &V);

    memCopyD2D<T_Scalar>(n, n, A.values(), A.ld(), V, n);

    cusolverStatus_t cusolverStatus =
    cusolverDnXsyevd(handle(),
                     params(),
                     cusolver::bool2cusolverEigMode(calcVectors),
                     cublas::cla3pUplo2cublasUplo(A.prop().uplo()),
                     n,
                     TypeTraits<T_Scalar>::cuda_type(), V, n,
                     TypeTraits<T_RScalar>::cuda_type(), W,
                     TypeTraits<T_Scalar>::cuda_type(),
                     deviceWork().data(), m_workspaceInBytesOnDevice,
                     hostWork().data(), m_workspaceInBytesOnHost,
                     info().data());

    err::check_cusolver(cusolverStatus);

    inputMeta() = InputMeta(A.nrows(), A.ncols(), TypeTraits<T_Scalar>::cuda_type());
}
/*-------------------------------------------------*/
template void CuSolverHandler::executeSyevd<dns::RdMatrix>(const dns::RdMatrix&, bool);
template void CuSolverHandler::executeSyevd<dns::RfMatrix>(const dns::RfMatrix&, bool);
template void CuSolverHandler::executeSyevd<dns::CdMatrix>(const dns::CdMatrix&, bool);
template void CuSolverHandler::executeSyevd<dns::CfMatrix>(const dns::CfMatrix&, bool);
/*-------------------------------------------------*/
template <typename T_Matrix>
void CuSolverHandler::reserveSyevdx(const T_Matrix& A, bool calcVectors, 
                                    eigRange_t range, 
                                    cuSolverInt il, 
                                    cuSolverInt iu, 
                                    typename TypeTraits<typename T_Matrix::value_type>::real_type vl, 
                                    typename TypeTraits<typename T_Matrix::value_type>::real_type vu)
{
    using T_Scalar = typename T_Matrix::value_type;
    using T_RScalar = typename TypeTraits<T_Scalar>::real_type;

    cuSolverInt n = A.nrows();
    cuSolverInt sizeW = n * sizeof(T_RScalar);
    cuSolverInt sizeV = n * n * sizeof(T_Scalar);

    customWork().reserve(sizeW + sizeV);

    T_RScalar *W = nullptr;
    T_Scalar *V = nullptr;
    syevdAssignInternalPointers<T_Scalar>(n, &W, &V);

    cuSolverInt eigsFound = 0;

    cusolverStatus_t cusolverStatus =
    cusolverDnXsyevdx_bufferSize(handle(),
                                 params(),
                                 cusolver::bool2cusolverEigMode(calcVectors),
                                 cusolver::cla3pEigRange2cusolverEigRange(range),
                                 cublas::cla3pUplo2cublasUplo(A.prop().uplo()),
                                 n,
                                 TypeTraits<T_Scalar>::cuda_type(), A.values(), A.ld(),
                                 &vl, &vu, il, iu, &eigsFound,
                                 TypeTraits<T_RScalar>::cuda_type(), W,
                                 TypeTraits<T_Scalar>::cuda_type(),
                                 &m_workspaceInBytesOnDevice,
                                 &m_workspaceInBytesOnHost);

    err::check_cusolver(cusolverStatus);

    info().reserve(1);
    deviceWork().reserve(m_workspaceInBytesOnDevice);
    hostWork().reserve(m_workspaceInBytesOnHost);
}
/*-------------------------------------------------*/
template void CuSolverHandler::reserveSyevdx<dns::RdMatrix>(const dns::RdMatrix&, bool, eigRange_t, cuSolverInt, cuSolverInt, real_t, real_t);
template void CuSolverHandler::reserveSyevdx<dns::RfMatrix>(const dns::RfMatrix&, bool, eigRange_t, cuSolverInt, cuSolverInt, real4_t, real4_t);
template void CuSolverHandler::reserveSyevdx<dns::CdMatrix>(const dns::CdMatrix&, bool, eigRange_t, cuSolverInt, cuSolverInt, real_t, real_t);
template void CuSolverHandler::reserveSyevdx<dns::CfMatrix>(const dns::CfMatrix&, bool, eigRange_t, cuSolverInt, cuSolverInt, real4_t, real4_t);
/*-------------------------------------------------*/
template <typename T_Matrix>
int_t CuSolverHandler::executeSyevdx(const T_Matrix& A, bool calcVectors, 
                                     eigRange_t range, 
                                     cuSolverInt il, 
                                     cuSolverInt iu, 
                                     typename TypeTraits<typename T_Matrix::value_type>::real_type vl, 
                                     typename TypeTraits<typename T_Matrix::value_type>::real_type vu)
{
    using T_Scalar = typename T_Matrix::value_type;
    using T_RScalar = typename TypeTraits<T_Scalar>::real_type;

    cuSolverInt n = A.nrows();

    T_RScalar *W = nullptr;
    T_Scalar *V = nullptr;
    syevdAssignInternalPointers<T_Scalar>(n, &W, &V);

    memCopyD2D<T_Scalar>(n, n, A.values(), A.ld(), V, n);

    cuSolverInt eigsFound = 0;

    cusolverStatus_t cusolverStatus =
    cusolverDnXsyevdx(handle(),
                      params(),
                      cusolver::bool2cusolverEigMode(calcVectors),
                      cusolver::cla3pEigRange2cusolverEigRange(range),
                      cublas::cla3pUplo2cublasUplo(A.prop().uplo()),
                      n,
                      TypeTraits<T_Scalar>::cuda_type(), V, n,
                      &vl, &vu, il, iu, &eigsFound,
                      TypeTraits<T_RScalar>::cuda_type(), W,
                      TypeTraits<T_Scalar>::cuda_type(),
                      deviceWork().data(), m_workspaceInBytesOnDevice,
                      hostWork().data(), m_workspaceInBytesOnHost,
                      info().data());

    err::check_cusolver(cusolverStatus);

    inputMeta() = InputMeta(A.nrows(), A.ncols(), TypeTraits<T_Scalar>::cuda_type());

    return eigsFound;
}
/*-------------------------------------------------*/
template int_t CuSolverHandler::executeSyevdx<dns::RdMatrix>(const dns::RdMatrix&, bool, eigRange_t, cuSolverInt, cuSolverInt, real_t, real_t);
template int_t CuSolverHandler::executeSyevdx<dns::RfMatrix>(const dns::RfMatrix&, bool, eigRange_t, cuSolverInt, cuSolverInt, real4_t, real4_t);
template int_t CuSolverHandler::executeSyevdx<dns::CdMatrix>(const dns::CdMatrix&, bool, eigRange_t, cuSolverInt, cuSolverInt, real_t, real_t);
template int_t CuSolverHandler::executeSyevdx<dns::CfMatrix>(const dns::CfMatrix&, bool, eigRange_t, cuSolverInt, cuSolverInt, real4_t, real4_t);
/*-------------------------------------------------*/
template <typename T_Matrix>
void CuSolverHandler::reserveGesvd(const T_Matrix& A, svdPolicy_t policyU, svdPolicy_t policyVT)
{
    using T_Scalar = typename T_Matrix::value_type;
    using T_RScalar = typename TypeTraits<T_Scalar>::real_type;

    cuSolverInt m = A.nrows();
    cuSolverInt n = A.ncols();
    cuSolverInt k = std::min(m, n);

    // TODO: If allocate with Real size, cuSOLVER returns 
    //       CUDA_ERROR_MISALIGNED_ADDRESS error (for complex types only)
    // cuSolverInt sizeS = k * sizeof(T_RScalar);
    cuSolverInt sizeS = k * sizeof(T_Scalar);
    cuSolverInt sizeA = m * n * sizeof(T_Scalar);
    cuSolverInt sizeU = svdVectorSize(policyU, m, k) * sizeof(T_Scalar);
    cuSolverInt sizeV = svdVectorSize(policyVT, n, k) * sizeof(T_Scalar);

    customWork().reserve(sizeS + sizeA + sizeU + sizeV);

    T_RScalar *S = nullptr;
    T_Scalar *U = nullptr;
    T_Scalar *VT = nullptr;
    gesvdAssignInternalPointers<T_Scalar>(policyU, policyVT, m, n, &S, nullptr, &U, &VT);

    cuSolverInt ldu = m;
    cuSolverInt ldvt = (policyVT == svdPolicy_t::Limited ? k : n);

    cusolverStatus_t cusolverStatus =
    cusolverDnXgesvd_bufferSize(handle(),
                                params(),
                                static_cast<signed char>(policyU),
                                static_cast<signed char>(policyVT),
                                m, n,
                                TypeTraits<T_Scalar>::cuda_type(), A.values(), A.ld(),
                                TypeTraits<T_RScalar>::cuda_type(), S,
                                TypeTraits<T_Scalar>::cuda_type(), U, ldu,
                                TypeTraits<T_Scalar>::cuda_type(), VT, ldvt,
                                TypeTraits<T_Scalar>::cuda_type(),
                                &m_workspaceInBytesOnDevice,
                                &m_workspaceInBytesOnHost);

    err::check_cusolver(cusolverStatus);

    info().reserve(1);
    deviceWork().reserve(m_workspaceInBytesOnDevice);
    hostWork().reserve(m_workspaceInBytesOnHost);
}
/*-------------------------------------------------*/
template void CuSolverHandler::reserveGesvd<dns::RdMatrix>(const dns::RdMatrix&, svdPolicy_t, svdPolicy_t);
template void CuSolverHandler::reserveGesvd<dns::RfMatrix>(const dns::RfMatrix&, svdPolicy_t, svdPolicy_t);
template void CuSolverHandler::reserveGesvd<dns::CdMatrix>(const dns::CdMatrix&, svdPolicy_t, svdPolicy_t);
template void CuSolverHandler::reserveGesvd<dns::CfMatrix>(const dns::CfMatrix&, svdPolicy_t, svdPolicy_t);
/*-------------------------------------------------*/
template <typename T_Matrix>
void CuSolverHandler::executeGesvd(const T_Matrix& A, svdPolicy_t policyU, svdPolicy_t policyVT)
{
    using T_Scalar = typename T_Matrix::value_type;
    using T_RScalar = typename TypeTraits<T_Scalar>::real_type;

    cuSolverInt m = A.nrows();
    cuSolverInt n = A.ncols();
    cuSolverInt k = std::min(m, n);

    T_RScalar *S = nullptr;
    T_Scalar *vA = nullptr;
    T_Scalar *U = nullptr;
    T_Scalar *VT = nullptr;
    gesvdAssignInternalPointers<T_Scalar>(policyU, policyVT, m, n, &S, &vA, &U, &VT);

    cuSolverInt lda = m;
    cuSolverInt ldu = m;
    cuSolverInt ldvt = (policyVT == svdPolicy_t::Limited ? k : n);

    memCopyD2D<T_Scalar>(m, n, A.values(), A.ld(), vA, lda);

    cusolverStatus_t cusolverStatus =
    cusolverDnXgesvd(handle(),
                     params(),
                     static_cast<signed char>(policyU),
                     static_cast<signed char>(policyVT),
                     m, n,
                     TypeTraits<T_Scalar>::cuda_type(), vA, lda,
                     TypeTraits<T_RScalar>::cuda_type(), S,
                     TypeTraits<T_Scalar>::cuda_type(), U, ldu,
                     TypeTraits<T_Scalar>::cuda_type(), VT, ldvt,
                     TypeTraits<T_Scalar>::cuda_type(),
                     deviceWork().data(), m_workspaceInBytesOnDevice,
                     hostWork().data(), m_workspaceInBytesOnHost,
                     info().data());

    err::check_cusolver(cusolverStatus);

    inputMeta() = InputMeta(A.nrows(), A.ncols(), TypeTraits<T_Scalar>::cuda_type());
}
/*-------------------------------------------------*/
template void CuSolverHandler::executeGesvd<dns::RdMatrix>(const dns::RdMatrix&, svdPolicy_t, svdPolicy_t);
template void CuSolverHandler::executeGesvd<dns::RfMatrix>(const dns::RfMatrix&, svdPolicy_t, svdPolicy_t);
template void CuSolverHandler::executeGesvd<dns::CdMatrix>(const dns::CdMatrix&, svdPolicy_t, svdPolicy_t);
template void CuSolverHandler::executeGesvd<dns::CfMatrix>(const dns::CfMatrix&, svdPolicy_t, svdPolicy_t);
/*-------------------------------------------------*/
template <typename T_Vector>
void CuSolverHandler::gesvdGetSingularValues(T_Vector& sigma)
{
    using T_Scalar = typename T_Vector::value_type;

    if(!TypeTraits<T_Scalar>::is_real()) {
        throw err::CudaException("Gesvd only returns real singular values.");
    }

    int_t m = inputMeta().nrows();
    int_t n = inputMeta().ncols();
    int_t k = std::min(m, n);

    if(!sigma) sigma = T_Vector(k);
    ::cla3p::similarity_dim_check(sigma.size(), k);

    T_Scalar *S = nullptr;
    gesvdAssignInternalPointers<T_Scalar>(svdPolicy_t::NoCalculation,
                                          svdPolicy_t::NoCalculation,
                                          m, n, &S, nullptr, nullptr, nullptr);

    memCopyD2D<T_Scalar>(k, S, sigma.values());
}
/*-------------------------------------------------*/
template void CuSolverHandler::gesvdGetSingularValues<dns::RdVector>(dns::RdVector&);
template void CuSolverHandler::gesvdGetSingularValues<dns::RfVector>(dns::RfVector&);
/*-------------------------------------------------*/
template <typename T_Matrix>
void CuSolverHandler::gesvdGetSingularVectors(svdPolicy_t policyU, svdPolicy_t policyVT, T_Matrix& U, T_Matrix& V, bool transposeVT)
{
    using T_Scalar = typename T_Matrix::value_type;

    int_t m = inputMeta().nrows();
    int_t n = inputMeta().ncols();
    int_t k = std::min(m, n);

    int_t nrowsU = 0;
    int_t ncolsU = 0;
    int_t nrowsVT = 0;
    int_t ncolsVT = 0;

    if(policyU != svdPolicy_t::NoCalculation) {
        nrowsU = m;
        ncolsU = (policyU == svdPolicy_t::Limited ? k : m);
        if(!U) U = T_Matrix(nrowsU, ncolsU);
        ::cla3p::similarity_dim_check(U.nrows(), nrowsU);
        ::cla3p::similarity_dim_check(U.ncols(), ncolsU);
    }

    if(policyVT != svdPolicy_t::NoCalculation) {
        nrowsVT = (policyVT == svdPolicy_t::Limited ? k : n);
        ncolsVT = n;
        int_t nrowsV = transposeVT ? ncolsVT : nrowsVT;
        int_t ncolsV = transposeVT ? nrowsVT : ncolsVT;
        if(!V) V = T_Matrix(nrowsV, ncolsV);
        ::cla3p::similarity_dim_check(V.nrows(), nrowsV);
        ::cla3p::similarity_dim_check(V.ncols(), ncolsV);
    }

    T_Scalar *storedU = nullptr;
    T_Scalar *storedVT = nullptr;
    gesvdAssignInternalPointers<T_Scalar>(policyU, policyVT, m, n, nullptr, nullptr, &storedU, &storedVT);

    if(policyU != svdPolicy_t::NoCalculation) {
        memCopyD2D<T_Scalar>(nrowsU, ncolsU, storedU, nrowsU, U.values(), U.ld());
    } // policyU

    if(policyVT != svdPolicy_t::NoCalculation) {
        if(transposeVT) {
            blk::dns::ctranspose<T_Scalar>(nrowsVT, ncolsVT, storedVT, nrowsVT, V.values(), V.ld());
        } else {
            memCopyD2D<T_Scalar>(nrowsVT, ncolsVT, storedVT, nrowsVT, V.values(), V.ld());
        }
    } // policyVT
}
/*-------------------------------------------------*/
template void CuSolverHandler::gesvdGetSingularVectors<dns::RdMatrix>(svdPolicy_t, svdPolicy_t, dns::RdMatrix&, dns::RdMatrix&, bool);
template void CuSolverHandler::gesvdGetSingularVectors<dns::RfMatrix>(svdPolicy_t, svdPolicy_t, dns::RfMatrix&, dns::RfMatrix&, bool);
template void CuSolverHandler::gesvdGetSingularVectors<dns::CdMatrix>(svdPolicy_t, svdPolicy_t, dns::CdMatrix&, dns::CdMatrix&, bool);
template void CuSolverHandler::gesvdGetSingularVectors<dns::CfMatrix>(svdPolicy_t, svdPolicy_t, dns::CfMatrix&, dns::CfMatrix&, bool);
/*-------------------------------------------------*/
CuSolverHandler::cuSolverInt CuSolverHandler::svdVectorSize(svdPolicy_t policy, cuSolverInt n, cuSolverInt k) const
{
    switch(policy) {
        case svdPolicy_t::Full         : return n * n;
        case svdPolicy_t::Limited      : return n * std::min(n, k);
        case svdPolicy_t::NoCalculation: return 0;
        default: return 0;
    }
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void CuSolverHandler::geevAssignInternalPointers(bool calcLeft, bool calcRight, int_t n, T_Scalar** W, T_Scalar** A, T_Scalar** VL, T_Scalar** VR)
{
    bool realCase = TypeTraits<T_Scalar>::is_real();

    if(W)  *W  = nullptr;
    if(A)  *A  = nullptr;
    if(VL) *VL = nullptr;
    if(VR) *VR = nullptr;

    char *charBuffer = static_cast<char*>(customWork().data());

    std::size_t sizeW  = (realCase ? 2 * n : n) * sizeof(T_Scalar);
    std::size_t sizeN2 = n * n * sizeof(T_Scalar);
    std::size_t sizeA  = sizeN2;
    std::size_t sizeVL = calcLeft ? sizeN2 : 0;
    std::size_t sizeVR = calcRight ? sizeN2 : 0;

    if(W)               { *W  = reinterpret_cast<T_Scalar*>(charBuffer); } charBuffer += sizeW ;
    if(A)               { *A  = reinterpret_cast<T_Scalar*>(charBuffer); } charBuffer += sizeA ;
    if(calcLeft  && VL) { *VL = reinterpret_cast<T_Scalar*>(charBuffer); } charBuffer += sizeVL;
    if(calcRight && VR) { *VR = reinterpret_cast<T_Scalar*>(charBuffer); } charBuffer += sizeVR;
}
/*-------------------------------------------------*/
template void CuSolverHandler::geevAssignInternalPointers<real_t>(bool, bool, int_t, real_t**, real_t**, real_t**, real_t**);
template void CuSolverHandler::geevAssignInternalPointers<real4_t>(bool, bool, int_t, real4_t**, real4_t**, real4_t**, real4_t**);
template void CuSolverHandler::geevAssignInternalPointers<complex_t>(bool, bool, int_t, complex_t**, complex_t**, complex_t**, complex_t**);
template void CuSolverHandler::geevAssignInternalPointers<complex8_t>(bool, bool, int_t, complex8_t**, complex8_t**, complex8_t**, complex8_t**);
/*-------------------------------------------------*/
template <typename T_Scalar>
void CuSolverHandler::syevdAssignInternalPointers(int_t n, typename TypeTraits<T_Scalar>::real_type** W, T_Scalar** V)
{
    using T_RScalar = typename TypeTraits<T_Scalar>::real_type;

    if(W) *W = nullptr;
    if(V) *V = nullptr;

    char *charBuffer = static_cast<char*>(customWork().data());

    std::size_t sizeW = n * sizeof(T_RScalar);
    std::size_t sizeV = n * n * sizeof(T_Scalar);

    if(W) { *W = reinterpret_cast<T_RScalar*>(charBuffer); } charBuffer += sizeW;
    if(V) { *V = reinterpret_cast<T_Scalar*>(charBuffer); } charBuffer += sizeV;
}
/*-------------------------------------------------*/
template void CuSolverHandler::syevdAssignInternalPointers<real_t>(int_t, real_t**, real_t**);
template void CuSolverHandler::syevdAssignInternalPointers<real4_t>(int_t, real4_t**, real4_t**);
template void CuSolverHandler::syevdAssignInternalPointers<complex_t>(int_t, real_t**, complex_t**);
template void CuSolverHandler::syevdAssignInternalPointers<complex8_t>(int_t, real4_t**, complex8_t**);
/*-------------------------------------------------*/
template <typename T_Scalar>
void CuSolverHandler::gesvdAssignInternalPointers(svdPolicy_t jobu, svdPolicy_t jobv, int_t m, int_t n, typename TypeTraits<T_Scalar>::real_type** S, T_Scalar** A, T_Scalar** U, T_Scalar** V)
{
    using T_RScalar = typename TypeTraits<T_Scalar>::real_type;

    int_t k = std::min(m, n);

    if(S) *S = nullptr;
    if(A) *A = nullptr;
    if(U) *U = nullptr;
    if(V) *V = nullptr;

    char *charBuffer = static_cast<char*>(customWork().data());

    // TODO: If allocate with Real size, cuSOLVER returns 
    //       CUDA_ERROR_MISALIGNED_ADDRESS error (for complex types only)
    // std::size_t sizeS = k * sizeof(T_RScalar);
    std::size_t sizeS = k * sizeof(T_Scalar);
    std::size_t sizeA = m * n * sizeof(T_Scalar);
    std::size_t sizeU = svdVectorSize(jobu, m, k) * sizeof(T_Scalar);
    std::size_t sizeV = svdVectorSize(jobv, n, k) * sizeof(T_Scalar);

    if(S) { *S = reinterpret_cast<T_RScalar*>(charBuffer); } charBuffer += sizeS;
    if(A) { *A = reinterpret_cast<T_Scalar *>(charBuffer); } charBuffer += sizeA;
    if(U) { *U = reinterpret_cast<T_Scalar *>(charBuffer); } charBuffer += sizeU;
    if(V) { *V = reinterpret_cast<T_Scalar *>(charBuffer); } charBuffer += sizeV;
}
/*-------------------------------------------------*/
template void CuSolverHandler::gesvdAssignInternalPointers<real_t>(svdPolicy_t, svdPolicy_t, int_t, int_t, real_t**, real_t**, real_t**, real_t**);
template void CuSolverHandler::gesvdAssignInternalPointers<real4_t>(svdPolicy_t, svdPolicy_t, int_t, int_t, real4_t**, real4_t**, real4_t**, real4_t**);
template void CuSolverHandler::gesvdAssignInternalPointers<complex_t>(svdPolicy_t, svdPolicy_t, int_t, int_t, real_t**, complex_t**, complex_t**, complex_t**);
template void CuSolverHandler::gesvdAssignInternalPointers<complex8_t>(svdPolicy_t, svdPolicy_t, int_t, int_t, real4_t**, complex8_t**, complex8_t**, complex8_t**);
/*-------------------------------------------------*/
} // namespace culite
/*-------------------------------------------------*/
