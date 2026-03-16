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

#ifndef CCULITE_CUBLAS_HANDLER_HPP_
#define CCULITE_CUBLAS_HANDLER_HPP_

/**
 * @file
 */

#include <culite/types/integer.hpp>
#include <culite/types/enums.hpp>
#include <culite/proxies/cublas_proxy.hpp>

/*-------------------------------------------------*/
namespace culite { 
/*-------------------------------------------------*/

class CuBlasHandler;

/*-------------------------------------------------*/

/**
 * @nosubgrouping
 * @brief RAII helper for managing cuBLAS pointer mode.
 * @details This class provides a scoped guard that temporarily sets a cuBLAS pointer mode
 *          and automatically restores the previous mode when the object goes out of scope.
 *          This ensures proper cleanup even in the presence of exceptions.
 */
class CuBlasSpm {
    public:

        // disable copying
        CuBlasSpm(const CuBlasSpm&) = delete;
        CuBlasSpm& operator=(const CuBlasSpm&) = delete;

        /**
         * @brief Constructor.
         * @details Sets the pointer mode for the given CuBlasHandler and stores the previous mode.
         * @param[in] cuBlasHandler The CuBlasHandler whose pointer mode will be modified.
         * @param[in] mode The new pointer mode to set.
         */
        CuBlasSpm(CuBlasHandler& cuBlasHandler, cublasPointerMode_t mode);
        
        /**
         * @brief Destructor.
         * @details Restores the pointer mode that was active before the constructor was called.
         */
        ~CuBlasSpm();

    private:
        CuBlasHandler& m_handler;
        cublasPointerMode_t m_oldMode;
};  

/*-------------------------------------------------*/

/**
 * @nosubgrouping
 * @brief The cuBlas handler class.
 * @details This class provides a wrapper around the cuBLAS library for performing
 *          Basic Linear Algebra Subprograms (BLAS) operations on GPU devices.
 *          It manages the cuBLAS handle and provides convenient interfaces for
 *          vector and matrix operations.
 */
class CuBlasHandler {

    public:

        // disable copying
        CuBlasHandler(const CuBlasHandler&) = delete;
        CuBlasHandler& operator=(const CuBlasHandler&) = delete;
        
        /**
         * @brief Constructor.
         * @details Initializes the cuBLAS handle and internal state.
         */
        CuBlasHandler();
        
        /**
         * @brief Destructor.
         * @details Destroys the cuBLAS handle and releases all allocated resources.
         */
        ~CuBlasHandler();

        /**
         * @brief Get the cuBlas handle.
         * @return The cuBlas handle.
         */
        cublasHandle_t handle();

        /**
         * @brief Set the pointer mode for cuBLAS operations.
         * @details Changes the pointer mode and returns the previous mode. The pointer mode
         *          determines whether scalar values (alpha, beta) are passed by reference
         *          on the host or device.
         * @param[in] mode The new pointer mode to set.
         * @return The pointer mode that was active before this call.
         */
        cublasPointerMode_t setPointerMode(cublasPointerMode_t mode);
        
        /**
         * @brief Get the current pointer mode.
         * @details Returns the current pointer mode setting for cuBLAS operations.
         * @return The current pointer mode.
         */
        cublasPointerMode_t pointerMode();

        /**
         * @brief Finds the index of the maximum absolute value element.
         * @details Computes the index of the element with the maximum absolute value in vector @p x.
         * @tparam T_Scalar The scalar type of the vector elements.
         * @param[in] n The number of elements in the vector.
         * @param[in] x Pointer to the device vector.
         * @param[in] incx The stride between consecutive elements of @p x.
         * @param[out] result Pointer to store the result index (1-based).
         */
        template <typename T_Scalar>
        void iamax(int_t n, const T_Scalar *x, int_t incx, int_t *result);

        /**
         * @brief Computes the sum of absolute values.
         * @details Computes @f$ \sum_{i=1}^{n} |x_i| @f$ for the elements of vector @p x.
         * @tparam T_Scalar The scalar type of the vector elements.
         * @param[in] n The number of elements in the vector.
         * @param[in] x Pointer to the device vector.
         * @param[in] incx The stride between consecutive elements of @p x.
         * @param[out] result Pointer to store the sum of absolute values.
         */
        template <typename T_Scalar>
        void asum(int_t n, const T_Scalar *x, int_t incx, typename TypeTraits<T_Scalar>::real_type *result);

        /**
         * @brief Computes the Euclidean norm (L2 norm).
         * @details Computes @f$ \sqrt{\sum_{i=1}^{n} |x_i|^2} @f$ for the elements of vector @p x.
         * @tparam T_Scalar The scalar type of the vector elements.
         * @param[in] n The number of elements in the vector.
         * @param[in] x Pointer to the device vector.
         * @param[in] incx The stride between consecutive elements of @p x.
         * @param[out] result Pointer to store the Euclidean norm.
         */
        template <typename T_Scalar>
        void nrm2(int_t n, const T_Scalar *x, int_t incx, typename TypeTraits<T_Scalar>::real_type *result);

        /**
         * @brief Copies a vector to another vector.
         * @details Copies vector @p x to vector @p y, computing @f$ y = x @f$.
         * @tparam T_Scalar The scalar type of the vector elements.
         * @param[in] n The number of elements in the vectors.
         * @param[in] x Pointer to the source device vector.
         * @param[in] incx The stride between consecutive elements of @p x.
         * @param[out] y Pointer to the destination device vector.
         * @param[in] incy The stride between consecutive elements of @p y.
         */
        template <typename T_Scalar>
        void copy(int_t n, const T_Scalar *x, int_t incx, T_Scalar *y, int_t incy);

        /**
         * @brief Scales a vector by a scalar.
         * @details Computes @f$ x = \alpha \cdot x @f$ for vector @p x.
         * @tparam T_Scalar The scalar type of the vector elements.
         * @param[in] n The number of elements in the vector.
         * @param[in] alpha Pointer to the scalar multiplier.
         * @param[in,out] x Pointer to the device vector (modified in-place).
         * @param[in] incx The stride between consecutive elements of @p x.
         */
        template <typename T_Scalar>
        void scal(int_t n, const T_Scalar *alpha, T_Scalar *x, int_t incx);

        /**
         * @brief Computes the dot product of two vectors.
         * @details Computes @f$ result = x^T \cdot y = \sum_{i=1}^{n} x_i \cdot y_i @f$.
         * @tparam T_Scalar The scalar type of the vector elements.
         * @param[in] n The number of elements in the vectors.
         * @param[in] x Pointer to the first vector on device.
         * @param[in] incx The stride between consecutive elements of @p x.
         * @param[in] y Pointer to the second vector on device.
         * @param[in] incy The stride between consecutive elements of @p y.
         * @param[out] result Pointer to store the dot product result.
         */
        template <typename T_Scalar>
        void dot(int_t n, const T_Scalar *x, int_t incx, const T_Scalar *y, int_t incy, T_Scalar *result);

        /**
         * @brief Computes the conjugate dot product of two vectors.
         * @details Computes @f$ result = x^H \cdot y = \sum_{i=1}^{n} \overline{x_i} \cdot y_i @f$,
         *          where @f$ \overline{x_i} @f$ denotes the complex conjugate of @f$ x_i @f$.
         * @tparam T_Scalar The scalar type of the vector elements.
         * @param[in] n The number of elements in the vectors.
         * @param[in] x Pointer to the first vector on device (will be conjugated).
         * @param[in] incx The stride between consecutive elements of @p x.
         * @param[in] y Pointer to the second vector on device.
         * @param[in] incy The stride between consecutive elements of @p y.
         * @param[out] result Pointer to store the conjugate dot product result.
         */
        template <typename T_Scalar>
        void dotc(int_t n, const T_Scalar *x, int_t incx, const T_Scalar *y, int_t incy, T_Scalar *result);

        /**
         * @brief Computes a vector plus scalar times a vector.
         * @details Computes @f$ y = \alpha \cdot x + y @f$, updating vector @p y in-place.
         * @tparam T_Scalar The scalar type of the vector elements.
         * @param[in] n The number of elements in the vectors.
         * @param[in] alpha Pointer to the scalar multiplier.
         * @param[in] x Pointer to the input vector on device.
         * @param[in] incx The stride between consecutive elements of @p x.
         * @param[in,out] y Pointer to the vector on device (modified in-place).
         * @param[in] incy The stride between consecutive elements of @p y.
         */
        template <typename T_Scalar>
        void axpy(int_t n, const T_Scalar *alpha, const T_Scalar *x, int_t incx, T_Scalar *y, int_t incy);

        /**
         * @brief Performs parametrized matrix addition.
         * @details Computes @f$ C = \alpha \cdot op_A(A) + \beta \cdot op_B(B) @f$, where
         *          @f$ op_A @f$ and @f$ op_B @f$ can independently be no-transpose, transpose,
         *          or conjugate transpose operations.
         * @tparam T_Scalar The scalar type of the matrix elements.
         * @param[in] opA The operation to apply to matrix @p a (no-transpose, transpose, or conjugate transpose).
         * @param[in] opB The operation to apply to matrix @p b (no-transpose, transpose, or conjugate transpose).
         * @param[in] m The number of rows in matrices @f$ op_A(a) @f$, @f$ op_B(b) @f$, and @p c.
         * @param[in] n The number of columns in matrices @f$ op_A(a) @f$, @f$ op_B(b) @f$, and @p c.
         * @param[in] alpha Pointer to the scalar multiplier for matrix @p a.
         * @param[in] a Pointer to matrix @p a on device.
         * @param[in] lda The leading dimension of matrix @p a.
         * @param[in] beta Pointer to the scalar multiplier for matrix @p b.
         * @param[in] b Pointer to matrix @p b on device.
         * @param[in] ldb The leading dimension of matrix @p b.
         * @param[out] c Pointer to matrix @p c on device (result).
         * @param[in] ldc The leading dimension of matrix @p c.
         */
        template <typename T_Scalar>
        void geam(op_t opA,
                  op_t opB,
                  int_t m, int_t n,
                  const T_Scalar *alpha, const T_Scalar *a, int_t lda,
                  const T_Scalar *beta,  const T_Scalar *b, int_t ldb,
                  T_Scalar *c, int_t ldc);

        /**
         * @brief Performs matrix-diagonal matrix multiplication.
         * @details Computes @f$ C = A \cdot diag(x) @f$ if @p side is @ref cla3p::side_t::Right,
         *          or @f$ C = diag(x) \cdot A @f$ if @p side is @ref cla3p::side_t::Left,
         *          where @f$ diag(x) @f$ is a diagonal matrix constructed from vector @p x.
         * @tparam T_Scalar The scalar type of the matrix and vector elements.
         * @param[in] side Specifies whether the diagonal matrix multiplies from the left or right.
         * @param[in] m The number of rows in matrices @p a and @p c.
         * @param[in] n The number of columns in matrices @p a and @p c.
         * @param[in] a Pointer to matrix @p a on device.
         * @param[in] lda The leading dimension of matrix @p a.
         * @param[in] x Pointer to the diagonal vector on device.
         * @param[in] incx The stride between consecutive elements of @p x.
         * @param[out] c Pointer to matrix @p c on device (result).
         * @param[in] ldc The leading dimension of matrix @p c.
         */
        template <typename T_Scalar>
        void dgmm(side_t side,
                  int_t m, int_t n,
                  const T_Scalar *a, int_t lda,
                  const T_Scalar *x, int_t incx,
                  T_Scalar *c, int_t ldc);

        /**
         * @brief Performs general rank-1 update (unconjugated).
         * @details Computes @f$ A = \alpha \cdot x \cdot y^T + A @f$.
         * @tparam T_Scalar The scalar type of the vector and matrix elements.
         * @param[in] m The number of rows in matrix @p a.
         * @param[in] n The number of columns in matrix @p a.
         * @param[in] alpha Pointer to the scalar multiplier.
         * @param[in] x Pointer to the first vector on device.
         * @param[in] incx The stride between consecutive elements of @p x.
         * @param[in] y Pointer to the second vector on device.
         * @param[in] incy The stride between consecutive elements of @p y.
         * @param[in,out] a Pointer to matrix @p a on device (modified in-place).
         * @param[in] lda The leading dimension of matrix @p a.
         */
        template <typename T_Scalar>
        void ger(int_t m, int_t n,
                 const T_Scalar *alpha,
                 const T_Scalar *x, int_t incx,
                 const T_Scalar *y, int_t incy,
                 T_Scalar *a, int_t lda);

        /**
         * @brief Performs general rank-1 update (conjugated).
         * @details Computes @f$ A = \alpha \cdot x \cdot y^H + A @f$, where @f$ y^H @f$ denotes
         *          the conjugate transpose of vector @p y.
         * @tparam T_Scalar The scalar type of the vector and matrix elements.
         * @param[in] m The number of rows in matrix @p a.
         * @param[in] n The number of columns in matrix @p a.
         * @param[in] alpha Pointer to the scalar multiplier.
         * @param[in] x Pointer to the first vector on device.
         * @param[in] incx The stride between consecutive elements of @p x.
         * @param[in] y Pointer to the second vector on device (will be conjugated).
         * @param[in] incy The stride between consecutive elements of @p y.
         * @param[in,out] a Pointer to matrix @p a on device (modified in-place).
         * @param[in] lda The leading dimension of matrix @p a.
         */
        template <typename T_Scalar>
        void gerc(int_t m, int_t n,
                  const T_Scalar *alpha,
                  const T_Scalar *x, int_t incx,
                  const T_Scalar *y, int_t incy,
                  T_Scalar *a, int_t lda);

        /**
         * @brief Performs symmetric rank-1 update.
         * @details Computes @f$ A = \alpha \cdot x \cdot x^T + A @f$, where @f$ A @f$ is symmetric.
         * @tparam T_Scalar The scalar type of the vector and matrix elements.
         * @param[in] uplo Specifies whether the upper or lower triangular part of @p a is referenced.
         * @param[in] n The number of rows and columns in matrix @p a.
         * @param[in] alpha Pointer to the scalar multiplier.
         * @param[in] x Pointer to the vector on device.
         * @param[in] incx The stride between consecutive elements of @p x.
         * @param[in,out] a Pointer to the symmetric matrix @p a on device (modified in-place).
         * @param[in] lda The leading dimension of matrix @p a.
         */
        template <typename T_Scalar>
        void syr(uplo_t uplo,
                 int_t n,
                 const T_Scalar *alpha,
                 const T_Scalar *x, int_t incx,
                 T_Scalar *a, int_t lda);

        /**
         * @brief Performs Hermitian rank-1 update.
         * @details Computes @f$ A = \alpha \cdot x \cdot x^H + A @f$, where @f$ A @f$ is Hermitian.
         * @tparam T_Scalar The scalar type of the vector and matrix elements.
         * @param[in] uplo Specifies whether the upper or lower triangular part of @p a is referenced.
         * @param[in] n The number of rows and columns in matrix @p a.
         * @param[in] alpha Pointer to the real scalar multiplier.
         * @param[in] x Pointer to the vector on device.
         * @param[in] incx The stride between consecutive elements of @p x.
         * @param[in,out] a Pointer to the Hermitian matrix @p a on device (modified in-place).
         * @param[in] lda The leading dimension of matrix @p a.
         */
        template <typename T_Scalar>
        void her(uplo_t uplo,
                 int_t n,
                 const typename TypeTraits<T_Scalar>::real_type *alpha,
                 const T_Scalar *x, int_t incx,
                 T_Scalar *a, int_t lda);

        /**
         * @brief Performs general matrix-vector multiplication.
         * @details Computes @f$ y = \alpha \cdot op(A) \cdot x + \beta \cdot y @f$.
         * @tparam T_Scalar The scalar type of the matrix and vector elements.
         * @param[in] op The operation to apply to matrix @p a (no-transpose, transpose, or conjugate transpose).
         * @param[in] m The number of rows in matrix @f$ op(A) @f$.
         * @param[in] n The number of columns in matrix @f$ op(A) @f$.
         * @param[in] alpha Pointer to the scalar multiplier.
         * @param[in] a Pointer to matrix @p a on device.
         * @param[in] lda The leading dimension of matrix @p a.
         * @param[in] x Pointer to the input vector on device.
         * @param[in] incx The stride between consecutive elements of @p x.
         * @param[in] beta Pointer to the scalar multiplier for @p y.
         * @param[in,out] y Pointer to the result vector on device (modified in-place).
         * @param[in] incy The stride between consecutive elements of @p y.
         */
        template <typename T_Scalar>
        void gemv(op_t op,
                  int_t m, int_t n, const T_Scalar* alpha,
                  const T_Scalar *a, int_t lda,
                  const T_Scalar *x, int_t incx,
                  const T_Scalar* beta, T_Scalar *y, int_t incy);

        /**
         * @brief Performs symmetric matrix-vector multiplication.
         * @details Computes @f$ y = \alpha \cdot A \cdot x + \beta \cdot y @f$, where @f$ A @f$ is symmetric.
         * @tparam T_Scalar The scalar type of the matrix and vector elements.
         * @param[in] uplo Specifies whether the upper or lower triangular part of @p a is referenced.
         * @param[in] n The number of rows and columns in matrix @p a.
         * @param[in] alpha Pointer to the scalar multiplier.
         * @param[in] a Pointer to the symmetric matrix @p a on device.
         * @param[in] lda The leading dimension of matrix @p a.
         * @param[in] x Pointer to the input vector on device.
         * @param[in] incx The stride between consecutive elements of @p x.
         * @param[in] beta Pointer to the scalar multiplier for @p y.
         * @param[in,out] y Pointer to the result vector on device (modified in-place).
         * @param[in] incy The stride between consecutive elements of @p y.
         */
        template <typename T_Scalar>
        void symv(uplo_t uplo,
                  int_t n, const T_Scalar* alpha,
                  const T_Scalar *a, int_t lda,
                  const T_Scalar *x, int_t incx,
                  const T_Scalar* beta, T_Scalar *y, int_t incy);

        /**
         * @brief Performs Hermitian matrix-vector multiplication.
         * @details Computes @f$ y = \alpha \cdot A \cdot x + \beta \cdot y @f$, where @f$ A @f$ is Hermitian.
         * @tparam T_Scalar The scalar type of the matrix and vector elements.
         * @param[in] uplo Specifies whether the upper or lower triangular part of @p a is referenced.
         * @param[in] n The number of rows and columns in matrix @p a.
         * @param[in] alpha Pointer to the scalar multiplier.
         * @param[in] a Pointer to the Hermitian matrix @p a on device.
         * @param[in] lda The leading dimension of matrix @p a.
         * @param[in] x Pointer to the input vector on device.
         * @param[in] incx The stride between consecutive elements of @p x.
         * @param[in] beta Pointer to the scalar multiplier for @p y.
         * @param[in,out] y Pointer to the result vector on device (modified in-place).
         * @param[in] incy The stride between consecutive elements of @p y.
         */
        template <typename T_Scalar>
        void hemv(uplo_t uplo,
                  int_t n, const T_Scalar* alpha,
                  const T_Scalar *a, int_t lda,
                  const T_Scalar *x, int_t incx,
                  const T_Scalar* beta, T_Scalar *y, int_t incy);

        /**
         * @brief Performs triangular matrix-vector multiplication.
         * @details Computes @f$ x = op(A) \cdot x @f$, where @f$ A @f$ is triangular.
         * @tparam T_Scalar The scalar type of the matrix and vector elements.
         * @param[in] uplo Specifies whether the upper or lower triangular part of @p a is used.
         * @param[in] trans The operation to apply to matrix @p a (no-transpose, transpose, or conjugate transpose).
         * @param[in] diag Specifies whether the matrix is unit triangular or not.
         * @param[in] n The number of rows and columns in matrix @p a.
         * @param[in] a Pointer to the triangular matrix @p a on device.
         * @param[in] lda The leading dimension of matrix @p a.
         * @param[in,out] x Pointer to the vector on device (modified in-place).
         * @param[in] incx The stride between consecutive elements of @p x.
         */
        template <typename T_Scalar>
        void trmv(uplo_t uplo,
                  op_t trans,
                  diag_t diag,
                  int_t n, const T_Scalar* a, int_t lda,
                  T_Scalar *x, int_t incx);

        /**
         * @brief Solves a triangular system of equations.
         * @details Solves @f$ op(A) \cdot x = b @f$, where @f$ A @f$ is triangular and @p x contains @f$ b @f$ on entry.
         * @tparam T_Scalar The scalar type of the matrix and vector elements.
         * @param[in] uplo Specifies whether the upper or lower triangular part of @p a is used.
         * @param[in] trans The operation to apply to matrix @p a (no-transpose, transpose, or conjugate transpose).
         * @param[in] diag Specifies whether the matrix is unit triangular or not.
         * @param[in] n The number of rows and columns in matrix @p a.
         * @param[in] a Pointer to the triangular matrix @p a on device.
         * @param[in] lda The leading dimension of matrix @p a.
         * @param[in,out] x Pointer to the vector on device containing @f$ b @f$ on entry and the solution on exit.
         * @param[in] incx The stride between consecutive elements of @p x.
         */
        template <typename T_Scalar>
        void trsv(uplo_t uplo,
                  op_t trans,
                  diag_t diag,
                  int_t n, const T_Scalar* a, int_t lda,
                  T_Scalar *x, int_t incx);

        /**
         * @brief Performs general matrix-matrix multiplication.
         * @details Computes @f$ C = \alpha \cdot op(A) \cdot op(B) + \beta \cdot C @f$.
         * @tparam T_Scalar The scalar type of the matrix elements.
         * @param[in] transa The operation to apply to matrix @p a (no-transpose, transpose, or conjugate transpose).
         * @param[in] transb The operation to apply to matrix @p b (no-transpose, transpose, or conjugate transpose).
         * @param[in] m The number of rows in matrices @f$ op(A) @f$ and @p c.
         * @param[in] n The number of columns in matrices @f$ op(B) @f$ and @p c.
         * @param[in] k The number of columns in @f$ op(A) @f$ and rows in @f$ op(B) @f$.
         * @param[in] alpha Pointer to the scalar multiplier.
         * @param[in] a Pointer to matrix @p a on device.
         * @param[in] lda The leading dimension of matrix @p a.
         * @param[in] b Pointer to matrix @p b on device.
         * @param[in] ldb The leading dimension of matrix @p b.
         * @param[in] beta Pointer to the scalar multiplier for @p c.
         * @param[in,out] c Pointer to matrix @p c on device (modified in-place).
         * @param[in] ldc The leading dimension of matrix @p c.
         */
        template <typename T_Scalar>
        void gemm(op_t transa,
                  op_t transb,
                  int_t m, int_t n, int_t k,
                  const T_Scalar* alpha,
                  const T_Scalar* a, int_t lda,
                  const T_Scalar* b, int_t ldb,
                  const T_Scalar* beta,
                  T_Scalar* c, int_t ldc);

        /**
         * @brief Performs symmetric matrix-matrix multiplication.
         * @details Computes @f$ C = \alpha \cdot A \cdot B + \beta \cdot C @f$ if @p side is Left,
         *          or @f$ C = \alpha \cdot B \cdot A + \beta \cdot C @f$ if @p side is Right, where @f$ A @f$ is symmetric.
         * @tparam T_Scalar The scalar type of the matrix elements.
         * @param[in] side Specifies whether the symmetric matrix @p a multiplies from the left or right.
         * @param[in] uplo Specifies whether the upper or lower triangular part of @p a is referenced.
         * @param[in] m The number of rows in matrices @p b and @p c.
         * @param[in] n The number of columns in matrices @p b and @p c.
         * @param[in] alpha Pointer to the scalar multiplier.
         * @param[in] a Pointer to the symmetric matrix @p a on device.
         * @param[in] lda The leading dimension of matrix @p a.
         * @param[in] b Pointer to matrix @p b on device.
         * @param[in] ldb The leading dimension of matrix @p b.
         * @param[in] beta Pointer to the scalar multiplier for @p c.
         * @param[in,out] c Pointer to matrix @p c on device (modified in-place).
         * @param[in] ldc The leading dimension of matrix @p c.
         */
        template <typename T_Scalar>
        void symm(side_t side,
                  uplo_t uplo,
                  int_t m, int_t n,
                  const T_Scalar* alpha,
                  const T_Scalar* a, int_t lda,
                  const T_Scalar* b, int_t ldb,
                  const T_Scalar* beta,
                  T_Scalar* c, int_t ldc);

        /**
         * @brief Performs Hermitian matrix-matrix multiplication.
         * @details Computes @f$ C = \alpha \cdot A \cdot B + \beta \cdot C @f$ if @p side is Left,
         *          or @f$ C = \alpha \cdot B \cdot A + \beta \cdot C @f$ if @p side is Right, where @f$ A @f$ is Hermitian.
         * @tparam T_Scalar The scalar type of the matrix elements.
         * @param[in] side Specifies whether the Hermitian matrix @p a multiplies from the left or right.
         * @param[in] uplo Specifies whether the upper or lower triangular part of @p a is referenced.
         * @param[in] m The number of rows in matrices @p b and @p c.
         * @param[in] n The number of columns in matrices @p b and @p c.
         * @param[in] alpha Pointer to the scalar multiplier.
         * @param[in] a Pointer to the Hermitian matrix @p a on device.
         * @param[in] lda The leading dimension of matrix @p a.
         * @param[in] b Pointer to matrix @p b on device.
         * @param[in] ldb The leading dimension of matrix @p b.
         * @param[in] beta Pointer to the scalar multiplier for @p c.
         * @param[in,out] c Pointer to matrix @p c on device (modified in-place).
         * @param[in] ldc The leading dimension of matrix @p c.
         */
        template <typename T_Scalar>
        void hemm(side_t side,
                  uplo_t uplo,
                  int_t m, int_t n,
                  const T_Scalar* alpha,
                  const T_Scalar* a, int_t lda,
                  const T_Scalar* b, int_t ldb,
                  const T_Scalar* beta,
                  T_Scalar* c, int_t ldc);

        /**
         * @brief Performs triangular matrix-matrix multiplication.
         * @details Computes @f$ C = \alpha \cdot op(A) \cdot B @f$ if @p side is Left,
         *          or @f$ C = \alpha \cdot B \cdot op(A) @f$ if @p side is Right, where @f$ A @f$ is triangular.
         * @tparam T_Scalar The scalar type of the matrix elements.
         * @param[in] side Specifies whether the triangular matrix @p a multiplies from the left or right.
         * @param[in] uplo Specifies whether the upper or lower triangular part of @p a is used.
         * @param[in] trans The operation to apply to matrix @p a (no-transpose, transpose, or conjugate transpose).
         * @param[in] diag Specifies whether the matrix is unit triangular or not.
         * @param[in] m The number of rows in matrices @p b and @p c.
         * @param[in] n The number of columns in matrices @p b and @p c.
         * @param[in] alpha Pointer to the scalar multiplier.
         * @param[in] a Pointer to the triangular matrix @p a on device.
         * @param[in] lda The leading dimension of matrix @p a.
         * @param[in] b Pointer to matrix @p b on device.
         * @param[in] ldb The leading dimension of matrix @p b.
         * @param[out] c Pointer to matrix @p c on device (result).
         * @param[in] ldc The leading dimension of matrix @p c.
         */
        template <typename T_Scalar>
        void trmm(side_t side,
                  uplo_t uplo,
                  op_t trans,
                  diag_t diag,
                  int_t m, int_t n,
                  const T_Scalar* alpha,
                  const T_Scalar* a, int_t lda,
                  const T_Scalar* b, int_t ldb,
                  T_Scalar* c, int_t ldc);

        /**
         * @brief Solves a triangular matrix equation.
         * @details Solves @f$ op(A) \cdot X = \alpha \cdot B @f$ if @p side is Left,
         *          or @f$ X \cdot op(A) = \alpha \cdot B @f$ if @p side is Right, where @f$ A @f$ is triangular.
         *          The solution @f$ X @f$ overwrites @p b on exit.
         * @tparam T_Scalar The scalar type of the matrix elements.
         * @param[in] side Specifies whether the triangular matrix @p a is on the left or right.
         * @param[in] uplo Specifies whether the upper or lower triangular part of @p a is used.
         * @param[in] trans The operation to apply to matrix @p a (no-transpose, transpose, or conjugate transpose).
         * @param[in] diag Specifies whether the matrix is unit triangular or not.
         * @param[in] m The number of rows in matrix @p b.
         * @param[in] n The number of columns in matrix @p b.
         * @param[in] alpha Pointer to the scalar multiplier.
         * @param[in] a Pointer to the triangular matrix @p a on device.
         * @param[in] lda The leading dimension of matrix @p a.
         * @param[in,out] b Pointer to matrix @p b on device (contains right-hand side on entry, solution on exit).
         * @param[in] ldb The leading dimension of matrix @p b.
         */
        template <typename T_Scalar>
        void trsm(side_t side,
                  uplo_t uplo,
                  op_t trans,
                  diag_t diag,
                  int_t m, int_t n,
                  const T_Scalar* alpha,
                  const T_Scalar* a, int_t lda,
                  T_Scalar* b, int_t ldb);

        /**
         * @brief Performs symmetric rank-k update.
         * @details Computes @f$ C = \alpha \cdot A \cdot A^T + \beta \cdot C @f$ if @p trans is NoTrans,
         *          or @f$ C = \alpha \cdot A^T \cdot A + \beta \cdot C @f$ if @p trans is Trans, where @f$ C @f$ is symmetric.
         * @tparam T_Scalar The scalar type of the matrix elements.
         * @param[in] uplo Specifies whether the upper or lower triangular part of @p c is updated.
         * @param[in] trans The operation to apply to matrix @p a (no-transpose or transpose).
         * @param[in] n The number of rows and columns in matrix @p c.
         * @param[in] k The inner dimension of the multiplication.
         * @param[in] alpha Pointer to the scalar multiplier.
         * @param[in] a Pointer to matrix @p a on device.
         * @param[in] lda The leading dimension of matrix @p a.
         * @param[in] beta Pointer to the scalar multiplier for @p c.
         * @param[in,out] c Pointer to the symmetric matrix @p c on device (modified in-place).
         * @param[in] ldc The leading dimension of matrix @p c.
         */
        template <typename T_Scalar>
        void syrk(uplo_t uplo,
                    op_t trans,
                    int_t n, int_t k,
                    const T_Scalar* alpha,
                    const T_Scalar* a, int_t lda,
                    const T_Scalar* beta,
                    T_Scalar* c, int_t ldc);

        /**
         * @brief Performs Hermitian rank-k update.
         * @details Computes @f$ C = \alpha \cdot A \cdot A^H + \beta \cdot C @f$ if @p trans is NoTrans,
         *          or @f$ C = \alpha \cdot A^H \cdot A + \beta \cdot C @f$ if @p trans is ConjTrans, where @f$ C @f$ is Hermitian.
         * @tparam T_Scalar The scalar type of the matrix elements.
         * @param[in] uplo Specifies whether the upper or lower triangular part of @p c is updated.
         * @param[in] trans The operation to apply to matrix @p a (no-transpose or conjugate transpose).
         * @param[in] n The number of rows and columns in matrix @p c.
         * @param[in] k The inner dimension of the multiplication.
         * @param[in] alpha Pointer to the real scalar multiplier.
         * @param[in] a Pointer to matrix @p a on device.
         * @param[in] lda The leading dimension of matrix @p a.
         * @param[in] beta Pointer to the real scalar multiplier for @p c.
         * @param[in,out] c Pointer to the Hermitian matrix @p c on device (modified in-place).
         * @param[in] ldc The leading dimension of matrix @p c.
         */
        template <typename T_Scalar>
        void herk(uplo_t uplo,
                  op_t trans,
                  int_t n, int_t k,
                  const typename TypeTraits<T_Scalar>::real_type* alpha,
                  const T_Scalar* a, int_t lda,
                  const typename TypeTraits<T_Scalar>::real_type* beta,
                  T_Scalar* c, int_t ldc);

        /**
         * @brief Performs symmetric rank-k update with two matrices.
         * @details Computes @f$ C = \alpha \cdot A \cdot B^T + \beta \cdot C @f$ if @p trans is NoTrans,
         *          or @f$ C = \alpha \cdot A^T \cdot B + \beta \cdot C @f$ if @p trans is Trans, where @f$ C @f$ is symmetric.
         * @tparam T_Scalar The scalar type of the matrix elements.
         * @param[in] uplo Specifies whether the upper or lower triangular part of @p c is updated.
         * @param[in] trans The operation to apply to matrices @p a and @p b (no-transpose or transpose).
         * @param[in] n The number of rows and columns in matrix @p c.
         * @param[in] k The inner dimension of the multiplication.
         * @param[in] alpha Pointer to the scalar multiplier.
         * @param[in] a Pointer to matrix @p a on device.
         * @param[in] lda The leading dimension of matrix @p a.
         * @param[in] b Pointer to matrix @p b on device.
         * @param[in] ldb The leading dimension of matrix @p b.
         * @param[in] beta Pointer to the scalar multiplier for @p c.
         * @param[in,out] c Pointer to the symmetric matrix @p c on device (modified in-place).
         * @param[in] ldc The leading dimension of matrix @p c.
         */
        template <typename T_Scalar>
        void syrkx(uplo_t uplo,
                    op_t trans,
                    int_t n, int_t k,
                    const T_Scalar* alpha,
                    const T_Scalar* a, int_t lda,
                    const T_Scalar* b, int_t ldb,
                    const T_Scalar* beta,
                    T_Scalar* c, int_t ldc);

        /**
         * @brief Performs Hermitian rank-k update with two matrices.
         * @details Computes @f$ C = \alpha \cdot A \cdot B^H + \beta \cdot C @f$ if @p trans is NoTrans,
         *          or @f$ C = \alpha \cdot A^H \cdot B + \beta \cdot C @f$ if @p trans is ConjTrans, where @f$ C @f$ is Hermitian.
         * @tparam T_Scalar The scalar type of the matrix elements.
         * @param[in] uplo Specifies whether the upper or lower triangular part of @p c is updated.
         * @param[in] trans The operation to apply to matrices @p a and @p b (no-transpose or conjugate transpose).
         * @param[in] n The number of rows and columns in matrix @p c.
         * @param[in] k The inner dimension of the multiplication.
         * @param[in] alpha Pointer to the scalar multiplier.
         * @param[in] a Pointer to matrix @p a on device.
         * @param[in] lda The leading dimension of matrix @p a.
         * @param[in] b Pointer to matrix @p b on device.
         * @param[in] ldb The leading dimension of matrix @p b.
         * @param[in] beta Pointer to the real scalar multiplier for @p c.
         * @param[in,out] c Pointer to the Hermitian matrix @p c on device (modified in-place).
         * @param[in] ldc The leading dimension of matrix @p c.
         */
        template <typename T_Scalar>
        void herkx(uplo_t uplo,
                   op_t trans,
                   int_t n, int_t k,
                   const T_Scalar* alpha,
                   const T_Scalar* a, int_t lda,
                   const T_Scalar* b, int_t ldb,
                   const typename TypeTraits<T_Scalar>::real_type* beta,
                   T_Scalar* c, int_t ldc);
                
    private:
        cublasHandle_t m_handle{nullptr};
};

/*-------------------------------------------------*/

/**
 * @ingroup culite_module_index_cublas_utils
 * @brief Returns the global cuBLAS handler instance.
 * @details Provides access to a singleton cuBLAS handler that can be used throughout
 *          the application for cuBLAS operations.
 * @return Reference to the global @ref CuBlasHandler instance.
 */
inline CuBlasHandler& globalCuBlasHandler()
{
    static CuBlasHandler handler;
    return handler;
}

/*-------------------------------------------------*/
} // namespace culite
/*-------------------------------------------------*/

#endif // CULITE_CUBLAS_HANDLER_HPP_
