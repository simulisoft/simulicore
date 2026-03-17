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

#ifndef CLA3P_ENUMS_HPP_
#define CLA3P_ENUMS_HPP_

/** 
 * @file
 */

#include <ostream>

/*-------------------------------------------------*/
namespace cla3p {
/*-------------------------------------------------*/

/**
 * @ingroup cla3p_module_index_datatypes
 * @enum prop_t
 * @brief Matrix algebraic property enumeration.
 * @details Defines the algebraic properties of matrices used to optimize computational operations.
 *          Each property type imposes specific structural constraints and enables specialized algorithms.
 */
enum class prop_t {
    Undefined   = 0, /**< No algebraic property defined. */
    General        , /**< General matrix with no special algebraic structure. */
    Symmetric      , /**< Symmetric matrix: @f$ A = A^T @f$. */
    Hermitian      , /**< Hermitian matrix: @f$ A = A^H @f$. */
    Triangular     , /**< Triangular or trapezoidal matrix structure. */
    Skew             /**< Skew-symmetric matrix: @f$ A = -A^T @f$. */
};

/**
 * @ingroup cla3p_module_index_datatypes
 * @enum uplo_t
 * @brief Matrix storage pattern enumeration.
 * @details Specifies which triangular portion of a matrix contains the active data.
 *          This enumeration is used to define storage efficiency for symmetric, Hermitian,
 *          and triangular matrices where only one triangular portion needs to be stored.
 */
enum class uplo_t : char {
    Full  = 'F', /**< Full matrix storage; both upper and lower triangular parts are stored and utilized. */
    Upper = 'U', /**< Upper triangular storage; only the upper triangular part is stored and utilized. */
    Lower = 'L'  /**< Lower triangular storage; only the lower triangular part is stored and utilized. */
};

/**
 * @ingroup cla3p_module_index_datatypes
 * @enum op_t
 * @brief Matrix operation enumeration.
 * @details Specifies virtual operations to be applied to matrices in computational routines.
 *          These operations are not explicitly performed; instead, the operation flag is passed
 *          to optimized algorithms that implicitly handle the transformation.
 * @par Example
 *          To compute @f$ y = A^T \cdot x @f$ without explicitly forming @f$ A^T @f$:@n
 @code
 cla3p::ops::mult(1, cla3p::op_t::T, A, x, 0, y);
 @endcode
 */
enum class op_t : char {
    N = 'N', /**< No operation: @f$ \text{op}(A) = A @f$. */
    T = 'T', /**< Transpose operation: @f$ \text{op}(A) = A^T @f$. */
    C = 'C'  /**< Conjugate transpose operation: @f$ \text{op}(A) = A^H @f$. */
};

/**
 * @ingroup cla3p_module_index_datatypes
 * @enum side_t
 * @brief Matrix operand positioning enumeration.
 * @details Specifies whether a matrix operand is positioned on the left or right side
 *          in binary matrix operations such as triangular solves and matrix multiplication.
 */
enum class side_t : char {
    Left  = 'L', /**< Matrix is positioned on the left side of the operation. */
    Right = 'R'  /**< Matrix is positioned on the right side of the operation. */
};

/**
 * @ingroup cla3p_module_index_datatypes
 * @enum diag_t
 * @brief Triangular matrix diagonal type enumeration.
 * @details Specifies whether the diagonal of a triangular matrix consists of unit elements
 *          or arbitrary values. This affects storage requirements and computational algorithms.
 */
enum class diag_t : char {
    NonUnit = 'N', /**< Diagonal elements are arbitrary (non-unit) values. */
    Unit    = 'U'  /**< Diagonal elements are implicitly one; not stored explicitly. */
};

/**
 * @ingroup cla3p_module_index_datatypes
 * @enum dup_t
 * @brief Duplicate entry resolution policy enumeration.
 * @details Defines the policy for handling duplicate entries when constructing sparse matrices
 *          from coordinate format. Multiple entries at the same matrix location are combined
 *          according to the specified policy.
 */
enum class dup_t {
    Sum  = 0, /**< Duplicated entries are summed: @f$ a_{ij} = \sum \text{values}_{ij} @f$. */
    Prod    , /**< Duplicated entries are multiplied: @f$ a_{ij} = \prod \text{values}_{ij} @f$. */
    Amax    , /**< The entry with maximum absolute value is retained: @f$ a_{ij} = \max|\text{values}_{ij}| @f$. */
    Amin      /**< The entry with minimum absolute value is retained: @f$ a_{ij} = \min|\text{values}_{ij}| @f$. */
};

/**
 * @ingroup cla3p_module_index_datatypes
 * @enum decomp_t
 * @brief Matrix decomposition method enumeration.
 * @details Specifies the factorization algorithm to be used for matrix decomposition.
 *          The choice of decomposition method depends on matrix properties and computational requirements.
 */
enum class decomp_t {
    Auto        = 0, /**< Automatic selection based on matrix properties. */
    LLT         = 1, /**< Cholesky decomposition: @f$ A = L \cdot L^T @f$ (for symmetric positive definite matrices). */
    LDLT        = 2, /**< @f$ LDL^T @f$ decomposition: @f$ A = L \cdot D \cdot L^T @f$ (for symmetric indefinite matrices). */
    LU          = 3, /**< LU decomposition with partial pivoting: @f$ P \cdot A = L \cdot U @f$. */
    CompleteLU  = 4, /**< Complete LU decomposition with full pivoting: @f$ P \cdot A \cdot Q = L \cdot U @f$. */
    SymmetricLU = 5  /**< Symmetric LU decomposition with Bunch-Kaufman pivoting. */
};

/**
 * @ingroup cla3p_module_index_datatypes
 * @enum eigRange_t
 * @brief Eigenvalue computation range enumeration.
 * @details Specifies which subset of eigenvalues to compute in eigenvalue decomposition routines.
 *          This allows selective computation to improve performance when only specific eigenvalues
 *          are needed, rather than computing the full spectrum.
 */
enum class eigRange_t : char {
    All   = 'A', /**< All eigenvalues are computed. */
    Value = 'V', /**< Eigenvalues in the half-open interval (vl, vu] are computed. */
    Index = 'I'  /**< Eigenvalues with indices il to iu are computed. */
};

/**
 * @ingroup cla3p_module_index_datatypes
 * @enum svdPolicy_t
 * @brief Singular vector computation policy enumeration.
 * @details Specifies the extent of singular vector computation in Singular Value Decomposition (SVD).
 *          Different policies trade off computational cost against the completeness of the decomposition.
 */
enum class svdPolicy_t : char {
    NoCalculation = 'N', /**< Singular vectors are not computed; only singular values are calculated. */
    Limited       = 'S', /**< @deprecated Use `Economy` instead. */
    Economy       = 'S', /**< Compute the first @f$ \min(m,n) @f$ singular vectors (thin/economy SVD). */
    Full          = 'A'  /**< Compute all @f$ m @f$ left and @f$ n @f$ right singular vectors (full SVD). */
};

/**
 * @ingroup cla3p_module_index_datatypes
 * @enum qrPolicy_t
 * @brief QR decomposition storage policy enumeration.
 * @details Specifies the representation format for the orthogonal matrix @f$ Q @f$ in QR decomposition.
 *          The choice affects memory usage and subsequent computational operations.
 */
enum class qrPolicy_t {
    Reflection = 0, /**< Matrix @f$ Q @f$ is stored implicitly as a product of elementary Householder reflectors. */
    Full            /**< Matrix @f$ Q @f$ is explicitly formed and stored in full matrix format. */
};

/**
 * @ingroup cla3p_module_index_datatypes
 * @enum qrtrans_t
 * @brief Input matrix transposition policy enumeration for QR decomposition.
 * @details Specifies whether the input matrix should be automatically transposed before QR decomposition
 *          to ensure optimal tall matrix orientation (number of rows @f$ \geq @f$ number of columns).
 * @par Behavior
 *          - NoTranspose: Decomposes @f$ A \rightarrow QR @f$ directly.
 *          - AutoTranspose: If @f$ m < n @f$, decomposes @f$ A^H \rightarrow QR @f$ instead.
 */
enum class qrtrans_t {
    NoTranspose   = 0, /**< The input matrix is decomposed without modification. */
    AutoTranspose      /**< The input matrix is conjugate-transposed if @f$ m < n @f$ to form a tall matrix. */
};

/**
 * @ingroup cla3p_module_index_datatypes
 * @enum lraMethod_t
 * @brief Low-rank approximation algorithm enumeration.
 * @details Specifies the computational method for matrix rank reduction in low-rank approximation (LRA) algorithms.
 *          Used as a parameter in the @ref cla3p::lra::RankModerator class to select the factorization strategy.
 */
enum class lraMethod_t {
    StandardSVD = 0, /**< Standard Singular Value Decomposition for rank determination (accurate but computationally expensive). */
    FastQR           /**< Hybrid QR-SVD method optimized for large matrices (faster approximation with acceptable accuracy). */
};

/*-------------------------------------------------*/
} // namespace cla3p
/*-------------------------------------------------*/

/**
 * @ingroup cla3p_module_index_stream_operators
 * @brief Outputs a property type enumeration to a stream.
 * @param[in,out] os The output stream to write to.
 * @param[in] prop The property type enumeration value.
 * @return A reference to the output stream.
 */
std::ostream& operator<<(std::ostream& os, const cla3p::prop_t& prop);

/**
 * @ingroup cla3p_module_index_stream_operators
 * @brief Outputs a storage pattern enumeration to a stream.
 * @param[in,out] os The output stream to write to.
 * @param[in] uplo The storage pattern enumeration value.
 * @return A reference to the output stream.
 */
std::ostream& operator<<(std::ostream& os, const cla3p::uplo_t& uplo);

/**
 * @ingroup cla3p_module_index_stream_operators
 * @brief Outputs a matrix operation enumeration to a stream.
 * @param[in,out] os The output stream to write to.
 * @param[in] op The operation type enumeration value.
 * @return A reference to the output stream.
 */
std::ostream& operator<<(std::ostream& os, const cla3p::op_t& op);

/**
 * @ingroup cla3p_module_index_stream_operators
 * @brief Outputs an operand positioning enumeration to a stream.
 * @param[in,out] os The output stream to write to.
 * @param[in] side The operand side enumeration value.
 * @return A reference to the output stream.
 */
std::ostream& operator<<(std::ostream& os, const cla3p::side_t& side);
/**
 * @ingroup cla3p_module_index_stream_operators
 * @brief Outputs a decomposition method enumeration to a stream.
 * @param[in,out] os The output stream to write to.
 * @param[in] decomp The decomposition method enumeration value.
 * @return A reference to the output stream.
 */
std::ostream& operator<<(std::ostream& os, const cla3p::decomp_t& decomp);

/**
 * @ingroup cla3p_module_index_stream_operators
 * @brief Outputs a low-rank approximation method enumeration to a stream.
 * @param[in,out] os The output stream to write to.
 * @param[in] method The LRA method enumeration value.
 * @return A reference to the output stream.
 */
std::ostream& operator<<(std::ostream& os, const cla3p::lraMethod_t& method);

#endif // CLA3P_ENUMS_HPP_
