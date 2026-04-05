# [v0.3.1] - 2026-04-05

## CLA3P Module

### New Features
- **LapackSyev Eigensolver**: Symmetric/Hermitian eigenvalue decomposition on CPU
  - Added `syev`, `heev`, and `geev` LAPACK proxies
  - Implemented `LapackSyev` class for symmetric and Hermitian eigenproblems
  - Added examples and tutorial/module documentation for `eigsol`
- **VirtualStrided Support**: Dense matrix diagonal can now be referenced via `VirtualStrided` class

### Bug Fixes
- Fixed multiplication dimension checks in `mult_dim_check()` for matrix-vector operations
- Fixed update operators for both dense and sparse cases

### Improvements
- **Header Include Restructuring**: All internal include paths reorganized into dedicated subdirectories (`algebra`, `bulk`, `checks`, `dense`, `error`, `generic`, `linsol`, `lra`, `perms`, `proxies`, `qr`, `sparse`, `support`, `svd`, `types`, `virtuals`); all includes now use angle-bracket form
- Renamed `dns?D` file family to `dns?d` for consistent naming conventions
- Replaced `cla3p::svdPolicy_t::Limited` with `Economy` throughout the codebase
- Replaced tabs with spaces in example source files

---

## cuLite Module

### New Features
- **Dense Object Fill**: Added `fill()` member function and `operator=(Scalar)` for dense vectors and matrices
- **Extended cuSOLVER Support**: Additional cuSOLVER operations exposed through the handler interface
- **SVD and GEEV Improvements**:
  - Re-implemented GESVD fix for S-size bug on CUDA
  - Stabilized `DefaultSVD` and `LapackGeev` constructors
  - Added additional members to `LapackGeev`
  - Fixed complex SVD implementation
  - Added examples and tutorial/module documentation for `geev` and `gesvd`
- **New Utility Class**: `GridXD` for GPU kernel launch configuration

### Bug Fixes
- Fixed complex SVD implementation for correctness
- Fixed `complex_t` class arithmetic operators (`+`, `-`)
- Fixed update operators for both dense and sparse cases

### Improvements
- **Implementation Split**: Moved `CuBlasHandler`, `CuSolverHandler`, `CuSparseHandler`, and cuSPARSE proxy template instantiations from header-only to `.cpp` files with explicit specializations, reducing compile times
- **Header Include Restructuring**: All internal include paths reorganized into dedicated subdirectories; all includes now use angle-bracket form
- Replaced tabs with spaces in example source files
- Updated README with SVD support documentation

---

# [v0.3.0] - 2026-03-26

## CLA3P Module

### New Features
- **Enhanced CSR Matrix Support**: Pardiso now accepts CSR matrices in Full and Upper storage formats
- **Template-Only Checks System**: Outer checks now use meta classes and template-only implementations for improved compile-time optimization

### Improvements
- **Code Maintenance**: General code cleanup and refactoring
  - Use type aliases throughout codebase
  - Removed obsolete comments
  - Enhanced virtual expression system maintenance
  - Added extra validation in outer product operations
  - Updated namespace usage in syntax
- **Examples**: Tidied up example programs for better clarity
- **Documentation**: 
  - Fixed documentation regarding the "bind" parameter
  - Updated syntax documentation

### Bug Fixes
- Improved dimension checking using meta classes in `mult_dim_check()`
- Enhanced property validation code

---

## cuLite Module

### New Features
- **Eigenvalue Decomposition (GEEV)**: Complete eigenvalue/eigenvector computation for general matrices
  - Added `geev()` support in CuSolverHandler
  - New eigenvalue solver directory and solver implementations
  - LapackGeev proxy functions for CPU-based eigenvalue decomposition
  - Dedicated cuSolver proxy file for eigenvalue operations
  - Support for both eigenvalues-only and full eigenpair computation
- **Virtual Expression System**: Lazy evaluation for GPU matrices and vectors
  - Virtual row vector class for transposed views
  - Virtual conjugate operations
  - Virtual transpose and conjugate transpose
  - Virtual arithmetic operations (+, -, scaling)
  - Outer product virtual expressions
  - Constructor support for creating objects from virtual expressions
  - Transpose and conjugate-transpose members for sparse matrices
- **Functional Algebra Operations**: High-level functional API for common linear algebra operations
  - Inner product (dot, dotc)
  - Vector addition and scaling (axpy)
  - Matrix-vector multiplication for symmetric/Hermitian matrices
  - Matrix-matrix multiplication for symmetric/Hermitian matrices
  - Outer product operations (ger, gerc, syr, her)
  - Rank-1 updates for symmetric/Hermitian matrices
  - Functional update operations for dense and sparse matrices
  - Diagonal matrix-vector multiplication for CSX matrices
- **Sparse Matrix Operations**:
  - CSR to CSC conversion via cuSPARSE proxy
  - Custom CSX addition kernels for sparse matrix arithmetic
  - Host-device copy operations for CSR matrices
  - Inner and outer algebra operators for sparse matrices
- **cuBLAS Handler Enhancements**: Added multiple new operations
  - `copy()` member for array copying
  - `dotc()` for conjugate dot product
  - `axpy()` for vector scaling and addition
  - `ger()`, `gerc()` for rank-1 updates
  - `syr()`, `her()` for symmetric/Hermitian rank-1 updates

### Improvements
- **Documentation**: 
  - Added comprehensive Doxygen documentation for modules
  - Enhanced documentation for update operations
  - Added documentation for LapackGeev functions
  - Fixed various documentation references
- **Code Organization**:
  - Extracted legacy cuSPARSE proxies to separate file for better maintainability
  - Improved visibility control for virtual expressions
  - Enhanced dns1D code structure
- **Utilities**:
  - Added `memSetZero()` utility function for GPU memory initialization
  - Added complex number arithmetic operators (+, -)

### Bug Fixes
- Fixed compilation error in `operator>>()` for stream operations
- Fixed CuSolver GEEV implementation issues
- Added const casts to fix compilation issues in cuSPARSE proxies
- Resolved virtual expression conflicts by reimplementing Virtual Expression system
- Fixed argument handling in update operations

### API Changes
- Modified `add()` function argument lists for consistency
- Enhanced operator support for update operations (+=, -=, *=)
- Added multiplication operators for matrix-vector and matrix-matrix operations

---

# [v0.2.0] - 2026-02-07

## CLA3P Module

### New Features
- **CSR (Compressed Sparse Row) Matrix Support (BETA)**: Added full CSR matrix class with comprehensive operations matching CSC functionality
  - CSR-specific bulk operations and math functions
  - Sparse proxy support for both MKL and Arm Performance Libraries
  - Introduced CSX base class to unify common CSR/CSC functionality
  - Enhanced sparse matrix conversion capabilities (COO to CSR/CSC)
- **Virtual Expression System**: Enhanced lazy evaluation capabilities
  - New `VirtualStrided` class for strided memory access patterns
  - Improved virtual transpose, scale, and plus/minus operations
  - Better support for complex algebraic expressions with deferred evaluation
- **Copy and Move Semantics**: Refined constructors and assignment operators for better memory management
- **Template-based Checks**: Converted basic validation checks to pure templates for improved compile-time optimization and cuLite integration

### Improvements
- **Documentation**: Comprehensive updates across all modules and tutorials
- **Memory Allocators**: Renamed allocator template functions (`i_?alloc<>()` → `i_?alloc_t<>()`) for consistency
- **Vector Notation**: Updated examples and documentation to use lowercase for vector variables (standard mathematical convention)
- **Code Quality**: Various refactoring and cleanup improvements
  - Extracted base classes from containers (XxContainerBase)
  - Improved const-correctness throughout the codebase
  - Enhanced hermitian coefficient handling

### Bug Fixes
- Fixed COO matrix to CSR conversion
- Fixed compilation issues on macOS
- Corrected various compiler warnings
- Fixed RankModerator documentation
- Corrected examples compilation issues

### API Changes
- Renamed internal allocator functions for better naming conventions
- Enhanced sparse matrix algebra operators to support CSR matrices
- Extended permutation matrix operations to work with new sparse formats

---

## cuLite Module

### New Features (Initial Release, Linux/Windows)
- **CUDA-Accelerated Dense Linear Algebra (BETA)**: Complete GPU-based dense matrix and vector operations
  - Dense vector and matrix classes with GPU memory management
  - Support for real and complex arithmetic (float, double, complex types)
  - Block operations and shallow copy capabilities
- **Dense Linear System Solvers**: GPU-accelerated linear equation solving
  - LU decomposition solver with cuSOLVER backend
  - Automatic memory management between host and device
- **Sparse Matrix Support**: GPU sparse matrix operations
  - CSR (Compressed Sparse Row) matrix class
  - CSC (Compressed Sparse Column) matrix class
  - Sparse matrix-vector and matrix-matrix multiplication via cuSPARSE
- **CUDA Library Integration**: Comprehensive wrappers for NVIDIA libraries
  - cuBLAS handler and proxies for dense linear algebra operations (GEMV, GEMM, TRSV, etc.)
  - cuSOLVER handler and proxies for decompositions and linear solvers
  - cuSPARSE handler and proxies for sparse operations
- **Utility Classes and Infrastructure**:
  - CUDA stream management for asynchronous operations
  - Scoped Pointer Mode (SPM) utility for cuBLAS and cuSPARSE pointer mode handling
  - Device memory allocators with automatic cleanup
  - Custom CUDA kernels for specialized operations
- **Matrix Operations**: Comprehensive set of GPU operations
  - Matrix norms (one-norm, infinity-norm, Frobenius, max-norm)
  - Matrix conjugate and real/imaginary extractors
  - Scaling, addition, and multiplication operations
  - 2D CUDA kernel launcher for custom operations

### Documentation & Examples
- Complete tutorial and API documentation
- Example programs demonstrating:
  - Dense vector and matrix creation
  - Data transfer between host and device
  - Linear system solving
- CMake integration and build scripts
- Comprehensive README with installation instructions

### Development Infrastructure
- Full CMake build system with CUDA support (nvcc compiler)
- Example builder and executor scripts
- Integration with CLA3P design patterns and allocators

### Architecture
- Consistent API design mirroring CLA3P module structure
- Bulk data operations (dns1D, dns2D) for low-level GPU memory operations
- Type traits and enums for type-safe GPU operations
- Exception handling with CUDA-specific error reporting

---

# [v0.1.0] - 2026-01-02

SimuliCore official launch.

