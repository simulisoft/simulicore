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

