<a name="readme-top"></a>




# Overview




## Table of contents

1. [About cuLite](#about-culite)
2. [How to use cuLite](#how-to-use-culite)
3. [Useful info](#useful-info)
4. [Documentation](#documentation)
5. [Third-Party Dependencies](#third-party-dependencies)
6. [Examples](#examples)
7. [License](#license)
8. [Contact](#contact)




<a name="about-culite"></a>
## About cuLite

**CUDA Lightweight Interface** is a modern, lightweight C++ wrapper library built upon NVIDIA CUDA libraries and **CLA3P**, distributed as a core component of the [**SimuliCore**](https://www.simulisoft.com/simulicore/index.html) suite.

With its minimalist design, **cuLite** provides an intuitive interface to GPU computing primitives, making it accessible for novice developers while offering advanced functionality required by experienced users to achieve optimal GPU performance. Built as a GPU acceleration layer on top of **CLA3P** (Compact Linear Algebra Parallel Portable Package), **cuLite** maintains full compatibility with **CLA3P** data structures and operations within the **SimuliCore** ecosystem, enabling seamless hybrid CPU-GPU computing workflows.

Utilizing a foundational framework of GPU memory management and accelerated computation, **cuLite** is designed for continuous evolution.

### Features Currently Implemented

**cuLite** is a dynamic and evolving library. Its currently supported capabilities and upcoming features include:

- [x] GPU memory management (device/host allocation and transfers)
- [x] Dense vector/matrix management
- [x] Sparse matrix management
- [x] Dense vector/matrix algebra
- [x] Sparse matrix algebra (limited)
- [x] Virtual operation layer
- [x] Dense linear system solvers (LU)
- [ ] Sparse linear system solvers (cuDSS)
- [ ] QR decomposition
- [x] Singular Value Decomposition (SVD)
- [x] Eigenvalue problems (limited support)
- [x] CUDA stream management (limited support)
- [x] cuBLAS integration (Basic Linear Algebra Subprograms on GPU - limited support)
- [x] cuSOLVER integration (Dense linear system solvers on GPU - limited support)
- [x] cuSPARSE integration (Sparse matrix operations on GPU - limited support)

All features described above are supported by a comprehensive suite of intuitive **handler classes** to facilitate rapid development, complemented by an advanced **proxy layer** for direct CUDA library access. This dual-layer approach enables users to select the interface that best aligns with their specific project requirements.

The library continues to expand, with new features integrated into each release. Technical inquiries or feature requests for future versions may be directed to us via the [contact](https://www.simulisoft.com/simulicore/index.html#contact) channels.

<p align="right">(<a href="#readme-top">back to top</a>)</p>




<a name="how-to-use-culite"></a>
## How to use cuLite

**cuLite** artifacts are located within the **SimuliCore** installation directory. To integrate the library into your project, perform the following steps:
* **Include Path**: Add `<simulicore_install>/include` to your project's include directories.
* **Library Link**: Link your project with the appropriate library found in `<simulicore_install>/lib`.
* **Dependencies**: Ensure your project is linked with the required CUDA libraries.

**cuLite** provides both 32-bit and 64-bit integer interfaces to accommodate different computational requirements:
* 32-bit Integer Interface: Link your project with the standard library (`libculite.so` for Linux or `culite.lib` for Windows).
* 64-bit Integer Interface:
  1. Add the `-DCULITE_I64` definition to your compilation flags.
  2. Link with the 64-bit specific library (`libculite_i64.so` for Linux or `culite_i64.lib` for Windows).

```cmake
#
# Sample Linux x86_64 CMake configuration (32-bit integers)
# cuLite needs CLA3P library and NVIDIA CUDA libraries
# CLA3P needs Intel MKL library
#

include(<simulicore_install>/cmake/3rd/mkl.lin.cmake)
include(<simulicore_install>/cmake/3rd/cuda.lin.cmake)
set(SIMULICORE_INC <simulicore_install>/include)
set(SIMULICORE_LIB -L<simulicore_install>/lib -lculite -lcla3p)
set(SIMULICORE_3RD_DEF ${INTEL_MKL_DEF})
set(SIMULICORE_3RD_INC ${INTEL_MKL_INC} ${NVIDIA_CUDA_INC} ${NVIDIA_CUDSS_INC})
set(SIMULICORE_3RD_LIB ${INTEL_MKL_LIB} ${NVIDIA_CUDA_LIB} ${NVIDIA_CUDSS_LIB})

add_executable(<target> main.cpp)
target_compile_definitions(<target> PRIVATE ${SIMULICORE_3RD_DEF})
target_include_directories(<target> PRIVATE  ${SIMULICORE_INC} ${SIMULICORE_3RD_INC})
target_link_libraries(<target> ${SIMULICORE_LIB} ${SIMULICORE_3RD_LIB})
```
See section [Third-Party Dependencies](#third-party-dependencies) for more information.

<p align="right">(<a href="#readme-top">back to top</a>)</p>




<a name="useful-info"></a>
## Useful info

**cuLite** provides a streamlined interface to NVIDIA CUDA libraries, abstracting away much of the boilerplate code typically required for GPU computing while maintaining high performance and flexibility.

### Memory Management

The library simplifies GPU memory operations through intuitive handler classes:

```cpp
#include <cla3p/dense.hpp>
#include <culite/dense.hpp>

// Transfer data from host to device
cla3p::dns::RdVector<cla3p::real_t> hostVec(1000);

// Allocate device memory
culite::dns::RdVector<culite::real_t> devVec1(1000);
culite::dns::RdVector<culite::real_t> devVec2;

// Transfer data from host to device
hostVec >> devVec1;
hostVec >> devVec2; // allocates memory for devVec2 internally

// Transfer data from device to host
devVec1 >> hostVec;
```

### GPU Linear Algebra Operations

**cuLite** exposes matrix-vector multiplication at three levels of abstraction, from low-level BLAS calls up to expressive C++ operators. The examples below demonstrate all three approaches for computing `y = alpha * A * x + beta * y` using dense GPU objects.

**1. Direct cuBLAS handler interface** — raw access to the underlying BLAS routine via `culite::CuBlasHandler::gemv`:

```cpp
culite::dns::RdMatrix A(m, n);
culite::dns::RdVector x(n);
culite::dns::RdVector y(m);

// Fill A and x with data...

// Compute y = alpha * A * x + beta * y through the global cuBLAS handler
culite::real_t alpha = 1.0, beta = 0.0;
culite::CuBlasHandler& cuBlasHandler = culite::globalCuBlasHandler();
cuBlasHandler.gemv(cla3p::op_t::N,
                   A.nrows(), A.ncols(),
                   &alpha, A.values(), A.ld(),
                   x.values(), 1,
                   &beta, y.values(), 1);
```

**2. Functional interface with a custom cuBLAS handler** — higher-level `culite::ops::mult` function accepting an explicit handler, useful when multiple independent BLAS contexts are required:

```cpp
culite::dns::RdMatrix A(m, n);
culite::dns::RdVector x(n);
culite::dns::RdVector y(m);

// Fill A and x with data...

// Create a dedicated cuBLAS handler and compute y = alpha * A * x + beta * y
culite::CuBlasHandler myHandler;
culite::real_t alpha = 1.0, beta = 0.0;
culite::ops::mult(alpha, cla3p::op_t::N, A, x, beta, y, myHandler);
```

**3. Operator interface** — the most concise form; the global cuBLAS handler is used internally:

```cpp
culite::dns::RdMatrix A(m, n);
culite::dns::RdVector x(n);

// Fill A and x with data...

// Compute y = A * x (global cuBLAS handler is used)
culite::dns::RdVector y = A * x;

// Or accumulate in-place: y += A * x
y += A * x;
```

### Stream Management

For advanced users requiring explicit control over GPU execution, **cuLite** provides stream management capabilities. The CUDA stream class is currently under development and offers limited functionality:

```cpp
#include <culite/generic/cuda_stream.hpp>

// Create and manage CUDA streams via cudaStream class
culite::CudaStream cudaStream;

// Get the underlying CUDA stream handle
cudaStream_t stream = cudaStream.stream();

// Synchronize the stream
cudaStream.sync();
```

The library is designed to work seamlessly with **CLA3P** for hybrid CPU-GPU computing workflows, enabling efficient data transfer and computation across different processing units.

<p align="right">(<a href="#readme-top">back to top</a>)</p>




<a name="documentation"></a>
## Documentation

You can find the latest **cuLite** version documentation [here](https://www.simulisoft.com/simulicore/v0.3.0/culite/html/index.html).

<p align="right">(<a href="#readme-top">back to top</a>)</p>




<a name="third-party-dependencies"></a>
## Third-Party Dependencies

To accelerate computations on GPU hardware, **cuLite** relies on NVIDIA CUDA libraries:
* **Linux/Windows (x86_64 with NVIDIA GPU)**:
  - [**NVIDIA CUDA Toolkit**](https://developer.nvidia.com/cuda-toolkit): The CUDA Toolkit includes cuBLAS, cuSOLVER, and cuSPARSE libraries required by **cuLite**.
  - [**cuDSS (CUDA Direct Sparse Solver)**](https://developer.nvidia.com/cudss): Required for advanced sparse linear system solving capabilities, provides high-performance direct solvers for sparse matrices on GPU.

Detailed instructions for linking with CUDA libraries are available in the provided linking [guide](https://github.com/simulisoft/simulicore/blob/main/cmake/3rd/cuda.md).

**Minimum Requirements**:
* CUDA Toolkit 13.0 or later
* NVIDIA GPU with compute capability 5.0 or higher
* Compatible NVIDIA GPU driver
* cuDSS library (for sparse solver functionality)

<p align="right">(<a href="#readme-top">back to top</a>)</p>




<a name="examples"></a>
## Examples

The directory `<simulicore_install>/examples/culite` contains comprehensive templates and source examples demonstrating the configuration and compilation of custom projects utilizing the **cuLite** library.

### Building Examples on Linux
Execute the following commands from a terminal:
```sh
cd <simulicore_install>/examples/culite
./example_builder.sh
```
Example executables are organized within the `i32/bin` and `i64/bin` directories, corresponding to their respective 32-bit and 64-bit integer interfaces.
To execute a demonstration, select the desired script using the format `ex<number>_<description>.sh` and run it from your terminal:
```sh
./i32/bin/ex01a_dense_vector_create.sh
```

### Building Examples on Windows

To build and execute the provided examples on Windows using Visual Studio, follow the procedure below:
* **Project Initialization**: Open the `<simulicore_install>/examples/culite` directory in Visual Studio to initialize the CMake project.
* **Configuration**: By default, examples are configured for the 32-bit integer interface. To utilize the 64-bit interface, set the CMake variable `-DCULITE_EXAMPLES_I64=ON`.
* **Compilation & Installation**: Build the solution to generate the example executables and perform installation.
* **Execution**: All binaries are located in the `ixx/bin` directory (where `ixx` corresponds to `i32` or `i64`, depending on your selection). To run an example, execute the desired batch file via the Visual Studio terminal:

```sh
./i32/bin/ex01a_dense_vector_create.bat
```
> <b><i>NOTE (Compliance Requirement): Build Configuration Alignment</i></b>  
> To ensure stable execution, verify that the build configuration (Debug or Release) of the application matches the configuration used to compile the library. Mixing configurations is strictly discouraged due to fundamental incompatibilities in the C Runtime Library (CRT).


<p align="right">(<a href="#readme-top">back to top</a>)</p>




<a name="license"></a>
## License

Distributed as a part of [**SimuliCore**](https://www.simulisoft.com/simulicore/index.html) licensed under the [Apache License, Version 2.0](https://www.apache.org/licenses/LICENSE-2.0).

<p align="right">(<a href="#readme-top">back to top</a>)</p>




<a name="contact"></a>
## Contact

As a recently released framework, **cuLite** is under active development and refinement. We welcome user feedback regarding the software's functionality, documentation, or implementation.

Please submit your insights, feature requests, documentation inquiries, or technical issue reports through the following channels:
* **Issue Tracking**: Submit a report via the SimuliCore GitHub Issues [page](https://github.com/simulisoft/simulicore/issues).
* **Direct Inquiry**: Contact the development team at info@simulisoft.com.

<p align="right">(<a href="#readme-top">back to top</a>)</p>
