<a name="readme-top"></a>




# Overview




## Table of contents

1. [About CLA3P](#about-cla3p)
2. [How to use CLA3P](#how-to-use-cla3p)
3. [Useful info](#useful-info)
4. [Documentation](#documentation)
5. [Third-Party Dependencies](#third-party-dependencies)
6. [Examples](#examples)
7. [License](#license)
8. [Contact](#contact)




<a name="about-cla3p"></a>
## About CLA3P

**Compact Linear Algebra Parallel Portable Package** is a modern, lightweight numerical linear algebra library distributed as a core component of the [**SimuliCore**](https://www.simulisoft.com/simulicore/index.html) suite.

With its minimalist design, **CLA3P** is ideally suited for novice developers, while its advanced API offers the sophisticated functionality required by experienced users to achieve optimal computational performance.

Utilizing a foundational framework of dense and sparse linear algebra, **CLA3P** is designed for continuous evolution.

### Features Currently Implemented

**CLA3P** is a dynamic and evolving library. Its currently supported capabilities and upcoming features include:

- [x] Dense Matrix algebra
- [x] Sparse Matrix algebra
- [x] Smart algebraic operator sets with symbolic operation layer
- [x] Property specification for matrices
- [x] Column & Block operations
- [x] Vector/Matrix content referencing
- [x] Dense linear system solvers
- [x] Dense QR decompositions
- [x] Dense eigensolvers
- [x] Sparse linear system solvers

All features described above are supported by a comprehensive suite of intuitive **operators** to facilitate rapid development, complemented by an advanced **functional mathematics interface** for optimal performance. This dual-layer approach enables users to select the interface that best aligns with their specific project requirements.

The library continues to expand, with new features integrated into each release. Technical inquiries or feature requests for future versions may be directed to us via the [contact](https://www.simulisoft.com/simulicore/index.html#contact) channels.

<p align="right">(<a href="#readme-top">back to top</a>)</p>




<a name="how-to-use-cla3p"></a>
## How to use CLA3P

**CLA3P** artifacts are located within the **SimuliCore** installation directory. To integrate the library into your project, perform the following steps:
* **Include Path**: Add `<simulicore_install>/include` to your project's include directories.
* **Library Link**: Link your project with the appropriate library found in `<simulicore_install>/lib`.
* **Dependencies**: Ensure your project is linked with the required third-party dependencies.

**CLA3P** provides both 32-bit and 64-bit integer interfaces to accommodate different computational requirements:
* 32-bit Integer Interface: Link your project with the standard library (`libcla3p.so` for Linux/macOS or `cla3p.lib` for Windows).
* 64-bit Integer Interface:
  1. Add the `-DCLA3P_I64` definition to your compilation flags.
  2. Link with the 64-bit specific library (`libcla3p_i64.so` for Linux/macOS or `cla3p_i64.lib` for Windows).

```cmake
#
# Sample Linux x86_64 CMake configuration (32-bit integers)
#

include(<simulicore_install>/cmake/3rd/mkl.lin.cmake)
set(CLA3P_INC <simulicore_install>/include)
set(CLA3P_LIB -L<simulicore_install>/lib -lcla3p)

add_executable(<target> main.cpp)
target_include_directories(<target> PRIVATE ${CLA3P_INC})
target_link_libraries(<target> ${CLA3P_LIB} ${INTEL_MKL_LIB})
```
See section [Third-Party Dependencies](#third-party-dependencies) for more information.

<p align="right">(<a href="#readme-top">back to top</a>)</p>




<a name="useful-info"></a>
## Useful info

**CLA3P** provides a dual-interface architecture designed for both ease of use and computational efficiency. Users may utilize the intuitive, operator-based interface for streamlined development or the advanced, functional-based interface to achieve maximum performance.
For example, a transposed matrix-vector multiplication may be executed as follows:

```cpp
Y = A.transpose() * X;
Z += A.transpose() * X;
```
or by utilizing the mult() function:
```cpp
cla3p::ops::mult(1., cla3p::op_t::T, A, X, 0., Y);
cla3p::ops::mult(1., cla3p::op_t::T, A, X, 1., Z);
```

Leveraging the **CLA3P Virtuals** symbolic operation layer, operations such as `A.transpose()` are not executed explicitly. Instead, a virtually transposed matrix object is generated and passed to the multiplication operator. This approach ensures that the calculation is performed without additional memory allocation or redundant computational overhead.

The **CLA3P Virtuals** framework enables the execution of complex, single-line operations while minimizing the creation of temporary objects. By utilizing template expressions, operations are decomposed into individual sub-operations that are stacked and lazily evaluated only upon assignment. This architectural design ensures that high-level, readable code maintains the performance characteristics of manually optimized, low-level routines.

```cpp
cla3p::dns::RdMatrix A = cla3p::dns::RdMatrix::random(3,3);
cla3p::dns::RdVector X = cla3p::dns::RdVector::random(3);
cla3p::dns::RdVector Y(3);

//
// (1) No additional space required
//
Y = 2. * A * X - A.transpose() * X;

//
// (2) Additional space required for evaluation of vector (- A * X + 2. * X)
//
Y += A.transpose() * ( - A * X + 2. * X );
```
For scenarios requiring optimal performance and explicit control, equivalent operations may be performed via the advanced functional interface. This interface provides direct access to high-performance computational kernels:
```cpp
// Operation (1)
cla3p::ops::mult( 2., cla3p::op_t::N, A, X, 0., Y); // Y = 2 * A * X
cla3p::ops::mult(-1., cla3p::op_t::T, A, X, 1., Y); // Y -= A' * X

// Operation (2)
cla3p::dns::RdVector tmp(3);
cla3p::ops::mult(-1., cla3p::op_t::N, A, X, 0., tmp); // tmp = - A * X
cla3p::ops::update(2., X, tmp); // tmp += 2 * X
cla3p::ops::mult(1., cla3p::op_t::T, A, tmp, 1., Y); // Y += A' * tmp
```

**CLA3P Virtuals** are currently available for all vector and matrix objects.

> <b><i>NOTE:</i></b>  
> In cases where a virtual matrix product is formed by a sparse matrix and a transposed dense matrix (e.g., `A * B.transpose()`), the library explicitly calculates the transpose of the dense matrix. If this specific operation is required within a computational loop, it is recommended to pre-calculate the transpose of matrix `B` outside the loop to avoid redundant explicit computations and minimize computational overhead.

<p align="right">(<a href="#readme-top">back to top</a>)</p>




<a name="documentation"></a>
## Documentation

You can find the latest **CLA3P** version documentation [here](https://www.simulisoft.com/simulicore/v0.1.0/cla3p/html/index.html).

<p align="right">(<a href="#readme-top">back to top</a>)</p>




<a name="third-party-dependencies"></a>
## Third-Party Dependencies

To accelerate computations, **CLA3P** relies on optimized third-party performance libraries tailored to specific hardware architectures:
* **Linux/Windows (x86_64)**: Requires [**Intel&reg; Math Kernel Library**](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html). Detailed instructions for linking with **Intel&reg; MKL** are available in the provided linking [guide](https://github.com/simulisoft/simulicore/blob/main/cmake/3rd/mkl.md).
* **macOS (arm64)**: Requires [**Arm Performance Libraries**](https://developer.arm.com/Tools%20and%20Software/Arm%20Performance%20Libraries). A dedicated guide for linking with **Arm Performance Libraries** is available [here](https://github.com/simulisoft/simulicore/blob/main/cmake/3rd/armpl.md).

<p align="right">(<a href="#readme-top">back to top</a>)</p>




<a name="examples"></a>
## Examples

The directory `<simulicore_install>/examples/cla3p` contains comprehensive templates and source examples demonstrating the configuration and compilation of custom projects utilizing the **CLA3P** library.

### Building Examples on Linux/macOS
Execute the following commands from a terminal:
```sh
cd <simulicore_install>/examples/cla3p
./example_builder.sh
```
Example executables are organized within the `i32/bin` and `i64/bin` directories, corresponding to their respective 32-bit and 64-bit integer interfaces.
To execute a demonstration, select the desired script using the format `ex<number>_<description>.sh` and run it from your terminal:
```sh
./i32/bin/ex01a_dense_vector_create.sh
```

### Building Examples on Windows

To build and execute the provided examples on Windows using Visual Studio, follow the procedure below:
* **Project Initialization**: Open the `<simulicore_install>/examples/cla3p` directory in Visual Studio to initialize the CMake project.
* **Configuration**: By default, examples are configured for the 32-bit integer interface. To utilize the 64-bit interface, set the CMake variable `-DCLA3P_EXAMPLES_I64=ON`.
* **Compilation & Installation**: Build the solution to generate the example executables and perform installation.
* **Execution**: All binaries are located in the `ixx/bin` directory (where `ixx` corresponds to `i32` or `i64`, depending on your selection). To run an example, execute the desired batch file via the Visual Studio terminal:

```sh
./i32/bin/ex01a_dense_vector_create.bat
```

<p align="right">(<a href="#readme-top">back to top</a>)</p>




<a name="license"></a>
## License

Distributed as a part of [**SimuliCore**](https://www.simulisoft.com/simulicore/index.html) licensed under the [Apache License, Version 2.0](https://www.apache.org/licenses/LICENSE-2.0).

<p align="right">(<a href="#readme-top">back to top</a>)</p>




<a name="contact"></a>
## Contact

As a recently released framework, **CLA3P** is under active development and refinement. We welcome user feedback regarding the software's functionality, documentation, or implementation.

Please submit your insights, feature requests, documentation inquiries, or technical issue reports through the following channels:
* **Issue Tracking**: Submit a report via the SimuliCore GitHub Issues [page](https://github.com/simulisoft/simulicore/issues).
* **Direct Inquiry**: Contact the development team at info@simulisoft.com.

<p align="right">(<a href="#readme-top">back to top</a>)</p>


