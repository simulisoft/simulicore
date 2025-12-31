<a name="readme-top"></a>




# Overview




## Table of contents

1. [About CMPI](#about-cmpi)
<!-- 
2. [How to use CMPI](#how-to-use-cmpi)
3. [Documentation](#documentation)
4. [License](#license)
5. [Contact](#contact)
-->




<a name="about-cmpi"></a>
## About CMPI

**Compact Message Passing Interface** is a modern C++ interface designed to streamline distributed computing and inter-node communication significantly. CMPI is in the final stages of development and is scheduled for release as a core component of the [**SimuliCore**](https://www.simulisoft.com/simulicore/index.html) suite.



<!--
<p align="right">(<a href="#readme-top">back to top</a>)</p>



<a name="how-to-use-cmpi"></a>
## How to use CMPI

**CMPI** artifacts are located in your **SimuliCore** installation directory. You can use them simply by:
* including directory `<simulicore_install>/include`
* linking your code with the shared libraries located in `<simulicore_install>/lib`
* linking with 3rd party dependencies

```cmake
#
# sample linux CMake setup
#

include(<simulicore_install>/3rd/mpi.lin.cmake)
set(CMPI_INC <simulicore_install>/include)
set(CMPI_LIB -L<simulicore_install>/lib -cmpi)

add_executable(<target> main.cpp)
target_include_directories(<target> PRIVATE ${CMPI_INC})
target_link_libraries(<target> ${CMPI_LIB} ${MPI_LIB})
```

<p align="right">(<a href="#readme-top">back to top</a>)</p>




<a name="documentation"></a>
## Documentation

You can find the latest **CMPI** version documentation [here](https://www.simulisoft.com/simulicore/v0.1.0/cmpi/html/index.html).

<p align="right">(<a href="#readme-top">back to top</a>)</p>




<a name="license"></a>
## License

Distributed as a part of [**SimuliCore**](https://www.simulisoft.com/simulicore/index.html) licensed under
the [Apache License, Version 2.0](https://www.apache.org/licenses/LICENSE-2.0).

<p align="right">(<a href="#readme-top">back to top</a>)</p>




<a name="contact"></a>
## Contact

**CMPI** is spanking new, so feel free to contact us about anything regarding this software.  
Please share your opinions, future requests and documentation clarifications, or report an issue [here](https://www.simulisoft.com/simulicore/index.html#contact).

<p align="right">(<a href="#readme-top">back to top</a>)</p>

-->
