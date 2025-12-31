<a name="readme-top"></a>




# Overview




## Table of contents

1. [About cuLite](#about-culite)
<!-- 
2. [How to use cuLite](#how-to-use-culite)
3. [Documentation](#documentation)
4. [License](#license)
5. [Contact](#contact)
-->




<a name="about-culite"></a>
## About cuLite

**cuLite** is a modern, lightweight C++ interface built upon CUDA libraries, specifically designed to streamline GPU computing. Currently in the final stages of development, **cuLite** is scheduled for integration and release as an integral component of the [**SimuliCore**](https://www.simulisoft.com/simulicore/index.html) suite in 2026.

<!--
<p align="right">(<a href="#readme-top">back to top</a>)</p>



<a name="how-to-use-culite"></a>
## How to use cuLite

**cuLite** artifacts are located in your **SimuliCore** installation directory. You can use them simply by:
* including directory `<simulicore_install>/include`
* linking your code with the shared libraries located in `<simulicore_install>/lib`
* linking with 3rd party dependencies

```cmake
#
# sample linux CMake setup
#

include(<simulicore_install>/3rd/mpi.lin.cmake)
set(CULITE_INC <simulicore_install>/include)
set(CULITE_LIB -L<simulicore_install>/lib -culite)

add_executable(<target> main.cpp)
target_include_directories(<target> PRIVATE ${CULITE_INC})
target_link_libraries(<target> ${CULITE_LIB} ${MPI_LIB})
```

<p align="right">(<a href="#readme-top">back to top</a>)</p>




<a name="documentation"></a>
## Documentation

You can find the latest **cuLite** version documentation [here](https://www.simulisoft.com/simulicore/v0.1.0/culite/html/index.html).

<p align="right">(<a href="#readme-top">back to top</a>)</p>




<a name="license"></a>
## License

Distributed as a part of [**SimuliCore**](https://www.simulisoft.com/simulicore/index.html) licensed under
the [Apache License, Version 2.0](https://www.apache.org/licenses/LICENSE-2.0).

<p align="right">(<a href="#readme-top">back to top</a>)</p>




<a name="contact"></a>
## Contact

**cuLite** is spanking new, so feel free to contact us about anything regarding this software.  
Please share your opinions, future requests and documentation clarifications, or report an issue [here](https://www.simulisoft.com/simulicore/index.html#contact).

<p align="right">(<a href="#readme-top">back to top</a>)</p>

-->
