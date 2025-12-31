## Intel&reg; MKL Linking for x86_64 Architectures

Certain modules within this project depend on **Intel&reg; Math Kernel Library** to accelerate computational operations. For instance, projects utilizing **CLA3P** require linkage to both **CLA3P** and **Intel&reg; MKL**.

Prior to compiling **SimuliCore**, it is necessary to edit the files `mkl.lin.cmake` and/or `mkl.win.cmake` and configure the CMake variables `INTEL_MKL_ROOT` and `INTEL_ICC_ROOT` to reference the installation directories of **Intel&reg; MKL** and **Intel&reg; Compiler**, respectively.

These configuration files implement a standardized setup of custom CMake variables to facilitate the compilation and linking processes. The default configuration employs dynamic linking with **Intel&reg; OpenMP** as the threading layer. Users may modify these files to implement alternative configurations, such as static linking with **Intel&reg; MKL** libraries or utilizing **GNU OpenMP** instead of **Intel&reg; OpenMP**. Additional information for creating customized configurations is available at the [**Intel&reg; MKL Link Line Advisor**](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-link-line-advisor.html).

Following the installation of **SimuliCore**, copies of the configured files `mkl.lin.cmake` and `mkl.win.cmake` will be located in the `<simulicore_install>/cmake/3rd` directory. These configuration files may be utilized to build custom projects using CMake.

To generate the MKL CMake variables, include the appropriate file in your CMakeLists.txt as follows:
```cmake
if(${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
  include(<simulicore_install>/cmake/3rd/mkl.win.cmake)
elseif(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
  include(<simulicore_install>/cmake/3rd/mkl.lin.cmake)
endif()

message(STATUS "IntelMKL include directory: ${INTEL_MKL_INC}")
message(STATUS "IntelMKL library list (32-bit integers): ${INTEL_MKL_LIB}")
message(STATUS "IntelMKL library list (64-bit integers): ${INTEL_MKL_I64_LIB}")
```
This configuration establishes the `INTEL_MKL_<...>` and `INTEL_MKL_I64_<...>` variables for subsequent use.

> <b><i>NOTE:</i></b>  
> **INTEL_MKL_LIB** specifies the library list required for building applications using 32-bit integers.  
> **INTEL_MKL_I64_LIB** specifies the library list required for building applications using 64-bit integers.  
> Only one library list should be selected per target.  
> **INTEL_MKL_DEF** defines the compile definitions for complex number representations and 32-bit integers.  
> **INTEL_MKL_I64_DEF** defines the compile definitions for complex number representations and 64-bit integers (includes MKL_ILP64 definition).  

To add the **Intel&reg; MKL** include directory to your CMake target, use:
```cmake
target_include_directories(<target> PRIVATE ${INTEL_MKL_INC})
```

To link your CMake target with **Intel&reg; MKL**, use the following methods:
```cmake
#
# For 32-bit integer lengths
#
target_link_libraries(<target> ${INTEL_MKL_LIB})

#
# For 64-bit integer lengths
#
target_compile_definitions(<target> PRIVATE MKL_ILP64)
target_link_libraries(<target> ${INTEL_MKL_I64_LIB})
```

The following is a complete enumeration of the CMake variables defined in `mkl.<os>.cmake`:
1) **INTEL_MKL_ROOT** (Intel&reg; MKL root directory)
2) **INTEL_ICC_ROOT** (Intel&reg; Compiler root directory)
3) **INTEL_MKL_DEF** (Compile definitions for LP64 and complex number representations)
4) **INTEL_MKL_I64_DEF** (Compile definitions for ILP64 and complex number representations)
5) **INTEL_MKL_INC** (Intel&reg; MKL include directory)
6) **INTEL_MKL_LIB_DIR** (Intel&reg; MKL library directory)
7) **INTEL_MKL_DLL_DIR** (Intel&reg; MKL redistributable library directory, Windows only)
8) **INTEL_ICC_LIB_DIR** (Intel&reg; Compiler library directory)
9) **INTEL_ICC_DLL_DIR** (Intel&reg; Compiler redistributable directory, Windows only)
10) **INTEL_MKL_LIB** (Complete Intel&reg; MKL and Intel&reg; Compiler linkage for 32-bit integers)
11) **INTEL_MKL_I64_LIB** (Complete Intel&reg; MKL and Intel&reg; Compiler linkage for 64-bit integers)
