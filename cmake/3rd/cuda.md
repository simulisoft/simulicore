## CUDA & cuDSS Linking for x86_64 Architectures with NVIDIA GPUs

Certain modules within this project depend on **NVIDIA CUDA Toolkit** and **NVIDIA cuDSS** to accelerate computational operations on GPU hardware. For instance, projects utilizing **cuLite** require linkage to both **cuLite**, **NVIDIA CUDA libraries** (cuBLAS, cuSOLVER, cuSPARSE), and **NVIDIA cuDSS**.

Prior to compiling **SimuliCore**, it is necessary to edit the files `cuda.lin.cmake` and/or `cuda.win.cmake` and configure the CMake variables `NVIDIA_CUDA_ROOT` and `NVIDIA_CUDSS_ROOT` to reference the installation directories of **NVIDIA CUDA Toolkit** and **NVIDIA cuDSS**, respectively.

These configuration files implement a standardized setup of custom CMake variables to facilitate the compilation and linking processes. The default configuration employs dynamic linking with CUDA runtime libraries. Users may modify these files to implement alternative configurations, such as static linking with CUDA libraries or specifying particular compute capabilities. Additional information for CUDA compilation and linking is available at the [**NVIDIA CUDA Toolkit Documentation**](https://docs.nvidia.com/cuda/).

Following the installation of **SimuliCore**, copies of the configured files `cuda.lin.cmake` and `cuda.win.cmake` will be located in the `<simulicore_install>/cmake/3rd` directory. These configuration files may be utilized to build custom projects using CMake.

To generate the CUDA and cuDSS CMake variables, include the appropriate file in your CMakeLists.txt as follows:
```cmake
if(${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
  include(<simulicore_install>/cmake/3rd/cuda.win.cmake)
elseif(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
  include(<simulicore_install>/cmake/3rd/cuda.lin.cmake)
endif()

message(STATUS "NVIDIA CUDA include directory: ${NVIDIA_CUDA_INC}")
message(STATUS "NVIDIA CUDA library list: ${NVIDIA_CUDA_LIB}")
message(STATUS "NVIDIA cuDSS include directory: ${NVIDIA_CUDSS_INC}")
message(STATUS "NVIDIA cuDSS library list: ${NVIDIA_CUDSS_LIB}")
```
This configuration establishes the `NVIDIA_CUDA_<...>` and `NVIDIA_CUDSS_<...>` variables for subsequent use.

To add the **NVIDIA CUDA** and **cuDSS** include directories to your CMake target, use:
```cmake
target_include_directories(<target> PRIVATE ${NVIDIA_CUDA_INC} ${NVIDIA_CUDSS_INC})
```

To link your CMake target with **NVIDIA CUDA** and **cuDSS** libraries, use the following methods:
```cmake

target_link_libraries(<target> ${NVIDIA_CUDA_LIB} ${NVIDIA_CUDSS_LIB})

```

The following is a complete enumeration of the CMake variables defined in `cuda.<os>.cmake`:
1) **NVIDIA_CUDA_ROOT** (NVIDIA CUDA Toolkit root directory)
2) **NVIDIA_CUDSS_ROOT** (NVIDIA cuDSS root directory)
3) **NVIDIA_CUDA_INC** (NVIDIA CUDA Toolkit include directory)
4) **NVIDIA_CUDSS_INC** (NVIDIA cuDSS include directory)
5) **NVIDIA_CUDA_LIB_DIR** (NVIDIA CUDA Toolkit library directory)
6) **NVIDIA_CUDA_DLL_DIR** (NVIDIA CUDA Toolkit binary/DLL directory)
7) **NVIDIA_CUDSS_LIB_DIR** (NVIDIA cuDSS library directory)
8) **NVIDIA_CUDSS_DLL_DIR** (NVIDIA cuDSS binary/DLL directory)
9) **NVIDIA_CUDA_LIB** (Complete NVIDIA CUDA library linkage, includes cuBLAS, cuSOLVER, cuSPARSE)
10) **NVIDIA_CUDSS_LIB** (NVIDIA cuDSS library linkage)

