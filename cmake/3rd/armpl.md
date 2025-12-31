## Arm Performance Libraries Linking for ARM64 Architectures

Certain modules within this project depend on **Arm Performance Libraries** to accelerate computational operations. For instance, projects utilizing **CLA3P** require linkage to both **CLA3P** and **Arm Performance Libraries**.

Prior to compiling **SimuliCore**, it is necessary to edit the file `armpl.mac.cmake` and configure the CMake variable `ARMPL_ROOT` to reference the installation directory of **Arm Performance Libraries**.

These configuration files implement a standardized setup of custom CMake variables to facilitate the compilation and linking processes. The default configuration employs dynamic linking with **OpenMP** support. Users may modify these files to implement alternative configurations, such as static linking. Additional information is available in the [**Arm Performance Libraries Reference Guide**](https://developer.arm.com/documentation/101004/2404/General-information/Access-Arm-Performance-Libraries).

Following the installation of **SimuliCore**, a copy of the configured `armpl.mac.cmake` file will be located in the `<simulicore_install>/cmake/3rd` directory. These configuration files may be utilized to build custom projects using CMake.

To generate the ArmPL CMake variables, include the appropriate file in your CMakeLists.txt as follows:
```cmake
if(APPLE)
  include(<simulicore_install>/cmake/3rd/armpl.mac.cmake)
endif()

message(STATUS "ArmPL include directory: ${ARMPL_INC}")
message(STATUS "ArmPL library list (32-bit integers): ${ARMPL_LIB}")
message(STATUS "ArmPL library list (64-bit integers): ${ARMPL_I64_LIB}")
```
This configuration establishes the `ARMPL_<...>` and `ARMPL_I64_<...>` variables for subsequent use.

> <b><i>NOTE:</i></b>  
> **ARMPL_LIB** specifies the library list required for building applications using 32-bit integers.  
> **ARMPL_I64_LIB** specifies the library list required for building applications using 64-bit integers.  
> Only one library list should be selected per target.  
> **ARMPL_DEF** defines the compile definitions for complex number representations and 32-bit integers.  
> **ARMPL_I64_DEF** defines the compile definitions for complex number representations and 64-bit integers (includes INTEGERP64 definition).  

To add the **Arm Performance Libraries** include directory to your CMake target, use:
```cmake
target_include_directories(<target> PRIVATE ${ARMPL_INC})
```

To link your CMake target with **Arm Performance Libraries**, use the following methods:
```cmake
#
# For 32-bit integer lengths
#
target_link_libraries(<target> ${ARMPL_LIB})

#
# For 64-bit integer lengths
#
target_compile_definitions(<target> PRIVATE INTEGER64)
target_link_libraries(<target> ${ARMPL_I64_LIB})
```

The following is a complete enumeration of the CMake variables defined in `armpl.<os>.cmake`:
1) **ARMPL_ROOT** (Arm Performance Libraries root directory)
2) **ARMPL_DEF** (Compile definitions for LP64 and complex number representations)
3) **ARMPL_I64_DEF** (Compile definitions for ILP64 and complex number representations)
4) **ARMPL_INC** (Arm Performance Libraries include directory)
5) **ARMPL_LIB_DIR** (Arm Performance Libraries library directory)
6) **ARMPL_LIB** (Complete Arm Performance Libraries linkage for 32-bit integers)
7) **ARMPL_I64_LIB** (Complete Arm Performance Libraries linkage for 64-bit integers)

<!-- 6) **ARMPL_DLL_DIR** (Arm Performance Libraries redist library directory, windows only) -->

