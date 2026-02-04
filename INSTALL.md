# SimuliCore Installation Guide

Before proceeding with the installation, please review the [Prerequisites](README.md#prerequisites) section to ensure your system meets all necessary requirements.

## Table of Contents

1. [Linux Installation](#linux-installation)
2. [Windows Installation](#windows-installation)
3. [macOS Installation](#macos-installation)
4. [Additional Configuration Arguments](#extra-configuration-arguments)




<a name="linux-installation"></a>
## Linux Installation

1. Open a terminal window.
  
2. Clone the repository:
   
   ```sh
   git clone git@github.com:simulisoft/simulicore.git
   cd simulicore
   ```
   
   Throughout this document, the top-level directory will be referred to as `<simulicore_root>`.

3. To configure Intel&reg; MKL, open the following file:
   
   ```sh
   <simulicore_root>/cmake/3rd/mkl.lin.cmake
   ```
   
   Configure the paths `INTEL_MKL_ROOT` and `INTEL_ICC_ROOT` to reference your installed copies of **Intel&reg; MKL** and **Intel&reg; Compiler**, respectively.
   
   For additional information, please consult the [Intel&reg; MKL Linking Guide](cmake/3rd/mkl.md).

4. To configure NVIDIA&reg; CUDA, open the following file:
   
   ```sh
   <simulicore_root>/cmake/3rd/cuda.lin.cmake
   ```
   
   Configure the paths `NVIDIA_CUDA_ROOT` and `NVIDIA_CUDSS_ROOT` to reference your installed copies of **NVIDIA&reg; CUDA Toolkit** and **NVIDIA&reg; cuDSS**, respectively.
   
   For additional information, please consult the [CUDA & cuDSS Linking Guide](cmake/3rd/cuda.md).

   <a name="linux-cmake-config"></a>
5. Select a build directory `build_dir` and execute the CMake command to configure the project:
   
   ```sh
   cmake -S <simulicore_root> -B <build_dir>
   ```
   
   For additional configuration options, please refer to the [Additional Configuration Arguments](#extra-configuration-arguments) section below.

6. Compile and install the project:
    
   ```sh
   cmake --build <build_dir>
   cmake --install <build_dir>
   ```




<a name="windows-installation"></a>
## Windows Installation

1. Clone the repository using your preferred cloning method:
   
   ```sh
   https: https://github.com/simulisoft/simulicore.git
   ssh: git@github.com:simulisoft/simulicore.git
   ```
   
   Throughout this document, the top-level directory will be referred to as `<simulicore_root>`.

2. To configure Intel&reg; MKL, open the following file:
   
   ```sh
   <simulicore_root>/cmake/3rd/mkl.win.cmake
   ```
   
   Configure the paths `INTEL_MKL_ROOT` and `INTEL_ICC_ROOT` to reference your installed copies of **Intel&reg; MKL** and **Intel&reg; Compiler**, respectively.

   For additional information, please consult the [Intel&reg; MKL Linking Guide](cmake/3rd/mkl.md).

3. To configure NVIDIA&reg; CUDA, open the following file:
   
   ```sh
   <simulicore_root>/cmake/3rd/cuda.win.cmake
   ```
   
   Configure the paths `NVIDIA_CUDA_ROOT` and `NVIDIA_CUDSS_ROOT` to reference your installed copies of **NVIDIA&reg; CUDA Toolkit** and **NVIDIA&reg; cuDSS**, respectively.

   For additional information, please consult the [CUDA & cuDSS Linking Guide](cmake/3rd/cuda.md).

   > **Warning:** Microsoft Visual Studio may fail to configure the project with CUDA 13. If you encounter configuration issues, please use the alternative BASH-based workflow below (steps 4-7) instead of the Visual Studio approach.

   **Option A (Visual Studio):** Launch Microsoft Visual Studio and open the `simulicore_root` directory.

   <a name="windows-cmake-config"></a>
   4a. Configure CMake within Microsoft Visual Studio by navigating to **Project > CMake Settings for SimuliCore**.
   
   For additional configuration options, please refer to the [Additional Configuration Arguments](#extra-configuration-arguments) section below.
   
   5a. Compile the project by selecting **Build > Build All**.
   
   6a. Install the project by selecting **Build > Install SimuliCore**.

   **Option B (Windows BASH - Alternative for CUDA 13 compatibility):**

   4b. Open a Windows BASH terminal (e.g., Git Bash).

   5b. Select a build directory `build_dir` and execute the CMake command to configure the project:
   
   ```sh
   cmake -S <simulicore_root> -B <build_dir>
   ```
   
   For additional configuration options, please refer to the [Additional Configuration Arguments](#extra-configuration-arguments) section below.

   6b. Compile and install the project:
    
   ```sh
   cmake --build <build_dir> --config Release
   cmake --install <build_dir> --config Release
   ```

7. The installation will be complete once the build and install commands finish successfully.




<a name="macos-installation"></a>
## macOS Installation

1. Open a terminal window.
  
2. Clone the repository:
   
   ```sh
   git clone git@github.com:simulisoft/simulicore.git
   cd simulicore
   ```
   
   Throughout this document, the top-level directory will be referred to as `<simulicore_root>`.

3. Open the following file:
   
   ```sh
   <simulicore_root>/cmake/3rd/armpl.mac.cmake
   ```
   
   Configure the path `ARMPL_ROOT` to reference your installed copy of **Arm Performance Libraries**.
   
   For additional information, please consult the [Arm Performance Libraries Linking Guide](cmake/3rd/armpl.md).

   <a name="macos-cmake-config"></a>
4. Select a build directory `build_dir` and execute the CMake command to configure the project:
   
   ```sh
   cmake -S <simulicore_root> -B <build_dir>
   ```
   
   For additional configuration options, please refer to the [Additional Configuration Arguments](#extra-configuration-arguments) section below.

5. Compile and install the project:
    
   ```sh
   cmake --build <build_dir>
   cmake --install <build_dir>
   ```



<a name="extra-configuration-arguments"></a>
## Additional Configuration Arguments

The following list presents configuration arguments applicable to the entire project.

- **Linux**: Append the desired options to the [CMake Configuration](#linux-cmake-config) command.
- **macOS**: Append the desired options to the [CMake Configuration](#macos-cmake-config) command.
- **Windows**: Add the desired options in the "CMake command arguments" section within the [CMake Settings](#windows-cmake-config).

#### <ins>Installation Directory</ins>

The default installation directories are as follows:
- Linux/macOS: `<simulicore_root>/install`
- Windows: `<simulicore_root>/out/install/<build-type>`

To specify a custom installation directory, use one of the following options:

```cmake
--install-prefix=<your_install_path>
```
or
```cmake
-DCMAKE_INSTALL_PREFIX=<your_install_path>
```



