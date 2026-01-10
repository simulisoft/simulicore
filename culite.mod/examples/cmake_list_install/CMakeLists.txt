cmake_minimum_required(VERSION 3.12.0)

project("CULITE Examples" VERSION "0.1.0" LANGUAGES CXX)

#-----------------------------------------------
# set global definitions
#-----------------------------------------------
option(CULITE_EXAMPLES_I64 "cuLite integer precision (32/64bits)" OFF)
set(SIMULICORE_INSTALL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../..")

#-----------------------------------------------
# build setup
#-----------------------------------------------
include(${SIMULICORE_INSTALL_DIR}/cmake/os.cmake)
include(${SIMULICORE_INSTALL_DIR}/cmake/compiler.cmake)
include(${SIMULICORE_INSTALL_DIR}/cmake/build_type.cmake)

message(STATUS "Detected OS: ${CMAKE_SYSTEM_NAME} (${CMAKE_SYSTEM_PROCESSOR})")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

#-----------------------------------------------
# installation setup
#-----------------------------------------------
if(${CULITE_EXAMPLES_I64})
	set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/i64" CACHE PATH "Installation directory" FORCE)
else()
	set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/i32" CACHE PATH "Installation directory" FORCE)
endif()
message(STATUS "cuLite Examples Install Directory: ${CMAKE_INSTALL_PREFIX}")

#-----------------------------------------------
# set openmp
#-----------------------------------------------
#include(${SIMULICORE_INSTALL_DIR}/cmake/openmp.cmake)

#-----------------------------------------------
# set 3rd parties (initialize with cuLite and CLA3P)
#-----------------------------------------------

set(CULITE_EXAMPLES_3RD_PARTY_DEF)
set(CULITE_EXAMPLES_3RD_PARTY_I64_DEF CULITE_I64 CLA3P_I64)
set(CULITE_EXAMPLES_3RD_PARTY_INC ${SIMULICORE_INSTALL_DIR}/include)

set(CULITE_EXAMPLES_3RD_PARTY_LIB_PATH ${CULITE_EXAMPLES_3RD_PARTY_LIB_PATH} ${SIMULICORE_INSTALL_DIR}/lib)
set(CULITE_EXAMPLES_3RD_PARTY_DLL_PATH ${CULITE_EXAMPLES_3RD_PARTY_DLL_PATH} ${SIMULICORE_INSTALL_DIR}/lib)
if(${SIMULICORE_SYSTEM_LINUX})
	set(CULITE_EXAMPLES_3RD_PARTY_LIB -L${CULITE_EXAMPLES_3RD_PARTY_LIB_PATH} -lculite -lcla3p)
	set(CULITE_EXAMPLES_3RD_PARTY_I64_LIB -L${CULITE_EXAMPLES_3RD_PARTY_LIB_PATH} -lculite_i64 -lcla3p_i64)
elseif(${SIMULICORE_SYSTEM_WIN32})
	set(CULITE_EXAMPLES_3RD_PARTY_LIB ${CULITE_EXAMPLES_3RD_PARTY_LIB_PATH}/culite.lib ${CULITE_EXAMPLES_3RD_PARTY_LIB_PATH}/cla3p.lib)
	set(CULITE_EXAMPLES_3RD_PARTY_I64_LIB ${CULITE_EXAMPLES_3RD_PARTY_LIB_PATH}/culite_i64.lib ${CULITE_EXAMPLES_3RD_PARTY_LIB_PATH}/cla3p_i64.lib)
elseif(${SIMULICORE_SYSTEM_MACOS})
	message(FATAL_ERROR "cuLite is not supported on macOS")
endif()

# 
# armpl/accelerate setup
#
#if((${SIMULICORE_SYSTEM_LINUX} OR ${SIMULICORE_SYSTEM_WIN32}) AND ${SIMULICORE_ARCH_ARM64})
#
#	if(${SIMULICORE_SYSTEM_LINUX})
#		include(${SIMULICORE_INSTALL_DIR}/cmake/3rd/armpl.lin.cmake)
#	elseif(${SIMULICORE_SYSTEM_WIN32})
#		include(${SIMULICORE_INSTALL_DIR}/cmake/3rd/armpl.win.cmake)
#	endif()
#	
#	# definitions
#	set(CULITE_EXAMPLES_3RD_PARTY_DEF ${CULITE_EXAMPLES_3RD_PARTY_DEF} ${ARMPL_DEF})
#	set(CULITE_EXAMPLES_3RD_PARTY_I64_DEF ${CULITE_EXAMPLES_3RD_PARTY_I64_DEF} ${ARMPL_I64_DEF})
#	# includes
#	set(CULITE_EXAMPLES_3RD_PARTY_INC ${CULITE_EXAMPLES_3RD_PARTY_INC} ${ARMPL_INC})
#	# libraries
#	set(CULITE_EXAMPLES_3RD_PARTY_LIB ${CULITE_EXAMPLES_3RD_PARTY_LIB} ${ARMPL_LIB})
#	set(CULITE_EXAMPLES_3RD_PARTY_I64_LIB ${CULITE_EXAMPLES_3RD_PARTY_I64_LIB} ${ARMPL_I64_LIB})
#	# library paths
#	set(CULITE_EXAMPLES_3RD_PARTY_LIB_PATH ${CULITE_EXAMPLES_3RD_PARTY_LIB_PATH} ${ARMPL_LIB_DIR})
#	set(CULITE_EXAMPLES_3RD_PARTY_DLL_PATH ${CULITE_EXAMPLES_3RD_PARTY_DLL_PATH} ${ARMPL_DLL_DIR})

#endif()

# 
# mkl setup
#
if((${SIMULICORE_SYSTEM_LINUX} OR ${SIMULICORE_SYSTEM_WIN32}) AND ${SIMULICORE_ARCH_X86_64})

	if(${SIMULICORE_SYSTEM_LINUX})
		include(${SIMULICORE_INSTALL_DIR}/cmake/3rd/mkl.lin.cmake)
	elseif(${SIMULICORE_SYSTEM_WIN32})
		include(${SIMULICORE_INSTALL_DIR}/cmake/3rd/mkl.win.cmake)
	endif()

	# definitions
	set(CULITE_EXAMPLES_3RD_PARTY_DEF ${CULITE_EXAMPLES_3RD_PARTY_DEF} ${INTEL_MKL_DEF})
	set(CULITE_EXAMPLES_3RD_PARTY_I64_DEF ${CULITE_EXAMPLES_3RD_PARTY_I64_DEF} ${INTEL_MKL_I64_DEF})
	# includes
	set(CULITE_EXAMPLES_3RD_PARTY_INC ${CULITE_EXAMPLES_3RD_PARTY_INC} ${INTEL_MKL_INC})
	# libraries
	set(CULITE_EXAMPLES_3RD_PARTY_LIB ${CULITE_EXAMPLES_3RD_PARTY_LIB} ${INTEL_MKL_LIB})
	set(CULITE_EXAMPLES_3RD_PARTY_I64_LIB  ${CULITE_EXAMPLES_3RD_PARTY_I64_LIB} ${INTEL_MKL_I64_LIB})
	# library paths
	set(CULITE_EXAMPLES_3RD_PARTY_LIB_PATH ${CULITE_EXAMPLES_3RD_PARTY_LIB_PATH} ${INTEL_MKL_LIB_DIR} ${INTEL_ICC_LIB_DIR})
	set(CULITE_EXAMPLES_3RD_PARTY_DLL_PATH ${CULITE_EXAMPLES_3RD_PARTY_DLL_PATH} ${INTEL_MKL_DLL_DIR} ${INTEL_ICC_DLL_DIR})

endif()

# 
# cuda setup
#
if((${SIMULICORE_SYSTEM_LINUX} OR ${SIMULICORE_SYSTEM_WIN32}) AND ${SIMULICORE_ARCH_X86_64})

	if(${SIMULICORE_SYSTEM_LINUX})
		include(${SIMULICORE_INSTALL_DIR}/cmake/3rd/cuda.lin.cmake)
	elseif(${SIMULICORE_SYSTEM_WIN32})
		include(${SIMULICORE_INSTALL_DIR}/cmake/3rd/cuda.win.cmake)
	endif()

	# definitions
	set(CULITE_EXAMPLES_3RD_PARTY_DEF ${CULITE_EXAMPLES_3RD_PARTY_DEF} ${NVIDIA_CUDA_DEF} ${NVIDIA_CUDSS_DEF})
	set(CULITE_EXAMPLES_3RD_PARTY_I64_DEF ${CULITE_EXAMPLES_3RD_PARTY_I64_DEF} ${NVIDIA_CUDA_DEF} ${NVIDIA_CUDSS_DEF})
	# includes
	set(CULITE_EXAMPLES_3RD_PARTY_INC ${CULITE_EXAMPLES_3RD_PARTY_INC} ${NVIDIA_CUDA_INC} ${NVIDIA_CUDSS_INC})
	# libraries
	set(CULITE_EXAMPLES_3RD_PARTY_LIB ${CULITE_EXAMPLES_3RD_PARTY_LIB} ${NVIDIA_CUDA_LIB} ${NVIDIA_CUDSS_LIB})
	set(CULITE_EXAMPLES_3RD_PARTY_I64_LIB  ${CULITE_EXAMPLES_3RD_PARTY_I64_LIB} ${NVIDIA_CUDA_LIB} ${NVIDIA_CUDSS_LIB})
	# library paths
	set(CULITE_EXAMPLES_3RD_PARTY_LIB_PATH ${CULITE_EXAMPLES_3RD_PARTY_LIB_PATH} ${NVIDIA_CUDA_LIB_DIR} ${NVIDIA_CUDSS_LIB_DIR})
	set(CULITE_EXAMPLES_3RD_PARTY_DLL_PATH ${CULITE_EXAMPLES_3RD_PARTY_DLL_PATH} ${NVIDIA_CUDA_DLL_DIR} ${NVIDIA_CUDSS_DLL_DIR})

endif()

#-----------------------------------------------
# function that creates wrap script
#-----------------------------------------------
if(${SIMULICORE_SYSTEM_LINUX})
	function(install_script_file)
		set(SCRIPT_FNAME "${CMAKE_CURRENT_BINARY_DIR}/${EXAMPLE}.sh")
		file(WRITE ${SCRIPT_FNAME} "#!/bin/bash\n")
		file(APPEND ${SCRIPT_FNAME} "ROOT_DIR=`dirname $0`\n")
		file(APPEND ${SCRIPT_FNAME} "if [[ -z $" "{LD_LIBRARY_PATH} ]]; then\n")
		file(APPEND ${SCRIPT_FNAME} "  export LD_LIBRARY_PATH=${CULITE_EXAMPLES_3RD_PARTY_LIB_PATH}\n")
		file(APPEND ${SCRIPT_FNAME} "else\n")
		file(APPEND ${SCRIPT_FNAME} "  export LD_LIBRARY_PATH=$" "{LD_LIBRARY_PATH}:${CULITE_EXAMPLES_3RD_PARTY_LIB_PATH}\n")
		file(APPEND ${SCRIPT_FNAME} "fi\n")
		file(APPEND ${SCRIPT_FNAME} "exec $" "{ROOT_DIR}/${EXAMPLE}")
		install(PROGRAMS ${SCRIPT_FNAME} DESTINATION bin)
	endfunction()
endif()

if(${SIMULICORE_SYSTEM_WIN32})
	function(install_script_file)
		set(SCRIPT_FNAME "${CMAKE_CURRENT_BINARY_DIR}/${EXAMPLE}.bat")
		file(WRITE  ${SCRIPT_FNAME} "@echo off\n")
		file(APPEND ${SCRIPT_FNAME} "setlocal\n")
		file(APPEND ${SCRIPT_FNAME} "set PATH=%PATH%;${CULITE_EXAMPLES_3RD_PARTY_DLL_PATH}\n")
		file(APPEND ${SCRIPT_FNAME} "\"%~dp0\\${EXAMPLE}\"")
		install(PROGRAMS ${SCRIPT_FNAME} DESTINATION bin)
	endfunction()
endif()

#-----------------------------------------------
# target/dependencies/installation setup
#-----------------------------------------------
if(${SIMULICORE_SYSTEM_LINUX})
	string(REPLACE ";" ":" CULITE_EXAMPLES_3RD_PARTY_LIB_PATH "${CULITE_EXAMPLES_3RD_PARTY_LIB_PATH}")
endif()

file(GLOB EXAMPLES "src/*.cpp")

foreach(EXAMPLE_CPP_FULL ${EXAMPLES})
	get_filename_component(EXAMPLE_CPP ${EXAMPLE_CPP_FULL} NAME)
	string(REPLACE ".cpp" "" EXAMPLE ${EXAMPLE_CPP})
	add_executable(${EXAMPLE} src/${EXAMPLE_CPP})
	if(${CULITE_EXAMPLES_I64})
		target_include_directories(${EXAMPLE} PRIVATE ${CULITE_EXAMPLES_3RD_PARTY_INC})
		target_compile_definitions(${EXAMPLE} PRIVATE ${CULITE_EXAMPLES_3RD_PARTY_I64_DEF})
		target_link_libraries(${EXAMPLE} ${CULITE_EXAMPLES_3RD_PARTY_I64_LIB})
	else()
		target_include_directories(${EXAMPLE} PRIVATE ${CULITE_EXAMPLES_3RD_PARTY_INC})
		target_compile_definitions(${EXAMPLE} PRIVATE ${CULITE_EXAMPLES_3RD_PARTY_DEF})
		target_link_libraries(${EXAMPLE} ${CULITE_EXAMPLES_3RD_PARTY_LIB})
	endif()
	install(TARGETS ${EXAMPLE} DESTINATION bin)
	install_script_file()
endforeach()

#-----------------------------------------------
# end
#-----------------------------------------------
