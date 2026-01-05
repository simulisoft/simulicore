cmake_minimum_required(VERSION 3.12.0)

project("CLA3P Examples" VERSION "0.1.0" LANGUAGES CXX)

#-----------------------------------------------
# set global definitions
#-----------------------------------------------
option(CLA3P_EXAMPLES_I64 "CLA3P integer precision (32/64bits)" OFF)
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
if(${CLA3P_EXAMPLES_I64})
	set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/i64" CACHE PATH "Installation directory" FORCE)
else()
	set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/i32" CACHE PATH "Installation directory" FORCE)
endif()
message(STATUS "CLA3P Examples Install Directory: ${CMAKE_INSTALL_PREFIX}")

#-----------------------------------------------
# set openmp
#-----------------------------------------------
include(${SIMULICORE_INSTALL_DIR}/cmake/openmp.cmake)

if(${SIMULICORE_SYSTEM_MACOS})
	set(CLA3P_EXAMPLES_3RD_PARTY_INC ${CLA3P_EXAMPLES_3RD_PARTY_INC} ${MACOS_OMP_INC})
	set(CLA3P_EXAMPLES_3RD_PARTY_LIB ${CLA3P_EXAMPLES_3RD_PARTY_LIB} ${MACOS_OMP_LIB})
	set(CLA3P_EXAMPLES_3RD_PARTY_I64_LIB  ${CLA3P_EXAMPLES_3RD_PARTY_I64_LIB} ${MACOS_OMP_LIB})
endif()

#-----------------------------------------------
# set 3rd parties (initialize with CLA3P)
#-----------------------------------------------

set(CLA3P_EXAMPLES_3RD_PARTY_DEF)
set(CLA3P_EXAMPLES_3RD_PARTY_I64_DEF CLA3P_I64)
set(CLA3P_EXAMPLES_3RD_PARTY_INC ${SIMULICORE_INSTALL_DIR}/include)

set(CLA3P_EXAMPLES_3RD_PARTY_LIB_PATH ${CLA3P_EXAMPLES_3RD_PARTY_LIB_PATH} ${SIMULICORE_INSTALL_DIR}/lib)
set(CLA3P_EXAMPLES_3RD_PARTY_DLL_PATH ${CLA3P_EXAMPLES_3RD_PARTY_DLL_PATH} ${SIMULICORE_INSTALL_DIR}/lib)
if(${SIMULICORE_SYSTEM_LINUX})
	set(CLA3P_EXAMPLES_3RD_PARTY_LIB -L${CLA3P_EXAMPLES_3RD_PARTY_LIB_PATH} -lcla3p)
	set(CLA3P_EXAMPLES_3RD_PARTY_I64_LIB -L${CLA3P_EXAMPLES_3RD_PARTY_LIB_PATH} -lcla3p_i64)
elseif(${SIMULICORE_SYSTEM_WIN32})
	set(CLA3P_EXAMPLES_3RD_PARTY_LIB ${CLA3P_EXAMPLES_3RD_PARTY_LIB_PATH}/cla3p.lib)
	set(CLA3P_EXAMPLES_3RD_PARTY_I64_LIB ${CLA3P_EXAMPLES_3RD_PARTY_LIB_PATH}/cla3p_i64.lib)
elseif(${SIMULICORE_SYSTEM_MACOS})
	set(CLA3P_EXAMPLES_3RD_PARTY_LIB -L${CLA3P_EXAMPLES_3RD_PARTY_LIB_PATH} -lcla3p)
	set(CLA3P_EXAMPLES_3RD_PARTY_I64_LIB -L${CLA3P_EXAMPLES_3RD_PARTY_LIB_PATH} -lcla3p_i64)
endif()

# 
# armpl/accelerate setup
#
if(${SIMULICORE_SYSTEM_MACOS})

	#set(CLA3P_EXAMPLES_3RD_PARTY_DEF ${CLA3P_EXAMPLES_3RD_PARTY_DEF} ACCELERATE_NEW_LAPACK)
	#set(CLA3P_EXAMPLES_3RD_PARTY_I64_DEF ${CLA3P_EXAMPLES_3RD_PARTY_I64_DEF} ACCELERATE_NEW_LAPACK ACCELERATE_LAPACK_ILP64)
	#set(CLA3P_EXAMPLES_3RD_PARTY_INC ${CLA3P_EXAMPLES_3RD_PARTY_INC} ${MACOS_OMP_INC})
	#set(CLA3P_EXAMPLES_3RD_PARTY_LIB ${CLA3P_EXAMPLES_3RD_PARTY_LIB} ${MACOS_OMP_LIB} "-framework Accelerate")
	#set(CLA3P_EXAMPLES_3RD_PARTY_I64_LIB  ${CLA3P_EXAMPLES_3RD_PARTY_I64_LIB} ${MACOS_OMP_LIB} "-framework Accelerate")

	include(${SIMULICORE_INSTALL_DIR}/cmake/3rd/armpl.mac.cmake)

	# definitions
	set(CLA3P_EXAMPLES_3RD_PARTY_DEF ${CLA3P_EXAMPLES_3RD_PARTY_DEF} ${ARMPL_DEF})
	set(CLA3P_EXAMPLES_3RD_PARTY_I64_DEF ${CLA3P_EXAMPLES_3RD_PARTY_I64_DEF} ${ARMPL_I64_DEF})
	# includes
	set(CLA3P_EXAMPLES_3RD_PARTY_INC ${CLA3P_EXAMPLES_3RD_PARTY_INC} ${ARMPL_INC})
	# libraries
	set(CLA3P_EXAMPLES_3RD_PARTY_LIB ${CLA3P_EXAMPLES_3RD_PARTY_LIB} ${ARMPL_LIB})
	set(CLA3P_EXAMPLES_3RD_PARTY_I64_LIB ${CLA3P_EXAMPLES_3RD_PARTY_I64_LIB} ${ARMPL_I64_LIB})
	# library paths
	set(CLA3P_EXAMPLES_3RD_PARTY_LIB_PATH ${CLA3P_EXAMPLES_3RD_PARTY_LIB_PATH} ${ARMPL_LIB_DIR})
	set(CLA3P_EXAMPLES_3RD_PARTY_DLL_PATH ${CLA3P_EXAMPLES_3RD_PARTY_DLL_PATH} ${ARMPL_DLL_DIR})

endif()

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
	set(CLA3P_EXAMPLES_3RD_PARTY_DEF ${CLA3P_EXAMPLES_3RD_PARTY_DEF} ${INTEL_MKL_DEF})
	set(CLA3P_EXAMPLES_3RD_PARTY_I64_DEF ${CLA3P_EXAMPLES_3RD_PARTY_I64_DEF} ${INTEL_MKL_I64_DEF})
	# includes
	set(CLA3P_EXAMPLES_3RD_PARTY_INC ${CLA3P_EXAMPLES_3RD_PARTY_INC} ${INTEL_MKL_INC})
	# libraries
	set(CLA3P_EXAMPLES_3RD_PARTY_LIB ${CLA3P_EXAMPLES_3RD_PARTY_LIB} ${INTEL_MKL_LIB})
	set(CLA3P_EXAMPLES_3RD_PARTY_I64_LIB  ${CLA3P_EXAMPLES_3RD_PARTY_I64_LIB} ${INTEL_MKL_I64_LIB})
	# library paths
	set(CLA3P_EXAMPLES_3RD_PARTY_LIB_PATH ${CLA3P_EXAMPLES_3RD_PARTY_LIB_PATH} ${INTEL_MKL_LIB_DIR} ${INTEL_ICC_LIB_DIR})
	set(CLA3P_EXAMPLES_3RD_PARTY_DLL_PATH ${CLA3P_EXAMPLES_3RD_PARTY_DLL_PATH} ${INTEL_MKL_DLL_DIR} ${INTEL_ICC_DLL_DIR})

endif()

#-----------------------------------------------
# function that creates wrap script
#-----------------------------------------------
if(${SIMULICORE_SYSTEM_MACOS})
	function(install_script_file)
		set(SCRIPT_FNAME "${CMAKE_CURRENT_BINARY_DIR}/${EXAMPLE}.sh")
		file(WRITE ${SCRIPT_FNAME} "#!/bin/bash\n")
		file(APPEND ${SCRIPT_FNAME} "ROOT_DIR=`dirname $0`\n")
		file(APPEND ${SCRIPT_FNAME} "if [[ -z $" "{DYLD_FALLBACK_LIBRARY_PATH} ]]; then\n")
		file(APPEND ${SCRIPT_FNAME} "  export DYLD_FALLBACK_LIBRARY_PATH=${CLA3P_EXAMPLES_3RD_PARTY_LIB_PATH}\n")
		file(APPEND ${SCRIPT_FNAME} "else\n")
		file(APPEND ${SCRIPT_FNAME} "  export DYLD_FALLBACK_LIBRARY_PATH=$" "{DYLD_FALLBACK_LIBRARY_PATH}:${CLA3P_EXAMPLES_3RD_PARTY_LIB_PATH}\n")
		file(APPEND ${SCRIPT_FNAME} "fi\n")
		file(APPEND ${SCRIPT_FNAME} "exec $" "{ROOT_DIR}/${EXAMPLE}")
		install(PROGRAMS ${SCRIPT_FNAME} DESTINATION bin)
	endfunction()
endif()

if(${SIMULICORE_SYSTEM_LINUX})
	function(install_script_file)
		set(SCRIPT_FNAME "${CMAKE_CURRENT_BINARY_DIR}/${EXAMPLE}.sh")
		file(WRITE ${SCRIPT_FNAME} "#!/bin/bash\n")
		file(APPEND ${SCRIPT_FNAME} "ROOT_DIR=`dirname $0`\n")
		file(APPEND ${SCRIPT_FNAME} "if [[ -z $" "{LD_LIBRARY_PATH} ]]; then\n")
		file(APPEND ${SCRIPT_FNAME} "  export LD_LIBRARY_PATH=${CLA3P_EXAMPLES_3RD_PARTY_LIB_PATH}\n")
		file(APPEND ${SCRIPT_FNAME} "else\n")
		file(APPEND ${SCRIPT_FNAME} "  export LD_LIBRARY_PATH=$" "{LD_LIBRARY_PATH}:${CLA3P_EXAMPLES_3RD_PARTY_LIB_PATH}\n")
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
		file(APPEND ${SCRIPT_FNAME} "set PATH=%PATH%;${CLA3P_EXAMPLES_3RD_PARTY_DLL_PATH}\n")
		file(APPEND ${SCRIPT_FNAME} "\"%~dp0\\${EXAMPLE}\"")
		install(PROGRAMS ${SCRIPT_FNAME} DESTINATION bin)
	endfunction()
endif()

#-----------------------------------------------
# target/dependencies/installation setup
#-----------------------------------------------
if(${SIMULICORE_SYSTEM_LINUX} OR ${SIMULICORE_SYSTEM_MACOS})
	string(REPLACE ";" ":" CLA3P_EXAMPLES_3RD_PARTY_LIB_PATH "${CLA3P_EXAMPLES_3RD_PARTY_LIB_PATH}")
endif()

file(GLOB EXAMPLES "src/*.cpp")

foreach(EXAMPLE_CPP_FULL ${EXAMPLES})
	get_filename_component(EXAMPLE_CPP ${EXAMPLE_CPP_FULL} NAME)
	string(REPLACE ".cpp" "" EXAMPLE ${EXAMPLE_CPP})
	add_executable(${EXAMPLE} src/${EXAMPLE_CPP})
	if(${CLA3P_EXAMPLES_I64})
		target_include_directories(${EXAMPLE} PRIVATE ${CLA3P_EXAMPLES_3RD_PARTY_INC})
		target_compile_definitions(${EXAMPLE} PRIVATE ${CLA3P_EXAMPLES_3RD_PARTY_I64_DEF})
		target_link_libraries(${EXAMPLE} ${CLA3P_EXAMPLES_3RD_PARTY_I64_LIB})
	else()
		target_include_directories(${EXAMPLE} PRIVATE ${CLA3P_EXAMPLES_3RD_PARTY_INC})
		target_compile_definitions(${EXAMPLE} PRIVATE ${CLA3P_EXAMPLES_3RD_PARTY_DEF})
		target_link_libraries(${EXAMPLE} ${CLA3P_EXAMPLES_3RD_PARTY_LIB})
	endif()
	install(TARGETS ${EXAMPLE} DESTINATION bin)
	install_script_file()
endforeach()

#-----------------------------------------------
# end
#-----------------------------------------------
