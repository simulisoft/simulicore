#-----------------------------------------------
# set 3rd parties
#-----------------------------------------------

set(CULITE_3RD_PARTY_DEF )
set(CULITE_3RD_PARTY_I64_DEF CULITE_I64)
set(CULITE_3RD_PARTY_INC )
set(CULITE_3RD_PARTY_LIB )
set(CULITE_3RD_PARTY_I64_LIB )

if(${SIMULICORE_ARCH_ARM64})
	message(FATAL "NVIDIA CUDA & CUDSS for arm64 not yet supported")
endif()

#-----------------------------------------------
# cuda/cudss setup
#-----------------------------------------------
set(CULITE_CMAKE_NVIDIA_CUDA false)
if((${SIMULICORE_SYSTEM_LINUX} OR ${SIMULICORE_SYSTEM_WIN32}) AND ${SIMULICORE_ARCH_X86_64})

	if(${SIMULICORE_SYSTEM_LINUX})
		include(${SIMULICORE_ROOT}/cmake/3rd/cuda.lin.cmake)
	elseif(${SIMULICORE_SYSTEM_WIN32})
		include(${SIMULICORE_ROOT}/cmake/3rd/cuda.win.cmake)
	endif()

	message(STATUS "Intel CUDA Root: ${NVIDIA_CUDA_ROOT}")
	message(STATUS "Intel CUDSS Root: ${NVIDIA_CUDSS_ROOT}")

	set(CULITE_CMAKE_NVIDIA_CUDA true)

	# definitions
	set(CULITE_3RD_PARTY_DEF ${CULITE_3RD_PARTY_DEF} CULITE_NVIDIA_CUDA CULITE_NVIDIA_CUDSS)
	set(CULITE_3RD_PARTY_I64_DEF ${CULITE_3RD_PARTY_I64_DEF} CULITE_NVIDIA_CUDA CULITE_NVIDIA_CUDSS)

	# includes
	set(CULITE_3RD_PARTY_INC ${CULITE_3RD_PARTY_INC} ${NVIDIA_CUDA_INC} ${NVIDIA_CUDSS_INC})

	# libraries
	set(CULITE_3RD_PARTY_LIB ${CULITE_3RD_PARTY_LIB} ${NVIDIA_CUDA_LIB} ${NVIDIA_CUDSS_LIB})
	set(CULITE_3RD_PARTY_I64_LIB ${CULITE_3RD_PARTY_I64_LIB} ${NVIDIA_CUDA_LIB} ${NVIDIA_CUDSS_LIB})

endif()

#-----------------------------------------------
# cla3p setup
#-----------------------------------------------
# definitions
#set(CULITE_3RD_PARTY_DEF ${CULITE_3RD_PARTY_DEF})
set(CULITE_3RD_PARTY_I64_DEF ${CULITE_3RD_PARTY_I64_DEF} CLA3P_I64)

# includes
set(CULITE_3RD_PARTY_INC ${CULITE_3RD_PARTY_INC} ${CLA3P_INC})

# libraries
set(CULITE_3RD_PARTY_LIB ${CULITE_3RD_PARTY_LIB} ${CLA3P_LIB})
set(CULITE_3RD_PARTY_I64_LIB ${CULITE_3RD_PARTY_I64_LIB} ${CLA3P_I64_LIB})

#-----------------------------------------------
# sub-directory setup
#-----------------------------------------------
add_subdirectory(culite)

#-----------------------------------------------
# end
#-----------------------------------------------
