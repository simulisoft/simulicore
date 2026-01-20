set(NVIDIA_CUDA_LICENSE_FILES 
    EULA.txt
	)

#-----------------------------------------------
# installation setup
#-----------------------------------------------
set(NVIDIA_CUDA_LICENSE_INSTALL license/nvidia_cuda)

install(FILES ${NVIDIA_CUDA_LICENSE_FILES} DESTINATION ${NVIDIA_CUDA_LICENSE_INSTALL})
