set(NVIDIA_CUDSS_LICENSE_FILES 
    LICENSE
	)

#-----------------------------------------------
# installation setup
#-----------------------------------------------
set(NVIDIA_CUDSS_LICENSE_INSTALL license/nvidia_cudss)

install(FILES ${NVIDIA_CUDSS_LICENSE_FILES} DESTINATION ${NVIDIA_CUDSS_LICENSE_INSTALL})
