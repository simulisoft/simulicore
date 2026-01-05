set(NVIDIA_CUDA_ROOT "C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v13.1")
set(NVIDIA_CUDSS_ROOT "C:/Program Files/NVIDIA cuDSS/v0.7")

#
# definitions
#

#
# include
#
set(NVIDIA_CUDA_INC ${NVIDIA_CUDA_ROOT}/include)
set(NVIDIA_CUDSS_INC ${NVIDIA_CUDSS_ROOT}/include)

#
# library paths
#
set(NVIDIA_CUDA_LIB_DIR ${NVIDIA_CUDA_ROOT}/lib/x64)
set(NVIDIA_CUDA_DLL_DIR ${NVIDIA_CUDA_ROOT}/bin/x64)
set(NVIDIA_CUDSS_LIB_DIR ${NVIDIA_CUDSS_ROOT}/lib/13)
set(NVIDIA_CUDSS_DLL_DIR ${NVIDIA_CUDSS_ROOT}/bin/13)

#
# libraries
#
set(NVIDIA_CUDA_LIB 
	${NVIDIA_CUDA_LIB_DIR}/cudart.lib 
	${NVIDIA_CUDA_LIB_DIR}/cublas.lib 
	${NVIDIA_CUDA_LIB_DIR}/cusolver.lib 
	${NVIDIA_CUDA_LIB_DIR}/cusparse.lib 
	${NVIDIA_CUDSS_LIB_DIR}/cudss.lib
	${NVIDIA_CUDSS_LIB_DIR}/cudss_mtlayer_vcomp140.lib)
