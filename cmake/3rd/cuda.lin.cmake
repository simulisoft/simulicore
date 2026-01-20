set(NVIDIA_CUDA_ROOT /opt/nvidia/cuda-toolkit/13.1.1)
set(NVIDIA_CUDSS_ROOT /opt/nvidia/cudss/0.7.1/cuda13)

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
set(NVIDIA_CUDA_LIB_DIR ${NVIDIA_CUDA_ROOT}/lib64)
set(NVIDIA_CUDSS_LIB_DIR ${NVIDIA_CUDSS_ROOT}/lib)

# 
# libraries
#
set(NVIDIA_CUDA_LIB 
	-L${NVIDIA_CUDA_LIB_DIR}
	-lcudart
	-lcublas
	-lcusolver 
	-lcusparse)

set(NVIDIA_CUDSS_LIB 
	-L${NVIDIA_CUDSS_LIB_DIR}
	-lcudss)
