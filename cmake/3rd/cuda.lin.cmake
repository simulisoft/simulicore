set(NVIDIA_CUDA_ROOT /usr/local/cuda-13.1)
set(NVIDIA_CUDSS_ROOT /usr/lib/x86_64-linux-gnu/libcudss)

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
set(NVIDIA_CUDSS_LIB_DIR ${NVIDIA_CUDSS_ROOT}/13)

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
