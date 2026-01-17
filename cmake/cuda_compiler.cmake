if(${SIMULICORE_SYSTEM_LINUX})
	set(CMAKE_CUDA_COMPILER "${NVIDIA_CUDA_ROOT}/bin/nvcc")
	if(CMAKE_CUDA_COMPILER_ID STREQUAL "NVIDIA")
		set(CUDA_FLAGS_LINUX "")
		if(${SIMULICORE_DEVELOPER_MODE})
			list(APPEND CUDA_FLAGS_LINUX -Xcompiler=-Werror
			                             -Xcompiler=-Wpedantic)
		endif()
		list(APPEND CUDA_FLAGS_LINUX -m64 
		                             -Xcompiler=-Wall
									 -Xcompiler=-Wextra
									 -Xcompiler=-Wshadow
									 -Xcompiler=-Wmissing-format-attribute
									 -Xcompiler=-fdiagnostics-color=always)	
		add_compile_options("$<$<COMPILE_LANGUAGE:CUDA>:${CUDA_FLAGS_LINUX}>")
	endif()
endif()

if(${SIMULICORE_SYSTEM_WINDOWS})
	set(CMAKE_CUDA_COMPILER "${NVIDIA_CUDA_ROOT}/bin/nvcc")
	if(CMAKE_CUDA_COMPILER_ID STREQUAL "NVIDIA")
		# TODO: add windows cuda compile options
	endif()
endif()

enable_language(CUDA)

set(CMAKE_CUDA_ARCHITECTURES 
    75      # Turing (RTX 20-series, T4)
    80 86   # Ampere (RTX 30-series, A100)
    89 90   # Ada Lovelace / Hopper (RTX 40-series, H100)
    100 120 # Blackwell (RTX 50-series, B100/B200)
    # 130?  # Rubin (Vera Rubin / R100 - check latest toolkit support)
)
