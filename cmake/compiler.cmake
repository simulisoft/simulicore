set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

if(${SIMULICORE_SYSTEM_LINUX})
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
		set(CPP_FLAGS_LINUX "")
        if(${SIMULICORE_DEVELOPER_MODE})
			list(APPEND CPP_FLAGS_LINUX -Werror 
			                            -Wpedantic)
        endif()
			list(APPEND CPP_FLAGS_LINUX -m64 
			                            -Wall 
										-Wextra 
										-Wshadow 
										-Wmissing-format-attribute 
										-fdiagnostics-color=always)
	endif()
	add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:${CPP_FLAGS_LINUX}>")
endif()

if(${SIMULICORE_SYSTEM_MACOS})
	execute_process(COMMAND brew --prefix
		OUTPUT_STRIP_TRAILING_WHITESPACE
		OUTPUT_VARIABLE BREW_PREFIX_DIR
		RESULTS_VARIABLE BREW_PREFIX_COMMAND_RET)
	if(${BREW_PREFIX_COMMAND_RET} EQUAL "0")
		message(STATUS "Homebrew directory: ${BREW_PREFIX_DIR}")
	else()
		unset(BREW_PREFIX_DIR)
		message(WARNING "Could not find Homebrew directory")
	endif()
endif()

# TODO: add options for all compilers
