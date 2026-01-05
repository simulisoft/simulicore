cmake_minimum_required(VERSION 3.12.0)
cmake_policy(SET CMP0074 NEW)

set(SIMULICORE_NAME "SimuliCore")
set(SIMULICORE_VERSION "0.1.0")
project(${SIMULICORE_NAME} VERSION ${SIMULICORE_VERSION} LANGUAGES CXX)

#-----------------------------------------------
# welcome message
#-----------------------------------------------
message(STATUS "Configuring ${SIMULICORE_NAME} Suite...")
if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.16.0")
	set(SIMULICORE_INDENTATION "    ")
	list(APPEND CMAKE_MESSAGE_INDENT ${SIMULICORE_INDENTATION})
endif()

#-----------------------------------------------
# set project root
#-----------------------------------------------
set(SIMULICORE_ROOT ${CMAKE_CURRENT_SOURCE_DIR})
message(STATUS "${SIMULICORE_NAME} Root Directory: ${SIMULICORE_ROOT}")

#-----------------------------------------------
# set global project definitions
#-----------------------------------------------
option(SIMULICORE_DEVELOPER_MODE "Install all ${SIMULICORE_NAME} featues for developers" OFF)

#-----------------------------------------------
# build & install setup
#-----------------------------------------------
include(${SIMULICORE_ROOT}/cmake/os.cmake)
include(${SIMULICORE_ROOT}/cmake/compiler.cmake)
include(${SIMULICORE_ROOT}/cmake/build_type.cmake)
include(${SIMULICORE_ROOT}/cmake/install.cmake)

message(STATUS "Detected OS: ${CMAKE_SYSTEM_NAME} (${CMAKE_SYSTEM_PROCESSOR})")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "${SIMULICORE_NAME} Install Directory: ${CMAKE_INSTALL_PREFIX}")

#-----------------------------------------------
# sub-module setup
#-----------------------------------------------
add_subdirectory(cmake)
add_subdirectory(licensing)
add_subdirectory(cla3p.mod)
add_subdirectory(cmpi.mod)
if(${SIMULICORE_SYSTEM_LINUX} OR ${SIMULICORE_SYSTEM_WIN32})
add_subdirectory(culite.mod)
endif()

if(${SIMULICORE_DEVELOPER_MODE})
	add_subdirectory(sample)
endif()

#-----------------------------------------------
# end
#-----------------------------------------------
