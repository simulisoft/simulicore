message(STATUS "Configuring CMPI...")
if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.16.0")
	list(APPEND CMAKE_MESSAGE_INDENT ${SIMULICORE_INDENTATION})
endif()

#-----------------------------------------------
# set module root
#-----------------------------------------------
set(CMPI_ROOT ${CMAKE_CURRENT_SOURCE_DIR})
message(STATUS "CMPI Root Directory: ${CMPI_ROOT}")

#-----------------------------------------------
# sub-directory setup
#-----------------------------------------------
#add_subdirectory(source)
#add_subdirectory(examples)

if(${SIMULICORE_DEVELOPER_MODE})
	add_subdirectory(documentation)
	add_subdirectory(images)
endif()

#-----------------------------------------------
# end
#-----------------------------------------------
