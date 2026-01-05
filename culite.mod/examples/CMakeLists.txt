#-----------------------------------------------
# sub-directory setup
#-----------------------------------------------
add_subdirectory(src)

#-----------------------------------------------
# installation setup
#-----------------------------------------------
set(CULITE_EXAMPLES_INSTALL examples/culite)

install(FILES cmake_list_install/CMakeLists.txt DESTINATION ${CULITE_EXAMPLES_INSTALL})

if(${SIMULICORE_SYSTEM_LINUX} OR ${SIMULICORE_SYSTEM_MACOS})
	install(PROGRAMS example_builder.sh DESTINATION ${CULITE_EXAMPLES_INSTALL})
	install(PROGRAMS example_builder_i32.sh DESTINATION ${CULITE_EXAMPLES_INSTALL})
	install(PROGRAMS example_builder_i64.sh DESTINATION ${CULITE_EXAMPLES_INSTALL})
	if(${SIMULICORE_DEVELOPER_MODE})
		install(PROGRAMS example_executor.sh DESTINATION ${CULITE_EXAMPLES_INSTALL})
		install(PROGRAMS example_executor_i32.sh DESTINATION ${CULITE_EXAMPLES_INSTALL})
		install(PROGRAMS example_executor_i64.sh DESTINATION ${CULITE_EXAMPLES_INSTALL})
	endif()
endif()

#-----------------------------------------------
# end
#-----------------------------------------------
