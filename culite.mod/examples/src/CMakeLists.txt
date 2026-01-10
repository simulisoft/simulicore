#-----------------------------------------------
# source setup
#-----------------------------------------------
set(CULITE_EXAMPLES_SRC
	ex01a_dense_vector_create.cpp
	)

#-----------------------------------------------
# installation setup
#-----------------------------------------------
set(CULITE_EXAMPLES_INSTALL examples/culite/src)

install(FILES ${CULITE_EXAMPLES_SRC} DESTINATION ${CULITE_EXAMPLES_INSTALL})

#-----------------------------------------------
# end
#-----------------------------------------------
