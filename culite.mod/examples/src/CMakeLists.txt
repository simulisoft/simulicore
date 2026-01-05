#-----------------------------------------------
# source setup
#-----------------------------------------------
set(CULITE_EXAMPLES_SRC
	ex01a_dense_vector_create.cpp
	ex01c_dense_vector_create_from_aux_data.cpp
	ex02a_dense_matrix_create.cpp
	ex02d_dense_matrix_create_from_aux_data.cpp
	ex05a_solving_linear_systems_lu.cpp
	)

#-----------------------------------------------
# installation setup
#-----------------------------------------------
set(CULITE_EXAMPLES_INSTALL examples/culite/src)

install(FILES ${CULITE_EXAMPLES_SRC} DESTINATION ${CULITE_EXAMPLES_INSTALL})

#-----------------------------------------------
# end
#-----------------------------------------------
