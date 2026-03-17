#-----------------------------------------------
# source setup
#-----------------------------------------------
set(CULITE_EXAMPLES_SRC
	ex01a_dense_vector_create.cpp
	ex01c_dense_vector_create_from_aux_data.cpp
	ex02a_dense_matrix_create.cpp
	ex02d_dense_matrix_create_from_aux_data.cpp
	ex02k_dense_matrix_algebra_scale.cpp
	ex02l_dense_matrix_algebra_add.cpp
	ex02m_dense_matrix_algebra_vmult.cpp
	ex02n_dense_matrix_algebra_mmult.cpp
	ex02o_dense_matrix_algebra_mtmult.cpp
	ex05a_solving_linear_systems_lu.cpp
	ex06a_sparse_matrix_create.cpp
	ex06d_sparse_matrix_create_from_aux_data.cpp
	ex06k_sparse_matrix_algebra_scale.cpp
	ex06l_sparse_matrix_algebra_add.cpp
	ex06m_sparse_matrix_algebra_vmult.cpp
    ex06n_sparse_matrix_algebra_mmult.cpp
	ex07a_eigenproblem_geev.cpp
	ex08a_sv_decomposition.cpp
	)

#-----------------------------------------------
# installation setup
#-----------------------------------------------
set(CULITE_EXAMPLES_INSTALL examples/culite/src)

install(FILES ${CULITE_EXAMPLES_SRC} DESTINATION ${CULITE_EXAMPLES_INSTALL})

#-----------------------------------------------
# end
#-----------------------------------------------
