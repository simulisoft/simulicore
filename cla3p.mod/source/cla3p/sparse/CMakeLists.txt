#-----------------------------------------------
# source setup
#-----------------------------------------------
set(CLA3P_SRC ${CLA3P_SRC}
	sparse/csc_xxmatrix.cpp
	sparse/coo_xxmatrix.cpp
	PARENT_SCOPE)

set(CLA3P_SPARSE_HPP 
	csc_xxcontainer_base.hpp
	csc_xxcontainer.hpp
	csc_xxmatrix.hpp
	coo_xxmatrix.hpp
	)

#-----------------------------------------------
# installation setup
#-----------------------------------------------
set(CLA3P_SPARSE_HPP_INSTALL include/cla3p/sparse)

install(FILES ${CLA3P_SPARSE_HPP} DESTINATION ${CLA3P_SPARSE_HPP_INSTALL})
#-----------------------------------------------
# end
#-----------------------------------------------
