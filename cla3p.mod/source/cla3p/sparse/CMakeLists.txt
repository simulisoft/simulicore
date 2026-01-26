#-----------------------------------------------
# source setup
#-----------------------------------------------
set(CLA3P_SRC ${CLA3P_SRC}
    sparse/csr_xxmatrix.cpp
	sparse/csc_xxmatrix.cpp
	sparse/coo_xxmatrix.cpp
	PARENT_SCOPE)

set(CLA3P_SPARSE_HPP 
    csx_xxcontainer_base.hpp
    csx_xxcontainer.hpp
    csr_xxmatrix.hpp
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
