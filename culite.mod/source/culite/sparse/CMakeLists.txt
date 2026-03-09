#-----------------------------------------------
# source setup
#-----------------------------------------------
set(CULITE_SRC ${CULITE_SRC}
    sparse/csr_xxmatrix.cpp
    sparse/csc_xxmatrix.cpp
    PARENT_SCOPE)

set(CULITE_SPARSE_HPP 
    csx_xxcontainer.hpp
    csr_xxmatrix.hpp
    csc_xxmatrix.hpp
    )

#-----------------------------------------------
# installation setup
#-----------------------------------------------
set(CULITE_SPARSE_HPP_INSTALL include/culite/sparse)
install(FILES ${CULITE_SPARSE_HPP} DESTINATION ${CULITE_SPARSE_HPP_INSTALL})

#-----------------------------------------------
# end
#-----------------------------------------------
