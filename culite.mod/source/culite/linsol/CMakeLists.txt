#-----------------------------------------------
# source setup
#-----------------------------------------------
set(CULITE_SRC ${CULITE_SRC}
	linsol/lapack_base.cpp
	PARENT_SCOPE)

set(CULITE_LINSOL_HPP 
	lapack_base.hpp
	lapack_lu.hpp
	)

#-----------------------------------------------
# installation setup
#-----------------------------------------------
set(CULITE_LINSOL_HPP_INSTALL include/culite/linsol)

install(FILES ${CULITE_LINSOL_HPP} DESTINATION ${CULITE_LINSOL_HPP_INSTALL})
#-----------------------------------------------
# end
#-----------------------------------------------
