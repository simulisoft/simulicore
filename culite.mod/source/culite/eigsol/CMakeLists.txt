#-----------------------------------------------
# source setup
#-----------------------------------------------
set(CULITE_SRC ${CULITE_SRC}
	eigsol/lapack_geev.cpp
	PARENT_SCOPE)

set(CULITE_EIGSOL_HPP 
	lapack_geev.hpp
	)

#-----------------------------------------------
# installation setup
#-----------------------------------------------
set(CULITE_EIGSOL_HPP_INSTALL include/culite/eigsol)

install(FILES ${CULITE_EIGSOL_HPP} DESTINATION ${CULITE_EIGSOL_HPP_INSTALL})
#-----------------------------------------------
# end
#-----------------------------------------------
