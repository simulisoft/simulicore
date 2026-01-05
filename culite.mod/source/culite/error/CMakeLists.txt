#-----------------------------------------------
# source setup
#-----------------------------------------------
set(CULITE_SRC ${CULITE_SRC}
	error/exceptions.cpp
	error/literals.cpp
	error/cuda.cpp
	PARENT_SCOPE)

set(CULITE_ERROR_HPP 
	exceptions.hpp
	cuda.hpp
	)

#-----------------------------------------------
# installation setup
#-----------------------------------------------
set(CULITE_ERROR_HPP_INSTALL include/culite/error)

install(FILES ${CULITE_ERROR_HPP} DESTINATION ${CULITE_ERROR_HPP_INSTALL})
#-----------------------------------------------
# end
#-----------------------------------------------
