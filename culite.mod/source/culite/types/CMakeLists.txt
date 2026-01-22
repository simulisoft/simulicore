#-----------------------------------------------
# source setup
#-----------------------------------------------
set(CULITE_SRC ${CULITE_SRC}
	types/integer.cpp
	types/scalar.cpp
	types/enums.cpp
	PARENT_SCOPE)

set(CULITE_TYPES_HPP 
	cuda_api.hpp
	traits.hpp
	integer.hpp
	complex.hpp
	scalar.hpp
	enums.hpp
	)

set(CULITE_CUSRC ${CULITE_CUSRC}
	PARENT_SCOPE)

set(CULITE_TYPES_CUH 
	)

#-----------------------------------------------
# installation setup
#-----------------------------------------------
set(CULITE_TYPES_HPP_INSTALL include/culite/types)

install(FILES ${CULITE_TYPES_HPP} DESTINATION ${CULITE_TYPES_HPP_INSTALL})
install(FILES ${CULITE_TYPES_CUH} DESTINATION ${CULITE_TYPES_HPP_INSTALL})
#-----------------------------------------------
# end
#-----------------------------------------------
