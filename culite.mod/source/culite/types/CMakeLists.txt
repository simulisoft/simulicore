#-----------------------------------------------
# source setup
#-----------------------------------------------
set(CULITE_SRC ${CULITE_SRC}
	#types/integer.cpp
	types/scalar.cpp
	#types/enums.cpp
	PARENT_SCOPE)

set(CULITE_TYPES_HPP 
	traits.hpp
	integer.hpp
	scalar.hpp
	#enums.hpp
	)

#-----------------------------------------------
# installation setup
#-----------------------------------------------
set(CULITE_TYPES_HPP_INSTALL include/culite/types)

install(FILES ${CULITE_TYPES_HPP} DESTINATION ${CULITE_TYPES_HPP_INSTALL})
#-----------------------------------------------
# end
#-----------------------------------------------
