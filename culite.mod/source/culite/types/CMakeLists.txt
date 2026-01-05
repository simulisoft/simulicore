#-----------------------------------------------
# source setup
#-----------------------------------------------
set(CULITE_SRC ${CULITE_SRC}
	#types/integer.cpp
	types/scalar.cpp
	#types/enums.cpp
	#types/operation.cpp
	#types/property.cpp
	#types/literals.cpp
	PARENT_SCOPE)

set(CULITE_TYPES_HPP 
	integer.hpp
	scalar.hpp
	#enums.hpp
	#operation.hpp
	#property.hpp
	#literals.hpp
	)

#-----------------------------------------------
# installation setup
#-----------------------------------------------
set(CULITE_TYPES_HPP_INSTALL include/culite/types)

install(FILES ${CULITE_TYPES_HPP} DESTINATION ${CULITE_TYPES_HPP_INSTALL})
#-----------------------------------------------
# end
#-----------------------------------------------
