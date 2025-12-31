#-----------------------------------------------
# source setup
#-----------------------------------------------
set(CLA3P_SRC ${CLA3P_SRC}
	generic/ownership.cpp
	PARENT_SCOPE)

set(CLA3P_GENERIC_HPP 
	meta1d.hpp
	meta2d.hpp
	matrix_meta.hpp
	coord.hpp
	tuple.hpp
	guard.hpp
	ownership.hpp
	type_traits.hpp
	)

#-----------------------------------------------
# installation setup
#-----------------------------------------------
set(CLA3P_GENERIC_HPP_INSTALL include/cla3p/generic)

install(FILES ${CLA3P_GENERIC_HPP} DESTINATION ${CLA3P_GENERIC_HPP_INSTALL})
#-----------------------------------------------
# end
#-----------------------------------------------
