#-----------------------------------------------
# source setup
#-----------------------------------------------
set(CLA3P_SRC ${CLA3P_SRC}
	lra/lra_xxmatrix.cpp
	lra/rank_moderator.cpp
	lra/lra_math.cpp
	PARENT_SCOPE)

set(CLA3P_LRA_HPP 
	lra_xxmatrix.hpp
	rank_moderator.hpp
	lra_math.hpp
	)

#-----------------------------------------------
# installation setup
#-----------------------------------------------
set(CLA3P_LRA_HPP_INSTALL include/cla3p/lra)

install(FILES ${CLA3P_LRA_HPP} DESTINATION ${CLA3P_LRA_HPP_INSTALL})
#-----------------------------------------------
# end
#-----------------------------------------------
