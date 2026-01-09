#-----------------------------------------------
# source setup
#-----------------------------------------------
set(CULITE_SRC ${CULITE_SRC}
	PARENT_SCOPE)

set(CULITE_GENERIC_HPP 
	)

#-----------------------------------------------
# installation setup
#-----------------------------------------------
set(CULITE_GENERIC_HPP_INSTALL include/culite/generic)

install(FILES ${CULITE_GENERIC_HPP} DESTINATION ${CULITE_GENERIC_HPP_INSTALL})
#-----------------------------------------------
# end
#-----------------------------------------------
