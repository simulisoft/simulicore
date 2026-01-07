#-----------------------------------------------
# source setup
#-----------------------------------------------
set(CULITE_SRC ${CULITE_SRC}
	support/imalloc.cpp
	PARENT_SCOPE)

set(CULITE_SUPPORT_HPP 
	imalloc.hpp
	)

#-----------------------------------------------
# installation setup
#-----------------------------------------------
set(CULITE_SUPPORT_HPP_INSTALL include/culite/support)

install(FILES ${CULITE_SUPPORT_HPP} DESTINATION ${CULITE_SUPPORT_HPP_INSTALL})
#-----------------------------------------------
# end
#-----------------------------------------------
