#-----------------------------------------------
# source setup
#-----------------------------------------------
set(CULITE_SRC ${CULITE_SRC}
    PARENT_SCOPE)

set(CULITE_CHECKS_HPP 
    cla3p_forwards.hpp
    )

#-----------------------------------------------
# installation setup
#-----------------------------------------------
set(CULITE_CHECKS_HPP_INSTALL include/culite/checks)

install(FILES ${CULITE_CHECKS_HPP} DESTINATION ${CULITE_CHECKS_HPP_INSTALL})
#-----------------------------------------------
# end
#-----------------------------------------------
