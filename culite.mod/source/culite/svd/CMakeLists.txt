#-----------------------------------------------
# source setup
#-----------------------------------------------
set(CULITE_SRC ${CULITE_SRC}
    svd/default_svd.cpp
    PARENT_SCOPE)

set(CULITE_SVD_HPP 
    default_svd.hpp
    )

#-----------------------------------------------
# installation setup
#-----------------------------------------------
set(CULITE_SVD_HPP_INSTALL include/culite/svd)

install(FILES ${CULITE_SVD_HPP} DESTINATION ${CULITE_SVD_HPP_INSTALL})
#-----------------------------------------------
# end
#-----------------------------------------------
