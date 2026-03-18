#-----------------------------------------------
# source setup
#-----------------------------------------------
set(CLA3P_SRC ${CLA3P_SRC}
    eigsol/lapack_syev.cpp
    eigsol/lapack_geev.cpp
    PARENT_SCOPE)

set(CLA3P_EIGSOL_HPP 
    lapack_syev.hpp
    lapack_geev.hpp
    )

#-----------------------------------------------
# installation setup
#-----------------------------------------------
set(CLA3P_EIGSOL_HPP_INSTALL include/cla3p/eigsol)

install(FILES ${CLA3P_EIGSOL_HPP} DESTINATION ${CLA3P_EIGSOL_HPP_INSTALL})
#-----------------------------------------------
# end
#-----------------------------------------------
