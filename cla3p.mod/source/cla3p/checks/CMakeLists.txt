#-----------------------------------------------
# source setup
#-----------------------------------------------
set(CLA3P_SRC ${CLA3P_SRC}
    checks/block_ops_checks.cpp
    checks/dns_checks.cpp
    checks/coo_checks.cpp
    checks/lapack_checks.cpp
    checks/matrix_math_checks.cpp
    checks/perm_checks.cpp
    checks/transp_checks.cpp
    PARENT_SCOPE)

set(CLA3P_CHECKS_HPP 
    basic_checks.hpp
    hermitian_coeff_checks.hpp
    block_ops_checks.hpp
    dns_checks.hpp
    csx_checks.hpp
    outer_checks.hpp
    matrix_math_checks.hpp
    decomp_xx_checks.hpp
    solve_checks.hpp
    )

#-----------------------------------------------
# installation setup
#-----------------------------------------------
set(CLA3P_CHECKS_HPP_INSTALL include/cla3p/checks)

install(FILES ${CLA3P_CHECKS_HPP} DESTINATION ${CLA3P_CHECKS_HPP_INSTALL})
#-----------------------------------------------
# end
#-----------------------------------------------
