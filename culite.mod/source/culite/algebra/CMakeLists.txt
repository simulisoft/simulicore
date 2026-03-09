#-----------------------------------------------
# source setup
#-----------------------------------------------
set(CULITE_SRC ${CULITE_SRC}
    algebra/functional_inner.cpp
    algebra/functional_outer.cpp
    algebra/functional_add.cpp
    algebra/functional_update.cpp
    algebra/functional_multmv.cpp
    algebra/functional_multmm.cpp
    PARENT_SCOPE)

set(CULITE_ALGEBRA_HPP 
    functional_inner.hpp
    functional_outer.hpp
    functional_add.hpp
    functional_update.hpp
    functional_multmv.hpp
    functional_multmm.hpp
    operators_mult_inner.hpp
    operators_mult_outer.hpp
    operators_scale.hpp
    operators_update.hpp
    operators_add.hpp
    operators_multxx.hpp
    operators_multmv.hpp
    operators_multmm.hpp
    )

#-----------------------------------------------
# installation setup
#-----------------------------------------------
set(CULITE_ALGEBRA_HPP_INSTALL include/culite/algebra)

install(FILES ${CULITE_ALGEBRA_HPP} DESTINATION ${CULITE_ALGEBRA_HPP_INSTALL})
#-----------------------------------------------
# end
#-----------------------------------------------
