#-----------------------------------------------
# source setup
#-----------------------------------------------
set(CLA3P_SRC ${CLA3P_SRC}
	virtuals/virtual_rowvec.cpp
	virtuals/virtual_transpose.cpp
	PARENT_SCOPE)

set(CLA3P_VIRTUALS_HPP 
	virtual_expression.hpp
	virtual_object.hpp
	virtual_rowvec.hpp
	virtual_transpose.hpp
	virtual_conjugate.hpp
	virtual_scale.hpp
	virtual_plus_minus.hpp
	virtual_outer.hpp
	virtual_product.hpp
	virtual_product_specs_impl.hpp
	)

#-----------------------------------------------
# installation setup
#-----------------------------------------------
set(CLA3P_VIRTUALS_HPP_INSTALL include/cla3p/virtuals)

install(FILES ${CLA3P_VIRTUALS_HPP} DESTINATION ${CLA3P_VIRTUALS_HPP_INSTALL})
#-----------------------------------------------
# end
#-----------------------------------------------
