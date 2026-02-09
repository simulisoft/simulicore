#-----------------------------------------------
# source setup
#-----------------------------------------------
set(CULITE_SRC ${CULITE_SRC}
	#virtuals/virtual_rowvec.cpp
    #virtuals/virtual_strided_block.cpp
	virtuals/virtual_transpose.cpp
	PARENT_SCOPE)

set(CULITE_VIRTUALS_HPP 
	virtual_expression.hpp
	virtual_object.hpp
	#virtual_rowvec.hpp
    #virtual_strided_block.hpp
	virtual_transpose.hpp
	#virtual_conjugate.hpp
	#virtual_scale.hpp
	#virtual_plus_minus.hpp
	#virtual_outer.hpp
	#virtual_product.hpp
	#virtual_product_specs_impl.hpp
	)

#-----------------------------------------------
# installation setup
#-----------------------------------------------
set(CULITE_VIRTUALS_HPP_INSTALL include/culite/virtuals)

install(FILES ${CULITE_VIRTUALS_HPP} DESTINATION ${CULITE_VIRTUALS_HPP_INSTALL})
#-----------------------------------------------
# end
#-----------------------------------------------
