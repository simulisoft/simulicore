#-----------------------------------------------
# source setup
#-----------------------------------------------
set(CULITE_SRC ${CULITE_SRC}
	dense/dns_xivector.cpp
	dense/dns_xxvector.cpp
	#dense/dns_cxvector.cpp
	#dense/dns_xxmatrix.cpp
	#dense/dns_cxmatrix.cpp
	PARENT_SCOPE)

set(CULITE_DENSE_HPP 
	dns_xivector.hpp
	dns_xxvector.hpp
	#dns_cxvector.hpp
	#dns_xxmatrix.hpp
	#dns_cxmatrix.hpp
	)

#-----------------------------------------------
# installation setup
#-----------------------------------------------
set(CULITE_DENSE_HPP_INSTALL include/culite/dense)

install(FILES ${CULITE_DENSE_HPP} DESTINATION ${CULITE_DENSE_HPP_INSTALL})
#-----------------------------------------------
# end
#-----------------------------------------------
