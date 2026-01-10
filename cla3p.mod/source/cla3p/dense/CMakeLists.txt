#-----------------------------------------------
# source setup
#-----------------------------------------------
set(CLA3P_SRC ${CLA3P_SRC}
	dense/dns_xivector.cpp
	dense/dns_xxvector.cpp
	dense/dns_cxvector.cpp
	dense/dns_xxmatrix.cpp
	dense/dns_cxmatrix.cpp
	PARENT_SCOPE)

set(CLA3P_DENSE_HPP 
	dns_xxcontainer_base.hpp
	dns_xxcontainer.hpp
	dns_xivector.hpp
	dns_xxvector.hpp
	dns_cxvector.hpp
	dns_xxmatrix.hpp
	dns_cxmatrix.hpp
	)

#-----------------------------------------------
# installation setup
#-----------------------------------------------
set(CLA3P_DENSE_HPP_INSTALL include/cla3p/dense)

install(FILES ${CLA3P_DENSE_HPP} DESTINATION ${CLA3P_DENSE_HPP_INSTALL})
#-----------------------------------------------
# end
#-----------------------------------------------
