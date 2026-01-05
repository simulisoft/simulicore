#-----------------------------------------------
# source setup
#-----------------------------------------------
set(CULITE_SRC ${CULITE_SRC}
	proxies/cublas_proxy.cpp
	PARENT_SCOPE)

set(CULITE_PROXY_HPP 
	cublas_proxy.hpp
	)

#-----------------------------------------------
# installation setup
#-----------------------------------------------
set(CULITE_PROXY_HPP_INSTALL include/culite/proxies)

install(FILES ${CULITE_PROXY_HPP} DESTINATION ${CULITE_PROXY_HPP_INSTALL})
#-----------------------------------------------
# end
#-----------------------------------------------
