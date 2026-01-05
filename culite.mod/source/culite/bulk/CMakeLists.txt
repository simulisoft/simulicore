#-----------------------------------------------
# source setup
#-----------------------------------------------
set(CULITE_SRC ${CULITE_SRC}
	PARENT_SCOPE)

set(CULITE_BULK_HPP 
	)

#-----------------------------------------------
# installation setup
#-----------------------------------------------
set(CULITE_BULK_HPP_INSTALL include/culite/bulk)

install(FILES ${CULITE_BULK_HPP} DESTINATION ${CULITE_BULK_HPP_INSTALL})
#-----------------------------------------------
# end
#-----------------------------------------------
