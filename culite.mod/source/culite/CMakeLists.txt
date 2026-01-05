#-----------------------------------------------
# set library name
#-----------------------------------------------
set(CULITE_INC "${CMAKE_CURRENT_SOURCE_DIR}/.." CACHE INTERNAL "The culite includes")
set(CULITE_LIB "culite" CACHE INTERNAL "The lp64 culite library")
set(CULITE_I64_LIB "${CULITE_LIB}_i64" CACHE INTERNAL "The ilp64 culite library")

#-----------------------------------------------
# source setup
#-----------------------------------------------
set(CULITE_SRC
	)

set(CULITE_HPP
#	types.hpp
#	error.hpp
#	support.hpp
#	dense.hpp
#	sparse.hpp
#	perms.hpp
#	virtuals.hpp
#	algebra.hpp
#	linsol.hpp
#	svd.hpp
#	qr.hpp
#	lra.hpp
	)

#-----------------------------------------------
# sub-directory setup
#-----------------------------------------------
#add_subdirectory(types)
#add_subdirectory(error)
#add_subdirectory(bulk)
#add_subdirectory(checks)
#add_subdirectory(proxies)
#add_subdirectory(support)
#add_subdirectory(generic)
#add_subdirectory(dense)
#add_subdirectory(sparse)
#add_subdirectory(perms)
#add_subdirectory(virtuals)
#add_subdirectory(algebra)
#add_subdirectory(linsol)
#add_subdirectory(svd)
#add_subdirectory(qr)
#add_subdirectory(lra)

#-----------------------------------------------
# target setup
#-----------------------------------------------
add_library(${CULITE_LIB} SHARED ${CULITE_SRC})
add_library(${CULITE_I64_LIB} SHARED ${CULITE_SRC})

#-----------------------------------------------
# dependencies setup
#-----------------------------------------------
target_compile_definitions(${CULITE_LIB} PRIVATE ${CULITE_3RD_PARTY_DEF})
target_compile_definitions(${CULITE_I64_LIB} PRIVATE ${CULITE_3RD_PARTY_I64_DEF})

target_include_directories(${CULITE_LIB} PRIVATE ${CULITE_INC} ${CULITE_3RD_PARTY_INC})
target_include_directories(${CULITE_I64_LIB} PRIVATE ${CULITE_INC} ${CULITE_3RD_PARTY_INC})

target_link_libraries(${CULITE_LIB} ${CULITE_3RD_PARTY_LIB})
target_link_libraries(${CULITE_I64_LIB} ${CULITE_3RD_PARTY_I64_LIB})

#-----------------------------------------------
# installation setup
#-----------------------------------------------
set(CULITE_HPP_INSTALL include/culite)
set(CULITE_LIB_INSTALL lib)

install(FILES ${CULITE_HPP} DESTINATION ${CULITE_HPP_INSTALL})
install(TARGETS ${CULITE_LIB} DESTINATION ${CULITE_LIB_INSTALL})
install(TARGETS ${CULITE_I64_LIB} DESTINATION ${CULITE_LIB_INSTALL})
#-----------------------------------------------
# end
#-----------------------------------------------
