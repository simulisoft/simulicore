#-----------------------------------------------
# source setup
#-----------------------------------------------
set(CULITE_SRC ${CULITE_SRC}
    PARENT_SCOPE)

set(CULITE_BULK_HPP 
    dns1d.hpp
    dns1d_impl.hpp
    dns2d.hpp
    dns2d_impl.hpp
    )

set(CULITE_CUSRC ${CULITE_CUSRC}
    bulk/dns1d_impl.cu
    bulk/dns2d_impl.cu
    bulk/csx_impl.cu
    PARENT_SCOPE)

set(CULITE_TYPES_CUH 
    )

#-----------------------------------------------
# installation setup
#-----------------------------------------------
set(CULITE_BULK_HPP_INSTALL include/culite/bulk)

install(FILES ${CULITE_BULK_HPP} DESTINATION ${CULITE_BULK_HPP_INSTALL})
#-----------------------------------------------
# end
#-----------------------------------------------
