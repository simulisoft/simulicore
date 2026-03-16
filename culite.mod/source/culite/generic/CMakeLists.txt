#-----------------------------------------------
# source setup
#-----------------------------------------------
set(CULITE_SRC ${CULITE_SRC}
    generic/cublas_handler.cpp
    generic/cuda_stream.cpp
    generic/cusolver_handler.cpp
    generic/cusparse_handler.cpp
    PARENT_SCOPE)

set(CULITE_GENERIC_HPP 
    guard.hpp
    meta1d.hpp
    meta2d.hpp
    matrix_meta.hpp
    cuda_stream.hpp
    cublas_handler.hpp
    cusolver_handler.hpp
    cusparse_handler.hpp
    )

#-----------------------------------------------
# installation setup
#-----------------------------------------------
set(CULITE_GENERIC_HPP_INSTALL include/culite/generic)

install(FILES ${CULITE_GENERIC_HPP} DESTINATION ${CULITE_GENERIC_HPP_INSTALL})
#-----------------------------------------------
# end
#-----------------------------------------------
