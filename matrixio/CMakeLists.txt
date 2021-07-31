set (MatrixIO_SRCS
    mtx_format.cpp
)
set (MatrixIO_HDRS
    csr_simple_format.h
    mtx_format.h
    csr_format.h
    coo_format.h
    mtx_format.h
    vec_format.h
)

add_library(matrixio ${MatrixIO_SRCS} ${MatrixIO_HDRS})

# INTERFACE include directories to customer target
target_include_directories(matrixio
          INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
          )
