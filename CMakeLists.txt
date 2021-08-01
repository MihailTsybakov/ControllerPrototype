cmake_minimum_required(VERSION 3.1)

project(ControllerPrototype)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake" ${CMAKE_MODULE_PATH})

set(ControllerPrototype_H 
	processtask.h
	processcontroller.h
	solvertest.h
	petsccgtest.h
)

set(ControllerPrototype_SRC 
	main.cpp
	processcontroller.cpp
	solvertest.cpp
	petsccgtest.cpp
)

find_package(MPILibraries REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

find_package(MKL REQUIRED)
include_directories(${MKL_INCLUDES})
link_directories(${MKL_LIBRARYDIR})
link_directories(${MKL_BINARYDIR})
link_directories(${OMP_LIBRARYDIR})

find_package(PETSc REQUIRED)

add_subdirectory(matrixio)

add_executable(${PROJECT_NAME} ${ControllerPrototype_H} ${ControllerPrototype_SRC})
target_include_directories(${PROJECT_NAME} PRIVATE ${PETSC_INCLUDES} ${PETSC_CONF_INCLUDES_64})
target_link_libraries(${PROJECT_NAME} LINK_PRIVATE ${PETSC_LIBRARY_64})
target_link_libraries(${PROJECT_NAME} PUBLIC matrixio ${MPI_C_LIBRARIES} ${MKL_LIBRARIES_ILP})
