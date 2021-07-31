if ((NOT PETSC_DIR) OR ("${PETSC_DIR}" STREQUAL ""))
set(PETSC_DIR $ENV{PETSC_DIR} CACHE TYPE STRING)
endif ()

# No testing!!!
#unset(PETSC_WITH_CUDA CACHE)
#set(PETSC_WITH_CUDA OFF)
option(PETSC_WITH_CUDA
    "PETSc compiled with CUDA" OFF)

if (WIN32)
option(PETSC_WITH_INTEL_COMPILER
    "PETSc compiled with the intel compiler" OFF)
else (WIN32)
option(PETSC_WITH_INTEL_COMPILER
    "PETSc compiled with the intel compiler" ON)
endif (WIN32)

if (("${CMAKE_CXX_COMPILER}" STREQUAL "icl") OR FIDESYS_USE_INTEL_COMPILER)
  set(PETSC_WITH_INTEL_COMPILER ON)
endif (("${CMAKE_CXX_COMPILER}" STREQUAL "icl") OR FIDESYS_USE_INTEL_COMPILER)

option(PETSC_DEBUG
    "PETSc compiled with debug info" OFF)

set(PETSC_ARCH_CMAKE "compiler_")

if(PETSC_WITH_INTEL_COMPILER)
	set(PETSC_ARCH_CMAKE ${PETSC_ARCH_CMAKE}intel-)
else(PETSC_WITH_INTEL_COMPILER)
	if (WIN32)
		set(PETSC_ARCH_CMAKE ${PETSC_ARCH_CMAKE}ms-)
	else (WIN32)
		set(PETSC_ARCH_CMAKE ${PETSC_ARCH_CMAKE}gnu-)
	endif (WIN32)
endif(PETSC_WITH_INTEL_COMPILER)

if (WIN32)
	string(FIND ${CMAKE_GENERATOR} "Visual Studio 11" is2012)
	if(${is2012} MATCHES 0)
		set(PETSC_ARCH_CMAKE ${PETSC_ARCH_CMAKE}msvc_2012-)
	endif(${is2012} MATCHES 0)
	string(FIND ${CMAKE_GENERATOR} "Visual Studio 12" is2013)
	if(${is2013} MATCHES 0)
		set(PETSC_ARCH_CMAKE ${PETSC_ARCH_CMAKE}msvc_2013-)
	endif(${is2013} MATCHES 0)
	string(FIND ${CMAKE_GENERATOR} "Visual Studio 14" is2015)
	if(${is2015} MATCHES 0)
		set(PETSC_ARCH_CMAKE ${PETSC_ARCH_CMAKE}msvc_2015-)
	endif(${is2015} MATCHES 0)
endif (WIN32)

if (MPI_Implementation STREQUAL "mpich")
	set(PETSC_ARCH_CMAKE ${PETSC_ARCH_CMAKE}mpich2)
endif (MPI_Implementation STREQUAL "mpich")
if (MPI_Implementation STREQUAL "intel")
	set(PETSC_ARCH_CMAKE ${PETSC_ARCH_CMAKE}intelmpi)
endif (MPI_Implementation STREQUAL "intel")
if (MPI_Implementation STREQUAL "openmpi")
	set(PETSC_ARCH_CMAKE ${PETSC_ARCH_CMAKE}openmpi)
endif (MPI_Implementation STREQUAL "openmpi")

if (PETSC_WITH_CUDA)
	set(PETSC_ARCH_CMAKE ${PETSC_ARCH_CMAKE}-cuda)
else (PETSC_WITH_CUDA)
	set(PETSC_ARCH_CMAKE ${PETSC_ARCH_CMAKE}-cudano)
endif (PETSC_WITH_CUDA)

if(FIDESYS_SOLVER_USE_MUMPS)
	set(PETSC_ARCH_CMAKE32 ${PETSC_ARCH_CMAKE}-mumps)
else(FIDESYS_SOLVER_USE_MUMPS)
	set(PETSC_ARCH_CMAKE32 ${PETSC_ARCH_CMAKE}-nomumps)
endif(FIDESYS_SOLVER_USE_MUMPS)

set(PETSC_ARCH_CMAKE64 ${PETSC_ARCH_CMAKE}-nomumps)

if((${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64") OR CMAKE_CL_64) # 64 bit
    set(PETSC_ARCH_CMAKE32 ${PETSC_ARCH_CMAKE32}-indexes_32)
    set(PETSC_ARCH_CMAKE64 ${PETSC_ARCH_CMAKE64}-indexes_64)
else((${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64") OR CMAKE_CL_64) # 32 bit
    set(PETSC_ARCH_CMAKE32 ${PETSC_ARCH_CMAKE32})
    set(PETSC_ARCH_CMAKE64 ${PETSC_ARCH_CMAKE32})
endif((${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64") OR CMAKE_CL_64)

unset(PETSC_ARCH_CMAKE)

if((${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64") OR CMAKE_CL_64) # 64 bit
    set(PETSC_ARCH_CMAKE32 ${PETSC_ARCH_CMAKE32}-mkl_lp64)
    set(PETSC_ARCH_CMAKE64 ${PETSC_ARCH_CMAKE64}-mkl_ilp64)
endif((${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64") OR CMAKE_CL_64) # 64 bit

if (PETSC_DEBUG)
    set(PETSC_ARCH_CMAKE32 ${PETSC_ARCH_CMAKE32}-debug)
    set(PETSC_ARCH_CMAKE64 ${PETSC_ARCH_CMAKE64}-debug)
endif (PETSC_DEBUG)

if(EXISTS ${PETSC_DIR}/${PETSC_ARCH_CMAKE32}/include)
unset(PETSC_CONF_INCLUDES_32 CACHE)
unset(PETSC_CONF_INCLUDES_64 CACHE)
unset(PETSC_LIBRARY_32 CACHE)
unset(PETSC_LIBRARY_64 CACHE)
endif(EXISTS ${PETSC_DIR}/${PETSC_ARCH_CMAKE32}/include)

find_path(PETSC_INCLUDES
    NAMES
    petsc.h
    PATHS
    ${PETSC_DIR}/include
)

### Discover PETSc version
if(PETSC_INCLUDES)
  file(READ  "${PETSC_INCLUDES}/petscversion.h" PETSC_VERSION_H_CONTENTS)
  string (REGEX MATCH "PETSC_VERSION_MAJOR[ \t]+([0-9]+)"    PETSC_VERSION_MAJOR    ${PETSC_VERSION_H_CONTENTS})
  string (REGEX MATCH "PETSC_VERSION_MINOR[ \t]+([0-9]+)"    PETSC_VERSION_MINOR    ${PETSC_VERSION_H_CONTENTS})
  string (REGEX MATCH "PETSC_VERSION_SUBMINOR[ \t]+([0-9]+)" PETSC_VERSION_SUBMINOR ${PETSC_VERSION_H_CONTENTS})
  string (REGEX MATCH "([0-9]+)" PETSC_VERSION_MAJOR    ${PETSC_VERSION_MAJOR})
  string (REGEX MATCH "([0-9]+)" PETSC_VERSION_MINOR    ${PETSC_VERSION_MINOR})
  string (REGEX MATCH "([0-9]+)" PETSC_VERSION_SUBMINOR ${PETSC_VERSION_SUBMINOR})
  message(STATUS "Found PETSc: " ${PETSC_VERSION_MAJOR}.${PETSC_VERSION_MINOR}.${PETSC_VERSION_SUBMINOR}) 

  if (PETSC_VERSION_MAJOR STREQUAL "")
    set(PETSc_VERSION 0)
    message(FATAL_ERROR "Found PETSc Library verion less than 3.4.3")
  else ()
    if (${PETSC_VERSION_MAJOR} EQUAL 3)
      if (${PETSC_VERSION_MINOR} LESS 4)
        if (${PETSC_VERSION_MINOR} LESS 3)
          message(FATAL_ERROR "Found PETSc Library verion less than 3.4.3")
        endif()
      endif()
    endif()
    if (${PETSC_VERSION_MAJOR} LESS 2)
      message(FATAL_ERROR "Found PETSc Library verion less than 3.4.3")
    endif()
  endif ()
else(PETSC_INCLUDES)
  message(FATAL_ERROR "PETSc library not found!")
endif(PETSC_INCLUDES)

find_path(PETSC_CONF_INCLUDES_32
    NAMES
    petscconf.h
    PATHS
    ${PETSC_DIR}/include
    ${PETSC_DIR}/${PETSC_ARCH_CMAKE32}/include
)

find_path(PETSC_CONF_INCLUDES_64
    NAMES
    petscconf.h
    PATHS
    ${PETSC_DIR}/include
    ${PETSC_DIR}/${PETSC_ARCH_CMAKE64}/include
)

if (WIN32)
    find_file(PETSC_LIBRARY_32
        libpetsc.lib
        PATHS
        ${PETSC_DIR}/${PETSC_ARCH_CMAKE32}/lib
    )
    find_file(PETSC_LIBRARY_64
        libpetsc.lib
        PATHS
        ${PETSC_DIR}/${PETSC_ARCH_CMAKE64}/lib
    )
else (WIN32)
    find_file(PETSC_LIBRARY_32
        libpetsc.so
        PATHS
        ${PETSC_DIR}/lib
        ${PETSC_DIR}/${PETSC_ARCH_CMAKE32}/lib
    )
    find_file(PETSC_LIBRARY_64
        libpetsc.so
        PATHS
        ${PETSC_DIR}/lib
        ${PETSC_DIR}/${PETSC_ARCH_CMAKE64}/lib
    )
endif (WIN32)

if (PETSC_WITH_CUDA)
    find_package(CUDA 5.0 REQUIRED)
    cuda_include_directories($ENV{CUSP_DIR})
    cuda_include_directories($ENV{THRUST_DIR})
    get_filename_component(CUDA_BIN_DIR ${CUDA_NVCC_EXECUTABLE} PATH CACHE)
    add_definitions(-DPETSC_WITH_CUDA)
endif (PETSC_WITH_CUDA)

include(FindPackageHandleStandardArgs)

message(STATUS "  PETSc directory: " ${PETSC_DIR})
message(STATUS "  PETSC_ARCH 32: " ${PETSC_ARCH_CMAKE32})
message(STATUS "  PETSC_ARCH 64: " ${PETSC_ARCH_CMAKE64})

find_package_handle_standard_args(PETSC DEFAULT_MSG PETSC_INCLUDES PETSC_CONF_INCLUDES_32 PETSC_LIBRARY_32 PETSC_ARCH_CMAKE32 PETSC_CONF_INCLUDES_64 PETSC_LIBRARY_64 PETSC_ARCH_CMAKE64)
