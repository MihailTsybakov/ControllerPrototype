# http://software.intel.com/en-us/articles/intel-mkl-link-line-advisor
# NEW MKL CMAKE
option(MKL_LINK_STATIC "MKL link only static libraries" ON)

if(WIN32)
  set (MKLROOT_TMP $ENV{MKLROOT})
  if(NOT ${MKLROOT_TMP} STREQUAL "")
    unset(MKLROOT CACHE)
    set (MKLROOT ${MKLROOT_TMP})
  endif(NOT ${MKLROOT_TMP} STREQUAL "")
  if((${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64") OR CMAKE_CL_64)
    find_path(MKL_LIBRARYDIR
      mkl_core.lib
      HINTS
      ${MKLROOT}/lib/intel64
      $ENV{MKLROOT}/lib/intel64
      ${LIB_INSTALL_DIR}
      $ENV{MKL_BINARYDIR}
      $ENV{MKLLIB}
      "C:/Program Files (x86)/Intel/oneAPI/mkl/latest/lib/intel64"
      "C:/Program Files (x86)/Intel/Composer XE 2015/mkl/lib/intel64"
      "C:/Program Files (x86)/Intel/Composer XE 2013 SP1/mkl/lib/intel64"
      "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2017/windows/mkl/lib/intel64"
      "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/mkl/lib/intel64"
    )
    find_path(MKL_BINARYDIR
      mkl_core.dll
      mkl_core.1.dll
      HINTS
      ${MKLROOT}/../redist/intel64/mkl
      $ENV{MKLROOT}/../redist/intel64/mkl
      ${LIB_INSTALL_DIR}
      $ENV{MKL_BINARYDIR}
      $ENV{MKLLIB}
	  "C:/Program Files (x86)/Intel/oneAPI/mkl/latest/redist/intel64"
      "C:/Program Files (x86)/Intel/Composer XE 2015/redist/intel64/mkl"
      "C:/Program Files (x86)/Intel/Composer XE 2013 SP1/redist/intel64/mkl"
      "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2017/windows/redist/intel64/mkl"
      "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/redist/intel64/mkl"
    )
    find_path(OMP_LIBRARYDIR
      libiomp5md.lib
      HINTS
	  "C:/Program Files (x86)/Intel/oneAPI/compiler/latest/windows/compiler/lib/intel64_win"
      "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/compiler/lib/intel64"
    )
  else((${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64") OR CMAKE_CL_64)
    find_path(MKL_LIBRARYDIR
      mkl_core.lib
      HINTS
      ${MKLROOT}/lib/ia32
      $ENV{MKLROOT}/lib/ia32
      ${LIB_INSTALL_DIR}
      $ENV{MKL_BINARYDIR}
      $ENV{MKLLIB}
      "C:/Program Files (x86)/Intel/Composer XE 2015/mkl/lib/ia32"
      "C:/Program Files (x86)/Intel/Composer XE 2013 SP1/mkl/lib/ia32"
      "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2017/windows/mkl/lib/ia32"
      "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/mkl/lib/ia32"
    )
    find_path(MKL_BINARYDIR
      mkl_core.dll
      HINTS
      ${MKLROOT}/../redist/ia32/mkl
      $ENV{MKLROOT}/../redist/ia32/mkl
      ${LIB_INSTALL_DIR}
      $ENV{MKL_BINARYDIR}
      $ENV{MKLLIB}
      "C:/Program Files (x86)/Intel/Composer XE 2015/redist/ia32/mkl"
      "C:/Program Files (x86)/Intel/Composer XE 2013 SP1/redist/ia32/mkl"
      "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2017/windows/mkl"
      "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/mkl"
    )
  endif((${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64") OR CMAKE_CL_64)
else(WIN32)
  if((${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64") OR CMAKE_CL_64) # 64 bit
    EXECUTE_PROCESS(
      COMMAND sh ${CMAKE_SOURCE_DIR}/CMake/IntelMKLLinuxEnv.sh intel64 lp64
      OUTPUT_VARIABLE MKLROOT
    )
  else((${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64") OR CMAKE_CL_64) # 64 bit
    EXECUTE_PROCESS(
      COMMAND sh ${CMAKE_SOURCE_DIR}/CMake/IntelMKLLinuxEnv.sh ia32 lp64
      OUTPUT_VARIABLE MKLROOT
    )
  endif((${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64") OR CMAKE_CL_64) # 64 bit
  STRING(REGEX REPLACE "\n" "" MKLROOT ${MKLROOT})
endif(WIN32)

find_path(MKL_INCLUDES
    mkl.h
  PATHS
    ${MKLROOT}/include
    "C:/Program Files (x86)/Intel/Composer XE 2015/mkl/include"
    "C:/Program Files (x86)/Intel/Composer XE 2013 SP1/mkl/include"
    "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2017/windows/mkl/include"
    "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/mkl/include"
    "C:/Program Files (x86)/Intel/oneAPI/mkl/latest/include"
)

### Discover MKL version
if(MKL_INCLUDES)
  file(STRINGS "${MKL_INCLUDES}/mkl.h" MKL_VERSION_TMP REGEX "^#define __INTEL_MKL_[A-Z_]+ [0-9]+$")
  if (MKL_VERSION_TMP STREQUAL "")
  file(STRINGS "${MKL_INCLUDES}/mkl_version.h" MKL_VERSION_TMP REGEX "^#define __INTEL_MKL_[A-Z_]+ [0-9]+$")
  endif (MKL_VERSION_TMP STREQUAL "")
  if (MKL_VERSION_TMP STREQUAL "")
    set(MKL_VERSION_MAJOR 0)
    message(ERROR 1 "Found Intel Math Kernel Library verion less than 11.1")
  else ()
    string(REGEX REPLACE ".*#define __INTEL_MKL__ ([0-9]+).*" "\\1" MKL_VERSION_MAJOR ${MKL_VERSION_TMP})
    string(REGEX REPLACE ".*#define __INTEL_MKL_MINOR__ ([0-9]+).*" "\\1" MKL_VERSION_MINOR ${MKL_VERSION_TMP})
    string(REGEX REPLACE ".*#define __INTEL_MKL_UPDATE__ ([0-9]+).*" "\\1" MKL_VERSION_UPDATE ${MKL_VERSION_TMP})
    if (${MKL_VERSION_MAJOR} EQUAL 11)
      if (${MKL_VERSION_MINOR} LESS 1)
        message(FATAL_ERROR " Found Intel Math Kernel Library verion less than 11.1")
      endif()
    endif()
    if (${MKL_VERSION_MAJOR} LESS 10)
      message(FATAL_ERROR " Found Intel Math Kernel Library verion less than 11.1")
    endif()
    message(STATUS "Found Intel Math Kernel Library version ${MKL_VERSION_MAJOR}.${MKL_VERSION_MINOR}.${MKL_VERSION_UPDATE}")
  endif ()
else(MKL_INCLUDES)
  message(FATAL_ERROR "Intel MKL libraries not found!")
endif(MKL_INCLUDES)

if(WIN32)
  if((${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64") OR CMAKE_CL_64)
    if(MKL_LINK_STATIC)
      set(MKL_LIBRARIES_LP
        mkl_scalapack_lp64.lib
        mkl_cdft_core.lib
        mkl_intel_lp64.lib
        mkl_core.lib
        mkl_intel_thread.lib
        mkl_blacs_intelmpi_lp64.lib
        libiomp5md.lib
      )
      set(MKL_LIBRARIES_ILP
        mkl_scalapack_ilp64.lib
        mkl_cdft_core.lib
        mkl_intel_ilp64.lib
        mkl_core.lib
        mkl_intel_thread.lib
        mkl_blacs_intelmpi_ilp64.lib
        libiomp5md.lib
      )
    else(MKL_LINK_STATIC)
      set(MKL_LIBRARIES_LP
        mkl_blacs_intelmpi_lp64.lib
        mkl_blas95_lp64.lib
        mkl_lapack95_lp64.lib
        mkl_scalapack_lp64_dll.lib
        mkl_cdft_core_dll.lib
        mkl_intel_lp64_dll.lib
        mkl_intel_thread_dll.lib
        mkl_core_dll.lib           
      )
      set(MKL_LIBRARIES_ILP
        mkl_blas95_ilp64.lib
        mkl_blacs_intelmpi_ilp64.lib
        mkl_lapack95_ilp64.lib
        mkl_scalapack_ilp64_dll.lib       
        mkl_intel_ilp64_dll.lib
        mkl_intel_thread_dll.lib
        mkl_cdft_core_dll.lib
        mkl_core_dll.lib              
      )
    endif(MKL_LINK_STATIC)
  else((${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64") OR CMAKE_CL_64)
    set(FIDESYS_MKL_USE_64_BIT OFF)
    if(MKL_LINK_STATIC)
      set(MKL_LIBRARIES_LP 
        libiomp5md.lib
        mkl_blas95.lib
        mkl_blacs_intelmpi.lib
        mkl_lapack95.lib
        mkl_scalapack_core.lib       
        mkl_intel_c.lib
        mkl_intel_thread.lib
        mkl_cdft_core.lib
        mkl_core.lib  
      )
    else(MKL_LINK_STATIC)
      set(MKL_LIBRARIES_LP         
        mkl_blas95.lib
        mkl_blacs_intelmpi.lib
        mkl_lapack95.lib
        mkl_scalapack_core.lib 
        mkl_intel_c.lib
        mkl_intel_thread.lib
        mkl_cdft_core.lib
        mkl_core.lib      
      )
    endif(MKL_LINK_STATIC)
    set(MKL_LIBRARIES_ILP ${MKL_LIBRARIES_LP})
  endif((${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64") OR CMAKE_CL_64)
else(WIN32)
  if((${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64") OR CMAKE_CL_64)
    set(MKL_LIBRARIES_LP 
      libmkl_blas95_lp64.a
      libmkl_lapack95_lp64.a
      libmkl_scalapack_lp64.a
      -Wl,--start-group
      libmkl_cdft_core.a   
      libmkl_intel_lp64.a
      libmkl_intel_thread.a
      libmkl_core.a   
      libmkl_blacs_intelmpi_lp64.a
      -Wl,--end-group
      iomp5
      pthread
      m
      dl
    )
    set(MKL_LIBRARIES_ILP
      libmkl_blas95_ilp64.a
      libmkl_lapack95_ilp64.a
      libmkl_scalapack_ilp64.a
      -Wl,--start-group
      libmkl_cdft_core.a   
      libmkl_intel_ilp64.a
      libmkl_intel_thread.a
      libmkl_core.a
      libmkl_blacs_intelmpi_ilp64.a
      -Wl,--end-group
      iomp5
      pthread
      m
      dl
    )
    set (MKL_LIBRARYDIR ${MKLROOT}/lib/intel64)
  else((${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64") OR CMAKE_CL_64)
    set(MKL_LIBRARIES_LP            
      libmkl_blas95.a
      libmkl_lapack95.a
      libmkl_scalapack_core.a
      libmkl_cdft_core.a
      libmkl_intel.a
      libmkl_intel_thread.a
      libmkl_core.a
      libmkl_blacs_intelmpi.a
      iomp5
      pthread
      m
      dl
    )
    set(MKL_LIBRARIES_ILP ${MKL_LIBRARIES_LP})
    set(FIDESYS_MKL_USE_64_BIT OFF)
    set (MKL_LIBRARYDIR ${MKLROOT}/lib/ia32)
  endif((${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64") OR CMAKE_CL_64)
  set (MKL_BINARYDIR ${MKL_LIBRARYDIR})
endif(WIN32)

if (FIDESYS_MKL_USE_64_BIT)
  add_definitions(-DMKL_ILP64)
else (FIDESYS_MKL_USE_64_BIT)
  set(MKL_LIBRARIES_ILP ${MKL_LIBRARIES_LP})
endif (FIDESYS_MKL_USE_64_BIT)

message(STATUS "Found MKL:")
message(STATUS "  MKL_ROOT: "  $ENV{MKL_ROOT})
message(STATUS "  MKL_INCLUDES: "  ${MKL_INCLUDES})
message(STATUS "  MKL_LIBRARIES_LP: "  ${MKL_LIBRARIES_LP})
message(STATUS "  MKL_LIBRARIES_ILP: " ${MKL_LIBRARIES_ILP})
message(STATUS "  MKL_LIBRARYDIR: " ${MKL_LIBRARYDIR})
message(STATUS "  MKL_BINARYDIR: " ${MKL_BINARYDIR})
message(STATUS "  OMP_LIBRARYDIR: " ${OMP_LIBRARYDIR})

unset(MKL_LIBRARIES CACHE)

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(MKL DEFAULT_MSG MKL_INCLUDES MKL_LIBRARYDIR MKL_BINARYDIR OMP_LIBRARYDIR MKL_LIBRARIES_LP MKL_LIBRARIES_ILP)

mark_as_advanced(MKL_INCLUDES MKL_LIBRARIES_LP MKL_LIBRARIES_ILP MKL_LIBRARYDIR MKL_BINARYDIR OMP_LIBRARYDIR)
