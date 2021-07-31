# NEW MPI CMAKE
if (WIN32)
if((${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64") OR CMAKE_CL_64)
	find_path(MPI_INTEL_INCLUDES
		NAMES
		mpi.h
		PATHS
		$ENV{I_MPI_ROOT}/intel64/include
		$ENV{PROGRAMFILES}/intel/MPI/*/em64t/include
		$ENV{I_MPI_ROOT}em64t/include
		$ENV{I_MPI_ROOT}/em64t/include
		$ENV{I_MPI_ONEAPI_ROOT}/include
	)
	find_path(MPI_BINARY_PATH
		NAMES
		impi.dll
		PATHS
		$ENV{I_MPI_ROOT}/intel64/bin
		$ENV{PROGRAMFILES}/intel/MPI/*/em64t/bin
		$ENV{I_MPI_ROOT}em64t/bin
		$ENV{I_MPI_ROOT}/em64t/bin
	)
	find_file(MPI_INTEL_LIB_C
		impi.lib
		PATHS
		$ENV{I_MPI_ROOT}/intel64/lib
		$ENV{PROGRAMFILES}/intel/MPI/*/em64t/lib
		$ENV{I_MPI_ROOT}/em64t/lib
		$ENV{I_MPI_ROOT}em64t/lib
		$ENV{I_MPI_ONEAPI_ROOT}/lib/release
	)
	find_file(MPI_INTEL_LIB_CXX
		impicxx.lib
		PATHS
		$ENV{I_MPI_ROOT}/intel64/lib
		$ENV{PROGRAMFILES}/intel/MPI/*/em64t/lib
		$ENV{I_MPI_ROOT}em64t/lib
		$ENV{I_MPI_ROOT}/em64t/lib
		$ENV{I_MPI_ONEAPI_ROOT}/lib/release
	)
	find_file(MPI_INTEL_LIB_Fortran
		impi.lib
		PATHS
		$ENV{I_MPI_ROOT}/intel64/lib
		$ENV{PROGRAMFILES}/intel/MPI/*/em64t/lib
		$ENV{I_MPI_ROOT}em64t/lib
		$ENV{I_MPI_ROOT}/em64t/lib
		$ENV{I_MPI_ONEAPI_ROOT}/lib/release
	)
	find_path(MPI_MPICH_INCLUDES
		NAMES
		mpi.h
		PATHS
		[HKEY_LOCAL_MACHINE\\SOFTWARE\\MPICH\\SMPD;binary]/../include
		[HKEY_LOCAL_MACHINE\\SOFTWARE\\MPICH2;Path]/include
		$ENV{ProgramW6432}/mpich2/include
	)
	find_file(MPI_MPICH_LIB_C
		mpi.lib
		PATHS
		$ENV{ProgramW6432}/mpich2/lib
		[HKEY_LOCAL_MACHINE\\SOFTWARE\\MPICH\\SMPD;binary]/../lib
		[HKEY_LOCAL_MACHINE\\SOFTWARE\\MPICH2;Path]/lib
	)
	find_file(MPI_MPICH_LIB_CXX
		cxx.lib
		PATHS
		$ENV{ProgramW6432}/mpich2/lib
		[HKEY_LOCAL_MACHINE\\SOFTWARE\\MPICH\\SMPD;binary]/../lib
		[HKEY_LOCAL_MACHINE\\SOFTWARE\\MPICH2;Path]/lib
	)
	find_file(MPI_MPICH_LIB_Fortran
		fmpich2.lib
		PATHS
		$ENV{ProgramW6432}/mpich2/lib
		[HKEY_LOCAL_MACHINE\\SOFTWARE\\MPICH\\SMPD;binary]/../lib
		[HKEY_LOCAL_MACHINE\\SOFTWARE\\MPICH2;Path]/lib
	)
else((${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64") OR CMAKE_CL_64)
	find_path(MPI_INTEL_INCLUDES
		NAMES
		mpi.h
		PATHS
		$ENV{PROGRAMFILES}/intel/MPI/*/ia32/include
		$ENV{I_MPI_ROOT}ia32/include
	)
	find_path(MPI_MPICH_INCLUDES
		NAMES
		mpi.h
		PATHS
		"$ENV{ProgramFiles(x86)}/mpich2/include"
		$ENV{ProgramFiles}/mpich2/include
		[HKEY_LOCAL_MACHINE\\SOFTWARE\\MPICH\\SMPD;binary]/../include
		[HKEY_LOCAL_MACHINE\\SOFTWARE\\MPICH2;Path]/include
	)
	find_file(MPI_INTEL_LIB_C
		impi.lib
		PATHS
		$ENV{PROGRAMFILES}/intel/MPI/*/ia32/lib
		$ENV{I_MPI_ROOT}ia32/lib
	)
	find_file(MPI_INTEL_LIB_CXX
		impicxx.lib
		PATHS
		$ENV{PROGRAMFILES}/intel/MPI/*/ia32/lib
		$ENV{I_MPI_ROOT}ia32/lib
	)
	find_file(MPI_INTEL_LIB_Fortran
		impi.lib
		PATHS
		$ENV{PROGRAMFILES}/intel/MPI/*/ia32/lib
		$ENV{I_MPI_ROOT}ia32/lib
	)
	find_file(MPI_MPICH_LIB_C
		mpi.lib
		PATHS
		"$ENV{ProgramFiles(x86)}/mpich2/lib"
		$ENV{ProgramFiles}/mpich2/lib
		[HKEY_LOCAL_MACHINE\\SOFTWARE\\MPICH\\SMPD;binary]/../lib
		[HKEY_LOCAL_MACHINE\\SOFTWARE\\MPICH2;Path]/lib
	)
	find_file(MPI_MPICH_LIB_CXX
		cxx.lib
		PATHS
		"$ENV{ProgramFiles(x86)}/mpich2/lib"
		$ENV{ProgramFiles}/mpich2/lib
		[HKEY_LOCAL_MACHINE\\SOFTWARE\\MPICH\\SMPD;binary]/../lib
		[HKEY_LOCAL_MACHINE\\SOFTWARE\\MPICH2;Path]/lib
	)
	find_file(MPI_MPICH_LIB_Fortran
		fmpich2.lib
		PATHS
		"$ENV{ProgramFiles(x86)}/mpich2/lib"
		$ENV{ProgramFiles}/mpich2/lib
		[HKEY_LOCAL_MACHINE\\SOFTWARE\\MPICH\\SMPD;binary]/../lib
		[HKEY_LOCAL_MACHINE\\SOFTWARE\\MPICH2;Path]/lib
	)
endif((${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64") OR CMAKE_CL_64)

else (WIN32)
if((${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64") OR CMAKE_CL_64)
	find_path(MPI_INTEL_INCLUDES
		NAMES
		mpi.h
		PATHS
		$ENV{I_MPI_ROOT}/intel64/include
	)
	find_file(MPI_INTEL_LIB_C
		libmpi.so
		PATHS
		$ENV{I_MPI_ROOT}/intel64/lib/release
	)
	find_file(MPI_INTEL_LIB_CXX
		libmpicxx.so
		PATHS
		$ENV{I_MPI_ROOT}/intel64/lib
	)
	find_file(MPI_INTEL_LIB_Fortran
		libmpifort.so
		PATHS
		$ENV{I_MPI_ROOT}/intel64/lib
	)

else((${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64") OR CMAKE_CL_64)
	find_path(MPI_INTEL_INCLUDES
		NAMES
		mpi.h
		PATHS
		$ENV{I_MPI_ROOT}/include32
	)
	find_file(MPI_INTEL_LIB_C
		libmpi.so
		PATHS
		$ENV{I_MPI_ROOT}/lib32
	)
	find_file(MPI_INTEL_LIB_CXX
		libmpi.so
		PATHS
		$ENV{I_MPI_ROOT}/lib32
	)
	find_file(MPI_INTEL_LIB_Fortran
		libmpigf.so
		PATHS
		$ENV{I_MPI_ROOT}/lib32
	)

endif((${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64") OR CMAKE_CL_64)

endif (WIN32)

if (MPI_INTEL_INCLUDES AND MPI_INTEL_LIB_C AND MPI_INTEL_LIB_CXX AND MPI_INTEL_LIB_Fortran)
	set(MPI_C_FOUND TRUE)
	set(MPI_CXX_FOUND TRUE)
	set(MPI_Fortran_FOUND TRUE)
	set(MPI_Implementation "intel" CACHE STRING "Intel MPI")
	set(MPI_Values "intel")
endif (MPI_INTEL_INCLUDES AND MPI_INTEL_LIB_C AND MPI_INTEL_LIB_CXX AND MPI_INTEL_LIB_Fortran)

if (MPI_MPICH_INCLUDES AND MPI_MPICH_LIB_C AND MPI_MPICH_LIB_CXX AND MPI_MPICH_LIB_Fortran)
	set(MPI_C_FOUND TRUE)
	set(MPI_CXX_FOUND TRUE)
	set(MPI_Fortran_FOUND TRUE)
	set(MPI_Implementation "mpich" CACHE STRING  "MPICH MPI")
	if (MPI_Default_LIBRARIES)
		set(MPI_Values ${MPI_Values} ";")
	endif (MPI_Default_LIBRARIES)
	set(MPI_Values ${MPI_Values} "mpich2")
endif (MPI_MPICH_INCLUDES AND MPI_MPICH_LIB_C AND MPI_MPICH_LIB_CXX AND MPI_MPICH_LIB_Fortran)

if (MPI_OPENMPI_INCLUDES AND MPI_OPENMPI_LIB_C AND MPI_OPENMPI_LIB_CXX AND MPI_OPENMPI_LIB_Fortran)
	set(MPI_C_FOUND TRUE)
	set(MPI_CXX_FOUND TRUE)
	set(MPI_Fortran_FOUND TRUE)
	set(MPI_Implementation "openmpi" CACHE STRING  "OPENMPI MPI")
	if (MPI_Default_LIBRARIES)
		set(MPI_Values ${MPI_Values} ";")
	endif (MPI_Default_LIBRARIES)
	set(MPI_Values ${MPI_Values} "openmpi")
endif (MPI_OPENMPI_INCLUDES AND MPI_OPENMPI_LIB_C AND MPI_OPENMPI_LIB_CXX AND MPI_OPENMPI_LIB_Fortran)

set_property(CACHE MPI_Implementation PROPERTY STRINGS ${MPI_Values})

if (MPI_Implementation STREQUAL "intel")
	set(MPI_BINARY_PATH ${MPI_INTEL_BINARY_PATH})
	set(MPI_INCLUDE_PATH ${MPI_INTEL_INCLUDES})
	set(MPI_CXX_INCLUDE_PATH ${MPI_INTEL_INCLUDES})
	set(MPI_Fortran_INCLUDE_PATH ${MPI_INTEL_INCLUDES})
	set(MPI_C_LIBRARIES ${MPI_INTEL_LIB_C})
	set(MPI_CXX_LIBRARIES ${MPI_INTEL_LIB_CXX})
	set(MPI_Fortran_LIBRARIES ${MPI_INTEL_LIB_Fortran})
endif (MPI_Implementation STREQUAL "intel")

if (MPI_Implementation STREQUAL "mpich2")
	set(MPI_INCLUDE_PATH ${MPI_MPICH_INCLUDES})
	set(MPI_Fortran_INCLUDE_PATH ${MPI_MPICH_INCLUDES})
	set(MPI_CXX_INCLUDE_PATH ${MPI_MPICH_INCLUDES})
	set(MPI_C_LIBRARIES ${MPI_MPICH_LIB_C})
	set(MPI_CXX_LIBRARIES ${MPI_MPICH_LIB_CXX})
	set(MPI_Fortran_LIBRARIES ${MPI_MPICH_LIB_Fortran})
endif (MPI_Implementation STREQUAL "mpich2")

if (MPI_Implementation STREQUAL "openmpi")
	set(MPI_INCLUDE_PATH ${MPI_OPENMPI_INCLUDES})
	set(MPI_CXX_INCLUDE_PATH ${MPI_OPENMPI_INCLUDES})
	set(MPI_Fortran_INCLUDE_PATH ${MPI_OPENMPI_INCLUDES})
	set(MPI_C_LIBRARIES ${MPI_OPENMPI_LIB_C})
	set(MPI_CXX_LIBRARIES ${MPI_OPENMPI_LIB_CXX})
	set(MPI_Fortran_LIBRARIES ${MPI_OPENMPI_LIB_Fortran})
endif (MPI_Implementation STREQUAL "openmpi")

mark_as_advanced(MPI_MPICH_INCLUDES MPI_MPICH_LIB_C MPI_MPICH_LIB_CXX MPI_MPICH_LIB_Fortran)
mark_as_advanced(MPI_INTEL_INCLUDES MPI_INTEL_LIB_C MPI_INTEL_LIB_CXX MPI_INTEL_LIB_Fortran)

message(STATUS "Found MPI:")
message(STATUS "  I_MPI_ROOT: " $ENV{I_MPI_ROOT}) 
message(STATUS "  MPI_Implementation: " ${MPI_Implementation})
message(STATUS "  MPI_INCLUDE_PATH: " ${MPI_INCLUDE_PATH})
message(STATUS "  MPI_CXX_INCLUDE_PATH: " ${MPI_CXX_INCLUDE_PATH})
message(STATUS "  MPI_C_LIBRARIES: " ${MPI_C_LIBRARIES})
message(STATUS "  MPI_CXX_LIBRARIES: " ${MPI_CXX_LIBRARIES})
message(STATUS "  MPI_Fortran_LIBRARIES: " ${MPI_Fortran_LIBRARIES})


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MPILibraries DEFAULT_MSG MPI_INCLUDE_PATH MPI_CXX_INCLUDE_PATH MPI_C_LIBRARIES MPI_CXX_LIBRARIES MPI_Fortran_LIBRARIES )
