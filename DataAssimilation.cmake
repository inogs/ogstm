# CMake project file for OGSTM
# original authors : E.Pascolo, S.Bna, L.Calori
# edited by: S.Spada (sspada@ogs.it)

# CMAKE setting
cmake_minimum_required (VERSION 2.6)
project (OGSTM)
set (CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/)
enable_language(Fortran C)
set (NETCDF_F90 "YES")
find_package(MPI REQUIRED)
find_package(NetCDF REQUIRED)
find_package(BFM REQUIRED)

# DATA ASSIMILATION PACKAGES
find_package(3DVAR REQUIRED)
find_package(PETSc REQUIRED)
find_package(PnetCDF REQUIRED)

# NEW EnsDA PACKAGES
set(BLA_VENDOR Intel10_64lp_seq)
find_package(LAPACK REQUIRED)

if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

get_filename_component (default_prefix "./bin" ABSOLUTE)
set (CMAKE_INSTALL_PREFIX ${default_prefix} CACHE STRING
      "Choose the installation directory; by default it installs in the OGSTM directory."
      FORCE)

# COMPILER

get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
add_definitions(-Dkey_trahdfcoef1d -Dkey_trahdfbilap -Dkey_trc_smolar)
add_definitions(-Dkey_trc_hdfbilap -Dkey_trc_dmp -Dkey_kef -Dkey_trc_sed )
add_definitions(-Dkey_mpp -Dkey_mpp_mpi)
IF (BFMv2)
    add_definitions(-DBFMv2)
ENDIF()

add_definitions(-DExecEns)
add_definitions(-DExecEnsDA)
add_definitions(-DExecEnsParams)

if (MPI_Fortran_COMPILER MATCHES "mpiifort.*")
  # mpiifort
  set (CMAKE_Fortran_FLAGS_RELEASE " -fno-math-errno -Ofast -ipo -xHost -qopt-report5 -fpp -align array64byte")
  set (CMAKE_Fortran_FLAGS_DEBUG   " -O0 -g -fpp -CB -fp-stack-check -check all -traceback -gen-interfaces -warn interfaces -extend_source") #-fpe0 removed due to dsyevr needing ieee exceptions
elseif (MPI_Fortran_COMPILER MATCHES "mpif90.*")
  # mpif90
  set (CMAKE_Fortran_FLAGS_RELEASE " -O2  -fimplicit-none -cpp -ffree-line-length-0")
  set (CMAKE_Fortran_FLAGS_DEBUG   " -O0 -g -Wall -Wextra -cpp -fbounds-check -fimplicit-none -ffpe-trap=invalid -pedantic -ffree-line-length-0 ") #-ffpe-trap=overflow removed because it may interfere with dsyevr
else ()
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif ()

#include
include_directories(${BFM_INCLUDES})
include_directories(${NETCDF_INCLUDES_C})
include_directories(${NETCDFF_INCLUDES_F90})

# DATA ASSIMILATION INCLUDE SECTION
include_directories(${DA_INCLUDES})
include_directories(${PETSC_INCLUDES})

# Search Fortran module to compile
file(GLOB_RECURSE TMP src/*)
list (APPEND FORTRAN_SOURCES ${TMP})

#building
add_library( ogstm_lib ${FORTRAN_SOURCES})
add_executable (ogstm.xx application/ogstm_main_caller.f90)
target_link_libraries( ogstm.xx ogstm_lib ${NETCDFF_LIBRARIES_F90} ${BFM_LIBRARIES} ${DA_LIBRARIES} ${PETSC_LIBRARIES} ${PNETCDF_LIBRARIES} ${LAPACK_LIBRARIES})
