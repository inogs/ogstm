# CMake project file for OGSTM
# author : E.Pascolo, S.Bna, L.Calori

# CMAKE setting
cmake_minimum_required (VERSION 3.18)
project (OGSTM)
set (CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/)
enable_language(Fortran C)
set (NETCDF_F90 "YES")
find_package(MPI REQUIRED)
find_package(NetCDF REQUIRED)
find_package(BFM REQUIRED)

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

if (MPI_Fortran_COMPILER MATCHES "mpiifort.*")
  # mpiifort
  set (CMAKE_Fortran_FLAGS_RELEASE " -fno-math-errno -O2 -xAVX -qopt-report5 -g -cpp -align array64byte") #-qopenmp
  set (CMAKE_Fortran_FLAGS_DEBUG   " -O0 -g -cpp -CB -fp-stack-check -check all -traceback -gen-interfaces -warn interfaces -fpe0 -extend_source") #-qopenmp
elseif (MPI_Fortran_COMPILER MATCHES "mpif90.*")
  # mpif90
#  set (CMAKE_Fortran_FLAGS_RELEASE " -O2  -fimplicit-none -cpp -ffree-line-length-0 ")
#  set (CMAKE_Fortran_FLAGS_DEBUG   " -O0 -g -Wall -Wextra -cpp -fbounds-check -fimplicit-none -ffpe-trap=invalid,overflow -pedantic")
  set (CMAKE_Fortran_FLAGS_RELEASE " -O2  -cpp -g -acc=gpu -gpu=cc70 ")
  set (CMAKE_Fortran_FLAGS_DEBUG   " -O0 -g -Wall -Wextra -cpp -fbounds-check -fimplicit-none -ffpe-trap=invalid,overflow -pedantic")
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

# Search Fortran module to compile
set( FOLDERS BIO  General  IO  MPI  namelists  PHYS BC)
  foreach(FOLDER ${FOLDERS})
  file(GLOB TMP src/${FOLDER}/*)
  list (APPEND FORTRAN_SOURCES ${TMP})
endforeach()

#building
add_library( ogstm_lib ${FORTRAN_SOURCES})
add_executable (ogstm.xx application/ogstm_main_caller.f90)
target_link_libraries( ogstm.xx ogstm_lib ${NETCDFF_LIBRARIES_F90} ${BFM_LIBRARIES})
