# CMake project file for OGSTM
# author : E.Pascolo, S.Bna, L.Calori

cmake_minimum_required (VERSION 2.6)
project (OGSTM)
set (CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/)
enable_language(Fortran)
find_package(MPI REQUIRED)
find_package(NetCDF REQUIRED)

message("Compiler :${MPI_Fortran_COMPILER}")
# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# default installation
get_filename_component (default_prefix "./bin" ABSOLUTE)
set (CMAKE_INSTALL_PREFIX ${default_prefix} CACHE STRING
      "Choose the installation directory; by default it installs in the OGSTM directory."
      FORCE)

# FFLAGS depend on the compiler

get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
if (MPI_Fortran_COMPILER MATCHES "mpiifort.*")
  # mpiifort
  set (CMAKE_Fortran_FLAGS_RELEASE "-fno-math-errno -O2 -g -cpp -unroll=3 -qopt-subscript-in-range -align all -heap-arrays")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -traceback -fp-stack-check -check bounds -fpe0")
elseif (MPI_Fortran_COMPILER MATCHES "mpif90.*")
  # mpif90
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2  -fimplicit-none -cpp  -ffixed-line-length-132")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -Wall -Wextra -fbounds-check -fimplicit-none -ffpe-trap=invalid,overflow -pedantic -ffixed-line-length-132")
else ()
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif ()

#include
include_directories(${CMAKE_SOURCE_DIR}/../bfm/include/)
include_directories(${NETCDF_INCLUDE_DIR})

# Fortran module
set ( FOLDERS BIO  General  IO  MPI  namelists  PHYS)
foreach (FOLDER ${FOLDERS})
file(GLOB TMP src/${FOLDER}/*)
#message(" tmp ->${TMP}<--" )
list (APPEND FORTRAN_SOURCES ${TMP})
endforeach()
message(" pippo ->${FORTRAN_SOURCES}<--" )

#set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)

add_library( ogstm_lib ${FORTRAN_SOURCES})
add_executable (ogstm.xx application/ogstm_main_caller.f90)

target_link_libraries( ogstm.xx ogstm_lib ${NETCDF_LIBRARY})
#set (EXECUTABLES "ogstm.xx" ${NMPROGRAMS})
#set (SCRIPTS "gen_pert.sh" "pert_multi_mode.sh")

#add_executable ("NORMA.exe" "NORMA.f")
#foreach (p ${NMPROGRAMS})
#  add_executable (${p} "${p}.f")
#endforeach (p)

# install executables and scripts
#install (TARGETS ${EXECUTABLES}
#         RUNTIME DESTINATION "bin")
#install (PROGRAMS ${SCRIPTS}
#         DESTINATION "bin")
