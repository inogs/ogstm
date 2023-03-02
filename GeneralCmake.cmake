# CMake project file for OGSTM
# author : E.Pascolo, S.Bna, L.Calori

# CMAKE setting
cmake_minimum_required (VERSION 3.20)
project (OGSTM Fortran C)
set (CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/)
enable_language(Fortran C)
set (NETCDF_F90 "YES")
set(CMAKE_NO_SYSTEM_FROM_IMPORTED TRUE)
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

add_definitions(-Dkey_trahdfcoef1d -Dkey_trahdfbilap -Dkey_trc_smolar)
add_definitions(-Dkey_trc_hdfbilap -Dkey_trc_dmp -Dkey_kef -Dkey_trc_sed )
add_definitions(-Dkey_mpp -Dkey_mpp_mpi)
IF (BFMv2)
    add_definitions(-DBFMv2)
ENDIF()

if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS_RELEASE " -fno-math-errno -O2 -xAVX -qopt-report5 -g -cpp -align array64byte") #-qopenmp
    set(CMAKE_Fortran_FLAGS_DEBUG   " -O0 -g -cpp -CB -fp-stack-check -check all -traceback -gen-interfaces -warn interfaces -fpe0 -extend-source")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS_RELEASE " -O2  -fimplicit-none -cpp -ffree-line-length-0 ")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-cpp -O0 -g3 -ggdb -Wall -Wextra -fno-omit-frame-pointer -fbounds-check -fimplicit-none -pedantic -ffpe-trap=invalid,zero,overflow")
    set(CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC|PGI")
    set(CMAKE_Fortran_FLAGS_RELEASE "-fast -Mipa=fast,inline -acc -gpu=cc70,pinned -Mextend -Mpreprocess")
    set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -traceback -acc -gpu=cc70,pinned,debug -Mbounds -Mextend -Mpreprocess -Minfo=all")
else()
    message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
    message ("Fortran compiler: " ${CMAKE_Fortran_COMPILER_ID})
    message ("No optimized Fortran compiler flags are known, we just try -O2...")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif()

#include
include_directories(${MPI_Fortran_INCLUDE_PATH})
link_directories(${MPI_Fortran_LIBRARIES})
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
target_link_libraries(ogstm.xx ogstm_lib ${NETCDFF_LIBRARIES_F90} ${BFM_LIBRARIES} MPI::MPI_Fortran)
