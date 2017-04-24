
# - Find PETSc
# Find the native PETSc includes and library
#
#  Petsc_has_interfaces   - True if NetCDF found including required interfaces
#

if (PETSC_INCLUDES AND PETSC_LIBRARIES)
  # Already in cache, be silent
  set (PETSC_FIND_QUIETLY TRUE)
endif (PETSC_INCLUDES AND PETSC_LIBRARIES)

find_path (PETSC_INCLUDES petscao.h HINTS $ENV{PETSC_INC})
message(STATUS "PETSC include =  ${PETSC_INCLUDES}  ")
find_library (PETSC_LIBRARIES  NAMES petsc HINTS $ENV{PETSC_LIB})
message(STATUS "PETSC library =  ${PETSC_LIBRARIES}  ")
mark_as_advanced(PETSC_LIBRARIES)

if (PETSC_INCLUDES AND PETSC_LIBRARIES)
  set (Petsc_has_interfaces "YES") # will be set to NO if we're missing any interfaces
else (PETSC_INCLUDES AND PETSC_LIBRARIES)
  set (Petsc_has_interfaces "NO")
  message (STATUS "Failed to find PETSc")
endif (PETSC_INCLUDES AND PETSC_LIBRARIES )

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (PETSc DEFAULT_MSG  Petsc_has_interfaces)
mark_as_advanced(PETSC_INCLUDES PETSC_LIBRARIES)
