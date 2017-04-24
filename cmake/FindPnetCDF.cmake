
# - Find PnetCDF
# Find the native PnetCDF includes and library
#
#  PnetCDF_has_interfaces   - True if NetCDF found including required interfaces
#

if (PNETCDF_INCLUDES AND PNETCDF_LIBRARIES)
  # Already in cache, be silent
  set (PNETCDF_FIND_QUIETLY TRUE)
endif (PNETCDF_INCLUDES AND PNETCDF_LIBRARIES)

find_path (PNETCDF_INCLUDES pnetcdf.h HINTS $ENV{PNETCDF_INC})
message(STATUS "PNETCDF include =  ${PNETCDF_INCLUDES}  ")
find_library (PNETCDF_LIBRARIES  NAMES pnetcdf)
message(STATUS "PNETCDF library =  ${PNETCDF_LIBRARIES}  ")
mark_as_advanced(NETCDF_LIBRARIES_C)

if (PNETCDF_INCLUDES AND PNETCDF_LIBRARIES)
  set (PnetCDF_has_interfaces "YES") # will be set to NO if we're missing any interfaces
  #message (STATUS "find NetCDF")
else (PNETCDF_INCLUDES AND PNETCDF_LIBRARIES)
  set (PnetCDF_has_interfaces "NO")
  message (STATUS "Failed to find PnetCDF")
endif (PNETCDF_INCLUDES AND PNETCDF_LIBRARIES)

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (PnetCDF DEFAULT_MSG  PnetCDF_has_interfaces)
mark_as_advanced(PNETCDF_INCLUDES PNETCDF_LIBRARIES)
