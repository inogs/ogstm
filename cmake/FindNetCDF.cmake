
# - Find NetCDF
# Find the native NetCDF includes and library
#
#  NETCDF_INCLUDES_LANG    - where to find netcdf.h, etc
#  NETCDF_LIBRARIES_LANG   - Link these libraries when using NetCDF
#  NetCDF_has_interfaces   - True if NetCDF found including required interfaces
#

if (NETCDF_INCLUDES AND NETCDF_LIBRARIES)
  # Already in cache, be silent
  set (NETCDF_FIND_QUIETLY TRUE)
endif (NETCDF_INCLUDES AND NETCDF_LIBRARIES)

find_path (NETCDF_INCLUDES_C netcdf.h HINTS NETCDF_DIR NETCDF_DIR)
message(STATUS "NETCDF   C include =  ${NETCDF_INCLUDES_C}  ")
find_library (NETCDF_LIBRARIES_C  NAMES netcdf)
message(STATUS "NETCDF   C library =  ${NETCDF_LIBRARIES_C}  ")
mark_as_advanced(NETCDF_LIBRARIES_C)

find_path (NETCDFF_INCLUDES_F90 netcdf.mod HINTS NETCDFF_DIR NETCDFF_DIR)
message(STATUS "NETCDF F90 include =  ${NETCDFF_INCLUDES_F90}  ")
find_library (NETCDFF_LIBRARIES_F90  NAMES netcdff HINTS NETCDFF_DIR )
message(STATUS "NETCDF F90 library=  ${NETCDFF_LIBRARIES_F90}  ")
mark_as_advanced(NETCDF_LIBRARIES_F90)

if (NETCDFF_INCLUDES_F90 AND NETCDFF_LIBRARIES_F90 AND NETCDF_INCLUDES_C AND NETCDF_LIBRARIES_C )
  set (NetCDF_has_interfaces "YES") # will be set to NO if we're missing any interfaces
  #message (STATUS "find NetCDF")
else (NETCDFF_INCLUDES_F90 AND NETCDFF_LIBRARIES_F90 AND NETCDF_INCLUDES_C AND NETCDF_LIBRARIES_C )
  set (NetCDF_has_interfaces "NO")
  message (STATUS "Failed to find NetCDF")
endif (NETCDFF_INCLUDES_F90 AND NETCDFF_LIBRARIES_F90 AND NETCDF_INCLUDES_C AND NETCDF_LIBRARIES_C )

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (NetCDF DEFAULT_MSG  NetCDF_has_interfaces)
mark_as_advanced (NETCDF_LIBRARIES NETCDF_INCLUDES)
