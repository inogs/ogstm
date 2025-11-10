# - Find FABM
# Find FABM package 


 find_library( FABM_LIBRARIES NAMES fabm libfabm HINTS $ENV{FABM_LIBRARY} )
 message(STATUS "FABM library =  ${FABM_LIBRARIES}  ")
 if (NOT "$ENV{FABM_LIBRARY}" STREQUAL "")
	 set(EnvironmentVariableName "$ENV{FABM_LIBRARY}" CACHE INTERNAL "Copied from environment variable")
endif()

message("EnvironmentVariableName = ${EnvironmentVariableName}")
 if (FABM_LIBRARIES)
    find_path (FABM_INCLUDES NAMES fabm.h HINTS $ENV{FABM_INCLUDE} NO_DEFAULT_PATH)
    add_definitions(-Dkey_trc_fabm -Dkey_INCLUDE_FABM_PELCO2)
    message(STATUS "FABM include =  ${FABM_INCLUDES}  ")
    set (FABM_has_interfaces "YES")
 else (FABM_LIBRARIES)
    set (FABM_has_interfaces "NO")
 endif (FABM_LIBRARIES)

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (FABM DEFAULT_MSG FABM_has_interfaces)
mark_as_advanced (FABM_LIBRARIES FABM_INCLUDES)
