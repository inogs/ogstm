# - Find OASIM_ATM
# Find OASIM_ATM package


find_library( OASIM_ATM_LIBRARIES NAMES oasim liboasim HINTS $ENV{OASIM_ATM_LIBRARY} )
if (NOT "$ENV{OASIM_ATM_LIBRARY}" STREQUAL "")
    set(EnvironmentVariableName2 "$ENV{OASIM_ATM_LIBRARY}" CACHE INTERNAL "Copied from environment variable")
endif()

message("EnvironmentVariableName2 = ${EnvironmentVariableName2}")

 message(STATUS "OASIM_ATM library =  ${OASIM_ATM_LIBRARIES}  ")
 if (OASIM_ATM_LIBRARIES)
	 find_path (OASIM_ATM_INCLUDES NAMES oasim.mod HINTS $ENV{OASIM_ATM_INCLUDE} NO_DEFAULT_PATH)
    add_definitions()
    message(STATUS "OASIM_ATM include =  ${OASIM_ATM_INCLUDES}  ")
    set (OASIM_ATM_has_interfaces "YES")
 else (OASIM_ATM_LIBRARIES)
    set (OASIM_ATM_has_interfaces "NO")
 endif (OASIM_ATM_LIBRARIES)

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (OASIM_ATM DEFAULT_MSG OASIM_ATM_has_interfaces)
mark_as_advanced (OASIM_ATM_LIBRARIES OASIM_ATM_INCLUDES)
