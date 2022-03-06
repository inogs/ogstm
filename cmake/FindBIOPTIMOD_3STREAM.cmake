# - Find BIOPTIMOD_3STREAM
# Find BIOPTIMOD_3STREAM package 


find_library( BIOPTIMOD_3STREAM_LIBRARIES NAMES adj libadj HINTS $ENV{BIOPTIMOD_3STREAM_LIBRARY} )
if (NOT "$ENV{BIOPTIMOD_3STREAM_LIBRARY}" STREQUAL "")
    set(EnvironmentVariableName2 "$ENV{BIOPTIMOD_3STREAM_LIBRARY}" CACHE INTERNAL "Copied from environment variable")
endif()

message("EnvironmentVariableName2 = ${EnvironmentVariableName2}")

 message(STATUS "BIOPTIMOD_3STREAM library =  ${BIOPTIMOD_3STREAM_LIBRARIES}  ")
 if (BIOPTIMOD_3STREAM_LIBRARIES)
	 find_path (BIOPTIMOD_3STREAM_INCLUDES NAMES adj_3stream.mod HINTS $ENV{BIOPTIMOD_3STREAM_INCLUDE} NO_DEFAULT_PATH)
    add_definitions()
    message(STATUS "BIOPTIMOD_3STREAM include =  ${BIOPTIMOD_3STREAM_INCLUDES}  ")
    set (BIOPTIMOD_3STREAM_has_interfaces "YES")
 else (BIOPTIMOD_3STREAM_LIBRARIES)
    set (BIOPTIMOD_3STREAM_has_interfaces "NO")
 endif (BIOPTIMOD_3STREAM_LIBRARIES)

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (BIOPTIMOD_3STREAM DEFAULT_MSG BIOPTIMOD_3STREAM_has_interfaces)
mark_as_advanced (BIOPTIMOD_3STREAM_LIBRARIES BIOPTIMOD_3STREAM_INCLUDES)
