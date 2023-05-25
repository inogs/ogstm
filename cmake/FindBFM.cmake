# - Find BFM
# Find BFM package 


 find_library( BFM_LIBRARIES NAMES bfm.a libbfm.a HINTS $ENV{BFM_LIBRARY} )
 message(STATUS "BFM library =  ${BFM_LIBRARIES}  ")
 if (BFM_LIBRARIES)
    find_path (BFM_INCLUDES NAMES BFM_var_list.h HINTS $ENV{BFM_INCLUDE} NO_DEFAULT_PATH)
    add_definitions(-Dkey_trc_bfm -Dkey_INCLUDE_BFM_PELCO2)
    message(STATUS "BFM include =  ${BFM_INCLUDES}  ")
    set (BFM_has_interfaces "YES")
 else (BFM_LIBRARIES)
    set (BFM_has_interfaces "NO")
 endif (BFM_LIBRARIES)

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (BFM DEFAULT_MSG BFM_has_interfaces)
mark_as_advanced (BFM_LIBRARIES BFM_INCLUDES)
