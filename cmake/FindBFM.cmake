# - Find NetCDF
# Find BFM package 


 find_library( BFM_LIBRARIES NAMES bfm libbfm HINTS $ENV{BFM_LIBRARY} )
 if (BFM_LIBRARIES)
    find_path (BFM_INCLUDES NAMES BFM_var_list.h HINTS $ENV{BFM_INCLUDE} NO_DEFAULT_PATH)
    add_definitions(-Dkey_trc_bfm -Dkey_INCLUDE_BFM_PELCO2)
 endif (BFM_LIBRARIES)

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (BFM DEFAULT_MSG BFM_LIBRARIES BFM_INCLUDES)