# - Find 3DVar
# Find 3DVar package 


 find_library( DA_LIBRARIES NAMES var_3d libvar_3d HINTS $ENV{DA_LIBRARY} )
 message(STATUS "DA library =  ${DA_LIBRARIES}  ")
 if (DA_LIBRARIES)
    find_path (DA_INCLUDES NAMES filenames.mod mpi_str.mod HINTS $ENV{DA_INCLUDE} NO_DEFAULT_PATH)
#    add_definitions(-Dkey_trc_bfm -Dkey_INCLUDE_BFM_PELCO2)
    message(STATUS "DA include =  ${DA_INCLUDES}  ")
    set (DA_has_interfaces "YES")
 else (DA_LIBRARIES)
    set (DA_has_interfaces "NO")
 endif (DA_LIBRARIES)

 find_path (EXIT_INCLUDES NAMES c_exit.h HINTS ${DA_INCLUDES}/libfexit NO_DEFAULT_PATH)
 message(STATUS "EXIT_INCLUDES = ${EXIT_INCLUDES}")
 find_library( EXIT_LIBRARIES NAMES f_exit HINTS ${EXIT_INCLUDES} )
 message(STATUS "EXIT library =  ${EXIT_LIBRARIES}  ")
 if (EXIT_LIBRARIES)
    set (EXIT_has_interfaces "YES")
 else (EXIT_LIBRARIES)
    set (EXIT_has_interfaces "NO")
 endif (EXIT_LIBRARIES)

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (libvar_3d DEFAULT_MSG DA_has_interfaces)
mark_as_advanced (DA_LIBRARIES DA_INCLUDES EXIT_LIBRARIES EXIT_INCLUDES)