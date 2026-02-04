# - Find FABM
# Find FABM package 


message("EnvironmentVariableName = ${EnvironmentVariableName}")
find_path (YAML_INCLUDES NAMES yaml.mod  yaml_settings.mod  yaml_types.mod HINTS $ENV{YAML_INCLUDE} NO_DEFAULT_PATH)
message(STATUS "YAML include =  ${YAML_INCLUDES}  ")
set (FABM_has_interfaces "YES")

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (YAML DEFAULT_MSG FABM_has_interfaces)
mark_as_advanced (YAML_INCLUDES)
