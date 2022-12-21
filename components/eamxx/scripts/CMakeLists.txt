# Generate the source file for the CF validator and build it.

set (CF_STANDARD_NAME_FILE "${PROJECT_SOURCE_DIR}/data/cf-standard-name-table.yaml"
     CACHE STRING "Location of the cf standard name yaml table")
set (CF_SCREAM_NAME_FILE "${PROJECT_SOURCE_DIR}/data/cf-scream-name-table.yaml"
     CACHE STRING "Location of the scream-specific cf name yaml table")

add_executable(query-cf-database query-cf-database.cpp)
target_compile_definitions(query-cf-database PUBLIC
                           CF_STANDARD_NAME_FILE=${CF_STANDARD_NAME_FILE}
                           CF_SCREAM_NAME_FILE=${CF_SCREAM_NAME_FILE})
target_link_libraries(query-cf-database ekat)
