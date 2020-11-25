# Generate the source file for the CF validator and build it.
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/query-cf-database.cpp.in
  ${CMAKE_CURRENT_BINARY_DIR}/query-cf-database.cpp
  @ONLY
)
add_executable(query-cf-database ${CMAKE_CURRENT_BINARY_DIR}/query-cf-database.cpp)
target_link_libraries(query-cf-database ekat)
