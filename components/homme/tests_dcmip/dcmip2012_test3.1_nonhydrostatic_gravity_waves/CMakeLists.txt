file(GLOB NCL_SCRIPTS    "*.ncl")                                       # get ncl-scripts
file(GLOB SHELL_SCRIPTS  "*.sh")                                        # get shell-scripts
file(GLOB NAMELISTS      "*.nl")                                        # get namelists
file(GLOB PYTHON_SCRIPTS "*.py")                                        # get python scripts
file(GLOB M_FILES "*.m")                                                # get matlab scripts

# Copy test files to build
install(PROGRAMS ${NCL_SCRIPTS} ${SHELL_SCRIPTS} ${NAMELISTS} ${PYTHON_SCRIPTS} ${M_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Create a symbolic link to the vertical coordinate directory
execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_SOURCE_DIR}/test/vcoord ${CMAKE_CURRENT_BINARY_DIR}/vcoord )

# Add subdirectory for results
add_subdirectory(movies)
