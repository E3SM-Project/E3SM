file(GLOB NCL_SCRIPTS    "*.ncl")                                       # get ncl-scripts
file(GLOB SHELL_SCRIPTS  "*.sh")                                        # get shell-scripts
file(GLOB NAMELISTS      "*.nl")                                        # get namelists
file(GLOB PYTHON_SCRIPTS "*.py")                                        # get python scripts
file(GLOB M_FILES        "*.m")                                         # get matlab scripts

# Copy test files to build
install(PROGRAMS ${NCL_SCRIPTS} ${SHELL_SCRIPTS} ${NAMELISTS} ${PYTHON_SCRIPTS} ${M_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_subdirectory(movies)
