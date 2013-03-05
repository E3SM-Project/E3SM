# Set the include directories
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}
                    ${SWEQX_INCLUDE_DIRS})

# There are only three separate executables that need to be compiled
#   separately
#                                     NP  NC PLEV NTRAC USE_PIO WITH_ENERGY
createTestExec(swtcB sweqx SWEQX_SRCS  4   4    1     1   FALSE      FALSE)
