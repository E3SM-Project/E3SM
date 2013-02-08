# Set the include directories
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}
                    ${SWEQX_INCLUDE_DIRS})

# There are only three separate executables that need to be compiled
#   separately
#                                     NP  PLEV  USE_PIO  WITH_ENERGY
createTestExec(swtcA sweqx SWEQX_SRCS  8    20    FALSE        TRUE)
