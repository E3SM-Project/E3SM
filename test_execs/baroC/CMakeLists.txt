# Set the include directories
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}
                    ${PREQX_INCLUDE_DIRS})

# There are only three separate executables that need to be compiled
#   separately
#                                     NP  PLEV  USE_PIO  WITH_ENERGY
createTestExec(baroC preqx PREQX_SRCS  4    26    FALSE        TRUE)
