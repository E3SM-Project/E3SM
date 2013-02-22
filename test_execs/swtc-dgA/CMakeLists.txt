# Set the include directories
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}
                    ${SWDGX_INCLUDE_DIRS})

# There are only three separate executables that need to be compiled
#   separately
#                                        NP  PLEV  USE_PIO  WITH_ENERGY
createTestExec(swtc-dgA swdgx SWDGX_SRCS  6     1    FALSE        FALSE)
