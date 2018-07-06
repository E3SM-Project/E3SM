
file(GLOB RRTMGP_SRCS "*.hpp" "*.cpp")

add_library(rrtmgp ${RRTMGP_SRCS})

add_subdirectory(tests)
