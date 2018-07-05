
file(GLOB COUPLER_SRCS "*.hpp" "*.cpp")

add_library(coupler ${COUPLER_SRCS})

add_subdirectory(tests)
