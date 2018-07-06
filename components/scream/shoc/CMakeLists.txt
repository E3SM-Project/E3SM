
file(GLOB SHOC_SRCS "*.hpp" "*.cpp")

add_library(shoc ${SHOC_SRCS})

add_subdirectory(tests)
