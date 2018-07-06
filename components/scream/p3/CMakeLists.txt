
file(GLOB P3_SRCS "*.hpp" "*.cpp")

add_library(p3 ${P3_SRCS})

add_subdirectory(tests)
