
set(P3_SRCS
  p3.cpp
)

set(P3_HEADERS
  p3.hpp
)

add_library(p3 ${P3_SRCS})

add_subdirectory(tests)
