
set(SHOC_SRCS
  shoc.cpp
)

set(SHOC_HEADERS
  shoc.hpp
)

add_library(shoc ${SHOC_SRCS})

add_subdirectory(tests)
