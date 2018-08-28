
set(COUPLER_SRCS
  coupler.cpp
)

set(COUPLER_HEADERS
  coupler.hpp
)

add_library(coupler ${COUPLER_SRCS})

add_subdirectory(tests)
