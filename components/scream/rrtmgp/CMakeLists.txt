
set(RRTMGP_SRCS
  rrtmgp.cpp
)

set(RRTMGP_HEADERS
  rrtmgp.hpp
)

add_library(rrtmgp ${RRTMGP_SRCS})

add_subdirectory(tests)
