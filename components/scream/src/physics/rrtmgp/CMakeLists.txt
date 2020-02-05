
set(RRTMGP_SRCS
  rrtmgp.cpp
)

set(RRTMGP_HEADERS
  rrtmgp.hpp
)

add_library(rrtmgp ${RRTMGP_SRCS})

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
