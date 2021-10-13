set(CLDFRAC_SRCS
  atmosphere_prescribed_aerosol.cpp
)

set(CLDFRAC_HEADERS
  atmosphere_prescribed_aerosol.hpp
  spa_functions.hpp
  spa_functions_impl.hpp
)

add_library(spa ${CLDFRAC_SRCS})
target_link_libraries(spa physics_share scream_share)

add_subdirectory(tests)
