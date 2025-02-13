# List cpp sources to be compiled
set(POMPEI_SRCS
  eamxx_pompei_process_interface.cpp
  pompei.cpp
)

set(POMPEI_LIBS "pompei")
add_library(pompei ${POMPEI_SRCS})

target_compile_definitions(pompei PUBLIC EAMXX_HAS_POMPEI)
target_link_libraries(pompei physics_share scream_share)
target_compile_options(pompei PUBLIC)

include (ScreamUtils)

# Add this library to eamxx_physics
target_link_libraries(eamxx_physics INTERFACE pompei)
