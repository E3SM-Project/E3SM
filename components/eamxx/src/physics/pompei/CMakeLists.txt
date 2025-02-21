
# List cpp sources to be compiled
# In this case we have the pompei model code in pompei.cpp
# and the interface code between pompei and eamxx.
set(POMPEI_SRCS
  eamxx_pompei_process_interface.cpp
  pompei.cpp
)

# Create a library of the pompei code.
set(POMPEI_LIBS "pompei")
add_library(pompei ${POMPEI_SRCS})
target_compile_definitions(pompei PUBLIC EAMXX_HAS_POMPEI)
target_link_libraries(pompei physics_share scream_share)
target_compile_options(pompei PUBLIC)

include (ScreamUtils)

# Add this library to eamxx_physics
target_link_libraries(eamxx_physics INTERFACE pompei)
