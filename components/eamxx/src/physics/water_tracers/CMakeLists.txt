set(WATER_TRACERS_SRCS
  eamxx_water_tracers_process_interface.cpp
  eamxx_water_isotopes_process_interface.cpp
)

add_library(water_tracers ${WATER_TRACERS_SRCS})

# Define public compile definitions so registration guards become active
target_compile_definitions(water_tracers PUBLIC
  EAMXX_HAS_WATER_TRACERS
  EAMXX_HAS_WATER_ISOTOPES
)

target_include_directories(water_tracers PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
)

# Link required EAMxx libraries
target_link_libraries(water_tracers scream_share)

# Conditionally link into eamxx_physics INTERFACE target (P3 pattern)
if (TARGET eamxx_physics)
  target_link_libraries(eamxx_physics INTERFACE water_tracers)
endif()
