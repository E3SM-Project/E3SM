add_library(eamxx_utils
  eamxx_bfbhash.cpp
  eamxx_time_stamp.cpp
  eamxx_timing.cpp
  eamxx_repro_sum_mod.F90
  eamxx_utils.cpp
)

set_target_properties(eamxx_utils PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)
target_include_directories (eamxx_utils PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_link_libraries (eamxx_utils PUBLIC
  eamxx_core
  csm_share # for shr_reprosum_mod
  ekat::Logging
  ekat::Pack)

# We need to link gptl, used by eamxx_timing.cpp, but it is a bit problematic.
# In standalone builds, it is built by scorpio, but its target properties are PRIVATE
# In CIME builds, it is attached to the spio target (which we link to our interface
# targets pioc/piof).
# Hence, in CIME builds, just link against pioc (sharedlib are already built anyways),
# while in standalone builds, link against gptl itself
if (SCREAM_CIME_BUILD)
  target_link_libraries (eamxx_utils PUBLIC pioc)
  target_include_directories(eamxx_utils PUBLIC ${INSTALL_SHAREDPATH}/include)
else ()
  target_link_libraries (eamxx_utils PUBLIC gptl)
  target_include_directories(eamxx_utils PUBLIC ${GPTL_PATH})
endif()

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
