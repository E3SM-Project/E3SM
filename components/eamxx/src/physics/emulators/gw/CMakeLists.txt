add_library(gw_emulator
  eamxx_gw_emulator_process.cpp
)

target_link_libraries(gw_emulator PUBLIC eamxx_atm_process)
target_compile_definitions(gw_emulator PUBLIC EAMXX_HAS_GW_EMULATOR)
