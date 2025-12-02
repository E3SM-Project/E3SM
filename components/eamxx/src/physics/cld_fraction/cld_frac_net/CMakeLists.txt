if (EAMXX_ENABLE_PYTHON)
  # Process wrapping ML emulator via python hooks
  add_library(cld_frac_net
    eamxx_cld_frac_net_process_interface.cpp
  )
  target_compile_definitions(cld_frac_net PUBLIC EAMXX_HAS_CLD_FRAC_NET)
  target_link_libraries (cld_frac_net PUBLIC eamxx_atm_process)
  if (TARGET eamxx_physics)
    # Add this library to eamxx_physics
    target_link_libraries(eamxx_physics INTERFACE cld_frac_net)
  endif()
endif()
