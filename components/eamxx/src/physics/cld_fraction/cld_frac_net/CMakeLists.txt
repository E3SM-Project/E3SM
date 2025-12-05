if (EAMXX_ENABLE_PYTHON)
  # Process wrapping ML emulator (either via python or LAPIS-generated c++ code)
  add_library(cld_frac_net
    cld_frac_net.cpp
    eamxx_cld_frac_net_process_interface.cpp
  )
  target_compile_definitions(cld_frac_net PUBLIC EAMXX_HAS_CLD_FRAC_NET)
  target_link_libraries (cld_frac_net PUBLIC eamxx_atm_process)
  if (TARGET eamxx_physics)
    # Add this library to eamxx_physics
    target_link_libraries(eamxx_physics INTERFACE cld_frac_net)
  endif()
endif()
