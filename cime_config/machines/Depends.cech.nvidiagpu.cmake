# Depends.cech.nvidiagpu.cmake
# Machine- and compiler-specific optimization overrides for cech with nvidiagpu compiler.

list(APPEND REDUCE_OPT_LIST
  homme/src/share/derivative_mod.F90
)

# -Mnovect avoids an internal compiler error in this file (with older nvidia compilers).
# NVHPC >= 25.9 has fixed the ICE; skip the flag for those versions.
if (NOT DEBUG AND CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "25.9")
  foreach(ITEM IN LISTS REDUCE_OPT_LIST)
    e3sm_add_flags("${ITEM}" " -Mnovect")
  endforeach()
endif()

# 1. Non-Atmosphere Support Components (Fortran & C):
# In atmosphere-centric runs (e.g. F2010-SCREAMv1), ELM and CICE are not compute-critical.
# nvfortran compiles them very slowly at -O2 due to extensive module dependencies and AST analysis.
# Lowering their default to -O1 shaves substantial build time without hurting simulation throughput.
if (NOT DEBUG AND (COMP_NAME STREQUAL "elm" OR COMP_NAME STREQUAL "cice"))
  if (DEFINED CMAKE_Fortran_FLAGS_RELEASE)
    string(REPLACE "-O2" "-O1" CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}")
  endif()
  if (DEFINED CMAKE_C_FLAGS_RELEASE)
    string(REPLACE "-O2" "-O1" CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}")
  endif()
  add_compile_options("-O1")
endif()

# 2. MPAS Framework Components (mpaso / mpassi / mali / common):
# On cech (nvidiagpu), nvfortran (fort2) OOM-kills (signal 9) or hangs for >70 minutes
# when optimizing large MPAS ocean/seaice files with deep derived-type pointer chains
# (e.g. mpas_ocn_particle_list.f90, mpas_ocn_rpn_calculator.f90,
# mpas_ocn_lagrangian_particle_tracking*.f90, mpas_ocn_diagnostics.f90).
# Note: build_mpas_model.cmake populates CORES and then calls add_subdirectory("mpas-framework/src")
# where Depends.cech.nvidiagpu.cmake is included with CORES defined.
# Since ATM is the focus on this machine, compiling MPAS at -O0 completely avoids
# nvfortran optimization hangs and OOMs while allowing all coupled/MPAS tests to build and pass.
if (NOT DEBUG AND (DEFINED CORES OR COMP_NAME STREQUAL "mpaso" OR COMP_NAME STREQUAL "mpassi" OR COMP_NAME STREQUAL "mali"))
  if (DEFINED CMAKE_Fortran_FLAGS_RELEASE)
    string(REPLACE "-O2" "-O0" CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}")
    string(REPLACE "-O1" "-O0" CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}")
  endif()
  if (DEFINED CMAKE_C_FLAGS_RELEASE)
    string(REPLACE "-O2" "-O0" CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}")
    string(REPLACE "-O1" "-O0" CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}")
  endif()
  add_compile_options("-O0")
endif()

# 3. SCREAM / EAMxx C++ View Utilities, ETI, and Process Interfaces:
# NVHPC template instantiation and optimization at -O2/-O3 is very slow for view sync,
# masking utilities, ETI template instantiations, and glue interfaces.
# Apply -O1 once during the final 'cpl' component pass (after build_eamxx() has defined
# the EAMxx targets) using TARGET_DIRECTORY so no changes to components/cmake/* are needed.
if (NOT DEBUG AND COMP_NAME STREQUAL "cpl")
  set(SCREAM_DEOPT_FIELD_FILES
    eamxx/src/share/field/field_sync.cpp
    eamxx/src/share/field/utils/compute_mask.cpp
    # ETI update masked
    eamxx/src/share/field/eti/field_update_masked_double_double.cpp
    eamxx/src/share/field/eti/field_update_masked_double_float.cpp
    eamxx/src/share/field/eti/field_update_masked_double_int.cpp
    eamxx/src/share/field/eti/field_update_masked_float_float.cpp
    eamxx/src/share/field/eti/field_update_masked_float_int.cpp
    eamxx/src/share/field/eti/field_update_masked_int_int.cpp
    # ETI update fill-aware
    eamxx/src/share/field/eti/field_update_fill_aware_double_double.cpp
    eamxx/src/share/field/eti/field_update_fill_aware_double_float.cpp
    eamxx/src/share/field/eti/field_update_fill_aware_double_int.cpp
    eamxx/src/share/field/eti/field_update_fill_aware_float_float.cpp
    eamxx/src/share/field/eti/field_update_fill_aware_float_int.cpp
    eamxx/src/share/field/eti/field_update_fill_aware_int_int.cpp
    # ETI update impl
    eamxx/src/share/field/eti/field_update_impl_double_double.cpp
    eamxx/src/share/field/eti/field_update_impl_double_float.cpp
    eamxx/src/share/field/eti/field_update_impl_double_int.cpp
    eamxx/src/share/field/eti/field_update_impl_float_float.cpp
    eamxx/src/share/field/eti/field_update_impl_float_int.cpp
    eamxx/src/share/field/eti/field_update_impl_int_int.cpp
    # ETI deep copy
    eamxx/src/share/field/eti/field_deep_copy_impl_double.cpp
    eamxx/src/share/field/eti/field_deep_copy_impl_float.cpp
    eamxx/src/share/field/eti/field_deep_copy_impl_int.cpp
    eamxx/src/share/field/eti/field_deep_copy_masked_double.cpp
    eamxx/src/share/field/eti/field_deep_copy_masked_float.cpp
    eamxx/src/share/field/eti/field_deep_copy_masked_int.cpp
    # ETI get views
    eamxx/src/share/field/eti/field_get_view_double_device.cpp
    eamxx/src/share/field/eti/field_get_view_float_device.cpp
    eamxx/src/share/field/eti/field_get_view_int_device.cpp
    eamxx/src/share/field/eti/field_get_view_double_host.cpp
    eamxx/src/share/field/eti/field_get_view_float_host.cpp
    eamxx/src/share/field/eti/field_get_view_int_host.cpp
    eamxx/src/share/field/eti/field_get_strided_view_double_device.cpp
    eamxx/src/share/field/eti/field_get_strided_view_float_device.cpp
    eamxx/src/share/field/eti/field_get_strided_view_int_device.cpp
    eamxx/src/share/field/eti/field_get_strided_view_double_host.cpp
    eamxx/src/share/field/eti/field_get_strided_view_float_host.cpp
    eamxx/src/share/field/eti/field_get_strided_view_int_host.cpp
  )

  set(SCREAM_DEOPT_RRTMGP_FILES
    eamxx/src/physics/rrtmgp/eamxx_rrtmgp_process_interface.cpp
  )

  set(SCREAM_DEOPT_ATM_FILES
    eamxx/src/mct_coupling/eamxx_cxx_f90_interface.cpp
  )

  set(SCREAM_DEOPT_REMAP_FILES
    eamxx/src/share/remap/horizontal_remapper.cpp
  )

  set(SCREAM_DEOPT_PROPERTY_CHECK_FILES
    eamxx/src/share/property_checks/field_nan_check.cpp
    eamxx/src/share/property_checks/field_within_interval_check.cpp
  )

  if (TARGET eamxx_field)
    foreach(ITEM IN LISTS SCREAM_DEOPT_FIELD_FILES)
      if (EXISTS "${PROJECT_SOURCE_DIR}/${ITEM}")
        set_property(SOURCE "${PROJECT_SOURCE_DIR}/${ITEM}" TARGET_DIRECTORY eamxx_field APPEND_STRING PROPERTY COMPILE_FLAGS " -O1")
      endif()
    endforeach()
  endif()

  if (TARGET scream_rrtmgp)
    foreach(ITEM IN LISTS SCREAM_DEOPT_RRTMGP_FILES)
      if (EXISTS "${PROJECT_SOURCE_DIR}/${ITEM}")
        set_property(SOURCE "${PROJECT_SOURCE_DIR}/${ITEM}" TARGET_DIRECTORY scream_rrtmgp APPEND_STRING PROPERTY COMPILE_FLAGS " -O1")
      endif()
    endforeach()
  endif()

  if (TARGET atm)
    foreach(ITEM IN LISTS SCREAM_DEOPT_ATM_FILES)
      if (EXISTS "${PROJECT_SOURCE_DIR}/${ITEM}")
        set_property(SOURCE "${PROJECT_SOURCE_DIR}/${ITEM}" TARGET_DIRECTORY atm APPEND_STRING PROPERTY COMPILE_FLAGS " -O1")
      endif()
    endforeach()
  endif()

  if (TARGET eamxx_remap)
    foreach(ITEM IN LISTS SCREAM_DEOPT_REMAP_FILES)
      if (EXISTS "${PROJECT_SOURCE_DIR}/${ITEM}")
        set_property(SOURCE "${PROJECT_SOURCE_DIR}/${ITEM}" TARGET_DIRECTORY eamxx_remap APPEND_STRING PROPERTY COMPILE_FLAGS " -O1")
      endif()
    endforeach()
  endif()

  if (TARGET eamxx_property_checks)
    foreach(ITEM IN LISTS SCREAM_DEOPT_PROPERTY_CHECK_FILES)
      if (EXISTS "${PROJECT_SOURCE_DIR}/${ITEM}")
        set_property(SOURCE "${PROJECT_SOURCE_DIR}/${ITEM}" TARGET_DIRECTORY eamxx_property_checks APPEND_STRING PROPERTY COMPILE_FLAGS " -O1")
      endif()
    endforeach()
  endif()
endif()
