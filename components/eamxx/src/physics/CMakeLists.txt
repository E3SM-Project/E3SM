add_library(eamxx_physics INTERFACE)

# When doing a CIME build, query namelist_scream.xml to determine which physics
# processes are active in this simulation, then only build those packages.
if (SCREAM_CIME_BUILD)
  execute_process(
    COMMAND ${CASEROOT}/atmquery --value --grep atm_procs_list
    WORKING_DIRECTORY ${CASEROOT}
    RESULT_VARIABLE EAMXX_QUERY_PROCS_STATUS
    OUTPUT_VARIABLE EAMXX_ACTIVE_ATM_PROCS_RAW
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
  )
  if (EAMXX_QUERY_PROCS_STATUS EQUAL 0)
    # Output is multiple lines of comma-separated proc names; flatten into one CMake list
    string(REPLACE "\n" "," EAMXX_ACTIVE_ATM_PROCS_CSV "${EAMXX_ACTIVE_ATM_PROCS_RAW}")
    string(REGEX REPLACE "[ \t]*,[ \t]*" ";" EAMXX_ACTIVE_ATM_PROCS "${EAMXX_ACTIVE_ATM_PROCS_CSV}")
    string(TOLOWER "${EAMXX_ACTIVE_ATM_PROCS}" EAMXX_ACTIVE_ATM_PROCS)
    list(REMOVE_DUPLICATES EAMXX_ACTIVE_ATM_PROCS)
    message(STATUS "EAMxx: Active atm procs: ${EAMXX_ACTIVE_ATM_PROCS}")
  else()
    message(WARNING "EAMxx: Could not query active atm procs from ${CASEROOT}/namelist_scream.xml; building all physics")
    set(EAMXX_ACTIVE_ATM_PROCS "")
  endif()
endif()

# Helper: returns TRUE if the physics described by proc_names should be built.
# In non-CIME builds or when the query failed (empty list), all physics are built.
# proc_names is a CMake list of exact proc names to look for in EAMXX_ACTIVE_ATM_PROCS.
function(eamxx_physics_is_active proc_names out_var)
  if (NOT SCREAM_CIME_BUILD OR NOT EAMXX_ACTIVE_ATM_PROCS)
    set(${out_var} TRUE PARENT_SCOPE)
    return()
  endif()
  foreach(_name IN LISTS proc_names)
    if (_name IN_LIST EAMXX_ACTIVE_ATM_PROCS)
      set(${out_var} TRUE PARENT_SCOPE)
      return()
    endif()
  endforeach()
  set(${out_var} FALSE PARENT_SCOPE)
endfunction()

# Add physics subdirectories, conditioned on active procs for CIME builds
eamxx_physics_is_active("p3" _build_p3)
if (_build_p3)
  add_subdirectory(p3)
endif()

if (SCREAM_DOUBLE_PRECISION)
  eamxx_physics_is_active("rrtmgp" _build_rrtmgp)
  if (_build_rrtmgp)
    add_subdirectory(rrtmgp)
  endif()

  eamxx_physics_is_active("cosp" _build_cosp)
  if (_build_cosp)
    add_subdirectory(cosp)
  endif()

  eamxx_physics_is_active("tms" _build_tms)
  if (_build_tms)
    add_subdirectory(tms)
  endif()

  eamxx_physics_is_active("iop_forcing" _build_iop_forcing)
  if (_build_iop_forcing)
    add_subdirectory(iop_forcing)
  endif()

  eamxx_physics_is_active("zm" _build_zm)
  if (_build_zm)
    add_subdirectory(zm)
  endif()
else()
  message(STATUS "WARNING: RRTMGP, COSP, TMS, ZM, and IOPForcing only supported for double precision builds; skipping")
endif()

eamxx_physics_is_active("shoc" _build_shoc)
if (_build_shoc)
  add_subdirectory(shoc)
endif()

eamxx_physics_is_active("cld_fraction;cld_frac_net" _build_cld_fraction)
if (_build_cld_fraction)
  add_subdirectory(cld_fraction)
endif()

eamxx_physics_is_active("spa" _build_spa)
if (_build_spa)
  add_subdirectory(spa)
endif()

eamxx_physics_is_active("spc" _build_spc)
if (_build_spc)
  add_subdirectory(spc)
endif()

eamxx_physics_is_active("nudging" _build_nudging)
if (_build_nudging)
  add_subdirectory(nudging)
endif()

if (SCREAM_ENABLE_MAM)
  set(_mam_procs mam4_aci;mam4_optics;mam4_wetscav;mam4_drydep;mam4_srf_online_emiss;mam4_aero_microphys;mam4_constituent_fluxes)
  eamxx_physics_is_active("${_mam_procs}" _build_mam)
  if (_build_mam)
    add_subdirectory(mam)
  endif()
endif()

eamxx_physics_is_active("gw" _build_gw)
if (_build_gw)
  add_subdirectory(gw)
endif()
