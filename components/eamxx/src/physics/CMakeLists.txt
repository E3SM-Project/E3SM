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

# Add a physics subdirectory if its process(es) are active.
# Usage: eamxx_add_physics_subdirectory(<dir> [proc_name1 proc_name2 ...])
# If no proc names are given, the directory name is used as the proc name.
macro(eamxx_add_physics_subdirectory dir)
  if (ARGN)
    eamxx_physics_is_active("${ARGN}" _eamxx_should_build)
  else()
    eamxx_physics_is_active("${dir}" _eamxx_should_build)
  endif()
  if (_eamxx_should_build)
    add_subdirectory(${dir})
  endif()
endmacro()

# Add physics subdirectories, conditioned on active procs for CIME builds
eamxx_add_physics_subdirectory(p3)

if (SCREAM_DOUBLE_PRECISION)
  eamxx_add_physics_subdirectory(rrtmgp)
  eamxx_add_physics_subdirectory(cosp)
  eamxx_add_physics_subdirectory(tms)
  eamxx_add_physics_subdirectory(iop_forcing)
  eamxx_add_physics_subdirectory(zm)
else()
  message(STATUS "WARNING: RRTMGP, COSP, TMS, ZM, and IOPForcing only supported for double precision builds; skipping")
endif()

eamxx_add_physics_subdirectory(shoc)
eamxx_add_physics_subdirectory(cld_fraction cld_fraction cld_frac_net)
eamxx_add_physics_subdirectory(spa)
eamxx_add_physics_subdirectory(spc)
eamxx_add_physics_subdirectory(nudging)

if (SCREAM_ENABLE_MAM)
  eamxx_add_physics_subdirectory(mam mam4_aci mam4_optics mam4_wetscav mam4_drydep mam4_srf_online_emiss mam4_aero_microphys mam4_constituent_fluxes)
endif()

eamxx_add_physics_subdirectory(gw)
