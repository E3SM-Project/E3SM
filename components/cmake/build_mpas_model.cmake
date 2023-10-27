function(build_mpas_models)

  # fix for new kokkos version; undo when mpas uses new kokkos build infrastructure
  set(USE_KOKKOS FALSE)

  file(GLOB MPASCONFS "${BUILDCONF}/mpas*conf" "${BUILDCONF}/maliconf")

  foreach(ITEM IN LISTS MPASCONFS)
    get_filename_component(MPASCONF ${ITEM} NAME)
    string(REPLACE "conf" "" COMP_NAME "${MPASCONF}")

    if (COMP_NAME STREQUAL "mpaso")
      list(APPEND CORES "ocean")
      set(COMP_CLASS "ocn")
    elseif (COMP_NAME STREQUAL "mpassi")
      list(APPEND CORES "seaice")
      set(COMP_CLASS "ice")
    elseif (COMP_NAME STREQUAL "mali")
      list(APPEND CORES "landice")
      set(COMP_CLASS "glc")
    else()
      message(FATAL_ERROR "Unrecognized MPAS model ${COMP_NAME}")
    endif()

    message("Found MPAS component ${COMP_CLASS} model '${COMP_NAME}'")
  endforeach()

  if (USE_ALBANY)
    set(ALBANY True)
  endif()

  if (CORES)
    add_subdirectory("mpas-framework/src")
  endif()

endfunction(build_mpas_models)
