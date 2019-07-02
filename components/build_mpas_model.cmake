function(build_mpas_models)

  file(GLOB MPASCONFS "${BUILDCONF}/mpas*conf" "${BUILDCONF}/maliconf")
  foreach(ITEM IN LISTS MPASCONFS)
    get_filename_component(MPASCONF ${ITEM} NAME)
    string(REPLACE "conf" "" MODEL "${MPASCONF}")

    if (MODEL STREQUAL "mpaso")
      list(APPEND CORES "ocean")
      list(APPEND MODELS "ocn")
    elseif (MODEL STREQUAL "mpassi")
      list(APPEND CORES "seaice")
      list(APPEND MODELS "ice")
    elseif (MODEL STREQUAL "mali")
      list(APPEND CORES "landice")
      list(APPEND MODELS "glc")
      if (USE_ALBANY)
        set(ALBANY True)
      endif()
    else()
      message(FATAL_ERROR "Unrecognized MPAS model ${MODEL_ARG}")
    endif()

    message("Found MPAS component model '${MODEL}'")
  endforeach()

  set(MODELS ${MODELS} PARENT_SCOPE)

  add_subdirectory("mpas-source/src")

endfunction(build_mpas_models)
