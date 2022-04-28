macro (CreateGPTLTarget)
  # If we didn't already parse this script, proceed
  if (NOT TARGET gptl)
    if (SCREAM_CIME_BUILD)
      # Some sanity checks
      if (NOT DEFINED INSTALL_SHAREDPATH)
        message (FATAL_ERROR "Error! The cmake variable 'INSTALL_SHAREDPATH' is not defined.")
      endif ()

      # Look for libgptl in INSTALL_SHAREDPATH/lib
      find_library(GPTL_LIB gptl REQUIRED PATHS ${INSTALL_SHAREDPATH}/lib)

      # Create the imported target that scream targets can link to
      add_library (gptl INTERFACE)
      target_link_libraries (gptl INTERFACE ${GPTL_LIB})
      target_include_directories (gptl INTERFACE ${INSTALL_SHAREDPATH}/include)
      if (NOT MPILIB STREQUAL "mpi-serial")
        target_compile_definitions (gptl INTERFACE HAVE_MPI)
      endif()
    endif ()
  endif ()
endmacro()
