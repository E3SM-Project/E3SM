# Need to add Fortran as a language, since that is what the shr source is in
enable_language(Fortran)

# Set variables needed for processing genf90 templates
set(CIMEROOT ${SCREAM_BASE_DIR}/../../cime)
list(APPEND CMAKE_MODULE_PATH ${CIMEROOT}/src/CMake)
set(GENF90 ${CIMEROOT}/src/externals/genf90/genf90.pl)
set(ENABLE_GENF90 True)
include(genf90_utils)
include(Sourcelist_utils)

# NOTE: I should be able to do this:
#    # Note that this sets the share_sources variable to the list of source files we need to build
#    add_subdirectory(${SCREAM_BASE_DIR}/../../share/util ${CMAKE_CURRENT_BINARY_DIR}/share_util_build)
# However, that did not seem to add the fortran sources derived from the genf90 templates to the targets.
# Instead, try explicitly adding sources and process the genf90 templates ourselves.
# TODO: Figure out why the CMakeLists.txt and the cmake macros there do not handle this correctly.
# Note also that our E3SM cmake build does not use the CMakeLists.txt in the share/util directory either
# apparently, opting to do this manually as well.
# GENF90_SOURCE lists source files we will need to run through the genf90 perl script
set (GENF90_SOURCE 
    ${SCREAM_BASE_DIR}/../../share/util/shr_infnan_mod.F90.in
    ${SCREAM_BASE_DIR}/../../share/util/shr_assert_mod.F90.in
)
# FORTRAN_SOURCE lists the source files we want to build that do NOT need to be run through the genf90
# perl script. We will append to this list below with our processed genf90 files.
set (FORTRAN_SOURCE
    shr_orb_mod_c.F90
    ${SCREAM_BASE_DIR}/../../share/util/shr_orb_mod.F90
    ${SCREAM_BASE_DIR}/../../share/util/shr_kind_mod.F90
    ${SCREAM_BASE_DIR}/../../share/util/shr_sys_mod.F90
    ${SCREAM_BASE_DIR}/../../share/util/shr_const_mod.F90
    ${SCREAM_BASE_DIR}/../../share/util/shr_log_mod.F90
    ${SCREAM_BASE_DIR}/../../share/util/shr_abort_mod.F90
    ${SCREAM_BASE_DIR}/../../share/util/shr_mpi_mod.F90
    ${SCREAM_BASE_DIR}/../../share/util/shr_strconvert_mod.F90
)
# Process genf90 template files. This adds a custom command (and hence target) for each f90 source
# that needs to be built from the genf90 template files listed in GENF90_SOURCE.
foreach (SRC_FILE ${GENF90_SOURCE})
    string(REPLACE ".in" "" SRC_FILE_STRIPPED ${SRC_FILE})
    get_filename_component(BASENAME ${SRC_FILE_STRIPPED} NAME)
    set(SRC_FILE_OUT "${CMAKE_CURRENT_BINARY_DIR}/${BASENAME}")
    add_custom_command (
        OUTPUT ${SRC_FILE_OUT}
        COMMAND ${CIMEROOT}/src/externals/genf90/genf90.pl
        ${SRC_FILE} > ${SRC_FILE_OUT}
        DEPENDS ${SRC_FILE} genf90)
    list(APPEND FORTRAN_SOURCE ${SRC_FILE_OUT})
endforeach ()
set(share_sources ${FORTRAN_SOURCE})

add_library(share_util ${share_sources})
