# Need to add Fortran as a language, since that is what the shr source is in
include(CheckLanguage)
check_language(Fortran)
if(CMAKE_Fortran_COMPILER)
    enable_language(Fortran)
else()
    message(STATUS "No Fortran support")
endif()

# Set variables needed for processing genf90 templates
set(CIMEROOT ${SCREAM_BASE_DIR}/../../cime)
list(APPEND CMAKE_MODULE_PATH ${CIMEROOT}/src/CMake)
set(GENF90 ${CIMEROOT}/src/externals/genf90/genf90.pl)
set(ENABLE_GENF90 True)
include(genf90_utils)
include(Sourcelist_utils)

# GENF90_SOURCE lists source files we will need to run through the genf90 perl script
set (GENF90_SOURCE 
    ${SCREAM_BASE_DIR}/../../share/util/shr_infnan_mod.F90.in
    ${SCREAM_BASE_DIR}/../../share/util/shr_assert_mod.F90.in
)
# FORTRAN_SOURCE lists the source files we want to build that do NOT need to be run through the genf90
# perl script. We will append to this list below with our processed genf90 files.
set (FORTRAN_SOURCE
    shr_orb_mod_c2f.F90
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

# These CPP macros are needed in shr_infnan_mod
target_compile_definitions(share_util PUBLIC
  $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CXX_COMPILER_ID:GNU>>:CPRGNU>
  $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CXX_COMPILER_ID:Intel>>:CPRINTEL>)
