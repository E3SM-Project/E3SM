#==============================================================================
#  GENERATE TARGET SOURCES
#==============================================================================

SET(GENERATED_SRCS pio_file_simple_tests.F90
                   pio_init_finalize.F90
                   pio_fail.F90 
                   pio_file_fail.F90
                   ncdf_simple_tests.F90
                   ncdf_get_put.F90
                   ncdf_fail.F90 
                   pio_decomp_tests.F90
                   pio_decomp_fillval.F90)

foreach (SRC_FILE IN LISTS GENERATED_SRCS)
    add_custom_command (
        OUTPUT ${SRC_FILE}
        COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/util/pio_tf_f90gen.pl
                --annotate-source
                --out=${CMAKE_CURRENT_BINARY_DIR}/${SRC_FILE}
                ${CMAKE_CURRENT_SOURCE_DIR}/${SRC_FILE}.in
        DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${SRC_FILE}.in)
endforeach ()

#==============================================================================
#  DEFINE THE TARGETS
#==============================================================================

#===== pio_init_finalize =====
add_executable (pio_init_finalize pio_init_finalize.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/util/pio_tutil.F90)
target_link_libraries (pio_init_finalize piof pioc)

#===== pio_file_simple_tests =====
add_executable (pio_file_simple_tests pio_file_simple_tests.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/util/pio_tutil.F90)
target_link_libraries (pio_file_simple_tests piof pioc)
#if (ENABLE_TIMING)
#    target_link_libraries (pio_file_simple_tests gptl)
#endif ()

#===== pio_file_fail =====
add_executable (pio_file_fail pio_file_fail.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/util/pio_tutil.F90)
target_link_libraries (pio_file_fail piof pioc)

#===== ncdf_simple_tests =====
add_executable (ncdf_simple_tests ncdf_simple_tests.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/util/pio_tutil.F90)
target_link_libraries (ncdf_simple_tests piof pioc)

#===== ncdf_get_put =====
add_executable (ncdf_get_put ncdf_get_put.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/util/pio_tutil.F90)
target_link_libraries (ncdf_get_put piof pioc)

#===== ncdf_fail =====
add_executable (ncdf_fail ncdf_fail.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/util/pio_tutil.F90)
target_link_libraries (ncdf_fail piof pioc)

#===== pio_decomp_tests =====
add_executable (pio_decomp_tests pio_decomp_tests.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/util/pio_tutil.F90)
target_link_libraries (pio_decomp_tests piof pioc)

#===== pio_decomp_fillval =====
add_executable (pio_decomp_fillval pio_decomp_fillval.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/util/pio_tutil.F90)
target_link_libraries (pio_decomp_fillval piof pioc)
