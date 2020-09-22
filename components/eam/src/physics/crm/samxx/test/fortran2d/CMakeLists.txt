
add_executable(fortran2d ../dmdf.F90 ../fortran_driver.F90 ${FORTRAN_SRC})
target_link_libraries(fortran2d ${NCFLAGS})
set_property(TARGET fortran2d APPEND PROPERTY COMPILE_FLAGS ${DEFS2D} )

