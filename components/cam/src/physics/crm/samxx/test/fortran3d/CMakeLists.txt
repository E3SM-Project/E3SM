
add_executable(fortran3d ../dmdf.F90 ../fortran_driver.F90 ${FORTRAN_SRC})
target_link_libraries(fortran3d ${NCFLAGS})
set_property(TARGET fortran3d APPEND PROPERTY COMPILE_FLAGS ${DEFS3D} )

