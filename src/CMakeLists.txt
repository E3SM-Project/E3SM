PROJECT(PIO C Fortran)
ENABLE_LANGUAGE(Fortran)
#INCLUDE(FortranCInterface)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)
SET(SRCS_C topology.c pio_nc.c pio_file.c pioc_support.c pio_lists.c pioc.c
        pioc_sc.c pio_spmd.c pio_rearrange.c pio_darray.c pio_put_nc.c 
        pio_get_nc.c bget.c)


SET(SRCS_F90 pio_nf.F90 pio.F90 pio_kinds.F90 pio_types.F90
             piolib_mod.F90 
             pio_support.F90                     )

SET(TEMPSRCF90    pionfatt_mod.F90
                                pionfput_mod.F90
                                pionfget_mod.F90
           	                piodarray.F90
        )

FOREACH(tempfile IN LISTS TEMPSRCF90)
ADD_CUSTOM_COMMAND(
	OUTPUT ${tempfile}
	COMMAND ${GENF90_PATH}/genf90.pl ${CMAKE_CURRENT_SOURCE_DIR}/${tempfile}.in > ${tempfile}
	DEPENDS  ${tempfile}.in
)
ENDFOREACH()
if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
   SET(CMAKE_Fortran_FLAGS	"${CMAKE_Fortran_FLAGS} -ffree-line-length-none")
endif()

if("${CMAKE_C_COMPILER_ID}" STREQUAL "PGI")
  SET(CMAKE_C_FLAGS	"${CMAKE_C_FLAGS} -c99")
else()
  if("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
    SET(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -std=c99")
  endif()
  if("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
    SET(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -std=c99")
  endif()
endif()

ADD_LIBRARY(pioc ${SRCS_C} )

ADD_LIBRARY(piof ${SRCS_F90} ${TEMPSRCF90})

MACRO(CAR var)
  SET(${var} ${ARGV1})
ENDMACRO(CAR)

if(PNETCDF_LIBRARY)
  MESSAGE("Adding pnetcdf  library")
  TARGET_LINK_LIBRARIES(pioc ${PNETCDF_LIBRARY})
ENDIF()
if(NETCDF_C_FOUND)
  MESSAGE("Adding netcdf library ")
#  if(DEFINED  ENV{LIB_NCAR}) 
#    TARGET_LINK_LIBRARIES(pioc netcdf)
#  else()
#    string(REPLACE "/usr/lib64/libm.so" "" netcdflib "${NETCDF_C_LIBRARY}") 
    TARGET_LINK_LIBRARIES(pioc ${NETCDF_C_LIBRARY})
#  endif()
#  LINK_DIRECTORIES(${NETCDF_LIB})
#  if(HDF5_HL_LIBRARIES)
#    if(NOT DEFINED ENV{LIB_NCAR}) 
#       string(REPLACE "/usr/lib64/libm.so" ""  HDF5_clean_C_LIBRARIES
#         "${HDF5_C_LIBRARIES}")
#       TARGET_LINK_LIBRARIES(pioc ${HDF5_HL_LIBRARIES} )
#       TARGET_LINK_LIBRARIES(pioc ${HDF5_clean_C_LIBRARIES})
#   endif()
#     LINK_DIRECTORIES(${HDF5_C_LIBRARIES})
#  endif()
ENDIF()


