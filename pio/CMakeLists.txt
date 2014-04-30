IF( NOT GENF90_PATH) 
  SET (GENF90_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin)
ENDIF()

PROJECT(PIO C Fortran)
ENABLE_LANGUAGE(Fortran)
#INCLUDE(FortranCInterface)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)

SET(SRCS_C topology.c pio_nc.c pio_file.c pioc_support.c pio_lists.c pioc.c
        pioc_sc.c pio_spmd.c pio_rearrange.c pio_darray.c pio_put_nc.c pio_get_nc.c)


SET(SRCS_F90 pio_nf.F90 pio.F90 pio_kinds.F90 pio_types.F90
             piolib_mod.F90 pio_mpi_utils.F90 
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
if("${CMAKE_C_COMPILER_ID}" STREQUAL "INTEL")
  SET(CMAKE_C_FLAGS	"${CMAKE_C_FLAGS} -O0 -traceback -ftrapuv -debug all  ")
endif()

SET(CMAKE_Fortran_FLAGS	"${CMAKE_Fortran_FLAGS} -g")

SET(CMAKE_C_FLAGS	"${CMAKE_C_FLAGS} -g")
if("${CMAKE_C_COMPILER_ID}" STREQUAL "PGI")
  SET(CMAKE_C_FLAGS	"${CMAKE_C_FLAGS} -c99")
else()
  SET(CMAKE_C_FLAGS	"${CMAKE_C_FLAGS} -std=c99")
endif()

ADD_LIBRARY(pioc ${SRCS_C} )

ADD_LIBRARY(piof ${SRCS_F90} ${TEMPSRCF90})


if(PNETCDF_LIBRARY)
  MESSAGE("Adding pnetcdf  library ${PNETCDF_LIBRARIES}")
  TARGET_LINK_LIBRARIES(pioc ${PNETCDF_LIBRARIES})
ENDIF()
if(NETCDF_FOUND)
  MESSAGE("Adding netcdf library ${NETCDF_C_LIBRARY}")
  TARGET_LINK_LIBRARIES(pioc ${NETCDF_C_LIBRARY})
ENDIF()
TARGET_LINK_LIBRARIES(pioc ${ADDITIONAL_LIBS})



