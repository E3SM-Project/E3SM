PROJECT(CPRNC C Fortran)

ENABLE_LANGUAGE(Fortran) 

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/CMake)

find_package(Netcdf REQUIRED)

ADD_CUSTOM_COMMAND(
  OUTPUT ${CMAKE_BINARY_DIR}/compare_vars_mod.F90
  COMMAND perl ${CMAKE_SOURCE_DIR}/genf90.pl
    ${CMAKE_SOURCE_DIR}/compare_vars_mod.F90.in > ${CMAKE_BINARY_DIR}/compare_vars_mod.F90
)

# Support extended line lengths
IF (${CMAKE_Fortran_COMPILER_ID} STREQUAL GNU)
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")
ELSEIF (CMAKE_COMPILER_ID STREQUAL PGI)
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Mextend")
ELSEIF (CMAKE_COMPILER_ID STREQUAL PathScale)
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -extend-source")
ENDIF ()

INCLUDE_DIRECTORIES(
  ${Netcdf_INCLUDE_DIR}
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_BINARY_DIR}
)                

SET (CPRNC_SRCS
  dtypes.h
  ${CMAKE_BINARY_DIR}/compare_vars_mod.F90
  filestruct.F90
  utils.F90
  prec.F90
  cprnc.F90
)

ADD_EXECUTABLE(cprnc ${CPRNC_SRCS})
ADD_DEPENDENCIES(cprnc ${CMAKE_BINARY_DIR}/compare_vars_mod.F90)

TARGET_LINK_LIBRARIES(cprnc ${Netcdf_LIBRARIES})

INSTALL(TARGETS cprnc
        RUNTIME DESTINATION bin)
