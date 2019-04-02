MESSAGE(" + CICE:")
MESSAGE("${ICE_DIR}")

# eventually the glob should be replaced with an add_subdirectory
#ADD_SUBDIRECTORY(${ICE_DIR}/src)
FILE(GLOB  CICE_SOURCE ${ICE_DIR}/src/*)

ADD_DEFINITIONS(-DNXGLOB -DNYGLOB -DBLCKX -DBLCKY -DMXBLCKS -DMYBLCKS -DNCAT -DNTR_AERO)
INCLUDE_DIRECTORIES (
          ${SHR_BINARY_DIR}
          ${STREAMS_BINARY_DIR}
          ${SEQ_BINARY_DIR}
          ${DRIVER_BINARY_DIR}
          ${ESMF_BINARY_DIR}
          ${GPTL_BINARY_DIR}
          ${MCT_BINARY_DIR}
          ${MPEU_BINARY_DIR}
          ${E3SM_PIO_LIBRARY_DIR}
          ${SHR_SYS_BINARY_DIR}
          ${SHR_CAL_BINARY_DIR}
)
ADD_LIBRARY(cice ${CICE_SOURCE})
SET_TARGET_PROPERTIES(cice PROPERTIES LINKER_LANGUAGE Fortran)
ADD_DEPENDENCIES(cice cam shr seq shrcal esmf shrsys  gptl  mct)
