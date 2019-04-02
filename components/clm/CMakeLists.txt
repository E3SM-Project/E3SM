MESSAGE(" + CLM:")
MESSAGE("${LND_DIR}")

INCLUDE_DIRECTORIES (
          ${SHR_BINARY_DIR}
          ${ESMF_BINARY_DIR}
          ${GPTL_BINARY_DIR}
          ${SEQ_BINARY_DIR}
          ${MCT_BINARY_DIR}
          ${MPEU_BINARY_DIR}
)
ADD_LIBRARY(clm ${LND_DIR}src)
SET_TARGET_PROPERTIES(clm PROPERTIES LINKER_LANGUAGE Fortran)
ADD_DEPENDENCIES(clm cice cam  driver shr seq shrcal esmf shrsys  gptl  mct)
