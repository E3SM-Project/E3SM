MESSAGE("   + RTM:")
FILE(GLOB RTM_SOURCE riverroute/*.F90 cpl/*.F90)

INCLUDE_DIRECTORIES (
      ${SHR_BINARY_DIR}
      ${SEQ_BINARY_DIR}
      ${GPTL_BINARY_DIR}
      ${ESMF_BINARY_DIR}
      ${ACME_SHR_DIR}
)

#ADD_DEFINITIONS(-ffree-form)
# 
ADD_LIBRARY (rtm ${RTM_SOURCE})
ADD_DEPENDENCIES(rtm  shr seq gptl esmf)
