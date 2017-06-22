MESSAGE("   + CAM:")

# Logic to choose which dynamics and which physics packages
set(CAM_DYNAMICS_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/dynamics/se)
set(CAM_PHYSICS_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/physics/cam)

####### CAM Setting
IF (ACME_CAM_PLEV)
  ADD_DEFINITIONS(-DPLEV=${ACME_CAM_PLEV})
  MESSAGE("     CAM setting PLEV=${ACME_CAM_PLEV}")
ELSE()
  MESSAGE(FATAL_ERROR "CAM ACME_CAM_PLEV variable to be set in CMake script")
ENDIF()

ADD_DEFINITIONS(-DPTRM=1 -DPTRN=1 -DPTRK=1)
MESSAGE("     CAM Setting PTRM=PTRN=PTRK=1")

####### dynamics/se configurations
IF (ACME_CAM_SE_NP)
  ADD_DEFINITIONS(-DNP=${ACME_CAM_SE_NP})
  MESSAGE("     CAM SE Share setting NP=${ACME_CAM_SE_NP}")
ELSE()
  MESSAGE(FATAL_ERROR "CAM SE Requires ACME_CAM_SE_NP variable to be set in CMake script")
ENDIF()

ADD_DEFINITIONS(-D_PRIM) #needed for element_state defs used in diffusion
ADD_DEFINITIONS(-DSPMD) #needed for spmd_utils to reconcile se/spmd_dyn
ADD_DEFINITIONS(-DCAM=1)
ADD_DEFINITIONS(-DNC=1)
MESSAGE("     CAM SE Share setting NC=1 (unused)")


####### physics/cam configurations
IF (ACME_CAM_PHYSICS_CAM_PCOLS)
  ADD_DEFINITIONS(-DPCOLS=${ACME_CAM_PHYSICS_CAM_PCOLS})
  MESSAGE("     CAM Physics CAM setting PCOLS=${ACME_CAM_PHYSICS_CAM_PCOLS}")
ELSE()
  MESSAGE(FATAL_ERROR "CAM Physics CAM Requires ACME_CAM_PHYSICS_CAM_PCOLS variable to be set in CMake script")
ENDIF()

IF (ACME_CAM_PHYSICS_CAM_PSUBCOLS)
  ADD_DEFINITIONS(-DPSUBCOLS=${ACME_CAM_PHYSICS_CAM_PSUBCOLS})
  MESSAGE("     CAM Physics CAM setting PSUBCOLS=${ACME_CAM_PHYSICS_CAM_PSUBCOLS}")
ELSE()
  MESSAGE(FATAL_ERROR "CAM Physics CAM Requires ACME_CAM_PHYSICS_CAM_PSUBCOLS variable to be set in CMake script")
ENDIF()

IF (ACME_CAM_PHYSICS_CAM_PCNST)
  ADD_DEFINITIONS(-DPCNST=${ACME_CAM_PHYSICS_CAM_PCNST})
  MESSAGE("     CAM Physics CAM setting PCNST=${ACME_CAM_PHYSICS_CAM_PCNST}")
ELSE()
  MESSAGE(FATAL_ERROR "CAM Physics CAM Requires ACME_CAM_PHYSICS_CAM_PCNST variable to be set in CMake script")
ENDIF()

ADD_DEFINITIONS(-DN_RAD_CNST=30)

####### chemistry configurations
ADD_DEFINITIONS(-DMODAL_AERO -DMODAL_AERO_3MODE)
###############

FILE(GLOB CAM_DYNAMICS_SOURCE ${CAM_DYNAMICS_SOURCE_DIR}/*F90
     ../../homme/src/share/*F90 ../../homme/src/preqx/*F90)
FILE(GLOB CAM_UTILS_SOURCE utils/*F90
                           #unitdrivers/*F90 ???
)

FILE(GLOB CAM_CPL_SOURCE cpl/*F90)
FILE(GLOB CAM_CONTROL_SOURCE control/*F90 control/mpishorthand.F)

FILE(GLOB CAM_PHYSICS_SOURCE 
          physics/cam/*F90
#          physics/rrtmg/*F90
#          physics/rrtmg/ext/rrtmg_mcica/*F90
#          physics/rrtmg/ext/rrtmg_lw/*F90
#          physics/rrtmg/ext/rrtmg_sw/*F90
)

FILE(GLOB CAM_CHEMISTRY_SOURCE
          chemistry/pp_trop_mam3/*F90
          chemistry/aerosol/*F90
          chemistry/modal_aero/*F90
          chemistry/mozart/*F90
          chemistry/utils/*F90
)

INCLUDE_DIRECTORIES (
          ${SHR_BINARY_DIR}
          ${ESMF_BINARY_DIR}
          ${MPEU_BINARY_DIR}
          ${GPTL_BINARY_DIR} #Need perf_mod.mod
          ${SEQ_BINARY_DIR}  #seq_cdata_mod
          ${ACME_SHR_DIR}    #shr_assert.h
)


ADD_LIBRARY (cam ${CAM_DYNAMICS_SOURCE} ${CAM_PHYSICS_SOURCE}
                 ${CAM_CHEMISTRY_SOURCE} ${CAM_CPL_SOURCE}
                 ${CAM_UTILS_SOURCE} ${CAM_CONTROL_SOURCE} )

ADD_DEPENDENCIES(cam  shrsys shr seq gptl esmf mpeu)

