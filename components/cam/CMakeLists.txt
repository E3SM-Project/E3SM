MESSAGE(" + CAM:")
MESSAGE(${ATM_DIR})
 
MESSAGE("FIND LOCATION OF HOMME:")
GET_FILENAME_COMPONENT(PARENT_DIR ${CIME_DIR} DIRECTORY)                                                     

SET(HOMME_DIR ${PARENT_DIR}/components/homme )
MESSAGE("HOMME_DIR SET TO:" ${HOMME_DIR})
# MESSAGE("PRJECT BINARY DIR:" ${PROJECT_BINARY_DIR})

# JIM FOUCAR:  IT'S UP TO YOU TO FIGURE OUT WHAT THE VARIABLES SHOULD BE SET TO, I JUST ADDED THEM SO THE CODE WOULD COMPILE
ADD_DEFINITIONS(-DPLON -DPLAT -DPLEV -DPCOLS -DPSUBCOLS -DPTRM -DPTRN -DPTRK -DPTRN -DPTRM -DCLUBB_REAL_TYPE=8 -DNP -DMPI_TAG_UB -DPCNST -DN_RAD_CNST -DDUST_DMT_GRD)


# here's how to set an ifdef variable
SET(WITH_PRIM "ON")
IF (WITH_PRIM)
MESSAGE(STATUS "IFDEF _PRIM = TRUE")
ADD_DEFINITIONS(-D_PRIM=1)
ENDIF()


# SHOULD PROBABLY REPLACE THE GLOBS WITH ADD_SUBDIRECTORIES
FILE(GLOB_RECURSE HOMME_SHR ${HOMME_DIR}/src/share/*)
FILE(GLOB_RECURSE HOMME_PREQX  ${HOMME_DIR}/src/preqx/* 
                             # ${HOMME_DIR}/src/sweqx/*
                             # ${HOMME_DIR}/src/theta-l/*
                              ${HOMME_DIR}/src/restart_io_mod.F90)
#FILE(GLOB CAM_MPI bld/mct-serial/mpi-serial/*)

SET(WITH_CHEM "OFF")
IF (WITH_CHEM)
MESSAGE(STATUS "IFDEF CHEM = FALSE")
ADD_DEFINITIONS(-D_CHEM=0)
ENDIF()

FILE(GLOB_RECURSE HOMME_SHR ${HOMME_DIR}/src/share/*)
FILE(GLOB_RECURSE HOMME_PREQX  ${HOMME_DIR}/src/preqx/*
                             # ${HOMME_DIR}/src/sweqx/*
                             # ${HOMME_DIR}/src/theta-l/*
                              ${HOMME_DIR}/src/restart_io_mod.F90)
#FILE(GLOB CAM_MPI bld/mct-serial/mpi-serial/*)

SET(WITH_CHEM "OFF")
IF (WITH_CHEM)
MESSAGE(STATUS "IFDEF CHEM = FALSE")
ADD_DEFINITIONS(-D_CHEM=0)
ENDIF()


FILE(GLOB PHYSICS_RESTART  src/physics/cam/restart_physics.F90 src/utils/cam_dom/*)

# DEPENDING ON THE CONFIGURATION, YOU'LL NEED TO REPEAT THE LOGIC STRUCTURE FOR THE BUILD HERE BASED ON THE 
# CAM-PERL-CONFIGURE FILE
FILE(GLOB_RECURSE CAM_SOURCE src/advection/*
                             src/cpl/*
                             src/control/*
                             src/dynamics/se/*
                      #       src/physics/cam/constituents.F90
                      #       src/physics/cam/ppgrid.F90
                      #       src/physics/cam/phys_control.F90
                      #       src/physics/cam/comsrf.F90
                      #       src/physics/cam/phys_grid.F90
                             src/physics/cam/*
                             src/chemistry/utils/*
                        #     src/physics/rrtmg/*
                             src/physics/camrt/*
                             src/chemistry/mozart/*
                             src/chemistry/pp_trop_mam3/*
                             src/chemistry/modal_aero/*
                             src/chemistry/aerosol/*
                             src/utils/*
                         #    src/physics/carma/*
                         #    src/physics/clubb/*

			                 #     src/physics/*
                        #

                        #      src/chemistry/pp_none/*
                        #     src/chemistry/aero/*
                        #     src/chemistry/bulk_aero/*
				#     src/chemistry/aerosol/*
                        #     src/chemistry/modal_aero/*
                        #     src/chemistry/pp_trop_mam3/*
                                      )

LIST(REMOVE_ITEM CAM_SOURCE ${PHYSICS_RESTART})
INCLUDE_DIRECTORIES (
          ${SHR_BINARY_DIR}
          ${ESMF_BINARY_DIR}
          ${GPTL_BINARY_DIR}
          ${SEQ_BINARY_DIR}
          ${MCT_BINARY_DIR}
          ${MPEU_BINARY_DIR}
          ${PIO_BINARY_DIR}
          ${DRIVER_BINARY_DIR}
          ${SHR_CAL_BINARY_DIR}
          ${SHR_SYS_BINARY_DIR}
          ${STREAMS_BINARY_DIR}
	    ${SHR_BINARY_DIR}
)

SET(WITH_CAM "ON")
IF (WITH_CAM)
MESSAGE(STATUS "IFDEF CAM = TRUE")
ADD_DEFINITIONS(-DCAM=1)
ENDIF()

SET(WITH_SPMD "ON")
IF (WITH_SPMD)
MESSAGE(STATUS "IFDEF SPMD = TRUE")
ADD_DEFINITIONS(-DSPMD=1)
ENDIF()
SET(WITH__MPI "ON")
IF(WITH__MPI)
MESSAGE(STATUS "IFDEF _MPI = TRUE")
ADD_DEFINITIONS(-D_MPI=1)
ENDIF()


# mam3 cpp defs
SET(WITH_MODAL_AERO_3MODE "ON")
IF(WITH_MODAL_AERO_3MODE)
MESSAGE(STATUS "IFDEF MODAL_AERO_3MODE = TRUE")
ADD_DEFINITIONS(-DMODAL_AERO_3MODE=1)
ENDIF()

SET(WITH_MODAL_AERO "ON")
IF(WITH_MODAL_AERO)
MESSAGE(STATUS "IFDEF MODAL_AERO = TRUE")
ADD_DEFINITIONS(-DMODAL_AERO=1)
ENDIF()



ADD_LIBRARY(cam ${CAM_SOURCE} ${CAM_SOURCE_MORE} ${HOMME_SHR} ${HOMME_PREQX})
ADD_DEPENDENCIES(cam shr seq shrcal esmf streams shrsys gptl mct pio)
SET_TARGET_PROPERTIES(cam PROPERTIES LINKER_LANGUAGE Fortran)
