
#TODO -  Glob all cam directories into single libcam

ADD_DEFINITIONS(-DCAM) #Needed in dynamics/se/shr

IF (ACME_CAM_PLEV)
  ADD_DEFINITIONS(-DPLEV=${ACME_CAM_PLEV})
  MESSAGE("++ CAM setting PLEV=${ACME_CAM_PLEV}")
ELSE()
  MESSAGE(FATAL_ERROR "CAM Requires ACME_CAM_PLEV variable to be set in CMake script")
ENDIF()

add_subdirectory(src)
