
set(PAM_DRIVER_SRC pam_driver.cpp
                   pam_driver.F90
                   crm_feedback.h)

add_library(pam_driver ${PAM_DRIVER_SRC})
target_include_directories(pam_driver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(pam_driver PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
include(${YAKL_HOME}/yakl_utils.cmake)
yakl_process_target(pam_driver)

