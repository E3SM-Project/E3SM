# Process the libyaml subdirectory, and make sure it is processed only once

DEFINE_PROPERTY(GLOBAL
                PROPERTY E3SM_YAMLCPP_ALREADY_BUILT
                BRIEF_DOCS "Wheter yaml-cpp subdir has already been processed"
                FULL_DOCS "This property is used by cmake to ensure that yaml-cpp
                           submodule directory is only included once (with add_subdirectory).")

GET_PROPERTY(IS_YAMLCPP_ALREADY_BUILT GLOBAL PROPERTY E3SM_YAMLCPP_ALREADY_BUILT SET)

IF (NOT IS_YAMLCPP_ALREADY_BUILT)
  SET (YAMLCPP_SRC_DIR ${SCREAM_SOURCE_DIR}/../../externals/yaml-cpp)
  SET (YAMLCPP_BIN_DIR ${CMAKE_BINARY_DIR}/externals/yaml-cpp)

  SET (BUILD_TESTING FALSE)
  ADD_SUBDIRECTORY (${YAMLCPP_SRC_DIR} ${YAMLCPP_BIN_DIR})

  set (SCREAM_TPL_INCLUDE_DIRS
    ${SCREAM_TPL_INCLUDE_DIRS}
    ${YAMLCPP_SRC_DIR}/include
  )

  set (SCREAM_TPL_LIBRARY_DIRS
    ${SCREAM_TPL_LIBRARY_DIRS}
    ${YAMLCPP_BIN_DIR}
  )

  set (SCREAM_TPL_LIBRARIES
    ${SCREAM_TPL_LIBRARIES}
    yaml-cpp
  )

  SET_PROPERTY(GLOBAL PROPERTY E3SM_YAMLCPP_ALREADY_BUILT TRUE)
ENDIF()
