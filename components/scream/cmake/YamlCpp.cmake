# Process the libyaml subdirectory, and make sure it is processed only once

define_property(GLOBAL
                PROPERTY E3SM_YAMLCPP_ALREADY_BUILT
                BRIEF_DOCS "Wheter yaml-cpp subdir has already been processed"
                FULL_DOCS "This property is used by cmake to ensure that yaml-cpp
                           submodule directory is only included once (with add_subdirectory).")

get_property(IS_YAMLCPP_ALREADY_BUILT GLOBAL PROPERTY E3SM_YAMLCPP_ALREADY_BUILT SET)

if (NOT IS_YAMLCPP_ALREADY_BUILT)
  set (YAMLCPP_SRC_DIR ${SCREAM_SOURCE_DIR}/../../externals/yaml-cpp)
  set (YAMLCPP_BIN_DIR ${CMAKE_BINARY_DIR}/externals/yaml-cpp)

  set (BUILD_TESTING FALSE)
  add_subdirectory (${YAMLCPP_SRC_DIR} ${YAMLCPP_BIN_DIR})

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

  set_property(GLOBAL PROPERTY E3SM_YAMLCPP_ALREADY_BUILT TRUE)
endif()
