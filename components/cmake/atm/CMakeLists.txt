build_model("${COMP_CLASS}" "${COMP_NAME}")

if (COMP_NAME STREQUAL "eatm")
  target_compile_features(${COMP_CLASS} PRIVATE cxx_std_17)
  target_include_directories(${COMP_CLASS} PRIVATE
    ${PROJECT_SOURCE_DIR}/components/emulators/eatm/src
    ${PROJECT_SOURCE_DIR}/components/emulators/common/src
  )
endif()
