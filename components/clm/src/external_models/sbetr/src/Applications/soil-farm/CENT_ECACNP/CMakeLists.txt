set(CENT_ECACNP_SOURCES
  BgcCentCnpIndexType.F90
  BgcCentCnpForcType.F90
  BgcCentCnpCompetType.F90
  BgcCentCnpDecompType.F90
  BgcCentCnpNitDenType.F90
  BgcCentCnpSOMType.F90
  BgcCentCnpType.F90
  SOMStateVarUpdateMod.F90
  BgcReactionsCentECACnpType.F90
  PlantSoilBgcCnpType.F90
)

if (BETR_BGC EQUAL 1)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_util)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_math)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_bgc)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_grid)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_dtype)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_core)
  include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/bgcfarm_util)
  include(add_betr_library)
  add_betr_library(cent_ecacnp ${CENT_ECACNP_SOURCES})

  set(BETR_LIBRARIES cent_ecacnp;${BETR_LIBRARIES} PARENT_SCOPE)
  set(BETR_LIBRARIES cent_ecacnp;${BETR_LIBRARIES})

#X#add_subdirectory(tests)

  if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
    install(TARGETS cent_ecacnp DESTINATION lib)
    file(GLOB HEADERS *.h)
    install(FILES ${HEADERS} DESTINATION include/soil-farm/cent_ecacnp)
  endif()
else()
  include_directories(${CMAKE_BINARY_DIR}/src/betr_share/betr_util)
  include_directories(${CMAKE_BINARY_DIR}/src/betr_share/betr_math)
  include_directories(${CMAKE_BINARY_DIR}/src/betr_share/betr_grid)
  include_directories(${CMAKE_BINARY_DIR}/src/betr_share/betr_dtype)
  include_directories(${CMAKE_BINARY_DIR}/src/betr_share/betr_core)
  include_directories(${CMAKE_BINARY_DIR}/src/bgcfarm/bgcfarm_util)

  include(add_bgcfarm_library)
  add_bgcfarm_library(CENT_ECACNP ${CENT_ECACNP_SOURCES})

  set(BGCFARM_LIBRARIES CENT_ECACNP;${BGCFARM_LIBRARIES} PARENT_SCOPE)
  set(BGCFARM_LIBRARIES CENT_ECACNP;${BGCFARM_LIBRARIES})

#add_subdirectory(tests)
  if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
    install(TARGETS CENT_ECACNP DESTINATION lib)
    file(GLOB HEADERS *.h)
    install(FILES ${HEADERS} DESTINATION include/bgcfarm/CENT_ECACNP)
  endif()
endif()
