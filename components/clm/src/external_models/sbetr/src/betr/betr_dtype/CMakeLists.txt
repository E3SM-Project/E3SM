set(BETR_DTYPE_SOURCES
  BeTR_aerocondType.F90
  BeTR_carbonfluxType.F90
  BeTR_nitrogenfluxType.F90
  BeTR_phosphorusfluxType.F90
  BeTR_carbonfluxRecvType.F90
  BeTR_nitrogenfluxRecvType.F90
  BeTR_phosphorusfluxRecvType.F90
  BeTR_carbonstateRecvType.F90
  BeTR_nitrogenstateRecvType.F90
  BeTR_phosphorusstateRecvType.F90
  BeTR_biogeophysInputType.F90
  BeTR_biogeoStateType.F90
  BeTR_biogeoFluxType.F90
  )


include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_util)


include(add_betr_library)
add_betr_library(betr_dtype ${BETR_DTYPE_SOURCES})

set(BETR_LIBRARIES betr_dtype;${BETR_LIBRARIES} PARENT_SCOPE)
set(BETR_LIBRARIES betr_dtype;${BETR_LIBRARIES})

#X#add_subdirectory(tests)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS betr_dtype DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/betr/dtype)
endif()
