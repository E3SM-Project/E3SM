
if (DEFINED ENV{MPAS_TOOL_DIR})
  message(STATUS "*** Using MPAS tools from $ENV{MPAS_TOOL_DIR} ***")
  file(COPY $ENV{MPAS_TOOL_DIR}/input_gen/namelist_gen DESTINATION .)
  file(COPY $ENV{MPAS_TOOL_DIR}/input_gen/streams_gen  DESTINATION .)
  file(COPY $ENV{MPAS_TOOL_DIR}/registry/parse         DESTINATION .)
else()
  message(STATUS "*** Building MPAS tools from source ***")
  # Make build tools, need to be compiled with serial compiler.
  set(CMAKE_C_COMPILER ${SCC})

  add_executable(streams_gen input_gen/streams_gen.c input_gen/test_functions.c ../external/ezxml/ezxml.c)
  add_executable(namelist_gen input_gen/namelist_gen.c input_gen/test_functions.c ../external/ezxml/ezxml.c)
  add_executable(parse registry/parse.c registry/dictionary.c registry/gen_inc.c registry/fortprintf.c registry/utility.c ../external/ezxml/ezxml.c)

  foreach(EXEITEM streams_gen namelist_gen parse)
    target_compile_definitions(${EXEITEM} PRIVATE ${CPPDEFS})
    target_compile_options(${EXEITEM} PRIVATE "-Uvector")
    target_include_directories(${EXEITEM} PRIVATE ${INCLUDES})
  endforeach()
endif()
