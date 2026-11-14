
if (DEFINED ENV{MPAS_TOOL_DIR})
  message(STATUS "*** Using MPAS tools from $ENV{MPAS_TOOL_DIR} ***")
  add_custom_target(namelist_gen)
  add_custom_command(
    TARGET namelist_gen PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $ENV{MPAS_TOOL_DIR}/input_gen/namelist_gen ${CMAKE_CURRENT_BINARY_DIR}/namelist_gen)
  add_custom_target(streams_gen)
  add_custom_command(
    TARGET streams_gen PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $ENV{MPAS_TOOL_DIR}/input_gen/streams_gen ${CMAKE_CURRENT_BINARY_DIR}/streams_gen)
  add_custom_target(parse)
  add_custom_command(
    TARGET parse PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $ENV{MPAS_TOOL_DIR}/registry/parse ${CMAKE_CURRENT_BINARY_DIR}/parse)

  set(MPAS_NAMELIST_GEN ${CMAKE_CURRENT_BINARY_DIR}/namelist_gen PARENT_SCOPE)
  set(MPAS_STREAMS_GEN ${CMAKE_CURRENT_BINARY_DIR}/streams_gen PARENT_SCOPE)
  set(MPAS_PARSE ${CMAKE_CURRENT_BINARY_DIR}/parse PARENT_SCOPE)
else()
  message(STATUS "*** Building MPAS tools from source ***")
  # The tools run on the build node during the build, so they are built with
  # the serial C compiler, as in the Make build. CMake does not officially
  # support changing CMAKE_C_COMPILER after project(), but the Makefile and
  # Ninja generators honor it for the targets in this directory. Where tools
  # built this way cannot run, set MPAS_TOOL_DIR to prebuilt tools instead.
  if (SCC)
    set(CMAKE_C_COMPILER ${SCC})
  endif()

  # ezxml is compiled once for all three tools
  add_library(ezxml_tools OBJECT ../external/ezxml/ezxml.c)
  add_executable(streams_gen input_gen/streams_gen.c input_gen/test_functions.c $<TARGET_OBJECTS:ezxml_tools>)
  add_executable(namelist_gen input_gen/namelist_gen.c input_gen/test_functions.c $<TARGET_OBJECTS:ezxml_tools>)
  add_executable(parse registry/parse.c registry/dictionary.c registry/gen_inc.c registry/fortprintf.c registry/utility.c $<TARGET_OBJECTS:ezxml_tools>)

  foreach(EXEITEM ezxml_tools streams_gen namelist_gen parse)
    target_compile_definitions(${EXEITEM} PRIVATE ${CPPDEFS})
    target_compile_options(${EXEITEM} PRIVATE "-Uvector")
    target_include_directories(${EXEITEM} PRIVATE ${INCLUDES})
  endforeach()

  set(MPAS_NAMELIST_GEN $<TARGET_FILE:namelist_gen> PARENT_SCOPE)
  set(MPAS_STREAMS_GEN $<TARGET_FILE:streams_gen> PARENT_SCOPE)
  set(MPAS_PARSE $<TARGET_FILE:parse> PARENT_SCOPE)
endif()
