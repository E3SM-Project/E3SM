# Set switch file on command line when running CMake
#set(SWITCH "" CACHE STRING "Switch file, either full path, relative path from location of top-level WW3/ dir, or a switch in model/bin")

function(parse_switches)
    # parse cache variable
    set(SWITCH "" CACHE STRING "Switch file")

    # path to the top level of WW3 submodule
    set(WW3_DIR "${PROJECT_SOURCE_DIR}/ww3/src/WW3")

    # Search for switch file as a full path or in model/bin
    if(EXISTS ${SWITCH})
      set(switch_file ${SWITCH})
    else()
      set(switch_file ${WW3_DIR}/model/bin/switch_${SWITCH})
      if(NOT EXISTS ${switch_file})
        message(FATAL_ERROR "Switch file '${switch_file}' does not exist, set switch with -DSWITCH=<switch>")
      endif()
    endif()

    # Copy switch file to build dir
    configure_file(${switch_file} ${CMAKE_BINARY_DIR}/cmake/wav/switch COPYONLY)

    # Open switch file and parse switches
    file(STRINGS ${CMAKE_BINARY_DIR}/cmake/wav/switch switch_strings)
    separate_arguments(switches UNIX_COMMAND ${switch_strings})

    # Include list of src files to make file more readable
    # defines variables "ftn_src", "pdlib_src", "scrip_src", and "scripnc_src"
    include(${WW3_DIR}/model/src/cmake/src_list.cmake)

    #-------------------------
    # Determine switch specific files
    #-------------------------
    include(${WW3_DIR}/model/src/cmake/check_switches.cmake)
    # make (temporary) directory strcuture that `check_switches` expects
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
    # Copy json file so that path checked by `check_switches` function is correct
    configure_file(${WW3_DIR}/model/src/cmake/switches.json ${CMAKE_CURRENT_SOURCE_DIR}/cmake COPYONLY)
    # use WW3 fucntion to parse the switches
    check_switches("${switches}" switch_files)
    # remove the temporary directory needed to make `check_switches` work
    file(REMOVE_RECURSE ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

    set(SWITCH_SOURCES)

    # add relative path to filenames
    list(APPEND srcfiles ${ftn_src} ${switch_files})
    foreach(filename ${srcfiles} )
      list(APPEND SWITCH_SOURCES ww3/src/WW3/model/src/${filename})
    endforeach()

    # if using switch fix the relative path to needed source files
    if("SCRIP" IN_LIST switches)
      foreach(filepath ${scrip_src})
        get_filename_component(BASENAME ${filepath} NAME)
        list(APPEND SWITCH_SOURCES ww3/src/WW3/model/src/SCRIP/${BASENAME})
      endforeach()
    endif()

    # if using switch fix the relative path to needed source files
    if("SCRIPNC" IN_LIST switches)
      foreach(filepath ${scripnc_src})
        get_filename_component(BASENAME ${filepath} NAME)
        list(APPEND SWITCH_SOURCES ww3/src/WW3/model/src/SCRIP/${BASENAME})
      endforeach()
    endif()
    # manually add NETCDF metadata module
    list(APPEND SWITCH_SOURCES ww3/src/WW3/model/src/w3ounfmetamd.F90)

    # return the master list of all files needed
    set(SWITCH_SOURCES ${SWITCH_SOURCES} PARENT_SCOPE)
    # and the list of switches used to determine master list
    set(switches ${switches} PARENT_SCOPE)
endfunction()

function(cull_sources_from_switches SOURCES GEN_F90_SOURCES)

  # parse the switch file and get lists of needed source files
  parse_switches()

  # first check that GEN_F90_SOURCES is an empty list
  list(LENGTH GEN_F90_SOURCES GEN_F90_LENGTH)
  if(GEN_F90_LENGTH GREATER 0)
      message(FATAL_ERROR
        "Culling of WW3 source files based on switches does not support"
        "generated files but, LENGTH(GEN_F90_SOURCES)=${GEN_F90_LENGTH}")
  endif()

  # create a copy before iterating to avoid issues with modifying list
  set(SOURCES_CULLED ${SOURCES})

  # loop over all source files found in SourceMods directories
  foreach(filepath ${SOURCES})
    # skip files in the cpl directory because we want all those.
    if(${filepath} MATCHES "ww3/src/cpl/.+\.[Ff]90")
      continue()
    endif()

    if(NOT ${filepath} IN_LIST SWITCH_SOURCES)
      list(REMOVE_ITEM SOURCES_CULLED ${filepath})
    endif()
  endforeach()

  list(LENGTH SOURCES ORIGINAL_LENGTH)
  list(LENGTH SOURCES_CULLED CULLED_LENGTH)

  if(${ORIGINAL_LENGTH} EQUAL ${CULLED_LENGTH})
    message(FATAL_ERROR "No culling based on switches occured. Something is incorrect")
  endif()

  # return list of switches
  set(switches ${switches} PARENT_SCOPE)
  set(SOURCES_CULLED ${SOURCES_CULLED} PARENT_SCOPE)
  set(GEN_F90_SOURCES_CULLED ${GEN_F90_SOURCES} PARENT_SCOPE)
endfunction()
