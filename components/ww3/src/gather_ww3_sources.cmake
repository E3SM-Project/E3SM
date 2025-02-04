# Set switch file on command line when running CMake
set(SWITCH "" CACHE STRING "Switch file, either full path, relative path from location of top-level WW3/ dir, or a switch in model/bin")

# Search for switch file as a full path or in model/bin
if(EXISTS ${SWITCH})
  set(switch_file ${SWITCH})
else()
    set(switch_file ${CMAKE_CURRENT_SOURCE_DIR}/../../ww3/src/WW3/model/bin/switch_${SWITCH})
  if(NOT EXISTS ${switch_file})
    message(FATAL_ERROR "Switch file '${switch_file}' does not exist, set switch with -DSWITCH=<switch>")
  endif()
endif()

message(STATUS "Build with switch: ${switch_file}")
# Copy switch file to build dir
configure_file(${switch_file} ${CMAKE_BINARY_DIR}/switch COPYONLY)

# Open switch file
file(STRINGS ${CMAKE_BINARY_DIR}/switch switch_strings)
separate_arguments(switches UNIX_COMMAND ${switch_strings})

# Include list of src files to make file more readable
# defines variables "ftn_src", "pdlib_src", "scrip_src", and "scripnc_src"
include(${CMAKE_CURRENT_SOURCE_DIR}/../../ww3/src/WW3/model/src/cmake/src_list.cmake)

#-------------------------
# Determine switch specific files
# Include check_switches as a function for less verbosity in this CMakeLists.txt
#-------------------------
#message(STATUS "switches are : ${switches}")
include(${CMAKE_CURRENT_SOURCE_DIR}/../../ww3/src/WW3/model/src/cmake/check_switches.cmake)
# Copy json file so that path checked by `check_switches` function is correct
file(MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../../ww3/src/WW3/model/src/cmake/switches.json ${CMAKE_CURRENT_SOURCE_DIR}/cmake COPYONLY)
check_switches("${switches}" switch_files)
file(REMOVE_RECURSE ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

message(STATUS "---")
message(STATUS "list of always source files is : ${scrip_src}")
#message(STATUS "list of switch files is : ${switch_files}")
message(STATUS "---")

set(BASENAME_SET)
set(SOURCES_RESULT)
set(GEN_F90_SOURCES_RESULT)

list(APPEND srcfiles ${ftn_src} ${switch_files})
foreach(file ${srcfiles} )
  list(APPEND SOURCES_RESULT ww3/src/WW3/model/src/${file})
  list(APPEND BASENAME_SET ${file})
endforeach()

# add the coupler files to source list
file(GLOB MATCHES RELATIVE "${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/ww3/src/cpl/*.[Ff]90")
foreach (MATCH IN LISTS MATCHES)
  get_filename_component(BASENAME ${MATCH} NAME)
  list(FIND BASENAME_SET ${BASENAME} BASENAME_WAS_FOUND)
  if (BASENAME_WAS_FOUND EQUAL -1)
    list(APPEND SOURCES_RESULT ${MATCH})
    list(APPEND BASENAME_SET ${BASENAME})
  else()
    message("Warning: Skipping repeated base filename ${BASENAME} for ${MATCH}")
  endif()
endforeach()

if("SCRIP" IN_LIST switches)
  foreach(filepath ${scrip_src})
    get_filename_component(BASENAME ${filepath} NAME)
    list(APPEND SOURCES_RESULT ww3/src/WW3/model/src/SCRIP/${BASENAME})
    list(APPEND BASENAME_SET ${BASENAME})
  endforeach()
endif()

if("SCRIPNC" IN_LIST switches)
  foreach(filepath ${scripnc_src})
    get_filename_component(BASENAME ${filepath} NAME)
    list(APPEND SOURCES_RESULT ww3/src/WW3/model/src/SCRIP/${BASENAME})
    list(APPEND BASENAME_SET ${BASENAME})
  endforeach()
endif()
