# Remove CMake cache 
SET(CMAKE_GENERATED ${CMAKE_BINARY_DIR}/CMakeCache.txt
                    ${CMAKE_BINARY_DIR}/cmake_install.cmake  
                    ${CMAKE_BINARY_DIR}/Makefile
                    ${CMAKE_BINARY_DIR}/CMakeFiles
)

FOREACH(genFile ${CMAKE_GENERATED})

  IF (EXISTS ${genFile})
     FILE(REMOVE_RECURSE ${genFile})
  ENDIF()

ENDFOREACH(genFile)


