# Remove CMake cache 
SET(CMAKE_GENERATED ${HOMME_BINARY_DIR}/CMakeCache.txt
                    ${HOMME_BINARY_DIR}/cmake_install.cmake  
                    ${HOMME_BINARY_DIR}/Makefile
                    ${HOMME_BINARY_DIR}/CMakeFiles
)

FOREACH(genFile ${CMAKE_GENERATED})

  IF (EXISTS ${genFile})
     FILE(REMOVE_RECURSE ${genFile})
  ENDIF()

ENDFOREACH(genFile)


