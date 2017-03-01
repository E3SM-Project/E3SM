###-------------------------------------------------------------------------###
### CMakeList.txt for examples using pio
###-------------------------------------------------------------------------###

# Include PIO include and lib directories
INCLUDE_DIRECTORIES(${PIO_INCLUDE_DIRS})
LINK_DIRECTORIES(${PIO_LIB_DIR})

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -O0")
# Compiler-specific compiler options
if ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c99")
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "PGI")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -c99")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -c99")
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c99")
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c99")
endif()
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O0")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0")

SET(SRC examplePio.c)
ADD_EXECUTABLE(examplePio_c ${SRC})
TARGET_LINK_LIBRARIES(examplePio_c pioc)

SET(SRC example1.c)
ADD_EXECUTABLE(example1 ${SRC})
TARGET_LINK_LIBRARIES(example1 pioc)

SET(SRC example2.c)
ADD_EXECUTABLE(example2 ${SRC})
TARGET_LINK_LIBRARIES(example2 pioc)

#===== MPE =====
find_package (MPE "2.4.8" COMPONENTS C)
if (MPE_C_FOUND)
    target_include_directories (example2
        PUBLIC ${MPE_C_INCLUDE_DIRS})
    target_compile_definitions (example2
        PUBLIC HAVE_MPE)
    target_link_libraries (example2 pioc lmpe mpe pthread m)
else ()
    target_compile_definitions (example2
        PUBLIC _NOMPE)
endif ()

