#
FILE(GLOB CPPSOURCES *.cpp)
FILE(GLOB CPPHEADERS *.hpp)

add_library(glimmercismcpp ${CPPSOURCES} ${CPPHEADERS})
include_directories ( ${GLIMMER_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}
                      ${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS})

