
SET(SRC check_mod.F90  gdecomp_mod.F90  kinds_mod.F90  namelist_mod.F90  
            testpio.F90  utils_mod.F90)

INCLUDE_DIRECTORIES(${PIO_INCLUDE_DIRS})
LINK_DIRECTORIES(${PIO_LIB_DIR})
ADD_EXECUTABLE(testpio ${SRC})
TARGET_LINK_LIBRARIES(testpio pio)
