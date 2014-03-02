
SET(SRC check_mod.F90  gdecomp_mod.F90  kinds_mod.F90  namelist_mod.F90  
            testpio.F90  utils_mod.F90)
SET(WSSRC wstest.c)

INCLUDE_DIRECTORIES(${PIO_INCLUDE_DIRS})
LINK_DIRECTORIES(${PIO_LIB_DIR})
ADD_EXECUTABLE(testpio ${SRC})
ADD_EXECUTABLE(wstest ${WSSRC})
TARGET_LINK_LIBRARIES(testpio piof pioc)
TARGET_LINK_LIBRARIES(wstest pioc)


