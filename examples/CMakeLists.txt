###-------------------------------------------------------------------------###
### CMakeList.txt for examples using pio
###-------------------------------------------------------------------------###

if (${PLATFORM} STREQUAL "goldbach" )
    ADD_SUBDIRECTORY(f03)
else()
    ADD_SUBDIRECTORY(f03)
    ADD_SUBDIRECTORY(c)
    ADD_SUBDIRECTORY(cxx)
endif()