###-------------------------------------------------------------------------###
### CMakeList.txt for examples using pio
###-------------------------------------------------------------------------###

# Adding PIO definitions - defined in the root directory
ADD_DEFINITIONS(${PIO_DEFINITIONS})

ADD_SUBDIRECTORY(c)
ADD_SUBDIRECTORY(cxx)
ADD_SUBDIRECTORY(f03)
