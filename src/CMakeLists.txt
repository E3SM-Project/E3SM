#==============================================================================
#  PRELIMINARIES
#==============================================================================

# Test for big-endian nature 
if (PIO_TEST_BIG_ENDIAN) 
  include (TestBigEndian)
  test_big_endian (PIO_BIG_ENDIAN_TEST_RESULT)
  if (PIO_BIG_ENDIAN_TEST_RESULT) 
    set (PIO_BIG_ENDIAN ON  CACHE BOOL "Whether machine is big endian")
  else ()
    set (PIO_BIG_ENDIAN OFF CACHE BOOL "Whether machine is big endian")
  endif ()
else ()
  set (PIO_BIG_ENDIAN OFF CACHE BOOL "Whether machine is big endian")
endif ()

#==============================================================================
#  INCLUDE SOURCE DIRECTORIES
#==============================================================================

add_subdirectory (gptl)
add_subdirectory (clib)
add_subdirectory (flib)
