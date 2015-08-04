if (ENABLE_TIMING)
    add_subdirectory(gptl EXCLUDE_FROM_ALL)
endif ()

add_subdirectory(c)

if (ENABLE_FORTRAN)
    add_subdirectory(fortran)
endif ()
