# This macro sets compiler specific flags
macro(set_up_compilers)
  #
  # General C compiler flags.
  #
  if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -pedantic-errors -Wextra")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Werror-implicit-function-declaration")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-sign-compare -Wno-unused-parameter")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-unused-but-set-variable")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-int-to-pointer-cast -Wno-pointer-to-int-cast")

    if (BUILD_SHARED_LIBS)
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fPIC")
    endif()

    if (LINUX EQUAL 1)
      # Counter some of GCC's more recent stinginess on Linux.
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_POSIX_C_SOURCE=200809L")
    endif()

  elseif (CMAKE_C_COMPILER_ID STREQUAL "Clang")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -pedantic-errors -Wextra")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Werror-implicit-function-declaration ")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-builtin")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-sign-compare -Wno-unused-parameter ")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-int-to-pointer-cast -Wno-pointer-to-int-cast")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-unused-function")

  elseif (CMAKE_C_COMPILER_ID STREQUAL "Intel")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall")
  endif()

  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${SYS_FLAGS}")

  #
  # Fortran compiler flags.
  #
  if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -W -Wall -std=gnu -pedantic")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DCPRGNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wno-unused-variable")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wno-unused-parameter")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Werror=use-without-only -ffree-line-length-none")

  elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DCPRINTEL")
    if(HOSTNAME MATCHES "edison")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mkl")
    elseif (HOSTNAME MATCHES "cori")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mkl")
    endif()
  endif()

endmacro()
