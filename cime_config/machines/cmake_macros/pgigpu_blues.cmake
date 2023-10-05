if (COMP_NAME STREQUAL gptl)
  string(APPEND CPPDEFS " -DHAVE_NANOTIME -DBIT64 -DHAVE_GETTIMEOFDAY")
endif()
string(APPEND LDFLAGS " -Minline -ta=tesla:ccall,fastmath,loadcache:L1,unroll,fma,managed,deepcopy,nonvvm -Mcuda -Minfo=accel")

