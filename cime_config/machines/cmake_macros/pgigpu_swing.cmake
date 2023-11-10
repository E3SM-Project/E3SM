if (COMP_NAME STREQUAL gptl)
  string(APPEND CPPDEFS " -DHAVE_NANOTIME -DBIT64 -DHAVE_GETTIMEOFDAY")
endif()
string(APPEND LDFLAGS " -gpu=cc70,cc60,deepcopy -Minfo=accel")
