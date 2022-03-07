if (DEBUG)
  string(APPEND FFLAGS " -qinitauto=7FF7FFFF -qflttrap=ov:zero:inv:en")
endif()
string(APPEND CXX_LIBS " -L/sw/ascent/gcc/8.1.1/lib64")
