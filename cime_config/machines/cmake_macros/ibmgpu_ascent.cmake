if (DEBUG)
  string(APPEND FFLAGS " -qinitauto=7FF7FFFF -qflttrap=ov:zero:inv:en")
endif()
