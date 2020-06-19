foreach(ITEM IN LISTS FILES_NEED_CUDA_FLAGS)
  e3sm_add_flags("${ITEM}" "-DUSE_OPENACC=1 -acc -ta=tesla,cc70,pinned -Minfo=accel")
endforeach()
