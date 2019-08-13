foreach(ITEM IN LISTS FILES_NEED_CUDA_FLAGS)
  set_property(SOURCE ${ITEM} APPEND_STRING PROPERTY COMPILE_FLAGS " -DUSE_OPENACC=1 -acc -ta=tesla,cc70,pinned -Minfo=accel ")
endforeach()
