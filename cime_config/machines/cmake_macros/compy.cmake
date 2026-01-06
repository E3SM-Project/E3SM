string(APPEND KOKKOS_OPTIONS " -DKokkos_ENABLE_OPENMP=On -DKokkos_ENABLE_AGGRESSIVE_VECTORIZATION=On -DKokkos_ENABLE_EXPLICIT_INSTANTIATION=Off -DKokkos_ENABLE_DEPRECATED_CODE_4=Off")
string(APPEND CMAKE_C_FLAGS " -I/share/apps/python/anaconda3-2020.02/include/python3.7m")
string(APPEND SLIBS " /share/apps/python/anaconda3-2020.02/lib/libpython3.7m.a -lutil -lpthread -ldl -lm")
