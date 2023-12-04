string(APPEND LDFLAGS " -gpu=cc70,cc60 -Minfo=accel")
string(APPEND SLIBS " -L$ENV{ESSL_PATH}/lib64 -lessl")
string(APPEND KOKKOS_OPTIONS " -DKokkos_ARCH_POWER9=On -DKokkos_ARCH_VOLTA70=On -DKokkos_ENABLE_CUDA=On -DKokkos_ENABLE_CUDA_LAMBDA=On")
