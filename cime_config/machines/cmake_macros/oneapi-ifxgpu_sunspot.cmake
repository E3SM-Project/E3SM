
set(CXX_LINKER "CXX")
set(BLA_VENDOR Intel10_64_dyn)
string(APPEND SLIBS " -fiopenmp -fopenmp-targets=spir64")
set(USE_SYCL "TRUE")
string(APPEND SYCL_FLAGS " -\-intel -fsycl -fsycl-targets=spir64_gen -Xsycl-target-backend \"-device 12.60.7\"")
#string(APPEND SYCL_FLAGS " -\-intel -fsycl")
string(APPEND CXX_LDFLAGS " -Wl,-\-defsym,main=MAIN_\_ -lifcore -\-intel -fsycl -lsycl -Xsycl-target-backend \"-device 12.60.7\"")
