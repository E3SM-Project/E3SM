
set(CXX_LINKER "CXX")
string(APPEND SLIBS " -fiopenmp -fopenmp-targets=spir64")
set(USE_SYCL "TRUE")
string(APPEND SYCL_FLAGS " -\-intel -fsycl -fsycl-targets=spir64_gen -Xsycl-target-backend \"-device 12.60.7\"") # for pvc node only
#string(APPEND SYCL_FLAGS " -\-intel -fsycl")
string(APPEND CXX_LDFLAGS " -Wl,-\-defsym,main=MAIN_\_ -lifcore -\-intel -fsycl -lsycl -Xsycl-target-backend \"-device 12.60.7\"") #for pvc node only
