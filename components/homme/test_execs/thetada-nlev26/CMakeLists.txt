#_______________________________________________________________________
# create executable needed by DCMIP2012_test3.0 (interpolated output)

#            DA, default is 0
thetal_setup(1)
#              name                   target   NP  NC PLEV USE_PIO  WITH_ENERGY QSIZE_D
createTestExec(thetada-nlev26         theta-l  4   4  26   FALSE    TRUE        4       )
