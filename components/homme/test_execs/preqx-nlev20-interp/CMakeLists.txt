#_______________________________________________________________________
# create executable needed by DCMIP2012_test3.0 (interpolated output)

preqx_setup()
#              name                 target NP  NC PLEV USE_PIO  WITH_ENERGY QSIZE_D
createTestExec(preqx-nlev20-interp  preqx  4   4  20   FALSE    FALSE       1       )
