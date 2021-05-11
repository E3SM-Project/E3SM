#_______________________________________________________________________
# create executable needed by DCMIP2016_test3 (interpolated output)

preqx_setup()
#              name                 target NP NC PLEV USE_PIO WITH_ENERGY QSIZE_D
createTestExec(preqx-nlev40-interp  preqx  4  4  40   FALSE   FALSE       4)

