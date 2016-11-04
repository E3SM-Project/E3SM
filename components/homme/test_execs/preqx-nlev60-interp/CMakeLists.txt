#_______________________________________________________________________
# create executable needed by DCMIP2012_test2.x (interpolated output)

preqx_setup()
#              name                 target NP NC PLEV USE_PIO WITH_ENERGY QSIZE_D
createTestExec(preqx-nlev60-interp  preqx  4  4  60   FALSE   FALSE       1)

