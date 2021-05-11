#_______________________________________________________________________
# create executable needed by DCMIP2012_test2.0 (interpolated output)

preqx_setup()
#              name                 target NP NC PLEV USE_PIO WITH_ENERGY QSIZE_D
createTestExec(preqx-nlev30-interp  preqx  4  4  30   FALSE   TRUE        6)

