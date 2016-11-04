#_______________________________________________________________________
# create executable needed by DCMIP2012_test2.x (native grid output)

preqx_setup()
#              name                 target NP NC PLEV USE_PIO WITH_ENERGY QSIZE_D
createTestExec(preqx-nlev60-native  preqx  4  4  60   TRUE    FALSE       1)
