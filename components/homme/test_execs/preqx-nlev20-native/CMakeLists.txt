#_______________________________________________________________________
# create executable needed by planar hgw (native) and DCMIP2012_test3.0 (native grid output)

preqx_setup()
#              name                 target NP  NC PLEV USE_PIO  WITH_ENERGY QSIZE_D
createTestExec(preqx-nlev20-native  preqx  4   4  20   TRUE    FALSE       1       )
