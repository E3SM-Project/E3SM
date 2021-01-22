#_______________________________________________________________________
# create executable needed by planar hydrostatic mountain wave (native output)

preqx_setup()
#              name                 target NP  NC PLEV USE_PIO  WITH_ENERGY QSIZE_D
createTestExec(preqx-nlev200-native  preqx  4   4  200   TRUE    FALSE       1       )
