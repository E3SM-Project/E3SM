#_______________________________________________________________________
# create executable needed by planar nonhydrostatic mountain wave (native output)

preqx_setup()
#              name                 target NP  NC PLEV USE_PIO  WITH_ENERGY QSIZE_D
createTestExec(preqx-nlev140-native  preqx  4   4  140   TRUE    FALSE       1       )
