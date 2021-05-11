#_______________________________________________________________________
# create executable needed by planar schar mountain wave (native output)

preqx_setup()
#              name                 target NP  NC PLEV USE_PIO  WITH_ENERGY QSIZE_D
createTestExec(preqx-nlev100-native  preqx  4   4  100   TRUE    FALSE       1       )
