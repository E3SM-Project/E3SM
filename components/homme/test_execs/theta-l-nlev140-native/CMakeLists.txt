#_______________________________________________________________________
# create executable needed by planar nonhydrostatic mountain wave (native output)

thetal_setup()
#              name                 target NP  NC PLEV USE_PIO  WITH_ENERGY QSIZE_D
createTestExec(theta-l-nlev140-native         theta-l  4   4  140   TRUE    TRUE       1       )
