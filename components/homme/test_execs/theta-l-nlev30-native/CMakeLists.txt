#_______________________________________________________________________
# create executable needed by planar baroclinic instability (native output)

thetal_setup()
#              name                 target NP  NC PLEV USE_PIO  WITH_ENERGY QSIZE_D
createTestExec(theta-l-nlev30-native         theta-l  4   4  30   TRUE    TRUE       1       )
