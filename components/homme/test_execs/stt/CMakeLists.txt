thetal_setup()

ADD_DEFINITIONS(-DHOMME_WITHOUT_PIOLIBRARY)

# Set the variables for this test executable
#                           NP  NC PLEV USE_PIO  WITH_ENERGY  QSIZE_D
createTestExec(stt theta-l 4   4   3   FALSE        TRUE     5)
