INCLUDE (ScreamUtils)

SET (p3LibDir ${CMAKE_BINARY_DIR}/src/physics/p3)
SET (NEED_LIBS ${p3LibDir} scream_control scream_share ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

# Test atmosphere processes
CreateUnitTest(p3_stand_alone "p3_stand_alone.cpp" "${NEED_LIBS}")

