INCLUDE (ScreamUtils)

SET (p3Lib p3)
SET (NEED_LIBS ${p3Lib} scream_control scream_share)

# Test atmosphere processes
CreateUnitTest(p3_stand_alone "p3_stand_alone.cpp" "${NEED_LIBS}")

