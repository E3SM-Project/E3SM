INCLUDE (ScreamUtils)
INCLUDE (ScreamUtils)

SET (shocLib shoc)
SET (NEED_LIBS ${shocLib} scream_control scream_share)

CreateUnitTest(shoc_tests shoc_tests.cpp "shoc;scream_share")
CreateUnitTest(shoc_stand_alone "shoc_stand_alone.cpp" "${NEED_LIBS}")
