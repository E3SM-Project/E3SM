include(ScreamUtils)

set(NEED_LIBS ekat scream_share)

# Add the fake test
CreateUnitTest(fake_test "fake_test.cpp" "${NEED_LIBS}")
