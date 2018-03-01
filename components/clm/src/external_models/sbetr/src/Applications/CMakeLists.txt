add_subdirectory(soil-farm)
add_subdirectory(app_util)
# FIXME(bja, 201603) is this manual list necessary with the subdir
# adding to parent scope?

set(BETR_LIBRARIES app_util;${BETR_LIBRARIES} PARENT_SCOPE)
