add_subdirectory(betr_util)
add_subdirectory(betr_math)
add_subdirectory(betr_grid)
add_subdirectory(betr_dtype)
add_subdirectory(betr_bgc)
add_subdirectory(betr_core)
add_subdirectory(betr_para)
add_subdirectory(betr_rxns)
add_subdirectory(betr_main)
# FIXME(bja, 201603) is this manual list necessary with the subdir
# adding to parent scope?

set(BETR_LIBRARIES betr_util;betr_math;betr_grid;betr_dtype;betr_bgc;betr_core;betr_para;betr_main;betr_rxns;${BETR_LIBRARIES} PARENT_SCOPE)
