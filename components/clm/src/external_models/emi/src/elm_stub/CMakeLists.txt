add_subdirectory(shr)
add_subdirectory(utils)
add_subdirectory(elm_subgrid_types)
add_subdirectory(elm_types_stub)

set(EMI_LIBRARIES emi_shr; emi_utils; emi_elm_subgrid_types; emi_elm_types_stub; ${EMI_LIBRARIES} PARENT_SCOPE)

