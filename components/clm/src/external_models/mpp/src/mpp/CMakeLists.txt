add_subdirectory(util)
add_subdirectory(dtypes)
add_subdirectory(vsfm)
add_subdirectory(thermal)
add_subdirectory(thermal-e)
add_subdirectory(th)

set(MPP_LIBRARIES mpp_dtypes; mpp_thermal; mpp_thermal_e; mpp_util; mpp_vsfm; mpp_th; ${MPP_LIBRARIES} PARENT_SCOPE)

