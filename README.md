
Revised nudging code for E3SMv1 
================================================================================

This is a branch created based on the branch jiansunpnnl/ms2019, which is a branch of the E3SM model repository (https://github.com/E3SM-Project/E3SM/) and contains the new nudging code for E3SMv1 (Sun et al., 2019). 

Nudging code 
--------------------------------------------------------------------------------
The original nudging code can be found at: 

https://github.com/E3SM-Project/E3SM/blob/master/components/cam/src/physics/cam/nudging.F90

The nudging code used for Sun et al. (2019) can be found at:
https://github.com/E3SM-Project/E3SM/blob/jiansunpnnl/ms2019/components/cam/src/physics/cam/nudging.F90

The revised nudging code can be found at: 

https://github.com/E3SM-Project/E3SM/blob/jiansunpnnl/ndg_loc/components/cam/src/physics/cam/nudging.F90

Note that the physpkg.F90 is also changed and can be found at:
https://github.com/E3SM-Project/E3SM/blob/jiansunpnnl/ndg_loc/components/cam/src/physics/cam/physpkg.F90

Code modifications
--------------------------------------------------------------------------------
Code modifications can be viewed at: 

https://github.com/E3SM-Project/E3SM/commit/267911c6d9deda95b81b1fc8850d880a8d3f6bb5

The modifications against the version used in Sun et al. (2019) mainly include:
  * Fixed the issue for restart/branch run.
  * The nudging tendency can be optionally calculated at the same location where the nudging data are output in the CLIM simulation.
  * The linear interpolation can interpolate the nudging data to the current or future (default) model time step.
  * Updated nudging code for FV dycore. Note that the location to apply the nudging tendency is changed in the physpkg.F90, because the implementation in the original nudging code is problematic. Such a change does not affect the simulations with the SE dycore but has a notable impact on those with the FV dycore. 
  * Add a namelist variable to explicitly specify the number of time slices per file.

Reference
--------------------------------------------------------------------------------
Sun, J., Zhang, K., Wan, H., Ma, P.-L., Tang, Q., Zhang, S. (2019), Impact of nudging strategy on the climate representativeness and hindcast skill of constrained EAMv1 simulations, Journal of Advances in Modeling Earth Systems, doi: 10.1029/2019MS001831.
