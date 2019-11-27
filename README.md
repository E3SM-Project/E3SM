
Revised nudging code for E3SMv1 
================================================================================

This is a branch of the E3SM model repository (https://github.com/E3SM-Project/E3SM/). It contains the new nudging code for E3SMv1 (Sun et al., 2019). 

Nudging code 
--------------------------------------------------------------------------------
The original nudging code can be found at: 

https://github.com/E3SM-Project/E3SM/blob/master/components/cam/src/physics/cam/nudging.F90

The revised nudging code can be found at: 

https://github.com/E3SM-Project/E3SM/blob/jiansunpnnl/ms2019/components/cam/src/physics/cam/nudging.F90

Code modifications
--------------------------------------------------------------------------------
Code modifications can be viewed at: 

https://github.com/E3SM-Project/E3SM/commit/c41728fdce93692e1480441dc4673a64f3f0ff72

The modifications mainly include:
  * The original code requires the nudging data to have only one time slice per file. The revised code can handle multiple time slices. 
  * The original code can only use the step-function nudging. The revised code can linearly interpolate the nudging data to the current model time step. 
  * Fixed bugs for the intermittent nudging configuration.

Further development
-------------------------------------------------------------------------------- 
The nudging code can be further improved as:
  * The modified nudging code is only tested for SE dycore, not for FV dycore.
  * The modified nudging code works as expected for start-up run, but not for restart/branch run.
  * The nudging tendency is calculated at the location different from where the nudging data are output in the CLIM simulation.

Such kind of improvements are available at:
https://github.com/E3SM-Project/E3SM/tree/jiansunpnnl/ndg_loc

Reference
--------------------------------------------------------------------------------
Sun, J., Zhang, K., Wan, H., Ma, P.-L., Tang, Q., Zhang, S. (2019), Impact of nudging strategy on the climate representativeness and hindcast skill of constrained EAMv1 simulations, journal of advances in modeling earth systems, doi: 10.1029/2019MS001831.
