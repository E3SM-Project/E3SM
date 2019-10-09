
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

https://github.com/E3SM-Project/E3SM/commit/1317ef4b411bf2c1c4363382ac0ac2c0ec5a5000

The modifications mainly include:
  * The original code requires the nudging data to have only one time slice per file. The revised code can handle multiple time slices. 
  * The original code can only use the step-function nudging. The revised code can linearly interpolate the nudging data to the current model time step. 
  * Fixed bugs for the intermittent nudging configuration. 


Reference
--------------------------------------------------------------------------------
Sun et al. (2019) Impact of nudging strategy on the climate representativeness and hindcast skill of constrained EAMv1 simulations. Under review for JAMES. 

