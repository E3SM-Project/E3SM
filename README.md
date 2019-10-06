
Revised nudging code for E3SMv1 
================================================================================

This is a branch of the E3SM model repository: 

https://github.com/E3SM-Project/E3SM/

The original nudging code can be found at: 

https://github.com/E3SM-Project/E3SM/blob/master/components/cam/src/physics/cam/nudging.F90

The revised nudging code can be found at: 

https://github.com/E3SM-Project/E3SM/blob/jiansunpnnl/ms2019/components/cam/src/physics/cam/nudging.F90

The modifications mainly include:
  * The original nudging code only works appropriately when there is one time slice per data file. The revised code can work well with multiple time slices per data file.
  * The original nudging code only applies the step-function nudging. The revised code linearly interpolates the nudging data to each model time step. The difference can be viewed in Sun et al.'s JAMES paper.
  * Fix bugs to perform a real intermittent nudged simulation.


Reference
--------------------------------------------------------------------------------
Sun et al. (2019) Impact of nudging strategy on the climate representativeness and hindcast skill of constrained EAMv1 simulations. Under review for JAMES. 

