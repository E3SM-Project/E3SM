
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

https://github.com/E3SM-Project/E3SM/commit/c41728fdce93692e1480441dc4673a64f3f0ff72

The modifications mainly include:
  * The original code requires the nudging data to have only one time slice per file. The revised code can handle multiple time slices. 
  * The original code can only use the step-function nudging. The revised code can linearly interpolate the nudging data to the current model time step. 
  * Fixed bugs for the intermittent nudging configuration.

Newer version
-------------------------------------------------------------------------------- 
The nudging code is under further development. The development branch can be found here: https://github.com/E3SM-Project/E3SM/tree/jiansunpnnl/ndg_loc. For details about the new code features and improvements, please contact Kai Zhang (kai.zhang@pnnl.gov) and Jian Sun (jian.sun@pnnl.gov).

Reference
--------------------------------------------------------------------------------
Sun, J., Zhang, K., Wan, H., Ma, P.-L., Tang, Q., Zhang, S. (2019), Impact of nudging strategy on the climate representativeness and hindcast skill of constrained EAMv1 simulations, Journal of Advances in Modeling Earth Systems, doi: 10.1029/2019MS001831. https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2019MS001831

