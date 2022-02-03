module constMod
! Define some constants

   use shr_kind_mod, only : r8 => shr_kind_r8
   implicit none
   save

   real(R8),parameter :: SHR_CONST_REARTH  = 6.37122e6_R8   ! radius of earth ~ m
   real(r8),parameter :: re_km = SHR_CONST_REARTH*0.001        ! radius of earth (km)

end module constMod
