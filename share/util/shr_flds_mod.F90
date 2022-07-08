module shr_flds_mod

   use shr_kind_mod      , only : CX => shr_kind_CX, CXX => shr_kind_CXX
   use shr_sys_mod       , only : shr_sys_abort

   implicit none
   public

   !----------------------------------------------------------------------------
   ! for the domain
   !----------------------------------------------------------------------------

   character(CXX) :: shr_flds_dom_coord
   character(CXX) :: shr_flds_dom_other

end module shr_flds_mod
