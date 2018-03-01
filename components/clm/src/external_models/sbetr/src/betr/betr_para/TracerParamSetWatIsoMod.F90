module TracerParamSetWatIsoMod
#include "bshr_assert.h"
  use bshr_kind_mod            , only : r8 => shr_kind_r8
 implicit none
  private

  public :: get_equi_lv_h2oiso_fractionation
  public :: get_equi_sv_h2oiso_fractionation
  public :: get_equi_sl_h2oiso_fractionation
  contains

   !------------------------------------------------------------------------
   function get_equi_lv_h2oiso_fractionation(trcid, temp, betrtracer_vars)result(ans)
   !
   ! DESCRIPTION
   ! get equilibrium isotopic fractionation of liquid against gaseous phase
   !
   !
   use BeTRTracerType        , only : betrtracer_type
   implicit none
   integer                   , intent(in) :: trcid
   real(r8)                  , intent(in) :: temp                !temperature
   type(betrtracer_type)     , intent(in) :: betrtracer_vars

   real(r8) :: ans

    ! remove compiler warnings for unused dummy args
    if (trcid > 0)                             continue
    if (temp > 0)                              continue
    if (len(betrtracer_vars%betr_simname) > 0) continue

   !now it is set to one, for O18 and H/D, pleasure refer to Braud et al. (2005, J. Hydrology)
   ans = 1._r8
   return
   end function get_equi_lv_h2oiso_fractionation

   !------------------------------------------------------------------------
   function get_equi_sv_h2oiso_fractionation(trcid, temp, betrtracer_vars)result(ans)
   !
   ! DESCRIPTION
   ! get equilibrium isotopic fractionation of ice against vapor phase
   !
   !
   use BeTRTracerType        , only : betrtracer_type
   implicit none
   integer                   , intent(in) :: trcid
   real(r8)                  , intent(in) :: temp                !temperature
   type(betrtracer_type)     , intent(in) :: betrtracer_vars

   real(r8) :: ans

    ! remove compiler warnings for unused dummy args
    if (trcid > 0)                             continue
    if (temp > 0)                              continue
    if (len(betrtracer_vars%betr_simname) > 0) continue

   !now it is set to one, for O18, Roche (2013, GMD) gives some information
   ans = 1._r8
   return
   end function get_equi_sv_h2oiso_fractionation


   !------------------------------------------------------------------------
   function get_equi_sl_h2oiso_fractionation(trcid, temp, betrtracer_vars)result(ans)
   !
   ! DESCRIPTION
   ! get equilibrium isotopic fractionation of ice against liquid water
   !
   !
   use BeTRTracerType        , only : betrtracer_type
   implicit none
   integer                   , intent(in) :: trcid
   real(r8)                  , intent(in) :: temp                !temperature
   type(betrtracer_type)     , intent(in) :: betrtracer_vars

   real(r8) :: ans

    ! remove compiler warnings for unused dummy args
    if (trcid > 0)                             continue
    if (temp > 0)                              continue
    if (len(betrtracer_vars%betr_simname) > 0) continue

   !now it is set to one, it is equal to alpha_sv*alpha_vl
   ans = 1._r8
   return
   end function get_equi_sl_h2oiso_fractionation

end module TracerParamSetWatIsoMod
