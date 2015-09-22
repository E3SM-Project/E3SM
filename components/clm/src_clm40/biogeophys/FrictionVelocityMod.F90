module FrictionVelocityMod

!------------------------------------------------------------------------------
!BOP
!
! !MODULE: FrictionVelocityMod
!
! !DESCRIPTION:
! Calculation of the friction velocity, relation for potential
! temperature and humidity profiles of surface boundary layer.
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: FrictionVelocity       ! Calculate friction velocity
  public :: MoninObukIni           ! Initialization of the Monin-Obukhov length
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: StabilityFunc1        ! Stability function for rib < 0.
  private :: StabilityFunc2        ! Stability function for rib < 0.
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------------

contains

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: FrictionVelocity
!
! !INTERFACE:
  subroutine FrictionVelocity(lbn, ubn, fn, filtern, &
                              displa, z0m, z0h, z0q, &
                              obu, iter, ur, um, ustar, &
                              temp1, temp2, temp12m, temp22m, fm, landunit_index)
!
! !DESCRIPTION:
! Calculation of the friction velocity, relation for potential
! temperature and humidity profiles of surface boundary layer.
! The scheme is based on the work of Zeng et al. (1998):
! Intercomparison of bulk aerodynamic algorithms for the computation
! of sea surface fluxes using TOGA CORE and TAO data. J. Climate,
! Vol. 11, 2628-2644.
!
! !USES:
   use clmtype
   use clm_atmlnd, only : clm_a2l
   use clm_varcon, only : vkc
   use clm_varctl, only : iulog
!
! !ARGUMENTS:
   implicit none
   integer , intent(in)  :: lbn, ubn         ! pft/landunit array bounds
   integer , intent(in)  :: fn               ! number of filtered pft/landunit elements
   integer , intent(in)  :: filtern(fn)      ! pft/landunit filter
   real(r8), intent(in)  :: displa(lbn:ubn)  ! displacement height (m)
   real(r8), intent(in)  :: z0m(lbn:ubn)     ! roughness length over vegetation, momentum [m]
   real(r8), intent(in)  :: z0h(lbn:ubn)     ! roughness length over vegetation, sensible heat [m]
   real(r8), intent(in)  :: z0q(lbn:ubn)     ! roughness length over vegetation, latent heat [m]
   real(r8), intent(in)  :: obu(lbn:ubn)     ! monin-obukhov length (m)
   integer,  intent(in)  :: iter             ! iteration number
   real(r8), intent(in)  :: ur(lbn:ubn)      ! wind speed at reference height [m/s]
   real(r8), intent(in)  :: um(lbn:ubn)      ! wind speed including the stablity effect [m/s]
   logical,  optional, intent(in)  :: landunit_index  ! optional argument that defines landunit or pft level
   real(r8), intent(out) :: ustar(lbn:ubn)   ! friction velocity [m/s]
   real(r8), intent(out) :: temp1(lbn:ubn)   ! relation for potential temperature profile
   real(r8), intent(out) :: temp12m(lbn:ubn) ! relation for potential temperature profile applied at 2-m
   real(r8), intent(out) :: temp2(lbn:ubn)   ! relation for specific humidity profile
   real(r8), intent(out) :: temp22m(lbn:ubn) ! relation for specific humidity profile applied at 2-m
   real(r8), intent(inout) :: fm(lbn:ubn)    ! diagnose 10m wind (DUST only)
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
! 12/19/01, Peter Thornton
! Added arguments to eliminate passing clm derived type into this function.
! Created by Mariana Vertenstein
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
   integer , pointer :: ngridcell(:)      !pft/landunit gridcell index
   real(r8), pointer :: forc_hgt_u_pft(:) !observational height of wind at pft level [m]
   real(r8), pointer :: forc_hgt_t_pft(:) !observational height of temperature at pft level [m]
   real(r8), pointer :: forc_hgt_q_pft(:) !observational height of specific humidity at pft level [m]
   integer , pointer :: pfti(:)           !beginning pfti index for landunit
   integer , pointer :: pftf(:)           !final pft index for landunit
!
! local pointers to implicit out arguments
!
   real(r8), pointer :: u10(:)         ! 10-m wind (m/s) (for dust model)
   real(r8), pointer :: fv(:)          ! friction velocity (m/s) (for dust model)
   real(r8), pointer :: vds(:)         ! dry deposition velocity term (m/s) (for SO4 NH4NO3)
   real(r8), pointer :: u10_clm(:)     ! 10-m wind (m/s)
   real(r8), pointer :: va(:)          ! atmospheric wind speed plus convective velocity (m/s)
!
!
! !OTHER LOCAL VARIABLES:
!EOP
!
   real(r8), parameter :: zetam = 1.574_r8 ! transition point of flux-gradient relation (wind profile)
   real(r8), parameter :: zetat = 0.465_r8 ! transition point of flux-gradient relation (temp. profile)
   integer :: f                            ! pft/landunit filter index
   integer :: n                            ! pft/landunit index
   integer :: g                            ! gridcell index
   integer :: pp                           ! pfti,pftf index
   real(r8):: zldis(lbn:ubn)               ! reference height "minus" zero displacement heght [m]
   real(r8):: zeta(lbn:ubn)                ! dimensionless height used in Monin-Obukhov theory
   real(r8) :: tmp1,tmp2,tmp3,tmp4         ! Used to diagnose the 10 meter wind
   real(r8) :: fmnew                       ! Used to diagnose the 10 meter wind
   real(r8) :: fm10                        ! Used to diagnose the 10 meter wind
   real(r8) :: zeta10                      ! Used to diagnose the 10 meter wind
   real(r8) :: vds_tmp                     ! Temporary for dry deposition velocity
!------------------------------------------------------------------------------

   ! Assign local pointers to derived type members (gridcell-level)

   if (present(landunit_index)) then
     ngridcell  => lun%gridcell
   else
     ngridcell  => pft%gridcell
   end if

   vds        => pps%vds
   u10        => pps%u10
   u10_clm    => pps%u10_clm
   va         => pps%va
   fv         => pps%fv

   ! Assign local pointers to derived type members (pft or landunit-level)

   pfti             => lun%pfti
   pftf             => lun%pftf

   ! Assign local pointers to derived type members (pft-level)

   forc_hgt_u_pft => pps%forc_hgt_u_pft
   forc_hgt_t_pft => pps%forc_hgt_t_pft
   forc_hgt_q_pft => pps%forc_hgt_q_pft

   ! Adjustment factors for unstable (moz < 0) or stable (moz > 0) conditions.

   do f = 1, fn
      n = filtern(f)
      g = ngridcell(n)

      ! Wind profile

      if (present(landunit_index)) then
        zldis(n) = forc_hgt_u_pft(pfti(n))-displa(n)
      else
        zldis(n) = forc_hgt_u_pft(n)-displa(n)
      end if
      zeta(n) = zldis(n)/obu(n)
      if (zeta(n) < -zetam) then
         ustar(n) = vkc*um(n)/(log(-zetam*obu(n)/z0m(n))&
              - StabilityFunc1(-zetam) &
              + StabilityFunc1(z0m(n)/obu(n)) &
              + 1.14_r8*((-zeta(n))**0.333_r8-(zetam)**0.333_r8))
      else if (zeta(n) < 0._r8) then
         ustar(n) = vkc*um(n)/(log(zldis(n)/z0m(n))&
              - StabilityFunc1(zeta(n))&
              + StabilityFunc1(z0m(n)/obu(n)))
      else if (zeta(n) <=  1._r8) then
         ustar(n) = vkc*um(n)/(log(zldis(n)/z0m(n)) + 5._r8*zeta(n) -5._r8*z0m(n)/obu(n))
      else
         ustar(n) = vkc*um(n)/(log(obu(n)/z0m(n))+5._r8-5._r8*z0m(n)/obu(n) &
              +(5._r8*log(zeta(n))+zeta(n)-1._r8))
      end if
      
      if (zeta(n) < 0._r8) then
         vds_tmp = 2.e-3_r8*ustar(n) * ( 1._r8 + (300._r8/(-obu(n)))**0.666_r8)
      else
         vds_tmp = 2.e-3_r8*ustar(n)
      endif

      if (present(landunit_index)) then
         do pp = pfti(n),pftf(n)
            vds(pp) = vds_tmp
         end do
      else
         vds(n) = vds_tmp
      end if

! Calculate a 10-m wind (10m + z0m + d)
! For now, this will not be the same as the 10-m wind calculated for the dust 
! model because the CLM stability functions are used here, not the LSM stability
! functions used in the dust model. We will eventually change the dust model to be 
! consistent with the following formulation.
! Note that the 10-m wind calculated this way could actually be larger than the
! atmospheric forcing wind because 1) this includes the convective velocity, 2)
! this includes the 1 m/s minimum wind threshold

! If forcing height is less than or equal to 10m, then set 10-m wind to um
      if (present(landunit_index)) then
        do pp = pfti(n),pftf(n)
          if (zldis(n)-z0m(n) .le. 10._r8) then
            u10_clm(pp) = um(n)
          else
            if (zeta(n) < -zetam) then
              u10_clm(pp) = um(n) - ( ustar(n)/vkc*(log(-zetam*obu(n)/(10._r8+z0m(n)))      &
                                      - StabilityFunc1(-zetam)                              &
                                      + StabilityFunc1((10._r8+z0m(n))/obu(n))              &
                                      + 1.14_r8*((-zeta(n))**0.333_r8-(zetam)**0.333_r8)) )
            else if (zeta(n) < 0._r8) then
              u10_clm(pp) = um(n) - ( ustar(n)/vkc*(log(zldis(n)/(10._r8+z0m(n)))           &
                                      - StabilityFunc1(zeta(n))                             &
                                      + StabilityFunc1((10._r8+z0m(n))/obu(n))) )
            else if (zeta(n) <=  1._r8) then
              u10_clm(pp) = um(n) - ( ustar(n)/vkc*(log(zldis(n)/(10._r8+z0m(n)))           &
                                      + 5._r8*zeta(n) - 5._r8*(10._r8+z0m(n))/obu(n)) )
            else
              u10_clm(pp) = um(n) - ( ustar(n)/vkc*(log(obu(n)/(10._r8+z0m(n)))             &
                                      + 5._r8 - 5._r8*(10._r8+z0m(n))/obu(n)                &
                                      + (5._r8*log(zeta(n))+zeta(n)-1._r8)) )

            end if
          end if
          va(pp) = um(n)
        end do
      else
        if (zldis(n)-z0m(n) .le. 10._r8) then
          u10_clm(n) = um(n)
        else
          if (zeta(n) < -zetam) then
            u10_clm(n) = um(n) - ( ustar(n)/vkc*(log(-zetam*obu(n)/(10._r8+z0m(n)))         &
                                   - StabilityFunc1(-zetam)                                 &
                                   + StabilityFunc1((10._r8+z0m(n))/obu(n))                 &
                                   + 1.14_r8*((-zeta(n))**0.333_r8-(zetam)**0.333_r8)) )
          else if (zeta(n) < 0._r8) then
            u10_clm(n) = um(n) - ( ustar(n)/vkc*(log(zldis(n)/(10._r8+z0m(n)))              &
                                   - StabilityFunc1(zeta(n))                                &
                                   + StabilityFunc1((10._r8+z0m(n))/obu(n))) )
          else if (zeta(n) <=  1._r8) then
            u10_clm(n) = um(n) - ( ustar(n)/vkc*(log(zldis(n)/(10._r8+z0m(n)))              &
                                   + 5._r8*zeta(n) - 5._r8*(10._r8+z0m(n))/obu(n)) )
          else
            u10_clm(n) = um(n) - ( ustar(n)/vkc*(log(obu(n)/(10._r8+z0m(n)))                &
                                   + 5._r8 - 5._r8*(10._r8+z0m(n))/obu(n)                   &
                                   + (5._r8*log(zeta(n))+zeta(n)-1._r8)) )
          end if
        end if
        va(n) = um(n)
      end if

      ! Temperature profile

      if (present(landunit_index)) then
        zldis(n) = forc_hgt_t_pft(pfti(n))-displa(n)
      else
        zldis(n) = forc_hgt_t_pft(n)-displa(n)
      end if
      zeta(n) = zldis(n)/obu(n)
      if (zeta(n) < -zetat) then
         temp1(n) = vkc/(log(-zetat*obu(n)/z0h(n))&
              - StabilityFunc2(-zetat) &
              + StabilityFunc2(z0h(n)/obu(n)) &
              + 0.8_r8*((zetat)**(-0.333_r8)-(-zeta(n))**(-0.333_r8)))
      else if (zeta(n) < 0._r8) then
         temp1(n) = vkc/(log(zldis(n)/z0h(n)) &
              - StabilityFunc2(zeta(n)) &
              + StabilityFunc2(z0h(n)/obu(n)))
      else if (zeta(n) <=  1._r8) then
         temp1(n) = vkc/(log(zldis(n)/z0h(n)) + 5._r8*zeta(n) - 5._r8*z0h(n)/obu(n))
      else
         temp1(n) = vkc/(log(obu(n)/z0h(n)) + 5._r8 - 5._r8*z0h(n)/obu(n) &
              + (5._r8*log(zeta(n))+zeta(n)-1._r8))
      end if

      ! Humidity profile

      if (present(landunit_index)) then
       if (forc_hgt_q_pft(pfti(n)) == forc_hgt_t_pft(pfti(n)) .and. z0q(n) == z0h(n)) then
         temp2(n) = temp1(n)
       else
         zldis(n) = forc_hgt_q_pft(pfti(n))-displa(n)
         zeta(n) = zldis(n)/obu(n)
         if (zeta(n) < -zetat) then
            temp2(n) = vkc/(log(-zetat*obu(n)/z0q(n)) &
                 - StabilityFunc2(-zetat) &
                 + StabilityFunc2(z0q(n)/obu(n)) &
                 + 0.8_r8*((zetat)**(-0.333_r8)-(-zeta(n))**(-0.333_r8)))
         else if (zeta(n) < 0._r8) then
            temp2(n) = vkc/(log(zldis(n)/z0q(n)) &
                 - StabilityFunc2(zeta(n)) &
                 + StabilityFunc2(z0q(n)/obu(n)))
         else if (zeta(n) <=  1._r8) then
            temp2(n) = vkc/(log(zldis(n)/z0q(n)) + 5._r8*zeta(n)-5._r8*z0q(n)/obu(n))
         else
            temp2(n) = vkc/(log(obu(n)/z0q(n)) + 5._r8 - 5._r8*z0q(n)/obu(n) &
                 + (5._r8*log(zeta(n))+zeta(n)-1._r8))
         end if
       end if
      else
       if (forc_hgt_q_pft(n) == forc_hgt_t_pft(n) .and. z0q(n) == z0h(n)) then
         temp2(n) = temp1(n)
       else
         zldis(n) = forc_hgt_q_pft(n)-displa(n)
         zeta(n) = zldis(n)/obu(n)
         if (zeta(n) < -zetat) then
            temp2(n) = vkc/(log(-zetat*obu(n)/z0q(n)) &
                 - StabilityFunc2(-zetat) &
                 + StabilityFunc2(z0q(n)/obu(n)) &
                 + 0.8_r8*((zetat)**(-0.333_r8)-(-zeta(n))**(-0.333_r8)))
         else if (zeta(n) < 0._r8) then
            temp2(n) = vkc/(log(zldis(n)/z0q(n)) &
                 - StabilityFunc2(zeta(n)) &
                 + StabilityFunc2(z0q(n)/obu(n)))
         else if (zeta(n) <=  1._r8) then
            temp2(n) = vkc/(log(zldis(n)/z0q(n)) + 5._r8*zeta(n)-5._r8*z0q(n)/obu(n))
         else
            temp2(n) = vkc/(log(obu(n)/z0q(n)) + 5._r8 - 5._r8*z0q(n)/obu(n) &
                 + (5._r8*log(zeta(n))+zeta(n)-1._r8))
         end if
       endif
      endif

      ! Temperature profile applied at 2-m

      zldis(n) = 2.0_r8 + z0h(n)
      zeta(n) = zldis(n)/obu(n)
      if (zeta(n) < -zetat) then
         temp12m(n) = vkc/(log(-zetat*obu(n)/z0h(n))&
              - StabilityFunc2(-zetat) &
              + StabilityFunc2(z0h(n)/obu(n)) &
              + 0.8_r8*((zetat)**(-0.333_r8)-(-zeta(n))**(-0.333_r8)))
      else if (zeta(n) < 0._r8) then
         temp12m(n) = vkc/(log(zldis(n)/z0h(n)) &
              - StabilityFunc2(zeta(n))  &
              + StabilityFunc2(z0h(n)/obu(n)))
      else if (zeta(n) <=  1._r8) then
         temp12m(n) = vkc/(log(zldis(n)/z0h(n)) + 5._r8*zeta(n) - 5._r8*z0h(n)/obu(n))
      else
         temp12m(n) = vkc/(log(obu(n)/z0h(n)) + 5._r8 - 5._r8*z0h(n)/obu(n) &
              + (5._r8*log(zeta(n))+zeta(n)-1._r8))
      end if

      ! Humidity profile applied at 2-m

      if (z0q(n) == z0h(n)) then
         temp22m(n) = temp12m(n)
      else
         zldis(n) = 2.0_r8 + z0q(n)
         zeta(n) = zldis(n)/obu(n)
         if (zeta(n) < -zetat) then
            temp22m(n) = vkc/(log(-zetat*obu(n)/z0q(n)) - &
                 StabilityFunc2(-zetat) + StabilityFunc2(z0q(n)/obu(n)) &
                 + 0.8_r8*((zetat)**(-0.333_r8)-(-zeta(n))**(-0.333_r8)))
         else if (zeta(n) < 0._r8) then
            temp22m(n) = vkc/(log(zldis(n)/z0q(n)) - &
                 StabilityFunc2(zeta(n))+StabilityFunc2(z0q(n)/obu(n)))
         else if (zeta(n) <=  1._r8) then
            temp22m(n) = vkc/(log(zldis(n)/z0q(n)) + 5._r8*zeta(n)-5._r8*z0q(n)/obu(n))
         else
            temp22m(n) = vkc/(log(obu(n)/z0q(n)) + 5._r8 - 5._r8*z0q(n)/obu(n) &
                 + (5._r8*log(zeta(n))+zeta(n)-1._r8))
         end if
      end if

      ! diagnose 10-m wind for dust model (dstmbl.F)
      ! Notes from C. Zender's dst.F:
      ! According to Bon96 p. 62, the displacement height d (here displa) is
      ! 0.0 <= d <= 0.34 m in dust source regions (i.e., regions w/o trees).
      ! Therefore d <= 0.034*z1 and may safely be neglected.
      ! Code from LSM routine SurfaceTemperature was used to obtain u10

      if (present(landunit_index)) then
        zldis(n) = forc_hgt_u_pft(pfti(n))-displa(n)
      else
        zldis(n) = forc_hgt_u_pft(n)-displa(n)
      end if
      zeta(n) = zldis(n)/obu(n)
      if (min(zeta(n), 1._r8) < 0._r8) then
         tmp1 = (1._r8 - 16._r8*min(zeta(n),1._r8))**0.25_r8
         tmp2 = log((1._r8+tmp1*tmp1)/2._r8)
         tmp3 = log((1._r8+tmp1)/2._r8)
         fmnew = 2._r8*tmp3 + tmp2 - 2._r8*atan(tmp1) + 1.5707963_r8
      else
         fmnew = -5._r8*min(zeta(n),1._r8)
      endif
      if (iter == 1) then
         fm(n) = fmnew
      else
         fm(n) = 0.5_r8 * (fm(n)+fmnew)
      end if
      zeta10 = min(10._r8/obu(n), 1._r8)
      if (zeta(n) == 0._r8) zeta10 = 0._r8
      if (zeta10 < 0._r8) then
         tmp1 = (1.0_r8 - 16.0_r8 * zeta10)**0.25_r8
         tmp2 = log((1.0_r8 + tmp1*tmp1)/2.0_r8)
         tmp3 = log((1.0_r8 + tmp1)/2.0_r8)
         fm10 = 2.0_r8*tmp3 + tmp2 - 2.0_r8*atan(tmp1) + 1.5707963_r8
      else                ! not stable
         fm10 = -5.0_r8 * zeta10
      end if
      if (present(landunit_index)) then
        tmp4 = log( max( 1.0_8, forc_hgt_u_pft(pfti(n)) / 10._r8) )
      else 
        tmp4 = log( max( 1.0_8, forc_hgt_u_pft(n) / 10._r8) )
      end if
      if (present(landunit_index)) then
        do pp = pfti(n),pftf(n)
          u10(pp) = ur(n) - ustar(n)/vkc * (tmp4 - fm(n) + fm10)
          fv(pp)  = ustar(n)
        end do 
      else
        u10(n) = ur(n) - ustar(n)/vkc * (tmp4 - fm(n) + fm10)
        fv(n)  = ustar(n)
      end if

   end do

   end subroutine FrictionVelocity

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: StabilityFunc
!
! !INTERFACE:
   real(r8) function StabilityFunc1(zeta)
!
! !DESCRIPTION:
! Stability function for rib < 0.
!
! !USES:
      use shr_const_mod, only: SHR_CONST_PI
!
! !ARGUMENTS:
      implicit none
      real(r8), intent(in) :: zeta  ! dimensionless height used in Monin-Obukhov theory
!
! !CALLED FROM:
! subroutine FrictionVelocity in this module
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
!
!
! !LOCAL VARIABLES:
!EOP
      real(r8) :: chik, chik2
!------------------------------------------------------------------------------

      chik2 = sqrt(1._r8-16._r8*zeta)
      chik = sqrt(chik2)
      StabilityFunc1 = 2._r8*log((1._r8+chik)*0.5_r8) &
           + log((1._r8+chik2)*0.5_r8)-2._r8*atan(chik)+SHR_CONST_PI*0.5_r8

    end function StabilityFunc1

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: StabilityFunc2
!
! !INTERFACE:
   real(r8) function StabilityFunc2(zeta)
!
! !DESCRIPTION:
! Stability function for rib < 0.
!
! !USES:
     use shr_const_mod, only: SHR_CONST_PI
!
! !ARGUMENTS:
     implicit none
     real(r8), intent(in) :: zeta  ! dimensionless height used in Monin-Obukhov theory
!
! !CALLED FROM:
! subroutine FrictionVelocity in this module
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
!
!
! !LOCAL VARIABLES:
!EOP
     real(r8) :: chik2
!------------------------------------------------------------------------------

     chik2 = sqrt(1._r8-16._r8*zeta)
     StabilityFunc2 = 2._r8*log((1._r8+chik2)*0.5_r8)

   end function StabilityFunc2

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: MoninObukIni
!
! !INTERFACE:
  subroutine MoninObukIni (ur, thv, dthv, zldis, z0m, um, obu)
!
! !DESCRIPTION:
! Initialization of the Monin-Obukhov length.
! The scheme is based on the work of Zeng et al. (1998):
! Intercomparison of bulk aerodynamic algorithms for the computation
! of sea surface fluxes using TOGA CORE and TAO data. J. Climate,
! Vol. 11, 2628-2644.
!
! !USES:
    use clm_varcon, only : grav
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(in)  :: ur    ! wind speed at reference height [m/s]
    real(r8), intent(in)  :: thv   ! virtual potential temperature (kelvin)
    real(r8), intent(in)  :: dthv  ! diff of vir. poten. temp. between ref. height and surface
    real(r8), intent(in)  :: zldis ! reference height "minus" zero displacement heght [m]
    real(r8), intent(in)  :: z0m   ! roughness length, momentum [m]
    real(r8), intent(out) :: um    ! wind speed including the stability effect [m/s]
    real(r8), intent(out) :: obu   ! monin-obukhov length (m)
!
! !CALLED FROM:
! subroutine BareGroundFluxes in module BareGroundFluxesMod.F90
! subroutine BiogeophysicsLake in module BiogeophysicsLakeMod.F90
! subroutine CanopyFluxes in module CanopyFluxesMod.F90
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
!
!
! !LOCAL VARIABLES:
!EOP
!
    real(r8) :: wc    ! convective velocity [m/s]
    real(r8) :: rib   ! bulk Richardson number
    real(r8) :: zeta  ! dimensionless height used in Monin-Obukhov theory
    real(r8) :: ustar ! friction velocity [m/s]
!-----------------------------------------------------------------------

    ! Initial values of u* and convective velocity

    ustar=0.06_r8
    wc=0.5_r8
    if (dthv >= 0._r8) then
       um=max(ur,0.1_r8)
    else
       um=sqrt(ur*ur+wc*wc)
    endif

    rib=grav*zldis*dthv/(thv*um*um)

    if (rib >= 0._r8) then      ! neutral or stable
       zeta = rib*log(zldis/z0m)/(1._r8-5._r8*min(rib,0.19_r8))
       zeta = min(2._r8,max(zeta,0.01_r8 ))
    else                     ! unstable
       zeta=rib*log(zldis/z0m)
       zeta = max(-100._r8,min(zeta,-0.01_r8 ))
    endif

    obu=zldis/zeta

  end subroutine MoninObukIni

end module FrictionVelocityMod
