module FrictionVelocityMod

!#py #include "shr_assert.h"

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculation of the friction velocity, relation for potential
  ! temperature and humidity profiles of surface boundary layer.
  !
  ! !USES:
  use shr_kind_mod         , only : r8 => shr_kind_r8
  !#py !#py use shr_log_mod          , only : errMsg => shr_log_errMsg
  use FrictionVelocityType , only : frictionvel_type
  use LandunitType         , only : lun_pp
  use VegetationType       , only : veg_pp
  !
  ! !PUBLIC TYPES:
  implicit none

  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: FrictionVelocity       ! Calculate friction velocity
  public :: MoninObukIni           ! Initialization of the Monin-Obukhov length
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: StabilityFunc1        ! Stability function for rib < 0.
  private :: StabilityFunc2        ! Stability function for rib < 0.
  !------------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------------
  subroutine FrictionVelocity(lbn, ubn, fn, filtern, &
       displa, z0m, z0h, z0q, &
       obu, iter, ur, um, ustar, &
       temp1, temp2, temp12m, temp22m, fm, frictionvel_vars, landunit_index)
    !$acc routine seq
    ! !DESCRIPTION:
    ! Calculation of the friction velocity, relation for potential
    ! temperature and humidity profiles of surface boundary layer.
    ! The scheme is based on the work of Zeng et al. (1998):
    ! Intercomparison of bulk aerodynamic algorithms for the computation
    ! of sea surface fluxes using TOGA CORE and TAO data. J. Climate,
    ! Vol. 11, 2628-2644.
    !
    ! !USES:
    use clm_varcon, only : vkc

    implicit none
    !
    ! !ARGUMENTS:
    integer  , intent(in)    :: lbn, ubn                 ! pft/landunit array bounds
    integer  , intent(in)    :: fn                       ! number of filtered pft/landunit elements
    integer  , intent(in)    :: filtern(fn)              ! pft/landunit filter
    real(r8) , intent(in)    :: displa  ( lbn: )         ! displacement height (m) [lbn:ubn]
    real(r8) , intent(in)    :: z0m     ( lbn: )         ! roughness length over vegetation, momentum [m] [lbn:ubn]
    real(r8) , intent(in)    :: z0h     ( lbn: )         ! roughness length over vegetation, sensible heat [m] [lbn:ubn]
    real(r8) , intent(in)    :: z0q     ( lbn: )         ! roughness length over vegetation, latent heat [m] [lbn:ubn]
    real(r8) , intent(in)    :: obu     ( lbn: )         ! monin-obukhov length (m) [lbn:ubn]
    integer  , intent(in)    :: iter                     ! iteration number
    real(r8) , intent(in)    :: ur      ( lbn: )         ! wind speed at reference height [m/s] [lbn:ubn]
    real(r8) , intent(in)    :: um      ( lbn: )         ! wind speed including the stablity effect [m/s] [lbn:ubn]
    real(r8) , intent(out)   :: ustar   ( lbn: )         ! friction velocity [m/s] [lbn:ubn]
    real(r8) , intent(out)   :: temp1   ( lbn: )         ! relation for potential temperature profile [lbn:ubn]
    real(r8) , intent(out)   :: temp12m ( lbn: )         ! relation for potential temperature profile applied at 2-m [lbn:ubn]
    real(r8) , intent(out)   :: temp2   ( lbn: )         ! relation for specific humidity profile [lbn:ubn]
    real(r8) , intent(out)   :: temp22m ( lbn: )         ! relation for specific humidity profile applied at 2-m [lbn:ubn]
    real(r8) , intent(inout) :: fm      ( lbn: )         ! diagnose 10m wind (DUST only) [lbn:ubn]
    type(frictionvel_type) , intent(inout) :: frictionvel_vars
    logical  , intent(in), optional :: landunit_index   ! optional argument that defines landunit or pft level
    !
    ! !LOCAL VARIABLES:
    real(r8), parameter :: zetam = 1.574_r8 ! transition point of flux-gradient relation (wind profile)
    real(r8), parameter :: zetat = 0.465_r8 ! transition point of flux-gradient relation (temp. profile)
    integer  :: f                           ! pft/landunit filter index
    integer  :: n                           ! pft/landunit index
    integer  :: g                           ! gridcell index
    integer  :: pp                          ! pfti,pftf index
    integer :: pfti, pftf
    real(r8) :: zldis(lbn:ubn)             ! reference height "minus" zero displacement heght [m]
    real(r8) :: zeta(lbn:ubn)                ! dimensionless height used in Monin-Obukhov theory
    real(r8) :: tmp1,tmp2,tmp3,tmp4         ! Used to diagnose the 10 meter wind
    real(r8) :: fmnew                       ! Used to diagnose the 10 meter wind
    real(r8) :: fm10                        ! Used to diagnose the 10 meter wind
    real(r8) :: zeta10                      ! Used to diagnose the 10 meter wind
    real(r8) :: vds_tmp                     ! Temporary for dry deposition velocity
    logical  :: lnd_index
    !------------------------------------------------------------------------------
    !Enforce expected array sizes



      if (present(landunit_index)) then
              lnd_index = .true.
      else
              lnd_index = .false.
      end if

      ! Adjustment factors for unstable (moz < 0) or stable (moz > 0) conditions.
      if(lnd_index) then

        do f = 1, fn
          n = filtern(f)
          g = lun_pp%gridcell(n)
          pfti = lun_pp%pfti(n)
          pftf = lun_pp%pftf(n)
          !Wind Profile
          zldis(n) = frictionvel_vars%forc_hgt_u_patch(pfti)-displa(n)
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

          do pp = pfti, pftf
             frictionvel_vars%vds_patch(pp) = vds_tmp
          end do

          do pp = pfti, pftf
             if (zldis(n)-z0m(n) <= 10._r8) then
                frictionvel_vars%u10_clm_patch(pp) = um(n)
             else
                if (zeta(n) < -zetam) then
                   frictionvel_vars%u10_clm_patch(pp) = um(n) - ( ustar(n)/vkc*(log(-zetam*obu(n)/(10._r8+z0m(n)))      &
                        - StabilityFunc1(-zetam)                              &
                        + StabilityFunc1((10._r8+z0m(n))/obu(n))              &
                        + 1.14_r8*((-zeta(n))**0.333_r8-(zetam)**0.333_r8)) )
                else if (zeta(n) < 0._r8) then
                   frictionvel_vars%u10_clm_patch(pp) = um(n) - ( ustar(n)/vkc*(log(zldis(n)/(10._r8+z0m(n)))           &
                        - StabilityFunc1(zeta(n))                             &
                        + StabilityFunc1((10._r8+z0m(n))/obu(n))) )
                else if (zeta(n) <=  1._r8) then
                   frictionvel_vars%u10_clm_patch(pp) = um(n) - ( ustar(n)/vkc*(log(zldis(n)/(10._r8+z0m(n)))           &
                        + 5._r8*zeta(n) - 5._r8*(10._r8+z0m(n))/obu(n)) )
                else
                   frictionvel_vars%u10_clm_patch(pp) = um(n) - ( ustar(n)/vkc*(log(obu(n)/(10._r8+z0m(n)))             &
                        + 5._r8 - 5._r8*(10._r8+z0m(n))/obu(n)                &
                        + (5._r8*log(zeta(n))+zeta(n)-1._r8)) )

                end if
             end if
             frictionvel_vars%va_patch(pp) = um(n)
          end do
          !===================!
          !Temperature Profile!
          !===================!
          zldis(n) = frictionvel_vars%forc_hgt_t_patch(pfti)-displa(n)
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

          !================!
          !Humidity profile!
          !================!
          if (frictionvel_vars%forc_hgt_q_patch(pfti) == frictionvel_vars%forc_hgt_t_patch(pfti) &
            .and. z0q(n) == z0h(n)) then
             temp2(n) = temp1(n)
          else
             zldis(n) = frictionvel_vars%forc_hgt_q_patch(pfti)-displa(n)
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


          zldis(n) = frictionvel_vars%forc_hgt_u_patch(pfti)-displa(n)

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
          tmp4 = log( max( 1.0_8, frictionvel_vars%forc_hgt_u_patch(pfti) / 10._r8) )

          do pp = pfti, pftf
              frictionvel_vars%u10_patch(pp) = ur(n) - ustar(n)/vkc * (tmp4 - fm(n) + fm10)
              frictionvel_vars%fv_patch(pp)  = ustar(n)
          end do

        end do

      end if
!!!!!========================================================================!!!!
!!!!========================================================================!!!!

      if(.not. lnd_index) then

        do f = 1, fn
          n = filtern(f)

          g = veg_pp%gridcell(n)
          !Wind Profile
          zldis(n) = frictionvel_vars%forc_hgt_u_patch(n)-displa(n)

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

          frictionvel_vars%vds_patch(n) = vds_tmp

          ! Calculate a 10-m wind (10m + z0m + d)
          ! For now, this will not be the same as the 10-m wind calculated for the dust
          ! model because the CLM stability functions are used here, not the LSM stability
          ! functions used in the dust model. We will eventually change the dust model to be
          ! consistent with the following formulation.
          ! Note that the 10-m wind calculated this way could actually be larger than the
          ! atmospheric forcing wind because 1) this includes the convective velocity, 2)
          ! this includes the 1 m/s minimum wind threshold

          ! If forcing height is less than or equal to 10m, then set 10-m wind to um

          if (zldis(n)-z0m(n) <= 10._r8) then
             frictionvel_vars%u10_clm_patch(n) = um(n)
          else
             if (zeta(n) < -zetam) then
                frictionvel_vars%u10_clm_patch(n) = um(n) - ( ustar(n)/vkc*(log(-zetam*obu(n)/(10._r8+z0m(n)))         &
                     - StabilityFunc1(-zetam)                                 &
                     + StabilityFunc1((10._r8+z0m(n))/obu(n))                 &
                     + 1.14_r8*((-zeta(n))**0.333_r8-(zetam)**0.333_r8)) )
             else if (zeta(n) < 0._r8) then
                frictionvel_vars%u10_clm_patch(n) = um(n) - ( ustar(n)/vkc*(log(zldis(n)/(10._r8+z0m(n)))              &
                     - StabilityFunc1(zeta(n))                                &
                     + StabilityFunc1((10._r8+z0m(n))/obu(n))) )
             else if (zeta(n) <=  1._r8) then
                frictionvel_vars%u10_clm_patch(n) = um(n) - ( ustar(n)/vkc*(log(zldis(n)/(10._r8+z0m(n)))              &
                     + 5._r8*zeta(n) - 5._r8*(10._r8+z0m(n))/obu(n)) )
             else
                frictionvel_vars%u10_clm_patch(n) = um(n) - ( ustar(n)/vkc*(log(obu(n)/(10._r8+z0m(n)))                &
                     + 5._r8 - 5._r8*(10._r8+z0m(n))/obu(n)                   &
                     + (5._r8*log(zeta(n))+zeta(n)-1._r8)) )
             end if
          end if
          frictionvel_vars%va_patch(n) = um(n)
          !===================!
          !Temperature Profile!
          !===================!
          zldis(n) = frictionvel_vars%forc_hgt_t_patch(n)-displa(n)
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
          !=================!
          !Humidity Profile !
          !=================!
          if (frictionvel_vars%forc_hgt_q_patch(n) == frictionvel_vars%forc_hgt_t_patch(n) .and. z0q(n) == z0h(n)) then
             temp2(n) = temp1(n)
          else
             zldis(n) = frictionvel_vars%forc_hgt_q_patch(n)-displa(n)
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


          zldis(n) = frictionvel_vars%forc_hgt_u_patch(n)-displa(n)

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

          tmp4 = log( max( 1.0_8, frictionvel_vars%forc_hgt_u_patch(n) / 10._r8) )

          frictionvel_vars%u10_patch(n) = ur(n) - ustar(n)/vkc * (tmp4 - fm(n) + fm10)
          frictionvel_vars%fv_patch(n)  = ustar(n)

        end do !! do loop of fn

      end if !! not landunit_index

  end subroutine FrictionVelocity

  !------------------------------------------------------------------------------
  real(r8) function StabilityFunc1(zeta)
    !$acc routine seq
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
    ! !LOCAL VARIABLES:
    real(r8) :: chik, chik2
    !------------------------------------------------------------------------------

    chik2 = sqrt(1._r8-16._r8*zeta)
    chik = sqrt(chik2)
    StabilityFunc1 = 2._r8*log((1._r8+chik)*0.5_r8) &
         + log((1._r8+chik2)*0.5_r8)-2._r8*atan(chik)+SHR_CONST_PI*0.5_r8

  end function StabilityFunc1

  !------------------------------------------------------------------------------
  real(r8) function StabilityFunc2(zeta)
    !$acc routine seq
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
    ! !LOCAL VARIABLES:
    real(r8) :: chik2
    !------------------------------------------------------------------------------

    chik2 = sqrt(1._r8-16._r8*zeta)
    StabilityFunc2 = 2._r8*log((1._r8+chik2)*0.5_r8)

  end function StabilityFunc2

  !-----------------------------------------------------------------------
  subroutine MoninObukIni (ur, thv, dthv, zldis, z0m, um, obu)
    !$acc routine seq
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
    ! !LOCAL VARIABLES:
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
