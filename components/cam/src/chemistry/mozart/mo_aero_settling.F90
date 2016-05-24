
!======================================================================
!
! ROUTINE
!   strat_aer_settling
!
!   Date...
!     8  November 1999
!
!   Programmed by...
!     Douglas E. Kinnison
!
!   Modified for WACCM2 on 3 September 2004 - Removed ICE
!   Modified for WACCM3 on 8 November  2004 - added NAT 7um mode.
!
! DESCRIPTION
!   This routine vertically redistributes condensed phase HNO3.
!
!   For each aerosol type the terminal velocity is calculated. This
!   quantity is dependent on: 1) the mass density of the aerosol;
!   2) the radius; 3) the dynamic viscosity; 4) shape; and the Cunningham
!   correction factor for spherical particles. See Fuchs, The Mechanics
!   of Aerosols Oxford, Pergmann Press, pp 27-31, 1964 and Kasten,
!   Falling Speed of Aerosol Particles, J. of Appl. Met., 7, 944-947,
!   1968 for details. For aerosol with a radius of 3 microns (e.g., NAT)
!   and 10 microns (e.g., ICE) the terminal velocity (cm sec-1) is  0.2
!   and 1.7 respectively.
!
!   The flux of condensed phase HNO3 is then derived using the
!   following equation:
!          Flux (molec cm-2 s-1) = V * C * exp (8*ln^2 sigma).
!
!          where: V is terminal velocity (cm sec-1)
!                 C is condensed phase conc (molecules cm-3)
!                 sigma is the width of the log normal distribution.
!
!   The approach of settling the entire aerosol size distribution is
!   based on work of Considine et al., JGR, 1999.
!
!   The routine is a straighforward approach, starting at the top zone
!   (highest altitude) and progressing towards the surface. The gross
!   settling of condensed phase HNO3 is simply the quanity:
!   Flux * dt/dz.
!
!
! ARGUMENTS
!
!      All of the following components are grid center quanitities.
!
!   INPUT:
!      ad            Air Density (molecules cm-3)
!      press         Pressure (Pascals)
!      timestep      Gross chemistry timestep (in seconds)
!      temp          Temperature (K)
!      hno3_cond     Condensed phase HNO3 (mole fraction)
!      radius_nat    Mean radius of NAT (cm)
!      zstar         log pressure altitude coordinate (km)
!
!   OUTPUT:
!      hno3_cond     vertically modified
!
!======================================================================

      module mo_aero_settling

      use shr_kind_mod, only : r8 => shr_kind_r8

      private
      public  :: strat_aer_settling
      public  :: strat_aer_settl_init

      contains

        subroutine strat_aer_settl_init

          use cam_history,  only : addfld
          use ppgrid,       only : pver

          implicit none

          call addfld( 'VEL_NAT1', (/ 'lev' /), 'I', 'cm/s', 'small nat settling velocity' )
          call addfld( 'VEL_NAT2', (/ 'lev' /), 'I', 'cm/s', 'large nat settling velocity' )

        end subroutine strat_aer_settl_init


      subroutine strat_aer_settling( ad, press, timestep, zstar, temp, &
                                     hno3_cond, radius_nat, ncol, lchnk, aero_ndx )

      use ppgrid,      only : pcols, pver
      use chem_mods,   only : adv_mass
      use physconst,   only : gravit
      use cam_history, only : outfld

!-----------------------------------------------------------------------
!     ... dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in)     ::  ncol                                    ! columns in chunk
      integer, intent(in)     ::  lchnk                                   ! chunk id
      integer, intent(in)     ::  aero_ndx                                ! aerosol index
      real(r8), intent(in)    ::  timestep                                ! model time step (s)
      real(r8), intent(in)    ::  ad(ncol,pver)                           ! Air density (molecules cm-3)
      real(r8), intent(in)    ::  radius_nat(ncol,pver)                   ! Mean radius of NAT (cm)
      real(r8), intent(in)    ::  zstar(ncol,pver)                        ! altitude (km)
      real(r8), intent(in)    ::  press(pcols,pver)                       ! Pressure (Pa)
      real(r8), intent(in)    ::  temp(pcols,pver)                        ! temperature (K)
      real(r8), intent(inout) ::  hno3_cond(ncol,pver)                    ! Condensed Phase HNO3 (VMR)

!-----------------------------------------------------------------------
!     ... local variables
!-----------------------------------------------------------------------
      real(r8), parameter :: avo_num       = 6.022e23_r8, &    ! molecules/mole
                             MW_air        = 28.8_r8, &        ! grams/mole air
                             nat_dens      = 1.6_r8, &         ! g/cm^3
                             shape_fac_nat = 1.0_r8, &         ! TBD
                             sigma_nat     = 1.6_r8, &         ! Width of distribution
                             av_const      = 2.117265e4_r8, &  ! (8*8.31448*1000 / PI)
                             km2cm         = 1.e5_r8, &        ! km to cm
                             m2cm          = 1.e2_r8, &        ! m to cm
                             c1            = 2._r8/9._r8

      integer  :: i, k, kp1
      real(r8) :: gravity                   ! gravity cm/s^2
      real(r8) :: Cc_nat                    ! Cunningham Correction Factor
      real(r8) :: dt_dz                     ! dt / dz, sec cm-1
      real(r8) :: flux_nat                  ! aerosol flux, molec cm-2 sec-1
      real(r8) :: mean_vel                  ! mean velocity, cm sec-1
      real(r8) :: mfp                       ! Mean Free Path
      real(r8) :: vel_nat                   ! Terminal velocity (cm sec-1)
      real(r8) :: rad_nat                   ! wrk radius nat (cm)
      real(r8) :: visc                      ! Dynamic viscosity of air
      real(r8) :: depos                     ! molecules/cm**3 deposited
      real(r8) :: atm_dens, atm_densa       ! total atm density and inverse
      real(r8) :: cond_hno3                 ! wrk variables
      real(r8) :: const_nat                 ! wrk variables
      real(r8) :: t                         ! working temperatue
      real(r8) :: velnat(ncol,pver)         ! holding variable for output
      logical  :: lon_mask(ncol)            ! longitude logic mask

      gravity       = gravit*m2cm       ! (cm/s^2)

      const_nat = exp( 8._r8*(log(sigma_nat))**2 )
      do k = 1,pver
	 velnat(:,k) = 0._r8
      end do


!----------------------------------------------------------------------
!     ... derive Aerosol Settling (explicit approach)
!----------------------------------------------------------------------
Level_loop : &
      do k = 2,pver-1
!----------------------------------------------------------------------
!     ... operate between 2.0hPa and 300hPa and only where nat exist
!----------------------------------------------------------------------
	 kp1 = k + 1
	 lon_mask(:) = press(:ncol,k) >= 2.e2_r8 .and. press(:ncol,k) <= 300.e2_r8 &
                       .and. (radius_nat(:,k) > 0._r8 )
	 if( any( lon_mask(:) ) ) then
Column_loop : &
            do i = 1,ncol
	       if( lon_mask(i) ) then
                  t = temp(i,k)
!----------------------------------------------------------------------
!     ... General Setup for NAT
!         Calculate the settling of the NAT Aerosol
!         NOTE: Index "k" is the box that is being calculated
!               Index "k-1" is flux from above into the k box
!               A positive NET {flux*dt/dz} adds to the "k" box
!----------------------------------------------------------------------
!     ... mean Molecular Velocity, cm sec-1
!----------------------------------------------------------------------
                  mean_vel  = sqrt( av_const*t/MW_air )*100._r8
!----------------------------------------------------------------------
!     ... dynamic Viscosity, g cm-1 sec-1
!----------------------------------------------------------------------
	          visc = (t*1.458e-6_r8)**1.5_r8 /(t + 110.4_r8)*10000._r8
!----------------------------------------------------------------------
!     ... mean Free Path, cm
!----------------------------------------------------------------------
		  atm_dens  = ad(i,k)
		  atm_densa = 1._r8/atm_dens
                  mfp = visc* avo_num / (.499_r8*atm_dens*MW_air*mean_vel)
!----------------------------------------------------------------------
!     ... dt / dz, sec cm-1; NOTE: zstar is in km
!----------------------------------------------------------------------
                  dt_dz = timestep / ((zstar(i,k-1) - zstar(i,k))*km2cm)
!----------------------------------------------------------------------
!     ... calculate NAT Aerosol Settling
!----------------------------------------------------------------------
		  rad_nat   = radius_nat(i,k)
                  cond_hno3 = hno3_cond(i,k)*atm_dens
!----------------------------------------------------------------------
!     ... Cunningham Correction Factor, Unitless
!----------------------------------------------------------------------
                  Cc_nat = 1._r8 + (mfp/rad_nat)*(1.246_r8 + .42_r8*exp( -.87_r8*rad_nat/mfp ))
!----------------------------------------------------------------------
!     ... terminal Velocity of Aerosol, cm sec-1
!----------------------------------------------------------------------
	          vel_nat = c1*rad_nat**2 * nat_dens*gravity*Cc_nat/visc*shape_fac_nat
		  velnat(i,k) = vel_nat
!----------------------------------------------------------------------
!     ... aerosol Flux, Cond-phase molecules cm-2 sec-1
!----------------------------------------------------------------------
         	  flux_nat = cond_hno3*vel_nat* const_nat
!----------------------------------------------------------------------
!     ... calculate NAT Aerosol Settling (i.e., HNO3 redistribution)
!----------------------------------------------------------------------
		  depos = min( cond_hno3,flux_nat*dt_dz )
!----------------------------------------------------------------------
!     ... modify the HNO3_cond in level "k"
!----------------------------------------------------------------------
	          hno3_cond(i,k) = (cond_hno3 - depos)*atm_densa
!----------------------------------------------------------------------
!     ... modify the HNO3_cond in level "k+1"
!----------------------------------------------------------------------
		  hno3_cond(i,kp1) = hno3_cond(i,kp1) + depos/ad(i,kp1)
               end if
            end do Column_loop
         end if
      end do Level_loop

!----------------------------------------------------------------------
!	... output nat velocity
!----------------------------------------------------------------------
      if( aero_ndx == 1 ) then
         call outfld( 'VEL_NAT1', velnat, ncol, lchnk )
      else if( aero_ndx == 2 ) then
         call outfld( 'VEL_NAT2', velnat, ncol, lchnk )
      end if

      end subroutine strat_aer_settling

      end module mo_aero_settling
