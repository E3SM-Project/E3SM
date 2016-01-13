!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_basal_traction.F90 - part of the Community Ice Sheet Model (CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2014
!   CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of CISM.
!
!   CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "glide_mask.inc"
#include "config.inc"

  module glissade_basal_traction

  !-----------------------------------------------------------------------------
  ! Compute or prescribe the basal traction coefficient 'beta' as required by
  ! the higher-order velocity solver.
  ! 
  ! Note that beta is assumed to be a positive constant.  In earlier versions of
  ! the code it was called 'betasquared'.
  !
  ! The units are Pa/(m/yr) if we assume a linear sliding law of the form
  !    taub_x = -beta * u, taub_y = -beta * v
  !
  ! However, the units are Pa if beta is treated as a till yield stress.
  !
  ! Current options are as follows:
  ! 
  ! [0] constant value of 10 Pa/(m/yr) (useful for debugging)
  ! [1] simple hard-coded pattern (useful for debugging)
  ! [2] treat beta value as a till yield stress (in Pa) using Picard iteration 
  ! [3] linear (inverse) function of basal water depth (bwat) 
  ! [4] very large value for beta to enforce no slip everywhere 
  ! [5] beta field passed in from .nc input file as part of standard i/o
  ! [6] no slip everywhere (using Dirichlet BC rather than large beta)
  ! [7] treat beta value as till yield stress (in Pa) using Newton-type iteration (in devel.)
  ! [8] set beta as prescribed for ISMIP-HOM test C (serial only)
  ! [9] power law that uses effective pressure
  ! [10] Coulomb friction law of Schoof (2005)
  ! [11] Coulomb friction law of Schoof (2005) with prescribed value of A = flwa at bed
  ! [12] Friction law of Tsai et al. (2015): minimum of power-law stress and Coulomb stress

  ! TODO - Renumber HO_BABC options so that, for example, the no-slip options have small numbers?
  !-----------------------------------------------------------------------------

  use glimmer_paramets, only : dp
  use glimmer_physcon,  only : scyr
  use glimmer_paramets, only : vel0, tau0
  use glimmer_log
  use glide_types
  use parallel,         only : staggered_parallel_halo  
  use glissade_grid_operators

  implicit none

!***********************************************************************

contains

!***********************************************************************

  subroutine calcbeta (whichbabc,                    &
                       dew,           dns,           &
                       ewn,           nsn,           &
                       thisvel,       othervel,      &
                       bwat,          beta_const,    &
                       mintauf,       basal_physics, &
                       flwa_basal,    thck,          &
                       mask,          beta_external, &
                       beta,                         &
                       f_ground,                     &
                       pmp_mask,                     &
                       beta_grounded_min)

  ! subroutine to calculate map of beta sliding parameter, based on 
  ! user input ("whichbabc" flag, from config file as "which_ho_babc").
   
  ! NOTE: Previously, the input arguments were assumed to be dimensionless
  ! and were rescaled in this routine.  Now the input arguments are
  ! assumed to have the units given below.
     
  use glimmer_paramets, only: len0
  use glimmer_physcon, only: gn
  use parallel, only: nhalo

  implicit none

  ! Input/output arguments

  integer, intent(in) :: whichbabc
  integer, intent(in) :: ewn, nsn

  real(dp), intent(in)                    :: dew, dns           ! m
  real(dp), intent(in), dimension(:,:)    :: thisvel, othervel  ! basal velocity components (m/yr)
  real(dp), intent(in), dimension(:,:)    :: bwat               ! basal water depth (m)
  real(dp), intent(in), dimension(:,:)    :: mintauf            ! till yield stress (Pa)
  real(dp), intent(in)                    :: beta_const         ! spatially uniform beta (Pa yr/m)
  type(glide_basal_physics), intent(in)   :: basal_physics      ! basal physics object
  real(dp), intent(in), dimension(:,:)    :: flwa_basal         ! flwa for the basal ice layer (Pa^{-3} yr^{-1})
  real(dp), intent(in), dimension(:,:)    :: thck               ! ice thickness
  integer,  intent(in), dimension(:,:)    :: mask               ! staggered grid mask
  real(dp), intent(in), dimension(:,:)    :: beta_external      ! fixed beta read from external file (Pa yr/m)
  real(dp), intent(inout), dimension(:,:) :: beta               ! basal traction coefficient (Pa yr/m)
                                                                ! Note: This is beta_internal in glissade
  real(dp), intent(in), dimension(:,:), optional :: f_ground    ! grounded ice fraction, 0 <= f_ground <= 1
  integer,  intent(in), dimension(:,:), optional :: pmp_mask    ! = 1 where bed is at pressure melting point, elsewhere = 0
  real(dp), intent(in), optional          :: beta_grounded_min  ! minimum beta value for grounded ice (Pa m^{-1} yr)

  ! Local variables

  real(dp) :: smallnum = 1.0d-2  ! m/yr

  ! SFP added for making beta a function of basal water flux 
  real(dp), dimension(:,:), allocatable :: unstagbeta
  real(dp) :: C, m

  real(dp) :: Ldomain   ! size of full domain
  real(dp) :: omega     ! frequency of beta field
  real(dp) :: dx, dy
  integer :: ilo, ihi, jlo, jhi  ! limits of beta field for ISHOM C case
  integer :: ew, ns

  ! variables for power law
  real(dp) :: powerlaw_p, powerlaw_q

  ! variables for Coulomb friction law
  real(dp) :: Coulomb_C   ! friction coefficient (unitless)
  real(dp) :: lambda_max  ! wavelength of bedrock bumps at subgrid scale (m)
  real(dp) :: m_max       ! maximum bed obstacle slope (unitless)
  real(dp), dimension(size(beta,1), size(beta,2)) :: big_lambda       ! bedrock characteristics
  real(dp), dimension(size(beta,1), size(beta,2)) :: speed            ! ice speed, sqrt(uvel^2 + vvel^2), m/yr 
  integer,  dimension(size(thck,1), size(thck,2)) :: imask            ! ice grid mask  1=ice, 0=no ice
  real(dp), dimension(size(beta,1), size(beta,2)) :: flwa_basal_stag  ! flwa for the basal ice layer on the staggered grid
                                                                      ! Note: Units are Pa^{-n} yr^{-1}
  ! variables for Tsai et al. parameterization
  real(dp) :: taub_powerlaw  ! basal shear stress given by a power law as in Tsai et al. (2015)
  real(dp) :: taub_Coulomb   ! basal shear stress given by Coulomb friction as in Tsai et al. (2015)

  character(len=300) :: message

  !WHL - debug - for diagnostic output
  integer, parameter :: itest = 450, jtest = 3
!  integer, parameter :: itest = 560, jtest = 3
!  integer, parameter :: itest = 1120, jtest = 3

  select case(whichbabc)

    case(HO_BABC_CONSTANT)  ! spatially uniform value; useful for debugging and test cases

!      beta(:,:) = 10.d0       ! This is the default value (Pa yr/m)
      beta(:,:) = beta_const   ! Pa yr/m

    case(HO_BABC_SIMPLE)    ! simple pattern; also useful for debugging and test cases
                            ! (here, a strip of weak bed surrounded by stronger bed to simulate an ice stream)

      beta(:,:) = 1.d4        ! Pa yr/m

      !TODO - Change this loop to work in parallel (set beta on the global grid and scatter to local)
      do ns = 5, nsn-5
         do ew = 1, ewn-1
            beta(ew,ns) = 100.d0      ! Pa yr/m
         end do
      end do

    case(HO_BABC_BETA_TPMP)     ! large value for frozen bed, lower value for bed at pressure melting point

       ! Set beta = beta_const wherever the bed is at the pressure melting point temperature.
       ! Elsewhere, set beta to a large value.

       if (present(pmp_mask)) then
          
          where(pmp_mask == 1)
             beta(:,:) = beta_const    ! constant that can be specified in config file; 10 Pa yr/m by default
          elsewhere 
             beta(:,:) = 1.d10         ! Pa yr/m
          endwhere
          
       else

          write(message,*) 'Must supply pressure-melting-point mask with HO_BABC_BETA_TPMP option'
          call write_log(trim(message), GM_FATAL)

       endif

    case(HO_BABC_YIELD_PICARD)  ! take input value for till yield stress and force beta to be implemented such
                                ! that plastic-till sliding behavior is enforced (see additional notes in documentation).

      !!! NOTE: Eventually, this option will provide the till yield stress as calculate from the basal processes
      !!! submodel. Currently, to enable sliding over plastic till, simple specify the value of "beta" as 
      !!! if it were the till yield stress (in units of Pascals).
      
      beta(:,:) = mintauf(:,:) &                                                         ! plastic yield stress (Pa)
                         / dsqrt( thisvel(:,:)**2 + othervel(:,:)**2 + (smallnum)**2 )   ! velocity components (m/yr)

      !!! since beta is updated here, communicate that info to halos
      call staggered_parallel_halo(beta)

    case(HO_BABC_BETA_BWAT)  ! set value of beta as proportional to value of bwat                                         

      !NOTE: This parameterization has not been scientifically tested.
      !TODO - Test option HO_BABC_BETA_BWAT
      !       Where do these constants come from?
      C = 10.d0   ! Does this assume that bwat is in units of m or dimensionless?
      m = 1.d0

      allocate(unstagbeta(ewn,nsn))

      unstagbeta(:,:) = 200.d0   ! Pa yr/m
                                 ! This setting ensures that the parameterization does nothing.  Remove it?

      where ( bwat > 0.d0 .and. unstagbeta > 200.d0 )
          unstagbeta = C / ( bwat**m )
      endwhere

      ! average beta from unstag grid onto stag grid
      beta = 0.5d0 * ( unstagbeta(1:ewn-1,:) + unstagbeta(2:ewn,:) )
      beta = 0.5d0 * ( unstagbeta(:,1:nsn-1) + unstagbeta(:,2:nsn) )
   
      deallocate(unstagbeta) 

    !Note: This is redundant in that it could be implemented by using HO_BETA_CONST with beta_const = 1.d10
    !      But keeping it for historical reasons since many config files use it

    case(HO_BABC_LARGE_BETA)      ! frozen (u=v=0) ice-bed interface

      beta(:,:) = 1.d10           ! Pa yr/m

    case(HO_BABC_ISHOMC)          ! prescribe according to ISMIP-HOM test C

       !Note: Ideally, beta would be read in from an external netCDF file.
       !      However, this is not possible given that the global velocity grid is smaller
       !       than the ice grid and hence not able to fit the full beta field.
       !      The following code sets beta on the full grid as prescribed by Pattyn et al. (2008).
       !NOTE: This works only in serial!

       Ldomain = (ewn-2*nhalo) * dew   ! size of full domain (must be square)
       omega = 2.d0*pi / Ldomain

       ilo = nhalo
       ihi = ewn-nhalo
       jlo = nhalo
       jhi = nsn-nhalo
       
       ! Prescribe beta as in Pattyn et al., The Cryosphere, 2008
       beta(:,:) = 0.d0
       do ns = jlo, jhi
          do ew = ilo, ihi
             dx = dew * (ew-ilo)
             dy = dns * (ns-jlo)
             beta(ew,ns) = 1000.d0 + 1000.d0 * sin(omega*dx) * sin(omega*dy)
          enddo
       enddo

    case(HO_BABC_EXTERNAL_BETA)   ! use value passed in externally from CISM

       ! set beta to the prescribed external value
       ! Note: This assumes that beta_external has units of Pa yr/m on input.
       beta(:,:) = beta_external(:,:)

       ! beta is initialized to a negative value; we can use that fact to check whether
       ! it has been read correctly from the file
       if (maxval(beta) < 0.d0) then
          call write_log('ERROR: Trying to use HO_BABC_EXTERNAL_BETA, but all beta values are < 0,')
          call write_log('which implies that beta could not be read from the input file.')
          call write_log('Make sure that beta is in the cism input file,')
          call write_log('or change which_ho_babc to a different option.')
          call write_log('Invalid value for beta. See log file for details.', GM_FATAL)
       end if

    case(HO_BABC_POWERLAW)   ! A power law that uses effective pressure
       ! See Cuffey & Paterson, Physics of Glaciers, 4th Ed. (2010), p. 240, eq. 7.17
       ! This is based on Weertman's classic sliding relation (1957) augmented by the bed-separation index described by Bindschadler (1983)
       !   ub = k taub^p N^-q
       ! rearranging for taub gives:
       !   taub = k^(-1/p) ub^(1/p) N^(q/p)

       ! p and q should be _positive_ exponents
       ! TODO: powerlaw_p and powerlaw_q could be turned into config parameters instead of hard-coded
       ! If p/=1, this is nonlinear in velocity
       ! Cuffey & Paterson recommend p=3 and q=1, and k dependent on thermal & mechanical properties of ice and inversely on bed roughness.   
       powerlaw_p = 3.0d0
       powerlaw_q = 1.0d0

       beta(:,:) = basal_physics%friction_powerlaw_k**(-1.0d0/powerlaw_p) &
            * basal_physics%effecpress_stag(:,:)**(powerlaw_q/powerlaw_p) &
            * dsqrt( thisvel(:,:)**2 + othervel(:,:)**2 )**(1.0d0/powerlaw_p-1.0d0)

    case(HO_BABC_COULOMB_FRICTION, HO_BABC_COULOMB_CONST_BASAL_FLWA)

      ! Basal stress representation using Coulomb friction law
      ! Coulomb sliding law: Schoof 2005 PRS, eqn. 6.2  (see also Pimentel, Flowers & Schoof 2010 JGR)

      ! Set up parameters needed for the friction law
      m_max = basal_physics%Coulomb_bump_max_slope       ! maximum bed obstacle slope(unitless)
      lambda_max = basal_physics%Coulomb_bump_wavelength ! wavelength of bedrock bumps (m)
      Coulomb_C = basal_physics%Coulomb_C                ! basal shear stress factor (Pa (m^-1 y)^1/3)

      ! Compute biglambda = wavelength of bedrock bumps [m] * flwa [Pa^-n yr^-1] / max bed obstacle slope [dimensionless]

      if (whichbabc == HO_BABC_COULOMB_FRICTION) then

         ! Need flwa of the basal layer on the staggered grid
         !TODO - Pass in ice_mask instead of computing imask here?
         !       (Small difference: ice_mask = 1 where thck > thklim rather than thck > 0)
         where (thck > 0.0)
            imask = 1
         elsewhere
            imask = 0
         end where
         call glissade_stagger(ewn,         nsn,               &
                               flwa_basal,  flwa_basal_stag,   &
                               imask,       stagger_margin_in = 1)
         ! TODO Not sure if a halo update is needed on flwa_basal_stag!  I don't think so if nhalo>=2.

         big_lambda(:,:) = (lambda_max / m_max) * flwa_basal_stag(:,:)

      else   ! which babc = HO_BABC_COULOMB_CONST_BASAL_FLWA; use a constant value of basal flwa
             ! NOTE: Units of flwa_basal are Pa{-n} yr{-1}

         big_lambda(:,:) = (lambda_max / m_max) * basal_physics%flwa_basal

      endif

      ! Note: For MISMIP3D, Coulomb_C is multiplied by a spatial factor (C_space_factor) which is
      !       read in during initialization. This factor is typically between 0 and 1.
      !       If this factor is not present in the input file, it is set to 1 everywhere.

      ! Compute beta
      ! gn = Glen's n from physcon module
      speed(:,:) = dsqrt(thisvel(:,:)**2 + othervel(:,:)**2 + smallnum**2)
      beta(:,:) = Coulomb_C * basal_physics%C_space_factor_stag(:,:) * &
           basal_physics%effecpress_stag(:,:) * speed(:,:)**(1.0d0/gn - 1.0d0) * &
           (speed(:,:) + basal_physics%effecpress_stag(:,:)**gn * big_lambda)**(-1.0d0/gn)

      ! Limit for numerical stability
      where (beta > 1.0d8)
         beta = 1.0d8
      end where

      !WHL - debug - Write values along a flowline
!      write(6,*) ' '
!      write(6,*) 'Apply Coulomb friction: i, j, speed, beta, taub:'
!      ns = jtest
!      do ew = itest, itest+15
!         write(6,*) ew, ns, speed(ew,ns), beta(ew,ns), beta(ew,ns)*speed(ew,ns)
!      enddo

    case(HO_BABC_COULOMB_POWERLAW_TSAI)

      ! Basal stress representation based on Tsai et al. (2015)
      ! The basal stress is the minimum of two values:
      ! (1) power law:          tau_b = powerlaw_C * |u_b|^(1/powerlaw_m)
      ! (2) Coulomb friction:   tau_b = Coulomb_C * N
      !                             N = effective pressure = rhoi*g*(H - H_f)
      !                           H_f = flotation thickness = (rhow/rhoi)*(eus-topg)
      ! This value of N is obtained by setting basal_water = BWATER_OCEAN_PENETRATION = 4 with p_ocean_penetration = 1.0 in the config file.
      ! The other parameters (powerlaw_C, powerlaw_m and Coulomb_C) can also be set in the config file. 

       !WHL - debug - write out basal stresses
!       write(6,*) ' '
!!       write(6,*) 'powerlaw_C, powerlaw_m, Coulomb_C =', basal_physics%powerlaw_C, basal_physics%powerlaw_m, basal_physics%Coulomb_C
!       write(6,*) 'Apply Tsai parameterization: i, j, speed, beta, taub, taub_powerlaw, taub_Coulomb, effecpress:'

       do ns = 1, nsn-1
          do ew = 1, ewn-1
             
             speed(ew,ns) = dsqrt(thisvel(ew,ns)**2 + othervel(ew,ns)**2 + smallnum**2)

             taub_powerlaw = basal_physics%powerlaw_C * speed(ew,ns)**(1.d0/basal_physics%powerlaw_m)
             taub_Coulomb  = basal_physics%Coulomb_C * basal_physics%effecpress_stag(ew,ns)

             if (taub_Coulomb <= taub_powerlaw) then   ! apply Coulomb stress, which is smaller
                beta(ew,ns) = taub_Coulomb / speed(ew,ns)
             else  ! apply power-law stress
                beta(ew,ns) = taub_powerlaw / speed(ew,ns)
             endif

!             !WHL - debug - Write values along a flowline
!             if (ns == jtest .and. ew >= itest .and. ew <= itest+15) then
!                write(6,*) ew, ns, speed(ew,ns), beta(ew,ns), speed(ew,ns)*beta(ew,ns), taub_powerlaw, taub_Coulomb, basal_physics%effecpress_stag(ew,ns)
!             endif

          enddo   ! ew
       enddo   ! ns

    case default
       ! do nothing

   end select

   ! If f_ground is passed in (as for Glissade), then multiply beta by f_ground (0 <= f_ground <= 1) 
   !  to reduce the basal traction in regions that are partially or totally floating.
   ! Note: With a GLP, f_ground will have values between 0 and 1 at vertices adjacent to the GL.
   !       Without a GLP, f_ground = 0 or 1 everywhere based on a flotation criterion.
   !       By convention, f_ground = 0 where no ice is present.
   !
   ! If f_ground in not passed in (as for Glam), then check for areas where the ice is floating
   !  and make sure beta in these regions is 0. 

   if (present(f_ground)) then   ! Multiply beta by grounded ice fraction

      beta(:,:) = beta(:,:) * f_ground(:,:)

   else    ! set beta = 0 where Glide mask says the ice is floating

      do ns = 1, nsn-1
         do ew = 1, ewn-1 
            if (GLIDE_IS_FLOAT(mask(ew,ns))) then
               beta(ew,ns) = 0.d0
            endif
         end do
      end do

   endif   ! present(f_ground)

   ! For beta close to 0 beneath grounded ice, it is possible to generate unrealistically fast flow.
   ! This could happen, for example, when reading beta from an external file.
   ! To prevent this, set beta to a minimum value beneath grounded ice.
   ! The default value of beta_grounded_min = 0.0, in which case this loop has no effect.
   ! However, beta_grounded_min can be set to a nonzero value in the config file.
  
   if (present(f_ground) .and. present(beta_grounded_min)) then

      do ns = 1, nsn-1
         do ew = 1, ewn-1
            if (f_ground(ew,ns) > 0.d0 .and. beta(ew,ns) < beta_grounded_min) then
               beta(ew,ns) = beta_grounded_min
!!               print*, 'Reset beta: ew, ns, f_ground, beta:', ew, ns, f_ground(ew,ns), beta(ew,ns)
            endif
         enddo
      enddo

   endif   ! present(f_ground) and present(beta_grounded_min)

   ! Bug check: Make sure beta >= 0
   ! This check will find negative values as well as NaNs
   do ns = 1, nsn-1
      do ew = 1, ewn-1 
         if (beta(ew,ns) >= 0.d0) then
            ! do nothing
         else
            write(message,*) 'Invalid beta value in calcbeta: ew, ns, beta:', ew, ns, beta(ew,ns)
            call write_log(trim(message), GM_FATAL)
         endif
      end do
   end do
   
 end subroutine calcbeta

!***********************************************************************

end module glissade_basal_traction

!***********************************************************************
