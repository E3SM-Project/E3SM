! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_stokes_drift
!
!> \brief  Routines related to calulationing Stokes drift profiles
!> \author Steven Brus 
!> \date   June 2021
!> \details
!>  This module contains the routines for calculating wave-induced 
!>  Stokes drift profiles from an active wave model.
!
!-----------------------------------------------------------------------

module ocn_stokes_drift

   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_log
   use ocn_constants
   use ocn_config
   use ocn_mesh
   use ocn_diagnostics_variables


   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_stokes_drift_reconstruct, &
             ocn_surface_stokes_drift, &
             ocn_stokes_drift_kpp_enhancement_factor, &
             ocn_stokes_drift_langmuir_number, &
             ocn_stokes_drift_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   real (kind=RKIND), dimension(:), allocatable :: stokesDriftWavenumbers


!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_stokes_drift_reconstruct
!
!> \brief   Reconstructs depth-dependent Stokes drift profile 
!> \author  Steven Brus
!> \date    June 2021 
!> \details
!>   This routine reconstructs the depth-dependent Stokes drift profile 
!>   using partitioned surface Stokes drift in 3, 4, or 6 wavenumber bins
!>   from an active wave model. 
!
!-----------------------------------------------------------------------

   subroutine ocn_stokes_drift_reconstruct(forcingPool)!{{{

      !-----------------------------------------------------------------
      ! Input variables
      !-----------------------------------------------------------------

      type(mpas_pool_type), intent(inout) :: &
         forcingPool

      !-----------------------------------------------------------------
      ! Output variables
      !-----------------------------------------------------------------


      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------

      integer :: iCell, k
      real(kind=RKIND), dimension(:,:), pointer :: stokesDriftZonal, stokesDriftMeridional
      real(kind=RKIND), dimension(:,:), pointer :: stokesDriftZonalWavenumber, stokesDriftMeridionalWavenumber

      !-----------------------------------------------------------------

      call mpas_pool_get_array(forcingPool, 'stokesDriftZonal', stokesDriftZonal)
      call mpas_pool_get_array(forcingPool, 'stokesDriftMeridional', stokesDriftMeridional)

      call mpas_pool_get_array(forcingPool, 'stokesDriftZonalWavenumber', stokesDriftZonalWavenumber)
      call mpas_pool_get_array(forcingPool, 'stokesDriftMeridionalWavenumber', stokesDriftMeridionalWavenumber)

      do iCell = 1,nCellsAll

         do k = minLevelCell(iCell), maxLevelCell(iCell)
 

            call ocn_stokes_drift_evaluate(zMid(k,iCell), &
                                           stokesDriftZonalWavenumber(:,iCell), stokesDriftMeridionalWavenumber(:,iCell), &
                                           stokesDriftZonal(k,iCell), stokesDriftMeridional(k,iCell))
         end do

      end do 

      !-----------------------------------------------------------------

   end subroutine ocn_stokes_drift_reconstruct!}}}

!***********************************************************************
!
!  routine ocn_surface_stokes_drift
!
!> \brief   Reconstructs surface Stokes drift
!> \author  Steven Brus
!> \date    December 2021 
!> \details
!>   This routine reconstructs the surface Stokes drift velocity
!>   using partitioned surface Stokes drift in 3, 4, or 6 wavenumber bins
!>   from an active wave model. 
!
!-----------------------------------------------------------------------

   subroutine ocn_surface_stokes_drift(forcingPool)!{{{

      !-----------------------------------------------------------------
      ! Input variables
      !-----------------------------------------------------------------

      type(mpas_pool_type), intent(inout) :: &
         forcingPool

      !-----------------------------------------------------------------
      ! Output variables
      !-----------------------------------------------------------------


      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------

      integer :: iCell
      real(kind=RKIND), dimension(:), pointer :: stokesDriftSurfaceZonal, stokesDriftSurfaceMeridional
      real(kind=RKIND), dimension(:,:), pointer :: stokesDriftZonalWavenumber
      real(kind=RKIND), dimension(:,:), pointer :: stokesDriftMeridionalWavenumber

      !-----------------------------------------------------------------

      call mpas_pool_get_array(forcingPool, 'stokesDriftSurfaceZonal', stokesDriftSurfaceZonal)
      call mpas_pool_get_array(forcingPool, 'stokesDriftSurfaceMeridional', stokesDriftSurfaceMeridional)

      call mpas_pool_get_array(forcingPool, 'stokesDriftZonalWavenumber', stokesDriftZonalWavenumber)
      call mpas_pool_get_array(forcingPool, 'stokesDriftMeridionalWavenumber', stokesDriftMeridionalWavenumber)

      do iCell = 1,nCellsAll

         call ocn_stokes_drift_evaluate(0.0_RKIND, &
                                        stokesDriftZonalWavenumber(:,iCell), stokesDriftMeridionalWavenumber(:,iCell), &
                                        stokesDriftSurfaceZonal(iCell), stokesDriftSurfaceMeridional(iCell))

      end do 

      !-----------------------------------------------------------------

   end subroutine ocn_surface_stokes_drift!}}}

!***********************************************************************
!
!  routine ocn_stokes_drift_langmuir_number
!
!> \brief Computes the langmuir number based on Van Roekel et al. 2012 
!> \author  Steven Brus
!> \date    June 2021 
!> \details
!>   This routine computes projected Langmuir number for misaligned wind and
!>   waves using the Stokes drift from an active wave model. 
!
!-----------------------------------------------------------------------

   subroutine ocn_stokes_drift_langmuir_number(windStressZonal, &
                                               windStressMeridional, &
                                               surfaceFrictionVelocity, &
                                               significantWaveHeight, &
                                               boundaryLayerDepth, &
                                               zonalWavenumberCoeff, &
                                               meridionalWavenumberCoeff, &
                                               alpha, &
                                               langmuirNumberMLAvgProjected)!{{{

      !-----------------------------------------------------------------
      ! Input variables
      !-----------------------------------------------------------------

      real(kind=RKIND), intent(in) :: windStressZonal
      real(kind=RKIND), intent(in) :: windStressMeridional
      real(kind=RKIND), intent(in) :: surfaceFrictionVelocity
      real(kind=RKIND), intent(in) :: significantWaveHeight 
      real(kind=RKIND), intent(in) :: boundaryLayerDepth
      real(kind=RKIND), dimension(:), intent(in) :: zonalWavenumberCoeff
      real(kind=RKIND), dimension(:), intent(in) :: meridionalWavenumberCoeff

      !-----------------------------------------------------------------
      ! Output variables
      !-----------------------------------------------------------------

      real(kind=RKIND), intent(out) :: alpha
      real(kind=RKIND), intent(out) :: langmuirNumberMLAvgProjected

      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------

      real(kind=RKIND) :: stokesDriftZonalSurface, stokesDriftMeridionalSurface
      real(kind=RKIND) :: stokesDriftMagnitudeSurface
      real(kind=RKIND) :: stokesDriftZonalMLAvg, stokesDriftMeridionalMLAvg
      real(kind=RKIND) :: stokesDriftMagnitudeMLAvg
      real(kind=RKIND) :: langmuirNumberMLAvg
      real(kind=RKIND) :: langmuirNumberTurb
      real(kind=RKIND) :: theta_ww
      real(kind=RKIND) :: denominator
      real(kind=RKIND) :: z0
      real(kind=RKIND), parameter :: small = 1e-10_RKIND
      real(kind=RKIND), parameter :: kappa = 0.4_RKIND

      !-----------------------------------------------------------------

      ! Compute Stokes drift from wavenumber partitioned surface values

      call ocn_stokes_drift_evaluate(0.0_RKIND, zonalWavenumberCoeff, meridionalWavenumberCoeff, &
                                     stokesDriftZonalSurface, stokesDriftMeridionalSurface)

      call ocn_stokes_drift_average(0.2_RKIND*boundaryLayerDepth, zonalWavenumberCoeff, meridionalWavenumberCoeff, &
                                    stokesDriftZonalMLAvg, stokesDriftMeridionalMLAvg)

      ! Calulate Stokes drift mangnitudes

      stokesDriftMagnitudeSurface = max(sqrt(stokesDriftZonalSurface**2 + stokesDriftMeridionalSurface**2),small)

      stokesDriftMagnitudeMLAvg = max(sqrt(stokesDriftZonalMLAvg**2 + stokesDriftMeridionalMLAvg**2),small)


      ! Calculate Langmuir numbers

      langmuirNumberTurb = sqrt(surfaceFrictionVelocity/stokesDriftMagnitudeSurface)

      langmuirNumberMLAvg = sqrt(surfaceFrictionVelocity/stokesDriftMagnitudeMLAvg)


      ! Calculate angle between wind and Langmuir cells

      theta_ww = atan2(stokesDriftMeridionalMLAvg, stokesDriftZonalMLAvg) - atan2(windStressMeridional, windStressZonal)

      z0 = 4.0_RKIND*max(significantWaveHeight,small)

      denominator = langmuirNumberTurb**2/kappa &
                  * log(max(abs(boundaryLayerDepth/z0),1.0_RKIND)) &
                  + cos(theta_ww)
      alpha = atan(sin(theta_ww)/denominator)


      ! Calculate surface layer averaged, projected Langmuir number

      langmuirNumberMLAvgProjected = langmuirNumberMLAvg*sqrt(abs(cos(alpha))/(abs(cos(theta_ww-alpha))+small))

      !-----------------------------------------------------------------

   end subroutine ocn_stokes_drift_langmuir_number!}}}

!***********************************************************************
!
!  routine ocn_stokes_drift_kpp_enhancement_factor
!
!> \brief Computes the KPP enhancement factor based on Van Roekel et al. 2012 
!> \author  Steven Brus
!> \date    June 2021 
!> \details
!>   This routine computes the KPP enhancment factor based on a projected 
!>   Langmuir number using the Stokes drift from an active wave model. 
!
!-----------------------------------------------------------------------

   subroutine ocn_stokes_drift_kpp_enhancement_factor(alpha, langmuirNumber, enhancementFactor)!{{{

      !-----------------------------------------------------------------
      ! Input variables
      !-----------------------------------------------------------------

      real(kind=RKIND), intent(in) :: alpha
      real(kind=RKIND), intent(in) :: langmuirNumber

      !-----------------------------------------------------------------
      ! Output variables
      !-----------------------------------------------------------------

      real(kind=RKIND), intent(out) :: enhancementFactor

      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------


      !-----------------------------------------------------------------

      enhancementFactor = min(2.0_RKIND,abs(cos(alpha)))    &
                        * sqrt(1.0_RKIND                    &
                        + (1.5_RKIND*langmuirNumber)**(-2)  &
                        + (5.4_RKIND*langmuirNumber)**(-4)) 

      !-----------------------------------------------------------------

   end subroutine ocn_stokes_drift_kpp_enhancement_factor!}}}

!***********************************************************************
!
!  routine ocn_stokes_drift_evaluate
!
!> \brief   Evaluate depth-dependent Stokes drift profile at a given depth
!> \author  Steven Brus
!> \date    June 2021 
!> \details
!>   This routine evaluates the depth-dependent Stokes drift profile at 
!>   a given depth using partitioned surface Stokes drift in 3, 4, or 6
!>   wavenumber bins from an active wave model. 
!
!-----------------------------------------------------------------------

   subroutine ocn_stokes_drift_evaluate(zDepth, zonalWavenumberCoeff, meridionalWavenumberCoeff, &
                                        stokesDriftZonal, stokesDriftMeridional)!{{{

      !-----------------------------------------------------------------
      ! Input variables
      !-----------------------------------------------------------------
      
      real(kind=RKIND), intent(in) :: zDepth
      real (kind=RKIND), dimension(:), intent(in) :: zonalWavenumberCoeff 
      real (kind=RKIND), dimension(:), intent(in) :: meridionalWavenumberCoeff

      !-----------------------------------------------------------------
      ! Output variables
      !-----------------------------------------------------------------

      real(kind=RKIND), intent(out) :: stokesDriftZonal
      real(kind=RKIND), intent(out) :: stokesDriftMeridional

      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------

      real(kind=RKIND) :: decayFactor
      integer :: i

      !-----------------------------------------------------------------

      stokesDriftZonal = 0.0_RKIND
      stokesDriftMeridional = 0.0_RKIND

      do i = 1,config_n_stokes_drift_wavenumber_partitions

         decayFactor = exp(2.0_RKIND*stokesDriftWavenumbers(i)*zDepth)
         stokesDriftZonal = stokesDriftZonal + zonalWavenumberCoeff(i)*decayFactor
         stokesDriftMeridional = stokesDriftMeridional + meridionalWavenumberCoeff(i)*decayFactor

      end do

      !-----------------------------------------------------------------

   end subroutine ocn_stokes_drift_evaluate!}}}


!***********************************************************************
!
!  routine ocn_stokes_drift_average
!
!> \brief   Averages Stokes drift profile over a specified depth
!> \author  Steven Brus
!> \date    June 2021 
!> \details
!>   This routine averages the Stokes drift profile over a specified 
!>   depth using partitioned surface Stokes drift in 3, 4, or 6 wavenumber
!>   bins from an active wave model. 
!
!-----------------------------------------------------------------------

   subroutine ocn_stokes_drift_average(zMLD, zonalWavenumberCoeff, meridionalWavenumberCoeff, &
                                       stokesDriftZonalMLAvg, stokesDriftMeridionalMLAvg)!{{{

      !-----------------------------------------------------------------
      ! Input variables
      !-----------------------------------------------------------------

      real (kind=RKIND), intent(in) :: zMLD
      real (kind=RKIND), dimension(:), intent(in) :: zonalWavenumberCoeff 
      real (kind=RKIND), dimension(:), intent(in) :: meridionalWavenumberCoeff

      !-----------------------------------------------------------------
      ! Output variables
      !-----------------------------------------------------------------

      real(kind=RKIND), intent(out) :: stokesDriftZonalMLAvg
      real(kind=RKIND), intent(out) :: stokesDriftMeridionalMLAvg

      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------

      real(kind=RKIND) :: integralFactor
      integer :: i

      !-----------------------------------------------------------------

      stokesDriftZonalMLAvg = 0.0_RKIND        
      stokesDriftMeridionalMLAvg = 0.0_RKIND

      do i = 1,config_n_stokes_drift_wavenumber_partitions
         integralFactor = (1.0_RKIND - exp(-2.0_RKIND*stokesDriftWavenumbers(i)*zMLD)) &
                        / (2.0_RKIND*stokesDriftWavenumbers(i))
         stokesDriftZonalMLAvg = stokesDriftZonalMLAvg + zonalWavenumberCoeff(i)*integralFactor
         stokesDriftMeridionalMLAvg = stokesDriftMeridionalMLAvg + meridionalWavenumberCoeff(i)*integralFactor
      end do

      stokesDriftZonalMLAvg = stokesDriftZonalMLAvg / zMLD
      stokesDriftMeridionalMLAvg = stokesDriftMeridionalMLAvg / zMLD

      !-----------------------------------------------------------------

   end subroutine ocn_stokes_drift_average!}}}

!***********************************************************************
!
!  routine ocn_stokes_drift_init
!
!> \brief   Initializes Stokes drift computations 
!> \author  Steven Brus
!> \date    June 2021
!> \details
!>  This routine initializes quantities related to Stokes drift-related
!>  computations. The wavenumber values are from a personal communication
!>  with Brandon Reichl at GFDL. This approach is used in Li et al. 2019,
!>  (https://doi.org/10.1029/2019MS001810) but these values are not
!>  documented there.
!
!----------------------------------------------------------------------

   subroutine ocn_stokes_drift_init(err)!{{{

      !-----------------------------------------------------------------
      ! Input variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! Output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err

      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------

      !*** initialize error flag

      err = 0

      if ((config_cvmix_kpp_use_active_wave .eqv. .true.) .and. (config_use_active_wave .eqv. .false.)) then
           call mpas_log_write(&
           "config_use_active_wave must be .true. for config_cvmix_kpp_use_active_wave = .true.")
           err = 1
      endif

      allocate(stokesDriftWavenumbers(config_n_stokes_drift_wavenumber_partitions))

      select case (config_n_stokes_drift_wavenumber_partitions)

         case(3)
           
           stokesDriftWavenumbers(1) = 0.04_RKIND
           stokesDriftWavenumbers(2) = 0.11_RKIND
           stokesDriftWavenumbers(3) = 0.3305_RKIND
           stokesDriftWavenumbers(4) = 0.0_RKIND
           stokesDriftWavenumbers(5) = 0.0_RKIND
           stokesDriftWavenumbers(6) = 0.0_RKIND

         case(4)

           stokesDriftWavenumbers(1) = 0.025_RKIND 
           stokesDriftWavenumbers(2) = 0.08_RKIND
           stokesDriftWavenumbers(3) = 0.16_RKIND
           stokesDriftWavenumbers(4) = 0.35_RKIND
           stokesDriftWavenumbers(5) = 0.0_RKIND
           stokesDriftWavenumbers(6) = 0.0_RKIND

         case(6)

           stokesDriftWavenumbers(1) = 0.01_RKIND
           stokesDriftWavenumbers(2) = 0.03_RKIND 
           stokesDriftWavenumbers(3) = 0.06_RKIND
           stokesDriftWavenumbers(4) = 0.10_RKIND
           stokesDriftWavenumbers(5) = 0.20_RKIND
           stokesDriftWavenumbers(6) = 0.35_RKIND

         case default
   
           call mpas_log_write(&
           "Invalid choice for config_n_stokes_drift_wavenumber_partitions. Choices are: 3, 4, and 6")
           err = 1

      end select

   !--------------------------------------------------------------------

   end subroutine ocn_stokes_drift_init!}}}


!***********************************************************************

end module ocn_stokes_drift

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
