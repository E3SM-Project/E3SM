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
             ocn_stokes_drift_average, &
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
!>   using partitioned surface Stokes drift in 6 wavenumber bins from 
!>   an active wave model. 
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

      integer :: iCell, k, i
      real(kind=RKIND), dimension(:,:), pointer :: stokesDriftZonal
      real(kind=RKIND), dimension(:,:), pointer :: stokesDriftMeridional
      real(kind=RKIND), dimension(:), pointer :: stokesDriftZonalWavenumber1, &
                                                 stokesDriftZonalWavenumber2, &
                                                 stokesDriftZonalWavenumber3, &
                                                 stokesDriftZonalWavenumber4, &
                                                 stokesDriftZonalWavenumber5, &
                                                 stokesDriftZonalWavenumber6
      real(kind=RKIND), dimension(:), pointer :: stokesDriftMeridionalWavenumber1, &
                                                 stokesDriftMeridionalWavenumber2, &
                                                 stokesDriftMeridionalWavenumber3, &
                                                 stokesDriftMeridionalWavenumber4, &
                                                 stokesDriftMeridionalWavenumber5, &
                                                 stokesDriftMeridionalWavenumber6
      real(kind=RKIND) :: decayFactor

      !-----------------------------------------------------------------

      call mpas_pool_get_array(forcingPool, 'stokesDriftZonal', stokesDriftZonal)
      call mpas_pool_get_array(forcingPool, 'stokesDriftMeridional', stokesDriftMeridional)
      call mpas_pool_get_array(forcingPool, 'stokesDriftZonalWavenumber1', stokesDriftZonalWavenumber1)
      call mpas_pool_get_array(forcingPool, 'stokesDriftMeridionalWavenumber1', stokesDriftMeridionalWavenumber1)
      call mpas_pool_get_array(forcingPool, 'stokesDriftZonalWavenumber2', stokesDriftZonalWavenumber2)
      call mpas_pool_get_array(forcingPool, 'stokesDriftMeridionalWavenumber2', stokesDriftMeridionalWavenumber2)
      call mpas_pool_get_array(forcingPool, 'stokesDriftZonalWavenumber3', stokesDriftZonalWavenumber3)
      call mpas_pool_get_array(forcingPool, 'stokesDriftMeridionalWavenumber3', stokesDriftMeridionalWavenumber3)
      call mpas_pool_get_array(forcingPool, 'stokesDriftZonalWavenumber4', stokesDriftZonalWavenumber4)
      call mpas_pool_get_array(forcingPool, 'stokesDriftMeridionalWavenumber4', stokesDriftMeridionalWavenumber4)
      call mpas_pool_get_array(forcingPool, 'stokesDriftZonalWavenumber5', stokesDriftZonalWavenumber5)
      call mpas_pool_get_array(forcingPool, 'stokesDriftZonalWavenumber5', stokesDriftZonalWavenumber5)
      call mpas_pool_get_array(forcingPool, 'stokesDriftMeridionalWavenumber6', stokesDriftMeridionalWavenumber6)
      call mpas_pool_get_array(forcingPool, 'stokesDriftMeridionalWavenumber6', stokesDriftMeridionalWavenumber6)

      do iCell = 1,nCellsAll

         do k = minLevelCell(iCell), maxLevelCell(iCell)
            stokesDriftZonal(k,iCell) = 0.0_RKIND
            stokesDriftMeridional(k,iCell) = 0.0_RKIND

            decayFactor = exp(2.0_RKIND*stokesDriftWavenumbers(1)*ZMid(k,iCell))
            stokesDriftZonal(k,iCell) = stokesDriftZonal(k,iCell) &
                                      + stokesDriftZonalWavenumber1(iCell)*decayFactor
            stokesDriftMeridional(k,iCell) = stokesDriftMeridional(k,iCell) &
                                           + stokesDriftMeridionalWavenumber1(iCell)*decayFactor

            decayFactor = exp(2.0_RKIND*stokesDriftWavenumbers(2)*ZMid(k,iCell))
            stokesDriftZonal(k,iCell) = stokesDriftZonal(k,iCell) &
                                      + stokesDriftZonalWavenumber2(iCell)*decayFactor
            stokesDriftMeridional(k,iCell) = stokesDriftMeridional(k,iCell) &
                                           + stokesDriftMeridionalWavenumber2(iCell)*decayFactor

            decayFactor = exp(2.0_RKIND*stokesDriftWavenumbers(3)*ZMid(k,iCell))
            stokesDriftZonal(k,iCell) = stokesDriftZonal(k,iCell) &
                                      + stokesDriftZonalWavenumber3(iCell)*decayFactor
            stokesDriftMeridional(k,iCell) = stokesDriftMeridional(k,iCell) &
                                           + stokesDriftMeridionalWavenumber3(iCell)*decayFactor

            if (config_cvmix_kpp_n_wavenumber_partitions > 3) then
               decayFactor = exp(2.0_RKIND*stokesDriftWavenumbers(4)*ZMid(k,iCell))
               stokesDriftZonal(k,iCell) = stokesDriftZonal(k,iCell) &
                                         + stokesDriftZonalWavenumber4(iCell)*decayFactor
               stokesDriftMeridional(k,iCell) = stokesDriftMeridional(k,iCell) &
                                              + stokesDriftMeridionalWavenumber4(iCell)*decayFactor
            end if

            if (config_cvmix_kpp_n_wavenumber_partitions > 4) then
               decayFactor = exp(2.0_RKIND*stokesDriftWavenumbers(5)*ZMid(k,iCell))
               stokesDriftZonal(k,iCell) = stokesDriftZonal(k,iCell) &
                                         + stokesDriftZonalWavenumber5(iCell)*decayFactor
               stokesDriftMeridional(k,iCell) = stokesDriftMeridional(k,iCell) &
                                              + stokesDriftMeridionalWavenumber5(iCell)*decayFactor

               decayFactor = exp(2.0_RKIND*stokesDriftWavenumbers(6)*ZMid(k,iCell))
               stokesDriftZonal(k,iCell) = stokesDriftZonal(k,iCell) &
                                         + stokesDriftZonalWavenumber6(iCell)*decayFactor
               stokesDriftMeridional(k,iCell) = stokesDriftMeridional(k,iCell) &
                                              + stokesDriftMeridionalWavenumber6(iCell)*decayFactor
            end if

         end do

      end do 

      !-----------------------------------------------------------------

   end subroutine ocn_stokes_drift_reconstruct!}}}

!***********************************************************************
!
!  routine ocn_stokes_drift_kpp_enhancement_factor
!
!> \brief Computes the Langmuir number based on Van Roekel et al. 2012 
!> \author  Steven Brus
!> \date    June 2021 
!> \details
!>   This routine computes the KPP enhancment factor based on a projected 
!>   Langmuir number using the Stokes drift from an active wave model. 
!
!-----------------------------------------------------------------------

   subroutine ocn_stokes_drift_kpp_enhancement_factor(forcingPool)!{{{

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

      real(kind=RKIND), dimension(:), pointer :: stokesDriftZonalMLAvg
      real(kind=RKIND), dimension(:), pointer :: stokesDriftMeridionalMLAvg
      real(kind=RKIND), dimension(:), pointer :: stokesDriftKPPEnhancementFactor
      real(kind=RKIND), dimension(:), pointer :: significantWaveHeight

      !-----------------------------------------------------------------

      call mpas_pool_get_array(forcingPool, 'stokesDriftZonalMLAvg', stokesDriftZonalMLAvg)
      call mpas_pool_get_array(forcingPool, 'stokesDriftMeridionalMLAvg', stokesDriftMeridionalMLAvg)
      call mpas_pool_get_array(forcingPool, 'stokesDriftKPPEnhancementFactor', stokesDriftKPPEnhancementFactor)
      call mpas_pool_get_array(forcingPool, 'significantWaveHeight', significantWaveHeight)


      !-----------------------------------------------------------------

   end subroutine ocn_stokes_drift_kpp_enhancement_factor!}}}

!***********************************************************************
!
!  routine ocn_stokes_drift_average
!
!> \brief   Averages Stokes drift profile over a specified depth
!> \author  Steven Brus
!> \date    June 2021 
!> \details
!>   This routine averages the Stokes drift profile over a specified 
!>   depth using partitioned surface Stokes drift in 6 wavenumber bins
!>   from an active wave model. 
!
!-----------------------------------------------------------------------

   subroutine ocn_stokes_drift_average(statePool,forcingPool,fractionOfMLD,timeLevelIn)!{{{

      !-----------------------------------------------------------------
      ! Input variables
      !-----------------------------------------------------------------

      type(mpas_pool_type), intent(in) :: &
         statePool

      type(mpas_pool_type), intent(inout) :: &
         forcingPool

      real, intent(in) :: fractionOfMLD

      integer, intent(in), optional :: timeLevelIn


      !-----------------------------------------------------------------
      ! Output variables
      !-----------------------------------------------------------------


      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------

      integer :: iCell, k
      integer :: timeLevel

      real(kind=RKIND), dimension(:), pointer :: ssh
      real(kind=RKIND), dimension(:), pointer :: stokesDriftZonalMLAvg
      real(kind=RKIND), dimension(:), pointer :: stokesDriftMeridionalMLAvg
      real(kind=RKIND), dimension(:), pointer :: stokesDriftZonalWavenumber1, &
                                                 stokesDriftZonalWavenumber2, &
                                                 stokesDriftZonalWavenumber3, &
                                                 stokesDriftZonalWavenumber4, &
                                                 stokesDriftZonalWavenumber5, &
                                                 stokesDriftZonalWavenumber6
      real(kind=RKIND), dimension(:), pointer :: stokesDriftMeridionalWavenumber1, &
                                                 stokesDriftMeridionalWavenumber2, &
                                                 stokesDriftMeridionalWavenumber3, &
                                                 stokesDriftMeridionalWavenumber4, &
                                                 stokesDriftMeridionalWavenumber5, &
                                                 stokesDriftMeridionalWavenumber6
      real(kind=RKIND) :: integralFactor
      real(kind=RKIND) :: zMLD

      !-----------------------------------------------------------------

      if (present(timeLevelIn)) then
         timeLevel = timeLevelIn
      else
         timeLevel = 1
      end if

      call mpas_pool_get_array(statePool, 'ssh', ssh, timeLevel)

      call mpas_pool_get_array(forcingPool, 'stokesDriftZonalMLAvg', stokesDriftZonalMLAvg)
      call mpas_pool_get_array(forcingPool, 'stokesDriftMeridionalMLAvg', stokesDriftMeridionalMLAvg)
      call mpas_pool_get_array(forcingPool, 'stokesDriftZonalWavenumber1', stokesDriftZonalWavenumber1)
      call mpas_pool_get_array(forcingPool, 'stokesDriftMeridionalWavenumber1', stokesDriftMeridionalWavenumber1)
      call mpas_pool_get_array(forcingPool, 'stokesDriftZonalWavenumber2', stokesDriftZonalWavenumber2)
      call mpas_pool_get_array(forcingPool, 'stokesDriftMeridionalWavenumber2', stokesDriftMeridionalWavenumber2)
      call mpas_pool_get_array(forcingPool, 'stokesDriftZonalWavenumber3', stokesDriftZonalWavenumber3)
      call mpas_pool_get_array(forcingPool, 'stokesDriftMeridionalWavenumber3', stokesDriftMeridionalWavenumber3)
      call mpas_pool_get_array(forcingPool, 'stokesDriftZonalWavenumber4', stokesDriftZonalWavenumber4)
      call mpas_pool_get_array(forcingPool, 'stokesDriftMeridionalWavenumber4', stokesDriftMeridionalWavenumber4)
      call mpas_pool_get_array(forcingPool, 'stokesDriftZonalWavenumber5', stokesDriftZonalWavenumber5)
      call mpas_pool_get_array(forcingPool, 'stokesDriftZonalWavenumber5', stokesDriftZonalWavenumber5)
      call mpas_pool_get_array(forcingPool, 'stokesDriftMeridionalWavenumber6', stokesDriftMeridionalWavenumber6)
      call mpas_pool_get_array(forcingPool, 'stokesDriftMeridionalWavenumber6', stokesDriftMeridionalWavenumber6)

      do iCell = 1,nCellsAll

        k =  max(minLevelCell(iCell), indMLD(iCell))
        zMLD = fractionOfMLD*(ssh(iCell) - zMid(k, iCell))

        stokesDriftZonalMLAvg(iCell) = 0.0_RKIND        
        stokesDriftMeridionalMLAvg(iCell) = 0.0_RKIND
        
        integralFactor = (1.0_RKIND - exp(-2.0_RKIND*stokesDriftWavenumbers(1)*zMLD)) &
                       / (2.0_RKIND*stokesDriftWavenumbers(1))
        stokesDriftZonalMLAvg(iCell) = stokesDriftZonalMLAvg(iCell) &
                                     + stokesDriftZonalWavenumber1(iCell)*integralFactor
        stokesDriftMeridionalMLAvg(iCell) = stokesDriftMeridionalMLAvg(iCell) &
                                          + stokesDriftMeridionalWavenumber1(iCell)*integralFactor
        
        integralFactor = (1.0_RKIND - exp(-2.0_RKIND*stokesDriftWavenumbers(2)*zMLD)) &
                       / (2.0_RKIND*stokesDriftWavenumbers(2))
        stokesDriftZonalMLAvg(iCell) = stokesDriftZonalMLAvg(iCell) &
                                     + stokesDriftZonalWavenumber2(iCell)*integralFactor
        stokesDriftMeridionalMLAvg(iCell) = stokesDriftMeridionalMLAvg(iCell) &
                                          + stokesDriftMeridionalWavenumber2(iCell)*integralFactor
        
        integralFactor = (1.0_RKIND - exp(-2.0_RKIND*stokesDriftWavenumbers(3)*zMLD)) &
                       / (2.0_RKIND*stokesDriftWavenumbers(3))
        stokesDriftZonalMLAvg(iCell) = stokesDriftZonalMLAvg(iCell) &
                                     + stokesDriftZonalWavenumber3(iCell)*integralFactor
        stokesDriftMeridionalMLAvg(iCell) = stokesDriftMeridionalMLAvg(iCell) &
                                          + stokesDriftMeridionalWavenumber3(iCell)*integralFactor
        
        if (config_cvmix_kpp_n_wavenumber_partitions > 3) then
           integralFactor = (1.0_RKIND - exp(-2.0_RKIND*stokesDriftWavenumbers(4)*zMLD)) &
                          / (2.0_RKIND*stokesDriftWavenumbers(4))
           stokesDriftZonalMLAvg(iCell) = stokesDriftZonalMLAvg(iCell) &
                                        + stokesDriftZonalWavenumber4(iCell)*integralFactor
           stokesDriftMeridionalMLAvg(iCell) = stokesDriftMeridionalMLAvg(iCell) &
                                             + stokesDriftMeridionalWavenumber4(iCell)*integralFactor
        end if
        
        if (config_cvmix_kpp_n_wavenumber_partitions > 4) then
           integralFactor = (1.0_RKIND - exp(-2.0_RKIND*stokesDriftWavenumbers(5)*zMLD)) &
                          / (2.0_RKIND*stokesDriftWavenumbers(5))
           stokesDriftZonalMLAvg(iCell) = stokesDriftZonalMLAvg(iCell) &
                                        + stokesDriftZonalWavenumber5(iCell)*integralFactor
           stokesDriftMeridionalMLAvg(iCell) = stokesDriftMeridionalMLAvg(iCell) &
                                             + stokesDriftMeridionalWavenumber5(iCell)*integralFactor
           
           integralFactor = (1.0_RKIND - exp(-2.0_RKIND*stokesDriftWavenumbers(6)*zMLD)) &
                          / (2.0_RKIND*stokesDriftWavenumbers(6))
           stokesDriftZonalMLAvg(iCell) = stokesDriftZonalMLAvg(iCell) &
                                        + stokesDriftZonalWavenumber6(iCell)*integralFactor
           stokesDriftMeridionalMLAvg(iCell) = stokesDriftMeridionalMLAvg(iCell) &
                                             + stokesDriftMeridionalWavenumber6(iCell)*integralFactor
        end if

        stokesDriftZonalMLAvg(iCell) = stokesDriftZonalMLAvg(iCell) / zMLD
        stokesDriftMeridionalMLAvg(iCell) = stokesDriftMeridionalMLAvg(iCell) / zMLD

      end do

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
!>  computations
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

      allocate(stokesDriftWavenumbers(config_cvmix_kpp_n_wavenumber_partitions))

      select case (config_cvmix_kpp_n_wavenumber_partitions)

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
           "Invalid choice for config_cvmix_kpp_n_wavenumber_partitions. Choices are: 3, 4, and 6")
           err = 1

      end select

   !--------------------------------------------------------------------

   end subroutine ocn_stokes_drift_init!}}}


!***********************************************************************

end module ocn_stokes_drift

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
