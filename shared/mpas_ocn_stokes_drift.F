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

   real (kind=RKIND), dimension(:), allocatable :: waveStokesDriftWavenumbers


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
      real(kind=RKIND), dimension(:,:), pointer :: waveStokesDriftZonal
      real(kind=RKIND), dimension(:,:), pointer :: waveStokesDriftMeridional
      real(kind=RKIND), dimension(:), pointer :: waveStokesDriftZonalWavenumber1, &
                                                 waveStokesDriftZonalWavenumber2, &
                                                 waveStokesDriftZonalWavenumber3, &
                                                 waveStokesDriftZonalWavenumber4, &
                                                 waveStokesDriftZonalWavenumber5, &
                                                 waveStokesDriftZonalWavenumber6
      real(kind=RKIND), dimension(:), pointer :: waveStokesDriftMeridionalWavenumber1, &
                                                 waveStokesDriftMeridionalWavenumber2, &
                                                 waveStokesDriftMeridionalWavenumber3, &
                                                 waveStokesDriftMeridionalWavenumber4, &
                                                 waveStokesDriftMeridionalWavenumber5, &
                                                 waveStokesDriftMeridionalWavenumber6
      real(kind=RKIND) :: decayFactor

      !-----------------------------------------------------------------

      call mpas_pool_get_array(forcingPool, 'waveStokesDriftZonal', waveStokesDriftZonal)
      call mpas_pool_get_array(forcingPool, 'waveStokesDriftMeridional', waveStokesDriftMeridional)
      call mpas_pool_get_array(forcingPool, 'waveStokesDriftZonalWavenumber1', waveStokesDriftZonalWavenumber1)
      call mpas_pool_get_array(forcingPool, 'waveStokesDriftMeridionalWavenumber1', waveStokesDriftMeridionalWavenumber1)
      call mpas_pool_get_array(forcingPool, 'waveStokesDriftZonalWavenumber2', waveStokesDriftZonalWavenumber2)
      call mpas_pool_get_array(forcingPool, 'waveStokesDriftMeridionalWavenumber2', waveStokesDriftMeridionalWavenumber2)
      call mpas_pool_get_array(forcingPool, 'waveStokesDriftZonalWavenumber3', waveStokesDriftZonalWavenumber3)
      call mpas_pool_get_array(forcingPool, 'waveStokesDriftMeridionalWavenumber3', waveStokesDriftMeridionalWavenumber3)
      call mpas_pool_get_array(forcingPool, 'waveStokesDriftZonalWavenumber4', waveStokesDriftZonalWavenumber4)
      call mpas_pool_get_array(forcingPool, 'waveStokesDriftMeridionalWavenumber4', waveStokesDriftMeridionalWavenumber4)
      call mpas_pool_get_array(forcingPool, 'waveStokesDriftZonalWavenumber5', waveStokesDriftZonalWavenumber5)
      call mpas_pool_get_array(forcingPool, 'waveStokesDriftZonalWavenumber5', waveStokesDriftZonalWavenumber5)
      call mpas_pool_get_array(forcingPool, 'waveStokesDriftMeridionalWavenumber6', waveStokesDriftMeridionalWavenumber6)
      call mpas_pool_get_array(forcingPool, 'waveStokesDriftMeridionalWavenumber6', waveStokesDriftMeridionalWavenumber6)

      do iCell = 1,nCellsAll

         do k = minLevelCell(iCell), maxLevelCell(iCell)
            waveStokesDriftZonal(k,iCell) = 0.0_RKIND
            waveStokesDriftMeridional(k,iCell) = 0.0_RKIND

            decayFactor = exp(2.0_RKIND*waveStokesDriftWavenumbers(1)*ZMid(k,iCell))
            waveStokesDriftZonal(k,iCell) = waveStokesDriftZonal(k,iCell) &
                                          + waveStokesDriftZonalWavenumber1(iCell)*decayFactor
            waveStokesDriftMeridional(k,iCell) = waveStokesDriftMeridional(k,iCell) &
                                               + waveStokesDriftMeridionalWavenumber1(iCell)*decayFactor

            decayFactor = exp(2.0_RKIND*waveStokesDriftWavenumbers(2)*ZMid(k,iCell))
            waveStokesDriftZonal(k,iCell) = waveStokesDriftZonal(k,iCell) &
                                          + waveStokesDriftZonalWavenumber2(iCell)*decayFactor
            waveStokesDriftMeridional(k,iCell) = waveStokesDriftMeridional(k,iCell) &
                                               + waveStokesDriftMeridionalWavenumber2(iCell)*decayFactor

            decayFactor = exp(2.0_RKIND*waveStokesDriftWavenumbers(3)*ZMid(k,iCell))
            waveStokesDriftZonal(k,iCell) = waveStokesDriftZonal(k,iCell) &
                                          + waveStokesDriftZonalWavenumber3(iCell)*decayFactor
            waveStokesDriftMeridional(k,iCell) = waveStokesDriftMeridional(k,iCell) &
                                               + waveStokesDriftMeridionalWavenumber3(iCell)*decayFactor

            if (config_cvmix_kpp_n_wavenumber_partitions > 3) then
               decayFactor = exp(2.0_RKIND*waveStokesDriftWavenumbers(4)*ZMid(k,iCell))
               waveStokesDriftZonal(k,iCell) = waveStokesDriftZonal(k,iCell) &
                                             + waveStokesDriftZonalWavenumber4(iCell)*decayFactor
               waveStokesDriftMeridional(k,iCell) = waveStokesDriftMeridional(k,iCell) &
                                                 + waveStokesDriftMeridionalWavenumber4(iCell)*decayFactor
            end if

            if (config_cvmix_kpp_n_wavenumber_partitions > 4) then
               decayFactor = exp(2.0_RKIND*waveStokesDriftWavenumbers(5)*ZMid(k,iCell))
               waveStokesDriftZonal(k,iCell) = waveStokesDriftZonal(k,iCell) &
                                             + waveStokesDriftZonalWavenumber5(iCell)*decayFactor
               waveStokesDriftMeridional(k,iCell) = waveStokesDriftMeridional(k,iCell) &
                                                  + waveStokesDriftMeridionalWavenumber5(iCell)*decayFactor

               decayFactor = exp(2.0_RKIND*waveStokesDriftWavenumbers(6)*ZMid(k,iCell))
               waveStokesDriftZonal(k,iCell) = waveStokesDriftZonal(k,iCell) &
                                             + waveStokesDriftZonalWavenumber6(iCell)*decayFactor
               waveStokesDriftMeridional(k,iCell) = waveStokesDriftMeridional(k,iCell) &
                                                  + waveStokesDriftMeridionalWavenumber6(iCell)*decayFactor
            end if

         end do

      end do 

      !-----------------------------------------------------------------

   end subroutine ocn_stokes_drift_reconstruct!}}}

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
      real(kind=RKIND), dimension(:), pointer :: waveStokesDriftZonalMLAvg
      real(kind=RKIND), dimension(:), pointer :: waveStokesDriftMeridionalMLAvg
      real(kind=RKIND), dimension(:), pointer :: waveStokesDriftZonalWavenumber1, &
                                                 waveStokesDriftZonalWavenumber2, &
                                                 waveStokesDriftZonalWavenumber3, &
                                                 waveStokesDriftZonalWavenumber4, &
                                                 waveStokesDriftZonalWavenumber5, &
                                                 waveStokesDriftZonalWavenumber6
      real(kind=RKIND), dimension(:), pointer :: waveStokesDriftMeridionalWavenumber1, &
                                                 waveStokesDriftMeridionalWavenumber2, &
                                                 waveStokesDriftMeridionalWavenumber3, &
                                                 waveStokesDriftMeridionalWavenumber4, &
                                                 waveStokesDriftMeridionalWavenumber5, &
                                                 waveStokesDriftMeridionalWavenumber6
      real(kind=RKIND) :: integralFactor
      real(kind=RKIND) :: zMLD

      !-----------------------------------------------------------------

      if (present(timeLevelIn)) then
         timeLevel = timeLevelIn
      else
         timeLevel = 1
      end if

      call mpas_pool_get_array(statePool, 'ssh', ssh, timeLevel)

      call mpas_pool_get_array(forcingPool, 'waveStokesDriftZonalMLAvg', waveStokesDriftZonalMLAvg)
      call mpas_pool_get_array(forcingPool, 'waveStokesDriftMeridionalMLAvg', waveStokesDriftMeridionalMLAvg)
      call mpas_pool_get_array(forcingPool, 'waveStokesDriftZonalWavenumber1', waveStokesDriftZonalWavenumber1)
      call mpas_pool_get_array(forcingPool, 'waveStokesDriftMeridionalWavenumber1', waveStokesDriftMeridionalWavenumber1)
      call mpas_pool_get_array(forcingPool, 'waveStokesDriftZonalWavenumber2', waveStokesDriftZonalWavenumber2)
      call mpas_pool_get_array(forcingPool, 'waveStokesDriftMeridionalWavenumber2', waveStokesDriftMeridionalWavenumber2)
      call mpas_pool_get_array(forcingPool, 'waveStokesDriftZonalWavenumber3', waveStokesDriftZonalWavenumber3)
      call mpas_pool_get_array(forcingPool, 'waveStokesDriftMeridionalWavenumber3', waveStokesDriftMeridionalWavenumber3)
      call mpas_pool_get_array(forcingPool, 'waveStokesDriftZonalWavenumber4', waveStokesDriftZonalWavenumber4)
      call mpas_pool_get_array(forcingPool, 'waveStokesDriftMeridionalWavenumber4', waveStokesDriftMeridionalWavenumber4)
      call mpas_pool_get_array(forcingPool, 'waveStokesDriftZonalWavenumber5', waveStokesDriftZonalWavenumber5)
      call mpas_pool_get_array(forcingPool, 'waveStokesDriftZonalWavenumber5', waveStokesDriftZonalWavenumber5)
      call mpas_pool_get_array(forcingPool, 'waveStokesDriftMeridionalWavenumber6', waveStokesDriftMeridionalWavenumber6)
      call mpas_pool_get_array(forcingPool, 'waveStokesDriftMeridionalWavenumber6', waveStokesDriftMeridionalWavenumber6)

      do iCell = 1,nCellsAll

        k =  max(minLevelCell(iCell), indMLD(iCell))
        zMLD = fractionOfMLD*(ssh(iCell) - zMid(k, iCell))

        waveStokesDriftZonalMLAvg(iCell) = 0.0_RKIND        
        waveStokesDriftMeridionalMLAvg(iCell) = 0.0_RKIND
        
        integralFactor = (1.0_RKIND - exp(-2.0_RKIND*waveStokesDriftWavenumbers(1)*zMLD)) &
                       / (2.0_RKIND*waveStokesDriftWavenumbers(1))
        waveStokesDriftZonalMLAvg(iCell) = waveStokesDriftZonalMLAvg(iCell) &
                                         + waveStokesDriftZonalWavenumber1(iCell)*integralFactor
        waveStokesDriftMeridionalMLAvg(iCell) = waveStokesDriftMeridionalMLAvg(iCell) &
                                              + waveStokesDriftMeridionalWavenumber1(iCell)*integralFactor
        
        integralFactor = (1.0_RKIND - exp(-2.0_RKIND*waveStokesDriftWavenumbers(2)*zMLD)) &
                       / (2.0_RKIND*waveStokesDriftWavenumbers(2))
        waveStokesDriftZonalMLAvg(iCell) = waveStokesDriftZonalMLAvg(iCell) &
                                         + waveStokesDriftZonalWavenumber2(iCell)*integralFactor
        waveStokesDriftMeridionalMLAvg(iCell) = waveStokesDriftMeridionalMLAvg(iCell) &
                                              + waveStokesDriftMeridionalWavenumber2(iCell)*integralFactor
        
        integralFactor = (1.0_RKIND - exp(-2.0_RKIND*waveStokesDriftWavenumbers(3)*zMLD)) &
                       / (2.0_RKIND*waveStokesDriftWavenumbers(3))
        waveStokesDriftZonalMLAvg(iCell) = waveStokesDriftZonalMLAvg(iCell) &
                                         + waveStokesDriftZonalWavenumber3(iCell)*integralFactor
        waveStokesDriftMeridionalMLAvg(iCell) = waveStokesDriftMeridionalMLAvg(iCell) &
                                              + waveStokesDriftMeridionalWavenumber3(iCell)*integralFactor
        
        if (config_cvmix_kpp_n_wavenumber_partitions > 3) then
           integralFactor = (1.0_RKIND - exp(-2.0_RKIND*waveStokesDriftWavenumbers(4)*zMLD)) &
                          / (2.0_RKIND*waveStokesDriftWavenumbers(4))
           waveStokesDriftZonalMLAvg(iCell) = waveStokesDriftZonalMLAvg(iCell) &
                                            + waveStokesDriftZonalWavenumber4(iCell)*integralFactor
           waveStokesDriftMeridionalMLAvg(iCell) = waveStokesDriftMeridionalMLAvg(iCell) &
                                                 + waveStokesDriftMeridionalWavenumber4(iCell)*integralFactor
        end if
        
        if (config_cvmix_kpp_n_wavenumber_partitions > 4) then
           integralFactor = (1.0_RKIND - exp(-2.0_RKIND*waveStokesDriftWavenumbers(5)*zMLD)) &
                          / (2.0_RKIND*waveStokesDriftWavenumbers(5))
           waveStokesDriftZonalMLAvg(iCell) = waveStokesDriftZonalMLAvg(iCell) &
                                            + waveStokesDriftZonalWavenumber5(iCell)*integralFactor
           waveStokesDriftMeridionalMLAvg(iCell) = waveStokesDriftMeridionalMLAvg(iCell) &
                                                 + waveStokesDriftMeridionalWavenumber5(iCell)*integralFactor
           
           integralFactor = (1.0_RKIND - exp(-2.0_RKIND*waveStokesDriftWavenumbers(6)*zMLD)) &
                          / (2.0_RKIND*waveStokesDriftWavenumbers(6))
           waveStokesDriftZonalMLAvg(iCell) = waveStokesDriftZonalMLAvg(iCell) &
                                            + waveStokesDriftZonalWavenumber6(iCell)*integralFactor
           waveStokesDriftMeridionalMLAvg(iCell) = waveStokesDriftMeridionalMLAvg(iCell) &
                                                 + waveStokesDriftMeridionalWavenumber6(iCell)*integralFactor
        end if

        waveStokesDriftZonalMLAvg(iCell) = waveStokesDriftZonalMLAvg(iCell) / zMLD
        waveStokesDriftMeridionalMLAvg(iCell) = waveStokesDriftMeridionalMLAvg(iCell) / zMLD

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

      allocate(waveStokesDriftWavenumbers(config_cvmix_kpp_n_wavenumber_partitions))

      select case (config_cvmix_kpp_n_wavenumber_partitions)

         case(3)
           
           waveStokesDriftWavenumbers(1) = 0.04_RKIND
           waveStokesDriftWavenumbers(2) = 0.11_RKIND
           waveStokesDriftWavenumbers(3) = 0.3305_RKIND
           waveStokesDriftWavenumbers(4) = 0.0_RKIND
           waveStokesDriftWavenumbers(5) = 0.0_RKIND
           waveStokesDriftWavenumbers(6) = 0.0_RKIND

         case(4)

           waveStokesDriftWavenumbers(1) = 0.025_RKIND 
           waveStokesDriftWavenumbers(2) = 0.08_RKIND
           waveStokesDriftWavenumbers(3) = 0.16_RKIND
           waveStokesDriftWavenumbers(4) = 0.35_RKIND
           waveStokesDriftWavenumbers(5) = 0.0_RKIND
           waveStokesDriftWavenumbers(6) = 0.0_RKIND

         case(6)

           waveStokesDriftWavenumbers(1) = 0.01_RKIND
           waveStokesDriftWavenumbers(2) = 0.03_RKIND 
           waveStokesDriftWavenumbers(3) = 0.06_RKIND
           waveStokesDriftWavenumbers(4) = 0.10_RKIND
           waveStokesDriftWavenumbers(5) = 0.20_RKIND
           waveStokesDriftWavenumbers(6) = 0.35_RKIND

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
