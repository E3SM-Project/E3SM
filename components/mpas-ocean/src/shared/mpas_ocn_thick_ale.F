! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_thick_ale
!
!> \brief MPAS ocean ALE thickness driver
!> \author Mark Petersen
!> \date   August 2013
!> \details
!>  This module contains the routines for computing ALE thickness.
!
!-----------------------------------------------------------------------

module ocn_thick_ale

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_timer

   use ocn_constants
   use ocn_config
   use ocn_mesh

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

   public :: ocn_ALE_thickness, &
             ocn_thick_ale_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   ! Choices for ALE thickness proportionality

   integer :: &
      ALEthickProportionality

   integer, parameter :: &
      ALEthickProportionThickTimesWgts = 1, &
      ALEthickProportionWgtsOnly       = 2

   ! Variables for adjusting thickness based on min, max parameters

   logical :: &
      useMinMaxThick ! flag for using min/max adjustment

   real (kind=RKIND) :: &
      minThick,     &! minimum thickness
      maxThickFact   ! factor for limiting max thickness

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_ALE_thickness
!
!> \brief   Computes desired ALE thickness at new time
!> \author  Mark Petersen
!> \date    August 2013
!> \details
!>  This routine computes the desired Arbitrary Lagrangian-Eulerian (ALE)
!>  thickness at the new time. It uses the ALE formulation, and includes
!>  contributions from SSH variations (z-star), high-frequency divergence
!>  (z-tilde), and imposes a minimum layer thickness.
!
!-----------------------------------------------------------------------

   subroutine ocn_ALE_thickness(verticalMeshPool, SSH, ALE_thickness, &
                                err, newHighFreqThickness)!{{{

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         verticalMeshPool     !< [in] vertical mesh information

      real (kind=RKIND), dimension(:), intent(in) :: &
         SSH                  !< [in] sea surface height

      real (kind=RKIND), dimension(:,:), intent(in), optional :: &
         newHighFreqThickness !< [in] high freq thickness alters ALE thickness.

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(out) :: &
         ALE_thickness            !< [out] desired thickness at new time

      integer, intent(out) :: err !< [out] error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      integer :: &
         iCell, k,    &! loop indices for cells, depth
         kMax, kMin,  &! indices for min, max active levels
         nCells        ! number of cells

      real (kind=RKIND) :: &
         weightSum,    &! sum of weights in vertical
         thicknessSum, &! total thickness
         remainder,    &! track remainder in mix/max alteration
         newThickness   ! temp used during min/max adjustment

      real (kind=RKIND), dimension(:), allocatable :: &
         prelim_ALE_thickness,   & ! ALE thickness at new time
         min_ALE_thickness_down, & ! thickness alteration in min/max adj
         min_ALE_thickness_up      ! thickness alteration in min/max adj

      real (kind=RKIND), dimension(:,:), pointer :: &
         restingThickness   ! Layer thickness when the ocean is at rest
                            ! i.e. without SSH or internal perturbations.

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      err = 0

      !TEMP: acc directives temporarily disabled here for bit-for-bit
      !      validation - will enable these on subsequent PR
      !      must also update host since ssh data on device - this
      !      will also be removed
      !$acc update host (ssh)

      ! Retrieve or allocate needed variables

      call mpas_pool_get_array(verticalMeshPool, 'restingThickness', &
                                                  restingThickness)
      !xacc enter data copyin(restingThickness) 

      nCells = nCellsHalo(1)

      !
      ! ALE thickness alteration due to SSH (z-star)
      !
      select case(ALEthickProportionality)

      ! restingThickness_times_weights
      case (ALEthickProportionThickTimesWgts)

#ifdef MPAS_OPENACC
         !xacc parallel loop &
         !xacc    present(ALE_thickness, SSH, restingThickness, &
         !xacc            minLevelCell, maxLevelCell, &
         !xacc            vertCoordMovementWeights) &
         !xacc    private(k, kMin, kMax, thicknessSum)
#else
         !$omp parallel
         !$omp do schedule(runtime) &
         !$omp    private(k, kMin, kMax, thicknessSum)
#endif
         do iCell = 1, nCells
            kMax = maxLevelCell(iCell)
            kMin = minLevelCell(iCell)
   
            thicknessSum = 1e-14_RKIND
            do k = kMin, kMax
               thicknessSum = thicknessSum &
                            + vertCoordMovementWeights(k) &
                            * restingThickness(k,iCell)
            end do
   
            ! Note that restingThickness is nonzero, and remaining
            ! terms are perturbations about zero.
            ! This is equation 4 and 6 in Petersen et al 2015,
            ! but with eqn 6
            do k = kMin, kMax
               ALE_thickness(k,iCell) = restingThickness(k,iCell) &
                  + (SSH(iCell)*vertCoordMovementWeights(k)* &
                     restingThickness(k,iCell) )/thicknessSum
            end do
         enddo
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif
   
      ! weights_only
      case (ALEthickProportionWgtsOnly)

#ifdef MPAS_OPENACC
         !xacc parallel loop &
         !xacc    present(ALE_thickness, restingThickness, ssh, &
         !xacc            vertCoordMovementWeights, &
         !xacc            minLevelCell, maxLevelCell) &
         !xacc    private(k, kMin, kMax, weightSum)
#else
         !$omp parallel
         !$omp do schedule(runtime) &
         !$omp    private(k, kMin, kMax, weightSum)
#endif
         do iCell = 1, nCells
            kMax = maxLevelCell(iCell)
            kMin = minLevelCell(iCell)
   
            weightSum = 1e-14_RKIND
            do k = kMin, kMax
               weightSum = weightSum + vertCoordMovementWeights(k) 
            end do
   
            do k = kMin, kMax
               ! Using this, we must require that the
               ! sum(restingThickness(k, iCell))
               ! summed over k is equal to bottomDepth.
               ! This is equation 4 and 6 in Petersen et al 2015,
               ! but with eqn 6 altered so only the W_k weights
               ! are used. The resting thickness shown in eqn 6 is
               ! not included here.
               ALE_thickness(k,iCell) = restingThickness(k,iCell) &
                                      + ssh(iCell)* &
                            vertCoordMovementWeights(k)/weightSum
            end do
         enddo
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif

      case default
         ! Checked on init
      end select ! ALE thickness proportionality

      ! thickness filtering active
      if (present(newHighFreqThickness)) then

#ifdef MPAS_OPENACC
         !xacc parallel loop &
         !xacc    present(ALE_thickness, newHighFreqThickness, &
         !xacc            minLevelCell, maxLevelCell) &
         !xacc    private(k, kMin, kMax)
#else
         !$omp parallel
         !$omp do schedule(runtime) private(k, kMin, kMax)
#endif
         do iCell = 1, nCells
            kMax = maxLevelCell(iCell)
            kMin = minLevelCell(iCell)

            do k=kMin,kMax
               ALE_thickness(k,iCell) = ALE_thickness(k,iCell) &
                                      + newHighFreqThickness(k,iCell)
            enddo
         enddo
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif

      endif

      !
      ! ALE thickness alteration due to minimum and maximum thickness
      !
      if (useMinMaxThick) then

         allocate(prelim_ALE_thickness  (nVertLevels), &
                  min_ALE_thickness_down(nVertLevels), &
                  min_ALE_thickness_up  (nVertLevels))
         !xacc enter data create(prelim_ALE_thickness, &
         !xacc                   min_ALE_thickness_down, &
         !xacc                   min_ALE_thickness_up)

#ifdef MPAS_OPENACC
         !xacc parallel loop &
         !xacc    present(ALE_thickness, restingThickness, &
         !xacc            minLevelCell, maxLevelCell) &
         !xacc    private(k, kMin, kMax, remainder, newThickness, &
         !xacc            prelim_ALE_thickness, &
         !xacc            min_ALE_thickness_down, min_ALE_thickness_up)
#else
         !$omp parallel
         !$omp do schedule(runtime) &
         !$omp    private(k, kMin, kMax, remainder, newThickness, &
         !$omp            prelim_ALE_thickness, &
         !$omp            min_ALE_thickness_down, min_ALE_thickness_up)
#endif
         do iCell = 1, nCells
            kMax = maxLevelCell(iCell)
            kMin = minLevelCell(iCell)

            ! go down the column:
            do k = kMin, kMax
               prelim_ALE_thickness(k) = ALE_thickness(k,iCell)
            end do

            remainder = 0.0_RKIND
            do k = kMin, kMax
               newThickness = max(minThick, &
                            min(prelim_ALE_thickness(k) + remainder, &
                                maxThickFact*restingThickness(k,iCell)))
               min_ALE_thickness_down(k) = newThickness &
                                         - prelim_ALE_thickness(k)
               remainder = remainder - min_ALE_thickness_down(k)
            end do

            ! go back up the column:
            min_ALE_thickness_up(kMax) = 0.0_RKIND
            do k=kMin,kMax
               prelim_ALE_thickness(k) = prelim_ALE_thickness(k) &
                                       + min_ALE_thickness_down(k)
            end do
            do k = kMax-1, kMin, -1
               newThickness = max(minThick, &
                         min(prelim_ALE_thickness(k) + remainder, &
                             maxThickFact * restingThickness(k,iCell)))
               min_ALE_thickness_up(k) = newThickness &
                                       - prelim_ALE_thickness(k)
               remainder = remainder - min_ALE_thickness_up(k)
            end do
            min_ALE_thickness_up(kMin) = min_ALE_thickness_up(kMin) &
                                       + remainder

            do k=kMin,kMax
               ALE_thickness(k,iCell) = ALE_thickness(k,iCell) &
                                      + min_ALE_thickness_down(k) &
                                      + min_ALE_thickness_up(k)

            enddo
         enddo
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif

         !xacc exit data delete(prelim_ALE_thickness, &
         !xacc                  min_ALE_thickness_down, &
         !xacc                  min_ALE_thickness_up)
         deallocate(prelim_ALE_thickness, &
                    min_ALE_thickness_down, &
                    min_ALE_thickness_up)

      endif ! useMinMaxThick

      !xacc exit data delete(restingThickness) 

      !TEMP: move result to device - this will be unnecessary once
      !      acc directives enabled
      !$acc update device(ALE_thickness)

   !--------------------------------------------------------------------

   end subroutine ocn_ALE_thickness!}}}

!***********************************************************************
!
!  routine ocn_thick_ale_init
!
!> \brief   Initializes flags used within ALE thickness computation.
!> \author  Mark Petersen
!> \date    August 2013
!> \details
!>  This routine initializes flags related to thickness computation in
!>  the Arbitrary Lagrangian-Eulerian (ALE) formulation.
!
!-----------------------------------------------------------------------

   subroutine ocn_thick_ale_init(err)!{{{

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: Error flag

      !-----------------------------------------------------------------

      err = 0 ! default error flag

      ! Determine ALE thickness proportionality type from user config
      select case (trim(config_ALE_thickness_proportionality))
      case ('restingThickness_times_weights')
          ALEthickProportionality = ALEthickProportionThickTimesWgts
      case ('weights_only')
          ALEthickProportionality = ALEthickProportionWgtsOnly
      case default
          call mpas_log_write( &
             ' Error: Invalid config_ALE_thickness_proportionality', &
             MPAS_LOG_CRIT)
      end select

      ! Determine whether to adjust thickness based on user input
      if (config_use_min_max_thickness) then
         useMinMaxThick = .true.
         minThick     = config_min_thickness
         maxThickFact = config_max_thickness_factor
      else
         useMinMaxThick = .false.
         minThick     = 0.0_RKIND
         maxThickFact = 0.0_RKIND
      endif

   !--------------------------------------------------------------------

   end subroutine ocn_thick_ale_init!}}}

!***********************************************************************

end module ocn_thick_ale

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
