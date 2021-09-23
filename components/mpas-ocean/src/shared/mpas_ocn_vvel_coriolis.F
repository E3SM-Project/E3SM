! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vel_hadv_coriolis
!
!> \brief MPAS ocean tendency of coriolis for vertical momentum equation
!> \author Luke Van Roekel
!> \date   December 2020
!> \details
!>  This module contains the routine for computing
!>  tendencies from the vertical momentum coriolis force.
!>
!
!-----------------------------------------------------------------------

module ocn_vvel_coriolis

   use mpas_timer
   use mpas_derived_types
   use mpas_pool_routines
   use ocn_constants
   use ocn_config
   use mpas_constants

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

   public :: ocn_vvel_coriolis_tend, &
             ocn_vvel_coriolis_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: vertCoriolisOn

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_vvel_coriolis_tend
!
!> \brief   Computes tendency term for coriolis force for vertical momentum equation
!> \author  Luke Van Roekel
!> \date    January 2021
!> \details

!>  This routine computes the coriolis tendencies for vertical momentum based on current state.

!
!-----------------------------------------------------------------------

   subroutine ocn_vvel_coriolis_tend(meshPool, velocityZonal, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         velocityZonal    !< Input: Horizontal velocity

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         tend          !< Input/Output: velocity tendency

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer, dimension(:), pointer :: maxLevelCell, nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnEdge, edgesOnCell, edgeMask

      integer :: j, k, i, iCell
      integer :: cell1, cell2, iEdge, nCells
      integer, pointer :: nVertLevels
      integer, dimension(:), pointer :: nCellsArray
      real (kind=RKIND) :: workVorticity, invLength, edgeWeight, r_tmp
      real, dimension(:), pointer :: latCell
      real (kind=RKIND), dimension(:), allocatable :: qArr

      err = 0

      if ( .not. vertCoriolisOn ) return

      call mpas_timer_start("vertical velocity coriolis")

      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

      call mpas_pool_get_array(meshPool, 'latCell', latCell)

      call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      nCells = nCellsArray( 2 )

      !$omp parallel
      !$omp do schedule(runtime) &
      !$omp private(cell1, cell2, invLength, k, qArr, j, eoe, edgeWeight, workVorticity)
      do iCell= 1, nCells

         do k = 2, maxLevelCell(iCell)
            tend(k, iCell) = tend(k, iCell) + 2.0_RKIND*omega*cos(latCell(iCell))* &
                    (0.5_RKIND*(velocityZonal(k,iCell) * velocityZonal(k-1,iCell)))
                    !FIXME: zonalvelocity may be a bad choice here as it's not updated during a timestep
                    ! but maybe not the worst assumption.
         end do

      end do
      !$omp end do
      !$omp end parallel

      call mpas_timer_stop("vertical velocity coriolis")

   !--------------------------------------------------------------------

 end subroutine ocn_vvel_coriolis_tend!}}}

!***********************************************************************
!
!  routine ocn_vvel_pgf_init
!
!> \brief   Initializes ocean vertical momentum PGF tendencies
!> \author  Luke Van Roekel
!> \date    January 2021
!
!-----------------------------------------------------------------------

   subroutine ocn_vvel_coriolis_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! Output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      err = 0

      vertCoriolisOn = .true.
      if ( config_disable_vvel_coriolis .or.  &
          .not. config_enable_nonhydrostatic_mode ) vertCoriolisOn = .false.

   !--------------------------------------------------------------------

 end subroutine ocn_vvel_coriolis_init!}}}

!***********************************************************************

end module ocn_vvel_coriolis

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
