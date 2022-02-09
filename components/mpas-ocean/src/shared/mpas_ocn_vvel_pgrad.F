! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vvel_pgrad
!
!> \brief MPAS ocean tendency of pressure gradient force for vertical momentum equation
!> \author Luke Van Roekel, Darren Engwirda, Sara Calandrini
!> \date   January 2021, updated November 2021
!> \details
!>  This module contains the routine for computing
!>  tendencies from the vertical momentum coriolis force.
!>
!
!-----------------------------------------------------------------------

module ocn_vvel_pgrad

   use mpas_timer
   use mpas_derived_types
   use mpas_pool_routines
   use ocn_constants
   use ocn_config
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

   public :: ocn_vvel_pgf_tend, &
             ocn_vvel_pgf_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: vertPGFOn

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_vvel_pgf_tend
!
!> \brief   Computes tendency term for pressure gradient in vertical momentum equation
!> \author  Luke Van Roekel, Darren Engwirda, Sara Calandrini
!> \date    January 2021, updated November 2021
!> \details

!>  This routine computes the PGF tendencies for the vertical momentum based on 
!>  current state.
!
!-----------------------------------------------------------------------

   subroutine ocn_vvel_pgf_tend(meshPool, layerThickness, nhPressure, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickness                     !< Input: layer thickness at cells

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         nhPressure                         !< Input: pressure (nonhydrostatic)        

      type (mpas_pool_type), intent(in) :: &
         meshPool           !< Input: mesh information

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

      integer :: j, k, i, iCell, nCells
      integer, pointer :: nVertLevels
      integer, dimension(:), pointer :: maxLevelCell, nCellsArray

      err = 0

      if ( .not. vertPGFOn ) return

      call mpas_timer_start("vertical velocity PGF")

      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      nCells = nCellsArray( 2 )

      !$omp parallel
      !$omp do schedule(runtime) &
      !$omp private(k)
      do iCell= 1, nCells
         tend(1,iCell) = tend(1,iCell) - & ! nhPressure = 0.0 at surface
            2.0_RKIND * (-nhPressure(1,iCell)) / layerThickness(1,iCell)
         do k = 2, maxLevelCell(iCell)
            tend(k,iCell) = tend(k,iCell) - &
               2.0_RKIND * (nhPressure(k-1,iCell) - nhPressure(k,iCell)) &
                  / (layerThickness(k-1,iCell) + layerThickness(k,iCell)) 
         end do

      end do
      !$omp end do
      !$omp end parallel

      call mpas_timer_stop("vertical velocity PGF")

   !--------------------------------------------------------------------

 end subroutine ocn_vvel_pgf_tend!}}}

!***********************************************************************
!
!  routine ocn_vvel_pgf_init
!
!> \brief   Initializes ocean momentum PGF tendencies
!> \author  Luke Van Roekel
!> \date    January 2021
!
!-----------------------------------------------------------------------

   subroutine ocn_vvel_pgf_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! Output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      err = 0

      vertPGFOn = .true.
      if ( .not. config_enable_nonhydrostatic_mode ) vertPGFOn = .false.

   !--------------------------------------------------------------------

 end subroutine ocn_vvel_pgf_init!}}}

!***********************************************************************

end module ocn_vvel_pgrad

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
