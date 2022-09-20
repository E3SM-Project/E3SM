! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_sediment_flux_index
!
!> \brief MPAS ocean analysis mode member: sediment_flux_index
!> \author Zhendong Cao and Phillip J. Wolfram
!> \date   2019/02/26
!> \details
!>  MPAS ocean analysis mode member: sediment_flux_index
!>
!-----------------------------------------------------------------------

module ocn_sediment_flux_index

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_timekeeping
   use mpas_stream_manager

   use ocn_constants
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

   public :: ocn_init_sediment_flux_index, &
             ocn_compute_sediment_flux_index, &
             ocn_restart_sediment_flux_index, &
             ocn_finalize_sediment_flux_index

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_sediment_flux_index
!
!> \brief   Initialize MPAS-Ocean analysis member
!> \author  Zhendong Cao and Phillip J. Wolfram
!> \date    2019/02/26
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_sediment_flux_index(domain, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

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

      err = 0

   end subroutine ocn_init_sediment_flux_index!}}}

!***********************************************************************
!
!  routine ocn_compute_sediment_flux_index
!
!> \brief   Compute MPAS-Ocean analysis member
!> \author  Zhendong Cao and Phillip J. Wolfram
!> \date    2019/02/26
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_sediment_flux_index(domain, timeLevel, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      integer, intent(in) :: timeLevel

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

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

      type (mpas_pool_type), pointer :: sedimentFluxIndexAMPool
      type (dm_info) :: dminfo
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: meshPool

      real (kind=RKIND), dimension(:,:), pointer :: velX, velY, velZ
      real (kind=RKIND), dimension(:), pointer :: posX, posY
      real (kind=RKIND), dimension(:), pointer :: sfiVAx, sfiVAy, sfiBx, sfiBy
      logical, pointer :: on_a_sphere, use_lat_lon_coords
      integer, pointer :: nCells, nVertLevels, nCellsSolve
      integer k, iCell, i
      err = 0

      dminfo = domain % dminfo

      block => domain % blocklist
      do while (associated(block))
         ! get dimensions
         call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
         call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)

         ! get pointers to pools
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'sedimentFluxIndexAM', sedimentFluxIndexAMPool)

         call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)
         call mpas_pool_get_config(ocnConfigs, 'config_AM_sedimentFluxIndex_use_lat_lon_coords', use_lat_lon_coords)

         if (.not. on_a_sphere) then
            use_lat_lon_coords = .false.
         end if

         if (use_lat_lon_coords) then
            call mpas_pool_get_array(meshPool, 'lonCell', posX)
            call mpas_pool_get_array(meshPool, 'latCell', posY)
            velX => velocityZonal
            velY => velocityMeridional
         else
            call mpas_pool_get_array(meshPool, 'xCell', posX)
            call mpas_pool_get_array(meshPool, 'yCell', posY)
            velX => velocityX
            velY => velocityY
         end if 

         call mpas_pool_get_array(sedimentFluxIndexAMPool, 'sedimentFluxIndexVAX', sfiVAx)
         call mpas_pool_get_array(sedimentFluxIndexAMPool, 'sedimentFluxIndexVAY', sfiVAy)
         call mpas_pool_get_array(sedimentFluxIndexAMPool, 'sedimentFluxIndexBX', sfiBx)
         call mpas_pool_get_array(sedimentFluxIndexAMPool, 'sedimentFluxIndexBY', sfiBy)
         ! Computations which are functions of nCells, nEdges, or nVertices
         ! must be placed within this block loop
         ! Here are some example loops
         do iCell = 1,nCellsSolve
            sfiVAx(iCell) = (sum(velX(:,iCell))/float(nVertLevels))**3.0_RKIND
            sfiVAy(iCell) = (sum(velY(:,iCell))/float(nVertLevels))**3.0_RKIND
            sfiBx(iCell) = velX(1,iCell)**3.0_RKIND
            sfiBy(iCell) = velY(1,iCell)**3.0_RKIND
         end do

         block => block % next
      end do

   end subroutine ocn_compute_sediment_flux_index!}}}

!***********************************************************************
!
!  routine ocn_restart_sediment_flux_index
!
!> \brief   Save restart for MPAS-Ocean analysis member
!> \author  Zhendong Cao and Phillip J. Wolfram
!> \date    2019/02/26
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_restart_sediment_flux_index(domain, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

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

      err = 0

   end subroutine ocn_restart_sediment_flux_index!}}}

!***********************************************************************
!
!  routine ocn_finalize_sediment_flux_index
!
!> \brief   Finalize MPAS-Ocean analysis member
!> \author  Zhendong Cao and Phillip J. Wolfram
!> \date    2019/02/26
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_finalize_sediment_flux_index(domain, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

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

      err = 0

   end subroutine ocn_finalize_sediment_flux_index!}}}

end module ocn_sediment_flux_index

! vim: foldmethod=marker
