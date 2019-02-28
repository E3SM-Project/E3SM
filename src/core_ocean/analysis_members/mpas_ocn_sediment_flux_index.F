! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
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
!>  In order to add a new analysis member, do the following:
!>  1. Copy these to your new analysis member name:
!>     cp mpas_ocn_TEMPLATE.F mpas_ocn_your_new_name.F
!>     cp Registry_TEMPLATE.xml Registry_your_new_name.xml
!>
!>  2. In those two new files, replace the following text:
!>     temPlate, TEM_PLATE, FILL_IN_AUTHOR, FILL_IN_DATE
!>     Typically temPlate uses camel case (variable names), like yourNewName,
!>     while TEM_PLATE uses underscores (subroutine names), like your_new_name.
!>     note: do not replace 'filename_template' in Registry_yourNewName.xml
!>
!>  3. Add a #include line for your registry to
!>     Registry_analysis_members.xml
!>
!>  4. In mpas_ocn_analysis_driver.F, add a use statement for your new analysis member.
!>     In addition, add lines for your analysis member, and replace TEM_PLATE
!>     and temPlate as described in step 2. There should be 5 places that need additions:
!>      - Adding the analysis member name to the analysis member list
!>      - Adding an init if test can subroutine call
!>      - Adding a compute if test can subroutine call
!>      - Adding a restart if test can subroutine call
!>      - Adding a finalize if test can subroutine call
!>
!>  5. In src/core_ocean/analysis_members/Makefile, add your
!>     new analysis member to the list of members. See another analysis member
!>     in that file for an example.
!>     NOTE: If your analysis member depends on other files, add a dependency
!>           line for the member and list them there. See okubo weiss for an example.
!>
!-----------------------------------------------------------------------

module ocn_sediment_flux_index

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_timekeeping
   use mpas_stream_manager

   use ocn_constants
   use ocn_diagnostics_routines

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
      type (mpas_pool_type), pointer :: diagnosticsPool

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
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(block % structs, 'sedimentFluxIndexAM', sedimentFluxIndexAMPool)

         call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)
         call mpas_pool_get_config(ocnConfigs, 'config_AM_sedimentFluxIndex_use_lat_lon_coords', use_lat_lon_coords)

         if (.not. on_a_sphere) then
            use_lat_lon_coords = .false.
         end if

         if (use_lat_lon_coords) then
            call mpas_pool_get_array(meshPool, 'lonCell', posX)
            call mpas_pool_get_array(meshPool, 'latCell', posY)
            call mpas_pool_get_array(diagnosticsPool, 'velocityZonal', velX)
            call mpas_pool_get_array(diagnosticsPool, 'velocityMeridional', velY)
         else
            call mpas_pool_get_array(meshPool, 'xCell', posX)
            call mpas_pool_get_array(meshPool, 'yCell', posY)
            call mpas_pool_get_array(diagnosticsPool, 'velocityX', velX)
            call mpas_pool_get_array(diagnosticsPool, 'velocityY', velY)
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
