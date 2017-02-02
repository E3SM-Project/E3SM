! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_eddy_product_variables
!
!> \brief MPAS ocean analysis mode member: eddy_product_variables
!> \author Mark Petersen
!> \date   Oct 2016
!> \details
!>  MPAS ocean analysis mode member: eddy_product_variables
!> computes:
!> SSH^2 SSHSquared
!> u^2 velocityZonalSquared
!> v^2  velocityMeridionalSquared
!> uT velocityZonalTimesTemperature
!> vT  velocityMeridionalTimesTemperature
!>
!> In future, could additionally compute:
!> uS velocityZonalTimesSalinity
!> vS  velocityMeridionalTimesSalinity
!> uw velocityZonalTimesVertVel
!> vw  velocityMeridionalTimesVerVel
!>
!> They you just stick these in the time-average AM.  We already time average u,v, and all tracers.
!> Then compute <v'T'> etc as a post processing step, same as POP, as follows:
!>
!> <v'T'> = <vT> - <v><T>
!> <(u')^2> = <u^2> - <u>^2
!>
!-----------------------------------------------------------------------

module ocn_eddy_product_variables

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

   public :: ocn_init_eddy_product_variables, &
             ocn_compute_eddy_product_variables, &
             ocn_restart_eddy_product_variables, &
             ocn_finalize_eddy_product_variables

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_eddy_product_variables
!
!> \brief   Initialize MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    Oct 2016
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_eddy_product_variables(domain, err)!{{{

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

   end subroutine ocn_init_eddy_product_variables!}}}

!***********************************************************************
!
!  routine ocn_compute_eddy_product_variables
!
!> \brief   Compute MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    Oct 2016
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_eddy_product_variables(domain, timeLevel, err)!{{{

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

      type (mpas_pool_type), pointer :: eddyProductVariablesAMPool
      type (dm_info) :: dminfo
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: tracersPool
      type (mpas_pool_type), pointer :: diagnosticsPool

      integer, pointer :: nVertLevels, nCellsSolve
      integer :: iTracer, k, iCell
      integer, pointer :: index_temperature, index_salinity
      integer, dimension(:), pointer :: maxLevelCell

      real (kind=RKIND), dimension(:), pointer :: ssh, SSHSquared
      real (kind=RKIND), dimension(:,:), pointer :: velocityZonal, velocityMeridional, &
           velocityZonalSquared, velocityMeridionalSquared, velocityZonalTimesTemperature, velocityMeridionalTimesTemperature
      real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers

      err = 0

      dminfo = domain % dminfo

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(block % structs, 'eddyProductVariablesAM', eddyProductVariablesAMPool)

         call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)
         call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)

         call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 1)
         call mpas_pool_get_dimension(tracersPool, 'index_temperature', index_temperature)

         call mpas_pool_get_array(statePool, 'ssh',ssh)
         call mpas_pool_get_array(diagnosticsPool, 'velocityZonal', velocityZonal)
         call mpas_pool_get_array(diagnosticsPool, 'velocityMeridional', velocityMeridional)

         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
         call mpas_pool_get_array(eddyProductVariablesAMPool, 'SSHSquared', SSHSquared)
         call mpas_pool_get_array(eddyProductVariablesAMPool, 'velocityZonalSquared', velocityZonalSquared)
         call mpas_pool_get_array(eddyProductVariablesAMPool, 'velocityMeridionalSquared', velocityMeridionalSquared)
         call mpas_pool_get_array(eddyProductVariablesAMPool, 'velocityZonalTimesTemperature', velocityZonalTimesTemperature)
         call mpas_pool_get_array(eddyProductVariablesAMPool, 'velocityMeridionalTimesTemperature', &
            velocityMeridionalTimesTemperature)

         do iCell = 1,nCellsSolve
            SSHSquared(iCell) = ssh(iCell)**2
            do k = 1, maxLevelCell(iCell)
               velocityZonalSquared(k,iCell) = velocityZonal(k,iCell)**2
               velocityMeridionalSquared(k,iCell) = velocityMeridional(k,iCell)**2
               velocityZonalTimesTemperature(k,iCell) = velocityZonal(k,iCell)*activeTracers(index_temperature,k,iCell)
               velocityMeridionalTimesTemperature(k,iCell) = velocityMeridional(k,iCell)*activeTracers(index_temperature,k,iCell)
            end do
         end do

         block => block % next
      end do

   end subroutine ocn_compute_eddy_product_variables!}}}

!***********************************************************************
!
!  routine ocn_restart_eddy_product_variables
!
!> \brief   Save restart for MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    Oct 2016
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_restart_eddy_product_variables(domain, err)!{{{

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

   end subroutine ocn_restart_eddy_product_variables!}}}

!***********************************************************************
!
!  routine ocn_finalize_eddy_product_variables
!
!> \brief   Finalize MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    Oct 2016
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_finalize_eddy_product_variables(domain, err)!{{{

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

   end subroutine ocn_finalize_eddy_product_variables!}}}

end module ocn_eddy_product_variables

! vim: foldmethod=marker
