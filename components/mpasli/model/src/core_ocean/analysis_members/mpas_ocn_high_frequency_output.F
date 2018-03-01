! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_high_frequency_output
!
!> \brief MPAS ocean analysis mode member: high_frequency_output
!> \author Todd Ringler
!> \date   2015/06/12
!> \details
!>  MPAS ocean analysis mode member: high_frequency_output
!>
!-----------------------------------------------------------------------

module ocn_high_frequency_output

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

   public :: ocn_init_high_frequency_output, &
             ocn_compute_high_frequency_output, &
             ocn_restart_high_frequency_output, &
             ocn_finalize_high_frequency_output

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_high_frequency_output
!
!> \brief   Initialize MPAS-Ocean analysis member
!> \author  Todd Ringler
!> \date    2015/06/12
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_high_frequency_output(domain, err)!{{{

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

   end subroutine ocn_init_high_frequency_output!}}}

!***********************************************************************
!
!  routine ocn_compute_high_frequency_output
!
!> \brief   Compute MPAS-Ocean analysis member
!> \author  Todd Ringler
!> \date    2015/06/12
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_high_frequency_output(domain, timeLevel, err)!{{{

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

      type (dm_info) :: dminfo
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: diagnosticsPool
      type (mpas_pool_type), pointer :: forcingPool
      type (mpas_pool_type), pointer :: highFrequencyOutputAMPool

      integer :: iLevel, iLevelTarget
      integer, pointer :: nVertLevels
      real (kind=RKIND), dimension(:), pointer :: refBottomDepth, kineticEnergyAt100m, relativeVorticityAt100m
      real (kind=RKIND), dimension(:,:), pointer :: kineticEnergyCell, relativeVorticityCell, tracersAtSurface
      real (kind=RKIND), dimension(:,:,:), pointer :: tracers

      err = 0

      dminfo = domain % dminfo

      block => domain % blocklist
      do while (associated(block))
         ! get dimensions
         call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)

         ! get pointers to pools
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
         call mpas_pool_get_subpool(domain % blocklist % structs, 'highFrequencyOutputAM', highFrequencyOutputAMPool)

         ! get static data from mesh pool
         call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)

         ! get arrays that will be 'sliced' and put into high frequency output
         call mpas_pool_get_array(diagnosticsPool, 'kineticEnergyCell', kineticEnergyCell, timeLevel)
         call mpas_pool_get_array(diagnosticsPool, 'relativeVorticityCell', relativeVorticityCell, timeLevel)
         call mpas_pool_get_array(statePool, 'tracers', tracers, timeLevel)

         ! get arrays that can be written to output at high freqency
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'kineticEnergyAt100m', kineticEnergyAt100m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'relativeVorticityAt100m', relativeVorticityAt100m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'tracersAtSurface', tracersAtSurface)

         !
         ! note for future build out
         ! HERE: interpolate vertically to target z-level or isopycnal or some other surface
         !

         ! for now, just get close enough
         iLevelTarget = 1
         do iLevel=2,nVertLevels
           if(refBottomDepth(iLevel) > 100.0_RKIND) then
              iLevelTarget = iLevel-1
              exit
           endif
         enddo

         ! copy data into high frequency output fields
         ! tracer data will be converted to new tracer infrastrcture (and this line removed) before June 23 2015.
         kineticEnergyAt100m(:) = kineticEnergyCell(iLevelTarget,:)
         relativeVorticityAt100m(:) = relativeVorticityCell(iLevelTarget,:)
         tracersAtSurface(1,:) = tracers(1,1,:)
         tracersAtSurface(2,:) = tracers(2,1,:)

         block => block % next
      end do

   end subroutine ocn_compute_high_frequency_output!}}}

!***********************************************************************
!
!  routine ocn_restart_high_frequency_output
!
!> \brief   Save restart for MPAS-Ocean analysis member
!> \author  Todd Ringler
!> \date    2015/06/12
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_restart_high_frequency_output(domain, err)!{{{

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

   end subroutine ocn_restart_high_frequency_output!}}}

!***********************************************************************
!
!  routine ocn_finalize_high_frequency_output
!
!> \brief   Finalize MPAS-Ocean analysis member
!> \author  Todd Ringler
!> \date    2015/06/12
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_finalize_high_frequency_output(domain, err)!{{{

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

   end subroutine ocn_finalize_high_frequency_output!}}}

end module ocn_high_frequency_output

! vim: foldmethod=marker
