! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_ice_present
!
!> \brief MPAS sea ice analysis mode member: ice_present
!> \author Adrian K. Turner
!> \date   14th November 2015
!> \details
!>  MPAS sea ice analysis mode member: ice_present
!>
!-----------------------------------------------------------------------

module seaice_ice_present

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_timekeeping
   use mpas_stream_manager

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

   public :: seaice_bootstrap_ice_present, &
             seaice_init_ice_present, &
             seaice_precompute_ice_present, &
             seaice_compute_ice_present, &
             seaice_restart_ice_present, &
             seaice_finalize_ice_present

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine seaice_bootstrap_ice_present
!
!> \brief   Bootstrap MPAS-Seaice analysis member
!> \author  Adrian K. Turner
!> \date    14th November 2015
!> \details
!>  This routine conducts all bootstraps required for the
!>  MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_bootstrap_ice_present(domain, instance, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

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

   end subroutine seaice_bootstrap_ice_present!}}}

!***********************************************************************
!
!  routine seaice_init_ice_present
!
!> \brief   Initialize MPAS-Seaice analysis member
!> \author  Adrian K. Turner
!> \date    14th November 2015
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_init_ice_present(domain, instance, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

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

   end subroutine seaice_init_ice_present!}}}

!***********************************************************************
!
!  routine seaice_precompute_ice_present
!
!> \brief   Precompute MPAS-Seaice analysis member
!> \author  Adrian K. Turner
!> \date    14th November 2015
!> \details
!>  This routine conducts all pre-computation required for this
!>  MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_precompute_ice_present(domain, instance, timeLevel, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

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

      err = 0

   end subroutine seaice_precompute_ice_present!}}}

!***********************************************************************
!
!  routine seaice_compute_ice_present
!
!> \brief   Compute MPAS-Seaice analysis member
!> \author  Adrian K. Turner
!> \date    14th November 2015
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_compute_ice_present(domain, instance, timeLevel, err)!{{{

      use seaice_constants, only: &
           seaicePuny

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

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

      type (block_type), pointer :: &
           block

      type (mpas_pool_type), pointer :: &
           icePresentAMPool, &
           tracersAggregatePool

      real(kind=RKIND), dimension(:), pointer :: &
           icePresent, &
           iceAreaCell

      integer, pointer :: &
           nCellsSolve

      integer :: &
           iCell

      err = 0

      block => domain % blocklist
      do while (associated(block))

         call mpas_pool_get_subpool(block % structs, 'icePresentAM', icePresentAMPool)
         call mpas_pool_get_subpool(block % structs, 'tracers_aggregate', tracersAggregatePool)

         call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)

         call mpas_pool_get_array(icePresentAMPool, 'icePresent', icePresent)
         call mpas_pool_get_array(tracersAggregatePool, 'iceAreaCell', iceAreaCell)

         do iCell = 1, nCellsSolve

            if (iceAreaCell(iCell) > seaicePuny) then

               icePresent(iCell) = 1.0_RKIND

            else

               icePresent(iCell) = 0.0_RKIND

            endif

         enddo ! iCell

         block => block % next
      end do

   end subroutine seaice_compute_ice_present!}}}

!***********************************************************************
!
!  routine seaice_restart_ice_present
!
!> \brief   Save restart for MPAS-Seaice analysis member
!> \author  Adrian K. Turner
!> \date    14th November 2015
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_restart_ice_present(domain, instance, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

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

   end subroutine seaice_restart_ice_present!}}}

!***********************************************************************
!
!  routine seaice_finalize_ice_present
!
!> \brief   Finalize MPAS-Seaice analysis member
!> \author  Adrian K. Turner
!> \date    14th November 2015
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_finalize_ice_present(domain, instance, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

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

   end subroutine seaice_finalize_ice_present!}}}

!-----------------------------------------------------------------------

end module seaice_ice_present

! vim: foldmethod=marker
