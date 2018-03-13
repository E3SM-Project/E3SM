! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_high_frequency_output
!
!> \brief MPAS-Seaice analysis mode member: high_frequency_output
!> \author Todd Ringler
!> \date   2015/06/12
!> \details
!>  MPAS-Seaice analysis mode member: high_frequency_output
!>
!-----------------------------------------------------------------------

module seaice_high_frequency_output

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_timekeeping
   use mpas_stream_manager
   use ice_constants_colpkg

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

   public :: seaice_bootstrap_high_frequency_output, &
             seaice_init_high_frequency_output, &
             seaice_precompute_high_frequency_output, &
             seaice_compute_high_frequency_output, &
             seaice_restart_high_frequency_output, &
             seaice_finalize_high_frequency_output

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine seaice_bootstrap_high_frequency_output
!
!> \brief   Bootstrap MPAS-Seaice analysis member
!> \author  Adrian Turner
!> \date    10th November 2015
!> \details
!>  This routine conducts all bootstraps required for the
!>  MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_bootstrap_high_frequency_output(domain, instance, err)!{{{

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

   end subroutine seaice_bootstrap_high_frequency_output!}}}

!***********************************************************************
!
!  routine seaice_init_high_frequency_output
!
!> \brief   Initialize MPAS-Seaice analysis member
!> \author  Todd Ringler
!> \date    2015/06/12
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_init_high_frequency_output(domain, instance, err)!{{{

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

   end subroutine seaice_init_high_frequency_output!}}}

!***********************************************************************
!
!  routine seaice_precompute_high_frequency_output
!
!> \brief   Precompute MPAS-Seaice analysis member
!> \author  Adrian K. Turner
!> \date    2015/09/09
!> \details
!>  This routine conducts all precomputation required for this
!>  MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_precompute_high_frequency_output(domain, instance, timeLevel, err)!{{{

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

   end subroutine seaice_precompute_high_frequency_output!}}}

!***********************************************************************
!
!  routine seaice_compute_high_frequency_output
!
!> \brief   Compute MPAS-Seaice analysis member
!> \author  Todd Ringler
!> \date    2015/06/12
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_compute_high_frequency_output(domain, instance, timeLevel, err)!{{{

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

      type (dm_info) :: dminfo
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: tracersPool
      type (mpas_pool_type), pointer :: highFrequencyOutputAMPool

      real (kind=RKIND), dimension(:), pointer :: iceAreaCellCategory1
      real (kind=RKIND), dimension(:,:,:), pointer :: iceAreaCategory

      err = 0

      dminfo = domain % dminfo

      block => domain % blocklist
      do while (associated(block))
         ! get dimensions

         ! get pointers to pools
         call mpas_pool_get_subpool(block % structs, 'tracers', tracersPool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(domain % blocklist % structs, 'highFrequencyOutputAM', highFrequencyOutputAMPool)

         ! get arrays that will be 'sliced' and put into high frequency output
         call mpas_pool_get_array(tracersPool, 'iceAreaCategory', iceAreaCategory, timeLevel)

         ! get arrays that can be written to output at high freqency
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'iceAreaCellCategory1', iceAreaCellCategory1)

         ! copy data into high frequency output fields
         iceAreaCellCategory1(:) = iceAreaCategory(1,1,:)

         block => block % next
      end do

   end subroutine seaice_compute_high_frequency_output!}}}

!***********************************************************************
!
!  routine seaice_restart_high_frequency_output
!
!> \brief   Save restart for MPAS-Seaice analysis member
!> \author  Todd Ringler
!> \date    2015/06/12
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_restart_high_frequency_output(domain, instance, err)!{{{

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

   end subroutine seaice_restart_high_frequency_output!}}}

!***********************************************************************
!
!  routine seaice_finalize_high_frequency_output
!
!> \brief   Finalize MPAS-Seaice analysis member
!> \author  Todd Ringler
!> \date    2015/06/12
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_finalize_high_frequency_output(domain, instance, err)!{{{

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

   end subroutine seaice_finalize_high_frequency_output!}}}

end module seaice_high_frequency_output

! vim: foldmethod=marker
