!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_ridging_diagnostics
!
!> \brief MPAS sea ice analysis mode member: ridging_diagnostics
!> \author Adrian K. Turner
!> \date   7th September 2015
!> \details
!>  MPAS sea ice analysis mode member: ridging_diagnostics
!>
!-----------------------------------------------------------------------

module seaice_ridging_diagnostics

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

   public :: seaice_bootstrap_ridging_diagnostics, &
             seaice_init_ridging_diagnostics, &
             seaice_precompute_ridging_diagnostics, &
             seaice_compute_ridging_diagnostics, &
             seaice_restart_ridging_diagnostics, &
             seaice_finalize_ridging_diagnostics

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine seaice_bootstrap_ridging_diagnostics
!
!> \brief   Bootstrap MPAS-Seaice analysis member
!> \author  Adrian K. Turner
!> \date    10th November 2015
!> \details
!>  This routine conducts all bootstraps required for the
!>  MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_bootstrap_ridging_diagnostics(domain, instance, err)!{{{

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

   end subroutine seaice_bootstrap_ridging_diagnostics!}}}

!***********************************************************************
!
!  routine seaice_init_ridging_diagnostics
!
!> \brief   Initialize MPAS-Seaice analysis member
!> \author  Adrian K. Turner
!> \date    7th September 2015
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_init_ridging_diagnostics(domain, instance, err)!{{{

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

   end subroutine seaice_init_ridging_diagnostics!}}}

!***********************************************************************
!
!  routine seaice_precompute_ridging_diagnostics
!
!> \brief   Precompute MPAS-Seaice analysis member
!> \author  Adrian K. Turner
!> \date    9th September 2015
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_precompute_ridging_diagnostics(domain, instance, timeLevel, err)!{{{

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

   end subroutine seaice_precompute_ridging_diagnostics!}}}

!***********************************************************************
!
!  routine seaice_compute_ridging_diagnostics
!
!> \brief   Compute MPAS-Seaice analysis member
!> \author  Adrian K. Turner
!> \date    7th September 2015
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_compute_ridging_diagnostics(domain, instance, timeLevel, err)!{{{

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

      type (dm_info) :: &
           dminfo

      type (block_type), pointer :: &
           block

      type (MPAS_pool_type), pointer :: &
           ridgingDiagnosticsAMPool, &
           tracersPool

      real(kind=RKIND), dimension(:), pointer :: &
           levelIceAreaAverage, &
           ridgedIceAreaAverage, &
           levelIceVolumeAverage, &
           ridgedIceVolumeAverage

      real(kind=RKIND), dimension(:,:,:), pointer :: &
           iceAreaCategory, &
           iceVolumeCategory, &
           levelIceArea, &
           levelIceVolume

      integer, pointer :: &
           nCellsSolve, &
           nCategories

      integer :: &
           iCell, &
           iCategory

      err = 0

      dminfo = domain % dminfo

      block => domain % blocklist
      do while (associated(block))

         call MPAS_pool_get_subpool(block % structs, 'ridgingDiagnosticsAM', ridgingDiagnosticsAMPool)
         call MPAS_pool_get_subpool(block % structs, 'tracers', tracersPool)

         call MPAS_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
         call MPAS_pool_get_dimension(block % dimensions, 'nCategories', nCategories)

         call MPAS_pool_get_array(ridgingDiagnosticsAMPool, "levelIceAreaAverage", levelIceAreaAverage)
         call MPAS_pool_get_array(ridgingDiagnosticsAMPool, "ridgedIceAreaAverage", ridgedIceAreaAverage)
         call MPAS_pool_get_array(ridgingDiagnosticsAMPool, "levelIceVolumeAverage", levelIceVolumeAverage)
         call MPAS_pool_get_array(ridgingDiagnosticsAMPool, "ridgedIceVolumeAverage", ridgedIceVolumeAverage)

         call MPAS_pool_get_array(tracersPool, "iceAreaCategory", iceAreaCategory, 1)
         call MPAS_pool_get_array(tracersPool, "iceVolumeCategory", iceVolumeCategory, 1)
         call MPAS_pool_get_array(tracersPool, "levelIceArea", levelIceArea, 1)
         call MPAS_pool_get_array(tracersPool, "levelIceVolume", levelIceVolume, 1)

         do iCell = 1, nCellsSolve

            levelIceAreaAverage(iCell)    = 0.0_RKIND
            ridgedIceAreaAverage(iCell)   = 0.0_RKIND

            levelIceVolumeAverage(iCell)  = 0.0_RKIND
            ridgedIceVolumeAverage(iCell) = 0.0_RKIND

            do iCategory = 1, nCategories

               levelIceAreaAverage(iCell) = &
                    levelIceAreaAverage(iCell) + iceAreaCategory(1,iCategory,iCell) * &
                    levelIceArea(1,iCategory,iCell)
               ridgedIceAreaAverage(iCell) = &
                    ridgedIceAreaAverage(iCell) + iceAreaCategory(1,iCategory,iCell) * &
                    (1.0_RKIND - levelIceArea(1,iCategory,iCell))

               levelIceVolumeAverage(iCell) = &
                    levelIceVolumeAverage(iCell) + iceVolumeCategory(1,iCategory,iCell) * &
                    levelIceVolume(1,iCategory,iCell)
               ridgedIceVolumeAverage(iCell) = &
                    ridgedIceVolumeAverage(iCell) + iceVolumeCategory(1,iCategory,iCell) * &
                    (1.0_RKIND - levelIceVolume(1,iCategory,iCell))

            enddo ! iCategory
         enddo ! iCell

         block => block % next
      enddo

   end subroutine seaice_compute_ridging_diagnostics!}}}

!***********************************************************************
!
!  routine seaice_restart_ridging_diagnostics
!
!> \brief   Save restart for MPAS-Seaice analysis member
!> \author  Adrian K. Turner
!> \date    7th September 2015
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_restart_ridging_diagnostics(domain, instance, err)!{{{

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

   end subroutine seaice_restart_ridging_diagnostics!}}}

!***********************************************************************
!
!  routine seaice_finalize_ridging_diagnostics
!
!> \brief   Finalize MPAS-Seaice analysis member
!> \author  Adrian K. Turner
!> \date    7th September 2015
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_finalize_ridging_diagnostics(domain, instance, err)!{{{

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

   end subroutine seaice_finalize_ridging_diagnostics!}}}

!-----------------------------------------------------------------------

end module seaice_ridging_diagnostics

! vim: foldmethod=marker
