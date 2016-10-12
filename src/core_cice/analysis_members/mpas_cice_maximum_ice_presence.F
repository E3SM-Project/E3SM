! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_maximum_ice_presence
!
!> \brief MPAS sea ice analysis mode member: maximum_ice_presence
!> \author Adrian K. Turner
!> \date   2nd December 2015
!> \details
!>  MPAS sea ice analysis mode member: maximum_ice_presence
!>
!-----------------------------------------------------------------------

module cice_maximum_ice_presence

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

   public :: cice_bootstrap_maximum_ice_presence, &
             cice_init_maximum_ice_presence, &
             cice_precompute_maximum_ice_presence, &
             cice_compute_maximum_ice_presence, &
             cice_restart_maximum_ice_presence, &
             cice_finalize_maximum_ice_presence

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine cice_bootstrap_maximum_ice_presence
!
!> \brief   Bootstrap MPAS-CICE analysis member
!> \author  Adrian K. Turner
!> \date    2nd December 2015
!> \details
!>  This routine conducts all bootstraps required for the
!>  MPAS-CICE analysis member.
!
!-----------------------------------------------------------------------

   subroutine cice_bootstrap_maximum_ice_presence(domain, instance, err)!{{{

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

   end subroutine cice_bootstrap_maximum_ice_presence!}}}

!***********************************************************************
!
!  routine cice_init_maximum_ice_presence
!
!> \brief   Initialize MPAS-CICE analysis member
!> \author  Adrian K. Turner
!> \date    2nd December 2015
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-CICE analysis member.
!
!-----------------------------------------------------------------------

   subroutine cice_init_maximum_ice_presence(domain, instance, err)!{{{

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

   end subroutine cice_init_maximum_ice_presence!}}}

!***********************************************************************
!
!  routine cice_precompute_maximum_ice_presence
!
!> \brief   Precompute MPAS-CICE analysis member
!> \author  Adrian K. Turner
!> \date    2nd December 2015
!> \details
!>  This routine conducts all pre-computation required for this
!>  MPAS-CICE analysis member.
!
!-----------------------------------------------------------------------

   subroutine cice_precompute_maximum_ice_presence(domain, instance, timeLevel, err)!{{{

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

   end subroutine cice_precompute_maximum_ice_presence!}}}

!***********************************************************************
!
!  routine cice_compute_maximum_ice_presence
!
!> \brief   Compute MPAS-CICE analysis member
!> \author  Adrian K. Turner
!> \date    2nd December 2015
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-CICE analysis member.
!
!-----------------------------------------------------------------------

   subroutine cice_compute_maximum_ice_presence(domain, instance, timeLevel, err)!{{{

     use cice_constants, only: &
          cicePuny

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

      type(block_type), pointer :: &
           block

      type(MPAS_pool_type), pointer :: &
           maximumIcePresenceAMPool, &
           tracersAggregatePool

      real(kind=RKIND), dimension(:), pointer :: &
           iceAreaCell

      integer, pointer :: &
           nCellsSolve, &
           nMonths

      integer :: &
           iCell, &
           iMonth, &
           ierrTime

      type(MPAS_Time_Type) :: &
           currentTime, &
           startTime

      integer, dimension(:,:), pointer :: &
           maximumIcePresence

      character(len=strKIND), pointer :: &
           config_AM_maximumIcePresence_start_time

      err = 0

      ! ice domain
      currentTime = mpas_get_clock_time(domain % clock, MPAS_NOW, ierrTime)
      call MPAS_get_time(currentTime, MM=iMonth)
      call MPAS_pool_get_config(domain % configs, "config_AM_maximumIcePresence_start_time", config_AM_maximumIcePresence_start_time)
      call MPAS_set_time(startTime, dateTimeString=config_AM_maximumIcePresence_start_time)

      if (currentTime > startTime) then

         block => domain % blocklist
         do while (associated(block))

            call MPAS_pool_get_subpool(block % structs, "maximumIcePresenceAM", maximumIcePresenceAMPool)
            call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracersAggregatePool)

            call MPAS_pool_get_array(maximumIcePresenceAMPool, "maximumIcePresence", maximumIcePresence)
            call MPAS_pool_get_array(tracersAggregatePool, "iceAreaCell", iceAreaCell)

            call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)
            call MPAS_pool_get_dimension(block % dimensions, "nMonths", nMonths)

            do iCell = 1, nCellsSolve

               if (iceAreaCell(iCell) > cicePuny) then

                  maximumIcePresence(iMonth, iCell) = 1
                  maximumIcePresence(nMonths,iCell) = 1

               endif

            enddo ! iCell

            block => block % next
         enddo

      endif

   end subroutine cice_compute_maximum_ice_presence!}}}

!***********************************************************************
!
!  routine cice_restart_maximum_ice_presence
!
!> \brief   Save restart for MPAS-CICE analysis member
!> \author  Adrian K. Turner
!> \date    2nd December 2015
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-CICE analysis member.
!
!-----------------------------------------------------------------------

   subroutine cice_restart_maximum_ice_presence(domain, instance, err)!{{{

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

   end subroutine cice_restart_maximum_ice_presence!}}}

!***********************************************************************
!
!  routine cice_finalize_maximum_ice_presence
!
!> \brief   Finalize MPAS-CICE analysis member
!> \author  Adrian K. Turner
!> \date    2nd December 2015
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-CICE analysis member.
!
!-----------------------------------------------------------------------

   subroutine cice_finalize_maximum_ice_presence(domain, instance, err)!{{{

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

   end subroutine cice_finalize_maximum_ice_presence!}}}

!-----------------------------------------------------------------------

end module cice_maximum_ice_presence

! vim: foldmethod=marker
