! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_load_balance
!
!> \brief MPAS sea ice analysis mode member: load_balance
!> \author Adrian K. Turner
!> \date   2nd December 2015
!> \details
!>  MPAS sea ice analysis mode member: load_balance
!>
!-----------------------------------------------------------------------

module seaice_load_balance

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_timekeeping
   use mpas_stream_manager
   use mpas_log, only: mpas_log_write

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

   public :: seaice_bootstrap_load_balance, &
             seaice_init_load_balance, &
             seaice_precompute_load_balance, &
             seaice_compute_load_balance, &
             seaice_restart_load_balance, &
             seaice_finalize_load_balance

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine seaice_bootstrap_load_balance
!
!> \brief   Bootstrap MPAS-Seaice analysis member
!> \author  Adrian K. Turner
!> \date    2nd December 2015
!> \details
!>  This routine conducts all bootstraps required for the
!>  MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_bootstrap_load_balance(domain, instance, err)!{{{

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

   end subroutine seaice_bootstrap_load_balance!}}}

!***********************************************************************
!
!  routine seaice_init_load_balance
!
!> \brief   Initialize MPAS-Seaice analysis member
!> \author  Adrian K. Turner
!> \date    2nd December 2015
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_init_load_balance(domain, instance, err)!{{{

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

      type(block_type), pointer :: &
           block

      type(MPAS_pool_type), pointer :: &
           loadBalanceAMPool

      integer, dimension(:), pointer :: &
           nCellsProc

      integer, pointer :: &
           nCellsSolve, &
           nProcs

      integer, dimension(:), allocatable :: &
           nCellsProcTmp

      err = 0

      call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nProcs", nProcs)
      if (nProcs /= domain % dminfo % nprocs) then
         call mpas_log_write('mpas_seaice_load_balance: wrong nProcs set', MPAS_LOG_CRIT)
      endif

      call MPAS_pool_get_subpool(domain % blocklist % structs, "loadBalanceAM", loadBalanceAMPool)
      call MPAS_pool_get_array(loadBalanceAMPool, "nCellsProc", nCellsProc)

      allocate(nCellsProcTmp(domain % dminfo % nprocs))

      nCellsProcTmp(:) = -1

      nCellsProcTmp(domain % dminfo % my_proc_id+1) = 0

      block => domain % blocklist
      do while (associated(block))

         call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)

         nCellsProcTmp(domain % dminfo % my_proc_id+1) = nCellsProcTmp(domain % dminfo % my_proc_id+1) + nCellsSolve

         block => block % next
      enddo

      call MPAS_dmpar_max_int_array(domain % dminfo, domain % dminfo % nprocs, nCellsProcTmp, nCellsProc)

      deallocate(nCellsProcTmp)

   end subroutine seaice_init_load_balance!}}}

!***********************************************************************
!
!  routine seaice_precompute_load_balance
!
!> \brief   Precompute MPAS-Seaice analysis member
!> \author  Adrian K. Turner
!> \date    2nd December 2015
!> \details
!>  This routine conducts all pre-computation required for this
!>  MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_precompute_load_balance(domain, instance, timeLevel, err)!{{{

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

   end subroutine seaice_precompute_load_balance!}}}

!***********************************************************************
!
!  routine seaice_compute_load_balance
!
!> \brief   Compute MPAS-Seaice analysis member
!> \author  Adrian K. Turner
!> \date    2nd December 2015
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_compute_load_balance(domain, instance, timeLevel, err)!{{{

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

      type(block_type), pointer :: &
           block

      type(MPAS_pool_type), pointer :: &
           loadBalanceAMPool, &
           tracersAggregatePool

      integer, dimension(:), pointer :: &
           nCellsProcWithSeaIce

      real(kind=RKIND), dimension(:), pointer :: &
           iceAreaCell

      integer, pointer :: &
           nCellsSolve

      integer, dimension(:), allocatable :: &
           nCellsProcWithSeaIceTmp

      integer :: &
           iCell

      err = 0

      ! cells per processor
      call MPAS_pool_get_subpool(domain % blocklist % structs, "loadBalanceAM", loadBalanceAMPool)
      call MPAS_pool_get_array(loadBalanceAMPool, "nCellsProcWithSeaIce", nCellsProcWithSeaIce)

      allocate(nCellsProcWithSeaIceTmp(domain % dminfo % nprocs))

      nCellsProcWithSeaIceTmp(:) = -1

      nCellsProcWithSeaIceTmp(domain % dminfo % my_proc_id+1) = 0

      block => domain % blocklist
      do while (associated(block))

         call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracersAggregatePool)
         call MPAS_pool_get_array(tracersAggregatePool, "iceAreaCell", iceAreaCell)

         call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)

         do iCell = 1, nCellsSolve

            if (iceAreaCell(iCell) > seaicePuny) then

               nCellsProcWithSeaIceTmp(domain % dminfo % my_proc_id+1) = &
                    nCellsProcWithSeaIceTmp(domain % dminfo % my_proc_id+1) + 1

            endif

         enddo ! iCell

         block => block % next
      enddo

      call MPAS_dmpar_max_int_array(domain % dminfo, domain % dminfo % nprocs, nCellsProcWithSeaIceTmp, nCellsProcWithSeaIce)

      deallocate(nCellsProcWithSeaIceTmp)

   end subroutine seaice_compute_load_balance!}}}

!***********************************************************************
!
!  routine seaice_restart_load_balance
!
!> \brief   Save restart for MPAS-Seaice analysis member
!> \author  Adrian K. Turner
!> \date    2nd December 2015
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_restart_load_balance(domain, instance, err)!{{{

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

   end subroutine seaice_restart_load_balance!}}}

!***********************************************************************
!
!  routine seaice_finalize_load_balance
!
!> \brief   Finalize MPAS-Seaice analysis member
!> \author  Adrian K. Turner
!> \date    2nd December 2015
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_finalize_load_balance(domain, instance, err)!{{{

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

   end subroutine seaice_finalize_load_balance!}}}

!-----------------------------------------------------------------------

end module seaice_load_balance

! vim: foldmethod=marker
