! !MODULE: seq_pauseresume_mod --- Module for managing pause/resume data for ESP components
!
! !DESCRIPTION:
!
!     A module to collect and distribute pause/resume information
!
!
! !INTERFACE: ------------------------------------------------------------------

module seq_pauseresume_mod

   ! !USES:

   use shr_kind_mod,       only: CL => SHR_KIND_CL
   use shr_sys_mod,        only: shr_sys_flush, shr_sys_abort
   use seq_comm_mct,       only: num_inst_driver

   implicit none

   private
#include <mpif.h>

   ! !PUBLIC INTERFACES:
   public :: seq_resume_store_comp ! Store component resume filenames
   public :: seq_resume_get_files  ! Retrieve pointer to resume filenames
   public :: seq_resume_free       ! Free resume filename storage
   public :: seq_resume_broadcast  ! Broadcast component filenames to all PEs

   ! Type to hold resume filenames
   type seq_resume_type
      character(len=CL), pointer :: atm_resume(:) => NULL() ! atm resume file(s)
      character(len=CL), pointer :: lnd_resume(:) => NULL() ! lnd resume file(s)
      character(len=CL), pointer :: ice_resume(:) => NULL() ! ice resume file(s)
      character(len=CL), pointer :: ocn_resume(:) => NULL() ! ocn resume file(s)
      character(len=CL), pointer :: glc_resume(:) => NULL() ! glc resume file(s)
      character(len=CL), pointer :: rof_resume(:) => NULL() ! rof resume file(s)
      character(len=CL), pointer :: wav_resume(:) => NULL() ! wav resume file(s)
      character(len=CL), pointer :: cpl_resume(:) => NULL() ! cpl resume file(s)
   end type seq_resume_type

   type(seq_resume_type), pointer :: resume => NULL() ! For storing pause/resume files

   private :: seq_resume_broadcast_array

CONTAINS

   !===========================================================================

   !BOP =======================================================================
   !
   ! !IROUTINE: seq_resume_store_comp - Allocate space & store resume filenames
   !
   ! !DESCRIPTION:
   !
   ! Allocate data for resume filenames from all component instances
   ! Store resume filenames for requested component
   !
   ! Assumptions about instance numbers:
   !   Multi-driver: num_inst_driver = total number of instances of each <comp>
   !                 num_inst_<comp> = 1
   !   Single-driver: num_inst_driver = 1
   !                  num_inst_<comp> = total number of <comp> instances
   !
   ! Assumption about resume names: All <comp> PEs should have the same value
   !                                for <comp>%cdata_cc%resume_filename
   !
   ! !INTERFACE: --------------------------------------------------------------

   subroutine seq_resume_store_comp(oid, filename, num_inst_comp, ninst, iamroot)
      character(len=1),      intent(in)    :: oid           ! 1 letter comp type
      character(len=*),      intent(in)    :: filename      ! resume filename
      integer,               intent(in)    :: num_inst_comp ! # comp instances
      integer,               intent(in)    :: ninst         ! comp  instance #
      logical,               intent(in)    :: iamroot       ! is comp root?

      integer                              :: num_inst      ! # store instances
      character(len=CL), pointer           :: fname_ptr(:)
      character(len=*),  parameter         :: subname = 'seq_resume_store_comp'

      nullify(fname_ptr)

      if (.not. associated(resume)) then
         allocate(resume)
      end if

      if (len_trim(filename) > 0) then
         num_inst = num_inst_comp * num_inst_driver
      else
         num_inst = 0
      end if

      ! Make sure each comp field is allocated correctly
      select case(oid)
      case ('a')
         if (associated(resume%atm_resume)) then
            if ((num_inst == 0) .or. (size(resume%atm_resume) /= num_inst)) then
               deallocate(resume%atm_resume)
               nullify(resume%atm_resume)
            end if
         end if
         if (num_inst > 0) then
            if (.not. associated(resume%atm_resume)) then
               allocate(resume%atm_resume(num_inst))
            end if
            fname_ptr => resume%atm_resume
         end if
      case ('l')
         if (associated(resume%lnd_resume)) then
            if ((num_inst == 0) .or. (size(resume%lnd_resume) /= num_inst)) then
               deallocate(resume%lnd_resume)
               nullify(resume%lnd_resume)
            end if
         end if
         if (num_inst > 0) then
            if (.not. associated(resume%lnd_resume)) then
               allocate(resume%lnd_resume(num_inst))
            end if
            fname_ptr => resume%lnd_resume
         end if
      case ('o')
         if (associated(resume%ocn_resume)) then
            if ((num_inst == 0) .or. (size(resume%ocn_resume) /= num_inst)) then
               deallocate(resume%ocn_resume)
               nullify(resume%ocn_resume)
            end if
         end if
         if (num_inst > 0) then
            if (.not. associated(resume%ocn_resume)) then
               allocate(resume%ocn_resume(num_inst))
            end if
            fname_ptr => resume%ocn_resume
         end if
      case ('i')
         if (associated(resume%ice_resume)) then
            if ((num_inst == 0) .or. (size(resume%ice_resume) /= num_inst)) then
               deallocate(resume%ice_resume)
               nullify(resume%ice_resume)
            end if
         end if
         if (num_inst > 0) then
            if (.not. associated(resume%ice_resume)) then
               allocate(resume%ice_resume(num_inst))
            end if
            fname_ptr => resume%ice_resume
         end if
      case ('r')
         if (associated(resume%rof_resume)) then
            if ((num_inst == 0) .or. (size(resume%rof_resume) /= num_inst)) then
               deallocate(resume%rof_resume)
               nullify(resume%rof_resume)
            end if
         end if
         if (num_inst > 0) then
            if (.not. associated(resume%rof_resume)) then
               allocate(resume%rof_resume(num_inst))
            end if
            fname_ptr => resume%rof_resume
         end if
      case ('g')
         if (associated(resume%glc_resume)) then
            if ((num_inst == 0) .or. (size(resume%glc_resume) /= num_inst)) then
               deallocate(resume%glc_resume)
               nullify(resume%glc_resume)
            end if
         end if
         if (num_inst > 0) then
            if (.not. associated(resume%glc_resume)) then
               allocate(resume%glc_resume(num_inst))
            end if
            fname_ptr => resume%glc_resume
         end if
      case ('w')
         if (associated(resume%wav_resume)) then
            if ((num_inst == 0) .or. (size(resume%wav_resume) /= num_inst)) then
               deallocate(resume%wav_resume)
               nullify(resume%wav_resume)
            end if
         end if
         if (num_inst > 0) then
            if (.not. associated(resume%wav_resume)) then
               allocate(resume%wav_resume(num_inst))
            end if
            fname_ptr => resume%wav_resume
         end if
      case ('x')
         if (associated(resume%cpl_resume)) then
            if ((num_inst == 0) .or. (size(resume%cpl_resume) /= num_inst)) then
               deallocate(resume%cpl_resume)
               nullify(resume%cpl_resume)
            end if
         end if
         if (num_inst > 0) then
            if (.not. associated(resume%cpl_resume)) then
               allocate(resume%cpl_resume(num_inst))
            end if
            fname_ptr => resume%cpl_resume
         end if
      case default
         call shr_sys_abort(subname//': Bad component id, '//oid)
      end select

      ! Copy in the resume filename if it exists
      if (associated(fname_ptr)) then
         fname_ptr(ninst) = filename
      end if

   end subroutine seq_resume_store_comp

   !===========================================================================
   !BOP =======================================================================
   !
   ! !IROUTINE: seq_resume_get_files -- Return resume filename info
   !
   ! !DESCRIPTION:
   !
   ! Return resume filename info
   !
   ! !INTERFACE: --------------------------------------------------------------
   subroutine seq_resume_get_files(oneletterid, files, bcast)
      character(len=1),           intent(in) :: oneletterid
      character(len=*), pointer              :: files(:)
      logical,          optional, intent(in) :: bcast

      character(len=*), parameter  :: subname = 'seq_resume_get_files'

      nullify(files)
      if (present(bcast)) then
         if (bcast) then
            call seq_resume_broadcast(oneletterid)
         end if
      ! No else: if not present, assume false
      end if
      select case(oneletterid)
      case ('a')
         if (associated(resume%atm_resume)) then
            files => resume%atm_resume
         end if
      case ('l')
         if (associated(resume%lnd_resume)) then
            files => resume%lnd_resume
         end if
      case ('o')
         if (associated(resume%ocn_resume)) then
            files => resume%ocn_resume
         end if
      case ('i')
         if (associated(resume%ice_resume)) then
            files => resume%ice_resume
         end if
      case ('r')
         if (associated(resume%rof_resume)) then
            files => resume%rof_resume
         end if
      case ('g')
         if (associated(resume%glc_resume)) then
            files => resume%glc_resume
         end if
      case ('w')
         if (associated(resume%wav_resume)) then
            files => resume%wav_resume
         end if
      case ('x')
         if (associated(resume%cpl_resume)) then
            files => resume%cpl_resume
         end if
      case default
         call shr_sys_abort(subname//': Bad component id, '//oneletterid)
      end select
   end subroutine seq_resume_get_files

   !===========================================================================
   !BOP =======================================================================
   !
   ! !IROUTINE: seq_resume_free -- Free space for resume filenames
   !
   ! !DESCRIPTION:
   !
   ! Free data for resume filenames from all component instances
   !
   ! !INTERFACE: --------------------------------------------------------------

   subroutine seq_resume_free()

      if (associated(resume)) then
         if (associated(resume%atm_resume)) then
            deallocate(resume%atm_resume)
            nullify(resume%atm_resume)
         end if

         if (associated(resume%lnd_resume)) then
            deallocate(resume%lnd_resume)
            nullify(resume%lnd_resume)
         end if

         if (associated(resume%ocn_resume)) then
            deallocate(resume%ocn_resume)
            nullify(resume%ocn_resume)
         end if

         if (associated(resume%ice_resume)) then
            deallocate(resume%ice_resume)
            nullify(resume%ice_resume)
         end if

         if (associated(resume%rof_resume)) then
            deallocate(resume%rof_resume)
            nullify(resume%rof_resume)
         end if

         if (associated(resume%glc_resume)) then
            deallocate(resume%glc_resume)
            nullify(resume%glc_resume)
         end if

         if (associated(resume%wav_resume)) then
            deallocate(resume%wav_resume)
            nullify(resume%wav_resume)
         end if

         if (associated(resume%cpl_resume)) then
            deallocate(resume%cpl_resume)
            nullify(resume%cpl_resume)
         end if
      end if
   end subroutine seq_resume_free

   !===========================================================================
   !BOP =======================================================================
   !
   ! !IROUTINE: seq_resume_broadcast
   !
   ! !DESCRIPTION:
   !
   ! Broadcast a component type's resume filenames to all PEs
   !
   ! !INTERFACE: --------------------------------------------------------------

   subroutine seq_resume_broadcast(oneletterid)
      character(len=1), intent(in) :: oneletterid

      character(len=CL), pointer   :: fname_ptr(:)
      character(len=*),  parameter :: subname = 'seq_resume_broadcast'

      ! This interface does a pointer dance. Because the array
      ! passed to seq_resume_broadcast_array may be NULL on input
      ! but allocated on output, we need to 'reconnect' it to the
      ! resume structure
      select case(oneletterid)
      case ('a')
         fname_ptr => resume%atm_resume
         call seq_resume_broadcast_array(fname_ptr)
         resume%atm_resume => fname_ptr
      case ('l')
         fname_ptr => resume%lnd_resume
         call seq_resume_broadcast_array(fname_ptr)
         resume%lnd_resume => fname_ptr
      case ('o')
         fname_ptr => resume%ocn_resume
         call seq_resume_broadcast_array(fname_ptr)
         resume%ocn_resume => fname_ptr
      case ('i')
         fname_ptr => resume%ice_resume
         call seq_resume_broadcast_array(fname_ptr)
         resume%ice_resume => fname_ptr
      case ('r')
         fname_ptr => resume%rof_resume
         call seq_resume_broadcast_array(fname_ptr)
         resume%rof_resume => fname_ptr
      case ('g')
         fname_ptr => resume%glc_resume
         call seq_resume_broadcast_array(fname_ptr)
         resume%glc_resume => fname_ptr
      case ('w')
         fname_ptr => resume%wav_resume
         call seq_resume_broadcast_array(fname_ptr)
         resume%wav_resume => fname_ptr
      case ('x')
         fname_ptr => resume%cpl_resume
         call seq_resume_broadcast_array(fname_ptr)
         resume%cpl_resume => fname_ptr
      case default
         call shr_sys_abort(subname//': Bad component id, '//oneletterid)
      end select
   end subroutine seq_resume_broadcast

   subroutine seq_resume_broadcast_array(filename_array)
     use shr_mpi_mod,  only: shr_mpi_bcast
     ! Used to bcast component filenames across multiple drivers
     use seq_comm_mct, only: global_comm

     character(len=CL), pointer   :: filename_array(:)

     integer, allocatable         :: active_entries(:)
     integer                      :: global_numpes
     integer                      :: num_entries
     integer                      :: my_entry
     integer                      :: index
     integer                      :: ierr
     character(len=128)           :: errmsg
     character(len=*),  parameter :: subname = "(fill_array_pes)"

     call MPI_comm_rank(global_comm, global_numpes, ierr)
     allocate(active_entries(global_numpes))

     ! Find filled array element (if any)
     ! Note, it is an error to find more than one.
     active_entries = 0
     my_entry = 0
     if (associated(filename_array)) then
        do index = 1, size(filename_array)
           if (len_trim(filename_array(index)) > 0) then
              if (my_entry > 0) then
                 write(errmsg, '(2(a,i0))') ': Bad entry, ', index,           &
                      ', already have ',my_entry
                 call shr_sys_abort(subname//trim(errmsg))
              end if
              my_entry = index
           end if
        end do
     end if
     ! Share my_entry with other PEs
     call MPI_allgather(my_entry, 1, MPI_INTEGER, active_entries, 1, MPI_INTEGER, global_comm, ierr)
     ! Allocate our array if needed
     num_entries = MAXVAL(active_entries)
     if ((num_entries > 0) .and. (.not. associated(filename_array))) then
        allocate(filename_array(num_entries))
     end if
     do index = 1, global_numpes
        my_entry = active_entries(index)
        if (my_entry > 0) then
           call shr_mpi_bcast(filename_array(my_entry), global_comm,          &
                subname//': bcast', pebcast=index-1)
        end if
     end do

     deallocate(active_entries)
  end subroutine seq_resume_broadcast_array

end module seq_pauseresume_mod
