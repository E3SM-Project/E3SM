module ocn_comp_mct

   ! shr mods
   use shr_kind_mod, only: &
      IN => SHR_KIND_IN, &
      CS => SHR_KIND_CS, &  ! short char
      CL => SHR_KIND_CX, &  ! long char
      CX => SHR_KIND_CX     ! extra-long char

   ! seq mods
   use seq_cdata_mod, only: seq_cdata, seq_cdata_setptrs
   use seq_infodata_mod, only: seq_infodata_type, seq_infodata_GetData
   use seq_timemgr_mod, only: seq_timemgr_EClockGetData

   ! toolkits mods
   use esmf, only: ESMF_Clock
   use mct_mod, only: mct_aVect, mct_gsMap, mct_gGrid

   ! MPI
   ! allow(use-all)
   use mpi

   implicit none
   save
   private

   !---------------------------------------------------------------------------
   ! Public interfaces
   !---------------------------------------------------------------------------
   public :: ocn_init_mct
   public :: ocn_run_mct
   public :: ocn_final_mct

   !--------------------------------------------------------------------------
   ! Private module data
   !--------------------------------------------------------------------------
   integer :: ocn_log_unit
   integer(IN), parameter :: master_task = 0       ! task number of master task

contains
   subroutine ocn_init_mct(EClock, cdata, x2o, o2x, NLFilename)
      use omega_f2cxx_mod, only: omega_ocn_init
      use, intrinsic :: iso_c_binding, only: c_null_char
      use, intrinsic :: iso_c_binding, only: c_ptr, c_loc, c_int, c_char

      use seq_comm_mct, only: seq_comm_suffix
      use seq_infodata_mod, only: &
         seq_infodata_start_type_start, &
         seq_infodata_start_type_cont, &
         seq_infodata_start_type_brnch
      use shr_sys_mod, only: shr_sys_flush
      use shr_cal_mod, only: shr_cal_noleap, shr_cal_gregorian
      use shr_file_mod, only: shr_file_getunit, shr_file_setIO

      ! !INPUT/OUTPUT PARAMETERS:
      type(ESMF_Clock), intent(inout) :: EClock
      type(seq_cdata), intent(inout) :: cdata
      type(mct_aVect), intent(inout) :: x2o, o2x
      character(len=*), optional, intent(in) :: NLFilename

      !--- local variables ---
      integer           :: mpicom_ocn  ! mpi communicator
      integer(IN)       :: my_task     ! my task in mpi communicator mpicom
      integer           :: inst_index  ! number of current instance (ie. 1)
      character(len=16), save :: &
         inst_suffix = ""  ! char str for instance (ie. "_0001" or "")
      integer(IN)       :: OCN_ID          ! mct comp id

      character(len=CS) :: start_type      ! "startup" | "hybrid" | "branch"
      character(len=CX) :: cpl_seq_option  ! coupling sequence option
      character(len=CX) :: ocn_log_fname   ! name of ocn log file
      character(len=CL) :: calendar        ! calendar string

      type(seq_infodata_type), pointer :: infodata
      type(mct_gsMap), pointer :: gsMap_ocn
      type(mct_gGrid), pointer :: dom_ocn
      integer(IN) :: ierr, mpi_ierr  ! error codes
      integer(IN) :: case_start_tod, case_start_ymd, cur_tod, cur_ymd
      integer(kind=c_int) :: start_type_c
      character(kind=c_char, len=CL), target :: calendar_c
      character(kind=c_char, len=CL), target :: ocn_log_fname_c

      ! set cdata pointers
      call seq_cdata_setptrs( &
         cdata, &
         ID=OCN_ID, &
         mpicom=mpicom_ocn, &
         gsMap=gsMap_ocn, &
         dom=dom_ocn, &
         infodata=infodata &
         )

      ! determine start type and coupling type
      call seq_infodata_GetData( &
         infodata, start_type=start_type, cpl_seq_option=cpl_seq_option &
         )

      ! Determine instance information
      inst_suffix = seq_comm_suffix(OCN_ID)

      ! Determine communicator group
      call mpi_comm_rank(mpicom_ocn, my_task, ierr)

      !------------------------------------------------------------------------
      ! Init ocn.log
      !------------------------------------------------------------------------
      if (my_task == master_task) then
         ocn_log_unit = shr_file_getUnit()
         call shr_file_setIO( &
            "ocn_modelio.nml"//trim(inst_suffix), ocn_log_unit &
            )
         inquire (unit=ocn_log_unit, name=ocn_log_fname)
      end if

      call mpi_bcast( &
         ocn_log_unit, 1, MPI_INTEGER, master_task, mpicom_ocn, mpi_ierr &
         )
      if (ierr /= 0) then
         print *, "[omega] ERROR broadcasting ocn.log unit"
         call mpi_abort(mpicom_ocn, ierr, mpi_ierr)
      end if

      call mpi_bcast( &
         ocn_log_fname, 256, MPI_CHARACTER, master_task, mpicom_ocn, ierr &
         )
      if (ierr /= 0) then
         print *, "[omega] ERROR broadcasting ocn.log file name"
         call mpi_abort(mpicom_ocn, ierr, mpi_ierr)
      end if

      !------------------------------------------------------------------------
      ! Initialize atm
      !------------------------------------------------------------------------
      if (trim(start_type) == trim(seq_infodata_start_type_start)) then
         start_type_c = 0
      else if (trim(start_type) == trim(seq_infodata_start_type_cont)) then
         start_type_c = 1
      else if (trim(start_type) == trim(seq_infodata_start_type_brnch)) then
         start_type_c = 2
      else
         print *, "[omega] ERROR! Unsupported run type: "//trim(start_type)
         call mpi_abort(mpicom_ocn, ierr, mpi_ierr)
      end if

      ! Get clock and calendar information
      call seq_timemgr_EClockGetData( &
         EClock, &
         calendar=calendar, &
         curr_ymd=cur_ymd, &
         curr_tod=cur_tod, &
         start_ymd=case_start_ymd, &
         start_tod=case_start_tod &
         )

      ! convert CIME calendar to a C string using Omega naming conventions
      if (trim(calendar) == trim(shr_cal_noleap)) then
         calendar_c = "No Leap"//c_null_char
      else if (trim(calendar) == trim(shr_cal_gregorian)) then
         calendar_c = "Gregorian"//c_null_char
      else
         print *, "[omega] ERROR! Unsupported calendar: "//trim(calendar)
         call mpi_abort(mpicom_ocn, ierr, mpi_ierr)
      end if

      call omega_ocn_init( &
         mpicom_ocn, &
         OCN_ID, &
         "omega.yml"//c_null_char, &
         trim(ocn_log_fname)//c_null_char, &
         calendar_c, &
         case_start_ymd, &
         case_start_tod &
         )

   end subroutine ocn_init_mct

   subroutine ocn_run_mct(EClock, cdata, x2o, o2x)
      ! !INPUT/OUTPUT PARAMETERS:
      type(ESMF_Clock), intent(inout) :: EClock
      type(seq_cdata), intent(inout) :: cdata
      type(mct_aVect), intent(inout) :: x2o, o2x
   end subroutine ocn_run_mct

   subroutine ocn_final_mct(EClock, cdata, x2o, o2x)
      ! !INPUT/OUTPUT PARAMETERS:
      type(ESMF_Clock), intent(inout) :: EClock
      type(seq_cdata), intent(inout) :: cdata
      type(mct_aVect), intent(inout) :: x2o, o2x
   end subroutine ocn_final_mct
end module ocn_comp_mct
