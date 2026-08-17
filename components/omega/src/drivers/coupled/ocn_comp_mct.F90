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
   use mct_mod, only: mct_aVect, mct_gsMap, mct_gGrid, mct_aVect_init

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
   integer(IN) :: my_task ! my task in mpi communicator mpicom
   integer(IN), parameter :: master_task = 0 ! task number of master task

contains
   subroutine ocn_init_mct(EClock, cdata, x2o, o2x, NLFilename)
      use, intrinsic :: iso_c_binding, only: c_null_char
      use, intrinsic :: iso_c_binding, only: c_ptr, c_loc, c_int, c_char

      use omega_f2cxx_mod, only: &
         omega_ocn_init1, &
         omega_ocn_init2
#ifdef HAVE_MOAB
      use omega_f2cxx_mod, only: omega_get_moab_pid
      use seq_comm_mct, only: mpoid
#endif

      use omega_cpl_indices, only: &
         num_coupler_imports, &
         num_coupler_exports, &
         import_field_names, &
         export_field_names, &
         import_field_indices, &
         export_field_indices, &
         cpl_x2o_field_names, &
         cpl_o2x_field_names, &
         omega_set_cpl_indices

      use mct_mod, only: mct_gsMap_lsize
      use seq_comm_mct, only: seq_comm_suffix
      use seq_flds_mod, only: seq_flds_o2x_fields, seq_flds_x2o_fields
      use seq_infodata_mod, only: &
         seq_infodata_PutData, &
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
      type(mct_gGrid), pointer :: gGrid_ocn
      integer :: lsize
      integer(IN) :: ierr, mpi_ierr  ! error codes
      integer(IN) :: &
         coupling_time_step, case_start_tod, case_start_ymd, cur_tod, cur_ymd
      integer(kind=c_int) :: start_type_c
      character(kind=c_char, len=CL), target :: calendar_c
      character(kind=c_char, len=CL), target :: ocn_log_fname_c

      ! set cdata pointers
      call seq_cdata_setptrs( &
         cdata, &
         ID=OCN_ID, &
         mpicom=mpicom_ocn, &
         gsMap=gsMap_ocn, &
         dom=gGrid_ocn, &
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
      ! Initialize ocn
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
         dtime=coupling_time_step, &
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

      ! populate the import/export field name and index arrays
      call omega_set_cpl_indices()

      call omega_ocn_init1( &
         mpicom_ocn, &
         OCN_ID, &
         "omega.yml"//c_null_char, &
         trim(ocn_log_fname)//c_null_char, &
         start_type_c, &
         calendar_c, &
         case_start_ymd, &
         case_start_tod, &
         coupling_time_step, &
         num_coupler_imports, &
         num_coupler_exports, &
         size(import_field_names), &
         size(export_field_names), &
         c_loc(import_field_names), &
         c_loc(export_field_names), &
         c_loc(import_field_indices), &
         c_loc(export_field_indices), &
         c_loc(cpl_x2o_field_names), &
         c_loc(cpl_o2x_field_names) &
         )

#ifdef HAVE_MOAB
      ! tell the coupler-side migration/mapping code (cplcomp_exchange_mod,
      ! prep_ocn_mod) which MOAB app id is this ocean instance's own mesh
      mpoid = omega_get_moab_pid()
#endif

      !-------------------------------------------------------------------------
      ! initialize MCT gsmap, domain, and attribute vectors
      !-------------------------------------------------------------------------
      call ocn_set_gsmap_mct(mpicom_ocn, ocn_id, gsMap_ocn)

      lsize = mct_gsMap_lsize(gsMap_ocn, mpicom_ocn)

      call ocn_set_domain_mct(lsize, gsMap_ocn, gGrid_ocn)

      ! Init import/export mct attribute vectors
      call mct_aVect_init(x2o, rList=seq_flds_x2o_fields, lsize=lsize)
      call mct_aVect_init(o2x, rList=seq_flds_o2x_fields, lsize=lsize)

      ! coupler needs Omega's decomposition before it can size x2o/o2x, so
      ! attach/export/import/halo-update must wait until they're allocated
      call seq_infodata_PutData( &
         infodata, &
         ocn_prognostic=.true., &
         ocnrof_prognostic=.true., &
         ocn_c2_glcshelf=.false. &
         )

      ! TODO: Get case config info and add as MetaData to Omega

      ! Under HAVE_MOAB, omega_ocn_init2 ignores these MCT attribute-vector
      ! pointers and attaches its own MOAB-backed buffers instead (see
      ! omega_cxx2f_interface.cpp); they're still passed here unconditionally
      ! since x2o/o2x are always allocated above regardless of driver.
      call omega_ocn_init2(c_loc(x2o%rAttr), c_loc(o2x%rAttr))

   end subroutine ocn_init_mct

   subroutine ocn_run_mct(EClock, cdata, x2o, o2x)

      use, intrinsic :: iso_c_binding, only: c_bool

      use omega_f2cxx_mod, only: omega_ocn_run
      use seq_timemgr_mod, only: seq_timemgr_RestartAlarmIsOn

      ! !INPUT/OUTPUT PARAMETERS:
      type(ESMF_Clock), intent(inout) :: EClock
      type(seq_cdata), intent(inout) :: cdata
      type(mct_aVect), intent(inout) :: x2o, o2x

      !--- local variables ---
      logical(kind=c_bool) :: write_restart

      ! check if coupler is requesting a restart write
      write_restart = logical( &
                      seq_timemgr_RestartAlarmIsOn(EClock), kind=c_bool)

      ! run omega for one coupling interval
      call omega_ocn_run(write_restart)

   end subroutine ocn_run_mct

   subroutine ocn_final_mct(EClock, cdata, x2o, o2x)
      use omega_f2cxx_mod, only: omega_ocn_finalize

      ! !INPUT/OUTPUT PARAMETERS:
      type(ESMF_Clock), intent(inout) :: EClock
      type(seq_cdata), intent(inout) :: cdata
      type(mct_aVect), intent(inout) :: x2o, o2x

      call omega_ocn_finalize()

   end subroutine ocn_final_mct

   subroutine ocn_set_gsmap_mct(mpicom_ocn, ocn_id, gsMap_ocn)
      use, intrinsic :: iso_c_binding, only: c_int
      use mct_mod, only: mct_gsMap_init
      use omega_f2cxx_mod, only: &
         omega_get_ncells_local, &
         omega_get_ncells_global, &
         omega_get_index_to_cell_id

      ! !INPUT/OUTPUT PARAMETERS:
      integer(IN), intent(in) :: mpicom_ocn
      integer(in), intent(in) :: OCN_ID
      type(mct_gsMap), intent(out) :: gsMap_ocn

      !--- local variables ---
      integer(kind=c_int), allocatable, target :: index_to_cell_id(:)
      integer(kind=c_int) :: ncells_local, ncells_global

      ! build the ocn cell numbering scheme for mct

      ncells_local = omega_get_ncells_local()
      ncells_global = omega_get_ncells_global()

      allocate (index_to_cell_id(ncells_local))

      call omega_get_index_to_cell_id(index_to_cell_id)

      call mct_gsMap_init( &
         gsmap_ocn, &
         index_to_cell_id, &
         mpicom_ocn, &
         ocn_id, &
         ncells_local, &
         ncells_global &
         )

      deallocate (index_to_cell_id)
   end subroutine ocn_set_gsmap_mct

   subroutine ocn_set_domain_mct(lsize, gsmap_ocn, ggrid_ocn)
      use, intrinsic :: iso_c_binding, only: c_int, c_double
      use mct_mod, only: &
         mct_gsMap, &
         mct_gGrid, &
         mct_gGrid_init, &
         mct_gGrid_importIAttr, &
         mct_gGrid_importRAttr, &
         mct_gsMap_orderedPoints
      use omega_f2cxx_mod, only: &
         omega_get_area_cell, &
         omega_get_lonlat_cell
      use seq_flds_mod, only: seq_flds_dom_coord, seq_flds_dom_other

      ! !INPUT/OUTPUT PARAMETERS:
      integer(IN), intent(in) :: lsize
      type(mct_gsMap), intent(in) :: gsmap_ocn
      type(mct_gGrid), intent(out) :: ggrid_ocn

      !--- local variables ---
      integer, pointer :: idata(:)
      real(kind=c_double), pointer :: data1(:), data2(:)

      allocate (data1(lsize))
      allocate (data2(lsize))

      ! initialize mct ocn domain
      call mct_gGrid_init( &
         GGrid=ggrid_ocn, &
         CoordChars=trim(seq_flds_dom_coord), &
         OtherChars=trim(seq_flds_dom_other), &
         lsize=lsize &
         )

      ! Initialize attribute vector with special value
      call mct_gsMap_orderedPoints(gsmap_ocn, my_task, idata)
      call mct_gGrid_importIAttr(ggrid_ocn, "GlobGridNum", idata, lsize)

      ! Fill in correct values for domain components
      call omega_get_lonlat_cell(data1, data2)
      call mct_gGrid_importRAttr(ggrid_ocn, "lon", data1, lsize)
      call mct_gGrid_importRAttr(ggrid_ocn, "lat", data2, lsize)

      call omega_get_area_cell(data1)
      call mct_gGrid_importRAttr(ggrid_ocn, "area", data1, lsize)

      ! mask and frac are both exactly 1, until landIceMask is suported
      data1(:) = real(1.0, kind=c_double)
      call mct_gGrid_importRAttr(ggrid_ocn, "mask", data1, lsize)
      call mct_gGrid_importRAttr(ggrid_ocn, "frac", data1, lsize)

      ! aream is computed by mct, so give invalid initial value
      data1(:) = real(-9999.0, kind=c_double)
      call mct_gGrid_importRAttr(ggrid_ocn, "aream", data1, lsize)

   end subroutine ocn_set_domain_mct
end module ocn_comp_mct
