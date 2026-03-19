module atm_comp_mct
   !------------------------------------------------------------------------
   ! Thin Fortran wrapper for the atmosphere emulator component.
   !
   ! Provides the MCT interface (atm_init_mct, atm_run_mct, atm_final_mct)
   ! and delegates all implementation to the C++ EmulatorAtm class via
   ! the emulatoratm_f2c interface module.
   !------------------------------------------------------------------------

   use esmf
   use mct_mod
   use seq_cdata_mod,    only: seq_cdata, seq_cdata_setptrs
   use seq_infodata_mod, only: seq_infodata_type, &
                               seq_infodata_getdata, &
                               seq_infodata_putdata
   use seq_timemgr_mod,  only: seq_timemgr_EClockGetData
   use seq_flds_mod,     only: seq_flds_a2x_fields, &
                               seq_flds_x2a_fields
   use seq_comm_mct,     only: seq_comm_inst, seq_comm_name, &
                               seq_comm_suffix
   use shr_kind_mod,     only: IN=>SHR_KIND_IN, R8=>SHR_KIND_R8, &
                               CL=>SHR_KIND_CL
   use shr_file_mod,     only: shr_file_getunit, shr_file_setIO, &
                               shr_file_getLogUnit, &
                               shr_file_setLogUnit
   use shr_sys_mod,      only: shr_sys_flush
   use iso_c_binding
   use emulator_f_api,  only : emulator_grid_desc, emulator_create_cfg, &
                               emulator_coupling_desc, &
                               create_config, create_grid_desc
   use emulator_f2c_api

   implicit none
   private

   public :: atm_init_mct
   public :: atm_run_mct
   public :: atm_final_mct

   integer                :: mpicom_atm
   integer(IN)            :: my_task
   integer                :: inst_index
   character(len=16)      :: inst_name
   character(len=16)      :: inst_suffix = ""
   integer(IN)            :: ATM_ID
   integer(IN), parameter :: master_task = 0
   integer                :: atm_log_unit

   type(c_ptr) :: emulators
CONTAINS

   !==========================================================================
   subroutine atm_init_mct(EClock, cdata, x2a, a2x, NLFilename)

      type(ESMF_Clock),           intent(inout) :: EClock
      type(seq_cdata),            intent(inout) :: cdata
      type(mct_aVect),            intent(inout) :: x2a, a2x
      character(len=*), optional, intent(in)    :: NLFilename

      type(seq_infodata_type), pointer :: infodata
      type(mct_gsMap),         pointer :: gsMap_atm
      type(mct_gGrid),         pointer :: dom_atm
      integer(IN) :: phase, ierr, lsize
      integer(IN) :: cur_ymd, cur_tod
      character(CL) :: run_type
      character(kind=c_char, len=256), target :: input_file_c
      character(kind=c_char, len=256), target :: log_file_c
      character(kind=c_char, len=256), target :: import_fields_c
      character(kind=c_char, len=256), target :: export_fields_c
      character(len=256) :: log_file_f
      integer(c_int) :: run_type_c
      integer :: shrlogunit
      type(emulator_create_cfg)    :: cfg
      type(emulator_coupling_desc) :: cpl

      call seq_cdata_setptrs(cdata, &
         id=ATM_ID, mpicom=mpicom_atm, &
         gsMap=gsMap_atm, dom=dom_atm, &
         infodata=infodata)

      call seq_infodata_getData(infodata, &
         atm_phase=phase, start_type=run_type)
      call seq_infodata_PutData(infodata, atm_prognostic=.true.)

      if (phase > 1) return

      inst_name   = seq_comm_name(ATM_ID)
      inst_index  = seq_comm_inst(ATM_ID)
      inst_suffix = seq_comm_suffix(ATM_ID)

      call MPI_Comm_rank(mpicom_atm, my_task, ierr)

      if (my_task == master_task) then
         atm_log_unit = shr_file_getunit()
         call shr_file_setIO( &
            'atm_modelio.nml'//trim(inst_suffix), atm_log_unit)
      else
         atm_log_unit = 6
      endif

      ! Redirect shared-lib logging to atm log file
      call shr_file_getLogUnit(shrlogunit)
      call shr_file_setLogUnit(atm_log_unit)

      if (my_task == master_task) then
         write(atm_log_unit,*) '(emulatoratm) ====================='
         write(atm_log_unit,*) '(emulatoratm) atm_init_mct starting'
         write(atm_log_unit,*) '(emulatoratm) ====================='
         call shr_sys_flush(atm_log_unit)
      endif

      if (trim(run_type) == 'startup') then
         run_type_c = 0
      else if (trim(run_type) == 'continue') then
         run_type_c = 1
      else if (trim(run_type) == 'branch') then
         run_type_c = 2
      else
         run_type_c = 0
      endif

      call seq_timemgr_EClockGetData(EClock, &
         curr_ymd=cur_ymd, curr_tod=cur_tod)

      input_file_c = "atm_in"//trim(inst_suffix)//C_NULL_CHAR

      log_file_c = C_NULL_CHAR
      if (my_task == master_task) then
         inquire(unit=atm_log_unit, name=log_file_f)
         log_file_c = trim(log_file_f)//C_NULL_CHAR
      endif

      cfg = create_config(f_comm=mpicom_atm, comp_id=ATM_ID, &
            run_type=run_type_c, start_ymd=cur_ymd, start_tod=cur_tod, &
            input_file=input_file_c, log_file=log_file_c)

      emulators = emulator_create(kind='atm'//C_NULL_CHAR, cfg=cfg)

      ! Setup MCT gsMap
      call atm_SetGSMap_mct(emulators, mpicom_atm, ATM_ID, gsMap_atm)
      lsize = mct_gsMap_lsize(gsMap_atm, mpicom_atm)

      ! Setup MCT domain
      call atm_domain_mct(emulators, lsize, gsMap_atm, dom_atm)

      ! Initialize MCT attribute vectors
      call mct_aVect_init(x2a, rList=seq_flds_x2a_fields, lsize=lsize)
      call mct_aVect_init(a2x, rList=seq_flds_a2x_fields, lsize=lsize)
      call mct_aVect_zero(x2a)
      call mct_aVect_zero(a2x)

      ! Build null-terminated field strings for C
      export_fields_c = trim(seq_flds_a2x_fields)//C_NULL_CHAR
      import_fields_c = trim(seq_flds_x2a_fields)//C_NULL_CHAR

      ! Initialize coupling indices in C++
      call emulator_init_coupling_indices(emulators, &
         c_loc(import_fields_c), c_loc(export_fields_c))

      ! Build coupling descriptor directly (x2a%rAttr is rank-2: nflds x lsize)
      cpl%import_data = c_loc(x2a%rAttr(1,1))
      cpl%export_data = c_loc(a2x%rAttr(1,1))
      cpl%num_imports = mct_aVect_nRattr(x2a)
      cpl%num_exports = mct_aVect_nRattr(a2x)
      cpl%field_size  = lsize

      ! Setup coupling buffers
      call emulator_setup_coupling(emulators, cpl)

      ! Initialize the emulator
      call emulator_init(emulators)

      call seq_infodata_PutData(infodata, &
         atm_nx=emulator_get_nx(emulators), &
         atm_ny=emulator_get_ny(emulators))

      if (my_task == master_task) then
         write(atm_log_unit,*) '(emulatoratm) atm_init_mct complete'
         write(atm_log_unit,*) '(emulatoratm)   nx  =', &
            emulator_get_nx(emulators)
         write(atm_log_unit,*) '(emulatoratm)   ny  =', &
            emulator_get_ny(emulators)
         write(atm_log_unit,*) '(emulatoratm)   nloc=', &
            emulator_get_num_local_cols(emulators)
         write(atm_log_unit,*) '(emulatoratm)   nglb=', &
            emulator_get_num_global_cols(emulators)
         call shr_sys_flush(atm_log_unit)
      endif

      ! Restore shared-lib log unit
      call shr_file_setLogUnit(shrlogunit)

   end subroutine atm_init_mct

   !==========================================================================
   subroutine atm_run_mct(EClock, cdata, x2a, a2x)

      type(ESMF_Clock), intent(inout) :: EClock
      type(seq_cdata),  intent(inout) :: cdata
      type(mct_aVect),  intent(inout) :: x2a, a2x

      type(seq_infodata_type), pointer :: infodata
      integer  :: dt, shrlogunit
      real(R8) :: nextsw_cday

      call seq_cdata_setptrs(cdata, infodata=infodata)
      call seq_timemgr_EClockGetData(EClock, &
         next_cday=nextsw_cday, dtime=dt)

      ! Redirect shared-lib logging to atm log file
      call shr_file_getLogUnit(shrlogunit)
      call shr_file_setLogUnit(atm_log_unit)

      call emulator_run(emulators, int(dt,c_int))

      call seq_infodata_PutData(infodata, nextsw_cday=nextsw_cday)

      ! Restore shared-lib log unit
      call shr_file_setLogUnit(shrlogunit)

   end subroutine atm_run_mct

   !==========================================================================
   subroutine atm_final_mct(EClock, cdata, x2a, a2x)

      type(ESMF_Clock), intent(inout) :: EClock
      type(seq_cdata),  intent(inout) :: cdata
      type(mct_aVect),  intent(inout) :: x2a, a2x

      if (my_task == master_task) then
         write(atm_log_unit,*) '(emulatoratm) atm_final_mct starting'
         call shr_sys_flush(atm_log_unit)
      endif

      call emulator_finalize(emulators)

      if (my_task == master_task) then
         write(atm_log_unit,*) '(emulatoratm) atm_final_mct complete'
         call shr_sys_flush(atm_log_unit)
      endif

   end subroutine atm_final_mct

   !==========================================================================
   subroutine atm_SetGSMap_mct(emu, mpicom_atm, ATMID, GSMap_atm)

      use, intrinsic :: iso_c_binding, only : c_ptr, c_int

      type(c_ptr),     intent(in)  :: emu
      integer,         intent(in)  :: mpicom_atm
      integer,         intent(in)  :: ATMID
      type(mct_gsMap), intent(out) :: GSMap_atm

      integer(c_int) :: num_local_cols, num_global_cols
      integer(c_int), allocatable, target :: col_gids(:)

      ! Query grid sizes from emulator
      num_local_cols  = emulator_get_num_local_cols(emu)
      num_global_cols = emulator_get_num_global_cols(emu)

      allocate(col_gids(num_local_cols))

      ! Get local global-IDs from emulator
      call emulator_get_local_col_gids(emu, col_gids)
      call mct_gsMap_init(GSMap_atm, col_gids, mpicom_atm, &
         ATMID, num_local_cols, num_global_cols)

      deallocate(col_gids)

   end subroutine atm_SetGSMap_mct

   !==========================================================================
   subroutine atm_domain_mct(emu, lsize, gsMap_atm, dom_atm)

      use, intrinsic :: iso_c_binding, only : c_ptr, c_int
      use shr_kind_mod, only : R8 => SHR_KIND_R8

      type(c_ptr),      intent(in)    :: emu
      integer,          intent(in)    :: lsize
      type(mct_gsMap),  intent(in)    :: gsMap_atm
      type(mct_gGrid),  intent(inout) :: dom_atm

      integer,  pointer :: idata(:)
      real(R8), pointer :: data1(:), data2(:)

      allocate(data1(lsize))
      allocate(data2(lsize))

      call mct_gGrid_init(GGrid=dom_atm, &
         CoordChars='lat:lon:hgt', &
         OtherChars='area:aream:mask:frac', &
         lsize=lsize)

      ! Global column ids for local points
      call mct_gsMap_orderedPoints(gsMap_atm, my_task, idata)
      call mct_gGrid_importIAttr(dom_atm, 'GlobGridNum', idata, lsize)

      data1(:) = -9999.0_R8
      call mct_gGrid_importRAttr(dom_atm, "lat",  data1, lsize)
      call mct_gGrid_importRAttr(dom_atm, "lon",  data1, lsize)
      call mct_gGrid_importRAttr(dom_atm, "area", data1, lsize)
      call mct_gGrid_importRAttr(dom_atm, "aream",data1, lsize)

      ! Get geometry from emulator
      call emulator_get_cols_latlon(emu, data1, data2)
      call mct_gGrid_importRAttr(dom_atm, "lat", data1, lsize)
      call mct_gGrid_importRAttr(dom_atm, "lon", data2, lsize)

      call emulator_get_cols_area(emu, data1)
      call mct_gGrid_importRAttr(dom_atm, "area",  data1, lsize)
      call mct_gGrid_importRAttr(dom_atm, "aream", data1, lsize)

      data1(:) = 1.0_R8
      call mct_gGrid_importRAttr(dom_atm, "mask", data1, lsize)
      call mct_gGrid_importRAttr(dom_atm, "frac", data1, lsize)

      deallocate(data1)
      deallocate(data2)
      if (associated(idata)) deallocate(idata)

   end subroutine atm_domain_mct

end module atm_comp_mct
