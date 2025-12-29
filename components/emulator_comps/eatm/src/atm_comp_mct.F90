module atm_comp_mct
   !----------------------------------------------------------------------------
   ! Thin Fortran wrapper for the atmosphere emulator component.
   !
   ! This module provides the MCT interface (atm_init_mct, atm_run_mct,
   ! atm_final_mct) and delegates all implementation to the C++ EmulatorAtm
   ! class via the emulator_f2c_mod interface.
   !
   ! This is part of the generalized Emulator Component Framework that supports
   ! atmosphere, ocean, ice, and other emulated E3SM components.
   !----------------------------------------------------------------------------

   use esmf
   use mct_mod
   use seq_cdata_mod,    only: seq_cdata, seq_cdata_setptrs
   use seq_infodata_mod, only: seq_infodata_type, seq_infodata_getdata, seq_infodata_putdata
   use seq_timemgr_mod,  only: seq_timemgr_EClockGetData
   use seq_flds_mod,     only: seq_flds_a2x_fields, seq_flds_x2a_fields
   use seq_comm_mct,     only: seq_comm_inst, seq_comm_name, seq_comm_suffix
   use shr_kind_mod,     only: IN=>SHR_KIND_IN, R8=>SHR_KIND_R8, CL=>SHR_KIND_CL
   use shr_file_mod,     only: shr_file_getunit, shr_file_setIO
   use iso_c_binding
   use emulator_atm_f2c_v2

   implicit none
   private

   !--------------------------------------------------------------------------
   ! Public interfaces
   !--------------------------------------------------------------------------
   public :: atm_init_mct
   public :: atm_run_mct
   public :: atm_final_mct

   !--------------------------------------------------------------------------
   ! Private module data
   !--------------------------------------------------------------------------
   integer                :: mpicom_atm          ! mpi communicator
   integer(IN)            :: my_task             ! my task in mpi communicator
   integer                :: inst_index          ! instance index
   character(len=16)      :: inst_name           ! instance name
   character(len=16)      :: inst_suffix = ""    ! instance suffix
   integer(IN)            :: ATM_ID              ! mct component id
   integer(IN), parameter :: master_task = 0     ! master task number
   integer                :: atm_log_unit

CONTAINS

   !===============================================================================
   subroutine atm_init_mct(EClock, cdata, x2a, a2x, NLFilename)
      !---------------------------------------------------------------------------
      ! Initialize the atmosphere emulator component
      !---------------------------------------------------------------------------

      ! Arguments
      type(ESMF_Clock),           intent(inout) :: EClock
      type(seq_cdata),            intent(inout) :: cdata
      type(mct_aVect),            intent(inout) :: x2a, a2x
      character(len=*), optional, intent(in)    :: NLFilename

      ! Local variables
      type(seq_infodata_type), pointer :: infodata
      type(mct_gsMap),         pointer :: gsMap_atm
      type(mct_gGrid),         pointer :: dom_atm
      integer(IN)                      :: phase
      integer(IN)                      :: ierr
      integer                          :: lsize
      integer(IN)                      :: cur_ymd, cur_tod
      integer(IN)                      :: case_start_ymd, case_start_tod
      character(CL)                    :: run_type
      character(kind=c_char, len=256), target :: input_file_c
      character(kind=c_char, len=256), target :: log_file_c
      character(len=256)               :: log_file_f
      integer(c_int)                   :: run_type_c
      !---------------------------------------------------------------------------

      ! Get data from cdata structure
      call seq_cdata_setptrs(cdata, &
         id=ATM_ID, &
         mpicom=mpicom_atm, &
         gsMap=gsMap_atm, &
         dom=dom_atm, &
         infodata=infodata)

      call seq_infodata_getData(infodata, atm_phase=phase, start_type=run_type)
      call seq_infodata_PutData(infodata, atm_prognostic=.true.)

      if (phase > 1) return

      ! Determine instance information
      inst_name   = seq_comm_name(ATM_ID)
      inst_index  = seq_comm_inst(ATM_ID)
      inst_suffix = seq_comm_suffix(ATM_ID)

      ! Get rank
      call MPI_Comm_rank(mpicom_atm, my_task, ierr)

      ! Setup log file
      if (my_task == master_task) then
         atm_log_unit = shr_file_getunit()
         call shr_file_setIO('atm_modelio.nml'//trim(inst_suffix), atm_log_unit)
      endif

      ! Determine run type
      if (trim(run_type) == 'startup') then
         run_type_c = 0
      else if (trim(run_type) == 'continue') then
         run_type_c = 1
      else if (trim(run_type) == 'branch') then
         run_type_c = 2
      else
         run_type_c = 0
      endif

      ! Get time info
      call seq_timemgr_EClockGetData(EClock, &
         curr_ymd=cur_ymd, curr_tod=cur_tod, &
         start_ymd=case_start_ymd, start_tod=case_start_tod)

      ! Create input file name (null-terminated for C++)
      ! This points to the YAML configuration file
      input_file_c = "atm_in"//trim(inst_suffix)//C_NULL_CHAR

      ! Get log file name (master task only)
      log_file_c = C_NULL_CHAR
      if (my_task == master_task) then
         inquire(unit=atm_log_unit, name=log_file_f)
         log_file_c = trim(log_file_f)//C_NULL_CHAR
      endif

      ! Call C++ to create the emulator instance
      call emulator_atm_create_instance( &
         mpicom_atm, ATM_ID, input_file_c, &
         log_file_c, &
         run_type_c, cur_ymd, cur_tod)

      ! Setup MCT gsMap using C++ grid info
      call atm_SetGSMap_mct(mpicom_atm, ATM_ID, gsMap_atm)
      lsize = mct_gsMap_lsize(gsMap_atm, mpicom_atm)

      ! Setup MCT domain
      call atm_domain_mct(lsize, gsMap_atm, dom_atm)

      ! Initialize MCT attribute vectors
      call mct_aVect_init(x2a, rList=seq_flds_x2a_fields, lsize=lsize)
      call mct_aVect_init(a2x, rList=seq_flds_a2x_fields, lsize=lsize)

      ! Zero-initialize the attribute vectors to prevent garbage values
      ! in fields we don't explicitly set
      call mct_aVect_zero(x2a)
      call mct_aVect_zero(a2x)

      ! Initialize coupling field indices in C++
      ! Pass the colon-separated field lists so C++ knows the index mapping
      call emulator_atm_init_coupling_indices( &
         trim(seq_flds_a2x_fields)//C_NULL_CHAR, &
         trim(seq_flds_x2a_fields)//C_NULL_CHAR)

      ! Setup coupling with C++ (pass pointers to mct data)
      call emulator_atm_setup_coupling( &
         c_loc(x2a%rAttr), c_loc(a2x%rAttr), &
         mct_aVect_nRattr(x2a), mct_aVect_nRattr(a2x), lsize)

      ! Initialize the emulator
      call emulator_atm_init()

      ! Set grid size info
      call seq_infodata_PutData(infodata, &
         atm_nx=emulator_atm_get_nx(), &
         atm_ny=emulator_atm_get_ny())

   end subroutine atm_init_mct

   !===============================================================================
   subroutine atm_run_mct(EClock, cdata, x2a, a2x)
      !---------------------------------------------------------------------------
      ! Run the atmosphere emulator for one timestep
      !---------------------------------------------------------------------------

      ! Arguments
      type(ESMF_Clock), intent(inout) :: EClock
      type(seq_cdata),  intent(inout) :: cdata
      type(mct_aVect),  intent(inout) :: x2a
      type(mct_aVect),  intent(inout) :: a2x

      ! Local variables
      type(seq_infodata_type), pointer :: infodata
      integer                          :: dt
      real(R8)                         :: nextsw_cday
      !---------------------------------------------------------------------------

      call seq_cdata_setptrs(cdata, infodata=infodata)
      call seq_timemgr_EClockGetData(EClock, next_cday=nextsw_cday, dtime=dt)

      ! Run the emulator
      call emulator_atm_run(dt)

      ! Set time of next radiation computation
      call seq_infodata_PutData(infodata, nextsw_cday=nextsw_cday)

   end subroutine atm_run_mct

   !===============================================================================
   subroutine atm_final_mct(EClock, cdata, x2a, a2x)
      !---------------------------------------------------------------------------
      ! Finalize the atmosphere emulator component
      !---------------------------------------------------------------------------

      ! Arguments
      type(ESMF_Clock), intent(inout) :: EClock
      type(seq_cdata),  intent(inout) :: cdata
      type(mct_aVect),  intent(inout) :: x2a
      type(mct_aVect),  intent(inout) :: a2x
      !---------------------------------------------------------------------------

      call emulator_atm_finalize()

   end subroutine atm_final_mct

   !===============================================================================
   subroutine atm_SetGSMap_mct(mpicom_atm, ATMID, GSMap_atm)
      !---------------------------------------------------------------------------
      ! Setup MCT global segment map using C++ grid info
      !---------------------------------------------------------------------------

      ! Arguments
      integer,         intent(in)  :: mpicom_atm
      integer,         intent(in)  :: ATMID
      type(mct_gsMap), intent(out) :: GSMap_atm

      ! Local variables
      integer(c_int) :: num_local_cols, num_global_cols
      integer(c_int), allocatable, target :: col_gids(:)
      !---------------------------------------------------------------------------

      num_local_cols  = emulator_atm_get_num_local_cols()
      num_global_cols = emulator_atm_get_num_global_cols()

      allocate(col_gids(num_local_cols))
      call emulator_atm_get_local_cols_gids(c_loc(col_gids))

      call mct_gsMap_init(GSMap_atm, col_gids, mpicom_atm, ATMID, &
         num_local_cols, num_global_cols)

      deallocate(col_gids)

   end subroutine atm_SetGSMap_mct

   !===============================================================================
   subroutine atm_domain_mct(lsize, gsMap_atm, dom_atm)
      !---------------------------------------------------------------------------
      ! Setup MCT domain using C++ grid info
      !---------------------------------------------------------------------------

      ! Arguments
      integer,         intent(in)    :: lsize
      type(mct_gsMap), intent(in)    :: gsMap_atm
      type(mct_gGrid), intent(inout) :: dom_atm

      ! Local variables
      integer,  pointer :: idata(:)
      real(R8), pointer :: data1(:), data2(:)
      !---------------------------------------------------------------------------

      allocate(data1(lsize))
      allocate(data2(lsize))

      ! Initialize MCT domain
      call mct_gGrid_init(GGrid=dom_atm, &
         CoordChars='lat:lon:hgt', &
         OtherChars='area:aream:mask:frac', &
         lsize=lsize)

      ! Get global IDs
      call mct_gsMap_orderedPoints(gsMap_atm, my_task, idata)
      call mct_gGrid_importIAttr(dom_atm, 'GlobGridNum', idata, lsize)

      ! Initialize with placeholder values
      data1(:) = -9999.0_R8
      call mct_gGrid_importRAttr(dom_atm, "lat", data1, lsize)
      call mct_gGrid_importRAttr(dom_atm, "lon", data1, lsize)
      call mct_gGrid_importRAttr(dom_atm, "area", data1, lsize)
      call mct_gGrid_importRAttr(dom_atm, "aream", data1, lsize)

      ! Get lat/lon from C++
      call emulator_atm_get_cols_latlon(c_loc(data1), c_loc(data2))
      call mct_gGrid_importRAttr(dom_atm, "lat", data1, lsize)
      call mct_gGrid_importRAttr(dom_atm, "lon", data2, lsize)

      ! Get area from C++
      call emulator_atm_get_cols_area(c_loc(data1))
      call mct_gGrid_importRAttr(dom_atm, "area", data1, lsize)
      call mct_gGrid_importRAttr(dom_atm, "aream", data1, lsize)

      ! Set mask and frac to 1
      ! TODO: get this from elsewhere?
      data1(:) = 1.0_R8
      call mct_gGrid_importRAttr(dom_atm, "mask", data1, lsize)
      call mct_gGrid_importRAttr(dom_atm, "frac", data1, lsize)

      deallocate(data1)
      deallocate(data2)
      deallocate(idata)

   end subroutine atm_domain_mct

end module atm_comp_mct
