module read_spa_data

  !====================================================================================
  ! These set of routines user tracer data codes to read and interpolate (space and
  ! time) spa (simplified prescribed aerosol) input data which is read in from
  ! an input file

  ! Author: Balwinder Singh
  !====================================================================================

  use shr_kind_mod, only: shr_kind_cl
  use tracer_data,  only: trfld, trfile ! data strutures for storing file and fields info
  use radconstants, only: nswbands, nlwbands ! short and long wave band numbers
  implicit none

  !declare all variables private by default
  private

  !list of public subroutines
  public :: read_spa_data_init, read_spa_data_adv, read_spa_data_register, spa_readnl

  !---------------------------------------------------------------------------
  !Module level data shared by more than one subroutine or protected variables
  !---------------------------------------------------------------------------

  logical, public, protected :: is_spa_active = .false.

  type(trfld), pointer :: spa_fields_type(:)     ! data struture to store spa fields related info (see tracer_data.F90 for details)
  type(trfile)         :: spa_file_type          ! data structure to store spa files related info (see tracer_data.F90 for details)

  !Following fields are read in via namelist but defaults are set here
  character(len=shr_kind_cl) :: filename = ' '
  character(len=shr_kind_cl) :: datapath = ' '
  character(len=32)          :: datatype = 'SERIAL'
  integer                    :: cycle_yr = 0


  integer, parameter :: N_CCN    = 1  !total number of CCN fields

  !nswbands is the number of sw bands and nlwbands is the number of lwbands
  integer, parameter :: N_FLDS   = nlwbands + (3*nswbands) + N_CCN !total number of fields to be read (3 SW arrays)

  !--------------------------------------------------
  !fields to be read from spa input file
  !--------------------------------------------------
  integer, parameter :: c_fld_len = 13 ! character field length (some compilers require fixed length characters)

  !Radiation LW and SW arrays
  character(len=c_fld_len), public, parameter :: aer_asm_sw_names(nswbands) = &
       ['AER_G_SW_0   ', 'AER_G_SW_1   ', 'AER_G_SW_2   ', 'AER_G_SW_3   ', 'AER_G_SW_4   ', &
       'AER_G_SW_5   ', 'AER_G_SW_6   ', 'AER_G_SW_7   ', 'AER_G_SW_8   ', 'AER_G_SW_9   ', &
       'AER_G_SW_10  ', 'AER_G_SW_11  ', 'AER_G_SW_12  ', 'AER_G_SW_13  ']
  character(len=c_fld_len), public, parameter :: aer_ssa_sw_names(nswbands) = &
       ['AER_SSA_SW_0 ', 'AER_SSA_SW_1 ', 'AER_SSA_SW_2 ', 'AER_SSA_SW_3 ', 'AER_SSA_SW_4 ', &
       'AER_SSA_SW_5 ', 'AER_SSA_SW_6 ', 'AER_SSA_SW_7 ', 'AER_SSA_SW_8 ', 'AER_SSA_SW_9 ', &
       'AER_SSA_SW_10', 'AER_SSA_SW_11', 'AER_SSA_SW_12', 'AER_SSA_SW_13']
  character(len=c_fld_len), public, parameter :: aer_tau_sw_names(nswbands) = &
       ['AER_TAU_SW_0 ', 'AER_TAU_SW_1 ', 'AER_TAU_SW_2 ', 'AER_TAU_SW_3 ', 'AER_TAU_SW_4 ', &
       'AER_TAU_SW_5 ', 'AER_TAU_SW_6 ', 'AER_TAU_SW_7 ', 'AER_TAU_SW_8 ', 'AER_TAU_SW_9 ', &
       'AER_TAU_SW_10', 'AER_TAU_SW_11', 'AER_TAU_SW_12','AER_TAU_SW_13']
  character(len=c_fld_len), public, parameter :: aer_tau_lw_names(nlwbands) = &
       ['AER_TAU_LW_0 ', 'AER_TAU_LW_1 ', 'AER_TAU_LW_2 ', 'AER_TAU_LW_3 ', 'AER_TAU_LW_4 ', &
       'AER_TAU_LW_5 ', 'AER_TAU_LW_6 ', 'AER_TAU_LW_7 ', 'AER_TAU_LW_8 ', 'AER_TAU_LW_9 ', &
       'AER_TAU_LW_10', 'AER_TAU_LW_11', 'AER_TAU_LW_12', 'AER_TAU_LW_13', 'AER_TAU_LW_14', &
       'AER_TAU_LW_15']

  !P3 array for CCN
  character(len=c_fld_len), public, parameter :: ccn_names (N_CCN) = ['CCN3         ']

  !combine all arrays into one array for reading this data using tracer data routines
  character(len=c_fld_len), parameter :: field_names(N_FLDS) = [ aer_asm_sw_names, aer_ssa_sw_names, &
       aer_tau_sw_names, aer_tau_lw_names, ccn_names]

contains

  !-------------------------------------------------------------------
  ! registers spa fields to the phys buffer
  !-------------------------------------------------------------------
  subroutine read_spa_data_register()

    use ppgrid,         only: pver,pcols
    use physics_buffer, only: pbuf_add_field, dtype_r8

    implicit none

    integer :: i, idx

    if(.not. is_spa_active) return

    !Add all the fields to physics buffer so that they can be retrieved later when needed
    do i = 1 , size(field_names)
       call pbuf_add_field(field_names(i),'physpkg',dtype_r8,(/pcols,pver/),idx)
    enddo

  endsubroutine read_spa_data_register

  !-------------------------------------------------------------------
  ! Read spa namelist
  !-------------------------------------------------------------------
  subroutine spa_readnl(nlfile)

    use shr_log_mod,     only: errMsg => shr_log_errMsg
    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use spmd_utils,      only: masterproc
    use mpishorthand,    only: mpiint, mpichar, mpicom
    use cam_abortutils,   only : endrun

    implicit none

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr

    character(len=shr_kind_cl) :: spa_file
    character(len=shr_kind_cl) :: spa_datapath
    character(len=32)          :: spa_type
    integer                    :: spa_cycle_yr

    !namelist fields to read
    namelist /spa_nl/ &
         spa_file,      & !file containing spa data
         spa_datapath,  & !path to the file
         spa_type,      & !Type of data(CYCLICAL,SERIAL,INTERP_MISSING_MONTHS,FIXED)
         spa_cycle_yr     !Year to cycle for the data

    !Initialize namelist variables from local module variables.
    spa_file     = filename
    spa_datapath = datapath
    spa_type     = datatype
    spa_cycle_yr = cycle_yr !BALLI replace cyle year with some other var name

    !Read namelist
    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'spa_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, spa_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun('ERROR reading spa namelist, '//errmsg(__FILE__,__LINE__))
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

#ifdef SPMD
   !Broadcast namelist variables
   call mpibcast(spa_file,     len(spa_file),     mpichar, 0, mpicom)
   call mpibcast(spa_datapath, len(spa_datapath), mpichar, 0, mpicom)
   call mpibcast(spa_type,     len(spa_type),     mpichar, 0, mpicom)
   call mpibcast(spa_cycle_yr, 1, mpiint,  0, mpicom)
#endif

   !Update module variables with user settings.
   filename   = spa_file
   datapath   = spa_datapath
   datatype   = spa_type
   cycle_yr   = spa_cycle_yr

   !Turn on spa if user has specified an input dataset.
   is_spa_active = len_trim(filename) > 0

  end subroutine spa_readnl

  !-------------------------------------------------------------------
  ! Initialize spa input data to be read from an input file
  !-------------------------------------------------------------------
  subroutine read_spa_data_init

    use tracer_data, only: trcdata_init

    implicit none

    !Local variables
    logical, parameter :: rmv_file = .false.                !if .true., local file will be removed after all data is read

    !Following variables are declared as parameters as this functionality is currently not implemented
    character(len=shr_kind_cl), parameter :: filelist = ' ' !not currently used
    integer, parameter :: fixed_ymd = 0                     !not currently used
    integer, parameter :: fixed_tod = 0                     !not currently used

    if(.not. is_spa_active) return

    !Allocate "in_pbuf" component of spa_file_type variable, so that tracer_init knows that these fields exist in PBUF
    allocate (spa_file_type%in_pbuf(size(field_names)))

    !set in_pbuf to true for all the fields since we added all the fields to pbuf during the "register" process
    spa_file_type%in_pbuf(:) = .true.

    !call the init routine so that tracer data routine can initialize all fields required to read and interpolated the data
    call trcdata_init( field_names, filename, filelist, datapath, spa_fields_type, spa_file_type, &
         rmv_file, cycle_yr, fixed_ymd, fixed_tod, datatype)


  end subroutine read_spa_data_init

  !-------------------------------------------------------------------
  ! Advance fields in time and interpolate both in space and time
  !-------------------------------------------------------------------
  subroutine read_spa_data_adv( state, pbuf2d )

    !advance fields in time and interpolate (space and time)
    use tracer_data,  only : advance_trcdata
    use physics_types,only : physics_state
    use ppgrid,       only : begchunk, endchunk

    use physics_buffer, only : physics_buffer_desc

    implicit none

    type(physics_state), intent(in)    :: state(begchunk:endchunk)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    if(.not. is_spa_active) return

    !interpolate in time and space
    call advance_trcdata( spa_fields_type, spa_file_type, state, pbuf2d )

  end subroutine read_spa_data_adv

end module read_spa_data
