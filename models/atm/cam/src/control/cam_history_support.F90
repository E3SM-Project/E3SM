module cam_history_support
  use shr_kind_mod, only: r8=>shr_kind_r8, i8=>shr_kind_i8, shr_kind_cl
  use pio,          only: var_desc_t, file_desc_t
  use abortutils,   only: endrun
  use cam_logfile,  only: iulog
  use spmd_utils,   only: masterproc

  implicit none
  private

  integer, parameter, public :: max_string_len = 256   ! Length of strings
  integer, parameter, public :: max_chars = shr_kind_cl         ! max chars for char variables
  integer, parameter, public :: fieldname_len = 24   ! max chars for field name
  integer, parameter, public :: fieldname_suffix_len =  3 ! length of field name suffix ("&IC")
  integer, parameter, public :: fieldname_lenp2      = fieldname_len + 2 ! allow for extra characters
  integer, parameter, public :: max_fieldname_len    = fieldname_len + fieldname_suffix_len ! max chars for field name (including suffix)
  integer, parameter         :: max_mdimname_len     = 16

  integer, parameter, public :: pflds = 750               ! max number of fields for namelist entries fincl and fexcl
                                                   ! also used in write restart
  integer, parameter, public :: ptapes    = 8             ! max number of tapes


  type, public :: column_info
    character(len=max_chars) :: lat_name ! latitude name for this column or columns
    character(len=max_chars) :: lon_name ! latitude name for this column or columns
    integer :: num_lats            ! number of lats in a group of contiguous columns
    integer :: num_lons            ! number of lons in a group of contiguous columns
    integer :: columnlat(2)       ! beginning and ending latitude (range) dimensioned by groups
    integer :: columnlon(2)       ! beginning and ending longitude (range) dimensioned by groups
    integer, pointer :: hmask(:,:)
    integer, pointer :: hmask_dyn(:,:)
    integer :: column_cnt
    integer :: htape             ! need to keep track of which tape we are writing for iodesc reuse
  end type column_info

  type, public :: field_info

    ! (i.e., how often "outfld" is called):  every other; only during LW/SW
    ! radiation calcs; etc.
    logical :: flag_xyfill                    ! non-applicable xy points flagged with fillvalue
    integer, pointer :: mdims(:)              ! indicies into hist_coords list

    real(r8) :: fillvalue                     ! fillvalue for this variable, set to default if not explicit in addfld

    integer :: numlev                         ! vertical dimension (.nc file and internal arr)

    integer :: begdim1                        ! on-node dim1 start index
    integer :: enddim1                        ! on-node dim1 end index

    integer :: begdim2                        ! on-node dim2 start index
    integer :: enddim2                        ! on-node dim2 end index

    integer :: begdim3                        ! on-node chunk or lat start index
    integer :: enddim3                        ! on-node chunk or lat end index

    integer :: decomp_type                    ! type of decomposition (physics or dynamics)

    integer :: vector_compliment              ! id for vector compliment for interpolation
    integer, pointer :: colperdim3(:)         ! number of valid elements per chunk or lat

    character(len=max_fieldname_len) :: name     ! field name
    character(len=max_chars) :: long_name        ! long name
    character(len=max_chars) :: units            ! units
    character(len=max_chars) :: sampling_seq     ! sampling sequence - if not every timestep, how often field is sampled
  end type field_info

  real(r8), parameter, public :: fillvalue = 1.e36_r8     ! fill value for netcdf fields


  !---------------------------------------------------------------------------
  !
  ! hentry: elements of an entry in the list of active fields on a single history file
  !
  !---------------------------------------------------------------------------
  type, public:: hentry
    type (field_info)     :: field            ! field information
    character(len=1)      :: avgflag          ! averaging flag
    character(len=max_chars) :: time_op          ! time operator (e.g. max, min, avg)
    character(len=max_chars),pointer :: field_column_name(:) ! names of column groups written to tape

    integer :: hwrt_prec                      ! history output precision
    real(r8), pointer :: hbuf(:,:,:)
    type(var_desc_t), pointer :: varid(:)      ! variable ids
    integer, pointer :: nacs(:,:)             ! accumulation counter
    type(var_desc_t) :: nacs_varid 
  end type hentry

  !---------------------------------------------------------------------------
  !
  !  active_entry: vehicle for producing a ragged array
  !
  !---------------------------------------------------------------------------
  type, public:: active_entry


    type(hentry), pointer :: hlist(:)

    type (column_info),pointer :: column   (:) ! array of history tape column entries
    type (column_info),pointer :: column_st(:) ! array of history tape column entries for staggered grid (FV)


    !
    ! PIO ids
    !

    type(file_desc_t) :: File            ! PIO file id

    type(var_desc_t) :: mdtid            ! var id for timestep
    type(var_desc_t) :: ndbaseid         ! var id for base day
    type(var_desc_t) :: nsbaseid         ! var id for base seconds of base day
    type(var_desc_t) :: nbdateid         ! var id for base date
    type(var_desc_t) :: nbsecid          ! var id for base seconds of base date
    type(var_desc_t) :: ndcurid          ! var id for current day
    type(var_desc_t) :: nscurid          ! var id for current seconds of current day
    type(var_desc_t) :: dateid           ! var id for current date
    type(var_desc_t) :: co2vmrid         ! var id for co2 volume mixing ratio
    type(var_desc_t) :: ch4vmrid         ! var id for ch4 volume mixing ratio
    type(var_desc_t) :: n2ovmrid         ! var id for n2o volume mixing ratio
    type(var_desc_t) :: f11vmrid         ! var id for f11 volume mixing ratio
    type(var_desc_t) :: f12vmrid         ! var id for f12 volume mixing ratio
    type(var_desc_t) :: sol_tsiid        ! var id for total solar irradiance (W/m2)
    type(var_desc_t) :: datesecid        ! var id for curent seconds of current date
#if ( defined BFB_CAM_SCAM_IOP )
    type(var_desc_t) :: bdateid         ! var id for base date
    type(var_desc_t) :: tsecid        ! var id for curent seconds of current date
#endif
    type(var_desc_t) :: nstephid         ! var id for current timestep
    type(var_desc_t) :: timeid           ! var id for time
    type(var_desc_t) :: tbndid           ! var id for time_bnds
    type(var_desc_t) :: date_writtenid   ! var id for date time sample written
    type(var_desc_t) :: time_writtenid   ! var id for time time sample written
    type(var_desc_t) :: nlonid           ! var id for number of longitudes
    type(var_desc_t) :: wnummaxid        ! var id for cutoff fourier wavenumber (reduced grid)
    type(var_desc_t) :: f107id           ! var id for f107
    type(var_desc_t) :: f107aid          ! var id for f107a
    type(var_desc_t) :: kpid             ! var id for kp
    type(var_desc_t) :: apid             ! var id for ap

  end type active_entry

  !---------------------------------------------------------------------------
  !
  !  formula_terms_t: Information for formula terms (CF convention) variables
  !                   Used to add a formula-terms variable to the history file
  !                   Also adds a string, '<name>: <var_name>' to the parent
  !                   mdim's 'formula_terms' attribute.
  !
  !---------------------------------------------------------------------------
  type, public :: formula_terms_t
    character(len=max_fieldname_len) :: a_name = ''   ! 'A' term variable name
    character(len=max_string_len)    :: a_long_name = '' ! 'A' long name
    real(r8), pointer                :: a_values(:) => null() ! 'A' variable values
    character(len=max_fieldname_len) :: b_name = ''   ! 'B' term variable name
    character(len=max_string_len)    :: b_long_name = '' ! 'B' long name
    real(r8), pointer                :: b_values(:) => null() ! 'B' variable values
    character(len=max_fieldname_len) :: p0_name = ''  ! 'p0' term variable name
    character(len=max_string_len)    :: p0_long_name = '' ! 'p0' long name
    character(len=max_chars)         :: p0_units = '' ! 'p0' variable units
    real(r8)                         :: p0_value = fillvalue ! 'p0' variable values
    character(len=max_fieldname_len) :: ps_name = ''  ! 'ps' term variable name
  end type formula_terms_t

  !---------------------------------------------------------------------------
  !
  !  hist_coord_t: Information for history variable dimension attributes
  !
  !---------------------------------------------------------------------------
  type, public :: hist_coord_t
    character(len=max_mdimname_len) :: name = ''  ! dimension name
    integer                  :: dimsize = 0       ! size of dimension
    character(len=max_chars) :: long_name = ''    ! 'long_name' attribute
    character(len=max_chars) :: units = ''        ! 'units' attribute
    character(len=max_chars) :: bounds_name = ''  ! 'bounds' attribute (& name of bounds variable)
    character(len=max_chars) :: standard_name = ''! 'standard_name' attribute
    character(len=4)         :: positive = ''     ! 'positive' attribute ('up' or 'down')
    integer,  pointer        :: integer_values(:) => null() ! dim values if integral
    real(r8), pointer        :: real_values(:) => null() ! dim values if real
    real(r8), pointer        :: bounds(:,:) => null() ! dim bounds
    type(formula_terms_t)    :: formula_terms     ! vars for formula terms
    logical                  :: integer_dim       ! .true. iff dim has integral values
  end type hist_coord_t

  integer, public                     :: registeredmdims=0, maxvarmdims=1
  character(len=9), parameter, public :: mdim_var_name = 'mdimnames'
  integer, parameter                  :: maxmdims = 25    ! arbitrary limit
  type(hist_coord_t), public, save    :: hist_coords(maxmdims)

  public             :: add_hist_coord, add_vert_coord
  public             :: write_hist_coord_attrs, write_hist_coord_vars
  public             :: lookup_hist_coord_indices
  public             :: get_hist_coord_index
  public             :: sec2hms, date2yyyymmdd

  interface add_hist_coord
    module procedure add_hist_coord_regonly
    module procedure add_hist_coord_int
    module procedure add_hist_coord_r8
  end interface

  interface assignment(=)
    module procedure field_copy
    module procedure formula_terms_copy
  end interface

  interface check_hist_coord
    ! NB: This is supposed to be a private interface
    ! check_hist_coord: returns 0 if <name> is not registered as an mdim
    !                   returns i if <name> is registered with compatible values
    !              calls endrun if <name> is registered with incompatible values
    ! Versions without the <name> argument return .true. or .false.
    module procedure check_hist_coord_char
    module procedure check_hist_coord_int
    module procedure check_hist_coord_int1
    module procedure check_hist_coord_r8
    module procedure check_hist_coord_r81
    module procedure check_hist_coord_r82
    module procedure check_hist_coord_ft
    module procedure check_hist_coord_all
  end interface

!!---------------------------------------------------------------------------

contains


  subroutine field_copy(f_out, f_in)
    type(field_info), intent(in) :: f_in
    type(field_info), intent(out) :: f_out

    f_out%flag_xyfill= f_in%flag_xyfill
    f_out%fillvalue= f_in%fillvalue
    f_out%numlev =  f_in%numlev                      ! vertical dimension (.nc file and internal arr)
    f_out%begdim1 = f_in%begdim1                     ! on-node dim1 start index
    f_out%enddim1 = f_in%enddim1                     ! on-node dim1 end index
    f_out%begdim2 = f_in%begdim2                     ! on-node dim2 start index
    f_out%enddim2 = f_in%enddim2                     ! on-node dim2 end index
    f_out%begdim3 = f_in%begdim3                     ! on-node chunk or lat start index
    f_out%enddim3 = f_in%enddim3                     ! on-node chunk or lat end index
    f_out%decomp_type = f_in%decomp_type             ! type of decomposition (physics or dynamics)

    f_out%vector_compliment = f_in%vector_compliment ! id for vector compliment for interpolation

    f_out%name = f_in%name                           ! field name
    f_out%long_name = f_in%long_name                 ! long name
    f_out%units = f_in%units                         ! units
    f_out%sampling_seq =  f_in%sampling_seq          ! sampling sequence - if not every timestep, how often field is sampled

    if(associated(f_in%mdims)) then
      f_out%mdims=>f_in%mdims
    else
      nullify(f_out%mdims)
    end if
    if(associated(f_in%colperdim3)) then
      f_out%colperdim3=>f_in%colperdim3
    else
      nullify(f_out%colperdim3)
    end if

  end subroutine field_copy

  subroutine formula_terms_copy(f_out, f_in)
    type(formula_terms_t), intent(in) :: f_in
    type(formula_terms_t), intent(out) :: f_out

    f_out%a_name = f_in%a_name
    f_out%a_long_name = f_in%a_long_name
    f_out%a_values => f_in%a_values
    f_out%b_name = f_in%b_name
    f_out%b_long_name = f_in%b_long_name
    f_out%b_values => f_in%b_values
    f_out%p0_name = f_in%p0_name
    f_out%p0_long_name = f_in%p0_long_name
    f_out%p0_units = f_in%p0_units
    f_out%p0_value = f_in%p0_value
    f_out%ps_name = f_in%ps_name
  end subroutine formula_terms_copy

  !!XXgoldyXX: I would like to put this in cam_pio_utils
  subroutine handle_pio_error(ierr, errorstr)
    use abortutils,   only: endrun
    use pio,          only: pio_noerr
    implicit none

    ! Input variables
    integer,          intent(in)  :: ierr
    character(len=*), intent(in)  :: errorstr

    ! Local variables
    character(len=256) :: errormsg

    if (ierr /= PIO_NOERR) then
      write(errormsg, *) trim(errorstr), ierr
      call endrun(errormsg)
    end if
    
  end subroutine handle_pio_error

  integer function get_hist_coord_index(mdimname)
    implicit none
    ! Input variables
    character(len=*), intent(in)            :: mdimname
    ! Local variables
    character(len=120) :: errormsg
    integer :: i
    
    get_hist_coord_index = -1
    do i=1,registeredmdims
      if(trim(mdimname) == trim(hist_coords(i)%name)) then
        get_hist_coord_index = i
        exit
      end if
    end do

  end function get_hist_coord_index

  ! Functions to check consistent term definition for hist coords
  logical function check_hist_coord_char(defined, input)
    implicit none

    ! Input variables
    character(len=*), intent(in)            :: defined
    character(len=*), intent(in), optional  :: input

    if (len_trim(defined) == 0) then
      ! In this case, we assume the current value is undefined so any input OK
      check_hist_coord_char = .true.
    else if (present(input)) then
      ! We have to match definitions
      check_hist_coord_char = (trim(input) == trim(defined))
    else
      ! Not sure here. We have a value and are redefining without one?
      check_hist_coord_char = .false.
    end if
  end function check_hist_coord_char

  logical function check_hist_coord_int(defined, input)
    implicit none

    ! Input variables
    integer, intent(in)            :: defined
    integer, intent(in), optional  :: input

    if (defined == 0) then
      ! In this case, we assume the current value is undefined so any input OK
      check_hist_coord_int = .true.
    else if (present(input)) then
      ! We have to match definitions
      check_hist_coord_int = (input == defined)
    else
      ! Not sure here. We have a value and are redefining without one?
      check_hist_coord_int = .false.
    end if
  end function check_hist_coord_int

  logical function check_hist_coord_int1(defined, input)
    implicit none

    ! Input variables
    integer,             pointer            :: defined(:)
    integer, intent(in),          optional  :: input(:)

    ! Local variables
    integer                                 :: i

    if (.not. associated(defined)) then
      ! In this case, we assume the current value is undefined so any input OK
      check_hist_coord_int1 = .true.
    else if (present(input)) then
      ! We have to match definitions
      check_hist_coord_int1 = (size(input) == size(defined))
    else
      ! Not sure here. We have a value and are redefining without one?
      check_hist_coord_int1 = .false.
    end if
    if (check_hist_coord_int1 .and. associated(defined)) then
      ! Need to check the values
      do i = 1, size(defined)
        if (defined(i) /= input(i)) then
          check_hist_coord_int1 = .false.
          exit
        end if
      end do
    end if
  end function check_hist_coord_int1

  logical function check_hist_coord_r8(defined, input)
    implicit none

    ! Input variables
    real(r8), intent(in)            :: defined
    real(r8), intent(in), optional  :: input

    if (defined == fillvalue) then
      ! In this case, we assume the current value is undefined so any input OK
      check_hist_coord_r8 = .true.
    else if (present(input)) then
      ! We have to match definitions
      check_hist_coord_r8 = (input == defined)
    else
      ! Not sure here. We have a value and are redefining without one?
      check_hist_coord_r8 = .false.
    end if
  end function check_hist_coord_r8

  logical function check_hist_coord_r81(defined, input)
    implicit none

    ! Input variables
    real(r8),             pointer            :: defined(:)
    real(r8), intent(in),          optional  :: input(:)

    ! Local variables
    integer                                  :: i

    if (.not. associated(defined)) then
      ! In this case, we assume the current value is undefined so any input OK
      check_hist_coord_r81 = .true.
    else if (present(input)) then
      ! We have to match definitions
      check_hist_coord_r81 = (size(input) == size(defined))
    else
      ! Not sure here. We have a value and are redefining without one?
      check_hist_coord_r81 = .false.
    end if
    if (check_hist_coord_r81 .and. associated(defined)) then
      ! Need to check the values
      do i = 1, size(defined)
        if (defined(i) /= input(i)) then
          check_hist_coord_r81 = .false.
          exit
        end if
      end do
    end if
  end function check_hist_coord_r81

  logical function check_hist_coord_r82(defined, input)
    implicit none

    ! Input variables
    real(r8),             pointer            :: defined(:,:)
    real(r8), intent(in),          optional  :: input(:,:)

    ! Local variables
    integer                                  :: i, j

    if (.not. associated(defined)) then
      ! In this case, we assume the current value is undefined so any input OK
      check_hist_coord_r82 = .true.
    else if (present(input)) then
      ! We have to match definitions
      check_hist_coord_r82 = ((size(input, 1) == size(defined, 1)) .and.    &
                              (size(input, 2) == size(defined, 2)))
    else
      ! Not sure here. We have a value and are redefining without one?
      check_hist_coord_r82 = .false.
    end if
    if (check_hist_coord_r82 .and. associated(defined)) then
      ! Need to check the values
      do j = 1, size(defined, 2)
        do i = 1, size(defined, 1)
          if (defined(i, j) /= input(i, j)) then
            check_hist_coord_r82 = .false.
            exit
          end if
        end do
      end do
    end if
  end function check_hist_coord_r82

  logical function check_hist_coord_ft(defined, input)
    implicit none

    ! Input variables
    type(formula_terms_t), intent(in)           :: defined
    type(formula_terms_t), intent(in), optional :: input

    ! We will assume that if formula_terms has been defined, a_name has a value
    if (len_trim(defined%a_name) == 0) then
      ! In this case, we assume the current value is undefined so any input OK
      check_hist_coord_ft = .true.
    else if (present(input)) then
      ! We have to match definitions
      ! Need to check the values
      check_hist_coord_ft =                                                   &
           check_hist_coord(defined%a_name,       input%a_name)         .and. &
           check_hist_coord(defined%a_long_name,  input%a_long_name)    .and. &
           check_hist_coord(defined%a_values,     input%a_values)       .and. &
           check_hist_coord(defined%b_name,       input%b_name)         .and. &
           check_hist_coord(defined%b_long_name,  input%b_long_name)    .and. &
           check_hist_coord(defined%b_values,     input%b_values)       .and. &
           check_hist_coord(defined%p0_name,      input%p0_name)        .and. &
           check_hist_coord(defined%p0_long_name, input%p0_long_name)   .and. &
           check_hist_coord(defined%p0_units,     input%p0_units)       .and. &
           check_hist_coord(defined%p0_value,     input%p0_value)       .and. &
           check_hist_coord(defined%ps_name,      input%ps_name)
    else
      ! Not sure here. We have a value and are redefining without one?
      check_hist_coord_ft = .false.
    end if
  end function check_hist_coord_ft

  ! check_hist_coord: returns 0 if <name> is not registered as a hist coord
  !                   returns i if <name> is registered with compatible values
  !                   calls endrun if <name> is registered with incompatible values
  integer function check_hist_coord_all(name, vlen, long_name, units, bounds, &
       i_values, r_values, bounds_name, positive, standard_name, formula_terms)
    implicit none

    ! Input variables
    character(len=*),      intent(in)            :: name
    integer,               intent(in)            :: vlen
    character(len=*),      intent(in),  optional :: long_name
    character(len=*),      intent(in),  optional :: units
    character(len=*),      intent(in),  optional :: bounds_name
    integer,               intent(in),  optional :: i_values(:)
    real(r8),              intent(in),  optional :: r_values(:)
    real(r8),              intent(in),  optional :: bounds(:,:)
    character(len=*),      intent(in),  optional :: positive
    character(len=*),      intent(in),  optional :: standard_name
    type(formula_terms_t), intent(in),  optional :: formula_terms

    ! Local variables
    character(len=120)                           :: errormsg
    integer                                      :: i

    i = get_hist_coord_index(trim(name))
    ! If i > 0, this mdim has already been registered
    if (i > 0) then
      check_hist_coord_all = i
      if (.not. check_hist_coord(hist_coords(i)%dimsize, vlen)) then
        write(errormsg, *) 'ERROR: Attempt to register dimension, '//trim(name)//' with incompatible size'
        call endrun(errormsg)
      end if
      if (.not. check_hist_coord(hist_coords(i)%long_name, long_name)) then
        write(errormsg, *) 'ERROR: Attempt to register dimension, ',trim(name),' with a different long_name'
        call endrun(errormsg)
      end if
      if (.not. check_hist_coord(hist_coords(i)%units, units)) then
        write(errormsg, *) 'ERROR: Attempt to register dimension, ',trim(name),' with different units'
        call endrun(errormsg)
      end if
      if (.not. check_hist_coord(hist_coords(i)%bounds_name, bounds_name)) then
        write(errormsg, *) 'ERROR: Attempt to register dimension, ',trim(name),' with a different bounds_name'
        call endrun(errormsg)
      end if
      if (.not. check_hist_coord(hist_coords(i)%standard_name, standard_name)) then
        write(errormsg, *) 'ERROR: Attempt to register dimension, ',trim(name),' with a different standard_name'
        call endrun(errormsg)
      end if
      if (.not. check_hist_coord(hist_coords(i)%positive, positive)) then
        write(errormsg, *) 'ERROR: Attempt to register dimension, ',trim(name),' with a different value of positive'
        call endrun(errormsg)
      end if
      ! Since the integer_dim defaults to .true., double check which to check
      if ((.not. hist_coords(i)%integer_dim) .or.                             &
           associated(hist_coords(i)%real_values)) then
        if (.not. check_hist_coord(hist_coords(i)%real_values, r_values)) then
          write(errormsg, *) 'ERROR: Attempt to register dimension, ',trim(name),' with different values'
          call endrun(errormsg)
        else if (present(i_values)) then
          write(errormsg, *) 'ERROR: Attempt to register integer values for real dimension'
          call endrun(errormsg)
        end if
      else
        if (.not. check_hist_coord(hist_coords(i)%integer_values, i_values)) then
          write(errormsg, *) 'ERROR: Attempt to register dimension, ',trim(name),' with different values'
          call endrun(errormsg)
        else if (present(i_values) .and. present(r_values)) then
          write(errormsg, *) 'ERROR: Attempt to register real values for integer dimension'
          call endrun(errormsg)
        end if
      end if
      if (.not. check_hist_coord(hist_coords(i)%bounds, bounds)) then
        write(errormsg, *) 'ERROR: Attempt to register dimension, ',trim(name),' with different bounds'
        call endrun(errormsg)
      end if
      if (.not. check_hist_coord(hist_coords(i)%formula_terms, formula_terms)) then
        write(errormsg, *) 'ERROR: Attempt to register dimension, ',trim(name),' with different formula_terms'
        call endrun(errormsg)
      end if
    else
      check_hist_coord_all = 0
    end if
  end function check_hist_coord_all

  subroutine add_hist_coord_regonly(name, index)
    implicit none

    ! Input variable
    character(len=*),  intent(in)    :: name
    integer, optional, intent(out)   :: index

    ! Local variables
    character(len=120)               :: errormsg
    integer                          :: i

    i = get_hist_coord_index(trim(name))
    ! If i > 0, this mdim has already been registered
    if (i <= 0) then
      registeredmdims = registeredmdims + 1
      if (registeredmdims > maxmdims) then
        call endrun('Too many dimensions in add_hist_coord.')
      end if
      if (len_trim(name) > max_mdimname_len) then
        write(errormsg,'(a,i3,a)') 'History coord name exceeds the ',         &
             max_mdimname_len, ' character length limit'
        call endrun(errormsg)
      end if
      hist_coords(registeredmdims)%name = trim(name)
      hist_coords(registeredmdims)%dimsize = 0
      hist_coords(registeredmdims)%long_name = ''
      hist_coords(registeredmdims)%units = ''
      hist_coords(registeredmdims)%integer_dim = .true.
      hist_coords(registeredmdims)%positive = ''
      hist_coords(registeredmdims)%standard_name = ''
      if (present(index)) then
        index = registeredmdims
      end if
    else
      if (present(index)) then
        index = i
      end if
    end if

  end subroutine add_hist_coord_regonly

  subroutine add_hist_coord_int(name, vlen, long_name, units, values,         &
       positive, standard_name)
    implicit none

    ! Input variables
    character(len=*), intent(in)                    :: name
    integer,          intent(in)                    :: vlen
    character(len=*), intent(in)                    :: long_name
    character(len=*), intent(in),          optional :: units
    integer,          intent(in),  target, optional :: values(:)
    character(len=*), intent(in),          optional :: positive
    character(len=*), intent(in),          optional :: standard_name

    ! Local variables
    integer                                         :: i

    ! First, check to see if it is OK to add this coord
    i = check_hist_coord(name, vlen=vlen, long_name=long_name, units=units,   &
         i_values=values, positive=positive, standard_name=standard_name)
    ! Register the name if necessary
    if (i == 0) then
      call add_hist_coord(trim(name), i)
      !  if(masterproc) write(iulog,*) 'Registering hist coord',name,'(',i,') with length: ',vlen
    end if

    ! Set the coord's values
    hist_coords(i)%dimsize = vlen
    if (len_trim(long_name) > max_chars) then
      if(masterproc) then
        write(iulog,*) 'WARNING: long_name for ',trim(name),' too long'
      end if
    end if
    hist_coords(i)%long_name = trim(long_name)
    if (present(units)) then
      hist_coords(i)%units = trim(units)
    else
      hist_coords(i)%units = ''
    end if
    hist_coords(i)%integer_dim = .true.
    if (present(values)) then
      hist_coords(i)%integer_values => values
    endif
    if (present(positive)) then
      hist_coords(i)%positive = trim(positive)
    end if
    if (present(standard_name)) then
      hist_coords(i)%standard_name = trim(standard_name)
    end if

  end subroutine add_hist_coord_int

  subroutine add_hist_coord_r8(name, vlen, long_name, units, values,         &
       bounds_name, bounds, positive, standard_name)
    implicit none

    ! Input variables
    character(len=*),      intent(in)                    :: name
    integer,               intent(in)                    :: vlen
    character(len=*),      intent(in)                    :: long_name
    character(len=*),      intent(in)                    :: units
    real(r8),              intent(in), target            :: values(:)
    character(len=*),      intent(in),          optional :: bounds_name
    real(r8),              intent(in), target,  optional :: bounds(:,:)
    character(len=*),      intent(in),          optional :: positive
    character(len=*),      intent(in),          optional :: standard_name

    ! Local variables
    character(len=120)                                   :: errormsg
    integer                                              :: i

    ! First, check to see if it is OK to add this coord
    i = check_hist_coord(name, vlen=vlen, long_name=long_name, units=units,   &
         r_values=values, positive=positive, standard_name=standard_name,     &
         bounds_name=bounds_name, bounds=bounds)
    ! Register the name if necessary
    if (i == 0) then
      call add_hist_coord(trim(name), i)
      !  if(masterproc) write(iulog,*) 'Registering hist coord',name,'(',i,') with length: ',vlen
    end if

    ! Set the coord's size
    hist_coords(i)%dimsize = vlen
    if (len_trim(long_name) > max_chars) then
      if(masterproc) then
        write(iulog,*) 'WARNING: long_name for ',trim(name),' too long'
      end if
    end if
    hist_coords(i)%long_name = trim(long_name)
    if (len_trim(units) > 0) then
      hist_coords(i)%units = trim(units)
    else
      hist_coords(i)%units = '1'
    end if
    hist_coords(i)%integer_dim = .false.
    hist_coords(i)%real_values => values
    if (present(positive)) then
      hist_coords(i)%positive = trim(positive)
    end if
    if (present(standard_name)) then
      hist_coords(i)%standard_name = trim(standard_name)
    end if
    if (present(bounds_name)) then
      hist_coords(i)%bounds_name = trim(bounds_name)
      if (.not. present(bounds)) then
        write(errormsg,*) 'bounds must be present for ',trim(bounds_name)
        call endrun(errormsg)
      end if
      hist_coords(i)%bounds => bounds
    else if (present(bounds)) then
      write(errormsg,*) 'bounds_name must be present for bounds values'
      call endrun(errormsg)
    else
      hist_coords(i)%bounds_name = ''
    end if

  end subroutine add_hist_coord_r8

  subroutine add_vert_coord(name, vlen, long_name, units, values,            &
       positive, standard_name, formula_terms)
    implicit none

    ! Input variables
    character(len=*),      intent(in)                    :: name
    integer,               intent(in)                    :: vlen
    character(len=*),      intent(in)                    :: long_name
    character(len=*),      intent(in)                    :: units
    real(r8),              intent(in), target            :: values(:)
    character(len=*),      intent(in),          optional :: positive
    character(len=*),      intent(in),          optional :: standard_name
    type(formula_terms_t), intent(in),          optional :: formula_terms

    ! Local variables
    character(len=120)                                   :: errormsg
    integer                                              :: i

    ! First, check to see if it is OK to add this coord
    i = check_hist_coord(name, vlen=vlen, long_name=long_name, units=units,   &
         r_values=values, positive=positive, standard_name=standard_name,     &
         formula_terms=formula_terms)
    ! Register the name and hist_coord values if necessary
    if (i == 0) then
      call add_hist_coord(trim(name), vlen, long_name, units, values,         &
           positive=positive, standard_name=standard_name)
      i = get_hist_coord_index(trim(name))
      !  if(masterproc) write(iulog,*) 'Registering hist coord',name,'(',i,') with length: ',vlen
    end if

    if (present(formula_terms)) then
      hist_coords(i)%formula_terms = formula_terms
    end if

  end subroutine add_vert_coord

  subroutine write_hist_coord_att(File, mdimind, boundsdim, dimonly, mdimid)
    use pio, only: file_desc_t, var_desc_t, pio_put_att, pio_noerr,           &
                   pio_int, pio_double, pio_inq_varid, pio_inq_dimid,         &
                   pio_def_dim, pio_def_var
    implicit none

    ! Input variables
    type(file_desc_t), intent(inout) :: File           ! PIO file Handle
    integer,           intent(in)    :: mdimind        ! Internal dim index
    integer,           intent(in)    :: boundsdim      ! Bounds dimension ID
    logical,           intent(in)    :: dimonly        ! No def_var if .true.
    integer, optional, intent(out)   :: mdimid

    ! Local variables
    integer                          :: dimid          ! PIO dimension ID
    type(var_desc_t)                 :: vardesc        ! PIO variable descriptor
    character(len=120)               :: errormsg
    character(len=max_chars)         :: formula_terms  ! Constructed string
    integer                          :: ierr

    ! Check to see if the dimension already exists in the file
    ierr = pio_inq_dimid(File, trim(hist_coords(mdimind)%name), dimid)
    if (ierr == PIO_NOERR) then
      write(errormsg, *) 'WARNING: A dimension already exists for ',          &
           trim(hist_coords(mdimind)%name)
      call endrun(errormsg)
    else
      ! Define the dimension
      ierr = pio_def_dim(File, trim(hist_coords(mdimind)%name),               &
           hist_coords(mdimind)%dimsize, dimid)
      call handle_pio_error(ierr, 'Unable to define dimension in write_hist_coord_att')
    end if
    ! If the caller wants to know the NetCDF dimension ID, set it here
    if (present(mdimid)) then
      mdimid = dimid
    end if
    if (.not. dimonly) then
      ! Check to see if the variable already exists in the file
      ierr = pio_inq_varid(File, hist_coords(mdimind)%name, vardesc)
      if (ierr == PIO_NOERR) then
        write(errormsg, *) 'WARNING: A variable already exists for ',         &
             trim(hist_coords(mdimind)%name)
        call endrun(errormsg)
      end if
      ! Time to define the variable
      if (hist_coords(mdimind)%integer_dim) then
        ierr = pio_def_var(File, trim(hist_coords(mdimind)%name), pio_int,    &
             (/dimid/), vardesc)
      else
        ierr = pio_def_var(File, trim(hist_coords(mdimind)%name), pio_double, &
             (/dimid/), vardesc)
      end if
      call handle_pio_error(ierr, 'Unable to define variable in write_hist_coord_att')
      ! long_name
      ierr=pio_put_att(File, vardesc, 'long_name', trim(hist_coords(mdimind)%long_name))
      call handle_pio_error(ierr, 'Error writing "long_name" attr in write_hist_coord_att')
      ! units
      if(len_trim(hist_coords(mdimind)%units) > 0) then
        ierr=pio_put_att(File, vardesc, 'units',                              &
             trim(hist_coords(mdimind)%units))
        call handle_pio_error(ierr, 'Error writing "units" attr in write_hist_coord_att')
      end if
      ! positive
      if(len_trim(hist_coords(mdimind)%positive) > 0) then
        ierr=pio_put_att(File, vardesc, 'positive',                           &
             trim(hist_coords(mdimind)%positive))
        call handle_pio_error(ierr, 'Error writing "positive" attr in write_hist_coord_att')
      end if
      ! standard_name
      if(len_trim(hist_coords(mdimind)%standard_name) > 0) then
        ierr=pio_put_att(File, vardesc, 'standard_name',                      &
             trim(hist_coords(mdimind)%standard_name))
        call handle_pio_error(ierr, 'Error writing "standard_name" attr in write_hist_coord_att')
      end if
      ! formula_terms
      if(len_trim(hist_coords(mdimind)%formula_terms%a_name) > 0) then
        write(formula_terms, '("a: ",a," b: ",a," p0: ",a," ps: ",a)')        &
             trim(hist_coords(mdimind)%formula_terms%a_name),                 &
             trim(hist_coords(mdimind)%formula_terms%b_name),                 &
             trim(hist_coords(mdimind)%formula_terms%p0_name),                &
             trim(hist_coords(mdimind)%formula_terms%ps_name)
        ierr=pio_put_att(File, vardesc, 'formula_terms', trim(formula_terms))
        call handle_pio_error(ierr, 'Error writing "formula_terms" attr in write_hist_coord_att')
      end if
      ! bounds
      if (associated(hist_coords(mdimind)%bounds)) then
        ! Write name of the bounds variable
        ierr=pio_put_att(File, vardesc, 'bounds', trim(hist_coords(mdimind)%bounds_name))
        call handle_pio_error(ierr, 'Error writing "bounds" attr in write_hist_coord_att')
      end if

      ! Now, we need to define and populate the associated bounds variable
      ! NB: Reusing vardesc, no longer assocated with main variable
      if (associated(hist_coords(mdimind)%bounds)) then
        if (size(hist_coords(mdimind)%bounds,2) /= hist_coords(mdimind)%dimsize) then
          ! If anyone hits this check, add a new dimension for this case
          write(errormsg, *) 'The bounds variable, ',                         &
               trim(hist_coords(mdimind)%bounds_name),                        &
               ', needs to have dimension (2,', hist_coords(mdimind)%dimsize
          call endrun(errormsg)
        end if
        ierr = pio_def_var(File, trim(hist_coords(mdimind)%bounds_name),      &
             pio_double, (/boundsdim,dimid/), vardesc)
        call handle_pio_error(ierr, 'Unable to define bounds variable in write_hist_coord_att')
      end if

      ! See if we have formula_terms variables to define
      ! Define the "a" variable name
      ! NB: Reusing vardesc, no longer assocated with previous variables
      if (associated(hist_coords(mdimind)%formula_terms%a_values)) then
        if (size(hist_coords(mdimind)%formula_terms%a_values) /= hist_coords(mdimind)%dimsize) then
          write(errormsg, *) 'The forumla_terms variable, ',                  &
               trim(hist_coords(mdimind)%formula_terms%a_name),               &
               ', needs to have dimension', hist_coords(mdimind)%dimsize
          call endrun(errormsg)
        end if
        ierr = pio_def_var(File, trim(hist_coords(mdimind)%formula_terms%a_name), &
             pio_double, (/dimid/), vardesc)
        call handle_pio_error(ierr, 'Unable to define "a" formula_terms variable in write_hist_coord_att')
        ierr = pio_put_att(File, vardesc, 'long_name', trim(hist_coords(mdimind)%formula_terms%a_long_name))
        call handle_pio_error(ierr, 'Error writing "long_name" attr for "a" formula_term in write_hist_coord_att')
      end if
      ! Define the "b" variable name
      ! NB: Reusing vardesc, no longer assocated with previous variables
      if (associated(hist_coords(mdimind)%formula_terms%b_values)) then
        if (size(hist_coords(mdimind)%formula_terms%b_values) /= hist_coords(mdimind)%dimsize) then
          write(errormsg, *) 'The forumla_terms variable, ',                  &
               trim(hist_coords(mdimind)%formula_terms%b_name),               &
               ', needs to have dimension', hist_coords(mdimind)%dimsize
          call endrun(errormsg)
        end if
        ierr = pio_def_var(File, trim(hist_coords(mdimind)%formula_terms%b_name), &
             pio_double, (/dimid/), vardesc)
        call handle_pio_error(ierr, 'Unable to define "b" formula_terms variable in write_hist_coord_att')
        ierr = pio_put_att(File, vardesc, 'long_name', trim(hist_coords(mdimind)%formula_terms%b_long_name))
        call handle_pio_error(ierr, 'Error writing "long_name" attr for "b" formula_term in write_hist_coord_att')
      end if
      ! Maybe define the p0 variable (this may be defined already which is OK)
      ! NB: Reusing vardesc, no longer assocated with previous variables
      if (hist_coords(mdimind)%formula_terms%p0_value /= fillvalue) then
        ierr = pio_inq_varid(File, trim(hist_coords(mdimind)%formula_terms%p0_name), vardesc)
        if (ierr /= PIO_NOERR) then
          ierr = pio_def_var(File, trim(hist_coords(mdimind)%formula_terms%p0_name), &
               pio_double, vardesc)
          call handle_pio_error(ierr, 'Unable to define "p0" formula_terms variable in write_hist_coord_att')
          ierr = pio_put_att(File, vardesc, 'long_name', trim(hist_coords(mdimind)%formula_terms%p0_long_name))
          call handle_pio_error(ierr, 'Error writing "long_name" attr for "p0" formula_term in write_hist_coord_att')
          ierr = pio_put_att(File, vardesc, 'units', trim(hist_coords(mdimind)%formula_terms%p0_units))
          call handle_pio_error(ierr, 'Error writing "units" attr for "p0" formula_term in write_hist_coord_att')
        end if
      end if
      ! PS is not our responsibility
    end if ! (.not. dimonly)

  end subroutine write_hist_coord_att

  subroutine write_hist_coord_attrs(File, boundsdim, mdimids, writemdims_in)
    use pio, only: file_desc_t, var_desc_t, pio_put_att, pio_def_var,         &
                   pio_bcast_error, pio_internal_error, pio_seterrorhandling, &
                   pio_char, pio_def_dim
    implicit none

    ! Input variables
    type(file_desc_t), intent(inout) :: File           ! PIO file Handle
    integer,           intent(in)    :: boundsdim      ! Bounds dimension ID
    integer, optional, allocatable, intent(out)   :: mdimids(:) ! NetCDF dim IDs
    logical, optional, intent(in)    :: writemdims_in  ! Write mdim variable

    ! Local variables
    integer                          :: i
    integer                          :: ierr
    integer                          :: dimids(2)      ! PIO dimension IDs
    logical                          :: writemdims     ! Define an mdim variable
    type(var_desc_t)                 :: vardesc        ! PIO variable descriptor

    if (present(mdimids)) then
      allocate(mdimids(registeredmdims))
    end if

    ! We will handle errors for this routine
    call pio_seterrorhandling(File, PIO_BCAST_ERROR)

    if (present(writemdims_in)) then
      writemdims = writemdims_in
    else
      writemdims = .false.
    end if

    ! NB: Currently, writemdims is for restart and we don't need to write
    ! these out in a history-restart file. This could change in the future.
    ! which would require a change to the function of the fourth argument
    ! Fill in the attribute information for each mdim
    do i = 1, registeredmdims
      if (present(mdimids)) then
        call write_hist_coord_att(File, i, boundsdim, writemdims, mdimids(i))
      else
        call write_hist_coord_att(File, i, boundsdim, writemdims)
      end if
    end do

    if (writemdims) then
      ierr = pio_def_dim(File, 'mdimslen', max_mdimname_len, dimids(1))
      call handle_pio_error(ierr, 'Unable to define len dimension for mdimnames in write_hist_coord_attrs')
      ierr = pio_def_dim(File, 'num_mdims', registeredmdims, dimids(2))
      call handle_pio_error(ierr, 'Unable to define dimension for mdimnames in write_hist_coord_attrs')
      ierr = pio_def_var(File, mdim_var_name, pio_char, dimids, vardesc)
      call handle_pio_error(ierr, 'Unable to define mdimnames in write_hist_coord_attrs')
      ierr = pio_put_att(File, vardesc, 'long_name', 'mdim dimension names')
      call handle_pio_error(ierr, 'Error writing "long_name" attr for mdimnames in write_hist_coord_attrs')
    end if

    ! Back to I/O or die trying
    call pio_seterrorhandling(File, PIO_INTERNAL_ERROR)

  end subroutine write_hist_coord_attrs

  subroutine write_hist_coord_var(File, mdimind)
    use pio, only: file_desc_t, var_desc_t, pio_put_var, pio_inq_varid, PIO_NOERR
    implicit none

    ! Input variables
    type(file_desc_t), intent(inout) :: File           ! PIO file Handle
    integer,           intent(in)    :: mdimind        ! Internal dim index

    ! Local variables
    type(var_desc_t)                 :: vardesc        ! PIO variable descriptor
    character(len=120)               :: errormsg
    integer                          :: ierr

    ! Check to make sure the variable already exists in the file
    ierr = pio_inq_varid(File, trim(hist_coords(mdimind)%name), vardesc)
    call handle_pio_error(ierr, 'Error writing values for nonexistent dimension variable write_hist_coord_var')
    ! Write out the values for this dimension variable
    if (hist_coords(mdimind)%integer_dim) then
      ierr = pio_put_var(File, vardesc, hist_coords(mdimind)%integer_values)
    else
      ierr = pio_put_var(File, vardesc, hist_coords(mdimind)%real_values)
    end if
    call handle_pio_error(ierr, 'Error writing variable values in write_hist_coord_var')

    ! Now, we need to possibly write values for the associated bounds variable
    if (associated(hist_coords(mdimind)%bounds)) then
      ! Check to make sure the variable already exists in the file
      ! NB: Reusing vardesc, no longer assocated with previous variables
      ierr = pio_inq_varid(File, trim(hist_coords(mdimind)%bounds_name), vardesc)
      call handle_pio_error(ierr, 'Error writing values for nonexistent bounds variable write_hist_coord_var')
    ! Write out the values for this bounds variable
      ierr = pio_put_var(File, vardesc, hist_coords(mdimind)%bounds)
      call handle_pio_error(ierr, 'Error writing bounds values in write_hist_coord_var')
    end if

    ! Write values for the "a" variable name
    if (associated(hist_coords(mdimind)%formula_terms%a_values)) then
      ! Check to make sure the variable already exists in the file
      ! NB: Reusing vardesc, no longer assocated with previous variables
      ierr = pio_inq_varid(File, trim(hist_coords(mdimind)%formula_terms%a_name), vardesc)
      call handle_pio_error(ierr, 'Error writing values for nonexistent "a" formula_terms variable write_hist_coord_var')
    ! Write out the values for this "a" formula_terms variable
      ierr = pio_put_var(File, vardesc, hist_coords(mdimind)%formula_terms%a_values)
      call handle_pio_error(ierr, 'Error writing "a" formula_terms values in write_hist_coord_var')
    end if
    ! Write values for the "b" variable name
    if (associated(hist_coords(mdimind)%formula_terms%b_values)) then
      ! Check to make sure the variable already exists in the file
      ! NB: Reusing vardesc, no longer assocated with previous variables
      ierr = pio_inq_varid(File, trim(hist_coords(mdimind)%formula_terms%b_name), vardesc)
      call handle_pio_error(ierr, 'Error writing values for nonexistent "b" formula_terms variable write_hist_coord_var')
    ! Write out the values for this "b" formula_terms variable
      ierr = pio_put_var(File, vardesc, hist_coords(mdimind)%formula_terms%b_values)
      call handle_pio_error(ierr, 'Error writing "b" formula_terms values in write_hist_coord_var')
    end if
    ! Write values for the "p0" variable name (this may be an overwrite, too bad)
    if (hist_coords(mdimind)%formula_terms%p0_value /= fillvalue) then
      ! Check to make sure the variable already exists in the file
      ! NB: Reusing vardesc, no longer assocated with previous variables
      ierr = pio_inq_varid(File, trim(hist_coords(mdimind)%formula_terms%p0_name), vardesc)
      call handle_pio_error(ierr, 'Error writing values for nonexistent "p0" formula_terms variable write_hist_coord_var')
    ! Write out the values for this "p0" formula_terms variable
      ierr = pio_put_var(File, vardesc, hist_coords(mdimind)%formula_terms%p0_value)
      call handle_pio_error(ierr, 'Error writing "p0" formula_terms values in write_hist_coord_var')
    end if

  end subroutine write_hist_coord_var

  subroutine write_hist_coord_vars(File, writemdims_in)
   use pio, only: file_desc_t, var_desc_t, pio_put_var,                       &
                  pio_bcast_error, pio_internal_error,                        &
                  pio_seterrorhandling, pio_inq_varid
    implicit none

    ! Input variables
    type(file_desc_t), intent(inout) :: File           ! PIO file Handle
    logical, optional, intent(in)    :: writemdims_in  ! Write mdim variable

    ! Local variables
    integer                          :: i
    integer                          :: ierr
    logical                          :: writemdims     ! Define an mdim variable
    type(var_desc_t)                 :: vardesc        ! PIO variable descriptor
    character(len=max_mdimname_len), allocatable :: mdimnames(:)

    ! We will handle errors for this routine
    call pio_seterrorhandling(File, PIO_BCAST_ERROR)

    if (present(writemdims_in)) then
      writemdims = writemdims_in
    else
      writemdims = .false.
    end if

    if (writemdims) then
      allocate(mdimnames(registeredmdims))
    end if

    ! Write out the variable values for each mdim
    do i = 1, registeredmdims
      if (.not. writemdims) then
        ! NB: Currently, writemdims is for restart and we don't need to write
        ! these out in a history-restart file. This could change in the future
        ! which is why it is a separate if block
        ! Fill in the attribute information for each mdim
        call write_hist_coord_var(File, i)
      end if
      if (writemdims) then
        mdimnames(i) = trim(hist_coords(i)%name)
      end if
    end do

    if (writemdims) then
      ierr = pio_inq_varid(File, mdim_var_name, vardesc)
      call handle_pio_error(ierr, 'Error writing values for nonexistent mdimnames variable in write_hist_coord_vars')
      ! Write out the values for mdim names
      ierr = pio_put_var(File, vardesc, mdimnames)
      call handle_pio_error(ierr, 'Error writing values for mdimnames variable in write_hist_coord_vars')
      deallocate(mdimnames)
    end if

    ! Back to I/O or die trying
    call pio_seterrorhandling(File, PIO_INTERNAL_ERROR)

  end subroutine write_hist_coord_vars

  subroutine lookup_hist_coord_indices(mdimnames, mdimindicies)
    character(len=*), intent(in) :: mdimnames(:)
    integer, intent(out) :: mdimindicies(:)

    integer :: i, j
    integer :: cnt
    character(len=120) :: errormsg
    character(len=16) :: name


    cnt = size(mdimnames)
    mdimindicies = -1


    do j=1,cnt
      name = mdimnames(j)
      do i=1,registeredmdims
        if(name .eq. hist_coords(i)%name) then
          mdimindicies(j)=i
        end if
      end do
    end do
    do j=1,cnt
      if(mdimindicies(j)<0) then
        do i=1,registeredmdims		
          print *,__FILE__,__LINE__,i,hist_coords(i)%name
        end do
        write(errormsg,*) 'Name ',mdimnames(j),' not in registered mdimnames'
        call endrun(errormsg)
      end if
    end do


  end subroutine lookup_hist_coord_indices

  !#######################################################################

  character(len=8) function sec2hms (seconds)

    ! Input arguments

    integer, intent(in) :: seconds

    ! Local workspace

    integer :: hours     ! hours of hh:mm:ss
    integer :: minutes   ! minutes of hh:mm:ss
    integer :: secs      ! seconds of hh:mm:ss

    if (seconds < 0 .or. seconds > 86400) then
      write(iulog,*)'SEC2HRS: bad input seconds:', seconds
      call endrun ()
    end if

    hours   = seconds / 3600
    minutes = (seconds - hours*3600) / 60
    secs    = (seconds - hours*3600 - minutes*60)

    if (minutes < 0 .or. minutes > 60) then
      write(iulog,*)'SEC2HRS: bad minutes = ',minutes
      call endrun ()
    end if

    if (secs < 0 .or. secs > 60) then
      write(iulog,*)'SEC2HRS: bad secs = ',secs
      call endrun ()
    end if

    write(sec2hms,80) hours, minutes, secs
80  format(i2.2,':',i2.2,':',i2.2)
    return
  end function sec2hms
  character(len=10) function date2yyyymmdd (date)

    ! Input arguments

    integer, intent(in) :: date

    ! Local workspace

    integer :: year    ! year of yyyy-mm-dd
    integer :: month   ! month of yyyy-mm-dd
    integer :: day     ! day of yyyy-mm-dd

    if (date < 0) then
      call endrun ('DATE2YYYYMMDD: negative date not allowed')
    end if

    year  = date / 10000
    month = (date - year*10000) / 100
    day   = date - year*10000 - month*100

    write(date2yyyymmdd,80) year, month, day
80  format(i4.4,'-',i2.2,'-',i2.2)
    return
  end function date2yyyymmdd

  !#######################################################################


end module cam_history_support
