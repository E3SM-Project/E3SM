module glb_verif_smry
!------------------------------------------------------------------------------------
! Description:

! This module provides utilities for getting statistical summaries for global fields.
! The module can be used to identify, for example, negative values in tracer 
! mixing ratios, or energy conservation errors exceeding a certain threshold.
!
! The basic idea is to build a list of fields for which such summaries will be provided 
! during  model integration.  A field is identified by its unique name that contains 
! information about both the physical quantity and piece of code in which the values 
! are monitored. 
!
! For each field, first identify violations within each 
! chunk of the physics grid; get a total count, and note down the extreme value and 
! its location (chunk/column index, lat, lon, and vertical level index if applicable).
! Then the total count and the extreme among all chunks on a single MPI process ("domain")
! are obtained. Lastly, the domain summaries are collected  by the master process to 
! provide a global summary.
! 
!
! History:
!
! First  version by Hui Wan (PNNL, 2017-05). 
!  - MPI_gather was used for global communication.
!
! Second version by Hui Wan (PNNL, 2017-07/08), 
! with inputs from Wuyin Lin (BNL)
!  - Replace MPI_gather by MPI_allreduce;
!  - Reduce computational cost related to string comparison;
!  - Distinguish different verbose levels;
!  - Add option to provide summary output infrequently.  
!------------------------------------------------------------------------------------

  use shr_kind_mod,   only: r8=>SHR_KIND_R8
  use shr_kind_mod,   only: longchar=>SHR_KIND_CL
  use cam_abortutils, only: endrun
  use cam_logfile,    only: iulog
  use physconst,      only: pi

  implicit none
  private

  public tp_stat_smry
  public add_smry_field
  public global_smry_init
  public timestep_smry_init
  public get_smry_field_idx
  public get_chunk_smry
  public get_global_smry

  !-------------------------------------------------------------------
  ! Namelist switches for users (see also phys_control.F90) 
  !-------------------------------------------------------------------
  ! Frequency of smry report: 
  !  Negative: unit is hours.
  !  Positive: unit is time steps.

  integer,public :: glb_verif_smry_frq =  -6   

  !-------------------------------------------------------------------
  ! Level of detail included in the summary messages.
  ! <0: no smry.
  !  0: provide a one-line summary for each registered field if 
  !     there is any value exceeding the corresponding threshold; 
  !     report on # of violations and the extreme values.
  !  1: in addition to 0, also report on the locations of extreme values.
  !  2: in addition to 1, print out summary for every chunk of the 
  !     physics grid. This is similar to the original implementation 
  !     in QNEG3 and QNEG4. 

  integer,public :: glb_verif_smry_level = 0

  !-------------------------------------------------------------------
  ! Print summary for all registered fields regardless of the total count of violations. 
  ! When set to .false., summary is printed only for fields with counts >0.

  logical,public :: l_print_smry_for_all_fields = .false.
                                
  !----------------------------------------------
  ! Types of comparison supported by this module

  integer,public,parameter :: SMALLER_THAN     = -1
  integer,public,parameter :: GREATER_EQ       =  1
  integer,public,parameter :: ABS_SMALLER_THAN = -2
  integer,public,parameter :: ABS_GREATER_EQ   =  2

  !----------------------------------------------
  ! Types of fixers supported by this module

  integer,public,parameter :: NO_FIX   = 0
  integer,public,parameter :: CLIPPING = 1

  !-------------------------------------------------------------------
  ! Define our own "short string" length. SHR_KIND_CS is very long (80)
  integer,private,parameter :: shortchar = 36

  ! Name of this module that will appear in error messages
  character(len=shortchar),private,parameter :: THIS_MODULE = 'glb_verif_smry'

  !----------------
  ! Constants

  integer, parameter :: INT_UNDEF = -999 
  real(r8),parameter :: FLT_UNDEF = -999._r8

  real(r8),parameter :: rad2deg = 180._r8/pi

  !----------------
  ! Data structure 

  type tp_stat_smry

    character(len=shortchar) :: field_name       ! the physical quantity to be evaluated
    character(len=shortchar) :: field_unit       ! unit of the field

    integer  :: cmpr_type  ! one of the comparison types defined above
    real(r8) :: threshold  ! threshold specified by developer/user
    integer  :: count = 0  ! total number of cells with values exceeding threshold
    integer  :: fixer = NO_FIX  ! by default, do not fix values exceeding threshold

    ! extreme value and its location

    real(r8) :: extreme_val  = FLT_UNDEF
    real(r8) :: extreme_lat  = FLT_UNDEF/rad2deg
    real(r8) :: extreme_lon  = FLT_UNDEF/rad2deg
    integer  :: extreme_chnk = INT_UNDEF
    integer  :: extreme_col  = INT_UNDEF
    integer  :: extreme_lev  = INT_UNDEF

  end type tp_stat_smry

  !-------------------------------
  ! Misc module variables

  integer,parameter       :: max_number_of_smry_fields = 1000
  integer,public          :: current_number_of_smry_fields = 0
  integer                 :: n_smry_fields_mpimax = 0   ! number of fields needing mpimax for mpi_allreduce
  integer                 :: n_smry_fields_mpimin = 0   ! number of fields needing mpimin for mpi_allreduce

  character(len=longchar) :: msg 
  logical                 :: l_smry_arrays_allocated = .false.

  logical                 :: timestep_smry_on = .true.  ! get smry for the current time step?
                                                        ! re-evaluated during phys_timestep_init.

  !-------------------------------------------------------------------
  ! List of registered fields 

  type(tp_stat_smry) :: global_smry_1d(max_number_of_smry_fields)

  !-------------------------------
  interface get_chunk_smry
    module procedure get_chunk_smry_1_lev_real   ! for fields that do not have a vertical distribution
    module procedure get_chunk_smry_m_lev_real   ! for fields with multiple vertical levels
  end interface get_chunk_smry

contains

  !--------------------------------------------------------------------------------------------
  ! Description: 
  !  The subroutine registers a new field for getting global summary. It is expected to be
  !  called during the initialization of various parameterizations.
  !--------------------------------------------------------------------------------------------
  subroutine add_smry_field( fldname, fldunit, cmprtype, threshold, fixer, fldidx )

    character(len=*), intent(in)   :: fldname
    character(len=*), intent(in)   :: fldunit
    integer                        :: cmprtype
    real(r8)                       :: threshold
    integer, intent(out), optional :: fldidx
    integer, intent(in),  optional :: fixer

    integer :: ii

    ! During the model initialization, we first register all the fields that we would
    ! like to get global summary for, then allocate memory for the chunk/domain summary
    ! arrays. No new fields can be added after that allocation.

    if (l_smry_arrays_allocated) then
        msg = trim(THIS_MODULE)//': subroutine add_smry_field should not be called'//&
              ' after global_smry_init has been called.'
        call endrun(trim(msg))
    end if

    ! Check if the field has already been registered.

    do ii = 1,current_number_of_smry_fields
       if (trim(global_smry_1d(ii)%field_name) == trim(fldname) ) then
           msg = trim(THIS_MODULE)//': field '//trim(fldname)// &
                 ' has already been added to list of statistical summary.'
           call endrun(trim(msg))
       end if
    end do

    ! Now add a new field

    current_number_of_smry_fields = current_number_of_smry_fields + 1

    if (current_number_of_smry_fields.gt.max_number_of_smry_fields) then
       msg = trim(THIS_MODULE)//': max. No. of fields exceeded when attempting to add '//trim(fldname)
       call endrun(trim(msg))
    end if

    ii = current_number_of_smry_fields
    global_smry_1d(ii)%field_name     = trim(fldname)
    global_smry_1d(ii)%field_unit     = trim(fldunit)
    global_smry_1d(ii)%threshold      = threshold
    global_smry_1d(ii)%cmpr_type      = cmprtype
    if (present(fixer)) &
    global_smry_1d(ii)%fixer          = fixer

    ! Increase the field count for mpimax/mpimin

    SELECT CASE (cmprtype)
    CASE (GREATER_EQ,ABS_GREATER_EQ)
       n_smry_fields_mpimax = n_smry_fields_mpimax + 1

    CASE (SMALLER_THAN,ABS_SMALLER_THAN)
       n_smry_fields_mpimin = n_smry_fields_mpimin + 1

    CASE DEFAULT
       msg = trim(THIS_MODULE)//': unknown comparison for field '//trim(fldname) 
       call endrun(trim(msg))
    END SELECT

    if (present(fldidx)) fldidx = ii 
 
  end subroutine add_smry_field

  !--------------------------------------------------------------------------------
  ! Description: 
  !   Allocate memory for the chunk/domain summary variables, and copy meta data 
  !   from global_smry_1d. This subroutine needs to be called during model
  !   after all "add_smry_field" calls. 
  !--------------------------------------------------------------------------------
  subroutine global_smry_init( chunk_smry_2d, domain_smry_1d, begchunk, endchunk )

    use spmd_utils,   only: masterproc

    type(tp_stat_smry), pointer ::  chunk_smry_2d(:,:)
    type(tp_stat_smry), pointer :: domain_smry_1d(:)
    integer, intent(in) :: begchunk, endchunk

    integer :: ierr, ichnk, ifld

    ! Sanity check

    if (l_smry_arrays_allocated) then
        msg = trim(THIS_MODULE)//': attempting to call global_smry_init multiple times.'
        call endrun(trim(msg))
    end if

    !---------------------------------------
    ! Initialize array for domain summaries 
    !---------------------------------------
    ! Allocate memory

    allocate(domain_smry_1d(current_number_of_smry_fields), stat=ierr)
    if( ierr /= 0 ) then
       write(msg,*) trim(THIS_MODULE)//': domain_smry allocation error = ',ierr
       call endrun(trim(msg))
    end if

    ! Copy metadata

    domain_smry_1d(1:current_number_of_smry_fields)%field_name = &
    global_smry_1d(1:current_number_of_smry_fields)%field_name

    domain_smry_1d(1:current_number_of_smry_fields)%field_unit = &
    global_smry_1d(1:current_number_of_smry_fields)%field_unit

    domain_smry_1d(1:current_number_of_smry_fields)%cmpr_type = &
    global_smry_1d(1:current_number_of_smry_fields)%cmpr_type

    domain_smry_1d(1:current_number_of_smry_fields)%threshold = &
    global_smry_1d(1:current_number_of_smry_fields)%threshold

    domain_smry_1d(1:current_number_of_smry_fields)%fixer     = &
    global_smry_1d(1:current_number_of_smry_fields)%fixer
    !--------------------------------------
    ! Initialize array for chunk summaries 
    !--------------------------------------
    ! Allocate memory

    allocate(chunk_smry_2d(begchunk:endchunk,current_number_of_smry_fields), stat=ierr)
    if( ierr /= 0 ) then
       write(msg,*) trim(THIS_MODULE)//': chunk_smry allocation error = ',ierr
       call endrun(trim(msg))
    end if

    ! Copy metadata

    do ichnk = begchunk,endchunk 

       chunk_smry_2d(ichnk,1:current_number_of_smry_fields)%field_name = &
      global_smry_1d(      1:current_number_of_smry_fields)%field_name

       chunk_smry_2d(ichnk,1:current_number_of_smry_fields)%field_unit = &
      global_smry_1d(      1:current_number_of_smry_fields)%field_unit

       chunk_smry_2d(ichnk,1:current_number_of_smry_fields)%cmpr_type = &
      global_smry_1d(      1:current_number_of_smry_fields)%cmpr_type

       chunk_smry_2d(ichnk,1:current_number_of_smry_fields)%threshold = &
      global_smry_1d(      1:current_number_of_smry_fields)%threshold

       chunk_smry_2d(ichnk,1:current_number_of_smry_fields)%fixer     = &
      global_smry_1d(      1:current_number_of_smry_fields)%fixer

       chunk_smry_2d(ichnk,1:current_number_of_smry_fields)%extreme_chnk = ichnk
    end do

    ! Set flag for sanity check later

    l_smry_arrays_allocated = .true.

    ! Send message to log file

    if (masterproc) then
       write(iulog,*)'****************************************************************************'
       write(iulog,*)'GLB_VERIF_SMRY: Initialization of glb_verif_smry done.'
       write(iulog,*)'GLB_VERIF_SMRY: current_number_of_smry_fields = ',current_number_of_smry_fields 
       write(iulog,*)'GLB_VERIF_SMRY: '

       write(iulog,'(a19,a8, a40,a12, a10,a11,a7)')  &
                     'GLB_VERIF_SMRY:','Index', 'Field','Unit', 'Compr.','Threshold','Fixer'

       do ifld = 1,current_number_of_smry_fields

          write(iulog,'(a19,i8, a40,a12, i10,e11.3,i7)')                &
                      'GLB_VERIF_SMRY:',ifld,                               &
                      trim(domain_smry_1d(ifld)%field_name),                &
                      trim(domain_smry_1d(ifld)%field_unit),                &
                      domain_smry_1d(ifld)%cmpr_type,                       &
                      domain_smry_1d(ifld)%threshold, &
                      domain_smry_1d(ifld)%fixer 
       end do

       write(iulog,*)'****************************************************************************'
    end if

  end subroutine global_smry_init

  !--------------------------------------------------------------------------------
  ! Description:
  ! Decide whether chunk/global summary should be collected at this time step
  !--------------------------------------------------------------------------------
  subroutine timestep_smry_init( nstep )

    integer, intent(in) :: nstep

    timestep_smry_on = (glb_verif_smry_level .ge. 0) .and. &
                        mod(nstep,glb_verif_smry_frq) == 0

  end subroutine timestep_smry_init

  !---------------------------------------------------------------------
  ! Description:
  !   Find a field on the global summary list and return the index.
  !---------------------------------------------------------------------
  subroutine get_smry_field_idx( fldname, fldidx )

    character(len=*), intent(in)   :: fldname
    integer, intent(out)           :: fldidx

    integer :: ii

    fldidx = INT_UNDEF
    do ii = 1,current_number_of_smry_fields 
       if (global_smry_1d(ii)%field_name == fldname ) then
          fldidx = ii
          exit
       end if
    end do

  end subroutine get_smry_field_idx

  !---------------------------------------------------------------------------------------
  ! Description:
  !   Identify values in a field that exceed a threshold. The total number of such values 
  !   and the location of the extreme are noted down for later use.
  !   This subroutine is meant to operate on all columns (or a subset of columns) 
  !   in a single chunk of CAM's physics grid. The particular incarnation of the 
  !   subroutine deals with fields that have multiple vertical levels. 
  !---------------------------------------------------------------------------------------
  subroutine get_chunk_smry_m_lev_real( fldname, ncol, nlev, array_in, &
                                      &  lat, lon, chunk_smry, ifld_out )

    character(len=*),  intent(in)    :: fldname
    integer,           intent(in)    :: ncol                 ! number of columns packed in array
    integer,           intent(in)    :: nlev                 ! number of vertical levels
    real(r8),          intent(inout) :: array_in(:,:)        ! input array of values to be checked
    real(r8),          intent(in)    :: lat(:)
    real(r8),          intent(in)    :: lon(:)
    type(tp_stat_smry),intent(inout) :: chunk_smry(:)
    integer, optional, intent(out)   :: ifld_out

    ! Local variables

    integer  :: ifld
    real(r8) :: array(ncol,nlev)  ! equals array_in or abs(array_in) 
    integer  :: iflag(ncol,nlev) 
    integer  :: idx(2)
    character(len=shortchar) :: cmpr_type_char
  
    !-------------------------------------------------------------------------
    if (.not.timestep_smry_on) return
 
    !--------------------------------
    ! Find field on the master list.

    call get_smry_field_idx( fldname, ifld )

    if (present(ifld_out)) ifld_out = ifld
    if (ifld.eq.INT_UNDEF) return

    !-------------------------------------------------------------------------
    ! Calculate the total number of grid cells with value exceeding threshold
    ! and identify location of the extremem value.
  
    iflag(:,:) = 0

    SELECT CASE (chunk_smry(ifld)%cmpr_type)
    CASE (GREATER_EQ)
      cmpr_type_char = '>='
      array = array_in
      where( array .ge. chunk_smry(ifld)%threshold ) iflag = 1
      idx = maxloc( array )

    CASE (SMALLER_THAN)
      cmpr_type_char = '<'
      array = array_in
      where( array .lt. chunk_smry(ifld)%threshold ) iflag = 1
      idx = minloc( array )

    CASE (ABS_GREATER_EQ)
      cmpr_type_char = 'ABS >='
      array = abs(array_in)
      WHERE( array .ge. chunk_smry(ifld)%threshold ) iflag = 1
      idx = maxloc( array )

    CASE (ABS_SMALLER_THAN)
      cmpr_type_char = 'ABS < '
      array = abs(array_in)
      WHERE( array .lt. chunk_smry(ifld)%threshold ) iflag = 1
      idx = minloc( array )

    END SELECT

    ! Total number of values exceeding threshold

    chunk_smry(ifld)%count = sum( iflag )

    ! The extreme value

    chunk_smry(ifld)%extreme_val  = array(idx(1),idx(2))
    chunk_smry(ifld)%extreme_col  =       idx(1)
    chunk_smry(ifld)%extreme_lev  =       idx(2)
    chunk_smry(ifld)%extreme_lat  =   lat(idx(1))
    chunk_smry(ifld)%extreme_lon  =   lon(idx(1))

    ! Clipping

    if (chunk_smry(ifld)%fixer.eq.CLIPPING) then
       where( iflag.eq.1)  array_in = chunk_smry(ifld)%threshold
    end if 
  
    ! Send message to log file
  
    if ( (glb_verif_smry_level.ge.2) .and. &
         (chunk_smry(ifld)%count.gt.0 .or. l_print_smry_for_all_fields) ) then

       write(iulog,'(2x,a,a40,a12,a2, i8,a,a7,e15.7, a,e15.7, a3,2(a,f7.2),(a,i4),(a,i10),(a,i4),(a,i2))') &
         'chunk_smry: ', &
         trim(chunk_smry(ifld)%field_name),'('//trim(chunk_smry(ifld)%field_unit)//')',': ', &
         chunk_smry(ifld)%count, ' values ',trim(cmpr_type_char), chunk_smry(ifld)%threshold, &
         ', extreme is ',chunk_smry(ifld)%extreme_val,' at ',&
         '  lat ',chunk_smry(ifld)%extreme_lat *rad2deg, &
         ', lon ',chunk_smry(ifld)%extreme_lon *rad2deg, &
         ', lev ',chunk_smry(ifld)%extreme_lev, &
         ', chnk ',chunk_smry(ifld)%extreme_chnk, &
         ', col ',chunk_smry(ifld)%extreme_col, &
         ', fixer = ',chunk_smry(ifld)%fixer  
    end if
  
  end subroutine get_chunk_smry_m_lev_real

  !---------------------------------------------------------------------------------------
  ! Description:
  !   Identify values in a field that exceed a threshold. The total number of such values 
  !   and the location of the extreme are noted down for later use.
  !   This subroutine is meant to operate on all columns (or a subset of columns) 
  !   in a single chunk of CAM's physics grid. The particular incarnation of the 
  !   subroutine deals with fields that do not have a vertical distribution (e.g. surface
  !   fluxes and vertical integrals). 
  !---------------------------------------------------------------------------------------
  subroutine get_chunk_smry_1_lev_real( fldname, ncol, array_in, &
                                        lat, lon, chunk_smry, ifld_out )

    character(len=*),  intent(in)    :: fldname
    integer,           intent(in)    :: ncol           ! number of columns packed in array
    real(r8),          intent(inout) :: array_in(:)    ! input array of values to be checked
    real(r8),          intent(in)    :: lat(:)
    real(r8),          intent(in)    :: lon(:)
    type(tp_stat_smry),intent(inout) :: chunk_smry(:)
    integer, optional, intent(out)   :: ifld_out

    ! Local variables

    integer  :: ifld
    real(r8) :: array(ncol)    ! equals array_in or abs(array_in) 
    integer  :: iflag(ncol) 
    integer  :: idx(1)
    character(len=shortchar) :: cmpr_type_char

    !-------------------------------------------------------------------------
    if (.not.timestep_smry_on) return

    !--------------------------------
    ! Find field on the master list

    call get_smry_field_idx( fldname, ifld )

    if (present(ifld_out)) ifld_out = ifld
    if (ifld.eq.INT_UNDEF) return
 
    !-----------------------------------------------------------------------
    ! Calculate the total number of columns with value exceeding threshold
    ! and identify location of the extremem value.
  
    iflag(:) = 0

    SELECT CASE (chunk_smry(ifld)%cmpr_type)
    CASE (GREATER_EQ)
      cmpr_type_char = '>='
      array = array_in
      where( array .ge. chunk_smry(ifld)%threshold ) iflag = 1
      idx = maxloc( array )

    CASE (SMALLER_THAN)
      cmpr_type_char = '<'
      array = array_in
      where( array .lt. chunk_smry(ifld)%threshold ) iflag = 1
      idx = minloc( array )

    CASE (ABS_GREATER_EQ)
      cmpr_type_char = 'ABS >='
      array = abs(array_in)
      WHERE( array .ge. chunk_smry(ifld)%threshold ) iflag = 1
      idx = maxloc( array )

    CASE (ABS_SMALLER_THAN)
      cmpr_type_char = 'ABS <'
      array = abs(array_in)
      WHERE( array .lt. chunk_smry(ifld)%threshold ) iflag = 1
      idx = minloc( array )

    END SELECT

    ! Total number of values exceeding threshold

    chunk_smry(ifld)%count = sum( iflag )

    ! The extreme value

    chunk_smry(ifld)%extreme_val  = array(idx(1))
    chunk_smry(ifld)%extreme_col  =       idx(1)
    chunk_smry(ifld)%extreme_lat  =   lat(idx(1))
    chunk_smry(ifld)%extreme_lon  =   lon(idx(1))

    ! Clipping

    if (chunk_smry(ifld)%fixer.eq.CLIPPING) then
       where( iflag.eq.1)  array_in = chunk_smry(ifld)%threshold
    end if 

    ! Send message to log file
  
    if ( (glb_verif_smry_level.ge.2) .and. &
         (chunk_smry(ifld)%count.gt.0 .or. l_print_smry_for_all_fields) ) then

       write(iulog,'(2x,a,a40,a12,a2, i8,a,a7,e15.7, a,e15.7, a3,2(a,f7.2),(a,i10),(a,i4),(a,i2))') &
         'chunk_smry: ', &
         trim(chunk_smry(ifld)%field_name),'('//trim(chunk_smry(ifld)%field_unit)//')',': ', &
         chunk_smry(ifld)%count, ' values ',trim(cmpr_type_char), chunk_smry(ifld)%threshold, &
         ', extreme is ',chunk_smry(ifld)%extreme_val,' at ',&
         '  lat ',chunk_smry(ifld)%extreme_lat *rad2deg, &
         ', lon ',chunk_smry(ifld)%extreme_lon *rad2deg, &
         ', chnk ',chunk_smry(ifld)%extreme_chnk, &
         ', col ',chunk_smry(ifld)%extreme_col, &
         ', fixer = ',chunk_smry(ifld)%fixer  
    end if
  
  end subroutine get_chunk_smry_1_lev_real

  !---------------------------------------------------------------------------------------
  ! Description:
  !   Assuming the values exceeding a threshold have been identified for each chunk 
  !   handled by the current MPI process, this subroutine gets a total count of 
  !   violations in the current MPI process ("domain"), and identify the extreme value 
  !   among all chunks. Each call of this subroutine handles one field.
  !---------------------------------------------------------------------------------------
  subroutine get_domain_smry( chunk_smry_of_all_chunks, domain_smry )

    type(tp_stat_smry), intent(inout) :: chunk_smry_of_all_chunks(:)  ! shape: (nchnk)
    type(tp_stat_smry), intent(inout) :: domain_smry

    ! Local variables

    integer  :: idx(1)
    integer  :: ichnk
    character(len=shortchar) :: cmpr_type_char
  
    !------------------------------------------------
    ! Get a total count of values exceeding threshold

    domain_smry%count = sum( chunk_smry_of_all_chunks(:)%count )

    !-------------------
    ! Locate the extreme 
  
    SELECT CASE (domain_smry%cmpr_type)
    CASE (GREATER_EQ)
      idx = maxloc( chunk_smry_of_all_chunks(:)%extreme_val )
      cmpr_type_char = '>='

    CASE (ABS_GREATER_EQ)
      idx = maxloc( abs(chunk_smry_of_all_chunks(:)%extreme_val) )
      cmpr_type_char = 'ABS >='

    CASE (SMALLER_THAN)
      idx = minloc( chunk_smry_of_all_chunks(:)%extreme_val )
      cmpr_type_char = '<'

    CASE (ABS_SMALLER_THAN)
      idx = minloc( abs(chunk_smry_of_all_chunks(:)%extreme_val) )
      cmpr_type_char = 'ABS <'
    END SELECT

    ichnk = idx(1)
    domain_smry%extreme_val  = chunk_smry_of_all_chunks(ichnk)%extreme_val
    domain_smry%extreme_lat  = chunk_smry_of_all_chunks(ichnk)%extreme_lat
    domain_smry%extreme_lon  = chunk_smry_of_all_chunks(ichnk)%extreme_lon 
    domain_smry%extreme_lev  = chunk_smry_of_all_chunks(ichnk)%extreme_lev 
    domain_smry%extreme_col  = chunk_smry_of_all_chunks(ichnk)%extreme_col
    domain_smry%extreme_chnk = chunk_smry_of_all_chunks(ichnk)%extreme_chnk 

    ! Send message to log file
  
    if ( (glb_verif_smry_level.ge.2) .and. &
         (domain_smry%count.gt.0 .or. l_print_smry_for_all_fields) ) then

       write(iulog,'(2x,a,a40,a12,a2, i8,a,a7,e15.7, a,e15.7, a3,2(a,f7.2),(a,i4),(a,i10),(a,i4),(a,i2))') &
         'domain_smry: ', &
         trim(domain_smry%field_name),'('//trim(domain_smry%field_unit)//')',': ', &
         domain_smry%count, ' values ',trim(cmpr_type_char), domain_smry%threshold, &
         ', extreme is ',domain_smry%extreme_val,' at ',&
         '  lat ',domain_smry%extreme_lat *rad2deg, &
         ', lon ',domain_smry%extreme_lon *rad2deg, &
         ', lev ',domain_smry%extreme_lev, &
         ', chnk ',domain_smry%extreme_chnk, &
         ', col ',domain_smry%extreme_col, &
         ', fixer = ',domain_smry%fixer  
    end if
  
  end subroutine get_domain_smry

  !------------------------------------------------------------------------------------------
  ! Description:
  !   Assuming the values exceeding a threshold have been identified for each chunk, 
  !   this subroutine first calls get_domain_smry to get a summary for each MPI process,
  !   then collect information from all MPI processes and get the global summary.
  !   This subroutine is meant to handle all fields on the summary list.
  !------------------------------------------------------------------------------------------
  subroutine get_global_smry( chunk_smry_2d, domain_smry_1d, nstep)

#ifdef SPMD
    use mpishorthand, only: mpir8, mpiint, mpicom
    use spmd_utils,   only: npes
#endif
    use spmd_utils,   only: masterproc
    use cam_logfile,  only: iulog

    type(tp_stat_smry)  ::  chunk_smry_2d(:,:)  ! shape: (nchunk,nfld)
    type(tp_stat_smry)  :: domain_smry_1d(:)    ! shape: (nfld)
    integer,intent(in)  :: nstep                ! model time step

    integer :: ii

#ifdef SPMD

    real(r8) :: snd_array_mpimax (n_smry_fields_mpimax)
    real(r8) :: rcv_array_mpimax (n_smry_fields_mpimax)

    real(r8) :: snd_array_mpimin (n_smry_fields_mpimin)
    real(r8) :: rcv_array_mpimin (n_smry_fields_mpimin)

    integer :: imax, imin

#endif
    character(len=shortchar) :: cmpr_type_char

    if (.not.timestep_smry_on) return
    if (current_number_of_smry_fields.lt.1) return

    !--------------------------------------------
    ! Get domain summaries for each MPI process
    !--------------------------------------------
    do ii = 1,current_number_of_smry_fields
       call get_domain_smry( chunk_smry_2d(:,ii), domain_smry_1d(ii) )
    end do

    !--------------------------------------------
    ! Get global summaries
    !--------------------------------------------
#ifdef SPMD
    ! Pack arrays for MPI communication

    imax = 0
    imin = 0

    do ii = 1,current_number_of_smry_fields

      SELECT CASE (domain_smry_1d(ii)%cmpr_type)
      CASE (GREATER_EQ,ABS_GREATER_EQ)
        imax = imax +1
        snd_array_mpimax(imax) = domain_smry_1d(ii)%extreme_val
      
      CASE (SMALLER_THAN,ABS_SMALLER_THAN)
        imin = imin +1
        snd_array_mpimin(imin) = domain_smry_1d(ii)%extreme_val
      END SELECT
    end do

    !- MPI communications ---
    !
    !  Find the global max values

    if (imax.ne.n_smry_fields_mpimax) then
       call endrun(trim(THIS_MODULE)//'imax .ne. n_smry_fields_mpimax!')
    end if

    if (imax.gt.0) then
       call mpiallmaxreal( snd_array_mpimax, rcv_array_mpimax, imax, mpir8,  mpicom )
    end if

    !  Find the global min values

    if (imin.ne.n_smry_fields_mpimin) then
       call endrun(trim(THIS_MODULE)//'imin .ne. n_smry_fields_mpimin!')
    end if

    if (imin.gt.0) then
       call mpiallminreal( snd_array_mpimin, rcv_array_mpimin, imin, mpir8,  mpicom )
    end if

    ! Sum up the violation counts

    if (current_number_of_smry_fields.gt.0) then
       call mpiallsumint( domain_smry_1d(:)%count, global_smry_1d(:)%count, current_number_of_smry_fields, mpicom)
    end if

    !- MPI communications done ---

    ! Unpack results after MPI communication

    imax = 0
    imin = 0

    do ii = 1,current_number_of_smry_fields

      SELECT CASE (domain_smry_1d(ii)%cmpr_type)
      CASE (GREATER_EQ,ABS_GREATER_EQ)
        imax = imax +1
        global_smry_1d(ii)%extreme_val = rcv_array_mpimax(imax)
  
      CASE (SMALLER_THAN,ABS_SMALLER_THAN)
        imin = imin +1
        global_smry_1d(ii)%extreme_val = rcv_array_mpimin(imin)
      END SELECT
    end do

#else

    global_smry_1d(:)%extreme_val = domain_smry_1d(:)%extreme_val 
    global_smry_1d(:)%count       = domain_smry_1d(:)%count

#endif

    !-------------------------------
    ! Print messages to log file
    !-------------------------------
    ! Header line

    if (masterproc .and.  &
        any(global_smry_1d(1:current_number_of_smry_fields)%count.gt.0) ) then

       write(iulog,*)
       write(iulog,'(a19,a8, a40,a12, a10,a11,a2,a8, a13,a7, a9,a9,a5,a10,a5)')  &
                   'GLB_VERIF_SMRY:','nstep',  &
                   'Field','Unit',      &
                   'Cmpr.','Threshold','','Count',  &
                   'Extreme', 'Fixer',              &
                   '(Lat','Lon','Lev','Chunk','Col)'
    end if

    ! Summary of each field

    do ii = 1,current_number_of_smry_fields

      SELECT CASE (global_smry_1d(ii)%cmpr_type)
      CASE( GREATER_EQ )
        cmpr_type_char = '>='

      CASE( ABS_GREATER_EQ)
        cmpr_type_char = 'ABS >='

      CASE( SMALLER_THAN)
        cmpr_type_char = '<'

      CASE( ABS_SMALLER_THAN)
        cmpr_type_char = 'ABS <'

      END SELECT

      ! Master proc prints out the global summary

      if (masterproc .and. (global_smry_1d(ii)%count.gt.0 .or. l_print_smry_for_all_fields) ) then

        write(iulog,'(a19,i8, a40,a12, a10,e11.3,a2, i8, e13.3,i7)') &
                    'GLB_VERIF_SMRY:',nstep, &
                    trim(global_smry_1d(ii)%field_name),             &
                    trim(global_smry_1d(ii)%field_unit),             &
                    trim(cmpr_type_char), global_smry_1d(ii)%threshold, ':', &
                    global_smry_1d(ii)%count,       &
                    global_smry_1d(ii)%extreme_val, &
                    global_smry_1d(ii)%fixer 
      end if

      ! Every proc compars its own extreme value with the global extreme.
      ! When equal, print location of extreme value to log file in verbose mode

      if ( (glb_verif_smry_level.ge.1) .and.    &! print location 
           (global_smry_1d(ii)%count.gt.0 .or. l_print_smry_for_all_fields) ) then

      if ( global_smry_1d(ii)%extreme_val.eq.domain_smry_1d(ii)%extreme_val ) then

         write(iulog,'(a19,i8, a40,a12, a10,e11.3,a2, i8, e13.3,i7, 2f9.2,i5,i10,i5)') &
                     'GLB_VERIF_SMRY_VB:',nstep, &
                     trim(domain_smry_1d(ii)%field_name),                     &
                     trim(domain_smry_1d(ii)%field_unit) ,                     &
                     trim(cmpr_type_char), domain_smry_1d(ii)%threshold, ':', &
                     domain_smry_1d(ii)%count,               &
                     domain_smry_1d(ii)%extreme_val,         &
                     domain_smry_1d(ii)%fixer,               &
                     domain_smry_1d(ii)%extreme_lat*rad2deg, &
                     domain_smry_1d(ii)%extreme_lon*rad2deg, &
                     domain_smry_1d(ii)%extreme_lev,         &
                     domain_smry_1d(ii)%extreme_chnk,        &
                     domain_smry_1d(ii)%extreme_col
      end if
      end if

    end do

  end subroutine get_global_smry

end module glb_verif_smry
