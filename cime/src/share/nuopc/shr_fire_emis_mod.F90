module shr_fire_emis_mod

  !================================================================================
  ! Coordinates carbon emissions fluxes from CLM fires for use as sources of
  ! chemical constituents in CAM
  !
  ! This module reads fire_emis_nl namelist which specifies the compound fluxes
  ! that are to be passed through the model coupler.
  !================================================================================

  use ESMF                , only : ESMF_VMGetCurrent, ESMF_VM, ESMF_VMGet
  use ESMF                , only : ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU, ESMF_SUCCESS
  use shr_kind_mod        , only : r8 => shr_kind_r8, cl=>shr_kind_cl, cx=>shr_kind_cx, cs=>shr_kind_cs
  use shr_sys_mod         , only : shr_sys_abort
  use shr_log_mod         , only : logunit => shr_log_Unit
  use shr_mpi_mod         , only : shr_mpi_bcast
  use shr_nl_mod          , only : shr_nl_find_group_name
  use shr_expr_parser_mod , only : shr_exp_parse, shr_exp_item_t, shr_exp_list_destroy

  implicit none
  private

  public :: shr_fire_emis_readnl           ! reads fire_emis_nl namelist
  public :: shr_fire_emis_mechcomps        ! points to an array of chemical compounds (in CAM-Chem mechanism) than have fire emissions
  public :: shr_fire_emis_mechcomps_n      ! number of unique compounds in the CAM chemical mechanism that have fire emissions
  public :: shr_fire_emis_comps_n          ! number of unique emissions components
  public :: shr_fire_emis_linkedlist       ! points to linked list of shr_fire_emis_comp_t objects
  public :: shr_fire_emis_elevated         ! elevated emissions in ATM
  public :: shr_fire_emis_comp_ptr         ! user defined type that points to fire emis data obj (shr_fire_emis_comp_t)
  public :: shr_fire_emis_comp_t           ! emission component data type
  public :: shr_fire_emis_mechcomp_t       ! data type for chemical compound in CAM mechanism than has fire emissions

  logical :: fire_emis_initialized  = .false. ! true => shr_fire_emis_readnl alreay called
  logical :: shr_fire_emis_elevated = .true.

  character(len=CL), public :: shr_fire_emis_factors_file = ''       ! a table of basic fire emissions compounds
  character(len=CS), public :: shr_fire_emis_ztop_token = 'Sl_fztop' ! token for emissions top of vertical distribution
  integer, parameter        :: name_len=16

  ! fire emissions component data structure (or user defined type)
  type shr_fire_emis_comp_t
     character(len=name_len)             :: name            ! emissions component name (in fire emissions input table)
     integer                             :: index
     real(r8), pointer                   :: emis_factors(:) ! function of plant-function-type (PFT)
     real(r8)                            :: coeff           ! emissions component coeffecient
     real(r8)                            :: molec_weight    ! molecular weight of the fire emissions compound (g/mole)
     type(shr_fire_emis_comp_t), pointer :: next_emiscomp   ! points to next member in the linked list
  endtype shr_fire_emis_comp_t

  type shr_fire_emis_comp_ptr
    type(shr_fire_emis_comp_t), pointer :: ptr ! points to fire emis data obj (shr_fire_emis_comp_t)
  endtype shr_fire_emis_comp_ptr

  ! chemical compound in CAM mechanism that has fire emissions
  type shr_fire_emis_mechcomp_t
     character(len=16)             :: name                  ! compound name
     type(shr_fire_emis_comp_ptr), pointer :: emis_comps(:) ! an array of pointers to fire emis components
     integer                       :: n_emis_comps          ! number of fire emis compounds that make up the emissions for this mechanis compound
  end type shr_fire_emis_mechcomp_t

  type(shr_fire_emis_mechcomp_t),  pointer :: shr_fire_emis_mechcomps(:) ! array of chemical compounds (in CAM mechanism) that have fire emissions
  type(shr_fire_emis_comp_t),      pointer :: shr_fire_emis_linkedlist   ! points to linked list top

  integer :: shr_fire_emis_comps_n = 0      ! number of unique fire components
  integer :: shr_fire_emis_mechcomps_n = 0  ! number of unique compounds in the CAM chemical mechanism that have fire emissions

!-------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------

  subroutine shr_fire_emis_readnl( NLFileName, emis_nflds )

    !-------------------------------------------------------------------------
    !
    ! This reads the fire_emis_nl namelist group in drv_flds_in and parses the
    ! namelist information for the driver, CLM, and CAM.
    !
    ! Namelist variables:
    !   fire_emis_specifier, fire_emis_factors_file, fire_emis_elevated
    !
    !   fire_emis_specifier (array of strings) -- Each array element specifies
    !     how CAM-Chem constituents are mapped to basic smoke compounds in
    !     the fire emissions factors table (fire_emis_factors_file).  Each
    !     chemistry constituent name (left of '=' sign) is mapped to one or more
    !     smoke compound (separated by + sign if more than one), which can be
    !     proceeded by a multiplication factor (separated by '*').
    !     Example:
    !       fire_emis_specifier = 'bc_a1 = BC','pom_a1 = 1.4*OC','SO2 = SO2'
    !
    !   fire_emis_factors_file (string) -- Input file that contains the table
    !     of basic compounds that make up the smoke from the CLM fires.  This is
    !     used in CLM module FireEmisFactorsMod.
    !
    !   fire_emis_elevated (locical) -- If true then CAM-Chem treats the fire
    !     emission sources as 3-D vertically distributed forcings for the
    !     corresponding chemical tracers.
    !
    !-------------------------------------------------------------------------
    
    ! input/output variables
    character(len=*), intent(in)  :: NLFileName  ! name of namelist file
    integer         , intent(out) :: emis_nflds

    ! local variables
    type(ESMF_VM)       :: vm
    integer             :: localPet
    integer             :: mpicom
    integer             :: rc
    integer             :: unitn            ! namelist unit number
    integer             :: ierr             ! error code
    logical             :: exists           ! if file exists or not
    integer, parameter  :: maxspc = 100
    character(len=2*CX) :: fire_emis_specifier(maxspc) = ' '
    character(len=CL)   :: fire_emis_factors_file = ' '
    logical             :: fire_emis_elevated = .true.
    integer             :: i, tmp(1)
    character(*),parameter :: F00   = "('(shr_fire_emis_readnl) ',2a)"
    character(len=*), parameter :: subname='(shr_fire_emis_readnl)'
    !------------------------------------------------------------------

    namelist /fire_emis_nl/ fire_emis_specifier, fire_emis_factors_file, fire_emis_elevated

    rc = ESMF_SUCCESS

    ! If other processes have already initialized megan - then the info will just be re-initialized

    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return 

    call ESMF_VMGet(vm, localPet=localPet, mpiCommunicator=mpicom, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return 

    ! Note the following still needs to be called on all processors since the mpi_bcast is a collective 
    ! call on all the pes of mpicom
    if (localPet==0) then
       inquire( file=trim(NLFileName), exist=exists)
       if ( exists ) then
          open(newunit=unitn, file=trim(NLFilename), status='old' )
          write(logunit,F00) 'Read in fire_emis_readnl namelist from: ', trim(NLFilename)
          call shr_nl_find_group_name(unitn, 'fire_emis_nl', status=ierr)
          ! If ierr /= 0, no namelist present.
          if (ierr == 0) then
             read(unitn, fire_emis_nl, iostat=ierr)
             if (ierr > 0) then
                call shr_sys_abort( 'problem on read of fire_emis_nl namelist in shr_fire_emis_readnl' )
             endif
          endif
          close( unitn )
       end if
    end if
    call shr_mpi_bcast(fire_emis_specifier    , mpicom)
    call shr_mpi_bcast(fire_emis_factors_file , mpicom)
    call shr_mpi_bcast(fire_emis_elevated     , mpicom)

    shr_fire_emis_factors_file = fire_emis_factors_file
    shr_fire_emis_elevated = fire_emis_elevated

    ! parse the namelist info and initialize the module data - only if it has not been initialized
    if (.not. fire_emis_initialized) then
       call shr_fire_emis_init( fire_emis_specifier )
    end if
    emis_nflds = shr_fire_emis_mechcomps_n

  end subroutine shr_fire_emis_readnl

!-------------------------------------------------------------------------
! private methods...
!-------------------------------------------------------------------------

  subroutine shr_fire_emis_init( specifier )

    !--------------------------------------------------
    ! module data initializer
    !--------------------------------------------------

    ! input/output variables
    character(len=*), intent(in) :: specifier(:)

    ! local variables
    integer :: n_entries
    integer :: i, j, k
    type(shr_exp_item_t), pointer :: items_list, item
    !------------------------------------------------------

    nullify(shr_fire_emis_linkedlist)

    items_list => shr_exp_parse( specifier, nitems=n_entries )

    allocate(shr_fire_emis_mechcomps(n_entries))
    shr_fire_emis_mechcomps(:)%n_emis_comps = 0

    item => items_list
    i = 1
    do while(associated(item))

       do k=1,shr_fire_emis_mechcomps_n
          if ( trim(shr_fire_emis_mechcomps(k)%name) == trim(item%name) ) then
             call shr_sys_abort( 'shr_fire_emis_init : multiple emissions definitions specified for : '//trim(item%name))
          endif
       enddo
       if (len_trim(item%name) .le. name_len) then
          shr_fire_emis_mechcomps(i)%name = item%name(1:name_len)
       else
          call shr_sys_abort("shr_file_emis_init : name too long for data structure :"//trim(item%name))
       endif
       shr_fire_emis_mechcomps(i)%n_emis_comps = item%n_terms
       allocate(shr_fire_emis_mechcomps(i)%emis_comps(item%n_terms))

       do j = 1,item%n_terms
          shr_fire_emis_mechcomps(i)%emis_comps(j)%ptr => add_emis_comp( item%vars(j), item%coeffs(j) )
       enddo
       shr_fire_emis_mechcomps_n = shr_fire_emis_mechcomps_n+1

       item => item%next_item
       i = i+1
    enddo
    if (associated(items_list)) call shr_exp_list_destroy(items_list)

    ! Need to explicitly add Fl_ based on naming convention

    fire_emis_initialized = .true.

  end subroutine shr_fire_emis_init

  !-------------------------------------------------------------------------

  function add_emis_comp( name, coeff ) result(emis_comp)

    character(len=*), intent(in) :: name
    real(r8),         intent(in) :: coeff
    type(shr_fire_emis_comp_t), pointer :: emis_comp

    emis_comp => get_emis_comp_by_name(shr_fire_emis_linkedlist, name)
    if(associated(emis_comp)) then
       ! already in the list so return...
       return
    endif

    ! create new emissions component and add it to the list
    allocate(emis_comp)

    !    element%index = lookup_element( name )
    !    element%emis_factors = get_factors( list_elem%index )

    emis_comp%index = shr_fire_emis_comps_n+1

    emis_comp%name = trim(name)
    emis_comp%coeff = coeff
    nullify(emis_comp%next_emiscomp)

    call add_emis_comp_to_list(emis_comp)

  end function add_emis_comp

  !-------------------------------------------------------------------------

  recursive function get_emis_comp_by_name(list_comp, name) result(emis_comp)

    type(shr_fire_emis_comp_t), pointer  :: list_comp
    character(len=*), intent(in) :: name  ! variable name
    type(shr_fire_emis_comp_t), pointer  :: emis_comp ! returned object

    if(associated(list_comp)) then
       if(list_comp%name .eq. name) then
          emis_comp => list_comp
       else
          emis_comp => get_emis_comp_by_name(list_comp%next_emiscomp, name)
       end if
    else
       nullify(emis_comp)
    end if

  end function get_emis_comp_by_name

  !-------------------------------------------------------------------------

  subroutine add_emis_comp_to_list( new_emis_comp )

    type(shr_fire_emis_comp_t), target, intent(in) :: new_emis_comp

    type(shr_fire_emis_comp_t), pointer :: list_comp

    if(associated(shr_fire_emis_linkedlist)) then
       list_comp => shr_fire_emis_linkedlist
       do while(associated(list_comp%next_emiscomp))
          list_comp => list_comp%next_emiscomp
       end do
       list_comp%next_emiscomp => new_emis_comp
    else
       shr_fire_emis_linkedlist => new_emis_comp
    end if

    shr_fire_emis_comps_n = shr_fire_emis_comps_n + 1

  end subroutine add_emis_comp_to_list

endmodule shr_fire_emis_mod
