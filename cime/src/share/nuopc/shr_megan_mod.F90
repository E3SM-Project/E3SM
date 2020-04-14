module shr_megan_mod

  !================================================================================
  ! Handles MEGAN VOC emissions metadata for CLM produced chemical emissions
  ! MEGAN = Model of Emissions of Gases and Aerosols from Nature
  !
  ! This reads the megan_emis_nl namelist in drv_flds_in and makes the relavent
  ! information available to CAM, CLM, and driver. 
  ! - The driver sets up CLM to CAM communication for the  VOC flux fields. 
  ! - CLM needs to know what specific VOC fluxes need to be passed to the coupler 
  !   and how to assemble the fluxes.
  ! - CAM needs to know what specific VOC fluxes to expect from CLM.
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

  public :: shr_megan_readnl           ! reads megan_emis_nl namelist
  public :: shr_megan_mechcomps        ! points to an array of chemical compounds (in CAM-Chem mechanism) that have MEGAN emissions
  public :: shr_megan_mechcomps_n      ! number of unique compounds in the CAM chemical mechanism that have MEGAN emissions
  public :: shr_megan_megcomps_n       ! number of unique MEGAN compounds
  public :: shr_megan_megcomp_t        ! MEGAN compound data type
  public :: shr_megan_mechcomp_t       ! data type for chemical compound in CAM mechanism that has MEGAN emissions
  public :: shr_megan_linkedlist       ! points to linked list of shr_megan_comp_t objects
  public :: shr_megan_mapped_emisfctrs ! switch to use mapped emission factors
  public :: shr_megan_comp_ptr

  logical          , public :: megan_initialized       = .false. ! true => shr_megan_readnl alreay called
  character(len=CL), public :: shr_megan_factors_file  = ''

  ! MEGAN compound data structure (or user defined type)
  type shr_megan_megcomp_t
     character(len=16)     :: name            ! MEGAN compound name (in MEGAN input table)
     integer               :: index
     real(r8), pointer     :: emis_factors(:) ! function of plant-function-type (PFT)
     integer               :: class_number    ! MEGAN class number
     real(r8)              :: coeff           ! emissions component coeffecient
     real(r8)              :: molec_weight    ! molecular weight of the MEGAN compound (g/mole)
     type(shr_megan_megcomp_t), pointer :: next_megcomp ! points to next member in the linked list
  endtype shr_megan_megcomp_t

  type shr_megan_comp_ptr
    type(shr_megan_megcomp_t), pointer :: ptr
  endtype shr_megan_comp_ptr

  ! chemical compound in CAM mechanism that has MEGAN emissions
  type shr_megan_mechcomp_t
     character(len=16)                 :: name           ! compound name
     type(shr_megan_comp_ptr), pointer :: megan_comps(:) ! an array of pointers to megan emis compounds
     integer                           :: n_megan_comps  ! number of megan emis compounds that make up the emissions for this mechanis compound
  end type shr_megan_mechcomp_t

  type(shr_megan_mechcomp_t), pointer :: shr_megan_mechcomps(:) ! array of chemical compounds (in CAM mechanism) that have MEGAN emissions
  type(shr_megan_megcomp_t),  pointer :: shr_megan_linkedlist   ! points to linked list top

  integer :: shr_megan_megcomps_n  = 0          ! number of unique megan compounds
  integer :: shr_megan_mechcomps_n = 0          ! number of unique compounds in the CAM chemical mechanism that have MEGAN emissions

  ! switch to use mapped emission factors
  logical :: shr_megan_mapped_emisfctrs = .false.

!--------------------------------------------------------
contains
!--------------------------------------------------------

  subroutine shr_megan_readnl( NLFileName, megan_nflds)

    !-------------------------------------------------------------------------
    !
    ! This reads the megan_emis_nl namelist group in drv_flds_in and parses the
    ! namelist information for the driver, CLM, and CAM.
    !
    ! Namelist variables:
    !   megan_specifier, megan_mapped_emisfctrs, megan_factors_file
    !
    ! megan_specifier is a series of strings where each string contains one
    !  CAM chemistry constituent name (left of = sign) and one or more MEGAN
    !  compound (separated by + sign if more than one).  Each MEGAN compound
    !  can be proceeded by a multiplication factor (separated by *).  The
    !  specification of the MEGAN compounds to the right of the = signs tells
    !  the MEGAN VOC model within CLM how to construct the VOC fluxes using
    !  the factors in megan_factors_file and land surface state.
    !
    ! megan_factors_file read by CLM contains valid MEGAN compound names,
    !  MEGAN class groupings and scalar emission factors
    !
    ! megan_mapped_emisfctrs switch is used to tell the MEGAN model to use
    !  mapped emission factors read in from the CLM surface data input file
    !  rather than the scalar factors from megan_factors_file
    !
    ! Example:
    ! &megan_emis_nl
    !  megan_specifier = 'ISOP = isoprene',
    !     'C10H16 = myrcene + sabinene + limonene + carene_3 + ocimene_t_b + pinene_b + ...',
    !     'CH3OH = methanol',
    !     'C2H5OH = ethanol',
    !     'CH2O = formaldehyde',
    !     'CH3CHO = acetaldehyde',
    ! ...
    !  megan_factors_file = '$datapath/megan_emis_factors.nc'
    ! /
    !-------------------------------------------------------------------------
    
    ! input/output variables
    character(len=*), intent(in)  :: NLFileName
    integer,          intent(out) :: megan_nflds

    ! local variables
    type(ESMF_VM)       :: vm
    integer             :: localPet
    integer             :: mpicom
    integer             :: unitn            ! namelist unit number
    integer             :: ierr             ! error code
    logical             :: exists           ! if file exists or not
    integer, parameter  :: maxspc = 100
    character(len=2*CX) :: megan_specifier(maxspc) = ' '
    logical             :: megan_mapped_emisfctrs = .false.
    character(len=CL)   :: megan_factors_file = ' '
    integer             :: rc
    integer             :: i, tmp(1)
    character(*), parameter :: F00   = "('(shr_megan_readnl) ',2a)"
    character(len=*), parameter :: subname='(shr_megan_readnl)'
    !--------------------------------------------------------------

    namelist /megan_emis_nl/ megan_specifier, megan_factors_file, megan_mapped_emisfctrs

    !--- Open and read namelist ---
    if ( len_trim(NLFilename) == 0 ) then
       call shr_sys_abort( subName//'ERROR: nlfilename not set' )
    end if

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
          write(logunit,F00) 'Read in megan_emis_readnl namelist from: ', trim(NLFilename)
          call shr_nl_find_group_name(unitn, 'megan_emis_nl', status=ierr)
          if (ierr == 0) then
             ! Note that ierr /= 0, no namelist is present.
             read (unitn, megan_emis_nl, iostat=ierr)
             if (ierr > 0) then
                call shr_sys_abort( 'problem on read of megan_emis_nl namelist in shr_megan_readnl' )
             endif
          endif
          close( unitn )
       end if
    end if
    call shr_mpi_bcast( megan_specifier        , mpicom )
    call shr_mpi_bcast( megan_factors_file     , mpicom )
    call shr_mpi_bcast( megan_mapped_emisfctrs , mpicom )

    shr_megan_factors_file = megan_factors_file
    shr_megan_mapped_emisfctrs = megan_mapped_emisfctrs

    ! parse the namelist info and initialize the module data - only if it has not been initialized
    if (.not. megan_initialized) then
       call shr_megan_init( megan_specifier )
    end if
    megan_nflds = shr_megan_mechcomps_n

  end subroutine shr_megan_readnl

!-------------------------------------------------------------------------
! private methods...
!-------------------------------------------------------------------------

  subroutine shr_megan_init( specifier)

    !-----------------------------------------
    ! Initialize module data
    !-----------------------------------------

    ! input/output variables
    character(len=*), intent(in) :: specifier(:)

    ! local variables
    integer                       :: n_entries
    integer                       :: i, j, k
    type(shr_exp_item_t), pointer :: items_list, item
    !--------------------------------------------------------------

    nullify(shr_megan_linkedlist)

    items_list => shr_exp_parse( specifier, nitems=n_entries )

    allocate(shr_megan_mechcomps(n_entries))
    shr_megan_mechcomps(:)%n_megan_comps = 0

    item => items_list
    i = 1
    do while(associated(item))

       do k=1,shr_megan_mechcomps_n
          if ( trim(shr_megan_mechcomps(k)%name) == trim(item%name) ) then
             call shr_sys_abort( 'shr_megan_init : duplicate compound names : '//trim(item%name))
          endif
       enddo
       if (len_trim(item%name) .le. len(shr_megan_mechcomps(i)%name)) then
          shr_megan_mechcomps(i)%name = item%name(1:len(shr_megan_mechcomps(i)%name))
       else
          call shr_sys_abort( 'shr_megan_init : name too long for data structure : '//trim(item%name))
       endif
       shr_megan_mechcomps(i)%n_megan_comps = item%n_terms
       allocate(shr_megan_mechcomps(i)%megan_comps(item%n_terms))

       do j = 1,item%n_terms
          shr_megan_mechcomps(i)%megan_comps(j)%ptr => add_megan_comp( item%vars(j), item%coeffs(j) )
       enddo
       shr_megan_mechcomps_n = shr_megan_mechcomps_n+1

       item => item%next_item
       i = i+1

    enddo
    if (associated(items_list)) call shr_exp_list_destroy(items_list)

    megan_initialized = .true.

  end subroutine shr_megan_init

  !-------------------------------------------------------------------------

  function add_megan_comp( name, coeff ) result(megan_comp)

    character(len=16), intent(in) :: name
    real(r8),          intent(in) :: coeff
    type(shr_megan_megcomp_t), pointer :: megan_comp

    megan_comp => get_megan_comp_by_name(shr_megan_linkedlist, name)
    if(associated(megan_comp)) then
       ! already in the list so return...
       return
    endif

    ! create new megan compound and add it to the list
    allocate(megan_comp)

    !    element%index = lookup_element( name )
    !    element%emis_factors = get_factors( list_elem%index )

    megan_comp%index = shr_megan_megcomps_n+1

    megan_comp%name = trim(name)
    megan_comp%coeff = coeff
    nullify(megan_comp%next_megcomp)

    call add_megan_comp_to_list(megan_comp)

  end function add_megan_comp

  !-------------------------------------------------------------------------

  recursive function get_megan_comp_by_name(list_comp, name) result(megan_comp)

    type(shr_megan_megcomp_t), pointer  :: list_comp
    character(len=*), intent(in) :: name  ! variable name
    type(shr_megan_megcomp_t), pointer  :: megan_comp ! returned object

    if(associated(list_comp)) then
       if(list_comp%name .eq. name) then
          megan_comp => list_comp
       else
          megan_comp => get_megan_comp_by_name(list_comp%next_megcomp, name)
       end if
    else
       nullify(megan_comp)
    end if

  end function get_megan_comp_by_name

  !-------------------------------------------------------------------------

  subroutine add_megan_comp_to_list( new_megan_comp )

    type(shr_megan_megcomp_t), target, intent(in) :: new_megan_comp

    type(shr_megan_megcomp_t), pointer :: list_comp

    if(associated(shr_megan_linkedlist)) then
       list_comp => shr_megan_linkedlist
       do while(associated(list_comp%next_megcomp))
          list_comp => list_comp%next_megcomp
       end do
       list_comp%next_megcomp => new_megan_comp
    else
       shr_megan_linkedlist => new_megan_comp
    end if

    shr_megan_megcomps_n = shr_megan_megcomps_n + 1

  end subroutine add_megan_comp_to_list

endmodule shr_megan_mod
