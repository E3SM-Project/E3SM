!================================================================================
! Handles MEGAN VOC emissions metadata for CLM produced chemical emissions
! MEGAN = Model of Emissions of Gases and Aerosols from Nature
!
! This reads the megan_emis_nl namelist in drv_flds_in and makes the relavent
! information available to CAM, CLM, and driver. The driver sets up CLM to CAM
! communication for the  VOC flux fields. CLM needs to know what specific VOC
! fluxes need to be passed to the coupler and how to assimble the fluxes.
! CAM needs to know what specific VOC fluxes to expect from CLM.
!
! Francis Vitt -- 26 Oct 2011
!================================================================================
module shr_megan_mod

  use shr_kind_mod,only : r8 => shr_kind_r8
  use shr_kind_mod,only : CL => SHR_KIND_CL, CX => SHR_KIND_CX, CS => SHR_KIND_CS
  use shr_sys_mod, only : shr_sys_abort
  use shr_log_mod, only : loglev  => shr_log_Level
  use shr_log_mod, only : logunit => shr_log_Unit

  implicit none
  save
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
  character(len=CS), public :: shr_megan_fields_token  = ''      ! First drydep fields token
  character(len=CL), public :: shr_megan_factors_file  = ''
  character(len=CX), public :: shr_megan_fields        = ''

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

contains

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
  subroutine shr_megan_readnl( NLFileName, megan_fields, megan_nflds )
    use ESMF, only : ESMF_VM, ESMF_VMGetCurrent, ESMF_VMBroadcast, ESMF_VMGet
    use shr_nl_mod,     only : shr_nl_find_group_name
    use shr_file_mod,   only : shr_file_getUnit, shr_file_freeUnit

    character(len=*), intent(in)  :: NLFileName
    character(len=*), intent(out) :: megan_fields
    integer,          intent(out) :: megan_nflds

    type(ESMF_VM)   :: vm
    integer :: localPet
    integer :: unitn            ! namelist unit number
    integer :: ierr             ! error code
    logical :: exists           ! if file exists or not
    integer, parameter  :: maxspc = 100
    character(len=2*CX) :: megan_specifier(maxspc) = ' '
    logical             :: megan_mapped_emisfctrs = .false.
    character(len=CL)   :: megan_factors_file = ' '
    integer :: rc
    integer :: i, tmp(1)
    character(*),parameter :: F00   = "('(shr_megan_readnl) ',2a)"

    namelist /megan_emis_nl/ megan_specifier, megan_factors_file, megan_mapped_emisfctrs

    ! If other processes have already initialized megan - then just return
    ! the megan_fields that have already been set
    if (megan_initialized) then
       megan_fields = trim(shr_megan_fields)
       megan_nflds = shr_megan_mechcomps_n
       return
    end if
    call ESMF_VMGetCurrent(vm, rc=rc)
    call ESMF_VMGet(vm, localpet=localpet, rc=rc)
    megan_nflds = 0
    if (localPet==0) then
       inquire( file=trim(NLFileName), exist=exists)
       if ( exists ) then
          unitn = shr_file_getUnit()
          open( unitn, file=trim(NLFilename), status='old' )
          if ( loglev > 0 ) write(logunit,F00) &
               'Read in megan_emis_readnl namelist from: ', trim(NLFilename)

          call shr_nl_find_group_name(unitn, 'megan_emis_nl', status=ierr)
          ! If ierr /= 0, no namelist present.

          if (ierr == 0) then
             read(unitn, megan_emis_nl, iostat=ierr)

             if (ierr > 0) then
                call shr_sys_abort( 'problem on read of megan_emis_nl namelist in shr_megan_readnl' )
             endif
          endif

          close( unitn )
          call shr_file_freeUnit( unitn )
          do i=1,maxspc
             if(len_trim(megan_specifier(i)) > 0) then
                megan_nflds=megan_nflds+1
             endif
          enddo
       end if
    end if
    tmp = megan_nflds
    call ESMF_VMBroadcast(vm, tmp, 1, 0, rc=rc)
    megan_nflds = tmp(1)
    if(megan_nflds > 0) then
       call ESMF_VMBroadcast(vm, megan_specifier, 2*CX*megan_nflds, 0, rc=rc)
       call ESMF_VMBroadcast(vm, megan_factors_file, CL, 0, rc=rc)
       tmp = 0
       if(megan_mapped_emisfctrs) tmp=1
       call ESMF_VMBroadcast(vm, tmp, 1, 0, rc=rc)
       if(tmp(1)==1) megan_mapped_emisfctrs=.true.
    endif

    shr_megan_factors_file = megan_factors_file
    shr_megan_mapped_emisfctrs = megan_mapped_emisfctrs

    ! parse the namelist info and initialize the module data
    call shr_megan_init( megan_specifier, megan_fields )
  end subroutine shr_megan_readnl

  !-------------------------------------------------------------------------
  ! module data initializer
  !-------------------------------------------------------------------------
  subroutine shr_megan_init( specifier, megan_fields )

    use shr_expr_parser_mod, only : shr_exp_parse, shr_exp_item_t, shr_exp_list_destroy

    character(len=*), intent(in) :: specifier(:)
    character(len=*), intent(out) :: megan_fields

    integer :: n_entries
    integer :: i, j, k

    type(shr_exp_item_t), pointer :: items_list, item
    character(len=12) :: token   ! megan field name to add

    nullify(shr_megan_linkedlist)

    items_list => shr_exp_parse( specifier, nitems=n_entries )

    allocate(shr_megan_mechcomps(n_entries))
    shr_megan_mechcomps(:)%n_megan_comps = 0

    megan_fields = ''

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

       write(token,333) shr_megan_mechcomps_n

       if ( shr_megan_mechcomps_n == 1 ) then
          ! do not prepend ":" to the string for the first token
          megan_fields = trim(token)
          shr_megan_fields_token = token
       else
          megan_fields = trim(megan_fields)//':'//trim(token)
       endif

       item => item%next_item
       i = i+1
    enddo
    if (associated(items_list)) call shr_exp_list_destroy(items_list)

    megan_initialized = .true.
    shr_megan_fields = trim(megan_fields)

    ! Need to explicitly add Fl_ based on naming convention
333 format ('Fall_voc',i3.3)

  end subroutine shr_megan_init

  !-------------------------------------------------------------------------
  ! private methods...

  !-------------------------------------------------------------------------
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
