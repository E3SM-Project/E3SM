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
  public :: shr_megan_mechcomps        ! points to an array of chemical compounds (in CAM-Chem mechanism) than have MEGAN emissions
  public :: shr_megan_mechcomps_n      ! number of unique compounds in the CAM chemical mechanism  than have MEGAN emissions
  public :: shr_megan_megcomps_n       ! number of unique MEGAN compounds
  public :: shr_megan_megcomp_t        ! MEGAN compound data type
  public :: shr_megan_mechcomp_t       ! data type for chemical compound in CAM mechanism than has MEGAN emissions
  public :: shr_megan_linkedlist       ! points to linked list of shr_megan_comp_t objects
  public :: shr_megan_mapped_emisfctrs ! switch to use mapped emission factors
  public :: shr_megan_comp_ptr

  character(len=CS), public :: shr_megan_fields_token = ''   ! First drydep fields token
  character(len=CL), public :: shr_megan_factors_file = ''

  ! MEGAN compound data structure (or user defined type)
  type shr_megan_megcomp_t
     character(len=16)     :: name            ! MEGAN compound name (in MEGAN input table)
     integer               :: index
     real(r8), pointer     :: emis_factors(:) ! function of plant-function-type (PFT)
     integer               :: class_number    ! MEGAN class number
     real(r8)              :: molec_weight    ! molecular weight of the MEGAN compound (g/mole)
     type(shr_megan_megcomp_t), pointer :: next_megcomp ! points to next member in the linked list
  endtype shr_megan_megcomp_t

  type shr_megan_comp_ptr
    type(shr_megan_megcomp_t), pointer :: ptr
  endtype shr_megan_comp_ptr

  ! chemical compound in CAM mechanism than has MEGAN emissions
  type shr_megan_mechcomp_t
     character(len=16)             :: name           ! compound name
     type(shr_megan_comp_ptr), pointer :: megan_comps(:) ! an array of pointers to megan emis compounds 
     integer                       :: n_megan_comps  ! number of megan emis compounds than make up the emissions for this mechanis compound
  end type shr_megan_mechcomp_t

  type(shr_megan_mechcomp_t), pointer :: shr_megan_mechcomps(:) ! array of chemical compounds (in CAM mechanism) than have MEGAN emissions
  type(shr_megan_megcomp_t),  pointer :: shr_megan_linkedlist   ! points to linked list top

  integer :: shr_megan_megcomps_n  = 0          ! number of unique megan compounds
  integer :: shr_megan_mechcomps_n = 0          ! number of unique compounds in the CAM chemical mechanism  than have MEGAN emissions

  ! switch to use mapped emission factors
  logical :: shr_megan_mapped_emisfctrs = .false.

  ! private data 
  type parser_items_t
     character(len=16),pointer :: megan_comp_names(:)
     character(len=16) :: mech_comp_name
     integer :: n_megan_comps
  end type parser_items_t

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
  !  compounds (seperated by + sign if more than one).  The specification of 
  !  the MEGAN compounds to the right of the = signs tells the MEGAN VOC 
  !  model within CLM how to construct the VOC fluxes using the factors in 
  !  megan_factors_file and land surface state.
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
  subroutine shr_megan_readnl( NLFileName, megan_fields )

    use shr_nl_mod,     only : shr_nl_find_group_name
    use shr_file_mod,   only : shr_file_getUnit, shr_file_freeUnit

    character(len=*), intent(in)  :: NLFileName
    character(len=*), intent(out) :: megan_fields	

    integer :: unitn            ! namelist unit number
    integer :: ierr             ! error code
    logical :: exists           ! if file exists or not

    integer, parameter :: maxspc = 100

    character(len=2*CX) :: megan_specifier(maxspc) = ' '
    logical           :: megan_mapped_emisfctrs = .false.
    character(len=CL) :: megan_factors_file = ' '

    character(*),parameter :: F00   = "('(seq_drydep_read) ',2a)" 

    namelist /megan_emis_nl/ megan_specifier, megan_factors_file, megan_mapped_emisfctrs

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

       shr_megan_factors_file = megan_factors_file
       shr_megan_mapped_emisfctrs = megan_mapped_emisfctrs

       ! parse the namelist info and initialize the module data
       call shr_megan_init( megan_specifier, megan_fields )

       close( unitn )
       call shr_file_freeUnit( unitn )

    end if

  end subroutine shr_megan_readnl

  !-------------------------------------------------------------------------
  ! module data initializer
  !-------------------------------------------------------------------------
  subroutine shr_megan_init( specifier, megan_fields )

    character(len=*), intent(in) :: specifier(:)
    character(len=*), intent(out) :: megan_fields	

    integer :: n_entries
    integer :: i, j, k
    integer :: spc_len
    type(parser_items_t), pointer :: items

    character(len=12) :: token   ! megan field name to add

    nullify(shr_megan_linkedlist)

    n_entries = size(specifier)
    allocate(shr_megan_mechcomps(n_entries))
    shr_megan_mechcomps(:)%n_megan_comps = 0

    megan_fields = ''

    do i = 1,n_entries
       spc_len=len_trim(specifier(i))
       if ( spc_len > 0 ) then

          items => get_parser_items( specifier(i) )

          do k=1,shr_megan_mechcomps_n
             if ( trim(shr_megan_mechcomps(k)%name) == trim(items%mech_comp_name) ) then
                call shr_sys_abort( 'shr_megan_init : duplicate compound names : '//trim(items%mech_comp_name))
             endif
          enddo

          shr_megan_mechcomps(i)%name = items%mech_comp_name
          shr_megan_mechcomps(i)%n_megan_comps = items%n_megan_comps
          allocate(shr_megan_mechcomps(i)%megan_comps(items%n_megan_comps))

          do j = 1,items%n_megan_comps
             shr_megan_mechcomps(i)%megan_comps(j)%ptr => add_megan_comp( items%megan_comp_names(j) )
          enddo
          shr_megan_mechcomps_n = shr_megan_mechcomps_n+1

          call destroy_parser_items( items )

          write(token,333) shr_megan_mechcomps_n

          if ( shr_megan_mechcomps_n == 1 ) then
             ! no not prepend ":" to the string for the first token
             megan_fields = trim(token)
             shr_megan_fields_token = token
          else
             megan_fields = trim(megan_fields)//':'//trim(token)                 
          endif

       endif

    enddo

    ! Need to explicitly add Fl_ based on naming convention
333 format ('Fall_voc',i3.3)

  end subroutine shr_megan_init

  !-------------------------------------------------------------------------
  ! private methods...

  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  function get_parser_items( spec_entry ) result(items)

    character(len=*), intent(in) :: spec_entry

    type(parser_items_t), pointer :: items ! items returned

    integer :: ndxs(512)
    integer :: nelem, j, i
    character(len=CL) :: tmp_str

    j = scan( spec_entry, '=' )

    nelem = 1
    ndxs(nelem) = j

    tmp_str = trim( spec_entry(j+1:) )
    j = scan( tmp_str, '+' )

    do while(j>0)

       nelem = nelem+1
       ndxs(nelem) = ndxs(nelem-1) + j

       tmp_str = tmp_str(j+1:)
       j = scan( tmp_str, '+' )

    enddo
    ndxs(nelem+1) = len(spec_entry)+1

    allocate(items)
    allocate(items%megan_comp_names(nelem))
    items%mech_comp_name = trim(adjustl( spec_entry(:ndxs(1)-1)))
    items%n_megan_comps = nelem
    do i = 1,nelem 
       items%megan_comp_names(i) = trim(adjustl( spec_entry(ndxs(i)+1:ndxs(i+1)-1)))
    enddo

  endfunction get_parser_items

  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  subroutine destroy_parser_items( items )
    type(parser_items_t), pointer :: items

    deallocate( items%megan_comp_names )
    deallocate( items )
    nullify( items )
  endsubroutine destroy_parser_items

  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  function add_megan_comp( name ) result(megan_comp)

    character(len=16), intent(in) :: name
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
