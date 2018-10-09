! This code is part of RRTM for GCM Applications - Parallel (RRTMGP)
!
! Contacts: Robert Pincus and Eli Mlawer
! email:  rrtmgp@aer.com
!
! Copyright 2015-2018,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
! -------------------------------------------------------------------------------------------------
! Encapsulates a collection of volume mixing ratios (concentrations) of gases.
!   Each concentration is associated with a name, normally the chemical formula.
!
! Values may be provided as scalars, 1-dimensional profiles (nlay), or 2-D fields (ncol,nlay).
!   (nlay and ncol are determined from the input arrays; self-consistency is enforced)
!   example:
!   error_msg = gas_concs%set_vmr('h2o', values(:,:))
!   error_msg = gas_concs%set_vmr('o3' , values(:)  )
!   error_msg = gas_concs%set_vmr('co2', value      )
!
! Values can be requested as profiles (valid only if there are no 2D fields present in the object)
!   or as 2D fields. Values for all columns are returned although the entire collection
!   can be subsetted in the column dimension
!
! Subsets can be extracted in the column dimension
!
! Functions return strings. Non-empty strings indicate an error.
!
! -------------------------------------------------------------------------------------------------

module mo_gas_concentrations
  use mo_rte_kind,    only: wp
  use mo_util_string, only: lower_case
  implicit none
  integer, parameter :: GAS_NOT_IN_LIST = -1

  type, private :: conc_field
    real(wp), dimension(:,:), allocatable :: conc
  end type conc_field

  type, public :: ty_gas_concs
    !
    ! Data
    !
    character(len=32), dimension(:), allocatable :: gas_name
    type(conc_field),  dimension(:), allocatable :: concs
    integer :: ncol = 0, nlay = 0
    contains
      !
      ! Procedures
      !
      procedure, private :: increase_list_size
      procedure, private :: find_gas
      procedure, private :: set_vmr_scalar
      procedure, private :: set_vmr_1d
      procedure, private :: set_vmr_2d
      procedure, private :: get_vmr_1d
      procedure, private :: get_vmr_2d
      procedure, private :: get_subset_range
      !
      ! public interface
      !
      procedure, public :: reset
      generic,   public :: set_vmr => set_vmr_scalar, &
                                      set_vmr_1d, &
                                      set_vmr_2d
      generic,   public :: get_vmr => get_vmr_1d, &
                                      get_vmr_2d
      generic,   public :: get_subset => get_subset_range
  end type ty_gas_concs
contains
  ! -------------------------------------------------------------------------------------
  !
  ! Set concentrations --- scalar, 1D, 2D
  !
  ! -------------------------------------------------------------------------------------
  function set_vmr_scalar(this, gas, w) result(error_msg)
    class(ty_gas_concs), intent(inout) :: this
    character(len=*),    intent(in   ) :: gas
    real(wp),            intent(in   ) :: w
    character(len=128)                 :: error_msg
    ! ---------
    integer :: igas
    ! ---------
    error_msg = ''
    if (w < 0._wp) then
      error_msg = 'ty_gas_concs%set_vmr: concentrations should be >= 0'
      return
    endif

    igas = this%find_gas(gas)
    if (igas == GAS_NOT_IN_LIST) then
      call this%increase_list_size()
      igas = size(this%gas_name)
    end if
    !
    ! Deallocate anything existing -- could be more efficient to test if it's already the correct size
    !
    if (allocated(this%concs(igas)%conc)) deallocate(this%concs(igas)%conc)
    allocate(this%concs(igas)%conc(1,1))
    this%concs(igas)%conc(:,:) = w
    this%gas_name(igas) = trim(gas)
  end function set_vmr_scalar
  ! -------------------------------------------------------------------------------------
  function set_vmr_1d(this, gas, w) result(error_msg)
    class(ty_gas_concs), intent(inout) :: this
    character(len=*),    intent(in   ) :: gas
    real(wp), dimension(:), &
                         intent(in   ) :: w
    character(len=128)                 :: error_msg
    ! ---------
    integer :: igas
    ! ---------
    error_msg = ''

    if (any(w < 0._wp)) then
      error_msg = 'ty_gas_concs%set_vmr: concentrations should be >= 0'
    endif
    if(this%nlay > 0) then
      if(size(w) /= this%nlay) error_msg = 'ty_gas_concs%set_vmr: different dimension (nlay)'
    else
      this%nlay = size(w)
    end if
    if(error_msg /= "") return

    igas = this%find_gas(gas)
    if (igas == GAS_NOT_IN_LIST) then
      call this%increase_list_size()
      igas = size(this%gas_name)
    end if
    !
    ! Deallocate anything existing -- could be more efficient to test if it's already the correct size
    !
    if (allocated(this%concs(igas)%conc)) deallocate(this%concs(igas)%conc)
    allocate(this%concs(igas)%conc(1,this%nlay))
    this%concs(igas)%conc(1,:) = w
    this%gas_name(igas) = trim(gas)
  end function set_vmr_1d
  ! --------------------
  function set_vmr_2d(this, gas, w) result(error_msg)
    class(ty_gas_concs), intent(inout) :: this
    character(len=*),    intent(in   ) :: gas
    real(wp), dimension(:,:),  &
                         intent(in   ) :: w
    character(len=128)                 :: error_msg
    ! ---------
    integer :: igas
    ! ---------
    error_msg = ''

    if (any(w < 0._wp)) then
      error_msg = 'ty_gas_concs%set_vmr: concentrations should be >= 0'
    endif
    if(this%ncol > 0 .and. size(w, 1) /= this%ncol) then
      error_msg = 'ty_gas_concs%set_vmr: different dimension (ncol)'
    else
      this%ncol = size(w, 1)
    end if
    if(this%nlay > 0 .and. size(w, 2) /= this%nlay) then
      error_msg = 'ty_gas_concs%set_vmr: different dimension (nlay)'
    else
      this%nlay = size(w, 2)
    end if
    if(error_msg /= "") return

    igas = this%find_gas(gas)
    if (igas == GAS_NOT_IN_LIST) then
      call this%increase_list_size()
      igas = size(this%gas_name)
    end if
    !
    ! Deallocate anything existing -- could be more efficient to test if it's already the correct size
    !
    if (allocated(this%concs(igas)%conc)) deallocate(this%concs(igas)%conc)
    allocate(this%concs(igas)%conc(this%ncol,this%nlay))
    this%concs(igas)%conc(:,:) = w
    this%gas_name(igas) = trim(gas)
  end function set_vmr_2d
  ! -------------------------------------------------------------------------------------
  !
  ! Return volume mixing ratio as 1D or 2D array
  !
  ! -------------------------------------------------------------------------------------
  !
  ! 1D array ( lay depdendence only)
  !
  function get_vmr_1d(this, gas, array) result(error_msg)
    class(ty_gas_concs) :: this
    character(len=*),         intent(in ) :: gas
    real(wp), dimension(:),   intent(out) :: array
    character(len=128) :: error_msg
    ! ---------------------
    integer :: igas
    ! ---------------------
    error_msg = ''

    igas = this%find_gas(gas)
    if (igas == GAS_NOT_IN_LIST) then
      error_msg = 'ty_gas_concs%get_vmr; gas ' // trim(gas) // ' not found'
      array(:) = 0._wp
    else if(size(this%concs(igas)%conc, 1) > 1) then ! Are we requesting a single profile when many are present?
      error_msg = 'ty_gas_concs%get_vmr; gas ' // trim(gas) // ' requesting single profile but many are available'
      array(:) = 0._wp
    end if
    if(this%nlay > 0 .and. this%nlay /= size(array)) then
      error_msg = 'ty_gas_concs%get_vmr; gas ' // trim(gas) // ' array is wrong size (nlay)'
      array(:) = 0._wp
    end if
    if(error_msg /= "") return

    if(size(this%concs(igas)%conc, 2) > 1) then
      array(:) = this%concs(igas)%conc(1,:)
    else
      array(:) = this%concs(igas)%conc(1,1)
    end if

  end function get_vmr_1d
  ! -------------------------------------------------------------------------------------
  !
  ! 2D array (col, lay)
  !
  function get_vmr_2d(this, gas, array) result(error_msg)
    class(ty_gas_concs) :: this
    character(len=*),         intent(in ) :: gas
    real(wp), dimension(:,:), intent(out) :: array
    character(len=128)                    :: error_msg
    ! ---------------------
    integer :: igas
    ! ---------------------
    error_msg = ''

    igas = this%find_gas(gas)
    if (igas == GAS_NOT_IN_LIST) then
      error_msg = 'ty_gas_concs%get_vmr; gas ' // trim(gas) // ' not found'
      array(:,:) = 0._wp
    end if
    !
    ! Is the requested array the correct size?
    !
    if(this%ncol > 0 .and. this%ncol /= size(array,1)) then
      error_msg = 'ty_gas_concs%get_vmr; gas ' // trim(gas) // ' array is wrong size (ncol)'
      array(:,:) = 0._wp
    end if
    if(this%nlay > 0 .and. this%nlay /= size(array,2)) then
      error_msg = 'ty_gas_concs%get_vmr; gas ' // trim(gas) // ' array is wrong size (nlay)'
      array(:,:) = 0._wp
    end if
    if(error_msg /= "") return

    if(size(this%concs(igas)%conc, 1) > 1) then      ! Concentration stored as 2D
      array(:,:) =        this%concs(igas)%conc(:,:)
    else if(size(this%concs(igas)%conc, 2) > 1) then ! Concentration stored as 1D
      array(:,:) = spread(this%concs(igas)%conc(1,:), dim=1, ncopies=max(this%ncol, size(array,1)))
    else                                                   ! Concentration stored as scalar
      array(:,:) =        this%concs(igas)%conc(1,1)
    end if

  end function get_vmr_2d
  ! -------------------------------------------------------------------------------------
  !
  ! Extract a subset of n columns starting with column 'start'
  !
  ! -------------------------------------------------------------------------------------
  function get_subset_range(this, start, n, subset) result(error_msg)
    class(ty_gas_concs),      intent(in   ) :: this
    integer,                  intent(in   ) :: start, n
    class(ty_gas_concs),      intent(inout) :: subset
    character(len=128)                      :: error_msg
    ! ---------------------
    integer :: i
    ! ---------------------
    error_msg = ''
    if(n <= 0) &
       error_msg = "gas_concs%get_vmr: Asking for 0 or fewer columns "
    if(start < 1 ) &
       error_msg = "gas_concs%get_vmr: Asking for columns outside range"
    if(this%ncol > 0 .and. start > this%ncol .or. start+n-1 > this%ncol ) &
       error_msg = "gas_concs%get_vmr: Asking for columns outside range"
    if(error_msg /= "") return

    call subset%reset()
    allocate(subset%gas_name(size(this%gas_name)), &
             subset%concs   (size(this%concs))) ! These two arrays should be the same length
    subset%nlay = this%nlay
    subset%ncol = merge(n, 0, this%ncol > 0)
    subset%gas_name(:)  = this%gas_name(:)

    do i = 1, size(this%gas_name)
      !
      ! Preserve scalar/1D/2D representation in subset,
      !   but need to ensure at least extent 1 in col dimension (ncol = 0 means no gas exploits this dimension)
      !
      allocate(subset%concs(i)%conc(min(max(subset%ncol,1), size(this%concs(i)%conc, 1)), &
                                    min(    subset%nlay,    size(this%concs(i)%conc, 2))))
      if(size(this%concs(i)%conc, 1) > 1) then      ! Concentration stored as 2D
        subset%concs(i)%conc(:,:) = this%concs(i)%conc(start:(start+n-1),:)
      else
        subset%concs(i)%conc(:,:) = this%concs(i)%conc(:,:)
      end if
    end do

  end function get_subset_range
  ! -------------------------------------------------------------------------------------
  !
  ! Deallocate memory
  !
  ! -------------------------------------------------------------------------------------
  subroutine reset(this)
    class(ty_gas_concs), intent(inout) :: this
    ! -----------------
    integer :: i
    ! -----------------
    this%nlay = 0
    this%ncol = 0
    if(allocated(this%gas_name)) deallocate(this%gas_name)
    if (allocated(this%concs)) then
      do i = 1, size(this%concs)
        if(allocated(this%concs(i)%conc)) deallocate(this%concs(i)%conc)
      end do
      deallocate(this%concs)
    end if
  end subroutine reset
  ! -------------------------------------------------------------------------------------
  !
  ! Private procedures
  !
  ! -------------------------------------------------------------------------------------
  !
  ! This routine is called when adding a new concentration if the
  !   the gas isn't in the list already
  !
  subroutine increase_list_size(this)
    class(ty_gas_concs), intent(inout) :: this
    ! -----------------
    character(len=32), dimension(:), allocatable :: new_names
    type(conc_field),  dimension(:), allocatable :: new_concs
    ! -----------------

    if(allocated(this%gas_name)) then
      allocate(new_names(size(this%gas_name)+1), new_concs(size(this%gas_name)+1))
      new_names(1:size(this%gas_name)) = this%gas_name(:)
      new_concs(1:size(this%gas_name)) = this%concs(:)
      call move_alloc(new_names, this%gas_name)
      call move_alloc(new_concs, this%concs)
    else
      allocate(this%gas_name(1))
      allocate(this%concs(1))
    end if
  end subroutine increase_list_size
  ! -------------------------------------------------------------------------------------
  !
  ! find gas in list; GAS_NOT_IN_LIST if not found
  !
  function find_gas(this, gas)
    character(len=*),   intent(in) :: gas
    class(ty_gas_concs), intent(in) :: this
    integer                        :: find_gas
    ! -----------------
    integer :: igas
    ! -----------------
    find_gas = GAS_NOT_IN_LIST
    if(.not. allocated(this%gas_name)) return
    do igas = 1, size(this%gas_name)
      if (lower_case(trim(this%gas_name(igas))) == lower_case(trim(gas))) then
        find_gas = igas
      end if
    end do
  end function
  ! -------------------------------------------------------------------------------------
end module mo_gas_concentrations
