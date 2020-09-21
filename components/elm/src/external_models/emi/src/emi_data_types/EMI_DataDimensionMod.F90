module EMI_DataDimensionMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module defines the dimensions of EMI data
  !

  use decompMod , only : bounds_type
  use elm_varpar, only : nlevgrnd
  use elm_varpar, only : nlevsoi
  use elm_varpar, only : nlevsno
  use elm_varpar, only : ndecomp_pools
  use elm_varpar, only : nlevdecomp_full
  use abortutils, only : endrun
  use clm_varctl, only : iulog

implicit none

    character(*), parameter :: dimname_begg              = 'begg'
    character(*), parameter :: dimname_endg              = 'endg'
    character(*), parameter :: dimname_begl              = 'begl'
    character(*), parameter :: dimname_endl              = 'endl'
    character(*), parameter :: dimname_begc              = 'begc'
    character(*), parameter :: dimname_endc              = 'endc'
    character(*), parameter :: dimname_begp              = 'begp'
    character(*), parameter :: dimname_endp              = 'endp'
    character(*), parameter :: dimname_nlevsno           = '-nlevsno'
    character(*), parameter :: dimname_nlevsno_plus_one  = '-nlevsno + 1'
    character(*), parameter :: dimname_nlevsoi           = 'nlevsoi'
    character(*), parameter :: dimname_nlevgrnd          = 'nlevgrnd'
    character(*), parameter :: dimname_zero              = 'zero'
    character(*), parameter :: dimname_one               = 'one'
    character(*), parameter :: dimname_two               = 'two'
    character(*), parameter :: dimname_col_one_based_idx = 'endc - begc + 1'
    character(*), parameter :: dimname_nlevdecomp_full   = 'nlevdecomp_full';
    character(*), parameter :: dimname_ndecomp_pools     = 'ndecomp_pools';

  type emi_data_dimension_type
     character(len=24) :: name ! String labelling this IO type

     type(emi_data_dimension_type), pointer :: next
   contains
     procedure, public :: SetName      => EMID_Dim_SetName
     procedure, public :: GetDimValue  => EMID_Dim_GetDimValue
  end type emi_data_dimension_type


  type emi_data_dimension_list_type
     integer                                :: num_dims
     class(emi_data_dimension_type), pointer :: first
     class(emi_data_dimension_type), pointer :: last
  contains
     procedure, public :: Init             => EMID_Dim_List_Init
     procedure, public :: AddDim           => EMID_Dim_List_AddDim
     procedure, public :: AddDimByName     => EMID_Dim_List_AddDim_ByName
     procedure, public :: GetDimValue      => EMID_Dim_List_GetDimValue
  end type emi_data_dimension_list_type

contains

  !------------------------------------------------------------------------
  subroutine EMID_Dim_SetName(this, name)
    !
    ! !DESCRIPTION:
    ! Initializes a EMI data dimension
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(emi_data_dimension_type) :: this
    character(len=*)               :: name

    this%name = trim(name)

  end subroutine EMID_Dim_SetName

  !------------------------------------------------------------------------
  subroutine EMID_Dim_GetDimValue(this, bounds_clump, dim_name, dim_value)
    !
    ! !DESCRIPTION:
    ! Returns value for EMI data dimension
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(emi_data_dimension_type) :: this
    type(bounds_type), intent (in) :: bounds_clump
    character(len=*)               :: dim_name
    integer, intent(out)           :: dim_value


    select case (trim(dim_name))
    case (dimname_begg)
       dim_value = bounds_clump%begg

    case (dimname_endg)
       dim_value = bounds_clump%endg

    case (dimname_begl)
       dim_value = bounds_clump%begl

    case (dimname_endl)
       dim_value = bounds_clump%endl

    case (dimname_begc)
       dim_value = bounds_clump%begc

    case (dimname_endc)
       dim_value = bounds_clump%endc

    case (dimname_begp)
       dim_value = bounds_clump%begp

    case (dimname_endp)
       dim_value = bounds_clump%endp

    case (dimname_nlevsno)
       dim_value = -nlevsno

    case (dimname_nlevsno_plus_one)
       dim_value = -nlevsno + 1

    case (dimname_nlevsoi)
       dim_value = nlevsoi

    case (dimname_nlevgrnd)
       dim_value = nlevgrnd

    case (dimname_zero)
       dim_value = 0

    case (dimname_one)
       dim_value = 1

    case (dimname_two)
       dim_value = 2

    case (dimname_col_one_based_idx)
       dim_value = bounds_clump%endc - bounds_clump%begc + 1

    case (dimname_ndecomp_pools)
       dim_value = ndecomp_pools;

    case (dimname_nlevdecomp_full)
       dim_value = nlevdecomp_full;

    case default
       write(iulog,*)'dim_name = ',dim_name
       call endrun(msg='Unknown dim_name while trying to get dimension value.')
    end select

  end subroutine EMID_Dim_GetDimValue

  !------------------------------------------------------------------------
  subroutine EMID_Dim_List_Init(this)
    !
    ! !DESCRIPTION:
    ! Initializes a EMI data dimension list
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(emi_data_dimension_list_type) :: this
    !
    class(emi_data_dimension_type), pointer :: new_dim

    this%num_dims = 0

    nullify(this%first)
    nullify(this%last)

    call this%AddDimByName(dimname_begg)
    call this%AddDimByName(dimname_endg)
    call this%AddDimByName(dimname_begl)
    call this%AddDimByName(dimname_endl)
    call this%AddDimByName(dimname_begc)
    call this%AddDimByName(dimname_endc)
    call this%AddDimByName(dimname_begp)
    call this%AddDimByName(dimname_endp)
    call this%AddDimByName(dimname_nlevsno)
    call this%AddDimByName(dimname_nlevsno_plus_one)
    call this%AddDimByName(dimname_nlevsoi)
    call this%AddDimByName(dimname_nlevgrnd)
    call this%AddDimByName(dimname_zero)
    call this%AddDimByName(dimname_one)
    call this%AddDimByName(dimname_two)
    call this%AddDimByName(dimname_col_one_based_idx)
    call this%AddDimByName(dimname_nlevdecomp_full)
    call this%AddDimByName(dimname_ndecomp_pools)

  end subroutine EMID_Dim_List_Init

  !------------------------------------------------------------------------
  subroutine EMID_Dim_List_AddDim(this, new_dim)
    !
    ! !DESCRIPTION:
    ! Adds a EMID dimension to a list
    !
    ! !ARGUMENTS:
    implicit none
    !
    class(emi_data_dimension_list_type)     :: this
    class(emi_data_dimension_type), pointer :: new_dim

    this%num_dims = this%num_dims + 1

    if (.not.associated(this%first)) then
       this%first => new_dim
    endif

    if (associated(this%last)) then
       this%last%next => new_dim
    endif

    this%last => new_dim

  end subroutine EMID_Dim_List_AddDim

  !------------------------------------------------------------------------
  subroutine EMID_Dim_List_AddDim_ByName(this, dim_name)
    !
    ! !DESCRIPTION:
    ! Adds a EMID dimension to a list
    !
    ! !ARGUMENTS:
    implicit none
    !
    class(emi_data_dimension_list_type) :: this
    character(len=*)                    :: dim_name

    class(emi_data_dimension_type), pointer :: new_dim

    allocate(new_dim)
    call new_dim%SetName(dim_name)
    call this%AddDim(new_dim)
    nullify(new_dim)

  end subroutine EMID_Dim_List_AddDim_ByName

  !------------------------------------------------------------------------
  subroutine EMID_Dim_List_GetDimValue(this, bounds_clump, dim_name, dim_value)
    !
    ! !DESCRIPTION:
    ! Adds a EMID dimension to a list
    !
    ! !ARGUMENTS:
    implicit none
    !
    class(emi_data_dimension_list_type) :: this
    type(bounds_type), intent (in)      :: bounds_clump
    character(len=*)                    :: dim_name
    integer, intent(out)           :: dim_value

    class(emi_data_dimension_type), pointer :: cur_dim
    logical                                 :: dim_found

    dim_found = .false.

    cur_dim => this%first
    do
       if (.not.associated(cur_dim)) exit

       if (trim(cur_dim%name) .eq. dim_name) then
          call cur_dim%GetDimValue(bounds_clump, dim_name, dim_value)
          dim_found = .true.
          exit
       endif

       cur_dim => cur_dim%next
    enddo

    if (.not. dim_found) then
       write(iulog,*)'dim_name = ',dim_name
       call endrun(msg='Dimension not found in the list.')
    endif

  end subroutine EMID_Dim_List_GetDimValue

end module EMI_DataDimensionMod
