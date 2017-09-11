module ExternalModelInterfaceDataMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module provides
  !
  use shr_kind_mod         , only : r8 => shr_kind_r8
  use shr_log_mod          , only : errMsg => shr_log_errMsg
  use abortutils           , only : endrun
  use clm_varctl           , only : iulog
  !
  implicit none
  !
  type, public :: emi_data

     !private

     integer             :: id                  ! ID
     character (len=32)  :: name                ! Short name
     character (len=128) :: long_name           ! Long name of data
     character (len=24)  :: units               ! Units

     character (len=1)   :: avgflag             ! Needed for output to history file

     logical             :: is_int_type         ! Is data of an integer type?
     logical             :: is_real_type        ! Is data of real type?
     logical             :: is_set              ! Is data set

     integer             :: num_em_stages       ! Number of EM stages in which the data is exchanged
     integer, pointer    :: em_stage_ids(:)     ! ID of EM stages in which the data is exchanged

     integer             :: ndim                ! Dimension of the data

     ! Dimension name of the data (e.g. 'column','levscpf', etc).
     character (len=10) :: dim1_name
     character (len=10) :: dim2_name           
     character (len=10) :: dim3_name
     character (len=10) :: dim4_name

     character (len=24) :: dim1_beg_name
     character (len=24) :: dim2_beg_name
     character (len=24) :: dim3_beg_name
     character (len=24) :: dim4_beg_name
     character (len=24) :: dim1_end_name
     character (len=24) :: dim2_end_name
     character (len=24) :: dim3_end_name
     character (len=24) :: dim4_end_name

     integer :: dim1_beg, dim1_end
     integer :: dim2_beg, dim2_end
     integer :: dim3_beg, dim3_end
     integer :: dim4_beg, dim4_end

     integer, pointer :: data_int_1d(:)
     integer, pointer :: data_int_2d(:,:)
     integer, pointer :: data_int_3d(:,:,:)

     real(r8), pointer :: data_real_1d(:)
     real(r8), pointer :: data_real_2d(:,:)
     real(r8), pointer :: data_real_3d(:,:,:)
     real(r8), pointer :: data_real_4d(:,:,:,:)

     type(emi_data), pointer :: next

   contains

     procedure, public :: Init              => EMIDInit
     procedure, public :: Copy              => EMIDCopy
     procedure, public :: Setup             => EMIDSetup
     procedure, public :: SetDimensions     => EMIDSetDimensions
     procedure, public :: SetNDimension     => EMIDSetNDimension
     procedure, public :: SetType           => EMIDSetType
     procedure, public :: SetID             => EMIDSetID
     procedure, public :: SetName           => EMIDSetName
     procedure, public :: SetDimBegEndNames => EMIDSetDimBegEndNames
     procedure, public :: SetUnits          => EMIDSetUnits
     procedure, public :: SetAvgFlag        => EMIDSetAvgFlag
     procedure, public :: SetLongName       => EMIDSetLongName
     procedure, public :: SetEMStages       => EMIDSetEMStages
     procedure, public :: AppendEMStages    => EMIDAppendEMStages
     procedure, public :: AllocateMemory    => EMIDAllocateMemory
     procedure, public :: Reset             => EMIDReset
     procedure, public :: Destroy           => EMIDDestroy

  end type emi_data

  type emi_data_ptr
     type(emi_data), pointer :: data
  end type emi_data_ptr

  type, public :: emi_data_list

     integer                     :: num_data

     type(emi_data)    , pointer :: first
     type(emi_data)    , pointer :: last
     type(emi_data_ptr), pointer :: data_ptr(:)

   contains

     procedure, public :: Init               => EMIDListInit
     procedure, public :: AddData            => EMIDListAddData
     procedure, public :: AddDataByID        => EMIDListAddDataByID
     procedure, public :: Copy               => EMIDListCopy
     procedure, public :: Destroy            => EMIDListDestroy
     procedure, public :: GetIntValue        => EMIDListGetIntValue
     procedure, public :: GetPointerToInt1D  => EMIDListGetPointerToInt1D
     procedure, public :: GetPointerToInt2D  => EMIDListGetPointerToInt2D
     procedure, public :: GetPointerToInt3D  => EMIDListGetPointerToInt3D
     procedure, public :: GetPointerToReal1D => EMIDListGetPointerToReal1D
     procedure, public :: GetPointerToReal2D => EMIDListGetPointerToReal2D
     procedure, public :: GetPointerToReal3D => EMIDListGetPointerToReal3D
     procedure, public :: GetPointerToReal4D => EMIDListGetPointerToReal4D
     procedure, public :: IsDataIDPresent    => EMIDIsDataIDPresent
     procedure, public :: AppendDataEMStages => EMIDListAppendDataEMStages

  end type emi_data_list

contains

  !------------------------------------------------------------------------
  subroutine EMIDInit(this)
    !
    ! !DESCRIPTION:
    ! Initializes a EMI data
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(emi_data) :: this

    this%id             = 0
    this%name           = ""
    this%long_name      = ""
    this%units          = ""
    this%avgflag        = ""

    this%is_int_type    = .false.
    this%is_real_type   = .false.
    this%is_set         = .false.

    this%num_em_stages  = 0
    nullify(this%em_stage_ids)

    this%dim1_name      = ""
    this%dim2_name      = ""
    this%dim3_name      = ""
    this%dim4_name      = ""

    this%dim1_beg  = 0
    this%dim2_beg  = 0
    this%dim3_beg  = 0
    this%dim4_beg  = 0
    
    this%dim1_end  = 0
    this%dim2_end  = 0
    this%dim3_end  = 0
    this%dim4_end  = 0

    nullify(this%data_int_1d)
    nullify(this%data_int_2d)
    nullify(this%data_int_3d)

    nullify(this%data_real_1d)
    nullify(this%data_real_2d)
    nullify(this%data_real_3d)
    nullify(this%data_real_4d)

    nullify(this%next)
    
  end subroutine EMIDInit

  !------------------------------------------------------------------------
  subroutine EMIDCopy(this, default_data)
    !
    ! !DESCRIPTION:
    ! Initializes a EMI data
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(emi_data), intent(inout) :: this
    class(emi_data), intent(in) :: default_data
    !
    ! !LOCAL VARIABLES:
    integer :: ier                    ! error code

    call this%Init()

    this%id             = default_data%id
    this%name           = trim(default_data%name)
    this%long_name      = trim(default_data%long_name)
    this%units          = trim(default_data%units)
    this%avgflag        = trim(default_data%avgflag)

    this%is_int_type    = default_data%is_int_type
    this%is_real_type   = default_data%is_int_type

    this%num_em_stages  = default_data%num_em_stages

    if (associated(default_data%em_stage_ids)) then
       allocate(this%em_stage_ids(this%num_em_stages))
       this%em_stage_ids(:) = default_data%em_stage_ids(:)
    endif

    this%dim1_name      = trim(default_data%dim1_name)
    this%dim2_name      = trim(default_data%dim2_name)
    this%dim3_name      = trim(default_data%dim3_name)
    this%dim4_name      = trim(default_data%dim4_name)

    this%dim1_beg  = default_data%dim1_beg
    this%dim2_beg  = default_data%dim2_beg
    this%dim3_beg  = default_data%dim3_beg
    this%dim4_beg  = default_data%dim4_beg
    
    this%dim1_end  = default_data%dim1_end
    this%dim2_end  = default_data%dim2_end
    this%dim3_end  = default_data%dim3_end
    this%dim4_end  = default_data%dim4_end

    if (associated(default_data%data_int_1d)) then
       call EMIDAllocateMemory_Int_1D(this)
       this%data_int_1d(:) = default_data%data_int_1d(:)
    endif

    if (associated(default_data%data_int_2d)) then
       call EMIDAllocateMemory_Int_2D(this)
       this%data_int_2d(:,:) = default_data%data_int_2d(:,:)
    endif

    if (associated(default_data%data_int_3d)) then
       call EMIDAllocateMemory_Int_3D(this)
       this%data_int_3d(:,:,:) = default_data%data_int_3d(:,:,:)
    endif

    if (associated(default_data%data_real_1d)) then
       call EMIDAllocateMemory_Real_1D(this)
       this%data_real_1d(:) = default_data%data_real_1d(:)
    endif

    if (associated(default_data%data_real_2d)) then
       call EMIDAllocateMemory_Real_2D(this)
       this%data_real_2d(:,:) = default_data%data_real_2d(:,:)
    endif

    if (associated(default_data%data_real_3d)) then
       call EMIDAllocateMemory_Real_3D(this)
       this%data_real_3d(:,:,:) = default_data%data_real_3d(:,:,:)
    endif

    if (associated(default_data%data_real_4d)) then
       call EMIDAllocateMemory_Real_4D(this)
       this%data_real_4d(:,:,:,:) = default_data%data_real_4d(:,:,:,:)
    endif

    nullify(this%next)
    
  end subroutine EMIDCopy

  !------------------------------------------------------------------------
  subroutine EMIDSetup(this, id, name, long_name, units, avgflag, &
       num_em_stages, em_stage_ids)
    !
    ! !DESCRIPTION:
    ! Set value to data members of EMID object
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(emi_data)                      , intent(inout) :: this
    integer          , optional          , intent(in)    :: id
    character(len=*) , optional          , intent(in)    :: name
    character(len=*) , optional          , intent(in)    :: long_name
    character(len=*) , optional          , intent(in)    :: units
    character(len=*) , optional          , intent(in)    :: avgflag
    integer          , optional          , intent(in)    :: num_em_stages
    integer          , optional, pointer , intent(in)    :: em_stage_ids(:)

    if (present(id))        call this%SetID(id)
    if (present(name))      call this%SetName(name)
    if (present(long_name)) call this%SetLongName(name)
    if (present(units))     call this%SetUnits(units)
    if (present(avgflag))   call this%SetAvgFlag(units)

    if (present(num_em_stages) .or. present(em_stage_ids)) then
       if (present(num_em_stages) .and. present(em_stage_ids)) then
          call this%SetEMStages(num_em_stages, em_stage_ids)
       else
          call endrun(msg='EMIDSetup: Number of EM stages AND their IDs required.')
       endif
    endif

  end subroutine EMIDSetup

  !------------------------------------------------------------------------
  subroutine EMIDSetID(this, id)
    !
    ! !DESCRIPTION:
    ! Sets ID
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(emi_data) , intent(inout) :: this
    integer         , intent(in)    :: id
    
    this%id = id

  end subroutine EMIDSetID

  !------------------------------------------------------------------------
  subroutine EMIDSetName(this, name)
    !
    ! !DESCRIPTION:
    ! Sets short name
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(emi_data) , intent(inout) :: this
    character(len=*), intent(in)    :: name

    this%name = trim(name)

  end subroutine EMIDSetName

  !------------------------------------------------------------------------
  subroutine EMIDSetDimBegEndNames(this, d1_beg_name, d1_end_name, &
       d2_beg_name, d2_end_name, d3_beg_name, d3_end_name, &
       d4_beg_name, d4_end_name)
    !
    ! !DESCRIPTION:
    ! Sets beg/end dimension names
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(emi_data) , intent(inout) :: this
    character(len=*), intent(in)    :: d1_beg_name
    character(len=*), intent(in)    :: d2_beg_name
    character(len=*), intent(in)    :: d3_beg_name
    character(len=*), intent(in)    :: d4_beg_name
    character(len=*), intent(in)    :: d1_end_name
    character(len=*), intent(in)    :: d2_end_name
    character(len=*), intent(in)    :: d3_end_name
    character(len=*), intent(in)    :: d4_end_name

    this%dim1_beg_name = trim(d1_beg_name)
    this%dim1_end_name = trim(d1_end_name)
    this%dim2_beg_name = trim(d2_beg_name)
    this%dim2_end_name = trim(d2_end_name)
    this%dim3_beg_name = trim(d3_beg_name)
    this%dim3_end_name = trim(d3_end_name)
    this%dim4_beg_name = trim(d4_beg_name)
    this%dim4_end_name = trim(d4_end_name)

  end subroutine EMIDSetDimBegEndNames

  !------------------------------------------------------------------------
  subroutine EMIDSetLongName(this, long_name)
    !
    ! !DESCRIPTION:
    ! Sets long name
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(emi_data) , intent(inout) :: this
    character(len=*), intent(in)    :: long_name

    this%long_name = trim(long_name)

  end subroutine EMIDSetLongName

  !------------------------------------------------------------------------
  subroutine EMIDSetUnits(this, units)
    !
    ! !DESCRIPTION:
    ! Sets units
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(emi_data) , intent(inout) :: this
    character(len=*), intent(in)    :: units

    this%units = trim(units)

  end subroutine EMIDSetUnits

  !------------------------------------------------------------------------
  subroutine EMIDSetAvgFlag(this, avgflag)
    !
    ! !DESCRIPTION:
    ! Sets averaging flag option
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(emi_data) , intent(inout) :: this
    character(len=*), intent(in)    :: avgflag

    this%avgflag = trim(avgflag)

  end subroutine EMIDSetAvgFlag

  !------------------------------------------------------------------------
  subroutine EMIDSetEMStages(this, num_em_stages, em_stages)
    !
    ! !DESCRIPTION:
    ! Set number of external model IDs and stages in which the data would be
    ! exchanged.
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(emi_data), intent(inout) :: this
    integer, intent (in)           :: num_em_stages
    integer, intent (in), pointer  :: em_stages(:)
    
    if (size(em_stages) /= num_em_stages) then
       call endrun(msg='Number of EM stage IDs /= Number of EM Stages.')
    endif

    allocate(this%em_stage_ids(num_em_stages))

    this%num_em_stages   = num_em_stages
    this%em_stage_ids(:) = em_stages(:)

  end subroutine EMIDSetEMStages

  !------------------------------------------------------------------------
  subroutine EMIDAppendEMStages(this, num_em_stages, em_stages)
    !
    ! !DESCRIPTION:
    ! Appends the number of external model IDs and stages in which the data would be
    ! exchanged.
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(emi_data), intent(inout) :: this
    integer, intent (in)           :: num_em_stages
    integer, intent (in), pointer  :: em_stages(:)
    !
    integer                        :: num_em_stages_combined
    integer, pointer               :: em_stages_combined(:)
    integer                        :: num_unique_em_stages
    integer, pointer               :: unique_em_stages(:)
    integer                        :: iem
    integer                        :: ii, jj

    if (size(em_stages) /= num_em_stages) then
       call endrun(msg='Number of EM stage IDs /= Number of EM Stages.')
    endif

    num_em_stages_combined = num_em_stages + this%num_em_stages

    allocate(em_stages_combined(num_em_stages_combined))
    allocate(unique_em_stages  (num_em_stages_combined))

    do iem = 1,this%num_em_stages
       em_stages_combined(iem) = this%em_stage_ids(iem)
    enddo
    do iem = this%num_em_stages+1,num_em_stages_combined
       em_stages_combined(iem) = em_stages(iem-this%num_em_stages)
    enddo

    unique_em_stages(:)  = 1
    num_unique_em_stages = num_em_stages_combined

    do ii = 1, num_em_stages_combined
       do jj = ii+1, num_em_stages_combined
          if (em_stages_combined(ii) == em_stages_combined(jj)) then
             unique_em_stages(jj) = 0
             num_unique_em_stages = num_unique_em_stages - 1
          endif
       enddo
    enddo

    deallocate(this%em_stage_ids)
    allocate(this%em_stage_ids(num_unique_em_stages))

    this%num_em_stages   = num_unique_em_stages
    jj = 0
    do ii = 1, num_em_stages_combined
       if (unique_em_stages(ii) == 1) then
          jj = jj + 1
          this%em_stage_ids(jj) = em_stages_combined(ii)
       endif
    enddo

    deallocate(em_stages_combined)
    deallocate(unique_em_stages  )

  end subroutine EMIDAppendEMStages

  !------------------------------------------------------------------------
  subroutine EMIDSetDimensions(this,        &
    dim1_beg, dim1_end, dim2_beg, dim2_end, &
    dim3_beg, dim3_end, dim4_beg, dim4_end)
    !
    ! !DESCRIPTION:
    ! Sets dimenions for the data and allocates memory
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(emi_data)  , intent(inout) :: this
    integer          , intent (in)   :: dim1_beg, dim1_end
    integer          , intent (in)   :: dim2_beg, dim2_end
    integer          , intent (in)   :: dim3_beg, dim3_end
    integer          , intent (in)   :: dim4_beg, dim4_end

    this%dim1_beg = dim1_beg
    this%dim2_beg = dim2_beg
    this%dim3_beg = dim3_beg
    this%dim4_beg = dim4_beg

    this%dim1_end = dim1_end
    this%dim2_end = dim2_end
    this%dim3_end = dim3_end
    this%dim4_end = dim4_end

  end subroutine EMIDSetDimensions

  !------------------------------------------------------------------------
  subroutine EMIDSetNDimension(this, ndim)
    !
    ! !DESCRIPTION:
    ! Sets the dimension of data
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(emi_data)  , intent(inout) :: this
    integer          , intent (in)   :: ndim

    this%ndim = ndim

  end subroutine EMIDSetNDimension

  !------------------------------------------------------------------------
  subroutine EMIDSetType(this, is_int, is_real)
    !
    ! !DESCRIPTION:
    ! Sets data type
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(emi_data) , intent(inout) :: this
    logical         , intent(in)    :: is_int
    logical         , intent(in)    :: is_real

    if (is_int .and. is_real) then
       call endrun(msg='EMIDType must be either integer or real.')
    endif

    if ((.not.is_int) .and. (.not.is_real)) then
       call endrun(msg='EMIDType must be either integer or real.')
    endif

    this%is_int_type  = is_int
    this%is_real_type = is_real

  end subroutine EMIDSetType

  !------------------------------------------------------------------------
  subroutine EMIDAllocateMemory(this)
    !
    ! !DESCRIPTION:
    ! Initializes a EMI data
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(emi_data) :: this

    if (this%is_int_type .and. this%is_real_type) then
       call endrun(msg='Data type is defined to be both int and real.')
    endif

    if ((.not.this%is_int_type) .and. (.not.this%is_real_type)) then
       call endrun(msg='Data type is not defined to be either int or real.')
    endif

    if (this%ndim == 0) return

    select case(this%ndim)
    case (1)
       if (this%is_real_type) then
          call EMIDAllocateMemory_Real_1D(this)
       else
          call EMIDAllocateMemory_Int_1D(this)
       endif

    case (2)
       if (this%is_real_type) then
          call EMIDAllocateMemory_Real_2D(this)
       else
          call EMIDAllocateMemory_Int_2D(this)
       endif

    case (3)
       if (this%is_real_type) then
          call EMIDAllocateMemory_Real_3D(this)
       else
          call EMIDAllocateMemory_Int_3D(this)
       endif

    case (4)
       if (this%is_real_type) then
          call EMIDAllocateMemory_Real_4D(this)
       else
          call endrun(msg='EMID of type integer for dimension=4 is not supported.')

       endif

    case default
       call endrun(msg='EMID dimension larger than 4 is not supported.')

    end select

  end subroutine EMIDAllocateMemory

  !------------------------------------------------------------------------
  subroutine EMIDAllocateMemory_Int_1D(this)
    !
    ! !DESCRIPTION:
    ! Allocate 1D integer data
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(emi_data) :: this
    !
    ! !LOCAL VARIABLES:
    integer         :: ier ! error code

    allocate(this%data_int_1d(this%dim1_beg:this%dim1_end), &
             stat=ier)

    if (ier /= 0) then
       write(iulog,*) 'ERROR: Allocation failure'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

  end subroutine EMIDAllocateMemory_Int_1D

  !------------------------------------------------------------------------
  subroutine EMIDAllocateMemory_Int_2D(this)
    !
    ! !DESCRIPTION:
    ! Allocate 2D integer data type
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(emi_data) :: this
    !
    ! !LOCAL VARIABLES:
    integer         :: ier ! error code

    allocate(this%data_int_2d(this%dim1_beg:this%dim1_end, &
                              this%dim2_beg:this%dim2_end  ), &
             stat=ier)

    if (ier /= 0) then
       write(iulog,*) 'ERROR: Allocation failure'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

  end subroutine EMIDAllocateMemory_Int_2D

  !------------------------------------------------------------------------
  subroutine EMIDAllocateMemory_Int_3D(this)
    !
    ! !DESCRIPTION:
    ! Allocate 3D integer data type
    !
    implicit none
    !
    class(emi_data) :: this
    !
    ! !LOCAL VARIABLES:
    integer         :: ier ! error code

    allocate(this%data_int_3d( this%dim1_beg:this%dim1_end, &
                               this%dim2_beg:this%dim2_end, &
                               this%dim3_beg:this%dim3_end  ), &
             stat=ier)

    if (ier /= 0) then
       write(iulog,*) 'ERROR: Allocation failure'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

  end subroutine EMIDAllocateMemory_Int_3D

  !------------------------------------------------------------------------
  subroutine EMIDAllocateMemory_Real_1D(this)
    !
    ! !DESCRIPTION:
    ! Allocate 1D real data
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(emi_data) :: this
    !
    ! !LOCAL VARIABLES:
    integer         :: ier ! error code

    allocate(this%data_real_1d(this%dim1_beg:this%dim1_end), &
             stat=ier)

    if (ier /= 0) then
       write(iulog,*) 'ERROR: Allocation failure'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

  end subroutine EMIDAllocateMemory_Real_1D

  !------------------------------------------------------------------------
  subroutine EMIDAllocateMemory_Real_2D(this)
    !
    ! !DESCRIPTION:
    ! Allocate 2D real data type
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(emi_data) :: this
    !
    ! !LOCAL VARIABLES:
    integer         :: ier ! error code

    allocate(this%data_real_2d(this%dim1_beg:this%dim1_end, &
                               this%dim2_beg:this%dim2_end  ), &
             stat=ier)

    if (ier /= 0) then
       write(iulog,*) 'ERROR: Allocation failure'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

  end subroutine EMIDAllocateMemory_Real_2D

  !------------------------------------------------------------------------
  subroutine EMIDAllocateMemory_Real_3D(this)
    !
    ! !DESCRIPTION:
    ! Allocate 3D real data type
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(emi_data) :: this
    !
    ! !LOCAL VARIABLES:
    integer         :: ier ! error code

    allocate(this%data_real_3d(this%dim1_beg:this%dim1_end, &
                               this%dim2_beg:this%dim2_end, &
                               this%dim3_beg:this%dim3_end  ), &
             stat=ier)

    if (ier /= 0) then
       write(iulog,*) 'ERROR: Allocation failure'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

  end subroutine EMIDAllocateMemory_Real_3D

  !------------------------------------------------------------------------
  subroutine EMIDAllocateMemory_Real_4D(this)
    !
    ! !DESCRIPTION:
    ! Allocates 4D real data type
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(emi_data) :: this
    !
    ! !LOCAL VARIABLES:
    integer         :: ier ! error code

    allocate(this%data_real_4d(this%dim1_beg:this%dim1_end, &
                               this%dim2_beg:this%dim2_end, &
                               this%dim3_beg:this%dim3_end, &
                               this%dim4_beg:this%dim4_end  ), &
             stat=ier)

    if (ier /= 0) then
       write(iulog,*) 'ERROR: Allocation failure'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

  end subroutine EMIDAllocateMemory_Real_4D

  !------------------------------------------------------------------------
  subroutine EMIDReset(this)
    !
    ! !DESCRIPTION:
    ! Resets values of a EMI data
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(emi_data) :: this

    if (this%is_int_type .and. this%is_real_type) then
       call endrun(msg='Data type is defined to be both int and real.')
    endif

    if ((.not.this%is_int_type) .and. (.not.this%is_real_type)) then
       call endrun(msg='Data type is not defined to be either int or real.')
    endif

    if (this%ndim == 0) return

    select case(this%ndim)
    case (1)
       if (this%is_real_type) then
          this%data_real_1d(:) = 0._r8
       else
          this%data_int_1d(:) = 0
       endif

    case (2)
       if (this%is_real_type) then
          this%data_real_2d(:,:) = 0._r8
       else
          this%data_int_2d(:,:) = 0
       endif

    case (3)
       if (this%is_real_type) then
          this%data_real_3d(:,:,:) = 0._r8
       else
          this%data_int_3d(:,:,:) = 0
       endif

    case (4)
       if (this%is_real_type) then
          this%data_real_4d(:,:,:,:) = 0._r8
       else
          call endrun(msg='EMID of type integer for dimension=4 is not supported.')

       endif

    case default
       call endrun(msg='EMID dimension larger than 4 is not supported.')

    end select

    this%is_set = .false.

  end subroutine EMIDReset

  !------------------------------------------------------------------------
  subroutine EMIDDestroy(this)
    !
    ! !DESCRIPTION:
    ! Destroys a EMI data
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(emi_data) :: this

    if (this%is_int_type .and. this%is_real_type) then
       call endrun(msg='Data type is defined to be both int and real.')
    endif

    if ((.not.this%is_int_type) .and. (.not.this%is_real_type)) then
       call endrun(msg='Data type is not defined to be either int or real.')
    endif

    if (this%ndim == 0) return

    select case(this%ndim)
    case (1)
       if (this%is_real_type) then
          deallocate(this%data_real_1d)
       else
          deallocate(this%data_int_1d)
       endif

    case (2)
       if (this%is_real_type) then
          deallocate(this%data_real_2d)
       else
          deallocate(this%data_int_2d)
       endif

    case (3)
       if (this%is_real_type) then
          deallocate(this%data_real_3d)
       else
          deallocate(this%data_int_3d)
       endif

    case (4)
       if (this%is_real_type) then
          deallocate(this%data_real_4d)
       else
          call endrun(msg='EMID of type integer for dimension=4 is not supported.')
       endif

    case default
       call endrun(msg='EMID dimension larger than 4 is not supported.')

    end select

  end subroutine EMIDDestroy

  !------------------------------------------------------------------------
  subroutine EMIDListInit(this)
    !
    ! !DESCRIPTION:
    ! Initializes a EMID list
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(emi_data_list) :: this

    this%num_data = 0

    nullify(this%first)
    nullify(this%last)
    nullify(this%data_ptr)

  end subroutine EMIDListInit

  !------------------------------------------------------------------------
  subroutine EMIDListAddData(this, new_data)
    !
    ! !DESCRIPTION:
    ! Add a EMID to a list
    !
    ! !ARGUMENTS:
    implicit none
    !
    class(emi_data_list)     :: this
    class(emi_data), pointer :: new_data

    this%num_data = this%num_data + 1

    if (.not.associated(this%first)) then
       this%first => new_data
    endif
    
    if (associated(this%last)) then
       this%last%next => new_data
    endif

    this%last => new_data

  end subroutine EMIDListAddData
  
  !------------------------------------------------------------------------
  subroutine EMIDIsDataIDPresent(this, data_id, data_present)
    !
    ! !DESCRIPTION:
    ! Determine if the EMID, given by data_id, is present in the list
    !
    ! !ARGUMENTS:
    implicit none
    !
    class(emi_data_list)     :: this
    integer, intent(in)      :: data_id
    logical, intent(inout)   :: data_present
    !
    class(emi_data), pointer :: cur_data

    data_present = .false.

    cur_data => this%first
    do
       if (.not.associated(cur_data)) exit

       if (cur_data%id == data_id) then
          data_present = .true.
          exit
       endif

       cur_data => cur_data%next
    enddo

  end subroutine EMIDIsDataIDPresent

  !------------------------------------------------------------------------
  subroutine EMIDListAppendDataEMStages(this, data_id, num_em_stages_val, em_stage_ids_val, index_of_new_data)
    !
    ! !DESCRIPTION:
    ! Append EM stages of a data
    !
    ! !ARGUMENTS:
    implicit none
    !
    class(emi_data_list)          :: this
    integer, intent(in)           :: data_id
    integer, intent(in)           :: num_em_stages_val
    integer, pointer , intent(in) :: em_stage_ids_val(:)
    integer, intent(out)          :: index_of_new_data
    !
    class(emi_data), pointer      :: cur_data
    integer                       :: index_of_data

    index_of_data = 0

    cur_data => this%first
    do
       if (.not.associated(cur_data)) exit

       index_of_data = index_of_data + 1
       if (cur_data%id == data_id) then
          call cur_data%AppendEMStages(num_em_stages_val, em_stage_ids_val)
          index_of_new_data = index_of_data
          exit
       endif

       cur_data => cur_data%next
    enddo

  end subroutine EMIDListAppendDataEMStages

  !------------------------------------------------------------------------
  subroutine EMIDListAddDataByID(this, data_id, num_em_stages_val, em_stage_ids_val, index_of_new_data)
    !
    ! !DESCRIPTION:
    ! Add a EMID to a list
    !
    ! !USES:
    use ExternalModelConstants    , only : L2E_STATE_TSOIL_NLEVGRND
    use ExternalModelConstants    , only : L2E_STATE_H2OSOI_LIQ_NLEVGRND
    use ExternalModelConstants    , only : L2E_STATE_H2OSOI_ICE_NLEVGRND
    use ExternalModelConstants    , only : L2E_STATE_WTD
    use ExternalModelConstants    , only : L2E_STATE_VSFM_PROGNOSTIC_SOILP
    use ExternalModelConstants    , only : L2E_STATE_FRAC_H2OSFC
    use ExternalModelConstants    , only : L2E_STATE_FRAC_INUNDATED
    use ExternalModelConstants    , only : L2E_STATE_H2OSOI_LIQ_VOL_NLEVSOI
    use ExternalModelConstants    , only : L2E_STATE_H2OSOI_ICE_VOL_NLEVSOI
    use ExternalModelConstants    , only : L2E_STATE_H2OSOI_VOL_NLEVSOI
    use ExternalModelConstants    , only : L2E_STATE_AIR_VOL_NLEVSOI
    use ExternalModelConstants    , only : L2E_STATE_RHO_VAP_NLEVSOI
    use ExternalModelConstants    , only : L2E_STATE_RHVAP_SOI_NLEVSOI
    use ExternalModelConstants    , only : L2E_STATE_SOIL_MATRIC_POTENTIAL_NLEVSOI
    use ExternalModelConstants    , only : L2E_STATE_H2OSOI_LIQ_NLEVSOI
    use ExternalModelConstants    , only : L2E_STATE_H2OSOI_ICE_NLEVSOI
    use ExternalModelConstants    , only : L2E_STATE_TSNOW
    use ExternalModelConstants    , only : L2E_STATE_TH2OSFC
    use ExternalModelConstants    , only : L2E_STATE_H2OSOI_LIQ_NLEVSNOW
    use ExternalModelConstants    , only : L2E_STATE_H2OSOI_ICE_NLEVSNOW
    use ExternalModelConstants    , only : L2E_STATE_H2OSNOW
    use ExternalModelConstants    , only : L2E_STATE_H2OSFC
    use ExternalModelConstants    , only : L2E_STATE_FRAC_SNOW_EFFECTIVE

    use ExternalModelConstants    , only : E2L_STATE_H2OSOI_LIQ
    use ExternalModelConstants    , only : E2L_STATE_H2OSOI_ICE
    use ExternalModelConstants    , only : E2L_STATE_SOIL_MATRIC_POTENTIAL
    use ExternalModelConstants    , only : E2L_STATE_WTD
    use ExternalModelConstants    , only : E2L_STATE_VSFM_PROGNOSTIC_SOILP
    use ExternalModelConstants    , only : E2L_STATE_FSUN
    use ExternalModelConstants    , only : E2L_STATE_LAISUN
    use ExternalModelConstants    , only : E2L_STATE_LAISHA
    use ExternalModelConstants    , only : E2L_STATE_TSOIL_NLEVGRND
    use ExternalModelConstants    , only : E2L_STATE_TSNOW_NLEVSNOW
    use ExternalModelConstants    , only : E2L_STATE_TH2OSFC

    use ExternalModelConstants    , only : L2E_FLUX_INFIL_MASS_FLUX
    use ExternalModelConstants    , only : L2E_FLUX_VERTICAL_ET_MASS_FLUX
    use ExternalModelConstants    , only : L2E_FLUX_DEW_MASS_FLUX
    use ExternalModelConstants    , only : L2E_FLUX_SNOW_SUBLIMATION_MASS_FLUX
    use ExternalModelConstants    , only : L2E_FLUX_SNOW_LYR_DISAPPERANCE_MASS_FLUX
    use ExternalModelConstants    , only : L2E_FLUX_RESTART_SNOW_LYR_DISAPPERANCE_MASS_FLUX
    use ExternalModelConstants    , only : L2E_FLUX_DRAINAGE_MASS_FLUX
    use ExternalModelConstants    , only : L2E_FLUX_SOLAR_DIRECT_RADDIATION
    use ExternalModelConstants    , only : L2E_FLUX_SOLAR_DIFFUSE_RADDIATION
    use ExternalModelConstants    , only : L2E_FLUX_ABSORBED_SOLAR_RADIATION
    use ExternalModelConstants    , only : L2E_FLUX_SOIL_HEAT_FLUX
    use ExternalModelConstants    , only : L2E_FLUX_SNOW_HEAT_FLUX
    use ExternalModelConstants    , only : L2E_FLUX_H2OSFC_HEAT_FLUX
    use ExternalModelConstants    , only : L2E_FLUX_DERIVATIVE_OF_HEAT_FLUX

    use ExternalModelConstants    , only : E2L_FLUX_AQUIFER_RECHARGE
    use ExternalModelConstants    , only : E2L_FLUX_SNOW_LYR_DISAPPERANCE_MASS_FLUX

    use ExternalModelConstants    , only : L2E_FILTER_HYDROLOGYC
    use ExternalModelConstants    , only : L2E_FILTER_NUM_HYDROLOGYC
    use ExternalModelConstants    , only : L2E_FILTER_NOLAKEC
    use ExternalModelConstants    , only : L2E_FILTER_NUM_NOLAKEC
    use ExternalModelConstants    , only : L2E_FILTER_NOLAKEC_AND_NOURBANC
    use ExternalModelConstants    , only : L2E_FILTER_NUM_NOLAKEC_AND_NOURBANC

    use ExternalModelConstants    , only : L2E_COLUMN_ACTIVE
    use ExternalModelConstants    , only : L2E_COLUMN_TYPE
    use ExternalModelConstants    , only : L2E_COLUMN_LANDUNIT_INDEX
    use ExternalModelConstants    , only : L2E_COLUMN_ZI
    use ExternalModelConstants    , only : L2E_COLUMN_DZ
    use ExternalModelConstants    , only : L2E_COLUMN_Z
    use ExternalModelConstants    , only : L2E_COLUMN_AREA
    use ExternalModelConstants    , only : L2E_COLUMN_GRIDCELL_INDEX
    use ExternalModelConstants    , only : L2E_COLUMN_PATCH_INDEX
    use ExternalModelConstants    , only : L2E_COLUMN_NUM_SNOW_LAYERS
    use ExternalModelConstants    , only : L2E_COLUMN_ZI_SNOW_AND_SOIL
    use ExternalModelConstants    , only : L2E_COLUMN_DZ_SNOW_AND_SOIL
    use ExternalModelConstants    , only : L2E_COLUMN_Z_SNOW_AND_SOIL

    use ExternalModelConstants    , only : L2E_LANDUNIT_TYPE
    use ExternalModelConstants    , only : L2E_LANDUNIT_LAKEPOINT
    use ExternalModelConstants    , only : L2E_LANDUNIT_URBANPOINT

    use ExternalModelConstants    , only : L2E_PARAMETER_WATSATC
    use ExternalModelConstants    , only : L2E_PARAMETER_HKSATC
    use ExternalModelConstants    , only : L2E_PARAMETER_BSWC
    use ExternalModelConstants    , only : L2E_PARAMETER_SUCSATC
    use ExternalModelConstants    , only : L2E_PARAMETER_EFFPOROSITYC
    use ExternalModelConstants    , only : L2E_PARAMETER_CSOL
    use ExternalModelConstants    , only : L2E_PARAMETER_TKMG
    use ExternalModelConstants    , only : L2E_PARAMETER_TKDRY

    use ExternalModelIntefaceDataDimensionMod, only : dimname_begg
    use ExternalModelIntefaceDataDimensionMod, only : dimname_endg
    use ExternalModelIntefaceDataDimensionMod, only : dimname_begl
    use ExternalModelIntefaceDataDimensionMod, only : dimname_endl
    use ExternalModelIntefaceDataDimensionMod, only : dimname_begc
    use ExternalModelIntefaceDataDimensionMod, only : dimname_endc
    use ExternalModelIntefaceDataDimensionMod, only : dimname_begp
    use ExternalModelIntefaceDataDimensionMod, only : dimname_endp
    use ExternalModelIntefaceDataDimensionMod, only : dimname_nlevsno
    use ExternalModelIntefaceDataDimensionMod, only : dimname_nlevsno_plus_one
    use ExternalModelIntefaceDataDimensionMod, only : dimname_nlevsoi
    use ExternalModelIntefaceDataDimensionMod, only : dimname_nlevgrnd
    use ExternalModelIntefaceDataDimensionMod, only : dimname_zero
    use ExternalModelIntefaceDataDimensionMod, only : dimname_one
    use ExternalModelIntefaceDataDimensionMod, only : dimname_two
    use ExternalModelIntefaceDataDimensionMod, only : dimname_col_one_based_idx

    !
    implicit none
    !
    ! !ARGUMENTS:
    class(emi_data_list)                   :: this
    integer         , intent(in)           :: data_id
    integer         , intent(in)           :: num_em_stages_val
    integer         , pointer , intent(in) :: em_stage_ids_val(:)
    integer         , intent(out)          :: index_of_new_data
    !
    class(emi_data) , pointer              :: data
    integer                                :: id_val
    integer                                :: ndim
    character (len=32)                     :: name_val
    character (len=128)                    :: long_name_val
    character (len=24)                     :: units_val
    character (len=24)                     :: dim1_beg_name
    character (len=24)                     :: dim2_beg_name
    character (len=24)                     :: dim3_beg_name
    character (len=24)                     :: dim4_beg_name
    character (len=24)                     :: dim1_end_name
    character (len=24)                     :: dim2_end_name
    character (len=24)                     :: dim3_end_name
    character (len=24)                     :: dim4_end_name
    logical                                :: is_int_type
    logical                                :: is_real_type
    logical                                :: data_present

    call this%IsDataIDPresent(data_id, data_present)
    if (data_present) then
       call this%AppendDataEMStages(data_id, num_em_stages_val, &
                  em_stage_ids_val, index_of_new_data)
       return
    endif

    is_int_type    = .false.
    is_real_type   = .false.
    dim1_beg_name  = ''
    dim2_beg_name  = ''
    dim3_beg_name  = ''
    dim4_beg_name  = ''
    dim1_end_name  = ''
    dim2_end_name  = ''
    dim3_end_name  = ''
    dim4_end_name  = ''

    select case(data_id)

       ! --------------------------------------------------------------
       ! ALM-to-EM: State variables
       ! --------------------------------------------------------------
    case (L2E_STATE_TSOIL_NLEVGRND)
       id_val        = L2E_STATE_TSOIL_NLEVGRND
       name_val      = 'Soil temperature'
       long_name_val = 'Soil temperature: ALM to External Model'
       units_val     = '[K]'
       is_real_type  = .true.
       ndim          = 2
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc
       dim2_beg_name = dimname_one
       dim2_end_name = dimname_nlevgrnd

    case (L2E_STATE_H2OSOI_LIQ_NLEVGRND)
       id_val        = L2E_STATE_H2OSOI_LIQ_NLEVGRND
       name_val      = 'Soil liquid water'
       long_name_val = 'Soil liquid water: ALM to External Model'
       units_val     = '[kg/m2]'
       is_real_type  = .true.
       ndim          = 2
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc
       dim2_beg_name = dimname_one
       dim2_end_name = dimname_nlevgrnd

    case (L2E_STATE_H2OSOI_ICE_NLEVGRND)
       id_val        = L2E_STATE_H2OSOI_ICE_NLEVGRND
       name_val      = 'Soil ice water'
       long_name_val = 'Soil ice water: ALM to External Model'
       units_val     = '[kg/m2]'
       is_real_type  = .true.
       ndim          = 2
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc
       dim2_beg_name = dimname_one
       dim2_end_name = dimname_nlevgrnd

    case (L2E_STATE_WTD)
       id_val        = L2E_STATE_WTD
       name_val      = 'Water table depth'
       long_name_val = 'Water table depth: ALM to External Model'
       units_val     = '[m]'
       is_real_type  = .true.
       ndim          = 1
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc

    case (L2E_STATE_VSFM_PROGNOSTIC_SOILP)
       id_val        = L2E_STATE_VSFM_PROGNOSTIC_SOILP
       name_val      = 'Soil liquid pressure'
       long_name_val = 'Soil liquid pressure: ALM to External Model'
       units_val     = '[Pa]'
       is_real_type  = .true.
       ndim          = 2
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc
       dim2_beg_name = dimname_one
       dim2_end_name = dimname_nlevgrnd

    case (L2E_STATE_FRAC_H2OSFC)
       id_val        = L2E_STATE_FRAC_H2OSFC
       name_val      = 'Fraction of standing water'
       long_name_val = 'Fraction of standing water: ALM to External Model'
       units_val     = '[-]'
       is_real_type  = .true.
       ndim          = 1
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc

    case (L2E_STATE_FRAC_INUNDATED)
       id_val        = L2E_STATE_FRAC_INUNDATED
       name_val      = 'Fraction of inundated surface'
       long_name_val = 'Fraction of inundated surface: ALM to External Model'
       units_val     = '[-]'
       is_real_type  = .true.
       ndim          = 1
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc

    case (L2E_STATE_H2OSOI_LIQ_VOL_NLEVSOI)
       id_val        = L2E_STATE_H2OSOI_LIQ_VOL_NLEVSOI
       name_val      = 'Volumetric liquid content'
       long_name_val = 'Volumetric liquid content: ALM to External Model'
       units_val     = '[m3/m3]'
       is_real_type  = .true.
       ndim          = 2
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc
       dim2_beg_name = dimname_one
       dim2_end_name = dimname_nlevsoi

    case (L2E_STATE_H2OSOI_ICE_VOL_NLEVSOI)
       id_val        = L2E_STATE_H2OSOI_ICE_VOL_NLEVSOI
       name_val      = 'Volumetric ice content'
       long_name_val = 'Volumetric ice content: ALM to External Model'
       units_val     = '[m3/m3]'
       is_real_type  = .true.
       ndim          = 2
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc
       dim2_beg_name = dimname_one
       dim2_end_name = dimname_nlevsoi

    case (L2E_STATE_H2OSOI_VOL_NLEVSOI)
       id_val        = L2E_STATE_H2OSOI_VOL_NLEVSOI
       name_val      = 'Volumetric soil water'
       long_name_val = 'Volumetric soil water: ALM to External Model'
       units_val     = '[m3/m3]'
       is_real_type  = .true.
       ndim          = 2
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc
       dim2_beg_name = dimname_one
       dim2_end_name = dimname_nlevsoi

    case (L2E_STATE_AIR_VOL_NLEVSOI)
       id_val        = L2E_STATE_AIR_VOL_NLEVSOI
       name_val      = 'Air filled porosity'
       long_name_val = 'Air filled porosity: ALM to External Model'
       units_val     = '[m3/m3]'
       is_real_type  = .true.
       ndim          = 2
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc
       dim2_beg_name = dimname_one
       dim2_end_name = dimname_nlevsoi

    case (L2E_STATE_RHO_VAP_NLEVSOI)
       id_val        = L2E_STATE_RHO_VAP_NLEVSOI
       name_val      = 'Water vapor pressure'
       long_name_val = 'Water vapor pressure: ALM to External Model'
       units_val     = '[Pa]'
       is_real_type  = .true.
       ndim          = 2
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc
       dim2_beg_name = dimname_one
       dim2_end_name = dimname_nlevsoi

    case (L2E_STATE_RHVAP_SOI_NLEVSOI)
       id_val        = L2E_STATE_RHVAP_SOI_NLEVSOI
       name_val      = 'Relative humidity'
       long_name_val = 'Relative humidity: ALM to External Model'
       units_val     = '[-]'
       is_real_type  = .true.
       ndim          = 2
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc
       dim2_beg_name = dimname_one
       dim2_end_name = dimname_nlevsoi

    case (L2E_STATE_SOIL_MATRIC_POTENTIAL_NLEVSOI)
       id_val        = L2E_STATE_SOIL_MATRIC_POTENTIAL_NLEVSOI
       name_val      = 'Matric potential in soil & snow lyrs'
       long_name_val = 'Matric potential in soil & snow lyrs: ALM to External Model'
       units_val     = '[mm]'
       is_real_type  = .true.
       ndim          = 2
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc
       dim2_beg_name = dimname_one
       dim2_end_name = dimname_nlevsoi

    case (L2E_STATE_H2OSOI_LIQ_NLEVSOI)
       id_val        = L2E_STATE_H2OSOI_LIQ_NLEVSOI
       name_val      = 'Soil liquid water'
       long_name_val = 'Soil liquid water: ALM to External Model'
       units_val     = '[kg/m2]'
       is_real_type  = .true.
       ndim          = 2
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc
       dim2_beg_name = dimname_one
       dim2_end_name = dimname_nlevsoi

    case (L2E_STATE_H2OSOI_ICE_NLEVSOI)
       id_val        = L2E_STATE_H2OSOI_ICE_NLEVSOI
       name_val      = 'Soil ice water'
       long_name_val = 'Soil ice water: ALM to External Model'
       units_val     = '[kg/m2]'
       is_real_type  = .true.
       ndim          = 2
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc
       dim2_beg_name = dimname_one
       dim2_end_name = dimname_nlevsoi

    case (L2E_STATE_TSNOW)
       id_val        = L2E_STATE_TSNOW
       name_val      = 'Temperature snow'
       long_name_val = 'Temperature snow: ALM to External Model'
       units_val     = '[K]'
       is_real_type  = .true.
       ndim          = 2
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc
       dim2_beg_name = dimname_nlevsno_plus_one
       dim2_end_name = dimname_zero

    case (L2E_STATE_TH2OSFC)
       id_val        = L2E_STATE_TH2OSFC
       name_val      = 'Temperature of standing water'
       long_name_val = 'Temperature of standing water: ALM to External Model'
       units_val     = '[K]'
       is_real_type  = .true.
       ndim          = 1
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc

    case (L2E_STATE_H2OSOI_LIQ_NLEVSNOW)
       id_val        = L2E_STATE_H2OSOI_LIQ_NLEVSNOW
       name_val      = 'Liquid water in snow'
       long_name_val = 'Liquid water in snow: ALM to External Model'
       units_val     = '[kg/m2]'
       is_real_type  = .true.
       ndim          = 2
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc
       dim2_beg_name = dimname_nlevsno_plus_one
       dim2_end_name = dimname_zero

    case (L2E_STATE_H2OSOI_ICE_NLEVSNOW)
       id_val        = L2E_STATE_H2OSOI_ICE_NLEVSNOW
       name_val      = 'Ice water in snow'
       long_name_val = 'Ice water in snow: ALM to External Model'
       units_val     = '[kg/m2]'
       is_real_type  = .true.
       ndim          = 2
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc
       dim2_beg_name = dimname_nlevsno_plus_one
       dim2_end_name = dimname_zero

    case (L2E_STATE_H2OSNOW)
       id_val        = L2E_STATE_H2OSNOW
       name_val      = 'Snow water'
       long_name_val = 'Snow water: ALM to External Model'
       units_val     = '[mm H2O]'
       is_real_type  = .true.
       ndim          = 1
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc

    case (L2E_STATE_H2OSFC)
       id_val        = L2E_STATE_H2OSFC
       name_val      = 'Standing water'
       long_name_val = 'Standing water: ALM to External Model'
       units_val     = '[mm H2O]'
       is_real_type  = .true.
       ndim          = 1
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc

    case (L2E_STATE_FRAC_SNOW_EFFECTIVE)
       id_val        = L2E_STATE_FRAC_SNOW_EFFECTIVE
       name_val      = 'Frac. of grnd covered with snow'
       long_name_val = 'Frac. of ground covered with snow: ALM to External Model'
       units_val     = '[-]'
       is_real_type  = .true.
       ndim          = 1
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc

       ! --------------------------------------------------------------
       ! EM-to-ALM: State variables
       ! --------------------------------------------------------------
    case (E2L_STATE_H2OSOI_LIQ)
       id_val        = E2L_STATE_H2OSOI_LIQ
       name_val      = 'Soil liquid water'
       long_name_val = 'Soil liquid water: External Model to ALM'
       units_val     = '[kg/m2]'
       is_real_type  = .true.
       ndim          = 2
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc
       dim2_beg_name = dimname_one
       dim2_end_name = dimname_nlevgrnd

    case (E2L_STATE_H2OSOI_ICE)
       id_val        = E2L_STATE_H2OSOI_ICE
       name_val      = 'Soil ice water'
       long_name_val = 'Soil ice water: External Model to ALM'
       units_val     = '[kg/m2]'
       is_real_type  = .true.
       ndim          = 2
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc
       dim2_beg_name = dimname_one
       dim2_end_name = dimname_nlevgrnd

    case (E2L_STATE_SOIL_MATRIC_POTENTIAL)
       id_val        = E2L_STATE_SOIL_MATRIC_POTENTIAL
       name_val      = 'Soil matric potential'
       long_name_val = 'Soil matric potential: External Model to ALM'
       units_val     = '[mm]'
       is_real_type  = .true.
       ndim          = 2
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc
       dim2_beg_name = dimname_one
       dim2_end_name = dimname_nlevgrnd

    case (E2L_STATE_WTD)
       id_val        = E2L_STATE_WTD
       name_val      = 'Water table depth'
       long_name_val = 'Water table depth: External Model to ALM'
       units_val     = '[m]'
       is_real_type  = .true.
       ndim          = 1
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc

    case (E2L_STATE_VSFM_PROGNOSTIC_SOILP)
       id_val        = E2L_STATE_VSFM_PROGNOSTIC_SOILP
       name_val      = 'Soil liquid pressure'
       long_name_val = 'Soil liquid pressure: External Model to ALM'
       units_val     = '[Pa]'
       is_real_type  = .true.
       ndim          = 2
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc
       dim2_beg_name = dimname_one
       dim2_end_name = dimname_nlevgrnd

    case (E2L_STATE_FSUN)
       id_val        = E2L_STATE_FSUN
       name_val      = 'Fraction of canopy sunlit'
       long_name_val = ': External Model to ALM'
       units_val     = '[-]'
       is_real_type  = .true.
       ndim          = 1
       dim1_beg_name = dimname_begp
       dim1_end_name = dimname_endp

    case (E2L_STATE_LAISUN)
       id_val        = E2L_STATE_LAISUN
       name_val      = 'Sunlit leaf area'
       long_name_val = 'Sunlit leaf area: External Model to ALM'
       units_val     = '[-]'
       is_real_type  = .true.
       ndim          = 1
       dim1_beg_name = dimname_begp
       dim1_end_name = dimname_endp

    case (E2L_STATE_LAISHA)
       id_val        = E2L_STATE_LAISHA
       name_val      = 'Shaded leaf area'
       long_name_val = 'Shaded leaf area: External Model to ALM'
       units_val     = '[-]'
       is_real_type  = .true.
       ndim          = 1
       dim1_beg_name = dimname_begp
       dim1_end_name = dimname_endp

    case (E2L_STATE_TSOIL_NLEVGRND)
       id_val        = E2L_STATE_TSOIL_NLEVGRND
       name_val      = 'Soil temperature'
       long_name_val = 'Soil temperature: External Model to ALM'
       units_val     = '[K]'
       is_real_type  = .true.
       ndim          = 2
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc
       dim2_beg_name = dimname_one
       dim2_end_name = dimname_nlevgrnd

    case (E2L_STATE_TSNOW_NLEVSNOW)
       id_val        = E2L_STATE_TSNOW_NLEVSNOW
       name_val      = 'Snow temperature'
       long_name_val = 'Snow temperature: External Model to ALM'
       units_val     = '[K]'
       is_real_type  = .true.
       ndim          = 2
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc
       dim2_beg_name = dimname_nlevsno_plus_one
       dim2_end_name = dimname_zero

    case (E2L_STATE_TH2OSFC)
       id_val        = E2L_STATE_TH2OSFC
       name_val      = 'Standing water temperature'
       long_name_val = 'Standing water temperature: External Model to ALM'
       units_val     = '[K]'
       is_real_type  = .true.
       ndim          = 1
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc

       ! --------------------------------------------------------------
       ! ALM-to-EM: Flux variables
       ! --------------------------------------------------------------
    case (L2E_FLUX_INFIL_MASS_FLUX)
       id_val        = L2E_FLUX_INFIL_MASS_FLUX
       name_val      = 'Soil infiltration source'
       long_name_val = 'Soil infiltration source: ALM to External Model'
       units_val     = '[kg/s]'
       is_real_type  = .true.
       ndim          = 1
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc

    case (L2E_FLUX_VERTICAL_ET_MASS_FLUX)
       id_val        = L2E_FLUX_VERTICAL_ET_MASS_FLUX
       name_val      = 'Evapotranspiration sink'
       long_name_val = 'Evapotranspiration sink: ALM to External Model'
       units_val     = '[kg/s]'
       is_real_type  = .true.
       ndim          = 2
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc
       dim2_beg_name = dimname_one
       dim2_end_name = dimname_nlevgrnd

    case (L2E_FLUX_DEW_MASS_FLUX)
       id_val        = L2E_FLUX_DEW_MASS_FLUX
       name_val      = 'Dew sink'
       long_name_val = 'Dew sink: ALM to External Model'
       units_val     = '[kg/s]'
       is_real_type  = .true.
       ndim          = 1
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc

    case (L2E_FLUX_SNOW_SUBLIMATION_MASS_FLUX)
       id_val        = L2E_FLUX_SNOW_SUBLIMATION_MASS_FLUX
       name_val      = 'Snow sublimation sink'
       long_name_val = 'Snow sublimation sink: ALM to External Model'
       units_val     = '[kg/s]'
       is_real_type  = .true.
       ndim          = 1
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc

    case (L2E_FLUX_SNOW_LYR_DISAPPERANCE_MASS_FLUX)
       id_val        = L2E_FLUX_SNOW_LYR_DISAPPERANCE_MASS_FLUX
       name_val      = 'Snow layer disappearance sink'
       long_name_val = 'Snow layer disappearance sink: ALM to External Model'
       units_val     = '[kg/s]'
       is_real_type  = .true.
       ndim          = 1
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc

    case (L2E_FLUX_RESTART_SNOW_LYR_DISAPPERANCE_MASS_FLUX)
       id_val        = L2E_FLUX_RESTART_SNOW_LYR_DISAPPERANCE_MASS_FLUX
       name_val      = 'Snow layer disappearance sink'
       long_name_val = 'Snow layer disappearance sink from restart: ALM to External Model'
       units_val     = '[kg/s]'
       is_real_type  = .true.
       ndim          = 1
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc

    case (L2E_FLUX_DRAINAGE_MASS_FLUX)
       id_val        = L2E_FLUX_DRAINAGE_MASS_FLUX
       name_val      = 'Drainage sink'
       long_name_val = 'Drainage sink: ALM to External Model'
       units_val     = '[kg/s]'
       is_real_type  = .true.
       ndim          = 2
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc
       dim2_beg_name = dimname_one
       dim2_end_name = dimname_nlevgrnd

    case (L2E_FLUX_SOLAR_DIRECT_RADDIATION)
       id_val        = L2E_FLUX_SOLAR_DIRECT_RADDIATION
       name_val      = 'Direct beam solar radiation'
       long_name_val = 'Direct beam solar radiation: ALM to External Model'
       units_val     = '[W/m2]'
       is_real_type  = .true.
       ndim          = 2
       dim1_beg_name = dimname_begg
       dim1_end_name = dimname_endg
       dim2_beg_name = dimname_one
       dim2_end_name = dimname_two

    case (L2E_FLUX_SOLAR_DIFFUSE_RADDIATION)
       id_val        = L2E_FLUX_SOLAR_DIFFUSE_RADDIATION
       name_val      = 'Diffuse beam solar radiation'
       long_name_val = 'Diffuse beam solar radiation: ALM to External Model'
       units_val     = '[W/m2]'
       is_real_type  = .true.
       ndim          = 2
       dim1_beg_name = dimname_begg
       dim1_end_name = dimname_endg
       dim2_beg_name = dimname_one
       dim2_end_name = dimname_two

    case (L2E_FLUX_ABSORBED_SOLAR_RADIATION)
       id_val        = L2E_FLUX_ABSORBED_SOLAR_RADIATION
       name_val      = 'Absorbed solar radiation'
       long_name_val = 'Absorbed beam solar radiation: ALM to External Model'
       units_val     = '[W/m2]'
       is_real_type  = .true.
       ndim          = 2
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc
       dim2_beg_name = dimname_nlevsno_plus_one
       dim2_end_name = dimname_one

    case (L2E_FLUX_SOIL_HEAT_FLUX)
       id_val        = L2E_FLUX_SOIL_HEAT_FLUX
       name_val      = 'Net heat flux into soil'
       long_name_val = 'Net heat flux into soil: ALM to External Model'
       units_val     = '[W/m2]'
       is_real_type  = .true.
       ndim          = 1
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc

    case (L2E_FLUX_SNOW_HEAT_FLUX)
       id_val        = L2E_FLUX_SNOW_HEAT_FLUX
       name_val      = 'Net heat flux into snow'
       long_name_val = 'Net heat flux into snow: ALM to External Model'
       units_val     = '[W/m2]'
       is_real_type  = .true.
       ndim          = 1
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc

    case (L2E_FLUX_H2OSFC_HEAT_FLUX)
       id_val        = L2E_FLUX_H2OSFC_HEAT_FLUX
       name_val      = 'Net heat flux into standing water'
       long_name_val = 'Net heat flux into standing water: ALM to External Model'
       units_val     = '[W/m2]'
       is_real_type  = .true.
       ndim          = 1
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc

    case (L2E_FLUX_DERIVATIVE_OF_HEAT_FLUX)
       id_val        = L2E_FLUX_DERIVATIVE_OF_HEAT_FLUX
       name_val      = 'Derivative of heat flux w.r.t. temperature'
       long_name_val = 'Derivative of heat flux w.r.t. temperature: ALM to External Model'
       units_val     = '[W/m2]'
       is_real_type  = .true.
       ndim          = 1
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc

       ! --------------------------------------------------------------
       ! EM-to-ALM: Flux variables
       ! --------------------------------------------------------------
    case (E2L_FLUX_AQUIFER_RECHARGE)
       id_val        = E2L_FLUX_AQUIFER_RECHARGE
       name_val      = 'Aquifer recharge rate'
       long_name_val = 'Aquifer recharge rate: External Model to ALM'
       units_val     = '[mm/s]'
       is_real_type  = .true.
       ndim          = 1
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc

    case (E2L_FLUX_SNOW_LYR_DISAPPERANCE_MASS_FLUX)
       id_val        = E2L_FLUX_SNOW_LYR_DISAPPERANCE_MASS_FLUX
       name_val      = 'Snow layer disappearance sink'
       long_name_val = 'Snow layer disappearance sink initial value: External Model to ALM'
       units_val     = '[kg/s]'
       is_real_type  = .true.
       ndim          = 1
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc

       ! --------------------------------------------------------------
       ! ALM-to-ELM: Filter variables
       ! --------------------------------------------------------------
    case (L2E_FILTER_HYDROLOGYC)
       id_val        = L2E_FILTER_HYDROLOGYC
       name_val      = 'Hydrology filter'
       long_name_val = 'Hydrology filter: ALM to External Model'
       units_val     = '[-]'
       is_int_type   = .true.
       ndim          = 1
       dim1_beg_name = dimname_one
       dim1_end_name = dimname_col_one_based_idx

    case (L2E_FILTER_NUM_HYDROLOGYC)
       id_val        = L2E_FILTER_NUM_HYDROLOGYC
       name_val      = 'Number of hydrology filter'
       long_name_val = 'Number of hydrology filter: ALM to External Model'
       units_val     = '[-]'
       is_int_type   = .true.
       ndim          = 1
       dim1_beg_name = dimname_one
       dim1_end_name = dimname_one

    case (L2E_FILTER_NOLAKEC)
       id_val        = L2E_FILTER_HYDROLOGYC
       name_val      = 'Non-lake filter'
       long_name_val = 'Non-lake filter: ALM to External Model'
       units_val     = '[-]'
       is_int_type   = .true.
       ndim          = 1
       dim1_beg_name = dimname_one
       dim1_end_name = dimname_col_one_based_idx

    case (L2E_FILTER_NUM_NOLAKEC)
       id_val        = L2E_FILTER_NUM_HYDROLOGYC
       name_val      = 'Number of non-lake filter'
       long_name_val = 'Number of non-lake filter: ALM to External Model'
       units_val     = '[-]'
       is_int_type   = .true.
       ndim          = 1
       dim1_beg_name = dimname_one
       dim1_end_name = dimname_one

    case (L2E_FILTER_NOLAKEC_AND_NOURBANC)
       id_val        = L2E_FILTER_NOLAKEC_AND_NOURBANC
       name_val      = 'Non-lake & non-urban filter'
       long_name_val = 'Non-lake & non-urban filter: ALM to External Model'
       units_val     = '[-]'
       is_int_type   = .true.
       ndim          = 1
       dim1_beg_name = dimname_one
       dim1_end_name = dimname_col_one_based_idx

    case (L2E_FILTER_NUM_NOLAKEC_AND_NOURBANC)
       id_val        = L2E_FILTER_NUM_NOLAKEC_AND_NOURBANC
       name_val      = 'Number of non-lake & non-urban filter'
       long_name_val = 'Number of non-lake & non-urban filter: ALM to External Model'
       units_val     = '[-]'
       is_int_type   = .true.
       ndim          = 1
       dim1_beg_name = dimname_one
       dim1_end_name = dimname_one

       ! --------------------------------------------------------------
       ! ALM-to-ELM: Column variables
       ! --------------------------------------------------------------
    case (L2E_COLUMN_ACTIVE)
       id_val        = L2E_COLUMN_ACTIVE
       name_val      = 'Column active'
       long_name_val = 'Column active: ALM to External Model'
       units_val     = '[-]'
       is_int_type   = .true.
       ndim          = 1
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc

    case (L2E_COLUMN_TYPE)
       id_val        = L2E_COLUMN_TYPE
       name_val      = 'Column type'
       long_name_val = 'Column type: ALM to External Model'
       units_val     = '[-]'
       is_int_type   = .true.
       ndim          = 1
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc

    case (L2E_COLUMN_LANDUNIT_INDEX)
       id_val        = L2E_COLUMN_LANDUNIT_INDEX
       name_val      = 'Column to landunit index'
       long_name_val = 'Column landunit index: ALM to External Model'
       units_val     = '[-]'
       is_int_type   = .true.
       ndim          = 1
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc

    case (L2E_COLUMN_ZI)
       id_val        = L2E_COLUMN_ZI
       name_val      = 'Column layer interface depth'
       long_name_val = 'Column layer interface depth: ALM to External Model'
       units_val     = '[m]'
       is_real_type  = .true.
       ndim          = 2
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc
       dim2_beg_name = dimname_zero
       dim2_end_name = dimname_nlevgrnd

    case (L2E_COLUMN_DZ)
       id_val        = L2E_COLUMN_DZ
       name_val      = 'Column layer thickness'
       long_name_val = 'Column layer thickness: ALM to External Model'
       units_val     = '[m]'
       is_real_type  = .true.
       ndim          = 2
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc
       dim2_beg_name = dimname_one
       dim2_end_name = dimname_nlevgrnd

    case (L2E_COLUMN_Z)
       id_val        = L2E_COLUMN_Z
       name_val      = 'Column layer centroid depth'
       long_name_val = 'Column layer centroid depth: ALM to External Model'
       units_val     = '[m]'
       is_real_type  = .true.
       ndim          = 2
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc
       dim2_beg_name = dimname_one
       dim2_end_name = dimname_nlevgrnd

    case (L2E_COLUMN_AREA)
       id_val        = L2E_COLUMN_AREA
       name_val      = 'Column surface area'
       long_name_val = 'Column surface area: ALM to External Model'
       units_val     = '[m2]'
       is_real_type  = .true.
       ndim          = 1
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc

    case (L2E_COLUMN_GRIDCELL_INDEX)
       id_val        = L2E_COLUMN_GRIDCELL_INDEX
       name_val      = 'Column to gridcell index'
       long_name_val = 'Column to gridcell index: ALM to External Model'
       units_val     = '[-]'
       is_int_type   = .true.
       ndim          = 1
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc

    case (L2E_COLUMN_PATCH_INDEX)
       id_val        = L2E_COLUMN_PATCH_INDEX
       name_val      = 'Column to patch index'
       long_name_val = 'Column to patch index: ALM to External Model'
       units_val     = '[-]'
       is_int_type   = .true.
       ndim          = 1
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc

    case (L2E_COLUMN_NUM_SNOW_LAYERS)
       id_val        = L2E_COLUMN_NUM_SNOW_LAYERS
       name_val      = 'Number of snow layers'
       long_name_val = 'Number of snow layers: ALM to External Model'
       units_val     = '[-]'
       is_int_type   = .true.
       ndim          = 1
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc

    case (L2E_COLUMN_ZI_SNOW_AND_SOIL)
       id_val        = L2E_COLUMN_ZI_SNOW_AND_SOIL
       name_val      = 'Column layer interface depth'
       long_name_val = 'Column layer interface depth: ALM to External Model'
       units_val     = '[m]'
       is_real_type  = .true.
       ndim          = 2
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc
       dim2_beg_name = dimname_nlevsno
       dim2_end_name = dimname_nlevgrnd

    case (L2E_COLUMN_DZ_SNOW_AND_SOIL)
       id_val        = L2E_COLUMN_DZ_SNOW_AND_SOIL
       name_val      = 'Column layer thickness'
       long_name_val = 'Column layer thickness: ALM to External Model'
       units_val     = '[m]'
       is_real_type  = .true.
       ndim          = 2
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc
       dim2_beg_name = dimname_nlevsno_plus_one
       dim2_end_name = dimname_nlevgrnd

    case (L2E_COLUMN_Z_SNOW_AND_SOIL)
       id_val        = L2E_COLUMN_Z_SNOW_AND_SOIL
       name_val      = 'Column layer centroid depth'
       long_name_val = 'Column layer centroid depth: ALM to External Model'
       units_val     = '[m]'
       is_real_type  = .true.
       ndim          = 2
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc
       dim2_beg_name = dimname_nlevsno_plus_one
       dim2_end_name = dimname_nlevgrnd

       ! --------------------------------------------------------------
       ! ALM-to-ELM: Landunit variables
       ! --------------------------------------------------------------
    case (L2E_LANDUNIT_TYPE)
       id_val        = L2E_LANDUNIT_TYPE
       name_val      = 'Landunit type'
       long_name_val = 'Landunit type: ALM to External Model'
       units_val     = '[-]'
       is_int_type   = .true.
       ndim          = 1
       dim1_beg_name = dimname_begl
       dim1_end_name = dimname_endl

    case (L2E_LANDUNIT_LAKEPOINT)
       id_val        = L2E_LANDUNIT_LAKEPOINT
       name_val      = 'Landunit lake point'
       long_name_val = 'Landunit lake point: ALM to External Model'
       units_val     = '[-]'
       is_int_type   = .true.
       ndim          = 1
       dim1_beg_name = dimname_begl
       dim1_end_name = dimname_endl

    case (L2E_LANDUNIT_URBANPOINT)
       id_val        = L2E_LANDUNIT_URBANPOINT
       name_val      = 'Landunit urban point'
       long_name_val = 'Landunit urban point: ALM to External Model'
       units_val     = '[-]'
       is_int_type   = .true.
       ndim          = 1
       dim1_beg_name = dimname_begl
       dim1_end_name = dimname_endl

       ! --------------------------------------------------------------
       ! ALM-to-ELM: Parameters variables
       ! --------------------------------------------------------------
    case (L2E_PARAMETER_WATSATC)
       id_val        = L2E_PARAMETER_WATSATC
       name_val      = 'Soil porosity'
       long_name_val = 'Soil porosity: ALM to External Model'
       units_val     = '[m^3/m^3]'
       is_real_type  = .true.
       ndim          = 2
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc
       dim2_beg_name = dimname_one
       dim2_end_name = dimname_nlevgrnd

    case (L2E_PARAMETER_HKSATC)
       id_val        = L2E_PARAMETER_HKSATC
       name_val      = 'Soil hydraulic conductivity'
       long_name_val = 'Soil hydraulic conductivity at saturation: ALM to External Model'
       units_val     = '[mm/s]'
       is_real_type  = .true.
       ndim          = 2
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc
       dim2_beg_name = dimname_one
       dim2_end_name = dimname_nlevgrnd

    case (L2E_PARAMETER_BSWC)
       id_val        = L2E_PARAMETER_BSWC
       name_val      = 'Clapp and Hornberger parameter'
       long_name_val = 'Clapp and Hornberger parameter: ALM to External Model'
       units_val     = '[-]'
       is_real_type  = .true.
       ndim          = 2
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc
       dim2_beg_name = dimname_one
       dim2_end_name = dimname_nlevgrnd

    case (L2E_PARAMETER_SUCSATC)
       id_val        = L2E_PARAMETER_SUCSATC
       name_val      = 'Minimum soil suction'
       long_name_val = 'Minimum soil suction: ALM to External Model'
       units_val     = '[mm]'
       is_real_type  = .true.
       ndim          = 2
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc
       dim2_beg_name = dimname_one
       dim2_end_name = dimname_nlevgrnd

    case (L2E_PARAMETER_EFFPOROSITYC)
       id_val        = L2E_PARAMETER_EFFPOROSITYC
       name_val      = 'Effective porosity'
       long_name_val = 'Effective porosity: ALM to External Model'
       units_val     = '[-]'
       is_real_type  = .true.
       ndim          = 2
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc
       dim2_beg_name = dimname_one
       dim2_end_name = dimname_nlevgrnd

    case (L2E_PARAMETER_CSOL)
       id_val        = L2E_PARAMETER_CSOL
       name_val      = 'Heat capacity'
       long_name_val = 'Heat capacity: ALM to External Model'
       units_val     = '[W/m/K]'
       is_real_type  = .true.
       ndim          = 2
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc
       dim2_beg_name = dimname_one
       dim2_end_name = dimname_nlevgrnd

    case (L2E_PARAMETER_TKMG)
       id_val        = L2E_PARAMETER_TKMG
       name_val      = 'Thermal conductivity minerals'
       long_name_val = 'Thermal conductivity minerals: ALM to External Model'
       units_val     = '[W/m/K]'
       is_real_type  = .true.
       ndim          = 2
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc
       dim2_beg_name = dimname_one
       dim2_end_name = dimname_nlevgrnd

    case (L2E_PARAMETER_TKDRY)
       id_val        = L2E_PARAMETER_TKDRY
       name_val      = 'Thermal conductivity dry soils'
       long_name_val = 'Thermal conductivity dry soils: ALM to External Model'
       units_val     = '[W/m/K]'
       is_real_type  = .true.
       ndim          = 2
       dim1_beg_name = dimname_begc
       dim1_end_name = dimname_endc
       dim2_beg_name = dimname_one
       dim2_end_name = dimname_nlevgrnd

    case default
       write(iulog,*)'Unknown EMID id = ',data_id
       call endrun(msg='EMIDListAddDataByID: Unknown EMID id.')
    end select

    allocate(data)
    call data%Init()
    call data%Setup(                        &
         id            = id_val,            &
         name          = name_val,          &
         long_name     = long_name_val,     &
         units         = units_val,         &
         num_em_stages = num_em_stages_val, &
         em_stage_ids  = em_stage_ids_val)
    call data%SetType(is_int_type, is_real_type)
    call data%SetNDimension(ndim)
    call data%SetDimBegEndNames( &
         dim1_beg_name, dim1_end_name, &
         dim2_beg_name, dim2_end_name, &
         dim3_beg_name, dim3_end_name, &
         dim4_beg_name, dim4_end_name)
    call this%AddData(data)
    index_of_new_data = this%num_data

    nullify(data)

  end subroutine EMIDListAddDataByID

  !------------------------------------------------------------------------
  subroutine EMIDListCopy(this, data_list)
    !
    ! !DESCRIPTION:
    ! Makes a copy of EMI data list
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(emi_data_list)              :: this
    class(emi_data_list) , intent(in) :: data_list
    !
    ! !LOCAL VARIABLES:
    class(emi_data)      , pointer    :: cur_data
    class(emi_data)      , pointer    :: new_data
    integer                           :: idata

    call this%Init()

    allocate(this%data_ptr(data_list%num_data))

    idata = 0
    cur_data => data_list%first
    do
       if (.not.associated(cur_data)) exit

       idata = idata + 1

       allocate(new_data)

       call new_data%Copy(cur_data)

       call this%AddData(new_data)

       this%data_ptr(idata)%data => new_data

       nullify(new_data)

       cur_data => cur_data%next
    enddo

  end subroutine EMIDListCopy

  !------------------------------------------------------------------------
  subroutine EMIDListGetIntValue(this, data_index, int_value)
    !
    ! !DESCRIPTION:
    ! Gives access to data pointer
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(emi_data_list)              :: this
    integer , intent(in) :: data_index
    integer, intent(out) :: int_value
    !
    ! !LOCAL VARIABLES:
    integer :: idx

    if (data_index > this%num_data) then
       call endrun(msg='EMIDListGetPointerToData: Attempting to access ' // &
          'data index that is larger than number of data in the list.')
    endif

    if (.not. associated(this%data_ptr(data_index)%data%data_int_1d)) then
       call endrun(msg='EMIDListGetPointerToData: Attempting to asscess ' // &
          'an unallocated data_int_1d.')
    else
       if (this%data_ptr(data_index)%data%dim1_end /= &
           this%data_ptr(data_index)%data%dim1_beg) then
          call endrun(msg='EMIDListGetIntValue: Only extracts values from data ' // &
             'that has 1 value.')
       endif
       idx = this%data_ptr(data_index)%data%dim1_beg
       int_value = this%data_ptr(data_index)%data%data_int_1d(idx)
    endif

  end subroutine EMIDListGetIntValue

  !------------------------------------------------------------------------
  subroutine EMIDListGetPointerToInt1D(this, data_index, int_1d)
    !
    ! !DESCRIPTION:
    ! Gives access to data pointer
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(emi_data_list)              :: this
    integer , intent(in) :: data_index
    integer, pointer :: int_1d(:)
    !
    ! !LOCAL VARIABLES:

    if (data_index > this%num_data) then
       call endrun(msg='EMIDListGetPointerToData: Attempting to access ' // &
          'data index that is larger than number of data in the list.')
    endif

    if (.not. associated(this%data_ptr(data_index)%data%data_int_1d)) then
       call endrun(msg='EMIDListGetPointerToData: Attempting to asscess ' // &
          'an unallocated data_int_1d.')
    else
       int_1d => this%data_ptr(data_index)%data%data_int_1d
    endif

  end subroutine EMIDListGetPointerToInt1D

  !------------------------------------------------------------------------
  subroutine EMIDListGetPointerToInt2D(this, data_index, int_2d)
    !
    ! !DESCRIPTION:
    ! Gives access to data pointer
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(emi_data_list)              :: this
    integer , intent(in) :: data_index
    integer, pointer :: int_2d(:,:)
    !
    ! !LOCAL VARIABLES:

    if (data_index > this%num_data) then
       call endrun(msg='EMIDListGetPointerToData: Attempting to access ' // &
          'data index that is larger than number of data in the list.')
    endif

    if (.not. associated(this%data_ptr(data_index)%data%data_int_2d)) then
       call endrun(msg='EMIDListGetPointerToData: Attempting to asscess ' // &
          'an unallocated data_int_2d.')
    else
       int_2d => this%data_ptr(data_index)%data%data_int_2d
    endif

  end subroutine EMIDListGetPointerToInt2D

  !------------------------------------------------------------------------
  subroutine EMIDListGetPointerToInt3D(this, data_index, int_3d)
    !
    ! !DESCRIPTION:
    ! Gives access to data pointer
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(emi_data_list)              :: this
    integer , intent(in) :: data_index
    integer, pointer :: int_3d(:,:,:)
    !
    ! !LOCAL VARIABLES:

    if (data_index > this%num_data) then
       call endrun(msg='EMIDListGetPointerToData: Attempting to access ' // &
          'data index that is larger than number of data in the list.')
    endif

    if (.not. associated(this%data_ptr(data_index)%data%data_int_3d)) then
       call endrun(msg='EMIDListGetPointerToData: Attempting to asscess ' // &
          'an unallocated data_int_3d.')
    else
       int_3d => this%data_ptr(data_index)%data%data_int_3d
    endif

  end subroutine EMIDListGetPointerToInt3D

  !------------------------------------------------------------------------
  subroutine EMIDListGetPointerToReal1D(this, data_index, real_1d)
    !
    ! !DESCRIPTION:
    ! Gives access to data pointer
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(emi_data_list)              :: this
    integer , intent(in) :: data_index
    real(r8), pointer :: real_1d(:)
    !
    ! !LOCAL VARIABLES:

    if (data_index > this%num_data) then
       call endrun(msg='EMIDListGetPointerToData: Attempting to access ' // &
          'data index that is larger than number of data in the list.')
    endif

    if (.not. associated(this%data_ptr(data_index)%data%data_real_1d)) then
       call endrun(msg='EMIDListGetPointerToData: Attempting to asscess ' // &
          'an unallocated data_real_1d.')
    else
       real_1d => this%data_ptr(data_index)%data%data_real_1d
    endif

  end subroutine EMIDListGetPointerToReal1D

  !------------------------------------------------------------------------
  subroutine EMIDListGetPointerToReal2D(this, data_index, real_2d)
    !
    ! !DESCRIPTION:
    ! Gives access to data pointer
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(emi_data_list)              :: this
    integer , intent(in) :: data_index
    real(r8), pointer :: real_2d(:,:)
    !
    ! !LOCAL VARIABLES:

    if (data_index > this%num_data) then
       call endrun(msg='EMIDListGetPointerToData: Attempting to access ' // &
          'data index that is larger than number of data in the list.')
    endif

    if (.not. associated(this%data_ptr(data_index)%data%data_real_2d)) then
       call endrun(msg='EMIDListGetPointerToData: Attempting to asscess ' // &
          'an unallocated data_real_2d.')
    else
       real_2d => this%data_ptr(data_index)%data%data_real_2d
    endif

  end subroutine EMIDListGetPointerToReal2D

  !------------------------------------------------------------------------
  subroutine EMIDListGetPointerToReal3D(this, data_index, real_3d)
    !
    ! !DESCRIPTION:
    ! Gives access to data pointer
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(emi_data_list)              :: this
    integer , intent(in) :: data_index
    real(r8), pointer :: real_3d(:,:,:)
    !
    ! !LOCAL VARIABLES:

    if (data_index > this%num_data) then
       call endrun(msg='EMIDListGetPointerToData: Attempting to access ' // &
          'data index that is larger than number of data in the list.')
    endif

    if (.not. associated(this%data_ptr(data_index)%data%data_real_3d)) then
       call endrun(msg='EMIDListGetPointerToData: Attempting to asscess ' // &
          'an unallocated data_real_3d.')
    else
       real_3d => this%data_ptr(data_index)%data%data_real_3d
    endif

  end subroutine EMIDListGetPointerToReal3D

  !------------------------------------------------------------------------
  subroutine EMIDListGetPointerToReal4D(this, data_index, real_4d)
    !
    ! !DESCRIPTION:
    ! Gives access to data pointer
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(emi_data_list)              :: this
    integer , intent(in) :: data_index
    real(r8), pointer :: real_4d(:,:,:,:)
    !
    ! !LOCAL VARIABLES:

    if (data_index > this%num_data) then
       call endrun(msg='EMIDListGetPointerToData: Attempting to access ' // &
          'data index that is larger than number of data in the list.')
    endif

    if (.not. associated(this%data_ptr(data_index)%data%data_real_4d)) then
       call endrun(msg='EMIDListGetPointerToData: Attempting to asscess ' // &
          'an unallocated data_real_4d.')
    else
       real_4d => this%data_ptr(data_index)%data%data_real_4d
    endif

  end subroutine EMIDListGetPointerToReal4D

  !------------------------------------------------------------------------
  subroutine EMIDListDestroy(this)
    !
    ! !DESCRIPTION:
    ! Destroys a EMID list
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(emi_data_list) :: this
    ! !LOCAL VARIABLES:
    class(emi_data)      , pointer    :: cur_data
    class(emi_data)      , pointer    :: old_data

    cur_data => this%first
    do
       if (.not.associated(cur_data)) exit

       old_data => cur_data
       cur_data => cur_data%next
       call old_data%Destroy()

    enddo

    this%num_data = 0

    nullify(this%first)
    nullify(this%last)
    nullify(this%data_ptr)

  end subroutine EMIDListDestroy

end module ExternalModelInterfaceDataMod
