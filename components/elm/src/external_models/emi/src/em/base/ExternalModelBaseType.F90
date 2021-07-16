module ExternalModelBaseType

  use abortutils                   , only : endrun
  use shr_kind_mod                 , only : r8 => shr_kind_r8
  use shr_log_mod                  , only : errMsg => shr_log_errMsg
  use decompMod                    , only : bounds_type
  use EMI_DataMod, only : emi_data_list, emi_data
  use elm_varctl                   , only : iulog

  type, public :: em_base_type
   contains
     procedure, public :: Populate_L2E_Init_List  => EMBase_Populate_L2E_Init_List
     procedure, public :: Populate_E2L_Init_List  => EMBase_Populate_E2L_Init_List
     procedure, public :: Populate_L2E_List       => EMBase_Populate_L2E_List
     procedure, public :: Populate_E2L_List       => EMBase_Populate_E2L_List
     procedure, public :: Init                    => EMBase_Init
     procedure, public :: Solve                   => EMBase_Solve
  end type em_base_type


contains

  !------------------------------------------------------------------------
  subroutine EMBase_Populate_L2E_Init_List(this, l2e_init_list)
    !
    ! !DESCRIPTION:
    ! Initialze an emi_list for exchanging data from land model to external
    ! model during initialization stage
    !
    implicit none
    !
    ! !ARGUMENTS
    class(em_base_type)                 :: this
    class(emi_data_list), intent(inout) :: l2e_init_list

    write(iulog,*)'EMBase_Populate_L2E_Init_List must be extended by a child class.'
    call endrun(msg=errMsg(__FILE__, __LINE__))

  end subroutine EMBase_Populate_L2E_Init_List

  !------------------------------------------------------------------------
  subroutine EMBase_Populate_E2L_Init_List(this, e2l_init_list)
    !
    ! !DESCRIPTION:
    ! Initialze an emi_list for exchanging data from external model to land
    ! model during initialization stage
    !
    implicit none
    !
    ! !ARGUMENTS
    class(em_base_type)                 :: this
    class(emi_data_list), intent(inout) :: e2l_init_list

    write(iulog,*)'EMBase_Populate_E2L_Init_List must be extended by a child class.'
    call endrun(msg=errMsg(__FILE__, __LINE__))

  end subroutine EMBase_Populate_E2L_Init_List

  !------------------------------------------------------------------------
  subroutine EMBase_Populate_L2E_List(this, l2e_list)
    !
    ! !DESCRIPTION:
    ! Initialze an emi_list for exchanging data from land model to external
    ! model during time integration stage
    !
    implicit none
    !
    ! !ARGUMENTS
    class(em_base_type)                 :: this
    class(emi_data_list), intent(inout) :: l2e_list

    write(iulog,*)'EMBase_Populate_L2E_List must be extended by a child class.'
    call endrun(msg=errMsg(__FILE__, __LINE__))

  end subroutine EMBase_Populate_L2E_List

  !------------------------------------------------------------------------
  subroutine EMBase_Populate_E2L_List(this, e2l_list)
    !
    ! !DESCRIPTION:
    ! Initialze an emi_list for exchanging data from external model to land
    ! model during time integration stage
    !
    implicit none
    !
    ! !ARGUMENTS
    class(em_base_type)                 :: this
    class(emi_data_list), intent(inout) :: e2l_list

    write(iulog,*)'EMBase_Populate_E2L_List must be extended by a child class.'
    call endrun(msg=errMsg(__FILE__, __LINE__))

  end subroutine EMBase_Populate_E2L_List

  !------------------------------------------------------------------------
  subroutine EMBase_Init(this, l2e_init_list, e2l_init_list, iam, bounds_clump)
    !
    ! !DESCRIPTION:
    ! Initialze an emi_list for exchanging data from land model to external
    ! model during time integration stage
    !
    implicit none
    !
    ! !ARGUMENTS
    class(em_base_type)                  :: this
    class(emi_data_list) , intent(in)    :: l2e_init_list
    class(emi_data_list) , intent(inout) :: e2l_init_list
    integer              , intent(in)    :: iam
    type(bounds_type)    , intent (in)   :: bounds_clump

    write(iulog,*)'EMBase_Init must be extended by a child class.'
    call endrun(msg=errMsg(__FILE__, __LINE__))

  end subroutine EMBase_Init

  !------------------------------------------------------------------------
  subroutine EMBase_Solve(this, em_stage, dt, nstep, clump_rank, l2e_list, e2l_list, &
       bounds_clump)
    !
    ! !DESCRIPTION:
    ! Initialze an emi_list for exchanging data from land model to external
    ! model during time integration stage
    !
    implicit none
    !
    ! !ARGUMENTS
    class(em_base_type)                  :: this
    integer              , intent(in)    :: em_stage
    real(r8)             , intent(in)    :: dt
    integer              , intent(in)    :: nstep
    integer              , intent(in)    :: clump_rank
    class(emi_data_list) , intent(in)    :: l2e_list
    class(emi_data_list) , intent(inout) :: e2l_list
    type(bounds_type)    , intent (in)   :: bounds_clump

    write(iulog,*)'EMBase_Solve must be extended by a child class.'
    call endrun(msg=errMsg(__FILE__, __LINE__))

  end subroutine EMBase_Solve

end module ExternalModelBaseType
