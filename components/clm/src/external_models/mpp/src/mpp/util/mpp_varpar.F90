module mpp_varpar
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  integer :: nlevsoi           = 10 ! number of hydrologically active soil layers
  integer :: nlevgrnd          = 15 ! number of ground layers
  integer :: nlevsno           =  5 ! maximum number of snow layers
  integer :: max_patch_per_col = 20
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public mpp_varpar_init          ! set parameters
  public mpp_varpar_set_nlevsoi
  public mpp_varpar_set_nlevgrnd
  public mpp_varpar_set_nlevsno
  !
  !-----------------------------------------------------------------------

contains

  !------------------------------------------------------------------------------
  subroutine mpp_varpar_init(nlevsoi_val, nlevgrnd_val, nlevsno_val, &
      max_patch_per_col_val)
    !
    ! !DESCRIPTION:
    ! Initialize module variables 
    !
    implicit none
    !
    ! !ARGUMENTS:
    integer, intent(in) :: nlevsoi_val
    integer, intent(in) :: nlevgrnd_val
    integer, intent(in) :: nlevsno_val
    integer, intent(in) :: max_patch_per_col_val

    nlevsoi           = nlevsoi_val
    nlevgrnd          = nlevgrnd_val
    nlevsno           = nlevsno_val
    max_patch_per_col = max_patch_per_col_val
    
  end subroutine mpp_varpar_init

  !------------------------------------------------------------------------------
  subroutine mpp_varpar_set_nlevsoi(nlev)
    !
    ! !DESCRIPTION:
    ! Set number of soil layers
    !
    ! !ARGUMENTS:
    implicit none
    !
    ! !LOCAL VARIABLES:
    !
    integer, intent (in) :: nlev

    nlevsoi  = nlev
    
  end subroutine mpp_varpar_set_nlevsoi

  !------------------------------------------------------------------------------
  subroutine mpp_varpar_set_nlevgrnd(nlev)
    !
    ! !DESCRIPTION:
    ! Set number of ground layers
    !
    ! !ARGUMENTS:
    implicit none
    !
    ! !LOCAL VARIABLES:
    !
    integer, intent (in) :: nlev

    nlevgrnd  = nlev
    
  end subroutine mpp_varpar_set_nlevgrnd

  !------------------------------------------------------------------------------
  subroutine mpp_varpar_set_nlevsno(nlev)
    !
    ! !DESCRIPTION:
    ! Set number of snow layers
    !
    ! !ARGUMENTS:
    implicit none
    !
    ! !LOCAL VARIABLES:
    !
    integer, intent (in) :: nlev

    nlevsno  = nlev
    
  end subroutine mpp_varpar_set_nlevsno

end module mpp_varpar
