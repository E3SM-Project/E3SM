module mpp_varpar
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  integer :: nlevsoi   ! number of hydrologically active soil layers
  integer :: nlevgrnd  ! number of ground layers
  integer :: nlevsno   ! maximum number of snow layers
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
  subroutine mpp_varpar_init()
    !
    ! !DESCRIPTION:
    ! Initialize module variables 
    !
    ! !ARGUMENTS:
    implicit none
    !
    ! !LOCAL VARIABLES:
    !
    nlevsoi  = 10
    nlevgrnd = 15
    nlevsno  = 5
    
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
