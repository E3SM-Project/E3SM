module iondrag
  !-------------------------------------------------------------------------------
  !  Dummy interface for waccm/iondrag module
  !-------------------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid       ,only: pver
  use physics_types,only: physics_state, physics_ptend
  use physics_buffer ,only: physics_buffer_desc

  implicit none

  save

  private                         ! Make default type private to the module

  !-------------------------------------------------------------------------------
  ! Public interfaces:
  !-------------------------------------------------------------------------------
  public :: iondrag_register         ! Register variables in pbuf physics buffer
  public :: iondrag_init             ! Initialization
  public :: iondrag_calc             ! ion drag tensors lxx,lyy,lxy,lyx
  public :: iondrag_readnl
  public :: do_waccm_ions

  interface iondrag_calc
     module procedure iondrag_calc_ions
     module procedure iondrag_calc_ghg
  end interface

  logical, parameter :: do_waccm_ions = .false.

contains

  !================================================================================================

  subroutine iondrag_readnl(nlfile)

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

  end subroutine iondrag_readnl

  !==============================================================================     

  subroutine iondrag_register

  end subroutine iondrag_register

  !================================================================================================

  subroutine iondrag_init( pref_mid )
   
    !-------------------------------------------------------------------------------
    ! dummy arguments
    !-------------------------------------------------------------------------------
    real(r8), intent(in) :: pref_mid(pver)

  end subroutine iondrag_init

  !================================================================================================
  subroutine iondrag_calc_ions( lchnk, ncol, state, ptend, pbuf, delt )

    !-------------------------------------------------------------------------------
    ! dummy arguments
    !-------------------------------------------------------------------------------
    integer,intent(in)   :: lchnk               ! current chunk index
    integer,intent(in)   :: ncol                ! number of atmospheric columns
    real(r8), intent(in) :: delt                ! time step (s)
    type(physics_state), intent(in), target    :: state ! Physics state variables
    type(physics_ptend), intent(out)   :: ptend   ! Physics tendencies
    type(physics_buffer_desc), pointer :: pbuf(:) ! physics buffer

  end subroutine iondrag_calc_ions

  !=========================================================================

  subroutine iondrag_calc_ghg (lchnk,ncol,state,ptend)

    !--------------------Input arguments------------------------------------

    integer, intent(in) :: lchnk                   ! chunk identifier
    integer, intent(in) :: ncol                    ! number of atmospheric columns

    type(physics_state), intent(in) :: state
    type(physics_ptend), intent(out):: ptend

  end subroutine iondrag_calc_ghg

  !===================================================================================

end module iondrag
