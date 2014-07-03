module iondrag
  !-------------------------------------------------------------------------------
  !  Dummy interface for waccm/iondrag module
  !-------------------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid       ,only: pver
  use physics_types,only: physics_state, physics_ptend
  use physics_buffer ,only: physics_buffer_desc
  use abortutils,   only: endrun

  implicit none

  save

  private                         ! Make default type private to the module

  !-------------------------------------------------------------------------------
  ! Public interfaces:
  !-------------------------------------------------------------------------------
  public :: iondrag_register         ! Register variables in pbuf physics buffer
  public :: iondrag_init             ! Initialization
  public :: iondrag_calc             ! ion drag tensors lxx,lyy,lxy,lyx
  public :: iondrag_defaultopts
  public :: iondrag_setopts
  public :: do_waccm_ions

  interface iondrag_calc
     module procedure iondrag_calc_ions
     module procedure iondrag_calc_ghg
  end interface

  logical :: do_waccm_ions

contains

!==============================================================================     

  subroutine iondrag_register

    call endrun('iondrag_register: dummy interface should not be called')

  end subroutine iondrag_register

!================================================================================================

  subroutine iondrag_defaultopts( &
       efield_lflux_file_out,     &
       efield_hflux_file_out,     &
       efield_wei96_file_out    )

    character(len=*), intent(out), optional :: efield_lflux_file_out
    character(len=*), intent(out), optional :: efield_hflux_file_out
    character(len=*), intent(out), optional :: efield_wei96_file_out

    call endrun('iondrag_defaultopts: dummy interface should not be called')

  end subroutine iondrag_defaultopts

  !================================================================================================

  subroutine iondrag_setopts(    &
       efield_lflux_file_out,     &
       efield_hflux_file_out,     &
       efield_wei96_file_out    )

    character(len=*), intent(out), optional :: efield_lflux_file_out
    character(len=*), intent(out), optional :: efield_hflux_file_out
    character(len=*), intent(out), optional :: efield_wei96_file_out

    call endrun('iondrag_setopts: dummy interface should not be called')

  end subroutine iondrag_setopts

  !================================================================================================

  subroutine iondrag_init( pref_mid )
   
    !-------------------------------------------------------------------------------
    ! dummy arguments
    !-------------------------------------------------------------------------------
    real(r8), intent(in) :: pref_mid(pver)

    call endrun('iondrag_init: dummy interface should not be called')

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

    call endrun('iondrag_calc_ions: dummy interface should not be called')

  end subroutine iondrag_calc_ions

  !=========================================================================

  subroutine iondrag_calc_ghg (lchnk,ncol,state,ptend, pbuf)

    !--------------------Input arguments------------------------------------

    integer, intent(in) :: lchnk                   ! chunk identifier
    integer, intent(in) :: ncol                    ! number of atmospheric columns

    type(physics_state), intent(in) :: state
    type(physics_ptend), intent(out):: ptend
    type(physics_buffer_desc), pointer :: pbuf(:) ! physics buffer

    call endrun('iondrag_calc_ghg: dummy interface should not be called')

  end subroutine iondrag_calc_ghg

  !===================================================================================

end module iondrag
