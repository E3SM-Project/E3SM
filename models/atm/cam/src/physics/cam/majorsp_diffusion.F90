module majorsp_diffusion

!--------------------------------------------------------------------------
! Dummy interface for waccmx/majorsp_diffusion module
!--------------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid,       only: pver
  use physics_types,only: physics_state, physics_ptend
  use abortutils,   only: endrun

  implicit none

  private          ! Make default type private to the module
  save
!-----------------------
! Public interfaces
!-----------------------
  public mspd_init   ! Initialization
  public mspd_intr   ! Full routine

contains

!===============================================================================
  subroutine mspd_init()
  
    call endrun('mspd_init: dummy interface should not be called')

  end subroutine mspd_init

!===============================================================================
  subroutine mspd_intr(ztodt    ,state    ,ptend)

    !------------------------------Arguments--------------------------------
    real(r8), intent(in) :: ztodt                  ! 2 delta-t
    type(physics_state), intent(in)     :: state   ! Physics state variables
    type(physics_ptend), intent(inout)  :: ptend   ! indivdual parameterization tendencies

    call endrun('mspd_intr: dummy interface should not be called')

  end subroutine mspd_intr

end module majorsp_diffusion

