
module rayleigh_friction

  !---------------------------------------------------------------------------------
  ! Module to apply rayleigh friction in region of model top.
  ! We specify a decay rate profile that is largest at the model top and
  ! drops off vertically using a hyperbolic tangent profile.
  ! We compute the tendencies in u and v using an Euler backward scheme.
  ! We then apply the negative of the kinetic energy tendency to "s", the dry
  ! static energy.
  !
  ! calling sequence:
  !
  !  rayleigh_friction_init          initializes rayleigh friction constants
  !  rayleigh_friction_tend          computes rayleigh friction tendencies
  !
  !---------------------------Code history--------------------------------
  ! This is a new routine written by Art Mirin in collaboration with Phil Rasch.
  ! Initial coding for this version:  Art Mirin, May 2007.
  !---------------------------------------------------------------------------------

  use shr_kind_mod,     only: r8 => shr_kind_r8
  use ppgrid,           only: pver
  use spmd_utils,       only: masterproc
  use cam_logfile,      only: iulog

  implicit none
  private          ! Make default type private to the module
  save
  
  !
  ! Public interfaces
  !
  public rayleigh_friction_init          ! Initialization
  public rayleigh_friction_tend          ! Computation of tendencies

  !
  ! Public data
  !
  integer, public   :: rayk0 = 2           ! vertical level at which rayleigh friction term is centered
  real (r8), public :: raykrange = 0._r8   ! range of rayleigh friction profile 
                                           ! if 0, range is set to satisfy x=2 (see below)
  real (r8), public :: raytau0 = huge(1._r8)! approximate value of decay time at model top (days)
                                           ! if 0., no rayleigh friction is applied (set by namelist)

  ! 
  ! Private data
  !
  real (r8) :: krange         ! range of rayleigh friction profile 
  real (r8) :: tau0           ! approximate value of decay time at model top
  real (r8) :: otau0          ! inverse of tau0
  real (r8) :: otau(pver)     ! inverse decay time versus vertical level

  ! We apply a profile of the form otau0 * [1 + tanh (x)] / 2 , where
  ! x = (k0 - k) / krange. The default is for x to equal 2 at k=1, meaning
  ! krange = (k0 - 1) / 2. The default is applied when raykrange is set to 0.
  ! If otau0 = 0, no term is applied.

contains

  !===============================================================================
  subroutine rayleigh_friction_init()
    !------------------------------Arguments--------------------------------

    !---------------------------Local storage-------------------------------
    real (r8) x
    integer k

    !-----------------------------------------------------------------------
    ! Compute tau array
    !-----------------------------------------------------------------------

    krange = raykrange
    if (raykrange .eq. 0._r8) krange = (rayk0 - 1) / 2._r8

    tau0 = (86400._r8) * raytau0   ! convert to seconds
    otau0 = 0._r8
    if (tau0 .ne. 0._r8) otau0 = 1._r8/tau0

    do k = 1, pver
       x = (rayk0 - k) / krange
       otau(k) = otau0 * (1 + tanh(x)) / (2._r8)
    enddo

    if (masterproc) then
       write (iulog,*) 'Rayleigh friction - rayk0 = ', rayk0
       write (iulog,*) 'Rayleigh friction - raykrange = ', raykrange
       write (iulog,*) 'Rayleigh friction - raytau0 = ', raytau0
       write (iulog,*) 'Rayleigh friction - krange = ', krange
       write (iulog,*) 'Rayleigh friction - otau0 = ', otau0
       write (iulog,*) 'Rayleigh friction decay rate profile'
       do k = 1, pver
          write (iulog,*) '   k = ', k, '   otau = ', otau(k)
       enddo
    end if

    return

  end subroutine rayleigh_friction_init
  
!=========================================================================================
  subroutine rayleigh_friction_tend(                                     &
       ztodt    ,state    ,ptend    )
    !-----------------------------------------------------------------------
    ! interface routine for rayleigh friction
    !-----------------------------------------------------------------------
    use physics_types, only: physics_state, physics_ptend, physics_ptend_init


    !------------------------------Arguments--------------------------------
    real(r8), intent(in) :: ztodt                  ! physics timestep
    type(physics_state), intent(in)  :: state      ! physics state variables
    
    type(physics_ptend), intent(out) :: ptend      ! individual parameterization tendencies
    !
    !---------------------------Local storage-------------------------------
    integer :: ncol                                ! number of atmospheric columns
    integer :: k                                   ! level
    real(r8) :: rztodt                             ! 1./ztodt
    real(r8) :: c1, c2, c3                         ! temporary variables
    !-----------------------------------------------------------------------

    call physics_ptend_init(ptend, state%psetcols, 'rayfri', ls=.true., lu=.true., lv=.true.)

    if (otau0 .eq. 0._r8) return

    rztodt = 1._r8/ztodt
    ncol  = state%ncol

    ! u, v and s are modified by rayleigh friction

    do k = 1, pver
       c2 = 1._r8 / (1._r8 + otau(k)*ztodt)
       c1 = -otau(k) * c2
       c3 = 0.5_r8 * (1._r8 - c2*c2) * rztodt
       ptend%u(:ncol,k) = c1 * state%u(:ncol,k)
       ptend%v(:ncol,k) = c1 * state%v(:ncol,k)
       ptend%s(:ncol,k) = c3 * (state%u(:ncol,k)**2 + state%v(:ncol,k)**2)
    enddo

    return
  end subroutine rayleigh_friction_tend

end module rayleigh_friction
