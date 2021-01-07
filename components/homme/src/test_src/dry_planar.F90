
#ifndef CAM
#include "config.h"

module dry_planar_tests

  use element_mod,          only: element_t
  use hybrid_mod,           only: hybrid_t
  use hybvcoord_mod,        only: hvcoord_t, set_layer_locations
  use parallel_mod,         only: abortmp
  use element_ops,          only: set_state, set_state_i, tests_finalize
  use kinds,                only: rl=>real_kind, iulog
  use physical_constants,   only : dd_pi
  use dcmip12_wrapper, only : set_tracers, get_evenly_spaced_z, get_evenly_spaced_p, set_hybrid_coefficients, pressure_thickness
  use dimensions_mod,       only: np, nlev, nlevp, qsize, qsize_d, nelemd
  use element_state,        only: nt=>timelevels
  use control_mod,          only: planar_slice

  implicit none

  real(rl), parameter :: Rd 	= 287.0d0,	&	! Ideal gas const dry air (J kg^-1 K^1)
				                        g	= 9.80616d0,	&	! Gravity (m s^2)
				                        cp	= 1004.5d0,	&	! Specific heat capacity (J kg^-1 K^1)
                                p0	= 100000.d0, &! reference pressure (Pa)
                                kappa   = Rd/cp

  real(rl):: zi(nlevp), zm(nlev)                                          ! z coordinates
  real(rl):: ddn_hyai(nlevp), ddn_hybi(nlevp)                             ! vertical derivativess of hybrid coefficients
  real(rl):: ztop

  contains

! planar hydrostatic gravity wave
subroutine planar_hydro_gravity_wave_init(elem,hybrid,hvcoord,nets,nete)


  type(element_t),    intent(inout), target :: elem(:)                  ! element array
  type(hybrid_t),     intent(in)            :: hybrid                   ! hybrid parallel structure
  type(hvcoord_t),    intent(inout)         :: hvcoord                  ! hybrid vertical coordinates
  integer,            intent(in)            :: nets,nete                ! start, end element index

  call general_gravity_wave_init(elem,hybrid,hvcoord,nets,nete,100.0D0 * 1000.0D0, 0.0D0) ! EVENTUALLY ALLOW A CONSTATNT CORIOLIS FORCE HERE...

end subroutine planar_hydro_gravity_wave_init

! planar nonhydrostatic gravity wave
subroutine planar_nonhydro_gravity_wave_init(elem,hybrid,hvcoord,nets,nete)


  type(element_t),    intent(inout), target :: elem(:)                  ! element array
  type(hybrid_t),     intent(in)            :: hybrid                   ! hybrid parallel structure
  type(hvcoord_t),    intent(inout)         :: hvcoord                  ! hybrid vertical coordinates
  integer,            intent(in)            :: nets,nete                ! start, end element index

  call general_gravity_wave_init(elem,hybrid,hvcoord,nets,nete,5.0D0 * 1000.0D0, 0.0D0)

end subroutine planar_nonhydro_gravity_wave_init



subroutine general_gravity_wave_init(elem,hybrid,hvcoord,nets,nete,d,f)

  type(element_t),    intent(inout), target :: elem(:)                  ! element array
  type(hybrid_t),     intent(in)            :: hybrid                   ! hybrid parallel structure
  type(hvcoord_t),    intent(inout)         :: hvcoord                  ! hybrid vertical coordinates
  integer,            intent(in)            :: nets,nete                ! start, end element index
  real(rl),           intent(in)            :: d                        ! radius of perturbation
  real(rl),           intent(in)            :: f                        ! (const) Coriolis force

  integer,  parameter :: zcoords = 0                                    ! we are not using z coords
  logical,  parameter :: use_eta = .true.                               ! we are using hybrid eta coords

  real(rl), parameter ::    &                                           ! parameters needed to get eta from z
    T0      = 300.d0,       &	! temperature (k)
    ztop    = 10000.d0,     & ! model top (m)
    N       = 0.01d0,       & ! Brunt-Vaisala frequency
    bigG    = (g*g)/(N*N*Cp)  ! temperature, isothermal

  integer :: i,j,k,ie                                                   ! loop indices
  real(rl):: x,y,hyam,hybm,hyai,hybi                                ! pointwise coordiantes
  real(rl):: p,z,phis,u,v,w,T,T_mean,phis_ps,ps,rho,rho_mean,q(1),dp    ! pointwise field values

  if (hybrid%masterthread) write(iulog,*) 'initializing planar gravity wave'

  ! set analytic vertical coordinates
  call get_evenly_spaced_z(zi,zm, 0.0_rl,ztop)                                   ! get evenly spaced z levels
  hvcoord%etai  = ( (bigG/T0)*(exp(-zi*N*N/g) -1 )+1 ) **(1.0/kappa)    ! set eta levels from z
  call set_hybrid_coefficients(hvcoord,hybrid,  hvcoord%etai(1), 1.0_rl)! set hybrid A and B from eta levels
  call set_layer_locations(hvcoord, .true., hybrid%masterthread)

  ! set initial conditions
  do ie = nets,nete
     do k=1,nlev; do j=1,np; do i=1,np
        call get_xycoordinates(x,y,hyam,hybm, i,j,k,elem(ie),hvcoord)
        call gravity_wave(x,y,p,z,zcoords,use_eta,hyam,hybm,d,u,v,w,T,T_mean,phis,ps,rho,rho_mean,q(1))
        dp = pressure_thickness(ps,k,hvcoord)
        call set_state(u,v,w,T,ps,phis,p,dp,zm(k),g, i,j,k,elem(ie),1,nt)
        call set_tracers(q,qsize, dp,i,j,k,y,x,elem(ie))
     enddo; enddo; enddo;
     do k=1,nlevp; do j=1,np; do i=1,np
        call get_xycoordinates(x,y,hyai,hybi, i,j,k,elem(ie),hvcoord)
        call gravity_wave(x,y,p,z,zcoords,use_eta,hyai,hybi,d,u,v,w,T,T_mean,phis,ps,rho,rho_mean,q(1))
        call set_state_i(u,v,w,T,ps,phis,p,zi(k),g, i,j,k,elem(ie),1,nt)
     enddo; enddo; enddo;
     elem(ie)%fcor(:,:) = f
     call tests_finalize(elem(ie),hvcoord)
  enddo

end subroutine general_gravity_wave_init

! planar hydrostatic mountain wave
subroutine planar_hydro_mountain_wave_init(elem,hybrid,hvcoord,nets,nete)


  type(element_t),    intent(inout), target :: elem(:)                  ! element array
  type(hybrid_t),     intent(in)            :: hybrid                   ! hybrid parallel structure
  type(hvcoord_t),    intent(inout)         :: hvcoord                  ! hybrid vertical coordinates
  integer,            intent(in)            :: nets,nete                ! start, end element index

  call abortmp('planar hydrostatic mountain wave not yet implemented')

end subroutine planar_hydro_mountain_wave_init

! planar nonhydrostatic mountain wave
subroutine planar_nonhydro_mountain_wave_init(elem,hybrid,hvcoord,nets,nete)


  type(element_t),    intent(inout), target :: elem(:)                  ! element array
  type(hybrid_t),     intent(in)            :: hybrid                   ! hybrid parallel structure
  type(hvcoord_t),    intent(inout)         :: hvcoord                  ! hybrid vertical coordinates
  integer,            intent(in)            :: nets,nete                ! start, end element index

  call abortmp('planar nonhydrostatic mountain wave not yet implemented')

end subroutine planar_nonhydro_mountain_wave_init

! planar Schar mountain wave
subroutine planar_schar_mountain_wave_init(elem,hybrid,hvcoord,nets,nete)


  type(element_t),    intent(inout), target :: elem(:)                  ! element array
  type(hybrid_t),     intent(in)            :: hybrid                   ! hybrid parallel structure
  type(hvcoord_t),    intent(inout)         :: hvcoord                  ! hybrid vertical coordinates
  integer,            intent(in)            :: nets,nete                ! start, end element index

  call abortmp('planar Schar mountain wave not yet implemented')

end subroutine planar_schar_mountain_wave_init


! planar density current
subroutine planar_density_current_init(elem,hybrid,hvcoord,nets,nete)


  type(element_t),    intent(inout), target :: elem(:)                  ! element array
  type(hybrid_t),     intent(in)            :: hybrid                   ! hybrid parallel structure
  type(hvcoord_t),    intent(inout)         :: hvcoord                  ! hybrid vertical coordinates
  integer,            intent(in)            :: nets,nete                ! start, end element index

  call abortmp('planar density current not yet implemented')

end subroutine planar_density_current_init

! planar rising bubble
subroutine planar_rising_bubble_init(elem,hybrid,hvcoord,nets,nete)


  type(element_t),    intent(inout), target :: elem(:)                  ! element array
  type(hybrid_t),     intent(in)            :: hybrid                   ! hybrid parallel structure
  type(hvcoord_t),    intent(inout)         :: hvcoord                  ! hybrid vertical coordinates
  integer,            intent(in)            :: nets,nete                ! start, end element index

  call abortmp('planar rising bubble not yet implemented')

end subroutine planar_rising_bubble_init


! planar baroclinic instability
subroutine planar_baroclinic_instab_init(elem,hybrid,hvcoord,nets,nete)


  type(element_t),    intent(inout), target :: elem(:)                  ! element array
  type(hybrid_t),     intent(in)            :: hybrid                   ! hybrid parallel structure
  type(hvcoord_t),    intent(inout)         :: hvcoord                  ! hybrid vertical coordinates
  integer,            intent(in)            :: nets,nete                ! start, end element index

  call abortmp('planar baroclinic instability not yet implemented')

end subroutine planar_baroclinic_instab_init



SUBROUTINE gravity_wave (x,y,p,z,zcoords,hybrid_eta,hyam,hybm,d,u,v,w,t,t_mean,phis,ps,rho,rho_mean,q)

IMPLICIT NONE
!-----------------------------------------------------------------------
!     input/output params parameters at given location
!-----------------------------------------------------------------------

	real(rl), intent(in)     :: x        ! x (m)
  real(rl), intent(in)     :: y        ! y (m)
  real(rl), intent(inout)  :: z          ! Height (m)
  real(rl), intent(in)     :: hyam       ! A coefficient for hybrid-eta coordinate, at model level midpoint
  real(rl), intent(in)     :: hybm       ! B coefficient for hybrid-eta coordinate, at model level midpoint
  logical, intent(in)     :: hybrid_eta ! flag to indicate whether the hybrid sigma-p (eta) coordinate is used
	real(rl), intent(inout)  :: p          ! Pressure  (Pa)
	integer, intent(in)     :: zcoords    ! 0 or 1 see below
  real(rl), intent(in)    :: d           ! Width for Pert
	real(rl), intent(out)    :: u          ! x-dir wind (m s^-1)
  real(rl), intent(out)    :: v          ! y-dir wind (m s^-1)
  real(rl), intent(out)    :: w          ! Vertical Velocity (m s^-1)
  real(rl), intent(out)    :: t          ! Temperature (K)
  real(rl), intent(out)    :: t_mean     ! Temperature (K)
  real(rl), intent(out)    :: phis       ! Surface Geopotential (m^2 s^-2)
  real(rl), intent(out)    :: ps         ! Surface Pressure (Pa)
  real(rl), intent(out)    :: rho        ! density (kg m^-3)
  real(rl), intent(out)    :: rho_mean   ! density (kg m^-3)
  real(rl), intent(out)    :: q          ! Specific Humidity (kg/kg)

	! if zcoords = 1, then we use z and output z
	! if zcoords = 0, then we use p

!-----------------------------------------------------------------------
!     test case parameters
!-----------------------------------------------------------------------
	real(rl), parameter :: &
    u0      = 20.d0,        &	! Reference Velocity
    Tref     = 300.d0,       &	! Reference Temperature + Potential Temperature
    Pref     = 100000.d0,		&	! Reference PS
    ztop    = 10000.d0,     &	! Model Top
    yc    = 0.d0,         & ! y-loc of Pert Center
    xc    = 0.d0,         & ! x-loc of Pert Center
    delta_theta = 1.d0,     & ! Max Amplitude of Pert
    Lz      = 20000.d0, 		& ! Vertical Wavelength of Pert
    N       = 0.01d0,       & ! Brunt-Vaisala frequency
    N2      = N*N,          &	! Brunt-Vaisala frequency Squared
    bigG    = (g*g)/(N2*cp)   ! Constant

  real(rl) :: height           ! Model level height
  real(rl) :: r2, s							! Shape of perturbation
  real(rl) :: Ts 							! Surface temperature
  real(rl) :: t_pert           ! temperature perturbation
  real(rl) :: theta_pert       ! Pot-temp perturbation

!-----------------------------------------------------------------------
!    THE VELOCITIES
!-----------------------------------------------------------------------

	! Zonal Velocity

	u = u0

	! Meridional Velocity

	v = 0.d0

	! Vertical Velocity = Vertical Pressure Velocity = 0

	w = 0.d0

!-----------------------------------------------------------------------
!    PHIS (surface geopotential)
!-----------------------------------------------------------------------

	phis = 0.d0

!-----------------------------------------------------------------------
!    SURFACE TEMPERATURE
!-----------------------------------------------------------------------
	Ts = Tref

!-----------------------------------------------------------------------
!    PS (surface pressure)
!-----------------------------------------------------------------------
	ps = Pref

!-----------------------------------------------------------------------
!    HEIGHT AND PRESSURE AND MEAN TEMPERATURE
!-----------------------------------------------------------------------
	if (zcoords .eq. 1) then

		height = z
		p = ps*( (bigG/Ts)*exp(-N2*height/g)+1.d0 - (bigG/Ts)  )**(cp/Rd)

	else

    if (hybrid_eta) p = hyam*p0 + hybm*ps
		height = (-g/N2)*log( (Ts/bigG)*( (p/ps)**(Rd/cp) - 1.d0  ) + 1.d0 )
    z      = height

	endif

	t_mean = bigG*(1.d0 - exp(N2*height/g))+ Ts*exp(N2*height/g)

!-----------------------------------------------------------------------
!    rho (density), unperturbed using the background temperature t_mean
!-----------------------------------------------------------------------

	rho_mean = p/(Rd*t_mean)

!-----------------------------------------------------------------------
!    POTENTIAL TEMPERATURE PERTURBATION,
!    here: converted to temperature and added to the temperature field
!    models with a prognostic potential temperature field can utilize
!    the potential temperature perturbation theta_pert directly and add it
!    to the background theta field (not included here)
!-----------------------------------------------------------------------

  if (planar_slice .eqv. .true.) then
    r2  = (x-xc)**2
  else
	r2  = (x-xc)**2 + (y-yc)**2
end if

	s = (d**2)/(d**2 + r2)

	theta_pert = delta_theta*s*sin(2.d0*DD_PI*height/Lz)

	t_pert = theta_pert*(p/p0)**(Rd/cp)
	t = t_mean + t_pert

  rho = p/(Rd*t)

!-----------------------------------------------------------------------
!     initialize Q, set to zero
!-----------------------------------------------------------------------

	q = 0.d0

END SUBROUTINE gravity_wave


subroutine get_xycoordinates(x,y,hyam,hybm, i,j,k,elem,hvcoord)

  ! get x,y, vertical coords at node(i,j,k)

  real(rl),         intent(out):: x,y,hyam,hybm
  integer,          intent(in) :: i,j,k
  type(element_t),  intent(in) :: elem
  type(hvcoord_t),  intent(in) :: hvcoord

  ! get horizontal coordinates at column i,j
  x  = elem%spherep(i,j)%lon
  y  = elem%spherep(i,j)%lat

  ! get hybrid coeffiecients at midpoint of vertical level k
  hyam = hvcoord%hyam(k)
  hybm = hvcoord%hybm(k)

end subroutine





end module dry_planar_tests
#endif
