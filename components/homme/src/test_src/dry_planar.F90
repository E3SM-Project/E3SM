
#ifndef CAM
#include "config.h"

module dry_planar_tests

  use element_mod,          only: element_t
  use hybrid_mod,           only: hybrid_t
  use hybvcoord_mod,        only: hvcoord_t, set_layer_locations
  use parallel_mod,         only: abortmp
  use element_ops,          only: set_state, set_state_i, tests_finalize
  use eos,                  only: phi_from_eos
  use kinds,                only: rl=>real_kind, iulog
  use physical_constants,   only : dd_pi, rgas, rwater_vapor
  use dcmip12_wrapper, only : set_tracers, get_evenly_spaced_z, get_evenly_spaced_p, set_hybrid_coefficients, pressure_thickness
  use dimensions_mod,       only: np, nlev, nlevp, qsize, qsize_d, nelemd
  use element_state,        only: nt=>timelevels
  use control_mod,          only: planar_slice

  implicit none

!OG take from physical const mod instead

  real(rl), parameter :: Rd 	= 287.0d0,	&	! Ideal gas const dry air (J kg^-1 K^1)
				                        g	= 9.80616d0,	&	! Gravity (m s^2)
				                        cp	= 1004.5d0,	&	! Specific heat capacity (J kg^-1 K^1)
                                p0	= 100000.d0, &! reference pressure (Pa)
                                kappa   = Rd/cp

  real(rl), parameter :: bubble_const1=3.8, bubble_const2=17.27, bubble_const3=273.0, bubble_const4=36.0


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

  call dry_bubble_init(elem,hybrid,hvcoord,nets,nete,.25d0,0.d0)

end subroutine planar_rising_bubble_init


subroutine dry_bubble_init(elem,hybrid,hvcoord,nets,nete,d,f)

  use control_mod, only: bubble_T0, bubble_dT, bubble_xycenter, bubble_zcenter, bubble_ztop, &
                         bubble_xyradius,bubble_zradius, bubble_cosine, &
                         bubble_moist, bubble_moist_dq, &
                         Lx, Ly, Sx, Sy

  type(element_t),    intent(inout), target :: elem(:)                  ! element array
  type(hybrid_t),     intent(in)            :: hybrid                   ! hybrid parallel structure
  type(hvcoord_t),    intent(inout)         :: hvcoord                  ! hybrid vertical coordinates
  integer,            intent(in)            :: nets,nete                ! start, end element index
  real(rl),           intent(in)            :: d                        ! radius of perturbation
  real(rl),           intent(in)            :: f                        ! (const) Coriolis force

  integer :: i,j,k,ie,ii
  real(rl):: x,y,one,two,offset
  real(rl):: pi(nlevp), dpm(nlev), th0(nlevp), th0m(nlev), ai(nlevp), bi(nlevp), rr, &
             qi_s(nlevp), qm_s(nlev), Ti(nlevp), qi(nlevp), ri(nlevp)

#ifdef MODEL_THETA_L

  if (qsize < 3 .and. bubble_moist) then
    call abortmp('planar moist bubble requires at least 3 tracers')
  endif

  if (hybrid%masterthread) then
     write(iulog,*) 'initializing hot bubble with'
     print *, 'Lx, Ly =', Lx, Ly
     print *, 'Sx, Sy =', Sx, Sy
     print *, 'bubble_T0',  bubble_T0
     print *, 'bubble_dT', bubble_dT
     print *, 'bubble_xycenter', bubble_xycenter
     print *, 'bubble_zcenter', bubble_zcenter
     print *, 'bubble_ztop', bubble_ztop
     print *, 'bubble_xyradius', bubble_xyradius
     print *, 'bubble_zradius', bubble_zradius
     print *, 'bubble_cosine', bubble_cosine
     print *, 'bubble_moist', bubble_moist
     print *, 'bubble_moist_dq', bubble_moist_dq
  endif

  one = 1.0
  two = 2.0

  !evenly spaced with reversed indexing, just like in homme eta coord
  call get_evenly_spaced_z(zi,zm, 0.0_rl,bubble_ztop)

  !get pressure on interfaces
  do k=1,nlevp
    !set interface pressure from const lapse rate, does not depend on x,y, or _dT
    !this formula assumes ideal hy balance and ideal thermodynamics, not exactly consistent with
    !the dycore, that is why we reset zi from EOS below
    Ti(k) = bubble_T0 - zi(k)*g/cp
    pi(k) = p0*( Ti(k)/bubble_T0  )**(one/kappa)
  enddo
  do k=1,nlev
    dpm(k)=(pi(k+1)-pi(k))
  enddo

  !create hybrid coords now from pi
  !note that this code should not depend on partitioning
  !it depends on ref pressure profile, whcih is init-ed in the same way for all elements/gll

#if 0
  !almost pure pressure
  ai(:) = 0.0; bi(:) = 0.0
  ai(1) = pi(1)/p0; bi(nlevp) = one

  do k=2,nlevp
    bi(k) = pi(k)/pi(nlevp)
  enddo
#else
  !old version, hybrid
  ai(:) = 0.0; bi(:) = 0.0
  ai(1) = pi(1)/p0;  bi(nlevp) = one

  do k=2,nlev
    bi(k) = 1.0 - zi(k)/zi(1)
    !restore ai frop given pressure
    ai(k)=(pi(k)-bi(k)*pi(nlevp))/p0
  enddo
#endif

  hvcoord%hyai = ai; hvcoord%hybi = bi

  !set : hyam hybm 
  hvcoord%hyam = 0.5_rl *(ai(2:nlev+1) + ai(1:nlev))
  hvcoord%hybm = 0.5_rl *(bi(2:nlev+1) + bi(1:nlev))

  !  call set_layer_locations: sets  etam, etai, dp0, checks that Am=ai/2+ai/2
  call set_layer_locations(hvcoord, .true., hybrid%masterthread)

  !build specific humidity at saturation qs as in dcmip2016-kessler
  do k=1,nlevp
     qi_s(k) = bubble_const1 / pi(k) * exp( bubble_const2 * (Ti(k) - bubble_const3) / ( Ti(k) - bubble_const4 ) )
!CE function     
!   real tc = temp - 273.15;
!   return 610.94 * exp( 17.625*tc / (243.04+tc) );
!   qi_s(k)  = 610.94 * exp( 17.625*(Ti(k) - 273.15) / (Ti(k)-273.15 + 243.04  ) )
  enddo

  do k=1,nlev
    qm_s(k)=(qi_s(k+1)+qi_s(k)) / two
  enddo

  !reset z to the discrete hydro balance
  !if one wants to keep z levels exactly evenly spaced instead
  !they would have to compute p (hydro) and dp3d from the discrete EOS
  !via iteration (?) 
  !whole init-ed element is needed to call phi from EOS
  do j=1,np; do i=1,np
    elem(1)%state%phis(i,j) = 0.0
    elem(1)%state%dp3d(i,j,1:nlev,1) = dpm(1:nlev)
    elem(1)%state%vtheta_dp(i,j,1:nlev,1) = bubble_T0*dpm(1:nlev)
    elem(1)%state%Q(i,j,1:nlev,1) = qm_s(1:nlev)
    elem(1)%state%qdp(i,j,1:nlev,1,1) = dpm(1:nlev) * qm_s(1:nlev)
  enddo; enddo

  call phi_from_eos(hvcoord,elem(1)%state%phis,elem(1)%state%vtheta_dp(:,:,:,1),&
       elem(1)%state%dp3d(:,:,:,1),elem(1)%state%phinh_i(:,:,:,1))

  !reset zi array
  zi(1:nlevp) = elem(1)%state%phinh_i(1,1,1:nlevp,1)/g !<- i,j,k,tl
  !if needed,
  !reset zm to be an average as it is consistent with Taylor2020 eqn (30)
  !but zm is not used below, so, not done


  ! set initial conditions
  do ie = nets,nete
    do j=1,np; do i=1,np

      ! get horizontal coordinates at column i,j
      x  = elem(ie)%spherep(i,j)%lon; y  = elem(ie)%spherep(i,j)%lat

      do k=1,nlevp

        if (planar_slice .eqv. .true.) then
        !no y dependence
          rr =sqrt( (x-bubble_xycenter)   *(x-bubble_xycenter) / bubble_xyradius / bubble_xyradius + &
                    (zi(k)-bubble_zcenter)*(zi(k)-bubble_zcenter) / bubble_zradius / bubble_zradius    )
        else
          rr =sqrt( (x-bubble_xycenter)   *(x-bubble_xycenter) / bubble_xyradius / bubble_xyradius + &
                    (y-bubble_xycenter)   *(y-bubble_xycenter) / bubble_xyradius / bubble_xyradius + &
                    (zi(k)-bubble_zcenter)*(zi(k)-bubble_zcenter) / bubble_zradius / bubble_zradius    )
        endif
      
        !set pot. temperature on interfaces
        if ( rr < one ) then 

          qi(k) = bubble_moist_dq ! qi_s in many forms did not work

          if (bubble_cosine) then
            offset = cos(rr*dd_pi / two)
            th0(k) = bubble_T0 + bubble_dT * offset

            ! q = Rel Humidity * qs
            ! relative humidity = offset, or offset*, say, 0.9?
            qi(k) = offset * qi(k)
          else
            ! qi is set above to const
            th0(k) = bubble_T0 + bubble_dT
          endif

        else

          !set to reference profile
          th0(k) = bubble_T0
          qi(k) = qi_s(k)  

        endif

        !R_star(:,:,k) =(Rgas + (Rwater_vapor - Rgas)*Q(:,:,k))
        ri(k) = Rgas + (Rwater_vapor - Rgas)*qi(k)

      enddo ! k loop

      elem(ie)%state%ps_v(i,j,:) = pi(nlevp)
      elem(ie)%state%v = 0.0; elem(ie)%state%w_i= 0.0

      do k=1,nlev
        !set pottemp, dp, other state vars on midlevels
        !th0m(k)=(th0(k)+th0(k+1))/ two
        th0m(k)=(th0(k)*ri(k)+th0(k+1)*ri(k+1))/ two /rgas

        elem(ie)%state%dp3d(i,j,k,:)   = dpm(k)
      
        elem(ie)%state%phis(i,j) = 0.0
        elem(ie)%state%phinh_i(i,j,k,:) = zi(k)*g

        elem(ie)%state%vtheta_dp(i,j,k,:) = dpm(k)*th0m(k)

        if (bubble_moist) then
        elem(ie)%state%Q(i,j,k,1) =   ( qi(k) + qi(k+1) ) / two        
        elem(ie)%state%Qdp(i,j,k,1,:) = dpm(k) * elem(ie)%state%Q(i,j,k,1)
        end if
      enddo !k loop

    enddo; enddo !i,j loop
  enddo !ie loop

!    real (kind=real_kind) :: Q   (np,np,nlev,qsize_d)                 ! Tracer concentration               6
!    real (kind=real_kind) :: Qdp (np,np,nlev,qsize_d,2) 
  if (bubble_moist) then
     ii=2
  else 
     ii=1
  endif

  do ie = nets,nete
     elem(ie)%fcor(:,:) = f

     !all but vapor
     elem(ie)%state%Q(:,:,:,ii:qsize) = 0.0
     elem(ie)%state%Qdp(:,:,:,ii:qsize,:) = 0.0

     !sets hydro phi from theta and pressure, checks for hydrostatic balance after that, saves a state
     !call tests_finalize(elem(ie),hvcoord)
  enddo

#else
  !not THETA_L
  !abort, no bubble for preqx
  call abortmp('planar rising bubble does not work with preqx')
#endif

end subroutine dry_bubble_init



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
    u0      = 20.d0,        &! Reference Velocity
    Tref     = 300.d0,      &! Reference Temperature + Potential Temperature
    Pref     = 100000.d0,   &! Reference PS
    ztop    = 10000.d0,     &! Model Top
    yc    = 0.d0,         & ! y-loc of Pert Center
    xc    = 0.d0,         & ! x-loc of Pert Center
    delta_theta = 1.d0,     & ! Max Amplitude of Pert
    Lz      = 20000.d0,   & ! Vertical Wavelength of Pert
    N       = 0.01d0,       & ! Brunt-Vaisala frequency
    N2      = N*N,          &! Brunt-Vaisala frequency Squared
    bigG    = (g*g)/(N2*cp)   ! Constant

  real(rl) :: height           ! Model level height
  real(rl) :: r2, s            ! Shape of perturbation
  real(rl) :: Ts               ! Surface temperature
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
    z = height

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
