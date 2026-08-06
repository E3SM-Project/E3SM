
#ifndef CAM
#include "config.h"

module dry_planar_tests

  use element_mod,          only: element_t
  use hybrid_mod,           only: hybrid_t
  use hybvcoord_mod,        only: hvcoord_t, set_layer_locations
  use parallel_mod,         only: abortmp
  use element_ops,          only: set_state, set_state_i, tests_finalize
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

!consts for kessler-defined qsat
  real(rl), parameter :: bubble_const1=3.8, bubble_const2=17.27, bubble_const3=273.0, bubble_const4=36.0
!consts for RJ-defined qsat
  real(rl), parameter :: bubble_t0_const=273.16, bubble_epsilo=Rgas/Rwater_vapor, bubble_e0=610.78
!this one is used in dcmip as 2.5e6 instead of the one in cam, 2.501e6
  real(rl), parameter :: bubble_latvap=2.5e6

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

  !last input is coriolis parameter
  call bubble_init(elem,hybrid,hvcoord,nets,nete,0.d0)

end subroutine planar_rising_bubble_init


! was not tested with preqx

! A note on relation T = \Pi * theta where \Pi = (p/p0)^dry_kappa that is used below:
! depending on assumptions for Gibbs potential, definition of \Pi with dry kapa
! may not hold. It is correct for theta-l and is used (though may be
! inconsistent) in preqx in other places.
subroutine bubble_init(elem,hybrid,hvcoord,nets,nete,f)

  use control_mod, only: bubble_T0, bubble_dT, bubble_xycenter, bubble_zcenter, bubble_ztop, &
                         bubble_xyradius,bubble_zradius, bubble_cosine, &
                         bubble_moist, bubble_moist_drh, bubble_prec_type, bubble_rh_background
  use physical_constants, only: Lx, Ly, Sx, Sy
  use element_ops, only: set_elem_state

  type(element_t),    intent(inout), target :: elem(:)                  ! element array
  type(hybrid_t),     intent(in)            :: hybrid                   ! hybrid parallel structure
  type(hvcoord_t),    intent(inout)         :: hvcoord                  ! hybrid vertical coordinates
  integer,            intent(in)            :: nets,nete                ! start, end element index
  real(rl),           intent(in)            :: f                        ! (const) Coriolis force

  integer :: i,j,k,ie,ii
  real(rl):: x,y,offset
  real(rl):: pi(nlevp), pm(nlev), dpm(nlev), th0(nlevp), th0m(nlev), ai(nlevp), bi(nlevp), rr, &
             qi_s(nlevp), qm_s(nlev), Ti(nlevp), Tm(nlev), qi(nlevp)

  real(rl):: zero_mid_init(np,np,nlev),zero_int_init(np,np,nlevp), &
             dp_init(np,np,nlev),ps_init(np,np), &
             phis_init(np,np,nlevp),t_init(np,np,nlev),p_init(np,np,nlev), &
             zi_init(np,np,nlevp), zm_init(np,np,nlev), &
             rh
             
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
     print *, 'bubble_moist_drh', bubble_moist_drh
     print *, 'bubble_rh_background', bubble_rh_background
     print *, 'bubble_prec_type (0 is Kessler (default), 1 is RJ)', bubble_prec_type
  endif

  !for the background state
  !evenly spaced with reversed indexing, just like in homme eta coord
  call get_evenly_spaced_z(zi,zm,0.0_rl,bubble_ztop)

  !for the background state
  do k=1,nlevp
    Ti(k) = bubble_T0 - zi(k)*g/cp
    pi(k) = p0*( Ti(k)/bubble_T0  )**(1.0/kappa)
  enddo
  do k=1,nlev
    Tm(k) = bubble_T0 - zm(k)*g/cp
    pm(k) = p0*( Tm(k)/bubble_T0  )**(1.0/kappa)
  enddo
  ! this T, z, pressure, and q1 (set below) together do not obey 
  ! homme's EOS, meaning if one to compute phi from EOS, it won't match z
  ! one could have zi replaced with output from phi_from_eos(),
  ! but it won't satisfy model-independent setup.
  ! instead, we use uniform zi as above.

  !create hybrid coords now from pi
  !note that this code should not depend on partitioning
  !it depends on ref pressure profile, whcih is init-ed in the same way for all elements/gll
#if 0
  !almost sigma
  ai(:) = 0.0; bi(:) = 0.0
  ai(1) = pi(1)/p0; bi(nlevp) = one

  do k=2,nlevp
    bi(k) = pi(k)/pi(nlevp)
  enddo
#else
  !old version, hybrid
  ai(:) = 0.0; bi(:) = 0.0
  ai(1) = pi(1)/p0;  bi(nlevp) = 1.0

  do k=2,nlev
    bi(k) = 1.0 - zi(k)/zi(1)
    !restore ai frop given pressure
    ai(k)=(pi(k)-bi(k)*pi(nlevp))/p0
  enddo
#endif

  hvcoord%hyai = ai; hvcoord%hybi = bi

  do k = 1,nlev
    dpm(k) = pressure_thickness(pi(nlevp),k,hvcoord)
  end do

  !set : hyam hybm 
  hvcoord%hyam = 0.5 *(ai(2:nlev+1) + ai(1:nlev))
  hvcoord%hybm = 0.5 *(bi(2:nlev+1) + bi(1:nlev))

  !call set_layer_locations: sets  etam, etai, dp0, checks that Am=ai/2+ai/2
  call set_layer_locations(hvcoord, .true., hybrid%masterthread)

  !for the background state
  if (bubble_moist) then
    !build specific humidity at saturation qs as in dcmip2016-kessler

    if(bubble_prec_type == 0) then
    !kessler
    do k=1,nlevp
       qi_s(k) = bubble_const1 / pi(k) * exp( bubble_const2 * (Ti(k) - bubble_const3) / ( Ti(k) - bubble_const4 ) )
    enddo
    elseif(bubble_prec_type == 1) then
    !RJ physics
    do k=1,nlevp
       !dcmip line
       !qsat = epsilo * e0 / p(k) * exp(-(latvap/rh2o) * ((one/t(k))-(one/T0)))
       qi_s(k) = bubble_epsilo * bubble_e0 / pi(k) * &
                 exp(-(bubble_latvap/Rwater_vapor) * ((1.0/Ti(k))-(1.0/bubble_t0_const)))
    enddo
    else
    call abortmp('planar moist bubble bubble_prec_type should be 0 or 1')
    endif   

  else
     qi_s(1:nlevp) = 0.0
  endif

  do k=1,nlev
    qm_s(k)=(qi_s(k+1)+qi_s(k)) / 2.0
  enddo

  !set some of element-size arrays
  zero_mid_init(:,:,:) = 0.0
  zero_int_init(:,:,:) = 0.0
  do j=1,np; do i=1,np;
    zi_init(i,j,1:nlevp) = zi(1:nlevp)
    zm_init(i,j,1:nlev)  = zm(1:nlev)
    ps_init(i,j) = pi(nlevp)
    p_init(i,j,1:nlev) = pm(1:nlev)
    dp_init(i,j,1:nlev) = dpm(1:nlev)
  enddo; enddo

  ! set the rest of conditions for each element
  do ie = nets,nete
    do j=1,np; do i=1,np

      ! get horizontal coordinates at column i,j
      x  = elem(ie)%spherep(i,j)%lon; y  = elem(ie)%spherep(i,j)%lat

      do k=1,nlevp

        !intermediate value
        rr=(x-bubble_xycenter)   *(x-bubble_xycenter) / bubble_xyradius / bubble_xyradius + &
           (zi(k)-bubble_zcenter)*(zi(k)-bubble_zcenter) / bubble_zradius / bubble_zradius 

        if (planar_slice .eqv. .true.) then
          !no y dependence
          rr=sqrt(rr)
        else
          rr =sqrt( rr + &
                   (y-bubble_xycenter) * (y-bubble_xycenter) / bubble_xyradius / bubble_xyradius )
        endif
      
        rh=bubble_rh_background
        if ( rr < 1.0 ) then
          if (bubble_cosine) then
            offset = cos(rr*dd_pi / 2.0 )
            th0(k) = bubble_T0 + bubble_dT * offset
            rh = rh + bubble_moist_drh * offset
          else
            !0/1 nonsmooth function
            th0(k) = bubble_T0 + bubble_dT
            rh = rh + bubble_moist_drh
          endif
        else
          !set to reference profile
          th0(k) = bubble_T0
        endif
        qi(k) = rh * qi_s(k)
      enddo ! k loop

      !set theta on midlevels and then T from theta, exner
      th0m(1:nlev) = (th0(1:nlev) + th0(2:nlevp) ) / 2.0
      t_init(i,j,1:nlev) = th0m(1:nlev) * ( pm(1:nlev)/p0 )**kappa
   
      !set Q before set_elem_state   
      if (bubble_moist) then
        do k=1,nlev
          elem(ie)%state%Q(i,j,k,1) =   ( qi(k) + qi(k+1) ) / 2.0
        enddo        
      else
        elem(ie)%state%Q(i,j,:,1) =   0.0
      end if

      !call set_elem_state(u,v,w,w_i,T,ps,phis,p,dp,zm,zi,g,elem,n0,n1,ntQ)
      call set_elem_state(zero_mid_init,zero_mid_init,zero_mid_init, &
                          zero_int_init,t_init,ps_init,zero_mid_init(:,:,1), &
                          p_init,dp_init,zm_init,zi_init,g,elem(ie),1,3,-1)

    enddo; enddo !i,j loop
  enddo !ie loop

  !indexing of Q, Qdp
  !Q   (np,np,nlev,qsize_d)   
  !Qdp (np,np,nlev,qsize_d,2) 
  if (bubble_moist) then
     ii=2
  else 
     ii=1
  endif

  do ie = nets,nete
     elem(ie)%fcor(:,:) = f

     !all but vapor
     elem(ie)%state%Q(:,:,:,ii:qsize) = 0.0

     !sets hydro phi from (perturbed) theta and pressure, checks for hydrostatic balance after that, saves a state
     !call tests_finalize(elem(ie),hvcoord)
  enddo

end subroutine bubble_init



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
