#ifndef CAM
#include "config.h"

module dcmip_tests

use control_mod,          only: test_case
use dcmip2012_test1_2_3,  only: test2_steady_state_mountain, test2_schaer_mountain ,test3_gravity_wave
use dimensions_mod,       only: np, nlev, qsize, nlevp, qsize_d, nelemd
use element_mod,          only: element_t, elem_state_t, derived_state_t, nt=>timelevels
use hybrid_mod,           only: hybrid_t
use hybvcoord_mod,        only: hvcoord_t, set_layer_locations
use kinds,                only: rl=>real_kind, iulog
use parallel_mod,         only: abortmp

implicit none

! physical constants used by dcmip2012 test 3-1
real(rl), parameter ::              &
  g       = 9.80616,                & ! grav const
  a       = 6371229.0,              & ! earth radius in meters
  Rd      = 287.0,                  & ! dry gas const
  cp      = 1004.5,                 & ! heat capacity const pressure
  kappa   = Rd/cp,                  &
  pi      = 3.141592654,            &
  p0      = 100000.0                  ! reference pressure

contains

!_____________________________________________________________________
function evenly_spaced_z(z_b,z_t) result(zi)

  real(rl), intent(in) :: z_b,z_t   ! top and bottom coordinates
  real(rl)             :: zi(nlevp) ! z at interfaces
  integer :: k

  forall(k=1:nlevp) zi(k) = z_t-(k-1)*(z_t-z_b)/(nlevp-1)

end function

!_____________________________________________________________________
subroutine set_hybrid_coefficients(hv, hybrid, eta_t, c)

  type(hvcoord_t),    intent(inout) :: hv       ! hybrid vertical coordinate stucture
  type(hybrid_t),     intent(in)    :: hybrid   ! hybrid parallal structure
  real(rl),           intent(in)    :: eta_t    ! top eta level
  real(rl),           intent(in)    :: c        ! exponent

  real(rl)  :: eta_c, tmp
  integer   :: k

  ! place cutoff halfway between bottom and top eta coordiantes
  eta_c = hv%etai(nlev/2)

  do k=1,nlevp
    tmp        = max( (hv%etai(k)-eta_c)/(1.0-eta_c), 0.0_rl)
    hv%hybi(k) = tmp**c
    hv%hyai(k) = hv%etai(k) - hv%hybi(k)
    if(hybrid%par%masterproc) print *,k,': etai = ',hv%etai(k),' Ai = ',hv%hyai(k),' Bi = ',hv%hybi(k);
  enddo

  hv%hyam = 0.5_rl *(hv%hyai(2:nlev+1) + hv%hyai(1:nlev))
  hv%hybm = 0.5_rl *(hv%hybi(2:nlev+1) + hv%hybi(1:nlev))
  hv%etam = hv%hyam + hv%hybm

end subroutine

!_____________________________________________________________________
subroutine get_coordinates(lat,lon,hyam,hybm, i,j,k,elem,hvcoord)

  real(rl),         intent(out):: lon,lat,hyam,hybm
  integer,          intent(in) :: i,j,k
  type(element_t),  intent(in) :: elem
  type(hvcoord_t),  intent(in) :: hvcoord

  ! get horizontal coordinates at column i,j
  lon  = elem%spherep(i,j)%lon
  lat  = elem%spherep(i,j)%lat

  ! get hybrid coeffiecients at midpoint of vertical level k
  hyam = hvcoord%hyam(k)
  hybm = hvcoord%hybm(k)

end subroutine

!_____________________________________________________________________
subroutine set_element_state(u,v,T,ps,phis,q, i,j,k,elem,hvcoord)

  real(rl),         intent(in)    :: u,v,T,ps,phis,q
  integer,          intent(in)    :: i,j,k
  type(element_t),  intent(inout) :: elem
  type(hvcoord_t),  intent(in)    :: hvcoord

  real(rl):: p_np1, p_n, dp                                             ! pressure thickness

  ! get pressure level thickness
  p_np1 = hvcoord%hyai(k+1)*p0 + hvcoord%hybi(k+1)*ps
  p_n   = hvcoord%hyai(k  )*p0 + hvcoord%hybi(k  )*ps
  dp = p_np1-p_n

  ! set prognostic state variables at level midpoints
  elem%state%v   (i,j,1,k,:) = u
  elem%state%v   (i,j,2,k,:) = v
  elem%state%T   (i,j,k,:)   = T
  elem%state%dp3d(i,j,k,:)   = dp
  elem%state%ps_v(i,j,:)     = ps
  elem%state%phis(i,j)       = phis
  elem%state%Q   (i,j,k,1)   = q
  elem%state%Qdp (i,j,k,1,:) = q*dp

  ! set some diagnostic variables
  elem%derived%omega_p(i,j,k)      = 0.0d0
  elem%derived%eta_dot_dpdn(i,j,k) = 0.0d0
  elem%derived%dp(i,j,k)           = dp

end subroutine

!_____________________________________________________________________
subroutine dcmip2012_test2_0(elem,hybrid,hvcoord,nets,nete)

  ! steady state atmosphere with orography

  type(element_t),    intent(inout), target :: elem(:)                  ! element array
  type(hybrid_t),     intent(in)            :: hybrid                   ! hybrid parallel structure
  type(hvcoord_t),    intent(inout)         :: hvcoord                  ! hybrid vertical coordinates
  integer,            intent(in)            :: nets,nete                ! start, end element index

  integer,  parameter :: zcoords = 0                                    ! we are not using z coords
  logical,  parameter :: use_eta = .true.                               ! we are using hybrid eta coords
   real(8), parameter ::      &
      T0      = 300.d0,       &                                         ! temperature (K)
      gamma   = 0.0065d0,     &                                         ! temperature lapse rate (K/m)
      ztop    = 12000.d0,     &                                         ! model top (m)
      exponent= g/(Rd*gamma)

  integer :: i,j,k,ie                                                   ! loop indices
  real(rl):: lon,lat,hyam,hybm                                          ! pointwise coordiantes
  real(rl):: p,z,phis,u,v,w,T,T_mean,phis_ps,ps,rho,rho_mean,q          ! pointwise field values
  real(rl):: zi(nlevp)                                                  ! z coords of interfaces

  if (hybrid%masterthread) write(iulog,*) 'initializing dcmip2012 test 2-0: steady state atmosphere with orography'

  ! set analytic vertical coordinates
  zi            = evenly_spaced_z(0.0_rl, ztop)                         ! get evenly spaced z levels
  hvcoord%etai  = (1.d0 - gamma/T0*zi)**exponent                        ! set eta levels from z in orography-free region
  call set_hybrid_coefficients(hvcoord,hybrid,  hvcoord%etai(1), 1.0_rl)! set hybrid A and B from eta levels
  call set_layer_locations(hvcoord, .true., hybrid%masterthread)

  ! set initial conditions
  do ie = nets,nete; do k=1,nlev; do j=1,np; do i=1,np
    call get_coordinates(lat,lon,hyam,hybm, i,j,k,elem(ie),hvcoord)
    call test2_steady_state_mountain(lon,lat,p,z,zcoords,use_eta,hyam,hybm,u,v,w,T,phis,ps,rho,q)
    call set_element_state(u,v,T,ps,phis,q, i,j,k,elem(ie),hvcoord)
  enddo; enddo; enddo; enddo

end subroutine

!_____________________________________________________________________
subroutine dcmip2012_test2_x(elem,hybrid,hvcoord,nets,nete,shear)

  ! nonhydrostatic orographic waves (with or without shear)

  type(element_t),    intent(inout), target :: elem(:)                  ! element array
  type(hybrid_t),     intent(in)            :: hybrid                   ! hybrid parallel structure
  type(hvcoord_t),    intent(inout)         :: hvcoord                  ! hybrid vertical coordinates
  integer,            intent(in)            :: nets,nete                ! start, end element index
  integer,            intent(in)            :: shear                    ! flag: 1=shear 0=no shear

  integer,  parameter :: zcoords = 0                                    ! we are not using z coords
  logical,  parameter :: use_eta = .true.                               ! we are using hybrid eta coords
  real(8),  parameter ::   &
      Teq     = 300.d0,    &                                            ! temperature at equator
      ztop    = 30000.d0,	 &                                            ! model top (m)
      H       = Rd*Teq/g                                                ! characteristic height scale

  integer :: i,j,k,ie                                                   ! loop indices
  real(rl):: lon,lat,hyam,hybm                                          ! pointwise coordiantes
  real(rl):: p,z,phis,u,v,w,T,T_mean,phis_ps,ps,rho,rho_mean,q          ! pointwise field values
  real(rl):: zi(nlevp)                                                  ! z coords of interfaces

  if (hybrid%masterthread) write(iulog,*) 'initializing dcmip2012 test 2-x: steady state atmosphere with orography'

  ! set analytic vertical coordinates
  zi            = evenly_spaced_z(0.0_rl, ztop)                         ! get evenly spaced z levels
  hvcoord%etai  = exp(-zi/H)                                            ! set eta levels from z in orography-free region
  call set_hybrid_coefficients(hvcoord,hybrid,  hvcoord%etai(1), 1.0_rl)! set hybrid A and B from eta levels
  call set_layer_locations(hvcoord, .true., hybrid%masterthread)

  ! set initial conditions
  do ie = nets,nete; do k=1,nlev; do j=1,np; do i=1,np
    call get_coordinates(lat,lon,hyam,hybm, i,j,k,elem(ie),hvcoord)
    call test2_schaer_mountain(lon,lat,p,z,zcoords,use_eta,hyam,hybm,shear,u,v,w,T,phis,ps,rho,q)
    call set_element_state(u,v,T,ps,phis,q, i,j,k,elem(ie),hvcoord)
  enddo; enddo; enddo; enddo

end subroutine

!_____________________________________________________________________
subroutine dcmip2012_test3(elem,hybrid,hvcoord,nets,nete)

  ! nonhydrostatic gravity waves

  type(element_t),    intent(inout), target :: elem(:)                  ! element array
  type(hybrid_t),     intent(in)            :: hybrid                   ! hybrid parallel structure
  type(hvcoord_t),    intent(inout)         :: hvcoord                  ! hybrid vertical coordinates
  integer,            intent(in)            :: nets,nete                ! start, end element index

  integer,  parameter :: zcoords = 0                                    ! we are not using z coords
  logical,  parameter :: use_eta = .true.                               ! we are using hybrid eta coords

  real(rl), parameter ::    &                                           ! parameters needed to get eta from z
    T0      = 300.d0,       &	! temperature (k)
    ztop    = 10000.d0,     & ! model top (m)
    N       = 0.01d0,       & ! Brunt-Vaisala frequency
    bigG    = (g*g)/(N*N*Cp)  ! temperature, isothermal

  integer :: i,j,k,ie                                                   ! loop indices
  real(rl):: lon,lat,hyam,hybm                                          ! pointwise coordiantes
  real(rl):: p,z,phis,u,v,w,T,T_mean,phis_ps,ps,rho,rho_mean,q          ! pointwise field values
  real(rl):: p_np1, p_n, dp                                             ! pressure thickness
  real(rl):: zi(nlevp)                                                  ! z coords of interfaces

  if (hybrid%masterthread) write(iulog,*) 'initializing dcmip2012 test 3-0: nonhydrostatic gravity waves'

  ! set analytic vertical coordinates
  zi            = evenly_spaced_z(0.0_rl, ztop)                         ! get evenly spaced z levels
  hvcoord%etai  = ( (bigG/T0)*(exp(-zi*N*N/g) -1 )+1 ) **(1.0/kappa)    ! set eta levels from z at equator
  call set_hybrid_coefficients(hvcoord,hybrid,  hvcoord%etai(1), 1.0_rl)! set hybrid A and B from eta levels
  call set_layer_locations(hvcoord, .true., hybrid%masterthread)

  ! set initial conditions
  do ie = nets,nete; do k=1,nlev; do j=1,np; do i=1,np
    call get_coordinates(lat,lon,hyam,hybm, i,j,k,elem(ie),hvcoord)
    call test3_gravity_wave(lon,lat,p,z,zcoords,use_eta,hyam,hybm,u,v,w,T,T_mean,phis,ps,rho,rho_mean,q)
    call set_element_state(u,v,T,ps,phis,q, i,j,k,elem(ie),hvcoord)
  enddo; enddo; enddo; enddo

end subroutine

end module dcmip_tests
#endif
