#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module test_mod

use control_mod,    only: test_case, sub_case
use derivative_mod, only: derivative_t, gradient_sphere
use element_mod,    only: element_t
use element_state,  only: elem_state_t, nt=>timelevels
use hybrid_mod,     only: hybrid_t
use hybvcoord_mod,  only: hvcoord_t
use kinds,          only: rl => real_kind, iulog
use parallel_mod,   only: abortmp
use time_mod,       only: timelevel_t, time_at

! test case routines
use asp_tests,            only: asp_tracer, asp_baroclinic, asp_rossby, asp_mountain, asp_gravity_wave
use baroclinic_inst_mod,  only: binst_init_state, jw_baroclinic
use dcmip_tests,          only: dcmip2012_test1_1, dcmip2012_test1_2, dcmip2012_test1_3,&
                                dcmip2012_test2_0, dcmip2012_test2_x, dcmip2012_test3
use held_suarez_mod,      only: hs0_init_state

implicit none

contains

!_______________________________________________________________________
subroutine set_test_initial_conditions(elem, deriv, hybrid, hvcoord, tl, nets, nete)

  ! set initial conditions for HOMME stand-alone tests cases

  implicit none
  type(element_t),    intent(inout), target :: elem(:)                  ! element array
  type (derivative_t),intent(in)            :: deriv
  type(hybrid_t),     intent(in)            :: hybrid                   ! hybrid parallel structure
  type(hvcoord_t),    intent(inout)         :: hvcoord                  ! hybrid vertical coordinates
  type(timelevel_t),  intent(in)            :: tl                       ! time level sctructure
  integer,            intent(in)            :: nets,nete                ! start, end element index

  select case(test_case)

    case('asp_baroclinic');     call asp_baroclinic   (elem,hybrid,hvcoord,nets,nete)
    case('asp_gravity_wave');   call asp_gravity_wave (elem,hybrid,hvcoord,nets,nete,sub_case)
    case('asp_mountain');       call asp_mountain     (elem,hybrid,hvcoord,nets,nete)
    case('asp_rossby');         call asp_rossby       (elem,hybrid,hvcoord,nets,nete)
    case('asp_tracer');         call asp_tracer       (elem,hybrid,hvcoord,nets,nete)
    case('baroclinic');         call binst_init_state (elem,hybrid, nets, nete, hvcoord)
    case('dcmip2012_test1_1');  call dcmip2012_test1_1(elem,hybrid,hvcoord,nets,nete,0.0d0,1,nt)
    case('dcmip2012_test1_2');  call dcmip2012_test1_2(elem,hybrid,hvcoord,nets,nete,0.0d0,1,nt)
    case('dcmip2012_test1_3');  call dcmip2012_test1_3(elem,hybrid,hvcoord,nets,nete,0.0d0,1,nt,deriv)
    case('dcmip2012_test2_0');  call dcmip2012_test2_0(elem,hybrid,hvcoord,nets,nete)
    case('dcmip2012_test2_1');  call dcmip2012_test2_x(elem,hybrid,hvcoord,nets,nete,0)
    case('dcmip2012_test2_2');  call dcmip2012_test2_x(elem,hybrid,hvcoord,nets,nete,1)
    case('dcmip2012_test3');    call dcmip2012_test3  (elem,hybrid,hvcoord,nets,nete)
    case('held_suarez0');       call hs0_init_state   (elem,hybrid,hvcoord,nets,nete,300.0_rl)
    case('jw_baroclinic');      call jw_baroclinic    (elem,hybrid,hvcoord,nets,nete)
    case default;               call abortmp('unrecognized test case')

  endselect

end subroutine

!_______________________________________________________________________
subroutine set_test_prescribed_wind(elem, deriv, hybrid, hvcoord, dt, tl, nets, nete)

  implicit none
  type(element_t),    intent(inout), target :: elem(:)                  ! element array
  type (derivative_t),intent(in)            :: deriv
  type(hybrid_t),     intent(in)            :: hybrid                   ! hybrid parallel structure
  type(hvcoord_t),    intent(inout)         :: hvcoord                  ! hybrid vertical coordinates
  real(rl),           intent(in)            :: dt                       ! fixed timestep size
  type(timelevel_t),  intent(in)            :: tl                       ! time level structure
  integer,            intent(in)            :: nets,nete                ! start, end element index

  integer :: n0,np1, ie
  real(rl):: time
  type(elem_state_t), pointer :: s

  time = tl%nstep*dt
  n0   = tl%n0
  np1  = tl%np1

  do ie=nets,nete
    ! duplicate state from last timestep
    s => elem(ie)%state
    s%v   (:,:,:,:,np1) = s%v   (:,:,:,:,n0)
    s%T   (:,:,:,  np1) = s%T   (:,:,:,  n0)
    s%dp3d(:,:,:,  np1) = s%dp3d(:,:,:,  n0)
    s%ps_v(:,:,    np1) = s%ps_v(:,:,    n0)
  enddo

  ! set prescribed quantities
  select case(test_case)
    case('dcmip2012_test1_1'); call dcmip2012_test1_1(elem,hybrid,hvcoord,nets,nete,time,np1,np1)
    case('dcmip2012_test1_2'); call dcmip2012_test1_2(elem,hybrid,hvcoord,nets,nete,time,np1,np1)
    case('dcmip2012_test1_3'); call dcmip2012_test1_3(elem,hybrid,hvcoord,nets,nete,time,np1,np1,deriv)
  endselect

end subroutine

!_______________________________________________________________________
subroutine compute_test_forcing(elem,hybrid,hvcoord,n,n_tracer,dt,nets,nete)

  ! apply forcing terms produced by HOMME stand-alone tests

  use dcmip_tests, only:  dcmip2012_test2_x_forcing
  use held_suarez_mod, only: hs_forcing
  implicit none
  type(element_t),  intent(inout) :: elem(:)                            ! element array
  type(hybrid_t),   intent(in)    :: hybrid                             ! hybrid parallel structure
  type(hvcoord_t),  intent(in)    :: hvcoord
  real(kind=rl),    intent(in)    :: dt
  integer,          intent(in)    :: n,nets,nete,n_tracer

  integer :: ie,q

  do ie=nets,nete
    elem(ie)%derived%FT = 0
    elem(ie)%derived%FM = 0
    elem(ie)%derived%FQ = 0
  enddo

  ! get forcing from test case
  select case(test_case)
    case('dcmip2012_test2_1');  call dcmip2012_test2_x_forcing(elem, hybrid,hvcoord,nets,nete,n,dt)
    case('dcmip2012_test2_2');  call dcmip2012_test2_x_forcing(elem, hybrid,hvcoord,nets,nete,n,dt)
    case('held_suarez0');       
       do ie=nets,nete
          call hs_forcing(elem(ie),hvcoord,n,n_tracer,dt)
       enddo
  endselect

end subroutine

end module
