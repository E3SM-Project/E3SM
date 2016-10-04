#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module test_mod

use control_mod,    only: test_case, sub_case
use element_mod,    only: element_t
use fvm_control_volume_mod, only: fvm_struct
use hybrid_mod,     only: hybrid_t
use hybvcoord_mod,  only: hvcoord_t
use kinds,          only: real_kind, iulog
use parallel_mod,   only: abortmp
use time_mod,       only: timelevel_t

implicit none

contains

!_____________________________________________________________________
subroutine set_test_initial_conditions(elem, hybrid, hvcoord, tl, nets, nete, fvm)

  use asp_tests,            only: asp_tracer, asp_baroclinic, asp_rossby, asp_mountain, asp_gravity_wave, dcmip2_schar
  use baroclinic_inst_mod,  only: binst_init_state, jw_baroclinic
  use dcmip_tests,          only: dcmip2012_test2_0, dcmip2012_test2_x, dcmip2012_test3
  use held_suarez_mod,      only: hs0_init_state

  ! set initial conditions for HOMME stand-alone tests cases

  implicit none
  type(element_t),    intent(inout), target :: elem(:)                  ! element array
  type(hybrid_t),     intent(in)            :: hybrid                   ! hybrid parallel structure
  type(hvcoord_t),    intent(inout)         :: hvcoord                  ! hybrid vertical coordinates
  type (timelevel_t), intent(in)            :: tl                       ! time level sctructure
  integer,            intent(in)            :: nets,nete                ! start, end element index
  type(fvm_struct),   intent(inout)         :: fvm(:)                   ! finite volume structure

  select case(test_case)

    case('asp_baroclinic');     call asp_baroclinic   (elem, hybrid,hvcoord,nets,nete,fvm)
    case('asp_gravity_wave');   call asp_gravity_wave (elem, hybrid,hvcoord,nets,nete, sub_case)
    case('asp_mountain');       call asp_mountain     (elem, hybrid,hvcoord,nets,nete)
    case('asp_rossby');         call asp_rossby       (elem, hybrid,hvcoord,nets,nete)
    case('asp_tracer');         call asp_tracer       (elem, hybrid,hvcoord,nets,nete)
    case('baroclinic');         call binst_init_state (elem, hybrid, nets, nete, hvcoord)
    case('dcmip2012_test2_0');  call dcmip2012_test2_0(elem, hybrid,hvcoord,nets,nete)
    case('dcmip2012_test2_1');  call dcmip2012_test2_x(elem, hybrid,hvcoord,nets,nete,0)
    case('dcmip2012_test2_2');  call dcmip2012_test2_x(elem, hybrid,hvcoord,nets,nete,1)
    case('dcmip2012_test3');    call dcmip2012_test3  (elem, hybrid,hvcoord,nets,nete)
    case('held_suarez0');       call hs0_init_state   (elem, hvcoord,nets,nete,300.0_real_kind)
    case('jw_baroclinic');      call jw_baroclinic    (elem, hybrid,hvcoord,nets,nete)
    case default;               call abortmp('unrecognized test case')

  endselect

end subroutine

#if 0
  ! todo: move initialization messages into their respective routines
  if (hybrid%masterthread) write(iulog,*) 'initializing Polvani-Scott-Thomas baroclinic instability test'
  if (hybrid%masterthread) write(iulog,*) 'initializing ASP gravity wave test'
  if (hybrid%masterthread) write(iulog,*) 'initializing ASP mountain Rossby test'
  if (hybrid%masterthread) write(iulog,*) 'initializing ASP Rossby Haurwitz test'
  if (hybrid%masterthread) write(iulog,*) 'initializing pure tracer advection tests'
  if (hybrid%masterthread) write(iulog,*) 'initializing Jablonowski and Williamson ASP baroclinic instability test'
  if (hybrid%masterthread) write(iulog,*) 'initializing Jablonowski and Williamson baroclinic instability test V1'
  if (hybrid%masterthread) write(iulog,*) 'initializing Held-Suarez primitive equations test'
  if (hybrid%masterthread) write(iulog,*) 'initializing DCMIP2 test 2-0'
#endif

!_____________________________________________________________________
subroutine apply_test_forcing(elem,fvm,hybrid,hvcoord,n,n_tracer,dt,nets,nete)

  use dcmip_tests, only:  dcmip2012_test2_x_forcing

  ! compute and apply forcing terms needed by HOMME stand-alone tests

  implicit none
  type (element_t),       intent(inout) :: elem(:)                      ! element array
  type(fvm_struct),       intent(inout) :: fvm(:)                       ! finite volume structure
  type(hybrid_t),         intent(in)    :: hybrid                       ! hybrid parallel structure
  type (hvcoord_t),       intent(in)    :: hvcoord
  real (kind=real_kind),  intent(in)    :: dt
  integer,                intent(in)    :: n,nets,nete,n_tracer

  integer :: ie

  do ie=nets,nete
    elem(ie)%derived%FT = 0
    elem(ie)%derived%FM = 0
  enddo

  ! get forcing terms from test case
  select case(test_case)
    case('dcmip2012_test2_1');  call dcmip2012_test2_x_forcing(elem, hybrid,hvcoord,nets,nete,n,dt)
    case('dcmip2012_test2_2');  call dcmip2012_test2_x_forcing(elem, hybrid,hvcoord,nets,nete,n,dt)
  endselect

  ! apply forcing to state variables
  do ie=nets,nete

    ! apply forcing to dynamic variables
    elem(ie)%state%T(:,:,:,  n) = elem(ie)%state%T(:,:,:,  n) + dt * elem(ie)%derived%FT(:,:,:,  1)
    elem(ie)%state%v(:,:,:,:,n) = elem(ie)%state%v(:,:,:,:,n) + dt * elem(ie)%derived%FM(:,:,:,:,1)

    ! apply forcing to tracers (todo)
  enddo

end subroutine



end module
