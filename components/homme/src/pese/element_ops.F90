#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

!  get, set functions provided by each model

module element_ops

  use control_mod,    only: statefreq
  use derivative_mod, only: derivative_t
  use dimensions_mod, only: np, nlev, nlevp, nelemd, qsize, qsize_d
  use edge_mod,       only: edgevpack, edgevunpack
  use edgetype_mod,   only: edgebuffer_t
  use element_mod,    only: element_t
  use element_state,  only: timelevels, elem_state_t, derived_state_t
  use hybrid_mod,     only: hybrid_t
  use hybvcoord_mod,  only: hvcoord_t
  use kinds,          only: rl => real_kind, real_kind, iulog
  use perf_mod,       only: t_startf, t_stopf, t_barrierf, t_adj_detailf ! _EXTERNAL
  use parallel_mod,   only: abortmp
  use physical_constants, only: kappa, p0, g
  use shr_const_mod,  only: unset => shr_const_spval
  use vertical_se,    only: eta_derivative, elem_height, vertical_dss, v_interpolate

  implicit none

contains

  !_____________________________________________________________________
  subroutine get_field(elem,name,field,hvcoord,nt,ntQ)

    ! get scalar field by name

    implicit none
    type (element_t),       intent(in) :: elem
    character(len=*),       intent(in) :: name
    real (kind=real_kind),  intent(out):: field(np,np,nlev)
    type (hvcoord_t),       intent(in) :: hvcoord
    integer,                intent(in) :: nt
    integer,                intent(in) :: ntQ

    integer :: i,j,k, qi
    real(rl) :: var(np,np,nlev)

    ! get prognostic variables

    select case(name)

      case('T'  ); var = elem%state%T   (:,:,:,nt)
      case('u'  ); var = elem%state%v   (:,:,1,:,nt)
      case('v'  ); var = elem%state%v   (:,:,2,:,nt)
      case('p'  ); var = elem%state%dp3d(:,:,:,nt)
      case('ps' ); forall(k=1:nlev) var(:,:,k) = elem%state%ps_v(:,:,nt)

      case('Q'  ); var = elem%state%Q(:,:,:,1)
      case('Q2' ); qi=2; if(qsize_d>1) var = elem%state%Q(:,:,:,qi)
      case('Q3' ); qi=3; if(qsize_d>2) var = elem%state%Q(:,:,:,qi)
      case('Q4' ); qi=4; if(qsize_d>3) var = elem%state%Q(:,:,:,qi)
      case('Q5' ); qi=5; if(qsize_d>4) var = elem%state%Q(:,:,:,qi)

      case('phi');  var = elem%derived%phi
      case('geo');  var = elem%derived%phi
      case('omega');var = elem%derived%omega

      case default
        print *,'name = ',trim(name)
        call abortmp('ERROR: get_field name not supported in this model')
    end select

    ! interpolate each column from vertical gll nodes to vertical ouput levels
    do i=1,np
      do j=1,np
        field(i,j,:) = v_interpolate(var(i,j,:),nlev)
      enddo
    enddo
  end subroutine

  !_____________________________________________________________________________
  function get_vector_field(elem, name, n0, ni) result(vfield)

    ! Get vector-field data by name

    character*(*),    intent(in)  :: name
    type(element_t),  intent(in)  :: elem
    integer,          intent(in)  :: n0       ! time level index
    integer,          intent(in)  :: ni       ! num vertical interp pts

    real(rl) :: vfield(np,np,2,ni)
    real(rl) :: var(np,np,2,nlev)
    integer  :: i,j
    select case(name)                                                   ! switch on short variable name

      case('v');
        var = elem%state%v(:,:,:,:,n0)

      case default; var = unset                                         ! assign special "missing" value
    endselect

    ! interpolate each column
    do i=1,np
      do j=1,np
        vfield(i,j,1,:) = v_interpolate(var(i,j,1,:),ni)
        vfield(i,j,2,:) = v_interpolate(var(i,j,2,:),ni)
      enddo
    enddo

  end function

  !_____________________________________________________________________
  subroutine set_state(u,v,w,T,ps,phis,p,dp,zm,g,i,j,k,elem,n0,n1)

    ! set state variables at node(i,j,k) at layer midpoints

    real(real_kind),  intent(in)    :: u,v,w,T,ps,phis,p,dp,zm,g
    integer,          intent(in)    :: i,j,k,n0,n1
    type(element_t),  intent(inout) :: elem

    ! set prognostic state variables at level midpoints
    elem%state%v   (i,j,1,k,n0:n1) = u
    elem%state%v   (i,j,2,k,n0:n1) = v
    elem%state%T   (i,j,k,n0:n1)   = T
    elem%state%dp3d(i,j,k,n0:n1)   = dp
    elem%state%ps_v(i,j,n0:n1)     = ps
    elem%state%phis(i,j)           = phis

    ! set some diagnostic variables
    elem%derived%dp(i,j,k)         = dp
    elem%derived%phi(i,j,k)        = g*zm

  end subroutine

  !_____________________________________________________________________
  subroutine set_thermostate(elem,temperature,hvcoord,n0,n0_q)

    implicit none
    
    type (element_t),       intent(inout) :: elem
    real (kind=real_kind),  intent(in)    :: temperature(np,np,nlev)
    type (hvcoord_t),       intent(in)    :: hvcoord                    ! hybrid vertical coordinate struct
    integer :: n0,n0_q

    elem%state%T(:,:,:,n0)=temperature(:,:,:)

  end subroutine set_thermostate

  !_____________________________________________________________________
  subroutine copy_state(elem,nin,nout)
    implicit none
    type (element_t), intent(inout) :: elem
    integer,          intent(in)    :: nin,nout

    elem%state%v   (:,:,:,:,nout)= elem%state%v   (:,:,:,:,nin)
    elem%state%T   (:,:,:,  nout)= elem%state%T   (:,:,:,  nin)
    elem%state%dp3d(:,:,:,  nout)= elem%state%dp3d(:,:,:,  nin)
    elem%state%ps_v(:,:,    nout)= elem%state%ps_v(:,:,    nin)

  end subroutine copy_state

  !_____________________________________________________________________
  subroutine apply_vertical_dss(elem, nt)

    ! apply vertical dss to all prognostics and tracers

    type (element_t), intent(inout), target :: elem
    integer, intent(in) :: nt
    integer::l

    call vertical_dss( elem%state%T (:,:,  :,nt) )
    call vertical_dss( elem%state%v (:,:,1,:,nt) )
    call vertical_dss( elem%state%v (:,:,2,:,nt) )

    do l=1,qsize_d
      call vertical_dss(elem%state%Qdp(:,:,:,l,nt))
    enddo

  end subroutine

  !_____________________________________________________________________
  subroutine pack_edge_data(buffer, elem, nt, ie)

    ! pack prognostics into edge buffer

    type (EdgeBuffer_t),  intent(inout)       :: buffer                 ! buffer for edge-data exchange
    type (element_t),     intent(in), target  :: elem                   ! array of element_t structures
    integer,              intent(in)          :: nt, ie                 ! time level and element index

    integer :: i
    i = 0
    call edgeVpack(buffer, elem%state%ps_v(:,:,    nt), 1     , i, ie); i=i+1
    call edgeVpack(buffer, elem%state%T   (:,:,:,  nt), nlev  , i, ie); i=i+nlev
    call edgeVpack(buffer, elem%state%v   (:,:,:,:,nt), 2*nlev, i, ie); i=i+2*nlev

  end subroutine

  !_____________________________________________________________________
  subroutine unpack_edge_data(buffer, elem, nt, ie)

    ! unpack prognostics from edge buffer

    type (EdgeBuffer_t),  intent(inout)       :: buffer                 ! buffer for edge-data exchange
    type (element_t),     intent(inout), target  :: elem                ! array of element_t structures
    integer,              intent(in)          :: nt, ie                 ! time level and element index

    integer :: i
    i = 0
    call edgeVunpack(buffer, elem%state%ps_v(:,:,    nt), 1     , i, ie); i=i+1
    call edgeVunpack(buffer, elem%state%T   (:,:,:,  nt), nlev  , i, ie); i=i+nlev
    call edgeVunpack(buffer, elem%state%v   (:,:,:,:,nt), 2*nlev, i, ie); i=i+2*nlev

  end subroutine

  !_____________________________________________________________________
  subroutine apply_map(M, elem, nt)

    ! apply pointwise 2d map to each prognostic

    real(rl),         intent(in)            :: M(np,np)
    type (element_t), intent(inout), target :: elem                     ! array of element_t structures
    integer,          intent(in)            :: nt                       ! time level

    integer :: k
    elem%state%ps_v(:,:, nt) = elem%state%ps_v(:,:,nt)*M
    do k=1,nlev
      elem%state%T(:,:,k,  nt) = elem%state%T(:,:,  k,nt)*M
      elem%state%v(:,:,1,k,nt) = elem%state%v(:,:,1,k,nt)*M
      elem%state%v(:,:,2,k,nt) = elem%state%v(:,:,2,k,nt)*M
    enddo

  end subroutine

  !_____________________________________________________________________
  subroutine display_max_and_min(elem, hybrid, ie, nt, count)

    ! display verbose diagnostics

    type (element_t),   intent(inout), target :: elem                 ! element
    type (hybrid_t),    intent(in)		:: hybrid												! mpi/omp data struct
    integer,            intent(in)    :: ie                           ! element number
    integer,            intent(in)    :: nt                           ! time level to display
    integer,            intent(in)    :: count                        ! iteration counter

    real(rl):: maxu,minu,maxT,minT,maxv,minv,maxps,minps

    if( mod(count, statefreq)==0 ) then

      maxu = maxval(elem%state%v   (:,:,1,:,nt))
      minu = minval(elem%state%v   (:,:,1,:,nt))
      maxv = maxval(elem%state%v   (:,:,2,:,nt))
      minv = minval(elem%state%v   (:,:,2,:,nt))
      maxps= maxval(elem%state%ps_v(:,:,    nt))
      minps= minval(elem%state%ps_v(:,:,    nt))
      maxT = maxval(elem%state%T   (:,:,:,  nt))
      minT = minval(elem%state%T   (:,:,:,  nt))

      if (hybrid%masterthread) then
        print *,"nt =",nt, "ie =",ie, " max u =",maxu," max v =",maxv," max T =",maxT," max ps=",maxps
      endif

      if (any(isnan(elem%state%v(:,:,1,:,nt))))  stop 'detected NaN in u'
      if (any(isnan(elem%state%v(:,:,2,:,nt))))  stop 'detected NaN in v'
      if (any(isnan(elem%state%T(:,:,:,nt))))    stop 'detected NaN in T'
      if (any(isnan(elem%state%ps_v(:,:,nt))))   stop 'detected NaN in ps_v'

    endif

  end subroutine

end module

