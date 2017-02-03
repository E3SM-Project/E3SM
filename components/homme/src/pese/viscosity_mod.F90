#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module viscosity_mod

use viscosity_base, only: compute_zeta_C0, compute_div_C0, compute_zeta_C0_contra, compute_div_C0_contra, make_c0, make_c0_vector, neighbor_minmax, biharmonic_wk_scalar, neighbor_minmax_start,neighbor_minmax_finish, smooth_phis

	use kinds,              only: rl => real_kind, dd => longdouble_kind
	use dimensions_mod,     only: np, nlev, nlevp, nelem, nelemd
	use hybrid_mod,         only: hybrid_t
	use element_mod,        only: element_t
  use element_state,      only: elem_state_t, derived_state_t
	use derivative_mod,     only: derivative_t, laplace_sphere_wk, vlaplace_sphere_wk
	use edge_mod,           only: initEdgeBuffer, edgevpack, edgevunpack
  use edgetype_mod,       only: EdgeBuffer_t
	use bndry_mod,          only: bndry_exchangev
	use hybvcoord_mod,      only: hvcoord_t
  use control_mod,        only: nu, nu_s, nu_p, nu_div, nu_top, hypervis_order, hypervis_subcycle, psurf_vis
  use element_ops,        only: pack_edge_data, unpack_edge_data, apply_vertical_dss, apply_map
  use physical_constants, only: Cp, cpwater_vapor, Rgas, kappa, p0, g
  use vertical_se,        only: vertical_dss

	implicit none

  CONTAINS

    !_____________________________________________________________________
    subroutine apply_laplacian(T,v,elem,edge_buffer,hybrid,deriv,nets,nete,apply_map)

      real(rl),             intent(inout) :: T(np,np,nlev,nelemd)
      real(rl),             intent(inout) :: v(np,np,2,nlev,nelemd)
      type (element_t),			intent(inout), target :: elem(:)							! array of element_t structures
      type (EdgeBuffer_t),  intent(inout) :: edge_buffer
      type (hybrid_t),			intent(in)		:: hybrid												! mpi/omp data struct
      type (derivative_t),	intent(in)		:: deriv												! horizontal derivative data struct
      integer,							intent(in)		:: nets,nete										! start and end element indices
      logical,              intent(in)    :: apply_map

      type (elem_state_t), pointer :: s                                   ! ptr to element state variables
      integer :: ie,i,k                                              ! loop indicies

      do ie=nets,nete
        s => elem(ie)%state

          ! apply laplace operator
          do k=1,nlev
              T(:,:,k,ie)   =  laplace_sphere_wk(T(:,:,k,ie)   ,deriv,elem(ie),var_coef=.false.)
              v(:,:,:,k,ie) = vlaplace_sphere_wk(v(:,:,:,k,ie),deriv,elem(ie),var_coef=.false.)
          enddo!k

          ! apply mass matrix
          if(apply_map) then
            do k=1,nlev
              T(:,:,k,ie)   = T(:,:,k,ie)  *elem(ie)%spheremp
              v(:,:,1,k,ie) = v(:,:,1,k,ie)*elem(ie)%spheremp
              v(:,:,2,k,ie) = v(:,:,2,k,ie)*elem(ie)%spheremp
            enddo!k
          endif

          ! pack data into edge buffer
          i = 0
          call edgeVpack(edge_buffer, T(:,:,:,ie)  , nlev  , i, ie); i=i+nlev
          call edgeVpack(edge_buffer, v(:,:,:,:,ie), 2*nlev, i, ie); i=i+2*nlev

        enddo!ie

        ! exchange edge data
        call bndry_exchangeV(hybrid, edge_buffer)                       ! exchange horizontal edge data and perform dss

        do ie=nets,nete
          s => elem(ie)%state

          ! unpack data
          i = 0
          call edgeVunpack(edge_buffer, T(:,:,:,ie)  , nlev  , i, ie); i=i+nlev
          call edgeVunpack(edge_buffer, v(:,:,:,:,ie), 2*nlev, i, ie); i=i+2*nlev

          ! apply inverse mass matrix
          do k=1,nlev
            T(:,:,k,ie)   = T(:,:,k,ie)  *elem(ie)%rspheremp
            v(:,:,1,k,ie) = v(:,:,1,k,ie)*elem(ie)%rspheremp
            v(:,:,2,k,ie) = v(:,:,2,k,ie)*elem(ie)%rspheremp
          enddo!k

          ! perform stiffness summation
          call vertical_dss( T(:,:,:,ie) )
          call vertical_dss( v(:,:,1,:,ie) )
          call vertical_dss( v(:,:,2,:,ie) )

        enddo!ie
  end subroutine

  !_____________________________________________________________________
  subroutine apply_hyperviscosity(elem,edge_buffer,hvcoord,hybrid,deriv,nt,nets,nete,dt2,eta_ave_w)

    ! apply artificial horizontal diffusion to combat spurrious energy build-up

    type (element_t),			intent(inout), target :: elem(:)							! array of element_t structures
    type (EdgeBuffer_t),  intent(inout) :: edge_buffer
    type (hvcoord_t),			intent(inout)	:: hvcoord											! hybrid vertical coord data struct
    type (hybrid_t),			intent(in)		:: hybrid												! mpi/omp data struct
    type (derivative_t),	intent(in)		:: deriv												! horizontal derivative data struct
    integer,							intent(in)		:: nt                           ! time index
    integer,							intent(in)		:: nets,nete										! start and end element indices
    real*8,               intent(in)    :: dt2
    real (rl),            intent(in)    :: eta_ave_w

    real(rl), dimension(np,np,nlev,nelemd)   :: T_diff
    real(rl), dimension(np,np,2,nlev,nelemd) :: v_diff

    type (elem_state_t), pointer :: s                                   ! ptr to element state variables
    real(rl):: dt
    integer :: ic,ie,i,k                                              ! loop indicies

    ! Exit if hypervis coefficients are all zero
    if (nu_s == 0 .and. nu == 0 .and. nu_p==0 ) return;
    ! get subcycled time-step
    dt=dt2/hypervis_subcycle

    ! subcycle hyperviscosity if needed
    do ic=1,hypervis_subcycle                                         ! subcycle if needed

      ! copy data to diff-fields
      do ie=nets,nete
        s => elem(ie)%state
        T_diff(:,:,:,ie)   = s%T(:,:,:,nt)
        v_diff(:,:,:,:,ie) = s%v(:,:,:,:,nt)
      enddo

      ! apply 1st-order hyperviscosity
      if (hypervis_order == 1) then

        ! apply laplacian once
        call apply_laplacian(T_diff,v_diff,elem,edge_buffer,hybrid,deriv,nets,nete,.true.)

        do ie=nets,nete
          s => elem(ie)%state
          s%T(:,:,:,nt)   =s%T (:,:,:,nt)  + dt * nu_s * T_diff(:,:,:,ie)
          s%v(:,:,:,:,nt) =s%v (:,:,:,:,nt)+ dt * nu   * v_diff(:,:,:,:,ie)
          call apply_vertical_dss(elem(ie), nt)
        enddo!ie

      endif

      ! apply 2nd-order hyperviscosity
      if (hypervis_order == 2) then

        ! apply laplacian twice
        call apply_laplacian(T_diff,v_diff,elem,edge_buffer, hybrid,deriv,nets,nete,.false.)
        call apply_laplacian(T_diff,v_diff,elem,edge_buffer, hybrid,deriv,nets,nete,.false.)

        do ie=nets,nete
          s => elem(ie)%state
          s%T(:,:,:,nt)   =s%T (:,:,:,nt)  - dt * nu_s * T_diff(:,:,:,ie)
          s%v(:,:,:,:,nt) =s%v (:,:,:,:,nt)- dt * nu   * v_diff(:,:,:,:,ie)
          call apply_vertical_dss(elem(ie), nt)
        enddo!ie
      endif
    enddo ! ic

  end subroutine

end module 



