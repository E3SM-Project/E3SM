
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module viscosity_mod
  use viscosity_mod_base, only: compute_zeta_C0, compute_div_C0, compute_zeta_C0_contra, compute_div_C0_contra, make_c0, make_c0_vector
  use viscosity_mod_base, only: biharmonic_wk_scalar,neighbor_minmax, neighbor_minmax_start,neighbor_minmax_finish
  use viscosity_preqx_base, only: biharmonic_wk_dp3d

  use thread_mod, only : omp_get_num_threads
  use kinds, only : real_kind, iulog
  use dimensions_mod, only : np, nlev,qsize,nelemd
  use hybrid_mod, only : hybrid_t, hybrid_create
  use parallel_mod, only : parallel_t
  use element_mod, only : element_t
  use edgetype_mod, only : EdgeBuffer_t, EdgeDescriptor_t
  use bndry_mod, only : bndry_exchangev, bndry_exchangeS, bndry_exchangeS_start,bndry_exchangeS_finish
  use control_mod, only : hypervis_scaling, nu, nu_div
  use perf_mod, only: t_startf, t_stopf
  implicit none
  private

  public :: compute_zeta_C0, compute_div_C0, compute_zeta_C0_contra, compute_div_C0_contra, make_c0, make_c0_vector
  public :: biharmonic_wk_scalar, neighbor_minmax, neighbor_minmax_start,neighbor_minmax_finish, biharmonic_wk_dp3d
  public :: biharmonic_wk_scalar_openacc
  public :: neighbor_minmax_openacc



contains

  subroutine biharmonic_wk_scalar_openacc(elem,qtens,grads,deriv,edgeq,hybrid,nets,nete)
    use hybrid_mod            , only: hybrid_t
    use element_mod           , only: element_t
    use edgetype_mod          , only: edgeBuffer_t
    use derivative_mod        , only: derivative_t
    use control_mod           , only: hypervis_scaling
    use perf_mod              , only: t_startf, t_stopf
    use derivative_mod, only: laplace_sphere_wk_openacc
    use edge_mod      , only: edgeVpack_openacc, edgeVunpack_openacc
    use bndry_mod     , only: bndry_exchangeV => bndry_exchangeV_simple_overlap
    implicit none
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! compute weak biharmonic operator
    !    input:  qtens = Q
    !    output: qtens = weak biharmonic of Q
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    type (element_t)     , intent(inout) :: elem(:)
    real (kind=real_kind), intent(inout) :: qtens(np,np,nlev,qsize,nelemd)
    real(kind=real_kind) , intent(inout) :: grads(np,np,2,nlev,qsize,nelemd)
    type (derivative_t)  , intent(in   ) :: deriv
    type (EdgeBuffer_t)  , intent(inout) :: edgeq
    type (hybrid_t)      , intent(in   ) :: hybrid
    integer              , intent(in   ) :: nets,nete
    ! local
    integer :: k,kptr,i,j,ie,ic,q
    logical :: var_coef1
    !if tensor hyperviscosity with tensor V is used, then biharmonic operator is (\grad\cdot V\grad) (\grad \cdot \grad) 
    !so tensor is only used on second call to laplace_sphere_wk
    var_coef1 = .true.
    if(hypervis_scaling > 0) var_coef1 = .false.
    !$omp barrier
    !$omp master
    call laplace_sphere_wk_openacc(qtens,grads,deriv,elem,var_coef1,qtens,nlev*qsize,nets,nete,1,1)
    call t_startf('biwksc_PEU')
    call edgeVpack_openacc(edgeq,qtens,qsize*nlev,0,elem(:),nets,nete,1,1)
    !$omp end master
    !$omp barrier

    call t_startf('biwksc_exch')
    call bndry_exchangeV(hybrid,edgeq)
    call t_stopf('biwksc_exch')
    
    !$omp barrier
    !$omp master
    call edgeVunpack_openacc(edgeq,qtens,qsize*nlev,0,elem(:),nets,nete,1,1)
    call t_stopf('biwksc_PEU')
    !$acc parallel loop gang vector collapse(5) present(qtens,elem(:))
    do ie = nets , nete
      ! apply inverse mass matrix, then apply laplace again
      do q = 1 , qsize      
        do k = 1 , nlev    !  Potential loop inversion (AAM)
          do j = 1 , np
            do i = 1 , np
              qtens(i,j,k,q,ie) = elem(ie)%rspheremp(i,j)*qtens(i,j,k,q,ie)
            enddo
          enddo
        enddo
      enddo
    enddo
    call laplace_sphere_wk_openacc(qtens,grads,deriv,elem,.true.,qtens,nlev*qsize,nets,nete,1,1)
    !$omp end master
    !$omp barrier
  end subroutine biharmonic_wk_scalar_openacc

  subroutine neighbor_minmax_openacc(elem,hybrid,edgeMinMax,nets,nete,min_neigh,max_neigh)
    use hybrid_mod       , only: hybrid_t
    use element_mod      , only: element_t
    use perf_mod         , only: t_startf, t_stopf
    use edgetype_mod     , only: edgeBuffer_t
    use edge_mod , only: edgeSpack_openacc, edgeSunpackMin_openacc, edgeSunpackMax_openacc
    use bndry_mod, only: bndry_exchangeS => bndry_exchangeS_simple_overlap
    implicit none
    ! compute Q min&max over the element and all its neighbors
    integer :: nets,nete
    type (element_t)     , intent(in   ) :: elem(:)
    type (hybrid_t)      , intent(in   ) :: hybrid
    type (EdgeBuffer_t)  , intent(inout) :: edgeMinMax
    real (kind=real_kind), intent(inout) :: min_neigh(nlev,qsize,nelemd)
    real (kind=real_kind), intent(inout) :: max_neigh(nlev,qsize,nelemd)
    ! local
    integer :: ie,k,q,j,i
    ! compute Qmin, Qmax
    !$omp barrier
    !$omp master
    call t_startf('nmm_PEU')
    call edgeSpack_openacc(edgeMinMax,min_neigh,nlev*qsize,0         ,elem(:),nets,nete,1,1)
    call edgeSpack_openacc(edgeMinMax,max_neigh,nlev*qsize,nlev*qsize,elem(:),nets,nete,1,1)
    !$omp end master
    !$omp barrier

    call t_startf('nmm_exch')
    call bndry_exchangeS(hybrid,edgeMinMax)
    call t_stopf('nmm_exch')
       
    !$omp barrier
    !$omp master
    call edgeSunpackMin_openacc(edgeMinMax,min_neigh,nlev*qsize,0         ,elem(:),nets,nete,1,1)
    call edgeSunpackMax_openacc(edgeMinMax,max_neigh,nlev*qsize,nlev*qsize,elem(:),nets,nete,1,1)
    call t_stopf('nmm_PEU')
    !$omp end master
    !$omp barrier
  end subroutine neighbor_minmax_openacc

end module viscosity_mod

