
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module viscosity_openacc_mod
  use kinds, only: real_kind
  use dimensions_mod, only: np,nlev,qsize,nelemd
  implicit none
  private
  real(kind=real_kind), allocatable :: Qmin_pack(:,:,:,:,:)
  real(kind=real_kind), allocatable :: Qmax_pack(:,:,:,:,:)
  logical :: tracer_pack_allocated = .false.

  public :: neighbor_minmax
  public :: biharmonic_wk_scalar
  public :: biharmonic_wk_scalar_minmax

contains

  subroutine allocate_tracer_pack_arrays()
    implicit none
    allocate(qmin_pack(np,np,nlev,qsize,nelemd))
    allocate(qmax_pack(np,np,nlev,qsize,nelemd))
    !$acc enter data pcreate(qmin_pack,qmax_pack)
    tracer_pack_allocated = .true.
  end subroutine allocate_tracer_pack_arrays

  subroutine biharmonic_wk_scalar(elem,qtens,grads,deriv,edgeq,hybrid,nets,nete)
    use hybrid_mod            , only: hybrid_t
    use element_mod           , only: element_t
    use edgetype_mod          , only: edgeBuffer_t
    use derivative_mod        , only: derivative_t
    use control_mod           , only: hypervis_scaling
    use perf_mod              , only: t_startf, t_stopf
    use derivative_openacc_mod, only: laplace_sphere_wk
    use edge_openacc_mod      , only: edgeVpack, edgeVunpack
    use bndry_openacc_mod     , only: bndry_exchangeV => bndry_exchangeV_finer_overlap
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
    if (.not. tracer_pack_allocated) call allocate_tracer_pack_arrays()
    call laplace_sphere_wk(qtens,grads,deriv,elem,var_coef1,qtens,nlev*qsize,nets,nete,1,1)
    call t_startf('biwksc_PEU')
    call edgeVpack(edgeq,qtens,qsize*nlev,0,elem(:),nets,nete,1,1)
    !$omp end master
    !$omp barrier

    call t_startf('biwksc_exch')
    call bndry_exchangeV(hybrid,edgeq)
    call t_stopf('biwksc_exch')
    
    !$omp barrier
    !$omp master
    call edgeVunpack(edgeq%buf,edgeq%nlyr,qtens,qsize*nlev,0,elem(:),nets,nete,1,1)
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
    call laplace_sphere_wk(qtens,grads,deriv,elem,.true.,qtens,nlev*qsize,nets,nete,1,1)
    !$omp end master
    !$omp barrier
  end subroutine biharmonic_wk_scalar

  subroutine biharmonic_wk_scalar_minmax(elem,qtens,grads,deriv,edgeq,hybrid,nets,nete,emin,emax)
    use hybrid_mod            , only: hybrid_t
    use element_mod           , only: element_t
    use derivative_mod        , only: derivative_t
    use control_mod           , only: hypervis_scaling, hypervis_power
    use perf_mod              , only: t_startf, t_stopf
    use derivative_openacc_mod, only: laplace_sphere_wk
    use edgetype_mod          , only: edgeBuffer_t
    use edge_openacc_mod      , only: edgeVpack, edgeVunpack, edgeVunpackMin, edgeVunpackMax
    use bndry_openacc_mod     , only: bndry_exchangeV => bndry_exchangeV_finer_overlap
    implicit none
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! compute weak biharmonic operator
    !    input:  qtens = Q
    !    output: qtens = weak biharmonic of Q and Q element min/max
    !
    !    note: emin/emax must be initialized with Q element min/max.  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    type (hybrid_t)      , intent(in   ) :: hybrid
    type (element_t)     , intent(in   ) :: elem(:)
    integer              , intent(in   ) :: nets,nete
    real (kind=real_kind), intent(inout) :: qtens(np,np,nlev,qsize,nelemd)
    real(kind=real_kind) , intent(inout) :: grads(np,np,2,nlev,qsize,nelemd)
    type (EdgeBuffer_t)  , intent(inout) :: edgeq
    type (derivative_t)  , intent(in   ) :: deriv
    real (kind=real_kind), intent(inout) :: emin(nlev,qsize,nelemd)
    real (kind=real_kind), intent(inout) :: emax(nlev,qsize,nelemd)
    ! local
    integer :: k,kptr,i,j,ie,ic,q
    real (kind=real_kind) :: lap_p(np,np)
    logical :: var_coef1
    !if tensor hyperviscosity with tensor V is used, then biharmonic operator is (\grad\cdot V\grad) (\grad \cdot \grad) 
    !so tensor is only used on second call to laplace_sphere_wk
    var_coef1 = .true.
    if(hypervis_scaling > 0)    var_coef1 = .false.
    !$omp barrier
    !$omp master
    if (.not. tracer_pack_allocated) call allocate_tracer_pack_arrays()
    call minmax_pack(qmin_pack,qmax_pack,emin,emax,nets,nete)
    call laplace_sphere_wk(qtens,grads,deriv,elem,var_coef1,qtens,nlev*qsize,nets,nete,1,1)
    call t_startf('biwkscmm_PEU')
    call edgeVpack(edgeq,    qtens,qsize*nlev,0           ,elem(:),nets,nete,1,1)
    call edgeVpack(edgeq,Qmin_pack,nlev*qsize,nlev*qsize  ,elem(:),nets,nete,1,1)
    call edgeVpack(edgeq,Qmax_pack,nlev*qsize,2*nlev*qsize,elem(:),nets,nete,1,1)
    !$omp end master
    !$omp barrier

    call t_startf('biwkscmm_exch')
    call bndry_exchangeV(hybrid,edgeq)
    call t_stopf('biwkscmm_exch')

    !$omp barrier
    !$omp master
    call edgeVunpack   (edgeq%buf,edgeq%nlyr,    qtens,qsize*nlev,0           ,elem(:),nets,nete,1,1)
    call edgeVunpackMin(edgeq%buf,edgeq%nlyr,Qmin_pack,qsize*nlev,qsize*nlev  ,elem(:),nets,nete,1,1)
    call edgeVunpackMax(edgeq%buf,edgeq%nlyr,Qmax_pack,qsize*nlev,2*qsize*nlev,elem(:),nets,nete,1,1)
    call t_stopf('biwkscmm_PEU')
    !$acc parallel loop gang vector collapse(5) present(qtens,elem(:))
    do ie = nets , nete
      do q = 1 , qsize      
        do k = 1 , nlev
          do j = 1 , np
            do i = 1 , np
              qtens(i,j,k,q,ie) = elem(ie)%rspheremp(i,j)*qtens(i,j,k,q,ie)  ! apply inverse mass matrix
            enddo
          enddo
        enddo
      enddo
    enddo
    call laplace_sphere_wk(qtens,grads,deriv,elem,.true.,qtens,nlev*qsize,nets,nete,1,1)
    call minmax_reduce_corners(qmin_pack,qmax_pack,emin,emax,nets,nete)
    !$omp end master
    !$omp barrier
  end subroutine biharmonic_wk_scalar_minmax

  subroutine neighbor_minmax(elem,hybrid,edgeMinMax,nets,nete,min_neigh,max_neigh)
    use hybrid_mod       , only: hybrid_t
    use element_mod      , only: element_t
    use perf_mod         , only: t_startf, t_stopf
    use edgetype_mod     , only: edgeBuffer_t
    use edge_openacc_mod , only: edgeVpack, edgeVunpackMin, edgeVunpackMax
    use bndry_openacc_mod, only: bndry_exchangeV => bndry_exchangeV_finer_overlap
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
    if (.not. tracer_pack_allocated) call allocate_tracer_pack_arrays()
    call minmax_pack(qmin_pack,qmax_pack,min_neigh,max_neigh,nets,nete)
    call t_startf('nmm_PEU')
    call edgeVpack(edgeMinMax,Qmin_pack,nlev*qsize,0         ,elem(:),nets,nete,1,1)
    call edgeVpack(edgeMinMax,Qmax_pack,nlev*qsize,nlev*qsize,elem(:),nets,nete,1,1)
    !$omp end master
    !$omp barrier

    call t_startf('nmm_exch')
    call bndry_exchangeV(hybrid,edgeMinMax)
    call t_stopf('nmm_exch')
       
    !$omp barrier
    !$omp master
    call edgeVunpackMin(edgeMinMax%buf,edgeMinMax%nlyr,Qmin_pack,nlev*qsize,0         ,elem(:),nets,nete,1,1)
    call edgeVunpackMax(edgeMinMax%buf,edgeMinMax%nlyr,Qmax_pack,nlev*qsize,nlev*qsize,elem(:),nets,nete,1,1)
    call t_stopf('nmm_PEU')
    call minmax_reduce_corners(qmin_pack,qmax_pack,min_neigh,max_neigh,nets,nete)
    !$omp end master
    !$omp barrier
  end subroutine neighbor_minmax

  subroutine minmax_pack(qmin_pack,qmax_pack,min_neigh,max_neigh,nets,nete)
    implicit none
    real(kind=real_kind), intent(  out) :: qmin_pack(np,np,nlev,qsize,nelemd)
    real(kind=real_kind), intent(  out) :: qmax_pack(np,np,nlev,qsize,nelemd)
    real(kind=real_kind), intent(in   ) :: min_neigh(nlev,qsize,nelemd)
    real(kind=real_kind), intent(in   ) :: max_neigh(nlev,qsize,nelemd)
    integer             , intent(in   ) :: nets , nete
    integer :: ie,q,k,j,i
    !$acc parallel loop gang vector collapse(5) present(min_neigh,max_neigh,qmin_pack,qmax_pack)
    do ie = nets , nete
      do q = 1 , qsize
        do k = 1 , nlev
          do j = 1 , np
            do i = 1 , np
              Qmin_pack(i,j,k,q,ie) = min_neigh(k,q,ie)
              Qmax_pack(i,j,k,q,ie) = max_neigh(k,q,ie)
            enddo
          enddo
        enddo
      enddo
    enddo
  end subroutine minmax_pack

  subroutine minmax_reduce_corners(qmin_pack,qmax_pack,min_neigh,max_neigh,nets,nete)
    implicit none
    real(kind=real_kind), intent(in   ) :: qmin_pack(np,np,nlev,qsize,nelemd)
    real(kind=real_kind), intent(in   ) :: qmax_pack(np,np,nlev,qsize,nelemd)
    real(kind=real_kind), intent(  out) :: min_neigh(nlev,qsize,nelemd)
    real(kind=real_kind), intent(  out) :: max_neigh(nlev,qsize,nelemd)
    integer             , intent(in   ) :: nets , nete
    integer :: ie,q,k
    !$acc parallel loop gang vector collapse(3) present(min_neigh,max_neigh,qmin_pack,qmax_pack)
    do ie = nets , nete
      do q = 1 , qsize
        do k = 1 , nlev
          ! note: only need to consider the corners, since the data we packed was
          ! constant within each element
          min_neigh(k,q,ie)=min(qmin_pack(1,1,k,q,ie),qmin_pack(1,np,k,q,ie),qmin_pack(np,1,k,q,ie),qmin_pack(np,np,k,q,ie))
          min_neigh(k,q,ie)=max(min_neigh(k,q,ie),0d0)
          max_neigh(k,q,ie)=max(qmax_pack(1,1,k,q,ie),qmax_pack(1,np,k,q,ie),qmax_pack(np,1,k,q,ie),qmax_pack(np,np,k,q,ie))
        enddo
      enddo
    enddo
  end subroutine minmax_reduce_corners

end module viscosity_openacc_mod

