#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#undef _DGEMV
module solver_mod
  use kinds, only : real_kind, int_kind
  use dimensions_mod, only : npsq, nlev
  use perf_mod, only: t_startf, t_stopf ! _EXTERNAL
  use parallel_mod, only : abortmp
  implicit none
  private

  character(len=8), private, parameter :: blkjac_storage = "inverse"
  !  character(len=8), private, parameter :: blkjac_storage = "LUfactor"

  type, public :: blkjac_t
     real (kind=real_kind), dimension(npsq,npsq,nlev) :: E
     integer(kind=int_kind),     dimension(npsq,nlev) :: ipvt
  end type blkjac_t


  public  :: pcg_solver
  public  :: blkjac_init

  interface pcg_solver
     module procedure pcg_solver_stag
     module procedure pcg_solver_nonstag
  end interface

contains

  function pcg_solver_stag(elem,  & 
       rhs,      &
       cg,       &
       red,      &
       edge2,    &   
       lambdasq, &   
       deriv,    &   
       nets,     & 
       nete,     &
       blkjac) result(x) 
    use dimensions_mod, only : nlev, np, npsq
    use element_mod, only : element_t
    use reduction_mod, only : reductionbuffer_ordered_1d_t
    use cg_mod, only : cg_t, congrad
    use edge_mod, only : edgebuffer_t, edgevpack, edgevunpack!,edgerotate
    use derivative_mod, only : derivative_t, gradient_wk, gradient, divergence, derivative_stag_t
    use control_mod, only : maxits, while_iter, tol, precon_method
    use physical_constants, only : rrearth
    use bndry_mod, only : bndry_exchangeV
    use linear_algebra_mod, only : matvec
    use parallel_mod, only : syncmp, haltmp


    type(element_t), intent(in), target :: elem(:)
    integer, intent(in)  :: nets,nete
    real (kind=real_kind), intent(in) :: rhs(np,np,nlev,nets:nete) ! right hand side of operator
    type (cg_t)                       :: cg             ! conjugate gradient    (private)
    type (ReductionBuffer_ordered_1d_t)  :: red         ! CG reduction buffer   (shared memory)
    type (EdgeBuffer_t)               :: edge2          ! Laplacian edge buffer (shared memory)
    real (kind=real_kind)             :: lambdasq(nlev) ! Helmholtz lengthscale (private)
    type (derivative_t)          :: deriv          ! Staggered derivative struct     (private)
    type (blkjac_t)		      :: blkjac(nets:nete)
    real (kind=real_kind)             :: x(np,np,nlev,nets:nete)     ! solution (result)

#ifdef CAM
    call haltmp('semi-implicit method not yet supported in cam')
#else

    ! ===========
    ! Local
    ! ===========

    real (kind=real_kind),pointer :: metinv(:,:,:,:), Dinv(:,:,:,:)
    real (kind=real_kind),pointer :: metdet(:,:)
    real (kind=real_kind),pointer :: rmp(:,:)
    real (kind=real_kind),pointer :: mp(:,:)

    real (kind=real_kind) :: gradp(np,np,2,nlev,nets:nete)
    real (kind=real_kind) :: div(np,np)
    real (kind=real_kind) :: p(np,np)
    real (kind=real_kind) :: r(npsq)
    real (kind=real_kind) :: z(npsq)

    real (kind=real_kind) :: gradp1
    real (kind=real_kind) :: gradp2

    integer :: ie
    integer :: i,j,k
    integer :: kptr
    integer :: iptr
    integer :: ieptr

    character(len=1) :: Trans
    real(kind=real_kind) :: one,zero
    integer(kind=int_kind) :: inc

    call t_startf('pcg_solver_stag')

    ! ========================================
    ! Load the rhs in the cg struct...
    ! ========================================

#ifndef CAM

    !DBG print *,'pcg_solver: point #1'
    do ie=nets,nete
       ieptr=ie-nets+1
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,iptr,i,j)
#endif
       do k=1,nlev
          iptr=1
          do j=1,np
             do i=1,np
                cg%state(ieptr)%r(iptr,k) = rhs(i,j,k,ie)
                iptr=iptr+1
             end do
          end do
       end do
    end do
    !DBG print *,'pcg_solver: point #2'

    Trans='N'
    one=1.0D0
    zero=0.0D0
    inc=1

    ! cg%debug_level = 3
    do while (congrad(cg,red,maxits,tol))
       call t_startf('timer_helm')

!       print *,'CG inter = ',cg%iter
       do ie=nets,nete
          ieptr=ie-nets+1
          Dinv   => elem(ie)%Dinv
          metinv => elem(ie)%metinv
          metdet => elem(ie)%metdet

#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,Trans,inc,deriv,i,j,gradp1,gradp2)
#endif
          do k=1,nlev
             if (.not.cg%converged(k)) then

                ! ========================================
                ! Apply preconditioner: wrk2 = (M^-1) * wrk1 
                ! ========================================

                if (precon_method == "block_jacobi") then

                   if (blkjac_storage == "LUfactor") then
                      call abortmp('dgesl needs linpack')
!                      call dgesl(cg%state(ieptr)%r(:,k), &
!                                 cg%state(ieptr)%z(:,k), &
!                                 blkjac(ie)%E(:,:,k),    &
!                                 blkjac(ie)%ipvt(:,k)) 
                   else if (blkjac_storage == "inverse") then
#ifdef _DGEMV
                      call dgemv(Trans,npsq,npsq,one,blkjac(ie)%E(:,:,k),npsq,cg%state(ieptr)%r(:,k),inc,zero,cg%state(ieptr)%z(:,k),inc)
#else
                      call matvec(cg%state(ieptr)%r(:,k),cg%state(ieptr)%z(:,k),blkjac(ie)%E(:,:,k),npsq)
#endif
                   end if

                else if (precon_method == "identity") then

                   do j=1,npsq
                      cg%state(ieptr)%z(j,k)=cg%state(ieptr)%r(j,k)
                   end do

                end if


                !JMD===========================================
                !JMD   2*np*np*(np + np) Flops 
   		!JMD  SR = (4*np*np + 2*np*np + np*np)*Ld
   		!JMD  SUM(WS) = (6*np*np + 2*np*np + np*np
                !JMD===========================================
#ifdef _WK_GRAD
#ifdef _NEWSTRUCT
                gradp(:,:,:,k,ie)=gradient_wk(cg%state(ieptr,k)%z(:),deriv)*rrearth
#else
                !JPE temp solution get rid of reshape!
                gradp(:,:,:,k,ie)=gradient_wk(reshape(cg%state(ieptr)%z(:,k),(/np,np/)),deriv)*rrearth
#endif
#else
#ifdef _NEWSTRUCT
                gradp(:,:,:,k,ie)=gradient(cg%state(ieptr,k)%z(:),deriv)*rrearth
#else
                !JPE temp solution get rid of reshape!
                gradp(:,:,:,k,ie)=gradient(reshape(cg%state(ieptr)%z(:,k),(/np,np/)),deriv)*rrearth
#endif
#endif

                ! =======================================
                ! rotate gradient to form contravariant
                !JMD  4*np*np Flops
                ! =======================================

                do j=1,np
                   do i=1,np
                      gradp1       = gradp(i,j,1,k,ie)
                      gradp2       = gradp(i,j,2,k,ie)
                      gradp(i,j,1,k,ie) = metdet(i,j)*(Dinv(1,1,i,j)*gradp1 + &
                           Dinv(2,1,i,j)*gradp2)
                      gradp(i,j,2,k,ie) = metdet(i,j)*(Dinv(1,2,i,j)*gradp1 + &
                           Dinv(2,2,i,j)*gradp2)
                   end do
                end do

             end if
          end do

          kptr=0
          call edgeVpack(edge2,gradp(1,1,1,1,ie),2*nlev,kptr,elem(ie)%desc)

          !DBG print *,'pcg_solver: point #11'
          kptr=0
          !call edgerotate(edge2,2*nlev,kptr,elem(ie)%desc)

       end do
       while_iter = while_iter + 1 

       call bndry_exchangeV(cg%hybrid,edge2)

       do ie=nets,nete
          ieptr=ie-nets+1

          rmp     => elem(ie)%rmp
          Dinv    => elem(ie)%Dinv
          metdet => elem(ie)%metdet
          mp      => elem(ie)%mp

          kptr=0
          call edgeVunpack(edge2, gradp(1,1,1,1,ie), 2*nlev, kptr, elem(ie)%desc)

#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,i,j,gradp1,gradp2,div,deriv,iptr)
#endif
          do k=1,nlev
             if (.not.cg%converged(k)) then

                ! ====================
                ! 2*np*np Flops
                ! ====================

                do j=1,np
                   do i=1,np
                      gradp1 = gradp(i,j,1,k,ie)
                      gradp2 = gradp(i,j,2,k,ie)
                      gradp(i,j,1,k,ie) = rmp(i,j)*&
                                   (Dinv(1,1,i,j)*gradp1+Dinv(1,2,i,j)*gradp2)
                      gradp(i,j,2,k,ie) = rmp(i,j)*&
                                   (Dinv(2,1,i,j)*gradp1+Dinv(2,2,i,j)*gradp2)
                   end do
                end do

                ! ================================================
                ! Compute  Pseudo Laplacian(p), store in div
                !JMD   2*np*np*(np + np) Flops 
                ! ================================================

                div(:,:) = divergence(gradp(:,:,:,k,ie),deriv)*rrearth

                ! ==========================================
                ! compute Helmholtz operator, store in wrk3
                !  4*np*np Flops
                ! ==========================================

                iptr=1
                do j=1,np
                   do i=1,np
                      cg%state(ieptr)%s(iptr,k) = mp(i,j)*(metdet(i,j)*cg%state(ieptr)%z(iptr,k) + lambdasq(k)*div(i,j))
                      iptr=iptr+1
                   end do
                end do

             end if
          end do
       end do
#ifdef DEBUGOMP
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
#endif
       call t_stopf('timer_helm')

    end do  ! CG solver while loop


    ! ===============================
    ! Converged! unpack wrk3 into x
    ! ===============================

    do ie=nets,nete
       ieptr=ie-nets+1
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,i,j,iptr)
#endif
       do k=1,nlev
          iptr=1
          do j=1,np
             do i=1,np
#ifdef _NEWSTRUCT
                x(i,j,k,ie)=cg%state(ieptr,k)%x(iptr)
#else
                x(i,j,k,ie)=cg%state(ieptr)%x(iptr,k)
#endif
                iptr=iptr+1
             end do
          end do
       end do
    end do

#endif
    call t_stopf('pcg_solver_stag')
! CAM
#endif 
  end function pcg_solver_stag

  ! ================================================
  ! pcg_solver_nonstag:
  !
  ! Preconditioned conjugate gradient solver on the
  ! Gauss-Lobatto nonstaggered grid (np = np).
  ! 
  ! ================================================

  function pcg_solver_nonstag(elem,       &
       rhs,        &
       cg,         &
       red,        &
       edge1,      &
       edge2,      &   
       lambdasq,   &   
       deriv,      &   
       nets,       & 
       nete,       &
       blkjac) result(x) 
    use dimensions_mod, only : nlev, np, npsq
    use element_mod, only : element_t
    use reduction_mod, only : reductionbuffer_ordered_1d_t
    use cg_mod, only : cg_t, congrad
    use edge_mod, only : edgebuffer_t, edgevpack, edgevunpack!, edgerotate
    use derivative_mod, only : derivative_t, gradient_wk, gradient, divergence
    use control_mod, only : maxits, while_iter, tol, precon_method
    use physical_constants, only : rrearth
    use bndry_mod, only : bndry_exchangeV
    use linear_algebra_mod, only : matvec
    use parallel_mod, only : haltmp

    integer, intent(in)  :: nets,nete
    type(element_t), intent(in), target :: elem(:)
    real (kind=real_kind), intent(in) :: rhs(np,np,nlev,nets:nete) ! right hand side of operator
    type (cg_t)                       :: cg             ! conjugate gradient    (private)
    type (ReductionBuffer_ordered_1d_t)  :: red         ! CG reduction buffer   (shared memory)
    type (EdgeBuffer_t)               :: edge1          ! Laplacian divergence edge buffer (shared memory)
    type (EdgeBuffer_t)               :: edge2          ! Laplacian gradient edge buffer (shared memory)
    real (kind=real_kind), intent(in) :: lambdasq(nlev) ! Helmholtz lengthscale (private)
    type (derivative_t)               :: deriv          ! non staggered derivative struct     (private)
    type (blkjac_t)		      :: blkjac(nets:nete)

    real (kind=real_kind)             :: x(np,np,nlev,nets:nete)     ! solution (result)

    ! ===========
    ! Local
    ! ===========
#ifdef CAM
    call haltmp('semi-implicit method not yet supported in cam')
#else

    real (kind=real_kind),pointer :: metinv(:,:,:,:), Dinv(:,:,:,:)
    real (kind=real_kind),pointer :: metdet(:,:)
    real (kind=real_kind),pointer :: rmetdet(:,:)
    real (kind=real_kind),pointer :: rmp(:,:)
    real (kind=real_kind),pointer :: mp(:,:)

    real (kind=real_kind) :: gradp(np,np,2,nlev,nets:nete)
    real (kind=real_kind) :: div(np,np,nlev,nets:nete)
    real (kind=real_kind) :: p(np,np)
    real (kind=real_kind) :: r(npsq)
    real (kind=real_kind) :: z(npsq)

    real (kind=real_kind) :: gradp1
    real (kind=real_kind) :: gradp2

    integer :: ie
    integer :: i,j,k
    integer :: kptr
    integer :: iptr
    integer :: ieptr

    ! ========================================
    ! Load the rhs in the cg struct...
    ! ========================================

    call t_startf('pcg_solver_nonstag')

    do ie=nets,nete
       ieptr=ie-nets+1
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,i,j,iptr)
#endif
       do k=1,nlev
          iptr=1
          do j=1,np
             do i=1,np
                cg%state(ieptr)%r(iptr,k) = rhs(i,j,k,ie)
                iptr=iptr+1
             end do
          end do
       end do
    end do

    !cg%debug_level=1
    do while (congrad(cg,red,maxits,tol))
       call t_startf('timer_helm')

       do ie=nets,nete
          ieptr=ie-nets+1
          Dinv => elem(ie)%Dinv
          metinv => elem(ie)%metinv
          metdet => elem(ie)%metdet
          rmetdet  => elem(ie)%rmetdet

#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,i,j,deriv,iptr,gradp1,gradp2)
#endif
          do k=1,nlev
             if (.not.cg%converged(k)) then

                ! ========================================
                ! Apply preconditioner: wrk2 = (M^-1) * wrk1 
                ! ========================================

                if (precon_method == "block_jacobi") then

                   if (blkjac_storage == "LUfactor") then
                      call abortmp( 'dgesl needs linpack')
!                      call dgesl(cg%state(ieptr)%r(:,k), &
!                                 cg%state(ieptr)%z(:,k), &
!                                 blkjac(ie)%E(:,:,k),    &
!                                 blkjac(ie)%ipvt(:,k),   &
!                                 npsq)
                   else if (blkjac_storage == "inverse") then
                      call matvec(cg%state(ieptr)%r(:,k),cg%state(ieptr)%z(:,k),blkjac(ie)%E(:,:,k),npsq)
                   end if

                   !                   iptr=1
                   !                   do j=1,np
                   !                      do i=1,np
                   !                         cg%wrk2(iptr,k,ieptr)=z(iptr)
                   !                         p(i,j) = cg%state(ieptr)%z(iptr,k)
                   !                         iptr=iptr+1
                   !                      end do
                   !                   end do

                else if (precon_method == "identity") then

                   iptr=1
                   do j=1,np
                      do i=1,np
                         cg%state(ieptr)%z(iptr,k) = cg%state(ieptr)%r(iptr,k)*rmetdet(i,j)
                         iptr=iptr+1
                      end do
                   end do


                end if

                !JMD===========================================
                !JMD   2*np*np*(np + np) Flops 
   		!JMD  SR = (4*np*np + 2*np*np + np*np)*Ld
   		!JMD  SUM(WS) = (6*np*np + 2*np*np + np*np
                !JMD===========================================

#ifdef _WK_GRAD
                gradp(:,:,:,k,ie)=gradient_wk(reshape(cg%state(ieptr)%z(:,k),(/np,np/)),deriv)*rrearth
#else
                gradp(:,:,:,k,ie)=gradient(cg%state(ieptr)%z(:,k),deriv)*rrearth
#endif


                ! =======================================
                ! rotate gradient to form contravariant
                !JMD  4*np*np Flops
                ! =======================================

                do j=1,np
                   do i=1,np
                      gradp1       = gradp(i,j,1,k,ie)
                      gradp2       = gradp(i,j,2,k,ie)
#if 1
                      gradp(i,j,1,k,ie) = Dinv(1,1,i,j)*gradp1 + &
                           Dinv(2,1,i,j)*gradp2
                      gradp(i,j,2,k,ie) = Dinv(1,2,i,j)*gradp1 + &
                           Dinv(2,2,i,j)*gradp2
#else
                      gradp(i,j,1,k,ie) = (metinv(1,1,i,j)*gradp1 + metinv(1,2,i,j)*gradp2)
                      gradp(i,j,2,k,ie) = (metinv(2,1,i,j)*gradp1 + metinv(2,2,i,j)*gradp2)
#endif
                   end do
                end do
             end if
          end do

          kptr=0
          call edgeVpack(edge2,gradp(1,1,1,1,ie),2*nlev,kptr,elem(ie)%desc)

       end do

       while_iter = while_iter + 1 

       call bndry_exchangeV(cg%hybrid,edge2)

       do ie=nets,nete
          ieptr=ie-nets+1
          
          rmp     => elem(ie)%rmp
          Dinv    => elem(ie)%Dinv
          metdet => elem(ie)%metdet
          mp      => elem(ie)%mp

          kptr=0
          call edgeVunpack(edge2, gradp(1,1,1,1,ie), 2*nlev, kptr, elem(ie)%desc)
#ifdef DEBUGOMP
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
#endif
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,i,j,gradp1,gradp2,deriv)
#endif
          do k=1,nlev
             if (.not.cg%converged(k)) then

                ! ====================
                ! 2*np*np Flops
                ! ====================

                do j=1,np
                   do i=1,np
                      gradp1 = gradp(i,j,1,k,ie)
                      gradp2 = gradp(i,j,2,k,ie)
                      gradp(i,j,1,k,ie) = metdet(i,j) * rmp(i,j) * &
                             (Dinv(1,1,i,j)*gradp1 + Dinv(1,2,i,j)*gradp2)
                      gradp(i,j,2,k,ie) = metdet(i,j) * rmp(i,j) * &
                             (Dinv(2,1,i,j)*gradp1 + Dinv(2,2,i,j)*gradp2)
                   end do
                end do

                ! ================================================
                ! Compute  Pseudo Laplacian(p), store in div
                !JMD   2*np*np*(np + np) Flops 
                ! ================================================

                div(:,:,k,ie) = divergence(gradp(:,:,:,k,ie),deriv)*rrearth

                do j=1,np
                   do i=1,np
                      div(i,j,k,ie) = mp(i,j)*div(i,j,k,ie)
                   end do
                end do
             end if
          end do

          kptr=0
          call edgeVpack(edge1, div(1,1,1,ie), nlev, kptr, elem(ie)%desc)

       end do


       call bndry_exchangeV(cg%hybrid,edge1)


       ! ==========================================
       ! compute Helmholtz operator, store in wrk3
       !  4*np*np Flops
       ! ==========================================

       do ie=nets,nete
          ieptr=ie-nets+1

          rmp      => elem(ie)%rmp
          rmetdet  => elem(ie)%rmetdet
          metdet   => elem(ie)%metdet

          kptr=0
          call edgeVunpack(edge1, div(1,1,1,ie), nlev, kptr, elem(ie)%desc)

#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,i,j,iptr)
#endif
          do k=1,nlev
             if (.not.cg%converged(k)) then

                iptr=1
                do j=1,np
                   do i=1,np
                      cg%state(ieptr)%s(iptr,k) = metdet(i,j)*cg%state(ieptr)%z(iptr,k)+lambdasq(k)*rmp(i,j)*div(i,j,k,ie)
                      iptr=iptr+1
                   end do
                end do

             end if
          end do
       end do
#ifdef DEBUGOMP
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
#endif
       call t_stopf('timer_helm')

    end do  ! CG solver while loop
!    print *,'CG inter = ',cg%iter
    ! ===============================
    ! Converged! unpack wrk3 into x
    ! ===============================

    do ie=nets,nete
       ieptr=ie-nets+1
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,i,j,iptr)
#endif
       do k=1,nlev
          iptr=1
          do j=1,np
             do i=1,np
                !                x(i,j,k,ie)=cg%wrk3(iptr,k,ieptr)
                x(i,j,k,ie)=cg%state(ieptr)%x(iptr,k)
                iptr=iptr+1
             end do
          end do
       end do
    end do
    call t_stopf('pcg_solver_nonstag')
#endif
  end function pcg_solver_nonstag

  subroutine blkjac_init(elem, deriv,lambdasq,nets,nete,blkjac)
    use element_mod, only : element_t
    use derivative_mod, only : derivative_t, gradient_wk, gradient, divergence, derivative_stag_t
    use physical_constants, only : rearth, rrearth
    use dimensions_mod, only : nlev, np
    use parallel_mod, only : haltmp
    type(element_t), intent(in), target :: elem(:)
    type (derivative_t)                  :: deriv
    real (kind=real_kind), intent(in)    :: lambdasq(nlev)
    integer(kind=int_kind), intent(in) :: nets
    integer(kind=int_kind), intent(in) :: nete
    type (blkjac_t)      , intent(inout) :: blkjac(nets:nete)
#if 0
    !JMD    real (kind=real_kind)             :: E(npsq,npsq,nlev,nets:nete)
    !JMD    integer                           :: ipvt(npsq,nlev,nets:nete)
#endif

    ! ===============
    ! Local Variables
    ! ===============

    integer :: ie,kk
    integer :: i,j,k
    integer :: iptr
    integer :: ieptr
    integer :: lwork
    integer :: info

    real (kind=real_kind) :: p(npsq)
    real (kind=real_kind) :: z(npsq)
    real (kind=real_kind) :: gradp(np,np,2) ! velocity buffer
    real (kind=real_kind) :: div(np,np)

    real (kind=real_kind), pointer :: mp(:, :)

    real (kind=real_kind), pointer :: rmp(:,:)
    real (kind=real_kind), pointer :: mv(:,:)

    real (kind=real_kind), pointer :: metdet(:,:)
    real (kind=real_kind), pointer :: rmetdet(:,:)
    real (kind=real_kind), pointer :: metinv(:,:,:,:)

    real (kind=real_kind) :: gradp1
    real (kind=real_kind) :: gradp2
    real (kind=real_kind) :: det(2)

    real (kind=real_kind) :: lsq

    ! =================================
    ! Begin executable code
    ! =================================    
#ifndef CAM
    !DBG print *,'blkjac_init: point #0 nets,nete',nets,nete
    !DBG print *,'blkjac_init: point #0.1 nets,nete',nets,nete

    lsq=rearth*rearth

    ! ===========================
    ! zero out E matrix
    ! ===========================
    !DBG print *,'blkjac_init: point #0.2 nets,nete',nets,nete
    do ie=nets,nete
       blkjac(ie)%E = 0.0d0
       blkjac(ie)%ipvt = 0
    enddo
    !DBG print *,'blkjac_init: point #0.3 nets,nete',nets,nete

    do k=1,nlev

       !JMDprint *,'blkjac_init: point #1'
       do ie=nets,nete

          !JMD print *,'blkjac_init: point #2 ie,nets,nete',ie,nets,nete
          !JMD print *,'blkjac_init: point #3 ie,nets,nete',ie,nets,nete


          metdet  => elem(ie)%metdet
          metinv  => elem(ie)%metinv
          rmp     => elem(ie)%rmp
          mp      => elem(ie)%mp          

#if (defined ELEMENT_OPENMP)
!$omp parallel do private(kk,p,gradp,gradp1,gradp2,deriv,i,j,div,iptr)
#endif
          do kk = 1, npsq     ! delta fn excitation index

             p(:) = 0.0d0
             p(kk)= 1.0d0 

             ! =========================================================
             ! Compute (weak form) pressure gradient(s) on velocity grid
             ! =========================================================

             gradp(:,:,:)= gradient_wk(reshape(p,(/np,np/)),deriv)*rrearth

             do j=1,np
                do i=1,np
                   gradp1       = gradp(i,j,1)
                   gradp2       = gradp(i,j,2)
                   gradp(i,j,1) = metdet(i,j)*(metinv(1,1,i,j)*gradp1 + metinv(1,2,i,j)*gradp2)
                   gradp(i,j,2) = metdet(i,j)*(metinv(2,1,i,j)*gradp1 + metinv(2,2,i,j)*gradp2)
                end do
             end do

             ! ================================================
             ! Apply inverse velocity mass matrix Mu^{-1}
             ! ================================================

             do j=2,np-1
                do i=2,np-1
                   gradp(i,j,1) = gradp(i,j,1)*rmp(i,j)
                   gradp(i,j,2) = gradp(i,j,2)*rmp(i,j)
                end do
             end do

             ! ================================================
             ! Zero (Neumann) pressure gradients along the boundaries
             !
             ! zero north/south edge
             ! ================================================

             do j=1,np
                gradp(1 ,j,1)=0.0D0
                gradp(1 ,j,2)=0.0D0
                gradp(np,j,1)=0.0D0
                gradp(np,j,2)=0.0D0
             end do

             ! ================================================
             ! zero east/west edge
             ! ================================================

             do i=1,np
                gradp(i,1 ,1)=0.0D0
                gradp(i,1 ,2)=0.0D0
                gradp(i,np,1)=0.0D0
                gradp(i,np,2)=0.0D0
             end do

             ! ==================================================
             ! Compute divergence on pressure grid
             ! ==================================================

             div = divergence(gradp,deriv)*rrearth

             ! ==================================================
             ! Compute Rt = Mp.(pp + dt^2.pmean.D.Mu^-1.D'.pp)
             ! ==================================================

             iptr=1
             do j=1,np
                do i=1,np
                   blkjac(ie)%E(iptr,kk,k) = metdet(i,j)*p(iptr) + lambdasq(k) * div(i,j)/(lsq*mv(i,j))
                   iptr=iptr+1
                end do
             end do

             !JMD print *,'blkjac_init: point #10'
          end do  ! end loop over kk (delta function exitation location)
          !DBG print *,'blkjac_init: point #11 ie,nets,nete',ie,nets,nete

          ! ======================================
          ! Lu factorize Block Jacobi matrix...
          ! ======================================
          !             print *, __FILE__,__LINE__,minval(blkjac(ie)%E(:,:,k)),maxval(blkjac(ie)%E(:,:,k))

          !          call dgefa(blkjac(ie)%E(1,1,k),npsq,npsq,blkjac(ie)%ipvt(1,k), info) !LINPACK
          call dgetrf(npsq,npsq,blkjac(ie)%E(1,1,k),npsq,blkjac(ie)%ipvt(1,k), info) !LAPACK

          !DBG print *,'blkjac_init: point #12 ie,nets,nete',ie,nets,nete

          if (info /= 0) then
             print *,"WARNING: singular Block Jacobi matrix detected during preconditioning setup"
          end if

          if (blkjac_storage == "inverse") then

             if (k==1 .and. ie == 1) then
                !debug                 print *,"storing Block Jacobi Inverse"
             end if

             ! =======================================
             ! invert E matrix on each processor
             ! =======================================


             !             call dgedi( blkjac(ie)%E(:,:,k), npsq, npsq, blkjac(ie)%ipvt(:,k), det, z(:), 1 ) !LINPACK
             lwork = npsq
             call dgetri(npsq, blkjac(ie)%E(1,1,k), npsq, blkjac(ie)%ipvt(1,k),z,lwork,info) !LAPACK


          else if (blkjac_storage == "LUfactor") then

             do i=1,npsq
                blkjac(ie)%E(i,i,k)=1.0D0/blkjac(ie)%E(i,i,k)
             end do

             if (k==1 .and. ie == 1) then
                print *,"LU Block Jacobi storage"
             end if

          else
             call haltmp("bad Block Jacobi storage option")
          end if

          !DBG print *,'blkjac_init: point #15'
       end do  ! end element loop
       !DBG print *,'blkjac_init: point #16'

    end do ! end level loop
    !DBG print *,'blkjac_init: point #17'
#endif
  end subroutine blkjac_init

end module solver_mod
