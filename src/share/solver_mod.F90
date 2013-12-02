#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#undef _DGEMV
module solver_mod
  use kinds, only : real_kind, int_kind
  use dimensions_mod, only : npsq, nlev
  use perf_mod, only: t_startf, t_stopf ! _EXTERNAL
  use parallel_mod, only : abortmp

  use ,intrinsic :: iso_c_binding


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
  public  :: solver_test

#ifdef TRILINOS
  public  :: solver_test_ml
#endif

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
    type (derivative_t), intent(in) :: deriv          ! Staggered derivative struct     (private)
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
#if (defined COLUMN_OPENMP)
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

#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,Trans,inc,i,j,gradp1,gradp2)
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
                      call dgemv(Trans,npsq,npsq,one,blkjac(ie)%E(:,:,k),npsq,&
                           cg%state(ieptr)%r(:,k),inc,zero,cg%state(ieptr)%z(:,k),inc)
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

#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,i,j,gradp1,gradp2,div,iptr)
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
#if (defined HORIZ_OPENMP)
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
#if (defined COLUMN_OPENMP)
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
    type (derivative_t), intent(in) :: deriv          ! non staggered derivative struct     (private)
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
#if (defined COLUMN_OPENMP)
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

#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,i,j,iptr,gradp1,gradp2)
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
#if (defined HORIZ_OPENMP)
!$OMP BARRIER
#endif
#endif
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,i,j,gradp1,gradp2)
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

#if (defined COLUMN_OPENMP)
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
#if (defined HORIZ_OPENMP)
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
#if (defined COLUMN_OPENMP)
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
    type (derivative_t), intent(in) :: deriv
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

#if (defined COLUMN_OPENMP)
!$omp parallel do private(kk,p,gradp,gradp1,gradp2,i,j,div,iptr)
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


#ifdef TRILINOS

  ! ================================================
  ! helm_graph:
  !
  ! ================================================
  subroutine helm_graph(N,global_index_mat,op_data_ptr) bind(C,name='homme_globalIDs')
  
    use dimensions_mod, only : nlev, np,npsq
    use derived_type_mod, only :derived_type
    implicit none
    integer(c_int), intent(in)      :: N
    integer(c_int), intent(inout)      :: global_index_mat(N)
    type(c_ptr)                      :: op_data_ptr
    type(derived_type),  pointer          :: fptr=>NULL()

    integer      :: temp_global_index_mat(N)


    integer  :: nets,nete

    ! ===========
    ! Local
    ! ===========
    integer :: ie
    integer :: i,j,k
    integer gindex

    
    call c_f_pointer(op_data_ptr,fptr) !convert c pointer to f pointer


    nets=fptr%nets
    nete=fptr%nete


    gindex=1
       do ie=nets,nete
          do k=1,nlev
             do j=1,np
                do i=1,np
                  global_index_mat(gindex)=fptr%base(ie)%gdofP(i,j) 
                  gindex=gindex+1;
                enddo
             enddo
          end do
       end do
  end subroutine helm_graph





  ! ================================================
  ! helm_mat:
  !
  !    L(x) = laplace_sphere_wk(x) = -< grad(PHI) dot grad(x) >
  !    <   > = spheremp weighted inner-product
  !    L is self-adjoint:  < L(x),y> = < x,L(y) >
  !
  ! return helmholtz matrix on element El
  !  form by applying to a unit vector x for each column:
  !     <PHI,x> + a*L(x) =  < PHI, rhs >        
  !     contant a ~ 10 dx^2 (typical scaling in semi-implicit solve)
  !
  ! 2D - but applied to every level  k=1,nlev
  !
  ! In matrix notation, following the convention in M.T. and A.L.'s 
  ! "implicit.pdf" notes:
  !   HelmOp*x=  D QQ^t (M + a*L ) x 
  ! with:
  !   M    = multiply by spheremp 
  !   QQ^t = pack, exchange, unpack
  !   D    = multiply by rspheremp  D = Q V^-1 Q^-L   
  !          where V = the SEM diagonal mass matrix acting on vectors with no
  !          duplicate degrees of freedom.  V does not appear in HOMME.
  !   L is self adjoint w.r.t. M:    L(x) M y = x M L(y)
  !
  !
  ! ================================================
  subroutine helm_mat(ie,ElMatSize,ElementMatOut,Indices,op_data_ptr) bind(C,name='helm_mat')
  
    use dimensions_mod, only : nlev, np,npsq
    use derivative_mod, only : derivative_t, laplace_sphere_wk
    use derived_type_mod, only :derived_type
    implicit none
    integer(c_int), intent(in),value :: ie 
    integer(c_int), intent(in),value :: ElMatSize
    real(c_double), intent(out)      :: ElementMatOut(ElMatSize,ElMatSize)
    integer(c_int), intent(out)      :: Indices(ElMatSize)
    type(c_ptr)                      :: op_data_ptr
    type(derived_type),  pointer          :: fptr=>NULL()

    integer  :: nets,nete
!    type (derivative_t)               :: deriv          ! non staggered derivative struct     (private)

    ! ===========
    ! Local
    ! ===========
    real (kind=real_kind), allocatable :: ElMatCol(:,:)
    real (kind=real_kind) :: x(np,np)
    real (kind=real_kind) :: alambda = 10*250e3**2      ! low res test, dx=250km grid spacing

    integer :: i,j,k
    integer ColIndex
    
    integer colindices(np,np)
    call c_f_pointer(op_data_ptr,fptr) !convert c pointer to f pointer

    ElementMatOut=0.0d0

    allocate(ElMatCol(np,np))


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Application Loop
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !do ColIndex=1,nlev*np**2
    do ColIndex=1,np**2
      do j=1,np
        do i=1,np
           x(i,j) = 0.0d0
           if(ColIndex == ((j-1)*np +i)) then 
                   x(i,j)=1 
           endif
        enddo
      enddo
             
             ! Apply x + laplace(x)
             ! weak laplace operator already includes mass
             ! so only multiply x by mass;

          !ElMatCol(:,:)= fptr%base(ie)%rspheremp(:,:)*fptr%base(ie)%spheremp(:,:)*x(:,:) 
          !ElMatCol(:,:)= fptr%base(ie)%spheremp(:,:)*x(:,:) 
          !ElMatCol(:,:)= x(:,:) 
          ElMatCol(:,:)= fptr%base(ie)%rspheremp(:,:)*(fptr%base(ie)%spheremp(:,:)*x(:,:) +&
            alambda*laplace_sphere_wk(x,fptr%deriv,fptr%base(ie),var_coef=.false.))
           



          ElementMatOut(:,ColIndex)=reshape(ElMatCol,(/np*np/))
           
!          ElementMatOut(1+(1-k)*np**2:1+k*np**2,ColIndex)=reshape(ElMatCol,(/np*np/))
    enddo

    k=1
    do j=1,np
      do i=1,np
         Indices(k)=fptr%base(ie)%gdofP(i,j)
         !Indices(k)=k
         k=k+1
      enddo
    enddo




!       write(6,*)'vecOut=[',vecOut,']'


!    if (hybrid%masterthread) print *,'applied helmholtz'
!    call flush(6)
       deallocate(ElMatCol)
  end subroutine helm_mat


  ! ================================================
  ! helm_rhs:
  !
  !    L(x) = laplace_sphere_wk(x) = -< grad(PHI) dot grad(x) >
  !    <   > = spheremp weighted inner-product
  !    L is self-adjoint:  < L(x),y> = < x,L(y) >
  !
  ! return helmholtz matrix on element El
  !  form by applying to a unit vector x for each column:
  !     <PHI,x> + a*L(x) =  < PHI, rhs >        
  !     contant a ~ 10 dx^2 (typical scaling in semi-implicit solve)
  !
  ! 2D - but applied to every level  k=1,nlev
  !
  ! In matrix notation, following the convention in M.T. and A.L.'s 
  ! "implicit.pdf" notes:
  !   HelmOp*x=  D QQ^t (M + a*L ) x 
  ! with:
  !   M    = multiply by spheremp 
  !   QQ^t = pack, exchange, unpack
  !   D    = multiply by rspheremp  D = Q V^-1 Q^-L   
  !          where V = the SEM diagonal mass matrix acting on vectors with no
  !          duplicate degrees of freedom.  V does not appear in HOMME.
  !   L is self adjoint w.r.t. M:    L(x) M y = x M L(y)
  !
  !
  ! ================================================
  subroutine helm_map(ie,ElMatSize,Indices,op_data_ptr) bind(C,name='helm_map')
  
 !~ elem,edge1,red,hybrid,deriv,nets,nete,vecIn)
    use dimensions_mod, only : nlev, np,npsq
!    use element_mod, only : element_t
!    use edge_mod, only : edgebuffer_t, edgevpack, edgevunpack!, edgerotate
!    use derivative_mod, only : derivative_t, laplace_sphere_wk
!    use control_mod, only : maxits, while_iter, tol, precon_method

!    use physical_constants, only : rrearth, dd_pi, rearth, omega
!    use bndry_mod, only : bndry_exchangeV
!    use linear_algebra_mod, only : matvec
!    use parallel_mod, only : haltmp
!    use hybrid_mod, only : hybrid_t
    use derived_type_mod, only :derived_type
    implicit none
    integer(c_int), intent(in),value :: ie 
    integer(c_int), intent(in),value :: ElMatSize
    !real(c_double), intent(out)      :: ElementMatOut(ElMatSize,ElMatSize)
    integer(c_int), intent(out)      :: Indices(ElMatSize)
    type(c_ptr)                      :: op_data_ptr
    type(derived_type),  pointer          :: fptr=>NULL()


    ! ===========
    ! Local
    ! ===========

    integer :: i,j,k
    
    call c_f_pointer(op_data_ptr,fptr) !convert c pointer to f pointer

     k=1
      do j=1,np
        do i=1,np
           Indices(k)=fptr%base(ie)%gdofP(i,j)
           k=k+1
        enddo
      enddo


  end subroutine helm_map




  ! ================================================
  ! solver_test_ml:
  !
  !    L(x) = laplace_sphere_wk(x) = -< grad(PHI) dot grad(x) >
  !    <   > = spheremp weighted inner-product
  !    L is self-adjoint:  < L(x),y> = < x,L(y) >
  !
  ! solve for x:
  !     <PHI,x> + a*L(x) =  < PHI, rhs >        
  !     contant a ~ 10 dx^2 (typical scaling in semi-implicit solve)
  !
  ! 2D solve - but applied to every level  k=1,nlev
  !
  ! In matrix notation, following the convention in M.T. and A.L.'s 
  ! "implicit.pdf" notes:
  !     D QQ^t (M + a*L ) x = rhs
  ! with:
  !   M    = multiply by spheremp 
  !   QQ^t = pack, exchange, unpack
  !   D    = multiply by rspheremp  D = Q V^-1 Q^-L   
  !          where V = the SEM diagonal mass matrix acting on vectors with no
  !          duplicate degrees of freedom.  V does not appear in HOMME.
  !   L is self adjoint w.r.t. M:    L(x) M y = x M L(y)
  !
  ! Note: if we solve laplace equation instead of Helmholz, we need to
  ! ensure < rhs,1>=0
  !
  ! ================================================
  subroutine solver_test_ml(elem,edge1,red,hybrid,deriv,nets,nete)
    use dimensions_mod, only : nlev, np,npsq
    use element_mod, only : element_t
    use reduction_mod, only : reductionbuffer_ordered_1d_t
    use cg_mod, only : cg_t, congrad, cg_create
    use edge_mod, only : edgebuffer_t, edgevpack, edgevunpack!, edgerotate
    use derivative_mod, only : derivative_t, laplace_sphere_wk
    use control_mod, only : maxits, while_iter, tol, precon_method
    use physical_constants, only : rrearth, dd_pi, rearth, omega
    use bndry_mod, only : bndry_exchangeV
    use linear_algebra_mod, only : matvec
    use parallel_mod, only : haltmp
    use hybrid_mod, only : hybrid_t
    use global_norms_mod, only : linf_snorm, l2_snorm
    use derived_type_mod, only : derived_type, initialize
    use time_mod, only : TimeLevel_t


    interface

     subroutine belosfinish()  bind(C,name='belosfinish')
       use ,intrinsic :: iso_c_binding ,only : c_double ,c_int ,c_ptr
     end subroutine belosfinish


     subroutine build_maps(vectorSize, nets,nete,np,nlev,datavector) &
                bind(C,name='BuildMaps')
       use ,intrinsic :: iso_c_binding ,only : c_double ,c_int ,c_ptr
            integer(c_int)                :: vectorSize
            integer(c_int)                :: nets,nete,np,nlev
            type(c_ptr)                   :: datavector
     end subroutine build_maps


     subroutine build_matrix(vectorSize, nets,nete,np,nlev,datavector) &
                bind(C,name='BuildMatrix')
       use ,intrinsic :: iso_c_binding ,only : c_double ,c_int ,c_ptr
            integer(c_int)                :: vectorSize
            integer(c_int)                :: nets,nete,np,nlev
            type(c_ptr)                   :: datavector
     end subroutine build_matrix



     subroutine build_rhs(nets,nete,np,nlev,rhs,datavector) &
                bind(C,name='BuildRHS')
       use ,intrinsic :: iso_c_binding ,only : c_double ,c_int ,c_ptr
            integer(c_int)                :: nets,nete,np,nlev
            real(c_double)  ,dimension(*) :: rhs
            type(c_ptr)                   :: datavector
     end subroutine build_rhs



     subroutine set_problem() &
                bind(C,name='SetProblem')
     end subroutine set_problem


     subroutine helmholtz_solve(vectorSize, lhs) &
                bind(C,name='HelmholtzSolve')
       use ,intrinsic :: iso_c_binding ,only : c_double ,c_int ,c_ptr
            integer(c_int)                :: vectorSize
            real(c_double)  ,dimension(*) :: lhs
     end subroutine helmholtz_solve

    end interface



    integer, intent(in)  :: nets,nete
    type(element_t), intent(in), target :: elem(:)
    type (ReductionBuffer_ordered_1d_t)  :: red         ! CG reduction buffer   (shared memory)
    type (derivative_t)               :: deriv          ! non staggered derivative struct     (private)
    type (hybrid_t)             :: hybrid
    type (EdgeBuffer_t)               :: edge1          ! Laplacian divergence edge buffer (shared memory)


    ! ===========
    ! Local
    ! ===========
    type (cg_t)                       :: cg             ! conjugate gradient    (private)
    real (kind=real_kind) :: LHS(np,np,nlev,nets:nete)
    real (kind=real_kind) :: RHS(np,np,nlev,nets:nete)
    real (kind=real_kind) :: sol(np,np,nlev,nets:nete)   ! exact solution
    real (kind=real_kind) :: solver_wts(npsq,nete-nets+1)
    real (kind=real_kind) :: x(np,np)
    real (kind=real_kind) :: alambda = 10*250e3**2      ! low res test, dx=250km grid spacing

    integer :: ie
    integer :: i,j,k
    integer :: kptr
    integer :: iptr
    integer :: ieptr
    real (kind=real_kind) :: snlat,cslat,cslon,snlon,xc,yc,zc, res, res_sol
    integer rhsVectorSize
    real(kind=real_kind):: rhsVector(np*np*nlev*(nete-nets+1))
    real(kind=real_kind):: solVector(np*np*nlev*(nete-nets+1))

    type(derived_type) ,target         :: helmholtzdata
    type(derived_type) ,pointer        :: f_ptr_helmholtzdata=>NULL()
    type(c_ptr)                        :: c_ptr_helmholtzdata

    integer linVecndx


    integer lenx
    real (kind=real_kind) pmean,dt
    type (TimeLevel_t) tl
    lenx=0
    pmean=0.0d0
    dt=0.0d0
    



!!  c_ptr_helmholtzdata should contain everything needed to evaluate the linear helmholtz operator

    call initialize(helmholtzdata, lenx, elem, pmean,edge1,edge1, edge1, &
            hybrid, deriv, dt, tl, nets, nete)


    f_ptr_helmholtzdata => helmholtzdata
    c_ptr_helmholtzdata =  c_loc(f_ptr_helmholtzdata)
!




    rhsVectorSize=np*np*nlev*(nete-nets+1)

    solVector=0.0d0



    if (hybrid%masterthread) print *,'creating manufactured solution'
    do ie=nets,nete
       iptr=1
       do j=1,np
          do i=1,np
             solver_wts(iptr,ie-nets+1) = elem(ie)%spheremp(i,j)
             iptr=iptr+1
          end do
       end do
    end do
    call cg_create(cg, npsq, nlev, nete-nets+1, hybrid, 0, solver_wts)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! make up an exact solution
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do ie=nets,nete
       ieptr=ie-nets+1
       do k=1,nlev
          iptr=1
          do j=1,np
             do i=1,np
                snlat = SIN(elem(ie)%spherep(i,j)%lat)
                cslat = COS(elem(ie)%spherep(i,j)%lat)
                snlon = SIN(elem(ie)%spherep(i,j)%lon)
                cslon = COS(elem(ie)%spherep(i,j)%lon)
  
                xc = cslat*cslon
                yc = cslat*snlon
                zc = snlat

                ! take a couple of low-freq spherical harmonics for the solution
                sol(i,j,k,ie) = 1*xc + 2*yc + 3*zc + 4*xc*yc + 5*xc*zc + 6*yc*zc + &
                     7*(xc*yc*zc) 

             end do
          end do
       end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! compute the RHS from our exact solution
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do ie=nets,nete
       do k=1,nlev
          RHS(:,:,k,ie)=elem(ie)%spheremp(:,:)*sol(:,:,k,ie) + &
               alambda*laplace_sphere_wk(sol(:,:,k,ie),deriv,elem(ie),var_coef=.false.)
          call edgeVpack(edge1, RHS(1,1,1,ie), nlev, 0, elem(ie)%desc)
       end do
    end do
    call bndry_exchangeV(cg%hybrid,edge1)

    linVecndx=1

    do ie=nets,nete
       ! unpack RHS
       call edgeVunpack(edge1, RHS(1,1,1,ie), nlev, 0, elem(ie)%desc)
       do k=1,nlev
          RHS(:,:,k,ie)=RHS(:,:,k,ie)*elem(ie)%rspheremp(:,:)
       enddo
       !
       !  Initialize CG solver:  set %r = residual from initial guess
       !  if initial guess = 0, then we take %r=RHS
       !
       ieptr=ie-nets+1
       do k=1,nlev
          iptr=1
          do j=1,np
             do i=1,np
                cg%state(ieptr)%r(iptr,k) = rhs(i,j,k,ie)
                iptr=iptr+1
                rhsVector(linVecndx)=rhs(i,j,k,ie)
                ! rhsVector(linVecndx)=1.0d0
                ! rhsVector(linVecndx)=sol(i,j,k,ie) 
                linVecndx=linVecndx+1
             enddo
          enddo
       enddo
    enddo




  ! construt maps to build trilinos Matrix and rhs vectors
  call build_maps(rhsVectorSize, nets,nete,np,nlev,c_ptr_helmholtzdata) 
  ! assemble trilinos RHS vector from HOMME rhs vector - unique dof's across processors
  call build_rhs(nets,nete,np,nlev,rhsVector,c_ptr_helmholtzdata) 
  ! assemble Trilinos Matrix - unique dof's across processors
  call build_matrix(rhsVectorSize, nets,nete,np,nlev,c_ptr_helmholtzdata) 
  ! configure linear solver with ML preconditioner for solving Helmholtz
  ! with Belos iterative solver
  call set_problem() 
  ! Solve Helmholtz equation with Trilinos return solution in solVector
  call helmholtz_solve(rhsVectorSize, solVector) 


! an example of building another rhs and calling solver again
!
!    linVecndx=1
!    do linVecndx=1,rhsVectorSize
!           rhsVector(linVecndx)=1.0*rhsVector(linVecndx)
!    enddo
!
!
!  for each rhs vector send to "build_rhs" and then call "helmholtz_solve" 
!   
!  call build_rhs(nets,nete,np,nlev,rhsVector,c_ptr_helmholtzdata) 
!       call t_startf('timer_trilinosmlhelm')
!  call helmholtz_solve(rhsVectorSize, solVector) 
!       call t_stopf('timer_trilinosmlhelm')

! call belosfinish to delete memory for the trilinos Matrix maps and Helmholtz
! solver
  !end second rhs test

  call belosfinish()

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Solver Loop
    ! in this version, we keep the residual C0 by DSSing the LHS
    ! and initializiont with a C0 RHS.  
    ! the update x, based on the last residual, will already be C0
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    cg%debug_level=1  ! 2=output every iterations
    maxits = 250
    tol=1d-10
    if (hybrid%masterthread) print *,'running solver V1 (C0 RHS) tol=',tol
    do while (congrad(cg,red,maxits,tol))
       do ie=nets,nete
          ieptr=ie-nets+1
          do k=1,nlev
             ! apply preconditioner here:
             cg%state(ieptr)%z(:,k) = cg%state(ieptr)%r(:,k)
             
             !reshape(cg%state(ieptr)%z(:,k),(/np,np/))
             iptr=1
             do j=1,np
                do i=1,np
                   x(i,j) = cg%state(ieptr)%z(iptr,k)
                   iptr=iptr+1
                enddo
             enddo
             
             ! solve x + laplace(x)
             ! weak laplace operator already includes mass
             ! so only multiply x by mass;
             !LHS(:,:,k,ie)=elem(ie)%spheremp(:,:)*x(:,:) + &
             !     alambda*laplace_sphere_wk(x,deriv,elem(ie),var_coef=.false.)

             LHS(:,:,k,ie)=elem(ie)%rspheremp(:,:)*(elem(ie)%spheremp(:,:)*x(:,:) + &
                  alambda*laplace_sphere_wk(x,deriv,elem(ie),var_coef=.false.))
             !LHS(:,:,k,ie)=elem(ie)%rspheremp(:,:)*(elem(ie)%spheremp(:,:)*x(:,:) )
             
          end do
          call edgeVpack(edge1, LHS(1,1,1,ie), nlev, 0, elem(ie)%desc)
       end do
       call bndry_exchangeV(cg%hybrid,edge1)
       do ie=nets,nete
          ! unpack LHS
          call edgeVunpack(edge1, LHS(1,1,1,ie), nlev, 0, elem(ie)%desc)

          ieptr=ie-nets+1
          do k=1,nlev
             iptr=1
             do j=1,np
                do i=1,np
                   cg%state(ieptr)%s(iptr,k) = LHS(i,j,k,ie)
                   iptr=iptr+1
                end do
             end do
          enddo
       enddo
       
    end do  ! CG solver while loop


    ! ===============================
    ! Converged! compute actual error (not residual computed in solver)
    ! ===============================
    do ie=nets,nete
       ieptr=ie-nets+1
       do k=1,nlev
          iptr=1
          do j=1,np
             do i=1,np
                LHS(i,j,k,ie) = cg%state(ieptr)%x(iptr,k)
                iptr=iptr+1
             enddo
          enddo
       enddo
    enddo
    res = l2_snorm(elem,LHS,sol,hybrid,np,nets,nete) 
    if (hybrid%masterthread) print *,'normalized l2 error= ',res
    res = linf_snorm(LHS,sol,hybrid,np,nets,nete) 
    if (hybrid%masterthread) print *,'normalized linf error= ',res

    tol=tol/2
!!PAL
tol=1.e-12
    if (hybrid%masterthread) print *,'running solver V2 (DG RHS) tol=',tol
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! VERSION 2.  SAVE 1 DSS
    ! (important since we need to get iterations down to about 5 to be competitive)
    ! compute the RHS from our exact solution 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! We should not have to apply DSS to RHS, since our equation:
    ! < PHI, lap(x) > = < PHI, DSS(RHS)>
    ! But since DSS(PHI)=PHI, and DSS is self adjoint,
    ! < PHI, DSS(RHS)>= < PHI, RHS>
    !
    ! in this version, we do not need to DSS the RHS (or LHS). 
    ! But the update x needs to be C0, so we DSS only x:
    !
    ! ONE ISSUE:  tolerence is based on <RHS,RHS>, which is not 
    ! computed correctly (and will always be larger than <DSS(RHS),DSS(RHS)>
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do ie=nets,nete
       ! note: weak from laplace operator includes mass, so remove it
       do k=1,nlev
          RHS(:,:,k,ie)=sol(:,:,k,ie) + &
               alambda*laplace_sphere_wk(sol(:,:,k,ie),deriv,elem(ie),var_coef=.false.)&
               / elem(ie)%spheremp(:,:)
       end do
       !
       !  Initialize CG solver:  set %r = residual from initial guess
       !  if initial guess = 0, then we take %r=RHS
       ieptr=ie-nets+1
       do k=1,nlev
          iptr=1
          do j=1,np
             do i=1,np
                cg%state(ieptr)%r(iptr,k) = rhs(i,j,k,ie)
                iptr=iptr+1
             enddo
          enddo
       enddo
    enddo
    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Solver Loop 2
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do while (congrad(cg,red,maxits,tol))
       do ie=nets,nete
          ieptr=ie-nets+1
          do k=1,nlev
             ! apply preconditioner here: (note: r is not C0)
             cg%state(ieptr)%z(:,k) = cg%state(ieptr)%r(:,k)
             
             !reshape(cg%state(ieptr)%z(:,k),(/np,np/))
             iptr=1
             do j=1,np
                do i=1,np
                   x(i,j) = cg%state(ieptr)%z(iptr,k)
                   iptr=iptr+1
                enddo
             enddo
             
             ! DSS x to make it C0.  use LHS for storage:
             LHS(:,:,k,ie)=x(:,:)*elem(ie)%spheremp(:,:)
             
          end do
          call edgeVpack(edge1, LHS(1,1,1,ie), nlev, 0, elem(ie)%desc)
       end do
       call bndry_exchangeV(cg%hybrid,edge1)
       do ie=nets,nete
          ! unpack LHS
          call edgeVunpack(edge1, LHS(1,1,1,ie), nlev, 0, elem(ie)%desc)
          do k=1,nlev
             LHS(:,:,k,ie)=LHS(:,:,k,ie)*elem(ie)%rspheremp(:,:)
          enddo

          ieptr=ie-nets+1
          do k=1,nlev
             x(:,:)=LHS(:,:,k,ie) ! x() is now C0

             ! compute LHS(x) = x + laplace(x)
             LHS(:,:,k,ie)=x(:,:) + &
                  alambda*laplace_sphere_wk(x,deriv,elem(ie),var_coef=.false.)&
                  /elem(ie)%spheremp(:,:)

             iptr=1
             do j=1,np
                do i=1,np
                   cg%state(ieptr)%s(iptr,k) = LHS(i,j,k,ie)    ! new LHS, DG
                   cg%state(ieptr)%z(iptr,k) = x(i,j)           ! z must be C0
                   iptr=iptr+1
                end do
             end do
          enddo
       enddo
       
    end do  ! CG solver while loop
    !print *,'solver test CG iter = ',cg%iter


    ! ===============================
    ! Converged! compute actual error (not residual computed in solver)
    ! ===============================

    if (hybrid%masterthread) print *,"Trilinos Solution" , "             ", "HOMME Solution","              ","Exact Solution"
    linVecndx=1
    do ie=nets,nete
       ieptr=ie-nets+1
       do k=1,nlev
          iptr=1
          do j=1,np
             do i=1,np
                LHS(i,j,k,ie) = cg%state(ieptr)%x(iptr,k)
                if (hybrid%masterthread) print *,solVector(linVecndx) , " ", LHS(i,j,k,ie)," ",sol(i,j,k,ie)
                !set LHS to the solVector to test the l2 norm
                LHS(i,j,k,ie)=solVector(linVecndx)
                linVecndx=linVecndx+1
                iptr=iptr+1
             enddo
          enddo
       enddo
    enddo

    res = l2_snorm(elem,LHS,sol,hybrid,np,nets,nete) 
    if (hybrid%masterthread) print *,'normalized l2 error= ',res
    res = linf_snorm(LHS,sol,hybrid,np,nets,nete) 
    if (hybrid%masterthread) print *,'normalized linf error= ',res

    deallocate(helmholtzdata%base)

  end subroutine solver_test_ml

#endif 


  ! ================================================
  ! solver_test:
  !
  !    L(x) = laplace_sphere_wk(x) = -< grad(PHI) dot grad(x) >
  !    <   > = spheremp weighted inner-product
  !    L is self-adjoint:  < L(x),y> = < x,L(y) >
  !
  ! solve for x:
  !     <PHI,x> + a*L(x) =  < PHI, rhs >        
  !     contant a ~ 10 dx^2 (typical scaling in semi-implicit solve)
  !
  ! 2D solve - but applied to every level  k=1,nlev
  !
  ! In matrix notation, following the convention in M.T. and A.L.'s 
  ! "implicit.pdf" notes:
  !     D QQ^t (M + a*L ) x = rhs
  ! with:
  !   M    = multiply by spheremp 
  !   QQ^t = pack, exchange, unpack
  !   D    = multiply by rspheremp  D = Q V^-1 Q^-L   
  !          where V = the SEM diagonal mass matrix acting on vectors with no
  !          duplicate degrees of freedom.  V does not appear in HOMME.
  !   L is self adjoint w.r.t. M:    L(x) M y = x M L(y)
  !
  ! Note: if we solve laplace equation instead of Helmholz, we need to
  ! ensure < rhs,1>=0
  !
  ! ================================================
  subroutine solver_test(elem,edge1,red,hybrid,deriv,nets,nete)
    use dimensions_mod, only : nlev, np,npsq
    use element_mod, only : element_t
    use reduction_mod, only : reductionbuffer_ordered_1d_t
    use cg_mod, only : cg_t, congrad, cg_create
    use edge_mod, only : edgebuffer_t, edgevpack, edgevunpack!, edgerotate
    use derivative_mod, only : derivative_t, laplace_sphere_wk
    use control_mod, only : maxits, while_iter, tol, precon_method
    use physical_constants, only : rrearth, dd_pi, rearth, omega
    use bndry_mod, only : bndry_exchangeV
    use linear_algebra_mod, only : matvec
    use parallel_mod, only : haltmp
    use hybrid_mod, only : hybrid_t
    use global_norms_mod, only : linf_snorm, l2_snorm

    integer, intent(in)  :: nets,nete
    type(element_t), intent(in), target :: elem(:)
    type (ReductionBuffer_ordered_1d_t)  :: red         ! CG reduction buffer   (shared memory)
    type (derivative_t), intent(in) :: deriv          ! non staggered derivative struct     (private)
    type (hybrid_t)             :: hybrid
    type (EdgeBuffer_t)               :: edge1          ! Laplacian divergence edge buffer (shared memory)


    ! ===========
    ! Local
    ! ===========
    type (cg_t)                       :: cg             ! conjugate gradient    (private)
    real (kind=real_kind) :: LHS(np,np,nlev,nets:nete)
    real (kind=real_kind) :: RHS(np,np,nlev,nets:nete)
    real (kind=real_kind) :: sol(np,np,nlev,nets:nete)   ! exact solution
    real (kind=real_kind) :: solver_wts(npsq,nete-nets+1)
    real (kind=real_kind) :: x(np,np)
    real (kind=real_kind) :: alambda = 10*250e3**2      ! low res test, dx=250km grid spacing

    integer :: ie
    integer :: i,j,k
    integer :: kptr
    integer :: iptr
    integer :: ieptr
    real (kind=real_kind) :: snlat,cslat,cslon,snlon,xc,yc,zc, res, res_sol

    if (hybrid%masterthread) print *,'creating manufactured solution'
    do ie=nets,nete
       iptr=1
       do j=1,np
          do i=1,np
             solver_wts(iptr,ie-nets+1) = elem(ie)%spheremp(i,j)
             iptr=iptr+1
          end do
       end do
    end do
    call cg_create(cg, npsq, nlev, nete-nets+1, hybrid, 0, solver_wts)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! make up an exact solution
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do ie=nets,nete
       ieptr=ie-nets+1
       do k=1,nlev
          iptr=1
          do j=1,np
             do i=1,np
                snlat = SIN(elem(ie)%spherep(i,j)%lat)
                cslat = COS(elem(ie)%spherep(i,j)%lat)
                snlon = SIN(elem(ie)%spherep(i,j)%lon)
                cslon = COS(elem(ie)%spherep(i,j)%lon)
  
                xc = cslat*cslon
                yc = cslat*snlon
                zc = snlat

                ! take a couple of low-freq spherical harmonics for the solution
                sol(i,j,k,ie) = 1*xc + 2*yc + 3*zc + 4*xc*yc + 5*xc*zc + 6*yc*zc + &
                     7*(xc*yc*zc) 

             end do
          end do
       end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! compute the RHS from our exact solution
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do ie=nets,nete
       do k=1,nlev
          RHS(:,:,k,ie)=elem(ie)%spheremp(:,:)*sol(:,:,k,ie) + &
               alambda*laplace_sphere_wk(sol(:,:,k,ie),deriv,elem(ie),var_coef=.false.)
          call edgeVpack(edge1, RHS(1,1,1,ie), nlev, 0, elem(ie)%desc)
       end do
    end do
    call bndry_exchangeV(cg%hybrid,edge1)
    do ie=nets,nete
       ! unpack RHS
       call edgeVunpack(edge1, RHS(1,1,1,ie), nlev, 0, elem(ie)%desc)
       do k=1,nlev
          RHS(:,:,k,ie)=RHS(:,:,k,ie)*elem(ie)%rspheremp(:,:)
       enddo
       !
       !  Initialize CG solver:  set %r = residual from initial guess
       !  if initial guess = 0, then we take %r=RHS
       !
       ieptr=ie-nets+1
       do k=1,nlev
          iptr=1
          do j=1,np
             do i=1,np
                cg%state(ieptr)%r(iptr,k) = rhs(i,j,k,ie)
                iptr=iptr+1
             enddo
          enddo
       enddo
    enddo
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Solver Loop
    ! in this version, we keep the residual C0 by DSSing the LHS
    ! and initializiont with a C0 RHS.  
    ! the update x, based on the last residual, will already be C0
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    cg%debug_level=1  ! 2=output every iterations
    maxits = 250
    tol=1d-7
    if (hybrid%masterthread) print *,'running solver V1 (C0 RHS) tol=',tol
    do while (congrad(cg,red,maxits,tol))
       do ie=nets,nete
          ieptr=ie-nets+1
          do k=1,nlev
             ! apply preconditioner here:
             cg%state(ieptr)%z(:,k) = cg%state(ieptr)%r(:,k)
             
             !reshape(cg%state(ieptr)%z(:,k),(/np,np/))
             iptr=1
             do j=1,np
                do i=1,np
                   x(i,j) = cg%state(ieptr)%z(iptr,k)
                   iptr=iptr+1
                enddo
             enddo
             
             ! solve x + laplace(x)
             ! weak laplace operator already includes mass
             ! so only multiply x by mass;
             LHS(:,:,k,ie)=elem(ie)%spheremp(:,:)*x(:,:) + &
                  alambda*laplace_sphere_wk(x,deriv,elem(ie),var_coef=.false.)
             
          end do
          call edgeVpack(edge1, LHS(1,1,1,ie), nlev, 0, elem(ie)%desc)
       end do
       call bndry_exchangeV(cg%hybrid,edge1)
       do ie=nets,nete
          ! unpack LHS
          call edgeVunpack(edge1, LHS(1,1,1,ie), nlev, 0, elem(ie)%desc)
          do k=1,nlev
             LHS(:,:,k,ie)=LHS(:,:,k,ie)*elem(ie)%rspheremp(:,:)
          enddo

          ieptr=ie-nets+1
          do k=1,nlev
             iptr=1
             do j=1,np
                do i=1,np
                   cg%state(ieptr)%s(iptr,k) = LHS(i,j,k,ie)
                   iptr=iptr+1
                end do
             end do
          enddo
       enddo
       
    end do  ! CG solver while loop


    ! ===============================
    ! Converged! compute actual error (not residual computed in solver)
    ! ===============================
    do ie=nets,nete
       ieptr=ie-nets+1
       do k=1,nlev
          iptr=1
          do j=1,np
             do i=1,np
                LHS(i,j,k,ie) = cg%state(ieptr)%x(iptr,k)
                iptr=iptr+1
             enddo
          enddo
       enddo
    enddo
    res = l2_snorm(elem,LHS,sol,hybrid,np,nets,nete) 
    if (hybrid%masterthread) print *,'normalized l2 error= ',res
    res = linf_snorm(LHS,sol,hybrid,np,nets,nete) 
    if (hybrid%masterthread) print *,'normalized linf error= ',res

    tol=tol/2
    if (hybrid%masterthread) print *,'running solver V2 (DG RHS) tol=',tol
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! VERSION 2.  SAVE 1 DSS
    ! (important since we need to get iterations down to about 5 to be competitive)
    ! compute the RHS from our exact solution 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! We should not have to apply DSS to RHS, since our equation:
    ! < PHI, lap(x) > = < PHI, DSS(RHS)>
    ! But since DSS(PHI)=PHI, and DSS is self adjoint,
    ! < PHI, DSS(RHS)>= < PHI, RHS>
    !
    ! in this version, we do not need to DSS the RHS (or LHS). 
    ! But the update x needs to be C0, so we DSS only x:
    !
    ! ONE ISSUE:  tolerence is based on <RHS,RHS>, which is not 
    ! computed correctly (and will always be larger than <DSS(RHS),DSS(RHS)>
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do ie=nets,nete
       ! note: weak from laplace operator includes mass, so remove it
       do k=1,nlev
          RHS(:,:,k,ie)=sol(:,:,k,ie) + &
               alambda*laplace_sphere_wk(sol(:,:,k,ie),deriv,elem(ie),var_coef=.false.)&
               / elem(ie)%spheremp(:,:)
       end do
       !
       !  Initialize CG solver:  set %r = residual from initial guess
       !  if initial guess = 0, then we take %r=RHS
       ieptr=ie-nets+1
       do k=1,nlev
          iptr=1
          do j=1,np
             do i=1,np
                cg%state(ieptr)%r(iptr,k) = rhs(i,j,k,ie)
                iptr=iptr+1
             enddo
          enddo
       enddo
    enddo
    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Solver Loop 2
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do while (congrad(cg,red,maxits,tol))
       do ie=nets,nete
          ieptr=ie-nets+1
          do k=1,nlev
             ! apply preconditioner here: (note: r is not C0)
             cg%state(ieptr)%z(:,k) = cg%state(ieptr)%r(:,k)
             
             !reshape(cg%state(ieptr)%z(:,k),(/np,np/))
             iptr=1
             do j=1,np
                do i=1,np
                   x(i,j) = cg%state(ieptr)%z(iptr,k)
                   iptr=iptr+1
                enddo
             enddo
             
             ! DSS x to make it C0.  use LHS for storage:
             LHS(:,:,k,ie)=x(:,:)*elem(ie)%spheremp(:,:)
             
          end do
          call edgeVpack(edge1, LHS(1,1,1,ie), nlev, 0, elem(ie)%desc)
       end do
       call bndry_exchangeV(cg%hybrid,edge1)
       do ie=nets,nete
          ! unpack LHS
          call edgeVunpack(edge1, LHS(1,1,1,ie), nlev, 0, elem(ie)%desc)
          do k=1,nlev
             LHS(:,:,k,ie)=LHS(:,:,k,ie)*elem(ie)%rspheremp(:,:)
          enddo

          ieptr=ie-nets+1
          do k=1,nlev
             x(:,:)=LHS(:,:,k,ie) ! x() is now C0

             ! compute LHS(x) = x + laplace(x)
             LHS(:,:,k,ie)=x(:,:) + &
                  alambda*laplace_sphere_wk(x,deriv,elem(ie),var_coef=.false.)&
                  /elem(ie)%spheremp(:,:)

             iptr=1
             do j=1,np
                do i=1,np
                   cg%state(ieptr)%s(iptr,k) = LHS(i,j,k,ie)    ! new LHS, DG
                   cg%state(ieptr)%z(iptr,k) = x(i,j)           ! z must be C0
                   iptr=iptr+1
                end do
             end do
          enddo
       enddo
       
    end do  ! CG solver while loop
    !print *,'solver test CG iter = ',cg%iter


    ! ===============================
    ! Converged! compute actual error (not residual computed in solver)
    ! ===============================
    do ie=nets,nete
       ieptr=ie-nets+1
       do k=1,nlev
          iptr=1
          do j=1,np
             do i=1,np
                LHS(i,j,k,ie) = cg%state(ieptr)%x(iptr,k)
                iptr=iptr+1
             enddo
          enddo
       enddo
    enddo
    res = l2_snorm(elem,LHS,sol,hybrid,np,nets,nete) 
    if (hybrid%masterthread) print *,'normalized l2 error= ',res
    res = linf_snorm(LHS,sol,hybrid,np,nets,nete) 
    if (hybrid%masterthread) print *,'normalized linf error= ',res


  end subroutine solver_test



end module solver_mod
