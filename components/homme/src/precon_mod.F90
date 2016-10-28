#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#undef _DGEMV
module precon_mod
! a preconditioner version of solver_mod to eventually be incorporated into
! solver_mod perhaps

  use kinds, only : real_kind, int_kind
  use dimensions_mod, only : npsq, nlev
  use perf_mod, only : t_startf, t_stopf
  use solver_mod, only : blkjac_t
  implicit none
  private

  character(len=8), private, parameter :: blkjac_storage = "inverse"
  !  character(len=8), private, parameter :: blkjac_storage = "LUfactor"

  public  :: pcg_presolver_nonstag

contains

  ! ================================================
  ! pcg_presolver_nonstag:
  !
  ! Preconditioned conjugate gradient solver on the
  ! Gauss-Lobatto nonstaggered grid
  ! 
  ! ================================================

  function pcg_presolver_nonstag(pptr, rhs) result(x) 
    use dimensions_mod, only : nlev, np, npsq, nelemd
    use element_mod, only : element_t
    use precon_type_mod, only : precon_type
    use reduction_mod, only : reductionbuffer_ordered_1d_t
    use cg_mod, only : cg_t, congrad
    use edge_mod, only : edgevpack, edgerotate, edgevunpack
    use edgetype_mod, only : edgebuffer_t
    use derivative_mod, only : derivative_t, gradient_wk, gradient, divergence
    use control_mod, only : maxits, tol, precon_method
    use physical_constants, only : rrearth
    use bndry_mod, only : bndry_exchangeV
    use linear_algebra_mod, only : matvec
    use solver_mod, only : blkjac_t

    integer  :: nets,nete

    type(element_t) :: elem(nelemd)

    real (kind=real_kind), intent(in) :: rhs(np,np,nlev,nelemd) ! right hand side of operator
    type (cg_t)                       :: cg             ! conjugate gradient    (private)
    type (ReductionBuffer_ordered_1d_t)  :: red         ! CG reduction buffer   (shared memory)
    type (EdgeBuffer_t)               :: edge1   ! Laplacian gradient edge buffer (shared memory)
    type (EdgeBuffer_t)               :: edge2   ! Laplacian gradient edge buffer (shared memory)
    real (kind=real_kind)             :: lambdasq(nlev) ! Helmholtz lengthscale (private)
    type (derivative_t)               :: deriv   ! non staggered derivative struct (private)
    type (blkjac_t)                   :: blkjac(nelemd)

    real (kind=real_kind)             :: x(np,np,nlev,nelemd)     ! solution (result)
    type(precon_type) ,pointer        :: pptr

    ! ===========
    ! Local
    ! ===========

    real (kind=real_kind), dimension(np,np,2,2) :: metinv
    real (kind=real_kind), dimension(np,np) :: metdet
    real (kind=real_kind), dimension(np,np) :: rmetdet
    real (kind=real_kind), dimension(np,np) :: rmp
    real (kind=real_kind), dimension(np,np) :: mp

    real (kind=real_kind) :: gradp(np,np,2,nlev,nelemd)
    real (kind=real_kind) :: div(np,np,nlev,nelemd)
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

! declarations of unpacked pptr pointer
    cg              =  pptr%cg 
    red             =  pptr%red 
    edge1           =  pptr%edge1 
    edge2           =  pptr%edge2 
   do k=1,nlev
    lambdasq(k)     =  pptr%lambdasq(k)
   end do 
    deriv           =  pptr%deriv
    nets            =  pptr%nets
    nete            =  pptr%nete
   do ie=nets,nete
    elem(ie)        =  pptr%base(ie)
   end do 
    blkjac      =  pptr%bc

    ! ========================================
    ! Load the rhs in the cg struct...
    ! ========================================

    call t_startf('pcg_presolver')

    do ie=nets,nete
       ieptr=ie-nets+1
       do k=1,nlev
          iptr=1
          do j=1,np
             do i=1,np
                pptr%cg%state(ieptr)%r(iptr,k) = rhs(i,j,k,ie)
                iptr=iptr+1
             end do
          end do
       end do
    end do

    do while (congrad(cg,red,maxits,tol))
       call t_startf('in_congrad')

       do ie=nets,nete
          ieptr=ie-nets+1
          metinv = elem(ie)%metinv
          metdet = elem(ie)%metdet
          rmetdet  = elem(ie)%rmetdet
          rmp     = elem(ie)%rmp
          mp      = elem(ie)%mp

          do k=1,nlev
             if (.not.cg%converged(k)) then

                ! ========================================
                ! Apply preconditioner: wrk2 = (M^-1) * wrk1 
                ! ========================================

                if (precon_method == "block_jacobi") then

                   if (blkjac_storage == "LUfactor") then
                      call dgesl(cg%state(ieptr)%r(:,k),cg%state(ieptr)%z(:,k),blkjac(ie)%E(:,:,k),blkjac(ie)%ipvt(:,k),npsq)
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
                      gradp(i,j,1,k,ie) = metdet(i,j)*(metinv(i,j,1,1)*gradp1 + &
                           metinv(i,j,1,2)*gradp2)
                      gradp(i,j,2,k,ie) = metdet(i,j)*(metinv(i,j,2,1)*gradp1 + &
                           metinv(i,j,2,2)*gradp2)
#else
                      gradp(i,j,1,k,ie) = (metinv(i,j,1,1)*gradp1 + metinv(i,j,1,2)*gradp2)
                      gradp(i,j,2,k,ie) = (metinv(i,j,2,1)*gradp1 + metinv(i,j,2,2)*gradp2)
#endif
                   end do
                end do
             end if
          end do

          kptr=0
          call edgeVpack(edge2,gradp(1,1,1,1,ie),2*nlev,kptr,ie)

          kptr=0
          call edgerotate(edge2,2*nlev,kptr,elem(ie)%desc)

       end do

       !$OMP BARRIER

       call bndry_exchangeV(cg%hybrid,edge2)

       !$OMP BARRIER

       do ie=nets,nete
          ieptr=ie-nets+1

          kptr=0
          call edgeVunpack(edge2, gradp(1,1,1,1,ie), 2*nlev, kptr, ie)

          do k=1,nlev
             if (.not.cg%converged(k)) then

                ! ====================
                ! 2*np*np Flops
                ! ====================

                do j=1,np
                   do i=1,np
                      gradp(i,j,1,k,ie) = rmp(i,j)*gradp(i,j,1,k,ie)
                      gradp(i,j,2,k,ie) = rmp(i,j)*gradp(i,j,2,k,ie)
                   ! gradp(i,j,1,k,ie) = metdet(i,j)*rmp(i,j)*gradp(i,j,1,k,ie)
                   ! gradp(i,j,2,k,ie) = metdet(i,j)*rmp(i,j)*gradp(i,j,2,k,ie)
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

          k=0
          call edgeVpack(edge1, div(1,1,1,ie), nlev, kptr, ie)

       end do


       !$OMP BARRIER

       call bndry_exchangeV(cg%hybrid,edge1)

       !$OMP BARRIER

       ! ==========================================
       ! compute Helmholtz operator, store in wrk3
       !  4*np*np Flops
       ! ==========================================

       do ie=nets,nete
          ieptr=ie-nets+1

          rmp      = elem(ie)%rmp
          rmetdet  = elem(ie)%rmetdet
          metdet   = elem(ie)%metdet

          kptr=0
          call edgeVunpack(edge1, div(1,1,1,ie), nlev, kptr, ie)

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
       call t_stopf('in_congrad')

    end do  ! CG solver while loop
    ! ===============================
    ! Converged! unpack wrk3 into x
    ! ===============================

    do ie=nets,nete
       ieptr=ie-nets+1
       do k=1,nlev
          iptr=1
          do j=1,np
             do i=1,np
!                !       x(i,j,k,ie)=cg%wrk3(iptr,k,ieptr)
                x(i,j,k,ie)=cg%state(ieptr)%x(iptr,k)
                iptr=iptr+1
             end do
          end do
       end do
    end do
    call t_stopf('pcg_presolver')

  end function pcg_presolver_nonstag

end module precon_mod
