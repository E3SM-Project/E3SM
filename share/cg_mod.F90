#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module cg_mod
  use kinds, only : real_kind
  use dimensions_mod, only : npsq, nlev, nelemd
  use global_norms_mod, only: wrap_repro_sum
  use parallel_mod, only: global_shared_buf, global_shared_sum
  use hybrid_mod, only : hybrid_t
  implicit none
  private

  integer, public, parameter :: CG_NO_DEBUG       = 0,  &
       CG_PR_FAST        = 1,  &
       CG_BIT_FOR_BIT    = 2,  &
       CG_PR_BIT_FOR_BIT = 3
  !CG-JMD
  type, public :: cg_state_t
     real(kind=real_kind) x(npsq,nlev)
     real(kind=real_kind) p(npsq,nlev)
     real(kind=real_kind) v(npsq,nlev)
     real(kind=real_kind) s(npsq,nlev)
     real(kind=real_kind) z(npsq,nlev)
     real(kind=real_kind) r(npsq,nlev)
  end type cg_state_t

  type, public :: cg_t

     integer len                                               ! length of domain
     integer ninst                                             ! number of instances
     integer nblks                                             ! number of blocks
     integer iter                                              ! internal iteration count
     integer ::   initalized,padding

     real (kind=real_kind), dimension(:,:), pointer :: wts           ! CG inner product weight

     type (cg_state_t), pointer  :: state(:)

     real (kind=real_kind), dimension(:),     pointer :: gammam1   ! gamma norm at iter-1
     real (kind=real_kind), dimension(:),     pointer :: sigmam1   ! gamma norm at iter-1
     real (kind=real_kind), dimension(:),     pointer :: rhs_norm  ! rhs norm <b,b>
     real (kind=real_kind), dimension(:),     pointer :: l2_norm   ! rhs norm <r_k,r_k>
     logical,               dimension(:),     pointer :: converged ! convergence mask
     integer debug_level                                       ! debug_level

     type (hybrid_t) :: hybrid
  end type cg_t

  public :: congrad
  public :: cg_create

  interface congrad
     module procedure congrad0
  end interface

contains

  ! ======================================================
  ! cg_create:
  !
  ! Create a conjugate gradient solver
  ! descriptor instantiation. This assumed to be 
  ! a thread private structure.
  !
  ! ninst        number of instances to solve
  ! len          problem vector length
  ! hybrid       general parallel descriptor
  ! debug_level  debugging level (optional)   0 = no debugging 
  !                                           1 = print fast residuals
  !                                           2 = use bit-for-bit reductions
  !                                           3 = print bit-for bit residuals
  !         
  ! =======================================================

  subroutine cg_create(cg,len,ninst,nblks,hybrid,debug_level,wts)
    use dimensions_mod, only : ne
    use hybrid_mod, only : hybrid_t
    use parallel_mod, only : abortmp
    type (cg_t)      , intent(inout)  :: cg
    integer          , intent(in)     :: len
    integer          , intent(in)     :: ninst
    integer          , intent(in)     :: nblks
    type (hybrid_t)  , intent(in)     :: hybrid
    integer, optional, intent(in)     :: debug_level
    real(kind=real_kind), optional, intent(in) :: wts(len,nblks)
    real(kind=real_kind) :: psum
    integer :: err, i
    cg%len         = len
    cg%ninst       = ninst
    cg%nblks       = nblks
    cg%iter        = 0

    if (0==ne) call abortmp('Error in cg_create: ne is zero')

    if (present(debug_level)) cg%debug_level=debug_level

    allocate(cg%converged(ninst),stat=err)
    if(err/=0) then
       print *, __FILE__,__LINE__,len,ninst,nblks,err
       call abortmp(' ')
    end if
    allocate(cg%gammam1(ninst),stat=err)
    if(err/=0) then
       print *, __FILE__,__LINE__,len,ninst,nblks,err
       call abortmp(' ')
    end if
    allocate(cg%sigmam1(ninst),stat=err)
    if(err/=0) then
       print *, __FILE__,__LINE__,len,ninst,nblks,err
       call abortmp(' ')
    end if
    allocate(cg%rhs_norm(ninst),stat=err)
    if(err/=0) then
       print *, __FILE__,__LINE__,len,ninst,nblks,err
       call abortmp(' ')
    end if
    allocate(cg%wts(len,nblks),stat=err)
    if(err/=0) then
       print *, __FILE__,__LINE__,len,ninst,nblks,err
       call abortmp(' ')
    end if

    allocate(cg%state(nblks),stat=err)

    if(err/=0) then
       print *, __FILE__,__LINE__,len,ninst,nblks,err
       call abortmp(' ')
    end if
    do i=1,nblks
       cg%state(i)%x=0.0_real_kind
       cg%state(i)%p=0.0_real_kind
       cg%state(i)%v=0.0_real_kind
       cg%state(i)%s=0.0_real_kind
       cg%state(i)%z=0.0_real_kind
       cg%state(i)%r=0.0_real_kind
    end do

    !CG-JMD

    if (present(wts)) then
       cg%wts(:,:) = wts(:,:)
    else
       cg%wts(:,:) = 1.0_real_kind
    endif

    cg%hybrid%par      = hybrid%par
    cg%hybrid%ithr     = hybrid%ithr
    cg%hybrid%NThreads = hybrid%NThreads
    cg%initalized      = 0

  end subroutine cg_create

  ! ======================================================
  ! congrad0:
  !
  ! preconditioned conjugate gradient solver core.
  ! operates by reverse communication. Assumes initial 
  ! guess (x) in zero.
  !
  ! =======================================================

  function congrad0(cg,red,maxits,tol) result(notdone)
    use reduction_mod, only :  reductionbuffer_ordered_1d_t
    type (cg_t)                       :: cg      ! private cg structure
    type (ReductionBuffer_ordered_1d_t)     :: red     ! shared data reduction data stucture
    integer, intent(in)               :: maxits  ! max iterations
    real (kind=real_kind), intent(in) :: tol     ! convergence tolerance

    logical                           :: notdone ! CG not done flag

    ! Local variables

    integer i,k,ie
    integer kptr

    real (kind=real_kind), dimension(cg%ninst)    :: eps
    real (kind=real_kind), dimension(3*cg%ninst)  :: redp

    real (kind=real_kind) :: gamma
    real (kind=real_kind) :: delta
    real (kind=real_kind) :: alpha
    real (kind=real_kind) :: beta
    real (kind=real_kind) :: sigma
    real (kind=real_kind) :: rhs_norm
    real (kind=real_kind) :: l2_norm


    if (cg%initalized .eq. 0) then

       cg%iter  = 0

       ! ============================
       ! Assumes that r=rhs and
       ! initial x1 = 0.0D0
       ! ============================

       cg%converged(:)=.false.

       notdone=.true.

       cg%initalized = 1
       return

    else if (cg%iter==0) then

       ! ==================================
       ! assumes that z_1 = (M^-1) * r_1
       ! assumes that v_1 = A * z_1
       ! ==================================

       do k=1,cg%ninst
          ! JPE This conditional should always be true - no need for it here
          !          if (.not.cg%converged(k)) then
          rhs_norm = 0.0D0
          gamma    = 0.0D0
          sigma    = 0.0D0

          do ie=1,cg%nblks

             global_shared_buf(ie,1:3) = 0.d0
             do i=1,cg%len
                cg%state(ie)%p(i,k)    = cg%state(ie)%z(i,k)                    ! p_1 = z_1
                cg%state(ie)%v(i,k)    = cg%state(ie)%s(i,k)                    ! p_1 = z_1
                global_shared_buf(ie,1)     = global_shared_buf(ie,1)     + cg%wts(i,ie)*cg%state(ie)%z(i,k)*cg%state(ie)%r(i,k) ! <z_1,r_1>
                global_shared_buf(ie,2)     = global_shared_buf(ie,2)     + cg%wts(i,ie)*cg%state(ie)%p(i,k)*cg%state(ie)%v(i,k) ! <p_1,v_1>
                global_shared_buf(ie,3)  = global_shared_buf(ie,3)  + cg%wts(i,ie)*cg%state(ie)%r(i,k)*cg%state(ie)%r(i,k) ! <b  ,b  >
             end do
          end do
          call wrap_repro_sum(nvars=3, comm=cg%hybrid%par%comm)
          
          red%buf(3*k-2,1) = global_shared_sum(1)
          red%buf(3*k-1,1) = global_shared_sum(2)
          red%buf(3*k  ,1) = global_shared_sum(3)

          !          end if
       end do

       kptr=0
       do k=1,cg%ninst
          ! JPE This conditional should always be true - no need for it here
          !          if (.not.cg%converged(k)) then
          gamma         = red%buf(3*kptr+1,1)     ! <z_1,r_1>
          sigma         = red%buf(3*kptr+2,1)     ! <p_1,v_1>
          cg%rhs_norm(k)= red%buf(3*kptr+3,1)     ! <b,b>
          cg%gammam1(k) = gamma
          cg%sigmam1(k) = sigma
#if 1
          eps(k)        = SQRT(cg%rhs_norm(k)/cg%rhs_norm(k))
#else
          eps(k) = ABS(gamma/cg%rhs_norm(k))
#endif             
          alpha = gamma/sigma
          do ie=1,cg%nblks
             do i=1,cg%len
                cg%state(ie)%x(i,k) =        alpha*cg%state(ie)%p(i,k)   ! x_2 = x_1 + alpha*p_1 , recall x1=0
                cg%state(ie)%r(i,k) = cg%state(ie)%r(i,k) - alpha*cg%state(ie)%v(i,k)   ! r_2 = r_1 - alpha*v_1
             end do
          end do
          kptr=kptr+1 
          !          end if
       end do

       if (cg%debug_level == 1 .or. cg%debug_level == 3) then
          if (cg%hybrid%par%masterproc .and. cg%hybrid%ithr == 0) then
             print *
             print *,"++++++++++++++++++++++++++++"
             do k=1,SIZE(eps)
                print *,"iter = ",cg%iter,"eps(",k,")=",eps(k)," rhs_norm(",k,")=",cg%rhs_norm(k)
             end do
          end if
       end if

       !CG-JMD cg%wrk1 => cg%r         ! assumes r_2 returns intact
       !CG-JMD cg%wrk2 => cg%z         ! request z_2 = (M^-1) * r_2
       !CG-JMD cg%wrk3 => cg%s         ! request s_2 = A*z_2

       notdone=.true.
       cg%iter=cg%iter+1

       return

    else !if (cg%iter>0) then

       ! ==================================
       ! assumes that z_k = (M^-1)*r_k
       ! and that     s_k =  A*z_k
       ! ==================================

       do k=1,cg%ninst
          if (.not.cg%converged(k)) then
             gamma=0.0D0
             delta=0.0D0
             l2_norm=0.0D0
             do ie=1,cg%nblks
                global_shared_buf(ie,1:3) = 0.d0
                do i=1,cg%len
                   global_shared_buf(ie,1)=global_shared_buf(ie,1) + cg%wts(i,ie)*cg%state(ie)%z(i,k)*cg%state(ie)%r(i,k)     ! gamma_k = <z_k,r_k>
                   global_shared_buf(ie,2)=global_shared_buf(ie,2) + cg%wts(i,ie)*cg%state(ie)%z(i,k)*cg%state(ie)%s(i,k)     ! delta_k = <z_k,s_k>
                   global_shared_buf(ie,3)=global_shared_buf(ie,3)  + cg%wts(i,ie)*cg%state(ie)%r(i,k)*cg%state(ie)%r(i,k)
                end do
             end do
             call wrap_repro_sum(nvars=3, comm=cg%hybrid%par%comm)
             red%buf(3*k-2,1) = global_shared_sum(1)
             red%buf(3*k-1,1) = global_shared_sum(2)
             red%buf(3*k  ,1) = global_shared_sum(3)
          end if
       end do

       kptr=0
       !
       !JPE: This initialization is required to assure that eps < tol when converged(k) is true.
       !
       eps=0.0 

       do k=1,cg%ninst
          if (.not.cg%converged(k)) then
             gamma  = red%buf(3*kptr+1,1)
             delta  = red%buf(3*kptr+2,1)
             l2_norm= red%buf(3*kptr+3,1)
#if 1
             eps(k) = SQRT(l2_norm/cg%rhs_norm(k))
#else
             eps(k) = ABS(gamma/cg%rhs_norm(k))
#endif             
             beta = gamma/cg%gammam1(k)
             sigma  = delta-beta*beta*cg%sigmam1(k)    ! eq(12) "method 1" of Lapack note 56
             alpha  = gamma/sigma
             cg%gammam1(k) = gamma
             cg%sigmam1(k) = sigma
             do ie=1,cg%nblks
                do i=1,cg%len
                   cg%state(ie)%p(i,k)=cg%state(ie)%z(i,k) +  beta*cg%state(ie)%p(i,k)
                   cg%state(ie)%v(i,k)=cg%state(ie)%s(i,k) +  beta*cg%state(ie)%v(i,k)
                   cg%state(ie)%x(i,k)=cg%state(ie)%x(i,k) + alpha*cg%state(ie)%p(i,k)
                   cg%state(ie)%r(i,k)=cg%state(ie)%r(i,k) - alpha*cg%state(ie)%v(i,k)
                end do
             end do
             if (eps(k) < tol) cg%converged(k)=.true.
             kptr=kptr+1
          end if
       end do

       if (cg%debug_level == 1 .or. cg%debug_level == 3) then
          if (cg%hybrid%par%masterproc .and. cg%hybrid%ithr == 0) then
             print *
             print *,"++++++++++++++++++++++++++++"
             do k=1,SIZE(eps)
                print *,"iter = ",cg%iter,"eps(",k,")=",eps(k)," rhs_norm(",k,")=",cg%rhs_norm(k)
             end do
          end if
       end if

       cg%iter=cg%iter+1
       if (cg%iter==maxits .or. ALL(eps(:) < tol)) then
	  !=================================================================
	  ! hopefully this will provide a slightly more reliable branch test 
	  !=================================================================
          cg%initalized=0
          notdone=.false.
       else
          notdone=.true.
       end if
       return               

    end if

10  format("cg:iter:",i4," residual=",e22.16," instance=",i4)

  end function congrad0

end module cg_mod



