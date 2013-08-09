#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module filter_mod
  use kinds, only : real_kind
  use dimensions_mod, only : np
  use perf_mod, only: t_startf, t_stopf, t_barrierf ! _EXTERNAL
  implicit none
  private

  real (kind=real_kind), parameter :: one=1.0D0, zero=0.0D0
  type, public :: filter_t
     real (kind=real_kind) :: FmatP(np,np)    ! Filter matrix
     character(len=8)      :: type    ! filter name: bv,fm, etc...
  end type filter_t

  public :: taylor_filter_create
  public :: fm_filter_create

  public :: filter_P
  public :: fm_transfer
  public :: bv_transfer
  public  :: bvsigma_test

  interface filter_matrix_create
     module procedure filter_matrix_create_inv
     module procedure filter_matrix_create_noinv
  end interface

  private :: bvsigma
  private :: filter_matrix_create
#ifdef _PRIM
  public :: preq_filter
  public :: prim_filter
#endif
contains
  ! ==================================
  ! taylor_filter_create:
  !
  ! Create a Taylor filter (requires
  ! communication) and use transfer
  ! function T...
  ! ==================================
  function taylor_filter_create(Tp, mu, gquadp) result(flt)
    use kinds, only : real_kind
    use quadrature_mod, only : quadrature_t, legendre, quad_norm
    real (kind=real_kind), dimension(:) :: Tp          ! transfer fn       
    real (kind=real_kind), intent(in)   :: mu             ! filter viscosity
    type (quadrature_t)  , intent(in)   :: gquadp  ! quadrature points and wts

    type (filter_t)                   :: flt

    ! Local variables 

    real (kind=real_kind), dimension(:), allocatable   :: gammaP ! normalization factor
    real (kind=real_kind), dimension(:,:), allocatable :: legP    ! legendre polynomials

    real (kind=real_kind), dimension(:,:), allocatable :: AP     ! spectral-physical 
    real (kind=real_kind), dimension(:,:), allocatable :: AinvP   ! physical-spectral


    integer                 :: npts
    integer                 :: n,j
    logical, parameter      :: Debug = .FALSE.

    call t_startf('taylor_filter_create')
    allocate(gammaP(np))

    allocate(legP(np,np))
    allocate(AP(np,np))

    allocate(AinvP(np,np))

    gammaP=quad_norm(gquadP,np)

    do j=1,np
       if(Debug) print *,'taylor_filter_create: gammaP(j) := ',gammaP(j)
       legP(:,j)=legendre(gquadP%points(j),np-1)
    end do

    ! -----------------------------------------
    ! compute T(k), the spectral coefficent 
    ! transfer function for the B-V filter
    ! -----------------------------------------

    AP = legP

    do n=1,np
       do j=1,np
          AinvP(j,n)=legP(n,j)*gquadP%weights(j)/gammaP(n)
       end do
    end do

    call filter_matrix_create(AP,Tp,AinvP,flt%FmatP,np,mu)

    flt%type = "bv"

    deallocate(gammaP)
    deallocate(legP)
    deallocate(AP)
    deallocate(AinvP)
    call t_stopf('taylor_filter_create')

  end function taylor_filter_create

  ! ==================================
  ! fm_filter_create:
  !
  ! Fischer-Mullen Filter Constructor
  ! with transfer function T...
  ! ==================================

  function fm_filter_create(Tp, mu, gquadP) result(flt)
    use quadrature_mod, only : quadrature_t, legendre

    real (kind=real_kind),   dimension(:) :: Tp    ! transfer fn
    real (kind=real_kind), intent(in) :: mu       ! filter viscosity
    type (quadrature_t)  , intent(in) :: gquadP    ! quadrature points/wts

    type (filter_t)                   :: flt

    ! Local variables 


    real (kind=real_kind), dimension(:,:), allocatable :: legP,legV    ! legendre polynomials

    real (kind=real_kind),   dimension(:,:), allocatable :: AP,AV      ! spectral-physical 

    integer                 :: k,n,j

    call t_startf('fm_filter_create')

    allocate(legP(np,np))

    allocate(AP(np,np))

    !JMD     allocate(flt%Fmat(npts,npts))
    !JMD     allocate(flt%tempt(npts,npts))

    do j=1,np
       legP(:,j)=legendre(gquadP%points(j),np-1)
       AP(1,j)=legP(1,j)
       AP(2,j)=legP(2,j)
       do n=3,np
          AP(n,j)=legP(n,j)-legP(n-2,j)
       end do
    end do


    call filter_matrix_create(AP,Tp,flt%FmatP,np,mu)

    flt%type = "fm"

    deallocate(legP)
    deallocate(AP)
    call t_stopf('fm_filter_create')

  end function fm_filter_create

  ! ==================================================
  ! This routing builds a Fischer-Mullen transfer 
  ! function that looks like:
  !
  !
  !     ^
  !  T  |
  !     |                 |
  !  1  |__________      _v_
  !     |          -_     
  !     |            \  wght
  !     |             \  ___
  !     |             |   ^
  !  0  |-------------|---|>
  !               ^   ^ 
  !     0         kc  npts  k-->
  !
  !  Where kc := npts-kcut is the point below which T = 1.
  ! ==================================================

  function fm_transfer(kcut,wght,npts) result(T)

    integer               :: kcut
    real (kind=real_kind) :: wght
    integer               :: npts
    real (kind=real_kind) :: T(npts)

    ! Local variables

    integer :: k
    integer :: kc
    real (kind=real_kind) :: cut
    real (kind=real_kind) :: amp

    call t_startf('fm_transfer')
    kc  = npts-kcut
    cut = kcut

    do k=1,kc
       T(k)=one
    end do

    do k=kc+1,npts
       amp = wght*(k-kc)*(k-kc)/(cut*cut)    ! quadratic growth
       T(k) = one - amp
    end do
    call t_stopf('fm_transfer')

  end function fm_transfer

  ! ===========================================
  ! bv_transfer:
  ! 
  ! compute Boyd-Vandeven transfer 
  ! function for a spectral element 
  ! with npts degrees of freedom...
  ! ===========================================

  function bv_transfer(p,s,npts) result(T)

    real (kind=real_kind), intent(in) :: p     ! order of sigmoid
    real (kind=real_kind), intent(in) :: s     ! scale of filter
    integer,               intent(in) :: npts  ! number of points
    real (kind=real_kind) :: T(npts)

    integer k
    real (kind=real_kind) :: arg
    real (kind=real_kind) :: rat

    call t_startf('bv_transfer')
    do k=1,npts
       rat = REAL(k,kind=real_kind)/REAL(npts,kind=real_kind)
       if (rat<s) then
          T(k)=one
       else
          arg = (rat-s)/(one-s)
          !
          !JPE: Need to assure that arg <= 1.0 
          !
          T(k)= bvsigma(p,min(arg,one))
       end if
    end do
    call t_stopf('bv_transfer')

  end function bv_transfer

  ! ===========================================
  ! bvsigma:
  !
  ! Boyd - Vandeven sigma function
  ! see p.98 of Taylor, Tribbia and Iskandarani
  ! ===========================================

  function bvsigma(p,x) result(sigma)

#ifdef CAM
    ! CAM code is required to use this function, to handle differences in
    ! Fortran 2008 support and compiler extensions.
    use shr_spfn_mod, only: erfc => shr_spfn_erfc
#endif

    real (kind=real_kind), intent(in) :: p
    real (kind=real_kind), intent(in) :: x
    real (kind=real_kind)             :: sigma

    ! Local variables

    real (kind=real_kind) :: om
    real (kind=real_kind) :: xfac
    real (kind=real_kind) :: arg
#if 0
    real (kind=real_kind) :: tmp
#endif

    ! For some old (pre-F2008) compilers, this is safer than "erfc", which
    ! is sometimes only single precision (not generic).
#ifndef CAM
    real*8  :: derfc
#endif
    call t_startf('bvsigma')

    om=ABS(x)-0.5D0

    if (x==zero) then
       sigma=one
    else if (x==one) then
       sigma=zero
    else
       xfac=one
       if (om /= zero) xfac = SQRT(-LOG((one-4.0D0*om**2))/(4.0D0*om**2))
       arg = 2.0D0*SQRT(p)*om*xfac
#if 0
       tmp = erfd(arg)
       print *, "x, erf", x, tmp
#endif
#ifdef CAM
       sigma = 0.50D0*erfc(arg)
#else
       sigma = 0.50D0*derfc(arg)
#endif
#if 0
       sigma = 0.50D0*(one - erfd(arg))
#endif
    end if
    call t_stopf('bvsigma')

  end function bvsigma

  ! ===========================================
  ! bvsigma_test:
  !
  ! Unit tester to reproduce figure 4, p 98 in 
  ! Taylor, Tribbia and Iskandarani
  ! ===========================================

  subroutine bvsigma_test(p)
    real (kind=real_kind), intent(in) :: p

    ! local variables

    integer, parameter    :: nvals = 20
    real (kind=real_kind) :: x,sigma
    integer               :: i

    write(6,10)p
10  format("Plotting Boyd-Vandeven transfer function of order ",f9.3)
    do i=1,nvals
       x=(real(i,kind=real_kind)-1.0_real_kind)/real(nvals,kind=real_kind)
       sigma=bvsigma(p,x)
       write(6,11)x,sigma
11     format(e22.15,3x,e22.15)
    end do
    write(6,*)""

  end subroutine bvsigma_test

  ! ===========================================
  ! filter:
  !
  ! Apply a Filter Matrix
  ! in both directions (x and y) of a 2-d domain.
  !
  ! ===========================================

  subroutine filter_P(p,flt) 
    use kinds, only : int_kind
    type (filter_t)                    :: flt
    real(kind=real_kind),intent(inout) :: p(np,np)

    ! Local

    integer  :: i,j,l
    integer(kind=int_kind), parameter :: unroll = 2
    real(kind=real_kind) :: temptp(np,np) 
    real(kind=real_kind) :: sumx00,sumx01
    real(kind=real_kind) :: sumx10,sumx11

    integer(kind=int_kind) :: rem,n2,npr,nprp1
    logical, parameter :: UseUnroll = .TRUE.

    !JMD if(MODULO(np,2) == 0 .and. UseUnroll) then 
    if(UseUnroll) then 
       rem = MODULO(np,unroll)
       npr  = np - rem 
       nprp1 = npr+1
       !JMD    print *,'Filter_P: np,npr,nprp1,rem:',np,npr,nprp1,rem

       do j=1,npr,2
          do l=1,npr,2
             sumx00=0.0d0
             sumx01=0.0d0
             sumx10=0.0d0
             sumx11=0.0d0

             do i=1,np
                sumx00 = sumx00 + flt%FmatP(i,l  )*p(i,j  )
                sumx01 = sumx01 + flt%FmatP(i,l+1)*p(i,j  )
                sumx10 = sumx10 + flt%FmatP(i,l  )*p(i,j+1)
                sumx11 = sumx11 + flt%FmatP(i,l+1)*p(i,j+1)
             end do

             temptP(j  ,l  ) = sumx00
             temptP(j  ,l+1) = sumx01
             temptP(j+1,l  ) = sumx10
             temptP(j+1,l+1) = sumx11

          end do
       end do
       if(rem>0) then 
          ! =======================
          ! Calculate the remainder
          ! =======================
          do l=1,npr
             sumx00=0.0d0
             sumx01=0.0d0
             do i=1,np
                sumx00 = sumx00 + flt%FmatP(i,l  )*p(i,np  )
                sumx01 = sumx01 + flt%FmatP(i,np  )*p(i,l  )
             enddo
             temptP(np  ,l  ) = sumx00	
             temptP(l  ,np  ) = sumx01	
          enddo
          sumx00=0.0d0
          do i=1,np
             sumx00 = sumx00 + flt%FmatP(i,np  )*p(i,np  )
          enddo
          temptP(np  ,np  ) = sumx00
       endif

       do j=1,npr,2
          do i=1,npr,2
             sumx00=0.0d0
             sumx01=0.0d0
             sumx10=0.0d0
             sumx11=0.0d0

             do l=1,np
                sumx00 = sumx00 +  flt%FmatP(l,j  )*temptP(l,i  )
                sumx01 = sumx01 +  flt%FmatP(l,j+1)*temptP(l,i  )
                sumx10 = sumx10 +  flt%FmatP(l,j  )*temptP(l,i+1)
                sumx11 = sumx11 +  flt%FmatP(l,j+1)*temptP(l,i+1)
             end do

             p(i  ,j  ) = sumx00
             p(i  ,j+1) = sumx01
             p(i+1,j  ) = sumx10
             p(i+1,j+1) = sumx11

          end do
       end do
       if(rem>0) then 
          ! =======================
          ! Calculate the remainder
          ! =======================
          do i=1,npr
             sumx00=0.0d0
             sumx01=0.0d0
             do l=1,np
                sumx00 = sumx00 +  flt%FmatP(l,i  )*temptP(l,np  )
                sumx01 = sumx01 +  flt%FmatP(l,np  )*temptP(l,i  )
             enddo
             p(np  ,i  ) = sumx00
             p(i  ,np  ) = sumx01
          enddo
          sumx00=0.0d0
          do l=1,np
             sumx00 = sumx00 +  flt%FmatP(l,np  )*temptP(l,np  )
          enddo
          p(np  ,np  ) = sumx00

       endif
    else
       do j=1,np
          do l=1,np
             sumx00=0.0d0
             do i=1,np
		sumx00 = sumx00 + flt%FmatP(i,l  )*p(i,j  )
             enddo
             temptP(j  ,l  ) = sumx00	
          enddo
       enddo
       do j=1,np
	  do i=1,np
	     sumx00=0.0d0
	     do l=1,np
		sumx00 = sumx00 +  flt%FmatP(l,j  )*temptP(l,i  )
	     enddo
	     p(i  ,j  ) = sumx00
          enddo
       enddo
    endif

  end subroutine filter_P
  

  ! =================================================================================
  !
  ! filter_matrix_create_inv:
  !
  ! This routing builds a 1D filter matrix, F, from a diagonal transfer function T 
  ! and the transform matrix A that maps a field u into spectral coefficients u_c,
  ! and its inverse A^-1, which is provided, via
  !    
  !   F = (1-mu)*I + mu*A.T.A^-1
  !
  !  where mu is the viscosity term, or weight of the filtering. 
  !
  ! ==================================================================================

  subroutine filter_matrix_create_inv(A,T,Ainv,F,npts,mu)

    integer                           :: npts      
    real (kind=real_kind), intent(in) :: A(npts,npts)
    real (kind=real_kind), intent(in) :: T(npts)
    real (kind=real_kind), intent(in) :: Ainv(npts,npts)
    real (kind=real_kind), intent(in) :: mu

    real (kind=real_kind) :: F(npts,npts)

    ! Local variables

    integer j,k,n
    real (kind=real_kind) :: fsum

    do k=1,npts  
       do j=1,npts
          fsum=zero
          do n=1,npts
             fsum = fsum + A(n,j)*T(n)*Ainv(k,n)
          end do
          if (j == k) then
             F(k,j) = (one-mu) + mu*fsum
          else
             F(k,j) = mu*fsum
          end if
       end do
    end do

  end subroutine filter_matrix_create_inv

  ! =================================================================================
  !
  ! filter_matrix_create_noinv:
  !
  ! This routing builds a 1D filter matrix, F, from a diagonal transfer function T 
  ! and the transform matrix A that maps a field u into spectral coefficients u_c,
  ! and its inverse A^-1, which is computed from A, via
  !    
  !   F = (1-mu)*I + mu*A.T.A^-1
  !
  !  where mu is the viscosity term, or weight of the filtering. 
  !
  ! =================================================================================

  subroutine filter_matrix_create_noinv(A,T,F,npts,mu)

    integer                           :: npts      
    real (kind=real_kind), intent(in) :: A(npts,npts)
    real (kind=real_kind), intent(in) :: T(npts)
    real (kind=real_kind), intent(in) :: mu

    real (kind=real_kind) :: F(npts,npts)

    ! Local variables

    integer               :: n,j,k
    integer               :: ipiv(npts)
    real (kind=real_kind) :: Ainv(npts,npts)
    real (kind=real_kind) :: fsum
    integer               :: ierr

    Ainv = A
    ierr=gaujordf(Ainv,npts,npts,ipiv)

    do k=1,npts  
       do j=1,npts
          fsum=zero
          do n=1,npts
             fsum = fsum + A(n,j)*T(n)*Ainv(k,n)
          end do
          if (j == k) then
             F(k,j) = (one-mu) + mu*fsum
          else
             F(k,j) = mu*fsum
          end if
       end do
    end do


  end subroutine filter_matrix_create_noinv

  ! ========================================================
  !
  ! Gauss-Jordan matrix inversion with full pivoting
  !
  ! Num. Rec. p. 30, 2nd Ed., Fortran
  !
  ! appropriated from code supplied by Paul Fischer and
  ! converted to F90 by RDL
  !
  ! a     is an m x n real matrix to be inverted in place
  ! ipiv  is the pivot vector
  !
  ! ========================================================

  function gaujordf(a,m,n,ipiv) result(ierr)

    integer, intent(in)   :: m,n
    real (kind=real_kind) :: a(m,n)
    integer               :: ipiv(n)
    integer               :: ierr

    integer :: indr(m),indc(n)
    integer :: i,j,k
    integer :: ir,jc

    real (kind=real_kind) :: amx,piv,tmp,work
    real (kind=real_kind) :: rmult(m)
    real (kind=real_kind) :: eps

    ierr = 0
    eps = 1.0D-14
    !
    do k=1,m
       ipiv(k)=0
    end do
    !
    do k=1,m
       amx=zero
       !
       !    Pivot search
       !
       do i=1,m
          if (ipiv(i) /= 1) then
             do j=1,m
                if (ipiv(j) == 0) then
                   if (abs(a(i,j)) >= amx) then
                      amx = ABS(a(i,j))
                      ir  = i
                      jc  = j
                   end if
                else if (ipiv(j) > 1) then
                   ierr = -ipiv(j)
                   return
                end if
             end do
          end if
       end do

       ipiv(jc) = ipiv(jc) + 1
       !
       !    Swap rows
       ! 
       if (ir /= jc) then
          do j=1,n
             tmp  = a(ir,j)
             a(ir,j) = a(jc,j)
             a(jc,j) = tmp
          end do
       end if

       indr(k)=ir
       indc(k)=jc

       !     write(6 ,*) k,' Piv:',jc,a(jc,jc)
       !       write(28,*) k,' Piv:',jc,a(jc,jc)

       if (abs(a(jc,jc)) < eps) then
          write(6 ,*) 'Gauss Jordan Pivot too small:',jc,a(jc,jc)
          !       write(28,*) 'small Gauss Jordan Piv:',jc,a(jc,jc)
          ierr = jc
          return
       end if

       piv = one/a(jc,jc)
       a(jc,jc)=one

       do j=1,n
          a(jc,j) = a(jc,j)*piv
       end do

       do j=1,n
          work    = a(jc,j)
          a(jc,j) = a(1 ,j)
          a(1 ,j) = work
       end do

       do i=2,m
          rmult(i) = a(i,jc)
          a(i,jc)  = zero
       end do

       do j=1,n
          do i=2,m
             a(i,j) = a(i,j) - rmult(i)*a(1,j)
          end do
       end do

       do j=1,n
          work    = a(jc,j)
          a(jc,j) = a(1 ,j)
          a(1 ,j) = work
       end do

    end do

    !
    !     Unscramble matrix
    !

    do j=m,1,-1
       if (indr(j) /= indc(j)) then
          do i=1,m
             tmp=a(i,indr(j))
             a(i,indr(j))=a(i,indc(j))
             a(i,indc(j))=tmp
          end do
       end if
    end do
  end function gaujordf

#ifdef _PRIM
  subroutine preq_filter(elem, edge3p1, flt, hybrid, nfilt, nets, nete)
    use element_mod, only : element_t
    use edge_mod, only  : Edgebuffer_t, edgevpack, edgevunpack
    use bndry_mod, only : bndry_exchangev
    use hybrid_mod, only  : hybrid_t
    use dimensions_mod, only : nlev
    type (element_t),    intent(inout)  :: elem(:)
    type (EdgeBuffer_t), intent(inout)  :: edge3p1
    type (filter_t),     intent(in)     :: flt
    integer,             intent(in)     :: nfilt
    integer,             intent(in)     :: nets,nete
    type (hybrid_t),     intent(in)     :: hybrid

    ! =====================
    ! Local variables
    ! =====================

    integer i,j,k,ie,kptr

    call t_barrierf('sync_preq_filter', hybrid%par%comm)
    call t_startf('preq_filter')
    if (flt%type == "bv") then
       do ie=nets,nete          
          elem(ie)%state%ps_v(:,:,nfilt) = EXP(elem(ie)%state%lnps(:,:,nfilt))
          call prim_filt_bv(elem(ie)%state%ps_v(:,:,nfilt),elem(ie)%state%T(:,:,:,nfilt), elem(ie)%state%v(:,:,:,:,nfilt), flt)

          do j=1,np
             do i=1,np
                elem(ie)%state%ps_v(i,j,nfilt) = elem(ie)%mp(i,j)*elem(ie)%state%ps_v(i,j,nfilt)
             end do
          end do
          do k=1,nlev
             do j=1,np
                do i=1,np
                   elem(ie)%state%T(i,j,k,nfilt)   = elem(ie)%mp(i,j)*elem(ie)%state%T(i,j,k,nfilt)
                   elem(ie)%state%v(i,j,1,k,nfilt) = elem(ie)%mp(i,j)*elem(ie)%state%v(i,j,1,k,nfilt)
                   elem(ie)%state%v(i,j,2,k,nfilt) = elem(ie)%mp(i,j)*elem(ie)%state%v(i,j,2,k,nfilt)
                end do
             end do
          end do

          kptr=0
          call edgeVpack(edge3p1,elem(ie)%state%v(:,:,1,1,nfilt),2*nlev,kptr,elem(ie)%desc)

          kptr=2*nlev
          call edgeVpack(edge3p1,elem(ie)%state%T(:,:,1,nfilt),nlev,kptr,elem(ie)%desc)

          kptr=3*nlev
          call edgeVpack(edge3p1,elem(ie)%state%ps_v(:,:,nfilt),1,kptr,elem(ie)%desc)

          !          kptr=0
          !          call edgerotate(edge3p1,2*nlev,kptr,elem(ie)%desc)

       end do

       call bndry_exchangeV(hybrid,edge3p1)

       do ie=nets,nete
          kptr=0
          call edgeVunpack(edge3p1, elem(ie)%state%v(:,:,1,1,nfilt), 2*nlev, kptr, elem(ie)%desc)

          kptr=2*nlev
          call edgeVunpack(edge3p1, elem(ie)%state%T(:,:,1,nfilt), nlev, kptr, elem(ie)%desc)

          kptr=3*nlev
          call edgeVunpack(edge3p1, elem(ie)%state%ps_v(:,:,nfilt), 1, kptr, elem(ie)%desc)

          do j=1,np
             do i=1,np
                elem(ie)%state%ps_v(i,j,nfilt) = elem(ie)%rmp(i,j)*elem(ie)%state%ps_v(i,j,nfilt)
             end do
          end do
          elem(ie)%state%lnps(:,:,nfilt) = LOG(elem(ie)%state%ps_v(:,:,nfilt))

          do k=1,nlev
             do j=1,np
                do i=1,np
                   elem(ie)%state%T(i,j,k,nfilt)   = elem(ie)%rmp(i,j)*elem(ie)%state%T(i,j,k,nfilt)
                   elem(ie)%state%v(i,j,1,k,nfilt) = elem(ie)%rmp(i,j)*elem(ie)%state%v(i,j,1,k,nfilt)
                   elem(ie)%state%v(i,j,2,k,nfilt) = elem(ie)%rmp(i,j)*elem(ie)%state%v(i,j,2,k,nfilt)
                end do
             end do
          end do

       end do
#ifdef DEBUGOMP
#if (defined HORIZ_OPENMP)
!$OMP BARRIER
#endif
#endif
    else if (flt%type == "fm") then
       do ie=nets,nete
          elem(ie)%state%ps_v(:,:,nfilt) = EXP(elem(ie)%state%lnps(:,:,nfilt))
          call prim_filt_fm(elem(ie)%state%ps_v(:,:,nfilt),          &
               elem(ie)%state%T(:,:,:,nfilt),   &
               elem(ie)%derived%div(:,:,:,nfilt), &
               elem(ie)%state%v(:,:,:,:,nfilt), &
               flt,                      &
               elem(ie)%D,                 &
               elem(ie)%Dinv)
          elem(ie)%state%lnps(:,:,nfilt) = LOG(elem(ie)%state%ps_v(:,:,nfilt))
       end do
    end if
    call t_stopf('preq_filter')

  end subroutine preq_filter

  ! =========================================
  ! Primitive Equations Boyd-Vandeven Filter
  ! =========================================

  subroutine prim_filt_bv(ps,T,v,flt)
    use dimensions_mod, only : nlev
    implicit none

    real(kind=real_kind), intent(inout) :: ps(np,np)
    real(kind=real_kind), intent(inout) :: T(np,np,nlev)
    real(kind=real_kind), intent(inout) :: v(np,np,2,nlev)

    type (filter_t), intent(in) :: flt

    integer :: k

    call filter_P(ps(:,:),flt)
    do k=1,nlev
       call filter_P(T(:,:,k),flt)
       call filter_P(v(:,:,1,k),flt)
       call filter_P(v(:,:,2,k),flt)
    end do

  end subroutine prim_filt_bv

  ! =========================================
  ! Primitive Equations Fisher-Mullen Filter
  ! =========================================

  subroutine prim_filt_fm(ps,T,div,v,flt,D,Dinv)
    use control_mod, only : integration
    use dimensions_mod, only : nlev
    implicit none

    real(kind=real_kind), intent(inout) :: ps(np,np)
    real(kind=real_kind), intent(inout) :: T(np,np,nlev)
    real(kind=real_kind), intent(inout) :: div(np,np,nlev)
    real(kind=real_kind), intent(inout) :: v(np,np,2,nlev)

    type (filter_t), intent(in)      :: flt
    real(kind=real_kind), intent(in) :: D(2,2,np,np)
    real(kind=real_kind), intent(in) :: Dinv(2,2,np,np)

    integer :: i,j,k
    real(kind=real_kind) :: v1,v2
    real(kind=real_kind) :: usph(np,np)
    real(kind=real_kind) :: vsph(np,np)

    call filter_P(ps(:,:),flt)
    do k=1,nlev
       call filter_P(T(:,:,k),flt)
       !     call filter_P(div(:,:,k),flt)

       ! ======================================
       ! Rotate cube velocities onto the sphere
       ! ======================================
       if( integration == "explicit" ) then
          do j=1,np
             do i=1,np
                v1 = v(i,j,1,k)
                v2 = v(i,j,2,k)
                usph(i,j)= v1*Dinv(1,1,i,j) + v2*Dinv(2,1,i,j)
                vsph(i,j)= v1*Dinv(1,2,i,j) + v2*Dinv(2,2,i,j)
             end do
          end do
       else
          do j=1,np
             do i=1,np
                v1 = v(i,j,1,k)
                v2 = v(i,j,2,k)
                usph(i,j)= v1*D(1,1,i,j) + v2*D(1,2,i,j)
                vsph(i,j)= v1*D(2,1,i,j) + v2*D(2,2,i,j)
             end do
          end do
       endif
       ! ======================================
       ! Filter the values on the sphere
       ! ======================================

       call filter_P(usph(:,:),flt)
       call filter_P(vsph(:,:),flt)

       ! ======================================
       ! Rotate sphere velocities onto the cube
       ! ======================================
       if( integration == "explicit" )  then
          do j=1,np
             do i=1,np
                v1 = usph(i,j)
                v2 = vsph(i,j)
                v(i,j,1,k)= v1*D(1,1,i,j) + v2*D(2,1,i,j)
                v(i,j,2,k)= v1*D(1,2,i,j) + v2*D(2,2,i,j)
             end do
          end do
       else
          do j=1,np
             do i=1,np
                v1 = usph(i,j)
                v2 = vsph(i,j)
                v(i,j,1,k)= v1*Dinv(1,1,i,j) + v2*Dinv(1,2,i,j)
                v(i,j,2,k)= v1*Dinv(2,1,i,j) + v2*Dinv(2,2,i,j)
             end do
          end do
       endif
    end do

  end subroutine prim_filt_fm

  subroutine prim_filter(elem, edge_3lp1,flt,nfilt,nets,nete,hybrid)
    use element_mod, only : element_t
    use parallel_mod, only : syncmp
    use dimensions_mod, only : nlev
    use edge_mod, only : EdgeBuffer_t, edgevpack, edgevunpack
    use bndry_mod, only : bndry_exchangev
    use hybrid_mod, only : hybrid_t
    use physical_constants, only : rearth,rrearth
    use control_mod, only : integration
    implicit none
    type (element_t), intent(inout), target :: elem(:)
    type (EdgeBuffer_t)           :: edge_3lp1
    type (filter_t) , intent(in)  :: flt
    integer         , intent(in)  :: nfilt
    integer         , intent(in)  :: nets,nete
    type (hybrid_t) , intent(in)  :: hybrid

    real (kind=real_kind) :: u(np,np)
    real (kind=real_kind) :: v(np,np)
    real (kind=real_kind) :: v1,v2
    real (kind=real_kind) :: tmp(np,np)

    integer i,j,k,ie
    integer kptr

    real (kind=real_kind), dimension(:,:), pointer :: mp,rmp

    call t_barrierf('sync_prim_filter', hybrid%par%comm)
    call t_startf('prim_filter')
 
    if (flt%type == "bv") then

#if (defined HORIZ_OPENMP)
       !$OMP BARRIER
#endif
       if (hybrid%ithr==0) call syncmp(hybrid%par)
#if (defined HORIZ_OPENMP)
       !$OMP BARRIER
#endif

       do ie=nets,nete

          mp => elem(ie)%mp

#ifdef _USE_VECTOR
          call vexp(elem(ie)%state%lnps(1,1,nfilt),elem(ie)%state%lnps(1,1,nfilt),np*np)
#else
          elem(ie)%state%lnps(:,:,nfilt)=EXP(elem(ie)%state%lnps(:,:,nfilt))
#endif
          call filter_P(elem(ie)%state%lnps(:,:,nfilt),flt)

          elem(ie)%state%lnps(:,:,nfilt) = mp*elem(ie)%state%lnps(:,:,nfilt)

#if (defined VERT_OPENMP)
!$omp parallel do private(k,i,j)
#endif
          do k=1,nlev
             call filter_P(elem(ie)%state%T(:,:,k,nfilt),flt)
             call filter_P(elem(ie)%state%v(:,:,1,k,nfilt),flt)
             call filter_P(elem(ie)%state%v(:,:,2,k,nfilt),flt)
             do j=1,np
                do i=1,np
                   elem(ie)%state%T(i,j,k,nfilt)   = mp(i,j)*elem(ie)%state%T(i,j,k,nfilt)
                   elem(ie)%state%v(i,j,1,k,nfilt) = mp(i,j)*elem(ie)%state%v(i,j,1,k,nfilt)
                   elem(ie)%state%v(i,j,2,k,nfilt) = mp(i,j)*elem(ie)%state%v(i,j,2,k,nfilt)
                end do
             end do
          end do

          kptr=0
          call edgeVpack(edge_3lp1,                 &
               elem(ie)%state%lnps(:,:,nfilt), &
               1,                         &
               kptr,                      &
               elem(ie)%desc)

          kptr=1
          call edgeVpack(edge_3lp1,                 &
               elem(ie)%state%T(:,:,1,nfilt),  &
               nlev,                      &
               kptr,                      &
               elem(ie)%desc)

          kptr=1+nlev
          call edgeVpack(edge_3lp1,                 &
               elem(ie)%state%v(:,:,1,1,nfilt),&
               2*nlev,                    &
               kptr,                      &
               elem(ie)%desc)


       end do

       call bndry_exchangeV(hybrid,edge_3lp1)

       do ie=nets,nete

          rmp => elem(ie)%rmp
          mp  => elem(ie)%mp

          kptr=0
          call edgeVunpack(edge_3lp1, elem(ie)%state%lnps(:,:,nfilt), 1, kptr, elem(ie)%desc)

          kptr=1
          call edgeVunpack(edge_3lp1, elem(ie)%state%T(:,:,1,nfilt), nlev, kptr, elem(ie)%desc)

          kptr=1+nlev
          call edgeVunpack(edge_3lp1, elem(ie)%state%v(:,:,1,1,nfilt), 2*nlev, kptr, elem(ie)%desc)

          elem(ie)%state%lnps(:,:,nfilt) = rmp*elem(ie)%state%lnps(:,:,nfilt)

#ifdef _USE_VECTOR
          call vlog(elem(ie)%state%lnps(1,1,nfilt),elem(ie)%state%lnps(1,1,nfilt),np*np)
#else
          elem(ie)%state%lnps(:,:,nfilt) = LOG(elem(ie)%state%lnps(:,:,nfilt))
#endif

#if (defined VERT_OPENMP)
!$omp parallel do private(k,i,j)
#endif
          do k=1,nlev
             do j=1,np
                do i=1,np
                   elem(ie)%state%T(i,j,k,nfilt)   = rmp(i,j)*elem(ie)%state%T(i,j,k,nfilt)
                   elem(ie)%state%v(i,j,1,k,nfilt) = rmp(i,j)*elem(ie)%state%v(i,j,1,k,nfilt)
                   elem(ie)%state%v(i,j,2,k,nfilt) = rmp(i,j)*elem(ie)%state%v(i,j,2,k,nfilt)
                end do
             end do
          end do

       end do
#ifdef DEBUGOMP
#if (defined HORIZ_OPENMP)
!$OMP BARRIER
#endif
#endif

    else if (flt%type == "fm") then

       do ie=nets,nete

#ifdef _USE_VECTOR
          call vexp(elem(ie)%state%lnps(1,1,nfilt),elem(ie)%state%lnps(1,1,nfilt),np*np)
#else
          elem(ie)%state%lnps(:,:,nfilt)=EXP(elem(ie)%state%lnps(:,:,nfilt))
#endif
          call filter_P(elem(ie)%state%lnps(:,:,nfilt),flt)
#ifdef _USE_VECTOR
          call vexp(elem(ie)%state%lnps(1,1,nfilt),elem(ie)%state%lnps(1,1,nfilt),np*np)
#else
          elem(ie)%state%lnps(:,:,nfilt) = LOG(elem(ie)%state%lnps(:,:,nfilt))
#endif
#if (defined VERT_OPENMP)
!$omp parallel do private(k,i,j,v1,v2,u,v)
#endif
          do k=1,nlev
             call filter_P(elem(ie)%state%T(:,:,k,nfilt),flt)
       if( integration == "explicit" ) then
                do j=1,np
                   do i=1,np
                      v1 = elem(ie)%state%v(i,j,1,k,nfilt)*rrearth
                      v2 = elem(ie)%state%v(i,j,2,k,nfilt)*rrearth
                      u(i,j)= v1*elem(ie)%Dinv(1,1,i,j) + v2*elem(ie)%Dinv(2,1,i,j)
                      v(i,j)= v1*elem(ie)%Dinv(1,2,i,j) + v2*elem(ie)%Dinv(2,2,i,j)
                   end do
                end do
             else
                do j=1,np
                   do i=1,np
                      v1 = elem(ie)%state%v(i,j,1,k,nfilt)*rrearth
                      v2 = elem(ie)%state%v(i,j,2,k,nfilt)*rrearth
                      u(i,j)= v1*elem(ie)%D(1,1,i,j) + v2*elem(ie)%D(1,2,i,j)
                      v(i,j)= v1*elem(ie)%D(2,1,i,j) + v2*elem(ie)%D(2,2,i,j)
                   end do
                end do
             endif

             call filter_P(u(:,:),flt)
             call filter_P(v(:,:),flt)

       if( integration == "explicit" ) then
                do j=1,np
                   do i=1,np
                      v1 = u(i,j)*rearth
                      v2 = v(i,j)*rearth
                      elem(ie)%state%v(i,j,1,k,nfilt)= v1*elem(ie)%D(1,1,i,j) + v2*elem(ie)%D(2,1,i,j)
                      elem(ie)%state%v(i,j,2,k,nfilt)= v1*elem(ie)%D(1,2,i,j) + v2*elem(ie)%D(2,2,i,j)
                   end do
                end do
             else
                do j=1,np
                   do i=1,np
                      v1 = u(i,j)*rearth
                      v2 = v(i,j)*rearth
                      elem(ie)%state%v(i,j,1,k,nfilt)= v1*elem(ie)%Dinv(1,1,i,j) + v2*elem(ie)%Dinv(1,2,i,j)
                      elem(ie)%state%v(i,j,2,k,nfilt)= v1*elem(ie)%Dinv(2,1,i,j) + v2*elem(ie)%Dinv(2,2,i,j)
                   end do
                end do
             endif
          end do

       end do

    end if
    call t_stopf('prim_filter')

  end subroutine prim_filter

#endif
end module filter_mod
