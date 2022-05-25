!-----------------------------------------------------------------------
! $Id$
!===============================================================================
module lapack_interfaces
    
  ! This module acts as an interface to Lapack routines. It's main purpose
  ! is to make interfaces available to clubb that can handle both 
  ! single and doulbe precision. This may be compiled along with Lapack
  ! source code, or along with a linked Lapack library such as MKL. 
    
  implicit none
  
  public :: lapack_gbsv, lapack_gbsvx, &
            lapack_gtsv, lapack_gtsvx, &
            lapack_isnan, lapack_potrf, &
            lapack_poequ, lapack_laqsy, &
            lapack_syev, lapack_trmv
            
  private :: &
    dgbsv_wrap, sgbsv_wrap, &
    dgbsvx_wrap, sgbsvx_wrap, &
    dgtsv_wrap, sgtsv_wrap, &
    dgtsvx_wrap, sgtsvx_wrap, &
    disnan_wrap, sisnan_wrap, &
    dpotrf_wrap, spotrf_wrap, &
    dpoequ_wrap, spoequ_wrap, &
    dlaqsy_wrap, slaqsy_wrap, &
    dsyev_wrap, ssyev_wrap, &
    dtrmv_wrap, strmv_wrap
    
  ! Interface for Lapack general band solver, single or double precision
  interface lapack_gbsv
      module procedure dgbsv_wrap
      module procedure sgbsv_wrap
  end interface

  ! Interface for Lapack general band solver, expert version, single or double precision
  interface lapack_gbsvx
      module procedure dgbsvx_wrap
      module procedure sgbsvx_wrap
  end interface
  
  ! Interface for Lapack tridiagonal matrix solver, single or double precision
  interface lapack_gtsv
      module procedure dgtsv_wrap
      module procedure sgtsv_wrap
  end interface
  
  ! Interface for Lapack tridiagonal matrix solver, expert version, single or double precision
  interface lapack_gtsvx
      module procedure dgtsvx_wrap
      module procedure sgtsvx_wrap
  end interface
  
  ! Interface for Lapack nan check, single or double precision
  interface lapack_isnan
      module procedure disnan_wrap
      module procedure sisnan_wrap
  end interface
  
  ! Interface for Lapack's Cholesky factorization of a real symmetric positive definite
  ! matrix, single or double precision
  interface lapack_potrf
      module procedure dpotrf_wrap
      module procedure spotrf_wrap
  end interface
  
  ! Interface for Lapack routine to compute row and column scalings intended to 
  ! equilibriate a symmetric positive definite matrix, single or doulbe precision
  interface lapack_poequ
      module procedure dpoequ_wrap
      module procedure spoequ_wrap
  end interface
  
  ! Interface for Lapack routine to equilibriate a symmetric matrix, single or double precision
  interface lapack_laqsy
      module procedure dlaqsy_wrap
      module procedure slaqsy_wrap
  end interface
  
  ! Interface for Lapack routine to compute all eigenvalues and, optionally, eigenvectors
  ! of a real symmetric matrix, single or double precision
  interface lapack_syev
      module procedure dsyev_wrap
      module procedure ssyev_wrap
  end interface
  
  ! Interface for Lapack routines to performe one of the following matrix-vector operations
  !     x := A*x,   or   x := A**T*x,
  ! where A is an upper or lower triangular matrix, single or double precision
  interface lapack_trmv
      module procedure dtrmv_wrap
      module procedure strmv_wrap
  end interface
      
      
  private ! Set Default Scope
  
  contains
      
  ! ==================== General Band Solver Wrappers ====================
  
  ! Double precision wrapper
  subroutine dgbsv_wrap( n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info )
      
      implicit none
      
      external :: dgbsv
      
      integer            info, kl, ku, ldab, ldb, n, nrhs
      integer            ipiv( * )
      double precision   ab( ldab, * ), b( ldb, * )
      
      call dgbsv( n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info )
      
  end subroutine dgbsv_wrap
  
  ! Single precision wrapper
  subroutine sgbsv_wrap( n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info )
      
      implicit none
      
      external :: sgbsv
      
      integer   info, kl, ku, ldab, ldb, n, nrhs
      integer   ipiv( * )
      real      ab( ldab, * ), b( ldb, * )
      
      call sgbsv( n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info )
      
  end subroutine sgbsv_wrap
  
  
  ! ==================== Band Solver Expert Wrappers ====================

  ! Double precision wrapper
  subroutine dgbsvx_wrap( fact, trans, n, kl, ku, nrhs, ab, ldab, afb, &
                          ldafb, ipiv, equed, r, c, b, ldb, x, ldx, &
                          rcond, ferr, berr, work, iwork, info )
    implicit none
    
    external :: dgbsvx
 
    character          equed, fact, trans
    integer            info, kl, ku, ldab, ldafb, ldb, ldx, n, nrhs
    double precision   rcond
    integer            ipiv( * ), iwork( * )
    double precision   ab( ldab, * ), afb( ldafb, * ), b( ldb, * ), &
                       berr( * ), c( * ), ferr( * ), r( * ), &
                       work( * ), x( ldx, * )
                       
    call dgbsvx( fact, trans, n, kl, ku, nrhs, ab, ldab, afb, &
                 ldafb, ipiv, equed, r, c, b, ldb, x, ldx, &
                 rcond, ferr, berr, work, iwork, info )
                 
  end subroutine dgbsvx_wrap
  
  ! Single precision wrapper
  subroutine sgbsvx_wrap( fact, trans, n, kl, ku, nrhs, ab, ldab, afb, &
                          ldafb, ipiv, equed, r, c, b, ldb, x, ldx, &
                          rcond, ferr, berr, work, iwork, info )                  
    implicit none
    
    external :: sgbsvx
 
    character   equed, fact, trans
    integer     info, kl, ku, ldab, ldafb, ldb, ldx, n, nrhs
    real        rcond
    integer     ipiv( * ), iwork( * )
    real        ab( ldab, * ), afb( ldafb, * ), b( ldb, * ), &
                berr( * ), c( * ), ferr( * ), r( * ), &
                work( * ), x( ldx, * )
                       
    call sgbsvx( fact, trans, n, kl, ku, nrhs, ab, ldab, afb, &
                 ldafb, ipiv, equed, r, c, b, ldb, x, ldx, &
                 rcond, ferr, berr, work, iwork, info )
                 
  end subroutine sgbsvx_wrap
  
  
  ! ==================== Tridiagonal Solver Wrappers ====================

  ! Double precision wrapper
  subroutine dgtsv_wrap( n, nrhs, dl, d, du, b, ldb, info )
      
    implicit none
    
    external :: dgtsv
      
    integer            info, ldb, n, nrhs
    double precision   b( ldb, * ), d( * ), dl( * ), du( * )
    
    call dgtsv( n, nrhs, dl, d, du, b, ldb, info )
    
  end subroutine dgtsv_wrap
    
  ! Single precision wrapper
  subroutine sgtsv_wrap( n, nrhs, dl, d, du, b, ldb, info )
      
    implicit none
    
    external :: sgtsv
      
    integer     info, ldb, n, nrhs
    real        b( ldb, * ), d( * ), dl( * ), du( * )
    
    call sgtsv( n, nrhs, dl, d, du, b, ldb, info )
    
  end subroutine sgtsv_wrap
  
  
  ! ==================== Tridiagonal Solver Expert Wrappers ====================

  ! Double precision wrapper
  subroutine dgtsvx_wrap( fact, trans, n, nrhs, dl, d, du, dlf, df, duf, &
                          du2, ipiv, b, ldb, x, ldx, rcond, ferr, berr, &
                          work, iwork, info )
    implicit none
    
    external :: dgtsvx
  
    character          fact, trans
    integer            info, ldb, ldx, n, nrhs
    double precision   rcond
    integer            ipiv( * ), iwork( * )
    double precision   b( ldb, * ), berr( * ), d( * ), df( * ), &
                       dl( * ), dlf( * ), du( * ), du2( * ), duf( * ), &
                       ferr( * ), work( * ), x( ldx, * )
                       
    call dgtsvx( fact, trans, n, nrhs, dl, d, du, dlf, df, duf, &
                 du2, ipiv, b, ldb, x, ldx, rcond, ferr, berr, &
                 work, iwork, info )
                 
  end subroutine dgtsvx_wrap
  
  ! Single precision wrapper
  subroutine sgtsvx_wrap( fact, trans, n, nrhs, dl, d, du, dlf, df, duf, &
                          du2, ipiv, b, ldb, x, ldx, rcond, ferr, berr, &
                          work, iwork, info )
    implicit none
    
    external :: sgtsvx
 
    character   fact, trans
    integer     info, ldb, ldx, n, nrhs
    real        rcond
    integer     ipiv( * ), iwork( * )
    real        b( ldb, * ), berr( * ), d( * ), df( * ), &
                dl( * ), dlf( * ), du( * ), du2( * ), duf( * ), &
                ferr( * ), work( * ), x( ldx, * )
                       
    call sgtsvx( fact, trans, n, nrhs, dl, d, du, dlf, df, duf, &
                 du2, ipiv, b, ldb, x, ldx, rcond, ferr, berr, &
                 work, iwork, info )
                 
  end subroutine sgtsvx_wrap
  
  
  ! ==================== NaN Check Wrappers ====================

  ! Double precision wrapper
  !-----------------------------------------------------------------------
  logical function disnan_wrap( ndim, nrhs, variable )

  ! Description:
  !   Check for NaN values in a variable using the LAPACK subroutines

  ! References:
  !   <http://www.netlib.org/lapack/single/sisnan.f>
  !   <http://www.netlib.org/lapack/double/disnan.f>
  !-----------------------------------------------------------------------
  
    implicit none
      
#ifdef NO_LAPACK_ISNAN /* Used for older LAPACK libraries that don't have sisnan/disnan */

    integer, intent(in) :: &
      ndim, & ! Size of variable
      nrhs    ! Number of right hand sides

    double precision, dimension(ndim,nrhs), intent(in) :: &
      variable ! Variable to check

      disnan_wrap = any( variable(:,1:nrhs) /= variable(:,1:nrhs) )
#else
  
    logical, external :: &
      disnan    ! Procedure

    integer, intent(in) :: &
      ndim, & ! Size of variable
      nrhs    ! Number of right hand sides

    double precision, dimension(ndim,nrhs), intent(in) :: &
      variable ! Variable to check

    integer :: k, j

    ! ---- Begin Code ----

    disnan_wrap = .false.

    do k = 1, ndim
      do j = 1, nrhs
            
        ! Lapack NaN check function, sisnan for single precision or disnan for double precision
        disnan_wrap = disnan( variable(k,j) )
          
        if ( disnan_wrap ) exit
      end do
      if ( disnan_wrap ) exit
    end do
      
#endif /* NO_LAPACK_ISNAN */

    return
  end function disnan_wrap
  
  ! Single precision wrapper
  !-----------------------------------------------------------------------
  logical function sisnan_wrap( ndim, nrhs, variable )

  ! Description:
  !   Check for NaN values in a variable using the LAPACK subroutines

  ! References:
  !   <http://www.netlib.org/lapack/single/sisnan.f>
  !   <http://www.netlib.org/lapack/double/disnan.f>
  !-----------------------------------------------------------------------

    implicit none
      
#ifdef NO_LAPACK_ISNAN /* Used for older LAPACK libraries that don't have sisnan/disnan */

    integer, intent(in) :: &
      ndim, & ! Size of variable
      nrhs    ! Number of right hand sides

    real, dimension(ndim,nrhs), intent(in) :: &
      variable ! Variable to check

      sisnan_wrap = any( variable(:,1:nrhs) /= variable(:,1:nrhs) )
#else
  
    logical, external :: &
      sisnan    ! Procedure

    integer, intent(in) :: &
      ndim, & ! Size of variable
      nrhs    ! Number of right hand sides

    real, dimension(ndim,nrhs), intent(in) :: &
      variable ! Variable to check

    integer :: k, j

    ! ---- Begin Code ----

    sisnan_wrap = .false.

    do k = 1, ndim
      do j = 1, nrhs
            
        ! Lapack NaN check function, sisnan for single precision or disnan for double precision
        sisnan_wrap = sisnan( variable(k,j) )
          
        if ( sisnan_wrap ) exit
      end do
      if ( sisnan_wrap ) exit
    end do
      
#endif /* NO_LAPACK_ISNAN */

    return
  end function sisnan_wrap
 
 
  ! ==================== Cholesky Factorization Wrappers ====================

  ! Double precision wrapper
  subroutine dpotrf_wrap( uplo, n, a, lda, info )
      
    implicit none
    
    external :: dpotrf
      
    character          uplo
    integer            info, lda, n
    double precision   a( lda, * )

    call dpotrf( uplo, n, a, lda, info )
    
  end subroutine dpotrf_wrap
  
  ! Single precision wrapper
  subroutine spotrf_wrap( uplo, n, a, lda, info )
      
    implicit none
    
    external :: spotrf
      
    character   uplo
    integer     info, lda, n
    real        a( lda, * )

    call spotrf( uplo, n, a, lda, info )
    
  end subroutine spotrf_wrap
  
  
  ! ==================== Equilibrium Scaling Calculation Wrappers ====================

  ! Double precision wrapper
  subroutine dpoequ_wrap( n, a, lda, s, scond, amax, info )
      
    implicit none
    
    external :: dpoequ
      
    integer   info, lda, n
    double precision      amax, scond
    double precision      a( lda, * ), s( * )

    call dpoequ( n, a, lda, s, scond, amax, info )
  
  end subroutine dpoequ_wrap
  
  ! Single precision wrapper
  subroutine spoequ_wrap( n, a, lda, s, scond, amax, info )
      
    implicit none
    
    external :: spoequ
      
    integer   info, lda, n
    real      amax, scond
    real      a( lda, * ), s( * )

    call spoequ( n, a, lda, s, scond, amax, info )
  
  end subroutine spoequ_wrap
  
  
  ! ==================== Matrix Equilibriator Wrappers ====================

  ! Double precision wrapper
  subroutine dlaqsy_wrap( uplo, n, a, lda, s, scond, amax, equed )
      
    implicit none
    
    external :: dlaqsy
      
    character          equed, uplo
    integer            lda, n
    double precision   amax, scond
    double precision   a( lda, * ), s( * )
  
    call dlaqsy( uplo, n, a, lda, s, scond, amax, equed )
  
  end subroutine dlaqsy_wrap
    
  ! Single precision wrapper
  subroutine slaqsy_wrap( uplo, n, a, lda, s, scond, amax, equed )
      
    implicit none
    
    external :: slaqsy
      
    character   equed, uplo
    integer     lda, n
    real        amax, scond
    real        a( lda, * ), s( * )
  
    call slaqsy( uplo, n, a, lda, s, scond, amax, equed )
    
  end subroutine slaqsy_wrap
  
  
  ! ==================== Eigenvalue/vector Calculation Wrappers ====================

  ! Double precision wrapper
  subroutine dsyev_wrap( jobz, uplo, n, a, lda, w, work, lwork, info )
      
    implicit none
    
    external :: dsyev
      
    character          jobz, uplo
    integer            info, lda, lwork, n
    double precision   a( lda, * ), w( * ), work( * )
    
    call dsyev( jobz, uplo, n, a, lda, w, work, lwork, info )
      
  end subroutine dsyev_wrap
  
  ! Single precision wrapper
  subroutine ssyev_wrap( jobz, uplo, n, a, lda, w, work, lwork, info )
     
    implicit none
    
    external :: ssyev
      
    character   jobz, uplo
    integer     info, lda, lwork, n
    real        a( lda, * ), w( * ), work( * )
    
    call ssyev( jobz, uplo, n, a, lda, w, work, lwork, info )
      
  end subroutine ssyev_wrap
  
  
  ! ==================== Matrix Operations Wrappers ====================
  
  ! Double precision wrapper
  subroutine dtrmv_wrap( uplo, trans, diag, n, a, lda, x, incx)
     
    implicit none
    
    external :: dtrmv
      
    integer             incx,lda,n
    character           diag,trans,uplo             
    double precision    a(lda,*),x(*)

    call dtrmv( uplo, trans, diag, n, a, lda, x, incx)
    
  end subroutine dtrmv_wrap
  
  ! Single precision wrapper
  subroutine strmv_wrap( uplo, trans, diag, n, a, lda, x, incx)
     
    implicit none
    
    external :: strmv
      
    integer     incx,lda,n
    character   diag,trans,uplo             
    real        a(lda,*),x(*)

    call strmv( uplo, trans, diag, n, a, lda, x, incx)
    
  end subroutine strmv_wrap
          
end module lapack_interfaces