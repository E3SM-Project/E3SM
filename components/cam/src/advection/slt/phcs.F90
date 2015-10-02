
subroutine phcs(pmn     ,hmn     ,ix      ,x1)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute associated Legendre functions of the first kind of order m and
! degree n, and the associated derivatives for arg x1.

! Method: 
! Compute associated Legendre functions of the first kind of order m and
! degree n, and the associated derivatives for arg x1.  The associated
! Legendre functions are evaluated using relationships contained in
! "Tables of Normalized Associated Legendre Polynomials",
! S. L. Belousov (1962).  Both the functions and their derivatives are
! ordered in a linear stored rectangular array (with a large enough
! domain to contain the particular wavenumber truncation defined in the
! pspect common block) by column. m = 0->ptrm, and  n = m->ptrn + m 
!                m
! The functions P (x) are normalized such that 
!                n
!                          /   m     2
!                          | [P  (x)] dx = 1/2
!                          /   n
!                             __
! and must be multiplied by  |2  to match Belousov tables.
!                           \|
!                  m
! The derivatives H (x) are defined as 
!                  n        m           2    m
!                          H (x) = -(1-x ) dP (x)/dx
!                           n                n
!
! and are evaluated using the recurrence relationship
!                          _________________________
!      m          m       |  2   2                     m
!     H (x) = nx P (x) -  |(n - m )(2n + 1)/(2n - 1)  P   (x)
!      n          n      \|                            n-1
!
! Modified 1/23/97 by Jim Rosinski to use real*16 arithmetic in order to 
! achieve (nearly) identical values on all machines.
! 
! Author: CCM1
! 
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8, i8 => shr_kind_i8
  use pspect
  implicit none

#ifdef NO_R16
   integer,parameter :: r16= selected_real_kind(12) ! 8 byte real
#else
   integer,parameter :: r16= selected_real_kind(24) ! 16 byte real
#endif

!------------------------------Arguments--------------------------------
  integer , intent(in)  :: ix       ! Dimension of Legendre funct arrays
  real(r8), intent(in)  :: x1       ! sin of latitude, [sin(phi), or mu]
  real(r8), intent(out) ::  pmn(ix) ! Legendre function array
  real(r8), intent(out) ::  hmn(ix) ! Derivative array
!-----------------------------------------------------------------------

!---------------------------Local variables-----------------------------
  integer jmax       ! Loop limit (N+1=> 2D wavenumber limit +1)
  integer nmax       ! Large enough n to envelope truncation
  integer(i8) n      ! 2-D wavenumber index (up/down column)
  integer ml         ! intermediate scratch variable
  integer k          ! counter on terms in trig series expansion
  integer(i8) n2     ! 2*n
  integer m          ! zonal wavenumber index
  integer nto        ! intermediate scratch variable
  integer mto        ! intermediate scratch variable
  integer j          ! 2-D wavenumber index in recurrence evaluation
  integer nmaxm      ! loop limit in recurrence evaluation

  real(r16) xtemp(3,pmmax+ptrn+1)  ! Workspace for evaluating recurrence
!                                    ! relation where xtemp(m-2,n) and
!                                    ! xtemp(m-1,n) contain Pmn's required 
!                                    ! to evaluate xtemp(m,n) (i.e.,always
!                                    ! contains three adjacent columns of
!                                    ! the Pmn data structure)
!
  real(r16) xx1        ! x1 in extended precision
  real(r16) xte        ! cosine latitude [cos(phi)]
  real(r16) teta       ! pi/2 - latitute (colatitude)
  real(r16) an         ! coefficient on trig. series expansion
  real(r16) sinpar     ! accumulator in trig. series expansion 
  real(r16) cospar     ! accumulator in trig. series expansion 
  real(r16) p          ! 2-D wavenumber (series expansion)
  real(r16) q          ! intermediate variable in series expansion 
  real(r16) r          ! zonal wavenumber (recurrence evaluation)
  real(r16) p2         ! intermediate variable in series expansion 
  real(r16) rr         ! twice the zonal wavenumber (recurrence)
  real(r16) sqp        ! intermediate variable in series expansion 
  real(r16) cosfak     ! coef. on cos term in series expansion
  real(r16) sinfak     ! coef. on sin term in series expansion
  real(r16) ateta      ! intermediate variable in series expansion 
  real(r16) costet     ! cos term in trigonometric series expansion
  real(r16) sintet     ! sin term in trigonometric series expansion
!
  real(r16) t          ! intermediate variable (recurrence evaluation)
  real(r16) wm2        ! intermediate variable (recurrence evaluation)
  real(r16) wmq2       ! intermediate variable (recurrence evaluation)
  real(r16) w          ! intermediate variable (recurrence evaluation)
  real(r16) wq         ! intermediate variable (recurrence evaluation)
  real(r16) q2         ! intermediate variable (recurrence evaluation)
  real(r16) wt         ! intermediate variable (recurrence evaluation)
  real(r16) q2d        ! intermediate variable (recurrence evaluation)
  real(r16) cmn        ! cmn  recurrence coefficient (see Belousov)
  real(r16) xdmn       ! dmn  recurrence coefficient (see Belousov)
  real(r16) emn        ! emn  recurrence coefficient (see Belousov)
  real(r16) n2m1       ! n2 - 1 in extended precision
  real(r16) n2m3       ! n2 - 3 in extended precision
  real(r16) n2p1nnm1   ! (n2+1)*(n*n-1) in extended precision
  real(r16) twopmq     ! p + p - q in extended precision
!-----------------------------------------------------------------------
!
! Begin procedure by evaluating the first two columns of the Legendre
! function matrix (i.e., all n for m=0,1) via a trigonometric series 
! expansion (see eqs. 19 and 21 in Belousov, 1962).  Note that indexing
! is offset by one (e.g., m index for wavenumber m=0 is 1 and so on)
! Setup first ...
!
  xx1  = x1
  jmax = ptrn + 1
  nmax = pmmax + jmax
  xte = (1._r16-xx1*xx1)**0.5_r16
  teta = acos(xx1)
  an = 1._r16
  xtemp(1,1) = 0.5_r16    ! P00
!
! begin loop over n (2D wavenumber, or degree of associated Legendre 
! function) beginning with n=1 (i.e., P00 was assigned above)
! note n odd/even distinction yielding 2 results per n cycle
!
  do n=2,nmax
     sinpar = 0._r16
     cospar = 0._r16
     ml = n
     p = n - 1
     p2 = p*p
     sqp = 1._r16/(p2+p)**0.5_r16
     an = an*(1._r16 - 1._r16/(4._r16*p2))**0.5_r16
     cosfak = 1._r16
     sinfak = p*sqp
     do k=1,ml,2
        q = k - 1
        twopmq = p + p - q
        ateta = (p-q)*teta
        costet = cos(ateta)
        sintet = sin(ateta)
        if (n==k) costet = costet*0.5_r16
        if (k/=1) then
           cosfak = (q-1._r16)/q*(twopmq+2._r16)/(twopmq+1._r16)*cosfak
           sinfak = cosfak*(p-q)*sqp
        end if
        cospar = cospar + costet*cosfak
        sinpar = sinpar + sintet*sinfak
     end do
     xtemp(1,n)   = an*cospar      ! P0n vector
     xtemp(2,n-1) = an*sinpar      ! P1n vector
  end do
!
! Assign Legendre functions and evaluate derivatives for all n and m=0,1
!
  pmn(1) = 0.5_r16
  pmn(1+jmax) = xtemp(2,1)
  hmn(1) = 0._r16
  hmn(1+jmax) = xx1*xtemp(2,1)
  do n=2,jmax
     pmn(n) = xtemp(1,n)
     pmn(n+jmax) = xtemp(2,n)
     n2 = n + n
     n2m1 = n2 - 1
     n2m3 = n2 - 3
     n2p1nnm1 = (n2+1)*(n*n-1)
     hmn(n) = (n-1)*(xx1*xtemp(1,n)-(n2m1/n2m3)**0.5_r16*xtemp(1,n-1))
     hmn(n+jmax) = n*xx1*xtemp(2,n)-(n2p1nnm1/n2m1)**0.5_r16*xtemp(2,n-1)
  end do
!
! Evaluate recurrence relationship for remaining Legendre functions
! (i.e., m=2 ... PTRM) and associated derivatives (see eq 17, Belousov)
!
  do m=3,pmmax
     r = m - 1
     rr = r + r
     xtemp(3,1) = (1._r16+1._r16/rr)**0.5_r16*xte*xtemp(2,1)
     nto = (m-1)*jmax
     pmn(nto+1) = xtemp(3,1)
     hmn(nto+1) = r*xx1*xtemp(3,1)
     nmaxm = nmax - m
!
! Loop over 2-D wavenumber (i.e., degree of Legendre function)
! Pmn's and Hmn's for current zonal wavenumber, r
!
     do j=2,nmaxm
        mto = nto + j
        t = j - 1
        q = rr + t - 1
        wm2 = q + t
        w = wm2 + 2
        wq = w*q
        q2 = q*q - 1
        wmq2 = wm2*q2
        wt = w*t
        q2d = q2 + q2
        cmn = ((wq*(q-2._r16))/(wmq2-q2d))**0.5_r16
        xdmn = ((wq*(t+1._r16))/wmq2)**0.5_r16
        emn = (wt/((q+1._r16)*wm2))**0.5_r16
        xtemp(3,j) = cmn*xtemp(1,j) - xx1*(xdmn*xtemp(1,j+1)-emn*xtemp(3,j-1))
        pmn(mto) = xtemp(3,j)
        hmn(mto) = (r+t)*xx1*xtemp(3,j) - (wt*(q+1._r16)/wm2)**0.5_r16*xtemp(3,j-1)
     end do
!
! shift Pmn's to left in workspace (setup for next recurrence pass)
!
!++pjr
! not initialized above
     xtemp(2,nmax) = 0._r16
     do j=nmaxm,nmax
        xtemp(3,j) = 0._r16
     end do
!--pjr
     do n=1,nmax
        xtemp(1,n) = xtemp(2,n)
        xtemp(2,n) = xtemp(3,n)
     end do
  end do

  return
end subroutine phcs

