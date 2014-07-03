subroutine settau(zdt)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Set time invariant hydrostatic matrices, which depend on the reference
! temperature and pressure in the semi-implicit time step. Note that
! this subroutine is actually called twice, because the effective time
! step changes between step 0 and step 1.
! 
! Method: 
! zdt = delta t for next semi-implicit time step.
! 
! Author: CCM1
! 
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  use pspect
  use commap
  use physconst,    only: cappa, rair, gravit
  use abortutils,   only: endrun
  use spmd_utils,   only: masterproc
  use hycoef,       only : hypi, hybi, hypd
  use sgexx,        only: dgeco, dgedi
  use cam_logfile,  only: iulog

  implicit none


!------------------------------Arguments--------------------------------
  real(r8), intent(in) :: zdt      ! time step (or dt/2 at time 0)
!---------------------------Local workspace-----------------------------
  real(r8) aq(plev,plev)
  real(r8) rcond,z(plev),det(2),work(plev)
  integer ipvt(plev)
  real(r8) zcr(plev)             ! gravity wave equivalent depth
  real(r8) zci(plev)             ! dummy, used to print phase speeds
  real(r8) zdt2                  ! zdt**2
  real(r8) factor                ! intermediate workspace
  real(r8) zdt0u                 ! vertical diff. of ref. temp (above)
  real(r8) zshu                  ! interface "sigma" (above)
  real(r8) zr2ds                 ! 1./(2.*hypd(k))
  real(r8) zdt0d                 ! vertical diff. of ref. temp (below)
  real(r8) zshd                  ! interface "sigma" (below)
  real(r8) ztd                   ! temporary accumulator
  real(r8) zcn                   ! sq(n)
  real(r8) zb(plev,plev)         ! semi-implicit matrix in d equation
  real(r8), save :: zdt_init=0   ! reinitialize if zdt <> zdt_init

  integer k,kk,kkk               ! level indices
  integer n                      ! n-wavenumber index
  integer nneg                   ! number of unstable mean temperatures
!-----------------------------------------------------------------------
!
  if (zdt == zdt_init) return

! save dt for which this code has performed the initialization
  zdt_init=zdt
	 	
  zdt2 = zdt*zdt
!
! Set mean temperature
! NOTE: Making t0 an actual function of height ***DOES NOT WORK***
!
  do k=1,plev
     t0(k) = 300._r8
  end do
!
! Calculate hydrostatic matrix tau
!
  zdt0u = 0._r8
  zshu = 0._r8
  do k=1,plev
     zr2ds = 1._r8/(2._r8*hypd(k))
     if (k < plev) then
        zdt0d = t0(k+1) - t0(k)
        zshd = hybi(k+1)
     else
        zdt0d = 0._r8
        zshd = 0._r8
     end if

     factor = ((zdt0u*zshu + zdt0d*zshd) - (zdt0d + zdt0u))*zr2ds
     do kk=1,k-1
        tau(kk,k) = factor*hypd(kk) + cappa*t0(k)*ecref(kk,k)
     end do

     factor = (zdt0u*zshu + zdt0d*zshd - zdt0d)*zr2ds
     tau(k,k) = factor*hypd(k) + cappa*t0(k)*ecref(k,k)

     factor = (zdt0u*zshu + zdt0d*zshd)*zr2ds
     do kk=k+1,plev
        tau(kk,k) = factor*hypd(kk)
     end do
     zdt0u = zdt0d
     zshu = zshd
  end do
!
! Vector for linear surface pressure term in divergence
! Pressure gradient and diagonal term of hydrostatic components
!
  do k=1,plev
     bps(k) = t0(k)
     bps(k) = bps(k)*rair
  end do
  do k=1,plev
     do kk=1,plev
        ztd = bps(k) * hypd(kk)/hypi(plevp)
        do kkk=1,plev
           ztd = ztd + href(kkk,k)*tau(kk,kkk)
        end do
        zb(kk,k) = ztd
        aq(kk,k) = ztd
     end do
  end do
!
! Compute and print gravity wave equivalent depths and phase speeds
!
  call qreig(zb      ,plev    ,zcr     )

  do k=1,plev
     zci(k) = sign(1._r8,zcr(k))*sqrt(abs(zcr(k)))
     zcr(k) = zcr(k) / gravit
  end do

  if (masterproc) then
     write(iulog,910) (t0(k),k=1,plev)
     write(iulog,920) (zci(k),k=1,plev)
     write(iulog,930) (zcr(k),k=1,plev)
  end if
!
! Test for unstable mean temperatures (negative phase speed and eqivalent
! depth) for at least one gravity wave.
!
  nneg = 0
  do k=1,plev
     if (zcr(k)<=0._r8) nneg = nneg + 1
  end do

  if (nneg/=0) then
     call endrun ('SETTAU: UNSTABLE MEAN TEMPERATURE.')
  end if
!
! Compute and invert matrix a(n)=(i+sq*b*delt**2)
!          
  do k=1,plev
     do kk=1,plev
        aq(kk,k) = aq(kk,k)*zdt2
        bm1(kk,k,1) = 0._r8
     end do
  end do
  do n=2,pnmax
     zcn = sq(n)
     do k=1,plev
        do kk=1,plev
           zb(kk,k) = zcn*aq(kk,k)
           if(kk.eq.k) zb(kk,k) = zb(kk,k) + 1._r8
        end do
     end do
!
! Use linpack routines to invert matrix
!  
     call dgeco(zb,plev,plev,ipvt,rcond,z)
     call dgedi(zb,plev,plev,ipvt,det,work,01)
     do k=1,plev
        do kk=1,plev
           bm1(kk,k,n) = zb(kk,k)
        end do
     end do
  end do

910 format(' REFERENCE TEMPERATURES FOR SEMI-IMPLICIT SCHEME = ',  /(1x,12f9.3))
920 format(' GRAVITY WAVE PHASE SPEEDS (M/S) FOR MEAN STATE = '    /(1x,12f9.3))
930 format(' GRAVITY WAVE EQUIVALENT DEPTHS (M) FOR MEAN STATE = ' /(1x,12f9.3))

  return
end subroutine settau

!============================================================================================

subroutine qreig(a       ,i       ,b       )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Create complex matrix P with real part = A and imaginary part = 0
! Find its eigenvalues and return their real parts.
! 
! Method: 
! This routine is of unknown lineage. It is only used to provide the
! equivalent depths of the reference atmosphere for a diagnostic print
! in SETTAU and has no effect on the model simulation. Therefore it can 
! be replaced at any time with a functionally equivalent, but more 
! understandable, procedure.  Consequently, the internal commenting has 
! not been brought up to CAM standards.                               
! 
! Author: 
! 
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  implicit none

!------------------------------Arguments--------------------------------
  real(r8), intent(in)  :: a(*)      ! Input real part
  integer , intent(in)  :: i
  real(r8), intent(out) :: b(*)
!-----------------------------------------------------------------------

!---------------------------Local variables-----------------------------
  complex(r8) p(plev*plev)    
  complex(r8) q(plev*plev)
  integer l,ij,ik               ! indicies
!-----------------------------------------------------------------------
!
!  l = 0
!  do ij=1,i
!     do ik=1,i
!        l = l + 1
!        p(l) = cmplx(a(l),0._r8,r8)
!     end do
!  end do

  do l = 1, i*i
     p(l) = cmplx( a(l), 0.0_r8, r8)
  end do

  call cmphes(p       ,i       ,1       ,i       )
  call cmplr(p       ,q       ,i)

  do ij=1,i
     b(ij) = real(q(ij),r8)
  end do

  return
end subroutine qreig

!============================================================================================

subroutine cmphes(ac      ,nac     ,k       ,l       )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Reduce complex matrix (ac) to upper Hessenburg matrix (ac)
! 
! Method: 
! This routine is of unknown lineage. It is only used to provide the
! equivalent depths of the reference atmosphere for a diagnostic print
! in SETTAU and has no effect on the model simulation. Therefore it can 
! be replaced at any time with a functionally equivalent, but more 
! understandable, procedure.  Consequently, the internal commenting has 
! not been brought up to CCM3 or CAM standards.                               
!
! Author: 
! 
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none

!------------------------------Arguments--------------------------------
  integer, intent(in) :: nac                ! Dimension of one side of matrix ac
  integer, intent(in) :: k,l                !
  complex(r8), intent(inout) :: ac(nac,nac) ! On input, complex matrix to be converted
                                            ! On output, upper Hessenburg matrix
!-----------------------------------------------------------------------

!---------------------------Local variables-----------------------------
  complex(r8) x        
  complex(r8) y         
  integer la            
  integer m1            
  integer i,m,j         ! Indices
  integer j1,i1         ! Loop limits
!-----------------------------------------------------------------------
!
  la = l - 1
  m1 = k + 1
  do m=m1,la
     i = m
     x = (0.0_r8,0.0_r8)
     do j=m,l
        if (abs(ac(j,m-1))>abs(x)) then
           x = ac(j,m-1)
           i = j
        end if
     end do
     if (i/=m) then
        j1 = m - 1
        do j=j1,nac
           y = ac(i,j)
           ac(i,j) = ac(m,j)
           ac(m,j) = y
        end do
        do j=1,l
           y = ac(j,i)
           ac(j,i) = ac(j,m)
           ac(j,m) = y
        end do
     end if
     if (x/=(0.0_r8,0.0_r8)) then
        i1 = m + 1
        do i=i1,l
           y = ac(i,m-1)
           if (y/=(0.0_r8,0.0_r8)) then
              y = y/x
              ac(i,m-1) = y
              do j=m,nac
                 ac(i,j) = ac(i,j) - y*ac(m,j)
              end do
              do j=1,l
                 ac(j,m) = ac(j,m) + y*ac(j,i)
              end do
           end if
        end do
     end if
  end do

  return
end subroutine cmphes

!============================================================================================

subroutine cmplr(hes     ,w       ,nc)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute w, eigenvalues of upper Hessenburg matrix hes
! 
! Method: 
! This routine is of unknown lineage. It is only used to provide the
! equivalent depths of the reference atmosphere for a diagnostic print
! in SETTAU and has no effect on the model simulation. Therefore it can 
! be replaced at any time with a functionally equivalent, but more 
! understandable, procedure.  Consequently, the internal commenting has 
! not been brought up to CCM3 or CAM standards.                               
!
! Author: 
! 
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use abortutils,   only: endrun
  use cam_logfile,  only: iulog

  implicit none

!------------------------------Arguments--------------------------------
  integer    , intent(in) :: nc            ! Dimension of input and output matrices
  complex(r8), intent(inout) :: hes(nc,nc) ! Upper hessenberg matrix from comhes
  complex(r8), intent(out):: w(nc)         ! Weights
!-----------------------------------------------------------------------

!---------------------------Local variables-----------------------------
  integer itest 
  integer nfail      ! Limit for number of iterations to convergence
  integer ntest 
  integer n,j,m 
  integer i          ! Eigenvalue
  integer its        ! Iteration counter
  integer l 
  integer l1,m1,n1,i1
  real(r8) a 
  real(r8) sr 
  real(r8) si 
  real(r8) tr 
  real(r8) ti 
  real(r8) xr 
  real(r8) yr 
  real(r8) zr 
  real(r8) xi 
  real(r8) yi 
  real(r8) areal 
  real(r8) eps
  complex(r8) s 
  complex(r8) t 
  complex(r8) x 
  complex(r8) y 
  complex(r8) z 
  complex(r8) u

  data itest/0/
  save a,eps,sr,itest
!-----------------------------------------------------------------------
!
  nfail = 30
  if (itest==0) then
     a = 1
5    continue
     eps = a
     sr = 1 + a
     a = a/2.0_r8
     if (sr/=1.0_r8) go to 5
     itest = 1
  end if
  if (nc.le.0) then
     write(iulog,*)'CMPLR: Entered with incorrect dimension  '
     write(iulog,*)'NC=',NC
     call endrun
  end if
  ntest = 10
  n = nc
  t = 0.0_r8
10 continue
  if (n==0) go to 300
  its = 0
20 continue
  if (n/=1) then
     do l1=2,n
        l = n + 2 - l1
        if (abs(hes(l,l-1)) <= eps*(abs(hes(l-1,l-1))+abs(hes(l,l)))) go to 50
     end do
  end if
  l = 1
50 continue
  if (l/=n) then
     if (its==nfail) then
        i = nc - n + 1
        write(iulog,*)'CMPLR: Failed to converge in ',nfail,' iterations'
        write(iulog,*)'Eigenvalue=',i
        call endrun
     end if
     if (its==ntest) then
        ntest = ntest + 10
        sr = hes(n,n-1)
        si = hes(n-1,n-2)
        sr = abs(sr)+abs(si)
        u = (0.0_r8,-1.0_r8)*hes(n,n-1)
        tr = u
        u = (0.0_r8,-1.0_r8)*hes(n-1,n-2)
        ti = u
        tr = abs(tr) + abs(ti)
        s = cmplx(sr,tr)
     else
        s = hes(n,n)
        x = hes(n-1,n)*hes(n,n-1)
        if (abs(x)/=0.0_r8) then
           y = 0.5_r8*(hes(n-1,n-1)-s)
           u = y*y + x
           z = sqrt(u)
           u = conjg(z)*y
           areal = u
           if (areal<0.0_r8) z = -z
           x = x/(y+z)
           s = s - x
        end if
     end if
     do i=1,n
        hes(i,i) = hes(i,i) - s
     end do
     t = t + s
     its = its + 1
     j = l + 1
     xr = abs(hes(n-1,n-1))
     yr = abs(hes(n,n-1))
     zr = abs(hes(n,n))
     n1 = n - 1
     if ((n1/=1).and.(n1>=j)) then
        do m1=j,n1
           m = n1 + j - m1
           yi = yr
           yr = abs(hes(m,m-1))
           xi = zr
           zr = xr
           xr = abs(hes(m-1,m-1))
           if (yr.le.eps*zr/yi*(zr+xr+xi)) go to 100
        end do
     end if
     m = l
100  continue
     m1 = m + 1
     do i=m1,n
        x = hes(i-1,i-1)
        y = hes(i,i-1)
        if (abs(x)<abs(y)) then
           i1 = i - 1
           do j=i1,n
              z = hes(i-1,j)
              hes(i-1,j) = hes(i,j)
              hes(i,j) = z
           end do
           z = x/y
           w(i) = 1.0_r8
        else
           z = y/x
           w(i) = -1.0_r8
        end if
        hes(i,i-1) = z
        do j=i,n
           hes(i,j) = hes(i,j) - z*hes(i-1,j)
        end do
     end do
     m1 = m + 1
     do j=m1,n
        x = hes(j,j-1)
        hes(j,j-1) = 0.0_r8
        areal = w(j)
        if (areal>0.0_r8) then
           do i=l,j
              z = hes(i,j-1)
              hes(i,j-1) = hes(i,j)
              hes(i,j) = z
           end do
        end if
        do i=l,j
           hes(i,j-1) = hes(i,j-1) + x*hes(i,j)
        end do
     end do
     go to 20
  end if
  w(n) = hes(n,n) + t
  n = n - 1
  go to 10
300 continue

  return
end subroutine cmplr

