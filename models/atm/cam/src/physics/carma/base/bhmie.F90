!! See Bohren and Huffman, "Absorption and Scattering of Light by
!! Small Particles", 1983, pg 480 (in Appendix A). 
!!
!! Subroutine bhmie calculates amplitude scattering matrix
!! elements and efficiencies for extinction, total scattering
!! and backscattering for a given size parameter and
!! relative refractive index.
!!
!! From the main program:
!!   refrel = cmplx(refre,refim) / refmed
!!
!! @author Chuck Bardeen
!! @version 2011
subroutine bhmie(carma, x, refrel, nang, s1, s2, Qext, Qsca, Qback, gfac, rc)

	! types
  use carma_precision_mod
  use carma_enums_mod, only     : RC_ERROR
  use carma_types_mod, only     : carma_type
  use carma_mod

  implicit none

  type(carma_type), intent(in)         :: carma         !! the carma object
  real(kind=f), intent(in)             :: x             !! radius / wavelength
  complex(kind=f), intent(in)          :: refrel        !! refractive index particle / reference refractive index
  integer, intent(in)                  :: nang          !! number of angles in s1 and s2
  complex(kind=f), intent(out)         :: s1(2*nang-1)  !! CORE RADIUS
  complex(kind=f), intent(out)         :: s2(2*nang-1)  !! REAL PART OF THE CORE INDEX OF REFRACTION
  real(kind=f), intent(out)            :: Qext          !! EFFICIENCY FACTOR FOR EXTINCTION
  real(kind=f), intent(out)            :: Qsca          !! EFFICIENCY FACTOR FOR SCATTERING
  real(kind=f), intent(out)            :: Qback         !! BACK SCATTER CROSS SECTION.
  real(kind=f), intent(out)            :: gfac          !! asymmetry factor
  integer, intent(inout)               :: rc            !! return code, negative indicates failure
  

  real(kind=f)    :: amu(100), theta(100), pi(100), tau(100), pi0(100), pi1(100)
  complex(kind=f) :: y, xi, xi0, xi1, an, bn
  complex(kind=f), allocatable ::   d(:)
  complex(kind=f) :: ccan, ccbn, anmi1, bnmi1
  real(kind=f)    :: psi0, psi1, psi, dn, dx, chi0, chi1, apsi0, apsi1, g1, g2
  real(kind=f)    :: dang, fn, ffn, apsi, chi, p, t, xstop, ymod
  integer         :: j, jj, n, nn, rn, nmx, nstop
      

  ! Mie x and y values.
  dx = x
  y  = x * refrel

  ! Series terminated after nstop terms
  xstop = x + 4._f * x**0.3333_f + 2.0_f
  nstop = xstop

  ! Will loop over nang angles.
  ymod = int(abs(y))
  nmx  = max(xstop, ymod) + 15
  dang = 1.570796327_f / float(nang - 1)
  allocate(d(nmx))
  
  do j = 1, nang
    theta(j) = (float(j) - 1._f) * dang
    amu(j)   = cos(theta(j))
  end do
  
  ! Logarithmic derivative d(j) calculated by downword
  ! recurrence beginning with initial value 0.0 + i*0.0
  ! at j = nmx
  d(nmx) = cmplx(0.0_f, 0.0_f)
  nn     = nmx-1
!  write(*,*) 'nmx=',nmx,' d(nmx)=',d(nmx), ' nn=',nn
      
  do n = 1, nn
    rn = nmx - n + 1
    d(nmx-n) = (rn/y) - (1._f / (d(nmx - n + 1) + rn / y))
    
!    write(*,*) 'n=',n,' rn=',rn,' y=', y,' d(nmx-n)=',d(nmx-n)
!    write(*,*) 'rn/y=',rn/y, 'd(nmx-n+1)=',d(nmx-n+1),'(d(nmx-n+1)+rn/y)',(d(nmx-n+1)+rn/y),'1./(d(nmx-n+1)+rn/y)',1./(d(nmx-n+1)+rn/y)
  end do
  
  pi0(1:nang) = 0.0_f
  pi1(1:nang) = 1.0_f
  
  nn = 2 * nang-1
  s1(1:nn) = cmplx(0.0_f, 0.0_f)
  s2(1:nn) = cmplx(0.0_f, 0.0_f)
  
  ! Riccati-Bessel functions with real argument x
  ! calculated by upward recurrence
  psi0  = cos(dx)
  psi1  = sin(dx)
  chi0  = -sin(x)
  chi1  = cos(x)
  apsi0 = psi0
  apsi1 = psi1
  xi0   = cmplx(apsi0,-chi0)
  xi1   = cmplx(apsi1,-chi1)
  Qsca  = 0.0_f
  g1    = 0.0_f
  g2    = 0.0_f
  n     = 1
  
  ! Loop over the terms n in the Mie series
  do while (.true.)
    dn   = n
    rn   = n
    fn   = (2._f * rn + 1._f) / (rn * (rn + 1._f))
    ffn  = (rn - 1._f) * (rn + 1._f) / rn
    psi  = (2._f * dn - 1._f) * psi1 / dx - psi0
    apsi = psi
    chi  = (2._f * rn - 1._f) * chi1 / x - chi0
    xi   = cmplx(apsi, -chi)
!    write(*,*) 'n=', n
!    write(*,*) 'd(n)=',d(n),' refrel=',refrel,' rn=',rn, ' x=',x,'apsi=',apsi,' apsi1=',apsi1

    an   = (d(n) / refrel + rn / x) * apsi - apsi1
!    write(*,*) 'an=',an,' xi=',xi,' xi1=',xi1

    an   = an / ((d(n) / refrel + rn / x) * xi - xi1)
    bn   = (refrel * d(n) + rn / x) * apsi - apsi1
    bn   = bn / ((refrel * d(n) + rn / x) * xi - xi1)
    ccan = conjg(an)
    ccbn = conjg(bn)
    g2   = g2 + fn * real(an * ccbn)
      
    if (n-1 > 0) then
      g1   = g1 + ffn * real(anmi1 * ccan + bnmi1 * ccbn)
    end if
    Qsca = Qsca + (2._f * rn + 1._f) * (abs(an) * abs(an) + abs(bn) * abs(bn))
      
    do j = 1, nang
      jj     = 2 * nang-j
      pi(j)  = pi1(j)
      tau(j) = rn * amu(j) * pi(j) - (rn + 1._f) * pi0(j)
      p      = (-1._f)**(n-1)
!      write(*,*) 'fn=',fn,' an=',an,' bn=',bn,' pi(j)=',pi(j),' tau(j)=',tau(j)
  
      s1(j)  = s1(j) + fn * (an * pi(j) + bn * tau(j))
      t      = (-1._f)**n
      s2(j)  = s2(j) + fn * (an * tau(j) + bn * pi(j))
      
      if (j.ne.jj) then
        s1(jj)=s1(jj) + fn*(an*pi(j)*p+bn*tau(j)*t)
        s2(jj)=s2(jj) + fn*(an*tau(j)*t+bn*pi(j)*p)
!        write(*,*) 'j=',j,' s1(j)=',s1(j),' s2(j)=',s2(j)
      end if
    end do
    
    psi0  = psi1
    psi1  = psi
    apsi1 = psi1
    chi0  = chi1
    chi1  = chi
    xi1   = cmplx(apsi1, -chi1)
    n     = n+1
    rn    = n
  
    do j = 1, nang
      pi1(j) = ((2._f * rn - 1._f) / (rn - 1._f)) * amu(j) * pi(j)
      pi1(j) = pi1(j) -rn * pi0(j) / (rn - 1._f)
      pi0(j) = pi(j)
    end do
  
    anmi1 = an
    bnmi1 = bn
    
    if (n - 1 - nstop >= 0) exit

  end do
  
  Qsca  = (2._f / (x * x)) * Qsca
  gfac  = (4._f / (x * x * Qsca)) * (g1+g2)
  Qext  = (4._f / (x * x)) * real(s1(1))
  Qback = (4._f / (x * x)) * abs(s1(2 * nang - 1)) * abs(s1(2 * nang - 1))
      
!  write(*,*) 'x',x,' s1(1)=',s1(1),' real(s1(1))=',real(s1(1))
!  write(*,*) 'Qsca=',Qsca,' gfac=',gfac,' Qext=',Qext,'Qback=',Qback

  deallocate(d)

  return
end subroutine bhmie
