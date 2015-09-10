program hybgrid
#ifdef NETCDF
      use netcdf_coord_file
#endif
implicit none

integer, parameter :: nlev = 18
integer, parameter :: nlevp = nlev+1

! ============================================
! ln(p) Layer thickness increments specified in
! Williamson, et al. April 1998, Table B-1
! (p. 1011) for 18 layer case.
! ============================================

real*8, parameter  :: inc_slyr=.015D0
real*8, parameter  :: inc_blyr=.015D0
real*8, parameter  :: inc_tropo=.0267D0
real*8, parameter  :: inc_lstrat=.175D0
real*8, parameter  :: inc_mstrat=.30D0

! ================================================
! reference pressure (assumed == ps)
! ================================================

real*8, parameter  :: p0       = 1000.0D0  ! one atmosphere

! ================================================
! pressure for top of boundary layer, troposphere
! lower stratosphere and middle stratosphere.
! ================================================

real*8, parameter  :: p_blyr   = 970.0D0
real*8, parameter  :: p_tropo  =  90.0D0
real*8, parameter  :: p_lstrat =  60.0D0
real*8, parameter  :: p_mstrat =   3.0D0

! ================================================
! maximum layer thicknesses of each region
! as specified by Williamson, et. al.
! ================================================

real*8, parameter  :: max_inc_blyr   = 1.0D0
real*8, parameter  :: max_inc_tropo  = 0.35D0
real*8, parameter  :: max_inc_lstrat = 1.0D0
real*8, parameter  :: max_inc_mstrat = 1.0D0

! ===================================================
! define the pressure interface layer number 
! (inclusive) below which pure pressure coordinates
! reign (npr) 
! and the pressure interface layer number above which 
! (inclusive) pure sigma coordinates reign (nsig)
! ===================================================

integer, parameter :: npr  = 5
integer, parameter :: nsig = nlev

real*8 lnp(nlevp)   

real*8 A(nlevp)   
real*8 B(nlevp)

real*8 Amid(nlev)
real*8 Bmid(nlev)

real*8 vinc(nlev),dlyr(nlev)

real*8 pmid(nlev)
real*8 pmid_ccm(nlev)

real*8 pint(nlevp)
real*8 pint_ccm(nlevp)

real*8 inc          ! working layer increment
real*8 p            ! pressure
character(len=3)  :: charlev
character(len=80) :: fname
integer k

inc=inc_slyr
lnp(nlevp) = 0.0D0
lnp(nlev)  = -inc
vinc(nlev)= inc

p   = p0

! =================================================================
! read in 18 level hybrid A and B coefficients for comparison...
! =================================================================

write(charlev,'(i3)') nlev
fname="habi-18ref.fbin"
open(unit=7,file=fname,status="old",access="sequential", form="unformatted")
do k=1,nlevp
   read(7)A(k)
   read(7)B(k)
   pint_ccm(k)=p0*(A(k) + B(k))
end do
close(7)

print *,"A,B: interface layer values"
do k=1,nlevp
   print *,k,A(k),B(k)
end do
print *

fname="habm-18ref.fbin"
open(unit=7,file=fname,status="old",access="sequential", form="unformatted")
do k=1,nlev
   read(7)Amid(k)
   read(7)Bmid(k)
   pmid_ccm(k)=p0*(Amid(k) + Bmid(k))
end do
close(7)

print *,"A,B: mid layer values"
do k=1,nlev
   print *,k,Amid(k),Bmid(k)
end do
print *

! =====================================
! Algorithm to generate first pass
! values of ln(p) at interfaces...
! =====================================

do k=nlev-1,1,-1
   p = p0*EXP(lnp(k+1))
   if (p>p_blyr) then
      inc = MIN(inc + inc_blyr,max_inc_mstrat)
   else if( p<p_blyr .and. p>p_tropo ) then
      inc = MIN(inc + inc_tropo,max_inc_tropo)
   else if( p<p_tropo .and. p>p_lstrat) then
      inc = MIN(inc + inc_lstrat,max_inc_lstrat)
   else if( p<p_lstrat .and. p>p_mstrat) then
      inc = MIN(inc + inc_mstrat,max_inc_mstrat)
   end if
   lnp(k)=lnp(k+1) - inc
   vinc(k)=inc
end do

print *
print *,"first pass interface pressures..."

do k=1,nlevp
   pint(k)=p0*EXP(lnp(k))
   print *,k,pint(k)
end do

! ============================================
! Generate pmid from pint using
! geometric mean p(k)=SQRT(p(k-1/2)*p(k+1/2))
! ============================================

print *
print *,"mid level pressures vs those in CCM file..."
do k=1,nlev
   pmid(k)=SQRT(pint(k)*pint(k+1))
   print *,k,pmid(k),pmid_ccm(k)
end do

! ================================================
! Back compute interface pressures starting at
! surface using eq. 34 in Williamson and Olson MWR
! 1994, page 1601
! ================================================

pint(nlevp)=p0
do k=nlev,1,-1
   pint(k)=2.0D0*pmid(k) - pint(k+1)
end do

print *
print *,"back computed interface pressures vs those in CCM file..."
do k=1,nlevp
   print *,k,pint(k),pint_ccm(k)
end do

print *
print *,"computed A and B coefficients..."

do k=1,npr
   A(k)=pint(k)/p0
end do

do k=npr,nsig
   A(k) = (pint(npr)/p0)*( (pint(nsig)-pint(k))/(pint(nsig)-pint(npr)))
end do

do k=nsig+1,nlevp
   A(k) =0.0D0
end do

do k=1,nlevp
   B(k) = pint(k)/p0 - A(k)
end do

print *
print *,"computed interface level A and B:"
do k=1,nlevp
   print *,k,A(k),B(k)
end do

do k=1,nlev
   Amid(k) = 0.5D0*(A(k+1) + A(k))
   Bmid(k) = 0.5D0*(B(k+1) + B(k))
end do

print *
print *,"computed mid level A and B:"
do k=1,nlev
   print *,k,Amid(k),Bmid(k)
end do
#ifdef NETCDF
fname="hyvert-"//TRIM(ADJUSTL(charlev))//".nc"
call write_netcdf_coord_file(fname,nlev,p0,a,b,amid,bmid)
#else
fname="habi-"//TRIM(ADJUSTL(charlev))//".fbin"
open(unit=7,file=fname,status="unknown",access="sequential", form="unformatted")
do k=1,nlevp
   write(7)A(k)
   write(7)B(k)
end do
close(7)

fname="habm-"//TRIM(ADJUSTL(charlev))//".fbin"
open(unit=7,file=fname,status="unknown",access="sequential", form="unformatted")
do k=1,nlev
   write(7)Amid(k)
   write(7)Bmid(k)
end do
close(7)
#endif


end program hybgrid


