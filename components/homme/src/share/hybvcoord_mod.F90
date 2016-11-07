#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module hybvcoord_mod
use kinds,              only: r8 => real_kind, iulog
use dimensions_mod,     only: plev => nlev, plevp => nlevp
use physical_constants, only: p0
use control_mod,        only: vanalytic

implicit none
private

!----------------------------------------------------------------------- 
! hvcoord_t: Hybrid level definitions: p = a*p0 + b*ps
!            interfaces   p(k) = hyai(k)*ps0 + hybi(k)*ps
!            midpoints    p(k) = hyam(k)*ps0 + hybm(k)*ps
!-----------------------------------------------------------------------
type, public :: hvcoord_t
  real(r8) ps0          ! base state surface-pressure for level definitions
  real(r8) hyai(plevp)  ! ps0 component of hybrid coordinate - interfaces
  real(r8) hyam(plev)   ! ps0 component of hybrid coordinate - midpoints
  real(r8) hybi(plevp)  ! ps  component of hybrid coordinate - interfaces
  real(r8) hybm(plev)   ! ps  component of hybrid coordinate - midpoints
  real(r8) hybd(plev)   ! difference in b (hybi) across layers
  real(r8) prsfac       ! log pressure extrapolation factor (time, space independent)
  real(r8) etam(plev)   ! eta-levels at midpoints
  real(r8) etai(plevp)  ! eta-levels at interfaces
  integer  nprlev       ! number of pure pressure levels at top  
  integer  pad
end type

public :: hvcoord_init, set_layer_locations
contains

  !_____________________________________________________________________
  function hvcoord_init(hvfile_mid, hvfile_int, lprint, masterproc, ierr) result(hvcoord)

    ! hvcoord_init: Initialize hybrid vertical coordinate system
    ! returns: hv coordinate structure, ierr=0 means success
    ! Modeled after Boville's hycoef

    character(len=*), intent(in) :: hvfile_mid      ! file containing mid levels vertical coordinate system
    character(len=*), intent(in) :: hvfile_int      ! file containing interface levels vertical coordinate system
    logical, intent(in)          :: lprint
    logical, intent(in)          :: masterproc
    integer, intent(out)         :: ierr
    type (hvcoord_t)             :: hvcoord
    integer k,ierr11, ierr12, plevp_in,plev_in, ln

    ierr=0            ! clear error state
    hvcoord%ps0 = p0  ! set base-state surface pressure

    hvcoord%hyai=0.0d0; hvcoord%hybi=0.0d0
    hvcoord%hyam=0.0d0; hvcoord%hybm=0.0d0

    if(vanalytic == 1) return ! exit if setting vertical coordinates analytically

#if (defined HORIZ_OPENMP)
!$OMP CRITICAL
#endif

    ln=len(trim(hvfile_int))
    if ( hvfile_int(ln-4:ln) == 'ascii' ) then

       ! Open ascii interface-level file
       open(UNIT=11,FILE=hvfile_int, form='formatted',status='old',iostat=ierr11)
       if (ierr11 /= 0) then
          write(iulog,*)'open() error:', __FILE__,__LINE__,hvfile_int
          ierr=ierr11
       end if

       ! Open ascii mid-level file
       open(UNIT=12,FILE=hvfile_mid, form='formatted',access='sequential',status='old',iostat=ierr12)
       if (ierr12 /= 0) then
          write(iulog,*)'open() error:', __FILE__,__LINE__,hvfile_mid
          ierr=ierr12
       end if

       if(ierr==0) then
          ! Read A interface coefficients from ascii file
          read(11,*) plevp_in
          if (plevp_in .ne. plevp) then
             write(iulog,*) 'Error: hyai input file and HOMME plevp do not match',plevp,plevp_in
             ierr=1
          endif
          read(11,*)hvcoord%hyai(1:plevp_in)

          ! Read B interface coefficients from ascii file
          read(11,*) plevp_in
          if (plevp_in .ne. plevp) then
             write(iulog,*) 'Error: hybi input file and HOMME plevp do not match',plevp,plevp_in
             ierr=1
          endif
          read(11,*)hvcoord%hybi(1:plevp_in)

          ! Read A midpoint coefficients from ascii file
          read(12,*) plev_in
          if (plev_in .ne. plev) then
             write(iulog,*) 'Error: hyam input file and HOMME plev do not match',plev,plev_in
             ierr=1
          endif
          read(12,*) hvcoord%hyam(1:plev_in)

          ! Read B midpoint coefficients from ascii file
          read(12,*) plev_in
          if (plev_in .ne. plev) then
             write(iulog,*) 'Error: hybm input file and HOMME plev do not match',plev,plev_in
             ierr=1
          endif
          read(12,*)hvcoord%hybm(1:plev_in)

          close(11)
          close(12)
       end if
    else

       ! Open binary interface-level file
       open(UNIT=11,FILE=hvfile_int, form='unformatted',status='old',iostat=ierr11)
       if (ierr11 /= 0) then
          write(iulog,*)'open() error:', __FILE__,__LINE__,hvfile_int
          ierr=ierr11
       end if

       ! Open ascii midpoint-level file
       open(UNIT=12,FILE=hvfile_mid, form='unformatted',access='sequential',status='old',iostat=ierr12)
       if (ierr12 /= 0) then
          write(iulog,*)'open() error:', __FILE__,__LINE__,hvfile_mid
          ierr=ierr12
       end if

       if(ierr==0) then
          ! Read A and B interface coefficients from binary file
          do k=1,plevp
             read(11)hvcoord%hyai(k)
             read(11)hvcoord%hybi(k)
          end do

          ! Read A and B midpoint coefficients from binary file
          do k=1,plev
             read(12)hvcoord%hyam(k)
             read(12)hvcoord%hybm(k)
          end do
          
          close(11)
          close(12)
       end if
    endif

#if (defined HORIZ_OPENMP)
!$OMP END CRITICAL
#endif

  if(ierr>0) return ! Exit if an error occured

  ! Set eta levels from A,B coefficients
  call set_layer_locations(hvcoord,lprint,masterproc)

end function

!_______________________________________________________________________
  subroutine set_layer_locations(hvcoord, lprint, masterproc)

  ! Set eta layers and other quantities derived from A,B vertical coefficients

  type (hvcoord_t), intent(inout) :: hvcoord
  logical,          intent(in)    :: lprint
  logical,          intent(in)    :: masterproc

  real(r8) :: amean,bmean,atest,btest,eps
  integer  :: k

  eps            = 1.D-05

  ! Get eta-levels from A,B
  forall(k=1:plev)  hvcoord%etam(k)=hvcoord%hyam(k)+hvcoord%hybm(k)
  forall(k=1:plevp) hvcoord%etai(k)=hvcoord%hyai(k)+hvcoord%hybi(k)

  ! Set nprlev to the lowest pure pressure interface
  hvcoord%nprlev = 0
  do k=1,plev; if(hvcoord%nprlev==0 .and. hvcoord%hybi(k).ne.0.0) hvcoord%nprlev = k - 1; enddo

  ! Set nprlev if all interfaces are pure pressure
  if (hvcoord%nprlev==0) hvcoord%nprlev = plev + 2

  ! Set delta-sigma part of layer thickness pressures
  forall(k=1:plev) hvcoord%hybd(k) = hvcoord%hybi(k+1)-hvcoord%hybi(k)

  ! Calculate the log pressure extrapolation factor
#if (PLEV>1)
  hvcoord%prsfac = log( hvcoord%hyam(plev) + hvcoord%hybm(plev)) / &
                   log((hvcoord%hyam(plev) + hvcoord%hybm(plev)) / (hvcoord%hyam(plev-1) + hvcoord%hybm(plev-1)))
#endif

  ! ======================================================================
  ! Test that midpoint A,B is mean of interface A,B
  ! ======================================================================
  do k = 1,plev

     amean = ( hvcoord%hyai(k+1) + hvcoord%hyai(k) )*0.5D0
     bmean = ( hvcoord%hybi(k+1) + hvcoord%hybi(k) )*0.5D0

     if(amean == 0. .and. hvcoord%hyam(k) == 0.) then
        atest = 0.
     else
        atest = abs( amean - hvcoord%hyam(k) )/ ( 0.5D0*( abs(amean + hvcoord%hyam(k)) ) )
     endif
     if(bmean == 0. .and. hvcoord%hybm(k) == 0.) then
        btest = 0.
     else
        btest = abs( bmean - hvcoord%hybm(k) )/ ( 0.5D0*( abs(bmean + hvcoord%hybm(k)) ) )
     endif

     if (atest > eps) then
        if (masterproc) then
           write(iulog,9850)
           write(iulog,*)'k,atest,eps=',k,atest,eps
           write(iulog,*)'hyai k,k+1:  ',hvcoord%hyai(k),hvcoord%hyai(k+1)
           write(iulog,*)'hyam k:  ',hvcoord%hyam(k)
        end if
     endif

     if (btest > eps) then
        if (masterproc) then
           write(iulog,9850)
           write(iulog,*)'k,btest,eps=',k,btest,eps
           write(iulog,*)'hybi k,k+1:  ',hvcoord%hybi(k),hvcoord%hybi(k+1)
           write(iulog,*)'hybm k:  ',hvcoord%hybm(k)
        end if
     endif
  end do

  if (masterproc) then
     if (lprint) then
        write(iulog,'(a)')'1 Layer Locations (*1000) '
        do k=1,plev
           write(iulog,9800)k,hvcoord%hyai(k),hvcoord%hybi(k),hvcoord%hyai(k)+hvcoord%hybi(k)
           write(iulog,9810) hvcoord%hyam(k), hvcoord%hybm(k), hvcoord%hyam(k)+hvcoord%hybm(k)
        end do
        write(iulog,9800)plevp,hvcoord%hyai(plevp),hvcoord%hybi(plevp),hvcoord%hyai(plevp)+hvcoord%hybi(plevp)
     else
        ! sanity check for endian problem with file:
        write(iulog,*)'min/max hybm() coordinates: ',minval(hvcoord%hybm(1:plev)),maxval(hvcoord%hybm(1:plev))
     endif
  end if
  
9800 format( 1x, i3, 3p, 3(f10.4,10x) )
9810 format( 1x, 3x, 3p, 3(10x,f10.4) )
9850 format('HYCOEF: A and/or B vertical level coefficients at full',/, &
            ' levels are not the arithmetic mean of half-level values')

  end subroutine

end module hybvcoord_mod
