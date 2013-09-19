#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module hybvcoord_mod

 use kinds,              only: r8 => real_kind, iulog
 use dimensions_mod,     only: plev => nlev, plevp => nlevp
 use physical_constants, only: p0

implicit none
private

!---------------------------------------------------------------------
! Purpose: Hybrid level definitions: p = A * p0 + B * ps
!          interfaces   p(k) = hyai(k)*ps0 + hybi(k)*ps
!          midpoints    p(k) = hyam(k)*ps0 + hybm(k)*ps
!---------------------------------------------------------------------
  type, public :: hvcoord_t
      real(r8) ps0          ! base state surface pressure for level definitions
      real(r8) hyai(plevp)  ! ps0 component of hybrid coordinate - interfaces
      real(r8) hyam(plev)   ! ps0 component of hybrid coordinate - midpoints
      real(r8) hybi(plevp)  ! ps  component of hybrid coordinate - interfaces
      real(r8) hybm(plev)   ! ps  component of hybrid coordinate - midpoints
      real(r8) hybd(plev)   ! difference in b (hybi) across layers
      real(r8) prsfac       ! log pressure extrapolation factor (time, space independent)
      real(r8) etam(plev)   ! eta - midpoint .  stored for conviencience
      real(r8) etai(plevp)  ! eta - interfaces
      integer  nprlev       ! number of pure pressure levels at top  
      integer  pad
  end type

  public :: hvcoord_init    ! function: initialize hvcoord from files

contains

  ! ==================================================================
  ! hvcoord_init: Initialize hybrid vertical coordinates.
  ! returns: hv coordinate structure, ierr=0 means success
  ! (After Boville's hycoef)
  ! ==================================================================

  function hvcoord_init(hvfile_mid, hvfile_int, lprint, masterproc, ierr) result(hvcoord)

    character(len=*), intent(in) :: hvfile_mid     ! file: mid-level vertical coords 
    character(len=*), intent(in) :: hvfile_int     ! file: interface-level vertical coords 
    logical,          intent(in) :: lprint         ! flag: print levels if true
    logical,          intent(in) :: masterproc     ! flag: is master process if true
    integer,          intent(out):: ierr           ! error code
    type (hvcoord_t)             :: hvcoord        ! output variable

    integer k, ierr11, ierr12, plevp_in, plev_in, ln
    real(r8) amean, bmean, atest, btest, eps

    ierr        = 0                                                   ! clear error flag
    hvcoord%ps0 = p0                                                  ! set base-state surface pressure (millibars)

#if (! defined ELEMENT_OPENMP)
!$OMP CRITICAL
#endif

    ln = len(trim(hvfile_int))                                        ! get length of file name
    if ( hvfile_int(ln-4:ln) == 'ascii' ) then                        ! check for .ascii suffix

    ! Read ascii format

       ! Open ascii interface file as unit 11
       open(UNIT=11,FILE=hvfile_int, form='formatted',status='old',iostat=ierr11)
       if (ierr11 /= 0) then
          write(iulog,*)'open() error:', __FILE__,__LINE__,hvfile_int
          ierr=ierr11
       end if

       ! Open ascii midpoint file as unit 12
       open(UNIT=12,FILE=hvfile_mid, form='formatted',access='sequential',status='old',iostat=ierr12)
       if (ierr12 /= 0) then
          write(iulog,*)'open() error:', __FILE__,__LINE__,hvfile_mid
          ierr=ierr12
       end if

       ! Read files if no error
       if(ierr==0) then

          ! Read A and B from ascii interface file

          read(11,*) plevp_in                                         ! read #A levels
          if (plevp_in .ne. plevp) then                               ! check for wrong # levels
             write(iulog,*) 'Error: hyai input file and HOMME plevp do not match',plevp,plevp_in
             ierr=1
          endif
          read(11,*) hvcoord%hyai(1:plevp)                            ! read component A

          read(11,*) plevp_in                                         ! read #B levels
          if (plevp_in .ne. plevp) then                               ! check for wrong # levels
             write(iulog,*) 'Error: hybi input file and HOMME plevp do not match',plevp,plevp_in
             ierr=1
          endif
          read(11,*)hvcoord%hybi(1:plevp)                             ! read component B

          ! Read A and B from ascii midpoint file

          read(12,*) plev_in                                          ! read #A levels
          if (plev_in .ne. plev) then                                 ! check for wrong # levels
             write(iulog,*) 'Error: hyam input file and HOMME plev do not match',plev,plev_in
             ierr=1
          endif
          read(12,*) hvcoord%hyam(1:plev)                             ! read component A

          read(12,*) plev_in                                          ! read #B levels
          if (plev_in .ne. plev) then                                 ! check for wrong # levels
             write(iulog,*) 'Error: hybi input file and HOMME plev do not match',plev,plev_in
             ierr=1
          endif
          read(12,*)hvcoord%hybm(1:plev)                              ! read component B

          close(11)                                                   ! close interface file
          close(12)                                                   ! close midpoint  file
       end if

    else

    ! Read binary format

       ! Open binary interface file as unit 11
       open(UNIT=11,FILE=hvfile_int, form='unformatted',status='old',iostat=ierr11)
       if (ierr11 /= 0) then
          write(iulog,*)'open() error:', __FILE__,__LINE__,hvfile_int
          ierr=ierr11
       end if

       ! Open binary midpoint file as unit 12
       open(UNIT=12,FILE=hvfile_mid, form='unformatted',access='sequential',status='old',iostat=ierr12)
       if (ierr12 /= 0) then
          write(iulog,*)'open() error:', __FILE__,__LINE__,hvfile_mid
          ierr=ierr12
       end if

       if(ierr==0) then	
          ! write(iulog,*)"plevp=",plevp; write(iulog,*)"Ai,Bi="

          ! Read A and B values for each interface level
          do k=1,plevp
             read(11) hvcoord%hyai(k)
             read(11) hvcoord%hybi(k)
             ! write(iulog,*)hvcoord%hyai(k),hvcoord%hybi(k)
          end do
          
          ! write(iulog,*)"plev=",plev; write(iulog,*)"Am,Bm="

          ! Read A and B values for each midpoint level
          do k=1,plev
             read(12) hvcoord%hyam(k)
             read(12) hvcoord%hybm(k)
             ! write(iulog,*)hvcoord%hyam(k),hvcoord%hybm(k)
          end do

          close(11)
          close(12)

       end if
    endif

#if (! defined ELEMENT_OPENMP)
!$OMP END CRITICAL
#endif
  if(ierr>0) return 

  eps            = 1.D-05
  hvcoord%nprlev = 0

  ! Compute eta = A + B for convenience
  do k=1,plev;  hvcoord%etam(k)=hvcoord%hyam(k)+hvcoord%hybm(k); enddo
  do k=1,plevp; hvcoord%etai(k)=hvcoord%hyai(k)+hvcoord%hybi(k); enddo

  ! ==================================
  ! Set layer locations
  ! ==================================

  ! ===============================================================
  ! Interfaces. Set nprlev to the interface above, the first time a 
  ! nonzero surface pressure contribution is found. "nprlev" 
  ! identifies the lowest pure pressure interface.
  ! ===============================================================

  do k=1,plev
     if (hvcoord%nprlev==0 .and. hvcoord%hybi(k).ne.0.0) hvcoord%nprlev = k - 1
  end do

  ! Set nprlev if no nonzero b's have been found. All interfaces are
  ! pure pressure. A pure pressure model requires other changes as well.
  
  if (hvcoord%nprlev==0) hvcoord%nprlev = plev + 2

  ! Set delta sigma part of layer thickness pressures
  do k=1,plev
     hvcoord%hybd(k) = hvcoord%hybi(k+1) - hvcoord%hybi(k)
  end do

  ! Calculate the log pressure extrapolation factor
#if (PLEV>1)
  hvcoord%prsfac = log( hvcoord%hyam(plev  ) + hvcoord%hybm(plev)) / &
                   log((hvcoord%hyam(plev  ) + hvcoord%hybm(plev)) / (hvcoord%hyam(plev-1) + hvcoord%hybm(plev-1)))
#endif

  ! Test that A, B at midpoints are means of A, B at interfaces

  do k = 1,plev
     amean = ( hvcoord%hyai(k+1) + hvcoord%hyai(k) )*0.5D0
     bmean = ( hvcoord%hybi(k+1) + hvcoord%hybi(k) )*0.5D0

     ! Compute error in A values
     if(amean == 0. .and. hvcoord%hyam(k) == 0.) then
        atest = 0.
     else
        atest = abs( amean - hvcoord%hyam(k) )/ ( 0.5D0*( abs(amean + hvcoord%hyam(k)) ) )
     endif

     ! Compute error in B values
     if(bmean == 0. .and. hvcoord%hybm(k) == 0.) then
        btest = 0.
     else
        btest = abs( bmean - hvcoord%hybm(k) )/ ( 0.5D0*( abs(bmean + hvcoord%hybm(k)) ) )
     endif

     ! Report midpoint A errors
     if (atest > eps) then
        if (masterproc) then
           write(iulog,9850)
           write(iulog,*)'k,atest,eps=',k,atest,eps
        end if
        ! call endrun
     endif

     ! Report midpoint B errors
     if (btest > eps) then
        if (masterproc) then
           write(iulog,9850)
           write(iulog,*)'k,btest,eps=',k,btest,eps
        end if
        ! call endrun
     endif
  end do

  ! Display layer locations if lprint = .true.
  if (masterproc) then
     if (lprint) then
        write(iulog,'(a)')'1 Layer Locations (*1000) A, B, A+B '
        do k=1,plev
           write(iulog,9800)k,hvcoord%hyai(k),hvcoord%hybi(k),hvcoord%hyai(k)+hvcoord%hybi(k)
           write(iulog,9810)  hvcoord%hyam(k),hvcoord%hybm(k),hvcoord%hyam(k)+hvcoord%hybm(k)
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

  end function hvcoord_init

end module hybvcoord_mod
