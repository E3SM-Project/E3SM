subroutine massborrow(subnam,lchnk,ncol,pcols,mbeg,mend,qmin,q,pdel) 

!!....................................................................... 
!! The mass borrower borrows tracer mass from an adjacent layer. 
!! It conserves the mass and can avoid negative tracers. 
!! 
!! At level k, it will first borrow the mass from the layer k+1 (lower level). 
!! If the mass is not sufficient in layer k+1, it will borrow mass from 
!! layer k+2. The borrower will proceed this process until the bottom layer. 
!! If the tracer mass in the bottom layer goes negative, it will repeat the 
!! process from the bottom to the top. In this way, the borrower works for 
!! any shape of mass profiles.
!! 
!! The code is adapted from the tracer mass borrower implemented in the 
!! global aerosol-climate model ECHAM-HAM (Feichter et al.,1996; 
!! Stier et al., 2005, Zhang et al., 2012).
!! 
!! Author : Kai Zhang (kai.zhang@pnnl.gov) 
!!....................................................................... 

  use shr_kind_mod,    only: r8 => shr_kind_r8
  use ppgrid,          only: pver
  use spmd_utils,      only: masterproc
  use cam_logfile,     only: iulog
  use phys_control,    only: print_fixer_message

  implicit none

!! interface 
!!....................................................................... 

  character*(*), intent(in) :: subnam                 ! name of calling routine
  integer, intent(in) :: lchnk                        ! chunk identifier
  integer, intent(in) :: ncol                         ! number of atmospheric columns
  integer, intent(in) :: pcols                        ! number of dim members
  integer, intent(in) :: mbeg                         ! first index 
  integer, intent(in) :: mend                         ! last index 
  real(r8), intent(in) :: qmin(mbeg:mend)             ! smallest value
  real(r8), intent(in) :: pdel(pcols,pver)            ! pressure thickness 
  real(r8), intent(inout) :: q(pcols,pver,mbeg:mend)  ! moisture/tracer field

!! local 
!!....................................................................... 

  integer :: i, k, m, j 
  integer :: ic(pcols) 
  real(r8):: nmass, zeps
  real(r8):: bmass(pcols)

  !! init
  !!....................................................................... 

  zeps = epsilon(1.0_r8)

  !! loop over tracers
  !!....................................................................... 

  do m = mbeg, mend

     ic(1:ncol) = 0 

     bmass(1:ncol) = 0.0_r8
     
     !! top to bottom
     !!....................................................................... 

     do k = 1, pver
        do i = 1, ncol

           !! new mass in the current layer
           !!....................................................................... 

           nmass = q(i,k,m) + bmass(i)/pdel(i,k)

           if ( nmass > qmin(m) ) then

              !! if new mass in the current layer is positive, don't borrow mass any more 
              !!....................................................................... 

              q(i,k,m) = nmass
              bmass(i) = 0.0_r8

           else

              !! set mass to qmin in the current layer, and save bmass
              !!....................................................................... 

              bmass(i) = (nmass - qmin(m)) * pdel(i,k)

              q(i,k,m) = qmin(m) 

              ic(i) = ic(i) + 1 

           end if !! nmass > 0.0_r8 

        end do !! i 
     end do !! k 

!!     do i = 1, ncol
!!
!!        if(print_fixer_message .and. ic(i).gt.0) then 
!!            write(iulog,*) '### mass borrower T2B ### tracer : ', m, ' column : ', i, ' chunk : ', lchnk  
!!        end if 
!!
!!        if(print_fixer_message .and. bmass(i) < 0._r8 ) then 
!!            write(iulog,*) '### mass borrower B2T ### tracer : ', m, ' column : ', i, ' chunk : ', lchnk  
!!        end if 
!!
!!     end do 

     !!....................................................................... 
     !! bottom to top
     !!....................................................................... 
     
     do k = pver, 1, -1 

        do i = 1, ncol

           !! if the surface layer still needs to borrow mass 
           !!....................................................................... 

           if (bmass(i) < 0._r8 ) then

              !! new mass in the current layer
              !!....................................................................... 

              nmass = q(i,k,m) + bmass(i)/pdel(i,k)

              if ( nmass > qmin(m) ) then

                 !! if new mass in the current layer is positive, don't borrow mass any more 
                 !!....................................................................... 

                 q(i,k,m) = nmass 
                 bmass(i) = 0.0_r8

              else

                 !! if new mass in the current layer is negative, continue to borrow mass
                 !!....................................................................... 

                 bmass(i) = (nmass - qmin(m))*pdel(i,k)
                 q(i,k,m) = qmin(m)

              end if !! nmass > 0.0_r8 

           end if !! bmass(i) < -zeps 

        end do !! i 

     end do !! k 

!!!     do k = 1, pver
!!!     do i = 1, ncol
!!!        if(print_fixer_message .and. q(i,k,m).lt.qmin(m)) then 
!!!            write(iulog,*) '### massborrow ### index : ', m, ' column : ', i, ' level : ', k, ' chunk : ', lchnk, ' tracer conc : ', q(i,k,m) 
!!!        end if 
!!!     end do 
!!!     end do 

  end do !! m

  return 
  end subroutine massborrow

