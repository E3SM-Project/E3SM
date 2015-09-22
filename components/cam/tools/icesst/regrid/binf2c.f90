subroutine binf2c (plonin,  nlatin,  lonin, latin, arrin, &
                   plonout, nlatout, nlon, rlon, latout, arrout, &
                   verbose)
!-----------------------------------------------------------------------
!
! Purpose: bin from a finer grid to a coarser one, taking account of missing or
! filled data points.  Configured to work correctly on a reduced grid.  Algorithm: 
! For each point on the fine grid with valid data, 1st find the nearest latitude 
! on the coarse mesh.  Then for that latitude find the nearest coarse longitude and 
! put the fine grid point into that coarse "bin".
!
!-----------------------------------------------------------------------
   use precision

   implicit none
!
! Arguments
!
   integer, intent(in) :: plonin                    ! longitude dimension of input
   integer, intent(in) :: nlatin                    ! latitude dimension of input
   real(r8), intent(in) :: lonin(plonin)            ! input longitudes
   real(r8), intent(in) :: latin(nlatin)            ! input latitudes
   real(r8), intent(in) :: arrin(plonin,nlatin)     ! input array

   integer, intent(in) :: plonout                   ! longitude dimension of output
   integer, intent(in) :: nlatout                   ! latitude dimension of output
   integer, intent(in) :: nlon(nlatout)             ! lons (deg.) at each lat (maybe on reduced grid)
   real(r8), intent(in) :: rlon(plonout,nlatout)    ! longitudes (deg.) at each latitude (maybe on reduced grid)
   real(r8), intent(in) :: latout(nlatout)          ! output latitudes
   real(r8), intent(out) :: arrout(plonout,nlatout) ! output array

   logical, intent(in) :: verbose                   ! added printout
!
! Local workspace
!  
   integer :: i,j
   integer :: ii,jj
   integer :: iiarr(1), jjarr(1)
   integer :: num
   integer :: bincount(plonout,nlatout)

   real(r16) :: arrloc(plonout,nlatout)           ! output array in real*16
   real(r8)  :: deltay(nlatout)
   real(r8)  :: deltax(plonout)

   if (nlatin < nlatout) then
      write(6,*)'Warning: input grid coarser than output grid'
   end if

   bincount(:,:) = 0
   arrloc(:,:) = 0.
   arrout(:,:) = 0.
!
! Find closest output grid lon index (ii) and lat index (jj) to input grid lon (i) and lat (j)
! Then sum arr into appropriate bin
!
   do j=1,nlatin
      deltay(:) = abs (latin(j) - latout(:))
      jjarr = minloc (deltay(:))                  ! lat index of closest output grid point
      jj = jjarr(1)
      do i=1,plonin
         num = nlon(jj)
         deltax(:num) = abs (lonin(i) - rlon(:num,jj))
         iiarr = minloc (deltax(:num))            ! lon index of closest output grid point
         ii = iiarr(1)
         arrloc(ii,jj) = arrloc(ii,jj) + arrin(i,j)
         bincount(ii,jj) = bincount(ii,jj) + 1
      end do
   end do
!
! Normalize by bin count
!
   do jj=1,nlatout
      do ii=1,nlon(jj)
         if (bincount(ii,jj) > 0) then
            arrout(ii,jj) = arrloc(ii,jj)/bincount(ii,jj)
         else
            write(6,*)'binf2c: Bincount(i=',ii,',j=',jj,') = 0: stopping'
            stop 999
         end if
      end do
   end do

   if (verbose) then
      write(6,*)'bincount:'
      do jj=1,nlatout
         write(6,'(1000i2)') (bincount(ii,jj),ii=1,nlon(jj))
      end do
   end if

   return
end subroutine binf2c
