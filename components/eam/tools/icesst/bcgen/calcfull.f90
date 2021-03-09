subroutine calcfull (j, nlon, nlat, nmon, ismf,           &
                     jsmf, iyr1rd, mon1rd, iyr1, mon1,    &
                     iyrn, monn, iyrnrd, monnrd, maxiter, &
                     tmin, tmax, dt, conv, bbmin,         &
                     obsclim, arr, centmon, a, c,         &
                     correl, oldttcalc)
   use prec
   use solver, only: solvmid

   implicit none

   integer, intent(in) :: j
   integer, intent(in) :: nlon
   integer, intent(in) :: nlat
   integer, intent(in) :: nmon
   integer, intent(in) :: monn
   integer, intent(in) :: iyrn
   integer, intent(in) :: maxiter
   integer, intent(in) :: iyr1
   integer, intent(in) :: mon1
   integer, intent(in) :: iyr1rd
   integer, intent(in) :: mon1rd
   integer, intent(in) :: iyrnrd
   integer, intent(in) :: monnrd

   real(r8), intent(in) :: tmin
   real(r8), intent(in) :: tmax
   real(r8), intent(in) :: dt
   real(r8), intent(in) :: conv
   real(r8), intent(in) :: bbmin

   real(r8), intent(inout) :: obsclim(nlon,nlat,12)
   real(r8), intent(inout) :: arr(nlon,nmon)
   real(r8), intent(in) :: centmon(nlon,nlat,12)
   real(r8), intent(in) :: a(nmon)
   real(r8), intent(in) :: c(nmon)
   real(r8), intent(in) :: correl(12)

   logical, intent(in) :: oldttcalc     ! for bfb agreement with original code

   integer :: ismf, jsmf
   integer :: i
   integer :: n, m, mm
   integer :: nprior, nafter

   real(r8) :: vecin(nmon)
   real(r8) :: vecout(nmon)
   real(r8) :: cc
   real(r8) :: tt

   write(6,*) 'Calculating mid-month values for latitude ', j

   do i=1,nlon

! Fill in some months prior to the 1st observed month and
! after the last observed month

      nprior = 12*(iyr1rd-iyr1) + mon1rd - mon1 
      if (nprior > 0) then
         do m=1,nprior
            if ((nprior+1-m) > 12) then
               cc = 0.0
            else
               cc = correl(nprior+1-m)
            end if
            mm = mod((mon1-2+m), 12) + 1  
            tt = arr(i,nprior+1)*cc
            if (obsclim(i,j,mm)+tt > tmax) tt = tmax-obsclim(i,j,mm)
!
!JR Bugfix from K. Taylor: do not zero tt.  Keep to allow bfb vs. original code
!
            if (oldttcalc) then
               if (obsclim(i,j,mm) >= tmax) tt = 0.0
            end if
            if (obsclim(i,j,mm)+tt  < tmin) tt = tmin-obsclim(i,j,mm)
            if (oldttcalc) then
               if (obsclim(i,j,mm) <= tmin) tt = 0.0
            end if
            arr(i,m) = tt
         end do
      end if

      nafter = 12*(iyrn-iyrnrd) + monn - monnrd 

      if (nafter > 0) then
         do m=1,nafter
            if (m > 12) then
               cc = 0.0
            else
               cc = correl(m)
            end if
            mm = mod((monnrd-1+m), 12) + 1
            tt = arr(i,nmon-nafter)*cc
            if ((obsclim(i,j,mm)+tt) > tmax) tt = tmax-obsclim(i,j,mm)
!
!JR Bugfix from K. Taylor: do not zero tt.  Keep to allow bfb vs. original code
!
            if (oldttcalc) then
               if (obsclim(i,j,mm)  >= tmax) tt = 0.0
            end if
            if ((obsclim(i,j,mm)+tt) < tmin) tt = tmin-obsclim(i,j,mm)
            if (oldttcalc) then
               if (obsclim(i,j,mm)  <= tmin) tt = 0.0
            end if
            arr(i,nmon-nafter+m) = tt
         end do
      end if

! Copy data into vecin, vecout

      do n=1,nmon
         mm = mod((mon1+n-2), 12) + 1
         vecin(n) = arr(i,n) + obsclim(i,j,mm)
         vecout(n) = vecin(n)
      end do

      call solvmid (i, j, nmon, conv, dt, &
                    tmin, tmax, bbmin, maxiter, a, &
                    c, vecin, vecout, ismf, jsmf)

      do n=1,nmon
         mm = mod((mon1+n-2), 12) + 1
         arr(i,n) = vecout(n) - centmon(i,j,mm)
      end do
   end do

   return
end subroutine calcfull
