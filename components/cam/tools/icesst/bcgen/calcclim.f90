subroutine calcclim (j, nlon, nlat, nmon, ismc,               &
                     jsmc, mon1clm, iyr1clm, iyr1, mon1,      &
                     iyr1rd, mon1rd, iyrnrd, monnrd, monnclm, &
                     iyrnclm, maxiter, tmin, tmax, dt,        &
                     conv, bbmin, isea, obsclim, arr,         &
                     centmon, ac, cc)

   use prec
   use solver, only: solvmid

   implicit none

   integer, intent(in) :: j
   integer, intent(in) :: nlon
   integer, intent(in) :: nlat
   integer, intent(in) :: nmon
   integer, intent(in) :: mon1clm
   integer, intent(in) :: iyr1clm
   integer, intent(in) :: monnclm
   integer, intent(in) :: iyrnclm
   integer, intent(in) :: maxiter
   integer, intent(inout) :: isea
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
   real(r8), intent(out)   :: centmon(nlon,nlat,12)
   real(r8), intent(in) :: ac(12)
   real(r8), intent(in) :: cc(12)

   integer :: ismc, jsmc
   integer :: i, k
   integer :: n, m, nn, mn, mm
   integer :: momax, momin, mmax, mmin, mmnn
   integer :: m1, m2, m3, m4

   real(r8) :: vecin(12)
   real(r8) :: vecout(12)
   real(r8) :: ocentmax
   real(r8) :: ocentmin
   real(r8) :: centmax
   real(r8) :: centmin
   real(r8) :: rmin, rmax
   real(r8) :: rrmin, rrmax

   mmnn = 0
   rrmin = 0.0
   rrmax = 0.0

   m1 = 12*(iyr1rd-iyr1) + mon1rd-mon1 + 1
   m2 = 12*(iyrnrd-iyr1) + monnrd-mon1 + 1

   m3 = 12*(iyr1clm-iyr1) + mon1clm-mon1 + 1
   m4 = 12*(iyrnclm-iyr1) + monnclm-mon1 + 1

   write(6,*) 'computing climatology for latitude ', j

   do i=1,nlon
      do n=1,12
         obsclim(i,j,n) = 0.0
      end do
      mn = 0
      rmin = 0.0
      rmax = 0.0
      do n=m1,m2
         if (arr(i,n) > tmax) then
            if (arr(i,n) > tmax+1.0e-2) then
               write(6,*) 'tmax exceeded ', i, j, n, arr(i,n)
               rmax = max (rmax,arr(i,n)-tmax)
               mn = mn + 1
            end if
            arr(i,n) = tmax
         else if (arr(i,n) < tmin) then
            if (arr(i,n) < tmin-dt) then
               if ((tmin-arr(i,n)) < 900.) then
                  rmin = min (rmin,arr(i,n)-tmin)
                  mn = mn + 1
               end if
            end if
            arr(i,n) = tmin                  
         end if
      end do
      
      if (mn > 0) then
         mmnn = mmnn + 1
         rrmin = min (rrmin, rmin)
         rrmax = max (rrmax, rmax)
         write(6,*) ' '
         write(6,*) ' WARNING -- observed value exceeds limits at ', mn, &
                    ' time points at latitude ', j
         write(6,*) 'and longitude ', i
         write(6,*) 'max error = ', rmax, '  min error = ', rmin
      end if

! Compute climatology  

      do n=1,12
         m = mod((mon1clm+n-2), 12) + 1
         obsclim(i,j,m)=0.0
         nn = 0
         do k = m3+n-1, m4, 12 
            nn = nn + 1
            obsclim(i,j,m) = obsclim(i,j,m) + arr(i,k)
         end do
         obsclim(i,j,m) = obsclim(i,j,m)/nn
      end do
      
! Remove climatology to generate anomalies

      do n=m1,m2
         mm = mod((mon1rd+n-m1-1), 12) + 1
         arr(i,n) = arr(i,n) - obsclim(i,j,mm)
      end do
   end do

   if (mmnn > 0) then 
      write(6,*) ' '
      write(6,*) ' WARNING -- observed value exceeds limits at'
      write(6,*) mmnn, ' grid cells'
      write(6,*) 'max error = ', rrmax, '  min error = ', rrmin
   end if

! Solve for climatological mid-month values

   write(6,*) 'Computing climatological mid-month values for latitude ', j
      
   do i=1,nlon
      do m=1,12
         vecin(m) = obsclim(i,j,m)
         vecout(m) = obsclim(i,j,m)
      end do

      call solvmid (i, j, 12, conv, dt, &
                    tmin, tmax, bbmin, maxiter, ac, &
                    cc, vecin, vecout, ismc, jsmc)

      do m=1,12
         centmon(i,j,m) = vecout(m)
      end do
            
      ocentmax = -1.e20
      ocentmin = 1.e20
      centmax = -1.e20
      centmin = 1.e20

      do m=1,12

! Find max and min values and months of max and min values

         if (obsclim(i,j,m) > ocentmax) then
            ocentmax = obsclim(i,j,m)
            momax = m
         end if
         
         if (obsclim(i,j,m) < ocentmin) then
            ocentmin = obsclim(i,j,m)
            momin = m
         end if
         
         if (centmon(i,j,m) > centmax) then
            centmax = centmon(i,j,m)
            mmax = m
         end if
         
         if (centmon(i,j,m) < centmin) then
            centmin = centmon(i,j,m)
            mmin = m
         end if
      end do
      isea = isea + 1   ! count ocean grid cells
   end do
   
   return
end subroutine calcclim
