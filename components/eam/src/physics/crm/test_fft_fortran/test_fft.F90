program test_MT_RNG
   use fftpack51D
   implicit none
   !----------------------------------------------------------------------------
   !----------------------------------------------------------------------------

   integer, parameter :: rknd = selected_real_kind(13)

   real(rknd),parameter :: pi=3.14159 
   integer,   parameter :: nx=16        ! length of the sequence
   integer,   parameter :: lensav=nx+15 ! must be at least N + INT(LOG(REAL(N))) + 4.

   real(rknd), dimension(nx) :: arr_in    ! input data
   real(rknd), dimension(nx) :: arr_out1  ! data after transforming forward and back
   real(rknd), dimension(nx) :: arr_out2  ! data after reconstructing with amp/phs
   real(rknd), dimension(nx) :: arr_out3  ! data after reconstructing with sin/cos
   real(rknd), dimension(nx) :: fft_out   ! FFT weights
   real(rknd), dimension(nx) :: wave_num  ! wave number
   real(rknd), dimension(nx) :: wave_len  ! wave length
   real(rknd), dimension(nx) :: amp       ! amplitude per wavenumber
   real(rknd), dimension(nx) :: phs       ! phase per wavenumber
   
   real(rknd) :: dx     ! ???
   real(rknd) :: angle  ! angle for creating input data
   real(rknd) :: phase  ! phase offset for input data
   real(rknd) :: freq   ! 
   real(rknd) :: diff1  ! diff for sanity check
   real(rknd) :: diff2  ! diff for sanity check
   real(rknd) :: diff3  ! diff for sanity check

   integer :: i,j       ! loop iterators
   integer :: ier       ! error return code
   integer :: nh        ! ???
   
   real(rknd), dimension(lensav) :: wsave  ! init output: prime factors of N and certain trigonometric values to be used in rfft1b or rfft1f.
   real(rknd), dimension(nx)     :: work

   !----------------------------------------------------------------------------
   ! create some data to transform
   !----------------------------------------------------------------------------
   do i = 1,nx
      angle = (2*pi)*(i-1)/nx
      phase = (2*pi)*0.4
      arr_in(i) = 100 
      arr_in(i) = arr_in(i) + 10*cos(1*angle)
      arr_in(i) = arr_in(i) + 10*sin(1*angle)
      ! arr_in(i) = arr_in(i) + 3*sin(2*angle-phase)
      ! initialize output array
      fft_out(i) = arr_in(i)
      amp(i) = 0
      phs(i) = 0
   end do

   !----------------------------------------------------------------------------
   ! Do the FFT
   !----------------------------------------------------------------------------

   ! initialization for FFT
   call rfft1i(nx,wsave,lensav,ier)
   if(ier /= 0) write(0,*) 'ERROR: rfftmi(): ESMT - FFT initialization error ',ier

   ! do the forward transform
   call rfft1f( nx, 1, fft_out(:), nx, wsave, lensav, work(:), nx, ier )
   if(ier /= 0) write(0,*) 'ERROR: rfftmf(): ESMT - Forward FFT error ',ier

   do i = 1,nx
      arr_out1(i) = fft_out(i)
      write(*,900) i,fft_out(i)
   end do
900 format(i6,'   ',f16.8)

   stop

   ! transform back
   call rfft1b( nx, 1, arr_out1(:), nx, wsave, lensav, work(:), nx, ier )
   if(ier /= 0) write(0,*) 'ERROR: rfftmb(): ESMT - backward FFT error ',ier

   !----------------------------------------------------------------------------
   !----------------------------------------------------------------------------

   ! ???
   if(mod(nx,2) == 0) then
      nh = nx/2 - 1
   else
      nh = (nx-1)/2
   endif

   do i = 1,nx
      ! phs(i) = atan2(fft_out(i+1),fft_out(i))
      ! amp(i) = sqrt( fft_out(i)**2. + fft_out(i+1)**2. )
      arr_out2(i) = 10
      arr_out3(i) = 10
   end do

   ! set wavenumbers
   dx = 1
   wave_num(1)=0.0
   do j = 1,nh
      ! wave_num(2*j)   = 2.*pi* (real(j,rknd)/real(nx,rknd)) /dx   ! cos
      ! wave_num(2*j+1) = 2.*pi* (real(j,rknd)/real(nx,rknd)) /dx   ! sin
      wave_num(2*j)   = nx/dx * (real(j,rknd)/real(nx,rknd))   ! cos
      wave_num(2*j+1) = nx/dx * (real(j,rknd)/real(nx,rknd))   ! sin
      wave_len(2*j)   = nx*dx * 1/wave_num(2*j)  
      wave_len(2*j+1) = nx*dx * 1/wave_num(2*j+1)
      phs(2*j) = atan2(fft_out(2*j+1),fft_out(2*j))
      amp(2*j) = sqrt( fft_out(2*j)**2. + fft_out(2*j+1)**2. )
      do i = 1,nx
         angle = (2*pi)*(i-1)/nx
         freq  = wave_num(2*j) !2.*pi* (real(j,rknd)/real(nx,rknd))
         arr_out2(i) = arr_out3(i) + amp(2*j)*cos(freq*angle-phs(2*j))
         arr_out3(i) = arr_out3(i) + fft_out(2*j)*cos(freq*angle) + fft_out(2*j+1)*sin(freq*angle)
      end do
   enddo

   ! ???
   if (mod(nx,2) == 0) wave_num(nx) = 2.*pi/(2.*dx)              ! nyquist wavelength for even n

   !----------------------------------------------------------------------------
   ! print sanity checks
   !----------------------------------------------------------------------------
   
   write(*,*)
   write(*,100) 'i','data in','data out 1','data out 2','data out 3','diff 1','diff 2','diff 3'
   do i = 1,nx
      diff1 = arr_in(i) - arr_out1(i)
      diff2 = arr_in(i) - arr_out2(i)
      diff3 = arr_in(i) - arr_out3(i)
      write(*,101) i, arr_in(i), arr_out1(i), arr_out2(i), arr_out3(i), diff1, diff2, diff3
   end do
100 format(A12,A16,A16,A16,A16,A16,A16,A16)
101 format(i12,f16.8,f16.8,f16.8,f16.8,f16.8,f16.8,f16.8)

   ! print FFT weights
   write(*,*)
   write(*,200) 'i','wave num' ,'wave len','FFT cos wgt','FFT sin wgt','amp','phase'
   do i = 1,nh
      write(*,201) i, wave_num(2*i), wave_len(2*i), fft_out(2*i), fft_out(2*i+1), amp(2*i), phs(2*i) 
   end do
200 format(A12,A16,A16,A16,A16,A16,A16)
201 format(i12,f16.8,f16.8,f16.8,f16.8,f16.8,f16.8)

   write(*,*)

   !----------------------------------------------------------------------------
   !----------------------------------------------------------------------------

end program test_MT_RNG