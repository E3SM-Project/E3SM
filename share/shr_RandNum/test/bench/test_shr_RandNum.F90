program test

! this program calls the available versions of the random generators

use shr_RandNum_mod, only: ShrRandGen, ShrIntrinsicRandGen, ShrKissRandGen, &
     ShrF95MtRandGen, ShrDsfmtRandGen, ShrMklMtRandGen

INTEGER, parameter :: r8 = selected_real_kind(12)

class(ShrRandGen), allocatable :: randStream

integer, parameter :: nstream = 16   ! number of streams of random numbers
integer, parameter :: length  = 1000 ! length of stream of random numbers
integer            :: ntrials = 50000

integer, dimension(nstream) :: seed = 7776578
integer, dimension(nstream,4) :: kiss_seed
integer, dimension(:,:), allocatable :: intrinsic_seed

real(r8), dimension(nstream,length) :: array

integer :: i, n, m, intrinsic_size
integer :: c1, c2, cr, cm
real(r8) :: dt, dt1,dt2

#ifdef INTEL_MKL
! intel math kernel library mersenne twister

  allocate(randStream, source=ShrMklMtRandGen(seed))

  call system_clock(c1, cr, cm)
  do m = 1,ntrials
     call randStream%finalize()
     deallocate(randStream)
     allocate(randStream, source=ShrMklMtRandGen(seed))
  enddo
  call system_clock(c2, cr, cm); dt1 = dble(c2-c1)/dble(cr)

  call system_clock(c1, cr, cm)
  do m = 1,ntrials
    call randStream%random(array)
  enddo
  call randStream%finalize()
  deallocate(randStream)
  call system_clock(c2, cr, cm); dt2 = dble(c2-c1)/dble(cr)
  dt = dt1+dt2
  print *, 'Init time   (SFMT_MKL): ',dt1
  print *, 'Gen  time   (SFMT_MKL): ',dt2
  print *, 'Total time  (SFMT_MKL): ',dt
  print *, 'MegaRNumbers (SFMT_MKL): ', 1.0e-6*dble(nstream*length*ntrials)/dt
  print *, 'Summation of Random Numbers: ', SUM(array)
  print *, '--------'; print *, ''
#endif

! keep it simple stupid random number

  call system_clock(c1, cr, cm)
  do m = 1,ntrials
    do n = 1,nstream
       do i = 1, 4
          kiss_seed(n,i) = seed(n)*i+n
       end do
    end do

    allocate(randStream, source=ShrKissRandGen(kiss_seed))

    call randStream%random(array)

    call randStream%finalize()
    deallocate(randStream)
  enddo
  call system_clock(c2, cr, cm); dt = dble(c2-c1)/dble(cr)

  print *, 'Total time   (KISSVEC): ',dt
  print *, 'MegaRNumbers (KISSVEC): ', 1.0e-6*dble(nstream*length*ntrials)/dt
  print *, 'Summation of Random Numbers: ', SUM(array)
  print *, '--------'; print *, ''

! fortran-95 implementation of merseene twister

  call system_clock(c1, cr, cm)
  do m = 1,ntrials
    allocate(randStream, source=ShrF95MtRandGen(seed))

    call randStream%random(array)

    call randStream%finalize()
    deallocate(randStream)
  enddo
  call system_clock(c2, cr, cm); dt = dble(c2-c1)/dble(cr)

  print *, 'Total time   (MT19937): ',dt
  print *, 'MegaRNumbers (MT19937): ', 1.0e-6*dble(nstream*length*ntrials)/dt
  print *, 'Summation of Random Numbers: ', SUM(array)
  print *, '--------'; print *, ''

! fortran-90 intrinsic pseudorandom number generator

  call system_clock(c1, cr, cm)
  call random_seed(size=intrinsic_size)
  allocate(intrinsic_seed(nstream,intrinsic_size))
  do m = 1,ntrials
     do n = 1, nstream
        do i = 1, intrinsic_size
           intrinsic_seed(n,i) = seed(n)*i+n
        end do
     end do
    allocate(randStream, source=ShrIntrinsicRandGen(intrinsic_seed))

    call randStream%random(array)

    call randStream%finalize()
    deallocate(randStream)
  enddo
  call system_clock(c2, cr, cm); dt = dble(c2-c1)/dble(cr)

  print *, 'Total time   (F90_INTRINSIC): ',dt
  print *, 'MegaRNumbers (F90_INTRINSIC): ', 1.0e-6*dble(nstream*length*ntrials)/dt
  print *, 'Summation of Random Numbers: ', SUM(array)
  print *, '--------'; print *, ''

! SIMD-orientated mersenne twister

  call system_clock(c1, cr, cm)
  do m = 1,ntrials
    allocate(randStream, source=ShrDsfmtRandGen(seed, length))

    call randStream%random(array)

    call randStream%finalize()
    deallocate(randStream)
  enddo
  call system_clock(c2, cr, cm); dt = dble(c2-c1)/dble(cr)

  print *, 'Total time   (DSFMT_F03): ',dt
  print *, 'MegaRNumbers (DSFMT_F03): ', 1.0e-6*dble(nstream*length*ntrials)/dt
  print *, 'Summation of Random Numbers: ', SUM(array)
  print *, '--------'; print *, ''

end program test
