program test

! this program calls the available versions of the random generators

use shr_RandNum_mod, only: ShrIntrinsicRandGen, ShrKissRandGen, &
     ShrF95MtRandGen, ShrDsfmtRandGen
#ifdef INTEL_MKL
use shr_RandNum_mod, only: ShrMklMtRandGen
#endif

INTEGER, parameter :: r8 = selected_real_kind(12)

#ifdef INTEL_MKL
type(ShrMklMtRandGen) :: mkl_gen
#endif
type(ShrKissRandGen) :: kiss_gen
type(ShrF95MtRandGen) :: f95_mt_gen
type(ShrIntrinsicRandGen) :: intrinsic_gen
type(ShrDsfmtRandGen) :: dsfmt_gen

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

  call system_clock(c1, cr, cm)
  do m = 1,ntrials
     mkl_gen = ShrMklMtRandGen(seed)
     call mkl_gen%finalize()
  enddo
  call system_clock(c2, cr, cm); dt1 = dble(c2-c1)/dble(cr)

  mkl_gen = ShrMklMtRandGen(seed)
  call system_clock(c1, cr, cm)
  do m = 1,ntrials
    call mkl_gen%random(array)
  enddo
  call mkl_gen%finalize()
  call system_clock(c2, cr, cm); dt2 = dble(c2-c1)/dble(cr)
  dt = dt1+dt2
  print *, 'Init/term time (SFMT_MKL): ',dt1
  print *, 'Gen  time      (SFMT_MKL): ',dt2
  print *, 'Total time     (SFMT_MKL): ',dt
  print *, 'MegaRNumbers   (SFMT_MKL): ', 1.0e-6*dble(nstream*length*ntrials)/dt
  print *, 'Summation of Random Numbers: ', SUM(array)
  print *, '--------'; print *, ''
#endif

! keep it simple stupid random number

  do n = 1,nstream
     do i = 1, 4
        kiss_seed(n,i) = seed(n)*i+n
     end do
  end do

  call system_clock(c1, cr, cm)
  do m = 1,ntrials
     kiss_gen = ShrKissRandGen(kiss_seed)
     call kiss_gen%finalize()
  enddo
  call system_clock(c2, cr, cm); dt1 = dble(c2-c1)/dble(cr)

  kiss_gen = ShrKissRandGen(kiss_seed)
  call system_clock(c1, cr, cm)
  do m = 1,ntrials
    call kiss_gen%random(array)
  enddo
  call kiss_gen%finalize()
  call system_clock(c2, cr, cm); dt2 = dble(c2-c1)/dble(cr)
  dt = dt1+dt2
  print *, 'Init/term time (KISSVEC): ',dt1
  print *, 'Gen  time      (KISSVEC): ',dt2
  print *, 'Total time     (KISSVEC): ',dt
  print *, 'MegaRNumbers   (KISSVEC): ', 1.0e-6*dble(nstream*length*ntrials)/dt
  print *, 'Summation of Random Numbers: ', SUM(array)
  print *, '--------'; print *, ''

! fortran-95 implementation of merseene twister

  call system_clock(c1, cr, cm)
  do m = 1,ntrials
     f95_mt_gen = ShrF95MtRandGen(seed)
     call f95_mt_gen%finalize()
  enddo
  call system_clock(c2, cr, cm); dt1 = dble(c2-c1)/dble(cr)

  f95_mt_gen = ShrF95MtRandGen(seed)
  call system_clock(c1, cr, cm)
  do m = 1,ntrials
    call f95_mt_gen%random(array)
  enddo
  call f95_mt_gen%finalize()
  call system_clock(c2, cr, cm); dt2 = dble(c2-c1)/dble(cr)
  dt = dt1+dt2
  print *, 'Init/term time (MT19937): ',dt1
  print *, 'Gen  time      (MT19937): ',dt2
  print *, 'Total time     (MT19937): ',dt
  print *, 'MegaRNumbers   (MT19937): ', 1.0e-6*dble(nstream*length*ntrials)/dt
  print *, 'Summation of Random Numbers: ', SUM(array)
  print *, '--------'; print *, ''

! fortran-90 intrinsic pseudorandom number generator

  call random_seed(size=intrinsic_size)
  allocate(intrinsic_seed(nstream,intrinsic_size))
  do n = 1, nstream
     do i = 1, intrinsic_size
        intrinsic_seed(n,i) = seed(n)*i+n
     end do
  end do

  call system_clock(c1, cr, cm)
  do m = 1,ntrials
     intrinsic_gen = ShrIntrinsicRandGen(intrinsic_seed)
     call intrinsic_gen%finalize()
  enddo
  call system_clock(c2, cr, cm); dt1 = dble(c2-c1)/dble(cr)

  intrinsic_gen = ShrIntrinsicRandGen(intrinsic_seed)
  call system_clock(c1, cr, cm)
  do m = 1,ntrials
     call intrinsic_gen%random(array)
  enddo
  call intrinsic_gen%finalize()
  call system_clock(c2, cr, cm); dt2 = dble(c2-c1)/dble(cr)
  dt = dt1+dt2
  print *, 'Init/term time (F90_INTRINSIC): ',dt1
  print *, 'Gen  time      (F90_INTRINSIC): ',dt2
  print *, 'Total time     (F90_INTRINSIC): ',dt
  print *, 'MegaRNumbers   (F90_INTRINSIC): ', 1.0e-6*dble(nstream*length*ntrials)/dt
  print *, 'Summation of Random Numbers: ', SUM(array)
  print *, '--------'; print *, ''

! SIMD-orientated mersenne twister

  call system_clock(c1, cr, cm)
  do m = 1,ntrials
     dsfmt_gen = ShrDsfmtRandGen(seed, length)
     call dsfmt_gen%finalize()
  enddo
  call system_clock(c2, cr, cm); dt1 = dble(c2-c1)/dble(cr)

  dsfmt_gen = ShrDsfmtRandGen(seed, length)
  call system_clock(c1, cr, cm)
  do m = 1,ntrials
    call dsfmt_gen%random(array)
  enddo
  call dsfmt_gen%finalize()
  call system_clock(c2, cr, cm); dt2 = dble(c2-c1)/dble(cr)
  dt = dt1+dt2
  print *, 'Init/term time (DSFMT_F03): ',dt1
  print *, 'Gen  time      (DSFMT_F03): ',dt2
  print *, 'Total time     (DSFMT_F03): ',dt
  print *, 'MegaRNumbers   (DSFMT_F03): ', 1.0e-6*dble(nstream*length*ntrials)/dt
  print *, 'Summation of Random Numbers: ', SUM(array)
  print *, '--------'; print *, ''

end program test
