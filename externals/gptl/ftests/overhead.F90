program overhead
  use gptl

  implicit none

  integer :: ret, iter, i
  real*8 :: wall1, usr1, sys1
  real*8 :: wall2, usr2, sys2
  integer :: cycles1, cycles2, cps
  character(len=16) :: msg
  integer :: handle
#ifdef THREADED_OMP
  integer, parameter :: maxthreads = 8
#else
  integer, parameter :: maxthreads = 1
#endif
  integer :: nthreads

#ifdef THREADED_OMP
  integer, external :: omp_set_num_threads
#endif

  open (unit=1, file='gettimeofday', form='formatted', status='replace')
  nthreads = 1
  do while (nthreads <= maxthreads)
    write(6,*)'overhead: running gettimeofday nthreads=', nthreads
    ret = gptlsetutr (gptlgettimeofday)
#ifdef THREADED_OMP
    ret = omp_set_num_threads (nthreads)
#endif
    ret = gptlsetoption (gptlabort_on_error, 1)
    ret = gptlinitialize ()
    ret = gptlstamp (wall1, usr1, sys1)
!$OMP PARALLEL DO PRIVATE (RET)
    do i=1,10000000
      ret = gptlstart ('loop')
      ret = gptlstop ('loop')
    end do
    ret = gptlstamp (wall2, usr2, sys2)
    write(1,'(i3,3f9.3)') nthreads, wall2-wall1, usr2-usr2, sys2-sys1
    nthreads = nthreads * 2
    ret = gptlpr (0)
    ret = gptlfinalize ()
  end do
  close (unit=1)

#ifdef _AIX
  open (unit=1, file='read_real_time', form='formatted', status='replace')
#else
  open (unit=1, file='nanotime', form='formatted', status='replace')
#endif
  nthreads = 1
  do while (nthreads <= maxthreads)
#ifdef _AIX
    write(6,*)'overhead: running read_real_time nthreads=', nthreads
    ret = gptlsetutr (gptlread_real_time)
#else
    write(6,*)'overhead: running nanotime nthreads=', nthreads
    ret = gptlsetutr (gptlnanotime)
#endif

#ifdef THREADED_OMP
    ret = omp_set_num_threads (nthreads)
#endif
    ret = gptlsetoption (gptlabort_on_error, 1)
    ret = gptlinitialize ()
    ret = gptlstamp (wall1, usr1, sys1)
!$OMP PARALLEL DO PRIVATE (RET)
    do i=1,10000000
      ret = gptlstart ('loop')
      ret = gptlstop ('loop')
    end do
    ret = gptlstamp (wall2, usr2, sys2)
    write(1,'(i3,3f9.3)') nthreads, wall2-wall1, usr2-usr2, sys2-sys1
    nthreads = nthreads * 2
    ret = gptlpr (0)
    ret = gptlfinalize ()
  end do
  close (unit=1)

  open (unit=1, file='no_wallclock', form='formatted', status='replace')
  nthreads = 1
  do while (nthreads <= maxthreads)
    write(6,*)'overhead: running no_wallclock nthreads=', nthreads
    ret = gptlsetutr (gptlgettimeofday)
#ifdef THREADED_OMP
    ret = omp_set_num_threads (nthreads)
#endif
    ret = gptlsetoption (gptlwall, 0)
    ret = gptlsetoption (gptlabort_on_error, 1)
    ret = gptlinitialize ()
    call system_clock (cycles1, cps)
!$OMP PARALLEL DO PRIVATE (RET)
    do i=1,10000000
      ret = gptlstart ('loop')
      ret = gptlstop ('loop')
    end do
    call system_clock (cycles2, cps)
    write(1,'(i3,f9.3)') nthreads, float((cycles2-cycles1))/cps
    nthreads = nthreads * 2
    ret = gptlpr (0)
    ret = gptlfinalize ()
  end do
  close (unit=1)

  open (unit=1, file='handle', form='formatted', status='replace')
  nthreads = 1
  do while (nthreads <= maxthreads)
    write(6,*)'overhead: running handle (no_wallclock) nthreads=', nthreads
    ret = gptlsetutr (gptlgettimeofday)
#ifdef THREADED_OMP
    ret = omp_set_num_threads (nthreads)
#endif
    ret = gptlsetoption (gptlwall, 0)
    ret = gptlsetoption (gptlabort_on_error, 1)
    ret = gptlinitialize ()
    handle = 0
    call system_clock (cycles1, cps)
!$OMP PARALLEL DO PRIVATE (RET) FIRSTPRIVATE (HANDLE)
    do i=1,10000000
      ret = gptlstart_handle ('loop', handle)
      ret = gptlstop_handle ('loop', handle)
    end do
    call system_clock (cycles2, cps)
    write(1,'(i3,f9.3)') nthreads, float((cycles2-cycles1))/cps
    nthreads = nthreads * 2
    ret = gptlpr (0)
    ret = gptlfinalize ()
  end do
  close (unit=1)

  open (unit=1, file='do_nothing', form='formatted', status='replace')
  nthreads = 1
  do while (nthreads <= maxthreads)
    write(6,*)'overhead: running do_nothing nthreads=', nthreads
    handle = 0
    call system_clock (cycles1, cps)
!$OMP PARALLEL DO PRIVATE (RET) FIRSTPRIVATE (HANDLE)
    do i=1,10000000
      call do_nothing1 ('string1', handle)
      call do_nothing2 ('string2', handle)
      call do_nothing1 ('string3', handle)
      call do_nothing2 ('string4', handle)
    end do
    call system_clock (cycles2, cps)
    write(1,'(i3,f9.3)') nthreads, float((cycles2-cycles1))/cps
    nthreads = nthreads * 2
  end do
  close (unit=1)

  stop 0
end program overhead

subroutine do_nothing1 (string, handle)
  implicit none

  character(len=*), intent(in) :: string
  integer, intent(in) :: handle

  if (string(1:1) == 'x') then
    write(6,*)'Bad string value'
  end if
end subroutine do_nothing1

subroutine do_nothing2 (string, handle)
  implicit none

  character(len=*), intent(in) :: string
  integer, intent(in) :: handle

  if (string(1:1) == 'x') then
    write(6,*)'Bad string value'
  end if
end subroutine do_nothing2
