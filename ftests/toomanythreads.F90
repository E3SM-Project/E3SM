program toomanythreads
  use gptl

  implicit none

  integer :: mythread
  integer :: ret(2)  ! for each of 2 threads
  integer :: n

  integer, external :: omp_get_thread_num
  
  write(6,*)'Testing setting of maxthreads...'
  if (gptlsetoption (gptlmaxthreads, 1) /= 0) then  ! only allow 1 thread
    write(6,*) 'Failure'
    call exit (1)
  end if
  write(6,*) 'Success'

  if (gptlinitialize () /= 0) then
    write(6,*) 'Failure from gptlinitialize'
    call exit (1)
  end if

  write(6,*) 'Testing using more threads than space was allocated for...'
  call omp_set_num_threads (2)
!$OMP PARALLEL DO PRIVATE (mythread)
  do n=1,2
    mythread = omp_get_thread_num ()
    ret(mythread+1) = gptlstart ('loop1')
  end do
!$OMP END PARALLEL DO

  if (ret(1) == 0 .and. ret(2) == 0) then
    write(6,*) 'Failure: Too many threads did NOT cause a GPTL failure'
    call exit (1)
  end if

  write(6,*) 'Success'
end program toomanythreads
