program handle
  use gptl
  implicit none

  integer :: handle1       ! Hash index
  integer :: n
  integer :: ret

  ret = gptlinitialize ()

  ret = gptlstart ('total') ! Time the entire code
! IMPORTANT: Start with a zero handle value so GPTLstart_handle knows to initialize
  handle1 = 0

!$OMP PARALLEL DO SHARED (handle1)
  do n=1,1000000
! First call the "_handle" versions of start and stop for the region
    ret = gptlstart_handle ('loop', handle1)
    ret = gptlstop_handle ('loop', handle1)
! Now call the standard start and stop functions for the same region
    ret = gptlstart ('loop')
    ret = gptlstop ('loop')
  end do
  ret = gptlstop ('total') ! Time the entire code

  ret = gptlpr (0)
  stop
end program handle
