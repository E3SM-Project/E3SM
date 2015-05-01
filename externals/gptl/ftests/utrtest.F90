program utrtest
  use gptl

  implicit none

  external :: sub

  double precision :: sum
  integer :: ret
  integer :: handle1
  integer :: handle2
  integer :: handle3
  integer :: handle4
  integer :: handle5
  integer :: handle6
  integer :: handle7
  integer :: handle8

  sum = 0.

  write(6,*) 'Purpose: estimate overhead of GPTL timing (UTR)'
  ret = gptlsetoption (gptlabort_on_error, 0)
  ret = gptlsetoption (gptlverbose, 1)
!  ret = gptlsetoption (gptltablesize, 111)
  
!  ret = gptlsetutr (gptlmpiwtime)
  ret = gptlsetutr (gptlread_real_time)
  ret = gptlsetutr (gptlclockgettime)
  ret = gptlsetutr (gptlgettimeofday)
  ret = gptlsetutr (gptlpapitime)
  ret = gptlsetutr (gptlnanotime)
  
  ret = gptlinitialize ()

  ret = gptlinit_handle ('1x1e7', handle1)
  ret = gptlinit_handle ('10x1e6', handle2)
  ret = gptlinit_handle ('100x1e5', handle3)
  ret = gptlinit_handle ('1000x1e4', handle4)
  ret = gptlinit_handle ('1e4x1000', handle5)
  ret = gptlinit_handle ('1e5x100', handle6)
  ret = gptlinit_handle ('1e6x10', handle7)
  ret = gptlinit_handle ('1e7x1', handle8)
  
  ret = gptlstart ('total')
  !      ret = GPTLdisable ()
  call sub (1, 10000000, "1x1e7", sum, handle1)
  call sub (10, 1000000, "10x1e6", sum, handle2)
  call sub (100, 100000, "100x1e5", sum, handle3)
  call sub (1000, 10000, "1000x1e4", sum, handle4)
  call sub (10000, 1000, "1e4x1000", sum, handle5)
  call sub (100000, 100, "1e5x100", sum, handle6)
  call sub (1000000, 10, "1e6x10", sum, handle7)
  call sub (10000000, 1, "1e7x1", sum, handle8)
  !      ret = gptlenable ()
  ret = gptlstop ("total")

  ret = gptlpr (0)
  stop 0
end program utrtest

subroutine sub (outer, inner, name, sum, handle)
  use gptl

  implicit none

  integer, intent(in) :: outer
  integer, intent(in) :: inner
  character(len=*), intent(in) :: name
  double precision, intent(inout) :: sum
  integer, intent(inout) :: handle
  
  integer :: i, j, ret

  do i=0,outer-1
    ret = gptlstart_handle (name, handle)
    do j=0,inner-1
      sum = sum + j
    end do
    ret = gptlstop_handle (name, handle)
  end do
  
  return
end subroutine sub
