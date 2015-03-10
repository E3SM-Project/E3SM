program testbacktrace
  use gptl

  implicit none

  integer :: n, ret

  ret = gptlsetoption (gptldopr_memusage, 1)
  ret = gptlinitialize ()

!$OMP PARALLEL DO
  do n=1,240
    call threaded_sub (n)
  end do
end program testbacktrace

subroutine threaded_sub (n)
  implicit none

  integer, intent(in) :: n

  call onemstack (n)
  return
end subroutine threaded_sub

subroutine onemstack (n)
  implicit none

  integer, intent(in) :: n
  character(len=1) :: chararr(1000000)
  
  chararr(:) = 'x'
  call testchars (chararr)
  return
end subroutine onemstack

subroutine testchars (chararr)
  implicit none

  character(len=1), intent(in) :: chararr(1000000)

  if (chararr(1000000) == 'y') then
    write(6,*)'charr(1000000) = y'
  end if
  return
end subroutine testchars
