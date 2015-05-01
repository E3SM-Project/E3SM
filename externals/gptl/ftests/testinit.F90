program testinit
  use gptl

  implicit none

  integer :: ret

  write(6,*)'testinit: Testing gptl_papilibraryinit...'
  ret = gptl_papilibraryinit ()

  if (ret == 0) then
    write(6,*)'Success'
  else
    write(6,*)'Failure'
    stop 999
  end if
end program testinit
