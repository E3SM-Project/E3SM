program nlreader
  use gptl

  implicit none

  integer :: ret

  ret = gptlsetoption (gptlverbose, 1)
  write(6,*)'nlreader: Testing gptlprocess_namelist...'
  call gptlprocess_namelist ('gptlnl', 1, ret)
  if (ret /= 0) then
    write(6,*)'Failure'
    call exit (1)
  end if
      
  ! Now turn off verbosity

  ret = gptlsetoption (gptlverbose, 0)
  ret = gptlinitialize ()
  ret = gptlstart ('main')
  ret = gptlstop ('main')
  ret = gptlpr (0)
  write(6,*)'Success'
end program nlreader
