program errtest
  use gptl

  implicit none

  integer ret, iter
  integer val

  write(6,*)'Purpose: test error conditions'

1 write(6,*)'Enter number for error test:'
  write(6,*)'0: bad option'
  write(6,*)'3: stop never started'
  write(6,*)'4: stop while already stopped'
  write(6,*)'5: instance not called'

  read (5,*) val
  if (val < 0 .or. val > 5) then
    write(6,*)'Val must be between 0 and 5'
    goto 1
  end if

  if (val == 0) then
    if (gptlsetoption (100, 1) < 0) write(6,*)'setoption failure'
    if (gptlinitialize () < 0) write(6,*)'initialize failure'
  else if (val == 3) then
    if (gptlinitialize () < 0) write(6,*)'initialize failure'
    if (gptlstop ('errtest') < 0) write(6,*)'stop failure'
  else if (val == 4) then
    if (gptlinitialize () < 0) write(6,*)'initialize failure'
    if (gptlstart ('errtest') < 0) write(6,*)'start failure'
    if (gptlstop ('errtest') < 0) write(6,*)'stop failure'
    if (gptlstop ('errtest') < 0) write(6,*)'stop failure'
  else if (val == 5) then
    if (gptlstart ('errtest') < 0) write(6,*)'start failure'
    if (gptlstop ('errtest') < 0) write(6,*)'stop failure'
    if (gptlpr (0) < 0) write(6,*)'stop failure'
  end if
  
  if (gptlfinalize () < 0) write(6,*)'gptlfinalize error'
  stop 0
end program errtest
