program outoforder
! Purpose: test behavior of imperfectly-nested regions
  use gptl
  implicit none

  integer :: n
  integer :: ret
  integer :: num_warn
  integer :: kount

  write(6,*) 'Testing out of order calls...'
  ret = gptlinitialize ()

  do n=1,10
    ret = gptlstart ("xxx")
    ret = gptlstart ("yyy")
    ret = gptlstart ("zzz")
    ret = gptlstop ("xxx")
    ret = gptlstop ("yyy")
    ret = gptlstop ("zzz")
  end do

  ret = gptlstart ("A")
  ret = gptlstart ("B")
  ret = gptlstart ("C")
  ret = gptlstop ("C")
  ret = gptlstop ("B")
  ret = gptlstop ("A")

  ret = gptlget_count ('xxx', 0, kount)
  if (kount /= 10) then
    write(6,*)'Failure: Got count=', kount, ' for xxx when expected 10'
    call exit (1)
  end if

  ret = gptlget_count ('yyy', 0, kount)
  if (kount /= 10) then
    write(6,*)'Failure: Got count=', kount, ' for yyy when expected 10'
    call exit (1)
  end if

  ret = gptlget_count ('zzz', 0, kount)
  if (kount /= 10) then
    write(6,*)'Failure: Got count=', kount, ' for zzz when expected 10'
    call exit (1)
  end if

  num_warn = gptlnum_warn ()
  if (num_warn > 0) then
    write(6,*)'Success: ', num_warn,' warnings were found'
  else
    write(6,*) 'Failure: no warnings were found when they should have been'
    call exit (1)
  end if
  ret = gptlpr (0)
end program outoforder
