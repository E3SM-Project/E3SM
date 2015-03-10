program testbasics
  use gptl

  implicit none

  integer :: ret
  integer :: nregions
  integer count, of       ! for gptlquery
  integer(8) :: pc        ! for gptlquery
  integer(8), allocatable :: iarr8(:)
  real(8) :: wc, usr, sys ! for gptlquery
  character(len=80) :: str
  integer :: size, rss, share, text, datastack
  
  write(6,*)'testbasics: Testing basic GPTL usage...'
  write(6,*)'Testing gptlinitialize...'
  if (gptlinitialize () /= 0) then
    write(6,*)'Failure in gptlinitialize'
    call exit (1)
  end if
  write(6,*)'Success'
  
  write(6,*)'Testing gptlstart...'
  if (gptlstart ('testbasics') /= 0) then
    write(6,*)'Failure in gptlstart'
    call exit(1)
  end if
  write(6,*)'Success'
  
  write(6,*)'Testing gptlstop...'
  if (gptlstop ('testbasics') /= 0) then
    write(6,*)'Failure in gptlstop'
    call exit(1)
  end if
  write(6,*)'Success'
  
  write(6,*)'Testing gptlpr...'
  if (gptlpr (0) /= 0) then
    write(6,*)'Failure in gptlpr(0)'
    call exit(1)
  end if
  write(6,*)'Success'
  
  write(6,*)'Testing gptlget_wallclock...'
  if (gptlget_wallclock ('testbasics', 0, wc) /= 0) then
    write(6,*)'Failure in gptlget_wallclock'
    call exit(1)
  end if
  write(6,*)'Success: wc=', wc
  
  write(6,*)'Testing gptlget_memusage...'
  if (gptlget_memusage (size, rss, share, text, datastack) /= 0) then
    write(6,*)'Failure in gptlget_memusage'
    call exit(1)
  end if
  write(6,*)'Success: size=', size, ' rss=', rss

  write(6,*)'Testing gptlprint_memusage...'
  str = 'testbasics before allocating 100 MB'
  if (gptlprint_memusage (trim(str)) /= 0) then
    write(6,*)'Failure in gptlprint_memusage'
    call exit(1)
  end if
  
  allocate (iarr8(13107200))
  iarr8(:) = 0
  
  str = 'testbasics after allocating 100 MB'
  if (gptlprint_memusage (trim(str)) /= 0) then
    write(6,*)'Failure in gptlprint_memusage'
    call exit(1)
  end if
  
  if (gptlprint_rusage (trim(str)) /= 0) then
    write(6,*)'Failure in gptlprint_rusage'
    call exit(1)
  end if
  
  write(6,*)'Success'
  
  write(6,*)'Testing gptldisable/gptlenable...'
  if (gptldisable () /= 0) then
    write(6,*)'Failure in gptldisable'
    call exit(1)
  end if
  if (gptlstart ('zzz') /= 0) then
    write(6,*)'Failure in disabled gptlstart'
    call exit(1)
  end if
  if (gptlstop ('zzz') /= 0) then
    write(6,*)'Failure in disabled gptlstop'
    call exit(1)
  end if
  if (gptlenable () /= 0) then
    write(6,*)'Failure in gptlenable'
    call exit(1)
  end if

  write(6,*)'Sub-testing gptlget_nregions...'
  if (gptlget_nregions (0, nregions) /= 0) then
    write(6,*)'Failure in gptlget_nregions'
    call exit(1)
  end if
  if (nregions /= 1) then
    write(6,*)'Failure: expected nregions=1 got', nregions
    call exit(1)
  end if
  write(6,*)'Success in gptlget_nregions'
  write(6,*)'Success'
  
  write(6,*)'Testing gptlquery...'
  ret = gptlquery ('testbasics', 0, count, of, wc, usr, sys, pc, 0)
  if (ret /= 0) then
    write(6,*)'Failure'
    call exit(1)
  end if
  if (count /= 1) then
    write(6,*)'Bad count value from gptlquery'
    call exit(1)
  end if
  if (of /= 0) then
    write(6,*)'Bad onflg value from gptlquery'
    call exit(1)
  end if
  write (6,*)'Success'
  
  write(6,*)'Testing gptlreset...'
  ret = gptlreset ()
  if (ret /= 0) then
    write(6,*)'Failure'
    call exit(1)
  end if
  ret = gptlquery ('testbasics', 0, count, of, wc, usr, sys, pc, 0)
  if (ret /= 0) then
    write(6,*)'Failure from gptlquery'
    call exit(1)
  end if
  if (count/=0 .or. of/=0 .or. wc/=0. .or. usr/=0.) then
    write(6,*)'Failure: one or more counts were not zeroed'
    call exit(1)
  end if
  write (6,*)'Success'
  
  write (6,*)'testbasics: All tests succeeded'
end program testbasics
