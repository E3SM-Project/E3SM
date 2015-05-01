program testpapi
  use gptl

  implicit none

  integer :: ret
  integer :: i, code
  integer(8) :: pc ! papi counters
  real(8) :: sum, val
  character(len=256) :: name
  character(len=2) :: tooshort
  
  write(6,*)'testpapi: Testing PAPI interface...'

  write(6,*)'Testing enabling gptlverbose...'
  if (gptlsetoption (gptlverbose, 1) /= 0) then
    write(6,*)'Failure'
    call exit(1)
  end if
  write(6,*)'Success'

  write(6,*)'Calling gptl_papilibraryinit...'
  if (gptl_papilibraryinit () /= 0) then
    write(6,*)'Failure from gptl_papilibraryinit...'
    call exit(1)
  end if
  write(6,*)'Success'
      
  write(6,*)'Testing gptlevent_name_to_code for PAPI_TOT_CYC...'
  if (gptlevent_name_to_code ('PAPI_TOT_CYC', code) /= 0) then
    write(6,*)'Failure from gptlevent_name_to_code'
    call exit(1)
  end if
  write(6,*)'Success: PAPI_TOT_CYC=',code

  write(6,*)'Testing passing PAPI_TOT_CYC to gptlsetoption...'
  if (gptlsetoption (code, 1) /= 0) then
    write(6,*)'Failure'
    call exit(1)
  end if
  write(6,*)'Success'
  
  write(6,*)'Testing duplicate enable PAPI_TOT_CYC...'
  if (gptlsetoption (code, 1) == 0) then
    write(6,*)'Failure to fail!'
    call exit(1)
  end if
  write(6,*)'Succeeded at failing!'
  
  write(6,*)'Testing turning off an already-on counter...'
  if (gptlsetoption (code, 0) == 0) then
    write(6,*)'Failure'
    call exit(1)
  end if
  write(6,*)'Succeeded at failing!'
  
  write(6,*)'Testing gptlevent_code_to_name for PAPI_TOT_CYC...'
  if (gptlevent_code_to_name (code, name) /= 0) then
    write(6,*)'Failure from gptlevent_code_to_name'
    call exit(1)
  end if

  if (trim(name) == 'PAPI_TOT_CYC') then
    write(6,*)'Success'
  else
    write(6,*)'Failure: got ',trim(name)
    write(6,*)'Expected PAPI_TOT_CYC'
    call exit(1)
  end if
  
  write(6,*)'Testing gptlevent_name_to_code for GPTL_CI...'
  if (gptlevent_name_to_code ('GPTL_CI', code) /= 0) then
    write(6,*)'Failure from gptlevent_name_to_code'
    call exit(1)
  end if
  write(6,*)'Success: GPTL_CI=',code
  
  write(6,*)'Testing too short var for gptlevent_code_to_name...'
  if (gptlevent_code_to_name (code, tooshort) == 0) then
    write(6,*)'Failure of gptlevent_code_to_name to fail'
    call exit(1)
  end if
  write(6,*)'Success at catching too short output var name'
  
  write(6,*)'Testing gptlevent_code_to_name for GPTL_CI...'
  if (gptlevent_code_to_name (code, name) /= 0) then
    write(6,*)'Failure from gptlevent_code_to_name'
    call exit(1)
  end if

  if (name == 'GPTL_CI') then
    write(6,*)'Success'
  else
    write(6,*)'Failure: got ',trim(name)
    write(6,*)'Expected GPTL_CI'
    call exit(1)
  end if

  write(6,*)'Testing bogus input to gptlevent_name_to_code...'
  if (gptlevent_name_to_code ('zzz', code) == 0) then
    write(6,*)'Failure of gptlevent_name_to_code to fail'
    call exit(1)
  end if
  write(6,*)'Success at catching bogus input name'
  
  write(6,*)'Testing bogus input to gptlevent_code_to_name...'
  code = -1
  if (gptlevent_code_to_name (code, name) == 0) then
    write(6,*)'Failure of gptlevent_code_to_name to fail'
    call exit(1)
  end if
  write(6,*)'Success at catching bogus input code'
  
  write(6,*)'Testing gptlinitialize'
  if (gptlinitialize () /= 0) then
    write(6,*)'Failure'
    call exit(1)
  end if
  write(6,*)'Success'
  
  ret = gptlstart ('sum')
  sum = 0.
  do i=1,1000000
    sum = sum + i
  end do
  ret = gptlstop ('sum')
  
  write(6,*)'Testing gptlquerycounters...'
  if (gptlquerycounters ('sum', 0, pc) /= 0) then
    write(6,*)'Failure'
    call exit(1)
  end if
  write(6,*)'Success: pc=', pc
  
  write(6,*)'Testing gptlget_eventvalue...'
  if (gptlget_eventvalue ('sum', 'PAPI_TOT_CYC', 0, val) /= 0) then
    write(6,*)'Failure'
    call exit(1)
  end if
  write(6,*)'Success: val=', val
  
  write(6,*)'sum,pc=',sum, pc
  if (pc < 1 .or. pc > 1.e9) then
    write(6,*)'Suspicious pc value=',pc
    call exit(1)
  else
    write(6,*)'Success'
  end if  
end program testpapi
