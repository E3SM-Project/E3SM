module perf_mod

implicit none
public

contains

  subroutine t_startf(timer_name)
    character(len=*) timer_name
    return
  end subroutine
  
  subroutine t_stopf(timer_name)
    character(len=*) timer_name
    return
  end subroutine

end module perf_mod
