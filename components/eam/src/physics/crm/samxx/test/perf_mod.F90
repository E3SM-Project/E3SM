
module perf_mod
  implicit none
  public

contains

  subroutine t_startf(str)
    character(len=*) :: str
  endsubroutine t_startf

  subroutine t_stopf(str)
    character(len=*) :: str
  endsubroutine t_stopf

end module perf_mod

