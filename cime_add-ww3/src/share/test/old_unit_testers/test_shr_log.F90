program test_shr_log
  use test_mod, only : test_init, test_final
  implicit none

  call test_init

  call test_shr_log_errMsg

  call test_final

contains

  subroutine test_shr_log_errMsg
    use shr_log_mod
    use test_mod

    implicit none

    character(len=256) :: my_result

    my_result = shr_log_errMsg('myfile.f90', 42)

    call test_is(my_result, "ERROR in myfile.f90 at line 42", "shr_log_errMsg: basic test")

  end subroutine test_shr_log_errMsg
end program test_shr_log


