module perfMod_GPU
  !This module will include intermediate subroutines
  !to guarantee compatibility with perfomance modules
  !within OpenACC code.
  !
  use perf_mod  , only : t_startf, t_stopf

contains

  subroutine  t_startGPU(event)

    character(len=64), intent(in) :: event

#ifndef _OPENACC
      call t_startf(trim(event))
#endif

  end subroutine t_startGPU

  subroutine  t_stopGPU(event)

    character(len=64), intent(in) :: event

#ifndef _OPENACC
      call t_stopf(trim(event))
#endif

end subroutine t_stopGPU


end module perfMod_GPU
