module control_mod
  use kinds, only : real_kind
  integer, parameter :: isrf_forc=0
  integer :: np, ntime, nplot
  real(kind=real_kind):: press_toa, srfpress, dt
contains
  subroutine readnl
    namelist /physicsutest_nl/ np, ntime, nplot, &
         press_toa,srfpress,dt
    
    open(unit=7,file='input.nl',status='old')
    read(unit=7,nml=physicsutest_nl)
    close(7)
          

  end subroutine readnl

end module control_mod
