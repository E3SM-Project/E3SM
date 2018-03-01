module betr_ctrl
!
!logical switches for global control of betr functions
implicit none
  character(len=*), private, parameter :: mod_filename = &
       __FILE__
  logical, public :: use_pH_data = .false.
  logical, public :: betr_use_cn =.false.
  integer, public :: biulog = 6        ! "stdout" log file unit number, default is 6
  logical, public :: do_betr_output = .true.
  integer, parameter, public :: max_betr_hist_type=400
  integer, parameter, public :: max_betr_rest_type=400
  logical, public :: betr_offline = .true.
end module betr_ctrl
