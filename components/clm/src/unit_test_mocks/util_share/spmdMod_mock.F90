module spmdMod
  ! Mock of spmdMod, which just defines masterproc

  implicit none
  save
  private

  logical, parameter, public :: masterproc = .true.
  integer, parameter, public :: iam = 0
end module spmdMod
