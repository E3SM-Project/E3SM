module BeTR_decompMod

  ! !PUBLIC TYPES:
  implicit none
  character(len=*), private, parameter :: mod_filename = &
       __FILE__
  type, public :: betr_bounds_type
    integer :: begg, endg       ! beginning and ending gridcell index
    integer :: begl, endl       ! beginning and ending landunit index
    integer :: begc, endc       ! beginning and ending column index
    integer :: begp, endp       ! beginning and ending pft index
    integer :: lbj, ubj
    integer :: begCohort, endCohort ! beginning and ending cohort indices

    integer :: level            ! whether defined on the proc or clump level
    integer :: clump_index      ! if defined on the clump level, this gives the clump index
  end type betr_bounds_type

end module BeTR_decompMod
