module ice_scam

  use ice_kinds_mod

  implicit none
  save

  ! single column control variables (only used for latlon grid)

  logical :: single_column       ! true => single column mode
  real (kind=dbl_kind) scmlat    ! single column latitude (degrees)
  real (kind=dbl_kind) scmlon    ! single column longitude (degrees)

end module ice_scam

	
