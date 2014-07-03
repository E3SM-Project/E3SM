module clm_varsur

  !-----------------------------------------------------------------------
  ! Module containing 2-d surface boundary data information
  ! surface boundary data, these are all "gdc" local 
  ! Note that some of these need to be pointers (as opposed to just allocatable arrays) to
  ! match the ncd_io interface; for consistency, we make them all pointers
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! weight of each landunit on the grid cell
  real(r8), pointer :: wt_lunit(:,:)     

  ! whether we have valid urban data in each grid cell
  logical , pointer :: urban_valid(:)

  ! for natural veg landunit, weight of each pft on the landunit (adds to 1.0 on the
  ! landunit for all all grid cells, even! those without any natural pft)
  ! (second dimension goes natpft_lb:natpft_ub)
  real(r8), pointer :: wt_nat_pft(:,:)   

  ! for crop landunit, weight of each cft on the landunit (adds to 1.0 on the
  ! landunit for all all grid cells, even  those without any crop)
  ! (second dimension goes cft_lb:cft_ub)
  real(r8), pointer :: wt_cft(:,:)       

  ! for glc_mec landunits, weight of glacier in each elevation class (adds to 1.0 on the
  ! landunit for all grid cells, even those without any glacier)
  real(r8), pointer :: wt_glc_mec(:,:)   

  ! subgrid glacier_mec sfc elevation
  real(r8), pointer :: topo_glc_mec(:,:) 
  !-----------------------------------------------------------------------

end module clm_varsur
