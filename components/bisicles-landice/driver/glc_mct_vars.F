module glc_mct_vars
  
  use seq_flds_mod
  use mct_mod

  implicit none

  save
  public

  integer :: MALI_MCT_GLCID
  type(mct_gsMap), pointer :: MALI_MCT_gsMap_g
  type(mct_gGrid), pointer :: MALI_MCT_dom_g
 
end module glc_mct_vars

