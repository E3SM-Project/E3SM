module mpassi_mct_vars
  
  use seq_flds_mod
  use mct_mod

  implicit none

  save
  public                               ! By default make data private

  integer :: MPASSI_MCT_ICEID
  type(mct_gsMap), pointer :: MPASSI_MCT_gsMap_i
  type(mct_gGrid), pointer :: MPASSI_MCT_dom_i
 
end module mpassi_mct_vars
