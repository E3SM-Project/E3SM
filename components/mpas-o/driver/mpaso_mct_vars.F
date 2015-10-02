module mpaso_mct_vars
  
  use seq_flds_mod
  use mct_mod

  implicit none

  save
  public                               ! By default make data private

  integer :: MPASO_MCT_OCNID
  type(mct_gsMap), pointer :: MPASO_MCT_gsMap_o
  type(mct_gGrid), pointer :: MPASO_MCT_dom_o
 
end module mpaso_mct_vars
