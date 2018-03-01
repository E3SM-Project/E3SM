module mpas_cice_mct_vars
  
  use seq_flds_mod
  use mct_mod

  implicit none

  save
  public                               ! By default make data private

  integer :: MPAS_CICE_MCT_ICEID
  type(mct_gsMap), pointer :: MPAS_CICE_MCT_gsMap_i
  type(mct_gGrid), pointer :: MPAS_CICE_MCT_dom_i
 
end module mpas_cice_mct_vars
