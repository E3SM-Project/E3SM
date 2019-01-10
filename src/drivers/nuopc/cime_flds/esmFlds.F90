module esmFlds

  use shr_nuopc_fldList_mod , only : shr_nuopc_fldList_type

  implicit none
  public

  !-----------------------------------------------
  ! Component and mapping array indices
  !-----------------------------------------------

  integer, parameter :: ncomps =8
  integer, parameter :: compmed=1
  integer, parameter :: compatm=2
  integer, parameter :: complnd=3
  integer, parameter :: compocn=4
  integer, parameter :: compice=5
  integer, parameter :: comprof=6
  integer, parameter :: compwav=7
  integer, parameter :: compglc=8

  character(len=*),parameter :: compname(ncomps) = (/'med','atm','lnd','ocn','ice','rof','wav','glc'/)

  type (shr_nuopc_fldList_type) :: fldListTo(ncomps) ! advertise fields to components
  type (shr_nuopc_fldList_type) :: fldListFr(ncomps) ! advertise fields from components

  type (shr_nuopc_fldList_type) :: fldListMed_aoflux 
  type (shr_nuopc_fldList_type) :: fldListMed_ocnalb 

end module esmFlds
