module FatesInterfaceMod

  ! ------------------------------------------------------------------------------------
  ! STUB STUB STUB
  
  ! This is the FATES public API
  ! A host land model has defined and allocated a structure "fates" as
  ! defined by fates_interface_type
  !
  ! It is also likely/possible that this type is defined as a vector
  ! which is allocated by thread
  ! ------------------------------------------------------------------------------------
  
  ! -------------------------------------------------------------------------------------
  ! Parameters that are dictated by FATES and known to be required knowledge
  !  needed by the HLMs
  ! -------------------------------------------------------------------------------------
  
  ! Variables mostly used for dimensioning host land model (HLM) array spaces

  implicit none

  public :: set_fates_global_elements
  
  integer, protected :: fates_maxElementsPerPatch
  integer, protected :: fates_maxElementsPerSite


  contains

  subroutine set_fates_global_elements(use_fates)
    implicit none
    
    logical,intent(in) :: use_fates    ! Is fates turned on?
    
    if (use_fates) then
       
       fates_maxElementsPerPatch = 1
       fates_maxElementsPerSite = 1
       
    else
       ! If we are not using FATES, the cohort dimension is still
       ! going to be initialized, lets set it to the smallest value
       ! possible so that the dimensioning info takes up little space
       
       fates_maxElementsPerPatch = 1
       
       fates_maxElementsPerSite = 1
         

      end if
  end module FatesInterfaceMod
