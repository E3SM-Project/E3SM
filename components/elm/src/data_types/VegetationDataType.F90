module VegetationDataType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Vegetation data type allocation and initialization
  ! --------------------------------------------------------
  !
  use VegetationStateTypes, only : vegetation_energy_state, vegetation_water_state, &
      vegetation_carbon_state, vegetation_nitrogen_state, vegetation_phosphorus_state
  use VegetationFluxTypes, only : vegetation_energy_flux, vegetation_water_flux, &
      vegetation_carbon_flux, vegetation_nitrogen_flux, vegetation_phosphorus_flux
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  public

  !-----------------------------------------------------------------------
  ! declare the public instances of vegetation-level data types
  !-----------------------------------------------------------------------
  type(vegetation_energy_state)          , public, target :: veg_es     ! vegetation energy state
  type(vegetation_water_state)           , public, target :: veg_ws     ! vegetation water state
  type(vegetation_carbon_state)          , public, target :: veg_cs     ! vegetation carbon state
  type(vegetation_carbon_state)          , public, target :: c13_veg_cs ! vegetation carbon state (C13)
  type(vegetation_carbon_state)          , public, target :: c14_veg_cs ! vegetation carbon state (C14)
  type(vegetation_nitrogen_state)        , public, target :: veg_ns     ! vegetation nitrogen state
  type(vegetation_phosphorus_state)      , public, target :: veg_ps     ! vegetation phosphorus state
  type(vegetation_energy_flux)           , public, target :: veg_ef     ! vegetation energy flux
  type(vegetation_water_flux)            , public, target :: veg_wf     ! vegetation water flux
  type(vegetation_carbon_flux)           , public, target :: veg_cf     ! vegetation carbon flux
  type(vegetation_carbon_flux)           , public, target :: c13_veg_cf ! vegetation carbon flux (C13)
  type(vegetation_carbon_flux)           , public, target :: c14_veg_cf ! vegetation carbon flux (C14)
  type(vegetation_nitrogen_flux)         , public, target :: veg_nf     ! vegetation nitrogen flux
  type(vegetation_phosphorus_flux)       , public, target :: veg_pf     ! vegetation phosphorus flux

   !$acc declare create(veg_es    )
   !$acc declare create(veg_ws    )
   !$acc declare create(veg_cs    )
   !$acc declare create(c13_veg_cs)
   !$acc declare create(c14_veg_cs)
   !$acc declare create(veg_ns    )
   !$acc declare create(veg_ps    )
   !$acc declare create(veg_ef    )
   !$acc declare create(veg_wf    )
   !$acc declare create(veg_cf    )
   !$acc declare create(c13_veg_cf)
   !$acc declare create(c14_veg_cf)
   !$acc declare create(veg_nf    )
   !$acc declare create(veg_pf    )

  !------------------------------------------------------------------------

    !------------------------------------------------------------------------

end module vegetationDataType
