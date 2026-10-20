module ColumnDataType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Column data type allocation and initialization
  ! --------------------------------------------------------
  !
  use ColumnType      , only : col_pp
  use ColumnStateTypes, only : column_energy_state, column_water_state, column_carbon_state, &
      column_nitrogen_state, column_phosphorus_state
  use ColumnFluxTypes, only : column_energy_flux, column_water_flux, column_carbon_flux, &
      column_nitrogen_flux, column_phosphorus_flux
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  public

  !-----------------------------------------------------------------------
  ! declare the public instances of column-level data types
  !-----------------------------------------------------------------------
  ! State types
  type(column_energy_state)          , public, target :: col_es     ! column energy state
  type(column_water_state)           , public, target :: col_ws     ! column water state
  type(column_carbon_state)          , public, target :: col_cs     ! column carbon state
  type(column_carbon_state)          , public, target :: c13_col_cs ! column carbon state (C13)
  type(column_carbon_state)          , public, target :: c14_col_cs ! column carbon state (C14)
  type(column_nitrogen_state)        , public, target :: col_ns     ! column nitrogen state
  type(column_phosphorus_state)      , public, target :: col_ps     ! column phosphorus state
  ! Flux types
  type(column_energy_flux)           , public, target :: col_ef     ! column energy flux
  type(column_water_flux)            , public, target :: col_wf     ! column water flux
  type(column_carbon_flux)           , public, target :: col_cf     ! column carbon flux
  type(column_carbon_flux)           , public, target :: c13_col_cf ! column carbon flux
  type(column_carbon_flux)           , public, target :: c14_col_cf ! column carbon flux
  type(column_nitrogen_flux)         , public, target :: col_nf     ! column nitrogen flux
  type(column_phosphorus_flux)       , public, target :: col_pf     ! column phosphorus flux

  !$acc declare create(col_es)
  !$acc declare create(col_ef)
  !$acc declare create(col_ws)
  !$acc declare create(col_wf)
  !$acc declare create(col_cs)
  !$acc declare create(c13_col_cs)
  !$acc declare create(c14_col_cs)
  !$acc declare create(col_ns    )
  !$acc declare create(col_ps    )

  !$acc declare create(col_cf    )
  !$acc declare create(c13_col_cf)
  !$acc declare create(c14_col_cf)
  !$acc declare create(col_nf    )
  !$acc declare create(col_pf    )
  !------------------------------------------------------------------------

    !------------------------------------------------------------------------

end module ColumnDataType
