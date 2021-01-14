module scream_cpl_indices

  use iso_c_binding, only: c_int, c_char

  implicit none
  private

  ! Focus only on the ones that scream imports/exports (subsets of x2a and a2x)
  integer, parameter, public :: num_required_imports = 21
  integer, parameter, public :: num_required_exports = 12
  integer, parameter, public :: num_optional_imports = 0
  integer, parameter, public :: num_optional_exports = 1
  integer, parameter, public :: num_imports = num_required_imports + num_optional_imports
  integer, parameter, public :: num_exports = num_required_exports + num_optional_exports

  integer(kind=c_int), public, allocatable, target :: index_x2a(:)
  integer(kind=c_int), public, allocatable, target :: index_a2x(:)

  ! Names used by the component coupler for import/export fields
  character(len=32,kind=c_char), public, allocatable, target :: cpl_names_a2x(:)
  character(len=32,kind=c_char), public, allocatable, target :: cpl_names_x2a(:)

  ! Names used by scream for import/export fields
  character(len=32,kind=c_char), public, allocatable, target :: scr_names_a2x(:)
  character(len=32,kind=c_char), public, allocatable, target :: scr_names_x2a(:)

  public :: scream_set_cpl_indices

contains

  subroutine scream_set_cpl_indices (x2a, a2x)
    use iso_c_binding,  only: C_NULL_CHAR
    use mct_mod,        only: mct_aVect, mct_avect_indexra
    !
    ! Input(s)
    !
    type(mct_avect), intent(in) :: x2a, a2x
    !
    ! Local(s)
    !
    integer :: i,idx

    allocate (index_x2a(num_imports))
    allocate (index_a2x(num_exports))

    allocate (cpl_names_x2a(num_imports))
    allocate (cpl_names_a2x(num_exports))

    allocate (scr_names_x2a(num_imports))
    allocate (scr_names_a2x(num_exports))

    ! Determine attribute vector indices

    ! List of cpl names of inputs that scream cares about

    !------------------------------------------------------------------------------------------
    !Following inputs are surface values so they are dimensioned (1:ncol) for each chunk.
    !"cam_in" derived type is populated using these inputs in atm_import_export.F90
    !using values from other model components
    !"cam_in" is then used by SCREAM model
    !------------------------------------------------------------------------------------------

    !Comments after the inputs below are organized as follows:
    !Long name [units] (cam_in member which captures the input) [List of parameterizations which are using this input currently]

    cpl_names_x2a(1)  = 'Faxx_evap' ! Surface water vapor flux    [kg/kg](cam_in%cflx(:,1)) [SHOC/check_energy_chng]
    cpl_names_x2a(2)  = 'Faxx_sen'  ! Surface sensible heat flux  [W/m2] (cam_in%shf)  [SHOC/check_energy_chng]
    cpl_names_x2a(3)  = 'Faxx_lat'  ! Surface latent heat flux    [W/m2] (cam_in%lhf)  [energy fixer qqflx_fixer/qneg4]
    cpl_names_x2a(4)  = 'Faxx_taux' ! Surface stress in X         [N/m2] (cam_in%wsx)  [SHOC]
    cpl_names_x2a(5)  = 'Faxx_tauy' ! Surface stress in Y         [N/m2] (cam_in%wsx)  [SHOC]
    cpl_names_x2a(6)  = 'Faxx_lwup' ! long wave up radiation flux [W/m2] (cam_in%lwup) [RRTMGP]
    cpl_names_x2a(7)  = 'Sx_avsdr'  ! short wave direct albedo    [no units] (cam_in%asdir)[RRTMGP]
    cpl_names_x2a(8)  = 'Sx_anidr'  ! long wave direct albedo     [no units] (cam_in%aldir)[RRTMGP]
    cpl_names_x2a(9)  = 'Sx_avsdf'  ! short wave difuse albedo    [no units] (cam_in%asdif)[RRTMGP]
    cpl_names_x2a(10) = 'Sx_anidf'  ! long wave difuse albedo     [no units] (cam_in%aldif)[RRTMGP]
    cpl_names_x2a(11) = 'Sx_t'      ! Surface temperature         [K]        (cam_in%ts)   [check_energy/output- not used anywhere else]
    cpl_names_x2a(12) = 'Sl_snowh'  ! Water equivalent snow depth [m]        (cam_in%snowhland) [SHOC]
    cpl_names_x2a(13) = 'Si_snowh'  ! Snow depth over ice         [m]        (cam_in%snowhice)  [***UNUSED***]
    cpl_names_x2a(14) = 'Sx_tref'   ! Reference height temperature[K]        (cam_in%tref)      [***UNUSED***]

    cpl_names_x2a(15) = 'Sx_qref'   ! Reference height humidity   [kg/kg]    (cam_in%qref)      [***UNUSED***]
    cpl_names_x2a(16) = 'Sx_u10'    ! 10m wind speed              [m/s]      (cam_in%u10)       [***UNUSED***]
    cpl_names_x2a(17) = 'Sf_ifrac'  ! Fraction of sfc area covered by sea-ice [no units] (cam_in%icefrac) [RRTMGP]

    !NOTE: Sf_ofrac (or ocean frac) is being used by aqua_planet and old schemes like vertical_diffusion,
    !Park stratiform_tend and macrophysics
    cpl_names_x2a(18) = 'Sf_ofrac'  ! Fraction of sfc area covered by ocean   [no units] (cam_in%ocnfrac) [***UNUSED***]
    cpl_names_x2a(19) = 'Sf_lfrac'  ! Fraction of sfc area covered by land    [no units] (cam_in%landfrac)[SHOC/RRTMGP/ZM]

    !NOTE:SHOC computes So_ustar (or ustar) internally
    cpl_names_x2a(20) = 'So_ustar'  ! Friction/shear velocity     [m/s]      (cam_in%ustar) [***UNUSED***]
    cpl_names_x2a(21) = 'So_re'     ! ???? (cam_in%re) [***UNUSED***]

    ! Names used by scream for the input fields above
    scr_names_x2a(1)  = 'surface_water_evaporation_flux'
    scr_names_x2a(2)  = 'Faxx_sen'
    scr_names_x2a(3)  = 'Faxx_lat'
    scr_names_x2a(4)  = 'Faxx_taux'
    scr_names_x2a(5)  = 'Faxx_tauy'
    scr_names_x2a(6)  = 'Faxx_lwup'
    scr_names_x2a(7)  = 'Sx_avsdr'
    scr_names_x2a(8)  = 'Sx_anidr'
    scr_names_x2a(9)  = 'Sx_avsdf'
    scr_names_x2a(10) = 'Sx_anidf'
    scr_names_x2a(11) = 'Sx_t'
    scr_names_x2a(12) = 'Sl_snowh'
    scr_names_x2a(13) = 'Si_snowh'
    scr_names_x2a(14) = 'Sx_tref'
    scr_names_x2a(15) = 'Sx_qref'
    scr_names_x2a(16) = 'Sx_u10'
    scr_names_x2a(17) = 'Sf_ifrac'
    scr_names_x2a(18) = 'Sf_ofrac'
    scr_names_x2a(19) = 'Sf_lfrac'
    scr_names_x2a(20) = 'So_ustar'
    scr_names_x2a(21) = 'So_re'

    ! List of cpl names of outputs that scream needs to pass back to cpl

    !Following outputs are computed in control/camsrfexch.F90 using internal SCREAM variables

    !comments after the outputs are organized as follows:
    !Long name [units] (cam_out derived type member) [optional comment about how it is computed]

    cpl_names_a2x(1)  = 'Sa_tbot'     ! Lowest model level temperature [K] (cam_out%tbot)
    cpl_names_a2x(2)  = 'Sa_ptem'     ! Potential temperature          [K] (cam_out%thbot)[Computed from temperature and exner function]
    cpl_names_a2x(3)  = 'Sa_z'        ! Geopotential height above surface at midpoints [m] (cam_out%zbot)
    cpl_names_a2x(4)  = 'Sa_u'        ! Zonal wind        [m/s]  (cam_out%ubot)
    cpl_names_a2x(5)  = 'Sa_v'        ! Meridional wind   [m/s]  (cam_out%vbot)
    cpl_names_a2x(6)  = 'Sa_pbot'     ! midpoint pressure [Pa]   (cam_out%pbot)
    cpl_names_a2x(7)  = 'Sa_dens'     ! Density           [kg/m3](cam_out%rho) [Computed as pbot/(rair*tbot)]
    cpl_names_a2x(8)  = 'Sa_shum'     ! Specific humidity [kg/kg](cam_out%qbot(i,1)[surface water vapor, i.e., state%q(1:ncol,pver,1)]

    !-------------------------------------------------------------------------------------------------
    !Important notes regarding following 4 cpl_names_a2x variables (for cpl_names_a2x indexed 9 to 12):
    !
    !1. All the prec* variables has units of m/s in the model but they are converted to mm/s when
    !they are assigned to the respective cam_out members in components/eam/src/cpl/atm_import_export.F90
    !
    !2. Convective precip variables (precc and precsc, definitions below) should be zero for SCREAM since
    !   convection schemes are turned off in SCREAM
    !   'precc'  is Convective precipitation rate (liq + ice)
    !   'precsc' is Convective snow rate (water equivalent)
    !
    !3. Large scale precip is carried in the following variables:
    !   'precl'  is Large-scale (stable) precipitation rate (liq + ice)
    !   'precsl' is Large-scale (stable) snow rate (water equivalent)
    !-------------------------------------------------------------------------------------------------

    !Faxa_rainc is (precc-precsc), therefore it is just the "liquid" part of the convective prec
    !cam_out variable corresponding to "Faxa_rainc" should be zero for SCREAM
    cpl_names_a2x(9)  = 'Faxa_rainc'  ! Liquid convective precip  [mm/s] (cam_out%precc-cam_out%precsc) [Obtained from Deep conv.]

    !Faxa_rainl is precl-precsl, therefore it is just the "liquid" part of the large scale prec
    cpl_names_a2x(10) = 'Faxa_rainl'  ! Liquid large-scale precip [mm/s] (cam_out%precl-cam_out%precsl) [obtained from P3]

    !cam_out variable corresponding to "Faxa_snowc" should be zero for SCREAM
    cpl_names_a2x(11) = 'Faxa_snowc'  ! Convective snow rate      [mm/s] (cam_out%precsc) [Obtained from Deep Conv.]
    cpl_names_a2x(12) = 'Faxa_snowl'  ! Large-scale (stable) snow rate [mm/s] (cam_out%precsl) [Obtained from P3]

    cpl_names_a2x(13)  = 'Sa_co2prog' ! Always 0.0_r8 as it is not computed by SCREAM (prognostic co2 is turned off)

    ! Names used by scream for the output fields above
    scr_names_a2x(1)  = 'surface_temperature'
    scr_names_a2x(2)  = 'Sa_ptem'
    scr_names_a2x(3)  = 'Sa_z'
    scr_names_a2x(4)  = 'Sa_u'
    scr_names_a2x(5)  = 'Sa_v'
    scr_names_a2x(6)  = 'Sa_pbot'
    scr_names_a2x(7)  = 'Sa_dens'
    scr_names_a2x(8)  = 'Sa_shum'
    scr_names_a2x(9)  = 'Faxa_rainc'
    scr_names_a2x(10) = 'Faxa_rainl'
    scr_names_a2x(11) = 'Faxa_snowc'
    scr_names_a2x(12) = 'Faxa_snowl'
    scr_names_a2x(13)  = 'Sa_co2prog'

    do i=1,num_required_imports
      index_x2a(i) = mct_avect_indexra(x2a,TRIM(cpl_names_x2a(i)))
      scr_names_x2a(i) = TRIM(scr_names_x2a(i)) // C_NULL_CHAR
    enddo
    do i=num_required_imports+1,num_imports
      index_x2a(i) = mct_avect_indexra(x2a,TRIM(cpl_names_x2a(i)),perrWith='quiet')
      scr_names_x2a(i) = TRIM(scr_names_x2a(i)) // C_NULL_CHAR
    enddo

    do i=1,num_required_exports
      index_a2x(i) = mct_avect_indexra(a2x,TRIM(cpl_names_a2x(i)))
      scr_names_a2x(i) = TRIM(scr_names_a2x(i)) // C_NULL_CHAR
    enddo
    do i=num_required_exports+1,num_exports
      index_a2x(i) = mct_avect_indexra(a2x,TRIM(cpl_names_a2x(i)),perrWith='quiet')
      scr_names_a2x(i) = TRIM(scr_names_a2x(i)) // C_NULL_CHAR
    enddo

    ! We no longer need the cpl names
    deallocate(cpl_names_a2x)
    deallocate(cpl_names_x2a)

  end subroutine scream_set_cpl_indices

end module scream_cpl_indices
