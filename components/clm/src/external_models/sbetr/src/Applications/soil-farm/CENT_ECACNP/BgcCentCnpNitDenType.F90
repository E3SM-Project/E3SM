module BgcCentCnpNitDenType

!
! do nitrification denitrification based on century's methods
  use bshr_kind_mod       , only : r8 => shr_kind_r8
implicit none

  private
  character(len=*), private, parameter :: mod_filename = &
       __FILE__

 !real(r8), parameter :: nitrif_n2o_loss_frac = 0.02_r8  ! fraction of N lost as N2O in nitrification (Parton et al., 2001)

  type, public :: century_nitden_type
    real(r8) :: d_con_g(3,2)    ! gas diffusivity constants (spp, #) (cm^2/s) (mult. by 10^-9)
    real(r8) :: d_con_w(3,3)    ! water diffusivity constants (spp, #)  (mult. by 10^-4)

   !parameters
    real(r8) :: organic_max
    real(r8) :: rij_kro_a
    real(r8) :: rij_kro_alpha
    real(r8) :: rij_kro_beta
    real(r8) :: rij_kro_gamma
    real(r8) :: rij_kro_delta
    real(r8) :: nitrif_n2o_loss_frac
    real(r8) :: surface_tension_water
  contains
    procedure, public :: init
    procedure, public :: calc_nitrif_denitrif_rate
    procedure, public :: get_nit_o2_scef
    procedure, public :: run_nitden
    procedure, public :: calc_pot_nitr
    procedure, private:: calc_anaerobic_frac
    procedure, private:: calc_cascade_matrix
    procedure, private:: InitPar
  end type century_nitden_type

  contains

  !-------------------------------------------------------------------------------
  subroutine init(this, biogeo_con)
  use BiogeoConType, only : BiogeoCon_type
  implicit none
  class(century_nitden_type) , intent(inout) :: this
  type(BiogeoCon_type)       , intent(in) :: biogeo_con

  call this%InitPar(biogeo_con)

  end subroutine init

  !-------------------------------------------------------------------------------
  subroutine InitPar(this, biogeo_con)
  use BiogeoConType, only : BiogeoCon_type
  implicit none
  class(century_nitden_type), intent(inout) :: this
  type(BiogeoCon_type),intent(in) :: biogeo_con

  this%organic_max  = biogeo_con%organic_max

  this%rij_kro_a = biogeo_con%rij_kro_a
  this%rij_kro_alpha = biogeo_con%rij_kro_alpha
  this%rij_kro_beta = biogeo_con%rij_kro_beta
  this%rij_kro_gamma = biogeo_con%rij_kro_gamma
  this%rij_kro_delta =  biogeo_con%rij_kro_delta
  this%surface_tension_water = biogeo_con%surface_tension_water

  this%d_con_g(1,:)=(/0.1875_r8, 0.0013_r8/) ! CH4
  this%d_con_g(2,:)=(/0.1759_r8, 0.00117_r8/) ! O2
  this%d_con_g(3,:)=(/0.1325_r8, 0.0009_r8/) ! CO2

  this%d_con_w(1,:)=(/0.9798_r8, 0.02986_r8, 0.0004381_r8/) ! CH4
  this%d_con_w(2,:)=(/1.172_r8, 0.03443_r8, 0.0005048_r8/) ! O2
  this%d_con_w(3,:)=(/0.939_r8, 0.02671_r8, 0.0004095_r8/) ! CO2

  this%nitrif_n2o_loss_frac = 6.e-4_r8
  end subroutine InitPar
  !-------------------------------------------------------------------------------

  subroutine calc_anaerobic_frac(this, o2b, centuryeca_forc, o2_decomp_depth, anaerobic_frac, diffus)
    !
    ! !DESCRIPTION:
    !
    ! calculate soil anoxia state for doing nitrification and denitrification
    ! Rewritten based on Charlie Koven's code by Jinyun Tang
    ! o2_decomp_depth = oxygen consumption from hr + potential nitrification
    ! !USES:
    use betr_varcon        , only : grav => bgrav
    use BgcCentCnpForcType , only : centuryeca_forc_type
    use MathfuncMod        , only : safe_div
    implicit none
    ! !ARGUMENTS:                               !indices
    class(century_nitden_type), intent(inout) :: this
    real(r8),                   intent(in)    :: o2b  ! mol /m3
    type(centuryeca_forc_type), intent(in)    :: centuryeca_forc
    real(r8),                   intent(in)    :: o2_decomp_depth !potential o2 consumption, as deduced from aerobic heteorotrophic decomposition, mol o2/m3/s
    real(r8),                   intent(out)   :: anaerobic_frac       !fraction of aerobic soil
    real(r8),                   intent(out)   :: diffus

    ! !LOCAL VARIABLES:
    real(r8), parameter :: rho_w  = 1.e3_r8        ![kg/m3]
    real(r8) :: f_a
    real(r8) :: eps
    real(r8) :: om_frac
    real(r8) :: organic_max
    real(r8) :: rij_kro_a
    real(r8) :: rij_kro_alpha
    real(r8) :: rij_kro_beta
    real(r8) :: rij_kro_gamma
    real(r8) :: rij_kro_delta
    real(r8) :: surface_tension_water
    real(r8) :: r_min_sat
    real(r8) :: r_psi_sat
    real(r8) :: r_max
    real(r8) :: r_min
    real(r8) :: r_psi
    real(r8) :: ratio_diffusivity_water_gas
    real(r8) :: o2g

    associate(                                        &
         temp       =>    centuryeca_forc%temp      , &
         h2osoi_vol =>    centuryeca_forc%h2osoi_vol, &
         air_vol    =>    centuryeca_forc%air_vol   , &
         watsat     =>    centuryeca_forc%watsat    , &
         watfc      =>    centuryeca_forc%watfc     , &
         bsw        =>    centuryeca_forc%bsw       , &
         soilpsi    =>    centuryeca_forc%soilpsi   , &
         o2_g2b     =>    centuryeca_forc%o2_g2b    , &
         cellorg    =>    centuryeca_forc%cellorg     &
      )
      o2g = o2b / o2_g2b
      surface_tension_water = this%surface_tension_water
      organic_max = this%organic_max
      ! Set parameters from simple-structure model to calculate anoxic fratction (Arah and Vinten 1995)
      rij_kro_a     = this%rij_kro_a
      rij_kro_alpha = this%rij_kro_alpha
      rij_kro_beta  = this%rij_kro_beta
      rij_kro_gamma = this%rij_kro_gamma
      rij_kro_delta = this%rij_kro_delta

      ! calculate gas diffusivity of soil at field capacity here
      ! use expression from methane code, but neglect OM for now
      f_a = 1._r8 - watfc / watsat
      eps =  watsat-watfc ! Air-filled fraction of total soil volume

      if (organic_max > 0._r8) then
        om_frac = min(cellorg/organic_max, 1._r8)
        ! Use first power, not square as in iniTimeConst
      else
        om_frac = 1._r8
      end if
      diffus = (this%d_con_g(2,1) + this%d_con_g(2,2)*temp) * 1.e-4_r8 * &
          (om_frac * f_a**(10._r8/3._r8) / watsat**2 + &
          (1._r8-om_frac) * eps**2 * f_a**(3._r8 / bsw) )

      ! calculate anoxic fraction of soils
      ! use rijtema and kroess model after Riley et al., 2000
      ! caclulated r_psi as a function of psi
      r_min = 2._r8 * surface_tension_water / (rho_w * grav * abs(soilpsi))
      r_max = 2._r8 * surface_tension_water / (rho_w * grav * 0.1_r8)
      r_psi = sqrt(r_min * r_max)
      ratio_diffusivity_water_gas = &
         (this%d_con_g(2,1) + this%d_con_g(2,2)*temp ) * 1.e-4_r8 / &
         ((this%d_con_w(2,1) + this%d_con_w(2,2)*temp + this%d_con_w(2,3)*temp**2) * 1.e-9_r8)

     anaerobic_frac = exp(-rij_kro_a * r_psi**(-rij_kro_alpha) * &
        o2_decomp_depth**(-rij_kro_beta) * o2g**rij_kro_gamma * (h2osoi_vol + &
            ratio_diffusivity_water_gas *  air_vol)**rij_kro_delta)
    end associate
  end subroutine calc_anaerobic_frac

  !-------------------------------------------------------------------------------
  subroutine calc_nitrif_denitrif_rate(this, centuryeca_forc, decompkf_eca,smin_nh4, smin_no3, &
      anaerobic_frac, diffus, pot_f_nit_mol_per_sec, pot_co2_hr, n2_n2o_ratio_denit, pot_f_nit, pot_f_denit)
    !
    ! !DESCRIPTION:
    ! calculate nitrification denitrification rate
    ! the actual nitrification rate will be f_nitr * [nh4]
    ! and the actual denitri rate will be of f_denit * [no3]
    !
    ! !USES:
    use betr_varcon        , only : rpi => brpi, secspday => bsecspday
    use MathfuncMod        , only : safe_div
    use bshr_const_mod     , only : SHR_CONST_TKFRZ
    use BgcCentCnpDecompType     , only : DecompCent_type
    use tracer_varcon      , only : catomw, natomw
    use BgcCentCnpForcType , only : centuryeca_forc_type
    implicit none
    ! !ARGUMENTS:
    class(century_nitden_type) , intent(inout) :: this
    type(centuryeca_forc_type) , intent(in) :: centuryeca_forc
    real(r8)                   , intent(in)  :: pot_f_nit_mol_per_sec
    real(r8)                   , intent(in)  :: pot_co2_hr            !potential aerobic heteotrophic respiration, mol CO2/m3/s
    real(r8)                   , intent(in)  :: anaerobic_frac        !fraction of anaerobic soil
    real(r8)                   , intent(in)  :: smin_nh4
    real(r8)                   , intent(in)  :: smin_no3              !soil no3 concentration [mol N/m3]
    type(DecompCent_type)      , intent(in)  :: decompkf_eca
    real(r8)                   , intent(in)  :: diffus
    real(r8)                   , intent(out) :: n2_n2o_ratio_denit    !ratio of n2 to n2o in denitrification
    real(r8)                   , intent(out) :: pot_f_nit             !nitrification rate of nh4
    real(r8)                   , intent(out) :: pot_f_denit           !potential denitrification rate of no3

    logical, parameter    :: no_frozen_nitrif_denitrif = .false.   !this is a testing parameter, just to make the model run

 ! !LOCAL VARIABLES:
    real(r8) :: soil_hr
    real(r8) :: k_nitr_t
    real(r8) :: k_nitr_ph
    real(r8) :: k_nitr_h2o
    real(r8) :: k_nitr
    real(r8) :: soil_bulkdensity
    real(r8) :: smin_no3_massdens
    real(r8) :: soil_co2_prod
    real(r8) :: fmax_denit_carbonsubstrate
    real(r8) :: fmax_denit_nitrate
    real(r8) :: f_denit_base
    real(r8) :: ratio_k1
    real(r8) :: ratio_no3_co2
    real(r8) :: wfps
    real(r8) :: fr_WFPS
    real(r8) :: co2diff_con(2) = (/0.1325_r8, 0.0009_r8/)
    real(r8) :: g_per_m3__to__ug_per_gsoil
    real(r8) :: g_per_m3_sec__to__ug_per_gsoil_day
    real(r8) :: k_nitr_max

    associate(                                         &
         bd         =>    centuryeca_forc%bd         , & !
         temp       =>    centuryeca_forc%temp       , &
         dz         =>    centuryeca_forc%dzsoi      , &
         watsat     =>    centuryeca_forc%watsat     , & !
         h2osoi_vol =>    centuryeca_forc%h2osoi_vol , & !
         h2osoi_liq =>    centuryeca_forc%h2osoi_liq , & !
         finundated =>    centuryeca_forc%finundated , & ! Input: [real(r8) (:)]
         t_scalar   =>    decompkf_eca%t_scalar      , & ! Input: [real(r8) (:,:)   ]  soil temperature scalar for decomp
         w_scalar   =>    decompkf_eca%w_scalar        & ! Input: [real(r8) (:,:)   ]  soil water scalar for decomp
         )

      pot_f_nit = pot_f_nit_mol_per_sec
      ! limit to oxic fraction of soils
      pot_f_nit  = pot_f_nit * (1._r8 - anaerobic_frac)   ![1/s]

      ! limit to non-frozen soil layers
      if ( temp <= SHR_CONST_TKFRZ .and. no_frozen_nitrif_denitrif) then
        pot_f_nit = 0._r8
      endif

     !---------------- denitrification
     ! first some input variables an unit conversions
     soil_hr = pot_co2_hr * catomw

     ! CENTURY papers give denitrification in units of per gram soil; need to convert from volumetric to mass-based units here
    soil_bulkdensity = bd + h2osoi_liq/dz

    g_per_m3__to__ug_per_gsoil = 1.e3_r8 / soil_bulkdensity

    g_per_m3_sec__to__ug_per_gsoil_day = g_per_m3__to__ug_per_gsoil * secspday

    smin_no3_massdens = max(smin_no3, 0._r8) * g_per_m3__to__ug_per_gsoil * natomw

    soil_hr = pot_co2_hr * catomw
    soil_co2_prod = (soil_hr * (g_per_m3_sec__to__ug_per_gsoil_day))

    !! maximum potential denitrification rates based on heterotrophic respiration rates or nitrate concentrations,
    !! from (del Grosso et al., 2000)
    fmax_denit_carbonsubstrate = (0.1_r8 * (soil_co2_prod**1.3_r8)) &
                 / g_per_m3_sec__to__ug_per_gsoil_day
    !
   fmax_denit_nitrate = (1.15_r8 * smin_no3_massdens**0.57_r8)  &
                 / g_per_m3_sec__to__ug_per_gsoil_day

   ! find limiting denitrification rate
   f_denit_base = max(min(fmax_denit_carbonsubstrate, fmax_denit_nitrate),0._r8)

   ! limit to non-frozen soil layers
   if ( temp <= SHR_CONST_TKFRZ .and. no_frozen_nitrif_denitrif ) then
     f_denit_base = 0._r8
   endif

   ! limit to anoxic fraction of soils
   pot_f_denit = f_denit_base * anaerobic_frac

   ! now calculate the ratio of N2O to N2 from denitrifictaion, following Del Grosso et al., 2000
   ! diffusivity constant (figure 6b)
   ratio_k1 = max(1.7_r8, 38.4_r8 - 350._r8 * diffus)

   ! ratio function (figure 7c)
   if ( soil_co2_prod > 0 ) then
     ratio_no3_co2 = smin_no3_massdens / soil_co2_prod
   else
     ! fucntion saturates at large no3/co2 ratios, so set as some nominally large number
     ratio_no3_co2 = 100._r8
   endif

   ! total water limitation function (Del Grosso et al., 2000, figure 7a)
   wfps = max(min(h2osoi_vol/watsat, 1._r8), 0._r8) * 100._r8
   fr_WFPS = max(0.1_r8, 0.015_r8 * wfps - 0.32_r8)

   ! final ratio expression
   n2_n2o_ratio_denit = max(0.16*ratio_k1, ratio_k1*exp(-0.8 * ratio_no3_co2)) * fr_WFPS

   end associate
  end subroutine calc_nitrif_denitrif_rate

  !-------------------------------------------------------------------------------
  subroutine calc_pot_nitr(this, smin_nh4, centuryeca_forc, decompkf_eca, pot_f_nit_mol_per_sec)
  !DESCRIPTION
  !calculate potential nitrificaiton rate

  use betr_varcon        , only : rpi => brpi, secspday => bsecspday
  use BgcCentCnpForcType , only : centuryeca_forc_type
  use BgcCentCnpDecompType     , only : DecompCent_type
  implicit none
  class(century_nitden_type) , intent(inout) :: this
  real(r8)                   , intent(in) :: smin_nh4
  type(centuryeca_forc_type) , intent(in) :: centuryeca_forc
  type(DecompCent_type)      , intent(in) :: decompkf_eca
  real(r8)                   , intent(out) :: pot_f_nit_mol_per_sec

  !local variables
  real(r8) :: k_nitr
  real(r8) :: k_nitr_max
  real(r8) :: k_nitr_t
  real(r8) :: k_nitr_ph
  real(r8) :: k_nitr_h2o

  associate(                                    &
     pH           =>    centuryeca_forc%pH    , &
     t_scalar     =>    decompkf_eca%t_scalar , & ! Input: [real(r8) (:,:)   ]  soil temperature scalar for decomp
     w_scalar     =>    decompkf_eca%w_scalar   & ! Input: [real(r8) (:,:)   ]  soil water scalar for decomp
  )

  ! Set maximum nitrification rate constant
  k_nitr_max =  0.1_r8 / secspday   ! [1/sec] 10%/day  Parton et al., 2001

  !---------------- nitrification
  ! follows CENTURY nitrification scheme (Parton et al., (2001, 1996))

  ! assume nitrification temp function equal to the HR scalar
  k_nitr_t = min(t_scalar, 1._r8)

  ! ph function from Parton et al., (2001, 1996)
  k_nitr_ph = 0.56_r8 + atan(rpi * 0.45_r8 * (-5._r8+ pH))/rpi

  ! moisture function-- assume the same moisture function as limits heterotrophic respiration
  ! Parton et al. base their nitrification- soil moisture rate constants based on heterotrophic rates-- can we do the same?
  k_nitr_h2o = w_scalar

  ! nitrification constant is a set scalar * temp, moisture, and ph scalars
  k_nitr = k_nitr_max * k_nitr_t * k_nitr_h2o * k_nitr_ph

  ! first-order decay of ammonium pool with scalar defined above
  pot_f_nit_mol_per_sec = max(k_nitr, 0._r8) * smin_nh4

  end associate
  end subroutine calc_pot_nitr

  !---------------------------------------------------------------------------------
  subroutine calc_cascade_matrix(this, centurybgc_index, n2_n2o_ratio_denit, cascade_matrix)

  use BgcCentCnpIndexType, only : centurybgc_index_type
  implicit none
  class(century_nitden_type)  , intent(inout) :: this
  type(centurybgc_index_type) , intent(in) :: centurybgc_index
  real(r8)                    , intent(in) :: n2_n2o_ratio_denit
  real(r8)                    , intent(inout)    :: cascade_matrix(centurybgc_index%nstvars, centurybgc_index%nreactions)

  integer :: reac

  associate(                                                &
    primvarid    => centurybgc_index%primvarid            , & !
    lid_nh4   => centurybgc_index%lid_nh4                 , & !
    lid_o2   => centurybgc_index%lid_o2                   , & !
    lid_n2   => centurybgc_index%lid_n2                   , & !
    lid_n2o   => centurybgc_index%lid_n2o                 , & !
    lid_no3   => centurybgc_index%lid_no3                 , & !
    lid_no3_den => centurybgc_index%lid_no3_den           , &
    lid_nh4_nit_reac => centurybgc_index%lid_nh4_nit_reac , & !
    lid_no3_den_reac => centurybgc_index%lid_no3_den_reac , & !
    lid_nh4_nit        => centurybgc_index%lid_nh4_nit    , & !
    lid_n2o_nit=> centurybgc_index%lid_n2o_nit              & !

  )
  !---------------------------------------------------------------------------------
  !reaction 8, nitrification
  reac = lid_nh4_nit_reac
  !NH4(+) + (2-f)O2 + (2-f)OH(-)-> (1-f)NO3(-) + (f/2)N2O + (3-f/2) H2O
  cascade_matrix(lid_nh4 ,reac) = -1._r8
  cascade_matrix(lid_o2  ,reac) = -(2._r8 - this%nitrif_n2o_loss_frac)
  cascade_matrix(lid_no3 ,reac) = 1._r8  - this%nitrif_n2o_loss_frac
  cascade_matrix(lid_n2o, reac) = 0.5_r8 * this%nitrif_n2o_loss_frac

  cascade_matrix(lid_nh4_nit,reac) = 1._r8
  cascade_matrix(lid_n2o_nit,reac) = this%nitrif_n2o_loss_frac

  !---------------------------------------------------------------------------------
  !reaction 9, denitrification
  reac = lid_no3_den_reac
  !NO3(-) -> 0.5*f N2O + 0.5* (1-f) N2, where f is a function determined from the century denitrification model
  cascade_matrix(lid_no3 ,reac)    = -1._r8
  cascade_matrix(lid_n2o ,reac)    = 0.5_r8 * 1._r8/(1._r8+n2_n2o_ratio_denit)
  cascade_matrix(lid_n2  ,reac)    = 0.5_r8 - cascade_matrix(lid_n2o ,reac)
  cascade_matrix(lid_no3_den,reac) = 1._r8

  end associate
  end subroutine calc_cascade_matrix

  !---------------------------------------------------------------------------------
  function get_nit_o2_scef(this)result(ans)
  !
  ! get the oxygen requirement for nitrification
  ! the following calculation is a crude estimation
  ! the full conversion of NH3 into nitrate costs 2 mol oxygen
  ! molecules. N2O may either produced from decomposition
  ! of hydroxylamine (NH2OH) or nitrifier detoxification process.
  ! the first part cost 0.5 mol oxygen for the formation of
  ! 1 mol NH2OH.
  implicit none
  class(century_nitden_type), intent(inout) :: this
  real(r8) :: ans

  ans = (2._r8 + 0.5_r8*this%nitrif_n2o_loss_frac)

  end function get_nit_o2_scef
  !---------------------------------------------------------------------------------
  subroutine run_nitden(this, centurybgc_index,centuryeca_forc, decompkf_eca,&
    smin_nh4, smin_no3, o2b, o2_decomp_depth, pot_f_nit_mol_per_sec, pot_co2_hr, &
    pot_f_nit, pot_f_denit, cascade_matrix)
  !this returns the vamx for nitrification and denitrification
  !also return is the stoichiometry matrix for nitden processes
  use BgcCentCnpDecompType      , only : DecompCent_type
  use BgcCentCnpForcType  , only : centuryeca_forc_type
  use BgcCentCnpIndexType , only : centurybgc_index_type
  implicit none
  class(century_nitden_type)  , intent(inout) :: this
  type(centurybgc_index_type) , intent(in) :: centurybgc_index
  type(centuryeca_forc_type)  , intent(in) :: centuryeca_forc
  type(DecompCent_type)       , intent(in) :: decompkf_eca
  real(r8)                    , intent(in) :: o2b
  real(r8)                    , intent(in) :: o2_decomp_depth
  real(r8)                    , intent(in) :: smin_nh4  ! mol Nh4/m3
  real(r8)                    , intent(in) :: smin_no3  ! mol NO3/m3
  real(r8)                    , intent(in) :: pot_f_nit_mol_per_sec ! potential nitrification mol N /s
  real(r8)                    , intent(in) :: pot_co2_hr  ! potential co2 emission from heterotrophic respiration
  real(r8)                    , intent(out) :: pot_f_nit  ! mol N /s
  real(r8)                    , intent(out) :: pot_f_denit  ! mol N/s
  real(r8)                    , intent(inout) :: cascade_matrix(centurybgc_index%nstvars, centurybgc_index%nreactions)

  !local variables
  real(r8) :: n2_n2o_ratio_denit
  real(r8) :: anaerobic_frac
  real(r8) :: diffus

  call this%calc_anaerobic_frac(o2b, centuryeca_forc, o2_decomp_depth, anaerobic_frac, diffus)

  !calculate potential denitrification and denitrification
  call this%calc_nitrif_denitrif_rate(centuryeca_forc, decompkf_eca, &
      smin_nh4, smin_no3, anaerobic_frac,diffus, pot_f_nit_mol_per_sec, pot_co2_hr,&
      n2_n2o_ratio_denit, pot_f_nit, pot_f_denit)

  !calcualte cascade matrix for nitrification denitrification
  call this%calc_cascade_matrix(centurybgc_index, n2_n2o_ratio_denit, cascade_matrix)
  end subroutine run_nitden

end module BgcCentCnpNitDenType
