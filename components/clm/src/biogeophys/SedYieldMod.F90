module SedYieldMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate the sediment flux caused by soil erosion based on the 
  ! improved Morgan model (Tan et al., 2017 & 2018) 
  !
  use shr_const_mod     , only : T0 => SHR_CONST_TKFRZ
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use decompMod         , only : bounds_type
  use clm_varcon        , only : grav, denh2o, rpi
  use clm_varpar        , only : mxpft, nlevsno, max_patch_per_col
  use atm2lndType       , only : atm2lnd_type
  use CanopyStateType   , only : CanopyState_type
  use EnergyFluxType    , only : energyflux_type
  use SoilHydrologyType , only : soilhydrology_type
  use SoilStateType     , only : soilstate_type
  use WaterfluxType     , only : waterflux_type
  use WaterStateType    , only : waterstate_type
  use TemperatureType   , only : temperature_type
  use ColumnType        , only : col_pp
  use LandunitType      , only : lun_pp
  use VegetationType    , only : veg_pp
  use SedFluxType       , only : sedflux_type 
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: SoilErosion          ! Calculate hillslope sediment flux
  public :: Landslide            ! Calculate landslide-driven sediment flux 
  ! !MODULE CONSTANTS
  ! a free parameter for soil cohesion
  real(r8), parameter :: lndsld_cs = 0.0196_r8
  ! a free parameter for root-induced soil shear strength
  real(r8), parameter :: lndsld_cr = 0.0031_r8
  ! a free parameter for power law of exponent for landslide
  real(r8), parameter :: lndsld_xi = 1.3529_r8
  ! a free parameter for landslide volume upper limit (m^3)
  real(r8), parameter :: lndsld_Vmax = 1.3344d+03
  ! landslide cut-off volume (m^3)
  real(r8), parameter :: lndsld_Vmin = 45._r8
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine SoilErosion (bounds, num_hydrologyc, filter_hydrologyc, &
    atm2lnd_vars, canopystate_vars, soilstate_vars, waterstate_vars, &
    waterflux_vars, sedflux_vars)
    !
    ! !DESCRIPTION:
    ! Calculate rainfall and runoff driven erosion 
    !
    ! !USES:
    use clm_time_manager, only : get_step_size
    use clm_varctl      , only : iulog
    use landunit_varcon , only : istcrop, istsoil
    use pftvarcon       , only : gcpsi, pftcc
    use spmdMod         , only : masterproc
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
    integer                  , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
    type(atm2lnd_type)       , intent(in)    :: atm2lnd_vars
    type(CanopyState_type)   , intent(in)    :: canopystate_vars
    type(soilstate_type)     , intent(in)    :: soilstate_vars
    type(waterstate_type)    , intent(in)    :: waterstate_vars
    type(waterflux_type)     , intent(in)    :: waterflux_vars
    type(sedflux_type)       , intent(inout) :: sedflux_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: c, fc, p, pi, l, j                         ! indices
    integer  :: dtime                                      ! timestep size [seconds]
    integer  :: nslp                                       ! slope level #
    real(r8) :: sinslp, fslp, factor_slp                   ! topo gradient
    real(r8) :: gndbare, vegcc, fbare                      ! veg cover factors
    real(r8) :: fungrvl                                    ! ground uncovered by gravel
    real(r8) :: K, COH                                     ! soil erodibility
    real(r8) :: Qs, Ptot, Ie, Dl                           ! water fluxes
    real(r8) :: Tc, Es_Q, Es_P, KE_DT, KE_LD               ! temporaries 
    real(r8) :: stxt(4)                                    ! soil texture including gravel
    character(len=32) :: subname = 'SoilErosion'           ! subroutine name
    !-----------------------------------------------------------------------

    associate(                                                        &
         forc_rain        =>    atm2lnd_vars%forc_rain_downscaled_col   , & ! Input: [real(r8) (:) ] rain rate [mm/s]
         forc_t           =>    atm2lnd_vars%forc_t_downscaled_col      , & ! Input: [real(r8) (:) ] atmospheric temperature (Kelvin) 

         tlai             =>    canopystate_vars%tlai_patch         , & ! Input: [real(r8) (:) ] LAI
         hbot             =>    canopystate_vars%hbot_patch         , & ! Input: [real(r8) (:) ] canopy bottom (m)
         htop             =>    canopystate_vars%htop_patch         , & ! Input: [real(r8) (:) ] canopy top (m) 

         hslp             =>    col_pp%hslp                         , & ! Input: [real(r8) (:) ]  hillslope gradient
         hslp_p10         =>    col_pp%hslp_p10                     , & ! Input: [real(r8) (:,:) ] hillslope gradient percentiles
    
         bd               =>    soilstate_vars%bd_col               , & ! Input: [real(r8) (:,:) ] soil bulk density (kg/m3)
         fsand            =>    soilstate_vars%cellsand_col         , & ! Input: [real(r8) (:,:) ] sand percentage
         fclay            =>    soilstate_vars%cellclay_col         , & ! Input: [real(r8) (:,:) ] clay percentage
         fgrvl            =>    soilstate_vars%cellgrvl_col         , & ! Input: [real(r8) (:,:) ] gravel percentage

         hsnow            =>    waterstate_vars%snow_depth_col      , & ! Input: [real(r8) (:) ] snow depth (m)

         prnsf            =>    sedflux_vars%prnsf_col              , & ! Input: [real(r8) (:) ] rainfall-driven erosion scaling factor
         qrnsf            =>    sedflux_vars%qrnsf_col              , & ! Input: [real(r8) (:) ] runoff-driven erosion scaling factor
         tcsf             =>    sedflux_vars%tcsf_col               , & ! Input: [real(r8) (:) ] transport capacity scaling factor 

         qflx_surf        =>    waterflux_vars%qflx_surf_col            , & ! Input: [real(r8) (:) ] surface runoff (mm/s)
         qflx_dirct_rain  =>    waterflux_vars%qflx_dirct_rain_patch    , & ! Input: [real(r8) (:) ] direct throughfall rain (mm/s)
         qflx_leafdrip    =>    waterflux_vars%qflx_leafdrip_patch      , & ! Input: [real(r8) (:) ] leaf rain drip (mm/s)

         flx_p_ero        =>    sedflux_vars%flx_p_ero_col          , & ! Output: [real(r8) (:) ] sed detached by rainfall (kg/m2/s)
         flx_q_ero        =>    sedflux_vars%flx_q_ero_col          , & ! Output: [real(r8) (:) ] sed detached by runoff (kg/m2/s)
         flx_sed_ero      =>    sedflux_vars%flx_sed_ero_col        , & ! Output: [real(r8) (:) ] total detachment (kg/m2/s)
         flx_sed_yld      =>    sedflux_vars%flx_sed_yld_col          & ! Output: [real(r8) (:) ] sed flux to inland waters (kg/m2/s)
         )

         ! Get time step
         dtime = get_step_size()

         ! nolakec or other col filters
         do fc = 1, num_hydrologyc
            c = filter_hydrologyc(fc)
            l = col_pp%landunit(c)

            ! initialization
            flx_p_ero(c)          = 0._r8
            flx_q_ero(c)          = 0._r8
            flx_sed_ero(c)        = 0._r8
            flx_sed_yld(c)        = 0._r8

            ! check landunit type and ground covered by snow/ice
            if ( lun_pp%itype(l)/=istsoil .and. lun_pp%itype(l)/=istcrop ) then
               cycle
            else if (hsnow(c)>0) then
               cycle
            end if

            ! soil detachment by rainfall
            stxt = (/fclay(c,1), 100.-fclay(c,1)-fsand(c,1), fsand(c,1), &
                fgrvl(c,1)/)  
            K = SoilDetachability(stxt)
            COH = SoilCohesion(stxt)
            Es_P = 0._r8    ! detachment by throughfall + leap drip
            if (forc_t(c)>T0 .and. forc_rain(c)>0._r8) then
               fungrvl = 1. - 0.01 * fgrvl(c,1)
               do pi = 1, max_patch_per_col
                  if ( pi<=col_pp%npfts(c) ) then
                     p = col_pp%pfti(c) + pi - 1
                     if (veg_pp%active(p) .and. veg_pp%wtcol(p)>0) then
                        ! throughfall power
                        Ptot = qflx_dirct_rain(p) * dtime   ! mm
                        Ie = 3.6d3 * forc_rain(c)           ! mm/hr
                        KE_DT = Ptot * fungrvl * max(0._r8,8.95+8.44*log10(Ie))
                        ! leaf drip power
                        Dl = max(0._r8, qflx_leafdrip(p)*dtime)      ! mm
                        KE_LD = max(0._r8,15.8*sqrt(0.5*(htop(p)+hbot(p)))-5.87) * &
                           fungrvl * Dl
                        Es_P = Es_P + prnsf(c) * veg_pp%wtcol(p) * K * (KE_DT+KE_LD)
                     end if
                  end if
               end do
            end if
            Es_P = 1d-7 / 8.64 * Es_P   ! kg/m2/s

            ! soil detachment by runoff
            gndbare = 0._r8
            vegcc = 0._r8
            do pi = 1, max_patch_per_col
               if ( pi<=col_pp%npfts(c) ) then
                  p = col_pp%pfti(c) + pi - 1
                  if (veg_pp%active(p)) then
                     fbare = exp( -gcpsi(veg_pp%itype(p)) * tlai(p) )
                     gndbare = gndbare + veg_pp%wtcol(p) * fbare
                     vegcc = vegcc + veg_pp%wtcol(p) * pftcc(veg_pp%itype(p))
                  end if
               end if
            end do
            
            Es_Q = 0._r8
            Tc = 0._r8
            if (qflx_surf(c)>0) then
               Qs = 8.64d4 * qflx_surf(c)  ! mm/d
               nslp = size(hslp_p10(c,:))
               factor_slp = 0._r8
               do j = 1, nslp
                  if (j==1) then
                     fslp = 1.5 / DBLE(nslp)
                  else if (j<nslp-1) then
                     fslp = 1.0 / DBLE(nslp)
                  else if (j==nslp-1) then
                     fslp = 1.3 / DBLE(nslp)
                  else
                     fslp = 0.2 / DBLE(nslp)
                  end if
                  sinslp = sin(atan(max(hslp_p10(c,j),1d-4)))
                  factor_slp = factor_slp + fslp * sinslp
               end do
               Es_Q = 19.1 * qrnsf(c) * 2.0 / COH * factor_slp * &
                  gndbare * Qs**1.5
               Tc = 19.1 * tcsf(c) * factor_slp * vegcc * Qs**2.0
            end if
            Es_Q = 1d-7 / 8.64 * Es_Q
            Tc = 1d-7 / 8.64 * Tc

            ! assign flux values
            flx_p_ero(c) = flx_p_ero(c) + Es_P
            flx_q_ero(c) = flx_q_ero(c) + Es_Q
            flx_sed_ero(c) = flx_sed_ero(c) + Es_P + Es_Q 
            flx_sed_yld(c) = flx_sed_yld(c) + min(Es_P+Es_Q, Tc)
         end do
    
    end associate

  end subroutine SoilErosion

  !-----------------------------------------------------------------------
  subroutine Landslide (bounds, num_hydrologyc, filter_hydrologyc, &
    soilstate_vars, waterstate_vars, sedflux_vars)
    !
    ! !DESCRIPTION:
    ! Calculate landslide-induced sediment detachment 
    !
    ! !USES:
    use clm_varctl      , only : iulog
    use clm_varpar      , only : nlevsoi, maxpatch_pft
    use landunit_varcon , only : istcrop, istsoil
    use spmdMod         , only : masterproc
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
    integer                  , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
    type(soilstate_type)     , intent(in)    :: soilstate_vars
    type(waterstate_type)    , intent(in)    :: waterstate_vars
    type(sedflux_type)       , intent(inout) :: sedflux_vars
    ! !LOCAL VARIABLES:
    integer  :: c, j, fc, p, pi, l                          ! indices
    integer  :: nslp                                        ! slope level #
    real(r8) :: FS                                          ! safety factor
    real(r8) :: hslp_max, omega, sin_omega, tan_omega       ! topo variables
    real(r8) :: tan_phi                                     ! soil internal friction
    real(r8) :: zsum, zn                                    ! soil depth (m)
    real(r8) :: ssum                                        ! total soil mass (kg/m2)
    real(r8) :: vwat, vwatsat, vwatmin                      ! cur, saturated and min soil moisture (v/v)
    real(r8) :: theta                                       ! effective saturation
    real(r8) :: Cs, Cr                                      ! intrinsic and root-induced soil cohesion (Pa)
    real(r8) :: ftree                                       ! tree fraction
    real(r8) :: Rwat, Rsed                                  ! unit weight (N/m3)
    real(r8) :: alpha                                       ! Van Genuchten parameter
    real(r8) :: rn                                          ! pore-shape distribution parameter
    real(r8) :: h                                           ! pressure head (m)
    real(r8) :: expn, vol                                   ! landslide run-out volume variables
    real(r8) :: fslp                                        ! slope level areal fraction
    real(r8) :: Es_L                                        ! landslide erosion (kg/m2/s)
    real(r8) :: stxt(4)                                     ! soil texture
    character(len=32) :: subname = 'Landslide'              ! subroutine name
    !-----------------------------------------------------------------------
  
    associate(                                                        &
         dz               =>    col_pp%dz                           , & ! Input: [real(r8) (:,:) ] layer thickness (m)                                 
         hslp_p10         =>    col_pp%hslp_p10                     , & ! Input: [real(r8) (:,:) ] hillslope gradient percentiles 
         znsoil           =>    col_pp%znsoil                       , & ! Input: [real(r8) (:)   ] soil depth (m)

         bd               =>    soilstate_vars%bd_col               , & ! Input: [real(r8) (:,:) ] soil bulk density (kg/m3)
         fsand            =>    soilstate_vars%cellsand_col         , & ! Input: [real(r8) (:,:) ] sand fraction
         fclay            =>    soilstate_vars%cellclay_col         , & ! Input: [real(r8) (:,:) ] clay fraction
         fgrvl            =>    soilstate_vars%cellgrvl_col         , & ! Input: [real(r8) (:,:) ] gravel fraction
         watmin           =>    soilstate_vars%watmin_col           , & ! Input: [real(r8) (:,:) ] minimum volumetric soil water (v/v)
         watsat           =>    soilstate_vars%watsat_col           , & ! Input: [real(r8) (:,:) ] soil porosity (v/v) 

         h2osoi_vol       =>    waterstate_vars%h2osoi_vol_col      , & ! Input: [real(r8) (:,:) ] volumetric soil water (v/v)

         flx_lndsld_ero   =>    sedflux_vars%flx_lndsld_ero_col     , & ! Output: [real(r8) (:) ] landslide sed detach (kg/m2/s)
         flx_lndsld_yld   =>    sedflux_vars%flx_lndsld_yld_col     , & ! Output: [real(r8) (:) ] landslide sed flux to rivers (kg/m2/s)
         flx_sed_yld      =>    sedflux_vars%flx_sed_yld_col          & ! Output: [real(r8) (:) ] total sed flux to rivers (kg/m2/s)
         )

         do fc = 1, num_hydrologyc
            c = filter_hydrologyc(fc)
            l = col_pp%landunit(c)

            ! initialization
            flx_lndsld_ero(c) = 0._r8
            flx_lndsld_yld(c) = 0._r8 

            ! check the landunit type
            if (lun_pp%itype(l)/=istsoil .and. lun_pp%itype(l)/=istcrop) then
               cycle
            end if

            ! no steep hillslopes
            tan_phi = tan(rpi*25./180.)
            hslp_max = maxval(hslp_p10(c,:))
            if (hslp_max<=tan_phi) then
               cycle
            end if

            zsum = 0._r8
            vwat = 0._r8
            vwatsat = 0._r8 
            vwatmin = 0._r8
            Cs = 0._r8
            Cr = 0._r8
            Rsed = 0._r8
            alpha = 0._r8
            rn = 0._r8
            ssum = 0._r8
            do j = 1, nlevsoi
               zsum = zsum + dz(c,j)
               vwat = vwat + h2osoi_vol(c,j)*dz(c,j)
               vwatsat = vwatsat + watsat(c,j)*dz(c,j) 
               vwatmin = vwatmin + watmin(c,j)*dz(c,j)
               Rsed = Rsed + (bd(c,j)+h2osoi_vol(c,j)*denh2o)*dz(c,j)
               stxt = (/fclay(c,j), 100.-fclay(c,j)-fsand(c,j), fsand(c,j), &
                  fgrvl(c,j)/)
               Cs = Cs + 1d3*lndsld_cs*SoilCohesion2(stxt)*dz(c,j)
               Cr = Cr + 5d3*lndsld_cr*dz(c,j)
               alpha = alpha + SoilVanGenuchten(stxt)*dz(c,j)
               rn = rn + SoilPoreShape(stxt)*dz(c,j)
               ssum = ssum + bd(c,j)*dz(c,j)
            end do

            ftree = 0._r8
            do pi = 1, max_patch_per_col
               if ( pi<=col_pp%npfts(c) ) then
                  p = col_pp%pfti(c) + pi - 1
                  if (veg_pp%active(p)) then
                     if (veg_pp%itype(p)>=1 .and. veg_pp%itype(p)<=8) then
                        ftree = ftree + veg_pp%wtcol(p)
                     end if
                  end if
               end if
            end do

            vwat = vwat / zsum
            vwatsat = vwatsat / zsum
            vwatmin = vwatmin / zsum
            Rwat = denh2o * grav
            Rsed = Rsed * grav / zsum
            Cs = Cs / zsum
            Cr = ftree * Cr / zsum
            alpha = alpha / zsum
            rn = rn / zsum
            theta = min( (vwat-vwatmin)/(vwatsat-vwatmin), 1.0 )
            h = -1._r8/alpha*(theta**(-rn/(rn-1.))-1.)**(1./rn)

            Es_L = 0._r8
            nslp = size(hslp_p10(c,:))
            do j = 1, nslp
               if (j==1) then
                  fslp = 1.5 / DBLE(nslp)
               else if (j<nslp-1) then
                  fslp = 1.0 / DBLE(nslp)
               else if (j==nslp-1) then
                  fslp = 1.3 / DBLE(nslp)
               else
                  fslp = 0.2 / DBLE(nslp)
               end if
               if (hslp_p10(c,j)>tan_phi) then
                  omega = atan(hslp_p10(c,j))
                  sin_omega = sin(omega)
                  tan_omega = tan(omega)
                  zn = znsoil(c) * cos(omega)
                  FS = (Cs+Cr)/(zn*Rsed*sin_omega) + tan_phi/tan_omega - &
                     (Rwat*h)/(Rsed*zn)*theta*tan_phi/sin_omega
                  if (FS<1.0) then
                     expn = 2. - lndsld_xi
                     if (expn==0.0) then
                        vol = log(1d3/lndsld_Vmin)
                     else
                        vol = (lndsld_Vmax**expn - lndsld_Vmin**expn) / expn
                     end if
                     Es_L = Es_L + 1d-10 / 8.64 * fslp * Rsed * vol
                  end if
               end if
            end do

            flx_lndsld_ero(c) = flx_lndsld_ero(c) + Es_L
            flx_lndsld_yld(c) = flx_lndsld_yld(c) + Es_L
            flx_sed_yld(c) = flx_sed_yld(c) + Es_L
         end do

    end associate

  end subroutine Landslide

  !------------------------------------------------------------------------------
  character(len=32) function SoilTextureType(stxt)
    !
    ! !DESCRIPTION:
    ! soil texture 
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: stxt(4)  ! [clay, silt, sand, gravel] in percentage
    !
    ! !LOCAL VARIABLES:
    real(r8) :: clay, silt, sand, tsum 
    !------------------------------------------------------------------------------

    tsum = sum(stxt(1:3))
    clay = stxt(1) / tsum 
    silt = stxt(2) / tsum
    sand = stxt(3) / tsum
    if (silt + 1.5*clay < 15) then
       SoilTextureType = 'sand'
    else if (silt + 2*clay < 30) then
       SoilTextureType = 'loamy sand'
    else if (((clay>=7 .and. clay<20) .and. sand>52) .or. &
            (clay<7 .and. silt<50)) then
       SoilTextureType = 'sandy loam'
    else if ((clay>=7 .and. clay<27) .and. (silt>=28 .and. silt<50) .and. &
            sand<=52) then
       SoilTextureType = 'loam'
    else if ((silt>=50 .and. (clay>=12 .and. clay<27)) .or. &
            ((silt>=50 .and. silt<80) .and. clay<12)) then
       SoilTextureType = 'silt loam'
    else if (silt>=80 .and. clay<12) then
       SoilTextureType = 'silt'
    else if ((clay>=20 .and. clay<35) .and. silt<28 .and. sand>45) then
       SoilTextureType = 'sandy clay loam'
    else if ((clay>=27 .and. clay<40) .and. (sand>20 .and. sand<=45)) then
       SoilTextureType = 'clay loam'
    else if ((clay>=27 .and. clay<40) .and. sand<=20) then
       SoilTextureType = 'silty clay loam'
    else if (clay>=35 .and. sand>45) then
       SoilTextureType = 'sandy clay'
    else if (clay>=40 .and. silt>=40) then
       SoilTextureType = 'silty clay'
    else if (clay>=40 .and. sand<=45 .and. silt<40) then
       SoilTextureType = 'clay'     
    end if

  end function SoilTextureType

  !------------------------------------------------------------------------------
  real(r8) function SoilDetachability(stxt)
    !
    ! !DESCRIPTION:
    ! soil detachability by rain drops (g J-1)
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: stxt(4)  ! [clay, silt, sand, gravel] in fraction 
    !
    ! !LOCAL VARIABLES:
    character(len=32) :: txttype
    !------------------------------------------------------------------------------

    txttype = SoilTextureType(stxt) 
    if (trim(txttype)=='clay') then
       SoilDetachability = 0.05
    else if (trim(txttype)=='silty clay') then
       SoilDetachability = 0.5
    else if (trim(txttype)=='sandy clay') then
       SoilDetachability = 0.3
    else if (trim(txttype)=='silty clay loam') then
       SoilDetachability = 0.8
    else if (trim(txttype)=='clay loam') then
       SoilDetachability = 0.7
    else if (trim(txttype)=='sandy clay loam') then
       SoilDetachability = 0.1
    else if (trim(txttype)=='silt') then
       SoilDetachability = 1.0
    else if (trim(txttype)=='silt loam') then
       SoilDetachability = 0.9
    else if (trim(txttype)=='loam') then
       SoilDetachability = 0.8
    else if (trim(txttype)=='sandy loam') then
       SoilDetachability = 0.7
    else if (trim(txttype)=='loamy sand') then
       SoilDetachability = 0.3
    else if (trim(txttype)=='sand') then
       SoilDetachability = 1.9
    end if

  end function SoilDetachability

  !------------------------------------------------------------------------------
  real(r8) function SoilCohesion(stxt)
    !
    ! !DESCRIPTION:
    ! soil cohesion against concentrated flow (kPa)
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: stxt(4)  ! [clay, silt, sand, gravel] in fraction 
    !
    ! !LOCAL VARIABLES:
    character(len=32) :: txttype
    !------------------------------------------------------------------------------

    txttype = SoilTextureType(stxt)
    if (trim(txttype)=='clay') then
       SoilCohesion = 12.0
    else if (trim(txttype)=='silty clay') then
       SoilCohesion = 10.0
    else if (trim(txttype)=='sandy clay') then
       SoilCohesion = 10.0
    else if (trim(txttype)=='silty clay loam') then
       SoilCohesion = 9.0
    else if (trim(txttype)=='clay loam') then
       SoilCohesion = 10.0
    else if (trim(txttype)=='sandy clay loam') then
       SoilCohesion = 3.0
    else if (trim(txttype)=='silt') then
       SoilCohesion = 3.0
    else if (trim(txttype)=='silt loam') then
       SoilCohesion = 3.0
    else if (trim(txttype)=='loam') then
       SoilCohesion = 3.0
    else if (trim(txttype)=='sandy loam') then
       SoilCohesion = 2.0
    else if (trim(txttype)=='loamy sand') then
       SoilCohesion = 2.0
    else if (trim(txttype)=='sand') then
       SoilCohesion = 2.0
    end if

  end function SoilCohesion

  !------------------------------------------------------------------------------
  real(r8) function SoilCohesion2(stxt)
    !
    ! !DESCRIPTION:
    ! soil cohesion against landslide (kPa)
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: stxt(4)  ! [clay, silt, sand, gravel] in fraction 
    !
    ! !LOCAL VARIABLES:
    character(len=32) :: txttype
    !------------------------------------------------------------------------------

    txttype = SoilTextureType(stxt)
    if (trim(txttype)=='clay') then
       SoilCohesion2 = 10.0
    else if (trim(txttype)=='silty clay') then
       SoilCohesion2 = 6.0
    else if (trim(txttype)=='sandy clay') then
       SoilCohesion2 = 6.0
    else if (trim(txttype)=='silty clay loam') then
       SoilCohesion2 = 8.0
    else if (trim(txttype)=='clay loam') then
       SoilCohesion2 = 8.0
    else if (trim(txttype)=='sandy clay loam') then
       SoilCohesion2 = 6.0
    else if (trim(txttype)=='silt') then
       SoilCohesion2 = 4.0
    else if (trim(txttype)=='silt loam') then
       SoilCohesion2 = 8.0
    else if (trim(txttype)=='loam') then
       SoilCohesion2 = 8.0
    else if (trim(txttype)=='sandy loam') then
       SoilCohesion2 = 6.0
    else if (trim(txttype)=='loamy sand') then
       SoilCohesion2 = 2.0
    else if (trim(txttype)=='sand') then
       SoilCohesion2 = 2.0
    end if

  end function SoilCohesion2

  !------------------------------------------------------------------------------
  real(r8) function SoilVanGenuchten(stxt)
    !
    ! !DESCRIPTION:
    ! soil van genuchten parameter (Schaap et al., 2001)
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: stxt(4)  ! [clay, silt, sand, gravel] in fraction 
    !
    ! !LOCAL VARIABLES:
    character(len=32) :: txttype
    !------------------------------------------------------------------------------

    txttype = SoilTextureType(stxt)
    if (trim(txttype)=='clay') then
       SoilVanGenuchten = 1.496
    else if (trim(txttype)=='silty clay') then
       SoilVanGenuchten = 1.622
    else if (trim(txttype)=='sandy clay') then
       SoilVanGenuchten = 3.342
    else if (trim(txttype)=='silty clay loam') then
       SoilVanGenuchten = 0.840
    else if (trim(txttype)=='clay loam') then
       SoilVanGenuchten = 1.581
    else if (trim(txttype)=='sandy clay loam') then
       SoilVanGenuchten = 2.109
    else if (trim(txttype)=='silt') then
       SoilVanGenuchten = 0.658
    else if (trim(txttype)=='silt loam') then
       SoilVanGenuchten = 0.506
    else if (trim(txttype)=='loam') then
       SoilVanGenuchten = 1.112
    else if (trim(txttype)=='sandy loam') then
       SoilVanGenuchten = 2.667
    else if (trim(txttype)=='loamy sand') then
       SoilVanGenuchten = 3.475
    else if (trim(txttype)=='sand') then
       SoilVanGenuchten = 3.524
    end if

  end function SoilVanGenuchten

  !------------------------------------------------------------------------------
  real(r8) function SoilPoreShape(stxt)
    !
    ! !DESCRIPTION:
    ! soil van genuchten parameter (Schaap et al., 2001)
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: stxt(4)  ! [clay, silt, sand, gravel] in fraction 
    !
    ! !LOCAL VARIABLES:
    character(len=32) :: txttype
    !------------------------------------------------------------------------------

    txttype = SoilTextureType(stxt)
    if (trim(txttype)=='clay') then
       SoilPoreShape = 1.253
    else if (trim(txttype)=='silty clay') then
       SoilPoreShape = 1.321
    else if (trim(txttype)=='sandy clay') then
       SoilPoreShape = 1.208
    else if (trim(txttype)=='silty clay loam') then
       SoilPoreShape = 1.521
    else if (trim(txttype)=='clay loam') then
       SoilPoreShape = 1.416
    else if (trim(txttype)=='sandy clay loam') then
       SoilPoreShape = 1.331
    else if (trim(txttype)=='silt') then
       SoilPoreShape = 1.679
    else if (trim(txttype)=='silt loam') then
       SoilPoreShape = 1.663
    else if (trim(txttype)=='loam') then
       SoilPoreShape = 1.472
    else if (trim(txttype)=='sandy loam') then
       SoilPoreShape = 1.449
    else if (trim(txttype)=='loamy sand') then
       SoilPoreShape = 1.746
    else if (trim(txttype)=='sand') then
       SoilPoreShape = 3.177
    end if

  end function SoilPoreShape

end module SedYieldMod
