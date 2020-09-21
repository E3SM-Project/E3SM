module SedYieldMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate the sediment flux caused by soil erosion based on the 
  ! improved Morgan model (Tan et al., 2017 & 2018) 
  !
  use shr_const_mod     , only : T0 => SHR_CONST_TKFRZ
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use abortutils        , only : endrun
  use decompMod         , only : bounds_type
  use elm_varcon        , only : grav, denh2o, rpi
  use clm_varpar        , only : mxpft, nlevsno, max_patch_per_col
  use clm_varpar        , only : nlevslp
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
  use TopounitDataType  , only : top_as, top_af ! Atmospheric state and flux variables
  use ColumnDataType    , only : col_ws, col_wf
  use VegetationDataType, only : veg_wf
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: SoilErosion          ! Calculate hillslope sediment flux
  ! !MODULE CONSTANTS
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine SoilErosion (bounds, num_soilc, filter_soilc, &
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
    use pftvarcon       , only : nc4_grass
    use spmdMod         , only : masterproc
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_soilc       ! number of column soil points in column filter
    integer                  , intent(in)    :: filter_soilc(:) ! column filter for soil points
    type(atm2lnd_type)       , intent(in)    :: atm2lnd_vars
    type(CanopyState_type)   , intent(in)    :: canopystate_vars
    type(soilstate_type)     , intent(in)    :: soilstate_vars
    type(waterstate_type)    , intent(in)    :: waterstate_vars
    type(waterflux_type)     , intent(in)    :: waterflux_vars
    type(sedflux_type)       , intent(inout) :: sedflux_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: c, fc, p, t, l, j                          ! indices
    integer  :: dtime                                      ! timestep size [seconds]
    real(r8) :: sinslp, fslp, factor_slp                   ! topo gradient
    real(r8) :: gndbare, vegcc, fbare                      ! veg cover factors
    real(r8) :: gndbare_crp                                ! crop veg cover factor
    real(r8) :: fungrvl                                    ! ground uncovered by gravel
    real(r8) :: K, COH                                     ! soil erodibility
    real(r8) :: Qs, Ptot, Ie, Dl                           ! water fluxes
    real(r8) :: Tc, Es_Q, Es_P, KE_DT, KE_LD               ! temporaries 
    real(r8) :: Es_Qcrp, Es_Pcrp                           ! cropland temporaries
    real(r8) :: stxt(4)                                    ! soil texture including gravel
    character(len=32) :: subname = 'SoilErosion'           ! subroutine name
    !-----------------------------------------------------------------------

    associate(                                                        &
         forc_rain        =>    top_af%rain                         , & ! Input: [real(r8) (:) ] rain rate [mm/s]
         forc_t           =>    top_as%tbot                         , & ! Input: [real(r8) (:) ] atmospheric temperature (Kelvin) 

         tlai             =>    canopystate_vars%tlai_patch         , & ! Input: [real(r8) (:) ] LAI
         hbot             =>    canopystate_vars%hbot_patch         , & ! Input: [real(r8) (:) ] canopy bottom (m)
         htop             =>    canopystate_vars%htop_patch         , & ! Input: [real(r8) (:) ] canopy top (m) 

         hslp_p10         =>    col_pp%hslp_p10                     , & ! Input: [real(r8) (:,:) ] hillslope gradient percentiles
    
         bd               =>    soilstate_vars%bd_col               , & ! Input: [real(r8) (:,:) ] soil bulk density (kg/m3)
         fsand            =>    soilstate_vars%cellsand_col         , & ! Input: [real(r8) (:,:) ] sand percentage
         fclay            =>    soilstate_vars%cellclay_col         , & ! Input: [real(r8) (:,:) ] clay percentage
         fgrvl            =>    soilstate_vars%cellgrvl_col         , & ! Input: [real(r8) (:,:) ] gravel percentage

         hsnow            =>    col_ws%snow_depth                   , & ! Input: [real(r8) (:) ] snow depth (m)

         pfactor          =>    sedflux_vars%pfactor_col            , & ! Input: [real(r8) (:) ] rainfall-driven erosion scaling factor
         qfactor          =>    sedflux_vars%qfactor_col            , & ! Input: [real(r8) (:) ] runoff-driven erosion scaling factor
         tfactor          =>    sedflux_vars%tfactor_col            , & ! Input: [real(r8) (:) ] transport capacity scaling factor 

         qflx_surf        =>    col_wf%qflx_surf                    , & ! Input: [real(r8) (:) ] surface runoff (mm/s)
         qflx_dirct_rain  =>    veg_wf%qflx_dirct_rain              , & ! Input: [real(r8) (:) ] direct throughfall rain (mm/s)
         qflx_leafdrip    =>    veg_wf%qflx_leafdrip                , & ! Input: [real(r8) (:) ] leaf rain drip (mm/s)

         flx_p_ero        =>    sedflux_vars%sed_p_ero_col          , & ! Output: [real(r8) (:) ] sed detached by rainfall (kg/m2/s)
         flx_q_ero        =>    sedflux_vars%sed_q_ero_col          , & ! Output: [real(r8) (:) ] sed detached by runoff (kg/m2/s)
         flx_sed_ero      =>    sedflux_vars%sed_ero_col            , & ! Output: [real(r8) (:) ] total detachment (kg/m2/s)
         flx_sed_crop_ero =>    sedflux_vars%sed_crop_ero_col       , & ! Output: [real(r8) (:) ] cropland detachment (kg/m2/s)
         flx_sed_yld      =>    sedflux_vars%sed_yld_col              & ! Output: [real(r8) (:) ] sed flux to inland waters (kg/m2/s)
         )

         ! Get time step
         dtime = get_step_size()

         ! nolakec or other col filters
         do fc = 1, num_soilc
            c = filter_soilc(fc)
            l = col_pp%landunit(c)
            t = col_pp%topounit(c)

            ! initialization
            flx_p_ero(c)          = 0._r8
            flx_q_ero(c)          = 0._r8
            flx_sed_ero(c)        = 0._r8
            flx_sed_crop_ero(c)   = 0._r8
            flx_sed_yld(c)        = 0._r8

            ! check landunit type and ground covered by snow/ice
            if ( lun_pp%itype(l)/=istsoil .and. lun_pp%itype(l)/=istcrop ) then
               cycle
            else if (hsnow(c)>0) then
               cycle
            end if

            ! soil detachment by rainfall
            stxt = (/fclay(c,1), 100._r8-fclay(c,1)-fsand(c,1), fsand(c,1), &
                fgrvl(c,1)/)  
            K = SoilDetachability(stxt)
            COH = SoilCohesion(stxt)
            Es_P = 0._r8    ! detachment by throughfall + leap drip
            Es_Pcrp = 0._r8 ! cropland detachment by throughfall + leap drip
            if (forc_t(t)>T0 .and. forc_rain(t)>0._r8) then
               fungrvl = 1._r8 - 0.01_r8 * fgrvl(c,1)
               do p = col_pp%pfti(c), col_pp%pftf(c)
                  if (veg_pp%active(p) .and. veg_pp%wtcol(p)>0._r8) then
                     ! throughfall power
                     Ptot = qflx_dirct_rain(p) * dtime   ! mm
                     Ie = 3.6e3_r8 * forc_rain(t)           ! mm/hr
                     KE_DT = Ptot * fungrvl * max(0._r8,8.95_r8+8.44_r8*log10(Ie))
                     ! leaf drip power
                     Dl = max(0._r8, qflx_leafdrip(p)*dtime)      ! mm
                     KE_LD = max(0._r8,15.8_r8*sqrt(0.5_r8*(htop(p)+hbot(p)))-5.87_r8) * &
                        fungrvl * Dl
                     Es_P = Es_P + pfactor(c) * veg_pp%wtcol(p) * K * (KE_DT+KE_LD)
                     ! For crop veg types
                     if( veg_pp%itype(p) > nc4_grass )then
                        Es_Pcrp = Es_Pcrp + pfactor(c) * veg_pp%wtcol(p) * K * (KE_DT+KE_LD)
                     end if 
                  end if
               end do
            end if
            Es_P = 1e-7_r8 / 8.64_r8 * Es_P        ! kg/m2/s
            Es_Pcrp = 1e-7_r8 / 8.64_r8 * Es_Pcrp  ! kg/m2/s

            ! soil detachment by runoff
            gndbare = 0._r8
            gndbare_crp = 0._r8
            vegcc = 0._r8
            do p = col_pp%pfti(c), col_pp%pftf(c)
               if (veg_pp%active(p) .and. veg_pp%wtcol(p)>0._r8) then
                  fbare = exp( -gcpsi(veg_pp%itype(p)) * tlai(p) )
                  gndbare = gndbare + veg_pp%wtcol(p) * fbare
                  if ( veg_pp%itype(p) > nc4_grass ) then
                     gndbare_crp = gndbare_crp + veg_pp%wtcol(p) * fbare
                  end if
                  vegcc = vegcc + veg_pp%wtcol(p) * pftcc(veg_pp%itype(p))
               end if
            end do
            
            Es_Q = 0._r8
            Es_Qcrp = 0._r8
            Tc = 0._r8
            if (qflx_surf(c)>0._r8) then
               Qs = 8.64e4_r8 * qflx_surf(c)  ! mm/d
               fslp = 1.0_r8 / DBLE(nlevslp-1)
               factor_slp = 0._r8
               do j = 1, nlevslp-1
                  sinslp = sin(atan(max(0.5_r8*(hslp_p10(c,j)+hslp_p10(c,j+1)),1e-4_r8)))
                  factor_slp = factor_slp + fslp * sinslp
               end do
               Es_Q = 19.1_r8 * qfactor(c) * 2.0_r8 / COH * factor_slp * &
                  gndbare * Qs**1.5_r8
               Es_Qcrp = 19.1_r8 * qfactor(c) * 2.0_r8 / COH * factor_slp * &
                  gndbare_crp * Qs**1.5_r8
               Tc = 19.1_r8 * tfactor(c) * factor_slp * vegcc * Qs**2.0_r8
            end if
            Es_Q = 1e-7_r8 / 8.64_r8 * Es_Q        ! kg/m2/s
            Es_Qcrp = 1e-7_r8 / 8.64_r8 * Es_Qcrp  ! kg/m2/s
            Tc = 1e-7_r8 / 8.64_r8 * Tc            ! kg/m2/s

            ! assign flux values
            flx_p_ero(c) = flx_p_ero(c) + Es_P
            flx_q_ero(c) = flx_q_ero(c) + Es_Q
            flx_sed_ero(c) = flx_sed_ero(c) + Es_P + Es_Q 
            flx_sed_crop_ero(c) = flx_sed_crop_ero(c) + Es_Pcrp + Es_Qcrp
            flx_sed_yld(c) = flx_sed_yld(c) + min(Es_P+Es_Q, Tc)
         end do
    
    end associate

  end subroutine SoilErosion

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
    if ( silt + 1.5_r8*clay < 15._r8 ) then
       SoilTextureType = 'sand'
    else if ( silt + 2.0_r8*clay < 30._r8 ) then
       SoilTextureType = 'loamy sand'
    else if ( ((clay>=7._r8 .and. clay<20._r8) .and. sand>52._r8) .or. &
            (clay<7._r8 .and. silt<50._r8) ) then
       SoilTextureType = 'sandy loam'
    else if ( (clay>=7._r8 .and. clay<27._r8) .and. & 
            (silt>=28._r8 .and. silt<50._r8) .and. sand<=52._r8 ) then
       SoilTextureType = 'loam'
    else if ( (silt>=50._r8 .and. (clay>=12._r8 .and. clay<27._r8)) .or. &
            ((silt>=50._r8 .and. silt<80._r8) .and. clay<12._r8) ) then
       SoilTextureType = 'silt loam'
    else if ( silt>=80._r8 .and. clay<12._r8 ) then
       SoilTextureType = 'silt'
    else if ( (clay>=20._r8 .and. clay<35._r8) .and. &
            silt<28._r8 .and. sand>45._r8 ) then
       SoilTextureType = 'sandy clay loam'
    else if ( (clay>=27._r8 .and. clay<40._r8) .and. &
            (sand>20._r8 .and. sand<=45._r8) ) then
       SoilTextureType = 'clay loam'
    else if ( (clay>=27._r8 .and. clay<40._r8) .and. sand<=20._r8 ) then
       SoilTextureType = 'silty clay loam'
    else if ( clay>=35._r8 .and. sand>45._r8 ) then
       SoilTextureType = 'sandy clay'
    else if ( clay>=40._r8 .and. silt>=40._r8 ) then
       SoilTextureType = 'silty clay'
    else if ( clay>=40._r8 .and. sand<=45._r8 .and. silt<40._r8 ) then
       SoilTextureType = 'clay'
    else
       call endrun(msg=' ERROR: soil clay, silt and sand are out of bounds.'//&
            errMsg(__FILE__, __LINE__))
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
       SoilDetachability = 0.05_r8
    else if (trim(txttype)=='silty clay') then
       SoilDetachability = 0.5_r8
    else if (trim(txttype)=='sandy clay') then
       SoilDetachability = 0.3_r8
    else if (trim(txttype)=='silty clay loam') then
       SoilDetachability = 0.8_r8
    else if (trim(txttype)=='clay loam') then
       SoilDetachability = 0.7_r8
    else if (trim(txttype)=='sandy clay loam') then
       SoilDetachability = 0.1_r8
    else if (trim(txttype)=='silt') then
       SoilDetachability = 1.0_r8
    else if (trim(txttype)=='silt loam') then
       SoilDetachability = 0.9_r8
    else if (trim(txttype)=='loam') then
       SoilDetachability = 0.8_r8
    else if (trim(txttype)=='sandy loam') then
       SoilDetachability = 0.7_r8
    else if (trim(txttype)=='loamy sand') then
       SoilDetachability = 0.3_r8
    else if (trim(txttype)=='sand') then
       SoilDetachability = 1.9_r8
    else
       call endrun(msg=' ERROR: no soil texture type is found.'//&
            errMsg(__FILE__, __LINE__)) 
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
       SoilCohesion = 12.0_r8
    else if (trim(txttype)=='silty clay') then
       SoilCohesion = 10.0_r8
    else if (trim(txttype)=='sandy clay') then
       SoilCohesion = 10.0_r8
    else if (trim(txttype)=='silty clay loam') then
       SoilCohesion = 9.0_r8
    else if (trim(txttype)=='clay loam') then
       SoilCohesion = 10.0_r8
    else if (trim(txttype)=='sandy clay loam') then
       SoilCohesion = 3.0_r8
    else if (trim(txttype)=='silt') then
       SoilCohesion = 3.0_r8
    else if (trim(txttype)=='silt loam') then
       SoilCohesion = 3.0_r8
    else if (trim(txttype)=='loam') then
       SoilCohesion = 3.0_r8
    else if (trim(txttype)=='sandy loam') then
       SoilCohesion = 2.0_r8
    else if (trim(txttype)=='loamy sand') then
       SoilCohesion = 2.0_r8
    else if (trim(txttype)=='sand') then
       SoilCohesion = 2.0_r8
    else
       call endrun(msg=' ERROR: no soil texture type is found.'//&
            errMsg(__FILE__, __LINE__))
    end if

  end function SoilCohesion

end module SedYieldMod
