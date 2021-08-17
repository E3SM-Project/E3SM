module UrbanRadiationMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate solar and longwave radiation, and turbulent fluxes for urban landunit
  !
  ! !USES:
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_sys_mod       , only : shr_sys_flush
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use decompMod         , only : bounds_type
  use elm_varpar        , only : numrad
  use elm_varcon        , only : isecspday, degpsec, namel
  use elm_varctl        , only : iulog
  use abortutils        , only : endrun
  use UrbanParamsType   , only : urbanparams_type
  use atm2lndType       , only : atm2lnd_type
  use SolarAbsorbedType , only : solarabs_type
  use SurfaceAlbedoType , only : surfalb_type
  use UrbanParamsType   , only : urbanparams_type
  use EnergyFluxType    , only : energyflux_type
  use TopounitDataType  , only : top_af
  use LandunitType      , only : lun_pp
  use ColumnType        , only : col_pp
  use ColumnDataType    , only : col_es, col_ws
  use VegetationType    , only : veg_pp
  use VegetationDataType, only : veg_es, veg_ef
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: UrbanRadiation    ! Urban physics - radiative fluxes
  !
  ! PRIVATE MEMBER FUNCTIONS
  private :: net_longwave     ! Net longwave radiation for road and both walls in urban canyon
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine UrbanRadiation ( &
       num_nourbanl, filter_nourbanl                                  , &
       num_urbanl, filter_urbanl                                      , &
       atm2lnd_vars, urbanparams_vars, &
       solarabs_vars, surfalb_vars, energyflux_vars)
    !
    ! !DESCRIPTION:
    ! Solar fluxes absorbed and reflected by roof and canyon (walls, road).
    ! Also net and upward longwave fluxes.

    ! !USES:
      !$acc routine seq
    use elm_varcon          , only : spval, sb, tfrz
    use column_varcon       , only : icol_road_perv, icol_road_imperv
    use column_varcon       , only : icol_roof, icol_sunwall, icol_shadewall
    !
    ! !ARGUMENTS:
    integer                , intent(in)    :: num_nourbanl       ! number of non-urban landunits in clump
    integer                , intent(in)    :: filter_nourbanl(:) ! non-urban landunit filter
    integer                , intent(in)    :: num_urbanl         ! number of urban landunits in clump
    integer                , intent(in)    :: filter_urbanl(:)   ! urban landunit filter
    type(atm2lnd_type)     , intent(in)    :: atm2lnd_vars
    type(urbanparams_type) , intent(in)    :: urbanparams_vars
    type(solarabs_type)    , intent(inout) :: solarabs_vars
    type(surfalb_type)     , intent(in)    :: surfalb_vars
    type(energyflux_type)  , intent(inout) :: energyflux_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: fp,fl,p,c,l,t,g     ! indices
    integer  :: local_secp1         ! seconds into current date in local time
    real(r8) :: dtime               ! land model time step (sec)
    integer  :: year,month,day      ! temporaries (not used)
    integer  :: secs                ! seconds into current date

    real(r8), parameter :: mpe    = 1.e-06_r8 ! prevents overflow for division by zero
    real(r8), parameter :: snoem  = 0.97_r8   ! snow emissivity (should use value from Biogeophysics1)

    real(r8) :: lwnet_roof     ! net (outgoing-incoming) longwave radiation (per unit ground area), roof (W/m**2)
    real(r8) :: lwnet_improad  ! net (outgoing-incoming) longwave radiation (per unit ground area), impervious road (W/m**2)
    real(r8) :: lwnet_perroad  ! net (outgoing-incoming) longwave radiation (per unit ground area), pervious road (W/m**2)
    real(r8) :: lwnet_sunwall  ! net (outgoing-incoming) longwave radiation (per unit wall area), sunlit wall (W/m**2)
    real(r8) :: lwnet_shadewall! net (outgoing-incoming) longwave radiation (per unit wall area), shaded wall (W/m**2)
    real(r8) :: lwnet_canyon   ! net (outgoing-incoming) longwave radiation for canyon, per unit ground area (W/m**2)
    real(r8) :: lwup_roof      ! upward longwave radiation (per unit ground area), roof (W/m**2)
    real(r8) :: lwup_improad   ! upward longwave radiation (per unit ground area), impervious road (W/m**2)
    real(r8) :: lwup_perroad   ! upward longwave radiation (per unit ground area), pervious road (W/m**2)
    real(r8) :: lwup_sunwall   ! upward longwave radiation, (per unit wall area), sunlit wall (W/m**2)
    real(r8) :: lwup_shadewall ! upward longwave radiation, (per unit wall area), shaded wall (W/m**2)
    real(r8) :: lwup_canyon    ! upward longwave radiation for canyon, per unit ground area (W/m**2)
    real(r8) :: t_roof         ! roof temperature (K)
    real(r8) :: t_improad      ! imppervious road temperature (K)
    real(r8) :: t_perroad      ! pervious road temperature (K)
    real(r8) :: t_sunwall      ! sunlit wall temperature (K)
    real(r8) :: t_shadewall    ! shaded wall temperature (K)
    real(r8) :: lwdown         ! atmospheric downward longwave radiation (W/m**2)
    real(r8) :: em_roof_s      ! roof emissivity with snow effects
    real(r8) :: em_improad_s   ! impervious road emissivity with snow effects
    real(r8) :: em_perroad_s   ! pervious road emissivity with snow effects
    !-----------------------------------------------------------------------

    associate(                                               &
         ctype              =>    col_pp%itype                   , & ! Input:  [integer (:)    ]  column type
         coli               =>    lun_pp%coli                    , & ! Input:  [integer (:)    ]  beginning column index for landunit
         colf               =>    lun_pp%colf                    , & ! Input:  [integer (:)    ]  ending column index for landunit
         canyon_hwr         =>    lun_pp%canyon_hwr              , & ! Input:  [real(r8) (:)   ]  ratio of building height to street width
         wtroad_perv        =>    lun_pp%wtroad_perv             , & ! Input:  [real(r8) (:)   ]  weight of pervious road wrt total road

         forc_solad         =>    top_af%solad                   , & ! Input:  [real(r8) (:,:) ]  direct beam radiation  (vis=forc_sols , nir=forc_soll ) (W/m**2)
         forc_solai         =>    top_af%solai                   , & ! Input:  [real(r8) (:,:) ]  diffuse beam radiation (vis=forc_sols , nir=forc_soll ) (W/m**2)
         forc_solar         =>    top_af%solar                   , & ! Input:  [real(r8) (:)   ]  incident solar radiation (W/m**2)
         forc_lwrad         =>    top_af%lwrad                   , & ! Input:  [real(r8) (:)   ]  downward infrared (longwave) radiation (W/m**2)

         frac_sno           =>    col_ws%frac_sno               , & ! Input:  [real(r8) (:)   ]  fraction of ground covered by snow (0 to 1)

         t_ref2m            =>    veg_es%t_ref2m             , & ! Input:  [real(r8) (:)   ]  2 m height surface air temperature (K)
         t_grnd             =>    col_es%t_grnd              , & ! Input:  [real(r8) (:)   ]  ground temperature (K)

         em_roof            =>    urbanparams_vars%em_roof                   , & ! Input:  [real(r8) (:)   ]  roof emissivity
         em_improad         =>    urbanparams_vars%em_improad                , & ! Input:  [real(r8) (:)   ]  impervious road emissivity
         em_perroad         =>    urbanparams_vars%em_perroad                , & ! Input:  [real(r8) (:)   ]  pervious road emissivity
         em_wall            =>    urbanparams_vars%em_wall                   , & ! Input:  [real(r8) (:)   ]  wall emissivity

         albd               =>    surfalb_vars%albd_patch                    , & ! Input:  [real(r8) (:,:) ] pft surface albedo (direct)
         albi               =>    surfalb_vars%albi_patch                    , & ! Input:  [real(r8) (:,:) ] pft surface albedo (diffuse)

         sabs_roof_dir      =>    solarabs_vars%sabs_roof_dir_lun            , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by roof per unit ground area per unit incident flux
         sabs_roof_dif      =>    solarabs_vars%sabs_roof_dif_lun            , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by roof per unit ground area per unit incident flux
         sabs_sunwall_dir   =>    solarabs_vars%sabs_sunwall_dir_lun         , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by sunwall per unit wall area per unit incident flux
         sabs_sunwall_dif   =>    solarabs_vars%sabs_sunwall_dif_lun         , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by sunwall per unit wall area per unit incident flux
         sabs_shadewall_dir =>    solarabs_vars%sabs_shadewall_dir_lun       , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by shadewall per unit wall area per unit incident flux
         sabs_shadewall_dif =>    solarabs_vars%sabs_shadewall_dif_lun       , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by shadewall per unit wall area per unit incident flux
         sabs_improad_dir   =>    solarabs_vars%sabs_improad_dir_lun         , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by impervious road per unit ground area per unit incident flux
         sabs_improad_dif   =>    solarabs_vars%sabs_improad_dif_lun         , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by impervious road per unit ground area per unit incident flux
         sabs_perroad_dir   =>    solarabs_vars%sabs_perroad_dir_lun         , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by pervious road per unit ground area per unit incident flux
         sabs_perroad_dif   =>    solarabs_vars%sabs_perroad_dif_lun         , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by pervious road per unit ground area per unit incident flux
         sabg               =>    solarabs_vars%sabg_patch                   , & ! Output: [real(r8) (:)   ]  solar radiation absorbed by ground (W/m**2)
         sabv               =>    solarabs_vars%sabv_patch                   , & ! Output: [real(r8) (:)   ]  solar radiation absorbed by vegetation (W/m**2)
         fsa                =>    solarabs_vars%fsa_patch                    , & ! Output: [real(r8) (:)   ]  solar radiation absorbed (total) (W/m**2)
         fsa_u              =>    solarabs_vars%fsa_u_patch                  , & ! Output: [real(r8) (:)   ]  urban solar radiation absorbed (total) (W/m**2)

         eflx_lwrad_out     =>    veg_ef%eflx_lwrad_out       , & ! Output: [real(r8) (:)   ]  emitted infrared (longwave) radiation (W/m**2)
         eflx_lwrad_net     =>    veg_ef%eflx_lwrad_net       , & ! Output: [real(r8) (:)   ]  net infrared (longwave) rad (W/m**2) [+ = to atm]
         eflx_lwrad_net_u   =>    veg_ef%eflx_lwrad_net_u      & ! Output: [real(r8) (:)   ]  urban net infrared (longwave) rad (W/m**2) [+ = to atm]

         )

      ! Define fields that appear on the restart file for non-urban landunits

      do fl = 1,num_nourbanl
         l = filter_nourbanl(fl)
         sabs_roof_dir(l,:)      = spval
         sabs_roof_dif(l,:)      = spval
         sabs_sunwall_dir(l,:)   = spval
         sabs_sunwall_dif(l,:)   = spval
         sabs_shadewall_dir(l,:) = spval
         sabs_shadewall_dif(l,:) = spval
         sabs_improad_dir(l,:)   = spval
         sabs_improad_dif(l,:)   = spval
         sabs_perroad_dir(l,:)   = spval
         sabs_perroad_dif(l,:)   = spval
      end do

      ! Set input forcing fields
      do fl = 1,num_urbanl
         l = filter_urbanl(fl)
         t = lun_pp%topounit(l)
         g = lun_pp%gridcell(l)

         ! Need to set the following temperatures to some defined value even if it
         ! does not appear in the urban landunit for the net_longwave computation

         t_roof      = 19._r8 + tfrz
         t_sunwall   = 19._r8 + tfrz
         t_shadewall = 19._r8 + tfrz
         t_improad   = 19._r8 + tfrz
         t_perroad   = 19._r8 + tfrz

         ! Initial assignment of emissivity
         em_roof_s    = em_roof(l)
         em_improad_s = em_improad(l)
         em_perroad_s = em_perroad(l)

         ! Set urban temperatures and emissivity including snow effects.
         do c = coli(l),colf(l)
            if (ctype(c) == icol_roof )  then
               t_roof     = t_grnd(c)
               em_roof_s = em_roof(l)*(1._r8-frac_sno(c)) + snoem*frac_sno(c)
            else if (ctype(c) == icol_road_imperv) then
               t_improad   = t_grnd(c)
               em_improad_s = em_improad(l)*(1._r8-frac_sno(c)) + snoem*frac_sno(c)
            else if (ctype(c) == icol_road_perv  ) then
               t_perroad   = t_grnd(c)
               em_perroad_s = em_perroad(l)*(1._r8-frac_sno(c)) + snoem*frac_sno(c)
            else if (ctype(c) == icol_sunwall    ) then
               t_sunwall   = t_grnd(c)
            else if (ctype(c) == icol_shadewall  ) then
               t_shadewall = t_grnd(c)
            end if
         end do
         lwdown = forc_lwrad(t)

         ! Net longwave radiation for road and both walls in urban canyon allowing for multiple re-emission
         ! NOTE: is it better to do surfaces independently in separate routines?
         ! This would cut local variables in half and make net_longwave more readable
         if (num_urbanl > 0) then
           call net_longwave (l ,     &
                canyon_hwr(l),      &
                wtroad_perv(l),     &
                lwdown,          &
                em_roof_s,       &
                em_improad_s,    &
                em_perroad_s,    &
                em_wall(l),      &
                t_roof,          &
                t_improad,       &
                t_perroad,       &
                t_sunwall,       &
                t_shadewall,     &
                lwnet_roof,      &
                lwnet_improad,   &
                lwnet_perroad,   &
                lwnet_sunwall,   &
                lwnet_shadewall, &
                lwnet_canyon,    &
                lwup_roof,       &
                lwup_improad,    &
                lwup_perroad,    &
                lwup_sunwall,    &
                lwup_shadewall,  &
                lwup_canyon,     &
                urbanparams_vars) !why pass entire urbanparams here?

         end if

      ! Determine variables needed for history output and communication with atm
      ! Loop over urban patches in clump
      !NOTE:  Rewrote this loop to have all landunit calcs within one overall loop.
      !       This saves local memory allocations.

      do p = lun_pp%pfti(l), lun_pp%pftf(l)
        !!are all patches on a active lanunit also active?
          if(.not. veg_pp%active(p)) then
            cycle
          end if

          g = veg_pp%gridcell(p)
          t = veg_pp%topounit(p)
          c = veg_pp%column(p)

         ! Solar absorbed and longwave out and net
         ! per unit ground area (roof, road) and per unit wall area (sunwall, shadewall)
         ! Each urban pft has its own column - this is used in the logic below

         if (ctype(c) == icol_roof) then
            eflx_lwrad_out(p) = lwup_roof
            eflx_lwrad_net(p) = lwnet_roof
            eflx_lwrad_net_u(p) = lwnet_roof
            sabg(p) = sabs_roof_dir(l,1)*forc_solad(t,1) + &
                 sabs_roof_dif(l,1)*forc_solai(t,1) + &
                 sabs_roof_dir(l,2)*forc_solad(t,2) + &
                 sabs_roof_dif(l,2)*forc_solai(t,2)

         else if (ctype(c) == icol_sunwall) then
            eflx_lwrad_out(p)   = lwup_sunwall
            eflx_lwrad_net(p)   = lwnet_sunwall
            eflx_lwrad_net_u(p) = lwnet_sunwall
            sabg(p) = sabs_sunwall_dir(l,1)*forc_solad(t,1) + &
                 sabs_sunwall_dif(l,1)*forc_solai(t,1) + &
                 sabs_sunwall_dir(l,2)*forc_solad(t,2) + &
                 sabs_sunwall_dif(l,2)*forc_solai(t,2)

         else if (ctype(c) == icol_shadewall) then
            eflx_lwrad_out(p)   = lwup_shadewall
            eflx_lwrad_net(p)   = lwnet_shadewall
            eflx_lwrad_net_u(p) = lwnet_shadewall
            sabg(p) = sabs_shadewall_dir(l,1)*forc_solad(t,1) + &
                 sabs_shadewall_dif(l,1)*forc_solai(t,1) + &
                 sabs_shadewall_dir(l,2)*forc_solad(t,2) + &
                 sabs_shadewall_dif(l,2)*forc_solai(t,2)

         else if (ctype(c) == icol_road_perv) then
            eflx_lwrad_out(p)   = lwup_perroad
            eflx_lwrad_net(p)   = lwnet_perroad
            eflx_lwrad_net_u(p) = lwnet_perroad
            sabg(p) = sabs_perroad_dir(l,1)*forc_solad(t,1) + &
                 sabs_perroad_dif(l,1)*forc_solai(t,1) + &
                 sabs_perroad_dir(l,2)*forc_solad(t,2) + &
                 sabs_perroad_dif(l,2)*forc_solai(t,2)

         else if (ctype(c) == icol_road_imperv) then
            eflx_lwrad_out(p)   = lwup_improad
            eflx_lwrad_net(p)   = lwnet_improad
            eflx_lwrad_net_u(p) = lwnet_improad
            sabg(p) = sabs_improad_dir(l,1)*forc_solad(t,1) + &
                 sabs_improad_dif(l,1)*forc_solai(t,1) + &
                 sabs_improad_dir(l,2)*forc_solad(t,2) + &
                 sabs_improad_dif(l,2)*forc_solai(t,2)
         end if

         sabv(p)   = 0._r8
         fsa(p)    = sabv(p) + sabg(p)
         fsa_u(p)  = fsa(p)

      end do ! end loop over urban patches

    end do ! end loop over urban landunit

    end associate

  end subroutine UrbanRadiation

  !-----------------------------------------------------------------------
  subroutine net_longwave (l, canyon_hwr, wtroad_perv                , &
       lwdown, em_roof, em_improad, em_perroad, em_wall              , &
       t_roof,  t_improad, t_perroad, t_sunwall, t_shadewall         , &
       lwnet_roof, lwnet_improad, lwnet_perroad, lwnet_sunwall, lwnet_shadewall, lwnet_canyon , &
       lwup_roof, lwup_improad, lwup_perroad, lwup_sunwall, lwup_shadewall, lwup_canyon, &
       urbanparams_vars)
    !
    ! !DESCRIPTION:
    ! Net longwave radiation for road and both walls in urban canyon allowing for
    ! multiple reflection. Also net longwave radiation for urban roof.
    !
    ! !USES:
      !$acc routine seq
    use elm_varcon , only : sb
    !
    ! !ARGUMENTS:
    integer, value, intent(in) :: l
    real(r8), intent(in)  :: canyon_hwr     ! ratio of building height to street width [landunit]
    real(r8), intent(in)  :: wtroad_perv    ! weight of pervious road wrt total road [landunit]

    real(r8), intent(in)  :: lwdown         ! atmospheric longwave radiation (W/m**2) [landunit]
    real(r8), intent(in)  :: em_roof        ! roof emissivity [landunit]
    real(r8), intent(in)  :: em_improad     ! impervious road emissivity [landunit]
    real(r8), intent(in)  :: em_perroad     ! pervious road emissivity [landunit]
    real(r8), intent(in)  :: em_wall        ! wall emissivity [landunit]

    real(r8), intent(in)  :: t_roof         ! roof temperature (K) [landunit]
    real(r8), intent(in)  :: t_improad      ! impervious road temperature (K) [landunit]
    real(r8), intent(in)  :: t_perroad      ! ervious road temperature (K) [landunit]
    real(r8), intent(in)  :: t_sunwall      ! sunlit wall temperature (K) [landunit]
    real(r8), intent(in)  :: t_shadewall    ! shaded wall temperature (K) [landunit]

    real(r8), intent(out) :: lwnet_roof     ! net (outgoing-incoming) longwave radiation, roof (W/m**2) [landunit]
    real(r8), intent(out) :: lwnet_improad  ! net (outgoing-incoming) longwave radiation, impervious road (W/m**2) [landunit]
    real(r8), intent(out) :: lwnet_perroad  ! net (outgoing-incoming) longwave radiation, pervious road (W/m**2) [landunit]
    real(r8), intent(out) :: lwnet_sunwall  ! net (outgoing-incoming) longwave radiation (per unit wall area), sunlit wall (W/m**2) [landunit]
    real(r8), intent(out) :: lwnet_shadewall! net (outgoing-incoming) longwave radiation (per unit wall area), shaded wall (W/m**2) [landunit]
    real(r8), intent(out) :: lwnet_canyon   ! net (outgoing-incoming) longwave radiation for canyon, per unit ground area (W/m**2) [landunit]

    real(r8), intent(out) :: lwup_roof      ! upward longwave radiation, roof (W/m**2) [landunit]
    real(r8), intent(out) :: lwup_improad   ! upward longwave radiation, impervious road (W/m**2) [landunit]
    real(r8), intent(out) :: lwup_perroad   ! upward longwave radiation, pervious road (W/m**2) [landunit]
    real(r8), intent(out) :: lwup_sunwall   ! upward longwave radiation (per unit wall area), sunlit wall (W/m**2) [landunit]
    real(r8), intent(out) :: lwup_shadewall ! upward longwave radiation (per unit wall area), shaded wall (W/m**2) [landunit]
    real(r8), intent(out) :: lwup_canyon    ! upward longwave radiation for canyon, per unit ground area (W/m**2) [landunit]
    !
    type(urbanparams_type) , intent(in) :: urbanparams_vars
    !
    ! !LOCAL VARIABLES:
    real(r8) :: lwdown_road          ! atmospheric longwave radiation for total road (W/m**2)
    real(r8) :: lwdown_sunwall       ! atmospheric longwave radiation (per unit wall area) for sunlit wall (W/m**2)
    real(r8) :: lwdown_shadewall     ! atmospheric longwave radiation (per unit wall area) for shaded wall (W/m**2)
    real(r8) :: lwtot                ! incoming longwave radiation (W/m**2)

    real(r8) :: improad_a            ! absorbed longwave for improad (W/m**2)
    real(r8) :: improad_r            ! reflected longwave for improad (W/m**2)
    real(r8) :: improad_r_sky        ! improad_r to sky (W/m**2)
    real(r8) :: improad_r_sunwall    ! improad_r to sunlit wall (W/m**2)
    real(r8) :: improad_r_shadewall  ! improad_r to shaded wall (W/m**2)
    real(r8) :: improad_e            ! emitted longwave for improad (W/m**2)
    real(r8) :: improad_e_sky        ! improad_e to sky (W/m**2)
    real(r8) :: improad_e_sunwall    ! improad_e to sunlit wall (W/m**2)
    real(r8) :: improad_e_shadewall  ! improad_e to shaded wall (W/m**2)

    real(r8) :: perroad_a            ! absorbed longwave for perroad (W/m**2)
    real(r8) :: perroad_r            ! reflected longwave for perroad (W/m**2)
    real(r8) :: perroad_r_sky        ! perroad_r to sky (W/m**2)
    real(r8) :: perroad_r_sunwall    ! perroad_r to sunlit wall (W/m**2)
    real(r8) :: perroad_r_shadewall  ! perroad_r to shaded wall (W/m**2)
    real(r8) :: perroad_e            ! emitted longwave for perroad (W/m**2)
    real(r8) :: perroad_e_sky        ! perroad_e to sky (W/m**2)
    real(r8) :: perroad_e_sunwall    ! perroad_e to sunlit wall (W/m**2)
    real(r8) :: perroad_e_shadewall  ! perroad_e to shaded wall (W/m**2)

    real(r8) :: road_a               ! absorbed longwave for total road (W/m**2)
    real(r8) :: road_r               ! reflected longwave for total road (W/m**2)
    real(r8) :: road_r_sky           ! total road_r to sky (W/m**2)
    real(r8) :: road_r_sunwall       ! total road_r to sunlit wall (W/m**2)
    real(r8) :: road_r_shadewall     ! total road_r to shaded wall (W/m**2)
    real(r8) :: road_e               ! emitted longwave for total road (W/m**2)
    real(r8) :: road_e_sky           ! total road_e to sky (W/m**2)
    real(r8) :: road_e_sunwall       ! total road_e to sunlit wall (W/m**2)
    real(r8) :: road_e_shadewall     ! total road_e to shaded wall (W/m**2)

    real(r8) :: sunwall_a            ! absorbed longwave (per unit wall area) for sunlit wall (W/m**2)
    real(r8) :: sunwall_r            ! reflected longwave (per unit wall area) for sunlit wall (W/m**2)
    real(r8) :: sunwall_r_sky        ! sunwall_r to sky (W/m**2)
    real(r8) :: sunwall_r_road       ! sunwall_r to road (W/m**2)
    real(r8) :: sunwall_r_shadewall  ! sunwall_r to opposing (shaded) wall (W/m**2)
    real(r8) :: sunwall_e            ! emitted longwave (per unit wall area) for sunlit wall (W/m**2)
    real(r8) :: sunwall_e_sky        ! sunwall_e to sky (W/m**2)
    real(r8) :: sunwall_e_road       ! sunwall_e to road (W/m**2)
    real(r8) :: sunwall_e_shadewall  ! sunwall_e to opposing (shaded) wall (W/m**2)

    real(r8) :: shadewall_a          ! absorbed longwave (per unit wall area) for shaded wall (W/m**2)
    real(r8) :: shadewall_r          ! reflected longwave (per unit wall area) for shaded wall (W/m**2)
    real(r8) :: shadewall_r_sky      ! shadewall_r to sky (W/m**2)
    real(r8) :: shadewall_r_road     ! shadewall_r to road (W/m**2)
    real(r8) :: shadewall_r_sunwall  ! shadewall_r to opposing (sunlit) wall (W/m**2)
    real(r8) :: shadewall_e          ! emitted longwave (per unit wall area) for shaded wall (W/m**2)
    real(r8) :: shadewall_e_sky      ! shadewall_e to sky (W/m**2)
    real(r8) :: shadewall_e_road     ! shadewall_e to road (W/m**2)
    real(r8) :: shadewall_e_sunwall  ! shadewall_e to opposing (sunlit) wall (W/m**2)
    integer  :: fl,iter            ! indices
    integer, parameter  :: n = 50    ! number of interations
    real(r8) :: crit                 ! convergence criterion (W/m**2)
    real(r8) :: err                  ! energy conservation error (W/m**2)
    real(r8) :: wtroad_imperv        ! weight of impervious road wrt total road
    !-----------------------------------------------------------------------

    associate(                             &
         vf_sr => urbanparams_vars%vf_sr , & ! Input:  [real(r8) (:)]  view factor of sky for road
         vf_wr => urbanparams_vars%vf_wr , & ! Input:  [real(r8) (:)]  view factor of one wall for road
         vf_sw => urbanparams_vars%vf_sw , & ! Input:  [real(r8) (:)]  view factor of sky for one wall
         vf_rw => urbanparams_vars%vf_rw , & ! Input:  [real(r8) (:)]  view factor of road for one wall
         vf_ww => urbanparams_vars%vf_ww   & ! Input:  [real(r8) (:)]  view factor of opposing wall for one wall
         )

      ! Calculate impervious road

      wtroad_imperv = 1._r8 - wtroad_perv

      ! atmospheric longwave radiation incident on walls and road in urban canyon.
      ! check for conservation (need to convert wall fluxes to ground area).
      ! lwdown (from atmosphere) = lwdown_road + (lwdown_sunwall + lwdown_shadewall)*canyon_hwr

      lwdown_road      = lwdown * vf_sr(l)
      lwdown_sunwall   = lwdown * vf_sw(l)
      lwdown_shadewall = lwdown * vf_sw(l)

      err = lwdown - (lwdown_road + (lwdown_shadewall + lwdown_sunwall)*canyon_hwr)
      if (abs(err) > 0.10_r8 ) then
#ifndef _OPENACC
         write(iulog,*) 'urban incident atmospheric longwave radiation balance error',err
         write(iulog,*) 'l          = ',l
         write(iulog,*) 'lwdown     = ',lwdown
         write(iulog,*) 'vf_sr      = ',vf_sr(l)
         write(iulog,*) 'vf_sw      = ',vf_sw(l)
         write(iulog,*) 'canyon_hwr = ',canyon_hwr
         write(iulog,*) 'elm model is stopping'
         call endrun(decomp_index=l, elmlevel=namel, msg=errmsg(__FILE__, __LINE__))
#endif
      endif

      ! initial absorption, reflection, and emission for road and both walls.
      ! distribute reflected and emitted radiation to sky, road, and walls according
      ! to appropriate view factor. radiation reflected to road and walls will
      ! undergo multiple reflections within the canyon.

      road_a              = 0.0_r8
      road_r              = 0.0_r8
      road_e              = 0.0_r8
      improad_a            =     em_improad   * lwdown_road
      improad_r            = (1._r8-em_improad ) * lwdown_road
      improad_r_sky        = improad_r  * vf_sr(l)
      improad_r_sunwall    = improad_r  * vf_wr(l)
      improad_r_shadewall  = improad_r  * vf_wr(l)
      improad_e            = em_improad  * sb * (t_improad **4)
      improad_e_sky        = improad_e  * vf_sr(l)
      improad_e_sunwall    = improad_e  * vf_wr(l)
      improad_e_shadewall  = improad_e  * vf_wr(l)
      road_a               = road_a  + improad_a *wtroad_imperv
      road_r               = road_r  + improad_r *wtroad_imperv
      road_e               = road_e  + improad_e *wtroad_imperv

      perroad_a            =     em_perroad   * lwdown_road
      perroad_r            = (1._r8-em_perroad ) * lwdown_road
      perroad_r_sky        = perroad_r  * vf_sr(l)
      perroad_r_sunwall    = perroad_r  * vf_wr(l)
      perroad_r_shadewall  = perroad_r  * vf_wr(l)
      perroad_e            = em_perroad * sb * (t_perroad **4)
      perroad_e_sky        = perroad_e  * vf_sr(l)
      perroad_e_sunwall    = perroad_e  * vf_wr(l)
      perroad_e_shadewall  = perroad_e  * vf_wr(l)
      road_a               = road_a  + perroad_a *wtroad_perv
      road_r               = road_r  + perroad_r *wtroad_perv
      road_e               = road_e  + perroad_e *wtroad_perv

      road_r_sky           = road_r  * vf_sr(l)
      road_r_sunwall       = road_r  * vf_wr(l)
      road_r_shadewall     = road_r  * vf_wr(l)
      road_e_sky           = road_e  * vf_sr(l)
      road_e_sunwall       = road_e  * vf_wr(l)
      road_e_shadewall     = road_e  * vf_wr(l)

      sunwall_a            = em_wall  * lwdown_sunwall
      sunwall_r            = (1._r8-em_wall ) * lwdown_sunwall
      sunwall_r_sky        = sunwall_r  * vf_sw(l)
      sunwall_r_road       = sunwall_r  * vf_rw(l)
      sunwall_r_shadewall  = sunwall_r  * vf_ww(l)
      sunwall_e            = em_wall  * sb * (t_sunwall **4)
      sunwall_e_sky        = sunwall_e  * vf_sw(l)
      sunwall_e_road       = sunwall_e  * vf_rw(l)
      sunwall_e_shadewall  = sunwall_e  * vf_ww(l)

      shadewall_a          = em_wall  * lwdown_shadewall
      shadewall_r          = (1._r8-em_wall ) * lwdown_shadewall
      shadewall_r_sky      = shadewall_r  * vf_sw(l)
      shadewall_r_road     = shadewall_r  * vf_rw(l)
      shadewall_r_sunwall  = shadewall_r  * vf_ww(l)
      shadewall_e          = em_wall  * sb * (t_shadewall **4)
      shadewall_e_sky      = shadewall_e  * vf_sw(l)
      shadewall_e_road     = shadewall_e  * vf_rw(l)
      shadewall_e_sunwall  = shadewall_e  * vf_ww(l)

      ! initialize sum of net and upward longwave radiation for road and both walls

      lwnet_improad    = improad_e    - improad_a
      lwnet_perroad    = perroad_e    - perroad_a
      lwnet_sunwall    = sunwall_e    - sunwall_a
      lwnet_shadewall  = shadewall_e  - shadewall_a

      lwup_improad    = improad_r_sky    + improad_e_sky
      lwup_perroad    = perroad_r_sky    + perroad_e_sky
      lwup_sunwall    = sunwall_r_sky    + sunwall_e_sky
      lwup_shadewall  = shadewall_r_sky  + shadewall_e_sky

      ! now account for absorption and reflection within canyon of fluxes from road and walls
      ! allowing for multiple reflections
      !
      ! (1) absorption and reflection. note: emission from road and walls absorbed by walls and roads
      !     only occurs in first iteration. zero out for later iterations.
      !
      !     road: fluxes from walls need to be projected to ground area
      !     wall: fluxes from road need to be projected to wall area
      !
      ! (2) add net longwave for ith reflection to total net longwave
      !
      ! (3) distribute reflected radiation to sky, road, and walls according to view factors
      !
      ! (4) add upward longwave radiation to sky from road and walls for ith reflection to total
      !
      ! (5) stop iteration when absorption for ith reflection is less than some nominal amount.
      !     small convergence criteria is required to ensure radiation is conserved

      do iter = 1, n
          ! step (1)

          lwtot  =  (sunwall_r_road  + sunwall_e_road   &
               + shadewall_r_road  + shadewall_e_road )*canyon_hwr
          road_a     = 0.0_r8
          road_r     = 0.0_r8
          improad_r  = (1._r8-em_improad ) * lwtot
          improad_a  =     em_improad   * lwtot
          road_a     = road_a  + improad_a *wtroad_imperv
          road_r     = road_r  + improad_r *wtroad_imperv
          perroad_r  = (1._r8-em_perroad ) * lwtot
          perroad_a  =     em_perroad   * lwtot
          road_a     = road_a  + perroad_a *wtroad_perv
          road_r     = road_r  + perroad_r *wtroad_perv

          lwtot  = (road_r_sunwall  + road_e_sunwall )/canyon_hwr  &
               + (shadewall_r_sunwall  + shadewall_e_sunwall )
          sunwall_a  =     em_wall   * lwtot
          sunwall_r  = (1._r8-em_wall ) * lwtot

          lwtot  = (road_r_shadewall  + road_e_shadewall )/canyon_hwr  &
               + (sunwall_r_shadewall  + sunwall_e_shadewall )
          shadewall_a  =     em_wall   * lwtot
          shadewall_r  = (1._r8-em_wall ) * lwtot

          sunwall_e_road       = 0._r8
          shadewall_e_road     = 0._r8
          road_e_sunwall       = 0._r8
          shadewall_e_sunwall  = 0._r8
          road_e_shadewall     = 0._r8
          sunwall_e_shadewall  = 0._r8

          ! step (2)

          lwnet_improad    = lwnet_improad    - improad_a
          lwnet_perroad    = lwnet_perroad    - perroad_a
          lwnet_sunwall    = lwnet_sunwall    - sunwall_a
          lwnet_shadewall  = lwnet_shadewall  - shadewall_a

          ! step (3)

          improad_r_sky        = improad_r  * vf_sr(l)
          improad_r_sunwall   = improad_r * vf_wr(l)
          improad_r_shadewall = improad_r * vf_wr(l)

          perroad_r_sky       = perroad_r * vf_sr(l)
          perroad_r_sunwall   = perroad_r * vf_wr(l)
          perroad_r_shadewall = perroad_r * vf_wr(l)

          road_r_sky          = road_r * vf_sr(l)
          road_r_sunwall      = road_r * vf_wr(l)
          road_r_shadewall    = road_r * vf_wr(l)

          sunwall_r_sky       = sunwall_r * vf_sw(l)
          sunwall_r_road      = sunwall_r * vf_rw(l)
          sunwall_r_shadewall = sunwall_r * vf_ww(l)

          shadewall_r_sky     = shadewall_r * vf_sw(l)
          shadewall_r_road    = shadewall_r * vf_rw(l)
          shadewall_r_sunwall = shadewall_r * vf_ww(l)

          ! step (4)

          lwup_improad   = lwup_improad   + improad_r_sky
          lwup_perroad   = lwup_perroad   + perroad_r_sky
          lwup_sunwall   = lwup_sunwall   + sunwall_r_sky
          lwup_shadewall = lwup_shadewall + shadewall_r_sky

          ! step (5)

          crit = max(road_a, sunwall_a, shadewall_a)
          if (crit < .001_r8) exit
        end do

        if (iter >= n) then
#ifndef _OPENACC
           write (iulog,*) 'urban net longwave radiation error: no convergence'
           write (iulog,*) 'elm model is stopping'
           call endrun(decomp_index=l, elmlevel=namel, msg=errmsg(__FILE__, __LINE__))
#endif
        endif

         ! total net longwave radiation for canyon. project wall fluxes to horizontal surface

         lwnet_canyon = 0.0_r8
         lwnet_canyon = lwnet_canyon + lwnet_improad*wtroad_imperv
         lwnet_canyon = lwnet_canyon + lwnet_perroad*wtroad_perv
         lwnet_canyon = lwnet_canyon + (lwnet_sunwall + lwnet_shadewall)*canyon_hwr

         ! total emitted longwave for canyon. project wall fluxes to horizontal

         lwup_canyon = 0.0_r8
         lwup_canyon = lwup_canyon + lwup_improad*wtroad_imperv
         lwup_canyon = lwup_canyon + lwup_perroad*wtroad_perv
         lwup_canyon = lwup_canyon + (lwup_sunwall + lwup_shadewall)*canyon_hwr

         ! conservation check. note: previous conservation check confirms partioning of incident
         ! atmospheric longwave radiation to road and walls is conserved as
         ! lwdown (from atmosphere) = lwdown_improad + lwdown_perroad + (lwdown_sunwall + lwdown_shadewall)*canyon_hwr

         err = lwnet_canyon - (lwup_canyon - lwdown)
         if (abs(err) > .10_r8 ) then
#ifndef _OPENACC
            write (iulog,*) 'urban net longwave radiation balance error',err
            write (iulog,*) 'elm model is stopping'
            call endrun(decomp_index=l, elmlevel=namel, msg=errmsg(__FILE__, __LINE__))
#endif
         end if

      ! Net longwave radiation for roof

      lwup_roof = em_roof*sb*(t_roof**4) + (1._r8-em_roof)*lwdown
      lwnet_roof = lwup_roof - lwdown

    end associate

  end subroutine net_longwave

end module UrbanRadiationMod
