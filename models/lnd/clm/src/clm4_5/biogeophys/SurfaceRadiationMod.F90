module SurfaceRadiationMod

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate solar fluxes absorbed by vegetation and ground surface
  !
  ! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use shr_log_mod , only: errMsg => shr_log_errMsg
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: SurfaceRadiation ! Solar fluxes absorbed by veg and ground surface
  !------------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------------
   subroutine SurfaceRadiation(bounds, num_nourbanp, filter_nourbanp)
     !
     ! !DESCRIPTION: 
     ! Solar fluxes absorbed by vegetation and ground surface
     ! Note possible problem when land is on different grid than atmosphere.
     ! Land may have sun above the horizon (coszen > 0) but atmosphere may
     ! have sun below the horizon (forc_solad = 0 and forc_solai = 0). This is okay
     ! because all fluxes (absorbed, reflected, transmitted) are multiplied
     ! by the incoming flux and all will equal zero.
     ! Atmosphere may have sun above horizon (forc_solad > 0 and forc_solai > 0) but
     ! land may have sun below horizon. This is okay because fabd, fabi,
     ! ftdd, ftid, and ftii all equal zero so that sabv=sabg=fsa=0. Also,
     ! albd and albi equal one so that fsr=forc_solad+forc_solai. In other words, all
     ! the radiation is reflected. NDVI should equal zero in this case.
     ! However, the way the code is currently implemented this is only true
     ! if (forc_solad+forc_solai)|vis = (forc_solad+forc_solai)|nir.
     ! Output variables are parsun,parsha,sabv,sabg,fsa,fsr,ndvi
     !
     ! !USES:
     use clmtype
     use clm_atmlnd      , only : clm_a2l
     use clm_varpar      , only : numrad, nlevsno
     use clm_varcon      , only : spval, istsoil, degpsec, isecspday, istcrop
     use clm_varctl      , only : subgridflag, use_snicar_frc, iulog
     use clm_time_manager, only : get_curr_date, get_step_size
     use SNICARMod       , only : DO_SNO_OC
     use abortutils      , only : endrun
     use decompMod       , only : bounds_type
     !
     ! !ARGUMENTS:
     implicit none
     type(bounds_type), intent(in) :: bounds   ! bounds
     integer, intent(in) :: num_nourbanp       ! number of pfts in non-urban points in pft filter
     integer, intent(in) :: filter_nourbanp(:) ! pft filter for non-urban points
     !
     ! !LOCAL VARIABLES:
     integer , parameter :: nband = numrad           ! number of solar radiation waveband classes
     real(r8), parameter :: mpe = 1.e-06_r8          ! prevents overflow for division by zero
     integer  :: fp                                  ! non-urban filter pft index
     integer  :: p                                   ! pft index
     integer  :: c                                   ! column index
     integer  :: l                                   ! landunit index
     integer  :: g                                   ! grid cell index
     integer  :: ib                                  ! waveband number (1=vis, 2=nir)
     integer  :: iv                                  ! canopy layer
     real(r8) :: absrad                              ! absorbed solar radiation (W/m**2)
     real(r8) :: rnir                                ! reflected solar radiation [nir] (W/m**2)
     real(r8) :: rvis                                ! reflected solar radiation [vis] (W/m**2)
     real(r8) :: trd(bounds%begp:bounds%endp,numrad) ! transmitted solar radiation: direct (W/m**2)
     real(r8) :: tri(bounds%begp:bounds%endp,numrad) ! transmitted solar radiation: diffuse (W/m**2)
     real(r8) :: cad(bounds%begp:bounds%endp,numrad) ! direct beam absorbed by canopy (W/m**2)
     real(r8) :: cai(bounds%begp:bounds%endp,numrad) ! diffuse radiation absorbed by canopy (W/m**2)
     integer  :: local_secp1                         ! seconds into current date in local time
     real(r8) :: dtime                               ! land model time step (sec)
     integer  :: year,month,day,secs                 ! calendar info for current time step
     integer  :: i                                   ! layer index [idx]
     real(r8) :: sabg_snl_sum                        ! temporary, absorbed energy in all active snow layers [W/m2]
     real(r8) :: absrad_pur                          ! temp: absorbed solar radiation by pure snow [W/m2]
     real(r8) :: absrad_bc                           ! temp: absorbed solar radiation without BC [W/m2]
     real(r8) :: absrad_oc                           ! temp: absorbed solar radiation without OC [W/m2]
     real(r8) :: absrad_dst                          ! temp: absorbed solar radiation without dust [W/m2]
     real(r8) :: sabg_pur(bounds%begp:bounds%endp)   ! solar radiation absorbed by ground with pure snow [W/m2]
     real(r8) :: sabg_bc(bounds%begp:bounds%endp)    ! solar radiation absorbed by ground without BC [W/m2]
     real(r8) :: sabg_oc(bounds%begp:bounds%endp)    ! solar radiation absorbed by ground without OC [W/m2]
     real(r8) :: sabg_dst(bounds%begp:bounds%endp)   ! solar radiation absorbed by ground without dust [W/m2]
     real(r8) :: parveg(bounds%begp:bounds%endp)     ! absorbed par by vegetation (W/m**2)
     integer  :: nstep                               ! time step number
     !------------------------------------------------------------------------------

   associate(& 
   forc_solad                =>    clm_a2l%forc_solad      , & ! Input:  [real(r8) (:,:)]  direct beam radiation (W/m**2)        
   forc_solai                =>    clm_a2l%forc_solai      , & ! Input:  [real(r8) (:,:)]  diffuse radiation (W/m**2)            
   albgrd                    =>    cps%albgrd              , & ! Input:  [real(r8) (:,:)]  ground albedo (direct)                
   albgri                    =>    cps%albgri              , & ! Input:  [real(r8) (:,:)]  ground albedo (diffuse)               
   coszen                    =>    cps%coszen              , & ! Input:  [real(r8) (:)]  cosine of solar zenith angle            
   albsod                    =>    cps%albsod              , & ! Input:  [real(r8) (:,:)]  direct-beam soil albedo (col,bnd) [frc]
   albsoi                    =>    cps%albsoi              , & ! Input:  [real(r8) (:,:)]  diffuse soil albedo (col,bnd) [frc]   
   sabg_soil                 =>    pef%sabg_soil           , & ! Input:  [real(r8) (:)]  solar radiation absorbed by soil (W/m**2)
   sabg_snow                 =>    pef%sabg_snow           , & ! Input:  [real(r8) (:)]  solar radiation absorbed by snow (W/m**2)
   elai                      =>    pps%elai                , & ! Input:  [real(r8) (:)]  one-sided leaf area index with burying by snow
   esai                      =>    pps%esai                , & ! Input:  [real(r8) (:)]  one-sided stem area index with burying by snow
   albd                      =>    pps%albd                , & ! Input:  [real(r8) (:,:)]  surface albedo (direct)               
   albi                      =>    pps%albi                , & ! Input:  [real(r8) (:,:)]  surface albedo (diffuse)              
   fabd                      =>    pps%fabd                , & ! Input:  [real(r8) (:,:)]  flux absorbed by canopy per unit direct flux
   fabd_sun                  =>    pps%fabd_sun            , & ! Input:  [real(r8) (:,:)]  flux absorbed by sunlit canopy per unit direct flux
   fabd_sha                  =>    pps%fabd_sha            , & ! Input:  [real(r8) (:,:)]  flux absorbed by shaded canopy per unit direct flux
   fabi                      =>    pps%fabi                , & ! Input:  [real(r8) (:,:)]  flux absorbed by canopy per unit diffuse flux
   fabi_sun                  =>    pps%fabi_sun            , & ! Input:  [real(r8) (:,:)]  flux absorbed by sunlit canopy per unit diffuse flux
   fabi_sha                  =>    pps%fabi_sha            , & ! Input:  [real(r8) (:,:)]  flux absorbed by shaded canopy per unit diffuse flux
   ftdd                      =>    pps%ftdd                , & ! Input:  [real(r8) (:,:)]  down direct flux below canopy per unit direct flux
   ftid                      =>    pps%ftid                , & ! Input:  [real(r8) (:,:)]  down diffuse flux below canopy per unit direct flux
   ftii                      =>    pps%ftii                , & ! Input:  [real(r8) (:,:)]  down diffuse flux below canopy per unit diffuse flux
   nrad                      =>    pps%nrad                , & ! Input:  [integer (:)]  number of canopy layers, above snow for radiative transfer
   fabd_sun_z                =>    pps%fabd_sun_z          , & ! Input:  [real(r8) (:,:)]  absorbed sunlit leaf direct  PAR (per unit lai+sai) for each canopy layer
   fabd_sha_z                =>    pps%fabd_sha_z          , & ! Input:  [real(r8) (:,:)]  absorbed shaded leaf direct  PAR (per unit lai+sai) for each canopy layer
   fabi_sun_z                =>    pps%fabi_sun_z          , & ! Input:  [real(r8) (:,:)]  absorbed sunlit leaf diffuse PAR (per unit lai+sai) for each canopy layer
   fabi_sha_z                =>    pps%fabi_sha_z          , & ! Input:  [real(r8) (:,:)]  absorbed shaded leaf diffuse PAR (per unit lai+sai) for each canopy layer
   fsun_z                    =>    pps%fsun_z              , & ! Input:  [real(r8) (:,:)]  sunlit fraction of canopy layer       
   tlai_z                    =>    pps%tlai_z              , & ! Input:  [real(r8) (:,:)]  tlai increment for canopy layer       
   tsai_z                    =>    pps%tsai_z              , & ! Input:  [real(r8) (:,:)]  tsai increment for canopy layer       
   laisun                    =>    pps%laisun              , & ! Output: [real(r8) (:)]  sunlit leaf area                        
   laisha                    =>    pps%laisha              , & ! Output: [real(r8) (:)]  shaded leaf area                        
   laisun_z                  =>    pps%laisun_z            , & ! Output: [real(r8) (:,:)]  sunlit leaf area for canopy layer     
   laisha_z                  =>    pps%laisha_z            , & ! Output: [real(r8) (:,:)]  shaded leaf area for canopy layer     
   fsun                      =>    pps%fsun                , & ! Output: [real(r8) (:)]  sunlit fraction of canopy               
   sabg                      =>    pef%sabg                , & ! Output: [real(r8) (:)]  solar radiation absorbed by ground (W/m**2)
   sabv                      =>    pef%sabv                , & ! Output: [real(r8) (:)]  solar radiation absorbed by vegetation (W/m**2)
   snow_depth                =>    cps%snow_depth          , & ! Output: [real(r8) (:)]  snow height (m)                         
   fsa                       =>    pef%fsa                 , & ! Output: [real(r8) (:)]  solar radiation absorbed (total) (W/m**2)
   fsa_r                     =>    pef%fsa_r               , & ! Output: [real(r8) (:)]  rural solar radiation absorbed (total) (W/m**2)
   fsr                       =>    pef%fsr                 , & ! Output: [real(r8) (:)]  solar radiation reflected (W/m**2)      
   parsun_z                  =>    pef%parsun_z            , & ! Output: [real(r8) (:,:)]  absorbed PAR for sunlit leaves in canopy layer
   parsha_z                  =>    pef%parsha_z            , & ! Output: [real(r8) (:,:)]  absorbed PAR for shaded leaves in canopy layer
   fsds_vis_d                =>    pef%fsds_vis_d          , & ! Output: [real(r8) (:)]  incident direct beam vis solar radiation (W/m**2)
   fsds_nir_d                =>    pef%fsds_nir_d          , & ! Output: [real(r8) (:)]  incident direct beam nir solar radiation (W/m**2)
   fsds_vis_i                =>    pef%fsds_vis_i          , & ! Output: [real(r8) (:)]  incident diffuse vis solar radiation (W/m**2)
   fsds_nir_i                =>    pef%fsds_nir_i          , & ! Output: [real(r8) (:)]  incident diffuse nir solar radiation (W/m**2)
   fsr_vis_d                 =>    pef%fsr_vis_d           , & ! Output: [real(r8) (:)]  reflected direct beam vis solar radiation (W/m**2)
   fsr_nir_d                 =>    pef%fsr_nir_d           , & ! Output: [real(r8) (:)]  reflected direct beam nir solar radiation (W/m**2)
   fsr_vis_i                 =>    pef%fsr_vis_i           , & ! Output: [real(r8) (:)]  reflected diffuse vis solar radiation (W/m**2)
   fsr_nir_i                 =>    pef%fsr_nir_i           , & ! Output: [real(r8) (:)]  reflected diffuse nir solar radiation (W/m**2)
   fsds_vis_d_ln             =>    pef%fsds_vis_d_ln       , & ! Output: [real(r8) (:)]  incident direct beam vis solar rad at local noon (W/m**2)
   fsds_nir_d_ln             =>    pef%fsds_nir_d_ln       , & ! Output: [real(r8) (:)]  incident direct beam nir solar rad at local noon (W/m**2)
   parveg_ln                 =>    pef%parveg_ln           , & ! Output: [real(r8) (:)]  absorbed par by vegetation at local noon (W/m**2)
   fsds_vis_i_ln             =>    pef%fsds_vis_i_ln       , & ! Output: [real(r8) (:)]  incident diffuse beam vis solar rad at local noon (W/m**2)
   fsr_vis_d_ln              =>    pef%fsr_vis_d_ln        , & ! Output: [real(r8) (:)]  reflected direct beam vis solar rad at local noon (W/m**2)
   fsr_nir_d_ln              =>    pef%fsr_nir_d_ln        , & ! Output: [real(r8) (:)]  reflected direct beam nir solar rad at local noon (W/m**2)
   frac_sno                  =>    cps%frac_sno            , & ! Output: [real(r8) (:)]  fraction of ground covered by snow (0 to 1)
   flx_absdv                 =>    cps%flx_absdv           , & ! Output: [real(r8) (:,:)]  direct flux absorption factor (col,lyr): VIS [frc]
   flx_absdn                 =>    cps%flx_absdn           , & ! Output: [real(r8) (:,:)]  direct flux absorption factor (col,lyr): NIR [frc]
   flx_absiv                 =>    cps%flx_absiv           , & ! Output: [real(r8) (:,:)]  diffuse flux absorption factor (col,lyr): VIS [frc]
   flx_absin                 =>    cps%flx_absin           , & ! Output: [real(r8) (:,:)]  diffuse flux absorption factor (col,lyr): NIR [frc]
   sabg_lyr                  =>    pef%sabg_lyr            , & ! Output: [real(r8) (:,:)]  absorbed radiative flux (pft,lyr) [W/m2]
   sabg_pen                  =>    pef%sabg_pen            , & ! Output: [real(r8) (:)]  (rural) shortwave radiation penetrating top soisno layer [W/m2]
   snl                       =>    cps%snl                 , & ! Output: [integer (:)]  negative number of snow layers [nbr]     
   sfc_frc_aer               =>    pef%sfc_frc_aer         , & ! Output: [real(r8) (:)]  surface forcing of snow with all aerosols (pft) [W/m2]
   sfc_frc_aer_sno           =>    pef%sfc_frc_aer_sno     , & ! Output: [real(r8) (:)]  surface forcing of snow with all aerosols, averaged only when snow is present (pft) [W/m2]
   albgrd_pur                =>    cps%albgrd_pur          , & ! Output: [real(r8) (:,:)]  pure snow ground albedo (direct)      
   albgri_pur                =>    cps%albgri_pur          , & ! Output: [real(r8) (:,:)]  pure snow ground albedo (diffuse)     
   sfc_frc_bc                =>    pef%sfc_frc_bc          , & ! Output: [real(r8) (:)]  surface forcing of snow with BC (pft) [W/m2]
   sfc_frc_bc_sno            =>    pef%sfc_frc_bc_sno      , & ! Output: [real(r8) (:)]  surface forcing of snow with BC, averaged only when snow is present (pft) [W/m2]
   albgrd_bc                 =>    cps%albgrd_bc           , & ! Output: [real(r8) (:,:)]  ground albedo without BC (direct) (col,bnd)
   albgri_bc                 =>    cps%albgri_bc           , & ! Output: [real(r8) (:,:)]  ground albedo without BC (diffuse) (col,bnd)
   sfc_frc_oc                =>    pef%sfc_frc_oc          , & ! Output: [real(r8) (:)]  surface forcing of snow with OC (pft) [W/m2]
   sfc_frc_oc_sno            =>    pef%sfc_frc_oc_sno      , & ! Output: [real(r8) (:)]  surface forcing of snow with OC, averaged only when snow is present (pft) [W/m2]
   albgrd_oc                 =>    cps%albgrd_oc           , & ! Output: [real(r8) (:,:)]  ground albedo without OC (direct) (col,bnd)
   albgri_oc                 =>    cps%albgri_oc           , & ! Output: [real(r8) (:,:)]  ground albedo without OC (diffuse) (col,bnd)
   sfc_frc_dst               =>    pef%sfc_frc_dst         , & ! Output: [real(r8) (:)]  surface forcing of snow with dust (pft) [W/m2]
   sfc_frc_dst_sno           =>    pef%sfc_frc_dst_sno     , & ! Output: [real(r8) (:)]  surface forcing of snow with dust, averaged only when snow is present (pft) [W/m2]
   albgrd_dst                =>    cps%albgrd_dst          , & ! Output: [real(r8) (:,:)]  ground albedo without dust (direct) (col,bnd)
   albgri_dst                =>    cps%albgri_dst          , & ! Output: [real(r8) (:,:)]  ground albedo without dust (diffuse) (col,bnd)
   albsnd_hst                =>    cps%albsnd_hst          , & ! Output: [real(r8) (:,:)]  snow albedo, direct, for history files (col,bnd) [frc]
   albsni_hst                =>    cps%albsni_hst          , & ! Output: [real(r8) (:,:)]  snow ground albedo, diffuse, for history files (col,bnd
   fsr_sno_vd                =>    pef%fsr_sno_vd          , & ! Output: [real(r8) (:)]  reflected visible, direct radiation from snow (for history files) (pft) [W/m2]
   fsr_sno_nd                =>    pef%fsr_sno_nd          , & ! Output: [real(r8) (:)]  reflected near-IR, direct radiation from snow (for history files) (pft) [W/m2]
   fsr_sno_vi                =>    pef%fsr_sno_vi          , & ! Output: [real(r8) (:)]  reflected visible, diffuse radiation from snow (for history files) (pft) [W/m2]
   fsr_sno_ni                =>    pef%fsr_sno_ni          , & ! Output: [real(r8) (:)]  reflected near-IR, diffuse radiation from snow (for history files) (pft) [W/m2]
   fsds_sno_vd               =>    pef%fsds_sno_vd         , & ! Output: [real(r8) (:)]  incident visible, direct radiation on snow (for history files) (pft) [W/m2]
   fsds_sno_nd               =>    pef%fsds_sno_nd         , & ! Output: [real(r8) (:)]  incident near-IR, direct radiation on snow (for history files) (pft) [W/m2]
   fsds_sno_vi               =>    pef%fsds_sno_vi         , & ! Output: [real(r8) (:)]  incident visible, diffuse radiation on snow (for history files) (pft) [W/m2]
   fsds_sno_ni               =>    pef%fsds_sno_ni         , & ! Output: [real(r8) (:)]  incident near-IR, diffuse radiation on snow (for history files) (pft) [W/m2]
   sub_surf_abs_SW           =>    cps%sub_surf_abs_SW       & ! Output: [real(r8) (:)]
   )

     ! Determine seconds off current time step
     
     dtime = get_step_size()
     call get_curr_date (year, month, day, secs)

     ! Initialize fluxes

     do fp = 1,num_nourbanp
        p = filter_nourbanp(fp)
        sabg_soil(p)  = 0._r8
        sabg_snow(p)  = 0._r8
        sabg(p)       = 0._r8
        sabv(p)       = 0._r8
        fsa(p)        = 0._r8
        l = pft%landunit(p)
        if (lun%itype(l)==istsoil .or. lun%itype(l)==istcrop) then
           fsa_r(p)      = 0._r8
        end if
        sabg_lyr(p,:) = 0._r8
        sabg_pur(p)   = 0._r8
        sabg_bc(p)    = 0._r8
        sabg_oc(p)    = 0._r8
        sabg_dst(p)   = 0._r8
        do iv = 1, nrad(p)
           parsun_z(p,iv) = 0._r8
           parsha_z(p,iv) = 0._r8
           laisun_z(p,iv) = 0._r8
           laisha_z(p,iv) = 0._r8
        end do
     end do 

     ! Loop over pfts to calculate laisun_z and laisha_z for each layer.
     ! Derive canopy laisun, laisha, and fsun from layer sums.
     ! If sun/shade big leaf code, nrad=1 and fsun_z(p,1) and tlai_z(p,1) from
     ! SurfaceAlbedo is canopy integrated so that layer value equals canopy value.

     do fp = 1,num_nourbanp
        p = filter_nourbanp(fp)
        laisun(p) = 0._r8
        laisha(p) = 0._r8
        do iv = 1, nrad(p)
           laisun_z(p,iv) = tlai_z(p,iv) * fsun_z(p,iv)
           laisha_z(p,iv) = tlai_z(p,iv) * (1._r8 - fsun_z(p,iv))
           laisun(p) = laisun(p) + laisun_z(p,iv) 
           laisha(p) = laisha(p) + laisha_z(p,iv) 
        end do
        if (elai(p) > 0._r8) then
           fsun(p) = laisun(p) / elai(p)
        else
           fsun(p) = 0._r8
        end if
     end do
        
     ! Loop over nband wavebands
     do ib = 1, nband
        do fp = 1,num_nourbanp
           p = filter_nourbanp(fp)
           c = pft%column(p)
           l = pft%landunit(p)
           g = pft%gridcell(p)

           ! Absorbed by canopy

           cad(p,ib) = forc_solad(g,ib)*fabd(p,ib)
           cai(p,ib) = forc_solai(g,ib)*fabi(p,ib)
           sabv(p) = sabv(p) + cad(p,ib) + cai(p,ib)
           fsa(p)  = fsa(p)  + cad(p,ib) + cai(p,ib)
           if (ib == 1) then
              parveg(p) = cad(p,ib) + cai(p,ib)
           end if
           if (lun%itype(l)==istsoil .or. lun%itype(l)==istcrop) then
              fsa_r(p)  = fsa_r(p)  + cad(p,ib) + cai(p,ib)
           end if

           ! Absorbed PAR profile through canopy
           ! If sun/shade big leaf code, nrad=1 and fluxes from SurfaceAlbedo
           ! are canopy integrated so that layer values equal big leaf values.

           if (ib == 1) then
              do iv = 1, nrad(p)
                 parsun_z(p,iv) = forc_solad(g,ib)*fabd_sun_z(p,iv) + forc_solai(g,ib)*fabi_sun_z(p,iv)
                 parsha_z(p,iv) = forc_solad(g,ib)*fabd_sha_z(p,iv) + forc_solai(g,ib)*fabi_sha_z(p,iv)
              end do
           end if

           ! Transmitted = solar fluxes incident on ground

           trd(p,ib) = forc_solad(g,ib)*ftdd(p,ib)
           tri(p,ib) = forc_solad(g,ib)*ftid(p,ib) + forc_solai(g,ib)*ftii(p,ib)

           ! Solar radiation absorbed by ground surface

           ! calculate absorbed solar by soil/snow separately
           absrad  = trd(p,ib)*(1._r8-albsod(c,ib)) + tri(p,ib)*(1._r8-albsoi(c,ib))
           sabg_soil(p) = sabg_soil(p) + absrad
           absrad  = trd(p,ib)*(1._r8-albsnd_hst(c,ib)) + tri(p,ib)*(1._r8-albsni_hst(c,ib))
           sabg_snow(p) = sabg_snow(p) + absrad
           absrad  = trd(p,ib)*(1._r8-albgrd(c,ib)) + tri(p,ib)*(1._r8-albgri(c,ib))
           sabg(p) = sabg(p) + absrad
           fsa(p)  = fsa(p)  + absrad
           if (lun%itype(l)==istsoil .or. lun%itype(l)==istcrop) then
              fsa_r(p)  = fsa_r(p)  + absrad
           end if
           if (snl(c) == 0) then
              sabg_snow(p) = sabg(p)
              sabg_soil(p) = sabg(p)
           endif
           ! if no subgrid fluxes, make sure to set both components equal to weighted average
           if (subgridflag == 0) then 
              sabg_snow(p) = sabg(p)
              sabg_soil(p) = sabg(p)
           endif

           if (use_snicar_frc) then
              ! Solar radiation absorbed by ground surface without BC
              absrad_bc = trd(p,ib)*(1._r8-albgrd_bc(c,ib)) + tri(p,ib)*(1._r8-albgri_bc(c,ib))
              sabg_bc(p) = sabg_bc(p) + absrad_bc

              ! Solar radiation absorbed by ground surface without OC
              absrad_oc = trd(p,ib)*(1._r8-albgrd_oc(c,ib)) + tri(p,ib)*(1._r8-albgri_oc(c,ib))
              sabg_oc(p) = sabg_oc(p) + absrad_oc

              ! Solar radiation absorbed by ground surface without dust
              absrad_dst = trd(p,ib)*(1._r8-albgrd_dst(c,ib)) + tri(p,ib)*(1._r8-albgri_dst(c,ib))
              sabg_dst(p) = sabg_dst(p) + absrad_dst

              ! Solar radiation absorbed by ground surface without any aerosols
              absrad_pur = trd(p,ib)*(1._r8-albgrd_pur(c,ib)) + tri(p,ib)*(1._r8-albgri_pur(c,ib))
              sabg_pur(p) = sabg_pur(p) + absrad_pur
           end if

        end do ! end of pft loop
     end do ! end nbands loop   

     !   compute absorbed flux in each snow layer and top soil layer,
     !   based on flux factors computed in the radiative transfer portion of SNICAR.

     do fp = 1,num_nourbanp
        p = filter_nourbanp(fp)
        c = pft%column(p)
        l = pft%landunit(p)
        sabg_snl_sum = 0._r8

        sub_surf_abs_SW(c) = 0._r8

        ! CASE1: No snow layers: all energy is absorbed in top soil layer
        if (snl(c) == 0) then
           sabg_lyr(p,:) = 0._r8
           sabg_lyr(p,1) = sabg(p)
           sabg_snl_sum  = sabg_lyr(p,1)

           ! CASE 2: Snow layers present: absorbed radiation is scaled according to 
           ! flux factors computed by SNICAR
        else
           do i = -nlevsno+1,1,1
              sabg_lyr(p,i) = flx_absdv(c,i)*trd(p,1) + flx_absdn(c,i)*trd(p,2) + &
                   flx_absiv(c,i)*tri(p,1) + flx_absin(c,i)*tri(p,2)
              ! summed radiation in active snow layers:
              if (i >= snl(c)+1) then
                 sabg_snl_sum = sabg_snl_sum + sabg_lyr(p,i)
              endif
              if (i > snl(c)+1) then ! if snow layer is below surface snow layer
                 !accumulate subsurface flux as a diagnostic for history file
                 sub_surf_abs_SW(c) = sub_surf_abs_SW(c) + sabg_lyr(p,i)
              endif
           enddo

           ! Divide absorbed by total, to get % absorbed in subsurface
           if (sabg_snl_sum /= 0._r8) then
              sub_surf_abs_SW(c) = sub_surf_abs_SW(c)/sabg_snl_sum
           else
              sub_surf_abs_SW(c) = 0._r8
           endif

           ! Error handling: The situation below can occur when solar radiation is 
           ! NOT computed every timestep.
           ! When the number of snow layers has changed in between computations of the 
           ! absorbed solar energy in each layer, we must redistribute the absorbed energy
           ! to avoid physically unrealistic conditions. The assumptions made below are 
           ! somewhat arbitrary, but this situation does not arise very frequently. 
           ! This error handling is implemented to accomodate any value of the
           ! radiation frequency.
           ! change condition to match sabg_snow isntead of sabg
           if (abs(sabg_snl_sum-sabg_snow(p)) > 0.00001_r8) then
              if (snl(c) == 0) then
                 sabg_lyr(p,-4:0) = 0._r8
                 sabg_lyr(p,1) = sabg(p)
              elseif (snl(c) == -1) then
                 sabg_lyr(p,-4:-1) = 0._r8
                 sabg_lyr(p,0) = sabg_snow(p)*0.6_r8
                 sabg_lyr(p,1) = sabg_snow(p)*0.4_r8
              else
                 sabg_lyr(p,:) = 0._r8
                 sabg_lyr(p,snl(c)+1) = sabg_snow(p)*0.75_r8
                 sabg_lyr(p,snl(c)+2) = sabg_snow(p)*0.25_r8
              endif
           endif

           ! If shallow snow depth, all solar radiation absorbed in top or top two snow layers
           ! to prevent unrealistic timestep soil warming 
           if (subgridflag == 0) then 
              if (snow_depth(c) < 0.10_r8) then
                 if (snl(c) == 0) then
                    sabg_lyr(p,-4:0) = 0._r8
                    sabg_lyr(p,1) = sabg(p)
                 elseif (snl(c) == -1) then
                    sabg_lyr(p,-4:-1) = 0._r8
                    sabg_lyr(p,0) = sabg(p)
                    sabg_lyr(p,1) = 0._r8
                 else
                    sabg_lyr(p,:) = 0._r8
                    sabg_lyr(p,snl(c)+1) = sabg(p)*0.75_r8
                    sabg_lyr(p,snl(c)+2) = sabg(p)*0.25_r8
                 endif
              endif
           endif
        endif

        ! This situation should not happen:
        if (abs(sum(sabg_lyr(p,:))-sabg_snow(p)) > 0.00001_r8) then
           write(iulog,*)"SNICAR ERROR: Absorbed ground radiation not equal to summed snow layer radiation"
           write(iulog,*)"Diff        = ",sum(sabg_lyr(p,:))-sabg_snow(p)
           write(iulog,*)"sabg_snow(p)= ",sabg_snow(p)
           write(iulog,*)"sabg_sum(p) = ",sum(sabg_lyr(p,:))
           write(iulog,*)"snl(c)      = ",snl(c)
           write(iulog,*)"flx_absdv1  = ",trd(p,1)*(1.-albgrd(c,1))
           write(iulog,*)"flx_absdv2  = ",sum(flx_absdv(c,:))*trd(p,1)
           write(iulog,*)"flx_absiv1  = ",tri(p,1)*(1.-albgri(c,1))
           write(iulog,*)"flx_absiv2  = ",sum(flx_absiv(c,:))*tri(p,1)
           write(iulog,*)"flx_absdn1  = ",trd(p,2)*(1.-albgrd(c,2))
           write(iulog,*)"flx_absdn2  = ",sum(flx_absdn(c,:))*trd(p,2)
           write(iulog,*)"flx_absin1  = ",tri(p,2)*(1.-albgri(c,2))
           write(iulog,*)"flx_absin2  = ",sum(flx_absin(c,:))*tri(p,2)
           write(iulog,*)"albgrd_nir  = ",albgrd(c,2)
           write(iulog,*)"coszen      = ",coszen(c)
           call endrun(decomp_index=c, clmlevel=namec, msg=errmsg(__FILE__, __LINE__))
        endif

        ! Diagnostic: shortwave penetrating ground (e.g. top layer)
        if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
           sabg_pen(p) = sabg(p) - sabg_lyr(p, snl(c)+1)
        end if

        if (use_snicar_frc) then
           
           ! BC aerosol forcing (pft-level):
           sfc_frc_bc(p) = sabg(p) - sabg_bc(p)

           ! OC aerosol forcing (pft-level):
           if (DO_SNO_OC) then
              sfc_frc_oc(p) = sabg(p) - sabg_oc(p)
           else
              sfc_frc_oc(p) = 0._r8
           endif

           ! dust aerosol forcing (pft-level):
           sfc_frc_dst(p) = sabg(p) - sabg_dst(p)

           ! all-aerosol forcing (pft-level):
           sfc_frc_aer(p) = sabg(p) - sabg_pur(p)        

           ! forcings averaged only over snow:
           if (frac_sno(c) > 0._r8) then
              sfc_frc_bc_sno(p)  = sfc_frc_bc(p)/frac_sno(c)
              sfc_frc_oc_sno(p)  = sfc_frc_oc(p)/frac_sno(c)
              sfc_frc_dst_sno(p) = sfc_frc_dst(p)/frac_sno(c)
              sfc_frc_aer_sno(p) = sfc_frc_aer(p)/frac_sno(c)
           else
              sfc_frc_bc_sno(p)  = spval
              sfc_frc_oc_sno(p)  = spval
              sfc_frc_dst_sno(p) = spval
              sfc_frc_aer_sno(p) = spval
           endif
        end if
     enddo

     ! Radiation diagnostics

     do fp = 1,num_nourbanp
        p = filter_nourbanp(fp)
        g = pft%gridcell(p)

        ! NDVI and reflected solar radiation

        rvis = albd(p,1)*forc_solad(g,1) + albi(p,1)*forc_solai(g,1)
        rnir = albd(p,2)*forc_solad(g,2) + albi(p,2)*forc_solai(g,2)
        fsr(p) = rvis + rnir

        fsds_vis_d(p) = forc_solad(g,1)
        fsds_nir_d(p) = forc_solad(g,2)
        fsds_vis_i(p) = forc_solai(g,1)
        fsds_nir_i(p) = forc_solai(g,2)
        fsr_vis_d(p)  = albd(p,1)*forc_solad(g,1)
        fsr_nir_d(p)  = albd(p,2)*forc_solad(g,2)
        fsr_vis_i(p)  = albi(p,1)*forc_solai(g,1)
        fsr_nir_i(p)  = albi(p,2)*forc_solai(g,2)

        local_secp1 = secs + nint((grc%londeg(g)/degpsec)/dtime)*dtime
        local_secp1 = mod(local_secp1,isecspday)
        if (local_secp1 == isecspday/2) then
           fsds_vis_d_ln(p) = forc_solad(g,1)
           fsds_nir_d_ln(p) = forc_solad(g,2)
           fsr_vis_d_ln(p) = albd(p,1)*forc_solad(g,1)
           fsr_nir_d_ln(p) = albd(p,2)*forc_solad(g,2)
           fsds_vis_i_ln(p) = forc_solai(g,1)
           parveg_ln(p)     = parveg(p)
        else
           fsds_vis_d_ln(p) = spval
           fsds_nir_d_ln(p) = spval
           fsr_vis_d_ln(p) = spval
           fsr_nir_d_ln(p) = spval
           fsds_vis_i_ln(p) = spval
           parveg_ln(p)     = spval
        end if

        ! diagnostic variables (downwelling and absorbed radiation partitioning) for history files
        ! (OPTIONAL)
        c = pft%column(p)
        if (snl(c) < 0) then
           fsds_sno_vd(p) = forc_solad(g,1)
           fsds_sno_nd(p) = forc_solad(g,2)
           fsds_sno_vi(p) = forc_solai(g,1)
           fsds_sno_ni(p) = forc_solai(g,2)

           fsr_sno_vd(p) = fsds_vis_d(p)*albsnd_hst(c,1)
           fsr_sno_nd(p) = fsds_nir_d(p)*albsnd_hst(c,2)
           fsr_sno_vi(p) = fsds_vis_i(p)*albsni_hst(c,1)
           fsr_sno_ni(p) = fsds_nir_i(p)*albsni_hst(c,2)
        else
           fsds_sno_vd(p) = spval
           fsds_sno_nd(p) = spval
           fsds_sno_vi(p) = spval
           fsds_sno_ni(p) = spval

           fsr_sno_vd(p) = spval
           fsr_sno_nd(p) = spval
           fsr_sno_vi(p) = spval
           fsr_sno_ni(p) = spval
        endif
     end do 

    end associate 
    end subroutine SurfaceRadiation

end module SurfaceRadiationMod
