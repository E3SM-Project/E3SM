module UrbanAlbedoMod

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
  use clm_varctl        , only : iulog
  use abortutils        , only : endrun  
  use UrbanParamsType   , only : urbanparams_type
  use WaterstateType    , only : waterstate_type
  use SolarAbsorbedType , only : solarabs_type
  use SurfaceAlbedoType , only : surfalb_type
  use LandunitType      , only : lun_pp                
  use ColumnType        , only : col_pp
  use ColumnDataType    , only : col_ws  
  use VegetationType    , only : veg_pp                
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: UrbanAlbedo       ! Urban physics - albedos  
  !
  ! PRIVATE MEMBER FUNCTIONS
  private :: SnowAlbedo       ! Snow albedos
  private :: incident_direct  ! Direct beam solar rad incident on walls and road in urban canyon 
  private :: incident_diffuse ! Diffuse solar rad incident on walls and road in urban canyon
  private :: net_solar        ! Solar radiation absorbed by road and both walls in urban canyon 
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine UrbanAlbedo (bounds, num_urbanl, filter_urbanl, &
       num_urbanc, filter_urbanc, num_urbanp, filter_urbanp, &
       waterstate_vars, urbanparams_vars, solarabs_vars, surfalb_vars) 
    !
    ! !DESCRIPTION: 
    ! Determine urban landunit component albedos
    !
    ! Note that this is called with the "inactive_and_active" version of the filters, because
    ! the variables computed here are needed over inactive points that might later become
    ! active (due to landuse change). Thus, this routine cannot depend on variables that are
    ! only computed over active points.
    !
    ! !USES:
    use shr_orb_mod   , only : shr_orb_decl, shr_orb_cosz
    use elm_varcon    , only : sb
    use column_varcon , only : icol_roof, icol_sunwall, icol_shadewall
    use column_varcon , only : icol_road_perv, icol_road_imperv
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds  
    integer                , intent(in)    :: num_urbanl       ! number of urban landunits in clump
    integer                , intent(in)    :: filter_urbanl(:) ! urban landunit filter
    integer                , intent(in)    :: num_urbanc       ! number of urban columns in clump
    integer                , intent(in)    :: filter_urbanc(:) ! urban column filter
    integer                , intent(in)    :: num_urbanp       ! number of urban patches in clump
    integer                , intent(in)    :: filter_urbanp(:) ! urban pft filter
    type(waterstate_type)  , intent(in)    :: waterstate_vars
    type(urbanparams_type) , intent(inout) :: urbanparams_vars
    type(solarabs_type)    , intent(inout) :: solarabs_vars
    type(surfalb_type)     , intent(inout) :: surfalb_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: fl,fp,fc,g,l,p,c,ib                                  ! indices
    integer  :: ic                                                   ! 0=unit incoming direct; 1=unit incoming diffuse
    integer  :: num_solar                                            ! counter
    real(r8) :: coszen             (bounds%begl:bounds%endl)         ! cosine solar zenith angle for next time step (landunit)
    real(r8) :: zen                (bounds%begl:bounds%endl)         ! solar zenith angle (radians)
    real(r8) :: sdir               (bounds%begl:bounds%endl, numrad) ! direct beam solar radiation on horizontal surface
    real(r8) :: sdif               (bounds%begl:bounds%endl, numrad) ! diffuse solar radiation on horizontal surface
    real(r8) :: sdir_road          (bounds%begl:bounds%endl, numrad) ! direct beam solar radiation incident on road
    real(r8) :: sdif_road          (bounds%begl:bounds%endl, numrad) ! diffuse solar radiation incident on road
    real(r8) :: sdir_sunwall       (bounds%begl:bounds%endl, numrad) ! direct beam solar radiation (per unit wall area) incident on sunlit wall per unit incident flux
    real(r8) :: sdif_sunwall       (bounds%begl:bounds%endl, numrad) ! diffuse solar radiation (per unit wall area) incident on sunlit wall per unit incident flux
    real(r8) :: sdir_shadewall     (bounds%begl:bounds%endl, numrad) ! direct beam solar radiation (per unit wall area) incident on shaded wall per unit incident flux
    real(r8) :: sdif_shadewall     (bounds%begl:bounds%endl, numrad) ! diffuse solar radiation (per unit wall area) incident on shaded wall per unit incident flux
    real(r8) :: albsnd_roof        (bounds%begl:bounds%endl, numrad) ! snow albedo for roof (direct)
    real(r8) :: albsni_roof        (bounds%begl:bounds%endl, numrad) ! snow albedo for roof (diffuse)
    real(r8) :: albsnd_improad     (bounds%begl:bounds%endl, numrad) ! snow albedo for impervious road (direct)
    real(r8) :: albsni_improad     (bounds%begl:bounds%endl, numrad) ! snow albedo for impervious road (diffuse)
    real(r8) :: albsnd_perroad     (bounds%begl:bounds%endl, numrad) ! snow albedo for pervious road (direct)
    real(r8) :: albsni_perroad     (bounds%begl:bounds%endl, numrad) ! snow albedo for pervious road (diffuse)
    real(r8) :: alb_roof_dir_s     (bounds%begl:bounds%endl, numrad) ! direct roof albedo with snow effects
    real(r8) :: alb_roof_dif_s     (bounds%begl:bounds%endl, numrad) ! diffuse roof albedo with snow effects
    real(r8) :: alb_improad_dir_s  (bounds%begl:bounds%endl, numrad) ! direct impervious road albedo with snow effects
    real(r8) :: alb_perroad_dir_s  (bounds%begl:bounds%endl, numrad) ! direct pervious road albedo with snow effects
    real(r8) :: alb_improad_dif_s  (bounds%begl:bounds%endl, numrad) ! diffuse impervious road albedo with snow effects
    real(r8) :: alb_perroad_dif_s  (bounds%begl:bounds%endl, numrad) ! diffuse pervious road albedo with snow effects
    real(r8) :: sref_roof_dir      (bounds%begl:bounds%endl, numrad) ! direct  solar reflected by roof per unit ground area per unit incident flux   
    real(r8) :: sref_roof_dif      (bounds%begl:bounds%endl, numrad) ! diffuse solar reflected by roof per unit ground area per unit incident flux   
    real(r8) :: sref_sunwall_dir   (bounds%begl:bounds%endl, numrad) ! direct  solar reflected by sunwall per unit wall area per unit incident flux  
    real(r8) :: sref_sunwall_dif   (bounds%begl:bounds%endl, numrad) ! diffuse solar reflected by sunwall per unit wall area per unit incident flux  
    real(r8) :: sref_shadewall_dir (bounds%begl:bounds%endl, numrad) ! direct  solar reflected by shadewall per unit wall area per unit incident flux  
    real(r8) :: sref_shadewall_dif (bounds%begl:bounds%endl, numrad) ! diffuse solar reflected by shadewall per unit wall area per unit incident flux  
    real(r8) :: sref_improad_dir   (bounds%begl:bounds%endl, numrad) ! direct  solar reflected by impervious road per unit ground area per unit incident flux   
    real(r8) :: sref_improad_dif   (bounds%begl:bounds%endl, numrad) ! diffuse solar reflected by impervious road per unit ground area per unit incident flux   
    real(r8) :: sref_perroad_dir   (bounds%begl:bounds%endl, numrad) ! direct  solar reflected by pervious road per unit ground area per unit incident flux   
    real(r8) :: sref_perroad_dif   (bounds%begl:bounds%endl, numrad) ! diffuse solar reflected by pervious road per unit ground area per unit incident flux   
    !-----------------------------------------------------------------------

    associate(                                                        &
         ctype              => col_pp%itype                            , & ! Input:  [integer (:)    ]  column type                                        
         coli               => lun_pp%coli                             , & ! Input:  [integer (:)    ]  beginning column index for landunit                
         canyon_hwr         => lun_pp%canyon_hwr                       , & ! Input:  [real(r8) (:)   ]  ratio of building height to street width          
         wtroad_perv        => lun_pp%wtroad_perv                      , & ! Input:  [real(r8) (:)   ]  weight of pervious road wrt total road            
         
         frac_sno           => col_ws%frac_sno         , & ! Input:  [real(r8) (:)   ]  fraction of ground covered by snow (0 to 1)       
         
         alb_roof_dir       => urbanparams_vars%alb_roof_dir        , & ! Output: [real(r8) (:,:) ]  direct roof albedo                              
         alb_roof_dif       => urbanparams_vars%alb_roof_dif        , & ! Output: [real(r8) (:,:) ]  diffuse roof albedo                             
         alb_improad_dir    => urbanparams_vars%alb_improad_dir     , & ! Output: [real(r8) (:,:) ]  direct impervious road albedo                   
         alb_improad_dif    => urbanparams_vars%alb_improad_dif     , & ! Output: [real(r8) (:,:) ]  diffuse imprevious road albedo                  
         alb_perroad_dir    => urbanparams_vars%alb_perroad_dir     , & ! Output: [real(r8) (:,:) ]  direct pervious road albedo                     
         alb_perroad_dif    => urbanparams_vars%alb_perroad_dif     , & ! Output: [real(r8) (:,:) ]  diffuse pervious road albedo                    
         alb_wall_dir       => urbanparams_vars%alb_wall_dir        , & ! Output: [real(r8) (:,:) ]  direct wall albedo                              
         alb_wall_dif       => urbanparams_vars%alb_wall_dif        , & ! Output: [real(r8) (:,:) ]  diffuse wall albedo                             
        
         sabs_roof_dir      => solarabs_vars%sabs_roof_dir_lun      , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by roof per unit ground area per unit incident flux
         sabs_roof_dif      => solarabs_vars%sabs_roof_dif_lun      , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by roof per unit ground area per unit incident flux
         sabs_sunwall_dir   => solarabs_vars%sabs_sunwall_dir_lun   , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by sunwall per unit wall area per unit incident flux
         sabs_sunwall_dif   => solarabs_vars%sabs_sunwall_dif_lun   , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by sunwall per unit wall area per unit incident flux
         sabs_shadewall_dir => solarabs_vars%sabs_shadewall_dir_lun , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by shadewall per unit wall area per unit incident flux
         sabs_shadewall_dif => solarabs_vars%sabs_shadewall_dif_lun , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by shadewall per unit wall area per unit incident flux
         sabs_improad_dir   => solarabs_vars%sabs_improad_dir_lun   , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by impervious road per unit ground area per unit incident flux
         sabs_improad_dif   => solarabs_vars%sabs_improad_dif_lun   , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by impervious road per unit ground area per unit incident flux
         sabs_perroad_dir   => solarabs_vars%sabs_perroad_dir_lun   , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by pervious road per unit ground area per unit incident flux
         sabs_perroad_dif   => solarabs_vars%sabs_perroad_dif_lun   , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by pervious road per unit ground area per unit incident flux
         
         fabd               => surfalb_vars%fabd_patch              , & ! Output:  [real(r8) (:,:) ]  flux absorbed by canopy per unit direct flux
         fabd_sun           => surfalb_vars%fabd_sun_patch          , & ! Output:  [real(r8) (:,:) ]  flux absorbed by sunlit canopy per unit direct flux
         fabd_sha           => surfalb_vars%fabd_sha_patch          , & ! Output:  [real(r8) (:,:) ]  flux absorbed by shaded canopy per unit direct flux
         fabi               => surfalb_vars%fabi_patch              , & ! Output:  [real(r8) (:,:) ]  flux absorbed by canopy per unit diffuse flux
         fabi_sun           => surfalb_vars%fabi_sun_patch          , & ! Output:  [real(r8) (:,:) ]  flux absorbed by sunlit canopy per unit diffuse flux
         fabi_sha           => surfalb_vars%fabi_sha_patch          , & ! Output:  [real(r8) (:,:) ]  flux absorbed by shaded canopy per unit diffuse flux
         ftdd               => surfalb_vars%ftdd_patch              , & ! Output:  [real(r8) (:,:) ]  down direct flux below canopy per unit direct flux
         ftid               => surfalb_vars%ftid_patch              , & ! Output:  [real(r8) (:,:) ]  down diffuse flux below canopy per unit direct flux
         ftii               => surfalb_vars%ftii_patch              , & ! Output:  [real(r8) (:,:) ]  down diffuse flux below canopy per unit diffuse flux
         albgrd             => surfalb_vars%albgrd_col              , & ! Output: [real(r8) (:,:) ]  urban col ground albedo (direct) 
         albgri             => surfalb_vars%albgri_col              , & ! Output: [real(r8) (:,:) ]  urban col ground albedo (diffuse)
         albd               => surfalb_vars%albd_patch              , & ! Output  [real(r8) (:,:) ]  urban pft surface albedo (direct)                         
         albi               => surfalb_vars%albi_patch              , & ! Output: [real(r8) (:,:) ]  urban pft surface albedo (diffuse)                        
         
         begl               => bounds%begl                          , &
         endl               => bounds%endl                            &
         )

      ! ----------------------------------------------------------------------------
      ! Solar declination and cosine solar zenith angle and zenith angle for 
      ! next time step
      ! ----------------------------------------------------------------------------
      
      do fl = 1,num_urbanl
         l = filter_urbanl(fl)
         g = lun_pp%gridcell(l)
         coszen(l) = surfalb_vars%coszen_col(coli(l))  ! Assumes coszen for each column are the same
         zen(l)    = acos(coszen(l))
      end do

       ! Initialize output because solar radiation only done if coszen > 0

      do ib = 1, numrad
         do fc = 1,num_urbanc
            c = filter_urbanc(fc)
            albgrd(c,ib) = 0._r8
            albgri(c,ib) = 0._r8
         end do

         do fp = 1,num_urbanp  
            p = filter_urbanp(fp)
            l = veg_pp%landunit(p)
            c = veg_pp%column(p)

            ! Initialize direct and diffuse albedo such that if the Sun is below
            ! the horizon, p2g scaling returns an albedo of 1.0.
            if (col_pp%itype(c) == icol_sunwall) then
               albd(p,ib) = 1._r8 / (3.0_r8 * lun_pp%canyon_hwr(l))
               albi(p,ib) = 1._r8 / (3.0_r8 * lun_pp%canyon_hwr(l))
            else if (col_pp%itype(c) == icol_shadewall) then
               albd(p,ib) = 1._r8 / (3.0_r8 * lun_pp%canyon_hwr(l))
               albi(p,ib) = 1._r8 / (3.0_r8 * lun_pp%canyon_hwr(l))
            else if (col_pp%itype(c) == icol_road_perv .or. col_pp%itype(c) == icol_road_imperv) then
               albd(p,ib) = 1._r8 / (3.0_r8)
               albi(p,ib) = 1._r8 / (3.0_r8)
            else if (col_pp%itype(c) == icol_roof) then
               albd(p,ib) = 1._r8
               albi(p,ib) = 1._r8
            endif

            fabd(p,ib)     = 0._r8
            fabd_sun(p,ib) = 0._r8
            fabd_sha(p,ib) = 0._r8
            fabi(p,ib)     = 0._r8
            fabi_sun(p,ib) = 0._r8
            fabi_sha(p,ib) = 0._r8
            if (coszen(l) > 0._r8) then
               ftdd(p,ib)  = 1._r8
            else
               ftdd(p,ib)  = 0._r8
            end if
            ftid(p,ib)     = 0._r8
            if (coszen(l) > 0._r8) then
               ftii(p,ib)  = 1._r8
            else
               ftii(p,ib)  = 0._r8
            end if
         end do
      end do

      ! ----------------------------------------------------------------------------
      ! Urban Code
      ! ----------------------------------------------------------------------------

      num_solar = 0
      do fl = 1,num_urbanl
         l = filter_urbanl(fl)
         if (coszen(l) > 0._r8) num_solar = num_solar + 1
      end do

      do ib = 1,numrad
         do fl = 1,num_urbanl
            l = filter_urbanl(fl)
            sabs_roof_dir(l,ib)      = 0._r8
            sabs_roof_dif(l,ib)      = 0._r8
            sabs_sunwall_dir(l,ib)   = 0._r8
            sabs_sunwall_dif(l,ib)   = 0._r8
            sabs_shadewall_dir(l,ib) = 0._r8
            sabs_shadewall_dif(l,ib) = 0._r8
            sabs_improad_dir(l,ib)   = 0._r8
            sabs_improad_dif(l,ib)   = 0._r8
            sabs_perroad_dir(l,ib)   = 0._r8
            sabs_perroad_dif(l,ib)   = 0._r8
            sref_roof_dir(l,ib)      = 1._r8
            sref_roof_dif(l,ib)      = 1._r8
            sref_sunwall_dir(l,ib)   = 1._r8
            sref_sunwall_dif(l,ib)   = 1._r8
            sref_shadewall_dir(l,ib) = 1._r8
            sref_shadewall_dif(l,ib) = 1._r8
            sref_improad_dir(l,ib)   = 1._r8
            sref_improad_dif(l,ib)   = 1._r8
            sref_perroad_dir(l,ib)   = 1._r8
            sref_perroad_dif(l,ib)   = 1._r8

            ! Initialize direct and diffuse albedos such that if the Sun is below
            ! the horizon, p2g scaling returns an albedo of 1.0.
            sref_sunwall_dir(l,ib)   = 1._r8 / (3.0_r8 * lun_pp%canyon_hwr(l))
            sref_sunwall_dif(l,ib)   = 1._r8 / (3.0_r8 * lun_pp%canyon_hwr(l))
            sref_shadewall_dir(l,ib) = 1._r8 / (3.0_r8 * lun_pp%canyon_hwr(l))
            sref_shadewall_dif(l,ib) = 1._r8 / (3.0_r8 * lun_pp%canyon_hwr(l))
            sref_improad_dir(l,ib)   = 1._r8 / (3.0_r8)
            sref_improad_dif(l,ib)   = 1._r8 / (3.0_r8)
            sref_perroad_dir(l,ib)   = 1._r8 / (3.0_r8)
            sref_perroad_dif(l,ib)   = 1._r8 / (3.0_r8)
            sref_roof_dir(l,ib)      = 1._r8
            sref_roof_dif(l,ib)      = 1._r8
         end do
      end do

      ! ----------------------------------------------------------------------------
      ! Only do the rest if all coszen are positive 
      ! ----------------------------------------------------------------------------

      if (num_solar > 0)then

         ! Set constants - solar fluxes are per unit incoming flux

         do ib = 1,numrad
            do fl = 1,num_urbanl
               l = filter_urbanl(fl)
               sdir(l,ib) = 1._r8
               sdif(l,ib) = 1._r8
            end do
         end do

         ! Incident direct beam radiation for 
         ! (a) roof and (b) road and both walls in urban canyon

         if (num_urbanl > 0) then
            call incident_direct (bounds, &
                 num_urbanl, filter_urbanl, &
                 canyon_hwr(begl:endl), &
                 coszen(begl:endl), &
                 zen(begl:endl), &
                 sdir(begl:endl, :), &
                 sdir_road(begl:endl, :), &
                 sdir_sunwall(begl:endl, :), &
                 sdir_shadewall(begl:endl, :))
         end if

         ! Incident diffuse radiation for 
         ! (a) roof and (b) road and both walls in urban canyon.

         if (num_urbanl > 0) then
            call incident_diffuse (bounds, &
                 num_urbanl, filter_urbanl, &
                 canyon_hwr(begl:endl), &
                 sdif(begl:endl, :), &
                 sdif_road(begl:endl, :), &
                 sdif_sunwall(begl:endl, :), &
                 sdif_shadewall(begl:endl, :), &
                 urbanparams_vars)
         end if

         ! Get snow albedos for roof and impervious and pervious road
         if (num_urbanl > 0) then
            ic = 0
            call SnowAlbedo(bounds, &
                 num_urbanc, filter_urbanc, &
                 coszen(begl:endl), &
                 ic, &
                 albsnd_roof(begl:endl, :), &
                 albsnd_improad(begl:endl, :), &
                 albsnd_perroad(begl:endl, :), &
                 waterstate_vars)

            ic = 1
            call SnowAlbedo(bounds, &
                 num_urbanc, filter_urbanc, &
                 coszen(begl:endl), &
                 ic, &
                 albsni_roof(begl:endl, :), &
                 albsni_improad(begl:endl, :), &
                 albsni_perroad(begl:endl, :), &
                 waterstate_vars)
         end if

         ! Combine snow-free and snow albedos
         do ib = 1,numrad
            do fc = 1,num_urbanc
               c = filter_urbanc(fc)
               l = col_pp%landunit(c)
               if (ctype(c) == icol_roof) then    
                  alb_roof_dir_s(l,ib) = alb_roof_dir(l,ib)*(1._r8-frac_sno(c))  &
                       + albsnd_roof(l,ib)*frac_sno(c)
                  alb_roof_dif_s(l,ib) = alb_roof_dif(l,ib)*(1._r8-frac_sno(c))  &
                       + albsni_roof(l,ib)*frac_sno(c)
               else if (ctype(c) == icol_road_imperv) then    
                  alb_improad_dir_s(l,ib) = alb_improad_dir(l,ib)*(1._r8-frac_sno(c))  &
                       + albsnd_improad(l,ib)*frac_sno(c)
                  alb_improad_dif_s(l,ib) = alb_improad_dif(l,ib)*(1._r8-frac_sno(c))  &
                       + albsni_improad(l,ib)*frac_sno(c)
               else if (ctype(c) == icol_road_perv) then    
                  alb_perroad_dir_s(l,ib) = alb_perroad_dir(l,ib)*(1._r8-frac_sno(c))  &
                       + albsnd_perroad(l,ib)*frac_sno(c)
                  alb_perroad_dif_s(l,ib) = alb_perroad_dif(l,ib)*(1._r8-frac_sno(c))  &
                       + albsni_perroad(l,ib)*frac_sno(c)
               end if
            end do
         end do

         ! Reflected and absorbed solar radiation per unit incident radiation 
         ! for road and both walls in urban canyon allowing for multiple reflection
         ! Reflected and absorbed solar radiation per unit incident radiation for roof

         if (num_urbanl > 0) then
            call net_solar (bounds, &
                 num_urbanl, filter_urbanl, &
                 coszen             (begl:endl), &
                 canyon_hwr         (begl:endl), &
                 wtroad_perv        (begl:endl), &
                 sdir               (begl:endl, :), &
                 sdif               (begl:endl, :), &
                 alb_improad_dir_s  (begl:endl, :), &
                 alb_perroad_dir_s  (begl:endl, :), &
                 alb_wall_dir       (begl:endl, :), &
                 alb_roof_dir_s     (begl:endl, :), &
                 alb_improad_dif_s  (begl:endl, :), &
                 alb_perroad_dif_s  (begl:endl, :), &
                 alb_wall_dif       (begl:endl, :), &
                 alb_roof_dif_s     (begl:endl, :), &
                 sdir_road          (begl:endl, :), &
                 sdir_sunwall       (begl:endl, :), &
                 sdir_shadewall     (begl:endl, :),  &
                 sdif_road          (begl:endl, :), &
                 sdif_sunwall       (begl:endl, :), &
                 sdif_shadewall     (begl:endl, :),  &
                 sref_improad_dir   (begl:endl, :), &
                 sref_perroad_dir   (begl:endl, :), &
                 sref_sunwall_dir   (begl:endl, :), &
                 sref_shadewall_dir (begl:endl, :), &
                 sref_roof_dir      (begl:endl, :), &
                 sref_improad_dif   (begl:endl, :), &
                 sref_perroad_dif   (begl:endl, :), &
                 sref_sunwall_dif   (begl:endl, :), &
                 sref_shadewall_dif (begl:endl, :), &
                 sref_roof_dif      (begl:endl, :), &
                 urbanparams_vars, solarabs_vars)
         end if

         ! ----------------------------------------------------------------------------
         ! Map urban output to surfalb_vars components 
         ! ----------------------------------------------------------------------------

         !  Set albgrd and albgri (ground albedos) and albd and albi (surface albedos)

         do ib = 1,numrad
            do fc = 1,num_urbanc
               c = filter_urbanc(fc)
               l = col_pp%landunit(c)
               if (ctype(c) == icol_roof) then    
                  albgrd(c,ib) = sref_roof_dir(l,ib) 
                  albgri(c,ib) = sref_roof_dif(l,ib) 
               else if (ctype(c) == icol_sunwall) then   
                  albgrd(c,ib) = sref_sunwall_dir(l,ib)
                  albgri(c,ib) = sref_sunwall_dif(l,ib)
               else if (ctype(c) == icol_shadewall) then 
                  albgrd(c,ib) = sref_shadewall_dir(l,ib)
                  albgri(c,ib) = sref_shadewall_dif(l,ib)
               else if (ctype(c) == icol_road_perv) then
                  albgrd(c,ib) = sref_perroad_dir(l,ib)
                  albgri(c,ib) = sref_perroad_dif(l,ib)
               else if (ctype(c) == icol_road_imperv) then
                  albgrd(c,ib) = sref_improad_dir(l,ib)
                  albgri(c,ib) = sref_improad_dif(l,ib)
               endif
            end do
            do fp = 1,num_urbanp
               p = filter_urbanp(fp)
               c = veg_pp%column(p)
               albd(p,ib) = albgrd(c,ib)
               albi(p,ib) = albgri(c,ib)
            end do
         end do
      end if

    end associate

  end subroutine UrbanAlbedo

  !-----------------------------------------------------------------------
  subroutine SnowAlbedo (bounds          , &
       num_urbanc, filter_urbanc, coszen, ind , &
       albsn_roof, albsn_improad, albsn_perroad, &
       waterstate_vars)
    !
    ! !DESCRIPTION:
    ! Determine urban snow albedos
    !
    ! !USES:
    use column_varcon, only : icol_roof, icol_road_perv, icol_road_imperv
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds                     
    integer , intent(in) :: num_urbanc                          ! number of urban columns in clump
    integer , intent(in) :: filter_urbanc(:)                    ! urban column filter
    integer , intent(in) :: ind                                 ! 0=direct beam, 1=diffuse radiation
    real(r8), intent(in) :: coszen        ( bounds%begl: )      ! cosine solar zenith angle [landunit]
    real(r8), intent(out):: albsn_roof    ( bounds%begl: , 1: ) ! roof snow albedo by waveband [landunit, numrad]
    real(r8), intent(out):: albsn_improad ( bounds%begl: , 1: ) ! impervious road snow albedo by waveband [landunit, numrad]
    real(r8), intent(out):: albsn_perroad ( bounds%begl: , 1: ) ! pervious road snow albedo by waveband [landunit, numrad]
    type(waterstate_type), intent(in) :: waterstate_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: fc,c,l              ! indices
    !
    ! These values are derived from Marshall (1989) assuming soot content of 1.5e-5 
    ! (three times what LSM uses globally). Note that snow age effects are ignored here.
    real(r8), parameter :: snal0 = 0.66_r8 ! vis albedo of urban snow
    real(r8), parameter :: snal1 = 0.56_r8 ! nir albedo of urban snow
    !-----------------------------------------------------------------------

    ! this code assumes that numrad = 2 , with the following
    ! index values: 1 = visible, 2 = NIR
    SHR_ASSERT_ALL(numrad == 2, errMsg(__FILE__, __LINE__))

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(coszen)        == (/bounds%endl/)),         errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(albsn_roof)    == (/bounds%endl, numrad/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(albsn_improad) == (/bounds%endl, numrad/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(albsn_perroad) == (/bounds%endl, numrad/)), errMsg(__FILE__, __LINE__))

    associate(                            & 
         h2osno =>  col_ws%h2osno & ! Input:  [real(r8) (:) ]  snow water (mm H2O)                               
         )
      
      do fc = 1,num_urbanc
         c = filter_urbanc(fc)
         l = col_pp%landunit(c)
         if (coszen(l) > 0._r8 .and. h2osno(c) > 0._r8) then
            if (col_pp%itype(c) == icol_roof) then
               albsn_roof(l,1) = snal0
               albsn_roof(l,2) = snal1
            else if (col_pp%itype(c) == icol_road_imperv) then
               albsn_improad(l,1) = snal0
               albsn_improad(l,2) = snal1
            else if (col_pp%itype(c) == icol_road_perv) then
               albsn_perroad(l,1) = snal0
               albsn_perroad(l,2) = snal1
            end if
         else
            if (col_pp%itype(c) == icol_roof) then
               albsn_roof(l,1) = 0._r8
               albsn_roof(l,2) = 0._r8
            else if (col_pp%itype(c) == icol_road_imperv) then
               albsn_improad(l,1) = 0._r8
               albsn_improad(l,2) = 0._r8
            else if (col_pp%itype(c) == icol_road_perv) then
               albsn_perroad(l,1) = 0._r8
               albsn_perroad(l,2) = 0._r8
            end if
         end if
      end do

    end associate

  end subroutine SnowAlbedo

  !-----------------------------------------------------------------------
  subroutine incident_direct (bounds                      ,                    &
       num_urbanl, filter_urbanl, canyon_hwr, coszen, zen ,                    &
       sdir, sdir_road, sdir_sunwall, sdir_shadewall)
    !
    ! !DESCRIPTION: 
    ! Direct beam solar radiation incident on walls and road in urban canyon
    !
    !                           Sun
    !                            /
    !             roof          /
    !            ------        /---            -
    !                 |       / |              |
    !    sunlit wall  |      /  | shaded wall  h
    !                 |     /   |              |
    !                 -----/-----              -
    !                    road
    !                 <--- w --->
    !
    ! Method:
    ! Road          = Horizontal surface. Account for shading by wall. Integrate over all canyon orientations
    ! Wall (sunlit) = Adjust horizontal radiation for 90 degree surface. Account for shading by opposing wall.
    !                 Integrate over all canyon orientations
    ! Wall (shaded) = 0
    !
    ! Conservation check: Total incoming direct beam (sdir) = sdir_road + (sdir_shadewall + sdir_sunwall)*canyon_hwr
    ! Multiplication by canyon_hwr scales wall fluxes (per unit wall area) to per unit ground area
    !
    ! Source: Masson, V. (2000) A physically-based scheme for the urban energy budget in 
    ! atmospheric models. Boundary-Layer Meteorology 94:357-397
    !
    ! This analytical solution from Masson (2000) agrees with the numerical solution to
    ! within 0.6 W/m**2 for sdir = 1000 W/m**2 and for all H/W from 0.1 to 10 by 0.1
    ! and all solar zenith angles from 1 to 90 deg by 1
    !
    ! !USES:
    use elm_varcon, only : rpi
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds                      
    integer , intent(in)  :: num_urbanl                          ! number of urban landunits
    integer , intent(in)  :: filter_urbanl(:)                    ! urban landunit filter
    real(r8), intent(in)  :: canyon_hwr( bounds%begl: )          ! ratio of building height to street width [landunit]
    real(r8), intent(in)  :: coszen( bounds%begl: )              ! cosine solar zenith angle [landunit]
    real(r8), intent(in)  :: zen( bounds%begl: )                 ! solar zenith angle (radians) [landunit]
    real(r8), intent(in)  :: sdir( bounds%begl: , 1: )           ! direct beam solar radiation incident on horizontal surface [landunit, numrad]
    real(r8), intent(out) :: sdir_road( bounds%begl: , 1: )      ! direct beam solar radiation incident on road per unit incident flux [landunit, numrad]
    real(r8), intent(out) :: sdir_sunwall( bounds%begl: , 1: )   ! direct beam solar radiation (per unit wall area) incident on sunlit wall per unit incident flux [landunit, numrad]
    real(r8), intent(out) :: sdir_shadewall( bounds%begl: , 1: ) ! direct beam solar radiation (per unit wall area) incident on shaded wall per unit incident flux [landunit, numrad]
    !
    ! !LOCAL VARIABLES:
    integer  :: fl,l,i,ib                       ! indices
    logical  :: numchk = .false.                ! true => perform numerical check of analytical solution
    real(r8) :: theta0(bounds%begl:bounds%endl) ! critical canyon orientation for which road is no longer illuminated
    real(r8) :: tanzen(bounds%begl:bounds%endl) ! tan(zenith angle)
    real(r8) :: swall_projected                 ! direct beam solar radiation (per unit ground area) incident on wall
    real(r8) :: err1(bounds%begl:bounds%endl)   ! energy conservation error
    real(r8) :: err2(bounds%begl:bounds%endl)   ! energy conservation error
    real(r8) :: err3(bounds%begl:bounds%endl)   ! energy conservation error
    real(r8) :: sumr                            ! sum of sroad for each orientation (0 <= theta <= pi/2)
    real(r8) :: sumw                            ! sum of swall for each orientation (0 <= theta <= pi/2)
    real(r8) :: num                             ! number of orientations
    real(r8) :: theta                           ! canyon orientation relative to sun (0 <= theta <= pi/2)
    real(r8) :: zen0                            ! critical solar zenith angle for which sun begins to illuminate road
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(canyon_hwr)     == (/bounds%endl/)),         errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(coszen)         == (/bounds%endl/)),         errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(zen)            == (/bounds%endl/)),         errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(sdir)           == (/bounds%endl, numrad/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(sdir_road)      == (/bounds%endl, numrad/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(sdir_sunwall)   == (/bounds%endl, numrad/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(sdir_shadewall) == (/bounds%endl, numrad/)), errMsg(__FILE__, __LINE__))

    do fl = 1,num_urbanl
       l = filter_urbanl(fl)
       if (coszen(l) > 0._r8) then
          theta0(l) = asin(min( (1._r8/(canyon_hwr(l)*tan(max(zen(l),0.000001_r8)))), 1._r8 ))
          tanzen(l) = tan(zen(l))
       end if
    end do

    do ib = 1,numrad

       do fl = 1,num_urbanl
          l = filter_urbanl(fl)
          if (coszen(l) > 0._r8) then
             sdir_shadewall(l,ib) = 0._r8

             ! incident solar radiation on wall and road integrated over all canyon orientations (0 <= theta <= pi/2)

             sdir_road(l,ib) = sdir(l,ib) *                                    &
                  (2._r8*theta0(l)/rpi - 2./rpi*canyon_hwr(l)*tanzen(l)*(1._r8-cos(theta0(l))))
             sdir_sunwall(l,ib) = 2._r8 * sdir(l,ib) * ((1._r8/canyon_hwr(l))* &
                  (0.5_r8-theta0(l)/rpi) + (1._r8/rpi)*tanzen(l)*(1._r8-cos(theta0(l))))

             ! conservation check for road and wall. need to use wall fluxes converted to ground area

             swall_projected = (sdir_shadewall(l,ib) + sdir_sunwall(l,ib)) * canyon_hwr(l)
             err1(l) = sdir(l,ib) - (sdir_road(l,ib) + swall_projected)
          else
             sdir_road(l,ib) = 0._r8
             sdir_sunwall(l,ib) = 0._r8
             sdir_shadewall(l,ib) = 0._r8
          endif
       end do

       do fl = 1,num_urbanl
          l = filter_urbanl(fl)
          if (coszen(l) > 0._r8) then
             if (abs(err1(l)) > 0.001_r8) then
                write (iulog,*) 'urban direct beam solar radiation balance error',err1(l)
                write (iulog,*) 'clm model is stopping'
                call endrun(decomp_index=l, clmlevel=namel, msg=errmsg(__FILE__, __LINE__))
             endif
          endif
       end do

       ! numerical check of analytical solution
       ! sum sroad and swall over all canyon orientations (0 <= theta <= pi/2)

       if (numchk) then
          do fl = 1,num_urbanl
             l = filter_urbanl(fl)
             if (coszen(l) > 0._r8) then
                sumr = 0._r8
                sumw = 0._r8
                num  = 0._r8
                do i = 1, 9000
                   theta = i/100._r8 * rpi/180._r8
                   zen0 = atan(1._r8/(canyon_hwr(l)*sin(theta)))
                   if (zen(l) >= zen0) then 
                      sumr = sumr + 0._r8
                      sumw = sumw + sdir(l,ib) / canyon_hwr(l)
                   else
                      sumr = sumr + sdir(l,ib) * (1._r8-canyon_hwr(l)*sin(theta)*tanzen(l))
                      sumw = sumw + sdir(l,ib) * sin(theta)*tanzen(l)
                   end if
                   num = num + 1._r8
                end do
                err2(l) = sumr/num - sdir_road(l,ib)
                err3(l) = sumw/num - sdir_sunwall(l,ib)
             endif
          end do
          do fl = 1,num_urbanl
             l = filter_urbanl(fl)
             if (coszen(l) > 0._r8) then
                if (abs(err2(l)) > 0.0006_r8 ) then
                   write (iulog,*) 'urban road incident direct beam solar radiation error',err2(l)
                   write (iulog,*) 'clm model is stopping'
                   call endrun(decomp_index=l, clmlevel=namel, msg=errmsg(__FILE__, __LINE__))
                endif
                if (abs(err3(l)) > 0.0006_r8 ) then
                   write (iulog,*) 'urban wall incident direct beam solar radiation error',err3(l)
                   write (iulog,*) 'clm model is stopping'
                   call endrun(decomp_index=l, clmlevel=namel, msg=errmsg(__FILE__, __LINE__))
                end if
             end if
          end do
       end if

    end do

  end subroutine incident_direct

  !-----------------------------------------------------------------------
  subroutine incident_diffuse (bounds, &
       num_urbanl, filter_urbanl, canyon_hwr, &
       sdif, sdif_road, sdif_sunwall, sdif_shadewall, &
       urbanparams_vars)
    !
    ! !DESCRIPTION: 
    ! Diffuse solar radiation incident on walls and road in urban canyon
    ! Conservation check: Total incoming diffuse 
    ! (sdif) = sdif_road + (sdif_shadewall + sdif_sunwall)*canyon_hwr
    ! Multiplication by canyon_hwr scales wall fluxes (per unit wall area) to per unit ground area
    !
    ! !ARGUMENTS:
    type(bounds_type)     , intent(in)  :: bounds                      
    integer               , intent(in)  :: num_urbanl                           ! number of urban landunits
    integer               , intent(in)  :: filter_urbanl(:)                     ! urban landunit filter
    real(r8)              , intent(in)  :: canyon_hwr     ( bounds%begl: )      ! ratio of building height to street width [landunit]
    real(r8)              , intent(in)  :: sdif           ( bounds%begl: , 1: ) ! diffuse solar radiation incident on horizontal surface [landunit, numrad]
    real(r8)              , intent(out) :: sdif_road      ( bounds%begl: , 1: ) ! diffuse solar radiation incident on road [landunit, numrad]
    real(r8)              , intent(out) :: sdif_sunwall   ( bounds%begl: , 1: ) ! diffuse solar radiation (per unit wall area) incident on sunlit wall [landunit, numrad]
    real(r8)              , intent(out) :: sdif_shadewall ( bounds%begl: , 1: ) ! diffuse solar radiation (per unit wall area) incident on shaded wall [landunit, numrad]
    type(urbanparams_type), intent(in)  :: urbanparams_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: l, fl, ib       ! indices      
    real(r8) :: err(bounds%begl:bounds%endl)    ! energy conservation error (W/m**2)
    real(r8) :: swall_projected ! diffuse solar radiation (per unit ground area) incident on wall (W/m**2)
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(canyon_hwr)     == (/bounds%endl/)),         errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(sdif)           == (/bounds%endl, numrad/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(sdif_road)      == (/bounds%endl, numrad/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(sdif_sunwall)   == (/bounds%endl, numrad/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(sdif_shadewall) == (/bounds%endl, numrad/)), errMsg(__FILE__, __LINE__))

    associate(                            & 
         vf_sr =>    urbanparams_vars%vf_sr , & ! Input:  [real(r8) (:) ]  view factor of sky for road                       
         vf_sw =>    urbanparams_vars%vf_sw   & ! Input:  [real(r8) (:) ]  view factor of sky for one wall                   
         )

      do ib = 1, numrad
         
         ! diffuse solar and conservation check. need to convert wall fluxes to ground area
         
         do fl = 1,num_urbanl
            l = filter_urbanl(fl)
            sdif_road(l,ib)      = sdif(l,ib) * vf_sr(l)
            sdif_sunwall(l,ib)   = sdif(l,ib) * vf_sw(l) 
            sdif_shadewall(l,ib) = sdif(l,ib) * vf_sw(l) 

            swall_projected = (sdif_shadewall(l,ib) + sdif_sunwall(l,ib)) * canyon_hwr(l)
            err(l) = sdif(l,ib) - (sdif_road(l,ib) + swall_projected)
         end do

         ! error check

         do fl = 1, num_urbanl
            l = filter_urbanl(fl)
            if (abs(err(l)) > 0.001_r8) then
               write (iulog,*) 'urban diffuse solar radiation balance error',err(l) 
               write (iulog,*) 'clm model is stopping'
               call endrun(decomp_index=l, clmlevel=namel, msg=errmsg(__FILE__, __LINE__))
            endif
         end do

      end do

    end associate

  end subroutine incident_diffuse

  !-----------------------------------------------------------------------
  subroutine net_solar (bounds                                                                 , &
       num_urbanl, filter_urbanl, coszen, canyon_hwr, wtroad_perv, sdir, sdif                  , &
       alb_improad_dir, alb_perroad_dir, alb_wall_dir, alb_roof_dir                            , &
       alb_improad_dif, alb_perroad_dif, alb_wall_dif, alb_roof_dif                            , &
       sdir_road, sdir_sunwall, sdir_shadewall,                                                  &
       sdif_road, sdif_sunwall, sdif_shadewall,                                                  &
       sref_improad_dir, sref_perroad_dir, sref_sunwall_dir, sref_shadewall_dir, sref_roof_dir , &
       sref_improad_dif, sref_perroad_dif, sref_sunwall_dif, sref_shadewall_dif, sref_roof_dif , &
       urbanparams_vars, solarabs_vars)
    !
    ! !DESCRIPTION: 
    ! Solar radiation absorbed by road and both walls in urban canyon allowing 
    ! for multiple reflection.
    !
    ! !ARGUMENTS:
    type (bounds_type), intent(in) :: bounds                            
    integer , intent(in)    :: num_urbanl                               ! number of urban landunits
    integer , intent(in)    :: filter_urbanl(:)                         ! urban landunit filter
    real(r8), intent(in)    :: coszen             ( bounds%begl: )      ! cosine solar zenith angle [landunit]
    real(r8), intent(in)    :: canyon_hwr         ( bounds%begl: )      ! ratio of building height to street width [landunit]
    real(r8), intent(in)    :: wtroad_perv        ( bounds%begl: )      ! weight of pervious road wrt total road [landunit]
    real(r8), intent(in)    :: sdir               ( bounds%begl: , 1: ) ! direct beam solar radiation incident on horizontal surface [landunit, numrad]
    real(r8), intent(in)    :: sdif               ( bounds%begl: , 1: ) ! diffuse solar radiation on horizontal surface [landunit, numrad]
    real(r8), intent(in)    :: alb_improad_dir    ( bounds%begl: , 1: ) ! direct impervious road albedo [landunit, numrad]
    real(r8), intent(in)    :: alb_perroad_dir    ( bounds%begl: , 1: ) ! direct pervious road albedo [landunit, numrad]
    real(r8), intent(in)    :: alb_wall_dir       ( bounds%begl: , 1: ) ! direct  wall albedo [landunit, numrad]
    real(r8), intent(in)    :: alb_roof_dir       ( bounds%begl: , 1: ) ! direct  roof albedo [landunit, numrad]
    real(r8), intent(in)    :: alb_improad_dif    ( bounds%begl: , 1: ) ! diffuse impervious road albedo [landunit, numrad]
    real(r8), intent(in)    :: alb_perroad_dif    ( bounds%begl: , 1: ) ! diffuse pervious road albedo [landunit, numrad]
    real(r8), intent(in)    :: alb_wall_dif       ( bounds%begl: , 1: ) ! diffuse wall albedo [landunit, numrad]
    real(r8), intent(in)    :: alb_roof_dif       ( bounds%begl: , 1: ) ! diffuse roof albedo [landunit, numrad]
    real(r8), intent(in)    :: sdir_road          ( bounds%begl: , 1: ) ! direct beam solar radiation incident on road per unit incident flux [landunit, numrad]
    real(r8), intent(in)    :: sdir_sunwall       ( bounds%begl: , 1: ) ! direct beam solar radiation (per unit wall area) incident on sunlit wall per unit incident flux [landunit, numrad]
    real(r8), intent(in)    :: sdir_shadewall     ( bounds%begl: , 1: ) ! direct beam solar radiation (per unit wall area) incident on shaded wall per unit incident flux [landunit, numrad]
    real(r8), intent(in)    :: sdif_road          ( bounds%begl: , 1: ) ! diffuse solar radiation incident on road per unit incident flux [landunit, numrad]
    real(r8), intent(in)    :: sdif_sunwall       ( bounds%begl: , 1: ) ! diffuse solar radiation (per unit wall area) incident on sunlit wall per unit incident flux [landunit, numrad]
    real(r8), intent(in)    :: sdif_shadewall     ( bounds%begl: , 1: ) ! diffuse solar radiation (per unit wall area) incident on shaded wall per unit incident flux [landunit, numrad]
    real(r8), intent(inout) :: sref_improad_dir   ( bounds%begl: , 1: ) ! direct  solar rad reflected by impervious road (per unit ground area) per unit incident flux [landunit, numrad]
    real(r8), intent(inout) :: sref_perroad_dir   ( bounds%begl: , 1: ) ! direct  solar rad reflected by pervious road (per unit ground area) per unit incident flux [landunit, numrad]
    real(r8), intent(inout) :: sref_improad_dif   ( bounds%begl: , 1: ) ! diffuse solar rad reflected by impervious road (per unit ground area) per unit incident flux [landunit, numrad]
    real(r8), intent(inout) :: sref_perroad_dif   ( bounds%begl: , 1: ) ! diffuse solar rad reflected by pervious road (per unit ground area) per unit incident flux [landunit, numrad]
    real(r8), intent(inout) :: sref_sunwall_dir   ( bounds%begl: , 1: ) ! direct solar  rad reflected by sunwall (per unit wall area) per unit incident flux [landunit, numrad]
    real(r8), intent(inout) :: sref_sunwall_dif   ( bounds%begl: , 1: ) ! diffuse solar rad reflected by sunwall (per unit wall area) per unit incident flux [landunit, numrad]
    real(r8), intent(inout) :: sref_shadewall_dir ( bounds%begl: , 1: ) ! direct solar  rad reflected by shadewall (per unit wall area) per unit incident flux [landunit, numrad]
    real(r8), intent(inout) :: sref_shadewall_dif ( bounds%begl: , 1: ) ! diffuse solar rad reflected by shadewall (per unit wall area) per unit incident flux [landunit, numrad]
    real(r8), intent(inout) :: sref_roof_dir      ( bounds%begl: , 1: ) ! direct  solar rad reflected by roof (per unit ground area) per unit incident flux [landunit, numrad]
    real(r8), intent(inout) :: sref_roof_dif      ( bounds%begl: , 1: ) ! diffuse solar rad reflected by roof (per unit ground area)  per unit incident flux [landunit, numrad]
    type(urbanparams_type), intent(in)    :: urbanparams_vars
    type(solarabs_type)   , intent(inout) :: solarabs_vars
    !
    ! !LOCAL VARIABLES
    real(r8) :: wtroad_imperv(bounds%begl:bounds%endl)           ! weight of impervious road wrt total road
    real(r8) :: sabs_canyon_dir(bounds%begl:bounds%endl)         ! direct solar rad absorbed by canyon per unit incident flux
    real(r8) :: sabs_canyon_dif(bounds%begl:bounds%endl)         ! diffuse solar rad absorbed by canyon per unit incident flux
    real(r8) :: sref_canyon_dir(bounds%begl:bounds%endl)         ! direct solar reflected by canyon per unit incident flux 
    real(r8) :: sref_canyon_dif(bounds%begl:bounds%endl)         ! diffuse solar reflected by canyon per unit incident flux

    real(r8) :: improad_a_dir(bounds%begl:bounds%endl)           ! absorbed direct solar for impervious road after "n" reflections per unit incident flux
    real(r8) :: improad_a_dif(bounds%begl:bounds%endl)           ! absorbed diffuse solar for impervious road after "n" reflections per unit incident flux
    real(r8) :: improad_r_dir(bounds%begl:bounds%endl)           ! reflected direct solar for impervious road after "n" reflections per unit incident flux
    real(r8) :: improad_r_dif(bounds%begl:bounds%endl)           ! reflected diffuse solar for impervious road after "n" reflections per unit incident flux
    real(r8) :: improad_r_sky_dir(bounds%begl:bounds%endl)       ! improad_r_dir to sky per unit incident flux
    real(r8) :: improad_r_sunwall_dir(bounds%begl:bounds%endl)   ! improad_r_dir to sunlit wall per unit incident flux
    real(r8) :: improad_r_shadewall_dir(bounds%begl:bounds%endl) ! improad_r_dir to shaded wall per unit incident flux
    real(r8) :: improad_r_sky_dif(bounds%begl:bounds%endl)       ! improad_r_dif to sky per unit incident flux
    real(r8) :: improad_r_sunwall_dif(bounds%begl:bounds%endl)   ! improad_r_dif to sunlit wall per unit incident flux
    real(r8) :: improad_r_shadewall_dif(bounds%begl:bounds%endl) ! improad_r_dif to shaded wall per unit incident flux

    real(r8) :: perroad_a_dir(bounds%begl:bounds%endl)           ! absorbed direct solar for pervious road after "n" reflections per unit incident flux
    real(r8) :: perroad_a_dif(bounds%begl:bounds%endl)           ! absorbed diffuse solar for pervious road after "n" reflections per unit incident flux
    real(r8) :: perroad_r_dir(bounds%begl:bounds%endl)           ! reflected direct solar for pervious road after "n" reflections per unit incident flux
    real(r8) :: perroad_r_dif(bounds%begl:bounds%endl)           ! reflected diffuse solar for pervious road after "n" reflections per unit incident flux
    real(r8) :: perroad_r_sky_dir(bounds%begl:bounds%endl)       ! perroad_r_dir to sky per unit incident flux 
    real(r8) :: perroad_r_sunwall_dir(bounds%begl:bounds%endl)   ! perroad_r_dir to sunlit wall per unit incident flux
    real(r8) :: perroad_r_shadewall_dir(bounds%begl:bounds%endl) ! perroad_r_dir to shaded wall per unit incident flux
    real(r8) :: perroad_r_sky_dif(bounds%begl:bounds%endl)       ! perroad_r_dif to sky per unit incident flux
    real(r8) :: perroad_r_sunwall_dif(bounds%begl:bounds%endl)   ! perroad_r_dif to sunlit wall per unit incident flux
    real(r8) :: perroad_r_shadewall_dif(bounds%begl:bounds%endl) ! perroad_r_dif to shaded wall per unit incident flux

    real(r8) :: road_a_dir(bounds%begl:bounds%endl)              ! absorbed direct solar for total road after "n" reflections per unit incident flux
    real(r8) :: road_a_dif(bounds%begl:bounds%endl)              ! absorbed diffuse solar for total road after "n" reflections per unit incident flux
    real(r8) :: road_r_dir(bounds%begl:bounds%endl)              ! reflected direct solar for total road after "n" reflections per unit incident flux
    real(r8) :: road_r_dif(bounds%begl:bounds%endl)              ! reflected diffuse solar for total road after "n" reflections per unit incident flux
    real(r8) :: road_r_sky_dir(bounds%begl:bounds%endl)          ! road_r_dir to sky per unit incident flux
    real(r8) :: road_r_sunwall_dir(bounds%begl:bounds%endl)      ! road_r_dir to sunlit wall per unit incident flux
    real(r8) :: road_r_shadewall_dir(bounds%begl:bounds%endl)    ! road_r_dir to shaded wall per unit incident flux
    real(r8) :: road_r_sky_dif(bounds%begl:bounds%endl)          ! road_r_dif to sky per unit incident flux
    real(r8) :: road_r_sunwall_dif(bounds%begl:bounds%endl)      ! road_r_dif to sunlit wall per unit incident flux
    real(r8) :: road_r_shadewall_dif(bounds%begl:bounds%endl)    ! road_r_dif to shaded wall per unit incident flux

    real(r8) :: sunwall_a_dir(bounds%begl:bounds%endl)           ! absorbed direct solar for sunlit wall (per unit wall area) after "n" reflections per unit incident flux
    real(r8) :: sunwall_a_dif(bounds%begl:bounds%endl)           ! absorbed diffuse solar for sunlit wall (per unit wall area) after "n" reflections per unit incident flux
    real(r8) :: sunwall_r_dir(bounds%begl:bounds%endl)           ! reflected direct solar for sunlit wall (per unit wall area) after "n" reflections per unit incident flux
    real(r8) :: sunwall_r_dif(bounds%begl:bounds%endl)           ! reflected diffuse solar for sunlit wall (per unit wall area) after "n" reflections per unit incident flux
    real(r8) :: sunwall_r_sky_dir(bounds%begl:bounds%endl)       ! sunwall_r_dir to sky per unit incident flux
    real(r8) :: sunwall_r_road_dir(bounds%begl:bounds%endl)      ! sunwall_r_dir to road per unit incident flux
    real(r8) :: sunwall_r_shadewall_dir(bounds%begl:bounds%endl) ! sunwall_r_dir to opposing (shaded) wall per unit incident flux
    real(r8) :: sunwall_r_sky_dif(bounds%begl:bounds%endl)       ! sunwall_r_dif to sky per unit incident flux
    real(r8) :: sunwall_r_road_dif(bounds%begl:bounds%endl)      ! sunwall_r_dif to road per unit incident flux
    real(r8) :: sunwall_r_shadewall_dif(bounds%begl:bounds%endl) ! sunwall_r_dif to opposing (shaded) wall per unit incident flux

    real(r8) :: shadewall_a_dir(bounds%begl:bounds%endl)         ! absorbed direct solar for shaded wall (per unit wall area) after "n" reflections per unit incident flux
    real(r8) :: shadewall_a_dif(bounds%begl:bounds%endl)         ! absorbed diffuse solar for shaded wall (per unit wall area) after "n" reflections per unit incident flux
    real(r8) :: shadewall_r_dir(bounds%begl:bounds%endl)         ! reflected direct solar for shaded wall (per unit wall area) after "n" reflections per unit incident flux
    real(r8) :: shadewall_r_dif(bounds%begl:bounds%endl)         ! reflected diffuse solar for shaded wall (per unit wall area) after "n" reflections per unit incident flux
    real(r8) :: shadewall_r_sky_dir(bounds%begl:bounds%endl)     ! shadewall_r_dir to sky per unit incident flux
    real(r8) :: shadewall_r_road_dir(bounds%begl:bounds%endl)    ! shadewall_r_dir to road per unit incident flux
    real(r8) :: shadewall_r_sunwall_dir(bounds%begl:bounds%endl) ! shadewall_r_dir to opposing (sunlit) wall per unit incident flux
    real(r8) :: shadewall_r_sky_dif(bounds%begl:bounds%endl)     ! shadewall_r_dif to sky per unit incident flux
    real(r8) :: shadewall_r_road_dif(bounds%begl:bounds%endl)    ! shadewall_r_dif to road per unit incident flux
    real(r8) :: shadewall_r_sunwall_dif(bounds%begl:bounds%endl) ! shadewall_r_dif to opposing (sunlit) wall per unit incident flux

    real(r8) :: canyon_alb_dir(bounds%begl:bounds%endl)          ! direct canyon albedo
    real(r8) :: canyon_alb_dif(bounds%begl:bounds%endl)          ! diffuse canyon albedo

    real(r8) :: stot(bounds%begl:bounds%endl)                    ! sum of radiative terms
    real(r8) :: stot_dir(bounds%begl:bounds%endl)                ! sum of direct radiative terms
    real(r8) :: stot_dif(bounds%begl:bounds%endl)                ! sum of diffuse radiative terms

    integer  :: l,fl,ib                          ! indices
    integer  :: iter_dir,iter_dif                ! iteration counter
    real(r8) :: crit                             ! convergence criterion
    real(r8) :: err                              ! energy conservation error
    integer  :: pass
    integer, parameter :: n = 50                 ! number of interations
    real(r8) :: sabs_road                        ! temporary for absorption over road
    real(r8) :: sref_road                        ! temporary for reflected over road
    real(r8), parameter :: errcrit  = .00001_r8  ! error criteria
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(coszen)             == (/bounds%endl/)),         errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(canyon_hwr)         == (/bounds%endl/)),         errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(wtroad_perv)        == (/bounds%endl/)),         errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(sdir)               == (/bounds%endl, numrad/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(sdif)               == (/bounds%endl, numrad/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(alb_improad_dir)    == (/bounds%endl, numrad/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(alb_perroad_dir)    == (/bounds%endl, numrad/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(alb_wall_dir)       == (/bounds%endl, numrad/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(alb_roof_dir)       == (/bounds%endl, numrad/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(alb_improad_dif)    == (/bounds%endl, numrad/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(alb_perroad_dif)    == (/bounds%endl, numrad/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(alb_wall_dif)       == (/bounds%endl, numrad/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(alb_roof_dif)       == (/bounds%endl, numrad/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(sdir_road)          == (/bounds%endl, numrad/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(sdir_sunwall)       == (/bounds%endl, numrad/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(sdir_shadewall)     == (/bounds%endl, numrad/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(sdif_road)          == (/bounds%endl, numrad/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(sdif_sunwall)       == (/bounds%endl, numrad/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(sdif_shadewall)     == (/bounds%endl, numrad/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(sref_improad_dir)   == (/bounds%endl, numrad/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(sref_perroad_dir)   == (/bounds%endl, numrad/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(sref_improad_dif)   == (/bounds%endl, numrad/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(sref_perroad_dif)   == (/bounds%endl, numrad/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(sref_sunwall_dir)   == (/bounds%endl, numrad/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(sref_sunwall_dif)   == (/bounds%endl, numrad/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(sref_shadewall_dir) == (/bounds%endl, numrad/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(sref_shadewall_dif) == (/bounds%endl, numrad/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(sref_roof_dir)      == (/bounds%endl, numrad/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(sref_roof_dif)      == (/bounds%endl, numrad/)), errMsg(__FILE__, __LINE__))

    associate(                                                           & 
         vf_sr              =>    urbanparams_vars%vf_sr               , & ! Input:  [real(r8) (:)   ]  view factor of sky for road                       
         vf_wr              =>    urbanparams_vars%vf_wr               , & ! Input:  [real(r8) (:)   ]  view factor of one wall for road                  
         vf_sw              =>    urbanparams_vars%vf_sw               , & ! Input:  [real(r8) (:)   ]  view factor of sky for one wall                   
         vf_rw              =>    urbanparams_vars%vf_rw               , & ! Input:  [real(r8) (:)   ]  view factor of road for one wall                  
         vf_ww              =>    urbanparams_vars%vf_ww               , & ! Input:  [real(r8) (:)   ]  view factor of opposing wall for one wall         

         sabs_roof_dir      =>    solarabs_vars%sabs_roof_dir_lun      , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by roof per unit ground area per unit incident flux
         sabs_roof_dif      =>    solarabs_vars%sabs_roof_dif_lun      , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by roof per unit ground area per unit incident flux
         sabs_sunwall_dir   =>    solarabs_vars%sabs_sunwall_dir_lun   , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by sunwall per unit wall area per unit incident flux
         sabs_sunwall_dif   =>    solarabs_vars%sabs_sunwall_dif_lun   , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by sunwall per unit wall area per unit incident flux
         sabs_shadewall_dir =>    solarabs_vars%sabs_shadewall_dir_lun , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by shadewall per unit wall area per unit incident flux
         sabs_shadewall_dif =>    solarabs_vars%sabs_shadewall_dif_lun , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by shadewall per unit wall area per unit incident flux
         sabs_improad_dir   =>    solarabs_vars%sabs_improad_dir_lun   , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by impervious road per unit ground area per unit incident flux
         sabs_improad_dif   =>    solarabs_vars%sabs_improad_dif_lun   , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by impervious road per unit ground area per unit incident flux
         sabs_perroad_dir   =>    solarabs_vars%sabs_perroad_dir_lun   , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by pervious road per unit ground area per unit incident flux
         sabs_perroad_dif   =>    solarabs_vars%sabs_perroad_dif_lun     & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by pervious road per unit ground area per unit incident flux
         )

      ! Calculate impervious road
      
      do fl = 1,num_urbanl 
         l = filter_urbanl(fl)
         wtroad_imperv(l) = 1._r8 - wtroad_perv(l)
      end do

      do ib = 1,numrad
         do fl = 1,num_urbanl
            l = filter_urbanl(fl)
            if (coszen(l) > 0._r8) then

               ! initial absorption and reflection for road and both walls. 
               ! distribute reflected radiation to sky, road, and walls 
               ! according to appropriate view factor. radiation reflected to 
               ! road and walls will undergo multiple reflections within the canyon. 
               ! do separately for direct beam and diffuse radiation.

               ! direct beam

               road_a_dir(l)              = 0.0_r8
               road_r_dir(l)              = 0.0_r8
               improad_a_dir(l)           = (1._r8-alb_improad_dir(l,ib)) * sdir_road(l,ib) 
               improad_r_dir(l)           =     alb_improad_dir(l,ib)  * sdir_road(l,ib) 
               improad_r_sky_dir(l)       = improad_r_dir(l) * vf_sr(l)
               improad_r_sunwall_dir(l)   = improad_r_dir(l) * vf_wr(l)
               improad_r_shadewall_dir(l) = improad_r_dir(l) * vf_wr(l)
               road_a_dir(l)              = road_a_dir(l) + improad_a_dir(l)*wtroad_imperv(l)
               road_r_dir(l)              = road_r_dir(l) + improad_r_dir(l)*wtroad_imperv(l)

               perroad_a_dir(l)           = (1._r8-alb_perroad_dir(l,ib)) * sdir_road(l,ib) 
               perroad_r_dir(l)           =     alb_perroad_dir(l,ib)  * sdir_road(l,ib) 
               perroad_r_sky_dir(l)       = perroad_r_dir(l) * vf_sr(l)
               perroad_r_sunwall_dir(l)   = perroad_r_dir(l) * vf_wr(l)
               perroad_r_shadewall_dir(l) = perroad_r_dir(l) * vf_wr(l)
               road_a_dir(l)              = road_a_dir(l) + perroad_a_dir(l)*wtroad_perv(l)
               road_r_dir(l)              = road_r_dir(l) + perroad_r_dir(l)*wtroad_perv(l)

               road_r_sky_dir(l)          = road_r_dir(l) * vf_sr(l)
               road_r_sunwall_dir(l)      = road_r_dir(l) * vf_wr(l)
               road_r_shadewall_dir(l)    = road_r_dir(l) * vf_wr(l)

               sunwall_a_dir(l)           = (1._r8-alb_wall_dir(l,ib)) * sdir_sunwall(l,ib)
               sunwall_r_dir(l)           =     alb_wall_dir(l,ib)  * sdir_sunwall(l,ib)
               sunwall_r_sky_dir(l)       = sunwall_r_dir(l) * vf_sw(l)
               sunwall_r_road_dir(l)      = sunwall_r_dir(l) * vf_rw(l)
               sunwall_r_shadewall_dir(l) = sunwall_r_dir(l) * vf_ww(l)

               shadewall_a_dir(l)         = (1._r8-alb_wall_dir(l,ib)) * sdir_shadewall(l,ib)
               shadewall_r_dir(l)         =     alb_wall_dir(l,ib)  * sdir_shadewall(l,ib)
               shadewall_r_sky_dir(l)     = shadewall_r_dir(l) * vf_sw(l)
               shadewall_r_road_dir(l)    = shadewall_r_dir(l) * vf_rw(l)
               shadewall_r_sunwall_dir(l) = shadewall_r_dir(l) * vf_ww(l)

               ! diffuse

               road_a_dif(l)              = 0.0_r8
               road_r_dif(l)              = 0.0_r8
               improad_a_dif(l)           = (1._r8-alb_improad_dif(l,ib)) * sdif_road(l,ib) 
               improad_r_dif(l)           =     alb_improad_dif(l,ib)  * sdif_road(l,ib) 
               improad_r_sky_dif(l)       = improad_r_dif(l) * vf_sr(l)
               improad_r_sunwall_dif(l)   = improad_r_dif(l) * vf_wr(l)
               improad_r_shadewall_dif(l) = improad_r_dif(l) * vf_wr(l)
               road_a_dif(l)              = road_a_dif(l) + improad_a_dif(l)*wtroad_imperv(l)
               road_r_dif(l)              = road_r_dif(l) + improad_r_dif(l)*wtroad_imperv(l)

               perroad_a_dif(l)           = (1._r8-alb_perroad_dif(l,ib)) * sdif_road(l,ib) 
               perroad_r_dif(l)           =     alb_perroad_dif(l,ib)  * sdif_road(l,ib) 
               perroad_r_sky_dif(l)       = perroad_r_dif(l) * vf_sr(l)
               perroad_r_sunwall_dif(l)   = perroad_r_dif(l) * vf_wr(l)
               perroad_r_shadewall_dif(l) = perroad_r_dif(l) * vf_wr(l)
               road_a_dif(l)              = road_a_dif(l) + perroad_a_dif(l)*wtroad_perv(l)
               road_r_dif(l)              = road_r_dif(l) + perroad_r_dif(l)*wtroad_perv(l)

               road_r_sky_dif(l)          = road_r_dif(l) * vf_sr(l)
               road_r_sunwall_dif(l)      = road_r_dif(l) * vf_wr(l)
               road_r_shadewall_dif(l)    = road_r_dif(l) * vf_wr(l)

               sunwall_a_dif(l)           = (1._r8-alb_wall_dif(l,ib)) * sdif_sunwall(l,ib)
               sunwall_r_dif(l)           =     alb_wall_dif(l,ib)  * sdif_sunwall(l,ib)
               sunwall_r_sky_dif(l)       = sunwall_r_dif(l) * vf_sw(l)
               sunwall_r_road_dif(l)      = sunwall_r_dif(l) * vf_rw(l)
               sunwall_r_shadewall_dif(l) = sunwall_r_dif(l) * vf_ww(l)

               shadewall_a_dif(l)         = (1._r8-alb_wall_dif(l,ib)) * sdif_shadewall(l,ib)
               shadewall_r_dif(l)         =     alb_wall_dif(l,ib)  * sdif_shadewall(l,ib)
               shadewall_r_sky_dif(l)     = shadewall_r_dif(l) * vf_sw(l)
               shadewall_r_road_dif(l)    = shadewall_r_dif(l) * vf_rw(l) 
               shadewall_r_sunwall_dif(l) = shadewall_r_dif(l) * vf_ww(l) 

               ! initialize sum of direct and diffuse solar absorption and reflection for road and both walls

               sabs_improad_dir(l,ib)   = improad_a_dir(l)
               sabs_perroad_dir(l,ib)   = perroad_a_dir(l)
               sabs_sunwall_dir(l,ib)   = sunwall_a_dir(l)
               sabs_shadewall_dir(l,ib) = shadewall_a_dir(l)

               sabs_improad_dif(l,ib)   = improad_a_dif(l)
               sabs_perroad_dif(l,ib)   = perroad_a_dif(l)
               sabs_sunwall_dif(l,ib)   = sunwall_a_dif(l)
               sabs_shadewall_dif(l,ib) = shadewall_a_dif(l)

               sref_improad_dir(l,ib)   = improad_r_sky_dir(l) 
               sref_perroad_dir(l,ib)   = perroad_r_sky_dir(l) 
               sref_sunwall_dir(l,ib)   = sunwall_r_sky_dir(l) 
               sref_shadewall_dir(l,ib) = shadewall_r_sky_dir(l) 

               sref_improad_dif(l,ib)   = improad_r_sky_dif(l)
               sref_perroad_dif(l,ib)   = perroad_r_sky_dif(l)
               sref_sunwall_dif(l,ib)   = sunwall_r_sky_dif(l)
               sref_shadewall_dif(l,ib) = shadewall_r_sky_dif(l)
            endif

         end do

         ! absorption and reflection for walls and road with multiple reflections
         ! (i.e., absorb and reflect initial reflection in canyon and allow for 
         ! subsequent scattering)
         !
         ! (1) absorption and reflection of scattered solar radiation
         !     road: reflected fluxes from walls need to be projected to ground area
         !     wall: reflected flux from road needs to be projected to wall area
         !
         ! (2) add absorbed radiation for ith reflection to total absorbed
         !
         ! (3) distribute reflected radiation to sky, road, and walls according to view factors
         !
         ! (4) add solar reflection to sky for ith reflection to total reflection
         !
         ! (5) stop iteration when absorption for ith reflection is less than some nominal amount. 
         !     small convergence criteria is required to ensure solar radiation is conserved
         !
         ! do separately for direct beam and diffuse

         do fl = 1,num_urbanl
            l = filter_urbanl(fl)
            if (coszen(l) > 0._r8) then

               ! reflected direct beam

               do iter_dir = 1, n
                  ! step (1)

                  stot(l) = (sunwall_r_road_dir(l) + shadewall_r_road_dir(l))*canyon_hwr(l)

                  road_a_dir(l) = 0.0_r8
                  road_r_dir(l) = 0.0_r8
                  improad_a_dir(l) = (1._r8-alb_improad_dir(l,ib)) * stot(l) 
                  improad_r_dir(l) =     alb_improad_dir(l,ib)  * stot(l) 
                  road_a_dir(l)    = road_a_dir(l) + improad_a_dir(l)*wtroad_imperv(l)
                  road_r_dir(l)    = road_r_dir(l) + improad_r_dir(l)*wtroad_imperv(l)
                  perroad_a_dir(l) = (1._r8-alb_perroad_dir(l,ib)) * stot(l) 
                  perroad_r_dir(l) =     alb_perroad_dir(l,ib)  * stot(l) 
                  road_a_dir(l)    = road_a_dir(l) + perroad_a_dir(l)*wtroad_perv(l)
                  road_r_dir(l)    = road_r_dir(l) + perroad_r_dir(l)*wtroad_perv(l)

                  stot(l) = road_r_sunwall_dir(l)/canyon_hwr(l) + shadewall_r_sunwall_dir(l)
                  sunwall_a_dir(l) = (1._r8-alb_wall_dir(l,ib)) * stot(l)
                  sunwall_r_dir(l) =     alb_wall_dir(l,ib)  * stot(l)

                  stot(l) = road_r_shadewall_dir(l)/canyon_hwr(l) + sunwall_r_shadewall_dir(l)
                  shadewall_a_dir(l) = (1._r8-alb_wall_dir(l,ib)) * stot(l)
                  shadewall_r_dir(l) =     alb_wall_dir(l,ib)  * stot(l)

                  ! step (2)

                  sabs_improad_dir(l,ib)   = sabs_improad_dir(l,ib)   + improad_a_dir(l)
                  sabs_perroad_dir(l,ib)   = sabs_perroad_dir(l,ib)   + perroad_a_dir(l)
                  sabs_sunwall_dir(l,ib)   = sabs_sunwall_dir(l,ib)   + sunwall_a_dir(l)
                  sabs_shadewall_dir(l,ib) = sabs_shadewall_dir(l,ib) + shadewall_a_dir(l)

                  ! step (3)

                  improad_r_sky_dir(l)       = improad_r_dir(l) * vf_sr(l)
                  improad_r_sunwall_dir(l)   = improad_r_dir(l) * vf_wr(l)
                  improad_r_shadewall_dir(l) = improad_r_dir(l) * vf_wr(l)

                  perroad_r_sky_dir(l)       = perroad_r_dir(l) * vf_sr(l)
                  perroad_r_sunwall_dir(l)   = perroad_r_dir(l) * vf_wr(l)
                  perroad_r_shadewall_dir(l) = perroad_r_dir(l) * vf_wr(l)

                  road_r_sky_dir(l)          = road_r_dir(l) * vf_sr(l)
                  road_r_sunwall_dir(l)      = road_r_dir(l) * vf_wr(l)
                  road_r_shadewall_dir(l)    = road_r_dir(l) * vf_wr(l)

                  sunwall_r_sky_dir(l)       = sunwall_r_dir(l) * vf_sw(l)
                  sunwall_r_road_dir(l)      = sunwall_r_dir(l) * vf_rw(l)
                  sunwall_r_shadewall_dir(l) = sunwall_r_dir(l) * vf_ww(l)

                  shadewall_r_sky_dir(l)     = shadewall_r_dir(l) * vf_sw(l)
                  shadewall_r_road_dir(l)    = shadewall_r_dir(l) * vf_rw(l)
                  shadewall_r_sunwall_dir(l) = shadewall_r_dir(l) * vf_ww(l)

                  ! step (4)

                  sref_improad_dir(l,ib)   = sref_improad_dir(l,ib) + improad_r_sky_dir(l)
                  sref_perroad_dir(l,ib)   = sref_perroad_dir(l,ib) + perroad_r_sky_dir(l)
                  sref_sunwall_dir(l,ib)   = sref_sunwall_dir(l,ib) + sunwall_r_sky_dir(l)
                  sref_shadewall_dir(l,ib) = sref_shadewall_dir(l,ib) + shadewall_r_sky_dir(l)

                  ! step (5)

                  crit = max(road_a_dir(l), sunwall_a_dir(l), shadewall_a_dir(l))
                  if (crit < errcrit) exit
               end do
               if (iter_dir >= n) then
                  write (iulog,*) 'urban net solar radiation error: no convergence, direct beam'
                  write (iulog,*) 'clm model is stopping'
                  call endrun(decomp_index=l, clmlevel=namel, msg=errmsg(__FILE__, __LINE__))
               endif

               ! reflected diffuse

               do iter_dif = 1, n
                  ! step (1)

                  stot(l) = (sunwall_r_road_dif(l) + shadewall_r_road_dif(l))*canyon_hwr(l)
                  road_a_dif(l)    = 0.0_r8
                  road_r_dif(l)    = 0.0_r8
                  improad_a_dif(l) = (1._r8-alb_improad_dif(l,ib)) * stot(l) 
                  improad_r_dif(l) =     alb_improad_dif(l,ib)  * stot(l) 
                  road_a_dif(l)    = road_a_dif(l) + improad_a_dif(l)*wtroad_imperv(l)
                  road_r_dif(l)    = road_r_dif(l) + improad_r_dif(l)*wtroad_imperv(l)
                  perroad_a_dif(l) = (1._r8-alb_perroad_dif(l,ib)) * stot(l) 
                  perroad_r_dif(l) =     alb_perroad_dif(l,ib)  * stot(l) 
                  road_a_dif(l)    = road_a_dif(l) + perroad_a_dif(l)*wtroad_perv(l)
                  road_r_dif(l)    = road_r_dif(l) + perroad_r_dif(l)*wtroad_perv(l)

                  stot(l) = road_r_sunwall_dif(l)/canyon_hwr(l) + shadewall_r_sunwall_dif(l)
                  sunwall_a_dif(l) = (1._r8-alb_wall_dif(l,ib)) * stot(l)
                  sunwall_r_dif(l) =     alb_wall_dif(l,ib)  * stot(l)

                  stot(l) = road_r_shadewall_dif(l)/canyon_hwr(l) + sunwall_r_shadewall_dif(l)
                  shadewall_a_dif(l) = (1._r8-alb_wall_dif(l,ib)) * stot(l)
                  shadewall_r_dif(l) =     alb_wall_dif(l,ib)  * stot(l)

                  ! step (2)

                  sabs_improad_dif(l,ib)   = sabs_improad_dif(l,ib)   + improad_a_dif(l)
                  sabs_perroad_dif(l,ib)   = sabs_perroad_dif(l,ib)   + perroad_a_dif(l)
                  sabs_sunwall_dif(l,ib)   = sabs_sunwall_dif(l,ib)   + sunwall_a_dif(l)
                  sabs_shadewall_dif(l,ib) = sabs_shadewall_dif(l,ib) + shadewall_a_dif(l)

                  ! step (3)

                  improad_r_sky_dif(l)       = improad_r_dif(l) * vf_sr(l)
                  improad_r_sunwall_dif(l)   = improad_r_dif(l) * vf_wr(l)
                  improad_r_shadewall_dif(l) = improad_r_dif(l) * vf_wr(l)

                  perroad_r_sky_dif(l)       = perroad_r_dif(l) * vf_sr(l)
                  perroad_r_sunwall_dif(l)   = perroad_r_dif(l) * vf_wr(l)
                  perroad_r_shadewall_dif(l) = perroad_r_dif(l) * vf_wr(l)

                  road_r_sky_dif(l)          = road_r_dif(l) * vf_sr(l)
                  road_r_sunwall_dif(l)      = road_r_dif(l) * vf_wr(l)
                  road_r_shadewall_dif(l)    = road_r_dif(l) * vf_wr(l)

                  sunwall_r_sky_dif(l)       = sunwall_r_dif(l) * vf_sw(l)
                  sunwall_r_road_dif(l)      = sunwall_r_dif(l) * vf_rw(l)
                  sunwall_r_shadewall_dif(l) = sunwall_r_dif(l) * vf_ww(l)

                  shadewall_r_sky_dif(l)     = shadewall_r_dif(l) * vf_sw(l)
                  shadewall_r_road_dif(l)    = shadewall_r_dif(l) * vf_rw(l)
                  shadewall_r_sunwall_dif(l) = shadewall_r_dif(l) * vf_ww(l)

                  ! step (4)

                  sref_improad_dif(l,ib)   = sref_improad_dif(l,ib)   + improad_r_sky_dif(l)
                  sref_perroad_dif(l,ib)   = sref_perroad_dif(l,ib)   + perroad_r_sky_dif(l)
                  sref_sunwall_dif(l,ib)   = sref_sunwall_dif(l,ib)   + sunwall_r_sky_dif(l)
                  sref_shadewall_dif(l,ib) = sref_shadewall_dif(l,ib) + shadewall_r_sky_dif(l)

                  ! step (5)

                  crit = max(road_a_dif(l), sunwall_a_dif(l), shadewall_a_dif(l))
                  if (crit < errcrit) exit
               end do
               if (iter_dif >= n) then
                  write (iulog,*) 'urban net solar radiation error: no convergence, diffuse'
                  write (iulog,*) 'clm model is stopping'
                  call endrun(decomp_index=l, clmlevel=namel, msg=errmsg(__FILE__, __LINE__))
               endif

               ! total reflected by canyon - sum of solar reflection to sky from canyon.
               ! project wall fluxes to horizontal surface

               sref_canyon_dir(l) = 0.0_r8
               sref_canyon_dif(l) = 0.0_r8
               sref_canyon_dir(l) = sref_canyon_dir(l) + sref_improad_dir(l,ib)*wtroad_imperv(l)
               sref_canyon_dif(l) = sref_canyon_dif(l) + sref_improad_dif(l,ib)*wtroad_imperv(l)
               sref_canyon_dir(l) = sref_canyon_dir(l) + sref_perroad_dir(l,ib)*wtroad_perv(l)
               sref_canyon_dif(l) = sref_canyon_dif(l) + sref_perroad_dif(l,ib)*wtroad_perv(l)
               sref_canyon_dir(l) = sref_canyon_dir(l) + (sref_sunwall_dir(l,ib) + sref_shadewall_dir(l,ib))*canyon_hwr(l)
               sref_canyon_dif(l) = sref_canyon_dif(l) + (sref_sunwall_dif(l,ib) + sref_shadewall_dif(l,ib))*canyon_hwr(l)

               ! total absorbed by canyon. project wall fluxes to horizontal surface

               sabs_canyon_dir(l) = 0.0_r8
               sabs_canyon_dif(l) = 0.0_r8
               sabs_canyon_dir(l) = sabs_canyon_dir(l) + sabs_improad_dir(l,ib)*wtroad_imperv(l)
               sabs_canyon_dif(l) = sabs_canyon_dif(l) + sabs_improad_dif(l,ib)*wtroad_imperv(l)
               sabs_canyon_dir(l) = sabs_canyon_dir(l) + sabs_perroad_dir(l,ib)*wtroad_perv(l)
               sabs_canyon_dif(l) = sabs_canyon_dif(l) + sabs_perroad_dif(l,ib)*wtroad_perv(l)
               sabs_canyon_dir(l) = sabs_canyon_dir(l) + (sabs_sunwall_dir(l,ib) + sabs_shadewall_dir(l,ib))*canyon_hwr(l)
               sabs_canyon_dif(l) = sabs_canyon_dif(l) + (sabs_sunwall_dif(l,ib) + sabs_shadewall_dif(l,ib))*canyon_hwr(l)

               ! conservation check. note: previous conservation checks confirm partioning of total direct
               ! beam and diffuse radiation from atmosphere to road and walls is conserved as
               !    sdir (from atmosphere) = sdir_road + (sdir_sunwall + sdir_shadewall)*canyon_hwr
               !    sdif (from atmosphere) = sdif_road + (sdif_sunwall + sdif_shadewall)*canyon_hwr

               stot_dir(l) = sdir_road(l,ib) + (sdir_sunwall(l,ib) + sdir_shadewall(l,ib))*canyon_hwr(l)
               stot_dif(l) = sdif_road(l,ib) + (sdif_sunwall(l,ib) + sdif_shadewall(l,ib))*canyon_hwr(l)

               err = stot_dir(l) + stot_dif(l) &
                    - (sabs_canyon_dir(l) + sabs_canyon_dif(l) + sref_canyon_dir(l) + sref_canyon_dif(l))
               if (abs(err) > 0.001_r8 ) then
                  write(iulog,*)'urban net solar radiation balance error for ib=',ib,' err= ',err
                  write(iulog,*)' l= ',l,' ib= ',ib 
                  write(iulog,*)' stot_dir        = ',stot_dir(l)
                  write(iulog,*)' stot_dif        = ',stot_dif(l)
                  write(iulog,*)' sabs_canyon_dir = ',sabs_canyon_dir(l)
                  write(iulog,*)' sabs_canyon_dif = ',sabs_canyon_dif(l)
                  write(iulog,*)' sref_canyon_dir = ',sref_canyon_dir(l)
                  write(iulog,*)' sref_canyon_dif = ',sref_canyon_dir(l)
                  write(iulog,*) 'clm model is stopping'
                  call endrun(decomp_index=l, clmlevel=namel, msg=errmsg(__FILE__, __LINE__))
               endif

               ! canyon albedo

               canyon_alb_dif(l) = sref_canyon_dif(l) / max(stot_dif(l), 1.e-06_r8)
               canyon_alb_dir(l) = sref_canyon_dir(l) / max(stot_dir(l), 1.e-06_r8)
            end if

         end do   ! end of landunit loop

         ! Refected and absorbed solar radiation per unit incident radiation for roof

         do fl = 1,num_urbanl
            l = filter_urbanl(fl)
            if (coszen(l) > 0._r8) then
               sref_roof_dir(l,ib) = alb_roof_dir(l,ib) * sdir(l,ib)
               sref_roof_dif(l,ib) = alb_roof_dif(l,ib) * sdif(l,ib)
               sabs_roof_dir(l,ib) = sdir(l,ib) - sref_roof_dir(l,ib)
               sabs_roof_dif(l,ib) = sdif(l,ib) - sref_roof_dif(l,ib)
            end if
         end do

      end do   ! end of radiation band loop

    end associate

  end subroutine net_solar

end module UrbanAlbedoMod

