module UrbanMod

#include "shr_assert.h"

  !----------------------------------------------------------------------- 
  ! !DESCRIPTION: 
  ! Calculate solar and longwave radiation, and turbulent fluxes for urban landunit
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use decompMod      , only : bounds_type
  use clm_varpar     , only : numrad
  use clm_varcon     , only : isecspday, degpsec
  use clm_varctl     , only : iulog
  use abortutils     , only : endrun  
  use shr_sys_mod    , only : shr_sys_flush 
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: UrbanAlbedo       ! Urban albedos  
  public :: UrbanSnowAlbedo   ! Urban snow albedos
  public :: UrbanParamInit    ! Initialization of urban parameter data structure
  public :: UrbanRadiation    ! Urban radiative fluxes
  public :: UrbanFluxes       ! Urban turbulent fluxes

  ! !Urban control variables
  character(len= *), parameter, public :: urban_hac_off = 'OFF'               ! 
  character(len= *), parameter, public :: urban_hac_on =  'ON'                ! 
  character(len= *), parameter, public :: urban_wasteheat_on = 'ON_WASTEHEAT' ! 
  character(len= 16), public :: urban_hac = urban_hac_off
  logical, public :: urban_traffic = .false.        ! urban traffic fluxes
  !
  ! PRIVATE MEMBER FUNCTIONS
  private :: incident_direct  ! Direct beam solar rad incident on walls and road in urban canyon 
  private :: incident_diffuse ! Diffuse solar rad incident on walls and road in urban canyon
  private :: net_solar        ! Solar radiation absorbed by road and both walls in urban canyon 
  private :: net_longwave     ! Net longwave radiation for road and both walls in urban canyon 
  !
  ! PRIVATE TYPES
  private
  type urban_params_t
     real(r8), allocatable :: wind_hgt_canyon(:)       ! height above road at which wind in canyon is to be computed (m)
     real(r8), allocatable :: em_roof(:)               ! roof emissivity
     real(r8), allocatable :: em_improad(:)            ! impervious road emissivity
     real(r8), allocatable :: em_perroad(:)            ! pervious road emissivity
     real(r8), allocatable :: em_wall(:)               ! wall emissivity
     real(r8), allocatable :: alb_roof_dir(:,:)        ! direct  roof albedo
     real(r8), allocatable :: alb_roof_dif(:,:)        ! diffuse roof albedo
     real(r8), allocatable :: alb_improad_dir(:,:)     ! direct  impervious road albedo
     real(r8), allocatable :: alb_improad_dif(:,:)     ! diffuse impervious road albedo
     real(r8), allocatable :: alb_perroad_dir(:,:)     ! direct  pervious road albedo
     real(r8), allocatable :: alb_perroad_dif(:,:)     ! diffuse pervious road albedo
     real(r8), allocatable :: alb_wall_dir(:,:)        ! direct  wall albedo
     real(r8), allocatable :: alb_wall_dif(:,:)        ! diffuse wall albedo
  end type urban_params_t

  type (urban_params_t), private :: urban_params  ! urban parameters

  integer,  private, parameter :: noonsec   = isecspday / 2 ! seconds at local noon
  !----------------------------------------------------------------------- 

contains

  !-----------------------------------------------------------------------
  subroutine UrbanAlbedo (bounds, &
                          num_urbanl, filter_urbanl, &
                          num_urbanc, filter_urbanc, &
                          num_urbanp, filter_urbanp)
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
    use clmtype
    use shr_orb_mod  , only : shr_orb_decl, shr_orb_cosz
    use clm_varcon   , only : icol_roof, icol_sunwall, icol_shadewall, icol_road_perv, icol_road_imperv, &
                              sb
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    integer , intent(in) :: num_urbanl       ! number of urban landunits in clump
    integer , intent(in) :: filter_urbanl(:) ! urban landunit filter
    integer , intent(in) :: num_urbanc       ! number of urban columns in clump
    integer , intent(in) :: filter_urbanc(:) ! urban column filter
    integer , intent(in) :: num_urbanp       ! number of urban pfts in clump
    integer , intent(in) :: filter_urbanp(:) ! urban pft filter
    !
    ! !LOCAL VARIABLES:
    real(r8), pointer :: czen(:)      ! cosine of solar zenith angle for each column
    !
    real(r8) :: coszen(bounds%begl:bounds%endl)                 ! cosine solar zenith angle
    real(r8) :: coszen_pft(bounds%begp:bounds%endp)             ! cosine solar zenith angle for next time step (pft level)
    real(r8) :: zen(bounds%begl:bounds%endl)                    ! solar zenith angle (radians)
    real(r8) :: sdir(bounds%begl:bounds%endl, numrad)           ! direct beam solar radiation on horizontal surface
    real(r8) :: sdif(bounds%begl:bounds%endl, numrad)           ! diffuse solar radiation on horizontal surface

    real(r8) :: sdir_road(bounds%begl:bounds%endl, numrad)      ! direct beam solar radiation incident on road
    real(r8) :: sdif_road(bounds%begl:bounds%endl, numrad)      ! diffuse solar radiation incident on road
    real(r8) :: sdir_sunwall(bounds%begl:bounds%endl, numrad)   ! direct beam solar radiation (per unit wall area) incident on sunlit wall per unit incident flux
    real(r8) :: sdif_sunwall(bounds%begl:bounds%endl, numrad)   ! diffuse solar radiation (per unit wall area) incident on sunlit wall per unit incident flux
    real(r8) :: sdir_shadewall(bounds%begl:bounds%endl, numrad) ! direct beam solar radiation (per unit wall area) incident on shaded wall per unit incident flux
    real(r8) :: sdif_shadewall(bounds%begl:bounds%endl, numrad) ! diffuse solar radiation (per unit wall area) incident on shaded wall per unit incident flux
    real(r8) :: albsnd_roof(bounds%begl:bounds%endl,numrad)     ! snow albedo for roof (direct)
    real(r8) :: albsni_roof(bounds%begl:bounds%endl,numrad)     ! snow albedo for roof (diffuse)
    real(r8) :: albsnd_improad(bounds%begl:bounds%endl,numrad)  ! snow albedo for impervious road (direct)
    real(r8) :: albsni_improad(bounds%begl:bounds%endl,numrad)  ! snow albedo for impervious road (diffuse)
    real(r8) :: albsnd_perroad(bounds%begl:bounds%endl,numrad)  ! snow albedo for pervious road (direct)
    real(r8) :: albsni_perroad(bounds%begl:bounds%endl,numrad)  ! snow albedo for pervious road (diffuse)
                                       
    integer  :: fl,fp,fc,g,l,p,c,ib                ! indices
    integer  :: ic                                 ! 0=unit incoming direct; 1=unit incoming diffuse
    integer  :: num_solar                          ! counter
    real(r8) :: alb_roof_dir_s(bounds%begl:bounds%endl,numrad)     ! direct roof albedo with snow effects
    real(r8) :: alb_roof_dif_s(bounds%begl:bounds%endl,numrad)     ! diffuse roof albedo with snow effects
    real(r8) :: alb_improad_dir_s(bounds%begl:bounds%endl,numrad)  ! direct impervious road albedo with snow effects
    real(r8) :: alb_perroad_dir_s(bounds%begl:bounds%endl,numrad)  ! direct pervious road albedo with snow effects
    real(r8) :: alb_improad_dif_s(bounds%begl:bounds%endl,numrad)  ! diffuse impervious road albedo with snow effects
    real(r8) :: alb_perroad_dif_s(bounds%begl:bounds%endl,numrad)  ! diffuse pervious road albedo with snow effects
    real(r8) :: sref_roof_dir(bounds%begl:bounds%endl,numrad)      ! direct  solar reflected by roof per unit ground area per unit incident flux   
    real(r8) :: sref_roof_dif(bounds%begl:bounds%endl,numrad)      ! diffuse solar reflected by roof per unit ground area per unit incident flux   
    real(r8) :: sref_sunwall_dir(bounds%begl:bounds%endl,numrad)   ! direct  solar reflected by sunwall per unit wall area per unit incident flux  
    real(r8) :: sref_sunwall_dif(bounds%begl:bounds%endl,numrad)   ! diffuse solar reflected by sunwall per unit wall area per unit incident flux  
    real(r8) :: sref_shadewall_dir(bounds%begl:bounds%endl,numrad) ! direct  solar reflected by shadewall per unit wall area per unit incident flux  
    real(r8) :: sref_shadewall_dif(bounds%begl:bounds%endl,numrad) ! diffuse solar reflected by shadewall per unit wall area per unit incident flux  
    real(r8) :: sref_improad_dir(bounds%begl:bounds%endl,numrad)   ! direct  solar reflected by impervious road per unit ground area per unit incident flux   
    real(r8) :: sref_improad_dif(bounds%begl:bounds%endl,numrad)   ! diffuse solar reflected by impervious road per unit ground area per unit incident flux   
    real(r8) :: sref_perroad_dir(bounds%begl:bounds%endl,numrad)   ! direct  solar reflected by pervious road per unit ground area per unit incident flux   
    real(r8) :: sref_perroad_dif(bounds%begl:bounds%endl,numrad)   ! diffuse solar reflected by pervious road per unit ground area per unit incident flux   
!-----------------------------------------------------------------------

    associate(&
   vf_sr              =>    lps%vf_sr                    , & ! Output: [real(r8) (:)   ]  view factor of sky for road                       
   vf_wr              =>    lps%vf_wr                    , & ! Output: [real(r8) (:)   ]  view factor of one wall for road                  
   vf_sw              =>    lps%vf_sw                    , & ! Output: [real(r8) (:)   ]  view factor of sky for one wall                   
   vf_rw              =>    lps%vf_rw                    , & ! Output: [real(r8) (:)   ]  view factor of road for one wall                  
   vf_ww              =>    lps%vf_ww                    , & ! Output: [real(r8) (:)   ]  view factor of opposing wall for one wall         
   sabs_roof_dir      =>    lps%sabs_roof_dir            , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by roof per unit ground area per unit incident flux
   sabs_roof_dif      =>    lps%sabs_roof_dif            , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by roof per unit ground area per unit incident flux
   sabs_sunwall_dir   =>    lps%sabs_sunwall_dir         , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by sunwall per unit wall area per unit incident flux
   sabs_sunwall_dif   =>    lps%sabs_sunwall_dif         , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by sunwall per unit wall area per unit incident flux
   sabs_shadewall_dir =>    lps%sabs_shadewall_dir       , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by shadewall per unit wall area per unit incident flux
   sabs_shadewall_dif =>    lps%sabs_shadewall_dif       , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by shadewall per unit wall area per unit incident flux
   sabs_improad_dir   =>    lps%sabs_improad_dir         , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by impervious road per unit ground area per unit incident flux
   sabs_improad_dif   =>    lps%sabs_improad_dif         , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by impervious road per unit ground area per unit incident flux
   sabs_perroad_dir   =>    lps%sabs_perroad_dir         , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by pervious road per unit ground area per unit incident flux
   sabs_perroad_dif   =>    lps%sabs_perroad_dif         , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by pervious road per unit ground area per unit incident flux
   albd               =>    pps%albd                     , & ! Output: [real(r8) (:,:) ]  surface albedo (direct)                         
   albi               =>    pps%albi                     , & ! Output: [real(r8) (:,:) ]  surface albedo (diffuse)                        
   fabd               =>    pps%fabd                     , & ! Output: [real(r8) (:,:) ]  flux absorbed by veg per unit direct  flux      
   fabd_sun           =>    pps%fabd_sun                 , & ! Output: [real(r8) (:,:) ]  flux absorbed by sunlit leaf per unit direct flux
   fabd_sha           =>    pps%fabd_sha                 , & ! Output: [real(r8) (:,:) ]  flux absorbed by shaded leaf per unit direct flux
   fabi               =>    pps%fabi                     , & ! Output: [real(r8) (:,:) ]  flux absorbed by veg per unit diffuse flux      
   fabi_sun           =>    pps%fabi_sun                 , & ! Output: [real(r8) (:,:) ]  flux absorbed by sunlit leaf per unit diffuse flux
   fabi_sha           =>    pps%fabi_sha                 , & ! Output: [real(r8) (:,:) ]  flux absorbed by shaded leaf per unit diffuse flux
   ftdd               =>    pps%ftdd                     , & ! Output: [real(r8) (:,:) ]  down direct  flux below veg per unit dir flx    
   ftid               =>    pps%ftid                     , & ! Output: [real(r8) (:,:) ]  down diffuse flux below veg per unit dir flx    
   ftii               =>    pps%ftii                     , & ! Output: [real(r8) (:,:) ]  down diffuse flux below veg per unit dif flx    
   fsun               =>    pps%fsun                     , & ! Output: [real(r8) (:)   ]  sunlit fraction of canopy                         
   lat                =>    grc%lat                      , & ! Input:  [real(r8) (:)   ]  latitude (radians)                                
   lon                =>    grc%lon                      , & ! Input:  [real(r8) (:)   ]  longitude (radians)                               
   lgridcell          =>   lun%gridcell                  , & ! Input:  [integer (:)    ]  gridcell of corresponding landunit                 
   coli               =>   lun%coli                      , & ! Input:  [integer (:)    ]  beginning column index for landunit                
   ctype              =>    col%itype                    , & ! Input:  [integer (:)    ]  column type                                        
   frac_sno           =>    cps%frac_sno                 , & ! Input:  [real(r8) (:)   ]  fraction of ground covered by snow (0 to 1)       
   clandunit          =>   col%landunit                  , & ! Input:  [integer (:)    ]  column's landunit                                  
   cgridcell          =>   col%gridcell                  , & ! Input:  [integer (:)    ]  gridcell of corresponding column                   
   pgridcell          =>   pft%gridcell                  , & ! Input:  [integer (:)    ]  gridcell of corresponding pft                      
   pcolumn            =>   pft%column                    , & ! Input:  [integer (:)    ]  column of corresponding pft                        
   canyon_hwr         =>    lun%canyon_hwr               , & ! Output: [real(r8) (:)   ]  ratio of building height to street width          
   wtroad_perv        =>    lun%wtroad_perv              , & ! Output: [real(r8) (:)   ]  weight of pervious road wrt total road            
   alb_roof_dir       =>    urban_params%alb_roof_dir    , & ! Output: [real(r8) (:,:) ]  direct roof albedo                              
   alb_roof_dif       =>    urban_params%alb_roof_dif    , & ! Output: [real(r8) (:,:) ]  diffuse roof albedo                             
   alb_improad_dir    =>    urban_params%alb_improad_dir , & ! Output: [real(r8) (:,:) ]  direct impervious road albedo                   
   alb_improad_dif    =>    urban_params%alb_improad_dif , & ! Output: [real(r8) (:,:) ]  diffuse imprevious road albedo                  
   alb_perroad_dir    =>    urban_params%alb_perroad_dir , & ! Output: [real(r8) (:,:) ]  direct pervious road albedo                     
   alb_perroad_dif    =>    urban_params%alb_perroad_dif , & ! Output: [real(r8) (:,:) ]  diffuse pervious road albedo                    
   alb_wall_dir       =>    urban_params%alb_wall_dir    , & ! Output: [real(r8) (:,:) ]  direct wall albedo                              
   alb_wall_dif       =>    urban_params%alb_wall_dif    , & ! Output: [real(r8) (:,:) ]  diffuse wall albedo                             
   albgrd             =>    cps%albgrd                   , & ! Output: [real(r8) (:,:) ]  ground albedo (direct)                          
   albgri             =>    cps%albgri                   , & ! Output: [real(r8) (:,:) ]  ground albedo (diffuse)                         
   begl               =>    bounds%begl                  , &
   endl               =>    bounds%endl                    &
   )

    !TODO: this must be kept as a pointer - putting it as an assoicate completely breaks things
    czen => cps%coszen
    !TODO

    ! ----------------------------------------------------------------------------
    ! Solar declination and cosine solar zenith angle and zenith angle for 
    ! next time step
    ! ----------------------------------------------------------------------------

    do fl = 1,num_urbanl
       l = filter_urbanl(fl)
       g = lgridcell(l)
       coszen(l) = czen(coli(l))               ! Assumes coszen for each column are the same
       zen(l)    = acos(coszen(l))
    end do

    do fp = 1,num_urbanp
       p = filter_urbanp(fp)
       g = pgridcell(p)
       c = pcolumn(p)
       coszen_pft(p) = czen(c)
    end do

    ! ----------------------------------------------------------------------------
    ! Initialize clmtype output since solar radiation is only done if coszen > 0
    ! ----------------------------------------------------------------------------
    
    do ib = 1,numrad
       do fc = 1,num_urbanc
          c = filter_urbanc(fc)

          albgrd(c,ib) = 0._r8
          albgri(c,ib) = 0._r8
       end do

       do fp = 1,num_urbanp
          p = filter_urbanp(fp)
          g = pgridcell(p)
          albd(p,ib)    = 1._r8
          albi(p,ib)    = 1._r8
          fabd(p,ib)    = 0._r8
          fabd_sun(p,ib) = 0._r8
          fabd_sha(p,ib) = 0._r8
          fabi(p,ib)    = 0._r8
          fabi_sun(p,ib) = 0._r8
          fabi_sha(p,ib) = 0._r8
          if (coszen_pft(p) > 0._r8) then
             ftdd(p,ib) = 1._r8
          else
             ftdd(p,ib) = 0._r8
          end if
          ftid(p,ib)    = 0._r8
          if (coszen_pft(p) > 0._r8) then
             ftii(p,ib) = 1._r8
          else
             ftii(p,ib) = 0._r8
          end if
          if (ib == 1) then
             fsun(p)    = 0._r8
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
          
       if (num_urbanl .gt. 0) then
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
       
       if (num_urbanl .gt. 0) then
          call incident_diffuse (bounds, &
               num_urbanl, filter_urbanl, &
               canyon_hwr(begl:endl), &
               sdif(begl:endl, :), &
               sdif_road(begl:endl, :), &
               sdif_sunwall(begl:endl, :), &
               sdif_shadewall(begl:endl, :))
       end if

       ! Get snow albedos for roof and impervious and pervious road
       if (num_urbanl .gt. 0) then
          ic = 0
          call UrbanSnowAlbedo(bounds, &
               num_urbanc, filter_urbanc, &
               coszen(begl:endl), &
               ic, &
               albsnd_roof(begl:endl, :), &
               albsnd_improad(begl:endl, :), &
               albsnd_perroad(begl:endl, :))

          ic = 1
          call UrbanSnowAlbedo(bounds, &
               num_urbanc, filter_urbanc, &
               coszen(begl:endl), &
               ic, &
               albsni_roof(begl:endl, :), &
               albsni_improad(begl:endl, :), &
               albsni_perroad(begl:endl, :))
       end if

       ! Combine snow-free and snow albedos
       do ib = 1,numrad
          do fc = 1,num_urbanc
             c = filter_urbanc(fc)
             l = clandunit(c)
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
       
       if (num_urbanl .gt. 0) then
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
               sref_roof_dif      (begl:endl, :))
       end if

       ! ----------------------------------------------------------------------------
       ! Map urban output to clmtype components 
       ! ----------------------------------------------------------------------------

       !  Set albgrd and albgri (ground albedos) and albd and albi (surface albedos)

       do ib = 1,numrad
          do fc = 1,num_urbanc
             c = filter_urbanc(fc)
             l = clandunit(c)
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
             c = pcolumn(p)
             albd(p,ib) = albgrd(c,ib)
             albi(p,ib) = albgri(c,ib)
          end do
       end do
    end if

  end associate
  end subroutine UrbanAlbedo

  !-----------------------------------------------------------------------
  subroutine UrbanSnowAlbedo (bounds, &
                              num_urbanc, filter_urbanc, coszen, ind, &
                              albsn_roof, albsn_improad, albsn_perroad)
    !
    ! !DESCRIPTION:
    ! Determine urban snow albedos
    !
    ! !USES:
    use clmtype
    use clm_varcon   , only : icol_roof, icol_road_perv, icol_road_imperv
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                    ! bounds
    integer , intent(in) :: num_urbanc                         ! number of urban columns in clump
    integer , intent(in) :: filter_urbanc(:)                   ! urban column filter
    integer , intent(in) :: ind                                ! 0=direct beam, 1=diffuse radiation
    real(r8), intent(in) :: coszen( bounds%begl: )             ! cosine solar zenith angle [landunit]
    real(r8), intent(out):: albsn_roof( bounds%begl: , 1: )    ! roof snow albedo by waveband [landunit, numrad]
    real(r8), intent(out):: albsn_improad( bounds%begl: , 1: ) ! impervious road snow albedo by waveband [landunit, numrad]
    real(r8), intent(out):: albsn_perroad( bounds%begl: , 1: ) ! pervious road snow albedo by waveband [landunit, numrad]
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
    SHR_ASSERT(numrad == 2, errMsg(__FILE__, __LINE__))

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(coszen)        == (/bounds%endl/)),         errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(albsn_roof)    == (/bounds%endl, numrad/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(albsn_improad) == (/bounds%endl, numrad/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(albsn_perroad) == (/bounds%endl, numrad/)), errMsg(__FILE__, __LINE__))

   associate(& 
   ctype                               =>    col%itype                                   , & ! Input:  [integer (:)]  column type                                        
   h2osno                              =>    cws%h2osno                                    & ! Input:  [real(r8) (:)]  snow water (mm H2O)                               
   )

    do fc = 1,num_urbanc
       c = filter_urbanc(fc)
       l = col%landunit(c)
       if (coszen(l) > 0._r8 .and. h2osno(c) > 0._r8) then
          if (ctype(c) == icol_roof) then
             albsn_roof(l,1) = snal0
             albsn_roof(l,2) = snal1
          else if (ctype(c) == icol_road_imperv) then
             albsn_improad(l,1) = snal0
             albsn_improad(l,2) = snal1
          else if (ctype(c) == icol_road_perv) then
             albsn_perroad(l,1) = snal0
             albsn_perroad(l,2) = snal1
          end if
       else
          if (ctype(c) == icol_roof) then
             albsn_roof(l,1) = 0._r8
             albsn_roof(l,2) = 0._r8
          else if (ctype(c) == icol_road_imperv) then
             albsn_improad(l,1) = 0._r8
             albsn_improad(l,2) = 0._r8
          else if (ctype(c) == icol_road_perv) then
             albsn_perroad(l,1) = 0._r8
             albsn_perroad(l,2) = 0._r8
          end if
       end if
    end do

    end associate 
  end subroutine UrbanSnowAlbedo

  !-----------------------------------------------------------------------
  subroutine UrbanRadiation (bounds, &
                             num_nourbanl, filter_nourbanl, &
                             num_urbanl, filter_urbanl, &
                             num_urbanc, filter_urbanc, &
                             num_urbanp, filter_urbanp)
    !
    ! !DESCRIPTION: 
    ! Solar fluxes absorbed and reflected by roof and canyon (walls, road).
    ! Also net and upward longwave fluxes.
    
    ! !USES:
    use clmtype
    use clm_varcon       , only : spval, icol_roof, icol_sunwall, icol_shadewall, &
                                  icol_road_perv, icol_road_imperv, sb
    use clm_varcon       , only : tfrz                ! To use new constant..
    use clm_time_manager , only : get_curr_date, get_step_size
    use clm_atmlnd       , only : clm_a2l, a2l_not_downscaled_gcell
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds    ! bounds
    integer , intent(in) :: num_nourbanl       ! number of non-urban landunits in clump
    integer , intent(in) :: filter_nourbanl(:) ! non-urban landunit filter
    integer , intent(in) :: num_urbanl         ! number of urban landunits in clump
    integer , intent(in) :: filter_urbanl(:)   ! urban landunit filter
    integer , intent(in) :: num_urbanc         ! number of urban columns in clump
    integer , intent(in) :: filter_urbanc(:)   ! urban column filter
    integer , intent(in) :: num_urbanp         ! number of urban pfts in clump
    integer , intent(in) :: filter_urbanp(:)   ! urban pft filter
    !
    ! !LOCAL VARIABLES:
    integer  :: fp,fl,p,c,l,g              ! indices
    integer  :: local_secp1                ! seconds into current date in local time
    real(r8) :: dtime                      ! land model time step (sec)
    integer  :: year,month,day             ! temporaries (not used)
    integer  :: secs                       ! seconds into current date

    real(r8), parameter :: mpe    = 1.e-06_r8 ! prevents overflow for division by zero
    real(r8), parameter :: snoem  = 0.97_r8   ! snow emissivity (should use value from Biogeophysics1)
                                 
    real(r8) :: lwnet_roof(bounds%begl:bounds%endl)     ! net (outgoing-incoming) longwave radiation (per unit ground area), roof (W/m**2)
    real(r8) :: lwnet_improad(bounds%begl:bounds%endl)  ! net (outgoing-incoming) longwave radiation (per unit ground area), impervious road (W/m**2)
    real(r8) :: lwnet_perroad(bounds%begl:bounds%endl)  ! net (outgoing-incoming) longwave radiation (per unit ground area), pervious road (W/m**2)
    real(r8) :: lwnet_sunwall(bounds%begl:bounds%endl)  ! net (outgoing-incoming) longwave radiation (per unit wall area), sunlit wall (W/m**2)
    real(r8) :: lwnet_shadewall(bounds%begl:bounds%endl)! net (outgoing-incoming) longwave radiation (per unit wall area), shaded wall (W/m**2)
    real(r8) :: lwnet_canyon(bounds%begl:bounds%endl)   ! net (outgoing-incoming) longwave radiation for canyon, per unit ground area (W/m**2)
    real(r8) :: lwup_roof(bounds%begl:bounds%endl)      ! upward longwave radiation (per unit ground area), roof (W/m**2)
    real(r8) :: lwup_improad(bounds%begl:bounds%endl)   ! upward longwave radiation (per unit ground area), impervious road (W/m**2)
    real(r8) :: lwup_perroad(bounds%begl:bounds%endl)   ! upward longwave radiation (per unit ground area), pervious road (W/m**2)
    real(r8) :: lwup_sunwall(bounds%begl:bounds%endl)   ! upward longwave radiation, (per unit wall area), sunlit wall (W/m**2)
    real(r8) :: lwup_shadewall(bounds%begl:bounds%endl) ! upward longwave radiation, (per unit wall area), shaded wall (W/m**2)
    real(r8) :: lwup_canyon(bounds%begl:bounds%endl)    ! upward longwave radiation for canyon, per unit ground area (W/m**2)
    real(r8) :: t_roof(bounds%begl:bounds%endl)         ! roof temperature (K)
    real(r8) :: t_improad(bounds%begl:bounds%endl)      ! imppervious road temperature (K)
    real(r8) :: t_perroad(bounds%begl:bounds%endl)      ! pervious road temperature (K)
    real(r8) :: t_sunwall(bounds%begl:bounds%endl)      ! sunlit wall temperature (K)
    real(r8) :: t_shadewall(bounds%begl:bounds%endl)    ! shaded wall temperature (K)
    real(r8) :: lwdown(bounds%begl:bounds%endl)         ! atmospheric downward longwave radiation (W/m**2)
    real(r8) :: em_roof_s(bounds%begl:bounds%endl)      ! roof emissivity with snow effects
    real(r8) :: em_improad_s(bounds%begl:bounds%endl)   ! impervious road emissivity with snow effects
    real(r8) :: em_perroad_s(bounds%begl:bounds%endl)   ! pervious road emissivity with snow effects
    !-----------------------------------------------------------------------

   associate(& 
   canyon_hwr                          =>    lun%canyon_hwr                              , & ! Input:  [real(r8) (:)]  ratio of building height to street width          
   wtroad_perv                         =>    lun%wtroad_perv                             , & ! Input:  [real(r8) (:)]  weight of pervious road wrt total road            
   em_roof                             =>    urban_params%em_roof                        , & ! Input:  [real(r8) (:)]  roof emissivity                                   
   em_improad                          =>    urban_params%em_improad                     , & ! Input:  [real(r8) (:)]  impervious road emissivity                        
   em_perroad                          =>    urban_params%em_perroad                     , & ! Input:  [real(r8) (:)]  pervious road emissivity                          
   em_wall                             =>    urban_params%em_wall                        , & ! Input:  [real(r8) (:)]  wall emissivity                                   
   londeg                              =>    grc%londeg                                  , & ! Input:  [real(r8) (:)]  longitude (degrees)                               
   forc_solad                          =>    clm_a2l%forc_solad                          , & ! Input:  [real(r8) (:,:)]  direct beam radiation  (vis=forc_sols , nir=forc_soll ) (W/m**2)
   forc_solai                          =>    clm_a2l%forc_solai                          , & ! Input:  [real(r8) (:,:)]  diffuse beam radiation (vis=forc_sols , nir=forc_soll ) (W/m**2)
   forc_solar                          =>    clm_a2l%forc_solar                          , & ! Input:  [real(r8) (:)]  incident solar radiation (W/m**2)                 
   forc_lwrad                          =>    a2l_not_downscaled_gcell%forc_lwrad           , & ! Input:  [real(r8) (:)]  downward infrared (longwave) radiation (W/m**2)   
   pfti                                =>   lun%pfti                                     , & ! Input:  [integer (:)]  beginning pfti index for landunit                  
   pftf                                =>   lun%pftf                                     , & ! Input:  [integer (:)]  ending pftf index for landunit                     
   coli                                =>   lun%coli                                     , & ! Input:  [integer (:)]  beginning column index for landunit                
   colf                                =>   lun%colf                                     , & ! Input:  [integer (:)]  ending column index for landunit                   
   lgridcell                           =>   lun%gridcell                                 , & ! Input:  [integer (:)]  gridcell of corresponding landunit                 
   sabs_roof_dir                       =>    lps%sabs_roof_dir                           , & ! Input:  [real(r8) (:,:)]  direct  solar absorbed  by roof per unit ground area per unit incident flux
   sabs_roof_dif                       =>    lps%sabs_roof_dif                           , & ! Input:  [real(r8) (:,:)]  diffuse solar absorbed  by roof per unit ground area per unit incident flux
   sabs_sunwall_dir                    =>    lps%sabs_sunwall_dir                        , & ! Input:  [real(r8) (:,:)]  direct  solar absorbed  by sunwall per unit wall area per unit incident flux
   sabs_sunwall_dif                    =>    lps%sabs_sunwall_dif                        , & ! Input:  [real(r8) (:,:)]  diffuse solar absorbed  by sunwall per unit wall area per unit incident flux
   sabs_shadewall_dir                  =>    lps%sabs_shadewall_dir                      , & ! Input:  [real(r8) (:,:)]  direct  solar absorbed  by shadewall per unit wall area per unit incident flux
   sabs_shadewall_dif                  =>    lps%sabs_shadewall_dif                      , & ! Input:  [real(r8) (:,:)]  diffuse solar absorbed  by shadewall per unit wall area per unit incident flux
   sabs_improad_dir                    =>    lps%sabs_improad_dir                        , & ! Input:  [real(r8) (:,:)]  direct  solar absorbed  by impervious road per unit ground area per unit incident flux
   sabs_improad_dif                    =>    lps%sabs_improad_dif                        , & ! Input:  [real(r8) (:,:)]  diffuse solar absorbed  by impervious road per unit ground area per unit incident flux
   sabs_perroad_dir                    =>    lps%sabs_perroad_dir                        , & ! Input:  [real(r8) (:,:)]  direct  solar absorbed  by pervious road per unit ground area per unit incident flux
   sabs_perroad_dif                    =>    lps%sabs_perroad_dif                        , & ! Input:  [real(r8) (:,:)]  diffuse solar absorbed  by pervious road per unit ground area per unit incident flux
   ctype                               =>    col%itype                                   , & ! Input:  [integer (:)]  column type                                        
   t_grnd                              =>    ces%t_grnd                                  , & ! Input:  [real(r8) (:)]  ground temperature (K)                            
   frac_sno                            =>    cps%frac_sno                                , & ! Input:  [real(r8) (:)]  fraction of ground covered by snow (0 to 1)       
   pgridcell                           =>   pft%gridcell                                 , & ! Input:  [integer (:)]  gridcell of corresponding pft                      
   plandunit                           =>   pft%landunit                                 , & ! Input:  [integer (:)]  landunit of corresponding pft                      
   pcolumn                             =>   pft%column                                   , & ! Input:  [integer (:)]  column of corresponding pft                        
   albd                                =>    pps%albd                                    , & ! Input:  [real(r8) (:,:)]  surface albedo (direct)                         
   albi                                =>    pps%albi                                    , & ! Input:  [real(r8) (:,:)]  surface albedo (diffuse)                        
   sabg                                =>    pef%sabg                                    , & ! Output: [real(r8) (:)]  solar radiation absorbed by ground (W/m**2)       
   sabv                                =>    pef%sabv                                    , & ! Output: [real(r8) (:)]  solar radiation absorbed by vegetation (W/m**2)   
   fsa                                 =>    pef%fsa                                     , & ! Output: [real(r8) (:)]  solar radiation absorbed (total) (W/m**2)         
   fsa_u                               =>    pef%fsa_u                                   , & ! Output: [real(r8) (:)]  urban solar radiation absorbed (total) (W/m**2)   
   fsr                                 =>    pef%fsr                                     , & ! Output: [real(r8) (:)]  solar radiation reflected (total) (W/m**2)        
   fsds_vis_d                          =>    pef%fsds_vis_d                              , & ! Output: [real(r8) (:)]  incident direct beam vis solar radiation (W/m**2) 
   fsds_nir_d                          =>    pef%fsds_nir_d                              , & ! Output: [real(r8) (:)]  incident direct beam nir solar radiation (W/m**2) 
   fsds_vis_i                          =>    pef%fsds_vis_i                              , & ! Output: [real(r8) (:)]  incident diffuse vis solar radiation (W/m**2)     
   fsds_nir_i                          =>    pef%fsds_nir_i                              , & ! Output: [real(r8) (:)]  incident diffuse nir solar radiation (W/m**2)     
   fsr_vis_d                           =>    pef%fsr_vis_d                               , & ! Output: [real(r8) (:)]  reflected direct beam vis solar radiation (W/m**2)
   fsr_nir_d                           =>    pef%fsr_nir_d                               , & ! Output: [real(r8) (:)]  reflected direct beam nir solar radiation (W/m**2)
   fsr_vis_i                           =>    pef%fsr_vis_i                               , & ! Output: [real(r8) (:)]  reflected diffuse vis solar radiation (W/m**2)    
   fsr_nir_i                           =>    pef%fsr_nir_i                               , & ! Output: [real(r8) (:)]  reflected diffuse nir solar radiation (W/m**2)    
   fsds_vis_d_ln                       =>    pef%fsds_vis_d_ln                           , & ! Output: [real(r8) (:)]  incident direct beam vis solar rad at local noon (W/m**2)
   fsds_nir_d_ln                       =>    pef%fsds_nir_d_ln                           , & ! Output: [real(r8) (:)]  incident direct beam nir solar rad at local noon (W/m**2)
   fsds_vis_i_ln                       =>    pef%fsds_vis_i_ln                           , & ! Output: [real(r8) (:)]  incident diffuse beam vis solar rad at local noon (W/m**2)
   parveg_ln                           =>    pef%parveg_ln                               , & ! Output: [real(r8) (:)]  absorbed par by vegetation at local noon (W/m**2) 
   fsr_vis_d_ln                        =>    pef%fsr_vis_d_ln                            , & ! Output: [real(r8) (:)]  reflected direct beam vis solar rad at local noon (W/m**2)
   fsr_nir_d_ln                        =>    pef%fsr_nir_d_ln                            , & ! Output: [real(r8) (:)]  reflected direct beam nir solar rad at local noon (W/m**2)
   eflx_lwrad_out                      =>    pef%eflx_lwrad_out                          , & ! Output: [real(r8) (:)]  emitted infrared (longwave) radiation (W/m**2)    
   eflx_lwrad_net                      =>    pef%eflx_lwrad_net                          , & ! Output: [real(r8) (:)]  net infrared (longwave) rad (W/m**2) [+ = to atm] 
   eflx_lwrad_net_u                    =>    pef%eflx_lwrad_net_u                        , & ! Output: [real(r8) (:)]  urban net infrared (longwave) rad (W/m**2) [+ = to atm]
   t_ref2m                             =>    pes%t_ref2m                                 , & ! Input:  [real(r8) (:)]  2 m height surface air temperature (K)            
   begl                                =>    bounds%begl                                 , &
   endl                                =>    bounds%endl                                   &
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
       g = lgridcell(l) 

       ! Need to set the following temperatures to some defined value even if it
       ! does not appear in the urban landunit for the net_longwave computation

       t_roof(l)      = 19._r8 + tfrz
       t_sunwall(l)   = 19._r8 + tfrz
       t_shadewall(l) = 19._r8 + tfrz
       t_improad(l)   = 19._r8 + tfrz
       t_perroad(l)   = 19._r8 + tfrz

       ! Initial assignment of emissivity
       em_roof_s(l)    = em_roof(l)
       em_improad_s(l) = em_improad(l)
       em_perroad_s(l) = em_perroad(l)

       ! Set urban temperatures and emissivity including snow effects.
       do c = coli(l),colf(l)
          if (ctype(c) == icol_roof       )  then
             t_roof(l)      = t_grnd(c)
             em_roof_s(l) = em_roof(l)*(1._r8-frac_sno(c)) + snoem*frac_sno(c)
          else if (ctype(c) == icol_road_imperv) then 
             t_improad(l)   = t_grnd(c)
             em_improad_s(l) = em_improad(l)*(1._r8-frac_sno(c)) + snoem*frac_sno(c)
          else if (ctype(c) == icol_road_perv  ) then
             t_perroad(l)   = t_grnd(c)
             em_perroad_s(l) = em_perroad(l)*(1._r8-frac_sno(c)) + snoem*frac_sno(c)
          else if (ctype(c) == icol_sunwall    ) then
             t_sunwall(l)   = t_grnd(c)
          else if (ctype(c) == icol_shadewall  ) then
             t_shadewall(l) = t_grnd(c)
          end if
       end do
       lwdown(l) = forc_lwrad(g)
    end do

    ! Net longwave radiation for road and both walls in urban canyon allowing for multiple re-emission
    
    if (num_urbanl .gt. 0) then
       call net_longwave (bounds,       &
            num_urbanl, filter_urbanl,  &
            canyon_hwr(begl:endl),      &
            wtroad_perv(begl:endl),     &
            lwdown(begl:endl),          &
            em_roof_s(begl:endl),       &
            em_improad_s(begl:endl),    &
            em_perroad_s(begl:endl),    &
            em_wall(begl:endl),         &
            t_roof(begl:endl),          &
            t_improad(begl:endl),       &
            t_perroad(begl:endl),       &
            t_sunwall(begl:endl),       &
            t_shadewall(begl:endl),     &
            lwnet_roof(begl:endl),      &
            lwnet_improad(begl:endl),   &
            lwnet_perroad(begl:endl),   &
            lwnet_sunwall(begl:endl),   &
            lwnet_shadewall(begl:endl), &
            lwnet_canyon(begl:endl),    &
            lwup_roof(begl:endl),       &
            lwup_improad(begl:endl),    &
            lwup_perroad(begl:endl),    &
            lwup_sunwall(begl:endl),    &
            lwup_shadewall(begl:endl),  &
            lwup_canyon(begl:endl))
    end if
    
    dtime = get_step_size()
    call get_curr_date (year, month, day, secs)
    
    ! Determine clmtype variables needed for history output and communication with atm
    ! Loop over urban pfts in clump

    do fp = 1,num_urbanp
       p = filter_urbanp(fp)
       g = pgridcell(p)
       
       local_secp1 = secs + nint((grc%londeg(g)/degpsec)/dtime)*dtime
       local_secp1 = mod(local_secp1,isecspday)

       ! Solar incident 
       
       fsds_vis_d(p) = forc_solad(g,1)
       fsds_nir_d(p) = forc_solad(g,2)    
       fsds_vis_i(p) = forc_solai(g,1)
       fsds_nir_i(p) = forc_solai(g,2)
       ! Determine local noon incident solar
       if (local_secp1 == noonsec) then
          fsds_vis_d_ln(p) = forc_solad(g,1)
          fsds_nir_d_ln(p) = forc_solad(g,2)
          fsds_vis_i_ln(p) = forc_solai(g,1)
          parveg_ln(p)     = 0._r8
       else
          fsds_vis_d_ln(p) = spval 
          fsds_nir_d_ln(p) = spval 
          fsds_vis_i_ln(p) = spval
          parveg_ln(p)     = spval
       endif
       
       ! Solar reflected 
       ! per unit ground area (roof, road) and per unit wall area (sunwall, shadewall)

       fsr_vis_d(p) = albd(p,1) * forc_solad(g,1)
       fsr_nir_d(p) = albd(p,2) * forc_solad(g,2)
       fsr_vis_i(p) = albi(p,1) * forc_solai(g,1)
       fsr_nir_i(p) = albi(p,2) * forc_solai(g,2)

       ! Determine local noon reflected solar
       if (local_secp1 == noonsec) then
          fsr_vis_d_ln(p) = fsr_vis_d(p)
          fsr_nir_d_ln(p) = fsr_nir_d(p)
       else
          fsr_vis_d_ln(p) = spval 
          fsr_nir_d_ln(p) = spval 
       endif
       fsr(p) = fsr_vis_d(p) + fsr_nir_d(p) + fsr_vis_i(p) + fsr_nir_i(p)  
    end do


    ! Loop over urban pfts in clump

    do fp = 1,num_urbanp
       p = filter_urbanp(fp)
       c = pcolumn(p)
       l = plandunit(p)
       g = pgridcell(p)

       ! Solar absorbed and longwave out and net
       ! per unit ground area (roof, road) and per unit wall area (sunwall, shadewall)
       ! Each urban pft has its own column - this is used in the logic below
 
       if (ctype(c) == icol_roof) then   
          eflx_lwrad_out(p) = lwup_roof(l)
          eflx_lwrad_net(p) = lwnet_roof(l)
          eflx_lwrad_net_u(p) = lwnet_roof(l)
          sabg(p) = sabs_roof_dir(l,1)*forc_solad(g,1) + &
                    sabs_roof_dif(l,1)*forc_solai(g,1) + &
                    sabs_roof_dir(l,2)*forc_solad(g,2) + &
                    sabs_roof_dif(l,2)*forc_solai(g,2) 
          
       else if (ctype(c) == icol_sunwall) then   
          eflx_lwrad_out(p) = lwup_sunwall(l)
          eflx_lwrad_net(p) = lwnet_sunwall(l)
          eflx_lwrad_net_u(p) = lwnet_sunwall(l)
          sabg(p) = sabs_sunwall_dir(l,1)*forc_solad(g,1) + &
                    sabs_sunwall_dif(l,1)*forc_solai(g,1) + &
                    sabs_sunwall_dir(l,2)*forc_solad(g,2) + &
                    sabs_sunwall_dif(l,2)*forc_solai(g,2) 
       else if (ctype(c) == icol_shadewall) then   
          eflx_lwrad_out(p) = lwup_shadewall(l)
          eflx_lwrad_net(p) = lwnet_shadewall(l)
          eflx_lwrad_net_u(p) = lwnet_shadewall(l)
          sabg(p) = sabs_shadewall_dir(l,1)*forc_solad(g,1) + &
                    sabs_shadewall_dif(l,1)*forc_solai(g,1) + &
                    sabs_shadewall_dir(l,2)*forc_solad(g,2) + &
                    sabs_shadewall_dif(l,2)*forc_solai(g,2) 
       else if (ctype(c) == icol_road_perv) then       
          eflx_lwrad_out(p) = lwup_perroad(l)
          eflx_lwrad_net(p) = lwnet_perroad(l)
          eflx_lwrad_net_u(p) = lwnet_perroad(l)
          sabg(p) = sabs_perroad_dir(l,1)*forc_solad(g,1) + &
                    sabs_perroad_dif(l,1)*forc_solai(g,1) + &
                    sabs_perroad_dir(l,2)*forc_solad(g,2) + &
                    sabs_perroad_dif(l,2)*forc_solai(g,2) 
       else if (ctype(c) == icol_road_imperv) then       
          eflx_lwrad_out(p) = lwup_improad(l)
          eflx_lwrad_net(p) = lwnet_improad(l)
          eflx_lwrad_net_u(p) = lwnet_improad(l)
          sabg(p) = sabs_improad_dir(l,1)*forc_solad(g,1) + &
                    sabs_improad_dif(l,1)*forc_solai(g,1) + &
                    sabs_improad_dir(l,2)*forc_solad(g,2) + &
                    sabs_improad_dif(l,2)*forc_solai(g,2) 
       end if
       sabv(p)   = 0._r8
       fsa(p)    = sabv(p) + sabg(p)
       fsa_u(p)  = fsa(p)
    end do ! end loop over urban pfts

    end associate 
  end subroutine UrbanRadiation

  !-----------------------------------------------------------------------
  subroutine incident_direct (bounds, &
       num_urbanl, filter_urbanl, canyon_hwr, coszen, zen, &
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
    use clmtype
    use clm_varcon, only : rpi
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                      ! bounds
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
    integer  :: fl,l,i,ib          ! indices
    logical  :: numchk = .false.   ! true => perform numerical check of analytical solution
    real(r8) :: theta0(bounds%begl:bounds%endl)    ! critical canyon orientation for which road is no longer illuminated
    real(r8) :: tanzen(bounds%begl:bounds%endl)    ! tan(zenith angle)
    real(r8) :: swall_projected    ! direct beam solar radiation (per unit ground area) incident on wall
    real(r8) :: err1(bounds%begl:bounds%endl)      ! energy conservation error
    real(r8) :: err2(bounds%begl:bounds%endl)      ! energy conservation error
    real(r8) :: err3(bounds%begl:bounds%endl)      ! energy conservation error
    real(r8) :: sumr               ! sum of sroad for each orientation (0 <= theta <= pi/2)
    real(r8) :: sumw               ! sum of swall for each orientation (0 <= theta <= pi/2)
    real(r8) :: num                ! number of orientations
    real(r8) :: theta              ! canyon orientation relative to sun (0 <= theta <= pi/2)
    real(r8) :: zen0               ! critical solar zenith angle for which sun begins to illuminate road
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
          
             sdir_road(l,ib) = sdir(l,ib) * &
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
       num_urbanl, filter_urbanl, canyon_hwr, sdif, sdif_road, sdif_sunwall, sdif_shadewall)
    !
    ! !DESCRIPTION: 
    ! Diffuse solar radiation incident on walls and road in urban canyon
    ! Conservation check: Total incoming diffuse 
    ! (sdif) = sdif_road + (sdif_shadewall + sdif_sunwall)*canyon_hwr
    ! Multiplication by canyon_hwr scales wall fluxes (per unit wall area) to per unit ground area
    !
    ! !USES:
    use clmtype
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                      ! bounds
    integer , intent(in)  :: num_urbanl                          ! number of urban landunits
    integer , intent(in)  :: filter_urbanl(:)                    ! urban landunit filter
    real(r8), intent(in)  :: canyon_hwr( bounds%begl: )          ! ratio of building height to street width [landunit]
    real(r8), intent(in)  :: sdif( bounds%begl: , 1: )           ! diffuse solar radiation incident on horizontal surface [landunit, numrad]
    real(r8), intent(out) :: sdif_road( bounds%begl: , 1: )      ! diffuse solar radiation incident on road [landunit, numrad]
    real(r8), intent(out) :: sdif_sunwall( bounds%begl: , 1: )   ! diffuse solar radiation (per unit wall area) incident on sunlit wall [landunit, numrad]
    real(r8), intent(out) :: sdif_shadewall( bounds%begl: , 1: ) ! diffuse solar radiation (per unit wall area) incident on shaded wall [landunit, numrad]
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

   associate(& 
   vf_sr =>    lps%vf_sr , & ! Input:  [real(r8) (:)]  view factor of sky for road                       
   vf_sw =>    lps%vf_sw   & ! Input:  [real(r8) (:)]  view factor of sky for one wall                   
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
  subroutine net_solar (bounds, &
       num_urbanl, filter_urbanl, coszen, canyon_hwr, wtroad_perv, sdir, sdif, &
       alb_improad_dir, alb_perroad_dir, alb_wall_dir, alb_roof_dir, &
       alb_improad_dif, alb_perroad_dif, alb_wall_dif, alb_roof_dif, &
       sdir_road, sdir_sunwall, sdir_shadewall,  &
       sdif_road, sdif_sunwall, sdif_shadewall,  &
       sref_improad_dir, sref_perroad_dir, sref_sunwall_dir, sref_shadewall_dir, sref_roof_dir, &
       sref_improad_dif, sref_perroad_dif, sref_sunwall_dif, sref_shadewall_dif, sref_roof_dif)
    !
    ! !DESCRIPTION: 
    ! Solar radiation absorbed by road and both walls in urban canyon allowing 
    ! for multiple reflection.
    !
    ! !USES:
    use clmtype
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                            ! bounds
    integer , intent(in)  :: num_urbanl                                ! number of urban landunits
    integer , intent(in)  :: filter_urbanl(:)                          ! urban landunit filter
    real(r8), intent(in)  :: coszen( bounds%begl: )                    ! cosine solar zenith angle [landunit]
    real(r8), intent(in)  :: canyon_hwr( bounds%begl: )                ! ratio of building height to street width [landunit]
    real(r8), intent(in)  :: wtroad_perv( bounds%begl: )               ! weight of pervious road wrt total road [landunit]
    real(r8), intent(in)  :: sdir( bounds%begl: , 1: )                 ! direct beam solar radiation incident on horizontal surface [landunit, numrad]
    real(r8), intent(in)  :: sdif( bounds%begl: , 1: )                 ! diffuse solar radiation on horizontal surface [landunit, numrad]
    real(r8), intent(in)  :: alb_improad_dir( bounds%begl: , 1: )      ! direct impervious road albedo [landunit, numrad]
    real(r8), intent(in)  :: alb_perroad_dir( bounds%begl: , 1: )      ! direct pervious road albedo [landunit, numrad]
    real(r8), intent(in)  :: alb_wall_dir( bounds%begl: , 1: )         ! direct  wall albedo [landunit, numrad]
    real(r8), intent(in)  :: alb_roof_dir( bounds%begl: , 1: )         ! direct  roof albedo [landunit, numrad]
    real(r8), intent(in)  :: alb_improad_dif( bounds%begl: , 1: )      ! diffuse impervious road albedo [landunit, numrad]
    real(r8), intent(in)  :: alb_perroad_dif( bounds%begl: , 1: )      ! diffuse pervious road albedo [landunit, numrad]
    real(r8), intent(in)  :: alb_wall_dif( bounds%begl: , 1: )         ! diffuse wall albedo [landunit, numrad]
    real(r8), intent(in)  :: alb_roof_dif( bounds%begl: , 1: )         ! diffuse roof albedo [landunit, numrad]
    real(r8), intent(in)  :: sdir_road( bounds%begl: , 1: )            ! direct beam solar radiation incident on road per unit incident flux [landunit, numrad]
    real(r8), intent(in)  :: sdir_sunwall( bounds%begl: , 1: )         ! direct beam solar radiation (per unit wall area) incident on sunlit wall per unit incident flux [landunit, numrad]
    real(r8), intent(in)  :: sdir_shadewall( bounds%begl: , 1: )       ! direct beam solar radiation (per unit wall area) incident on shaded wall per unit incident flux [landunit, numrad]
    real(r8), intent(in)  :: sdif_road( bounds%begl: , 1: )            ! diffuse solar radiation incident on road per unit incident flux [landunit, numrad]
    real(r8), intent(in)  :: sdif_sunwall( bounds%begl: , 1: )         ! diffuse solar radiation (per unit wall area) incident on sunlit wall per unit incident flux [landunit, numrad]
    real(r8), intent(in)  :: sdif_shadewall( bounds%begl: , 1: )       ! diffuse solar radiation (per unit wall area) incident on shaded wall per unit incident flux [landunit, numrad]
    real(r8), intent(inout) :: sref_improad_dir( bounds%begl: , 1: )   ! direct  solar rad reflected by impervious road (per unit ground area) per unit incident flux [landunit, numrad]
    real(r8), intent(inout) :: sref_perroad_dir( bounds%begl: , 1: )   ! direct  solar rad reflected by pervious road (per unit ground area) per unit incident flux [landunit, numrad]
    real(r8), intent(inout) :: sref_improad_dif( bounds%begl: , 1: )   ! diffuse solar rad reflected by impervious road (per unit ground area) per unit incident flux [landunit, numrad]
    real(r8), intent(inout) :: sref_perroad_dif( bounds%begl: , 1: )   ! diffuse solar rad reflected by pervious road (per unit ground area) per unit incident flux [landunit, numrad]
    real(r8), intent(inout) :: sref_sunwall_dir( bounds%begl: , 1: )   ! direct solar  rad reflected by sunwall (per unit wall area) per unit incident flux [landunit, numrad]
    real(r8), intent(inout) :: sref_sunwall_dif( bounds%begl: , 1: )   ! diffuse solar rad reflected by sunwall (per unit wall area) per unit incident flux [landunit, numrad]
    real(r8), intent(inout) :: sref_shadewall_dir( bounds%begl: , 1: ) ! direct solar  rad reflected by shadewall (per unit wall area) per unit incident flux [landunit, numrad]
    real(r8), intent(inout) :: sref_shadewall_dif( bounds%begl: , 1: ) ! diffuse solar rad reflected by shadewall (per unit wall area) per unit incident flux [landunit, numrad]
    real(r8), intent(inout) :: sref_roof_dir( bounds%begl: , 1: )      ! direct  solar rad reflected by roof (per unit ground area) per unit incident flux [landunit, numrad]
    real(r8), intent(inout) :: sref_roof_dif( bounds%begl: , 1: )      ! diffuse solar rad reflected by roof (per unit ground area)  per unit incident flux [landunit, numrad]
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

   associate(& 
   vf_sr                               =>    lps%vf_sr                                   , & ! Input:  [real(r8) (:)]  view factor of sky for road                       
   vf_wr                               =>    lps%vf_wr                                   , & ! Input:  [real(r8) (:)]  view factor of one wall for road                  
   vf_sw                               =>    lps%vf_sw                                   , & ! Input:  [real(r8) (:)]  view factor of sky for one wall                   
   vf_rw                               =>    lps%vf_rw                                   , & ! Input:  [real(r8) (:)]  view factor of road for one wall                  
   vf_ww                               =>    lps%vf_ww                                   , & ! Input:  [real(r8) (:)]  view factor of opposing wall for one wall         
   sabs_roof_dir                       =>    lps%sabs_roof_dir                           , & ! Input:  [real(r8) (:,:)]  direct  solar absorbed  by roof per unit ground area per unit incident flux
   sabs_roof_dif                       =>    lps%sabs_roof_dif                           , & ! Input:  [real(r8) (:,:)]  diffuse solar absorbed  by roof per unit ground area per unit incident flux
   sabs_sunwall_dir                    =>    lps%sabs_sunwall_dir                        , & ! Input:  [real(r8) (:,:)]  direct  solar absorbed  by sunwall per unit wall area per unit incident flux
   sabs_sunwall_dif                    =>    lps%sabs_sunwall_dif                        , & ! Input:  [real(r8) (:,:)]  diffuse solar absorbed  by sunwall per unit wall area per unit incident flux
   sabs_shadewall_dir                  =>    lps%sabs_shadewall_dir                      , & ! Input:  [real(r8) (:,:)]  direct  solar absorbed  by shadewall per unit wall area per unit incident flux
   sabs_shadewall_dif                  =>    lps%sabs_shadewall_dif                      , & ! Input:  [real(r8) (:,:)]  diffuse solar absorbed  by shadewall per unit wall area per unit incident flux
   sabs_improad_dir                    =>    lps%sabs_improad_dir                        , & ! Input:  [real(r8) (:,:)]  direct  solar absorbed  by impervious road per unit ground area per unit incident flux
   sabs_improad_dif                    =>    lps%sabs_improad_dif                        , & ! Input:  [real(r8) (:,:)]  diffuse solar absorbed  by impervious road per unit ground area per unit incident flux
   sabs_perroad_dir                    =>    lps%sabs_perroad_dir                        , & ! Input:  [real(r8) (:,:)]  direct  solar absorbed  by pervious road per unit ground area per unit incident flux
   sabs_perroad_dif                    =>    lps%sabs_perroad_dif                          & ! Input:  [real(r8) (:,:)]  diffuse solar absorbed  by pervious road per unit ground area per unit incident flux
   )

    ! Calculate impervious road

    do fl = 1,num_urbanl 
       l = filter_urbanl(fl)
       wtroad_imperv(l) = 1._r8 - wtroad_perv(l)
    end do

    do ib = 1,numrad
       do fl = 1,num_urbanl
          l = filter_urbanl(fl)
          if (coszen(l) .gt. 0._r8) then

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
          if (coszen(l) .gt. 0._r8) then

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
          if (coszen(l) .gt. 0._r8) then
             sref_roof_dir(l,ib) = alb_roof_dir(l,ib) * sdir(l,ib)
             sref_roof_dif(l,ib) = alb_roof_dif(l,ib) * sdif(l,ib)
             sabs_roof_dir(l,ib) = sdir(l,ib) - sref_roof_dir(l,ib)
             sabs_roof_dif(l,ib) = sdif(l,ib) - sref_roof_dif(l,ib)
          end if
       end do

    end do   ! end of radiation band loop

    end associate 
  end subroutine net_solar

  !-----------------------------------------------------------------------
  subroutine net_longwave (bounds, &
       num_urbanl, filter_urbanl, canyon_hwr, wtroad_perv, &
       lwdown, em_roof, em_improad, em_perroad, em_wall, &
       t_roof,  t_improad, t_perroad, t_sunwall, t_shadewall, &
       lwnet_roof, lwnet_improad, lwnet_perroad, lwnet_sunwall, lwnet_shadewall, lwnet_canyon, &
       lwup_roof, lwup_improad, lwup_perroad, lwup_sunwall, lwup_shadewall, lwup_canyon)
    !
    ! !DESCRIPTION: 
    ! Net longwave radiation for road and both walls in urban canyon allowing for 
    ! multiple reflection. Also net longwave radiation for urban roof. 
    !
    ! !USES:
    use clm_varcon   , only : sb
    use clmtype
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                  ! bounds
    integer , intent(in)  :: num_urbanl                      ! number of urban landunits
    integer , intent(in)  :: filter_urbanl(:)                ! urban landunit filter
    real(r8), intent(in)  :: canyon_hwr( bounds%begl: )      ! ratio of building height to street width [landunit]
    real(r8), intent(in)  :: wtroad_perv( bounds%begl: )     ! weight of pervious road wrt total road [landunit]

    real(r8), intent(in)  :: lwdown( bounds%begl: )          ! atmospheric longwave radiation (W/m**2) [landunit]
    real(r8), intent(in)  :: em_roof( bounds%begl: )         ! roof emissivity [landunit]
    real(r8), intent(in)  :: em_improad( bounds%begl: )      ! impervious road emissivity [landunit]
    real(r8), intent(in)  :: em_perroad( bounds%begl: )      ! pervious road emissivity [landunit]
    real(r8), intent(in)  :: em_wall( bounds%begl: )         ! wall emissivity [landunit]

    real(r8), intent(in)  :: t_roof( bounds%begl: )          ! roof temperature (K) [landunit]
    real(r8), intent(in)  :: t_improad( bounds%begl: )       ! impervious road temperature (K) [landunit]
    real(r8), intent(in)  :: t_perroad( bounds%begl: )       ! ervious road temperature (K) [landunit]
    real(r8), intent(in)  :: t_sunwall( bounds%begl: )       ! sunlit wall temperature (K) [landunit]
    real(r8), intent(in)  :: t_shadewall( bounds%begl: )     ! shaded wall temperature (K) [landunit]

    real(r8), intent(out) :: lwnet_roof( bounds%begl: )      ! net (outgoing-incoming) longwave radiation, roof (W/m**2) [landunit]
    real(r8), intent(out) :: lwnet_improad( bounds%begl: )   ! net (outgoing-incoming) longwave radiation, impervious road (W/m**2) [landunit]
    real(r8), intent(out) :: lwnet_perroad( bounds%begl: )   ! net (outgoing-incoming) longwave radiation, pervious road (W/m**2) [landunit]
    real(r8), intent(out) :: lwnet_sunwall( bounds%begl: )   ! net (outgoing-incoming) longwave radiation (per unit wall area), sunlit wall (W/m**2) [landunit]
    real(r8), intent(out) :: lwnet_shadewall( bounds%begl: ) ! net (outgoing-incoming) longwave radiation (per unit wall area), shaded wall (W/m**2) [landunit]
    real(r8), intent(out) :: lwnet_canyon( bounds%begl: )    ! net (outgoing-incoming) longwave radiation for canyon, per unit ground area (W/m**2) [landunit]

    real(r8), intent(out) :: lwup_roof( bounds%begl: )       ! upward longwave radiation, roof (W/m**2) [landunit]
    real(r8), intent(out) :: lwup_improad( bounds%begl: )    ! upward longwave radiation, impervious road (W/m**2) [landunit]
    real(r8), intent(out) :: lwup_perroad( bounds%begl: )    ! upward longwave radiation, pervious road (W/m**2) [landunit]
    real(r8), intent(out) :: lwup_sunwall( bounds%begl: )    ! upward longwave radiation (per unit wall area), sunlit wall (W/m**2) [landunit]
    real(r8), intent(out) :: lwup_shadewall( bounds%begl: )  ! upward longwave radiation (per unit wall area), shaded wall (W/m**2) [landunit]
    real(r8), intent(out) :: lwup_canyon( bounds%begl: )     ! upward longwave radiation for canyon, per unit ground area (W/m**2) [landunit]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: lwdown_road(bounds%begl:bounds%endl)         ! atmospheric longwave radiation for total road (W/m**2)
    real(r8) :: lwdown_sunwall(bounds%begl:bounds%endl)      ! atmospheric longwave radiation (per unit wall area) for sunlit wall (W/m**2)
    real(r8) :: lwdown_shadewall(bounds%begl:bounds%endl)    ! atmospheric longwave radiation (per unit wall area) for shaded wall (W/m**2)
    real(r8) :: lwtot(bounds%begl:bounds%endl)               ! incoming longwave radiation (W/m**2)

    real(r8) :: improad_a(bounds%begl:bounds%endl)           ! absorbed longwave for improad (W/m**2)
    real(r8) :: improad_r(bounds%begl:bounds%endl)           ! reflected longwave for improad (W/m**2)
    real(r8) :: improad_r_sky(bounds%begl:bounds%endl)       ! improad_r to sky (W/m**2)
    real(r8) :: improad_r_sunwall(bounds%begl:bounds%endl)   ! improad_r to sunlit wall (W/m**2)
    real(r8) :: improad_r_shadewall(bounds%begl:bounds%endl) ! improad_r to shaded wall (W/m**2)
    real(r8) :: improad_e(bounds%begl:bounds%endl)           ! emitted longwave for improad (W/m**2)
    real(r8) :: improad_e_sky(bounds%begl:bounds%endl)       ! improad_e to sky (W/m**2)
    real(r8) :: improad_e_sunwall(bounds%begl:bounds%endl)   ! improad_e to sunlit wall (W/m**2)
    real(r8) :: improad_e_shadewall(bounds%begl:bounds%endl) ! improad_e to shaded wall (W/m**2)

    real(r8) :: perroad_a(bounds%begl:bounds%endl)           ! absorbed longwave for perroad (W/m**2)
    real(r8) :: perroad_r(bounds%begl:bounds%endl)           ! reflected longwave for perroad (W/m**2)
    real(r8) :: perroad_r_sky(bounds%begl:bounds%endl)       ! perroad_r to sky (W/m**2)
    real(r8) :: perroad_r_sunwall(bounds%begl:bounds%endl)   ! perroad_r to sunlit wall (W/m**2)
    real(r8) :: perroad_r_shadewall(bounds%begl:bounds%endl) ! perroad_r to shaded wall (W/m**2)
    real(r8) :: perroad_e(bounds%begl:bounds%endl)           ! emitted longwave for perroad (W/m**2)
    real(r8) :: perroad_e_sky(bounds%begl:bounds%endl)       ! perroad_e to sky (W/m**2)
    real(r8) :: perroad_e_sunwall(bounds%begl:bounds%endl)   ! perroad_e to sunlit wall (W/m**2)
    real(r8) :: perroad_e_shadewall(bounds%begl:bounds%endl) ! perroad_e to shaded wall (W/m**2)

    real(r8) :: road_a(bounds%begl:bounds%endl)              ! absorbed longwave for total road (W/m**2)
    real(r8) :: road_r(bounds%begl:bounds%endl)              ! reflected longwave for total road (W/m**2)
    real(r8) :: road_r_sky(bounds%begl:bounds%endl)          ! total road_r to sky (W/m**2)
    real(r8) :: road_r_sunwall(bounds%begl:bounds%endl)      ! total road_r to sunlit wall (W/m**2)
    real(r8) :: road_r_shadewall(bounds%begl:bounds%endl)    ! total road_r to shaded wall (W/m**2)
    real(r8) :: road_e(bounds%begl:bounds%endl)              ! emitted longwave for total road (W/m**2)
    real(r8) :: road_e_sky(bounds%begl:bounds%endl)          ! total road_e to sky (W/m**2)
    real(r8) :: road_e_sunwall(bounds%begl:bounds%endl)      ! total road_e to sunlit wall (W/m**2)
    real(r8) :: road_e_shadewall(bounds%begl:bounds%endl)    ! total road_e to shaded wall (W/m**2)

    real(r8) :: sunwall_a(bounds%begl:bounds%endl)           ! absorbed longwave (per unit wall area) for sunlit wall (W/m**2)
    real(r8) :: sunwall_r(bounds%begl:bounds%endl)           ! reflected longwave (per unit wall area) for sunlit wall (W/m**2)
    real(r8) :: sunwall_r_sky(bounds%begl:bounds%endl)       ! sunwall_r to sky (W/m**2)
    real(r8) :: sunwall_r_road(bounds%begl:bounds%endl)      ! sunwall_r to road (W/m**2)
    real(r8) :: sunwall_r_shadewall(bounds%begl:bounds%endl) ! sunwall_r to opposing (shaded) wall (W/m**2)
    real(r8) :: sunwall_e(bounds%begl:bounds%endl)           ! emitted longwave (per unit wall area) for sunlit wall (W/m**2)
    real(r8) :: sunwall_e_sky(bounds%begl:bounds%endl)       ! sunwall_e to sky (W/m**2)
    real(r8) :: sunwall_e_road(bounds%begl:bounds%endl)      ! sunwall_e to road (W/m**2)
    real(r8) :: sunwall_e_shadewall(bounds%begl:bounds%endl) ! sunwall_e to opposing (shaded) wall (W/m**2)

    real(r8) :: shadewall_a(bounds%begl:bounds%endl)         ! absorbed longwave (per unit wall area) for shaded wall (W/m**2)
    real(r8) :: shadewall_r(bounds%begl:bounds%endl)         ! reflected longwave (per unit wall area) for shaded wall (W/m**2)
    real(r8) :: shadewall_r_sky(bounds%begl:bounds%endl)     ! shadewall_r to sky (W/m**2)
    real(r8) :: shadewall_r_road(bounds%begl:bounds%endl)    ! shadewall_r to road (W/m**2)
    real(r8) :: shadewall_r_sunwall(bounds%begl:bounds%endl) ! shadewall_r to opposing (sunlit) wall (W/m**2)
    real(r8) :: shadewall_e(bounds%begl:bounds%endl)         ! emitted longwave (per unit wall area) for shaded wall (W/m**2)
    real(r8) :: shadewall_e_sky(bounds%begl:bounds%endl)     ! shadewall_e to sky (W/m**2)
    real(r8) :: shadewall_e_road(bounds%begl:bounds%endl)    ! shadewall_e to road (W/m**2)
    real(r8) :: shadewall_e_sunwall(bounds%begl:bounds%endl) ! shadewall_e to opposing (sunlit) wall (W/m**2)
    integer  :: l,fl,iter                    ! indices
    integer, parameter  :: n = 50            ! number of interations
    real(r8) :: crit                         ! convergence criterion (W/m**2)
    real(r8) :: err                          ! energy conservation error (W/m**2)
    real(r8) :: wtroad_imperv(bounds%begl:bounds%endl)       ! weight of impervious road wrt total road
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(canyon_hwr)      == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(wtroad_perv)     == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(lwdown)          == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(em_roof)         == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(em_improad)      == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(em_perroad)      == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(em_wall)         == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(t_roof)          == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(t_improad)       == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(t_perroad)       == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(t_sunwall)       == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(t_shadewall)     == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(lwnet_roof)      == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(lwnet_improad)   == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(lwnet_perroad)   == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(lwnet_sunwall)   == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(lwnet_shadewall) == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(lwnet_canyon)    == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(lwup_roof)       == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(lwup_improad)    == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(lwup_perroad)    == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(lwup_sunwall)    == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(lwup_shadewall)  == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(lwup_canyon)     == (/bounds%endl/)), errMsg(__FILE__, __LINE__))

   associate(& 
   vf_sr =>    lps%vf_sr , & ! Input:  [real(r8) (:)]  view factor of sky for road                       
   vf_wr =>    lps%vf_wr , & ! Input:  [real(r8) (:)]  view factor of one wall for road                  
   vf_sw =>    lps%vf_sw , & ! Input:  [real(r8) (:)]  view factor of sky for one wall                   
   vf_rw =>    lps%vf_rw , & ! Input:  [real(r8) (:)]  view factor of road for one wall                  
   vf_ww =>    lps%vf_ww   & ! Input:  [real(r8) (:)]  view factor of opposing wall for one wall         
   )
    ! Calculate impervious road

    do fl = 1,num_urbanl 
       l = filter_urbanl(fl)
       wtroad_imperv(l) = 1._r8 - wtroad_perv(l)
    end do

    do fl = 1,num_urbanl
       l = filter_urbanl(fl)
       ! atmospheric longwave radiation incident on walls and road in urban canyon.
       ! check for conservation (need to convert wall fluxes to ground area).
       ! lwdown (from atmosphere) = lwdown_road + (lwdown_sunwall + lwdown_shadewall)*canyon_hwr
       
       lwdown_road(l)      = lwdown(l) * vf_sr(l) 
       lwdown_sunwall(l)   = lwdown(l) * vf_sw(l)
       lwdown_shadewall(l) = lwdown(l) * vf_sw(l) 

       err = lwdown(l) - (lwdown_road(l) + (lwdown_shadewall(l) + lwdown_sunwall(l))*canyon_hwr(l))
       if (abs(err) > 0.10_r8 ) then
          write(iulog,*) 'urban incident atmospheric longwave radiation balance error',err
          write(iulog,*) 'l          = ',l
          write(iulog,*) 'lwdown     = ',lwdown(l)
          write(iulog,*) 'vf_sr      = ',vf_sr(l)
          write(iulog,*) 'vf_sw      = ',vf_sw(l)
          write(iulog,*) 'canyon_hwr = ',canyon_hwr(l)
          write(iulog,*) 'clm model is stopping'
          call endrun(decomp_index=l, clmlevel=namel, msg=errmsg(__FILE__, __LINE__))
       endif
    end do

    do fl = 1,num_urbanl
       l = filter_urbanl(fl)

       ! initial absorption, reflection, and emission for road and both walls. 
       ! distribute reflected and emitted radiation to sky, road, and walls according 
       ! to appropriate view factor. radiation reflected to road and walls will
       ! undergo multiple reflections within the canyon.

       road_a(l)              = 0.0_r8
       road_r(l)              = 0.0_r8
       road_e(l)              = 0.0_r8
       improad_a(l)           =     em_improad(l)  * lwdown_road(l) 
       improad_r(l)           = (1._r8-em_improad(l)) * lwdown_road(l) 
       improad_r_sky(l)       = improad_r(l) * vf_sr(l)
       improad_r_sunwall(l)   = improad_r(l) * vf_wr(l)
       improad_r_shadewall(l) = improad_r(l) * vf_wr(l)
       improad_e(l)           = em_improad(l) * sb * (t_improad(l)**4) 
       improad_e_sky(l)       = improad_e(l) * vf_sr(l)
       improad_e_sunwall(l)   = improad_e(l) * vf_wr(l)
       improad_e_shadewall(l) = improad_e(l) * vf_wr(l)
       road_a(l)              = road_a(l) + improad_a(l)*wtroad_imperv(l)
       road_r(l)              = road_r(l) + improad_r(l)*wtroad_imperv(l)
       road_e(l)              = road_e(l) + improad_e(l)*wtroad_imperv(l)

       perroad_a(l)           =     em_perroad(l)  * lwdown_road(l)
       perroad_r(l)           = (1._r8-em_perroad(l)) * lwdown_road(l)
       perroad_r_sky(l)       = perroad_r(l) * vf_sr(l)
       perroad_r_sunwall(l)   = perroad_r(l) * vf_wr(l)
       perroad_r_shadewall(l) = perroad_r(l) * vf_wr(l)
       perroad_e(l)           = em_perroad(l) * sb * (t_perroad(l)**4) 
       perroad_e_sky(l)       = perroad_e(l) * vf_sr(l)
       perroad_e_sunwall(l)   = perroad_e(l) * vf_wr(l)
       perroad_e_shadewall(l) = perroad_e(l) * vf_wr(l)
       road_a(l)              = road_a(l) + perroad_a(l)*wtroad_perv(l)
       road_r(l)              = road_r(l) + perroad_r(l)*wtroad_perv(l)
       road_e(l)              = road_e(l) + perroad_e(l)*wtroad_perv(l)

       road_r_sky(l)          = road_r(l) * vf_sr(l)
       road_r_sunwall(l)      = road_r(l) * vf_wr(l)
       road_r_shadewall(l)    = road_r(l) * vf_wr(l)
       road_e_sky(l)          = road_e(l) * vf_sr(l)
       road_e_sunwall(l)      = road_e(l) * vf_wr(l)
       road_e_shadewall(l)    = road_e(l) * vf_wr(l)

       sunwall_a(l)           = em_wall(l) * lwdown_sunwall(l)
       sunwall_r(l)           = (1._r8-em_wall(l)) * lwdown_sunwall(l)
       sunwall_r_sky(l)       = sunwall_r(l) * vf_sw(l)
       sunwall_r_road(l)      = sunwall_r(l) * vf_rw(l)
       sunwall_r_shadewall(l) = sunwall_r(l) * vf_ww(l)
       sunwall_e(l)           = em_wall(l) * sb * (t_sunwall(l)**4) 
       sunwall_e_sky(l)       = sunwall_e(l) * vf_sw(l)
       sunwall_e_road(l)      = sunwall_e(l) * vf_rw(l)
       sunwall_e_shadewall(l) = sunwall_e(l) * vf_ww(l)

       shadewall_a(l)         = em_wall(l) * lwdown_shadewall(l)
       shadewall_r(l)         = (1._r8-em_wall(l)) * lwdown_shadewall(l)
       shadewall_r_sky(l)     = shadewall_r(l) * vf_sw(l)
       shadewall_r_road(l)    = shadewall_r(l) * vf_rw(l)
       shadewall_r_sunwall(l) = shadewall_r(l) * vf_ww(l)
       shadewall_e(l)         = em_wall(l) * sb * (t_shadewall(l)**4) 
       shadewall_e_sky(l)     = shadewall_e(l) * vf_sw(l)
       shadewall_e_road(l)    = shadewall_e(l) * vf_rw(l)
       shadewall_e_sunwall(l) = shadewall_e(l) * vf_ww(l)

       ! initialize sum of net and upward longwave radiation for road and both walls

       lwnet_improad(l)   = improad_e(l)   - improad_a(l)
       lwnet_perroad(l)   = perroad_e(l)   - perroad_a(l)
       lwnet_sunwall(l)   = sunwall_e(l)   - sunwall_a(l)
       lwnet_shadewall(l) = shadewall_e(l) - shadewall_a(l)

       lwup_improad(l)   = improad_r_sky(l)   + improad_e_sky(l)
       lwup_perroad(l)   = perroad_r_sky(l)   + perroad_e_sky(l)
       lwup_sunwall(l)   = sunwall_r_sky(l)   + sunwall_e_sky(l)
       lwup_shadewall(l) = shadewall_r_sky(l) + shadewall_e_sky(l)

    end do

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

    do fl = 1,num_urbanl
       l = filter_urbanl(fl)

       do iter = 1, n
          ! step (1)

          lwtot(l) =  (sunwall_r_road(l) + sunwall_e_road(l)  &
                       + shadewall_r_road(l) + shadewall_e_road(l))*canyon_hwr(l)
          road_a(l)    = 0.0_r8
          road_r(l)    = 0.0_r8
          improad_r(l) = (1._r8-em_improad(l)) * lwtot(l)
          improad_a(l) =     em_improad(l)  * lwtot(l)
          road_a(l)    = road_a(l) + improad_a(l)*wtroad_imperv(l)
          road_r(l)    = road_r(l) + improad_r(l)*wtroad_imperv(l)
          perroad_r(l) = (1._r8-em_perroad(l)) * lwtot(l) 
          perroad_a(l) =     em_perroad(l)  * lwtot(l) 
          road_a(l)    = road_a(l) + perroad_a(l)*wtroad_perv(l)
          road_r(l)    = road_r(l) + perroad_r(l)*wtroad_perv(l)

          lwtot(l) = (road_r_sunwall(l) + road_e_sunwall(l))/canyon_hwr(l) &
                      + (shadewall_r_sunwall(l) + shadewall_e_sunwall(l))
          sunwall_a(l) =     em_wall(l)  * lwtot(l)
          sunwall_r(l) = (1._r8-em_wall(l)) * lwtot(l)

          lwtot(l) = (road_r_shadewall(l) + road_e_shadewall(l))/canyon_hwr(l) &
                      + (sunwall_r_shadewall(l) + sunwall_e_shadewall(l))
          shadewall_a(l) =     em_wall(l)  * lwtot(l)
          shadewall_r(l) = (1._r8-em_wall(l)) * lwtot(l)

          sunwall_e_road(l)      = 0._r8
          shadewall_e_road(l)    = 0._r8
          road_e_sunwall(l)      = 0._r8
          shadewall_e_sunwall(l) = 0._r8
          road_e_shadewall(l)    = 0._r8
          sunwall_e_shadewall(l) = 0._r8

          ! step (2)

          lwnet_improad(l)   = lwnet_improad(l)   - improad_a(l)
          lwnet_perroad(l)   = lwnet_perroad(l)   - perroad_a(l)
          lwnet_sunwall(l)   = lwnet_sunwall(l)   - sunwall_a(l)
          lwnet_shadewall(l) = lwnet_shadewall(l) - shadewall_a(l)

          ! step (3)

          improad_r_sky(l)       = improad_r(l) * vf_sr(l)
          improad_r_sunwall(l)   = improad_r(l) * vf_wr(l)
          improad_r_shadewall(l) = improad_r(l) * vf_wr(l)

          perroad_r_sky(l)       = perroad_r(l) * vf_sr(l)
          perroad_r_sunwall(l)   = perroad_r(l) * vf_wr(l)
          perroad_r_shadewall(l) = perroad_r(l) * vf_wr(l)

          road_r_sky(l)          = road_r(l) * vf_sr(l)
          road_r_sunwall(l)      = road_r(l) * vf_wr(l)
          road_r_shadewall(l)    = road_r(l) * vf_wr(l)

          sunwall_r_sky(l)       = sunwall_r(l) * vf_sw(l)
          sunwall_r_road(l)      = sunwall_r(l) * vf_rw(l)
          sunwall_r_shadewall(l) = sunwall_r(l) * vf_ww(l)

          shadewall_r_sky(l)     = shadewall_r(l) * vf_sw(l)
          shadewall_r_road(l)    = shadewall_r(l) * vf_rw(l)
          shadewall_r_sunwall(l) = shadewall_r(l) * vf_ww(l)

          ! step (4)

          lwup_improad(l)   = lwup_improad(l)   + improad_r_sky(l)
          lwup_perroad(l)   = lwup_perroad(l)   + perroad_r_sky(l)
          lwup_sunwall(l)   = lwup_sunwall(l)   + sunwall_r_sky(l)
          lwup_shadewall(l) = lwup_shadewall(l) + shadewall_r_sky(l)

          ! step (5)

          crit = max(road_a(l), sunwall_a(l), shadewall_a(l))
          if (crit < .001_r8) exit
       end do
       if (iter >= n) then
          write (iulog,*) 'urban net longwave radiation error: no convergence'
          write (iulog,*) 'clm model is stopping'
          call endrun(decomp_index=l, clmlevel=namel, msg=errmsg(__FILE__, __LINE__))
       endif

       ! total net longwave radiation for canyon. project wall fluxes to horizontal surface

       lwnet_canyon(l) = 0.0_r8
       lwnet_canyon(l) = lwnet_canyon(l) + lwnet_improad(l)*wtroad_imperv(l)
       lwnet_canyon(l) = lwnet_canyon(l) + lwnet_perroad(l)*wtroad_perv(l)
       lwnet_canyon(l) = lwnet_canyon(l) + (lwnet_sunwall(l) + lwnet_shadewall(l))*canyon_hwr(l)

       ! total emitted longwave for canyon. project wall fluxes to horizontal

       lwup_canyon(l) = 0.0_r8
       lwup_canyon(l) = lwup_canyon(l) + lwup_improad(l)*wtroad_imperv(l)
       lwup_canyon(l) = lwup_canyon(l) + lwup_perroad(l)*wtroad_perv(l)
       lwup_canyon(l) = lwup_canyon(l) + (lwup_sunwall(l) + lwup_shadewall(l))*canyon_hwr(l)

       ! conservation check. note: previous conservation check confirms partioning of incident
       ! atmospheric longwave radiation to road and walls is conserved as
       ! lwdown (from atmosphere) = lwdown_improad + lwdown_perroad + (lwdown_sunwall + lwdown_shadewall)*canyon_hwr

       err = lwnet_canyon(l) - (lwup_canyon(l) - lwdown(l))
       if (abs(err) > .10_r8 ) then
          write (iulog,*) 'urban net longwave radiation balance error',err
          write (iulog,*) 'clm model is stopping'
          call endrun(decomp_index=l, clmlevel=namel, msg=errmsg(__FILE__, __LINE__))
       end if

    end do

    ! Net longwave radiation for roof
    
    do fl = 1,num_urbanl
       l = filter_urbanl(fl)
       lwup_roof(l) = em_roof(l)*sb*(t_roof(l)**4) + (1._r8-em_roof(l))*lwdown(l)
       lwnet_roof(l) = lwup_roof(l) - lwdown(l)
    end do

    end associate 
  end subroutine net_longwave

  !-----------------------------------------------------------------------
  subroutine UrbanParamInit(bounds)
    !
    ! !DESCRIPTION: 
    ! Initialize urban surface dataset variables
    !
    ! This should NOT be called from within a threaded region, because of its array allocations
    !
    ! !USES:
    use clmtype
    use shr_infnan_mod, only : nan => shr_infnan_nan, assignment(=)
    use clm_varcon    , only : isturb_MIN
    use UrbanInputMod , only : urbinp
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: ib,l,g                    ! indices
    integer  :: ier                       ! error status
    integer  :: dindx                     ! urban density type index
    !-----------------------------------------------------------------------

    associate(&
    begl => bounds%begl, &
    endl => bounds%endl  &
    )

    ! Allocate memory for urban parameter components

    allocate( urban_params%wind_hgt_canyon   (begl:endl),        &
              urban_params%em_roof           (begl:endl),        &
              urban_params%em_improad        (begl:endl),        &
              urban_params%em_perroad        (begl:endl),        &
              urban_params%em_wall           (begl:endl),        &
              urban_params%alb_roof_dir      (begl:endl,numrad), &
              urban_params%alb_roof_dif      (begl:endl,numrad), &        
              urban_params%alb_improad_dir   (begl:endl,numrad), &        
              urban_params%alb_perroad_dir   (begl:endl,numrad), &        
              urban_params%alb_improad_dif   (begl:endl,numrad), &        
              urban_params%alb_perroad_dif   (begl:endl,numrad), &        
              urban_params%alb_wall_dir      (begl:endl,numrad), &        
              urban_params%alb_wall_dif      (begl:endl,numrad), &
              stat=ier)
    if (ier /= 0) then
       write(iulog,*)'UrbanParamInit: allocation error for urban derived type'
       call endrun(msg=errmsg(__FILE__, __LINE__))
    endif

    urban_params%wind_hgt_canyon(begl:endl)   = nan
    urban_params%em_roof(begl:endl)           = nan
    urban_params%em_improad(begl:endl)        = nan
    urban_params%em_perroad(begl:endl)        = nan
    urban_params%em_wall(begl:endl)           = nan
    urban_params%alb_roof_dir(begl:endl,:)    = nan
    urban_params%alb_roof_dif(begl:endl,:)    = nan
    urban_params%alb_improad_dir(begl:endl,:) = nan
    urban_params%alb_perroad_dir(begl:endl,:) = nan
    urban_params%alb_improad_dif(begl:endl,:) = nan
    urban_params%alb_perroad_dif(begl:endl,:) = nan
    urban_params%alb_wall_dir(begl:endl,:)    = nan
    urban_params%alb_wall_dif(begl:endl,:)    = nan

    
    ! Set constants in derived type
    
    do l = begl, endl
       if (lun%urbpoi(l)) then
          g =lun%gridcell(l)
          dindx = lun%itype(l) - isturb_MIN + 1
          urban_params%wind_hgt_canyon(l) = urbinp%wind_hgt_canyon(g,dindx)
          do ib = 1,numrad
             urban_params%alb_roof_dir   (l,ib) = urbinp%alb_roof_dir   (g,dindx,ib)
             urban_params%alb_roof_dif   (l,ib) = urbinp%alb_roof_dif   (g,dindx,ib)
             urban_params%alb_improad_dir(l,ib) = urbinp%alb_improad_dir(g,dindx,ib)
             urban_params%alb_perroad_dir(l,ib) = urbinp%alb_perroad_dir(g,dindx,ib)
             urban_params%alb_improad_dif(l,ib) = urbinp%alb_improad_dif(g,dindx,ib)
             urban_params%alb_perroad_dif(l,ib) = urbinp%alb_perroad_dif(g,dindx,ib)
             urban_params%alb_wall_dir   (l,ib) = urbinp%alb_wall_dir   (g,dindx,ib)
             urban_params%alb_wall_dif   (l,ib) = urbinp%alb_wall_dif   (g,dindx,ib)
          end do
          urban_params%em_roof   (l) = urbinp%em_roof   (g,dindx)
          urban_params%em_improad(l) = urbinp%em_improad(g,dindx)
          urban_params%em_perroad(l) = urbinp%em_perroad(g,dindx)
          urban_params%em_wall   (l) = urbinp%em_wall   (g,dindx)
       end if
    end do

    end associate

  end subroutine UrbanParamInit

  !-----------------------------------------------------------------------
  subroutine UrbanFluxes (bounds, &
       num_nourbanl, filter_nourbanl, &
       num_urbanl, filter_urbanl, &
       num_urbanc, filter_urbanc, &
       num_urbanp, filter_urbanp)
    !
    ! !DESCRIPTION: 
    ! Turbulent and momentum fluxes from urban canyon (consisting of roof, sunwall, 
    ! shadewall, pervious and impervious road).
    
    ! !USES:
    use clmtype
    use clm_varcon         , only : cpair, vkc, spval, icol_roof, icol_sunwall, &
                                    icol_shadewall, icol_road_perv, icol_road_imperv, &
                                    grav, pondmx_urban, rpi, rgas, &
                                    ht_wasteheat_factor, ac_wasteheat_factor, &
                                    wasteheat_limit
    use filterMod          , only : filter
    use FrictionVelocityMod, only : FrictionVelocity, MoninObukIni
    use QSatMod            , only : QSat
    use clm_varpar         , only : maxpatch_urb, nlevurb, nlevgrnd
    use clm_time_manager   , only : get_curr_date, get_step_size, get_nstep
    use clm_atmlnd         , only : clm_a2l, a2l_not_downscaled_gcell

    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds    ! bounds
    integer , intent(in) :: num_nourbanl       ! number of non-urban landunits in clump
    integer , intent(in) :: filter_nourbanl(:) ! non-urban landunit filter
    integer , intent(in) :: num_urbanl         ! number of urban landunits in clump
    integer , intent(in) :: filter_urbanl(:)   ! urban landunit filter
    integer , intent(in) :: num_urbanc         ! number of urban columns in clump
    integer , intent(in) :: filter_urbanc(:)   ! urban column filter
    integer , intent(in) :: num_urbanp         ! number of urban pfts in clump
    integer , intent(in) :: filter_urbanp(:)   ! urban pft filter
    !
    ! !LOCAL VARIABLES:
    character(len=*), parameter :: sub="UrbanFluxes"
    integer  :: fp,fc,fl,f,p,c,l,g,j,pi,i     ! indices
                                 
    real(r8) :: canyontop_wind(bounds%begl:bounds%endl)    ! wind at canyon top (m/s) 
    real(r8) :: canyon_u_wind(bounds%begl:bounds%endl)     ! u-component of wind speed inside canyon (m/s)
    real(r8) :: canyon_wind(bounds%begl:bounds%endl)       ! net wind speed inside canyon (m/s)
    real(r8) :: canyon_resistance(bounds%begl:bounds%endl) ! resistance to heat and moisture transfer from canyon road/walls to canyon air (s/m)

    real(r8) :: ur(bounds%begl:bounds%endl)        ! wind speed at reference height (m/s)
    real(r8) :: ustar(bounds%begl:bounds%endl)     ! friction velocity (m/s)
    real(r8) :: ramu(bounds%begl:bounds%endl)      ! aerodynamic resistance (s/m)
    real(r8) :: rahu(bounds%begl:bounds%endl)      ! thermal resistance (s/m)
    real(r8) :: rawu(bounds%begl:bounds%endl)      ! moisture resistance (s/m)
    real(r8) :: temp1(bounds%begl:bounds%endl)     ! relation for potential temperature profile
    real(r8) :: temp12m(bounds%begl:bounds%endl)   ! relation for potential temperature profile applied at 2-m
    real(r8) :: temp2(bounds%begl:bounds%endl)     ! relation for specific humidity profile
    real(r8) :: temp22m(bounds%begl:bounds%endl)   ! relation for specific humidity profile applied at 2-m
    real(r8) :: thm_g(bounds%begl:bounds%endl)     ! intermediate variable (forc_t+0.0098*forc_hgt_t)
    real(r8) :: thv_g(bounds%begl:bounds%endl)     ! virtual potential temperature (K)
    real(r8) :: dth(bounds%begl:bounds%endl)       ! diff of virtual temp. between ref. height and surface
    real(r8) :: dqh(bounds%begl:bounds%endl)       ! diff of humidity between ref. height and surface
    real(r8) :: zldis(bounds%begl:bounds%endl)     ! reference height "minus" zero displacement height (m)
    real(r8) :: um(bounds%begl:bounds%endl)        ! wind speed including the stablity effect (m/s)
    real(r8) :: obu(bounds%begl:bounds%endl)       ! Monin-Obukhov length (m)
    real(r8) :: taf_numer(bounds%begl:bounds%endl) ! numerator of taf equation (K m/s)
    real(r8) :: taf_denom(bounds%begl:bounds%endl) ! denominator of taf equation (m/s)
    real(r8) :: qaf_numer(bounds%begl:bounds%endl) ! numerator of qaf equation (kg m/kg s)
    real(r8) :: qaf_denom(bounds%begl:bounds%endl) ! denominator of qaf equation (m/s)
    real(r8) :: wtas(bounds%begl:bounds%endl)      ! sensible heat conductance for urban air to atmospheric air (m/s)
    real(r8) :: wtaq(bounds%begl:bounds%endl)      ! latent heat conductance for urban air to atmospheric air (m/s)
    real(r8) :: wts_sum(bounds%begl:bounds%endl)   ! sum of wtas, wtus_roof, wtus_road_perv, wtus_road_imperv, wtus_sunwall, wtus_shadewall
    real(r8) :: wtq_sum(bounds%begl:bounds%endl)   ! sum of wtaq, wtuq_roof, wtuq_road_perv, wtuq_road_imperv, wtuq_sunwall, wtuq_shadewall
    real(r8) :: beta(bounds%begl:bounds%endl)      ! coefficient of convective velocity
    real(r8) :: zii(bounds%begl:bounds%endl)       ! convective boundary layer height (m)

    real(r8) :: fm(bounds%begl:bounds%endl)        ! needed for BGC only to diagnose 10m wind speed

    real(r8) :: wtus(bounds%begc:bounds%endc)      ! sensible heat conductance for urban columns (scaled) (m/s)
    real(r8) :: wtuq(bounds%begc:bounds%endc)      ! latent heat conductance for urban columns (scaled) (m/s)

    integer  :: iter               ! iteration index
    real(r8) :: dthv               ! diff of vir. poten. temp. between ref. height and surface
    real(r8) :: tstar              ! temperature scaling parameter
    real(r8) :: qstar              ! moisture scaling parameter
    real(r8) :: thvstar            ! virtual potential temperature scaling parameter
    real(r8) :: wtus_roof(bounds%begl:bounds%endl)                ! sensible heat conductance for roof (scaled) (m/s)
    real(r8) :: wtuq_roof(bounds%begl:bounds%endl)                ! latent heat conductance for roof (scaled) (m/s)
    real(r8) :: wtus_road_perv(bounds%begl:bounds%endl)           ! sensible heat conductance for pervious road (scaled) (m/s)
    real(r8) :: wtuq_road_perv(bounds%begl:bounds%endl)           ! latent heat conductance for pervious road (scaled) (m/s)
    real(r8) :: wtus_road_imperv(bounds%begl:bounds%endl)         ! sensible heat conductance for impervious road (scaled) (m/s)
    real(r8) :: wtuq_road_imperv(bounds%begl:bounds%endl)         ! latent heat conductance for impervious road (scaled) (m/s)
    real(r8) :: wtus_sunwall(bounds%begl:bounds%endl)             ! sensible heat conductance for sunwall (scaled) (m/s)
    real(r8) :: wtuq_sunwall(bounds%begl:bounds%endl)             ! latent heat conductance for sunwall (scaled) (m/s)
    real(r8) :: wtus_shadewall(bounds%begl:bounds%endl)           ! sensible heat conductance for shadewall (scaled) (m/s)
    real(r8) :: wtuq_shadewall(bounds%begl:bounds%endl)           ! latent heat conductance for shadewall (scaled) (m/s)
    real(r8) :: wtus_roof_unscl(bounds%begl:bounds%endl)          ! sensible heat conductance for roof (not scaled) (m/s)
    real(r8) :: wtuq_roof_unscl(bounds%begl:bounds%endl)          ! latent heat conductance for roof (not scaled) (m/s)
    real(r8) :: wtus_road_perv_unscl(bounds%begl:bounds%endl)     ! sensible heat conductance for pervious road (not scaled) (m/s)
    real(r8) :: wtuq_road_perv_unscl(bounds%begl:bounds%endl)     ! latent heat conductance for pervious road (not scaled) (m/s)
    real(r8) :: wtus_road_imperv_unscl(bounds%begl:bounds%endl)   ! sensible heat conductance for impervious road (not scaled) (m/s)
    real(r8) :: wtuq_road_imperv_unscl(bounds%begl:bounds%endl)   ! latent heat conductance for impervious road (not scaled) (m/s)
    real(r8) :: wtus_sunwall_unscl(bounds%begl:bounds%endl)       ! sensible heat conductance for sunwall (not scaled) (m/s)
    real(r8) :: wtuq_sunwall_unscl(bounds%begl:bounds%endl)       ! latent heat conductance for sunwall (not scaled) (m/s)
    real(r8) :: wtus_shadewall_unscl(bounds%begl:bounds%endl)     ! sensible heat conductance for shadewall (not scaled) (m/s)
    real(r8) :: wtuq_shadewall_unscl(bounds%begl:bounds%endl)     ! latent heat conductance for shadewall (not scaled) (m/s)
    real(r8) :: t_sunwall_innerl(bounds%begl:bounds%endl)         ! temperature of inner layer of sunwall (K)
    real(r8) :: t_shadewall_innerl(bounds%begl:bounds%endl)       ! temperature of inner layer of shadewall (K)
    real(r8) :: t_roof_innerl(bounds%begl:bounds%endl)            ! temperature of inner layer of roof (K)
    real(r8) :: lngth_roof                                        ! length of roof (m)
    real(r8) :: wc                                                ! convective velocity (m/s)
    real(r8) :: zeta                                              ! dimensionless height used in Monin-Obukhov theory
    real(r8) :: eflx_sh_grnd_scale(bounds%begp:bounds%endp)       ! scaled sensible heat flux from ground (W/m**2) [+ to atm] 
    real(r8) :: qflx_evap_soi_scale(bounds%begp:bounds%endp)      ! scaled soil evaporation (mm H2O/s) (+ = to atm) 
    real(r8) :: eflx_wasteheat_roof(bounds%begl:bounds%endl)      ! sensible heat flux from urban heating/cooling sources of waste heat for roof (W/m**2)
    real(r8) :: eflx_wasteheat_sunwall(bounds%begl:bounds%endl)   ! sensible heat flux from urban heating/cooling sources of waste heat for sunwall (W/m**2)
    real(r8) :: eflx_wasteheat_shadewall(bounds%begl:bounds%endl) ! sensible heat flux from urban heating/cooling sources of waste heat for shadewall (W/m**2)
    real(r8) :: eflx_heat_from_ac_roof(bounds%begl:bounds%endl)   ! sensible heat flux put back into canyon due to heat removal by AC for roof (W/m**2)
    real(r8) :: eflx_heat_from_ac_sunwall(bounds%begl:bounds%endl)! sensible heat flux put back into canyon due to heat removal by AC for sunwall (W/m**2)
    real(r8) :: eflx_heat_from_ac_shadewall(bounds%begl:bounds%endl) ! sensible heat flux put back into canyon due to heat removal by AC for shadewall (W/m**2)
    real(r8) :: eflx(bounds%begl:bounds%endl)                     ! total sensible heat flux for error check (W/m**2)
    real(r8) :: qflx(bounds%begl:bounds%endl)                     ! total water vapor flux for error check (kg/m**2/s)
    real(r8) :: eflx_scale(bounds%begl:bounds%endl)               ! sum of scaled sensible heat fluxes for urban columns for error check (W/m**2)
    real(r8) :: qflx_scale(bounds%begl:bounds%endl)               ! sum of scaled water vapor fluxes for urban columns for error check (kg/m**2/s)
    real(r8) :: eflx_err(bounds%begl:bounds%endl)                 ! sensible heat flux error (W/m**2)
    real(r8) :: qflx_err(bounds%begl:bounds%endl)                 ! water vapor flux error (kg/m**2/s)
    real(r8) :: fwet_roof                                         ! fraction of roof surface that is wet (-)
    real(r8) :: fwet_road_imperv                                  ! fraction of impervious road surface that is wet (-)

    integer, parameter  :: niters = 3  ! maximum number of iterations for surface temperature
    integer  :: local_secp1(bounds%begl:bounds%endl)   ! seconds into current date in local time (sec)
    real(r8) :: dtime                  ! land model time step (sec)
    integer  :: year,month,day,secs    ! calendar info for current time step
    logical  :: found                  ! flag in search loop
    integer  :: indexl                 ! index of first found in search loop
    integer  :: nstep                  ! time step number
    real(r8), parameter :: lapse_rate = 0.0098_r8     ! Dry adiabatic lapse rate (K/m)
    real(r8) :: e_ref2m                ! 2 m height surface saturated vapor pressure [Pa]
    real(r8) :: de2mdT                 ! derivative of 2 m height surface saturated vapor pressure on t_ref2m
    real(r8) :: qsat_ref2m             ! 2 m height surface saturated specific humidity [kg/kg]
    real(r8) :: dqsat2mdT              ! derivative of 2 m height surface saturated specific humidity on t_ref2m

!-----------------------------------------------------------------------

   associate(& 
   ht_roof                             =>    lun%ht_roof                                 , & ! Input:  [real(r8) (:)]  height of urban roof (m)                          
   wtlunit_roof                        =>    lun%wtlunit_roof                            , & ! Input:  [real(r8) (:)]  weight of roof with respect to landunit           
   canyon_hwr                          =>    lun%canyon_hwr                              , & ! Input:  [real(r8) (:)]  ratio of building height to street width          
   wtroad_perv                         =>    lun%wtroad_perv                             , & ! Input:  [real(r8) (:)]  weight of pervious road wrt total road            
   wind_hgt_canyon                     =>    urban_params%wind_hgt_canyon                , & ! Input:  [real(r8) (:)]  height above road at which wind in canyon is to be computed (m)
   forc_t                              =>    a2l_not_downscaled_gcell%forc_t               , & ! Input:  [real(r8) (:)]  atmospheric temperature (K)                       
   forc_th                             =>    a2l_not_downscaled_gcell%forc_th              , & ! Input:  [real(r8) (:)]  atmospheric potential temperature (K)             
   forc_u                              =>    clm_a2l%forc_u                              , & ! Input:  [real(r8) (:)]  atmospheric wind speed in east direction (m/s)    
   forc_v                              =>    clm_a2l%forc_v                              , & ! Input:  [real(r8) (:)]  atmospheric wind speed in north direction (m/s)   
   forc_rho                            =>    a2l_not_downscaled_gcell%forc_rho             , & ! Input:  [real(r8) (:)]  density (kg/m**3)                                 
   forc_q                              =>    a2l_not_downscaled_gcell%forc_q               , & ! Input:  [real(r8) (:)]  atmospheric specific humidity (kg/kg)             
   forc_pbot                           =>    a2l_not_downscaled_gcell%forc_pbot            , & ! Input:  [real(r8) (:)]  atmospheric pressure (Pa)                         
   londeg                              =>    grc%londeg                                  , & ! Input:  [real(r8) (:)]  longitude (degrees)                               
   pfti                                =>   lun%pfti                                     , & ! Input:  [integer (:)]  beginning pft index for landunit                   
   pftf                                =>   lun%pftf                                     , & ! Input:  [integer (:)]  ending pft index for landunit                      
   coli                                =>   lun%coli                                     , & ! Input:  [integer (:)]  beginning column index for landunit                
   colf                                =>   lun%colf                                     , & ! Input:  [integer (:)]  ending column index for landunit                   
   lgridcell                           =>   lun%gridcell                                 , & ! Input:  [integer (:)]  gridcell of corresponding landunit                 
   z_0_town                            =>   lun%z_0_town                                 , & ! Input:  [real(r8) (:)]  momentum roughness length of urban landunit (m)   
   z_d_town                            =>   lun%z_d_town                                 , & ! Input:  [real(r8) (:)]  displacement height of urban landunit (m)         
   taf                                 =>    lps%taf                                     , & ! Input:  [real(r8) (:)]  urban canopy air temperature (K)                  
   qaf                                 =>    lps%qaf                                     , & ! Input:  [real(r8) (:)]  urban canopy air specific humidity (kg/kg)        
   npfts                               =>   lun%npfts                                    , & ! Input:  [integer (:)]  landunit's number of pfts (columns)                
   eflx_traffic                        =>    lef%eflx_traffic                            , & ! Input:  [real(r8) (:)]  traffic sensible heat flux (W/m**2)               
   eflx_traffic_factor                 =>    lef%eflx_traffic_factor                     , & ! Input:  [real(r8) (:)]  multiplicative urban traffic factor for sensible heat flux
   eflx_wasteheat                      =>    lef%eflx_wasteheat                          , & ! Input:  [real(r8) (:)]  sensible heat flux from urban heating/cooling sources of waste heat (W/m**2)
   eflx_heat_from_ac                   =>    lef%eflx_heat_from_ac                       , & ! Input:  [real(r8) (:)]  sensible heat flux put back into canyon due to removal by AC (W/m**2)
   t_building                          =>    lps%t_building                              , & ! Output: [real(r8) (:)]  internal building temperature (K)                 
   ctype                               =>    col%itype                                   , & ! Input:  [integer (:)]  column type                                        
   t_grnd                              =>    ces%t_grnd                                  , & ! Input:  [real(r8) (:)]  ground surface temperature (K)                    
   qg                                  =>    cws%qg                                      , & ! Input:  [real(r8) (:)]  specific humidity at ground surface (kg/kg)       
   htvp                                =>    cps%htvp                                    , & ! Input:  [real(r8) (:)]  latent heat of evaporation (/sublimation) (J/kg)  
   dqgdT                               =>    cws%dqgdT                                   , & ! Input:  [real(r8) (:)]  temperature derivative of "qg"                    
   t_soisno                            =>    ces%t_soisno                                , & ! Input:  [real(r8) (:,:)]  soil temperature (K)                            
   eflx_urban_ac                       =>    cef%eflx_urban_ac                           , & ! Input:  [real(r8) (:)]  urban air conditioning flux (W/m**2)              
   eflx_urban_heat                     =>    cef%eflx_urban_heat                         , & ! Input:  [real(r8) (:)]  urban heating flux (W/m**2)                       
   h2osoi_ice                          =>    cws%h2osoi_ice                              , & ! Input:  [real(r8) (:,:)]  ice lens (kg/m2)                                
   h2osoi_liq                          =>    cws%h2osoi_liq                              , & ! Input:  [real(r8) (:,:)]  liquid water (kg/m2)                            
   frac_sno                            =>    cps%frac_sno                                , & ! Input:  [real(r8) (:)]  fraction of ground covered by snow (0 to 1)       
   snow_depth                          =>    cps%snow_depth                              , & ! Input:  [real(r8) (:)]  snow height (m)                                   
   h2osno                              =>    cws%h2osno                                  , & ! Input:  [real(r8) (:)]  snow water (mm H2O)                               
   snl                                 =>    cps%snl                                     , & ! Input:  [integer (:)]  number of snow layers                              
   rootr_road_perv                     =>    cps%rootr_road_perv                         , & ! Input:  [real(r8) (:,:)]  effective fraction of roots in each soil layer for urban pervious road
   soilalpha_u                         =>    cws%soilalpha_u                             , & ! Input:  [real(r8) (:)]  Urban factor that reduces ground saturated specific humidity (-)
   pgridcell                           =>   pft%gridcell                                 , & ! Input:  [integer (:)]  gridcell of corresponding pft                      
   pcolumn                             =>   pft%column                                   , & ! Input:  [integer (:)]  column of corresponding pft                        
   plandunit                           =>   pft%landunit                                 , & ! Input:  [integer (:)]  pft's landunit index                               
   ram1                                =>    pps%ram1                                    , & ! Output: [real(r8) (:)]  aerodynamical resistance (s/m)                    
   dlrad                               =>    pef%dlrad                                   , & ! Output: [real(r8) (:)]  downward longwave radiation below the canopy (W/m**2)
   ulrad                               =>    pef%ulrad                                   , & ! Output: [real(r8) (:)]  upward longwave radiation above the canopy (W/m**2)
   cgrnds                              =>    pef%cgrnds                                  , & ! Output: [real(r8) (:)]  deriv, of soil sensible heat flux wrt soil temp (W/m**2/K)
   cgrndl                              =>    pef%cgrndl                                  , & ! Output: [real(r8) (:)]  deriv of soil latent heat flux wrt soil temp (W/m**2/K)
   cgrnd                               =>    pef%cgrnd                                   , & ! Output: [real(r8) (:)]  deriv. of soil energy flux wrt to soil temp (W/m**2/K)
   taux                                =>    pmf%taux                                    , & ! Output: [real(r8) (:)]  wind (shear) stress: e-w (kg/m/s**2)              
   tauy                                =>    pmf%tauy                                    , & ! Output: [real(r8) (:)]  wind (shear) stress: n-s (kg/m/s**2)              
   eflx_sh_grnd                        =>    pef%eflx_sh_grnd                            , & ! Output: [real(r8) (:)]  sensible heat flux from ground (W/m**2) [+ to atm]
   eflx_sh_tot                         =>    pef%eflx_sh_tot                             , & ! Output: [real(r8) (:)]  total sensible heat flux (W/m**2) [+ to atm]      
   eflx_sh_tot_u                       =>    pef%eflx_sh_tot_u                           , & ! Output: [real(r8) (:)]  urban total sensible heat flux (W/m**2) [+ to atm]
   qflx_evap_soi                       =>    pwf%qflx_evap_soi                           , & ! Output: [real(r8) (:)]  soil evaporation (mm H2O/s) (+ = to atm)          
   qflx_tran_veg                       =>    pwf%qflx_tran_veg                           , & ! Output: [real(r8) (:)]  vegetation transpiration (mm H2O/s) (+ = to atm)  
   qflx_evap_veg                       =>    pwf%qflx_evap_veg                           , & ! Output: [real(r8) (:)]  vegetation evaporation (mm H2O/s) (+ = to atm)    
   qflx_evap_tot                       =>    pwf%qflx_evap_tot                           , & ! Output: [real(r8) (:)]  qflx_evap_soi + qflx_evap_can + qflx_tran_veg     
   t_ref2m                             =>    pes%t_ref2m                                 , & ! Output: [real(r8) (:)]  2 m height surface air temperature (K)            
   q_ref2m                             =>    pes%q_ref2m                                 , & ! Output: [real(r8) (:)]  2 m height surface specific humidity (kg/kg)      
   t_ref2m_u                           =>    pes%t_ref2m_u                               , & ! Output: [real(r8) (:)]  Urban 2 m height surface air temperature (K)      
   t_veg                               =>    pes%t_veg                                   , & ! Output: [real(r8) (:)]  vegetation temperature (K)                        
   rootr                               =>    pps%rootr                                   , & ! Output: [real(r8) (:,:)]  effective fraction of roots in each soil layer  
   psnsun                              =>    pcf%psnsun                                  , & ! Output: [real(r8) (:)]  sunlit leaf photosynthesis (umol CO2 /m**2/ s)    
   psnsha                              =>    pcf%psnsha                                  , & ! Output: [real(r8) (:)]  shaded leaf photosynthesis (umol CO2 /m**2/ s)    
   forc_hgt_u_pft                      =>    pps%forc_hgt_u_pft                          , & ! Input:  [real(r8) (:)]  observational height of wind at pft-level (m)     
   forc_hgt_t_pft                      =>    pps%forc_hgt_t_pft                          , & ! Input:  [real(r8) (:)]  observational height of temperature at pft-level (m)
   rh_ref2m                            =>    pes%rh_ref2m                                , & ! Output: [real(r8) (:)]  2 m height surface relative humidity (%)          
   rh_ref2m_u                          =>    pes%rh_ref2m_u                              , & ! Output: [real(r8) (:)]  Urban 2 m height surface relative humidity (%)    
   eflx_sh_snow                        =>    pef%eflx_sh_snow                            , & ! Output: [real(r8) (:)]  sensible heat flux from snow (W/m**2) [+ to atm]  
   eflx_sh_soil                        =>    pef%eflx_sh_soil                            , & ! Output: [real(r8) (:)]  sensible heat flux from soil (W/m**2) [+ to atm]  
   eflx_sh_h2osfc                      =>    pef%eflx_sh_h2osfc                          , & ! Output: [real(r8) (:)]  sensible heat flux from soil (W/m**2) [+ to atm]  
   begl                                =>    bounds%begl                                 , &
   endl                                =>    bounds%endl                                   &
   )
    ! Define fields that appear on the restart file for non-urban landunits 

    do fl = 1,num_nourbanl
       l = filter_nourbanl(fl)
       taf(l) = spval
       qaf(l) = spval
    end do

    ! Get time step
    nstep = get_nstep()

    ! Set constants (same as in Biogeophysics1Mod)
    beta(begl:endl) = 1._r8             ! Should be set to the same values as in Biogeophysics1Mod
    zii(begl:endl)  = 1000._r8          ! Should be set to the same values as in Biogeophysics1Mod

    ! Get current date
    dtime = get_step_size()
    call get_curr_date (year, month, day, secs)
    
    ! Compute canyontop wind using Masson (2000)

    do fl = 1, num_urbanl
       l = filter_urbanl(fl)
       g = lgridcell(l)

       local_secp1(l)        = secs + nint((grc%londeg(g)/degpsec)/dtime)*dtime
       local_secp1(l)        = mod(local_secp1(l),isecspday)

       ! Error checks

       if (ht_roof(l) - z_d_town(l) <= z_0_town(l)) then
          write (iulog,*) 'aerodynamic parameter error in UrbanFluxes'
          write (iulog,*) 'h_r - z_d <= z_0'
          write (iulog,*) 'ht_roof, z_d_town, z_0_town: ', ht_roof(l), z_d_town(l), &
                       z_0_town(l)
          write (iulog,*) 'clm model is stopping'
          call endrun(decomp_index=l, clmlevel=namel, msg=errmsg(__FILE__, __LINE__))
       end if
       if (forc_hgt_u_pft(pfti(l)) - z_d_town(l) <= z_0_town(l)) then
          write (iulog,*) 'aerodynamic parameter error in UrbanFluxes'
          write (iulog,*) 'h_u - z_d <= z_0'
          write (iulog,*) 'forc_hgt_u_pft, z_d_town, z_0_town: ', forc_hgt_u_pft(pfti(l)), z_d_town(l), &
                       z_0_town(l)
          write (iulog,*) 'clm model is stopping'
          call endrun(decomp_index=l, clmlevel=namel, msg=errmsg(__FILE__, __LINE__))
       end if

       ! Magnitude of atmospheric wind

       ur(l) = max(1.0_r8,sqrt(forc_u(g)*forc_u(g)+forc_v(g)*forc_v(g)))

       ! Canyon top wind

       canyontop_wind(l) = ur(l) * &
                           log( (ht_roof(l)-z_d_town(l)) / z_0_town(l) ) / &
                           log( (forc_hgt_u_pft(pfti(l))-z_d_town(l)) / z_0_town(l) )

       ! U component of canyon wind 

       if (canyon_hwr(l) < 0.5_r8) then  ! isolated roughness flow
         canyon_u_wind(l) = canyontop_wind(l) * exp( -0.5_r8*canyon_hwr(l)* &
                            (1._r8-(wind_hgt_canyon(l)/ht_roof(l))) )
       else if (canyon_hwr(l) < 1.0_r8) then ! wake interference flow
         canyon_u_wind(l) = canyontop_wind(l) * (1._r8+2._r8*(2._r8/rpi - 1._r8)* &
                            (ht_roof(l)/(ht_roof(l)/canyon_hwr(l)) - 0.5_r8)) * &
                            exp(-0.5_r8*canyon_hwr(l)*(1._r8-(wind_hgt_canyon(l)/ht_roof(l))))
       else  ! skimming flow
         canyon_u_wind(l) = canyontop_wind(l) * (2._r8/rpi) * &
                            exp(-0.5_r8*canyon_hwr(l)*(1._r8-(wind_hgt_canyon(l)/ht_roof(l))))
       end if

    end do

! Compute fluxes - Follows CLM approach for bare soils (Oleson et al 2004)

    do fl = 1, num_urbanl
       l = filter_urbanl(fl)
       g = lgridcell(l)

       thm_g(l) = forc_t(g) + lapse_rate*forc_hgt_t_pft(pfti(l))
       thv_g(l) = forc_th(g)*(1._r8+0.61_r8*forc_q(g))
       dth(l)   = thm_g(l)-taf(l)
       dqh(l)   = forc_q(g)-qaf(l)
       dthv     = dth(l)*(1._r8+0.61_r8*forc_q(g))+0.61_r8*forc_th(g)*dqh(l)
       zldis(l) = forc_hgt_u_pft(pfti(l)) - z_d_town(l)

       ! Initialize Monin-Obukhov length and wind speed including convective velocity

       call MoninObukIni(ur(l), thv_g(l), dthv, zldis(l), z_0_town(l), um(l), obu(l))

    end do

    ! Initialize conductances
    wtus_roof(begl:endl)        = 0._r8
    wtus_road_perv(begl:endl)   = 0._r8
    wtus_road_imperv(begl:endl) = 0._r8
    wtus_sunwall(begl:endl)     = 0._r8
    wtus_shadewall(begl:endl)   = 0._r8
    wtuq_roof(begl:endl)        = 0._r8
    wtuq_road_perv(begl:endl)   = 0._r8
    wtuq_road_imperv(begl:endl) = 0._r8
    wtuq_sunwall(begl:endl)     = 0._r8
    wtuq_shadewall(begl:endl)   = 0._r8
    wtus_roof_unscl(begl:endl)        = 0._r8
    wtus_road_perv_unscl(begl:endl)   = 0._r8
    wtus_road_imperv_unscl(begl:endl) = 0._r8
    wtus_sunwall_unscl(begl:endl)     = 0._r8
    wtus_shadewall_unscl(begl:endl)   = 0._r8
    wtuq_roof_unscl(begl:endl)        = 0._r8
    wtuq_road_perv_unscl(begl:endl)   = 0._r8
    wtuq_road_imperv_unscl(begl:endl) = 0._r8
    wtuq_sunwall_unscl(begl:endl)     = 0._r8
    wtuq_shadewall_unscl(begl:endl)   = 0._r8

    ! Start stability iteration

    do iter = 1,niters

     ! Get friction velocity, relation for potential
     ! temperature and humidity profiles of surface boundary layer.

      if (num_urbanl .gt. 0) then
         call FrictionVelocity(begl, endl, &
              num_urbanl, filter_urbanl, &
              z_d_town(begl:endl), z_0_town(begl:endl), z_0_town(begl:endl), z_0_town(begl:endl), &
              obu(begl:endl), iter, ur(begl:endl), um(begl:endl), ustar(begl:endl), &
              temp1(begl:endl), temp2(begl:endl), temp12m(begl:endl), temp22m(begl:endl), fm(begl:endl), &
              landunit_index=.true.)
      end if

      do fl = 1, num_urbanl
         l = filter_urbanl(fl)
         g = lgridcell(l)

         ! Determine aerodynamic resistance to fluxes from urban canopy air to
         ! atmosphere

         ramu(l) = 1._r8/(ustar(l)*ustar(l)/um(l))
         rahu(l) = 1._r8/(temp1(l)*ustar(l))
         rawu(l) = 1._r8/(temp2(l)*ustar(l))

         ! Determine magnitude of canyon wind by using horizontal wind determined
         ! previously and vertical wind from friction velocity (Masson 2000)

         canyon_wind(l) = sqrt(canyon_u_wind(l)**2._r8 + ustar(l)**2._r8)

         ! Determine canyon_resistance (currently this single resistance determines the
         ! resistance from urban surfaces (roof, pervious and impervious road, sunlit and
         ! shaded walls) to urban canopy air, since it is only dependent on wind speed
         ! Also from Masson 2000.

         canyon_resistance(l) = cpair * forc_rho(g) / (11.8_r8 + 4.2_r8*canyon_wind(l))

      end do

      ! This is the first term in the equation solutions for urban canopy air temperature
      ! and specific humidity (numerator) and is a landunit quantity
      do fl = 1, num_urbanl
         l = filter_urbanl(fl)
         g = lgridcell(l)

         taf_numer(l) = thm_g(l)/rahu(l)
         taf_denom(l) = 1._r8/rahu(l)
         qaf_numer(l) = forc_q(g)/rawu(l)
         qaf_denom(l) = 1._r8/rawu(l)

         ! First term needed for derivative of heat fluxes
         wtas(l) = 1._r8/rahu(l)
         wtaq(l) = 1._r8/rawu(l)

      end do


      ! Gather other terms for other urban columns for numerator and denominator of
      ! equations for urban canopy air temperature and specific humidity

      do fc = 1,num_urbanc
         c = filter_urbanc(fc)
         l = col%landunit(c)

         if (ctype(c) == icol_roof) then

            ! scaled sensible heat conductance
            wtus(c) = wtlunit_roof(l)/canyon_resistance(l)
            wtus_roof(l) = wtus(c)
            ! unscaled sensible heat conductance
            wtus_roof_unscl(l) = 1._r8/canyon_resistance(l)

            if (snow_depth(c) > 0._r8) then
               fwet_roof = min(snow_depth(c)/0.05_r8, 1._r8)
            else
               fwet_roof = (max(0._r8, h2osoi_liq(c,1)+h2osoi_ice(c,1))/pondmx_urban)**0.666666666666_r8
               fwet_roof = min(fwet_roof,1._r8)
            end if
            if (qaf(l) > qg(c)) then 
               fwet_roof = 1._r8
            end if
            ! scaled latent heat conductance
            wtuq(c) = fwet_roof*(wtlunit_roof(l)/canyon_resistance(l))
            wtuq_roof(l) = wtuq(c)
            ! unscaled latent heat conductance
            wtuq_roof_unscl(l) = fwet_roof*(1._r8/canyon_resistance(l))

            ! wasteheat from heating/cooling
            if (trim(urban_hac) == urban_wasteheat_on) then
               eflx_wasteheat_roof(l) = ac_wasteheat_factor * eflx_urban_ac(c) + &
                    ht_wasteheat_factor * eflx_urban_heat(c)
            else
               eflx_wasteheat_roof(l) = 0._r8
            end if

            ! If air conditioning on, always replace heat removed with heat into canyon
            if (trim(urban_hac) == urban_hac_on .or. trim(urban_hac) == urban_wasteheat_on) then
               eflx_heat_from_ac_roof(l) = abs(eflx_urban_ac(c))
            else
               eflx_heat_from_ac_roof(l) = 0._r8
            end if

         else if (ctype(c) == icol_road_perv) then

            ! scaled sensible heat conductance
            wtus(c) = wtroad_perv(l)*(1._r8-wtlunit_roof(l))/canyon_resistance(l)
            wtus_road_perv(l) = wtus(c)
            ! unscaled sensible heat conductance
            wtus_road_perv_unscl(l) = 1._r8/canyon_resistance(l)

            ! scaled latent heat conductance
            wtuq(c) = wtroad_perv(l)*(1._r8-wtlunit_roof(l))/canyon_resistance(l)
            wtuq_road_perv(l) = wtuq(c)
            ! unscaled latent heat conductance
            wtuq_road_perv_unscl(l) = 1._r8/canyon_resistance(l)

         else if (ctype(c) == icol_road_imperv) then

            ! scaled sensible heat conductance
            wtus(c) = (1._r8-wtroad_perv(l))*(1._r8-wtlunit_roof(l))/canyon_resistance(l)
            wtus_road_imperv(l) = wtus(c)
            ! unscaled sensible heat conductance
            wtus_road_imperv_unscl(l) = 1._r8/canyon_resistance(l)

            if (snow_depth(c) > 0._r8) then
               fwet_road_imperv = min(snow_depth(c)/0.05_r8, 1._r8)
            else
               fwet_road_imperv = (max(0._r8, h2osoi_liq(c,1)+h2osoi_ice(c,1))/pondmx_urban)**0.666666666666_r8
               fwet_road_imperv = min(fwet_road_imperv,1._r8)
            end if
            if (qaf(l) > qg(c)) then 
               fwet_road_imperv = 1._r8
            end if
            ! scaled latent heat conductance
            wtuq(c) = fwet_road_imperv*(1._r8-wtroad_perv(l))*(1._r8-wtlunit_roof(l))/canyon_resistance(l)
            wtuq_road_imperv(l) = wtuq(c)
            ! unscaled latent heat conductance
            wtuq_road_imperv_unscl(l) = fwet_road_imperv*(1._r8/canyon_resistance(l))

         else if (ctype(c) == icol_sunwall) then

            ! scaled sensible heat conductance
            wtus(c) = canyon_hwr(l)*(1._r8-wtlunit_roof(l))/canyon_resistance(l)
            wtus_sunwall(l) = wtus(c)
            ! unscaled sensible heat conductance
            wtus_sunwall_unscl(l) = 1._r8/canyon_resistance(l)

            ! scaled latent heat conductance
            wtuq(c) = 0._r8
            wtuq_sunwall(l) = wtuq(c)
            ! unscaled latent heat conductance
            wtuq_sunwall_unscl(l) = 0._r8

            ! wasteheat from heating/cooling
            if (trim(urban_hac) == urban_wasteheat_on) then
               eflx_wasteheat_sunwall(l) = ac_wasteheat_factor * eflx_urban_ac(c) + &
                    ht_wasteheat_factor * eflx_urban_heat(c)
            else
               eflx_wasteheat_sunwall(l) = 0._r8
            end if

            ! If air conditioning on, always replace heat removed with heat into canyon
            if (trim(urban_hac) == urban_hac_on .or. trim(urban_hac) == urban_wasteheat_on) then
               eflx_heat_from_ac_sunwall(l) = abs(eflx_urban_ac(c))
            else
               eflx_heat_from_ac_sunwall(l) = 0._r8
            end if

         else if (ctype(c) == icol_shadewall) then

            ! scaled sensible heat conductance
            wtus(c) = canyon_hwr(l)*(1._r8-wtlunit_roof(l))/canyon_resistance(l)
            wtus_shadewall(l) = wtus(c)
            ! unscaled sensible heat conductance
            wtus_shadewall_unscl(l) = 1._r8/canyon_resistance(l)

            ! scaled latent heat conductance
            wtuq(c) = 0._r8
            wtuq_shadewall(l) = wtuq(c)
            ! unscaled latent heat conductance
            wtuq_shadewall_unscl(l) = 0._r8

            ! wasteheat from heating/cooling
            if (trim(urban_hac) == urban_wasteheat_on) then
               eflx_wasteheat_shadewall(l) = ac_wasteheat_factor * eflx_urban_ac(c) + &
                    ht_wasteheat_factor * eflx_urban_heat(c)
            else
               eflx_wasteheat_shadewall(l) = 0._r8
            end if

            ! If air conditioning on, always replace heat removed with heat into canyon
            if (trim(urban_hac) == urban_hac_on .or. trim(urban_hac) == urban_wasteheat_on) then
               eflx_heat_from_ac_shadewall(l) = abs(eflx_urban_ac(c))
            else
               eflx_heat_from_ac_shadewall(l) = 0._r8
            end if
         else
            write(iulog,*) 'c, ctype, pi = ', c, ctype(c), pi
            write(iulog,*) 'Column indices for: shadewall, sunwall, road_imperv, road_perv, roof: '
            write(iulog,*) icol_shadewall, icol_sunwall, icol_road_imperv, icol_road_perv, icol_roof
            call endrun(decomp_index=l, clmlevel=namel, msg="ERROR, ctype out of range"//errmsg(__FILE__, __LINE__))
         end if

         taf_numer(l) = taf_numer(l) + t_grnd(c)*wtus(c)
         taf_denom(l) = taf_denom(l) + wtus(c)
         qaf_numer(l) = qaf_numer(l) + qg(c)*wtuq(c)
         qaf_denom(l) = qaf_denom(l) + wtuq(c)

      end do

      ! Calculate new urban canopy air temperature and specific humidity

      do fl = 1, num_urbanl
         l = filter_urbanl(fl)
         g = lgridcell(l)

         ! Total waste heat and heat from AC is sum of heat for walls and roofs
         ! accounting for different surface areas
         eflx_wasteheat(l) = wtlunit_roof(l)*eflx_wasteheat_roof(l) + &
            (1._r8-wtlunit_roof(l))*(canyon_hwr(l)*(eflx_wasteheat_sunwall(l) + &
            eflx_wasteheat_shadewall(l)))

         ! Limit wasteheat to ensure that we don't get any unrealistically strong 
         ! positive feedbacks due to AC in a warmer climate
         eflx_wasteheat(l) = min(eflx_wasteheat(l),wasteheat_limit)

         eflx_heat_from_ac(l) = wtlunit_roof(l)*eflx_heat_from_ac_roof(l) + &
            (1._r8-wtlunit_roof(l))*(canyon_hwr(l)*(eflx_heat_from_ac_sunwall(l) + &
            eflx_heat_from_ac_shadewall(l)))

         ! Calculate traffic heat flux
         ! Only comes from impervious road
         eflx_traffic(l) = (1._r8-wtlunit_roof(l))*(1._r8-wtroad_perv(l))* &
                           eflx_traffic_factor(l)

         taf(l) = taf_numer(l)/taf_denom(l)
         qaf(l) = qaf_numer(l)/qaf_denom(l)

         wts_sum(l) = wtas(l) + wtus_roof(l) + wtus_road_perv(l) + &
                      wtus_road_imperv(l) + wtus_sunwall(l) + wtus_shadewall(l)

         wtq_sum(l) = wtaq(l) + wtuq_roof(l) + wtuq_road_perv(l) + &
                      wtuq_road_imperv(l) + wtuq_sunwall(l) + wtuq_shadewall(l)

      end do

      ! This section of code is not required if niters = 1
      ! Determine stability using new taf and qaf
      ! TODO: Some of these constants replicate what is in FrictionVelocity and BareGround fluxes should consildate. EBK
      do fl = 1, num_urbanl
         l = filter_urbanl(fl)
         g = lgridcell(l)

         dth(l) = thm_g(l)-taf(l)
         dqh(l) = forc_q(g)-qaf(l)
         tstar = temp1(l)*dth(l)
         qstar = temp2(l)*dqh(l)
         thvstar = tstar*(1._r8+0.61_r8*forc_q(g)) + 0.61_r8*forc_th(g)*qstar
         zeta = zldis(l)*vkc*grav*thvstar/(ustar(l)**2*thv_g(l))

         if (zeta >= 0._r8) then                   !stable
            zeta = min(2._r8,max(zeta,0.01_r8))
            um(l) = max(ur(l),0.1_r8)
         else                                      !unstable
            zeta = max(-100._r8,min(zeta,-0.01_r8))
            wc = beta(l)*(-grav*ustar(l)*thvstar*zii(l)/thv_g(l))**0.333_r8
            um(l) = sqrt(ur(l)*ur(l) + wc*wc)
         end if

         obu(l) = zldis(l)/zeta
      end do

    end do   ! end iteration

! Determine fluxes from canyon surfaces

    ! the following initializations are needed to ensure that the values are 0 over non-
    ! active urban PFTs
    eflx_sh_grnd_scale(bounds%begp : bounds%endp) = 0._r8
    qflx_evap_soi_scale(bounds%begp : bounds%endp) = 0._r8

    do f = 1, num_urbanp

       p = filter_urbanp(f)
       c = pcolumn(p)
       g = pgridcell(p)
       l = plandunit(p)

       ram1(p) = ramu(l)  !pass value to global variable

       ! Upward and downward canopy longwave are zero

       ulrad(p)  = 0._r8
       dlrad(p)  = 0._r8

       ! Derivative of sensible and latent heat fluxes with respect to 
       ! ground temperature

       if (ctype(c) == icol_roof) then
         cgrnds(p) = forc_rho(g) * cpair * (wtas(l) + wtus_road_perv(l) +  &
                     wtus_road_imperv(l) + wtus_sunwall(l) + wtus_shadewall(l)) * &
                     (wtus_roof_unscl(l)/wts_sum(l))
         cgrndl(p) = forc_rho(g) * (wtaq(l) + wtuq_road_perv(l) +  &
                     wtuq_road_imperv(l) + wtuq_sunwall(l) + wtuq_shadewall(l)) * &
                     (wtuq_roof_unscl(l)/wtq_sum(l))*dqgdT(c)
       else if (ctype(c) == icol_road_perv) then
         cgrnds(p) = forc_rho(g) * cpair * (wtas(l) + wtus_roof(l) +  &
                     wtus_road_imperv(l) + wtus_sunwall(l) + wtus_shadewall(l)) * &
                     (wtus_road_perv_unscl(l)/wts_sum(l))
         cgrndl(p) = forc_rho(g) * (wtaq(l) + wtuq_roof(l) +  &
                     wtuq_road_imperv(l) + wtuq_sunwall(l) + wtuq_shadewall(l)) * &
                     (wtuq_road_perv_unscl(l)/wtq_sum(l))*dqgdT(c)
       else if (ctype(c) == icol_road_imperv) then
         cgrnds(p) = forc_rho(g) * cpair * (wtas(l) + wtus_roof(l) +  &
                     wtus_road_perv(l) + wtus_sunwall(l) + wtus_shadewall(l)) * &
                     (wtus_road_imperv_unscl(l)/wts_sum(l))
         cgrndl(p) = forc_rho(g) * (wtaq(l) + wtuq_roof(l) +  &
                     wtuq_road_perv(l) + wtuq_sunwall(l) + wtuq_shadewall(l)) * &
                     (wtuq_road_imperv_unscl(l)/wtq_sum(l))*dqgdT(c)
       else if (ctype(c) == icol_sunwall) then
         cgrnds(p) = forc_rho(g) * cpair * (wtas(l) + wtus_roof(l) +  &
                     wtus_road_perv(l) + wtus_road_imperv(l) + wtus_shadewall(l)) * &
                     (wtus_sunwall_unscl(l)/wts_sum(l))
         cgrndl(p) = 0._r8
       else if (ctype(c) == icol_shadewall) then
         cgrnds(p) = forc_rho(g) * cpair * (wtas(l) + wtus_roof(l) +  &
                     wtus_road_perv(l) + wtus_road_imperv(l) + wtus_sunwall(l)) * &
                     (wtus_shadewall_unscl(l)/wts_sum(l))
         cgrndl(p) = 0._r8
       end if
       cgrnd(p)  = cgrnds(p) + cgrndl(p)*htvp(c)

       ! Surface fluxes of momentum, sensible and latent heat

       taux(p)          = -forc_rho(g)*forc_u(g)/ramu(l)
       tauy(p)          = -forc_rho(g)*forc_v(g)/ramu(l)

       ! Use new canopy air temperature
       dth(l) = taf(l) - t_grnd(c)

       if (ctype(c) == icol_roof) then
         eflx_sh_grnd(p)  = -forc_rho(g)*cpair*wtus_roof_unscl(l)*dth(l)
         eflx_sh_snow(p)  = 0._r8
         eflx_sh_soil(p)  = 0._r8
         eflx_sh_h2osfc(p)= 0._r8
       else if (ctype(c) == icol_road_perv) then
         eflx_sh_grnd(p)  = -forc_rho(g)*cpair*wtus_road_perv_unscl(l)*dth(l)
         eflx_sh_snow(p)  = 0._r8
         eflx_sh_soil(p)  = 0._r8
         eflx_sh_h2osfc(p)= 0._r8
       else if (ctype(c) == icol_road_imperv) then
         eflx_sh_grnd(p)  = -forc_rho(g)*cpair*wtus_road_imperv_unscl(l)*dth(l)
         eflx_sh_snow(p)  = 0._r8
         eflx_sh_soil(p)  = 0._r8
         eflx_sh_h2osfc(p)= 0._r8
       else if (ctype(c) == icol_sunwall) then
         eflx_sh_grnd(p)  = -forc_rho(g)*cpair*wtus_sunwall_unscl(l)*dth(l)
         eflx_sh_snow(p)  = 0._r8
         eflx_sh_soil(p)  = 0._r8
         eflx_sh_h2osfc(p)= 0._r8
       else if (ctype(c) == icol_shadewall) then
         eflx_sh_grnd(p)  = -forc_rho(g)*cpair*wtus_shadewall_unscl(l)*dth(l)
         eflx_sh_snow(p)  = 0._r8
         eflx_sh_soil(p)  = 0._r8
         eflx_sh_h2osfc(p)= 0._r8
       end if

       eflx_sh_tot(p)   = eflx_sh_grnd(p)
       eflx_sh_tot_u(p) = eflx_sh_tot(p)

       dqh(l) = qaf(l) - qg(c)

       if (ctype(c) == icol_roof) then
         qflx_evap_soi(p) = -forc_rho(g)*wtuq_roof_unscl(l)*dqh(l)
       else if (ctype(c) == icol_road_perv) then
         ! Evaporation assigned to soil term if dew or snow
         ! or if no liquid water available in soil column
         if (dqh(l) > 0._r8 .or. frac_sno(c) > 0._r8 .or. soilalpha_u(c) .le. 0._r8) then
            qflx_evap_soi(p) = -forc_rho(g)*wtuq_road_perv_unscl(l)*dqh(l)
            qflx_tran_veg(p) = 0._r8
         ! Otherwise, evaporation assigned to transpiration term
         else
            qflx_evap_soi(p) = 0._r8
            qflx_tran_veg(p) = -forc_rho(g)*wtuq_road_perv_unscl(l)*dqh(l)
         end if
         qflx_evap_veg(p) = qflx_tran_veg(p)
       else if (ctype(c) == icol_road_imperv) then
         qflx_evap_soi(p) = -forc_rho(g)*wtuq_road_imperv_unscl(l)*dqh(l)
       else if (ctype(c) == icol_sunwall) then
         qflx_evap_soi(p) = 0._r8
       else if (ctype(c) == icol_shadewall) then
         qflx_evap_soi(p) = 0._r8
       end if

       ! SCALED sensible and latent heat flux for error check
       eflx_sh_grnd_scale(p)  = -forc_rho(g)*cpair*wtus(c)*dth(l)
       qflx_evap_soi_scale(p) = -forc_rho(g)*wtuq(c)*dqh(l)

    end do

    ! Check to see that total sensible and latent heat equal the sum of
    ! the scaled heat fluxes above
    do fl = 1, num_urbanl
       l = filter_urbanl(fl)
       g = lgridcell(l)
       eflx(l)       = -(forc_rho(g)*cpair/rahu(l))*(thm_g(l) - taf(l))
       qflx(l)       = -(forc_rho(g)/rawu(l))*(forc_q(g) - qaf(l))
       eflx_scale(l) = sum(eflx_sh_grnd_scale(pfti(l):pftf(l)))
       qflx_scale(l) = sum(qflx_evap_soi_scale(pfti(l):pftf(l)))
       eflx_err(l)   = eflx_scale(l) - eflx(l)
       qflx_err(l)   = qflx_scale(l) - qflx(l)
    end do 

    found = .false.
    do fl = 1, num_urbanl
       l = filter_urbanl(fl)
       if (abs(eflx_err(l)) > 0.01_r8) then
          found = .true.
          indexl = l
          exit
       end if
    end do
    if ( found ) then
       write(iulog,*)'WARNING:  Total sensible heat does not equal sum of scaled heat fluxes for urban columns ',&
            ' nstep = ',nstep,' indexl= ',indexl,' eflx_err= ',eflx_err(indexl)
       if (abs(eflx_err(indexl)) > .01_r8) then
          write(iulog,*)'clm model is stopping - error is greater than .01 W/m**2'
          write(iulog,*)'eflx_scale    = ',eflx_scale(indexl)
          write(iulog,*)'eflx_sh_grnd_scale: ',eflx_sh_grnd_scale(pfti(indexl):pftf(indexl))
          write(iulog,*)'eflx          = ',eflx(indexl)
          call endrun(decomp_index=indexl, clmlevel=namel, msg=errmsg(__FILE__, __LINE__))
       end if
    end if

    found = .false.
    do fl = 1, num_urbanl
       l = filter_urbanl(fl)
       ! 4.e-9 kg/m**2/s = 0.01 W/m**2
       if (abs(qflx_err(l)) > 4.e-9_r8) then
          found = .true.
          indexl = l
          exit
       end if
    end do
    if ( found ) then
       write(iulog,*)'WARNING:  Total water vapor flux does not equal sum of scaled water vapor fluxes for urban columns ',&
            ' nstep = ',nstep,' indexl= ',indexl,' qflx_err= ',qflx_err(indexl)
       if (abs(qflx_err(indexl)) > 4.e-9_r8) then
          write(iulog,*)'clm model is stopping - error is greater than 4.e-9 kg/m**2/s'
          write(iulog,*)'qflx_scale    = ',qflx_scale(indexl)
          write(iulog,*)'qflx          = ',qflx(indexl)
          call endrun(decomp_index=indexl, clmlevel=namel, msg=errmsg(__FILE__, __LINE__))
       end if
    end if

    ! Gather terms required to determine internal building temperature

    do fc = 1,num_urbanc
       c = filter_urbanc(fc)
       l = col%landunit(c)

       if (ctype(c) == icol_roof) then
          t_roof_innerl(l) = t_soisno(c,nlevurb)
       else if (ctype(c) == icol_sunwall) then
          t_sunwall_innerl(l) = t_soisno(c,nlevurb)
       else if (ctype(c) == icol_shadewall) then
          t_shadewall_innerl(l) = t_soisno(c,nlevurb)
       end if

    end do

    ! Calculate internal building temperature
    do fl = 1, num_urbanl
       l = filter_urbanl(fl)
     
       lngth_roof = (ht_roof(l)/canyon_hwr(l))*wtlunit_roof(l)/(1._r8-wtlunit_roof(l))
       t_building(l) = (ht_roof(l)*(t_shadewall_innerl(l) + t_sunwall_innerl(l)) &
                        +lngth_roof*t_roof_innerl(l))/(2._r8*ht_roof(l)+lngth_roof)
    end do

    ! No roots for urban except for pervious road

    do j = 1, nlevgrnd
      do f = 1, num_urbanp
         p = filter_urbanp(f)
         c = pcolumn(p)
         if (ctype(c) == icol_road_perv) then
            rootr(p,j) = rootr_road_perv(c,j)
         else
            rootr(p,j) = 0._r8
         end if
       end do
    end do

    do f = 1, num_urbanp

       p = filter_urbanp(f)
       c = pcolumn(p)
       g = pgridcell(p)
       l = plandunit(p)

       ! Use urban canopy air temperature and specific humidity to represent 
       ! 2-m temperature and humidity

       t_ref2m(p) = taf(l)
       q_ref2m(p) = qaf(l)
       t_ref2m_u(p) = taf(l)

       ! 2 m height relative humidity

       call QSat(t_ref2m(p), forc_pbot(g), e_ref2m, de2mdT, qsat_ref2m, dqsat2mdT)
       rh_ref2m(p) = min(100._r8, q_ref2m(p) / qsat_ref2m * 100._r8)
       rh_ref2m_u(p) = rh_ref2m(p)

       ! Variables needed by history tape

       t_veg(p) = forc_t(g)

       ! Add the following to avoid NaN

       psnsun(p)                   = 0._r8
       psnsha(p)                   = 0._r8

    end do

    end associate 
  end subroutine UrbanFluxes

end module UrbanMod
