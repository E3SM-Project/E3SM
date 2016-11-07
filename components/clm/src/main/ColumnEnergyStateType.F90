module ColumnEnergyStateType

!   moved from TemperatureType.F90

#include "shr_assert.h"

  !------------------------------------------------------------------------------
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use decompMod       , only : bounds_type
  use abortutils      , only : endrun
  use clm_varctl      , only : use_ed, use_cndv, iulog
  use clm_varpar      , only : nlevsno, nlevgrnd, nlevlak, nlevlak, nlevurb, crop_prog 
  use clm_varcon      , only : spval
  use GridcellType    , only : grc
  use LandunitType    , only : lun                
  use ColumnType      , only : col                
  use PatchType       , only : pft                
  !
  implicit none
  save
  private
  !
  type, public :: column_energy_state

     ! Temperatures
!     real(r8), pointer :: t_veg_patch              (:)   ! patch vegetation temperature (Kelvin)
     real(r8), pointer :: t_h2osfc_col             (:)   ! col surface water temperature
     real(r8), pointer :: t_h2osfc_bef_col         (:)   ! col surface water temperature from time-step before  
     real(r8), pointer :: t_ssbef_col              (:,:) ! col soil/snow temperature before update (-nlevsno+1:nlevgrnd) 
     real(r8), pointer :: t_soisno_col             (:,:) ! col soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd) 
     real(r8), pointer :: t_soi10cm_col            (:)   ! col soil temperature in top 10cm of soil (Kelvin)
     real(r8), pointer :: t_soi17cm_col            (:)   ! col soil temperature in top 17cm of soil (Kelvin)
     real(r8), pointer :: t_lake_col               (:,:) ! col lake temperature (Kelvin)  (1:nlevlak)          
     real(r8), pointer :: t_grnd_col               (:)   ! col ground temperature (Kelvin)
     real(r8), pointer :: t_grnd_r_col             (:)   ! col rural ground temperature (Kelvin)
     real(r8), pointer :: t_grnd_u_col             (:)   ! col urban ground temperature (Kelvin) (needed by Hydrology2Mod)
!     real(r8), pointer :: t_building_lun           (:)   ! lun internal building temperature (K)
     real(r8), pointer :: snot_top_col             (:)   ! col temperature of top snow layer [K]
     real(r8), pointer :: dTdz_top_col             (:)   ! col temperature gradient in top layer  [K m-1]
!     real(r8), pointer :: dt_veg_patch             (:)   ! patch change in t_veg, last iteration (Kelvin)

     real(r8), pointer :: dt_grnd_col              (:)   ! col change in t_grnd, last iteration (Kelvin)
     real(r8), pointer :: thv_col                  (:)   ! col virtual potential temperature (kelvin)
!     real(r8), pointer :: thm_patch                (:)   ! patch intermediate variable (forc_t+0.0098*forc_hgt_t_patch)
!     real(r8), pointer :: t_a10_patch              (:)   ! patch 10-day running mean of the 2 m temperature (K)
!     real(r8), pointer :: t_a10min_patch           (:)   ! patch 10-day running mean of min 2-m temperature
!     real(r8), pointer :: t_a5min_patch            (:)   ! patch 5-day running mean of min 2-m temperature

!     real(r8), pointer :: taf_lun                  (:)   ! lun urban canopy air temperature (K)

!     real(r8), pointer :: t_ref2m_patch            (:)   ! patch 2 m height surface air temperature (Kelvin)
!     real(r8), pointer :: t_ref2m_r_patch          (:)   ! patch rural 2 m height surface air temperature (Kelvin)
!     real(r8), pointer :: t_ref2m_u_patch          (:)   ! patch urban 2 m height surface air temperature (Kelvin)
!     real(r8), pointer :: t_ref2m_min_patch        (:)   ! patch daily minimum of average 2 m height surface air temperature (K)
!    real(r8), pointer :: t_ref2m_min_r_patch      (:)   ! patch daily minimum of average 2 m height surface air temperature - rural(K)
!     real(r8), pointer :: t_ref2m_min_u_patch      (:)   ! patch daily minimum of average 2 m height surface air temperature - urban (K)
!    real(r8), pointer :: t_ref2m_max_patch        (:)   ! patch daily maximum of average 2 m height surface air temperature (K)
!     real(r8), pointer :: t_ref2m_max_r_patch      (:)   ! patch daily maximum of average 2 m height surface air temperature - rural(K)
!     real(r8), pointer :: t_ref2m_max_u_patch      (:)   ! patch daily maximum of average 2 m height surface air temperature - urban (K)
!     real(r8), pointer :: t_ref2m_min_inst_patch   (:)   ! patch instantaneous daily min of average 2 m height surface air temp (K)
!     real(r8), pointer :: t_ref2m_min_inst_r_patch (:)   ! patch instantaneous daily min of average 2 m height surface air temp - rural (K)
!     real(r8), pointer :: t_ref2m_min_inst_u_patch (:)   ! patch instantaneous daily min of average 2 m height surface air temp - urban (K)
!     real(r8), pointer :: t_ref2m_max_inst_patch   (:)   ! patch instantaneous daily max of average 2 m height surface air temp (K)
!    real(r8), pointer :: t_ref2m_max_inst_r_patch (:)   ! patch instantaneous daily max of average 2 m height surface air temp - rural (K)
!     real(r8), pointer :: t_ref2m_max_inst_u_patch (:)   ! patch instantaneous daily max of average 2 m height surface air temp - urban (K)

     ! Accumulated quantities
     !
     ! TODO(wjs, 2014-08-05) Move these to the module(s) where they are used, to improve
     ! modularity. In cases where they are used by two completely different modules,
     ! which only use the same variable out of convenience, introduce a duplicate (point
     ! being: that way one parameterization is free to change the exact meaning of its
     ! accumulator without affecting the other).
     !
!     real(r8), pointer :: t_veg24_patch           (:)   ! patch 24hr average vegetation temperature (K)
!     real(r8), pointer :: t_veg240_patch          (:)   ! patch 240hr average vegetation temperature (Kelvin)
!     real(r8), pointer :: gdd0_patch              (:)   ! patch growing degree-days base  0C from planting  (ddays)
!     real(r8), pointer :: gdd8_patch              (:)   ! patch growing degree-days base  8C from planting  (ddays)
!     real(r8), pointer :: gdd10_patch             (:)   ! patch growing degree-days base 10C from planting  (ddays)
!     real(r8), pointer :: gdd020_patch            (:)   ! patch 20-year average of gdd0                     (ddays)
!     real(r8), pointer :: gdd820_patch            (:)   ! patch 20-year average of gdd8                     (ddays)
!     real(r8), pointer :: gdd1020_patch           (:)   ! patch 20-year average of gdd10                    (ddays)

     ! Heat content
     real(r8), pointer :: beta_col                 (:)   ! coefficient of convective velocity [-]
     real(r8), pointer :: hc_soi_col               (:)   ! col soil heat content (MJ/m2)
     real(r8), pointer :: hc_soisno_col            (:)   ! col soil plus snow heat content (MJ/m2)
!     real(r8), pointer :: heat1_grc                (:)   ! grc initial gridcell total heat content
!     real(r8), pointer :: heat2_grc                (:)   ! grc post land cover change total heat content

     ! Flags
     integer , pointer :: imelt_col                (:,:) ! flag for melting (=1), freezing (=2), Not=0 (-nlevsno+1:nlevgrnd) 

     ! Emissivities
!     real(r8), pointer :: emv_patch                (:)   ! patch vegetation emissivity 
     real(r8), pointer :: emg_col                  (:)   ! col ground emissivity

     ! Misc
     real(r8), pointer    :: xmf_col               (:)   ! total latent heat of phase change of ground water
     real(r8), pointer    :: xmf_h2osfc_col        (:)   ! latent heat of phase change of surface water
     real(r8), pointer    :: fact_col              (:,:) ! used in computing tridiagonal matrix
     real(r8), pointer    :: c_h2osfc_col          (:)   ! heat capacity of surface water

     ! For VSFM model
     real(r8), pointer :: t_soil_col_1d            (:)   ! 1D temperature of soil layers (Kelvin)

   contains

     procedure, public  :: Init  => init_col_es     
     procedure, public  :: Restart  => restart_col_es    
     procedure, private :: InitAllocate => initallocat_col_es
     procedure, private :: InitHistory  => inithistory_col_es
     procedure, private :: InitCold     =initcold_col_es
!     procedure, public  :: InitAccBuffer
!     procedure, public  :: InitAccVars
!     procedure, public  :: UpdateAccVars

  end type column_energy_state


contains

  !------------------------------------------------------------------------
  subroutine init_col_es(this, bounds, &
       em_roof_lun,  em_wall_lun, em_improad_lun, em_perroad_lun)

!  subroutine init_col_es(this, bounds)

    class(temperature_type)        :: this
    type(bounds_type) , intent(in) :: bounds  
!    real(r8)          , intent(in) :: em_roof_lun(bounds%begl:)
!    real(r8)          , intent(in) :: em_wall_lun(bounds%begl:)
!    real(r8)          , intent(in) :: em_improad_lun(bounds%begl:)
!    real(r8)          , intent(in) :: em_perroad_lun(bounds%begl:)

    call this%InitAllocate ( bounds )
    call this%InitHistory ( bounds )
    call this%InitCold ( bounds,                  &
         em_roof_lun(bounds%begl:bounds%endl),    &
         em_wall_lun(bounds%begl:bounds%endl),    &
         em_improad_lun(bounds%begl:bounds%endl), &
         em_perroad_lun(bounds%begl:bounds%endl))

  end subroutine init_col_es

 subroutine initallocate_col_es(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    class(temperature_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    integer :: begl, endl
    integer :: begg, endg
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    begl = bounds%begl; endl= bounds%endl
    begg = bounds%begg; endg= bounds%endg

    ! Temperatures

    allocate(this%t_h2osfc_col             (begc:endc))                      ; this%t_h2osfc_col             (:)   = nan
    allocate(this%t_h2osfc_bef_col         (begc:endc))                      ; this%t_h2osfc_bef_col         (:)   = nan
    allocate(this%t_ssbef_col              (begc:endc,-nlevsno+1:nlevgrnd))  ; this%t_ssbef_col              (:,:) = nan
    allocate(this%t_soisno_col             (begc:endc,-nlevsno+1:nlevgrnd))  ; this%t_soisno_col             (:,:) = nan
    allocate(this%t_lake_col               (begc:endc,1:nlevlak))            ; this%t_lake_col               (:,:) = nan
    allocate(this%t_grnd_col               (begc:endc))                      ; this%t_grnd_col               (:)   = nan
    allocate(this%t_grnd_r_col             (begc:endc))                      ; this%t_grnd_r_col             (:)   = nan
    allocate(this%t_grnd_u_col             (begc:endc))                      ; this%t_grnd_u_col             (:)   = nan
    allocate(this%snot_top_col             (begc:endc))                      ; this%snot_top_col             (:)   = nan
    allocate(this%dTdz_top_col             (begc:endc))                      ; this%dTdz_top_col             (:)   = nan

    allocate(this%t_soi10cm_col            (begc:endc))                      ; this%t_soi10cm_col            (:)   = nan
    allocate(this%t_soi17cm_col            (begc:endc))                      ; this%t_soi17cm_col            (:)   = spval
    allocate(this%dt_grnd_col              (begc:endc))                      ; this%dt_grnd_col              (:)   = nan
    allocate(this%thv_col                  (begc:endc))                      ; this%thv_col                  (:)   = nan

    ! Heat content
    allocate(this%beta_col                 (begc:endc))                      ; this%beta_col                 (:)   = nan
    allocate(this%hc_soi_col               (begc:endc))                      ; this%hc_soi_col               (:)   = nan
    allocate(this%hc_soisno_col            (begc:endc))                      ; this%hc_soisno_col            (:)   = nan

    ! flags
    allocate(this%imelt_col                (begc:endc,-nlevsno+1:nlevgrnd))  ; this%imelt_col                (:,:) = huge(1)

    ! emissivities
    allocate(this%emg_col                  (begc:endc))                      ; this%emg_col                  (:)   = nan

    allocate(this%xmf_col                  (begc:endc))                      ; this%xmf_col                  (:)   = nan
    allocate(this%xmf_h2osfc_col           (begc:endc))                      ; this%xmf_h2osfc_col           (:)   = nan
    allocate(this%fact_col                 (begc:endc, -nlevsno+1:nlevgrnd)) ; this%fact_col                 (:,:) = nan
    allocate(this%c_h2osfc_col             (begc:endc))                      ; this%c_h2osfc_col             (:)   = nan

    ! For VSFM model
    allocate(this%t_soil_col_1d            ((endc-begc+1)*nlevgrnd))         ; this%t_soil_col_1d            (:) = nan

  end subroutine initallocate_col_es

 !------------------------------------------------------------------------
  subroutine inithistory_col_es(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    use clm_varctl     , only : use_cn, use_cndv
    use histFileMod    , only : hist_addfld1d, hist_addfld2d, no_snow_normal
    !
    ! !ARGUMENTS:
    class(temperature_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer           :: begp, endp
    integer           :: begc, endc
    integer           :: begl, endl
    integer           :: begg, endg
    character(10)     :: active
    real(r8), pointer :: data2dptr(:,:), data1dptr(:) ! temp. pointers for slicing larger arrays
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    begl = bounds%begl; endl= bounds%endl
    begg = bounds%begg; endg= bounds%endg

    this%t_h2osfc_col(begc:endc) = spval
    call hist_addfld1d (fname='TH2OSFC',  units='K',  &
         avgflag='A', long_name='surface water temperature', &
         ptr_col=this%t_h2osfc_col)

    this%t_grnd_u_col(begc:endc) = spval
    call hist_addfld1d (fname='TG_U', units='K',  &
         avgflag='A', long_name='Urban ground temperature', &
         ptr_col=this%t_grnd_u_col, set_nourb=spval, c2l_scale_type='urbans')

    this%t_lake_col(begc:endc,:) = spval
    call hist_addfld2d (fname='TLAKE',  units='K', type2d='levlak', &
         avgflag='A', long_name='lake temperature', &
         ptr_col=this%t_lake_col)

    this%t_soisno_col(begc:endc,-nlevsno+1:0) = spval
    data2dptr => this%t_soisno_col(:,-nlevsno+1:0)
    call hist_addfld2d (fname='SNO_T', units='K', type2d='levsno',  &
         avgflag='A', long_name='Snow temperatures', &
         ptr_col=data2dptr, no_snow_behavior=no_snow_normal, default='inactive')

    this%t_grnd_col(begc:endc) = spval
    call hist_addfld1d (fname='TG',  units='K',  &
         avgflag='A', long_name='ground temperature', &
         ptr_col=this%t_grnd_col, c2l_scale_type='urbans')

    this%t_grnd_r_col(begc:endc) = spval
    call hist_addfld1d (fname='TG_R', units='K',  &
         avgflag='A', long_name='Rural ground temperature', &
         ptr_col=this%t_grnd_r_col, set_spec=spval)

    this%t_soisno_col(begc:endc,:) = spval
    call hist_addfld2d (fname='TSOI',  units='K', type2d='levgrnd', &
         avgflag='A', long_name='soil temperature (vegetated landunits only)', &
         ptr_col=this%t_soisno_col, l2g_scale_type='veg')

    this%t_soisno_col(begc:endc,:) = spval
    call hist_addfld2d (fname='TSOI_ICE',  units='K', type2d='levgrnd', &
         avgflag='A', long_name='soil temperature (ice landunits only)', &
         ptr_col=this%t_soisno_col, l2g_scale_type='ice')

    this%t_soi10cm_col(begc:endc) = spval
    call hist_addfld1d (fname='TSOI_10CM',  units='K', &
         avgflag='A', long_name='soil temperature in top 10cm of soil', &
         ptr_col=this%t_soi10cm_col, set_urb=spval)


    this%hc_soi_col(begc:endc) = spval
    call hist_addfld1d (fname='HCSOI',  units='MJ/m2',  &
         avgflag='A', long_name='soil heat content', &
         ptr_col=this%hc_soi_col, set_lake=spval, set_urb=spval, l2g_scale_type='veg')

    this%hc_soisno_col(begc:endc) = spval
    call hist_addfld1d (fname='HC',  units='MJ/m2',  &
         avgflag='A', long_name='heat content of soil/snow/lake', &
         ptr_col=this%hc_soisno_col, set_urb=spval) 

    this%snot_top_col(begc:endc) = spval 
    call hist_addfld1d (fname='SNOTTOPL', units='K/m', &
         avgflag='A', long_name='snow temperature (top layer)', &
         ptr_col=this%snot_top_col, set_urb=spval, default='inactive')

    this%dTdz_top_col(begc:endc) = spval 
    call hist_addfld1d (fname='SNOdTdzL', units='K/m', &
         avgflag='A', long_name='top snow layer temperature gradient (land)', &
         ptr_col=this%dTdz_top_col, set_urb=spval, default='inactive')

    if (use_cn) then
       this%emg_col(begc:endc) = spval
       call hist_addfld1d (fname='EMG', units='proportion', &
            avgflag='A', long_name='ground emissivity', &
            ptr_col=this%emg_col, default='inactive')
    end if

    if (use_cn) then
       this%beta_col(begc:endc) = spval
       call hist_addfld1d (fname='BETA', units='none', &
            avgflag='A', long_name='coefficient of convective velocity', &
            ptr_col=this%beta_col, default='inactive')
    end if


  end subroutine inithistory_col_es


  !-----------------------------------------------------------------------
  subroutine initcold_col_es(this, bounds, &
       em_roof_lun,  em_wall_lun, em_improad_lun, em_perroad_lun)
    !
    ! !DESCRIPTION:
    ! Initialize cold start conditions for module variables
    !
    ! !USES:
    use shr_kind_mod   , only : r8 => shr_kind_r8
    use shr_const_mod  , only : SHR_CONST_TKFRZ
    use clm_varcon     , only : denice, denh2o, sb
    use landunit_varcon, only : istice, istwet, istsoil, istdlak, istice_mec
    use column_varcon  , only : icol_road_imperv, icol_roof, icol_sunwall
    use column_varcon  , only : icol_shadewall, icol_road_perv
    use clm_varctl     , only : iulog, use_vancouver, use_mexicocity
    !
    ! !ARGUMENTS:
    class(temperature_type)        :: this
    type(bounds_type) , intent(in) :: bounds  
    real(r8)          , intent(in) :: em_roof_lun(bounds%begl:bounds%endl)
    real(r8)          , intent(in) :: em_wall_lun(bounds%begl:bounds%endl)
    real(r8)          , intent(in) :: em_improad_lun(bounds%begl:bounds%endl)
    real(r8)          , intent(in) :: em_perroad_lun(bounds%begl:bounds%endl)
    !
    ! !LOCAL VARIABLES:
    integer  :: j,l,c,p ! indices
    integer  :: nlevs   ! number of levels
    real(r8) :: snowbd  ! temporary calculation of snow bulk density (kg/m3)
    real(r8) :: fmelt   ! snowbd/100
    integer  :: lev
    !-----------------------------------------------------------------------


    !!!!!!!!   ********DW********* need to handle this later

    SHR_ASSERT_ALL((ubound(em_roof_lun)    == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(em_wall_lun)    == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(em_improad_lun) == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(em_perroad_lun) == (/bounds%endl/)), errMsg(__FILE__, __LINE__))

    associate(snl => col%snl) ! Output: [integer (:)    ]  number of snow layers   

      ! Set snow/soil temperature
      ! t_lake only has valid values over non-lake   
      ! t_soisno, t_grnd and t_veg have valid values over all land 

      do c = bounds%begc,bounds%endc
         l = col%landunit(c)

         this%t_soisno_col(c,-nlevsno+1:nlevgrnd) = spval

         ! Snow level temperatures - all land points
         if (snl(c) < 0) then
            do j = snl(c)+1, 0
               this%t_soisno_col(c,j) = 250._r8
            end do
         end if

         ! Below snow temperatures - nonlake points (lake points are set below)
         if (.not. lun%lakpoi(l)) then 

            if (lun%itype(l)==istice .or. lun%itype(l)==istice_mec) then
               this%t_soisno_col(c,1:nlevgrnd) = 250._r8

            else if (lun%itype(l) == istwet) then
               this%t_soisno_col(c,1:nlevgrnd) = 277._r8

            else if (lun%urbpoi(l)) then
               if (use_vancouver) then
                  if (col%itype(c) == icol_road_perv .or. col%itype(c) == icol_road_imperv) then 
                     ! Set road top layer to initial air temperature and interpolate other
                     ! layers down to 20C in bottom layer
                     do j = 1, nlevgrnd
                        this%t_soisno_col(c,j) = 297.56 - (j-1) * ((297.56-293.16)/(nlevgrnd-1)) 
                     end do
                     ! Set wall and roof layers to initial air temperature
                  else if (col%itype(c) == icol_sunwall .or. col%itype(c) == icol_shadewall .or. col%itype(c) == icol_roof) then
                     this%t_soisno_col(c,1:nlevurb) = 297.56
                  else
                     this%t_soisno_col(c,1:nlevgrnd) = 283._r8
                  end if
               else if (use_mexicocity) then
                  if (col%itype(c) == icol_road_perv .or. col%itype(c) == icol_road_imperv) then 
                     ! Set road top layer to initial air temperature and interpolate other
                     ! layers down to 22C in bottom layer
                     do j = 1, nlevgrnd
                        this%t_soisno_col(c,j) = 289.46 - (j-1) * ((289.46-295.16)/(nlevgrnd-1)) 
                     end do
                  else if (col%itype(c) == icol_sunwall .or. col%itype(c) == icol_shadewall .or. col%itype(c) == icol_roof) then
                     ! Set wall and roof layers to initial air temperature
                     this%t_soisno_col(c,1:nlevurb) = 289.46
                  else
                     this%t_soisno_col(c,1:nlevgrnd) = 283._r8
                  end if
               else
                  if (col%itype(c) == icol_road_perv .or. col%itype(c) == icol_road_imperv) then 
                     this%t_soisno_col(c,1:nlevgrnd) = 274._r8
                  else if (col%itype(c) == icol_sunwall .or. col%itype(c) == icol_shadewall &
                       .or. col%itype(c) == icol_roof) then
                     ! Set sunwall, shadewall, roof to fairly high temperature to avoid initialization
                     ! shock from large heating/air conditioning flux
                     this%t_soisno_col(c,1:nlevurb) = 292._r8
                  end if
               end if
            else
               this%t_soisno_col(c,1:nlevgrnd) = 274._r8
            endif
         endif
      end do

      ! Set Ground temperatures

      do c = bounds%begc,bounds%endc
         l = col%landunit(c)

         if (lun%lakpoi(l)) then 
            this%t_grnd_col(c) = 277._r8
         else
            this%t_grnd_col(c) = this%t_soisno_col(c,snl(c)+1)
         end if
         this%t_soi17cm_col(c) = this%t_grnd_col(c)
      end do

      do c = bounds%begc,bounds%endc
         l = col%landunit(c)
         if (lun%lakpoi(l)) then ! lake
            this%t_lake_col(c,1:nlevlak) = this%t_grnd_col(c)
            this%t_soisno_col(c,1:nlevgrnd) = this%t_grnd_col(c)
         end if
      end do

      ! Set t_h2osfc_col

      this%t_h2osfc_col(bounds%begc:bounds%endc)  = 274._r8

    end associate

    do c = bounds%begc,bounds%endc
       l = col%landunit(c)

       if (col%itype(c) == icol_roof       ) this%emg_col(c) = em_roof_lun(l)
       if (col%itype(c) == icol_sunwall    ) this%emg_col(c) = em_wall_lun(l)
       if (col%itype(c) == icol_shadewall  ) this%emg_col(c) = em_wall_lun(l)
       if (col%itype(c) == icol_road_imperv) this%emg_col(c) = em_improad_lun(l)
       if (col%itype(c) == icol_road_perv  ) this%emg_col(c) = em_perroad_lun(l)
    end do

  end subroutine initcold_col_es


end module ColumnEnergyStateType