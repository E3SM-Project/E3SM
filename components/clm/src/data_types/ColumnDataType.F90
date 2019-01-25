module ColumnDataType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Column data type allocation and initialization
  ! -------------------------------------------------------- 
  !
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use clm_varpar     , only : nlevsno, nlevgrnd, nlevlak, nlevurb
  use clm_varcon     , only : spval, ispval
  use histFileMod    , only : hist_addfld1d, hist_addfld2d, no_snow_normal
  use ncdio_pio      , only : file_desc_t, ncd_double
  use decompMod      , only : bounds_type
  use restUtilMod
  use ColumnType     , only : col_pp
  use LandunitType   , only : lun_pp
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds energy state information at the column level.
  !-----------------------------------------------------------------------
  type, public :: column_energy_state
    ! temperature variables
    real(r8), pointer :: t_soisno      (:,:) => null() ! soil temperature (K)  (-nlevsno+1:nlevgrnd) 
    real(r8), pointer :: t_ssbef       (:,:) => null() ! soil/snow temperature before update (K) (-nlevsno+1:nlevgrnd) 
    real(r8), pointer :: t_h2osfc      (:)   => null() ! surface water temperature (K)
    real(r8), pointer :: t_h2osfc_bef  (:)   => null() ! surface water temperature at start of time step (K)
    real(r8), pointer :: t_soi10cm     (:)   => null() ! soil temperature in top 10cm of soil (K)
    real(r8), pointer :: t_soi17cm     (:)   => null() ! soil temperature in top 17cm of soil (K)
    real(r8), pointer :: t_grnd        (:)   => null() ! ground temperature (K)
    real(r8), pointer :: t_lake        (:,:) => null() ! lake temperature (K)  (1:nlevlak)          
    real(r8), pointer :: t_grnd_r      (:)   => null() ! rural ground temperature (K)
    real(r8), pointer :: t_grnd_u      (:)   => null() ! urban ground temperature (K) (needed by Hydrology2Mod)
    real(r8), pointer :: snot_top      (:)   => null() ! temperature of top snow layer (K)
    real(r8), pointer :: dTdz_top      (:)   => null() ! temperature gradient in top layer  (K m-1)
    real(r8), pointer :: thv           (:)   => null() ! virtual potential temperature (K)
    ! heat content variables (diagnostic)
    real(r8), pointer :: hc_soi        (:)   => null() ! soil heat content (MJ/m2)
    real(r8), pointer :: hc_soisno     (:)   => null() ! soil plus snow heat content (MJ/m2)
    ! other quantities related to energy state
    real(r8), pointer :: emg           (:)   => null() ! ground emissivity (unitless)
    real(r8), pointer :: fact          (:,:) => null() ! factors used in computing tridiagonal matrix
    real(r8), pointer :: c_h2osfc      (:)   => null() ! heat capacity of surface water (J/K)
    ! For coupling with pflotran TH
    real(r8), pointer :: t_nearsurf    (:)   => null() ! near-surface air temperature averaged over bare-veg (K)
    
  contains
    procedure, public :: Init    => col_es_init
    procedure, public :: Restart => col_es_restart
    procedure, public :: Clean   => col_es_clean
  end type column_energy_state
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds water state information at the column level.
  !-----------------------------------------------------------------------
  type, public :: column_water_state
    real(r8), pointer :: xxx      (:) => null() ! xxx (xxx)
  contains
    procedure, public :: Init  => init_col_ws
    procedure, public :: Clean => clean_col_ws
  end type column_water_state
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds carbon state information at the column level.
  !-----------------------------------------------------------------------
  type, public :: column_carbon_state
    real(r8), pointer :: xxx      (:) => null() ! xxx (xxx)
  contains
    procedure, public :: Init  => init_col_cs
    procedure, public :: Clean => clean_col_cs
  end type column_carbon_state
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds nitrogen state information at the column level.
  !-----------------------------------------------------------------------
  type, public :: column_nitrogen_state
    real(r8), pointer :: xxx      (:) => null() ! xxx (xxx)
  contains
    procedure, public :: Init  => init_col_ns
    procedure, public :: Clean => clean_col_ns
  end type column_nitrogen_state
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds phosphorus state information at the column level.
  !-----------------------------------------------------------------------
  type, public :: column_phosphorus_state
    real(r8), pointer :: xxx      (:) => null() ! xxx (xxx)
  contains
    procedure, public :: Init  => init_col_ps
    procedure, public :: Clean => clean_col_ps
  end type column_phosphorus_state

  !-----------------------------------------------------------------------
  ! Define the data structure that holds energy flux information at the column level.
  !-----------------------------------------------------------------------
  type, public :: column_energy_flux
    real(r8), pointer :: eflx_h2osfc_to_snow     (:)   => null() ! snow melt to h2osfc heat flux (W/m2)
    real(r8), pointer :: eflx_snomelt            (:)   => null() ! snow melt heat flux (W/m**2)
    real(r8), pointer :: eflx_snomelt_r          (:)   => null() ! rural snow melt heat flux (W/m2)
    real(r8), pointer :: eflx_snomelt_u          (:)   => null() ! urban snow melt heat flux (W/m2)
    real(r8), pointer :: eflx_bot                (:)   => null() ! heat flux from beneath the soil or ice column (W/m2)
    real(r8), pointer :: eflx_fgr12              (:)   => null() ! ground heat flux between soil layers 1 and 2 (W/m2)
    real(r8), pointer :: eflx_fgr                (:,:) => null() ! (rural) soil downward heat flux (W/m2) (1:nlevgrnd)  (pos upward; usually eflx_bot >= 0)
    real(r8), pointer :: eflx_building_heat      (:)   => null() ! heat flux from urban building interior to urban walls, roof (W/m2)
    real(r8), pointer :: eflx_urban_ac           (:)   => null() ! urban air conditioning flux (W/m2)
    real(r8), pointer :: eflx_urban_heat         (:)   => null() ! urban heating flux (W/m**2)
    real(r8), pointer :: eflx_hs_h2osfc          (:)   => null() ! heat flux on standing water (W/m2)
    real(r8), pointer :: eflx_hs_top_snow        (:)   => null() ! heat flux on top snow layer (W/m2)
    real(r8), pointer :: eflx_hs_soil            (:)   => null() ! heat flux on soil [W/m2
    real(r8), pointer :: eflx_sabg_lyr           (:,:) => null() ! absorbed solar radiation (col,lyr) (W/m2)
    ! Derivatives of energy fluxes                      
    real(r8), pointer :: eflx_dhsdT              (:)   => null() ! deriv. of energy flux into surface layer wrt temp (W/m2/K)
    ! Latent heat terms
    real(r8), pointer :: htvp                    (:)   => null() ! latent heat of vapor of water (or sublimation) [j/kg]
    real(r8), pointer :: xmf                     (:)   => null() ! total latent heat of phase change of ground water
    real(r8), pointer :: xmf_h2osfc              (:)   => null() ! latent heat of phase change of surface water
    integer , pointer :: imelt                   (:,:) ! flag for melting (=1), freezing (=2), Not=0 (-nlevsno+1:nlevgrnd) 
    ! for couplig with pflotran                         
    real(r8), pointer :: eflx_soil_grnd          (:)   => null() ! integrated soil ground heat flux (W/m2)  [+ = into ground]
    real(r8), pointer :: eflx_rnet_soil          (:)   => null() ! soil net (sw+lw) radiation flux (W/m2) [+ = into soil]
    real(r8), pointer :: eflx_fgr0_soil          (:)   => null() ! soil-air heat flux (W/m2) [+ = into soil]
    real(r8), pointer :: eflx_fgr0_snow          (:)   => null() ! soil-snow heat flux (W/m2) [+ = into soil]
    real(r8), pointer :: eflx_fgr0_h2osfc        (:)   => null() ! soil-surfacewater heat flux (W/m2) [+ = into soil]
    ! Balance Checks                                    
    real(r8), pointer :: errsoi                  (:)   => null() ! soil/lake energy conservation error   (W/m2)
    real(r8), pointer :: errseb                  (:)   => null() ! surface energy conservation error     (W/m2)
    real(r8), pointer :: errsol                  (:)   => null() ! solar radiation conservation error    (W/m2)
    real(r8), pointer :: errlon                  (:)   => null() ! longwave radiation conservation error (W/m2)

  contains
    procedure, public :: Init    => col_ef_init
    procedure, public :: Restart => col_ef_restart
    procedure, public :: Clean   => col_ef_clean
  end type column_energy_flux
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds water flux information at the column level.
  !-----------------------------------------------------------------------
  type, public :: column_water_flux
    real(r8), pointer :: xxx      (:) => null() ! xxx (xxx)
  contains
    procedure, public :: Init  => init_col_wf
    procedure, public :: Clean => clean_col_wf
  end type column_water_flux
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds carbon flux information at the column level.
  !-----------------------------------------------------------------------
  type, public :: column_carbon_flux
    real(r8), pointer :: xxx      (:) => null() ! xxx (xxx)
  contains
    procedure, public :: Init  => init_col_cf
    procedure, public :: Clean => clean_col_cf
  end type column_carbon_flux
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds nitrogen flux information at the column level.
  !-----------------------------------------------------------------------
  type, public :: column_nitrogen_flux
    real(r8), pointer :: xxx      (:) => null() ! xxx (xxx)
  contains
    procedure, public :: Init  => init_col_nf
    procedure, public :: Clean => clean_col_nf
  end type column_nitrogen_flux
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds phosphorus flux information at the column level.
  !-----------------------------------------------------------------------
  type, public :: column_phosphorus_flux
    real(r8), pointer :: xxx      (:) => null() ! xxx (xxx)
  contains
    procedure, public :: Init  => init_col_pf
    procedure, public :: Clean => clean_col_pf
  end type column_phosphorus_flux

   
  !-----------------------------------------------------------------------
  ! declare the public instances of column-level data types
  !-----------------------------------------------------------------------
  type(column_energy_state)          , public, target :: col_es    ! column energy state
  type(column_water_state)           , public, target :: col_ws    ! column water state
  type(column_carbon_state)          , public, target :: col_cs    ! column carbon state
  type(column_nitrogen_state)        , public, target :: col_ns    ! column nitrogen state
  type(column_phosphorus_state)      , public, target :: col_ps    ! column phosphorus state
  type(column_energy_flux)           , public, target :: col_ef    ! column energy flux
  type(column_water_flux)            , public, target :: col_wf    ! column water flux
  type(column_carbon_flux)           , public, target :: col_cf    ! column carbon flux
  type(column_nitrogen_flux)         , public, target :: col_nf    ! column nitrogen flux
  type(column_phosphorus_flux)       , public, target :: col_pf    ! column phosphorus flux

  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean column energy state data structure
  !------------------------------------------------------------------------
  subroutine col_es_init(this, begc, endc)
    !
    ! !USES:
    use landunit_varcon, only : istice, istwet, istsoil, istdlak, istice_mec
    use clm_varctl     , only : iulog, use_cn, use_vancouver, use_mexicocity
    use column_varcon  , only : icol_road_perv, icol_road_imperv, icol_roof, icol_sunwall, icol_shadewall
    use UrbanParamsType, only : urbanparams_vars
    !
    ! !ARGUMENTS:
    class(column_energy_state) :: this
    integer, intent(in) :: begc,endc
    !------------------------------------------------------------------------
    !
    ! !LOCAL VARIABLES:
    integer           :: c,l,j                        ! indices
    real(r8), pointer :: data2dptr(:,:), data1dptr(:) ! temp. pointers for slicing larger arrays

    !------------------------------------------------------------------------------
    associate(snl => col_pp%snl) ! Output: [integer (:)    ]  number of snow layers   

    !-----------------------------------------------------------------------
    ! allocate for each member of col_es
    !-----------------------------------------------------------------------
    allocate(this%t_soisno         (begc:endc,-nlevsno+1:nlevgrnd)) ; this%t_soisno           (:,:) = nan
    allocate(this%t_ssbef          (begc:endc,-nlevsno+1:nlevgrnd)) ; this%t_ssbef            (:,:) = nan
    allocate(this%t_h2osfc         (begc:endc))                     ; this%t_h2osfc           (:)   = nan
    allocate(this%t_h2osfc_bef     (begc:endc))                     ; this%t_h2osfc_bef       (:)   = nan
    allocate(this%t_soi10cm        (begc:endc))                     ; this%t_soi10cm          (:)   = nan
    allocate(this%t_soi17cm        (begc:endc))                     ; this%t_soi17cm          (:)   = spval
    allocate(this%t_grnd           (begc:endc))                     ; this%t_grnd             (:)   = nan
    allocate(this%t_lake           (begc:endc,1:nlevlak))           ; this%t_lake             (:,:) = nan
    allocate(this%t_grnd_r         (begc:endc))                     ; this%t_grnd_r           (:)   = nan
    allocate(this%t_grnd_u         (begc:endc))                     ; this%t_grnd_u           (:)   = nan
    allocate(this%snot_top         (begc:endc))                     ; this%snot_top           (:)   = nan
    allocate(this%dTdz_top         (begc:endc))                     ; this%dTdz_top           (:)   = nan
    allocate(this%thv              (begc:endc))                     ; this%thv                (:)   = nan
    allocate(this%hc_soi           (begc:endc))                     ; this%hc_soi             (:)   = nan
    allocate(this%hc_soisno        (begc:endc))                     ; this%hc_soisno          (:)   = nan
    allocate(this%emg              (begc:endc))                     ; this%emg                (:)   = nan
    allocate(this%fact             (begc:endc, -nlevsno+1:nlevgrnd)); this%fact               (:,:) = nan
    allocate(this%c_h2osfc         (begc:endc))                     ; this%c_h2osfc           (:)   = nan
    allocate(this%t_nearsurf       (begc:endc))                     ; this%t_nearsurf         (:)   = nan

    !-----------------------------------------------------------------------
    ! initialize history fields for select members of col_es
    !-----------------------------------------------------------------------
    this%t_soisno(begc:endc,-nlevsno+1:0) = spval
    data2dptr => this%t_soisno(:,-nlevsno+1:0)
    call hist_addfld2d (fname='SNO_T', units='K', type2d='levsno',  &
         avgflag='A', long_name='Snow temperatures', &
         ptr_col=data2dptr, no_snow_behavior=no_snow_normal, default='inactive')

    this%t_soisno(begc:endc,:) = spval
    call hist_addfld2d (fname='TSOI',  units='K', type2d='levgrnd', &
         avgflag='A', long_name='soil temperature (vegetated landunits only)', &
         ptr_col=this%t_soisno, l2g_scale_type='veg')

    this%t_soisno(begc:endc,:) = spval
    call hist_addfld2d (fname='TSOI_ICE',  units='K', type2d='levgrnd', &
         avgflag='A', long_name='soil temperature (ice landunits only)', &
         ptr_col=this%t_soisno, l2g_scale_type='ice')

    this%t_h2osfc(begc:endc) = spval
    call hist_addfld1d (fname='TH2OSFC',  units='K',  &
         avgflag='A', long_name='surface water temperature', &
         ptr_col=this%t_h2osfc)

    this%t_soi10cm(begc:endc) = spval
    call hist_addfld1d (fname='TSOI_10CM',  units='K', &
         avgflag='A', long_name='soil temperature in top 10cm of soil', &
         ptr_col=this%t_soi10cm, set_urb=spval)

    this%t_grnd(begc:endc) = spval
    call hist_addfld1d (fname='TG',  units='K',  &
         avgflag='A', long_name='ground temperature', &
         ptr_col=this%t_grnd, c2l_scale_type='urbans')

    this%t_lake(begc:endc,:) = spval
    call hist_addfld2d (fname='TLAKE',  units='K', type2d='levlak', &
         avgflag='A', long_name='lake temperature', &
         ptr_col=this%t_lake)

    this%t_grnd_r(begc:endc) = spval
    call hist_addfld1d (fname='TG_R', units='K',  &
         avgflag='A', long_name='Rural ground temperature', &
         ptr_col=this%t_grnd_r, set_spec=spval)

    this%t_grnd_u(begc:endc) = spval
    call hist_addfld1d (fname='TG_U', units='K',  &
         avgflag='A', long_name='Urban ground temperature', &
         ptr_col=this%t_grnd_u, set_nourb=spval, c2l_scale_type='urbans')

    this%snot_top(begc:endc) = spval 
    call hist_addfld1d (fname='SNOTTOPL', units='K', &
         avgflag='A', long_name='snow temperature (top layer)', &
         ptr_col=this%snot_top, set_urb=spval, default='inactive')

    this%dTdz_top(begc:endc) = spval 
    call hist_addfld1d (fname='SNOdTdzL', units='K/m', &
         avgflag='A', long_name='top snow layer temperature gradient (land)', &
         ptr_col=this%dTdz_top, set_urb=spval, default='inactive')

    this%hc_soi(begc:endc) = spval
    call hist_addfld1d (fname='HCSOI',  units='MJ/m2',  &
         avgflag='A', long_name='soil heat content', &
         ptr_col=this%hc_soi, set_lake=spval, set_urb=spval, l2g_scale_type='veg')

    this%hc_soisno(begc:endc) = spval
    call hist_addfld1d (fname='HC',  units='MJ/m2',  &
         avgflag='A', long_name='heat content of soil/snow/lake', &
         ptr_col=this%hc_soisno, set_urb=spval)

    if (use_cn) then
       this%emg(begc:endc) = spval
       call hist_addfld1d (fname='EMG', units='proportion', &
            avgflag='A', long_name='ground emissivity', &
            ptr_col=this%emg, default='inactive')
    end if

    !-----------------------------------------------------------------------
    ! set cold-start initial values for select members of col_es
    !-----------------------------------------------------------------------
    
    ! Initialize soil+snow temperatures
    do c = begc,endc
       l = col_pp%landunit(c)

       ! Snow level temperatures - all land points
       if (snl(c) < 0) then
          do j = snl(c)+1, 0
             this%t_soisno(c,j) = 250._r8
          end do
       end if

       ! Below snow temperatures - nonlake points (lake points are set below)
       if (.not. lun_pp%lakpoi(l)) then 

          if (lun_pp%itype(l)==istice .or. lun_pp%itype(l)==istice_mec) then
             this%t_soisno(c,1:nlevgrnd) = 250._r8

          else if (lun_pp%itype(l) == istwet) then
             this%t_soisno(c,1:nlevgrnd) = 277._r8

          else if (lun_pp%urbpoi(l)) then
             if (use_vancouver) then
                if (col_pp%itype(c) == icol_road_perv .or. col_pp%itype(c) == icol_road_imperv) then 
                   ! Set road top layer to initial air temperature and interpolate other
                   ! layers down to 20C in bottom layer
                   do j = 1, nlevgrnd
                      this%t_soisno(c,j) = 297.56 - (j-1) * ((297.56-293.16)/(nlevgrnd-1)) 
                   end do
                   ! Set wall and roof layers to initial air temperature
                else if (col_pp%itype(c) == icol_sunwall .or. col_pp%itype(c) == icol_shadewall &
                     .or. col_pp%itype(c) == icol_roof) then
                   this%t_soisno(c,1:nlevurb) = 297.56
                else
                   this%t_soisno(c,1:nlevgrnd) = 283._r8
                end if
             else if (use_mexicocity) then
                if (col_pp%itype(c) == icol_road_perv .or. col_pp%itype(c) == icol_road_imperv) then 
                   ! Set road top layer to initial air temperature and interpolate other
                   ! layers down to 22C in bottom layer
                   do j = 1, nlevgrnd
                      this%t_soisno(c,j) = 289.46 - (j-1) * ((289.46-295.16)/(nlevgrnd-1)) 
                   end do
                else if (col_pp%itype(c) == icol_sunwall .or. col_pp%itype(c) == icol_shadewall &
                     .or. col_pp%itype(c) == icol_roof) then
                   ! Set wall and roof layers to initial air temperature
                   this%t_soisno(c,1:nlevurb) = 289.46
                else
                   this%t_soisno(c,1:nlevgrnd) = 283._r8
                end if
             else
                if (col_pp%itype(c) == icol_road_perv .or. col_pp%itype(c) == icol_road_imperv) then 
                   this%t_soisno(c,1:nlevgrnd) = 274._r8
                else if (col_pp%itype(c) == icol_sunwall .or. col_pp%itype(c) == icol_shadewall &
                     .or. col_pp%itype(c) == icol_roof) then
                   ! Set sunwall, shadewall, roof to fairly high temperature to avoid initialization
                   ! shock from large heating/air conditioning flux
                   this%t_soisno(c,1:nlevurb) = 292._r8
                end if
             end if
          else
             this%t_soisno(c,1:nlevgrnd) = 274._r8
          endif
          this%t_grnd(c) = this%t_soisno(c,snl(c)+1)
       endif
       
       if (lun_pp%lakpoi(l)) then ! special handling for lake points
          this%t_soisno(c,1:nlevgrnd) = 277._r8
          this%t_lake(c,1:nlevlak) = 277._r8
          this%t_grnd(c) = 277._r8
       end if

       this%t_soi17cm(c) = this%t_grnd(c)
       
       ! set emissivity for urban columns
       if (col_pp%itype(c) == icol_roof       ) this%emg(c) = urbanparams_vars%em_roof(l)
       if (col_pp%itype(c) == icol_sunwall    ) this%emg(c) = urbanparams_vars%em_wall(l)
       if (col_pp%itype(c) == icol_shadewall  ) this%emg(c) = urbanparams_vars%em_wall(l)
       if (col_pp%itype(c) == icol_road_imperv) this%emg(c) = urbanparams_vars%em_improad(l)
       if (col_pp%itype(c) == icol_road_perv  ) this%emg(c) = urbanparams_vars%em_perroad(l)

    end do ! columns loop
    
    

    ! Initialize surface water temperatures
    this%t_h2osfc = 274._r8
    
    end associate

  end subroutine col_es_init
    
  !------------------------------------------------------------------------
  subroutine col_es_restart(this, bounds, ncid, flag)
    ! 
    ! !DESCRIPTION:
    ! Read/Write column energy state information to/from restart file.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(column_energy_state) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   
    character(len=*) , intent(in)    :: flag   
    !
    ! !LOCAL VARIABLES:
    logical :: readvar   ! determine if variable is on initial file
    !-----------------------------------------------------------------------

    call restartvar(ncid=ncid, flag=flag, varname='T_SOISNO', xtype=ncd_double,   &
         dim1name='column', dim2name='levtot', switchdim=.true., &
         long_name='soil-snow temperature', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_soisno)

    call restartvar(ncid=ncid, flag=flag, varname='TH2OSFC', xtype=ncd_double,  &
         dim1name='column', &
         long_name='surface water temperature', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_h2osfc)
    if (flag=='read' .and. .not. readvar) then
       this%t_h2osfc(bounds%begc:bounds%endc) = 274.0_r8
    end if

    call restartvar(ncid=ncid, flag=flag, varname='T_GRND', xtype=ncd_double,  &
         dim1name='column', &
         long_name='ground temperature', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_grnd)

    call restartvar(ncid=ncid, flag=flag, varname='T_LAKE', xtype=ncd_double,  &
         dim1name='column', dim2name='levlak', switchdim=.true., &
         long_name='lake temperature', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_lake)

    call restartvar(ncid=ncid, flag=flag, varname='T_GRND_R', xtype=ncd_double,  &
         dim1name='column', &
         long_name='rural ground temperature', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_grnd_r)
         
    call restartvar(ncid=ncid, flag=flag, varname='T_GRND_U', xtype=ncd_double, &
         dim1name='column',                    &
         long_name='urban ground temperature', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_grnd_u)

  end subroutine col_es_restart

  !------------------------------------------------------------------------
  subroutine col_es_clean(this)
    !
    ! !ARGUMENTS:
    class(column_energy_state) :: this
    !------------------------------------------------------------------------
    deallocate(this%t_h2osfc)
  end subroutine col_es_clean
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean column water state data structure
  !------------------------------------------------------------------------
  subroutine init_col_ws(this, begc, endc)
    !
    ! !ARGUMENTS:
    class(column_water_state) :: this
    integer, intent(in) :: begc,endc
    !------------------------------------------------------------------------
    
  end subroutine init_col_ws
    
  !------------------------------------------------------------------------
  subroutine clean_col_ws(this)
    !
    ! !ARGUMENTS:
    class(column_water_state) :: this
    !------------------------------------------------------------------------
    
  end subroutine clean_col_ws
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean column carbon state data structure
  !------------------------------------------------------------------------
  subroutine init_col_cs(this, begc, endc)
    !
    ! !ARGUMENTS:
    class(column_carbon_state) :: this
    integer, intent(in) :: begc,endc
    !------------------------------------------------------------------------
    
  end subroutine init_col_cs
    
  !------------------------------------------------------------------------
  subroutine clean_col_cs(this)
    !
    ! !ARGUMENTS:
    class(column_carbon_state) :: this
    !------------------------------------------------------------------------
    
  end subroutine clean_col_cs
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean column nitrogen state data structure
  !------------------------------------------------------------------------
  subroutine init_col_ns(this, begc, endc)
    !
    ! !ARGUMENTS:
    class(column_nitrogen_state) :: this
    integer, intent(in) :: begc,endc
    !------------------------------------------------------------------------
    
  end subroutine init_col_ns
    
  !------------------------------------------------------------------------
  subroutine clean_col_ns(this)
    !
    ! !ARGUMENTS:
    class(column_nitrogen_state) :: this
    !------------------------------------------------------------------------
    
  end subroutine clean_col_ns
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean column phosphorus state data structure
  !------------------------------------------------------------------------
  subroutine init_col_ps(this, begc, endc)
    !
    ! !ARGUMENTS:
    class(column_phosphorus_state) :: this
    integer, intent(in) :: begc,endc
    !------------------------------------------------------------------------
    
  end subroutine init_col_ps
    
  !------------------------------------------------------------------------
  subroutine clean_col_ps(this)
    !
    ! !ARGUMENTS:
    class(column_phosphorus_state) :: this
    !------------------------------------------------------------------------
    
  end subroutine clean_col_ps

  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean column energy flux data structure
  !------------------------------------------------------------------------
  subroutine col_ef_init(this, begc, endc)
    !
    ! !ARGUMENTS:
    class(column_energy_flux) :: this
    integer, intent(in) :: begc,endc
    ! !LOCAL VARIABLES:
    integer  :: l,c
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! allocate for each member of col_ef
    !-----------------------------------------------------------------------
    allocate(this%eflx_h2osfc_to_snow  (begc:endc))              ; this%eflx_h2osfc_to_snow  (:)   = nan
    allocate(this%eflx_snomelt         (begc:endc))              ; this%eflx_snomelt         (:)   = nan
    allocate(this%eflx_snomelt_r       (begc:endc))              ; this%eflx_snomelt_r       (:)   = nan
    allocate(this%eflx_snomelt_u       (begc:endc))              ; this%eflx_snomelt_u       (:)   = nan
    allocate(this%eflx_bot             (begc:endc))              ; this%eflx_bot             (:)   = nan
    allocate(this%eflx_fgr12           (begc:endc))              ; this%eflx_fgr12           (:)   = nan
    allocate(this%eflx_fgr             (begc:endc, 1:nlevgrnd))  ; this%eflx_fgr             (:,:) = nan
    allocate(this%eflx_building_heat   (begc:endc))              ; this%eflx_building_heat   (:)   = nan
    allocate(this%eflx_urban_ac        (begc:endc))              ; this%eflx_urban_ac        (:)   = nan
    allocate(this%eflx_urban_heat      (begc:endc))              ; this%eflx_urban_heat      (:)   = nan
    allocate(this%eflx_hs_h2osfc       (begc:endc))              ; this%eflx_hs_h2osfc       (:)   = nan
    allocate(this%eflx_hs_top_snow     (begc:endc))              ; this%eflx_hs_top_snow     (:)   = nan
    allocate(this%eflx_hs_soil         (begc:endc))              ; this%eflx_hs_soil         (:)   = nan
    allocate(this%eflx_sabg_lyr        (begc:endc, -nlevsno+1:1)); this%eflx_sabg_lyr        (:,:) = nan
    allocate(this%eflx_dhsdT           (begc:endc))              ; this%eflx_dhsdT           (:)   = nan
    allocate(this%htvp                 (begc:endc))              ; this%htvp                 (:)   = nan
    allocate(this%xmf                  (begc:endc))              ; this%xmf                  (:)   = nan
    allocate(this%xmf_h2osfc           (begc:endc))              ; this%xmf_h2osfc           (:)   = nan
    allocate(this%imelt                (begc:endc,-nlevsno+1:nlevgrnd))  ; this%imelt        (:,:) = huge(1)
    allocate(this%eflx_soil_grnd       (begc:endc))              ; this%eflx_soil_grnd       (:)   = nan
    allocate(this%eflx_rnet_soil       (begc:endc))              ; this%eflx_rnet_soil       (:)   = nan
    allocate(this%eflx_fgr0_soil       (begc:endc))              ; this%eflx_fgr0_soil       (:)   = nan
    allocate(this%eflx_fgr0_snow       (begc:endc))              ; this%eflx_fgr0_snow       (:)   = nan
    allocate(this%eflx_fgr0_h2osfc     (begc:endc))              ; this%eflx_fgr0_h2osfc     (:)   = nan
    allocate(this%errsoi               (begc:endc))              ; this%errsoi               (:)   = nan
    allocate(this%errseb               (begc:endc))              ; this%errseb               (:)   = nan
    allocate(this%errsol               (begc:endc))              ; this%errsol               (:)   = nan
    allocate(this%errlon               (begc:endc))              ; this%errlon               (:)   = nan

    !-----------------------------------------------------------------------
    ! initialize history fields for select members of col_ef
    !-----------------------------------------------------------------------
    this%eflx_snomelt(begc:endc) = spval
    call hist_addfld1d (fname='FSM',  units='W/m^2',  &
         avgflag='A', long_name='snow melt heat flux', &
         ptr_col=this%eflx_snomelt, c2l_scale_type='urbanf')

    this%eflx_snomelt_r(begc:endc) = spval
    call hist_addfld1d (fname='FSM_R',  units='W/m^2',  &
         avgflag='A', long_name='Rural snow melt heat flux', &
         ptr_col=this%eflx_snomelt_r, set_spec=spval)

    this%eflx_snomelt_u(begc:endc) = spval
    call hist_addfld1d (fname='FSM_U',  units='W/m^2',  &
         avgflag='A', long_name='Urban snow melt heat flux', &
         ptr_col=this%eflx_snomelt_u, c2l_scale_type='urbanf', set_nourb=spval)

    this%eflx_building_heat(begc:endc) = spval
    call hist_addfld1d (fname='BUILDHEAT', units='W/m^2',  &
         avgflag='A', long_name='heat flux from urban building interior to walls and roof', &
         ptr_col=this%eflx_building_heat, set_nourb=0._r8, c2l_scale_type='urbanf')

    this%eflx_urban_ac(begc:endc) = spval
    call hist_addfld1d (fname='URBAN_AC', units='W/m^2',  &
         avgflag='A', long_name='urban air conditioning flux', &
         ptr_col=this%eflx_urban_ac, set_nourb=0._r8, c2l_scale_type='urbanf')

    this%eflx_urban_heat(begc:endc) = spval
    call hist_addfld1d (fname='URBAN_HEAT', units='W/m^2',  &
         avgflag='A', long_name='urban heating flux', &
         ptr_col=this%eflx_urban_heat, set_nourb=0._r8, c2l_scale_type='urbanf')

    this%eflx_fgr12(begc:endc) = spval
    call hist_addfld1d (fname='FGR12',  units='W/m^2',  &
         avgflag='A', long_name='heat flux between soil layers 1 and 2', &
         ptr_col=this%eflx_fgr12, set_lake=spval)

    this%eflx_fgr(begc:endc,:) = spval
    call hist_addfld2d (fname='FGR_SOIL_R', units='watt/m^2', type2d='levgrnd', &
         avgflag='A', long_name='Rural downward heat flux at interface below each soil layer', &
         ptr_col=this%eflx_fgr, set_spec=spval, default='inactive')

    this%errsoi(begc:endc) = spval
    call hist_addfld1d (fname='ERRSOI',  units='W/m^2',  &
         avgflag='A', long_name='soil/lake energy conservation error', &
         ptr_col=this%errsoi)

    !-----------------------------------------------------------------------
    ! set cold-start initial values for select members of col_ef
    !-----------------------------------------------------------------------
    do c = begc, endc
       l = col_pp%landunit(c)
       if (lun_pp%urbpoi(l)) then
          this%eflx_building_heat(c) = 0._r8
          this%eflx_urban_ac(c)      = 0._r8
          this%eflx_urban_heat(c)    = 0._r8
       else
          this%eflx_building_heat(c) = 0._r8
          this%eflx_urban_ac(c)      = 0._r8
          this%eflx_urban_heat(c)    = 0._r8
       end if
    end do
    
  end subroutine col_ef_init
    
  !------------------------------------------------------------------------
  subroutine col_ef_restart(this, bounds, ncid, flag)
    ! 
    ! !DESCRIPTION:
    ! Read/Write column energy state information to/from restart file.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(column_energy_flux) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   
    character(len=*) , intent(in)    :: flag   
    !
    ! !LOCAL VARIABLES:
    logical :: readvar   ! determine if variable is on initial file
    !-----------------------------------------------------------------------

    call restartvar(ncid=ncid, flag=flag, varname='URBAN_AC', xtype=ncd_double,  dim1name='column', &
         long_name='urban air conditioning flux', units='watt/m^2', &
         interpinic_flag='interp', readvar=readvar, data=this%eflx_urban_ac)

    call restartvar(ncid=ncid, flag=flag, varname='URBAN_HEAT', xtype=ncd_double, dim1name='column', &
         long_name='urban heating flux', units='watt/m^2', &
         interpinic_flag='interp', readvar=readvar, data=this%eflx_urban_heat)

  end subroutine col_ef_restart
  
  !------------------------------------------------------------------------
  subroutine col_ef_clean(this)
    !
    ! !ARGUMENTS:
    class(column_energy_flux) :: this
    !------------------------------------------------------------------------
    
  end subroutine col_ef_clean
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean column water flux data structure
  !------------------------------------------------------------------------
  subroutine init_col_wf(this, begc, endc)
    !
    ! !ARGUMENTS:
    class(column_water_flux) :: this
    integer, intent(in) :: begc,endc
    !------------------------------------------------------------------------
    
  end subroutine init_col_wf
    
  !------------------------------------------------------------------------
  subroutine clean_col_wf(this)
    !
    ! !ARGUMENTS:
    class(column_water_flux) :: this
    !------------------------------------------------------------------------
    
  end subroutine clean_col_wf
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean column carbon flux data structure
  !------------------------------------------------------------------------
  subroutine init_col_cf(this, begc, endc)
    !
    ! !ARGUMENTS:
    class(column_carbon_flux) :: this
    integer, intent(in) :: begc,endc
    !------------------------------------------------------------------------
    
  end subroutine init_col_cf
    
  !------------------------------------------------------------------------
  subroutine clean_col_cf(this)
    !
    ! !ARGUMENTS:
    class(column_carbon_flux) :: this
    !------------------------------------------------------------------------
    
  end subroutine clean_col_cf
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean column nitrogen flux data structure
  !------------------------------------------------------------------------
  subroutine init_col_nf(this, begc, endc)
    !
    ! !ARGUMENTS:
    class(column_nitrogen_flux) :: this
    integer, intent(in) :: begc,endc
    !------------------------------------------------------------------------
    
  end subroutine init_col_nf
    
  !------------------------------------------------------------------------
  subroutine clean_col_nf(this)
    !
    ! !ARGUMENTS:
    class(column_nitrogen_flux) :: this
    !------------------------------------------------------------------------
    
  end subroutine clean_col_nf
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean column phosphorus flux data structure
  !------------------------------------------------------------------------
  subroutine init_col_pf(this, begc, endc)
    !
    ! !ARGUMENTS:
    class(column_phosphorus_flux) :: this
    integer, intent(in) :: begc,endc
    !------------------------------------------------------------------------
    
  end subroutine init_col_pf
    
  !------------------------------------------------------------------------
  subroutine clean_col_pf(this)
    !
    ! !ARGUMENTS:
    class(column_phosphorus_flux) :: this
    !------------------------------------------------------------------------
    
  end subroutine clean_col_pf
  
    !------------------------------------------------------------------------
    
end module ColumnDataType

