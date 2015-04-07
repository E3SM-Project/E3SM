module AerosolMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  use shr_kind_mod     , only : r8 => shr_kind_r8
  use shr_log_mod      , only : errMsg => shr_log_errMsg
  use shr_infnan_mod   , only : nan => shr_infnan_nan, assignment(=)
  use decompMod        , only : bounds_type
  use clm_varpar       , only : nlevsno, nlevgrnd 
  use clm_time_manager , only : get_step_size
  use atm2lndType      , only : atm2lnd_type
  use WaterfluxType    , only : waterflux_type
  use WaterstateType   , only : waterstate_type
  use ColumnType       , only : col               
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: AerosolMasses
  public :: AerosolFluxes
  !
  ! !PUBLIC DATA MEMBERS:
  real(r8), public, parameter :: snw_rds_min = 54.526_r8       ! minimum allowed snow effective radius (also "fresh snow" value) [microns
  !
  type, public :: aerosol_type
     real(r8), pointer, public  :: mss_bcpho_col(:,:)      ! mass of hydrophobic BC in snow (col,lyr)     [kg]
     real(r8), pointer, public  :: mss_bcphi_col(:,:)      ! mass of hydrophillic BC in snow (col,lyr)    [kg]
     real(r8), pointer, public  :: mss_bctot_col(:,:)      ! total mass of BC in snow (pho+phi) (col,lyr) [kg]
     real(r8), pointer, public  :: mss_bc_col_col(:)       ! column-integrated mass of total BC           [kg]
     real(r8), pointer, public  :: mss_bc_top_col(:)       ! top-layer mass of total BC                   [kg]

     real(r8), pointer, public  :: mss_ocpho_col(:,:)      ! mass of hydrophobic OC in snow (col,lyr)     [kg]
     real(r8), pointer, public  :: mss_ocphi_col(:,:)      ! mass of hydrophillic OC in snow (col,lyr)    [kg]
     real(r8), pointer, public  :: mss_octot_col(:,:)      ! total mass of OC in snow (pho+phi) (col,lyr) [kg]
     real(r8), pointer, public  :: mss_oc_col_col(:)       ! column-integrated mass of total OC           [kg]
     real(r8), pointer, public  :: mss_oc_top_col(:)       ! top-layer mass of total OC                   [kg]

     real(r8), pointer, public  :: mss_dst1_col(:,:)       ! mass of dust species 1 in snow (col,lyr)     [kg]
     real(r8), pointer, public  :: mss_dst2_col(:,:)       ! mass of dust species 2 in snow (col,lyr)     [kg]
     real(r8), pointer, public  :: mss_dst3_col(:,:)       ! mass of dust species 3 in snow (col,lyr)     [kg]
     real(r8), pointer, public  :: mss_dst4_col(:,:)       ! mass of dust species 4 in snow (col,lyr)     [kg]
     real(r8), pointer, public  :: mss_dsttot_col(:,:)     ! total mass of dust in snow (col,lyr)         [kg]
     real(r8), pointer, public  :: mss_dst_col_col(:)      ! column-integrated mass of dust in snow       [kg]
     real(r8), pointer, public  :: mss_dst_top_col(:)      ! top-layer mass of dust in snow               [kg]

     real(r8), pointer, public  :: mss_cnc_bcphi_col(:,:)  ! mass concentration of hydrophilic BC in snow (col,lyr) [kg/kg]
     real(r8), pointer, public  :: mss_cnc_bcpho_col(:,:)  ! mass concentration of hydrophilic BC in snow (col,lyr) [kg/kg]
     real(r8), pointer, public  :: mss_cnc_ocphi_col(:,:)  ! mass concentration of hydrophilic OC in snow (col,lyr) [kg/kg]
     real(r8), pointer, public  :: mss_cnc_ocpho_col(:,:)  ! mass concentration of hydrophilic OC in snow (col,lyr) [kg/kg]
     real(r8), pointer, public  :: mss_cnc_dst1_col(:,:)   ! mass concentration of dust species 1 in snow (col,lyr) [kg/kg]
     real(r8), pointer, public  :: mss_cnc_dst2_col(:,:)   ! mass concentration of dust species 2 in snow (col,lyr) [kg/kg]
     real(r8), pointer, public  :: mss_cnc_dst3_col(:,:)   ! mass concentration of dust species 3 in snow (col,lyr) [kg/kg]
     real(r8), pointer, public  :: mss_cnc_dst4_col(:,:)   ! mass concentration of dust species 4 in snow (col,lyr) [kg/kg]

     real(r8), pointer, private :: flx_dst_dep_dry1_col(:) ! dust species 1 dry   deposition on ground (positive definite) [kg/s]
     real(r8), pointer, private :: flx_dst_dep_wet1_col(:) ! dust species 1 wet   deposition on ground (positive definite) [kg/s]
     real(r8), pointer, private :: flx_dst_dep_dry2_col(:) ! dust species 2 dry   deposition on ground (positive definite) [kg/s]
     real(r8), pointer, private :: flx_dst_dep_wet2_col(:) ! dust species 2 wet   deposition on ground (positive definite) [kg/s]
     real(r8), pointer, private :: flx_dst_dep_dry3_col(:) ! dust species 3 dry   deposition on ground (positive definite) [kg/s]
     real(r8), pointer, private :: flx_dst_dep_wet3_col(:) ! dust species 3 wet   deposition on ground (positive definite) [kg/s]
     real(r8), pointer, private :: flx_dst_dep_dry4_col(:) ! dust species 4 dry   deposition on ground (positive definite) [kg/s]
     real(r8), pointer, private :: flx_dst_dep_wet4_col(:) ! dust species 4 wet   deposition on ground (positive definite) [kg/s]
     real(r8), pointer, private :: flx_dst_dep_col(:)      ! total (dry+wet) dust deposition on ground (positive definite) [kg/s]

     real(r8), pointer, private :: flx_bc_dep_dry_col(:)   ! dry (BCPHO+BCPHI) BC deposition on ground (positive definite) [kg/s]
     real(r8), pointer, private :: flx_bc_dep_wet_col(:)   ! wet (BCPHI) BC       deposition on ground (positive definite) [kg/s]
     real(r8), pointer, private :: flx_bc_dep_pho_col(:)   ! hydrophobic BC       deposition on ground (positive definite) [kg/s]
     real(r8), pointer, private :: flx_bc_dep_phi_col(:)   ! hydrophillic BC      deposition on ground (positive definite) [kg/s]
     real(r8), pointer, private :: flx_bc_dep_col(:)       ! total (dry+wet) BC   deposition on ground (positive definite) [kg/s]

     real(r8), pointer, private :: flx_oc_dep_dry_col(:)   ! dry (OCPHO+OCPHI) OC deposition on ground (positive definite) [kg/s]
     real(r8), pointer, private :: flx_oc_dep_wet_col(:)   ! wet (OCPHI) OC       deposition on ground (positive definite) [kg/s]
     real(r8), pointer, private :: flx_oc_dep_pho_col(:)   ! hydrophobic OC       deposition on ground (positive definite) [kg/s]
     real(r8), pointer, private :: flx_oc_dep_phi_col(:)   ! hydrophillic OC      deposition on ground (positive definite) [kg/s]
     real(r8), pointer, private :: flx_oc_dep_col(:)       ! total (dry+wet) OC   deposition on ground (positive definite) [kg/s]

   contains

     ! Public procedures
     procedure, public  :: Init         
     procedure, public  :: Restart      
     procedure, public  :: Reset 

     ! Private procedures
     procedure, private :: InitAllocate 
     procedure, private :: InitHistory  
     procedure, private :: InitCold     
       
  end type aerosol_type
  !-----------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(aerosol_type) :: this
    type(bounds_type), intent(in) :: bounds  

    call this%InitAllocate(bounds)
    call this%InitHistory(bounds)
    call this%InitCold(bounds)

  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !ARGUMENTS:
    class(aerosol_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begc, endc
    !---------------------------------------------------------------------

    begc = bounds%begc; endc= bounds%endc

    allocate(this%flx_dst_dep_dry1_col (begc:endc))              ; this%flx_dst_dep_dry1_col (:)   = nan
    allocate(this%flx_dst_dep_wet1_col (begc:endc))              ; this%flx_dst_dep_wet1_col (:)   = nan
    allocate(this%flx_dst_dep_dry2_col (begc:endc))              ; this%flx_dst_dep_dry2_col (:)   = nan
    allocate(this%flx_dst_dep_wet2_col (begc:endc))              ; this%flx_dst_dep_wet2_col (:)   = nan
    allocate(this%flx_dst_dep_dry3_col (begc:endc))              ; this%flx_dst_dep_dry3_col (:)   = nan
    allocate(this%flx_dst_dep_wet3_col (begc:endc))              ; this%flx_dst_dep_wet3_col (:)   = nan
    allocate(this%flx_dst_dep_dry4_col (begc:endc))              ; this%flx_dst_dep_dry4_col (:)   = nan
    allocate(this%flx_dst_dep_wet4_col (begc:endc))              ; this%flx_dst_dep_wet4_col (:)   = nan
    allocate(this%flx_dst_dep_col      (begc:endc))              ; this%flx_dst_dep_col      (:)   = nan

    allocate(this%flx_bc_dep_dry_col   (begc:endc))              ; this%flx_bc_dep_dry_col   (:)   = nan
    allocate(this%flx_bc_dep_wet_col   (begc:endc))              ; this%flx_bc_dep_wet_col   (:)   = nan
    allocate(this%flx_bc_dep_pho_col   (begc:endc))              ; this%flx_bc_dep_pho_col   (:)   = nan
    allocate(this%flx_bc_dep_phi_col   (begc:endc))              ; this%flx_bc_dep_phi_col   (:)   = nan
    allocate(this%flx_bc_dep_col       (begc:endc))              ; this%flx_bc_dep_col       (:)   = nan

    allocate(this%flx_oc_dep_dry_col   (begc:endc))              ; this%flx_oc_dep_dry_col   (:)   = nan
    allocate(this%flx_oc_dep_wet_col   (begc:endc))              ; this%flx_oc_dep_wet_col   (:)   = nan
    allocate(this%flx_oc_dep_pho_col   (begc:endc))              ; this%flx_oc_dep_pho_col   (:)   = nan
    allocate(this%flx_oc_dep_phi_col   (begc:endc))              ; this%flx_oc_dep_phi_col   (:)   = nan
    allocate(this%flx_oc_dep_col       (begc:endc))              ; this%flx_oc_dep_col       (:)   = nan

    allocate(this%mss_bcpho_col        (begc:endc,-nlevsno+1:0)) ; this%mss_bcpho_col        (:,:) = nan
    allocate(this%mss_bcphi_col        (begc:endc,-nlevsno+1:0)) ; this%mss_bcphi_col        (:,:) = nan
    allocate(this%mss_bctot_col        (begc:endc,-nlevsno+1:0)) ; this%mss_bctot_col        (:,:) = nan
    allocate(this%mss_bc_col_col       (begc:endc))              ; this%mss_bc_col_col       (:)   = nan
    allocate(this%mss_bc_top_col       (begc:endc))              ; this%mss_bc_top_col       (:)   = nan

    allocate(this%mss_ocpho_col        (begc:endc,-nlevsno+1:0)) ; this%mss_ocpho_col        (:,:) = nan
    allocate(this%mss_ocphi_col        (begc:endc,-nlevsno+1:0)) ; this%mss_ocphi_col        (:,:) = nan
    allocate(this%mss_octot_col        (begc:endc,-nlevsno+1:0)) ; this%mss_octot_col        (:,:) = nan
    allocate(this%mss_oc_col_col       (begc:endc))              ; this%mss_oc_col_col       (:)   = nan
    allocate(this%mss_oc_top_col       (begc:endc))              ; this%mss_oc_top_col       (:)   = nan

    allocate(this%mss_dst1_col         (begc:endc,-nlevsno+1:0)) ; this%mss_dst1_col         (:,:) = nan
    allocate(this%mss_dst2_col         (begc:endc,-nlevsno+1:0)) ; this%mss_dst2_col         (:,:) = nan
    allocate(this%mss_dst3_col         (begc:endc,-nlevsno+1:0)) ; this%mss_dst3_col         (:,:) = nan
    allocate(this%mss_dst4_col         (begc:endc,-nlevsno+1:0)) ; this%mss_dst4_col         (:,:) = nan
    allocate(this%mss_dsttot_col       (begc:endc,-nlevsno+1:0)) ; this%mss_dsttot_col       (:,:) = nan
    allocate(this%mss_dst_col_col      (begc:endc))              ; this%mss_dst_col_col      (:)   = nan
    allocate(this%mss_dst_top_col      (begc:endc))              ; this%mss_dst_top_col      (:)   = nan

    allocate(this%mss_cnc_bcphi_col    (begc:endc,-nlevsno+1:0)) ; this%mss_cnc_bcphi_col    (:,:) = nan
    allocate(this%mss_cnc_bcpho_col    (begc:endc,-nlevsno+1:0)) ; this%mss_cnc_bcpho_col    (:,:) = nan
    allocate(this%mss_cnc_ocphi_col    (begc:endc,-nlevsno+1:0)) ; this%mss_cnc_ocphi_col    (:,:) = nan
    allocate(this%mss_cnc_ocpho_col    (begc:endc,-nlevsno+1:0)) ; this%mss_cnc_ocpho_col    (:,:) = nan
    allocate(this%mss_cnc_dst1_col     (begc:endc,-nlevsno+1:0)) ; this%mss_cnc_dst1_col     (:,:) = nan
    allocate(this%mss_cnc_dst2_col     (begc:endc,-nlevsno+1:0)) ; this%mss_cnc_dst2_col     (:,:) = nan
    allocate(this%mss_cnc_dst3_col     (begc:endc,-nlevsno+1:0)) ; this%mss_cnc_dst3_col     (:,:) = nan
    allocate(this%mss_cnc_dst4_col     (begc:endc,-nlevsno+1:0)) ; this%mss_cnc_dst4_col     (:,:) = nan

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! History fields initialization
    !
    ! !USES:
    use shr_infnan_mod, only: nan => shr_infnan_nan, assignment(=)
    use clm_varcon    , only: spval
    use clm_varpar    , only: nlevsno 
    use histFileMod   , only: hist_addfld1d, hist_addfld2d
    use histFileMod   , only: no_snow_normal, no_snow_zero
    !
    ! !ARGUMENTS:
    class(aerosol_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begc, endc
    real(r8), pointer :: data2dptr(:,:) ! temp. pointers for slicing larger arrays
    !---------------------------------------------------------------------

    begc = bounds%begc; endc= bounds%endc

    this%flx_dst_dep_col(begc:endc) = spval
    call hist_addfld1d (fname='DSTDEP', units='kg/m^2/s', &
         avgflag='A', long_name='total dust deposition (dry+wet) from atmosphere', &
         ptr_col=this%flx_dst_dep_col, set_urb=spval)

    this%flx_bc_dep_col(begc:endc) = spval
    call hist_addfld1d (fname='BCDEP', units='kg/m^2/s', &
         avgflag='A', long_name='total BC deposition (dry+wet) from atmosphere', &
         ptr_col=this%flx_bc_dep_col, set_urb=spval)

    this%flx_oc_dep_col(begc:endc) = spval    
    call hist_addfld1d (fname='OCDEP', units='kg/m^2/s', &
         avgflag='A', long_name='total OC deposition (dry+wet) from atmosphere', &
         ptr_col=this%flx_oc_dep_col, set_urb=spval)

    this%mss_bc_col_col(begc:endc) = spval
    call hist_addfld1d (fname='SNOBCMCL', units='kg/m2', &
         avgflag='A', long_name='mass of BC in snow column', &
         ptr_col=this%mss_bc_col_col, set_urb=spval)

    this%mss_bc_top_col(begc:endc) = spval
    call hist_addfld1d (fname='SNOBCMSL', units='kg/m2', &
         avgflag='A', long_name='mass of BC in top snow layer', &
         ptr_col=this%mss_bc_top_col, set_urb=spval)

    this%mss_oc_col_col(begc:endc) = spval
    call hist_addfld1d (fname='SNOOCMCL', units='kg/m2', &
         avgflag='A', long_name='mass of OC in snow column', &
         ptr_col=this%mss_oc_col_col, set_urb=spval)

    this%mss_oc_top_col(begc:endc) = spval
    call hist_addfld1d (fname='SNOOCMSL', units='kg/m2', &
         avgflag='A', long_name='mass of OC in top snow layer', &
         ptr_col=this%mss_oc_top_col, set_urb=spval)

    this%mss_dst_col_col(begc:endc) = spval
    call hist_addfld1d (fname='SNODSTMCL', units='kg/m2', &
         avgflag='A', long_name='mass of dust in snow column', &
         ptr_col=this%mss_dst_col_col, set_urb=spval)

    this%mss_dst_top_col(begc:endc) = spval
    call hist_addfld1d (fname='SNODSTMSL', units='kg/m2', &
         avgflag='A', long_name='mass of dust in top snow layer', &
         ptr_col=this%mss_dst_top_col, set_urb=spval)

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(aerosol_type) :: this
    type(bounds_type) , intent(in) :: bounds                                   
    !
    ! !LOCAL VARIABLES:
    integer  :: c       ! index
    !-----------------------------------------------------------------------

    do c = bounds%begc,bounds%endc
       this%mss_cnc_bcphi_col(c,:) = 0._r8
       this%mss_cnc_bcpho_col(c,:) = 0._r8
       this%mss_cnc_ocphi_col(c,:) = 0._r8
       this%mss_cnc_ocpho_col(c,:) = 0._r8
       this%mss_cnc_dst1_col(c,:)  = 0._r8
       this%mss_cnc_dst2_col(c,:)  = 0._r8
       this%mss_cnc_dst3_col(c,:)  = 0._r8
       this%mss_cnc_dst4_col(c,:)  = 0._r8

       this%mss_bctot_col(c,:)     = 0._r8
       this%mss_bcpho_col(c,:)     = 0._r8
       this%mss_bcphi_col(c,:)     = 0._r8

       this%mss_octot_col(c,:)     = 0._r8
       this%mss_ocpho_col(c,:)     = 0._r8
       this%mss_ocphi_col(c,:)     = 0._r8

       this%mss_dst1_col(c,:)      = 0._r8
       this%mss_dst2_col(c,:)      = 0._r8
       this%mss_dst3_col(c,:)      = 0._r8
       this%mss_dst4_col(c,:)      = 0._r8
       this%mss_dsttot_col(c,:)    = 0._r8
    end do

  end subroutine InitCold

  !------------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag, &
       h2osoi_ice_col, h2osoi_liq_col)
    ! 
    ! !DESCRIPTION:
    ! Read/Write module information to/from restart file.
    !
    ! !USES:
    use clm_varpar , only : nlevsno, nlevsoi
    use clm_varcon , only : spval
    use clm_varctl , only : iulog  
    use clm_varpar , only : nlevsno
    use spmdMod    , only : masterproc
    use ncdio_pio  , only : file_desc_t, ncd_defvar, ncd_io, ncd_double, ncd_int, ncd_inqvdlen
    use restUtilMod
    !
    ! !ARGUMENTS:
    class(aerosol_type) :: this
    type(bounds_type)   , intent(in)    :: bounds                                   
    type(file_desc_t)   , intent(inout) :: ncid                                         ! netcdf id
    character(len=*)    , intent(in)    :: flag                                         ! 'read' or 'write'
    real(r8)            , intent(in)    :: h2osoi_ice_col( bounds%begc: , -nlevsno+1: ) ! ice content (col,lyr) [kg/m2]
    real(r8)            , intent(in)    :: h2osoi_liq_col( bounds%begc: , -nlevsno+1: ) ! liquid water content (col,lyr) [kg/m2]
    !
    ! !LOCAL VARIABLES:
    integer :: j,c ! indices
    logical :: readvar      ! determine if variable is on initial file
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(h2osoi_ice_col) == (/bounds%endc,nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(h2osoi_liq_col) == (/bounds%endc,nlevgrnd/)), errMsg(__FILE__, __LINE__))

    call restartvar(ncid=ncid, flag=flag, varname='mss_bcpho', xtype=ncd_double,  &
         dim1name='column', dim2name='levsno', switchdim=.true., lowerb2=-nlevsno+1, upperb2=0, &
         long_name='snow layer hydrophobic black carbon mass', units='kg m-2', &
         interpinic_flag='interp', readvar=readvar, data=this%mss_bcpho_col)
    if (flag == 'read' .and. .not. readvar) then
       ! initial run, not restart: initialize mss_bcpho to zero
       this%mss_bcpho_col(bounds%begc:bounds%endc,-nlevsno+1:0) = 0._r8
    end if

    call restartvar(ncid=ncid, flag=flag, varname='mss_bcphi', xtype=ncd_double,  &
         dim1name='column', dim2name='levsno', switchdim=.true., lowerb2=-nlevsno+1, upperb2=0, &
         long_name='snow layer hydrophilic black carbon mass', units='kg m-2', &
         interpinic_flag='interp', readvar=readvar, data=this%mss_bcphi_col)
    if (flag == 'read' .and. .not. readvar) then
       ! initial run, not restart: initialize mss_bcphi to zero
       this%mss_bcphi_col(bounds%begc:bounds%endc,-nlevsno+1:0) = 0._r8
    end if

    call restartvar(ncid=ncid, flag=flag, varname='mss_ocpho', xtype=ncd_double,  &
         dim1name='column', dim2name='levsno', switchdim=.true., lowerb2=-nlevsno+1, upperb2=0, &
         long_name='snow layer hydrophobic organic carbon mass', units='kg m-2', &
         interpinic_flag='interp', readvar=readvar, data=this%mss_ocpho_col)
    if (flag == 'read' .and. .not. readvar) then
       ! initial run, not restart: initialize mss_ocpho to zero
       this%mss_ocpho_col(bounds%begc:bounds%endc,-nlevsno+1:0) = 0._r8
    end if

    call restartvar(ncid=ncid, flag=flag, varname='mss_ocphi', xtype=ncd_double,  &
         dim1name='column', dim2name='levsno', switchdim=.true., lowerb2=-nlevsno+1, upperb2=0, &
         long_name='snow layer hydrophilic organic carbon mass', units='kg m-2', &
         interpinic_flag='interp', readvar=readvar, data=this%mss_ocphi_col)
    if (flag == 'read' .and. .not. readvar) then
       ! initial run, not restart: initialize mss_ocphi to zero
       this%mss_ocphi_col(bounds%begc:bounds%endc,-nlevsno+1:0) = 0._r8
    end if

    call restartvar(ncid=ncid, flag=flag, varname='mss_dst1', xtype=ncd_double,  &
         dim1name='column', dim2name='levsno', switchdim=.true., lowerb2=-nlevsno+1, upperb2=0, &
         long_name='snow layer dust species 1 mass', units='kg m-2', &
         interpinic_flag='interp', readvar=readvar, data=this%mss_dst1_col)
    if (flag == 'read' .and. .not. readvar) then
       ! initial run, not restart: initialize mss_dst1 to zero
       this%mss_dst1_col(bounds%begc:bounds%endc,-nlevsno+1:0) = 0._r8
    end if
    
    call restartvar(ncid=ncid, flag=flag, varname='mss_dst2', xtype=ncd_double,  &
         dim1name='column', dim2name='levsno', switchdim=.true., lowerb2=-nlevsno+1, upperb2=0, &
         long_name='snow layer dust species 2 mass', units='kg m-2', &
         interpinic_flag='interp', readvar=readvar, data=this%mss_dst2_col)
    if (flag == 'read' .and. .not. readvar) then
       ! initial run, not restart: initialize mss_dst2 to zero
       this%mss_dst2_col(bounds%begc:bounds%endc,-nlevsno+1:0) = 0._r8
    endif

    call restartvar(ncid=ncid, flag=flag, varname='mss_dst3', xtype=ncd_double,  &
         dim1name='column', dim2name='levsno', switchdim=.true., lowerb2=-nlevsno+1, upperb2=0, &
         long_name='snow layer dust species 3 mass', units='kg m-2', &
         interpinic_flag='interp', readvar=readvar,  data=this%mss_dst3_col)
    if (flag == 'read' .and. .not. readvar) then
       ! initial run, not restart: initialize mss_dst3 to zero
       this%mss_dst3_col(bounds%begc:bounds%endc,-nlevsno+1:0) = 0._r8
    endif

    call restartvar(ncid=ncid, flag=flag, varname='mss_dst4', xtype=ncd_double,  &
         dim1name='column', dim2name='levsno', switchdim=.true., lowerb2=-nlevsno+1, upperb2=0, &
         long_name='snow layer dust species 4 mass', units='kg m-2', &
         interpinic_flag='interp', readvar=readvar, data=this%mss_dst4_col)
    if (flag == 'read' .and. .not. readvar) then
       ! initial run, not restart: initialize mss_dst4 to zero
       this%mss_dst4_col(bounds%begc:bounds%endc,-nlevsno+1:0) = 0._r8
    end if

    ! initialize other variables that are derived from those stored in the restart buffer (SNICAR variables) 
    if (flag == 'read' ) then
       do j = -nlevsno+1,0
          do c = bounds%begc, bounds%endc
             ! mass concentrations of aerosols in snow
             if (h2osoi_ice_col(c,j) + h2osoi_liq_col(c,j) > 0._r8) then
                this%mss_cnc_bcpho_col(c,j) = this%mss_bcpho_col(c,j) / (h2osoi_ice_col(c,j)+h2osoi_liq_col(c,j))
                this%mss_cnc_bcphi_col(c,j) = this%mss_bcphi_col(c,j) / (h2osoi_ice_col(c,j)+h2osoi_liq_col(c,j))
                this%mss_cnc_ocpho_col(c,j) = this%mss_ocpho_col(c,j) / (h2osoi_ice_col(c,j)+h2osoi_liq_col(c,j))
                this%mss_cnc_ocphi_col(c,j) = this%mss_ocphi_col(c,j) / (h2osoi_ice_col(c,j)+h2osoi_liq_col(c,j))

                this%mss_cnc_dst1_col(c,j) = this%mss_dst1_col(c,j) / (h2osoi_ice_col(c,j)+h2osoi_liq_col(c,j))
                this%mss_cnc_dst2_col(c,j) = this%mss_dst2_col(c,j) / (h2osoi_ice_col(c,j)+h2osoi_liq_col(c,j))
                this%mss_cnc_dst3_col(c,j) = this%mss_dst3_col(c,j) / (h2osoi_ice_col(c,j)+h2osoi_liq_col(c,j))
                this%mss_cnc_dst4_col(c,j) = this%mss_dst4_col(c,j) / (h2osoi_ice_col(c,j)+h2osoi_liq_col(c,j))
             else
                this%mss_cnc_bcpho_col(c,j) = 0._r8
                this%mss_cnc_bcphi_col(c,j) = 0._r8
                this%mss_cnc_ocpho_col(c,j) = 0._r8
                this%mss_cnc_ocphi_col(c,j) = 0._r8

                this%mss_cnc_dst1_col(c,j) = 0._r8
                this%mss_cnc_dst2_col(c,j) = 0._r8
                this%mss_cnc_dst3_col(c,j) = 0._r8
                this%mss_cnc_dst4_col(c,j) = 0._r8
             endif
          enddo
       enddo
    endif

  end subroutine Restart

  !-----------------------------------------------------------------------
  subroutine Reset(this, column)
    !
    ! !DESCRIPTION:
    ! Intitialize SNICAR variables for fresh snow column
    !
    ! !ARGUMENTS:
    class(aerosol_type)             :: this
    integer           , intent(in)  :: column     ! column index
    !-----------------------------------------------------------------------

    this%mss_bcpho_col(column,:)  = 0._r8
    this%mss_bcphi_col(column,:)  = 0._r8
    this%mss_bctot_col(column,:)  = 0._r8
    this%mss_bc_col_col(column)   = 0._r8
    this%mss_bc_top_col(column)   = 0._r8

    this%mss_ocpho_col(column,:)  = 0._r8
    this%mss_ocphi_col(column,:)  = 0._r8
    this%mss_octot_col(column,:)  = 0._r8
    this%mss_oc_col_col(column)   = 0._r8
    this%mss_oc_top_col(column)   = 0._r8

    this%mss_dst1_col(column,:)   = 0._r8
    this%mss_dst2_col(column,:)   = 0._r8
    this%mss_dst3_col(column,:)   = 0._r8
    this%mss_dst4_col(column,:)   = 0._r8
    this%mss_dsttot_col(column,:) = 0._r8
    this%mss_dst_col_col(column)  = 0._r8
    this%mss_dst_top_col(column)  = 0._r8

  end subroutine Reset

  !-----------------------------------------------------------------------
  subroutine AerosolMasses(bounds, num_on, filter_on, num_off, filter_off, &
       waterflux_inst, waterstate_inst, aerosol_inst)
    !
    ! !DESCRIPTION:
    ! Calculate column-integrated aerosol masses, and
    ! mass concentrations for radiative calculations and output
    ! (based on new snow level state, after SnowFilter is rebuilt.
    ! NEEDS TO BE AFTER SnowFiler is rebuilt in Hydrology2, otherwise there 
    ! can be zero snow layers but an active column in filter)
    !
    ! !ARGUMENTS:
    type(bounds_type)     , intent(in )   :: bounds
    integer               , intent(in)    :: num_on         ! number of column filter-ON points
    integer               , intent(in)    :: filter_on(:)   ! column filter for filter-ON points
    integer               , intent(in)    :: num_off        ! number of column non filter-OFF points
    integer               , intent(in)    :: filter_off(:)  ! column filter for filter-OFF points
    type(waterflux_type)  , intent(in)    :: waterflux_inst 
    type(waterstate_type) , intent(inout) :: waterstate_inst
    type(aerosol_type)    , intent(inout) :: aerosol_inst
    !
    ! !LOCAL VARIABLES:
    real(r8) :: dtime           ! land model time step (sec)
    integer  :: g,l,c,j,fc      ! indices
    real(r8) :: snowmass        ! liquid+ice snow mass in a layer [kg/m2]
    real(r8) :: snowcap_scl_fct ! temporary factor used to correct for snow capping
    !-----------------------------------------------------------------------

    associate(                                                & 
         snl           => col%snl                           , & ! Input:  [integer  (:)   ]  number of snow layers                    

         do_capsnow    => waterstate_inst%do_capsnow_col    , & ! Input:  [logical  (:)   ]  true => do snow capping                  
         h2osoi_ice    => waterstate_inst%h2osoi_ice_col    , & ! Input:  [real(r8) (:,:) ]  ice lens (kg/m2)                      
         h2osoi_liq    => waterstate_inst%h2osoi_liq_col    , & ! Input:  [real(r8) (:,:) ]  liquid water (kg/m2)                  
         qflx_snwcp_ice=> waterflux_inst%qflx_snwcp_ice_col , & ! Input:  [real(r8) (:)   ]  excess snowfall due to snow capping (mm H2O /s) [+]          

         h2osno_top    => waterstate_inst%h2osno_top_col    , & ! Output: [real(r8) (:)   |  top-layer mass of snow  [kg]
         snw_rds       => waterstate_inst%snw_rds_col       , & ! Output: [real(r8) (:,:) ]  effective snow grain radius (col,lyr) [microns, m^-6] 

         mss_bcpho     => aerosol_inst%mss_bcpho_col        , & ! Output: [real(r8) (:,:) ]  mass of hydrophobic BC in snow (col,lyr) [kg]
         mss_bcphi     => aerosol_inst%mss_bcphi_col        , & ! Output: [real(r8) (:,:) ]  mass of hydrophillic BC in snow (col,lyr) [kg]
         mss_bctot     => aerosol_inst%mss_bctot_col        , & ! Output: [real(r8) (:,:) ]  total mass of BC (pho+phi) (col,lyr) [kg]
         mss_bc_col    => aerosol_inst%mss_bc_col_col       , & ! Output: [real(r8) (:)   ]  total mass of BC in snow column (col) [kg]
         mss_bc_top    => aerosol_inst%mss_bc_top_col       , & ! Output: [real(r8) (:)   ]  total mass of BC in top snow layer (col) [kg]
         mss_ocpho     => aerosol_inst%mss_ocpho_col        , & ! Output: [real(r8) (:,:) ]  mass of hydrophobic OC in snow (col,lyr) [kg]
         mss_ocphi     => aerosol_inst%mss_ocphi_col        , & ! Output: [real(r8) (:,:) ]  mass of hydrophillic OC in snow (col,lyr) [kg]
         mss_octot     => aerosol_inst%mss_octot_col        , & ! Output: [real(r8) (:,:) ]  total mass of OC (pho+phi) (col,lyr) [kg]
         mss_oc_col    => aerosol_inst%mss_oc_col_col       , & ! Output: [real(r8) (:)   ]  total mass of OC in snow column (col) [kg]
         mss_oc_top    => aerosol_inst%mss_oc_top_col       , & ! Output: [real(r8) (:)   ]  total mass of OC in top snow layer (col) [kg]
         mss_dst1      => aerosol_inst%mss_dst1_col         , & ! Output: [real(r8) (:,:) ]  mass of dust species 1 in snow (col,lyr) [kg]
         mss_dst2      => aerosol_inst%mss_dst2_col         , & ! Output: [real(r8) (:,:) ]  mass of dust species 2 in snow (col,lyr) [kg]
         mss_dst3      => aerosol_inst%mss_dst3_col         , & ! Output: [real(r8) (:,:) ]  mass of dust species 3 in snow (col,lyr) [kg]
         mss_dst4      => aerosol_inst%mss_dst4_col         , & ! Output: [real(r8) (:,:) ]  mass of dust species 4 in snow (col,lyr) [kg]
         mss_dsttot    => aerosol_inst%mss_dsttot_col       , & ! Output: [real(r8) (:,:) ]  total mass of dust in snow (col,lyr) [kg]
         mss_dst_col   => aerosol_inst%mss_dst_col_col      , & ! Output: [real(r8) (:)   ]  total mass of dust in snow column (col) [kg]
         mss_dst_top   => aerosol_inst%mss_dst_top_col      , & ! Output: [real(r8) (:)   ]  total mass of dust in top snow layer (col) [kg]
         mss_cnc_bcphi => aerosol_inst%mss_cnc_bcphi_col    , & ! Output: [real(r8) (:,:) ]  mass concentration of BC species 1 (col,lyr) [kg/kg]
         mss_cnc_bcpho => aerosol_inst%mss_cnc_bcpho_col    , & ! Output: [real(r8) (:,:) ]  mass concentration of BC species 2 (col,lyr) [kg/kg]
         mss_cnc_ocphi => aerosol_inst%mss_cnc_ocphi_col    , & ! Output: [real(r8) (:,:) ]  mass concentration of OC species 1 (col,lyr) [kg/kg]
         mss_cnc_ocpho => aerosol_inst%mss_cnc_ocpho_col    , & ! Output: [real(r8) (:,:) ]  mass concentration of OC species 2 (col,lyr) [kg/kg]
         mss_cnc_dst1  => aerosol_inst%mss_cnc_dst1_col     , & ! Output: [real(r8) (:,:) ]  mass concentration of dust species 1 (col,lyr) [kg/kg]
         mss_cnc_dst2  => aerosol_inst%mss_cnc_dst2_col     , & ! Output: [real(r8) (:,:) ]  mass concentration of dust species 2 (col,lyr) [kg/kg]
         mss_cnc_dst3  => aerosol_inst%mss_cnc_dst3_col     , & ! Output: [real(r8) (:,:) ]  mass concentration of dust species 3 (col,lyr) [kg/kg]
         mss_cnc_dst4  => aerosol_inst%mss_cnc_dst4_col       & ! Output: [real(r8) (:,:) ]  mass concentration of dust species 4 (col,lyr) [kg/kg]
         )

      dtime = get_step_size()

      do fc = 1, num_on
         c = filter_on(fc)

         ! Zero column-integrated aerosol mass before summation
         mss_bc_col(c)  = 0._r8
         mss_oc_col(c)  = 0._r8
         mss_dst_col(c) = 0._r8

         do j = -nlevsno+1, 0

            ! layer mass of snow:
            snowmass = h2osoi_ice(c,j) + h2osoi_liq(c,j)

            ! Correct the top layer aerosol mass to account for snow capping. 
            ! This approach conserves the aerosol mass concentration
            ! (but not the aerosol amss) when snow-capping is invoked

            if (j == snl(c)+1) then
               if (do_capsnow(c)) then 

                  snowcap_scl_fct = snowmass / (snowmass + (qflx_snwcp_ice(c)*dtime))

                  mss_bcpho(c,j) = mss_bcpho(c,j)*snowcap_scl_fct
                  mss_bcphi(c,j) = mss_bcphi(c,j)*snowcap_scl_fct
                  mss_ocpho(c,j) = mss_ocpho(c,j)*snowcap_scl_fct
                  mss_ocphi(c,j) = mss_ocphi(c,j)*snowcap_scl_fct

                  mss_dst1(c,j)  = mss_dst1(c,j)*snowcap_scl_fct
                  mss_dst2(c,j)  = mss_dst2(c,j)*snowcap_scl_fct
                  mss_dst3(c,j)  = mss_dst3(c,j)*snowcap_scl_fct
                  mss_dst4(c,j)  = mss_dst4(c,j)*snowcap_scl_fct
               endif
            endif

            if (j >= snl(c)+1) then

               mss_bctot(c,j)     = mss_bcpho(c,j) + mss_bcphi(c,j)
               mss_bc_col(c)      = mss_bc_col(c)  + mss_bctot(c,j)
               mss_cnc_bcphi(c,j) = mss_bcphi(c,j) / snowmass
               mss_cnc_bcpho(c,j) = mss_bcpho(c,j) / snowmass

               mss_octot(c,j)     = mss_ocpho(c,j) + mss_ocphi(c,j)
               mss_oc_col(c)      = mss_oc_col(c)  + mss_octot(c,j)
               mss_cnc_ocphi(c,j) = mss_ocphi(c,j) / snowmass
               mss_cnc_ocpho(c,j) = mss_ocpho(c,j) / snowmass

               mss_dsttot(c,j)    = mss_dst1(c,j)  + mss_dst2(c,j) + mss_dst3(c,j) + mss_dst4(c,j)
               mss_dst_col(c)     = mss_dst_col(c) + mss_dsttot(c,j)
               mss_cnc_dst1(c,j)  = mss_dst1(c,j)  / snowmass
               mss_cnc_dst2(c,j)  = mss_dst2(c,j)  / snowmass
               mss_cnc_dst3(c,j)  = mss_dst3(c,j)  / snowmass
               mss_cnc_dst4(c,j)  = mss_dst4(c,j)  / snowmass

            else
               !set variables of empty snow layers to zero
               snw_rds(c,j)       = 0._r8

               mss_bcpho(c,j)     = 0._r8
               mss_bcphi(c,j)     = 0._r8
               mss_bctot(c,j)     = 0._r8
               mss_cnc_bcphi(c,j) = 0._r8
               mss_cnc_bcpho(c,j) = 0._r8

               mss_ocpho(c,j)     = 0._r8
               mss_ocphi(c,j)     = 0._r8
               mss_octot(c,j)     = 0._r8
               mss_cnc_ocphi(c,j) = 0._r8
               mss_cnc_ocpho(c,j) = 0._r8

               mss_dst1(c,j)      = 0._r8
               mss_dst2(c,j)      = 0._r8
               mss_dst3(c,j)      = 0._r8
               mss_dst4(c,j)      = 0._r8
               mss_dsttot(c,j)    = 0._r8
               mss_cnc_dst1(c,j)  = 0._r8
               mss_cnc_dst2(c,j)  = 0._r8
               mss_cnc_dst3(c,j)  = 0._r8
               mss_cnc_dst4(c,j)  = 0._r8
            endif
         enddo

         ! top-layer diagnostics
         h2osno_top(c)  = h2osoi_ice(c,snl(c)+1) + h2osoi_liq(c,snl(c)+1) !TODO MV - is this correct to be placed here???
         mss_bc_top(c)  = mss_bctot(c,snl(c)+1)
         mss_oc_top(c)  = mss_octot(c,snl(c)+1)
         mss_dst_top(c) = mss_dsttot(c,snl(c)+1)
      enddo

      ! Zero mass variables in columns without snow

      do fc = 1, num_off
         c = filter_off(fc)

         mss_bc_top(c)      = 0._r8
         mss_bc_col(c)      = 0._r8    
         mss_bcpho(c,:)     = 0._r8
         mss_bcphi(c,:)     = 0._r8
         mss_bctot(c,:)     = 0._r8
         mss_cnc_bcphi(c,:) = 0._r8
         mss_cnc_bcpho(c,:) = 0._r8

         mss_oc_top(c)      = 0._r8
         mss_oc_col(c)      = 0._r8    
         mss_ocpho(c,:)     = 0._r8
         mss_ocphi(c,:)     = 0._r8
         mss_octot(c,:)     = 0._r8
         mss_cnc_ocphi(c,:) = 0._r8
         mss_cnc_ocpho(c,:) = 0._r8

         mss_dst_top(c)     = 0._r8
         mss_dst_col(c)     = 0._r8
         mss_dst1(c,:)      = 0._r8
         mss_dst2(c,:)      = 0._r8
         mss_dst3(c,:)      = 0._r8
         mss_dst4(c,:)      = 0._r8
         mss_dsttot(c,:)    = 0._r8

         mss_cnc_dst1(c,:)  = 0._r8
         mss_cnc_dst2(c,:)  = 0._r8
         mss_cnc_dst3(c,:)  = 0._r8
         mss_cnc_dst4(c,:)  = 0._r8

      enddo

    end associate

  end subroutine AerosolMasses

  !-----------------------------------------------------------------------
  subroutine AerosolFluxes(bounds, num_snowc, filter_snowc, &
       atm2lnd_inst, aerosol_inst)
    !
    ! !DESCRIPTION:
    ! Compute aerosol fluxes through snowpack and aerosol deposition fluxes into top layere
    !
    ! !ARGUMENTS:
    type(bounds_type)  , intent(in)    :: bounds
    integer            , intent(in)    :: num_snowc       ! number of snow points in column filter
    integer            , intent(in)    :: filter_snowc(:) ! column filter for snow points
    type(atm2lnd_type) , intent(in)    :: atm2lnd_inst
    type(aerosol_type) , intent(inout) :: aerosol_inst
    !
    ! !LOCAL VARIABLES:
    real(r8) :: dtime      ! land model time step (sec)
    integer  :: c,g,j,fc
    !-----------------------------------------------------------------------

    associate(                                                   & 
         snl              => col%snl                           , & ! Input:  [integer  (:)   ] number of snow layers                     

         forc_aer         => atm2lnd_inst%forc_aer_grc         , & ! Input:  [real(r8) (:,:) ] aerosol deposition from atmosphere model (grd,aer) [kg m-1 s-1]

         mss_bcphi        => aerosol_inst%mss_bcphi_col        , & ! Output: [real(r8) (:,:) ] hydrophillic BC mass in snow (col,lyr) [kg]
         mss_bcpho        => aerosol_inst%mss_bcpho_col        , & ! Output: [real(r8) (:,:) ] hydrophobic  BC mass in snow (col,lyr) [kg]
         mss_ocphi        => aerosol_inst%mss_ocphi_col        , & ! Output: [real(r8) (:,:) ] hydrophillic OC mass in snow (col,lyr) [kg]
         mss_ocpho        => aerosol_inst%mss_ocpho_col        , & ! Output: [real(r8) (:,:) ] hydrophobic  OC mass in snow (col,lyr) [kg]
         mss_dst1         => aerosol_inst%mss_dst1_col         , & ! Output: [real(r8) (:,:) ] mass of dust species 1 in snow (col,lyr) [kg]
         mss_dst2         => aerosol_inst%mss_dst2_col         , & ! Output: [real(r8) (:,:) ] mass of dust species 2 in snow (col,lyr) [kg]
         mss_dst3         => aerosol_inst%mss_dst3_col         , & ! Output: [real(r8) (:,:) ] mass of dust species 3 in snow (col,lyr) [kg]
         mss_dst4         => aerosol_inst%mss_dst4_col         , & ! Output: [real(r8) (:,:) ] mass of dust species 4 in snow (col,lyr) [kg]

         flx_bc_dep       => aerosol_inst%flx_bc_dep_col       , & ! Output: [real(r8) (:)   ] total BC deposition (col) [kg m-2 s-1]  
         flx_bc_dep_wet   => aerosol_inst%flx_bc_dep_wet_col   , & ! Output: [real(r8) (:)   ] wet BC deposition (col) [kg m-2 s-1]    
         flx_bc_dep_dry   => aerosol_inst%flx_bc_dep_dry_col   , & ! Output: [real(r8) (:)   ] dry BC deposition (col) [kg m-2 s-1]    
         flx_bc_dep_phi   => aerosol_inst%flx_bc_dep_phi_col   , & ! Output: [real(r8) (:)   ] hydrophillic BC deposition (col) [kg m-1 s-1]
         flx_bc_dep_pho   => aerosol_inst%flx_bc_dep_pho_col   , & ! Output: [real(r8) (:)   ] hydrophobic BC deposition (col) [kg m-1 s-1]
         flx_oc_dep       => aerosol_inst%flx_oc_dep_col       , & ! Output: [real(r8) (:)   ] total OC deposition (col) [kg m-2 s-1]  
         flx_oc_dep_wet   => aerosol_inst%flx_oc_dep_wet_col   , & ! Output: [real(r8) (:)   ] wet OC deposition (col) [kg m-2 s-1]    
         flx_oc_dep_dry   => aerosol_inst%flx_oc_dep_dry_col   , & ! Output: [real(r8) (:)   ] dry OC deposition (col) [kg m-2 s-1]    
         flx_oc_dep_phi   => aerosol_inst%flx_oc_dep_phi_col   , & ! Output: [real(r8) (:)   ] hydrophillic OC deposition (col) [kg m-1 s-1]
         flx_oc_dep_pho   => aerosol_inst%flx_oc_dep_pho_col   , & ! Output: [real(r8) (:)   ] hydrophobic OC deposition (col) [kg m-1 s-1]
         flx_dst_dep      => aerosol_inst%flx_dst_dep_col      , & ! Output: [real(r8) (:)   ] total dust deposition (col) [kg m-2 s-1]
         flx_dst_dep_wet1 => aerosol_inst%flx_dst_dep_wet1_col , & ! Output: [real(r8) (:)   ] wet dust (species 1) deposition (col) [kg m-2 s-1]
         flx_dst_dep_dry1 => aerosol_inst%flx_dst_dep_dry1_col , & ! Output: [real(r8) (:)   ] dry dust (species 1) deposition (col) [kg m-2 s-1]
         flx_dst_dep_wet2 => aerosol_inst%flx_dst_dep_wet2_col , & ! Output: [real(r8) (:)   ] wet dust (species 2) deposition (col) [kg m-2 s-1]
         flx_dst_dep_dry2 => aerosol_inst%flx_dst_dep_dry2_col , & ! Output: [real(r8) (:)   ] dry dust (species 2) deposition (col) [kg m-2 s-1]
         flx_dst_dep_wet3 => aerosol_inst%flx_dst_dep_wet3_col , & ! Output: [real(r8) (:)   ] wet dust (species 3) deposition (col) [kg m-2 s-1]
         flx_dst_dep_dry3 => aerosol_inst%flx_dst_dep_dry3_col , & ! Output: [real(r8) (:)   ] dry dust (species 3) deposition (col) [kg m-2 s-1]
         flx_dst_dep_wet4 => aerosol_inst%flx_dst_dep_wet4_col , & ! Output: [real(r8) (:)   ] wet dust (species 4) deposition (col) [kg m-2 s-1]
         flx_dst_dep_dry4 => aerosol_inst%flx_dst_dep_dry4_col   & ! Output: [real(r8) (:)   ] dry dust (species 4) deposition (col) [kg m-2 s-1]
         )

      !  set aerosol deposition fluxes from forcing array
      !  The forcing array is either set from an external file 
      !  or from fluxes received from the atmosphere model

      do c = bounds%begc,bounds%endc
         g = col%gridcell(c)

         flx_bc_dep_dry(c)   = forc_aer(g,1) + forc_aer(g,2)
         flx_bc_dep_wet(c)   = forc_aer(g,3)
         flx_bc_dep_phi(c)   = forc_aer(g,1) + forc_aer(g,3)
         flx_bc_dep_pho(c)   = forc_aer(g,2)
         flx_bc_dep(c)       = forc_aer(g,1) + forc_aer(g,2) + forc_aer(g,3)

         flx_oc_dep_dry(c)   = forc_aer(g,4) + forc_aer(g,5)
         flx_oc_dep_wet(c)   = forc_aer(g,6)
         flx_oc_dep_phi(c)   = forc_aer(g,4) + forc_aer(g,6)
         flx_oc_dep_pho(c)   = forc_aer(g,5)
         flx_oc_dep(c)       = forc_aer(g,4) + forc_aer(g,5) + forc_aer(g,6)

         flx_dst_dep_wet1(c) = forc_aer(g,7)
         flx_dst_dep_dry1(c) = forc_aer(g,8)
         flx_dst_dep_wet2(c) = forc_aer(g,9)
         flx_dst_dep_dry2(c) = forc_aer(g,10)
         flx_dst_dep_wet3(c) = forc_aer(g,11)
         flx_dst_dep_dry3(c) = forc_aer(g,12)
         flx_dst_dep_wet4(c) = forc_aer(g,13)
         flx_dst_dep_dry4(c) = forc_aer(g,14)
         flx_dst_dep(c)      = forc_aer(g,7) + forc_aer(g,8) + forc_aer(g,9) + &
                               forc_aer(g,10) + forc_aer(g,11) + forc_aer(g,12) + &
                               forc_aer(g,13) + forc_aer(g,14)
      end do

      ! aerosol deposition fluxes into top layer
      ! This is done after the inter-layer fluxes so that some aerosol
      ! is in the top layer after deposition, and is not immediately
      ! washed out before radiative calculations are done

      dtime = get_step_size()

      do fc = 1, num_snowc
         c = filter_snowc(fc)
         mss_bcphi(c,snl(c)+1) = mss_bcphi(c,snl(c)+1) + (flx_bc_dep_phi(c)*dtime)
         mss_bcpho(c,snl(c)+1) = mss_bcpho(c,snl(c)+1) + (flx_bc_dep_pho(c)*dtime)
         mss_ocphi(c,snl(c)+1) = mss_ocphi(c,snl(c)+1) + (flx_oc_dep_phi(c)*dtime)
         mss_ocpho(c,snl(c)+1) = mss_ocpho(c,snl(c)+1) + (flx_oc_dep_pho(c)*dtime)

         mss_dst1(c,snl(c)+1) = mss_dst1(c,snl(c)+1) + (flx_dst_dep_dry1(c) + flx_dst_dep_wet1(c))*dtime
         mss_dst2(c,snl(c)+1) = mss_dst2(c,snl(c)+1) + (flx_dst_dep_dry2(c) + flx_dst_dep_wet2(c))*dtime
         mss_dst3(c,snl(c)+1) = mss_dst3(c,snl(c)+1) + (flx_dst_dep_dry3(c) + flx_dst_dep_wet3(c))*dtime
         mss_dst4(c,snl(c)+1) = mss_dst4(c,snl(c)+1) + (flx_dst_dep_dry4(c) + flx_dst_dep_wet4(c))*dtime
      end do

    end associate

  end subroutine AerosolFluxes

end module AerosolMod
