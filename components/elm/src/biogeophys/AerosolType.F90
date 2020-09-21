module AerosolType

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use elm_varpar     , only : nlevsno, nlevgrnd
  use decompMod      , only : bounds_type
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  ! !PUBLIC DATA:
  !
  type, public :: aerosol_type
     real(r8), pointer :: mss_bcpho_col(:,:)      ! col mass of hydrophobic BC in snow (col,lyr)     [kg]
     real(r8), pointer :: mss_bcphi_col(:,:)      ! col mass of hydrophillic BC in snow (col,lyr)    [kg]
     real(r8), pointer :: mss_bctot_col(:,:)      ! col total mass of BC in snow (pho+phi) (col,lyr) [kg]
     real(r8), pointer :: mss_bc_col_col(:)       ! col column-integrated mass of total BC           [kg]
     real(r8), pointer :: mss_bc_top_col(:)       ! col top-layer mass of total BC                   [kg]

     real(r8), pointer :: mss_ocpho_col(:,:)      ! col mass of hydrophobic OC in snow (col,lyr)     [kg]
     real(r8), pointer :: mss_ocphi_col(:,:)      ! col mass of hydrophillic OC in snow (col,lyr)    [kg]
     real(r8), pointer :: mss_octot_col(:,:)      ! col total mass of OC in snow (pho+phi) (col,lyr) [kg]
     real(r8), pointer :: mss_oc_col_col(:)       ! col column-integrated mass of total OC           [kg]
     real(r8), pointer :: mss_oc_top_col(:)       ! col top-layer mass of total OC                   [kg]

     real(r8), pointer :: mss_dst1_col(:,:)       ! col mass of dust species 1 in snow (col,lyr)     [kg]
     real(r8), pointer :: mss_dst2_col(:,:)       ! col mass of dust species 2 in snow (col,lyr)     [kg]
     real(r8), pointer :: mss_dst3_col(:,:)       ! col mass of dust species 3 in snow (col,lyr)     [kg]
     real(r8), pointer :: mss_dst4_col(:,:)       ! col mass of dust species 4 in snow (col,lyr)     [kg]
     real(r8), pointer :: mss_dsttot_col(:,:)     ! col total mass of dust in snow (col,lyr)         [kg]
     real(r8), pointer :: mss_dst_col_col(:)      ! col column-integrated mass of dust in snow       [kg]
     real(r8), pointer :: mss_dst_top_col(:)      ! col top-layer mass of dust in snow               [kg]

     real(r8), pointer :: mss_cnc_bcphi_col(:,:)  ! col mass concentration of hydrophilic BC in snow (col,lyr) [kg/kg]
     real(r8), pointer :: mss_cnc_bcpho_col(:,:)  ! col mass concentration of hydrophilic BC in snow (col,lyr) [kg/kg]
     real(r8), pointer :: mss_cnc_ocphi_col(:,:)  ! col mass concentration of hydrophilic OC in snow (col,lyr) [kg/kg]
     real(r8), pointer :: mss_cnc_ocpho_col(:,:)  ! col mass concentration of hydrophilic OC in snow (col,lyr) [kg/kg]
     real(r8), pointer :: mss_cnc_dst1_col(:,:)   ! col mass concentration of dust species 1 in snow (col,lyr) [kg/kg]
     real(r8), pointer :: mss_cnc_dst2_col(:,:)   ! col mass concentration of dust species 2 in snow (col,lyr) [kg/kg]
     real(r8), pointer :: mss_cnc_dst3_col(:,:)   ! col mass concentration of dust species 3 in snow (col,lyr) [kg/kg]
     real(r8), pointer :: mss_cnc_dst4_col(:,:)   ! col mass concentration of dust species 4 in snow (col,lyr) [kg/kg]

     real(r8), pointer :: flx_dst_dep_dry1_col(:) ! col dust species 1 dry   deposition on ground (positive definite) [kg/s]
     real(r8), pointer :: flx_dst_dep_wet1_col(:) ! col dust species 1 wet   deposition on ground (positive definite) [kg/s]
     real(r8), pointer :: flx_dst_dep_dry2_col(:) ! col dust species 2 dry   deposition on ground (positive definite) [kg/s]
     real(r8), pointer :: flx_dst_dep_wet2_col(:) ! col dust species 2 wet   deposition on ground (positive definite) [kg/s]
     real(r8), pointer :: flx_dst_dep_dry3_col(:) ! col dust species 3 dry   deposition on ground (positive definite) [kg/s]
     real(r8), pointer :: flx_dst_dep_wet3_col(:) ! col dust species 3 wet   deposition on ground (positive definite) [kg/s]
     real(r8), pointer :: flx_dst_dep_dry4_col(:) ! col dust species 4 dry   deposition on ground (positive definite) [kg/s]
     real(r8), pointer :: flx_dst_dep_wet4_col(:) ! col dust species 4 wet   deposition on ground (positive definite) [kg/s]
     real(r8), pointer :: flx_dst_dep_col(:)      ! col total (dry+wet) dust deposition on ground (positive definite) [kg/s]

     real(r8), pointer :: flx_bc_dep_dry_col(:)   ! col dry (BCPHO+BCPHI) BC deposition on ground (positive definite) [kg/s]
     real(r8), pointer :: flx_bc_dep_wet_col(:)   ! col wet (BCPHI) BC       deposition on ground (positive definite) [kg/s]
     real(r8), pointer :: flx_bc_dep_pho_col(:)   ! col hydrophobic BC       deposition on ground (positive definite) [kg/s]
     real(r8), pointer :: flx_bc_dep_phi_col(:)   ! col hydrophillic BC      deposition on ground (positive definite) [kg/s]
     real(r8), pointer :: flx_bc_dep_col(:)       ! col total (dry+wet) BC   deposition on ground (positive definite) [kg/s]

     real(r8), pointer :: flx_oc_dep_dry_col(:)   ! col dry (OCPHO+OCPHI) OC deposition on ground (positive definite) [kg/s]
     real(r8), pointer :: flx_oc_dep_wet_col(:)   ! col wet (OCPHI) OC       deposition on ground (positive definite) [kg/s]
     real(r8), pointer :: flx_oc_dep_pho_col(:)   ! col hydrophobic OC       deposition on ground (positive definite) [kg/s]
     real(r8), pointer :: flx_oc_dep_phi_col(:)   ! col hydrophillic OC      deposition on ground (positive definite) [kg/s]
     real(r8), pointer :: flx_oc_dep_col(:)       ! col total (dry+wet) OC   deposition on ground (positive definite) [kg/s]

   contains

     procedure, public  :: Init         
     procedure, public  :: Restart      
     procedure, public  :: Reset 
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
    use elm_varcon    , only: spval
    use elm_varpar    , only: nlevsno 
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
    use elm_varpar , only : nlevsno, nlevsoi
    use elm_varcon , only : spval
    use clm_varctl , only : iulog  
    use elm_varpar , only : nlevsno
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

end module AerosolType
