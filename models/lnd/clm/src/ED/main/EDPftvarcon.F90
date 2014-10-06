module EDPftvarcon

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing vegetation constants and method to
  ! read and initialize vegetation (PFT) constants.
  !
  ! !USES:
  use clm_varpar  , only : mxpft
  use shr_kind_mod, only : r8 => shr_kind_r8

  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private

  !ED specific variables. 
  type, public ::  EDPftvarcon_type
     real(r8) :: max_dbh            (0:mxpft) ! maximum dbh at which height growth ceases...
     real(r8) :: freezetol          (0:mxpft) ! minimum temperature tolerance...
     real(r8) :: wood_density       (0:mxpft) ! wood density  g cm^-3  ...
     real(r8) :: alpha_stem         (0:mxpft) ! live stem turnover rate. y-1
     real(r8) :: hgt_min            (0:mxpft) ! sapling height m
     real(r8) :: cushion            (0:mxpft) ! labile carbon storage target as multiple of leaf pool.
     real(r8) :: leaf_stor_priority (0:mxpft) ! leaf turnover vs labile carbon use prioritisation. (1 = lose  leaves, 0 = use store).
     real(r8) :: leafwatermax       (0:mxpft) ! degree to which respiration is limited by btran if btran = 0
     real(r8) :: rootresist         (0:mxpft)
     real(r8) :: soilbeta           (0:mxpft)
     real(r8) :: crown              (0:mxpft)
     real(r8) :: bark_scaler        (0:mxpft)
     real(r8) :: crown_kill         (0:mxpft)
     real(r8) :: initd              (0:mxpft)
     real(r8) :: sd_mort            (0:mxpft)
     real(r8) :: seed_rain          (0:mxpft)
     real(r8) :: BB_slope           (0:mxpft)
     real(r8) :: root_long          (0:mxpft) ! root longevity (yrs)
     real(r8) :: clone_alloc        (0:mxpft) ! fraction of carbon balance allocated to clonal reproduction.
     real(r8) :: seed_alloc         (0:mxpft) ! fraction of carbon balance allocated to seeds.
     real(r8) :: sapwood_ratio      (0:mxpft) ! amount of sapwood per unit leaf carbon and m of height. gC/gC/m
  end type EDPftvarcon_type

  type(EDPftvarcon_type), public :: EDPftvarcon_inst

  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: EDpftconrd ! Read and initialize vegetation (PFT) constants
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine EDpftconrd( ncid )
    !
    ! !DESCRIPTION:
    ! Read and initialize vegetation (PFT) constants
    !
    ! !USES:
    use ncdio_pio , only : file_desc_t, ncd_io
    use abortutils  , only : endrun
    !
    ! !ARGUMENTS:
    implicit none
    !
    type(file_desc_t), intent(inout) :: ncid   ! pio netCDF file id

    ! !LOCAL VARIABLES:

    logical :: readv            ! read variable in or not
    character(len=32) :: subname = 'EDpftconrd'              ! subroutine name

    call ncd_io('max_dbh',EDPftvarcon_inst%max_dbh, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )

    call ncd_io('freezetol',EDPftvarcon_inst%freezetol, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )

    call ncd_io('wood_density',EDPftvarcon_inst%wood_density, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )

    call ncd_io('alpha_stem',EDPftvarcon_inst%alpha_stem, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )

    call ncd_io('hgt_min',EDPftvarcon_inst%hgt_min, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )

    call ncd_io('cushion',EDPftvarcon_inst%cushion, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )

    call ncd_io('leaf_stor_priority',EDPftvarcon_inst%leaf_stor_priority, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )

    call ncd_io('leafwatermax',EDPftvarcon_inst%leafwatermax, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )

    call ncd_io('rootresist',EDPftvarcon_inst%rootresist,'read',     ncid, readvar=readv)
    if  ( .not. readv ) call endrun( trim(subname)//' ERROR: error in reading in pft data' )

    call ncd_io('soilbeta',EDPftvarcon_inst%soilbeta,'read',         ncid, readvar=readv)
    if   ( .not. readv) call endrun(trim(subname)// ' ERROR : error in reading in pft data')

    call ncd_io('crown',EDPftvarcon_inst%crown,'read',         ncid, readvar=readv)
    if   ( .not. readv) call endrun(trim(subname)// ' ERROR : error in reading in pft data')

    call ncd_io('bark_scaler',EDPftvarcon_inst%bark_scaler,'read',         ncid, readvar=readv)
    if   ( .not. readv) call endrun(trim(subname)// ' ERROR : error in reading in pft data')

    call ncd_io('crown_kill',EDPftvarcon_inst%crown_kill,'read',         ncid, readvar=readv)
    if   ( .not. readv) call endrun(trim(subname)// ' ERROR : error in reading in pft data')

    call ncd_io('initd',EDPftvarcon_inst%initd,'read',         ncid, readvar=readv)
    if   ( .not. readv) call endrun(trim(subname)// ' ERROR : error in reading in pft data')

    call ncd_io('sd_mort',EDPftvarcon_inst%sd_mort,'read',         ncid, readvar=readv)
    if   ( .not. readv) call endrun(trim(subname)// ' ERROR : error in reading in pft data')

    call ncd_io('seed_rain',EDPftvarcon_inst%seed_rain,'read',         ncid, readvar=readv)
    if   ( .not. readv) call endrun(trim(subname)// ' ERROR : error in reading in pft data')

    call ncd_io('BB_slope',EDPftvarcon_inst%BB_slope,'read',         ncid, readvar=readv)
    if   ( .not. readv) call endrun(trim(subname)// ' ERROR : error in reading in pft data')

    call ncd_io('root_long',EDPftvarcon_inst%root_long, 'read', ncid,  readvar=readv)
    if   ( .not. readv) call endrun(trim(subname)// ' ERROR : error in reading in pft data')

    call ncd_io('seed_alloc',EDPftvarcon_inst%seed_alloc, 'read', ncid,  readvar=readv)
    if   ( .not. readv) call endrun(trim(subname)// ' ERROR : error in reading in pft data')

    call ncd_io('clone_alloc',EDPftvarcon_inst%clone_alloc, 'read', ncid,  readvar=readv)
    if   ( .not. readv) call endrun(trim(subname)// ' ERROR : error in reading in pft data')

    call ncd_io('sapwood_ratio',EDPftvarcon_inst%sapwood_ratio, 'read', ncid,  readvar=readv)
    if   ( .not. readv) call endrun(trim(subname)// ' ERROR : error in reading in pft data')

  end subroutine EDpftconrd

end module EDPftvarcon

