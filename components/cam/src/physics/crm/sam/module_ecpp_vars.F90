module module_ecpp_vars
#ifdef ECPP

  use params, only: crm_rknd

  implicit none

  public

  integer, public, parameter :: nupdraft     = 1  ! Number of updraft class
  integer, public, parameter :: ndndraft     = 1  ! Number of dndraft class
  integer, public, parameter :: ncls_ecpp_in = 3  ! Number of total number of ecpp transport class
  integer, public, parameter :: ncc_in       = 2  ! number of clear/cloudy sub-calsses
  integer, public, parameter :: nprcp_in     = 2  ! Number of non-precipitating/precipitating sub-classes.

  integer, public, parameter :: QUI       = 1, &  !Quiescent class
  UP1       = 2     !First index for upward classes

  integer, public  :: DN1, & !First index of downward classes
  NCLASS_TR !Num. of transport classes
  !Both initialized based on
  !runtime settings

  integer, public :: NCLASS_CL = ncc_in, &  !Number of cloud classes
  CLR = 1, &        !Clear sub-class
  CLD = 2           !Cloudy sub-class

  integer, public :: NCLASS_PR = nprcp_in, &  !Number of precipitaion classes
  PRN = 1,       &  !Not precipitating sub-class
  PRY = 2           !Is precipitating sub-class

  real(crm_rknd),dimension(:,:,:,:)   , allocatable :: qlsink
  real(crm_rknd),dimension(:,:,:,:)   , allocatable :: precr
  real(crm_rknd),dimension(:,:,:,:)   , allocatable :: precsolid
  real(crm_rknd),dimension(:,:,:,:)   , allocatable :: rh
  real(crm_rknd),dimension(:,:,:,:)   , allocatable :: qlsink_bf
  real(crm_rknd),dimension(:,:,:,:)   , allocatable :: prain
  real(crm_rknd),dimension(:,:,:,:)   , allocatable :: qcloud_bf
  real(crm_rknd),dimension(:,:,:,:)   , allocatable :: qvs
  real(crm_rknd),dimension(:,:,:,:)   , allocatable :: qcloudsum1
  real(crm_rknd),dimension(:,:,:,:)   , allocatable :: qcloud_bfsum1
  real(crm_rknd),dimension(:,:,:,:)   , allocatable :: qrainsum1
  real(crm_rknd),dimension(:,:,:,:)   , allocatable :: qicesum1
  real(crm_rknd),dimension(:,:,:,:)   , allocatable :: qsnowsum1
  real(crm_rknd),dimension(:,:,:,:)   , allocatable :: qgraupsum1
  real(crm_rknd),dimension(:,:,:,:)   , allocatable :: qlsinksum1
  real(crm_rknd),dimension(:,:,:,:)   , allocatable :: qlsink_bfsum1
  real(crm_rknd),dimension(:,:,:,:)   , allocatable :: prainsum1
  real(crm_rknd),dimension(:,:,:,:)   , allocatable :: precrsum1
  real(crm_rknd),dimension(:,:,:,:)   , allocatable :: precsolidsum1
  real(crm_rknd),dimension(:,:,:,:)   , allocatable :: precallsum1
  real(crm_rknd),dimension(:,:,:,:)   , allocatable :: altsum1
  real(crm_rknd),dimension(:,:,:,:)   , allocatable :: rhsum1
  real(crm_rknd),dimension(:,:,:,:)   , allocatable :: cf3dsum1
  real(crm_rknd),dimension(:,:,:,:)   , allocatable :: wwsum1
  real(crm_rknd),dimension(:,:,:,:)   , allocatable :: tkesgssum1
  real(crm_rknd),dimension(:,:,:,:)   , allocatable :: qvssum1
  ! dim1 = z
  real(crm_rknd),dimension(:,:)       , allocatable :: xkhvsum
  real(crm_rknd),dimension(:,:)       , allocatable :: wup_thresh
  real(crm_rknd),dimension(:,:)       , allocatable :: wdown_thresh
  real(crm_rknd),dimension(:,:)       , allocatable :: wwqui_cen_sum
  real(crm_rknd),dimension(:,:)       , allocatable :: wwqui_bnd_sum
  real(crm_rknd),dimension(:,:)       , allocatable :: wwqui_cloudy_cen_sum
  real(crm_rknd),dimension(:,:)       , allocatable :: wwqui_cloudy_bnd_sum
  ! dims = (z, cloud sub-class, transport-class, precip sub-class)
  real(crm_rknd), dimension(:,:,:,:,:), allocatable :: area_bnd_final
  real(crm_rknd), dimension(:,:,:,:,:), allocatable :: area_bnd_sum
  real(crm_rknd), dimension(:,:,:,:,:), allocatable :: area_cen_final
  real(crm_rknd), dimension(:,:,:,:,:), allocatable :: area_cen_sum
  real(crm_rknd), dimension(:,:,:,:,:), allocatable :: mass_bnd_final
  real(crm_rknd), dimension(:,:,:,:,:), allocatable :: mass_bnd_sum
  real(crm_rknd), dimension(:,:,:,:,:), allocatable :: mass_cen_final
  real(crm_rknd), dimension(:,:,:,:,:), allocatable :: mass_cen_sum
  real(crm_rknd), dimension(:,:,:,:,:), allocatable :: ent_bnd_sum
  real(crm_rknd), dimension(:,:,:,:,:), allocatable :: rh_cen_sum
  real(crm_rknd), dimension(:,:,:,:,:), allocatable :: qcloud_cen_sum
  real(crm_rknd), dimension(:,:,:,:,:), allocatable :: qcloud_bf_cen_sum
  real(crm_rknd), dimension(:,:,:,:,:), allocatable :: qrain_cen_sum
  real(crm_rknd), dimension(:,:,:,:,:), allocatable :: qice_cen_sum
  real(crm_rknd), dimension(:,:,:,:,:), allocatable :: qsnow_cen_sum
  real(crm_rknd), dimension(:,:,:,:,:), allocatable :: qgraup_cen_sum
  real(crm_rknd), dimension(:,:,:,:,:), allocatable :: qlsink_cen_sum
  real(crm_rknd), dimension(:,:,:,:,:), allocatable :: precr_cen_sum
  real(crm_rknd), dimension(:,:,:,:,:), allocatable :: precsolid_cen_sum
  real(crm_rknd), dimension(:,:,:,:,:), allocatable :: precall_cen_sum
  real(crm_rknd), dimension(:,:,:,:,:), allocatable :: qlsink_bf_cen_sum
  real(crm_rknd), dimension(:,:,:,:,:), allocatable :: prain_cen_sum
  real(crm_rknd), dimension(:,:,:,:,:), allocatable :: qlsink_avg_cen_sum

contains

  !-----------------------------------------------------------------------------
  subroutine allocate_ecpp_vars(ncrms, ndraft)
    ! use openacc_utils
    implicit none
    integer, intent(in) :: ncrms

    allocate( qlsink(   nx, ny, nzm, ncrms) )
    allocate( precr(    nx, ny, nzm, ncrms) )
    allocate( precsolid(nx, ny, nzm, ncrms) )
    allocate( rh(       nx, ny, nzm, ncrms) )
    allocate( qvs(      nx, ny, nzm, ncrms) )

    allocate( qlsink_bf(nx, ny, nzm, ncrms) )
    allocate( prain(    nx, ny, nzm, ncrms) )
    allocate( qcloud_bf(nx, ny, nzm, ncrms) )

    allocate( qcloudsum1(   nx, ny, nzm,    ncrms) )
    allocate( qcloud_bfsum1(nx, ny, nzm,    ncrms) )
    allocate( qrainsum1(    nx, ny, nzm,    ncrms) )
    allocate( qicesum1(     nx, ny, nzm,    ncrms) )
    allocate( qsnowsum1(    nx, ny, nzm,    ncrms) )
    allocate( qgraupsum1(   nx, ny, nzm,    ncrms) )
    allocate( qlsinksum1(   nx, ny, nzm,    ncrms) ) 
    allocate( precrsum1(    nx, ny, nzm,    ncrms) )
    allocate( precsolidsum1(nx, ny, nzm,    ncrms) )
    allocate( precallsum1(  nx, ny, nzm,    ncrms) )
    allocate( altsum1(      nx, ny, nzm,    ncrms) )
    allocate( rhsum1(       nx, ny, nzm,    ncrms) )
    allocate( cf3dsum1(     nx, ny, nzm,    ncrms) )
    allocate( wwsum1(       nx, ny, nzstag, ncrms) )
    allocate( tkesgssum1(   nx, ny, nzm,    ncrms) )
    allocate( qlsink_bfsum1(nx, ny, nzm,    ncrms) )
    allocate( prainsum1(    nx, ny, nzm,    ncrms) )
    allocate( qvssum1(      nx, ny, nzm,    ncrms) )

    allocate( xkhvsum(nzm,ncrms) )

    allocate( wwqui_cen_sum(       nzm,  ncrms) )
    allocate( wwqui_bnd_sum(       nzm+1,ncrms) )
    allocate( wwqui_cloudy_cen_sum(nzm,  ncrms) )
    allocate( wwqui_cloudy_bnd_sum(nzm+1,ncrms) )

    allocate( wup_thresh(  nzm+1,ncrms) )
    allocate( wdown_thresh(nzm+1,ncrms) )

    allocate( area_bnd_final(    nzstag,NCLASS_CL, ndraft, NCLASS_PR, ncrms) )
    allocate( area_bnd_sum(      nzstag,NCLASS_CL, ndraft, NCLASS_PR, ncrms) )
    allocate( area_cen_final(    nzm   ,NCLASS_CL, ndraft, NCLASS_PR, ncrms) )
    allocate( area_cen_sum(      nzm   ,NCLASS_CL, ndraft, NCLASS_PR, ncrms) )
    allocate( mass_bnd_final(    nzstag,NCLASS_CL, ndraft, NCLASS_PR, ncrms) )
    allocate( mass_bnd_sum(      nzstag,NCLASS_CL, ndraft, NCLASS_PR, ncrms) )
    allocate( mass_cen_final(    nzm   ,NCLASS_CL, ndraft, NCLASS_PR, ncrms) )
    allocate( mass_cen_sum(      nzm   ,NCLASS_CL, ndraft, NCLASS_PR, ncrms) )
    allocate( ent_bnd_sum(       nzstag,NCLASS_CL, ndraft, NCLASS_PR, ncrms) )
    allocate( rh_cen_sum(        nzm   ,NCLASS_CL, ndraft, NCLASS_PR, ncrms) )
    allocate( qcloud_cen_sum(    nzm   ,NCLASS_CL, ndraft, NCLASS_PR, ncrms) )
    allocate( qcloud_bf_cen_sum( nzm   ,NCLASS_CL, ndraft, NCLASS_PR, ncrms) )
    allocate( qrain_cen_sum(     nzm   ,NCLASS_CL, ndraft, NCLASS_PR, ncrms) )
    allocate( qice_cen_sum(      nzm   ,NCLASS_CL, ndraft, NCLASS_PR, ncrms) )
    allocate( qsnow_cen_sum(     nzm   ,NCLASS_CL, ndraft, NCLASS_PR, ncrms) )
    allocate( qgraup_cen_sum(    nzm   ,NCLASS_CL, ndraft, NCLASS_PR, ncrms) )
    allocate( qlsink_cen_sum(    nzm   ,NCLASS_CL, ndraft, NCLASS_PR, ncrms) )
    allocate( precr_cen_sum(     nzm   ,NCLASS_CL, ndraft, NCLASS_PR, ncrms) )
    allocate( precsolid_cen_sum( nzm   ,NCLASS_CL, ndraft, NCLASS_PR, ncrms) )
    allocate( precall_cen_sum(   nzm   ,NCLASS_CL, ndraft, NCLASS_PR, ncrms) )
    allocate( qlsink_bf_cen_sum( nzm   ,NCLASS_CL, ndraft, NCLASS_PR, ncrms) )
    allocate( qlsink_avg_cen_sum(nzm   ,NCLASS_CL, ndraft, NCLASS_PR, ncrms) )
    allocate( prain_cen_sum(     nzm   ,NCLASS_CL, ndraft, NCLASS_PR, ncrms) )

  end subroutine allocate_ecpp_vars
  !-----------------------------------------------------------------------------
  subroutine deallocate_ecpp_vars()
    implicit none

    deallocate( qlsink )
    deallocate( precr )
    deallocate( precsolid )
    deallocate( rh )
    deallocate( qvs )
    deallocate( qlsink_bf )
    deallocate( prain )
    deallocate( cloud_bf )

    deallocate( qcloudsum1 )
    deallocate( qcloud_bfsum1 )
    deallocate( qrainsum1 )
    deallocate( qicesum1 )
    deallocate( qsnowsum1 )
    deallocate( qgraupsum1 )
    deallocate( qlsinksum1 )
    deallocate( precrsum1 )
    deallocate( precsolidsum1 )
    deallocate( precallsum1 )
    deallocate( altsum1 )
    deallocate( rhsum1 )
    deallocate( cf3dsum1 )
    deallocate( wwsum1 )
    deallocate( tkesgssum1 )
    deallocate( qlsink_bfsum1 )
    deallocate( prainsum1 )
    deallocate( qvssum1 )

    deallocate( xkhvsum )
    deallocate( wup_thresh )
    deallocate( wdown_thresh )

    deallocate( wwqui_cen_sum )
    deallocate( wwqui_bnd_sum )
    deallocate( wwqui_cloudy_cen_sum )
    deallocate( wwqui_cloudy_bnd_sum )

    deallocate( area_bnd_final )
    deallocate( area_bnd_sum )
    deallocate( area_cen_final )
    deallocate( area_cen_sum )
    deallocate( mass_bnd_final )
    deallocate( mass_bnd_sum )
    deallocate( mass_cen_final )
    deallocate( mass_cen_sum )
    deallocate( ent_bnd_sum )
    deallocate( rh_cen_sum )
    deallocate( qcloud_cen_sum )
    deallocate( qcloud_bf_cen_sum )
    deallocate( qrain_cen_sum )
    deallocate( qice_cen_sum )
    deallocate( qsnow_cen_sum )
    deallocate( qgraup_cen_sum )
    deallocate( qlsink_cen_sum )
    deallocate( qlsink_bf_cen_sum )
    deallocate( precr_cen_sum )
    deallocate( precsolid_cen_sum )
    deallocate( precall_cen_sum )
    deallocate( qlsink_avg_cen_sum )
    deallocate( prain_cen_sum  )

  end subroutine deallocate_ecpp_vars
  !-----------------------------------------------------------------------------
#endif /*ECPP*/
end module module_ecpp_vars
