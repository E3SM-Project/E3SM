module ecppvars
#ifdef ECPP

  use params, only: crm_rknd

  implicit none

  public

  integer, public, parameter :: nupdraft_in  = 1  ! Number of updraft class
  integer, public, parameter :: ndndraft_in  = 1  ! Number of dndraft class
  integer, public, parameter :: ncls_ecpp_in = 3  ! Number of total number of ecpp transport class
  ! = nupdraft_in+1+ndndraft_in
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
  real(crm_rknd),dimension(:,:,:,:)   , allocatable :: wwsqsum1
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
#endif /*ECPP*/
end module ecppvars
