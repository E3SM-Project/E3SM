module  module_ecpp_crm_driver
#ifdef ECPP
  !------------------------------------------------------------------------
  ! module to prepare CRM output for ECPP module in the MMF model.
  !
  ! This code was written originally by William Gustafson, and is adopted into
  ! the MMF model by Minghuai Wang (minghuai.wang@pnl.gov), November, 2009.
  !
  ! Differences between the methodology here and in Ferret:
  !  - When calculating wup_bar and wdown_bar, points with w==0 are ignored
  !    here and were included in wup in Ferret.
  !  - Clear fluxes are no longer chopped off at the cloud top.
  !  - When calculating the std. dev. in and below the cloud, the level
  !    just above the cloud top is now included so we include w out the
  !    cloud top.
  !  - When determining "cloudyother" in Ferret the cloud above the
  !    interface was used. Now, the average of the cloud above and below
  !    is used.
  !
  !----------------------------------------------------------------------------------------
  use shr_kind_mod, only: r8 => shr_kind_r8
  use params, only: crm_rknd
  use ecppvars
  use ecppvars,  only: QUI, UP1, DN1, NCLASS_TR, NCLASS_CL, CLR, CLD, NCLASS_PR, PRN, PRY
  use cam_abortutils,  only: endrun

  public ecpp_crm_stat
  public ecpp_crm_init
  public ecpp_crm_cleanup

  private
  save

  integer ::  nxstag, nystag, nzstag
  integer :: ntavg1_ss  ! # of seconds to average when computing categories
  integer :: ntavg2_ss  ! # of seconds to average between outputs.
  integer :: ntavg1, ntavg2  ! number of CRM steps in ntavg[12]_ss period
  integer :: itavg1, itavg2  ! level-1 and level-2 counters

  real(crm_rknd) :: cloudthresh, prcpthresh, downthresh, upthresh

  real(crm_rknd) :: cloudthresh_trans        !  the threshold total cloud water for updraft or downdraft
  real(crm_rknd) :: precthresh_trans         !  the threshold total rain, snow and graupel for clear, updraft or downdraft

  integer, dimension(:,:), allocatable :: updraftbase, updrafttop, dndrafttop, dndraftbase
  integer :: nupdraft, ndndraft, ndraft_max

contains

  subroutine ecpp_crm_init(ncrms,dt_gl)
    use grid, only: nx, ny, nzm
    use module_ecpp_stats, only: zero_out_sums1, zero_out_sums2
    use module_ecpp_ppdriver2, only: nupdraft_in, ndndraft_in, ncls_ecpp_in
    implicit none
    real(r8), intent(in) :: dt_gl  ! global model's time step
    integer , intent(in) :: ncrms
    integer :: kbase, ktop
    integer :: m
    integer :: nup, ndn, icrm
    character(len=100) :: msg

    nxstag = nx+1
    nystag = ny+1
    nzstag = nzm+1

    upthresh    = 1.    ! Multiples of std. dev. to classify as updraft
    downthresh  = 1.    ! Multiples of std. dev. to classify as downdraft
    cloudthresh = 1e-6  ! Cloud mixing ratio beyond which cell is "cloudy(liquid)" (kg/kg)
    prcpthresh  = 1e-6  ! Preciptation rate (precr) beyond which cell is raining (kg/m2/s)
                        ! this is used to classify precip vs. nonprecip for wet scavenging.

    ! high thresholds are used to classify transport classes (following Xu et al., 2002, Q.J.R.M.S.
    cloudthresh_trans = 1e-5  !Cloud mixing ratio beyond which cell is "cloudy" to classify transport classes (kg/kg)
    
    ! the maxium of cloudthres_trans and 0.01*qvs is used to classify transport class
    precthresh_trans  = 1e-4  !Preciptation mixing ratio beyond which cell is raining to classify transport classes (kg/kg)

    !---------------------------------------------------------------------------
    ! determine number of updrafts and downdrafts
    !---------------------------------------------------------------------------
    ! updraft kbase & ktop definition:
    !    ww(i,j,k ) >  wup_thresh for k=kbase+1:ktop
    !               <= wup_thresh at  k=kbase and k=ktop+1
    ! they identify the "T-points" which enclose the updraft "W-points"
    !    and are affected by the subgrid transport of this updraft
    !
    ! downdraft kbase & ktop definition:
    !    ww(i,j,k ) <  wdown_thresh for k=kbase+1:ktop
    !               >= wdown_thresh at  k=kbase and k=ktop+1
    ! they identify the "T-points" which enclose the downdraft "W-points"
    !    and are affected by the subgrid transport of this downdraft
    !
    ! for both updrafts and downdrafts,
    !    1 <= kbase < ktop < nzstag
  
    nupdraft = 1
    ndndraft = 1

    DN1 = nupdraft + 2 !Setup index of first downdraft class
    NCLASS_TR = nupdraft + ndndraft + 1

    ndraft_max = 1 + nupdraft + ndndraft

    if(NCLASS_TR.ne.ncls_ecpp_in) then
      call endrun('NCLASS_TR should be equal to ncls_ecpp_in')
    end if
    if((nupdraft.ne.nupdraft_in) .or. (ndndraft.ne.ndndraft_in)) then
      call endrun('nupdraft or ndndraft is not set correctly')
    end if

    allocate (updraftbase(nupdraft,ncrms), updrafttop( nupdraft,ncrms) )
    allocate (dndraftbase(ndndraft,ncrms), dndrafttop( ndndraft,ncrms) )

    do icrm = 1 , ncrms  
      updraftbase(1,icrm)=1
      updrafttop( 1,icrm)=nzm
      dndrafttop( 1,icrm)=nzm
      dndraftbase(1,icrm)=1
    enddo

    !---------------------------------------------------------------------------
    ! Allocate arrays
    !---------------------------------------------------------------------------
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
    allocate( wwsqsum1(     nx, ny, nzstag, ncrms) )
    allocate( tkesgssum1(   nx, ny, nzm,    ncrms) )
    allocate( qlsink_bfsum1(nx, ny, nzm,    ncrms) )
    allocate( prainsum1(    nx, ny, nzm,    ncrms) )
    allocate( qvssum1(      nx, ny, nzm,    ncrms) )

    allocate( xkhvsum(nzm,ncrms) )

    allocate( wwqui_cen_sum(       nzm,  ncrms) )
    allocate( wwqui_bnd_sum(       nzm+1,ncrms)
    allocate( wwqui_cloudy_cen_sum(nzm,  ncrms) )
    allocate( wwqui_cloudy_bnd_sum(nzm+1,ncrms) )

    allocate( wup_thresh(  nzm+1,ncrms) )
    allocate( wdown_thresh(nzm+1,ncrms) )

    allocate( area_bnd_final(    nzstag,NCLASS_CL, ndraft_max, NCLASS_PR, ncrms) )
    allocate( area_bnd_sum(      nzstag,NCLASS_CL, ndraft_max, NCLASS_PR, ncrms) )
    allocate( area_cen_final(    nzm   ,NCLASS_CL, ndraft_max, NCLASS_PR, ncrms) )
    allocate( area_cen_sum(      nzm   ,NCLASS_CL, ndraft_max, NCLASS_PR, ncrms) )
    allocate( mass_bnd_final(    nzstag,NCLASS_CL, ndraft_max, NCLASS_PR, ncrms) )
    allocate( mass_bnd_sum(      nzstag,NCLASS_CL, ndraft_max, NCLASS_PR, ncrms) )
    allocate( mass_cen_final(    nzm   ,NCLASS_CL, ndraft_max, NCLASS_PR, ncrms) )
    allocate( mass_cen_sum(      nzm   ,NCLASS_CL, ndraft_max, NCLASS_PR, ncrms) )
    allocate( ent_bnd_sum(       nzstag,NCLASS_CL, ndraft_max, NCLASS_PR, ncrms) )
    allocate( rh_cen_sum(        nzm   ,NCLASS_CL, ndraft_max, NCLASS_PR, ncrms) )
    allocate( qcloud_cen_sum(    nzm   ,NCLASS_CL, ndraft_max, NCLASS_PR, ncrms) )
    allocate( qcloud_bf_cen_sum( nzm   ,NCLASS_CL, ndraft_max, NCLASS_PR, ncrms) )
    allocate( qrain_cen_sum(     nzm   ,NCLASS_CL, ndraft_max, NCLASS_PR, ncrms) )
    allocate( qice_cen_sum(      nzm   ,NCLASS_CL, ndraft_max, NCLASS_PR, ncrms) )
    allocate( qsnow_cen_sum(     nzm   ,NCLASS_CL, ndraft_max, NCLASS_PR, ncrms) )
    allocate( qgraup_cen_sum(    nzm   ,NCLASS_CL, ndraft_max, NCLASS_PR, ncrms) )
    allocate( qlsink_cen_sum(    nzm   ,NCLASS_CL, ndraft_max, NCLASS_PR, ncrms) )
    allocate( precr_cen_sum(     nzm   ,NCLASS_CL, ndraft_max, NCLASS_PR, ncrms) )
    allocate( precsolid_cen_sum( nzm   ,NCLASS_CL, ndraft_max, NCLASS_PR, ncrms) )
    allocate( precall_cen_sum(   nzm   ,NCLASS_CL, ndraft_max, NCLASS_PR, ncrms) )
    allocate( qlsink_bf_cen_sum( nzm   ,NCLASS_CL, ndraft_max, NCLASS_PR, ncrms) )
    allocate( qlsink_avg_cen_sum(nzm   ,NCLASS_CL, ndraft_max, NCLASS_PR, ncrms) )
    allocate( prain_cen_sum(     nzm   ,NCLASS_CL, ndraft_max, NCLASS_PR, ncrms) )

    !---------------------------------------------------------------------------
    ! Initialize the running sums.
    !---------------------------------------------------------------------------
    do icrm = 1 , ncrms
      call zero_out_sums1( qcloudsum1(:,:,:,icrm), qcloud_bfsum1(:,:,:,icrm), qrainsum1(:,:,:,icrm), &
      qicesum1(:,:,:,icrm), qsnowsum1(:,:,:,icrm), qgraupsum1(:,:,:,icrm), &
      qlsinksum1(:,:,:,icrm), precrsum1(:,:,:,icrm), &
      precsolidsum1(:,:,:,icrm), precallsum1(:,:,:,icrm), &
      altsum1(:,:,:,icrm), rhsum1(:,:,:,icrm), cf3dsum1(:,:,:,icrm), &
      wwsum1(:,:,:,icrm), wwsqsum1(:,:,:,icrm), tkesgssum1(:,:,:,icrm), &
      qlsink_bfsum1(:,:,:,icrm), prainsum1(:,:,:,icrm), qvssum1(:,:,:,icrm) )
      ndn = ndndraft ; nup = nupdraft
      call zero_out_sums2( &
      xkhvsum(:,icrm), &
      wwqui_cen_sum(:,icrm), wwqui_bnd_sum(:,icrm), wwqui_cloudy_cen_sum(:,icrm), wwqui_cloudy_bnd_sum(:,icrm),  &
      area_bnd_final(:,:,1:1+nup+ndn,:,icrm), area_bnd_sum(:,:,1:1+nup+ndn,:,icrm), &
      area_cen_final(:,:,1:1+nup+ndn,:,icrm), area_cen_sum(:,:,1:1+nup+ndn,:,icrm), &
      mass_bnd_final(:,:,1:1+nup+ndn,:,icrm), mass_bnd_sum(:,:,1:1+nup+ndn,:,icrm), &
      mass_cen_final(:,:,1:1+nup+ndn,:,icrm), mass_cen_sum(:,:,1:1+nup+ndn,:,icrm), &
      ent_bnd_sum(:,:,1:1+nup+ndn,:,icrm), &
      rh_cen_sum(:,:,1:1+nup+ndn,:,icrm), &
      qcloud_cen_sum(:,:,1:1+nup+ndn,:,icrm), qcloud_bf_cen_sum(:,:,1:1+nup+ndn,:,icrm), qrain_cen_sum(:,:,1:1+nup+ndn,:,icrm), &
      qice_cen_sum(:,:,1:1+nup+ndn,:,icrm),  qsnow_cen_sum(:,:,1:1+nup+ndn,:,icrm), &
      qgraup_cen_sum(:,:,1:1+nup+ndn,:,icrm), &
      qlsink_cen_sum(:,:,1:1+nup+ndn,:,icrm), precr_cen_sum(:,:,1:1+nup+ndn,:,icrm), &
      precsolid_cen_sum(:,:,1:1+nup+ndn,:,icrm), precall_cen_sum(:,:,1:1+nup+ndn,:,icrm), &
      qlsink_bf_cen_sum(:,:,1:1+nup+ndn,:,icrm), qlsink_avg_cen_sum(:,:,1:1+nup+ndn,:,icrm), &
      prain_cen_sum(:,:,1:1+nup+ndn,:,icrm)  )

      wup_thresh(:,icrm) = 0.0
      wdown_thresh(:,icrm) = 0.0
    enddo
    itavg1 = 0
    itavg2 = 0

    !---------------------------------------------------------------------------
    ! set ntavg[12] and initialize itavg[12] counters
    !---------------------------------------------------------------------------
    call ecpp_set_ntavg(dt_gl)

  end subroutine ecpp_crm_init
  !---------------------------------------------------------------------------------------

  !=======================================================================================
  subroutine ecpp_crm_cleanup ()
    ! deallocate variables
    deallocate( updraftbase )
    deallocate( updrafttop )
    deallocate( dndraftbase )
    deallocate( dndrafttop )

    deallocate( qlsink )
    deallocate( precr )
    deallocate( precsolid )
    deallocate( rh )
    deallocate( qvs )
    deallocate( qlsink_bf )
    deallocate( prain, )
    deallocate( cloud_bf)

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
    deallocate( wwsqsum1 )
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

  end subroutine ecpp_crm_cleanup
  !---------------------------------------------------------------------------------------

  !========================================================================================
  subroutine ecpp_crm_stat(ncrms)
    use module_ecpp_stats
    use module_data_ecpp1, only: afrac_cut
    use grid,  only: nx, ny, nzm, pres
    use vars,  only: w, tabs, p, CF3D
    use sgs, only: tke, tk
    use microphysics, only: micro_field, iqv, iqci, iqr, iqs, iqg, cloudliq
    use module_mp_GRAUPEL, only: POLYSVP

    implicit none
    integer, intent(in) :: ncrms
    integer :: i, ierr, i_tidx, j, ncnt1, ncnt2, icrm
    integer :: nup, ndn
    integer :: kbase, ktop, m
    integer :: ii, jj, kk
    integer :: icl, icls, ipr
    real(crm_rknd), dimension(nx, ny, nzm   , ncrms) :: qcloud, qrain, qice, qsnow, qgraup, precall, alt, xkhv, tketmp
    real(crm_rknd), dimension(nx, ny, nzstag, ncrms) :: ww, wwsq
    real(crm_rknd) :: EVS

    !------------------------------------------------------------------------
    ! Main code section...
    !------------------------------------------------------------------------

    ndn = ndndraft ; nup = nupdraft
    itavg1 = itavg1 + 1
    itavg2 = itavg2 + 1
    ndn = ndndraft ; nup = nupdraft

    ! Get values from SAM cloud fields
    do icrm = 1 , ncrms
      qcloud(1:nx,1:ny,1:nzm,icrm) = cloudliq(icrm,1:nx,1:ny,1:nzm)
      qrain (1:nx,1:ny,1:nzm,icrm) = micro_field(icrm,1:nx,1:ny,1:nzm,iqr )
      qice  (1:nx,1:ny,1:nzm,icrm) = micro_field(icrm,1:nx,1:ny,1:nzm,iqci)
      qsnow (1:nx,1:ny,1:nzm,icrm) = micro_field(icrm,1:nx,1:ny,1:nzm,iqs )
      qgraup(1:nx,1:ny,1:nzm,icrm) = micro_field(icrm,1:nx,1:ny,1:nzm,iqg )

      precall(:,:,:,icrm)= precr(:,:,:,icrm) + precsolid(:,:,:,icrm)

      do ii=1, nx
        do jj=1, ny
          do kk=1, nzm
            EVS = POLYSVP(tabs(icrm,ii,jj,kk),0)   ! saturation water vapor pressure (PA)
            qvs(ii,jj,kk,icrm) = .622*EVS/(pres(icrm,kk)*100.-EVS)  ! pres(icrm,kk) with unit of hPa
            alt(ii,jj,kk,icrm) =  287.*tabs(icrm,ii,jj,kk)/(100.*pres(icrm,kk))
          end do
        end do
      end do

      ww(:,:,:,icrm)     = w(icrm,1:nx,1:ny,1:nzstag)
      wwsq(:,:,:,icrm)   = 0.  ! subgrid vertical velocity is not used in the current version of ECPP.
      xkhv(:,:,:,icrm)   = tk(icrm,1:nx,1:ny,1:nzm)  ! eddy viscosity m2/s
    enddo

    ! Increment the 3-D running sums for averaging period 1.
    do icrm = 1 , ncrms
      tketmp(:,:,:,icrm) = tke(icrm,1:nx,1:ny,:)
    enddo
    call rsums1( ncrms, &
                 qcloud,    qcloudsum1,    &
                 qcloud_bf, qcloud_bfsum1, &
                 qrain,     qrainsum1,     &
                 qice,      qicesum1,      &
                 qsnow,     qsnowsum1,     &
                 qgraup,    qgraupsum1,    &
                 qlsink,    qlsinksum1,    &
                 precr,     precrsum1,     &
                 precsolid, precsolidsum1, &
                 precall,   precallsum1,   &
                 alt,       altsum1,       &
                 rh,        rhsum1,        &
                 CF3D,      cf3dsum1,      &
                 ww,        wwsum1,        &
                 wwsq,      wwsqsum1,      &
                 tketmp,    tkesgssum1,    &
                 qlsink_bf, qlsink_bfsum1, &
                 prain,     prainsum1,     &
                 qvs,       qvssum1   )

    ! Increment the running sums for the level two variables that are not
    ! already incremented. Consolidate from 3-D to 1-D columns.
    call rsums2(ncrms, xkhv, xkhvsum )

    ! Check if we have reached the end of the level 1 time averaging period.
    if( mod(itavg1,ntavg1) == 0 ) then

      ! Turn the running sums into averages.
      if( itavg1 /= 0 ) then
        ncnt1 = ntavg1
      else
        ncnt1 = 1
      end if
      call rsums1ToAvg( ncrms, ncnt1, qcloudsum1, qcloud_bfsum1, qrainsum1, &
      qicesum1, qsnowsum1, qgraupsum1, qlsinksum1, precrsum1, precsolidsum1, precallsum1, &
      altsum1, rhsum1, cf3dsum1, wwsum1, wwsqsum1, tkesgssum1, qlsink_bfsum1, prainsum1, qvssum1  )

      ! Determine draft categories and get running sums of them.
      do icrm = 1 , ncrms
      call categorization_stats( .true., &
      nx, ny, nzm, nupdraft, ndndraft, ndraft_max, &
      upthresh, downthresh, &
      cloudthresh, prcpthresh, &
      cloudthresh_trans, precthresh_trans,  &
      qvssum1(:,:,:,icrm),          &
      qcloudsum1(:,:,:,icrm), qcloud_bfsum1(:,:,:,icrm), qrainsum1(:,:,:,icrm), &
      qicesum1(:,:,:,icrm), qsnowsum1(:,:,:,icrm), qgraupsum1(:,:,:,icrm), &
      qlsinksum1(:,:,:,icrm), precrsum1(:,:,:,icrm), &
      precsolidsum1(:,:,:,icrm), precallsum1(:,:,:,icrm), &
      altsum1(:,:,:,icrm), rhsum1(:,:,:,icrm), cf3dsum1(:,:,:,icrm), &
      wwsum1(:,:,:,icrm), wwsqsum1(:,:,:,icrm), tkesgssum1(:,:,:,icrm),  &
      qlsink_bfsum1(:,:,:,icrm), prainsum1(:,:,:,icrm), &
      area_bnd_final(:,:,1:1+ndn+nup,:,icrm), area_cen_final(:,:,1:1+ndn+nup,:,icrm), &
      area_bnd_sum(:,:,1:1+ndn+nup,:,icrm), area_cen_sum(:,:,1:1+ndn+nup,:,icrm), &
      ent_bnd_sum(:,:,1:1+ndn+nup,:,icrm), mass_bnd_sum(:,:,1:1+ndn+nup,:,icrm), &
      rh_cen_sum(:,:,1:1+ndn+nup,:,icrm), &
      qcloud_cen_sum(:,:,1:1+ndn+nup,:,icrm), qcloud_bf_cen_sum(:,:,1:1+ndn+nup,:,icrm), qrain_cen_sum(:,:,1:1+ndn+nup,:,icrm), &
      qice_cen_sum(:,:,1:1+ndn+nup,:,icrm), qsnow_cen_sum(:,:,1:1+ndn+nup,:,icrm), &
      qgraup_cen_sum(:,:,1:1+ndn+nup,:,icrm), &
      qlsink_cen_sum(:,:,1:1+ndn+nup,:,icrm), precr_cen_sum(:,:,1:1+ndn+nup,:,icrm), &
      precsolid_cen_sum(:,:,1:1+nup+ndn,:,icrm), precall_cen_sum(:,:,1:1+nup+ndn,:,icrm),  &
      qlsink_bf_cen_sum(:,:,1:1+nup+ndn,:,icrm), prain_cen_sum(:,:,1:1+nup+ndn,:,icrm),  &
      wwqui_cen_sum(:,icrm), wwqui_bnd_sum(:,icrm), wwqui_cloudy_cen_sum(:,icrm), wwqui_cloudy_bnd_sum(:,icrm), &
      wup_thresh(:,icrm), wdown_thresh(:,icrm) )
      enddo

      ! If we want final area categories based on the last avg1 period in each
      ! avg2 then we need to zero out the running sum just created for the areas
      ! if it is not the last block of time in ntavg2
      if( .not. mod(itavg2,ntavg2)==0 ) then
        do icrm = 1 , ncrms
        call zero_out_areas( &
        area_bnd_final(:,:,1:1+ndn+nup,:,icrm), &
        area_cen_final(:,:,1:1+ndn+nup,:,icrm) )
        enddo
      end if

      ! Done with time level one averages so zero them out for next period.
      do icrm = 1 , ncrms
      call zero_out_sums1( qcloudsum1(:,:,:,icrm), qcloud_bfsum1(:,:,:,icrm), qrainsum1(:,:,:,icrm), &
      qicesum1(:,:,:,icrm), qsnowsum1(:,:,:,icrm), qgraupsum1(:,:,:,icrm), &
      qlsinksum1(:,:,:,icrm), precrsum1(:,:,:,icrm), &
      precsolidsum1(:,:,:,icrm), precallsum1(:,:,:,icrm), &
      altsum1(:,:,:,icrm), rhsum1(:,:,:,icrm), cf3dsum1(:,:,:,icrm), &
      wwsum1(:,:,:,icrm), wwsqsum1(:,:,:,icrm), tkesgssum1(:,:,:,icrm), &
      qlsink_bfsum1(:,:,:,icrm), prainsum1(:,:,:,icrm), qvssum1(:,:,:,icrm) )
      enddo

    end if !End of time level one averaging period

    ! Check if we have reached the end of a level 2 averaging period.
    if( mod(itavg2,ntavg2) == 0 ) then

      ! Turn the running sums into averages. ncnt1 in this case is the number
      ! of calls to categorization_stats during the level 2 averaging period,
      ! which increment the bnd/cen arrays.
      if( itavg2 /= 0 ) then
        ncnt1 = ntavg2_ss/ntavg1_ss
        ncnt2 = ntavg2
      else
        ncnt1 = 1
        ncnt2 = 1
      end if

      do icrm = 1 , ncrms
      call rsums2ToAvg( nx, ny, ncnt1, ncnt2, &
      xkhvsum(:,icrm), &
      wwqui_cen_sum(:,icrm), wwqui_bnd_sum(:,icrm), wwqui_cloudy_cen_sum(:,icrm), wwqui_cloudy_bnd_sum(:,icrm),  &
      area_bnd_final(:,:,1:1+ndn+nup,:,icrm), area_bnd_sum(:,:,1:1+ndn+nup,:,icrm), &
      area_cen_final(:,:,1:1+ndn+nup,:,icrm), area_cen_sum(:,:,1:1+ndn+nup,:,icrm), &
      mass_bnd_final(:,:,1:1+ndn+nup,:,icrm), mass_bnd_sum(:,:,1:1+ndn+nup,:,icrm), &
      mass_cen_final(:,:,1:1+ndn+nup,:,icrm), mass_cen_sum(:,:,1:1+ndn+nup,:,icrm), &
      ent_bnd_sum(:,:,1:1+ndn+nup,:,icrm), &
      rh_cen_sum(:,:,1:1+ndn+nup,:,icrm), &
      qcloud_cen_sum(:,:,1:1+ndn+nup,:,icrm), qcloud_bf_cen_sum(:,:,1:1+ndn+nup,:,icrm), qrain_cen_sum(:,:,1:1+ndn+nup,:,icrm), &
      qice_cen_sum(:,:,1:1+ndn+nup,:,icrm), qsnow_cen_sum(:,:,1:1+ndn+nup,:,icrm), &
      qgraup_cen_sum(:,:,1:1+ndn+nup,:,icrm), &
      qlsink_cen_sum(:,:,1:1+ndn+nup,:,icrm), precr_cen_sum(:,:,1:1+ndn+nup,:,icrm), &
      precsolid_cen_sum(:,:,1:1+ndn+nup,:,icrm), precall_cen_sum(:,:,1:1+ndn+nup,:,icrm), &
      qlsink_bf_cen_sum(:,:,1:1+ndn+nup,:,icrm), prain_cen_sum(:,:,1:1+ndn+nup,:,icrm) )
      enddo

      ! get in-cloud value for rh, qcloud, qrain, qice, qsnow, qgraup,
      ! percr, precsolid, and precall. (qlsink is already in-cloud values)
      do icrm = 1 , ncrms
      do kk=1, nzm
        do icl=1, NCLASS_CL
          do icls=1, ncls_ecpp_in
            do ipr=1, NCLASS_PR
              if(area_cen_sum(kk, icl, icls, ipr,icrm).gt.afrac_cut) then
                rh_cen_sum(kk,icl,icls,ipr,icrm) = rh_cen_sum(kk,icl,icls,ipr,icrm)/area_cen_sum(kk,icl,icls,ipr,icrm)
                qcloud_cen_sum(kk,icl,icls,ipr,icrm) = qcloud_cen_sum(kk,icl,icls,ipr,icrm)/area_cen_sum(kk,icl,icls,ipr,icrm)
                qcloud_bf_cen_sum(kk,icl,icls,ipr,icrm) = qcloud_bf_cen_sum(kk,icl,icls,ipr,icrm)/area_cen_sum(kk,icl,icls,ipr,icrm)
                qrain_cen_sum(kk,icl,icls,ipr,icrm) = qrain_cen_sum(kk,icl,icls,ipr,icrm)/area_cen_sum(kk,icl,icls,ipr,icrm)
                qice_cen_sum(kk,icl,icls,ipr,icrm) = qice_cen_sum(kk,icl,icls,ipr,icrm)/area_cen_sum(kk,icl,icls,ipr,icrm)
                qsnow_cen_sum(kk,icl,icls,ipr,icrm) = qsnow_cen_sum(kk,icl,icls,ipr,icrm)/area_cen_sum(kk,icl,icls,ipr,icrm)
                qgraup_cen_sum(kk,icl,icls,ipr,icrm) = qgraup_cen_sum(kk,icl,icls,ipr,icrm)/area_cen_sum(kk,icl,icls,ipr,icrm)
                precr_cen_sum(kk,icl,icls,ipr,icrm) = precr_cen_sum(kk,icl,icls,ipr,icrm)/area_cen_sum(kk,icl,icls,ipr,icrm)
                precsolid_cen_sum(kk,icl,icls,ipr,icrm) = precsolid_cen_sum(kk,icl,icls,ipr,icrm)/area_cen_sum(kk,icl,icls,ipr,icrm)
                precall_cen_sum(kk,icl,icls,ipr,icrm) = precall_cen_sum(kk,icl,icls,ipr,icrm)/area_cen_sum(kk,icl,icls,ipr,icrm)
                prain_cen_sum(kk,icl,icls,ipr,icrm) = prain_cen_sum(kk,icl,icls,ipr,icrm)/area_cen_sum(kk,icl,icls,ipr,icrm)
                if(qcloud_bf_cen_sum(kk,icl,icls,ipr,icrm).gt.1.0e-10) then
                  qlsink_avg_cen_sum(kk,icl,icls,ipr,icrm) = min(1.0/ntavg2_ss, prain_cen_sum(kk,icl,icls,ipr,icrm)/qcloud_bf_cen_sum(kk,icl,icls,ipr,icrm))
                else
                  qlsink_avg_cen_sum(kk,icl,icls,ipr,icrm) = 0.0
                end if
                qlsink_bf_cen_sum(kk,icl,icls,ipr,icrm) = min(1.0/ntavg2_ss, qlsink_bf_cen_sum(kk,icl,icls,ipr,icrm))
                qlsink_cen_sum(kk,icl,icls,ipr,icrm) = min(1.0/ntavg2_ss, qlsink_cen_sum(kk,icl,icls,ipr,icrm))
              else
                rh_cen_sum(kk,icl,icls,ipr,icrm) = 0.0
                qcloud_cen_sum(kk,icl,icls,ipr,icrm) = 0.0
                qcloud_bf_cen_sum(kk,icl,icls,ipr,icrm) = 0.0
                qrain_cen_sum(kk,icl,icls,ipr,icrm) = 0.0
                qice_cen_sum(kk,icl,icls,ipr,icrm) = 0.0
                qsnow_cen_sum(kk,icl,icls,ipr,icrm) = 0.0
                qgraup_cen_sum(kk,icl,icls,ipr,icrm) = 0.0
                precr_cen_sum(kk,icl,icls,ipr,icrm) = 0.0
                precsolid_cen_sum(kk,icl,icls,ipr,icrm) = 0.0
                precall_cen_sum(kk,icl,icls,ipr,icrm) = 0.0
                qlsink_bf_cen_sum(kk,icl,icls,ipr,icrm) = 0.0
                prain_cen_sum(kk,icl,icls,ipr,icrm) = 0.0
                qlsink_avg_cen_sum(kk,icl,icls,ipr,icrm) = 0.0
                qlsink_bf_cen_sum(kk,icl,icls,ipr,icrm) = 0.0
                qlsink_cen_sum(kk,icl,icls,ipr,icrm) = 0.0
              end if
            end do
          end do
        end do
        !
        ! calculate vertical velocity variance for quiescent class
        if(sum(area_cen_sum(kk,1:NCLASS_CL, QUI, 1:NCLASS_PR,icrm)).gt.afrac_cut) then
          wwqui_cen_sum(kk,icrm) = wwqui_cen_sum(kk,icrm) / sum(area_cen_sum(kk,1:NCLASS_CL, QUI, 1:NCLASS_PR,icrm))
        else
          wwqui_cen_sum(kk,icrm) = 0.0
        end if
        if(sum(area_cen_sum(kk,CLD, QUI, 1:NCLASS_PR,icrm)).gt.afrac_cut) then
          wwqui_cloudy_cen_sum(kk,icrm) = wwqui_cloudy_cen_sum(kk,icrm) / sum(area_cen_sum(kk, CLD, QUI, 1:NCLASS_PR,icrm))
        else
          wwqui_cloudy_cen_sum(kk,icrm) = 0.0
        end if

      end do  ! kk
      !
      ! calcualte vertical velocity variance for quiescent class at layer boundary
      do kk=1, nzm+1
        if(sum(area_bnd_sum(kk,1:NCLASS_CL, QUI, 1:NCLASS_PR,icrm)).gt.afrac_cut) then
          wwqui_bnd_sum(kk,icrm) = wwqui_bnd_sum(kk,icrm) / sum(area_bnd_sum(kk,1:NCLASS_CL, QUI, 1:NCLASS_PR,icrm))
        else
          wwqui_bnd_sum(kk,icrm) = 0.0
        end if
        if(sum(area_bnd_sum(kk,CLD, QUI, 1:NCLASS_PR,icrm)).gt.afrac_cut) then
          wwqui_cloudy_bnd_sum(kk,icrm) = wwqui_cloudy_bnd_sum(kk,icrm) / sum(area_bnd_sum(kk, CLD, QUI, 1:NCLASS_PR,icrm))
        else
          wwqui_cloudy_bnd_sum(kk,icrm) = 0.0
        end if
      end do
      enddo

    end if !End of level two time averaging period

  end subroutine ecpp_crm_stat

  subroutine ecpp_set_ntavg(dt_gl)
  !----------------------------------------------------------------------------
  ! Sets ntavg1_ss and ntavg2_ss periods for "level-1" and "level-2" averages,
  ! respectively. Also sets ntavg1 and ntavg2 which are the number of CRM
  ! steps in each averaging period.
  !----------------------------------------------------------------------------
    use grid, only: dt    ! CRM timestep (calling frequency of ecpp_crm_stat)
    implicit none
    real(r8), intent(in) :: dt_gl  ! global model's time step

    ! set level-1 and level-2 averaging periods for ECPP
    ntavg1_ss = min(600._r8, dt_gl)  ! lesser of 10 minutes or the GCM timestep
    ntavg2_ss = dt_gl                ! level-2 averaging period is GCM timestep

    ! Current implementation of ECPP requires that ntavg2_ss is a multiple of
    ! ntavg1_ss. Adjust ntavg1_ss upward from 10 minutes necessary.
    ! Ex. 1: ntavg2_ss = 1200 => ntavg1_ss = 1200/ (1200/600) = 1200 / 2 = 600
    ! Ex. 2: ntavg2_ss = 900 => ntavg1_ss = 900 / (900/600) = 900 / 1 = 900
    ntavg1_ss = ntavg2_ss / (ntavg2_ss / ntavg1_ss)

    ! calculate number of steps assuming dt evenly divides ntavg[12]_ss
    ! (alaways the case under current usage)
    ntavg1 = ntavg1_ss / dt
    ntavg2 = ntavg2_ss / dt
  end subroutine ecpp_set_ntavg

#endif /*ECPP*/
end module  module_ecpp_crm_driver
