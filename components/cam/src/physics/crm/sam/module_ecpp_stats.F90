!------------------------------------------------------------------------
! F90 module to calculate cloud-model stats needed as innput into ECPP.
!
! Routines in this module:
!   boundary_inout
!   categorization_stats
!   cloud_prcp_check
!   determine_transport_thresh
!   rsums1
!   rsums1ToAvg
!   rsums2
!   rsums2ToAvg
!   setup_class_masks
!   xyrsumof2d
!   xyrsumof3d
!   zero_out_areas
!   zero_out_sums1
!   zero_out_sums2
!
! William.Gustafson@pnl.gov; 20-Jul-2006
! Last modified: 16-Apr-2009, William.Gustafson@pnl.gov
!------------------------------------------------------------------------
module module_ecpp_stats
#ifdef ECPP

  use ecppvars, only: QUI, UP1, DN1, NCLASS_TR, NCLASS_CL, CLR, CLD, NCLASS_PR, PRN, PRY
  !==Guangxing Lin
  !use abortutils, only: endrun
  use cam_abortutils, only: endrun
  use params, only: crm_rknd
  implicit none

contains

  !------------------------------------------------------------------------
  subroutine boundary_inout( &
    nx, ny, nz, &
    uu, vv, &
    u_insum, u_outsum, v_insum, v_outsum )
    ! Calculates the average in/out-flow velocities and increments the
    ! running sum of the results.
    ! William.Gustafson@pnl.gov; 25-Jul-2006
    !------------------------------------------------------------------------
    integer, intent(in) :: nx, ny, nz
    real(crm_rknd), dimension(:,:,:), intent(in) :: uu, vv
    real(crm_rknd), dimension(:), intent(inout) :: u_insum, u_outsum, v_insum, v_outsum

    integer :: i, j, k, nxstag, nystag
    real(crm_rknd) :: spd_in, spd_out

    nxstag = nx+1
    nystag = ny+1
    !
    ! Running sum of inflow/outflow horizontal velocities...
    !
    ! 02-nov-2006 r.easter
    !    calculate separate in/outflow along x and y boundaries
    !       because of possibility of fixed boundary conditions
    !       and non-square domains
    !    for u_in & u_out, we want the "lineal" average along
    !       the west and east boundaries, so divide by ny
    !    for v_in & v_out, we want the "lineal" average along
    !       the south and north boundaries, so divide by nx
    !    previous code version divided by "nin" and "nout"
    !       which is incorrect
    !
    do k=1,nz

      spd_in  = 0.;  spd_out = 0.
      do j=1,ny
        ! Western boundary
        if( uu(1,j,k) >= 0. ) then
          spd_in  = spd_in + uu(1,j,k)
        else
          spd_out = spd_out - uu(1,j,k)
        end if

        ! Eastern boundary
        if( uu(nxstag,j,k) <= 0. ) then
          spd_in  = spd_in - uu(nxstag,j,k)
        else
          spd_out = spd_out + uu(nxstag,j,k)
        end if
      end do !j=ny
      u_insum(k)  = u_insum(k)  + spd_in /real(ny,crm_rknd)
      u_outsum(k) = u_outsum(k) + spd_out/real(ny,crm_rknd)

      spd_in  = 0.;  spd_out = 0.
      do i=1,nx
        ! Southern boundary
        if( vv(i,1,k) >= 0. ) then
          spd_in  = spd_in + vv(i,1,k)
        else
          spd_out = spd_out - vv(i,1,k)
        end if

        ! Northern boundary
        if( vv(i,nystag,k) <= 0. ) then
          spd_in  = spd_in - vv(i,nystag,k)
        else
          spd_out = spd_out + vv(i,nystag,k)
        end if
      end do !i=nx
      v_insum(k)  = v_insum(k)  + spd_in /real(nx,crm_rknd)
      v_outsum(k) = v_outsum(k) + spd_out/real(nx,crm_rknd)

    end do !k=nz
  end subroutine boundary_inout

  !------------------------------------------------------------------------
  subroutine rsums1( ncrms, qcloud,    qcloudsum1,    &
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
    cf3d,      cf3dsum1,      &
    ww,        wwsum1,        &
    wwsq,      wwsqsum1,      &
    tkesgs,    tkesgssum1,    &
    qlsink_bf, qlsink_bfsum1, &
    prain,     prainsum1,     &
    qvs,       qvssum1       )

    ! Increments 3-D running sums for the variables averaged every
    ! ntavg1_mm minutes.
    ! William.Gustafson@pnl.gov; 20-Jul-2006
    ! Last modified: William.Gustafson@pnl.gof; 25-Nov-2008
    !------------------------------------------------------------------------
    integer, intent(in) :: ncrms
    real(crm_rknd), dimension(:,:,:,:), intent(in) :: &
    qcloud, qcloud_bf, qrain, qice, qsnow, qgraup, &
    qlsink, precr, precsolid, precall, &
    alt, rh, cf3d, ww, wwsq, tkesgs, qlsink_bf, prain, qvs
    real(crm_rknd), dimension(:,:,:,:), intent(inout) :: &
    qcloudsum1, qcloud_bfsum1, qrainsum1, &
    qicesum1, qsnowsum1, qgraupsum1, &
    qlsinksum1, precrsum1, precsolidsum1, precallsum1, &
    altsum1, rhsum1, cf3dsum1, wwsum1, wwsqsum1, tkesgssum1, &
    qlsink_bfsum1, prainsum1, qvssum1
    integer :: icrm

    do icrm = 1 , ncrms
      qcloudsum1   (:,:,:,icrm) = qcloudsum1   (:,:,:,icrm) + qcloud(:,:,:,icrm)
      qcloud_bfsum1(:,:,:,icrm) = qcloud_bfsum1(:,:,:,icrm) + qcloud_bf(:,:,:,icrm)
      qrainsum1    (:,:,:,icrm) = qrainsum1    (:,:,:,icrm) + qrain(:,:,:,icrm)
      qicesum1     (:,:,:,icrm) = qicesum1     (:,:,:,icrm) + qice(:,:,:,icrm)
      qsnowsum1    (:,:,:,icrm) = qsnowsum1    (:,:,:,icrm) + qsnow(:,:,:,icrm)
      qgraupsum1   (:,:,:,icrm) = qgraupsum1   (:,:,:,icrm) + qgraup(:,:,:,icrm)
      qlsinksum1   (:,:,:,icrm) = qlsinksum1   (:,:,:,icrm) + qlsink(:,:,:,icrm)*qcloud(:,:,:,icrm)  ! Note this is converted back in rsum2ToAvg
      precrsum1    (:,:,:,icrm) = precrsum1    (:,:,:,icrm) + precr(:,:,:,icrm)
      precsolidsum1(:,:,:,icrm) = precsolidsum1(:,:,:,icrm) + precsolid(:,:,:,icrm)
      precallsum1  (:,:,:,icrm) = precallsum1  (:,:,:,icrm) + precall(:,:,:,icrm)
      altsum1      (:,:,:,icrm) = altsum1      (:,:,:,icrm) + alt(:,:,:,icrm)
      rhsum1       (:,:,:,icrm) = rhsum1       (:,:,:,icrm) + rh(:,:,:,icrm)
      cf3dsum1     (:,:,:,icrm) = cf3dsum1     (:,:,:,icrm) + cf3d(icrm,:,:,:)
      wwsum1       (:,:,:,icrm) = wwsum1       (:,:,:,icrm) + ww(:,:,:,icrm)
      wwsqsum1     (:,:,:,icrm) = wwsqsum1     (:,:,:,icrm) + wwsq(:,:,:,icrm)
      tkesgssum1   (:,:,:,icrm) = tkesgssum1   (:,:,:,icrm) + tkesgs(:,:,:,icrm)
      qlsink_bfsum1(:,:,:,icrm) = qlsink_bfsum1(:,:,:,icrm) + qlsink_bf(:,:,:,icrm)*qcloud_bf(:,:,:,icrm)  ! Note this is converted back in rsum2ToAvg
      prainsum1    (:,:,:,icrm) = prainsum1    (:,:,:,icrm) + prain(:,:,:,icrm)
      qvssum1      (:,:,:,icrm) = qvssum1      (:,:,:,icrm) + qvs(:,:,:,icrm)
    enddo
  end subroutine rsums1


  !------------------------------------------------------------------------
  subroutine rsums1ToAvg( ncrms, nt, qcloudsum, qcloud_bfsum, qrainsum, &
    qicesum, qsnowsum, qgraupsum, &
    qlsinksum, precrsum, precsolidsum, precallsum, &
    altsum, rhsum, cf3dsum, wwsum, wwsqsum, tkesgssum, qlsink_bfsum, prainsum, qvssum )
    ! Turns the columns of running sums into averages for the level one time
    ! period.
    ! William.Gustafson@pnl.gov; 20-Jul-2006
    ! Last modified: William.Gustafson@pnl.gov; 25-Nov-2008
    !------------------------------------------------------------------------
    integer, intent(in) :: nt, ncrms
    real(crm_rknd), dimension(:,:,:,:), intent(inout) :: &
    qcloudsum, qcloud_bfsum, qrainsum, qicesum, qsnowsum, qgraupsum, &
    qlsinksum, precrsum, precsolidsum, precallsum, &
    altsum, rhsum, cf3dsum, wwsum, wwsqsum, tkesgssum, qlsink_bfsum, prainsum, qvssum
    real(crm_rknd) :: ncount
    integer :: icrm

    ncount = real(nt,crm_rknd)

    do icrm = 1 , ncrms
      qcloudsum   (:,:,:,icrm) = qcloudsum   (:,:,:,icrm)/ncount
      qcloud_bfsum(:,:,:,icrm) = qcloud_bfsum(:,:,:,icrm)/ncount
      qrainsum    (:,:,:,icrm) = qrainsum    (:,:,:,icrm)/ncount
      qicesum     (:,:,:,icrm) = qicesum     (:,:,:,icrm)/ncount
      qsnowsum    (:,:,:,icrm) = qsnowsum    (:,:,:,icrm)/ncount
      qgraupsum   (:,:,:,icrm) = qgraupsum   (:,:,:,icrm)/ncount
      qlsinksum   (:,:,:,icrm) = qlsinksum   (:,:,:,icrm)/ncount
      precrsum    (:,:,:,icrm) = precrsum    (:,:,:,icrm)/ncount
      precsolidsum(:,:,:,icrm) = precsolidsum(:,:,:,icrm)/ncount
      precallsum  (:,:,:,icrm) = precallsum  (:,:,:,icrm)/ncount
      altsum      (:,:,:,icrm) = altsum      (:,:,:,icrm)/ncount
      rhsum       (:,:,:,icrm) = rhsum       (:,:,:,icrm)/ncount
      cf3dsum     (:,:,:,icrm) = cf3dsum     (:,:,:,icrm)/ncount
      wwsum       (:,:,:,icrm) = wwsum       (:,:,:,icrm)/ncount
      wwsqsum     (:,:,:,icrm) = wwsqsum     (:,:,:,icrm)/ncount
      tkesgssum   (:,:,:,icrm) = tkesgssum   (:,:,:,icrm)/ncount
      qlsink_bfsum(:,:,:,icrm) = qlsink_bfsum(:,:,:,icrm)/ncount
      prainsum    (:,:,:,icrm) = prainsum    (:,:,:,icrm)/ncount
      qvssum      (:,:,:,icrm) = qvssum      (:,:,:,icrm)/ncount
    enddo
  end subroutine rsums1ToAvg

  !------------------------------------------------------------------------
  subroutine rsums2(ncrms, xkhv, xkhvsum )
    ! Increment the running sums for the level 2 time averaging period for
    ! variables that are not already incremented (i.e. not the area and mass
    ! flux categories and in/out-flow speed that are already done). The 3-D
    ! variables are collapsed to 1-D columns.
    ! William.Gustafson@pnl.gov; 20-Jul-2006
    ! Last modified: William.Gustafson@pnl.gov; 25-Nov-2008
    !------------------------------------------------------------------------
    integer, intent(in) :: ncrms
    real(crm_rknd), dimension(:,:,:,:), intent(in) :: xkhv
    real(crm_rknd), dimension(:,:), intent(inout) :: xkhvsum
    ! Running sums of the simple variables that will be averaged...
    call xyrsumof3d(ncrms,xkhv,xkhvsum)
  end subroutine rsums2


  !------------------------------------------------------------------------
  subroutine rsums2ToAvg( areaavgtype, nx, ny, nt1, nt2, &
    xkhvsum, &
    wwqui_cen_sum, wwqui_bnd_sum, wwqui_cloudy_cen_sum, wwqui_cloudy_bnd_sum,  &
    area_bnd_final, area_bnd_sum, &
    area_cen_final, area_cen_sum, &
    mass_bnd_final, mass_bnd_sum, &
    mass_cen_final, mass_cen_sum, &
    ent_bnd_sum, &
    rh_cen_sum, &
    qcloud_cen_sum, qcloud_bf_cen_sum, qrain_cen_sum, &
    qice_cen_sum, qsnow_cen_sum, qgraup_cen_sum, &
    qlsink_cen_sum, precr_cen_sum, &
    precsolid_cen_sum, precall_cen_sum, &
    qlsink_bf_cen_sum, prain_cen_sum )

    ! Turns the columns of level two time period running sums into averages.
    ! Note that variables that the statistics variables use a different
    ! number of times.
    !
    ! nt1 = time length of average for area and mass for areaavgtype=2
    ! nt2 = time length of average for 2nd averaging period (the whole time)
    !
    ! William.Gustafson@pnl.gov; 20-Jul-2006
    ! Last modified: 16-Apr-2009, wig
    !------------------------------------------------------------------------
    integer, intent(in) :: areaavgtype, nx, ny, nt1, nt2
    real(crm_rknd), dimension(:), intent(inout) :: &
    xkhvsum, wwqui_cen_sum, wwqui_bnd_sum, wwqui_cloudy_cen_sum, wwqui_cloudy_bnd_sum
    real(crm_rknd), dimension(:,:,:,:), intent(inout) :: &
    area_bnd_final, area_bnd_sum, &
    area_cen_final, area_cen_sum, &
    mass_bnd_final, mass_bnd_sum, &
    mass_cen_final, mass_cen_sum, &
    ent_bnd_sum, rh_cen_sum, &
    qcloud_cen_sum, qcloud_bf_cen_sum, qrain_cen_sum, &
    qice_cen_sum, qsnow_cen_sum, qgraup_cen_sum, &
    qlsink_cen_sum, precr_cen_sum, &
    precsolid_cen_sum, precall_cen_sum, &
    qlsink_bf_cen_sum, prain_cen_sum
    integer :: i, k
    real(crm_rknd) :: ncount2, ncountwind, thesum

    !  print*,"...end of level two averaging period."

    ncount2    = real(nx*ny*nt2,crm_rknd)
    ncountwind = real((nx+1)*ny*nt2,crm_rknd)

    xkhvsum      = xkhvsum/ncount2

    ! Only touch final areas if doing averages over ntavg2
    if( areaavgtype == 2 ) then
      area_bnd_final = area_bnd_final/real(nt1,crm_rknd)
      area_cen_final = area_cen_final/real(nt1,crm_rknd)
    end if

    area_bnd_sum       = area_bnd_sum   /real(nt1,crm_rknd)
    area_cen_sum       = area_cen_sum   /real(nt1,crm_rknd)
    ent_bnd_sum        = ent_bnd_sum    /real(nt1,crm_rknd)
    mass_bnd_sum       = mass_bnd_sum   /real(nt1,crm_rknd)
    mass_cen_sum       = mass_cen_sum   /real(nt1,crm_rknd)
    rh_cen_sum         = rh_cen_sum     /real(nt1,crm_rknd)
    qcloud_cen_sum     = qcloud_cen_sum /real(nt1,crm_rknd)
    qcloud_bf_cen_sum  = qcloud_bf_cen_sum/real(nt1,crm_rknd)
    qrain_cen_sum      = qrain_cen_sum  /real(nt1,crm_rknd)
    qice_cen_sum       = qice_cen_sum   /real(nt1,crm_rknd)
    qsnow_cen_sum      = qsnow_cen_sum  /real(nt1,crm_rknd)
    qgraup_cen_sum     = qgraup_cen_sum /real(nt1,crm_rknd)
    do k=1,size(qlsink_cen_sum,1) !Note: must be after qcloud_cen_sum is turned into an avg
      ! see rsums1 where qlsink=qlsink*qcloud
      thesum = sum(qcloud_cen_sum(k,:,:,:))
      if( thesum > 1e-25 ) then
        qlsink_cen_sum(k,:,:,:) = qlsink_cen_sum(k,:,:,:)/thesum/real(nt1,crm_rknd)
      else
        qlsink_cen_sum(k,:,:,:) = 0.
      end if
    end do
    precr_cen_sum      = precr_cen_sum/real(nt1,crm_rknd)
    precsolid_cen_sum  = precsolid_cen_sum/real(nt1,crm_rknd)
    precall_cen_sum    = precall_cen_sum/real(nt1,crm_rknd)
    do k=1,size(qlsink_bf_cen_sum,1) !Note: must be after qcloud_bf_cen_sum is turned into an avg
      ! see rsums1 where qlsink=qlsink*qcloud
      thesum = sum(qcloud_bf_cen_sum(k,:,:,:))
      if( thesum > 1e-25 ) then
        qlsink_bf_cen_sum(k,:,:,:) = qlsink_bf_cen_sum(k,:,:,:)/thesum/real(nt1,crm_rknd)
      else
        qlsink_bf_cen_sum(k,:,:,:) = 0.
      end if
    end do

    prain_cen_sum   = prain_cen_sum/real(nt1,crm_rknd)
    wwqui_cen_sum =  wwqui_cen_sum / real(nt1,crm_rknd)
    wwqui_bnd_sum =  wwqui_bnd_sum / real(nt1,crm_rknd)
    wwqui_cloudy_cen_sum = wwqui_cloudy_cen_sum / real(nt1,crm_rknd)
    wwqui_cloudy_bnd_sum = wwqui_cloudy_bnd_sum / real(nt1,crm_rknd)

  end subroutine rsums2ToAvg


  !------------------------------------------------------------------------
  subroutine xyrsumof2d(xin,sumout)
    ! For a 2-D intput variable (x,y), the x & y dimensions are summed and
    ! added to a running sum.
    ! William.Gustafson@pnl.gov; 25-Apr-2006
    !------------------------------------------------------------------------
    real(crm_rknd), dimension(:,:), intent(in) :: xin
    real(crm_rknd), intent(out) :: sumout

    sumout = sumout + sum(xin(:,:))
  end subroutine xyrsumof2d


  !------------------------------------------------------------------------
  subroutine xyrsumof3d(ncrms,xin,sumout)
    ! For a 3-D intput variable (x,y,z), the x & y dimensions are summed and
    ! added to a column  to return a running sum.
    ! William.Gustafson@pnl.gov; 26-Jun-2006
    !------------------------------------------------------------------------
    integer, intent(in) :: ncrms
    real(crm_rknd), dimension(:,:,:,:), intent(in) :: xin
    real(crm_rknd), dimension(:,:), intent(out) :: sumout
    integer :: k, icrm
    do icrm = 1 , ncrms
      do k=1,ubound(sumout,1)
        sumout(k,icrm) = sumout(k,icrm) + sum(xin(:,:,k,icrm))
      end do
    end do
  end subroutine xyrsumof3d


  !------------------------------------------------------------------------
  subroutine zero_out_areas( &
    area_bnd_final, area_cen_final )
    ! Zeros out the running sums of final area categories.
    ! William.Gustafson@pnl.gov; 19-Nov-2008
    !------------------------------------------------------------------------
    real(crm_rknd), dimension(:,:,:,:), intent(out) :: &
    area_bnd_final, area_cen_final

    area_bnd_final=0.
    area_cen_final=0.
  end subroutine zero_out_areas


  !------------------------------------------------------------------------
  subroutine zero_out_sums1( qcloudsum, qcloud_bfsum, qrainsum,                         &
    qicesum, qsnowsum, qgraupsum,                &
    qlsink, precr, precsolid, precall,           &
    altsum, rhsum, cf3dsum, wwsum, wwsqsum, tkesgssum,    &
    qlsink_bfsum, prainsum, qvssum     )
    ! Zeros out running sum arrays that are averaged every ntavg1_mm minutes.
    ! William.Gustafson@pnl.gov; 20-Jul-2006
    ! Last modified: William.Gustafson@pnl.gov; 25-Nov-2008
    !------------------------------------------------------------------------
    real(crm_rknd),dimension(:,:,:), intent(out) :: &
    qcloudsum, qcloud_bfsum, qrainsum, qicesum, qsnowsum, qgraupsum, &
    qlsink, precr, precsolid, precall, &
    altsum, rhsum, cf3dsum, wwsum, wwsqsum, tkesgssum, qlsink_bfsum, prainsum, qvssum

    qcloudsum=0.
    qcloud_bfsum=0.
    qrainsum=0.
    qicesum=0.
    qsnowsum=0.
    qgraupsum=0.
    qlsink=0.
    precr=0.
    precsolid=0.
    precall=0.
    altsum=0.
    rhsum=0.
    cf3dsum=0.
    wwsum=0.
    wwsqsum=0.
    tkesgssum=0.
    qlsink_bfsum=0.0
    prainsum=0.0
    qvssum=0.0
  end subroutine zero_out_sums1


  !------------------------------------------------------------------------
  subroutine zero_out_sums2( &
    xkhvsum, &
    wwqui_cen_sum, wwqui_bnd_sum, wwqui_cloudy_cen_sum, wwqui_cloudy_bnd_sum,  &
    area_bnd_final, area_bnd_sum, area_cen_final, area_cen_sum, &
    mass_bnd_final, mass_bnd_sum, mass_cen_final, mass_cen_sum, &
    ent_bnd_sum, &
    rh_cen_sum, &
    qcloud_cen_sum, qcloud_bf_cen_sum, qrain_cen_sum, &
    qice_cen_sum, qsnow_cen_sum, qgraup_cen_sum, &
    qlsink_cen_sum, &
    precr_cen_sum, precsolid_cen_sum, precall_cen_sum, &
    qlsink_bf_cen_sum, qlsink_avg_cen_sum, prain_cen_sum )
    ! Zeros out running sum arrays that are averaged every ntavg2_mm minutes.
    ! William.Gustafson@pnl.gov; 20-Jul-2006
    ! Last modified: 25-Nov-2008, wig
    !------------------------------------------------------------------------
    real(crm_rknd),dimension(:), intent(out) :: &
    xkhvsum, wwqui_cen_sum, wwqui_bnd_sum, wwqui_cloudy_cen_sum, wwqui_cloudy_bnd_sum
    real(crm_rknd),dimension(:,:,:,:), intent(out) :: &
    area_bnd_final, area_bnd_sum, area_cen_final, area_cen_sum, &
    mass_bnd_final, mass_bnd_sum, mass_cen_final, mass_cen_sum, &
    ent_bnd_sum, rh_cen_sum, &
    qcloud_cen_sum, qcloud_bf_cen_sum, qrain_cen_sum, &
    qice_cen_sum, qsnow_cen_sum, qgraup_cen_sum, &
    qlsink_cen_sum, &
    precr_cen_sum, precsolid_cen_sum, precall_cen_sum, &
    qlsink_bf_cen_sum, qlsink_avg_cen_sum, prain_cen_sum

    xkhvsum=0.
    wwqui_cen_sum=0.
    wwqui_bnd_sum=0.
    wwqui_cloudy_cen_sum=0.
    wwqui_cloudy_bnd_sum=0.
    area_bnd_final=0.
    area_bnd_sum=0.
    area_cen_final=0.
    area_cen_sum=0.
    mass_bnd_final=0.
    mass_bnd_sum=0.
    mass_cen_final=0.
    mass_cen_sum=0.
    ent_bnd_sum=0.
    rh_cen_sum=0.
    qcloud_cen_sum=0.
    qcloud_bf_cen_sum=0.
    qrain_cen_sum=0.
    qice_cen_sum=0.
    qsnow_cen_sum=0.
    qgraup_cen_sum=0.
    qlsink_cen_sum=0.
    precr_cen_sum=0.
    precsolid_cen_sum=0.
    precall_cen_sum=0.
    qlsink_bf_cen_sum=0.
    qlsink_avg_cen_sum=0.
    prain_cen_sum=0.
  end subroutine zero_out_sums2


  !------------------------------------------------------------------------
  subroutine categorization_stats( domass, &
    nx, ny, nz, nupdraft, ndndraft, ndraft_max, &
    mode_updnthresh, upthresh, downthresh, &
    upthresh2, downthresh2, cloudthresh, prcpthresh, &
    cloudthresh_trans, precthresh_trans,  &
    qvs,                    &
    plumetype, allcomb, &
    qcloud, qcloud_bf, qrain, qice, qsnow, qgraup, &
    qlsink, precr, precsolid, precall, &
    alt, rh, cf3d, ww, wwsq, tkesgs, &
    qlsink_bf, prain, &
    area_bnd_final, area_cen_final, &
    area_bnd_sum, area_cen_sum, ent_bnd_sum, mass_bnd_sum, &
    rh_cen_sum, &
    qcloud_cen_sum, qcloud_bf_cen_sum, qrain_cen_sum, &
    qice_cen_sum, qsnow_cen_sum, qgraup_cen_sum, &
    qlsink_cen_sum, precr_cen_sum, &
    precsolid_cen_sum, precall_cen_sum, &
    qlsink_bf_cen_sum, prain_cen_sum,  &
    wwqui_cen_sum, wwqui_bnd_sum, wwqui_cloudy_cen_sum, wwqui_cloudy_bnd_sum, &
    wup_thresh, wdown_thresh )
    !
    ! William.Gustafson@pnl.gov; 25-Nov-2008
    ! Last modified: William.Gustafson@pnl.gov; 16-Apr-2009
    !------------------------------------------------------------------------
    use module_data_ecpp1, only: a_quiescn_minaa
    !
    ! Subroutine arguments...
    !
    logical, intent(in) :: domass !calculate mass fluxes? T/F
    integer, intent(in) :: nx, ny, nz, nupdraft, ndndraft, ndraft_max, &
    mode_updnthresh, plumetype
    logical, intent(in) :: allcomb
    real(crm_rknd), intent(in) :: &
    cloudthresh, prcpthresh, &
    downthresh, upthresh, &
    downthresh2, upthresh2
    real(crm_rknd), intent(in) :: cloudthresh_trans,  precthresh_trans
    real(crm_rknd), dimension(:,:,:), intent(in) :: &
    qcloud, qcloud_bf, qrain, qice, qsnow, qgraup, &
    qlsink, precr, precsolid, precall, &
    alt, rh, cf3d, ww, wwsq, tkesgs, qlsink_bf, prain, qvs
    real(crm_rknd), dimension(:,:,:,:), intent(inout) :: &
    area_bnd_final, area_cen_final, &
    area_bnd_sum, area_cen_sum, ent_bnd_sum, mass_bnd_sum, &
    rh_cen_sum, &
    qcloud_cen_sum, qcloud_bf_cen_sum, qrain_cen_sum, &
    qice_cen_sum, qsnow_cen_sum, qgraup_cen_sum, &
    qlsink_cen_sum, precr_cen_sum, &
    precsolid_cen_sum, precall_cen_sum, qlsink_bf_cen_sum, prain_cen_sum

    real(crm_rknd), dimension(:), intent(inout) :: wwqui_cen_sum, wwqui_bnd_sum, wwqui_cloudy_cen_sum, wwqui_cloudy_bnd_sum
    real(crm_rknd), dimension(nz+1), intent(out)  :: wdown_thresh, wup_thresh
    !
    ! Local vars...
    !
    real(crm_rknd), dimension(nx,ny,nz+1,NCLASS_CL,ndraft_max,NCLASS_PR) :: mask_bnd
    real(crm_rknd), dimension(nx,ny,nz,NCLASS_CL,ndraft_max,NCLASS_PR) :: mask_cen
    real(crm_rknd), dimension(nz+1,2) :: wdown_thresh_k, wup_thresh_k
    real(crm_rknd), dimension(nx,ny,nz) :: cloudmixr, cloudmixr_total, precmixr_total
    integer, dimension(nx,ny)  :: cloudtop
    real(crm_rknd), dimension(nz+1)  :: wup_rms_k, wup_bar_k, wup_stddev_k  &
    , wdown_rms_k, wdown_bar_k, wdown_stddev_k
    integer  :: kup_top, kdown_top  ! defined as the maximum level that allows updraft and downdraft
    real(crm_rknd) :: mask, wwrho_k, wwrho_km1
    real(crm_rknd), dimension(nz+1) :: rhoair      ! layer-averaged air density
    real(crm_rknd) :: wlarge = 1.0e10   ! m/s
    real(crm_rknd) :: tmpa, tmpb
    real(crm_rknd), dimension(nz) :: thresh_factorbb_up, thresh_factorbb_down
    real(crm_rknd) :: acen_quiesc, acen_up, acen_down, abnd_quiesc, abnd_up, abnd_down
    real(crm_rknd) :: acen_quiesc_minaa
    real(crm_rknd) :: wwqui_bar_cen(nz), wwqui_bar_bnd(nz+1), wwqui_cloudy_bar_cen(nz), wwqui_cloudy_bar_bnd(nz+1)

    integer :: i, icl, ipr, itr, j, k, km0, km1, km2, nxy, nzstag
    integer :: iter

    logical :: thresh_calc_not_done

    acen_quiesc_minaa = a_quiescn_minaa + 0.01

    nxy = nx*ny
    nzstag = nz+1

    ! Transport classification is based on total condensate (cloudmixr_total), and
    ! cloudy (liquid) and clear (non-liquid) classification is based on liquid water,
    ! because wet deposition, aqueous chemistry, and droplet activaton, all are for liquid clouds.
    !
    ! Minghuai Wang, 2010-04
    !
    cloudmixr = qcloud
    cloudmixr_total = qcloud + qice

    ! total hydrometer (rain, snow, and graupel)
    precmixr_total = qrain+qsnow+qgraup

    rhoair(:) = 0.0
    do j=1,ny
      do i=1,nx
        !
        ! Get cloud top height
        ! Cloud top height is used to determine whether there is updraft/downdraft. No updraft and
        ! downdraft is allowed above the condensate level (both liquid and ice).
        cloudtop(i,j) = 1 !Default to bottom level if no cloud in column.
        do k=nz,1,-1
          if( cloudmixr_total(i,j,k) >= cloudthresh_trans ) then
            !
            ! 0.01*qvs may be too large at low level.
            !           if( cloudmixr_total(i,j,k) >= max(0.01*qvs(i,j,k), cloudthresh_trans) ) then
            cloudtop(i,j) = k
            exit
          end if
        end do !k
        !
        ! Get layer-averaged air density
        do k=1, nzstag
          km0 = min(nz,k)
          km1 = max(1,k-1)
          rhoair(k) = rhoair(k)+0.5*(1.0/alt(i,j,km1) + 1.0/alt(i,j,km0))/real(nxy,crm_rknd)
        end do
      end do !i
    end do !j

    call determine_transport_thresh( &
    nx, ny, nz, &
    mode_updnthresh, upthresh, downthresh, &
    upthresh2, downthresh2, cloudthresh, &
    ww, rhoair, &
    wdown_thresh_k, wup_thresh_k    &
    , cloudtop                       &
    , wup_rms_k, wup_bar_k, wup_stddev_k  &
    , wdown_rms_k, wdown_bar_k, wdown_stddev_k   &
    , kup_top, kdown_top )

    wdown_thresh(:) = wdown_thresh_k(:,1)
    wup_thresh(:) = wup_thresh_k(:,1)

    if ((nupdraft > 1) .or. (ndndraft > 1)) then
      call endrun('*** code for thresh_factorbb_up/down needs nup/dndraft = 1')
    end if
    thresh_factorbb_up(:) = 1.0 ; thresh_factorbb_down(:) = 1.0
    thresh_calc_not_done = .true.

    iter = 0
    thresh_calc_loop: &
    do while ( thresh_calc_not_done )

      iter = iter + 1
      ! if quiescent class area was too small on previous iteration,
      !     then thresh_factor_acen_quiesc will be > 1.0
      ! multiply wup/down_thresh_k by this factor to reduce the
      !     up/downdraft areas and increase the quiescent area
      do k = 1, nzstag
        if (k == 1) then
          tmpa = thresh_factorbb_up(k)
          tmpb = thresh_factorbb_down(k)
        else if (k == nzstag) then
          tmpa = thresh_factorbb_up(k-1)
          tmpb = thresh_factorbb_down(k-1)
        else
          tmpa = maxval( thresh_factorbb_up(k-1:k) )
          tmpb = maxval( thresh_factorbb_down(k-1:k) )
        end if
        wup_thresh_k(  k,:) = wup_thresh_k(  k,:) * tmpa
        wdown_thresh_k(k,:) = wdown_thresh_k(k,:) * tmpb
      end do ! k

      do k=1, max(1, kup_top-1)
        wup_thresh(k) =  wup_thresh_k(k,1)
      end do
      do k=1, max(1, kdown_top-1)
        wdown_thresh(k) = wdown_thresh_k(k,1)
      end do

      do k=1, nzstag
        if(wup_thresh(k).lt.0.05) then
          write(0,*) 'erros in wup_thresh', k, wup_thresh_k(:,1), thresh_factorbb_up(:)
          call endrun('wup_thresh errors in ecpp_stat')
        end if
      end do
      !
      !  fix a bug in the WRF_ECPP, Minghuai Wang, 2009-12.
      !  set wdown_thresh_k and wup_thresh_k to be an extreme value
      !  above updraft (kup_top) and downdraft top(kdown_top).
      !  This will make sure there is no updraft or downdraft above kup_top and kdown_top
      !
      do k=kup_top, nz+1
        wup_thresh_k(k, :) = wlarge
      end do
      do k=kdown_top, nz+1
        wdown_thresh_k(k,:) = -1. * wlarge
      end do

      call setup_class_masks( &
      nx, ny, nz, nupdraft, ndndraft, ndraft_max, &
      cloudmixr, cf3d, precall, ww, &
      wdown_thresh_k, wup_thresh_k, &
      cloudthresh, prcpthresh, &
      mask_bnd, mask_cen,  &
      cloudmixr_total, cloudthresh_trans, precthresh_trans,  &
      qvs, precmixr_total  )

      !
      ! ( code added on 14-dec-2009 to guarantee quiescent class
      !    area > acen_quiesc_minaa )
      ! at each level
      !    calculate total fractional area for quiescent class
      !       using the current level-1 averages
      !    if (acen_quiesc < acen_quiesc_minaa), increase the
      !       thresh_factorbb_up/down(k) by factor of 1.5 or 1.2
      !    (also, if acen_down > acen_up, increase thresh_factorbb_up by less
      !
      thresh_calc_not_done = .false.
      do k = 1,nz
        acen_quiesc = sum( mask_cen( 1:nx, 1:ny, k, 1:NCLASS_CL, QUI, 1:NCLASS_PR) )
        acen_quiesc = max( acen_quiesc/real(nxy,crm_rknd), real(0.0,crm_rknd) )
        acen_up     = sum( mask_cen( 1:nx, 1:ny, k, 1:NCLASS_CL, UP1, 1:NCLASS_PR) )
        acen_up     = max( acen_up/real(nxy,crm_rknd)   , real(0.0,crm_rknd) )
        acen_down   = max( (1.0 - acen_quiesc - acen_up), real(0.0,crm_rknd) )

        abnd_quiesc = sum( mask_bnd( 1:nx, 1:ny, k, 1:NCLASS_CL, QUI, 1:NCLASS_PR) )
        abnd_quiesc = max( abnd_quiesc/real(nxy,crm_rknd), real(0.0,crm_rknd) )
        abnd_up     = sum( mask_bnd( 1:nx, 1:ny, k, 1:NCLASS_CL, UP1, 1:NCLASS_PR) )
        abnd_up     = max( abnd_up/real(nxy,crm_rknd)   , real(0.0,crm_rknd) )
        abnd_down   = max( (1.0 - abnd_quiesc - abnd_up), real(0.0,crm_rknd) )

        if (min(acen_quiesc, abnd_quiesc) < acen_quiesc_minaa) then
          thresh_calc_not_done = .true.
          if (acen_down > acen_up ) then
            tmpa = acen_up/acen_down
          else if (abnd_down > abnd_up ) then
            tmpa = abnd_up/abnd_down
          else
            tmpa = real(1.0,crm_rknd)
          end if
          if (min(acen_quiesc,abnd_quiesc) < 0.5*acen_quiesc_minaa) then
            thresh_factorbb_down(k) = thresh_factorbb_down(k)*1.5
            thresh_factorbb_up(k) = thresh_factorbb_up(k)*max(1.5*tmpa, real(1.25,crm_rknd) )
          else
            thresh_factorbb_down(k) = thresh_factorbb_down(k)*1.25
            thresh_factorbb_up(k) = thresh_factorbb_up(k)*max(1.25*tmpa, real(1.125,crm_rknd) )
          end if
          if(iter.gt.5) then
            write(0, *) 'warning: The number of iteration is larger than 5 in ecpp_stat', 'iter=', iter ,   &
            'acen_quiesc=', acen_quiesc, 'acen_up=', acen_up, 'k=', k,    &
            'wthreshdown=', wdown_thresh_k(k,1), 'wthreshup=', wup_thresh_k(k,1)
            !           call endrun('The number of iteration is larger than 10 in ecpp_stat')
          end if
        end if
      end do ! k

      !  thresh_calc_not_done = .false.   ! not use this iteration method  +++mhwang

    end do thresh_calc_loop

    wwqui_bar_cen(:) = 0.0
    wwqui_cloudy_bar_cen(:) = 0.0
    wwqui_bar_bnd(:) = 0.0
    wwqui_cloudy_bar_bnd(:) = 0.0

    XYCLASSLOOPS: do j = 1,ny
      do i = 1,nx
        do ipr = 1,NCLASS_PR
          do itr = 1,ndraft_max
            do icl = 1,NCLASS_CL
              !
              ! We now have enough information to aggregate the variables into domain
              ! averages by class. Do this first for the cell centers...
              !
              do k = 1,nz
                mask = mask_cen(i,j,k,icl,itr,ipr)/real(nxy,crm_rknd)

                area_cen_final(k,icl,itr,ipr) = area_cen_final(k,icl,itr,ipr) + mask

                if( domass ) then
                  area_cen_sum(k,icl,itr,ipr) = area_cen_sum(k,icl,itr,ipr) + mask
                  rh_cen_sum(k,icl,itr,ipr) = rh_cen_sum(k,icl,itr,ipr) + rh(i,j,k)*mask
                  qcloud_cen_sum(k,icl,itr,ipr) = qcloud_cen_sum(k,icl,itr,ipr) + qcloud(i,j,k)*mask
                  qcloud_bf_cen_sum(k,icl,itr,ipr) = qcloud_bf_cen_sum(k,icl,itr,ipr) + qcloud_bf(i,j,k)*mask
                  qrain_cen_sum(k,icl,itr,ipr) = qrain_cen_sum(k,icl,itr,ipr) + qrain(i,j,k)*mask
                  qice_cen_sum(k,icl,itr,ipr) = qice_cen_sum(k,icl,itr,ipr) + qice(i,j,k)*mask
                  qsnow_cen_sum(k,icl,itr,ipr) = qsnow_cen_sum(k,icl,itr,ipr) + qsnow(i,j,k)*mask
                  qgraup_cen_sum(k,icl,itr,ipr) = qgraup_cen_sum(k,icl,itr,ipr) + qgraup(i,j,k)*mask
                  qlsink_cen_sum(k,icl,itr,ipr) = qlsink_cen_sum(k,icl,itr,ipr) + qlsink(i,j,k)*mask
                  precr_cen_sum(k,icl,itr,ipr) = precr_cen_sum(k,icl,itr,ipr) + precr(i,j,k)*mask
                  precsolid_cen_sum(k,icl,itr,ipr) = precsolid_cen_sum(k,icl,itr,ipr) + precsolid(i,j,k)*mask
                  precall_cen_sum(k,icl,itr,ipr) = precall_cen_sum(k,icl,itr,ipr) + precall(i,j,k)*mask
                  qlsink_bf_cen_sum(k,icl,itr,ipr) = qlsink_bf_cen_sum(k,icl,itr,ipr) + qlsink_bf(i,j,k)*mask
                  prain_cen_sum(k,icl,itr,ipr) = prain_cen_sum(k,icl,itr,ipr) + prain(i,j,k)*mask
                  !
                  ! calculate the mean vertical velocity over the quiescent class  +++mhwang
                  !
                  if(itr.eq.QUI) then
                    wwqui_bar_cen(k) = wwqui_bar_cen(k)+(ww(i,j,k)+ww(i,j,k+1))*0.5*mask
                    if(icl.eq.CLD) then
                      wwqui_cloudy_bar_cen(k)=wwqui_cloudy_bar_cen(k)+(ww(i,j,k)+ww(i,j,k+1))*0.5*mask
                    end if
                  end if

                end if
              end do !k
              !
              ! Now, we can do a similar aggregation for the cell boundaries. Here, we
              ! will also calculate the mass flux and entrainment.
              !
              do k = 1,nzstag
                mask = mask_bnd(i,j,k,icl,itr,ipr)/real(nxy,crm_rknd)

                area_bnd_final(k,icl,itr,ipr) = area_bnd_final(k,icl,itr,ipr) + mask

                if( domass ) then
                  !NOTE: technically we should interpolate and not do a simple
                  !      average to get density at the cell interface
                  km0 = min(nz,k)
                  km1 = max(1,k-1)
                  km2 = max(1,k-2)
                  wwrho_k   = 0.5*(1.0/alt(i,j,km1) + 1.0/alt(i,j,km0))*ww(i,j,k)
                  wwrho_km1 = 0.5*(1.0/alt(i,j,km2) + 1.0/alt(i,j,km1))*ww(i,j,km1)

                  area_bnd_sum(k,icl,itr,ipr) = area_bnd_sum(k,icl,itr,ipr) + mask
                  mass_bnd_sum(k,icl,itr,ipr) = mass_bnd_sum(k,icl,itr,ipr) + wwrho_k*mask
                  ent_bnd_sum(k,icl,itr,ipr) = ent_bnd_sum(k,icl,itr,ipr) + max(real(0.,crm_rknd), wwrho_k-wwrho_km1)*mask

                  !
                  ! calculate the mean vertical velocity over the quiescent class  +++mhwang
                  !
                  if(itr.eq.QUI) then
                    wwqui_bar_bnd(k) = wwqui_bar_bnd(k)+ww(i,j,k)*mask
                    if(icl.eq.CLD) then
                      wwqui_cloudy_bar_bnd(k)=wwqui_cloudy_bar_bnd(k)+ww(i,j,k)*mask
                    end if
                  end if

                end if
              end do !k

            end do !icl
          end do !itr
        end do !pr
      end do !i
    end do XYCLASSLOOPS !j

    !
    ! calcualte vertical velocity variance for quiescent class (total and cloudy)  +++mhwang
    !
    do k=1, nz
      if(sum(mask_cen(1:nx, 1:ny, k, 1:NCLASS_CL, QUI, 1:NCLASS_PR)).ge.0.5) then
        wwqui_bar_cen(k) = wwqui_bar_cen(k)* real(nxy,crm_rknd) /sum(mask_cen(1:nx, 1:ny, k, 1:NCLASS_CL, QUI, 1:NCLASS_PR))
      else
        wwqui_bar_cen(k) = 0.0
      end if
      if(sum(mask_cen(1:nx, 1:ny, k, CLD, QUI, 1:NCLASS_PR)).ge.0.5) then
        wwqui_cloudy_bar_cen(k) = wwqui_cloudy_bar_cen(k)* real(nxy,crm_rknd) /sum(mask_cen(1:nx, 1:ny, k, CLD, QUI, 1:NCLASS_PR))
      else
        wwqui_cloudy_bar_cen(k) = 0.0
      end if
    end do
    do k=1, nzstag
      if(sum(mask_bnd(1:nx, 1:ny, k, 1:NCLASS_CL, QUI, 1:NCLASS_PR)).ge.0.5) then
        wwqui_bar_bnd(k) = wwqui_bar_bnd(k)* real(nxy,crm_rknd) /sum(mask_bnd(1:nx, 1:ny, k, 1:NCLASS_CL, QUI, 1:NCLASS_PR))
      else
        wwqui_bar_bnd(k) = 0.0
      end if
      if(sum(mask_bnd(1:nx, 1:ny, k, CLD, QUI, 1:NCLASS_PR)).ge.0.5) then
        wwqui_cloudy_bar_bnd(k) = wwqui_cloudy_bar_bnd(k)* real(nxy,crm_rknd) /sum(mask_bnd(1:nx, 1:ny, k, CLD, QUI, 1:NCLASS_PR))
      else
        wwqui_cloudy_bar_bnd(k) = 0.0
      end if
    end do

    QUIELOOPS: do j = 1,ny
      do i = 1,nx
        do ipr = 1,NCLASS_PR
          do icl = 1,NCLASS_CL

            do k = 1,nz
              mask = mask_cen(i,j,k,icl,QUI,ipr)/real(nxy,crm_rknd)

              !
              ! calculate the vertical velocity variance over the quiescent class  +++mhwang
              ! wwqui_bar_cen is used in for both all sky and cloudy sky.
              ! when wwqui_cloudy_bar_cen was used for cloudy sky, wwqui_cloudy_cen_sum will be smaller than wwqui_cen_sum.
              !
#ifdef CLUBB_CRM
              wwqui_cen_sum(k) = wwqui_cen_sum(k)+mask * ((ww(i,j,k)+ww(i,j,k+1))*0.5-wwqui_bar_cen(k))**2 + mask * (wwsq(i,j,k)+wwsq(i,j,k+1))**2/4.
#else
              wwqui_cen_sum(k) = wwqui_cen_sum(k)+mask * ((ww(i,j,k)+ww(i,j,k+1))*0.5-wwqui_bar_cen(k))**2 + mask * tkesgs(i,j,k)/3.
#endif
              if(icl.eq.CLD) then
#ifdef CLUBB_CRM
                wwqui_cloudy_cen_sum(k)=wwqui_cloudy_cen_sum(k)+mask * ((ww(i,j,k)+ww(i,j,k+1))*0.5-wwqui_bar_cen(k))**2 + mask * (wwsq(i,j,k)+wwsq(i,j,k+1))**2/4.
#else
                wwqui_cloudy_cen_sum(k)=wwqui_cloudy_cen_sum(k)+mask * ((ww(i,j,k)+ww(i,j,k+1))*0.5-wwqui_bar_cen(k))**2 + mask * tkesgs(i,j,k)/3.
#endif
              end if
            end do !k

            !
            ! Now, we can do a similar aggregation for the cell boundaries.
            !
            do k = 1,nzstag
              mask = mask_bnd(i,j,k,icl,QUI,ipr)/real(nxy,crm_rknd)

              !NOTE: technically we should interpolate and not do a simple
              !      average to get density at the cell interface
              km0 = min(nz,k)
              km1 = max(1,k-1)
              !
              ! calculate the mean vertical velocity over the quiescent class  +++mhwang
              ! wwqui_bar_bnd is used in both all sky and cloudy sky.
              ! when wwqui_cloudy_bar_bnd was used for cloudy sky, wwqui_cloudy_bnd_sum will be smaller than wwqui_bnd_sum.
              !
#ifdef CLUBB_CRM
              wwqui_bnd_sum(k) = wwqui_bnd_sum(k)+mask * (ww(i,j,k)-wwqui_bar_bnd(k))**2 + mask * wwsq(i,j,k)**2
#else
              wwqui_bnd_sum(k) = wwqui_bnd_sum(k)+mask * (ww(i,j,k)-wwqui_bar_bnd(k))**2 + mask * (tkesgs(i,j,km0)+tkesgs(i,j,km1)) * 0.5/3.
#endif
              if(icl.eq.CLD) then
#ifdef CLUBB_CRM
                wwqui_cloudy_bnd_sum(k)=wwqui_cloudy_bnd_sum(k)+mask * (ww(i,j,k)-wwqui_bar_bnd(k))**2 + mask * wwsq(i,j,k)**2
#else
                wwqui_cloudy_bnd_sum(k)=wwqui_cloudy_bnd_sum(k)+mask * (ww(i,j,k)-wwqui_bar_bnd(k))**2 + mask * (tkesgs(i,j,km0)+tkesgs(i,j,km1)) * 0.5/3.
#endif
              end if

            end do !k

          end do !icl
        end do !pr
      end do !i
    end do QUIELOOPS !j

    ! testing small queiscent fraction +++mhwang
    do k=1, nz
      if(sum(area_cen_final(k,:,1,:)).lt.1.0e-3) then
        write(0, *) 'ecpp, area_cen_final, quiescent', sum(area_cen_final(k,:,1,:)), k, area_cen_final(k,:,1,:), wdown_thresh_k(k,1), wup_thresh_k(k,1)
        write(0, *)  'ecpp, area_cen_final, quiescent, wwk', ww(:,:,k), i, wup_rms_k(k), wup_bar_k(k), wup_stddev_k(k)
        write(0, *) 'ecpp, area_cen_final, quiescent, wwk+1', ww(:,:,k+1), i, wup_rms_k(k+1), wup_bar_k(k+1), wup_stddev_k(k+1)
        !      call endrun('area_cen_final less  then 1.0-e3')
      end if
    end do
    ! ---mhwang
  end subroutine categorization_stats

  !------------------------------------------------------------------------
  subroutine determine_transport_thresh( &
    nx, ny, nz, &
    mode_updnthresh, upthresh, downthresh, &
    upthresh2, downthresh2, cloudthresh, &
    !     ctime, &
    ww, rhoair, &
    wdown_thresh_k, wup_thresh_k         &
    , cloudtop                           &
    , wup_rms_k, wup_bar_k, wup_stddev_k  &
    , wdown_rms_k, wdown_bar_k, wdown_stddev_k  &
    , kup_top, kdown_top)
    !
    ! Deterines the velocity thresholds used to indicate whether a cell's
    ! motion is up, down, or quiescent. This is down for two threshold values
    ! in each direction by level. A dozen options are available on how this
    ! is done as documented below and at the top of postproc_wrfout.
    !
    ! William.Gustafosn@pnl.gov; 11-Sep-2008
    ! Modified: William.Gustafosn@pnl.gov; 14-Apr-2009
    !------------------------------------------------------------------------
    !  use timeroutines
    !
    ! Soubroutine arguments...
    !
    integer, intent(in) :: nx, ny, nz, mode_updnthresh
    real(crm_rknd), intent(in) :: &
    cloudthresh, &
    downthresh, upthresh, &
    downthresh2, upthresh2
    !  type(time), intent(in) :: ctime
    real(crm_rknd), dimension(:,:,:), intent(in) :: &
    ww
    real(crm_rknd), dimension(nz+1), intent(in) :: rhoair
    real(crm_rknd), dimension(nz+1,2), intent(out) :: wdown_thresh_k, wup_thresh_k
    integer, dimension(nx,ny), intent(in) :: cloudtop
    real(crm_rknd), dimension(nz+1), intent(out) :: wup_rms_k, wup_bar_k, wup_stddev_k, wdown_bar_k, wdown_rms_k, wdown_stddev_k
    integer, intent(out) :: kup_top, kdown_top  ! defined as the maximum level that allows updraft and downdraft
    !
    ! Local vars...
    !
    real(crm_rknd), dimension(nz+1) :: &
    tmpveca, tmpvecb, &
    !       wdown_bar_k, wdown_rms_k, wdown_stddev_k, &
    !       wup_bar_k, wup_rms_k, wup_stddev_k, &
    wup_rms_ksmo, wdown_rms_ksmo
    real(crm_rknd) :: tmpsuma, tmpsumb, tmpw, tmpw_minval, &
    wdown_bar, wdown_rms, wdown_stddev, &
    wup_bar, wup_rms, wup_stddev
    integer, dimension(nx,ny) ::  &
    cloudtop_upaa, cloudtop_upbb, cloudtop_downaa, cloudtop_downbb
    integer, dimension(nz+1) :: nup_k, ndown_k
    integer :: i, ib, ic, &
    j, jb, jc, &
    k, kk, kup_center, kdown_center
    integer :: ndown, nup
    integer :: ijdel, ijdel_cur, ijdel_upaa, ijdel_upbb, ijdel_downaa, ijdel_downbb

    ! Calc cloudtop_upaa(i,j) = max( cloudtop(i-del:i+del,j-del:j+del) )
    ! and similar for cloudtop_upbb, cloudtop_downaa/bb
    ! (assume periodic BC here)
    ijdel_upaa = 0 ; ijdel_downaa = 0
    ijdel_upbb = 0 ; ijdel_downbb = 0
    if ((mode_updnthresh == 12) .or. (mode_updnthresh == 13)) then
      !    ijdel_... = 1 corresponds to 3x3 stencil
      ijdel_upaa = 1 ; ijdel_downaa = 1
      ijdel_upbb = 1 ; ijdel_downbb = 1
    end if
    ijdel = max( ijdel_upaa, ijdel_upbb, ijdel_downaa, ijdel_downbb )

    if (ijdel > 0) then
      do j = 1, ny
        do i = 1, nx
          cloudtop_upaa(i,j) = cloudtop(i,j)
          cloudtop_downaa(i,j) = cloudtop(i,j)
          cloudtop_upbb(i,j) = cloudtop(i,j)
          cloudtop_downbb(i,j) = cloudtop(i,j)
          do jb = j-ijdel, j+ijdel
            jc = jb
            if (jc <  1) jc = jc + ny
            if (jc > ny) jc = jc - ny
            do ib = i-ijdel, i+ijdel
              ic = ib
              if (ic <  1) ic = ic + nx
              if (ic > nx) ic = ic - nx
              ijdel_cur = max( iabs(ib-i), iabs(jb-j) )
              ! cloudtop_downaa calculated over a (2*ijdel_downaa+1)**2 stencil
              if (ijdel_cur <= ijdel_downaa) &
              cloudtop_downaa(i,j) = max( cloudtop_downaa(i,j), cloudtop(ic,jc) )
              ! cloudtop_upaa calculated over a (2*ijdel_upaa+1)**2 stencil
              if (ijdel_cur <= ijdel_upaa) &
              cloudtop_upaa(i,j) = max( cloudtop_upaa(i,j), cloudtop(ic,jc) )
              ! cloudtop_downbb, cloudtop_upbb similarly
              if (ijdel_cur <= ijdel_downbb) &
              cloudtop_downbb(i,j) = max( cloudtop_downbb(i,j), cloudtop(ic,jc) )
              if (ijdel_cur <= ijdel_upbb) &
              cloudtop_upbb(i,j) = max( cloudtop_upbb(i,j), cloudtop(ic,jc) )
            end do   ! ib
          end do      ! jb
          ! add on 1 level as a "margin of error"
          cloudtop_upaa(  i,j) = min( cloudtop_upaa(  i,j)+1, nz )
          cloudtop_downaa(i,j) = min( cloudtop_downaa(i,j)+1, nz )
          cloudtop_upbb(  i,j) = min( cloudtop_upbb(  i,j)+1, nz )
          cloudtop_downbb(i,j) = min( cloudtop_downbb(i,j)+1, nz )
        end do   ! i
      end do   ! j
    end if   ! (ijdel > 0)

    ! new coding here and below
    !   cloudtop_up/downaa - only grid cells with k<=cloudtop_up/downaa
    !                        are used for calc of wup_rms and wdn_rms
    !   cloudtop_up/downbb - only grid cells with k<=cloudtop_up/downbb
    !                        can be classified as up/downdraft
    if ((mode_updnthresh == 12) .or. (mode_updnthresh == 13)) then
      ! mode_updnthresh >= 12 is a newer, more consistent usage of cloudtop info
      ! the cloudtop_upaa/upbb/downaa/downbb values are identical,
      !    and they correspond to the max cloudtop(i,j) over a 3x3 stencil
      ! only grid cells with k <= this "local" cloudtop can be up/downdraft grids
      continue
    else
      ! mode_updnthresh /= 12,13 corresponds to pre 11-jan-2008 versions of preprocessor
      !   where only grid cells with k <= cloudtop(i,j) are used for calc of wup/dn_rms,
      !   but any grid cells can be up/dn [even those with k >> cloudtop(i,j)]
      cloudtop_upaa(:,:)   = cloudtop(:,:)
      cloudtop_downaa(:,:) = cloudtop(:,:)
      cloudtop_upbb(:,:)   = nz
      cloudtop_downbb(:,:) = nz
    end if

    !
    ! Get standard deviation of up and down vertical velocity below the
    ! cloud tops. For now, each cell is treated equally. We may want to
    ! consider weighting each cell by its volume or mass.
    !
    ! Get the mean values first for wup and wdown
    ndown     = 0;   nup       = 0
    wdown_bar = 0.;  wup_bar   = 0.
    ndown_k(:)     = 0;   nup_k(:)       = 0
    wdown_bar_k(:) = 0.;  wup_bar_k(:)   = 0.
    kup_top = 1; kdown_top= 1
    do j=1,ny
      do i=1,nx
        do k=1,cloudtop_upaa(i,j)+1 !Plus 1 is so we get w across top of cloud.
          !It is dimmensionally ok since w is dimmed nz+1
          !We intentially ignore when w==0 as to not bias one direction
          !over the other for the count. This differs from the Ferret code which
          !assigns w=0 to up values.
          if( ww(i,j,k) > 0. ) then
            nup       = nup + 1
            wup_bar   = wup_bar + ww(i,j,k)
            nup_k(k)       = nup_k(k) + 1
            wup_bar_k(k)   = wup_bar_k(k) + ww(i,j,k)
            kup_top   = max(kup_top, k)
          end if
        end do
        do k=1,cloudtop_downaa(i,j)+1
          if( ww(i,j,k) < 0. ) then
            ndown     = ndown + 1
            wdown_bar = wdown_bar + ww(i,j,k)
            ndown_k(k)     = ndown_k(k) + 1
            wdown_bar_k(k) = wdown_bar_k(k) + ww(i,j,k)
            kdown_top  = max(kdown_top, k)
          end if
        end do

      end do
    end do
    if( nup > 0 )   wup_bar   = wup_bar / nup
    if( ndown > 0 ) wdown_bar = wdown_bar / ndown
    do k = 1, nz+1
      if( nup_k(k) > 0 )   wup_bar_k(k)   = wup_bar_k(k) / nup_k(k)
      if( ndown_k(k) > 0 ) wdown_bar_k(k) = wdown_bar_k(k) / ndown_k(k)
    end do

    !Now, we can get the std. dev. of wup and wdown.
    wdown_stddev      = 0.;  wup_stddev      = 0.
    wdown_stddev_k(:) = 0.;  wup_stddev_k(:) = 0.
    do j=1,ny
      do i=1,nx
        do k=1,cloudtop_upaa(i,j)+1 !Plus 1 is so we get w across top of cloud.
          !We intentionally ignore when w==0 as to not bias one direction
          !over the other.
          if( ww(i,j,k) > 0. ) then
            wup_stddev   = wup_stddev + (wup_bar-ww(i,j,k))**2
            wup_stddev_k(k)   = wup_stddev_k(k) + (wup_bar_k(k)-ww(i,j,k))**2
          end if
        end do
        do k=1,cloudtop_downaa(i,j)+1
          if( ww(i,j,k) < 0. ) then
            wdown_stddev = wdown_stddev + (wdown_bar-ww(i,j,k))**2
            wdown_stddev_k(k) = wdown_stddev_k(k) + (wdown_bar_k(k)-ww(i,j,k))**2
          end if
        end do
      end do
    end do
    if( nup > 0 )   wup_stddev   = sqrt(wup_stddev / nup)
    if( ndown > 0 ) wdown_stddev = sqrt(wdown_stddev / ndown)
    wup_rms = sqrt( wup_bar**2 + wup_stddev**2 )
    wdown_rms = sqrt( wdown_bar**2 + wdown_stddev**2 )
    do k = 1, nz+1
      if( nup_k(k) > 0 )   wup_stddev_k(k)   = sqrt(wup_stddev_k(k) / nup_k(k))
      if( ndown_k(k) > 0 ) wdown_stddev_k(k) = sqrt(wdown_stddev_k(k) / ndown_k(k))
      wup_rms_k(k) = sqrt( wup_bar_k(k)**2 + wup_stddev_k(k)**2 )
      wdown_rms_k(k) = sqrt( wdown_bar_k(k)**2 + wdown_stddev_k(k)**2 )
    end do

    ! calculated smoothed (3-point) wup/down_rms
    tmpveca(:) = wup_rms_k(  :)
    tmpvecb(:) = wdown_rms_k(:)
    do k = 2, nz
      wup_rms_ksmo(  k) = 0.0
      wdown_rms_ksmo(k) = 0.0
      tmpsuma = 0.0
      do kk = max(k-1,2), min(k+1,nz)
        wup_rms_ksmo(  k) = wup_rms_ksmo(  k) + tmpveca(kk)
        wdown_rms_ksmo(k) = wdown_rms_ksmo(k) + tmpvecb(kk)
        tmpsuma = tmpsuma + 1.0
      end do
      tmpsuma = max(tmpsuma,real(1.0,crm_rknd))
      wup_rms_ksmo(  k) = wup_rms_ksmo(  k)/tmpsuma
      wdown_rms_ksmo(k) = wdown_rms_ksmo(k)/tmpsuma
    end do
    wup_rms_ksmo(  1) = wup_rms_ksmo(  2)
    wdown_rms_ksmo(1) = wdown_rms_ksmo(2)
    wup_rms_ksmo(  nz+1) = wup_rms_ksmo(  nz)
    wdown_rms_ksmo(nz+1) = wdown_rms_ksmo(nz)

    !  print "(2a,2(2x,3f8.4))", &
    !     " ...wup_bar,std,rms;  wdown_bar,std,rms  ",   &
    !     wup_bar, wup_stddev, wup_rms, wdown_bar, wdown_stddev, wdown_rms
    !  if (mode_updnthresh >= 5) then
    !     print "(a/(15f7.3))", &
    !     " ...  wup_rms_k(2:nz)", (wup_rms_k(k), k=2,nz)
    !     print "(a/(15f7.3))", &
    !     " ...wdown_rms_k(2:nz)", (wdown_rms_k(k), k=2,nz)
    !  end if

    !
    ! Get masks to determine (cloud vs. clear) (up vs. down vs. other) categories.
    ! Vertical velocities are checked on the cell vertical interfaces to determine
    ! if they pass the threshold criteria. Clouds below the interface are then
    ! used for updrafts and above the int. for downdrafts. Quiescent (other)
    ! drafts use an average of the cloud above and below the interface to
    ! determine cloudiness.
    !
    select case ( mode_updnthresh )
    case ( 1 )
      wup_thresh_k(  :,1) =  wup_stddev*abs(upthresh)
      wdown_thresh_k(:,1) = -wdown_stddev*abs(downthresh)
      wup_thresh_k(  :,2) =  wup_stddev*abs(upthresh2)
      wdown_thresh_k(:,2) = -wdown_stddev*abs(downthresh2)
    case ( 2 )
      wup_thresh_k(  :,1) = wup_bar   + wup_stddev*abs(upthresh)
      wdown_thresh_k(:,1) = wdown_bar - wdown_stddev*abs(downthresh)
      wup_thresh_k(  :,2) = wup_bar   + wup_stddev*abs(upthresh2)
      wdown_thresh_k(:,2) = wdown_bar - wdown_stddev*abs(downthresh2)
    case ( 3 )
      wup_thresh_k(  :,1) =  abs(upthresh)
      wdown_thresh_k(:,1) = -abs(downthresh)
      wup_thresh_k(  :,2) =  abs(upthresh2)
      wdown_thresh_k(:,2) = -abs(downthresh2)
    case ( 4 )
      wup_thresh_k(  :,1) =  (wup_rms  )*abs(upthresh)
      wdown_thresh_k(:,1) = -(wdown_rms)*abs(downthresh)
      wup_thresh_k(  :,2) =  (wup_rms  )*abs(upthresh2)
      wdown_thresh_k(:,2) = -(wdown_rms)*abs(downthresh2)

    case ( 5 )
      ! For mode_updnthresh = 5, use a weighted average of wup_rms & wup_rms_ksmo(k)
      ! because wup_rms_ksmo will be zero (or close to it) at many levels
      wup_thresh_k(  :,1) =  (0.25*wup_rms  +0.75*wup_rms_ksmo(  :))*abs(upthresh)
      wdown_thresh_k(:,1) = -(0.25*wdown_rms+0.75*wdown_rms_ksmo(:))*abs(downthresh)
      wup_thresh_k(  :,2) =  (0.25*wup_rms  +0.75*wup_rms_ksmo(  :))*abs(upthresh2)
      wdown_thresh_k(:,2) = -(0.25*wdown_rms+0.75*wdown_rms_ksmo(:))*abs(downthresh2)

    case ( 6, 7 )
      ! For mode_updnthresh = 6 & 7, like case 4 except when k <= "updraft center k",
      ! use minimum of wup_rms and wup_rms_k for updraft threshold
      wup_thresh_k(  :,1) =  (wup_rms  )*abs(upthresh)
      wdown_thresh_k(:,1) = -(wdown_rms)*abs(downthresh)
      wup_thresh_k(  :,2) =  (wup_rms  )*abs(upthresh2)
      wdown_thresh_k(:,2) = -(wdown_rms)*abs(downthresh2)

      tmpsuma = 0.0 ; tmpsumb = 1.0e-30
      do k = 1, nz
        tmpw = wup_rms_k(k)
        if (mode_updnthresh == 7) tmpw = wup_rms_ksmo(k)
        tmpw = max(real(1.0e-4,crm_rknd),tmpw)
        tmpw = tmpw * rhoair(k)
        tmpsuma = tmpsuma + tmpw*k ; tmpsumb = tmpsumb + tmpw
      end do
      kup_center = nint(tmpsuma/tmpsumb)
      tmpw_minval = 0.10
      do k = 1, kup_center
        tmpw = wup_rms_k(k)
        if (mode_updnthresh == 7) tmpw = wup_rms_ksmo(k)
        tmpw = max( tmpw, tmpw_minval )
        tmpw = min( tmpw, wup_rms )
        wup_thresh_k(k,1) = tmpw*abs(upthresh)
        wup_thresh_k(k,2) = tmpw*abs(upthresh2)
      end do

    case ( 8, 9 )
      ! For mode_updnthresh = 8 & 9, like case 6, 7 except that updraft and
      ! downdraft are treated similarly.  So when k >= "downdraft center k",
      ! use minimum of wdown_rms and wdown_rms_k for downdraft threshold
      wup_thresh_k(  :,1) =  (wup_rms  )*abs(upthresh)
      wdown_thresh_k(:,1) = -(wdown_rms)*abs(downthresh)
      wup_thresh_k(  :,2) =  (wup_rms  )*abs(upthresh2)
      wdown_thresh_k(:,2) = -(wdown_rms)*abs(downthresh2)

      tmpsuma = 0.0 ; tmpsumb = 1.0e-30
      do k = 1, nz+1
        tmpw = wup_rms_k(k)
        if (mode_updnthresh == 9) tmpw = wup_rms_ksmo(k)
        tmpw = max(real(1.0e-4,crm_rknd),tmpw)
        tmpw = tmpw * rhoair(k)
        tmpsuma = tmpsuma + tmpw*k ; tmpsumb = tmpsumb + tmpw
      end do
      kup_center = nint(tmpsuma/tmpsumb)
      tmpw_minval = 0.10
      do k = 1, kup_center
        tmpw = wup_rms_k(k)
        if (mode_updnthresh == 9) tmpw = wup_rms_ksmo(k)
        tmpw = max( tmpw, tmpw_minval )
        tmpw = min( tmpw, wup_rms )
        wup_thresh_k(k,1) = tmpw*abs(upthresh)
        wup_thresh_k(k,2) = tmpw*abs(upthresh2)
      end do

      tmpsuma = 0.0 ; tmpsumb = 1.0e-30
      do k = 1, nz+1
        tmpw = wdown_rms_k(k)
        if (mode_updnthresh == 9) tmpw = wdown_rms_ksmo(k)
        tmpw = max(real(1.0e-4,crm_rknd),tmpw)
        tmpw = tmpw * rhoair(k)
        tmpsuma = tmpsuma + tmpw*k ; tmpsumb = tmpsumb + tmpw
      end do
      kdown_center = nint(tmpsuma/tmpsumb)
      tmpw_minval = 0.10
      do k = kdown_center, nz+1
        tmpw = wdown_rms_k(k)
        if (mode_updnthresh == 9) tmpw = wdown_rms_ksmo(k)
        tmpw = max( tmpw, tmpw_minval )
        tmpw = min( tmpw, wdown_rms )
        wdown_thresh_k(k,1) = -tmpw*abs(downthresh)
        wdown_thresh_k(k,2) = -tmpw*abs(downthresh2)
      end do

    case ( 14, 15 )
      ! case 14 & 15 -- added on 10-dec-2009
      !    updraft   and k  > "updraft   center k",  wup_rms
      !    updraft   and k <= "updraft   center k",  use min( wup_rms_k, wup_rms )
      !    downdraft and k  > "downdraft center k",  wdown_rms
      !    downdraft and k <= "downdraft center k",  min( use wdown_rms_k, wdown_rms )
      ! The idea is to have a higher threshold in upper troposphere to
      ! filter out gravity waves motions
      tmpsuma = 0.0 ; tmpsumb = 1.0e-30
      do k = 1, nz+1
        tmpw = wup_rms_k(k)
        if (mode_updnthresh == 15) tmpw = wup_rms_ksmo(k)
        tmpw = max(real(1.0e-4,crm_rknd),tmpw)
        tmpw = tmpw * rhoair(k)
        tmpsuma = tmpsuma + tmpw*k ; tmpsumb = tmpsumb + tmpw
      end do
      kup_center = nint(tmpsuma/tmpsumb)
      tmpw_minval = 0.10
      do k = 1, nz+1
        tmpw = wup_rms_k(k)
        if (mode_updnthresh == 15) tmpw = wup_rms_ksmo(k)
        if (k > kup_center) then
          tmpw = wup_rms
        else
          tmpw = min( tmpw, wup_rms )
        end if
        tmpw = max( tmpw, tmpw_minval )
        wup_thresh_k(k,1) = tmpw*abs(upthresh)
        wup_thresh_k(k,2) = tmpw*abs(upthresh2)
      end do

      tmpsuma = 0.0 ; tmpsumb = 1.0e-30
      do k = 1, nz+1
        tmpw = wdown_rms_k(k)
        if (mode_updnthresh == 15) tmpw = wdown_rms_ksmo(k)
        tmpw = max(real(1.0e-4,crm_rknd),tmpw)
        tmpw = tmpw * rhoair(k)
        tmpsuma = tmpsuma + tmpw*k ; tmpsumb = tmpsumb + tmpw
      end do
      kdown_center = nint(tmpsuma/tmpsumb)
      tmpw_minval = 0.10
      do k = 1, nz+1
        tmpw = wdown_rms_k(k)
        if (mode_updnthresh == 15) tmpw = wdown_rms_ksmo(k)
        if (k > kdown_center) then
          tmpw = wdown_rms
        else
          tmpw = min( tmpw, wdown_rms )
        end if
        tmpw = max( tmpw, tmpw_minval )
        wdown_thresh_k(k,1) = -tmpw*abs(downthresh)
        wdown_thresh_k(k,2) = -tmpw*abs(downthresh2)
      end do

    case ( 16, 17 )
      ! case 16 & 17 -- added on 10-dec-2009
      !    updraft   and k  > "updraft   center k",  use max( wup_rms_k, wup_rms )
      !    updraft   and k <= "updraft   center k",  use wup_rms_k
      !    downdraft and k  > "downdraft center k",  use max( wdown_rms_k, wdown_rms )
      !    downdraft and k <= "downdraft center k",  use wdown_rms_k
      ! The idea is to have a higher threshold in upper troposphere to
      ! filter out gravity waves motions
      tmpsuma = 0.0 ; tmpsumb = 1.0e-30
      do k = 1, nz+1
        tmpw = wup_rms_k(k)
        if (mode_updnthresh == 17) tmpw = wup_rms_ksmo(k)
        tmpw = max(real(1.0e-4,crm_rknd),tmpw)
        tmpw = tmpw * rhoair(k)
        tmpsuma = tmpsuma + tmpw*k ; tmpsumb = tmpsumb + tmpw
      end do
      kup_center = nint(tmpsuma/tmpsumb)
      tmpw_minval = 0.10
      do k = 1, nz+1
        tmpw = wup_rms_k(k)
        if (mode_updnthresh == 17) tmpw = wup_rms_ksmo(k)
        if (k > kup_center) tmpw = max( tmpw, wup_rms )
        tmpw = max( tmpw, tmpw_minval )
        wup_thresh_k(k,1) = tmpw*abs(upthresh)
        wup_thresh_k(k,2) = tmpw*abs(upthresh2)
      end do

      tmpsuma = 0.0 ; tmpsumb = 1.0e-30
      do k = 1, nz+1
        tmpw = wdown_rms_k(k)
        if (mode_updnthresh == 17) tmpw = wdown_rms_ksmo(k)
        tmpw = max(real(1.0e-4,crm_rknd),tmpw)
        tmpw = tmpw * rhoair(k)
        tmpsuma = tmpsuma + tmpw*k ; tmpsumb = tmpsumb + tmpw
      end do
      kdown_center = nint(tmpsuma/tmpsumb)
      tmpw_minval = 0.10
      do k = 1, nz+1
        tmpw = wdown_rms_k(k)
        if (mode_updnthresh == 17) tmpw = wdown_rms_ksmo(k)
        if (k > kdown_center) tmpw = max( tmpw, wdown_rms )
        tmpw = max( tmpw, tmpw_minval )
        wdown_thresh_k(k,1) = -tmpw*abs(downthresh)
        wdown_thresh_k(k,2) = -tmpw*abs(downthresh2)
      end do

    case ( 10, 11, 12, 13 )
      ! For mode_updnthresh = 10, 11, use wup_rms_k and wdown_rms_k at all
      ! levels (or the w---_rms_ksmo)
      tmpw_minval = 0.10
      do k = 1, nz+1
        tmpw = wup_rms_k(k)
        if (mode_updnthresh == 11) tmpw = wup_rms_ksmo(k)
        if (mode_updnthresh == 13) tmpw = wup_rms_ksmo(k)
        tmpw = max( tmpw, tmpw_minval )
        wup_thresh_k(k,1) = tmpw*abs(upthresh)
        wup_thresh_k(k,2) = tmpw*abs(upthresh2)
      end do
      tmpw_minval = 0.10
      do k = 1, nz+1
        tmpw = wdown_rms_k(k)
        if (mode_updnthresh == 11) tmpw = wdown_rms_ksmo(k)
        if (mode_updnthresh == 13) tmpw = wdown_rms_ksmo(k)
        tmpw = max( tmpw, tmpw_minval )
        wdown_thresh_k(k,1) = -tmpw*abs(downthresh)
        wdown_thresh_k(k,2) = -tmpw*abs(downthresh2)
      end do

    case default
      call endrun('determine_transport_thresh error - must have   1 <= mode_updnthresh <= 11')
    end select

  end subroutine determine_transport_thresh


  !------------------------------------------------------------------------
  subroutine setup_class_masks( &
    nx, ny, nz, nupdraft, ndndraft, ndraft_max, &
    cloudmixr, cf3d, precall, ww, &
    wdown_thresh_k, wup_thresh_k, &
    cloudthresh, prcpthresh, &
    mask_bnd, mask_cen,  &
    cloudmixr_total, cloudthresh_trans, precthresh_trans, &
    qvs, precmixr_total )
    !
    ! Sets up the masks used for determining quiescent/up/down, clear/cloudy,
    ! and non-precipitatin/precipitating classes.
    !
    ! William.Gustafosn@pnl.gov; 20-Nov-2008
    ! Last modified: William.Gustafson@pnl.gov; 16-Apr-2009

    ! Modification by Minghuai Wang (Minghuai.Wang@pnl.gov), April 23, 2010
    ! use total condensate (liquid+ice),  different condensate and precipitating thresholds
    ! to classify transport classes.
    ! See Xu et al., 2002, Q.J.R.M.S.
    !

    !------------------------------------------------------------------------
    !
    ! Soubroutine arguments...
    !
    integer, intent(in) :: nx, ny, nz, nupdraft, ndndraft, ndraft_max
    real(crm_rknd), dimension(:,:,:), intent(in) :: &
    cloudmixr, cf3d, precall, ww
    real(crm_rknd), dimension(nz+1,2), intent(in) :: wdown_thresh_k, wup_thresh_k
    real(crm_rknd), intent(in) :: cloudthresh, prcpthresh
    real(crm_rknd), dimension(nx,ny,nz+1,NCLASS_CL,ndraft_max,NCLASS_PR), &
    intent(out) :: mask_bnd
    real(crm_rknd), dimension(nx,ny,nz,NCLASS_CL,ndraft_max,NCLASS_PR), &
    intent(out) :: mask_cen
    real(crm_rknd), dimension( :, :, :), intent(in) :: cloudmixr_total   ! total condensate (liquid+ice)
    real(crm_rknd), intent(in) :: cloudthresh_trans, precthresh_trans  ! threshold for transport classes
    real(crm_rknd), dimension( :, :, :), intent(in) :: qvs, precmixr_total
    !
    ! Local vars...
    !
    integer, dimension(nz+1,nupdraft) :: maskup
    integer, dimension(nz+1,ndndraft) :: maskdn
    integer, dimension(nz+1) :: maskqu, &
    maskcld_bnd, maskclr_bnd, maskpry_bnd, maskprn_bnd
    integer, dimension(nz) :: maskcld, maskclr, maskpry, maskprn
    integer :: i, itr, icl, ipr, j, k, m, nzstag
    real(crm_rknd)  :: cloudthresh_trans_temp, precthresh_trans_temp

    nzstag = nz+1
    !
    ! Initialize the masks to zero and then we will accumulate values into
    ! them as we identify the various classes.
    !
    mask_bnd = 0.
    mask_cen = 0.
    !
    ! Loop over the horizontal dimensions...
    !
    XYLOOPS : do j = 1,ny
      do i=1,nx
        !
        ! Set initial mask values for the vertical cell boundaries...
        !
        maskup = 0
        maskdn = 0
        maskqu = 0
        maskcld = 0
        maskclr = 0
        maskcld_bnd = 0
        maskclr_bnd = 0
        maskpry = 0
        maskprn = 0
        maskpry_bnd = 0
        maskprn_bnd = 0

        if( nupdraft > 2 .or. ndndraft > 2 ) then
          call endrun('OOPS. Cannot have more than 2 updraft or 2 downdraft categories right now.')
        end if

        do k = 1,nzstag

          !Transport upward at cell boundaries...
          !We have to take into account the possibility of multiple
          !updraft categories. At this point, we handle only the
          !cases of one or two categories. We do not yet handle the
          !allcomb option.
          !
          ! updraft only exist in cloudy area or precipitating clear area ++++mhwang
          cloudthresh_trans_temp = cloudthresh_trans
          !           cloudthresh_trans_temp = max(cloudthresh_trans, 0.01 * (qvs(i,j,max(k-1,1))+qvs(i,j,min(k,nz)))*0.5)
          if( (cloudmixr_total(i,j,max(k-1,1))+cloudmixr_total(i,j,min(k,nz)))*0.5 > cloudthresh_trans_temp  &
          !               .or. (precall(i,j,max(k-1,1))+precall(i,j,min(k,nz)))*0.5 > prcpthresh_trans) then   !+++mhwang
          .or. (precmixr_total(i,j,max(k-1,1))+precmixr_total(i,j,min(k,nz)))*0.5 > precthresh_trans) then   !+++mhwang
          select case (nupdraft)
          case (1) !Only one threshold
            if( ww(i,j,k) > wup_thresh_k(k,1) ) then
              maskup(k,1) = 1
            end if
          case (2) !Two thresholds, assumes 1st is stronger wind
            if( ww(i,j,k) > wup_thresh_k(k,1) ) then
              maskup(k,1) = 1
            else if( ww(i,j,k) > wup_thresh_k(k,2) &
              .and. ww(i,j,k) <= wup_thresh_k(k,1) ) then
              maskup(k,2) = 1
            end if
          end select
        end if  ! end cloudmixr_total    +++mhwang

        !Transport downward at cell boundaries...
        !
        ! downdraft only exist in cloudy area or precipitating clear area   +++mhwang
        if( (cloudmixr_total(i,j,max(k-1,1))+cloudmixr_total(i,j,min(k,nz)))*0.5 > cloudthresh_trans_temp   &
        !               .or. (precall(i,j,max(k-1,1))+precall(i,j,min(k,nz)))*0.5 > prcpthresh_trans) then   !+++mhwang
        .or. (precmixr_total(i,j,max(k-1,1))+precmixr_total(i,j,min(k,nz)))*0.5 > precthresh_trans) then   !+++mhwang
        select case (ndndraft)
        case (1) !Only one threshold
          if( ww(i,j,k) < wdown_thresh_k(k,1) ) then
            maskdn(k,1) = 1
          end if
        case (2) !Two thresholds, assumes 1st is stronger wind
          if( ww(i,j,k) < wdown_thresh_k(k,1) ) then
            maskdn(k,1) = 1
          else if( ww(i,j,k) < wdown_thresh_k(k,2) &
            .and. ww(i,j,k) >= wdown_thresh_k(k,1) ) then
            maskdn(k,2) = 1
          end if
        end select
      end if  ! end cloudmixr_total, and precall   +++mhwang

      !Transport quiescent at cell boundaries if neither up or
      !down triggered...
      if( sum(maskup(k,:))+sum(maskdn(k,:)) < 1 ) then
        maskqu(k) = 1
      end if

      ! Cloudy or clear at cell boundaries...
      if( (cloudmixr(i,j,max(k-1,1))+cloudmixr(i,j,min(k,nz)))*0.5 > cloudthresh ) then
        maskcld_bnd(k) = 1
      else
        maskclr_bnd(k) = 1
      end if

      ! Raining or not at cell boundaries...
      if( (precall(i,j,max(k-1,1))+precall(i,j,min(k,nz)))*0.5 > prcpthresh ) then
        maskpry_bnd(k) = 1
      else
        maskprn_bnd(k) = 1
      end if

    end do !k
    do k = 1,nz

      ! Cloudy or clear at cell centers...
      if( cloudmixr(i,j,k) > cloudthresh ) then
        maskcld(k) = 1
      else
        maskclr(k) = 1
      end if

      ! Raining or not at cell centers...
      if( precall(i,j,k) > prcpthresh ) then
        maskpry(k) = 1
      else
        maskprn(k) = 1
      end if

    end do !k
    !
    ! Now, use the initial boundary masks by class to generate a combined
    ! mask for the cell boundaries.
    !
    do k = 1,nzstag

      !Upward, or at least upward quiescent
      if( sum(maskup(k,:)) > 0 .or. &
      (maskqu(k) > 0 .and. ww(i,j,k) > 0) ) then

      !Are we are here because of maskup? If so, then we need to
      !parse the correct updraft category.
      if( maskqu(k) < 1 ) then
        itr = UP1 + maxloc(maskup(k,:),1)-1
      else
        itr = QUI
      end if

      !For upward motion, determine cloud and precip characteristics
      !based on the cell-center values below the boundary.
      if( k==1 ) then
        icl = CLR
        ipr = PRN
      else
        call cloud_prcp_check(maskcld, CLD, maskclr, CLR, k-1, icl, &
        "setup_class_masks: bnd cloud up")
        call cloud_prcp_check(maskpry, PRY, maskprn, PRN, k-1, ipr, &
        "setup_class_masks: bnd prcp up")
      end if

      !Downward, or at least downward quiescent
    else if( sum(maskdn(k,:)) > 0 .or. &
      (maskqu(k) > 0 .and. ww(i,j,k) < 0) ) then

      !Are we here because of maskdn? If so, then we need to
      !parse the correct downdraft category.
      if( maskqu(k) < 1 ) then
        itr = DN1 + maxloc(maskdn(k,:),1)-1
      else
        itr = QUI
      end if

      !For downward motion, determine cloud and precip characteristics
      !based on the cell-center values above the boundary.
      if( k==nzstag ) then
        icl = CLR
        ipr = PRN
      else
        call cloud_prcp_check(maskcld, CLD, maskclr, CLR, k, icl, &
        "setup_class_masks: bnd cloud down")
        call cloud_prcp_check(maskpry, PRY, maskprn, PRN, k, ipr, &
        "setup_class_masks: bnd prcp down")
      end if

      !Quiescent with w=0. Use the cell-center values averaged
      !surrounding the boundary for the cloud/prcp states.
    else
      itr = QUI
      call cloud_prcp_check(maskcld_bnd, CLD, maskclr_bnd, CLR, k, icl, &
      "setup_class_masks: bnd cloud quiescent")
      call cloud_prcp_check(maskpry_bnd, PRY, maskprn_bnd, PRN, k, ipr, &
      "setup_class_masks: bnd prcp quiescent")
    end if

    ! +++mhwang
    ! Total condensate and different thresholds are used to classify transport classes. So the following change
    ! is not needed anymore. Minghuai Wang, 2010-04-23.
    !
    ! In the clear, and non-precipitating class, it is classified as quiescent class in the MMF simulation.
    ! If this is classed as updraft or downdraft in mode 16, this would lead to too much upraft and downdraft mass fluxes.
    ! Minghuai Wang, 2010-01-18 (Minghuai.Wang@pnl.gov)
    !           if(icl.eq.CLR .and. ipr.eq.PRN) then
    !             itr = QUI
    !           end if
    !---mhwang

    !We have all the class indices determined so now we can set
    !the correct mask location to 1.
    !           mask_bnd(i,j,k,icl,itr,ipr) = 1.
    ! use fractioal cloudiness in SAM
    if(icl.eq.CLR) then
      mask_bnd(i,j,k,icl,itr,ipr) = 1.
    else if(icl.eq.CLD) then
      mask_bnd(i,j,k,CLD,itr,ipr) = (cf3d(i,j,max(k-1,1))+cf3d(i,j,min(k, nz)))*0.5
      mask_bnd(i,j,k,CLR,itr,ipr) = 1. - (cf3d(i,j,max(k-1,1))+cf3d(i,j,min(k, nz)))*0.5
    end if


  end do !k-loop mask for boundaries
  !
  ! Now, use the initial boundary masks by class to generate a combined
  ! mask for the cell centers. We determine the transport class based on
  ! splitting the cell conceptually in half with the upper boundary
  ! influencing the top half of the cell and the bottom boundary the bottom
  ! half. Each contributes either 0 or 0.5 of the total contribution of the
  ! cell's transport. e.g. if both boundaries are upward, then the cell is
  ! fully an "up" transport cell. If the two boundaries are opposite, then
  ! the cell is weighted half in each direction for the masking.
  !
  do k = 1,nz

    !Get the cloud/prcp characteristics at cell center.
    call cloud_prcp_check(maskcld, CLD, maskclr, CLR, k, icl)
    call cloud_prcp_check(maskpry, PRY, maskprn, PRN, k, ipr)

    !Look at the bottom boundary first and determine it's
    !contribution to the cell center transport class.
    if( sum(maskup(k,:)) > 0 ) then
      itr = UP1 + maxloc(maskup(k,:),1)-1
    else if( sum(maskdn(k,:)) > 0 ) then
      itr = DN1 + maxloc(maskdn(k,:),1)-1
    else if( maskqu(k) > 0 ) then
      itr = QUI
    else
      call endrun("ERROR: setup_class_masks: We should not be in this place for cell bottoms.")
      stop
    end if

    ! +++mhwang
    ! ! Total condensate and different thresholds are used to classify transport classes. So the following change
    ! is not needed anymore. Minghuai Wang, 2010-04-23.

    ! In the clear, and non-precipitating class, it is classified as quiescent class in the MMF simulation.
    ! If this is classed as updraft or downdraft in mode 16, this would lead to too much upraft and downdraft mass fluxes.
    ! Minghuai Wang, 2010-01-18 (Minghuai.Wang@pnl.gov)
    !           if(icl.eq.CLR .and. ipr.eq.PRN) then
    !             itr = QUI
    !           end if
    !---mhwang

    !We have what we need for the cell bottom classes so increment
    !the center mask for the bottom half...
    !           mask_cen(i,j,k,icl,itr,ipr) = mask_cen(i,j,k,icl,itr,ipr) + 0.5
    ! Use fractional cloudiness at SAM
    if(icl.eq.CLR) then
      mask_cen(i,j,k,icl,itr,ipr) = mask_cen(i,j,k,icl,itr,ipr) + 0.5
    else if(icl.eq.CLD) then
      mask_cen(i,j,k,CLD,itr,ipr) =  mask_cen(i,j,k,CLD,itr,ipr) + (cf3d(i,j,k))*0.5
      mask_cen(i,j,k,CLR,itr,ipr) =  mask_cen(i,j,k,CLR,itr,ipr) + (1. - cf3d(i,j,k)) * 0.5
    end if

    !Next, look at the top boundary and determine it's
    !contribution to the cell center transport class.
    if( sum(maskup(k+1,:)) > 0 ) then
      itr = UP1 + maxloc(maskup(k+1,:),1)-1
    else if( sum(maskdn(k+1,:)) > 0 ) then
      itr = DN1 + maxloc(maskdn(k+1,:),1)-1
    else if( maskqu(k+1) > 0 ) then
      itr = QUI
    else
      call endrun("ERROR: setup_class_masks: We should not be in this place for cell tops.")
    end if

    ! +++mhwang
    ! In the clear, and non-precipitating class, it is classified as quiescent class in the MMF simulation.
    ! If this is classed as updraft or downdraft in mode 16, this would lead to too much upraft and downdraft mass fluxes.
    ! Minghuai Wang, 2010-01-18 (Minghuai.Wang@pnl.gov)
    !           if(icl.eq.CLR .and. ipr.eq.PRN) then
    !             itr = QUI
    !           end if
    !---mhwang

    !We have what we need for the cell top classes so increment
    !the center mask for the top half...
    !           mask_cen(i,j,k,icl,itr,ipr) = mask_cen(i,j,k,icl,itr,ipr) + 0.5
    ! use fractional cloudiness in SAM
    if(icl.eq.CLR) then
      mask_cen(i,j,k,icl,itr,ipr) = mask_cen(i,j,k,icl,itr,ipr) + 0.5
    else if(icl.eq.CLD) then
      mask_cen(i,j,k,CLD,itr,ipr) =  mask_cen(i,j,k,CLD,itr,ipr) + (cf3d(i,j,k))*0.5
      mask_cen(i,j,k,CLR,itr,ipr) =  mask_cen(i,j,k,CLR,itr,ipr) + (1. - cf3d(i,j,k)) * 0.5
    end if

  end do !k-loop mask for centers

end do
end do XYLOOPS
end subroutine setup_class_masks


!------------------------------------------------------------------------
subroutine cloud_prcp_check(mask1, flag1, mask2, flag2, k, iout, msg)
  !
  ! Assigns the flag associated with the mask value that is true to the
  ! output index. The masks are assumed to be 1-D arrays and k is the
  ! position in the array to check.
  ! William.Gustafson@pnl.gov; 11-Sep-2008
  !------------------------------------------------------------------------
  !
  ! Soubroutine arguments...
  !
  integer, dimension(:), intent(in) :: mask1, mask2
  integer, intent(in) :: flag1, flag2, k
  integer, intent(out) :: iout
  character(len=*), optional :: msg
  !
  ! Local var...
  !
  integer :: n
  !
  ! Sanity check
  !
  n = ubound(mask1,1)
  if( k < 1 .or. k > n) then
    write(0, *) 'cloud_prcp_check', 'k =',k, ' n =',n
    call endrun('ERROR: k out of bounds in cloud_prcp_check')
  end if
  !
  ! Whichever mask has the value 1 has the associated flag put into iout
  !
  if( mask1(k) > 0 .and. mask2(k) < 1 ) then
    iout = flag1
  else if( mask2(k) > 0 .and. mask1(k) < 1) then
    iout = flag2
  else
    write(0, *) 'cloud_prcp_check', 'k =', k
    call endrun("ERROR: neither mask dominates in cloud_prcp_check")
  end if

end subroutine cloud_prcp_check

#endif /*ECPP*/
end module module_ecpp_stats
