module precip_proc_mod
  use params, only: asyncid
  implicit none

contains

  subroutine precip_proc(ncrms,qpsrc,qpevp,q,qp,qn)

    use vars
    use micro_params
    use params
    use sat_mod

    implicit none
    integer, intent(in) :: ncrms
    real(crm_rknd) :: q(ncrms,dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    real(crm_rknd) :: qp(ncrms,dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    real(crm_rknd) qn(ncrms,nx,ny,nzm)  ! cloud condensate (liquid + ice)
    real(crm_rknd) qpsrc(ncrms,nz)  ! source of precipitation microphysical processes
    real(crm_rknd) qpevp(ncrms,nz)  ! sink of precipitating water due to evaporation

    integer i,j,k,icrm
    real(crm_rknd) autor, autos, accrr, accris, accrcs, accrig, accrcg
    real(crm_rknd) dq, omn, omp, omg, qsatt
    real(crm_rknd) pows1, pows2, powg1, powg2, powr1, powr2, tmp
    real(crm_rknd) qii, qcc, qrr, qss, qgg

    powr1 = (3 + b_rain) / 4.
    powr2 = (5 + b_rain) / 8.
    pows1 = (3 + b_snow) / 4.
    pows2 = (5 + b_snow) / 8.
    powg1 = (3 + b_grau) / 4.
    powg2 = (5 + b_grau) / 8.

    !$acc parallel loop collapse(2) async(asyncid)
    do k=1,nzm
      do icrm = 1 , ncrms
        qpsrc(icrm,k)=0.
        qpevp(icrm,k)=0.
      enddo
    enddo

    !$acc parallel loop collapse(4) async(asyncid)
    do k=1,nzm
      do j=1,ny
        do i=1,nx
          do icrm = 1 , ncrms

            !-------     Autoconversion/accretion

            if(qn(icrm,i,j,k)+qp(icrm,i,j,k).gt.0.) then


              omn = max(real(0.,crm_rknd),min(real(1.,crm_rknd),(tabs(icrm,i,j,k)-tbgmin)*a_bg))
              omp = max(real(0.,crm_rknd),min(real(1.,crm_rknd),(tabs(icrm,i,j,k)-tprmin)*a_pr))
              omg = max(real(0.,crm_rknd),min(real(1.,crm_rknd),(tabs(icrm,i,j,k)-tgrmin)*a_gr))

              if(qn(icrm,i,j,k).gt.0.) then

                qcc = qn(icrm,i,j,k) * omn
                qii = qn(icrm,i,j,k) * (1.-omn)

                if(qcc .gt. qcw0) then
                  autor = alphaelq
                else
                  autor = 0.
                endif

                if(qii .gt. qci0) then
                  autos = betaelq*coefice(icrm,k)
                else
                  autos = 0.
                endif

                accrr = 0.
                if(omp.gt.0.001) then
                  qrr = qp(icrm,i,j,k) * omp
                  accrr = accrrc(icrm,k) * qrr ** powr1
                endif
                accrcs = 0.
                accris = 0.
                if(omp.lt.0.999.and.omg.lt.0.999) then
                  qss = qp(icrm,i,j,k) * (1.-omp)*(1.-omg)
                  tmp = qss ** pows1
                  accrcs = accrsc(icrm,k) * tmp
                  accris = accrsi(icrm,k) * tmp
                endif
                accrcg = 0.
                accrig = 0.
                if(omp.lt.0.999.and.omg.gt.0.001) then
                  qgg = qp(icrm,i,j,k) * (1.-omp)*omg
                  tmp = qgg ** powg1
                  accrcg = accrgc(icrm,k) * tmp
                  accrig = accrgi(icrm,k) * tmp
                endif
                qcc = (qcc+dtn*autor*qcw0)/(1.+dtn*(accrr+accrcs+accrcg+autor))
                qii = (qii+dtn*autos*qci0)/(1.+dtn*(accris+accrig+autos))
                dq = dtn *(accrr*qcc + autor*(qcc-qcw0)+(accris+accrig)*qii + (accrcs+accrcg)*qcc + autos*(qii-qci0))
                dq = min(dq,qn(icrm,i,j,k))
                qp(icrm,i,j,k) = qp(icrm,i,j,k) + dq
                q(icrm,i,j,k) = q(icrm,i,j,k) - dq
                qn(icrm,i,j,k) = qn(icrm,i,j,k) - dq
                !$acc atomic update
                qpsrc(icrm,k) = qpsrc(icrm,k) + dq

              elseif(qp(icrm,i,j,k).gt.qp_threshold.and.qn(icrm,i,j,k).eq.0.) then

                qsatt = 0.
                if(omn.gt.0.001) qsatt = qsatt + omn*qsatw_crm(tabs(icrm,i,j,k),pres(icrm,k))
                if(omn.lt.0.999) qsatt = qsatt + (1.-omn)*qsati_crm(tabs(icrm,i,j,k),pres(icrm,k))
                dq = 0.
                if(omp.gt.0.001) then
                  qrr = qp(icrm,i,j,k) * omp
                  dq = dq + evapr1(icrm,k)*sqrt(qrr) + evapr2(icrm,k)*qrr**powr2
                endif
                if(omp.lt.0.999.and.omg.lt.0.999) then
                  qss = qp(icrm,i,j,k) * (1.-omp)*(1.-omg)
                  dq = dq + evaps1(icrm,k)*sqrt(qss) + evaps2(icrm,k)*qss**pows2
                endif
                if(omp.lt.0.999.and.omg.gt.0.001) then
                  qgg = qp(icrm,i,j,k) * (1.-omp)*omg
                  dq = dq + evapg1(icrm,k)*sqrt(qgg) + evapg2(icrm,k)*qgg**powg2
                endif
                dq = dq * dtn * (q(icrm,i,j,k) /qsatt-1.)
                dq = max(-0.5*qp(icrm,i,j,k),dq)
                qp(icrm,i,j,k) = qp(icrm,i,j,k) + dq
                q(icrm,i,j,k) = q(icrm,i,j,k) - dq
                !$acc atomic update
                qpevp(icrm,k) = qpevp(icrm,k) + dq

              else

                q(icrm,i,j,k) = q(icrm,i,j,k) + qp(icrm,i,j,k)
                !$acc atomic update
                qpevp(icrm,k) = qpevp(icrm,k) - qp(icrm,i,j,k)
                qp(icrm,i,j,k) = 0.

              endif

            endif

            dq = qp(icrm,i,j,k)
            qp(icrm,i,j,k)=max(real(0.,crm_rknd),qp(icrm,i,j,k))
            q(icrm,i,j,k) = q(icrm,i,j,k) + (dq-qp(icrm,i,j,k))

          enddo
        enddo
      enddo
    enddo

  end subroutine precip_proc

end module precip_proc_mod
