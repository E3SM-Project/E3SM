module cloud_mod
  use params, only: asyncid
  implicit none

contains

  subroutine cloud(ncrms,q,qp,qn)

    !  Condensation of cloud water/cloud ice.

    use vars
    use micro_params
    use params
    use sat_mod

    implicit none
    integer, intent(in) :: ncrms
    real(crm_rknd) :: q(ncrms,dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    real(crm_rknd) :: qp(ncrms,dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    real(crm_rknd) qn(ncrms,nx,ny,nzm)  ! cloud condensate (liquid + ice)

    integer i,j,k, kb, kc,icrm
    real(crm_rknd) dtabs, tabs1, an, bn, ap, bp, om, ag, omp
    real(crm_rknd) fac1,fac2
    real(crm_rknd) fff,dfff,qsatt,dqsat
    real(crm_rknd) lstarn,dlstarn,lstarp,dlstarp
    integer niter

    an = 1./(tbgmax-tbgmin)
    bn = tbgmin * an
    ap = 1./(tprmax-tprmin)
    bp = tprmin * ap
    fac1 = fac_cond+(1+bp)*fac_fus
    fac2 = fac_fus*ap
    ag = 1./(tgrmax-tgrmin)

    !$acc parallel loop collapse(4) async(asyncid)
    do k = 1, nzm
      do j = 1, ny
        do i = 1, nx
          do icrm = 1 , ncrms
            q(icrm,i,j,k)=max(real(0.,crm_rknd),q(icrm,i,j,k))
            ! Initial guess for temperature assuming no cloud water/ice:
            tabs(icrm,i,j,k) = t(icrm,i,j,k)-gamaz(icrm,k)
            tabs1=(tabs(icrm,i,j,k)+fac1*qp(icrm,i,j,k))/(1.+fac2*qp(icrm,i,j,k))
            ! Warm cloud:
            if(tabs1.ge.tbgmax) then
              tabs1=tabs(icrm,i,j,k)+fac_cond*qp(icrm,i,j,k)
              qsatt = qsatw_crm(tabs1,pres(icrm,k))
              ! Ice cloud:
            elseif(tabs1.le.tbgmin) then
              tabs1=tabs(icrm,i,j,k)+fac_sub*qp(icrm,i,j,k)
              qsatt = qsati_crm(tabs1,pres(icrm,k))
              ! Mixed-phase cloud:
            else
              om = an*tabs1-bn
              qsatt = om*qsatw_crm(tabs1,pres(icrm,k))+(1.-om)*qsati_crm(tabs1,pres(icrm,k))
            endif
            !  Test if condensation is possible:
            if(q(icrm,i,j,k).gt.qsatt) then
              niter=0
              dtabs = 100.
              do while(abs(dtabs).gt.0.01.and.niter.lt.10)
                if(tabs1.ge.tbgmax) then
                  om=1.
                  lstarn=fac_cond
                  dlstarn=0.
                  qsatt=qsatw_crm(tabs1,pres(icrm,k))
                  dqsat=dtqsatw_crm(tabs1,pres(icrm,k))
                else if(tabs1.le.tbgmin) then
                  om=0.
                  lstarn=fac_sub
                  dlstarn=0.
                  qsatt=qsati_crm(tabs1,pres(icrm,k))
                  dqsat=dtqsati_crm(tabs1,pres(icrm,k))
                else
                  om=an*tabs1-bn
                  lstarn=fac_cond+(1.-om)*fac_fus
                  dlstarn=an*fac_fus
                  qsatt=om*qsatw_crm(tabs1,pres(icrm,k))+(1.-om)*qsati_crm(tabs1,pres(icrm,k))
                  dqsat=om*dtqsatw_crm(tabs1,pres(icrm,k))+(1.-om)*dtqsati_crm(tabs1,pres(icrm,k))
                endif
                if(tabs1.ge.tprmax) then
                  omp=1.
                  lstarp=fac_cond
                  dlstarp=0.
                else if(tabs1.le.tprmin) then
                  omp=0.
                  lstarp=fac_sub
                  dlstarp=0.
                else
                  omp=ap*tabs1-bp
                  lstarp=fac_cond+(1.-omp)*fac_fus
                  dlstarp=ap*fac_fus
                endif
                fff = tabs(icrm,i,j,k)-tabs1+lstarn*(q(icrm,i,j,k)-qsatt)+lstarp*qp(icrm,i,j,k)
                dfff=dlstarn*(q(icrm,i,j,k)-qsatt)+dlstarp*qp(icrm,i,j,k)-lstarn*dqsat-1.
                dtabs=-fff/dfff
                niter=niter+1
                tabs1=tabs1+dtabs
              end do
              qsatt = qsatt + dqsat * dtabs
              qn(icrm,i,j,k) = max(real(0.,crm_rknd),q(icrm,i,j,k)-qsatt)
            else
              qn(icrm,i,j,k) = 0.
            endif
            tabs(icrm,i,j,k) = tabs1
            qp(icrm,i,j,k) = max(real(0.,crm_rknd),qp(icrm,i,j,k)) ! just in case
          end do
        end do
      end do
    end do

  end subroutine cloud

end module cloud_mod
