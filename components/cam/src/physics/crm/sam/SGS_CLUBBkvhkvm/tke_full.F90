
subroutine tke_full

!	this subroutine solves the TKE equation

use vars
use sgs
use params
implicit none

real def2(nx,ny,nzm)	
real grd,betdz,Ck,Ce,Ces,Ce1,Ce2,smix,Pr,Cee,Cs
real buoy_sgs,ratio,a_prod_sh,a_prod_bu,a_diss
real lstarn, lstarp, bbb, omn, omp
real qsatt,dqsat
integer i,j,k,kc,kb

!call t_startf('tke_full')

!Cs = 0.1944
Cs = 0.15
Ck=0.1
Ce=Ck**3/Cs**4
Ces=Ce/0.7*3.0	

if(RUN3D) then
  call shear_prod3D(def2)
else
  call shear_prod2D(def2)
endif

do k=1,nzm      
  kb=k-1
  kc=k+1

  grd=dz*adz(k)

  betdz=bet(k)/dz/(adzw(kc)+adzw(k))
  Ce1=Ce/0.7*0.19
  Ce2=Ce/0.7*0.51
  if(k.eq.1) then
    kb=1
    kc=2
    betdz=bet(k)/dz/adzw(kc)
    Ce1=Ces/0.7*0.19
    Ce2=Ces/0.7*0.51
  end if
  if(k.eq.nzm) then
    kb=nzm-1
    kc=nzm
    betdz=bet(k)/dz/adzw(k)
    Ce1=Ces/0.7*0.19
    Ce2=Ces/0.7*0.51
  end if
  tkelediss(k) = 0.
  tkesbdiss(k) = 0.
  tkesbshear(k)= 0.
  tkesbbuoy(k) = 0.
  do j=1,ny
  do i=1,nx
!  SGS buoyancy flux

!bloss: removed temperature diagnostics for omn.
!         - use mass weighted qsat, dqsat and latent heat for cloud 
!         - separate buoyancy contributions for precipitating water and ice.
   
   
     if(qcl(i,j,k)+qci(i,j,k) .gt. 0.) then
      
       omn = qcl(i,j,k)/(qcl(i,j,k)+qci(i,j,k)+1.e-20)
       lstarn = fac_cond+(1.-omn)*fac_fus
      
       dqsat = omn*dtqsatw_crm(tabs(i,j,k),pres(k))+ &
                             (1.-omn)*dtqsati_crm(tabs(i,j,k),pres(k))
       qsatt = omn*qsatw_crm(tabs(i,j,k),pres(k))+(1.-omn)*qsati_crm(tabs(i,j,k),pres(k))
       bbb = 1. + epsv*qsatt-qcl(i,j,k)-qci(i,j,k) -qpl(i,j,k)-qpi(i,j,k)+1.61*tabs(i,j,k)*dqsat
       bbb = bbb / (1.+lstarn*dqsat)
       buoy_sgs=betdz*(bbb*(t(i,j,kc)-t(i,j,kb)) &
         +(bbb*lstarn - (1.+lstarn*dqsat)*tabs(i,j,k))* &
             (qv(i,j,kc)+qcl(i,j,kc)+qci(i,j,kc)-qv(i,j,kb)-qcl(i,j,kb)-qci(i,j,kb)) & 
         + (bbb*fac_cond - (1.+fac_cond*dqsat)*tabs(i,j,k))*(qpl(i,j,kc)-qpl(i,j,kb))  &
         + (bbb*fac_sub  - (1.+fac_sub *dqsat)*tabs(i,j,k))*(qpi(i,j,kc)-qpi(i,j,kb)) )
!bloss   +(bbb*lstarp - (1.+lstarp*dqsat)*tabs(i,j,k))* &
!bloss            (qpl(i,j,kc)+qpi(i,j,kc)-qpl(i,j,kb)-qpi(i,j,kb)) )
     else
      
        bbb = 1.+epsv*qv(i,j,k)-qpl(i,j,k)-qpi(i,j,k)
        buoy_sgs=betdz*( bbb*(t(i,j,kc)-t(i,j,kb)) &
          +epsv*tabs(i,j,k)* &
           (qv(i,j,kc)+qcl(i,j,kc)+qci(i,j,kc)-qv(i,j,kb)-qcl(i,j,kb)-qci(i,j,kb)) &
          +(bbb*fac_cond-tabs(i,j,k))*(qpl(i,j,kc)-qpl(i,j,kb)) &
          +(bbb*fac_sub -tabs(i,j,k))*(qpi(i,j,kc)-qpi(i,j,kb)) )
!bloss    +(bbb*lstarp-tabs(i,j,k))* &
!bloss         (qpl(i,j,kc)+qpi(i,j,kc)-qpl(i,j,kb)-qpi(i,j,kb)) )
      end if

   if(buoy_sgs.le.0.) then
     smix=grd
   else
     smix=min(grd,max(0.1*grd, sqrt(0.76*tk(i,j,k)/Ck/sqrt(buoy_sgs+1.e-10))))
   end if


   ratio=smix/grd
   Pr=1. 
!   Pr=1. +2.*ratio
   Cee=Ce1+Ce2*ratio

   if(dosmagor) then

     tk(i,j,k)=sqrt(Ck**3/Cee*max(0.,def2(i,j,k)-Pr*buoy_sgs))*smix**2
     tke(i,j,k) = (tk(i,j,k)/(Ck*smix))**2
     a_prod_sh=(tk(i,j,k)+0.001)*def2(i,j,k)
     a_prod_bu=-(tk(i,j,k)+0.001)*Pr*buoy_sgs
     a_diss=a_prod_sh+a_prod_bu

   else

     tke(i,j,k)=max(0.,tke(i,j,k))
     a_prod_sh=(tk(i,j,k)+0.001)*def2(i,j,k)
     a_prod_bu=-(tk(i,j,k)+0.001)*Pr*buoy_sgs
     a_diss=min(tke(i,j,k)/(4.*dt),Cee/smix*tke(i,j,k)**1.5) ! cap the diss rate (useful for large time steps
     tke(i,j,k)=max(0.,tke(i,j,k)+dtn*(max(0.,a_prod_sh+a_prod_bu)-a_diss))
     tk(i,j,k)=Ck*smix*sqrt(tke(i,j,k))

   end if
	
   tkh(i,j,k)=Pr*tk(i,j,k)

   tkelediss(k) = tkelediss(k) - a_prod_sh
   tkesbdiss(k) = tkesbdiss(k) + a_diss
   tkesbshear(k)= tkesbshear(k)+ a_prod_sh
   tkesbbuoy(k) = tkesbbuoy(k) + a_prod_bu

  end do ! i
  end do ! j

  tkelediss(k) = tkelediss(k)/float(nx*ny)


end do ! k

!call t_stopf('tke_full')

end


