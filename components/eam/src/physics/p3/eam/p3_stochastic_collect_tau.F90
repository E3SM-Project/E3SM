module p3_stochastic_collect_tau
! Gettelman modified 2024 for P3 in E3SM
! From Morrison (Lebo, originally TAU bin code)
! Gettelman and Chen 2018
!the subroutines take in air density, air temperature, and the bin mass boundaries, and 
!output the mass and number mixing ratio tendencies in each bin directly.
!this is then wrapped for CAM. 

! note, this is now coded locally. Want the CAM interface to be over i,k I think.


use shr_kind_mod,      only: r8=>shr_kind_r8
use cam_history,       only: addfld
use physconst,         only: pi, rhoh2o
use micro_p3_utils,    only: qsmall 
use cam_logfile,       only: iulog
! if using ACC directives for GPUs, will need VLENS Here too

implicit none
private
save

! Subroutines
public :: p3_stochastic_kernel_init, p3_stochastic_collect_tau_tend

!In the module top, declare the following so that these can be used throughout the module:

integer, parameter, public  :: ncd = 35
integer, parameter, public  :: ncdp = ncd + 1
integer, parameter, public  :: ncdl = ncd
integer, parameter, public  :: ncdpl = ncdl+1

! for Zach's collision-coalescence code

real(r8), private :: knn(ncd,ncd)

real(r8), private :: mmean(ncd), diammean(ncd)       ! kg & m at bin mid-points
real(r8), private :: medge(ncdp), diamedge(ncdp)     ! kg & m at bin edges 
integer, private  :: cutoff_id                       ! cutoff between cloud water and rain drop, D = 40 microns

! Assume 6 microns for drop radius (if no number)...
real(r8), parameter :: m1 = 4._r8/3._r8*pi*rhoh2o*(6.e-6_r8)**3
! And 200 microns for rain
real(r8), parameter :: r1 = 4._r8/3._r8*pi*rhoh2o*(200.e-6_r8)**3

!$acc declare create(knn,cutoff_id,mmean,diammean,medge,diamedge)

!===============================================================================
contains
!===============================================================================

subroutine calc_bins    

  real(r8) :: DIAM(ncdp)
  real(r8) :: X(ncdp)
  real(r8) :: radsl(ncdp)
  real(r8) :: radl(ncd)
  integer  :: L, lcl  
  real(r8) :: kkfac

!Then before doing any calculations you'll need to calculate the bin mass grid 
! (note this code could be cleaned up, I'm just taking it as it's used in our bin scheme). 
! This only needs to be done once, since we'll use the same bin mass grid for all calculations. 

! use mass doubling bins from Graham Feingold (note cgs units)

  DIAM(1)=1.5625*2.E-04_r8                ! cm
  X(1)=PI/6._r8*DIAM(1)**3*rhoh2o/1000._r8  ! rhoh2o kg/m3 --> g/cm3 
  radsl(1) = X(1)                         ! grams 

  DO l=2,ncdp
     X(l)=2._r8*X(l-1)
     DIAM(l)=(6._r8/pi*X(l)*1000._r8/rhoh2o)**(1._r8/3._r8)  ! cm
     radsl(l)=X(l)             
  ENDDO

! now get bin mid-points

  do l=1,ncd
     radl(l)=(radsl(l)+radsl(l+1))/2._r8         ! grams   
     diammean(l) = (6._r8/pi*radl(l)*1000._r8/rhoh2o)**(1._r8/3._r8) ! cm
  end do

! set bin grid for method of moments

  ! for method of moments

  do lcl = 1,ncd+1
     medge(lcl) = radsl(lcl)               ! grams
     diamedge(lcl) = DIAM(lcl)             ! cm
  enddo

  do lcl = 1,ncd
     mmean(lcl) = radl(lcl)  
     diammean(lcl) = diammean(lcl)
  enddo

  do lcl = ncdp,1,-1
     if( diamedge(lcl).ge.40.e-4_r8 ) cutoff_id = lcl
  end do  

end subroutine calc_bins

subroutine p3_stochastic_kernel_init(lookup_file_dir,kernel_filename)

  character(len=*), intent(in) :: lookup_file_dir  ! Directory of P3 lookup tables (eventual home of sc kernel)

  character(len=*), intent(in) :: kernel_filename

  integer :: iunit ! unit number of opened file for collection kernel code from a lookup table.

  integer :: idd, jdd
  real(r8) :: kkfac

  ! Also assumed to be copied into run_directory. 
  ! Eventually/Ideally (1) use lookup_file_dir + kernel_filename and/or (2) Make namelist

  call calc_bins

! Read in the collection kernel code from a lookup table. Again, this only needs to be done once.
! use kernel from Zach (who got it from Jerry)

  KNN(:,:)=0._r8 ! initialize values
  kkfac=1.5_r8   ! from Zach

  open(newunit=iunit,file=kernel_filename,status='old')

  do idd=1,ncd
     do jdd=1,idd
        READ(iunit,941) KNN(IDD,JDD)
941     FORMAT(2X,E12.5)

        KNN(IDD,JDD)=(mmean(IDD)*kkfac+mmean(JDD)*kkfac)*KNN(IDD,JDD)
        if (knn(idd,jdd) < 0._r8) knn(idd,jdd)=0._r8
     end do
  end do

  !$acc update device(knn,cutoff_id,mmean,diammean,medge,diamedge)

end subroutine p3_stochastic_kernel_init

!main driver routine
!needs to pull in i,k fields (so might need dimensions here too)

subroutine p3_stochastic_collect_tau_tend(deltatin,t,rho,qcin,     &
                        ncin,qrin,nrin,mu_c,lambda_c,     &
                        n0r,lambda_r,qcin_new,ncin_new,qrin_new,nrin_new,   &
                        qctend_TAU,nctend_TAU, qrtend_TAU,nrtend_TAU,       &
                        scale_qc,scale_nc,scale_qr,   &
                        scale_nr,amk_c,ank_c,amk_r,ank_r,amk,ank,amk_out,   &
                        ank_out,gmnnn_lmnnn_TAU)

  !inputs: t,rho,qcin,ncin,qrin,nrin
  !outputs: qctend,nctend,qrtend,nrtend
  !not sure if we want to output bins (extra dimension). Good for testing?  
  
  real(r8), intent(in) :: deltatin
  real(r8), intent(in) :: t
  real(r8), intent(in) :: rho
  real(r8), intent(in) :: qcin
  real(r8), intent(in) :: ncin
  real(r8), intent(in) :: qrin
  real(r8), intent(in) :: nrin

  real(r8), intent(in) :: mu_c
  real(r8), intent(in) :: lambda_c
  real(r8), intent(in) :: lambda_r
  real(r8), intent(in) :: n0r

  real(r8), intent(out) :: qcin_new
  real(r8), intent(out) :: ncin_new
  real(r8), intent(out) :: qrin_new
  real(r8), intent(out) :: nrin_new

  real(r8), intent(out) :: qctend_TAU
  real(r8), intent(out) :: nctend_TAU
  real(r8), intent(out) :: qrtend_TAU
  real(r8), intent(out) :: nrtend_TAU
  
  real(r8), intent(out) :: scale_qc
  real(r8), intent(out) :: scale_nc
  real(r8), intent(out) :: scale_qr
  real(r8), intent(out) :: scale_nr
  
  real(r8), intent(out) :: amk_c(ncd)
  real(r8), intent(out) :: ank_c(ncd)
  real(r8), intent(out) :: amk_r(ncd)
  real(r8), intent(out) :: ank_r(ncd)
  real(r8), intent(out) :: amk(ncd)
  real(r8), intent(out) :: ank(ncd)
  real(r8), intent(out) :: amk_out(ncd)
  real(r8), intent(out) :: ank_out(ncd)
  
  real(r8), intent(out) :: gmnnn_lmnnn_TAU

  ! Local variables
  
  integer :: i,k,n,lcl
  integer :: cutoff_amk,cutoff
  
  real(r8) :: all_gmnnn,all_lmnnn
  real(r8) :: qscl
  
  real(r8) :: qcin_old
  real(r8) :: ncin_old
  real(r8) :: qrin_old
  real(r8) :: nrin_old
  
  real(r8) :: amk0(ncd)
  real(r8) :: ank0(ncd)
  real(r8) :: gnnnn(ncd)
  real(r8) :: gmnnn(ncd)
  real(r8) :: lnnnn(ncd)
  real(r8) :: lmnnn(ncd)
  real(r8) :: gnnnn0(ncd)
  real(r8) :: gmnnn0(ncd)
  real(r8) :: lnnnn0(ncd)
  real(r8) :: lmnnn0(ncd)
  
  integer, parameter :: sub_step = 60

  !$acc data create (cutoff_amk,cutoff,qcin_old,ncin_old,qrin_old, &
  !$acc              nrin_old,amk0,ank0,gnnnn,gmnnn,lnnnn,lmnnn,   &
  !$acc              gnnnn0,gmnnn0,lnnnn0,lmnnn0)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)  
  !do k=1,nlev
  !   do i=1,mgncol
       cutoff = cutoff_id - 1
  !   end do
  !end do
  !$acc end parallel

  ! First make bins from cam size distribution (bins are diagnostic)

  call cam_bin_distribute(qcin,ncin,qrin,nrin,mu_c,lambda_c, &
                          lambda_r,n0r,scale_qc, &
                          scale_nc,scale_qr,scale_nr,amk_c,ank_c,  &
                          amk_r,ank_r,amk,ank,cutoff_amk)


   if ( cutoff_amk > 0 ) then
      cutoff = cutoff_amk
   end if

!Then call the subroutines that actually do the calculations. The inputs/ouputs are described in comments below. 

!This part of the code needs to be called each time for each process rate calculation 
! (i.e., for each sampled cloud/rain gamma distribution):

! note: variables passed to compute_column_params are all inputs,
! outputs from this subroutine are stored as global variables

! inputs: t --> input air temperature (K)
!         rho --> input air density (kg/m^3)
!         medge --> bin mass boundary (g) 
!         amk --> array of bin mass mixing ratio, i.e., the input drop mass distribution (kg/kg)
!         ank --> array of bin number mixing ratio, i.e., the input drop number distribution (kg^-1)

! inputs: medge --> bin mass boundary (g), same as above

! outputs: gnnnn --> bin number mass mixing tendency gain, array in bins (#/cm^3/s)
!          lnnnn --> bin number mass mixing tendency loss, array in bins (#/cm^3/s)
!          gmnnn --> bin mass mixing ratio tendency gain, array in bins (g/cm^3/s) 
!          lmnnn --> bin mass mixing ratio tendency loss, array in bins (g/cm^3/s)


! Call Kernel

   qcin_new = 0._r8
   ncin_new = 0._r8
   qrin_new = 0._r8
   nrin_new = 0._r8
        
   qcin_old = 0._r8
   ncin_old = 0._r8
   qrin_old = 0._r8
   nrin_old = 0._r8
        
   qctend_TAU = 0._r8
   nctend_TAU = 0._r8
   qrtend_TAU = 0._r8
   nrtend_TAU = 0._r8


  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(1)
  do lcl=1,ncd
      gnnnn(lcl) = 0._r8
      gmnnn(lcl) = 0._r8
      lnnnn(lcl) = 0._r8
      lmnnn(lcl) = 0._r8
  end do
  !$acc end parallel

! update qc, nc, qr, nr

   !$acc loop seq
    do lcl=1,ncd
      amk0(lcl) = amk(lcl)
      ank0(lcl) = ank(lcl)
   end do
   ! substep bin code
   !$acc loop seq
   do n=1,sub_step
      call compute_coll_params(rho,medge,amk0(1:ncd),ank0(1:ncd),gnnnn0(1:ncd),gmnnn0(1:ncd),lnnnn0(1:ncd),lmnnn0(1:ncd))

      all_gmnnn=0._r8
      all_lmnnn=0._r8
      !scaling gmnnn, lmnnn
      !$acc loop seq
      do lcl=1,ncd
         all_gmnnn = all_gmnnn+gmnnn0(lcl)
         all_lmnnn = all_lmnnn+lmnnn0(lcl)
      end do
 
      if ( (all_gmnnn == 0._r8) .or. (all_lmnnn == 0._r8) ) then
         !$acc loop seq
         do lcl=1,ncd
            gmnnn0(lcl) = 0._r8
            lmnnn0(lcl) = 0._r8
         end do
      else
         !$acc loop seq
         do lcl=1,ncd
            lmnnn0(lcl) = lmnnn0(lcl)*(all_gmnnn/all_lmnnn)
         end do
      end if

      !$acc loop seq
      do lcl=1,ncd
         amk0(lcl) = amk0(lcl)+(gmnnn0(lcl)-lmnnn0(lcl))*1.e3_r8/ &
                     rho*deltatin/dble(sub_step)
         ank0(lcl) = ank0(lcl)+(gnnnn0(lcl)-lnnnn0(lcl))*1.e6_r8/ &
                     rho*deltatin/dble(sub_step)
         gmnnn(lcl) = gmnnn(lcl)+gmnnn0(lcl)/sub_step
         gnnnn(lcl) = gnnnn(lcl)+gnnnn0(lcl)/sub_step
         lmnnn(lcl) = lmnnn(lcl)+lmnnn0(lcl)/sub_step
         lnnnn(lcl) = lnnnn(lcl)+lnnnn0(lcl)/sub_step
      end do
   end do ! end of loop "sub_step"

   ! cloud water
   !$acc loop seq
   do lcl=1,cutoff
      qcin_old = qcin_old+amk(lcl)
      ncin_old = ncin_old+ank(lcl)
      qcin_new = qcin_new+(gmnnn(lcl)-lmnnn(lcl))*1.e3_r8/rho*deltatin
      ncin_new = ncin_new+(gnnnn(lcl)-lnnnn(lcl))*1.e6_r8/rho*deltatin
      qctend_TAU = qctend_TAU+(amk0(lcl)-amk(lcl))/deltatin
      nctend_TAU = nctend_TAU+(ank0(lcl)-ank(lcl))/deltatin
      gmnnn_lmnnn_TAU = gmnnn_lmnnn_TAU+gmnnn(lcl)-lmnnn(lcl)
   end do

   ! rain
   !$acc loop seq
   do lcl=cutoff+1,ncd
      qrin_old = qrin_old+amk(lcl)
      nrin_old = nrin_old+ank(lcl)
      qrin_new = qrin_new+(gmnnn(lcl)-lmnnn(lcl))*1.e3_r8/rho*deltatin
      nrin_new = nrin_new+(gnnnn(lcl)-lnnnn(lcl))*1.e6_r8/rho*deltatin
      qrtend_TAU = qrtend_TAU+(amk0(lcl)-amk(lcl))/deltatin
      nrtend_TAU = nrtend_TAU+(ank0(lcl)-ank(lcl))/deltatin
      gmnnn_lmnnn_TAU = gmnnn_lmnnn_TAU+gmnnn(lcl)-lmnnn(lcl)
   end do

   !$acc loop seq     
   do lcl=1,ncd
      amk_out(lcl) = amk(lcl) + (gmnnn(lcl)-lmnnn(lcl))*1.e3_r8/rho*deltatin
      ank_out(lcl) = ank(lcl) + (gnnnn(lcl)-lnnnn(lcl))*1.e6_r8/rho*deltatin
   end do
     
   qcin_new = qcin_new+qcin_old
   ncin_new = ncin_new+ncin_old
   qrin_new = qrin_new+qrin_old
   nrin_new = nrin_new+nrin_old
  

! Conservation checks 
! AG: Added May 2023

   ! First make sure all not negative
   qcin_new=max(qcin_new,0._r8)
   ncin_new=max(ncin_new,0._r8)
   qrin_new=max(qrin_new,0._r8)
   nrin_new=max(nrin_new,0._r8)

! Now adjust so that sign is correct. qc_new,nc_new <= input, qr_new >= input
! NOte that due to self collection nr can be larger or smaller than input....
! Makes above check redundant I think.

   qcin_new=min(qcin_new,qcin)
   ncin_new=min(ncin_new,ncin)
   qrin_new=max(qrin_new,qrin)

! Next scale mass...so output qc+qr is the same as input

   if ( (qcin_new+qrin_new) > 0._r8 ) then
      qscl = (qcin+qrin)/(qcin_new+qrin_new)
   else
      qscl = 0._r8
   end if
   qcin_new = qcin_new * qscl
   qrin_new = qrin_new * qscl

! Now zero nr,nc if either small or no mass?

   if ( qcin_new < qsmall ) then
      ncin_new = 0._r8
   end if
        
   if ( qrin_new < qsmall ) then
      nrin_new = 0._r8
   end if

!Finally add number if mass but no (or super small) number

   if ( qcin_new > qsmall .and. ncin_new < qsmall ) then
      ncin_new = qcin_new/m1
   end if
        
   if ( qrin_new > qsmall .and. nrin_new < qsmall) then
      nrin_new = qrin_new/r1
   end if

! Set a maximum on new number in case the adjustment above fails due to edge cases

   nrin_new=min(nrin_new,1.e9_r8)

! Then recalculate tendencies based on difference
! Clip tendencies for cloud (qc,nc) to be <= 0. 
! Qrtend is not used in pumas (-qctend is used) but clip that too). 
! Nr tend can be muliply signed. 

   qctend_TAU= min((qcin_new - qcin) / deltatin,0._r8)
   nctend_TAU= min((ncin_new - ncin) / deltatin,0._r8)
   qrtend_TAU= max((qrin_new - qrin) / deltatin,0._r8)
   nrtend_TAU= (nrin_new - nrin) / deltatin

  !$acc end data

end subroutine p3_stochastic_collect_tau_tend

subroutine cam_bin_distribute(qc,nc,qr,nr,mu_c,lambda_c, &
                              lambda_r,n0r,scale_qc, &
                              scale_nc,scale_qr,scale_nr,amk_c,ank_c,  &
                              amk_r,ank_r,amk,ank,cutoff_amk)

  implicit none

  real(r8), intent(in) :: qc,nc,qr,nr,mu_c, &
                          lambda_c,lambda_r,n0r
  real(r8), dimension(ncd), intent(out) :: amk_c,ank_c,amk_r,ank_r,amk,ank
  real(r8), intent(out) :: scale_nc,scale_qc,scale_nr,scale_qr 
  integer, intent(out) :: cutoff_amk

  ! Local variables

  integer  :: j
  real(r8) :: phi
  integer  :: id_max_qc, id_max_qr
  real(r8) :: max_qc, max_qr, min_amk 

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(1)
  do j=1,ncd
      ank_c(j) = 0._r8
      amk_c(j) = 0._r8
      ank_r(j) = 0._r8
      amk_r(j) = 0._r8
      ank(j) = 0._r8
      amk(j) = 0._r8
  end do
  !$acc end parallel

   scale_nc = 0._r8
   scale_qc = 0._r8
   scale_nr = 0._r8
   scale_qr = 0._r8
   cutoff_amk = 0

   id_max_qc = 0
   id_max_qr = 0
   max_qc = 0._r8
   max_qr = 0._r8

   ! cloud water, nc in #/m3 --> #/cm3
   if (qc > qsmall) then
      !$acc loop seq
      do j=1,ncd
         phi = nc*lambda_c**(mu_c+1._r8)/ &
               gamma(mu_c+1._r8)*(diammean(j)*1.e-2_r8)**mu_c* &
               exp(-lambda_c*diammean(j)*1.e-2_r8)  ! D cm --> m
         ank_c(j) = phi*(diamedge(j+1)-diamedge(j))*1.e-2_r8   ! D cm --> m
         amk_c(j) = phi*(diamedge(j+1)-diamedge(j))*1.e-2_r8*mmean(j)*1.e-3_r8  ! mass in bin g --> kg
         scale_nc = scale_nc+ank_c(j)
         scale_qc = scale_qc+amk_c(j) 
      end do
      scale_nc = scale_nc/nc
      scale_qc = scale_qc/qc

      !$acc loop seq
      do j=1,ncd
         ank_c(j) = ank_c(j)/scale_nc
         amk_c(j) = amk_c(j)/scale_qc
         if ( amk_c(j) > max_qc ) then
            id_max_qc = j
            max_qc = amk_c(j)
         end if
      end do
   end if

   ! rain drop
   if (qr > qsmall) then
      !$acc loop seq
      do j=1,ncd
         phi = n0r*exp(-lambda_r*diammean(j)*1.e-2_r8)   ! D cm --> m
         ank_r(j) = phi*(diamedge(j+1)-diamedge(j))*1.e-2_r8   ! D cm --> m  
         amk_r(j) = phi*(diamedge(j+1)-diamedge(j))*1.e-2_r8*mmean(j)*1.e-3_r8
         scale_nr = scale_nr + ank_r(j)
         scale_qr = scale_qr + amk_r(j)
      end do
      scale_nr = scale_nr/nr
      scale_qr = scale_qr/qr

      !$acc loop seq
      do j=1,ncd
         ank_r(j) = ank_r(j)/scale_nr
         amk_r(j) = amk_r(j)/scale_qr
         if ( amk_r(j) > max_qr ) then
            id_max_qr = j
            max_qr = amk_r(j)
         end if
      end do
   end if

   !$acc loop seq
   do j=1,ncd
      amk(j) = amk_c(j) + amk_r(j)
      ank(j) = ank_c(j) + ank_r(j)
   end do

   if ( (id_max_qc > 0) .and. (id_max_qr > 0) ) then
      if ( (max_qc/max_qr < 10._r8) .or. (max_qc/max_qr > 0.1_r8) ) then
         min_amk = amk(id_max_qc)
              !$acc loop seq
         do j=id_max_qc,id_max_qr
            if ( amk(j) <= min_amk ) then
               cutoff_amk = j
               min_amk = amk(j)
            end if
         end do
      end if
   end if

!input: qc,nc,qr,nr, medge (bin edges). May also need # bins?
!output: amk, ank (mixing ratio and number in each bin)

!this part will take a bit of thinking about.
!use size distribution parameters (mu, lambda) to generate the values at discrete size points
!need to also ensure mass conservation  

end subroutine cam_bin_distribute

! here are the subroutines called above that actually do the collision-coalescence calculations:

! The Kernel is from Jerry from many moons ago (included)

! I read in the file data and multiply by the summed mass of the individual bins 
! (with a factor of 1.5 so that the values represent the middle of the bin

! 941 FORMAT(2X,E12.5)
!     READ(iunit,941) KNN(IDD,JDD)
!     KNN(IDD,JDD)=(XK_GR(IDD)*kkfac+XK_GR(JDD)*kkfac)*KNN(IDD,JDD)

!where idd and jdd are the indexes for the bins and xk_gr is the mass of drops in a bin in grams
!

!************************************************************************************
! Setup variables needed for collection
! Either pass in or define globally the following variables
! tbase(height) - temperature in K as a function of height
! rhon(height) - air density as a function of height in kg/m^3
! xk_gr(bins) - mass of single drop in each bin in grams
! lsmall - small number
! QC - mass mixing ratio in kg/kg
! QN - number mixing ratio in #/kg
! All parameters are defined to be global in my version so that they are readily available throughout the code:
! SMN0,SNN0,SMCN,APN,AMN2,AMN3,PSIN,FN,FPSIN,XPSIN,HPSIN,FN2,XXPSIN (all arrays of drop bins)
!************************************************************************************

!AG: Global arrays need to be passed around I think? Right now at the module level. Is that okay?

SUBROUTINE COMPUTE_COLL_PARAMS(rhon,xk_gr,qc,qn,gnnnn,gmnnn,lnnnn,lmnnn)

  !$acc routine seq

  IMPLICIT NONE

! variable declarations (added by hm, 020118)
! note: vertical array looping is stripped out, this subroutine operates
! only on LOCAL values

  real(r8), dimension(ncd) :: qc,qn
  real(r8), dimension(ncdp) :: xk_gr
  real(r8) :: tbase,rhon
  integer :: lk
  integer :: l
  real(r8), parameter :: lsmall = 1.e-12_r8
  real(r8), dimension(ncd) :: smn0,snn0,smcn,amn2,amn3,psin,fn,fpsin, &
                               xpsin,hpsin,fn2,xxpsin
  real(r8) :: apn

  real(r8), dimension(ncd) :: gnnnn,gmnnn,lnnnn,lmnnn
  integer :: lm1,ll

  lk=ncd

  DO L=1,LK
     SMN0(L)=QC(L)*RHON/1.E3_r8
     SNN0(L)=QN(L)*RHON/1.E6_r8

     IF(SMN0(L).LT.lsmall.OR.SNN0(L).LT.lsmall)THEN
        SMN0(L)=0.0_r8
        SNN0(L)=0.0_r8
     ENDIF
  ENDDO

  DO L=1,LK
     IF(SMN0(L) .gt. 0._r8.AND.SNN0(L) .gt. 0._r8)THEN
        SMCN(L)=SMN0(L)/SNN0(L)
        IF((SMCN(L) .GT. 2._r8*XK_GR(L)))THEN
           SMCN(L) = (2._r8*XK_GR(L))
        ENDIF
        IF((SMCN(L) .LT. XK_GR(L)))THEN
           SMCN(L) = XK_GR(L)
        ENDIF
     ELSE
        SMCN(L)=0._r8
     ENDIF
     IF (SMCN(L).LT.XK_GR(L).OR.SMCN(L).GT.(2._r8*XK_GR(L)).OR.SMCN(L).EQ.0.0_r8)THEN
        APN=1.0_r8
     ELSE
        APN=0.5_r8*(1._r8+3._r8*(XK_GR(L)/SMCN(L))-2._r8*((XK_GR(L)/SMCN(L))**2._r8))
     ENDIF

     IF(SNN0(L) .GT. LSMALL)THEN
        AMN2(L)=APN*SMN0(L)*SMN0(L)/SNN0(L)
        AMN3(L)=APN*APN*APN*SMN0(L)*SMN0(L)*SMN0(L)/(SNN0(L)*SNN0(L))
     ELSE
        AMN2(L)=0._r8
        AMN3(L)=0._r8
     ENDIF

     IF(SMCN(L).LT.XK_GR(L))THEN
        PSIN(L)=0.0_r8
        FN(L)=2._r8*SNN0(L)/XK_GR(L)
     ELSE
        IF(SMCN(L).GT.(2._r8*XK_GR(L)))THEN
           FN(L)=0.0_r8
           PSIN(L)=2._r8*SNN0(L)/XK_GR(L)
        ELSE
           PSIN(L)=2._r8/XK_GR(L)*(SMN0(L)/XK_GR(L)-SNN0(L))
           FN(L)=2._r8/XK_GR(L)*(2._r8*SNN0(L)-SMN0(L)/XK_GR(L))
        ENDIF
     ENDIF

     IF(SNN0(L).LT.LSMALL.OR.SMN0(L).LT.LSMALL)THEN
        PSIN(L)=0.0_r8
        FN(L)=0.0_r8
     ENDIF

     FPSIN(L)=0.5_r8/XK_GR(L)*(PSIN(L)-FN(L))
     XPSIN(L)=2._r8*XK_GR(L)*PSIN(L)
     HPSIN(L)=PSIN(L)-0.5_r8*FN(L)
     FN2(L)=FN(L)/2._r8

     IF(L.GT.1)THEN
        XXPSIN(L)=XK_GR(L)*PSIN(L-1)
     ENDIF
  ENDDO

!************************************************************************************
! Compute collision coalescence
! Either pass in or define globally the following variables
! Gain terms begin with G, loss terms begin with L
! Second letter defines mass (M) or number (N)
! Third and fourth letters define the types of particles colling, i.e., NN means drops colliding with drops
! Last letter defines the category the new particles go into, in this case just N for liquid drops
! The resulting rates are in units of #/cm^3/s and g/cm^3/s
! Relies on predefined kernel array KNN(bins,bins) - see top of this file
!************************************************************************************

   GMNNN = 0._r8
   GNNNN = 0._r8
   LMNNN = 0._r8
   LNNNN = 0._r8
! remove verical array index, calculate gain/loss terms locally

  DO L=3,LK-1
     LM1=L-1
     DO LL=1,L-2
        GNNNN(L)=GNNNN(L)+(PSIN(LM1)*SMN0(LL)-FPSIN(LM1)*AMN2(LL))*KNN(LM1,LL)
        GMNNN(L)=GMNNN(L)+(XK_GR(L)*PSIN(LM1)*SMN0(LL)+FN2(LM1)*AMN2(LL)-FPSIN(LM1)*AMN3(LL))*KNN(LM1,LL)
     ENDDO
  ENDDO

  DO L=2,LK-1
     LM1=L-1
     GNNNN(L)=GNNNN(L)+0.5_r8*SNN0(LM1)*SNN0(LM1)*KNN(LM1,LM1)
     GMNNN(L)=GMNNN(L)+0.5_r8*(SNN0(LM1)*SMN0(LM1)+SMN0(LM1)*SNN0(LM1))*KNN(LM1,LM1)
     DO LL=1,L-1
        LNNNN(L)=LNNNN(L)+(PSIN(L)*SMN0(LL)-FPSIN(L)*AMN2(LL))*KNN(L,LL)
        GMNNN(L)=GMNNN(L)+(SMN0(LL)*SNN0(L)-PSIN(L)*AMN2(LL)+FPSIN(L)*AMN3(LL))*KNN(L,LL)
        LMNNN(L)=LMNNN(L)+(XPSIN(L)*SMN0(LL)-HPSIN(L)*AMN2(LL))*KNN(L,LL)
     ENDDO
  ENDDO

  DO L=1,LK-1
     DO LL=L,LK-1
        LNNNN(L)=LNNNN(L)+(SNN0(LL)*SNN0(L))*KNN(LL,L)
        LMNNN(L)=LMNNN(L)+(SNN0(LL)*SMN0(L))*KNN(LL,L)
     ENDDO
  ENDDO

END SUBROUTINE COMPUTE_COLL_PARAMS


end module p3_stochastic_collect_tau


