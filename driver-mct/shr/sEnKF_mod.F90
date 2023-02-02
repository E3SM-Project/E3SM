!BOP
!   !MODULE: sEnKF_mod
!   !INTERFACE:
module sEnKF_mod
!   !DESCRIPTION:
!   ! Random subgrouping ensemble kalman filter
!   ! two steps of ensemble kalman filter: EAKF and EnKF with perturbed observations
!   !First step: calculate the increment at observation space
!   !Second step: Project the analysis increment @ observation space to the model space.

!   !REVISION HISTORY:
!    July 25 2018 Y. Liu <liu6@tamu.edu> initial release
!   !USES:
use shr_kind_mod, only : r8=> shr_kind_r8
!use mod_kinds        , only : r8
    implicit none
    public :: obs_increment,update_from_obs_inc,obs_inflation, &
              obs_increment_sEAKF, update_increment_sEAKF, &
              obs_increment_bia,update_from_obs_inc_bia,&
              EnKF_squen, update_EnKF, &
              sEnKF_squen, update_sEnKF, &
              get_ens_mean,get_ens_mean3,get_ens_corr,&
              inflation_RTPS,inflation_RTPP, inflation_corvariance,&
              cal_lcwt1d,cal_lcwt2d,cal_lcwt_latlon
    private
        real(r8), parameter :: cor_cutoff = 0.0

    contains

!==============================================================================
!BOP
    !ROUTAINE: obs_increment
!   !INTERFACE: 

    subroutine obs_increment(ens, obs,ens_size, obs_var, obs_inc)
!   !DESCRIPTION:
       ! sequential EAKF step one:
       ! calculate analysis increment @ observation space using EAKF
       ! Anderson (2013)doi:10.1175/1520-0493(2003)131,0634:ALLSFF.2.0.CO;2.

       !ens     : forecast observation
       !obs     : observation value
       !obs_var : observation error variance 

        implicit none
        integer, intent(in) :: ens_size
        real(r8), intent(in) :: ens(ens_size), obs, obs_var
        real(r8), intent(out) :: obs_inc(ens_size)
!EOP

!BOC
        real(r8) :: af_rate, cov
        real(r8) :: mean, new_cov, new_mean
        integer :: ie

        mean = sum(ens) / ens_size

        cov = 0.0
    
        do ie = 1, ens_size
           cov = cov + (ens(ie)-mean)*(ens(ie)-mean)
        end do
        
        cov = cov / (ens_size-1)

        if(cov >= 1.e-8) then                ! ignore  too tiny corvariace 
           new_cov = (cov *obs_var)/(cov + obs_var)
           new_mean = new_cov * (obs_var*mean + cov*obs)/(cov*obs_var)
           af_rate = sqrt(new_cov/cov)    ! the std rate between analysis and forecast
           obs_inc = af_rate * (ens - mean) + new_mean - ens
        else
          obs_inc = 0.0
        endif

    end subroutine obs_increment

!==============================================================================
!BOP
    !ROUTAINE: obs_inflation  @YL need revisit for exactly values
    !  
!   !INTERFACE: 

    subroutine obs_inflation(ens, obs,ens_size, obs_var0,  obs_var, distance )
!   !DESCRIPTION:
      !inflation observation error for weakly constrain because model could crash when constrained by an obs far from model (mainly caused by model error) 
      ! model need reject an obs too far away from the forecast values
      ! here we follow (Parrish and Derber 1992; Dee 1995) for observation variance
      ! D*D=B+R

       !ens     : forecast observation
       !obs     : observation value
       !obs_var0 : default observation error variance 
       !obs_var : inflated observation error variance 


        implicit none
        integer, intent(in) :: ens_size
        real(r8), intent(in) :: ens(ens_size), obs, obs_var0
        real(r8), intent(out) :: obs_var
        real(r8), optional,intent(out) :: distance
!EOP

!BOC
        real(r8) :: cov, mean, min_obs_var
        integer :: ie

        mean = sum(ens) / ens_size
        distance = abs(mean - obs)  ! the distance between ensemble mean and observation
        
        cov = 0.0
        do ie = 1, ens_size
           cov = cov + (ens(ie)-mean)*(ens(ie)-mean)
        end do
        cov = cov / (ens_size-1)        
        min_obs_var= distance*distance-cov;

        obs_var= max(min_obs_var, obs_var0)

        return               

    end subroutine obs_inflation


!==============================================================================
!BOP
    !ROUTAINE: obs_increment_sEAKF
!   !INTERFACE: 
    subroutine obs_increment_sEAKF(fobs,zob,EFsize,ens_size,err, &
             obs_inc,rll)
!   !DESCRIPTION:
        ! sequential sEAKF step two:
        ! calculate analysis increment @ observation space using sEAKF
        ! Liu 2014 https://search.proquest.com/docview/1667464678?pq-origsite=gscholar

        ! EFsize      : total ensemble size
        ! ens_size    : sub-ensmble size
        ! rll         : ensemble member order for subgrouping, need consistant with next step, project to model space
        ! fobs        : forecast observation
        ! zob         : observation
        implicit none
        integer, intent(in) :: EFsize, ens_size
        integer, intent(in) :: rll(EFsize)
        real(r8), intent(in) :: fobs(EFsize),zob,err
        real(r8), intent(out) :: obs_inc(EFsize)
!EOP

!BOC
        real(r8) ::  obs, obs_var
        real(r8) ::  ens(ens_size),gobs_inc(ens_size)
        integer :: i, j,ie, sub_size, kk(ens_size), grp

        obs=zob
        obs_var=err*err
        sub_size=ens_size
        grp = EFsize/ens_size
        do i=1,grp
           do j=1,ens_size
               kk(j)=rll(j+(i-1)*ens_size)
               ens(j)=fobs(kk(j))
           enddo
           call obs_increment(ens, obs,sub_size, obs_var, gobs_inc) 
           do j=1,ens_size
                obs_inc(kk(j))=gobs_inc(j)
           enddo
        enddo
!EOC
    end subroutine obs_increment_sEAKF


!==============================================================================
!BOP
    !ROUTAINE: obs_increment_bia
!   !INTERFACE: 
    subroutine obs_increment_bia(ens, obs,ens_size, obs_var, obs_inc, bia_var)
!   !DESCRIPTION:
        ! calculate analysis increment @ observation space using EAKF with biased  forecast uncertainty (bia_var)
        ! Y. Liu 2012 unfinished test 

         implicit none
         integer, intent(in) :: ens_size
         real(r8), intent(in) :: ens(ens_size), obs, obs_var, bia_var
         real(r8), intent(out) :: obs_inc(ens_size)
!EOP

!BOC
         real(r8) :: af_rate, cov
         real(r8) :: mean, new_cov, new_mean
         integer :: ie
         mean = sum(ens) / ens_size

         cov = 0.0

         do ie = 1, ens_size
           cov = cov + (ens(ie)-mean)*(ens(ie)-mean)
         end do

         cov = cov / (ens_size-1)

         if(cov >= 1.e-4) then
           new_cov = ((cov+bia_var) *obs_var)/(cov +bia_var+ obs_var)
           new_mean = new_cov * (obs_var*mean + (cov+bia_var)*obs)/((cov+bia_var)*obs_var)
           af_rate = sqrt(new_cov/cov)
           obs_inc = af_rate * (ens - mean) + new_mean - ens
         else
           obs_inc = 0.0
         endif
!EOC

    end subroutine obs_increment_bia

!==============================================================================
!BOP
    !ROUTAINE: update_from_obs_inc
!   !INTERFACE: 

    subroutine update_from_obs_inc(fobs, obs_inc, state, ens_size, state_inc, &
               wt_factor)
!   !DESCRIPTION:
       ! sequential EAKF step two:
       ! calculate model state analysis increment based on the analysis increment @ observation space
       ! Anderson (2013)doi:10.1175/1520-0493(2003)131,0634:ALLSFF.2.0.CO;2.

       !fobs      : forecast observation
       !obs_inc   : analysis increment at observation space
       !state     : state ensemble at model space
       !ens_size  : ensemble size
       !state_inc :  analysis state @ model space
       !wt_factor:  weighting effect for localization

       implicit none

       integer, intent(in) :: ens_size
       real(r8), intent(in) :: fobs(ens_size), obs_inc(ens_size)
       real(r8), intent(in) :: state(ens_size), wt_factor
       real(r8), intent(out) :: state_inc(ens_size)

!EOP
!BOC
       real(r8) :: mean_s, mean_o, cv_s, cv_o, cv, cr
       real(r8) ::s_prim(ens_size),o_prim(ens_size)
       integer :: ie

       mean_s = sum(state) / ens_size
       s_prim=state-mean_s
       mean_o = sum(fobs) / ens_size
       o_prim=fobs-mean_o

       cv_s = 0.0
       cv_o = 0.0
       cv   = 0.0

       do ie = 1, ens_size
        !  cv_s = cv_s + (state(ie) - mean_s)*(state(ie) - mean_s)
        !  cv_o = cv_o + (fobs(ie) - mean_o)*(fobs(ie) - mean_o)
        !  cv   = cv + (state(ie) - mean_s)*(fobs(ie) - mean_o)
          cv_s = cv_s + s_prim(ie)*s_prim(ie)
          cv_o = cv_o +o_prim(ie)*o_prim(ie)
          cv   = cv + o_prim(ie)*s_prim(ie)
       end do

       cv_s = cv_s / (ens_size-1)
       cv_o = cv_o / (ens_size-1)
       cv   = cv / (ens_size-1)
       if(cv_o > 1.e-8 .and. cv_s > 1.e-8) then
          cr   = cv /(sqrt(cv_s)*sqrt(cv_o))
       else
          cr = 0.0
       endif

       if(abs(cr) > cor_cutoff .and. cv_o > 1.e-8)then
          state_inc = (wt_factor * cv/cv_o) * obs_inc
       else
          state_inc = 0.0
       endif
!EOC
    end subroutine update_from_obs_inc

!==============================================================================
!BOP
    !ROUTAINE: update_from_obs_inc
!   !INTERFACE: 

    subroutine update_increment_sEAKF(fobs,obs_inc,state,state_inc, &
               EFsize,ens_size,rll,sprt,wt_factor)
!   !DESCRIPTION:
       ! sequential sEAKF step two:
       ! calculate model state analysis increment based on the analysis increment @ observation space
       ! Liu 2014 https://search.proquest.com/docview/1667464678?pq-origsite=gscholar
       ! EAKF,Anderson (2013)doi:10.1175/1520-0493(2003)131,0634:ALLSFF.2.0.CO;2.

       !fobs      : forecast observation
       !obs_inc   : analysis increment at observation space
       !state     : state ensemble at model space
       !state_inc :  analysis state @ model space
       !EFsize    : full ensemble size
       !ens_size  : subgroup ensemble size
       !rll       : ensemble order for subgrouping
       !sprt      : switch for updating using 0 (full ensemble corvariace), 1(subgroup ensemble corvariaces) 
       !wt_factor :  weighting effect for localization

       implicit none

       integer, intent(in) :: EFsize, sprt,ens_size
       integer, intent(in) :: rll(EFsize)
       real(r8), intent(in) :: fobs(EFsize), obs_inc(EFsize)
       real(r8), intent(in) :: state(EFsize), wt_factor
       real(r8), intent(out) :: state_inc(EFsize)
!EOP
!BOC
       real(r8) :: obs(ens_size), gobs_inc(ens_size),cov_factor, &
               gstate(ens_size), gstate_inc(ens_size)
       integer :: grp , i ,j ,kk(ens_size)

       grp = EFsize/ens_size

       if( grp==1 .OR. sprt ==0) then
           call update_from_obs_inc(fobs, obs_inc, state, EFsize, &
                state_inc,   wt_factor)
       else

          do i=1,grp
             do j=1,ens_size
                kk(j)=rll(j+(i-1)*ens_size)
                obs(j)=fobs(kk(j))
                gobs_inc(j)=obs_inc(kk(j))
                gstate(j)=state(kk(j))
             enddo
             call update_from_obs_inc(obs, gobs_inc, gstate, ens_size, &
                  gstate_inc, wt_factor)
             do j=1,ens_size
                  state_inc(kk(j))=gstate_inc(j)
             enddo
          enddo
       endif
!EOC
    end subroutine update_increment_sEAKF
!==============================================================================
!BOP
    !ROUTAINE: update_from_obs_inc_bia
!   !INTERFACE: 

    subroutine update_from_obs_inc_bia(fobs, obs_inc, state, ens_size, state_inc, &
               wt_factor,obs_bia,state_bia,r)
!   !DESCRIPTION:
       ! sequential EAKF with negtive biased forecast uncertianty step two:
       ! calculate model state analysis increment based on the analysis increment @ observation space
       ! Y. Liu 2012, unfinished test 
       ! EAKF,Anderson (2013)doi:10.1175/1520-0493(2003)131,0634:ALLSFF.2.0.CO;2.

       !fobs      : forecast observation
       !obs_inc   : analysis increment at observation space
       !state     : state ensemble at model space
       !ens_size  : ensemble size
       !state_inc : analysis state @ model space
       !wt_factor : weighting effect for localization
       !obs_bia   : missing error variance of forecast observation
       !state_bia : missing error variance of forecast state
       !r         ; correlation corefficent of obs_bia and state_bia

       implicit none

       integer, intent(in) :: ens_size
       real(r8), intent(in) :: fobs(ens_size), obs_inc(ens_size)
       real(r8), intent(in) :: state(ens_size), wt_factor
       real(r8), intent(out) :: state_inc(ens_size)
       real(r8), intent(in) :: obs_bia,state_bia,r

       real(r8) :: mean_s, mean_o, cv_s, cv_o, cv, cr,cv_b
       integer :: ie

       mean_s = sum(state) / ens_size
       mean_o = sum(fobs) / ens_size

       cv_s = 0.0
       cv_o = 0.0
       cv   = 0.0

       do ie = 1, ens_size
         cv_s = cv_s + (state(ie) - mean_s)*(state(ie) - mean_s)
         cv_o = cv_o + (fobs(ie) - mean_o)*(fobs(ie) - mean_o)
         cv   = cv + (state(ie) - mean_s)*(fobs(ie) - mean_o)
       end do
       cv_s = cv_s / (ens_size-1)
       cv_o = cv_o / (ens_size-1)
       cv   = cv / (ens_size-1)
       cr   = cv /(sqrt(cv_s)*sqrt(cv_o))

       cv_b=r*sqrt(obs_bia*state_bia)

       if(abs(cr) > cor_cutoff .and. cv_o > 1.e-8)then
          state_inc=(wt_factor *(cv+cv_b)/(cv_o+obs_bia)) *obs_inc
       else
          state_inc = 0.0
       endif
!EOC
    end subroutine update_from_obs_inc_bia



!==============================================================================
!BOP
    !ROUTAINE: EnKF_squen
!   !INTERFACE: 
    subroutine EnKF_squen(ens,obs,ens_size,err,Kinv,Iens)
!   !DESCRIPTION:
       ! a sequential EnKF (pertuberd observation) with two steps setting
       use ran_mod, only : normal
       integer, intent(in) :: ens_size
       real(r8), intent(in) :: ens(ens_size), err,obs
       real(r8), intent(out) :: Kinv, Iens(ens_size)
!EOP
!BOC
       real(r8) :: cov, mean, pert_obs(ens_size)
       integer :: ie

       mean = sum(ens) / ens_size
       cov = 0.0

       do ie = 1, ens_size
         cov = cov + (ens(ie)-mean)*(ens(ie)-mean)
         pert_obs(ie)= normal() * err    ! random noice for observation
       end do
       cov = cov / (ens_size-1)
       pert_obs=pert_obs-sum(pert_obs)/ens_size+obs !  removed the peturbation mean
       if(cov >= 1.e-8) then
           Kinv= 1.0/(cov+err*err)
           Iens=Iens-sum(Iens)/ens_size+obs-ens
           Iens=pert_obs-ens
       else
           Kinv=0
           Iens=0
       endif
!EOC
    end subroutine EnKF_squen


!==============================================================================
!BOP
    !ROUTAINE: sEnKF_squen
!   !INTERFACE: 
    subroutine sEnKF_squen(ens,obs,EFsize,ens_size,err,rll,spd,Iens)
!   !DESCRIPTION:
       ! a sequential subgrouping EnKF (pertuberd observation) with two steps setting
       ! spd is not one, subgrouping method from Houtekamer and Mitchell (1998)
       use ran_mod, only : normal
       integer, intent(in) :: ens_size,spd,EFsize
       real(r8), intent(in) :: ens(EFsize), err,obs
       real(r8), intent(out) :: Iens(EFsize)
       integer ,intent(in) :: rll(EFsize)
!EOP
!BOC
       real(r8) :: cov, mean,fens(ens_size),IIens(ens_size),Kinv
       integer :: ie,i,k,grp,kk,spd1

       grp=EFsize/ens_size
       spd1=spd
       if(grp==1)spd1=1
       do i=1,grp
          do k=1,ens_size
             kk=ens_size*(i-1)+k
             fens(k)=ens(rll(kk))
          enddo

          if (spd1==1 )then
             mean = sum(fens) / ens_size
             cov = sum((fens-mean)**2)
             cov=cov/(ens_size-1)
          else
!            Houtekamer and Mitchell
             mean = (sum(ens)-sum(fens)) / (EFsize-ens_size)
             cov = sum((ens-mean)**2)-sum((fens-mean)**2) 
             cov=cov/(EFsize-ens_size-1)
          endif

          do ie = 1, ens_size
            IIens(ie)=normal()*err
          end do
          Kinv= 1.0/(cov+err*err)
          IIens=IIens-sum(IIens)/ens_size+obs-fens
          IIens=Kinv*IIens
          do k=1,ens_size
             kk=ens_size*(i-1)+k
             Iens(rll(kk))=IIens(k)
          enddo
       enddo
!EOC
    end subroutine sEnKF_squen

!==============================================================================
!BOP
    !ROUTAINE: update_EnKF
!   !INTERFACE: 
    subroutine update_EnKF(fens,Kinv,Iens,state,ens_size,cov_factor)
!   !DESCRIPTION:
       ! a sequential EnKF (pertuberd observation) with two steps setting
       ! Step 2 : update model state, ** here we are not calculate the increment only, wich is different than the EAKF
       ! fens     : forecast observation
       ! state    : model state

!
       integer, intent(in) :: ens_size
       real(r8), intent(in) :: fens(ens_size),Kinv, Iens(ens_size),cov_factor
       real(r8), intent(inout) :: state(ens_size)
!EOP
!BOC
       real(r8) :: cov, mean1,mean2
       integer :: ie,j,k

       mean1=sum(fens)/ens_size
       mean2=sum(state)/ens_size
       cov=0.0;
       do ie=1,ens_size
         cov = cov + (fens(ie)-mean1)*(state(ie)-mean2)
       enddo
       cov=cov/(ens_size-1)
       state=state+(cov*Kinv*cov_factor)*Iens
!EOC
    end subroutine update_EnKF

!==============================================================================
!BOP
    !ROUTAINE: update_sEnKF
!   !INTERFACE: 
    subroutine update_sEnKF(ens,Iens,state,EFsize,ens_size,rll,spd,wt)
!   !DESCRIPTION:
       ! a sequential subgrouping EnKF (pertuberd observation) with two steps setting
       ! Step 2 : update model state, ** here we are not calculate the increment only, wich is different than the EAKF
       ! ens     : forecast observation
       ! state    : model state
       integer, intent(in) :: EFsize,ens_size,spd
       real(r8), intent(in) :: ens(EFsize), Iens(EFsize),wt
       integer, intent(in) :: rll(EFsize)
       real(r8), intent(inout) :: state(EFsize)
!EOP
!BOC
       real(r8) :: cov, mean1,mean2,IIens(ens_size),fstate(ens_size),fens(ens_size)
       integer :: ie,i,k,j,grp,kk,spd1

       spd1=spd;
       grp=EFsize/ens_size
       if(grp==1)spd1=1
       do i=1,grp
          do k=1,ens_size
              kk=ens_size*(i-1)+k
              fstate(k)=state(rll(kk))
              fens(k)=ens(rll(kk))
              IIens(k)=Iens(rll(kk))
          enddo

          if (spd1==1)then
              mean1 = sum(fens) / ens_size
              mean2=sum(fstate)/ ens_size
              cov = sum((fens-mean1)*(fstate-mean2))/(ens_size-1)
          else
              mean1 = (sum(ens)-sum(fens)) / (EFsize-ens_size)
              mean2 = (sum(state)-sum(fstate)) / (EFsize-ens_size)
              cov = sum((ens-mean1)*(state-mean2)) &
                    -sum((fens-mean1)*(fstate-mean2))
              cov=cov/(EFsize-ens_size-1)
          endif

          fstate=fstate+(cov*wt)*IIens
          do k=1,ens_size
              kk=ens_size*(i-1)+k
              state(rll(kk))=fstate(k)
          enddo

       enddo
!EOC
    end subroutine update_sEnKF
!==============================================================================
!BOP
    !ROUTAINE: inflation_RTPP
!   !INTERFACE: 
    subroutine inflation_RTPP(anal,prior,alpha,ens_size)
!   !DESCRIPTION:
       ! the inflation scheme of relaxiation to prior (Fuqing Zhang's method) 
       !doi:10.1175/1520-0493(2004)132<1238:IOIEAO>2.0.CO;2
       ! alfa    : < 1.0, the partitial using prior spread
       implicit none
       integer, intent(in) :: ens_size
       real(r8) , intent(in) :: prior(ens_size), alpha
       real(r8), intent(inout) :: anal(ens_size) 
!EOP
!BOC
       real(r8) :: umean,xmean,xprm(ens_size),uprm(ens_size) 

       xmean=sum(prior)/ens_size
       umean=sum(anal)/ens_size
       xprm=prior-xmean
       uprm=anal-umean
       anal = uprm*(1-alpha)+xprm*alpha+umean
!EOC
    end subroutine inflation_RTPP

!==============================================================================
!BOP
    !ROUTAINE: inflation_RTPS
!   !INTERFACE: 
    subroutine inflation_RTPS(anal,prior,alpha,ens_size)
!   !DESCRIPTION:
      !Whitaker and Hamill 2012
      ! relaxiation to prior spread
      !x_u'= x_a'*(1+alpha*(sigma_p/sigma_a-1))
      implicit none
      integer, intent(in) :: ens_size
      real(r8) , intent(in) :: prior(ens_size),alpha
      real(r8), intent(inout) :: anal(ens_size)
!EOP
!BOC
      real(r8) :: prior_mean(2),anal_mean(2),anal_anom(ens_size)
      ! *_mean(1) is mean, *_mean(2) is std

      call get_ens_mean(anal,anal_mean,ens_size)
      call get_ens_mean(prior,prior_mean,ens_size)
      if(prior_mean(2) <= anal_mean(2) )  return  ! the observation did not affect the state
      anal_anom=anal-anal_mean(1)
      anal=anal_anom*(1+alpha*(prior_mean(2)/anal_mean(2)-1)) + anal_mean(1)
      return
!EOC
    end subroutine inflation_RTPS

!==============================================================================
!BOP
    !ROUTAINE: inflation_corvariance
!   !INTERFACE: 
    subroutine inflation_corvariance(x,rt,k,ens_size)
!   !DESCRIPTION:
       ! the multiplive  corvariace  inflation scheme of Anderson and Anderson (1999)

       integer, intent(in) ::  k,ens_size 
       real(r8), intent (in) :: rt 
       real(r8), intent(inout) :: x(k,ens_size)
!EOP
!BOC
       real(r8) :: xmean,xprm(ens_size)
       integer :: i

       do i=1,k
          xmean=sum(x(k,:))/ens_size
          xprm=x(k,:)-xmean
          x(k,:)=xprm*rt+xmean
       enddo
!EOC
    end subroutine inflation_corvariance
!==============================================================================
!BOP
    !ROUTAINE: cal_lcwt2d
!   !INTERFACE: 
    subroutine cal_lcwt2d(wt,lengthx,lengthy,R,xdim,ydim)
!   !DESCRIPTION:
       ! calculate corvariance localization weight in 2D based on Gaspari and Cohn (1999)
       ! doi:10.1002/qj.49712555417

       ! R: localization radus.
       ! wt: covariance weight coefficeint
       ! assume the observation location(0,0), the horizontal wt pattern is symmetry 
       implicit none
       integer, intent(in) :: xdim,ydim
       real(r8), intent(in) :: lengthx,lengthy,R
       real(r8), intent(inout):: wt(0:xdim,0:ydim)
!EOP

!BOC
       real(r8) :: dist
       integer :: i,k,m,n
       wt=0.0
       if (R < tiny(R)) then
         wt(0,0)=1.0
         return
       endif

       do k=0,ydim
         do i=0,xdim
            dist=sqrt((i*lengthx)**2+(k*lengthy)**2)/R
            if(dist>2.0) CYCLE

            if(dist<=1.0)then
               wt(i,k)=1+(((0.5-0.25*dist)*dist+0.625)*dist-1.6667)*dist*dist
            else
               wt(i,k)=((((0.0833*dist-0.5)*dist+0.625)*dist+1.6667)*dist-5.0)*dist &
                          + 4-0.6667/dist
            endif
         enddo
       enddo
!EOC
    end subroutine cal_lcwt2d

!==============================================================================
!BOP
    !ROUTAINE: cal_lcwt1d
!   !INTERFACE: 
    subroutine cal_lcwt1d(wt,length,R,xdim)
!   !DESCRIPTION:
       ! calculate corvariance localization weight in 1D based on Gaspari and Cohn (1999)
       ! doi:10.1002/qj.49712555417

       ! R: localization radus.
       ! wt: covariance weight coefficeint
       ! assume the observation location at 0, the  wt pattern is symmetry 

       implicit none
       integer, intent(in) :: xdim
       real(r8), intent(in) :: length,R
       real(r8), intent(inout):: wt(0:xdim)
!EOP

!BOC
       real(r8) :: dist
       integer :: i,k,m,n
       wt=0.0

       do i=0,xdim
          dist=(i*length)/R
          if(dist>=2.0) CYCLE
          if(dist<=1.0)then
               wt(i)=1+(((0.5-0.25*dist)*dist+0.625)*dist-1.6667)*dist*dist
          else
               wt(i)=((((0.0833*dist-0.5)*dist+0.625)*dist+1.6667)*dist-5.0)*dist &
                 + 4-0.6667/dist
          endif
       enddo
!EOC
    end subroutine cal_lcwt1d

!==============================================================================
!BOP
    !ROUTAINE: cal_covfactor
!   !INTERFACE: 
subroutine cal_lcwt_latlon(lon1,lat1,lon2,lat2,RR,wgt, ifdegree)
!   !DESCRIPTION:
       ! calculate corvariance localization weight on lat lon grid  based on Gaspari and Cohn (1999)
       ! doi:10.1002/qj.49712555417

       ! use the distance from (lon1,lat1) to (lon2,lat2),
       !  haversine formula
       ! lon1,lat1:  first grid in dgree/radius
       ! lon2,lat2:  second grid in degree/radius
       ! RR      :localization radus.
       ! wgt: covariance weight coefficeint
       implicit none
       real(r8),intent(in) :: lon1,lat1,lon2,lat2
       real(r8),intent(in) :: RR
       real(r8),intent(out) :: wgt
       logical, optional, intent(in) :: ifdegree  
!EOP
!BOC
       real(r8) :: lonb,latb,lone,late
       real(r8) :: RE,PI,Y,dist

       RE = 6430.0
       PI = 3.1415926
       if ( present(ifdegree) .and. ifdegree) then
           latb = lat1/180*PI;
           late = lat2/180*PI;
           lonb = lon1/180*PI;
           lone = lon2/180*PI;
       else
           latb = lat1
           late = lat2
           lonb = lon1
           lone = lon2
       endif

       Y = 2*asin(sqrt(sin((late-latb)/2)**2+cos(latb)*cos(late)*sin((lonb-lone)/2)**2));
       ! print *,'Y=',Y
       Y = Y*RE
       dist = Y/RR
       if (Y<=RR) then
            wgt = 1.0-1.0/4.0*dist**5+1.0/2.0*dist**4+5.0/8.0*dist**3-5.0/3.0*dist**2
       else
            if(Y>=2*RR) then
                wgt = 0.0
            else
                wgt = 1.0/12.0*(dist**5)-1.0/2.0*(dist**4)+5.0/8.0*(dist**3)+5.0/3.0*(dist**2)-5.0*dist+4.0-2.0/3.0/dist
           endif
       endif
!EOF
    end subroutine cal_lcwt_latlon

!==============================================================================
!BOP
    !ROUTAINE: get_ens_mean
!   !INTERFACE: 
    subroutine get_ens_mean(X,Y,ens_size)
!   !DESCRIPTION:
       ! calculate ensemble mean and standard derivation
        
       integer, intent(in) :: ens_size 
       real(r8) , intent (in) :: X(ens_size)
       real(r8), intent (inout) :: Y(2)
!EOP
!BOC
       real(r8) :: xmean, xprm(ens_size)
       xmean=sum(X)/ens_size
       Y(1)=xmean
       xprm=(X-xmean)**2
       Y(2)=sqrt(sum(xprm)/(ens_size-1))
!EOC
    end subroutine get_ens_mean


!==============================================================================
!BOP
    !ROUTAINE: get_ens_mean3
!   !INTERFACE: 
    subroutine get_ens_mean3(X,Y,ens_size)
!   !DESCRIPTION:
       ! calculate ensemble mean and standard derivation,kurtosis
       integer, intent(in) :: ens_size
       real(r8) , intent (in) :: X(ens_size)
       real(r8), intent (out) :: Y(3)
!EOP
!BOC
       real(r8) :: xmean, xprm(ens_size),yy(3)
       xmean=sum(X)/ens_size
       Y(1)=xmean
       xprm=(x-xmean)**2
       Y(2)=sqrt(sum(xprm)/(ens_size-1))
       Y(3)=sum(xprm**2)/(y(2)**4)/(ens_size-1)
!EOC
    end subroutine get_ens_mean3


!==============================================================================
!BOP
    !ROUTAINE: get_ens_corr
!   !INTERFACE: 
    subroutine get_ens_corr(X,X1,Y,ens_size)
!   !DESCRIPTION:
       ! calculate ensemble correlation
       integer, intent(in) :: ens_size
       real(r8) , intent (in) :: X(ens_size),X1(ens_size)
       real(r8), intent (inout) :: Y
!EOP
!EOC
       !local
       real(r8) :: xmean, xprm(ens_size),x1mean,x1prm(ens_size)

       xmean=sum(X)/ens_size
       x1mean=sum(X1)/ens_size
       xprm=(x-xmean)
       x1prm=(X1-x1mean)
       Y=sum(xprm*x1prm)/sqrt(sum(xprm*xprm)*sum(x1prm*x1prm))
!EOC
    end subroutine get_ens_corr

end module sEnKF_mod
