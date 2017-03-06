module tfilt_massfix

contains

subroutine tfilt_massfixrun (ztodt,         lat,    u3m1,   u3,     &
                          v3m1,   v3,    t3m1,   t3,     q3m1,   &
                          q3,     psm1,  ps,             alpha,  &
                          etamid, qfcst, vort,   div,    vortm2, &
                          divm2,         qminus, psm2,   um2,    &
                          vm2,    tm2,   qm2,    vortm1, divm1,  &
                          omga,   dpsl,  dpsm,   beta,   hadv ,  &
                          nlon,   pdeldry, pdelm1dry, pdelm2dry)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Time filter (second half of filter for vorticity and divergence only)
! 
! Method: 
! 
! Author: 
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use cam_history,  only: outfld
   use eul_control_mod, only : fixmas,eps
   use pmgrid,       only: plon, plev, plevp, plat
   use pspect
   use commap
   use constituents, only: pcnst, qmin, cnst_cam_outfld, &
                           tottnam, tendnam, cnst_get_type_byind, fixcnam, &
                           hadvnam, vadvnam
   use time_manager, only: get_nstep
   use physconst,    only: cpair, gravit
   use scamMod
   use phys_control, only: phys_getopts
   
#if ( defined BFB_CAM_SCAM_IOP )
   use iop
   use constituents, only: cnst_get_ind, cnst_name
#endif
   implicit none

!
! Input arguments
!
   real(r8), intent(in) :: ztodt                  ! two delta t (unless nstep<2)

   real(r8), intent(inout) :: qfcst(plon,plev,pcnst)! slt moisture forecast
   real(r8), intent(in) :: vort(plon,plev)
   real(r8), intent(in) :: div(plon,plev)
   real(r8), intent(inout) :: vortm2(plon,plev)
   real(r8), intent(inout) :: divm2(plon,plev)
   real(r8), intent(in) :: qminus(plon,plev,pcnst)
   real(r8), intent(inout) :: psm2(plon)
   real(r8), intent(inout) :: um2(plon,plev)
   real(r8), intent(inout) :: vm2(plon,plev)
   real(r8), intent(inout) :: tm2(plon,plev)
   real(r8), intent(inout) :: qm2(plon,plev,pcnst)
   real(r8), intent(inout) :: omga(plon,plev)
   real(r8), intent(in) :: dpsl(plon)
   real(r8), intent(in) :: dpsm(plon)
   real(r8), intent(in) :: beta                   ! energy fixer coefficient
   real(r8), intent(in) :: hadv(plon,plev,pcnst)  ! horizonal q advection tendency
   real(r8), intent(in) :: alpha(pcnst)
   real(r8), intent(in) :: etamid(plev)           ! vertical coords at midpoints 
   real(r8), intent(in) :: u3(plon,plev)
   real(r8), intent(in) :: v3(plon,plev)
   real(r8), intent(inout) :: t3(plon,plev)
   real(r8), intent(inout) :: pdeldry(:,:)      ! dry pressure difference at time n3
   real(r8), intent(inout) :: pdelm1dry(:,:)    ! dry pressure difference at time n3m1
   real(r8), intent(in) :: pdelm2dry(:,:)       ! dry pressure difference at time n3m2


   integer, intent(in) :: lat
   integer, intent(in) :: nlon

! Input/Output arguments

   real(r8), intent(inout) :: q3(plon,plev,pcnst)
   real(r8), intent(inout) :: ps(plon)
   real(r8), intent(inout) :: vortm1(plon,plev)
   real(r8), intent(inout) :: psm1(plon)
   real(r8), intent(inout) :: u3m1(plon,plev)
   real(r8), intent(inout) :: v3m1(plon,plev)
   real(r8), intent(inout) :: t3m1(plon,plev)
   real(r8), intent(inout) :: divm1(plon,plev)
   real(r8), intent(inout) :: q3m1(plon,plev,pcnst)
!
! Local workspace
!
   integer ifcnt                   ! Counter
   integer :: nstep                ! current timestep number
   integer :: timefiltstep         ! 

   real(r8) tfix    (plon)        ! T correction
   real(r8) engycorr(plon,plev)   ! energy equivalent to T correction
   real(r8) rpmid(plon,plev)      ! 1./pmid
   real(r8) pdel(plon,plev)       ! pdel(k)   = pint  (k+1)-pint  (k)
   real(r8) pint(plon,plevp)      ! pressure at model interfaces (n  )
   real(r8) pmid(plon,plev)       ! pressure at model levels (time n)
   real(r8) utend(plon,plev)      ! du/dt
   real(r8) vtend(plon,plev)      ! dv/dt
   real(r8) ttend(plon,plev)      ! dT/dt
   real(r8) qtend(plon,plev,pcnst)! dq/dt
   real(r8) pstend(plon)          ! d(ps)/dt
   real(r8) vadv  (plon,plev,pcnst) ! vertical q advection tendency
   real(r8) pintm1(plon,plevp)    ! pressure at model interfaces (n-1)
   real(r8) pmidm1(plon,plev)     ! pressure at model levels (time n-1)
   real(r8) pdelm1(plon,plev)     ! pdelm1(k) = pintm1(k+1)-pintm1(k)
   real(r8) om2eps
   real(r8) corm
   real(r8) wm
   real(r8) absf
   real(r8) worst
   logical lfixlim               ! flag to turn on fixer limiter

   real(r8) ta(plon,plev,pcnst)   ! total advection of constituents
   real(r8) dqfx3(plon,plev,pcnst)! q tendency due to mass adjustment
   real(r8) coslat                 ! cosine(latitude)
   real(r8) rcoslat(plon)         ! 1./cosine(latitude)
!   real(r8) engt                   ! Thermal   energy integral
!   real(r8) engk                   ! Kinetic   energy integral
!   real(r8) engp                   ! Potential energy integral
   integer i, k, m,j,ixcldliq,ixcldice,ixnumliq,ixnumice
#if ( defined BFB_CAM_SCAM_IOP )
   real(r8) :: t3forecast(plon,plev),delta_t3(plon,plev)
   real(r8) :: q3forecast(plon,plev,pcnst),delta_q3(plon,plev,pcnst)
#endif
   real(r8) fixmas_plon(plon)
   real(r8) beta_plon(plon) 
   real(r8) clat_plon(plon)
   real(r8) alpha_plon(plon)

!-----------------------------------------------------------------------
   nstep = get_nstep()
#if ( defined BFB_CAM_SCAM_IOP )
!
! Calculate 3d dynamics term
!
   do k=1,plev
      do i=1,nlon
         divt3dsav(i,k,lat)=(t3(i,k)-tm2(i,k))/ztodt -t2sav(i,k,lat)
         t3forecast(i,k)=tm2(i,k)+ztodt*t2sav(i,k,lat)+ztodt*divt3dsav(i,k,lat)
      end do
   end do
   do i=1,nlon
      do m=1,pcnst
         do k=1,plev
            divq3dsav(i,k,m,lat)= (qfcst(i,k,m)-qminus(i,k,m))/ztodt
            q3forecast(i,k,m)=qminus(i,k,m)+divq3dsav(i,k,m,lat)*ztodt
         end do
      end do
   end do


   q3(:nlon,:,:)=q3forecast(:nlon,:,:)
   t3(:nlon,:)=t3forecast(:nlon,:)
   qfcst(:nlon,:,:)=q3(:nlon,:,:)

!
! outflds for iop history tape - to get bit for bit with scam
! the n-1 values are put out.  After the fields are written out
! the current time level of info will be buffered for output next
! timestep
!
   call outfld('t',t3  ,plon   ,lat     )
   call outfld('q',q3  ,plon   ,lat     )
   call outfld('Ps',ps ,plon   ,lat     )
   call outfld('u',u3  ,plon   ,lat     )
   call outfld('v',v3  ,plon   ,lat     )
!
! read single values into plon arrays for output to history tape
! it would be nice if history tape supported 1 dimensional array variables
!
   fixmas_plon(:)=fixmas
   beta_plon(:)=beta
   clat_plon(:)=clat(lat)

   call outfld('fixmas',fixmas_plon,plon   ,lat     )
   call outfld('beta',beta_plon  ,plon   ,lat     )
   call outfld('CLAT    ',clat_plon  ,plon   ,lat     )
   call outfld('divT3d',divt3dsav(1,1,lat)  ,plon   ,lat     )
   do m =1,pcnst
      call outfld(trim(cnst_name(m))//'_dten',divq3dsav(1,1,m,lat)  ,plon   ,lat     )
   end do
#endif


   coslat = cos(clat(lat))
   do i=1,nlon
     rcoslat(i) = 1._r8/coslat
   enddo
   lfixlim = .true.


!
! Set average dry mass to specified constant preserving horizontal
! gradients of ln(ps). Proportionality factor was calculated in STEPON
! for nstep=0 or SCAN2 otherwise from integrals calculated in INIDAT
! and SCAN2 respectively.
! Set p*.
!
   do i=1,nlon
      ps(i) = ps(i)*fixmas
   end do
!
! Set current time pressure arrays for model levels etc.
!
   call plevs0(nlon    ,plon   ,plev    ,ps      ,pint    ,pmid    ,pdel)
!
!$OMP PARALLEL DO PRIVATE (K, I)
   do k=1,plev
      do i=1,nlon
         rpmid(i,k) = 1._r8/pmid(i,k)
      enddo
   enddo
!
! Add temperature correction for energy conservation
!
!$OMP PARALLEL DO PRIVATE (K, I)
   do k=1,plev
      do i=1,nlon
         engycorr(i,k) = (cpair/gravit)*beta*pdel(i,k)/ztodt
         t3      (i,k) = t3(i,k) + beta
      end do
   end do
   do i=1,nlon
      tfix(i) = beta/ztodt
   end do
!
! Output Energy correction term
!
! using do loop and select in order to enable functional parallelism with OpenMP
!$OMP PARALLEL DO PRIVATE (I)
   do i=1,2
      select case (i)
      case (1)
         call outfld ('ENGYCORR',engycorr ,plon   ,lat     )
      case (2)
         call outfld ('TFIX    ',tfix     ,plon   ,lat     )
      end select
   end do

!
! Compute q tendency due to mass adjustment
! If LFIXLIM = .T., then:
! Check to see if fixer is exceeding a desired fractional limit of the
! constituent mixing ratio ("corm").  If so, then limit the fixer to
! that specified limit.
!
   do m=1,pcnst
      if (cnst_get_type_byind(m).eq.'dry' ) then
         corm    = 1.e36_r8
      else
         corm    = 0.1_r8
      end if

!$OMP PARALLEL DO PRIVATE (K, I, IFCNT, WORST, WM, ABSF)
      do k=1,plev
         do i=1,nlon
            if (single_column) then
               dqfx3(i,k,m) = dqfxcam(i,k,m)
            else
               dqfx3(i,k,m) = alpha(m)*etamid(k)*abs(qfcst(i,k,m) - qminus(i,k,m))
#if ( defined BFB_CAM_SCAM_IOP )
               dqfx3sav(i,k,m,lat) = dqfx3(i,k,m)
#endif
            endif
         end do
         if (lfixlim) then
            ifcnt = 0
            worst = 0._r8
            wm    = 0._r8
            do i = 1,nlon
               absf = abs(dqfx3(i,k,m))
               if (absf.gt.corm) then
                  ifcnt = ifcnt + 1
                  worst = max(absf,worst)
                  wm = wm + absf
                  dqfx3(i,k,m) = sign(corm,dqfx3(i,k,m))
               endif
            end do
            if (ifcnt.gt.0) then
               wm = wm/real(ifcnt,r8)

! TBH:  Commented out as of CAM CRB meeting on 6/20/03
!              write(iulog,1000) m,corm,ifcnt,k,lat,wm,worst

            endif
         endif
         do i=1,nlon
            dqfx3(i,k,m) = qfcst(i,k,m)*dqfx3(i,k,m)/ztodt
            q3   (i,k,m) = qfcst(i,k,m) + ztodt*dqfx3(i,k,m)
            ta   (i,k,m) = (q3    (i,k,m) - qminus(i,k,m))/ztodt
            vadv (i,k,m) = (qfcst(i,k,m) - qminus(i,k,m))/ztodt - hadv(i,k,m)
         end do
      end do
   end do

!$OMP PARALLEL DO PRIVATE (K, I)
   do k=1,plev
      do i=1,nlon
         pdeldry(i,k) = pdel(i,k)*(1._r8-q3(i,k,1))
      end do ! i
   end do ! k


#if ( defined BFB_CAM_SCAM_IOP )
   do m=1,pcnst
      alpha_plon(:)= alpha(m)
      call outfld(trim(cnst_name(m))//'_alph',alpha_plon ,plon   ,lat     )
      call outfld(trim(cnst_name(m))//'_dqfx',dqfx3sav(1,1,m,lat) ,plon   ,lat     )
   end do
#endif
!
! Check for and correct invalid constituents
!
   call qneg3 ('TFILT_MASSFIX',lat   ,nlon    ,plon   ,plev    , &
               1, pcnst, qmin ,q3(1,1,1),.True.)
!
! Send slt tendencies to the history tape
!
!$OMP PARALLEL DO PRIVATE (M)
   do m=1,pcnst
      if ( cnst_cam_outfld(m) ) then
         call outfld(tottnam(m),ta(1,1,m),plon   ,lat     )
      end if
   end do
   if (.not. single_column) then
!
! Calculate vertical motion field
!
      call omcalc (rcoslat ,div     ,u3      ,v3      ,dpsl    ,  &
                dpsm    ,pmid    ,pdel    ,rpmid   ,pint(1,plevp), &
                omga    ,nlon    )

   endif

!   write(iulog,*)'tfilt: lat=',lat
!   write(iulog,*)'omga=',omga
!
! Time filter (second half of filter for vorticity and divergence only)
!
!   if(lat.eq.2) then
!      write(iulog,*)'tfilt: ps=',psm2(13),psm1(13),ps(13)
!      write(iulog,*)'tfilt: u=',um2(13,18),u3m1(13,18),u3(13,18)
!      write(iulog,*)'tfilt: t=',tm2(13,18),t3m1(13,18),t3(13,18)
!      write(iulog,*)'tfilt: water=',qm2(13,18,1),q3m1(13,18,1),q3(13,18,1)
!      write(iulog,*)'tfilt: cwat=',qm2(13,18,2),q3m1(13,18,2),q3(13,18,2)
!      write(iulog,*)'tfilt: vort=',vortm2(13,18),vortm1(13,18),vort(13,18)
!      write(iulog,*)'tfilt: div=',divm2(13,18),divm1(13,18),div(13,18)
!   end if

   om2eps = 1._r8 - 2._r8*eps

   if (nstep.ge.2) then
!$OMP PARALLEL DO PRIVATE (K, I, M)
      do k=1,plev
         do i=1,nlon
            u3m1(i,k) = om2eps*u3m1(i,k) + eps*um2(i,k) + eps*u3(i,k)
            v3m1(i,k) = om2eps*v3m1(i,k) + eps*vm2(i,k) + eps*v3(i,k)
            t3m1(i,k) = om2eps*t3m1(i,k) + eps*tm2(i,k) + eps*t3(i,k)
            q3m1(i,k,1) = om2eps*q3m1(i,k,1) + eps*qm2(i,k,1) + eps*q3(i,k,1)
            vortm1(i,k) = om2eps*vortm1(i,k) + eps*vortm2(i,k) + eps*vort(i,k)
            divm1(i,k) = om2eps*divm1(i,k) + eps*divm2(i,k) + eps*div(i,k)
         end do
         do m=2,pcnst
            if (cnst_get_type_byind(m) .eq. 'wet') then 
               do i=1,nlon
                  q3m1(i,k,m) = om2eps*q3m1(i,k,m) + eps*qm2(i,k,m) + eps*q3(i,k,m)
               end do
            endif 
         end do
         do m=2,pcnst
            if (cnst_get_type_byind(m) .eq. 'dry') then 
               do i=1,nlon  ! calculate numerator (timefiltered mass * pdeldry)
                  q3m1(i,k,m) = (om2eps*pdelm1dry(i,k)*q3m1(i,k,m) + &
                       eps*pdelm2dry(i,k)*qm2(i,k,m) + &
                       eps*pdeldry(i,k)*q3(i,k,m))
               end do !i
            endif  !dry
         end do !m
         do i=1,nlon     ! calculate time filtered value of pdeldry
               pdelm1dry(i,k) = om2eps*pdelm1dry(i,k) + &
                    eps*pdelm2dry(i,k) + eps*pdeldry(i,k)
         end do !i
         ! divide time filtered mass*pdeldry by timefiltered pdeldry
         do m=2,pcnst
            if (cnst_get_type_byind(m) == 'dry') then
               do i=1,nlon
                  q3m1(i,k,m) = q3m1(i,k,m)/pdelm1dry(i,k)
               end do !i
            endif ! dry
         end do !m

      end do
      do i=1,nlon
         psm1(i) = om2eps*psm1(i) + eps*psm2(i) + eps*ps(i)
      end do
   end if

   call plevs0 (nlon    ,plon   ,plev    ,psm1    ,pintm1  ,pmidm1  ,pdelm1)
!
! Compute time tendencies:comment out since currently not on h-t
!
!$OMP PARALLEL DO PRIVATE (K, I)
   do k=1,plev
      do i=1,nlon
         ttend(i,k) = (t3(i,k)-tm2(i,k))/ztodt
         utend(i,k) = (u3(i,k)-um2(i,k))/ztodt
         vtend(i,k) = (v3(i,k)-vm2(i,k))/ztodt
      end do
   end do

!$OMP PARALLEL DO PRIVATE (M, K, I)
   do m=1,pcnst
      do k=1,plev
         do i=1,nlon
            qtend(i,k,m) = (q3(i,k,m) - qm2(i,k,m))/ztodt
         end do
      end do
   end do

   do i=1,nlon
      pstend(i) = (ps(i) - psm2(i))/ztodt
   end do

!$OMP PARALLEL DO PRIVATE (M)
   do m=1,pcnst
      if ( cnst_cam_outfld(m) ) then
         call outfld (tendnam(m),qtend(1,1,m),plon,lat)
         call outfld (fixcnam(m),dqfx3(1,1,m),plon,lat)
         call outfld (hadvnam(m),hadv (1,1,m),plon,lat)
         call outfld (vadvnam(m),vadv (1,1,m),plon,lat)
      end if
   end do

! using do loop and select in order to enable functional parallelism with OpenMP
!$OMP PARALLEL DO PRIVATE (I)
   do i=1,4
      select case (i)
      case (1)
         call outfld ('UTEND   ',utend,plon,lat)
      case (2)
         call outfld ('VTEND   ',vtend,plon,lat)
      case (3)
         call outfld ('TTEND   ',ttend,plon,lat)
      case (4)
         call outfld ('LPSTEN  ',pstend,plon,lat)
      end select
   end do

   return
1000 format(' TIMEFILTER: WARNING: fixer for tracer ',i3,' exceeded ', &
      f8.5,' for ',i5,' points at k,lat = ',2i4, &
      ' Avg/Worst = ',1p2e10.2)

end subroutine  tfilt_massfixrun

end module tfilt_massfix
