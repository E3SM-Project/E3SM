
subroutine tfilt_massfix (ztodt   ,lat     ,u3m1    ,v3m1    ,t3m1    , &
                          q3      ,q3m1    ,ps      ,cwava   ,alpha   , &
                          etamid  ,qfcst   ,div     ,phis    ,omga    , &
                          dpsl    ,dpsm    ,nlon    ,um1     ,vm1     , &
                          tm1     ,qm1     ,psm1    ,beta    )
!-----------------------------------------------------------------------
!
! Purpose:
! Atmosphere and constituent mass fixer
!
! Author:  J. Olson
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  use pspect
  use commap
  use cam_history,  only: outfld
  use constituents, only: pcnst, qmin, cnst_cam_outfld, &
                          tottnam, tendnam, cnst_get_type_byind, fixcnam
  use physconst,    only: cpair, gravit
  use sld_control_mod, only : fixmas,eps

  implicit none

!
!------------------------------Arguments--------------------------------
!
  real(r8), intent(in)   :: ztodt                          ! timestep
  integer , intent(in)   :: lat                            ! latitude index
  real(r8), intent(in)   :: u3m1    (plon,plev)            ! u (n+1)
  real(r8), intent(in)   :: v3m1    (plon,plev)            ! v (n+1)
  real(r8), intent(inout):: t3m1    (plon,plev)            ! T (n+1)
  real(r8), intent(inout):: q3      (plon,plev,pcnst)      ! q + consts (time level n, after physics)
  real(r8), intent(inout):: q3m1    (plon,plev,pcnst)      ! q + consts (time level n+1)
  real(r8), intent(inout):: ps      (plon)                 ! Ps (n+1)
  real(r8), intent(in)   :: cwava                          ! weight for global integrals
  real(r8), intent(in)   :: alpha (pcnst)                  ! slt fixer coefficient
  real(r8), intent(in)   :: etamid(plev)                   ! vertical coords at midpoints 
  real(r8), intent(in)   :: qfcst (plon,plev,pcnst)        ! slt moisture forecast
  real(r8), intent(in)   :: div   (plon,plev)              ! divergence
  real(r8), intent(in)   :: phis  (plon)                   ! Geopotential field
  real(r8), intent(out)  :: omga  (plon,plev)              ! vertical motion
  real(r8), intent(in)   :: dpsl  (plon)                   ! long comp of grad ln(ps)
  real(r8), intent(in)   :: dpsm  (plon)                   ! lat comp  of grad ln(ps)
  integer , intent(in)   :: nlon                           ! number of longitudes
  real(r8), intent(in)   :: um1   (plon,plev)              ! u  (n)
  real(r8), intent(in)   :: vm1   (plon,plev)              ! v  (n)
  real(r8), intent(in)   :: tm1   (plon,plev)              ! temperature  (n)
  real(r8), intent(inout):: qm1   (plon,plev,pcnst)        ! q + consts input:   time level n, saved from last timestep;
                                                           !                     used to compute total tendency)
                                                           !            output:  time level n+1; saved for next time
  real(r8), intent(in)   :: psm1  (plon)                   ! Ps  (n)
  real(r8), intent(in)   :: beta                           ! energy fixer coefficient
!
!---------------------------Local workspace-----------------------------
!
  integer ifcnt                     ! Counter
  real(r8) tfix    (plon)           ! T correction
  real(r8) engycorr(plon,plev)      ! energy equivalent to T correction
  real(r8) rpmid (plon,plev)        ! 1./pmid
  real(r8) pdel  (plon,plev)        ! pdel(k)   = pint  (k+1)-pint  (k)
  real(r8) pint  (plon,plevp)       ! pressure at model interfaces (n  )
  real(r8) pmid  (plon,plev)        ! pressure at model levels (time n)
  real(r8) utend (plon,plev)        ! du/dt
  real(r8) vtend (plon,plev)        ! dv/dt
  real(r8) ttend (plon,plev)        ! dT/dt
  real(r8) qtend (plon,plev,pcnst)  ! dq/dt
  real(r8) pstend(plon)             ! d(ps)/dt
  real(r8) psl   (plon)             ! sea level pressure
  real(r8) corm                     ! fixer limit
  real(r8) wm                       ! accumulator 
  real(r8) absf                     ! absolute value of fixer
  real(r8) worst                    ! largest fixer contribution at each model level
  logical lfixlim                   ! flag to turn on fixer limiter

  real(r8) ta    (plon,plev,pcnst)  ! total advection of constituents
  real(r8) dqfx3 (plon,plev,pcnst)  ! q tendency due to mass adjustment
  real(r8) coslat                   ! cosine(latitude)
  real(r8) rcoslat(plon)            ! 1./cosine(latitude)

  integer i, k, m                   ! indices
!
!-----------------------------------------------------------------------
  coslat  = cos(clat(lat))
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
  rpmid(:nlon,:plev) = 1._r8/pmid(:nlon,:plev)
!
! Add temperature correction for energy conservation
!
  do k=1,plev
     do i=1,nlon
        engycorr(i,k) = (cpair/gravit)*beta*pdel(i,k)/ztodt
        t3m1    (i,k) = t3m1(i,k) + beta
     end do
  end do

  tfix(:nlon) = beta/ztodt
!
! Output Energy correction term
!
  call outfld ('ENGYCORR',engycorr ,plon   ,lat     )
  call outfld ('TFIX    ',tfix     ,plon   ,lat     )
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

     do k=1,plev

        do i=1,nlon
           dqfx3(i,k,m) = alpha(m)*etamid(k)*abs(qfcst(i,k,m) - q3(i,k,m))
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
!             write(iulog,1000) m,corm,ifcnt,k,lat,wm,worst

           endif
        endif

        do i=1,nlon
           dqfx3(i,k,m) = qfcst(i,k,m)*dqfx3(i,k,m)/ztodt
           q3m1(i,k,m) = qfcst(i,k,m) + ztodt*dqfx3(i,k,m)
           ta(i,k,m) = (q3m1(i,k,m) - q3(i,k,m))/ztodt
        end do
     end do
  end do

!
! stuff current water vapor into unused field so we can use it for 
! fixer next time step- This is required for dry mixing ratio conservation
!
  q3(:nlon,:,1) = q3m1(:nlon,:,1)

!
! Check for and correct invalid constituents
!
  call qneg3('TFILT_MASSFIX',lat   ,nlon    ,plon    ,plev    , &
             1, pcnst, qmin  ,q3m1 )
!
! Send slt tendencies to the history tape
!
  do m=1,pcnst
     if ( cnst_cam_outfld(m) ) then
        call outfld(tottnam(m),ta(1,1,m)   ,plon   ,lat     )
     end if
  end do
!
! Calculate vertical motion field
!
  call omcalc(rcoslat ,div     ,u3m1    ,v3m1    ,dpsl    , &
              dpsm    ,pmid    ,pdel    ,rpmid   ,pint(1,plevp), &
              omga    ,nlon    )

  call plevs0(nlon    ,plon   ,plev    ,ps      ,pint    ,pmid    ,pdel)
!
! Compute time tendencies
!
  do k=1,plev
     do i=1,nlon
        ttend(i,k) = (t3m1(i,k)-tm1(i,k))/ztodt
        utend(i,k) = (u3m1(i,k)-um1(i,k))/ztodt
        vtend(i,k) = (v3m1(i,k)-vm1(i,k))/ztodt
     end do
  end do

  do m=1,pcnst
     do k=1,plev
        do i=1,nlon
           qtend(i,k,m) = (q3m1(i,k,m) - qm1(i,k,m))/ztodt
           qm1  (i,k,m) =  q3m1(i,k,m)
        end do
     end do
  end do

  do i=1,nlon
     pstend(i) = (ps(i) - psm1(i))/ztodt
  end do

  do m=1,pcnst
     if ( cnst_cam_outfld(m) ) then
        call outfld (tendnam(m),qtend(1,1,m),plon,lat)
        call outfld (fixcnam(m),dqfx3(1,1,m),plon,lat)
     end if
  end do
  call outfld('UTEND   ',utend   ,plon   ,lat     )
  call outfld('VTEND   ',vtend   ,plon   ,lat     )
  call outfld('TTEND   ',ttend   ,plon   ,lat     )
  call outfld('LPSTEN  ',pstend  ,plon   ,lat     )

  call plevs0(nlon    ,plon   ,plev    ,ps      ,pint    ,pmid    ,pdel)

  return
1000 format(' TIMEFILTER: WARNING: fixer for tracer ',i3,' exceeded ', &
            f8.5,' for ',i5,' points at k,lat = ',2i4,' Avg/Worst = ',1p2e10.2)
end subroutine tfilt_massfix
