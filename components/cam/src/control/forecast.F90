subroutine forecast(lat, psm1, psm2,ps, &
                   u3, u3m1, u3m2, &
                   v3, v3m1, v3m2, &
                   t3, t3m1, t3m2, &
                   q3, q3m1, q3m2, ztodt, t2, &
                   fu, fv, qfcst,etamid, &
                   qminus, nlon)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Eularian forecast of t, u, and v.   Advection terms are also converted
! to flux form and integrated to check conservation
! 
! Author: 
! Original version:
!
!-----------------------------------------------------------------------

   use shr_kind_mod,   only: r8 => shr_kind_r8, i8 => shr_kind_i8
   use pmgrid
   use constituents,   only: pcnst, cnst_get_ind
   use time_manager,   only: is_first_step
   use pspect
   use physconst,      only: rair,cpair,gravit,rga
   use dycore,         only: dycore_is
   use cam_logfile,    only: iulog
   use scamMod
   use cam_history,    only: outfld
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
!
! Input arguments
!
   real(r8), intent(inout) :: t2(plev)         ! temp tendency
   real(r8), intent(inout) :: fu(plev)         ! u wind tendency
   real(r8), intent(inout) :: fv(plev)         ! v wind tendency
   real(r8), intent(inout) :: ps            ! surface pressure (time n)
   real(r8), intent(in) :: psm1          ! surface pressure (time n)
   real(r8), intent(in) :: psm2          ! surface pressure (time n-1)
   real(r8), intent(out) :: u3(plev)   ! u-wind (time n)
   real(r8), intent(in) :: u3m1(plev)   ! u-wind (time n)
   real(r8), intent(in) :: u3m2(plev) ! u-wind (time n-1)
   real(r8), intent(out) :: v3(plev)   ! u-wind (time n)
   real(r8), intent(in) :: v3m1(plev)   ! v-wind (time n)
   real(r8), intent(in) :: qminus(plon,plev,pcnst)
   real(r8), intent(in) :: v3m2(plev) ! v-wind (time n-1)
   real(r8), intent(out) :: t3(plev)   ! u-wind (time n)
   real(r8), intent(in) :: t3m1(plev)   ! temperature (time n)
   real(r8), intent(in) :: t3m2(plev)   ! temperature (time n)
   real(r8), intent(inout) :: q3(plev,pcnst)   ! constituent conc(tim
   real(r8), intent(inout) :: q3m1(plev,pcnst)   ! constituent conc(tim
   real(r8), intent(inout) :: q3m2(plev,pcnst)   ! constituent conc(time n: h2o first)
   real(r8), intent(in) :: etamid(plev)       ! vertical coords at midpoints
   real(r8), intent(inout) :: qfcst(plon,plev,pcnst)

   real(r8), intent(in) :: ztodt                       ! twice time step unless nstep=0
   integer lat               ! latitude index for S->N storage
   integer nlon
!
!---------------------------Local workspace-----------------------------
!
   integer jcen                ! lat index (extended grid) of forecast
   integer iter                ! number of iterations for
   integer itermx  ! number of iterations to be used in departure
!                     ! point calculation for nstep = 0 and 1
   integer itermn  ! number of iterations to be used in departure
!                     ! point calculation for nstep > 1
   parameter(itermx=4,itermn=1)
   real(r8) pmidm1(plev)  ! pressure at model levels (time n)
   real(r8) pintm1(plevp) ! pressure at model interfaces (n  )
   real(r8) pdelm1(plev)  ! pdel(k)   = pint  (k+1)-pint  (k)
   real(r8) pmidm2(plev)  ! pressure at model levels (time n)
   real(r8) pintm2(plevp) ! pressure at model interfaces (n  )
   real(r8) pdelm2(plev)  ! pdel(k)   = pint  (k+1)-pint  (k)
   real(r8) pmidm1f(plev)  ! pressure at model levels (time n)
   real(r8) pintm1f(plevp) ! pressure at model interfaces (n  )
   real(r8) pdelm1f(plev)  ! pdel(k)   = pint  (k+1)-pint  (k)
   real(r8) pdelb(plon,plev)  ! pressure diff bet intfcs (press defined using the "B" part 
   real(r8) pdela(plon,plev)
   real(r8) weight,fac
   real(r8) psfcst
   real(r8) tfcst(plev)
   real(r8) ufcst(plev)
   real(r8) vfcst(plev)
   real(r8) tdwdp(plev)
   real(r8) vdwdp(plev)
   real(r8) udwdp(plev)
   real(r8) qdwdp(plev,pcnst)
   real(r8) wfldint(plevp)     ! midpoint values of eta (a+b)
   real(r8) tfmod(plev)
   real(r8) ufmod(plev)
   real(r8) vfmod(plev)
   real(r8) qfmod(plev,pcnst)
   real(r8) alpha(pcnst)
   real(r8) sum
   real(r8) d_qdw
   real(r8) d_qdwdp(plev)
   real(r8) d_dqfx(plev)
   real(r8) d_qdv(plev)
   real(r8) d_qtd
   real(r8) d_qtv
   real(r8) d_qdvt
   real(r8) d_dqfxt
   real(r8) dqv(plev,pcnst)   ! constituent diffusion tendency
   save dqv
   real(r8) qphys(plev,pcnst)   ! constituent diffusion tendency

   real(r8) dqfx3m1(plev,pcnst) ! q tendency due to mass adjustment
!
   real(r8) qmassb(pcnst)     ! constituent mass integral before advection
   real(r8) hwava (pcnst)     ! temporary variable for mass fixer
   real(r8) ptb               ! potential temperature before advection
   real(r8) ptf               ! potential temperature after advection
   real(r8) dotproda           ! dot product
   real(r8) dotprodb           ! dot product
   integer i,k,m           ! longitude, level, constituent indices
!
!     Below are Variables Used in the Advection Diagnostics
!
   integer mplot
   parameter ( mplot = 1 ) ! The tracer for which all Advection Diagnostic
!                               ! are to be plotted 1 = q, 2 = tr01 etc...
!     
!     dummy arguments for outfld calls  in SCM
   integer  dummy

!
!  variables for relaxation addition
!
   real(r8) dist
   real(r8) denom
   real(r8) rtau(plev)
   real(r8) relaxt(plev)
   real(r8) relaxq(plev)
   logical relax
!
!  diagnostic variables for estimating vertical advection terms
!
   real(r8) tvadv(plev)       !estimate of vertical advection on T
   real(r8) qvadv(plev,pcnst)!estimate of vertical advection on q
   real(r8) qvadv1(plev,pcnst)!estimate of vertical advection on q
!
!  diagnostic variables for maintaining n-1 values of observed T and q
!
   real(r8) tobsm1(plev)
   real(r8) qobsm1(plev)
   save qobsm1, tobsm1

   real(r8) q3forecast,t3forecast
   real(r8) forecastdiff,bestforecastdiff
   real(r8) qmassf
   integer  j,icldliq,icldice

   l_conv  = .true.       ! .f. doesn't use divT and divq
   l_divtr = .false.      ! .t. includes some div of condensate
!     
   if(use_iop) then
      l_uvadvect = .false.
      l_uvphys   = .false.
   else
      l_uvadvect = .false.
      l_uvphys   = .false.
   end if

!
   call plevs0(nlon    ,plon   ,plev    ,psm1   ,pintm1  ,pmidm1 ,pdelm1)
   call plevs0(nlon    ,plon   ,plev    ,psm2   ,pintm2  ,pmidm2 ,pdelm2)
!
! Build interface vector for the specified omega profile
! (weighted average in pressure of specified level values)
!
   wfldint(1) = 0.0_r8

   do k=2,plev
      weight = (pintm1(k) - pmidm1(k-1))/(pmidm1(k) - pmidm1(k-1))
      wfldint(k) = (1.0_r8 - weight)*wfld(k-1) + weight*wfld(k)
   end do

   wfldint(plevp) = 0.0_r8

   if (use_3dfrc .and. use_iop) then

!  Complete a very simple forecast using supplied 3-dimensional forcing
!  by the large scale.  Obviates the need for any kind of vertical 
!  advection calculation.  Skip to diagnostic estimates of vertical term.
      i=1
      do k=1,plev
         tfcst(k) = t3m2(k) + ztodt*t2(k) + ztodt*divt3d(k)
      end do
      do m=1,pcnst
         do k=1,plev
            qfcst(1,k,m) = qminus(1,k,m) +  divq3d(k,m)*ztodt
         end do
      enddo

      go to 1000

   end if

!
!  provide an eulerian forecast.  First check to ensure that 2d forcing
!  is available.  If not and it is required for the forecast then calculate
!  it as a residule of the 3d forcing.  The gui will guarentte that the
!  appropriate 2d and/or 3d forcing is available so there is no need to
!  place software checks here to guard agains missing data.
!

      if((.not. (have_divt .and. have_divq)) .and. use_iop) then
!
!---ESTIMATE VERTICAL ADVECTION TENDENCY FOR T,q TO EVALUATE---
!---      HORIZONTAL ADVECTION COMPONENTS AS RESIDUALS      ---
!   using eulerian form for evaluating advection ... close enough!
!
         do k=2,plev-1
            fac = 1.0_r8/(2.0_r8*pdelm1(k))
            tvadv(k) =  - fac*(wfldint(k+1)*(t3m1(k+1) - t3m1(k)) &
                + wfldint(k)*(t3m1(k) - t3m1(k-1))) &
                + wfld(k)*t3m1(k)*rair/(cpair*pmidm1(k))
            do m=1,pcnst
               qvadv(k,m) =  (qfcst(1,k,m)-qminus(1,k,m))/ztodt
            end do
         end do
!     
!   - top and bottom levels next -
!
         k = 1
         fac = 1.0_r8/(2.0_r8*pdelm1(k))
         tvadv(k) = - fac*(wfldint(k+1)*(t3m1(k+1) - t3m1(k))) &
                      + wfld(k)*t3m1(k)*rair/(cpair*pmidm1(k))
         do m=1,pcnst
            qvadv(k,m) =   (qfcst(1,k,m)-qminus(1,k,m))/ztodt
         end do
!     
         k = plev
         fac = 1.0_r8/(2.0_r8*pdelm1(k))
         tvadv(k) = - fac*(wfldint(k)*(t3m1(k) - t3m1(k-1))) &
                      + wfld(k)*t3m1(k)*rair/(cpair*pmidm1(k))
         do m=1,pcnst
            qvadv(k,m) = (qfcst(1,k,m)-qminus(1,k,m))/ztodt
         end do
!
!     here's where the residuals are evaluated
!
         do k=1,plev
            divt(k) = divt3d(k) - tvadv(k)
            do m=1,pcnst
               divq(k,m) = divq3d(k,m) - qvadv(k,m)
            end do
         end do
!
      end if
!
! TIME FOR VERTICAL ADVECTION STEP
!
!
!  Eularian forecast for u,v and t
!

   if (dycore_is('EUL')) then 

     do k=2,plev-1
       fac = ztodt/(2.0_r8*pdelm1(k))
       tfcst(k) = t3m2(k) &
           - fac*(wfldint(k+1)*(t3m1(k+1) - t3m1(k)) &
           + wfldint(k)*(t3m1(k) - t3m1(k-1)))
       vfcst(k) = v3m2(k) &
           - fac*(wfldint(k+1)*(v3m1(k+1) - v3m1(k)) &
           + wfldint(k)*(v3m1(k) - v3m1(k-1)))
       ufcst(k) = u3m2(k) &
           - fac*(wfldint(k+1)*(u3m1(k+1) - u3m1(k)) &
           + wfldint(k)*(u3m1(k) - u3m1(k-1)))

     end do

!     
!     - top and bottom levels next -
!     

     k = 1
     fac = ztodt/(2.0_r8*pdelm1(k))
     tfcst(k) = t3m2(k) - fac*(wfldint(k+1)*(t3m1(k+1) - t3m1(k)))
     vfcst(k) = v3m2(k) - fac*(wfldint(k+1)*(v3m1(k+1) - v3m1(k)))
     ufcst(k) = u3m2(k) - fac*(wfldint(k+1)*(u3m1(k+1) - u3m1(k)))

     k = plev
     fac = ztodt/(2.0_r8*pdelm1(plev))
     tfcst(k) = t3m2(k) - fac*(wfldint(k)*(t3m1(k) - t3m1(k-1)))
     vfcst(k) = v3m2(k) - fac*(wfldint(k)*(v3m1(k) - v3m1(k-1)))
     ufcst(k) = u3m2(k) - fac*(wfldint(k)*(u3m1(k) - u3m1(k-1)))

!
!  SLT is used for constituents only
!  so that a centered approximation is used for T, U and V, and Q
!  check to see if we should be using a forward approximation for 
!  constituents
     do k=1,plev
       tdwdp(k) = t3m1(k)*(wfldint(k+1)-wfldint(k))/pdelm1(k)
       udwdp(k) = u3m1(k)*(wfldint(k+1)-wfldint(k))/pdelm1(k)
       vdwdp(k) = v3m1(k)*(wfldint(k+1)-wfldint(k))/pdelm1(k)
       do m=1,pcnst
         qdwdp(k,m) = qminus(1,k,m)*(wfldint(k+1)-wfldint(k))/pdelm2(k)
       end do
     end do
   
   else if (dycore_is('SE')) then
   
     tfcst(:) = t3m2(:)
     qfcst(1,:,:) = q3m2(:,:)
     ufcst(:) = u3m2(:)
     vfcst(:) = v3m2(:)
   
   endif

if (.not.use_iop) then
!
!
!  Modify advection forecast to properly enforce conservation
!  These terms are removed after conservation procedure has been applied
!
   do k=1,plev
      tfmod(k)      = - ztodt*tdwdp(k) + ztodt*wfld(k)*t3m1(k)*rair/(cpair*pmidm1(k))
      vfmod(k)      = - ztodt*vdwdp(k)
      ufmod(k)      = - ztodt*udwdp(k)
      do m=1,pcnst
         qfmod(k,m) = - ztodt*qdwdp(k,m)
      end do
!
      tfcst(k)   = tfcst(k)   + tfmod(k)
      vfcst(k)   = vfcst(k)   + vfmod(k)
      ufcst(k)   = ufcst(k)   + ufmod(k)
      do m=1,pcnst
         qfcst(1,k,m) = qfcst(1,k,m) + qfmod(k,m)
      end do
   end do

   call plevs0(nlon    ,plon   ,plev    ,psm1   ,pintm1f  ,pmidm1f ,pdelm1f)

!
! Place 1st set of Jims Diagnostics Here if desired
!
   if (l_diag) then !=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
!c
!        enthalpy conservation calculation
!
      ptb = 0.0_r8
      ptf = 0.0_r8
      do k=1,plev
         ptb = ptb + (t3m1(k)*((100000.0_r8/pmidm1(k))**.28571_r8)) &
            *(pdelm1(k)/(psm1 - pintm1(1)))
         ptf = ptf + (tfcst(k)*((100000.0_r8/pmidm1f(k))**.28571_r8)) &
            *(pdelm1f(k)/(psfcst - pintm1f(1)))
      end do
!
!        water vapor conservation
!
      qmassf = 0.0_r8
      do k=1,plev
         qmassf = qmassf + pdelm1f(k)*qfcst(1,k,1)/gravit
      end do
!
!        print t & q forecast information (before/after conservative advection)
!
      write(iulog,986)
986   format (' conservative advection characteristics')
      write(iulog,987) (t3m1(k), tfcst(k), &
         (tfcst(k)-t3m1(k)), &
         q3m1(k,1), qfcst(1,k,1), &
         (qfcst(1,k,1)-q3m1(k,1)), &
         864.0_r8*wfld(k), 0.01_r8*pdelm1(k), k=1,plev)
987   format (1x, 0p, 3f11.4, 3p, 3f11.4, 0p, 2f11.4)
!
!        print water vapor correction required for conservation
!
!
      write(iulog,1105) qmassb(1), qmassf, ptb, ptf
1105  format (' qmassb, qmassf; ptb, ptf =>',1p,2e12.3,'; ',3x,2e14.5)
!
   endif !=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

!
!
!  *** Remove flux correction term from advection forecast (after "fixer") ***
!  This is something the user should carefully consider, since in the
!  absence of specified or calculated horizontal advection tendencies
!  the advective form of the vertical transport term results in an
!  anomolous or implied source/sink for the respective equations
!
   do k=1,plev
      tfcst(k)   = tfcst(k)   - tfmod(k)
      vfcst(k)   = vfcst(k)   - vfmod(k)
      ufcst(k)   = ufcst(k)   - ufmod(k)
      do m=1,pcnst
         qfcst(1,k,m) = qfcst(1,k,m) - qfmod(k,m)
      end do
   end do

!
!     *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
!
! Place 2nd set of Jims Diagnostics here
!
   if (l_diag) then   ! check similar stuff as before w/o flux formalism
!
!        enthalpy conservation
!
      ptf = 0.0_r8
      do k=1,plev
         ptf = ptf + (tfcst(k)*((100000.0_r8/pmidm1f(k))**.28571_r8)) &
            *(pdelm1f(k)/(psfcst - pintm1f(1)))
      end do
!
!        water vapor conservation
!
      qmassf = 0.0_r8
      do k=1,plev
         qmassf = qmassf + pdelm1f(k)*qfcst(1,k,1)/gravit
      end do
!
!        print t & q forecast information (before/after advection)
!
      write(iulog,985)
985   format (' non-conservative advection characteristics')
      write(iulog,987) (t3m1(k), tfcst(k), &
         (tfcst(k)-t3m1(k)), &
         q3m1(k,1), qfcst(1,k,1), &
         (qfcst(1,k,1)-q3m1(k,1)), &
         864.0_r8*wfld(k), 0.01_r8*pdelm1(k), k=1,plev)
!
      write(iulog,1105) qmassb(1), qmassf, ptb, ptf
!
   endif                     !=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
! End of 2nd set
end if

!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!  *** MAKE THE FORECAST ***
!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
! Include energy conversion term in thermodynamic energy equation 
! Include all physics tendency information passed up from linems
! Include flux divergence information for T and q if available
! Code assumes that the flux divergence info is in tendency units
! -- Update temperature
! -- Update moisture
! -- Update momentum
!
!     Zero Convergence terms if l_conv is false
!
   if (.not.l_conv.or..not.use_iop) then
      do k=1,plev
         divt(k)   = 0.0_r8
         divq(k,1) = 0.0_r8
      enddo
   endif

!
!  Note: if including relaxation as part of the forward forecast step
!        add it here to t2 and dqv
!

   do k=1,plev
     tfcst(k) = tfcst(k) + ztodt*wfld(k)*t3m1(k)*rair/(cpair*pmidm1(k)) &
         + ztodt*(t2(k) + divt(k))
     do m=1,pcnst
       qfcst(1,k,m) = qfcst(1,k,m) + ztodt*divq(k,m)
     end do
   enddo

!     
!---ESTIMATE VERTICAL ADVECTION TENDENCY FOR T,q (DIAGNOSTIC)------
!   using eulerian form for evaluating advection (can actually
!   do this more accurately as residual before forecast step, 
!   but won't work if applying "revealed forcing" to model which
!   includes both horizontal and vertical large-scale forcing terms.
!   This is close enough for now!
!

1000 continue
   do k=2,plev-1
      fac = 1.0_r8/(2.0_r8*pdelm1(k))
      tvadv(k) =  - fac*(wfldint(k+1)*(t3m1(k+1) - t3m1(k)) &
         + wfldint(k)*(t3m1(k) - t3m1(k-1))) &
         + wfld(k)*t3m1(k)*rair/(cpair*pmidm1(k))
      do m=1,pcnst
         qvadv(k,m) =  - fac*(wfldint(k+1)*(q3m1(k+1,m) - q3m1(k,m)) &
            + wfldint(k)*(q3m1(k,m) - q3m1(k-1,m)))
      end do
   end do
!
!   - top and bottom levels next -
!
   k = 1
   fac = 1.0_r8/(2.0_r8*pdelm1(k))
   tvadv(k) = - fac*(wfldint(k+1)*(t3m1(k+1) - t3m1(k))) &
      + wfld(k)*t3m1(k)*rair/(cpair*pmidm1(k))
   do m=1,pcnst
      qvadv(k,m) = - fac*(wfldint(k+1)*(q3m1(k+1,m) - q3m1(k,m)))
   end do
!     
   k = plev
   fac = 1.0_r8/(2.0_r8*pdelm1(k))
   tvadv(k) = - fac*(wfldint(k)*(t3m1(k) - t3m1(k-1))) &
      + wfld(k)*t3m1(k)*rair/(cpair*pmidm1(k))
   do m=1,pcnst
      qvadv(k,m) = - fac*(wfldint(k)*(q3m1(k,m) - q3m1(k-1,m)))
   end do

!!$   call outfld('TVADV'   ,tvadv,plon,lat)
!!$   call outfld('QVADV'   ,qvadv,plon,lat)
!
!---end diagnostic estimates of vertical advection term----------
!
!     Using New Logicals for controlling changes to u,v
!
   if (.not.l_uvphys) then
      do k=1,plev
         fu(k) = 0.0_r8
         fv(k) = 0.0_r8
      enddo
   endif
!
   if(.not.l_uvadvect) then
      if (use_iop .and. have_v .and. have_u) then
         do k=1,plev
            ufcst(k) = uobs(k)
            vfcst(k) = vobs(k)
         enddo
!
      else
!
         do k=1,plev
            ufcst(k) = u3m2(k)
            vfcst(k) = v3m2(k)
         enddo
!
      endif      ! from  if (use_iop .and. have_v .and. have_u) 
!      
   else
!
      do k=1,plev
         ufcst(k) = ufcst(k) + ztodt*(fu(k) + divu(k))
         vfcst(k) = vfcst(k) + ztodt*(fv(k) + divv(k))
      enddo
   endif

!
! Copy fields from SLT/Eulerian forecast location to appropriate location in q3
!
   q3(:,:pcnst)=qfcst(1,:,:pcnst)
   t3(:)=tfcst(:)
   u3(:)=ufcst(:)
   v3(:)=vfcst(:)

   if (scm_relaxation) then
!
!    THIS IS WHERE WE RELAX THE SOLUTION IF REQUESTED
!    The relaxation can be thought of as a part of the "adjustment" physics
!
!    Another way to do this is to estimate the error at t3m2, q3m2 and
!    include it in the prediction equations (e.g., sum it with the t2
!    term from the tendency physics).  This is numerically stable, but
!    can not provide a "hard relaxation" because the adjustment physics 
!    then operates on the forecast value.  In order to minimize changes
!    to the code we move the outfld calls for the relaxed variables
!    (in this case T and q) from linemsbc into this routine after the
!    relaxation terms have been applied.
!
      do k=1,plev
         relaxt(k) = 0.0_r8
         relaxq(k) = 0.0_r8
      end do
!
      do k=1,plev
           
        if (pmidm1(k) .le. scm_relaxation_low*100._r8 .and. &
          pmidm1(k) .ge. scm_relaxation_high*100._r8) then

          rtau(k)   = 10800._r8          ! 3-hr adj. time scale
          rtau(k)   = max(ztodt,rtau(k))
          relaxt(k) = -(t3(k)   - tobs(k))/rtau(k)
          relaxq(k) = -(q3(k,1) - qobs(k))/rtau(k)
!
          t3(k)     = t3(k)   + relaxt(k)*ztodt
          q3(k,1)   = q3(k,1) + relaxq(k)*ztodt
        
        endif

      end do
!
         call outfld('TRELAX',relaxt,plon,lat )
         call outfld('QRELAX',relaxq,plon,lat )
         call outfld('TAURELAX',rtau,plon,lat )
!      end if
   end if

!     
!  evaluate the difference in state information from observed
!
   do k = 1, plev
      tdiff(k) = t3(k)   - tobs(k)
      qdiff(k) = q3(k,1) - qobs(k)
      udiff(k) = u3(k)   - uobs(k)
      vdiff(k) = v3(k)   - vobs(k)
   end do

!
! Copy observations into time n-1 storage (has diagnostics utility only)
!
   tobsm1(:)=tobs(:)
   qobsm1(:)=qobs(:)
!
!===============================================================
!
!  outfld calls moved from linemsbc

   call outfld('TOBS',tobs,plon,lat)
   call outfld('QOBS',qobs,plon,lat)
   call outfld('TDIFF',tdiff,plon,lat)
   call outfld('QDIFF',qdiff,plon,lat)
   if( use_iop ) then
      call outfld('DIVQ',divq,plon,lat)
      call outfld('DIVT',divt,plon,lat)
      call outfld('DIVQ3D',divq3d,plon,lat)
      call outfld('DIVT3D',divt3d,plon,lat)
!!$      call outfld('DIVU',divu,plon,lat)
!!$      call outfld('DIVV',divv,plon,lat)
      call outfld('PRECOBS',precobs,plon,lat )
      call outfld('LHFLXOBS',lhflxobs,plon,lat )
      call outfld('SHFLXOBS',shflxobs,plon,lat )
!!$      call outfld('Q1OBS',q1obs,plon,lat )
!!$      call outfld('Q2OBS',q2obs,plon,lat )
   endif

!
! Diagnose pressure arrays needed by DIFCOR
!
end subroutine forecast


!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------


