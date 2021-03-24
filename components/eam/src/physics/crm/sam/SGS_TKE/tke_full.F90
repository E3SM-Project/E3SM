module tke_full_mod

  use shear_prod2D_mod
  use shear_prod3D_mod
  use sat_mod

  implicit none

contains

subroutine tke_full(ncrms,dimx1_d, dimx2_d, dimy1_d, dimy2_d,   &
                    grdf_x, grdf_y, grdf_z, dosmagor,     &
                    tkesbdiss, tkesbshear, tkesbbuoy,     &
                    tke, tk, tkh)
  !-----------------------------------------------------------------------
  ! Purpose: solve the TKE equation
  !-----------------------------------------------------------------------

  use vars
  use params
  use openacc_utils
  implicit none
  integer, intent(in) :: ncrms
  !-----------------------------------------------------------------------
  !!! Interface Arguments
  integer       , intent(in)                 :: dimx1_d     ! scalar dimension parameter
  integer       , intent(in)                 :: dimx2_d     ! scalar dimension parameter
  integer       , intent(in)                 :: dimy1_d     ! scalar dimension parameter
  integer       , intent(in)                 :: dimy2_d     ! scalar dimension parameter
  real(crm_rknd), intent(in), dimension(ncrms,nzm) :: grdf_x      ! grid length in x direction
  real(crm_rknd), intent(in), dimension(ncrms,nzm) :: grdf_y      ! grid length in y direction
  real(crm_rknd), intent(in), dimension(ncrms,nzm) :: grdf_z      ! grid length in z direction
  logical       , intent(in)                 :: dosmagor    ! flag for diagnostic smagorinsky scheme

  real(crm_rknd), intent(out), dimension(ncrms,nz) :: tkesbdiss   ! TKE dissipation
  real(crm_rknd), intent(out), dimension(ncrms,nz) :: tkesbshear  ! TKE production by shear
  real(crm_rknd), intent(out), dimension(ncrms,nz) :: tkesbbuoy   ! TKE production by buoyancy

  real(crm_rknd), intent(  out) :: tke(ncrms, dimx1_s:dimx2_s , dimy1_s:dimy2_s , nzm )   ! SGS TKE
  real(crm_rknd), intent(inout) :: tk (ncrms, dimx1_d:dimx2_d , dimy1_d:dimy2_d , nzm )   ! SGS eddy viscosity
  real(crm_rknd), intent(inout) :: tkh(ncrms, dimx1_d:dimx2_d , dimy1_d:dimy2_d , nzm )   ! SGS eddy conductivity

  !-----------------------------------------------------------------------
  !!! Local Variables
  real(crm_rknd), allocatable :: def2          (:,:,:,:)
  real(crm_rknd), allocatable :: buoy_sgs_vert (:,:,:,:)
  real(crm_rknd), allocatable :: a_prod_bu_vert(:,:,:,:)
  real(crm_rknd) :: grd           !
  real(crm_rknd) :: betdz         !
  real(crm_rknd) :: Ck            !
  real(crm_rknd) :: Ce            !
  real(crm_rknd) :: Ces           !
  real(crm_rknd) :: Ce1           !
  real(crm_rknd) :: Ce2           !
  real(crm_rknd) :: smix          !
  real(crm_rknd) :: Pr            ! Prandtl number
  real(crm_rknd) :: Cee           !
  real(crm_rknd) :: Cs            !
  real(crm_rknd) :: buoy_sgs      !
  real(crm_rknd) :: ratio         !
  real(crm_rknd) :: a_prod_sh     ! shear production of TKE
  real(crm_rknd) :: a_prod_bu     ! buoyant production of TKE
  real(crm_rknd) :: a_diss        ! TKE dissipation
  real(crm_rknd) :: lstarn        !
  real(crm_rknd) :: lstarp        !
  real(crm_rknd) :: bbb           !
  real(crm_rknd) :: omn           !
  real(crm_rknd) :: omp           !
  real(crm_rknd) :: qsatt         !
  real(crm_rknd) :: dqsat         !
  real(crm_rknd) :: cx            ! correction factor for eddy visc CFL criteria
  real(crm_rknd) :: cy            ! correction factor for eddy visc CFL criteria
  real(crm_rknd) :: cz            ! correction factor for eddy visc CFL criteria
  real(crm_rknd) :: tkmax         ! Maximum TKE (CFL limiter)

  integer :: i,j,k,icrm
  integer :: kc      ! = k+1
  integer :: kb      ! = k-1

  real(crm_rknd) :: tabs_interface    ! tabs interpolated to interfaces
  real(crm_rknd) :: qp_interface      ! qp   interpolated to interfaces
  real(crm_rknd) :: qtot_interface    ! qtot interpolated to interfaces
  real(crm_rknd) :: qsat_check        ! used to check for cloud
  real(crm_rknd) :: qctot             ! total condensate

  real(crm_rknd) :: tk_min_value      ! min value for eddy viscosity (TK)
  real(crm_rknd) :: tk_min_depth      ! near-surface depth to apply tk_min (meters)
  real(crm_rknd) :: tmp

  allocate( def2          (ncrms,nx,ny,nzm  ) )
  allocate( buoy_sgs_vert (ncrms,nx,ny,0:nzm) )
  allocate( a_prod_bu_vert(ncrms,nx,ny,0:nzm) )
  call prefetch( def2           )
  call prefetch( buoy_sgs_vert  )
  call prefetch( a_prod_bu_vert )

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------

  tk_min_value = 0.05D0
  tk_min_depth = 500.D0

  Cs  = 0.15D0
  Ck  = 0.1D0
  Ce  = Ck**3.D0/Cs**4.D0
  Ces = Ce/0.7D0*3.0D0
  Pr  = 1.D0

  if(RUN3D) then
    call shear_prod3D(ncrms,def2)
  else
    call shear_prod2D(ncrms,def2)
  endif

  !!! initialize surface and top buoyancy flux to zero
  !$acc parallel loop collapse(3) async(asyncid)
  do j = 1 , ny
    do i = 1 , nx
      do icrm = 1 , ncrms
        a_prod_bu_vert(icrm,i,j,0) = 0.D0
        buoy_sgs_vert(icrm,i,j,0) = 0.D0
        a_prod_bu_vert(icrm,i,j,nzm) = 0.D0
        buoy_sgs_vert(icrm,i,j,nzm) = 0.D0
      enddo
    enddo
  enddo

  !-----------------------------------------------------------------------
  ! compute SGS buoyancy flux at w-levels and SGS quantities in interior of domain
  !-----------------------------------------------------------------------
  !!! compute subgrid buoyancy flux assuming clear conditions
  !!! we will over-write this later if conditions are cloudy
  !$acc parallel loop collapse(4) async(asyncid)
  do k = 1,nzm-1
    do j = 1,ny
      do i = 1,nx
        do icrm = 1 , ncrms
          if (k.lt.nzm) then
            kc = k+1
            kb = k
          else
            kc = k
            kb = k-1
          end if
          !!! first compute subgrid buoyancy flux at interface above this level.
          !!! average betdz to w-levels
          betdz = 0.5D0*(bet(icrm,kc)+bet(icrm,kb))/dz(icrm)/adzw(icrm,k+1)

          !!! compute temperature of mixture between two grid levels if all cloud
          !!!   were evaporated and sublimated
          tabs_interface = &
               0.5D0*( tabs(icrm,i,j,kc) + fac_cond*qcl(icrm,i,j,kc) + fac_sub*qci(icrm,i,j,kc) &
                     + tabs(icrm,i,j,kb) + fac_cond*qcl(icrm,i,j,kb) + fac_sub*qci(icrm,i,j,kb) )

          !!! similarly for water vapor if all cloud evaporated/sublimated
          qtot_interface = &
               0.5D0*( qv(icrm,i,j,kc) + qcl(icrm,i,j,kc) + qci(icrm,i,j,kc) &
                     + qv(icrm,i,j,kb) + qcl(icrm,i,j,kb) + qci(icrm,i,j,kb) )

          qp_interface = 0.5D0*( qpl(icrm,i,j,kc) + qpi(icrm,i,j,kc) + qpl(icrm,i,j,kb) + qpi(icrm,i,j,kb) )

          bbb = 1.D0+epsv*qtot_interface - qp_interface
          buoy_sgs=betdz*( bbb*(t(icrm,i,j,kc)-t(icrm,i,j,kb)) &
               +epsv*tabs_interface* &
               (qv(icrm,i,j,kc)+qcl(icrm,i,j,kc)+qci(icrm,i,j,kc)-qv(icrm,i,j,kb)-qcl(icrm,i,j,kb)-qci(icrm,i,j,kb)) &
               +(bbb*fac_cond-tabs_interface)*(qpl(icrm,i,j,kc)-qpl(icrm,i,j,kb)) &
               +(bbb*fac_sub -tabs_interface)*(qpi(icrm,i,j,kc)-qpi(icrm,i,j,kb)) )

          buoy_sgs_vert(icrm,i,j,k) = buoy_sgs
          a_prod_bu_vert(icrm,i,j,k) = -0.5D0*(tkh(icrm,i,j,kc)+tkh(icrm,i,j,kb)+0.002D0)*buoy_sgs

          !-----------------------------------------------------------------------
          ! now go back and check for cloud
          !-----------------------------------------------------------------------

          !!! if there's any cloud in the grid cells above or below, check to see if
          !!! the mixture between the two levels is also cloudy
          qctot = qcl(icrm,i,j,kc)+qci(icrm,i,j,kc)+qcl(icrm,i,j,kb)+qci(icrm,i,j,kb)
          if(qctot .gt. 0.D0) then

            !!! figure out the fraction of condensate that's liquid
            omn = (qcl(icrm,i,j,kc)+qcl(icrm,i,j,kb))/(qctot+1.D-20)

            !!! compute temperature of mixture between two grid levels
            !!! if all cloud were evaporated and sublimated
            tabs_interface = &
                 0.5D0*( tabs(icrm,i,j,kc) + fac_cond*qcl(icrm,i,j,kc) + fac_sub*qci(icrm,i,j,kc) &
                       + tabs(icrm,i,j,kb) + fac_cond*qcl(icrm,i,j,kb) + fac_sub*qci(icrm,i,j,kb) )

            !!! similarly for total water (vapor + cloud) mixing ratio
            qtot_interface = &
                 0.5D0*( qv(icrm,i,j,kc) + qcl(icrm,i,j,kc) + qci(icrm,i,j,kc) &
                       + qv(icrm,i,j,kb) + qcl(icrm,i,j,kb) + qci(icrm,i,j,kb) )

            !!! compute saturation mixing ratio at this temperature
            qsat_check =       omn *qsatw_crm(tabs_interface,presi(icrm,k+1)) &
                        +(1.D0-omn)*qsati_crm(tabs_interface,presi(icrm,k+1))

            !!! check to see if the total water exceeds this saturation mixing ratio.
            !!! if so, apply the cloudy relations for subgrid buoyancy flux
            if(qtot_interface.gt.qsat_check) then
              !!! apply cloudy relations for buoyancy flux, use the liquid-ice breakdown computed above.
              lstarn = fac_cond+(1.D0-omn)*fac_fus
              !!! use the average values of T from the two levels to compute qsat, dqsat
              !!! and the multipliers for the subgrid buoyancy fluxes.  Note that the
              !!! interface is halfway between neighboring levels, so that the potential
              !!! energy cancels out.  This is approximate and neglects the effects of
              !!! evaporation/condensation with mixing.  Hopefully good enough.
              tabs_interface = 0.5D0*( tabs(icrm,i,j,kc) + tabs(icrm,i,j,kb) )

              qp_interface = 0.5D0*( qpl(icrm,i,j,kc) + qpi(icrm,i,j,kc) + qpl(icrm,i,j,kb) + qpi(icrm,i,j,kb) )

              dqsat =      omn *dtqsatw_crm(tabs_interface,presi(icrm,k+1)) + &
                     (1.D0-omn)*dtqsati_crm(tabs_interface,presi(icrm,k+1))
              qsatt =        omn *qsatw_crm(tabs_interface,presi(icrm,k+1)) + &
                       (1.D0-omn)*qsati_crm(tabs_interface,presi(icrm,k+1))

              !!! condensate loading term
              bbb = 1.D0 + epsv*qsatt &
                   + qsatt - qtot_interface - qp_interface &
                   +1.61D0*tabs_interface*dqsat
              bbb = bbb / (1.D0+lstarn*dqsat)

              buoy_sgs = betdz*(bbb*(t(icrm,i,j,kc)-t(icrm,i,j,kb)) &
                   +(bbb*lstarn - (1.D0+lstarn*dqsat)*tabs_interface)* &
                   (qv(icrm,i,j,kc)+qcl(icrm,i,j,kc)+qci(icrm,i,j,kc)-qv(icrm,i,j,kb)-qcl(icrm,i,j,kb)-qci(icrm,i,j,kb)) &
                   + ( bbb*fac_cond-(1.D0+fac_cond*dqsat)*tabs(icrm,i,j,k) ) * ( qpl(icrm,i,j,kc)-qpl(icrm,i,j,kb) )  &
                   + ( bbb*fac_sub -(1.D0+fac_sub *dqsat)*tabs(icrm,i,j,k) ) * ( qpi(icrm,i,j,kc)-qpi(icrm,i,j,kb) ) )

              buoy_sgs_vert(icrm,i,j,k) = buoy_sgs
              a_prod_bu_vert(icrm,i,j,k) = -0.5D0*(tkh(icrm,i,j,kc)+tkh(icrm,i,j,kb)+0.002D0)*buoy_sgs
            end if ! if saturated at interface
          end if ! if either layer is cloudy
        end do ! i
      end do ! j
    enddo !k
  enddo !icrm

  !$acc parallel loop collapse(2) async(asyncid)
  do k = 1,nzm-1
    do icrm = 1 , ncrms
      tkelediss(icrm,k)  = 0.D0
      tkesbdiss(icrm,k)  = 0.D0
      tkesbshear(icrm,k) = 0.D0
      tkesbbuoy(icrm,k)  = 0.D0
    enddo
  enddo

  !$acc parallel loop collapse(4) async(asyncid)
  do k = 1,nzm-1
    do j = 1,ny
      do i = 1,nx
        do icrm = 1 , ncrms
          grd = dz(icrm)*adz(icrm,k)
          Ce1 = Ce/0.7D0*0.19D0
          Ce2 = Ce/0.7D0*0.51D0
          !!! compute correction factors for eddy visc/cond not to acceed 3D stability
          cx = dx**2.D0/dt/grdf_x(icrm,k)
          cy = dy**2.D0/dt/grdf_y(icrm,k)
          cz = (dz(icrm)*min(adzw(icrm,k),adzw(icrm,k+1)))**2.D0/dt/grdf_z(icrm,k)
          !!! maximum value of eddy visc/cond
          tkmax = 0.09D0/(1.D0/cx+1.D0/cy+1.D0/cz)
          buoy_sgs = 0.5D0*( buoy_sgs_vert(icrm,i,j,k-1) + buoy_sgs_vert(icrm,i,j,k) )
          if(buoy_sgs.le.0.) then
            smix = grd
          else
            smix = min(grd,max(0.1D0*grd, sqrt(0.76D0*tk(icrm,i,j,k)/Ck/sqrt(buoy_sgs+1.D-10))))
          end if
          ratio = smix/grd
          Cee = Ce1+Ce2*ratio
          if(dosmagor) then
            tk(icrm,i,j,k) = sqrt(Ck**3/Cee*max(0._crm_rknd,def2(icrm,i,j,k)-Pr*buoy_sgs))*smix**2
            tke(icrm,i,j,k) = (tk(icrm,i,j,k)/(Ck*smix))**2
            a_prod_sh = (tk(icrm,i,j,k)+0.001D0)*def2(icrm,i,j,k)
            ! a_prod_bu=-(tk(icrm,i,j,k)+0.001)*Pr*buoy_sgs
            a_prod_bu = 0.5D0*( a_prod_bu_vert(icrm,i,j,k-1) + a_prod_bu_vert(icrm,i,j,k) )
            a_diss = a_prod_sh+a_prod_bu
          else
            tke(icrm,i,j,k) = max(real(0.,crm_rknd),tke(icrm,i,j,k))
            a_prod_sh = (tk(icrm,i,j,k)+0.001D0)*def2(icrm,i,j,k)
            a_prod_bu = 0.5D0*( a_prod_bu_vert(icrm,i,j,k-1) + a_prod_bu_vert(icrm,i,j,k) )
            !!! cap the diss rate (useful for large time steps)
            a_diss = min(tke(icrm,i,j,k)/(4.D0*dt),Cee/smix*tke(icrm,i,j,k)**1.5D0)
            tke(icrm,i,j,k) = max(real(0.,crm_rknd),tke(icrm,i,j,k)+dtn*(max(0._crm_rknd,a_prod_sh+a_prod_bu)-a_diss))
            tk(icrm,i,j,k)  = Ck*smix*sqrt(tke(icrm,i,j,k))
          end if
          tk(icrm,i,j,k)  = min(tk(icrm,i,j,k),tkmax)
          tkh(icrm,i,j,k) = Pr*tk(icrm,i,j,k)

          tmp = a_prod_sh/dble(nx*ny)
          !$acc atomic update
          tkelediss(icrm,k)  = tkelediss(icrm,k) - tmp
          !$acc atomic update
          tkesbdiss(icrm,k)  = tkesbdiss(icrm,k) + a_diss
          !$acc atomic update
          tkesbshear(icrm,k) = tkesbshear(icrm,k)+ a_prod_sh
          !$acc atomic update
          tkesbbuoy(icrm,k)  = tkesbbuoy(icrm,k) + a_prod_bu
        end do ! i
      end do ! j
    end do ! k
  enddo !icrm

  deallocate( def2           )
  deallocate( buoy_sgs_vert  )
  deallocate( a_prod_bu_vert )

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------

end subroutine tke_full

end module tke_full_mod
