module shoc

!--------------------------------------------------------------
! SHOC parameterization
!   SHOC = Simplified Higher Order Closure 
!   reference, Bogenschutz and Krueger 2013
!
! Parameterization for low clouds and turbulence   
! email: bogenschutz1@llnl.gov
!--------------------------------------------------------------

use shr_kind_mod,  only: r8=>shr_kind_r8

implicit none

!=========================================================
! Constants set in initialization
!=========================================================

real(r8) :: ggr 	! gravity
real(r8) :: rgas	! dry air gas constant
real(r8) :: rv
real(r8) :: cp
real(r8) :: lcond

!=========================================================
! Private module parameters
!=========================================================
real(r8), parameter :: thl2tune=1.0_r8
real(r8), parameter :: qw2tune=1.0_r8
real(r8), parameter :: qwthl2tune=1.0_r8

logical, parameter :: dothetal_skew = .false.

real, parameter :: sqrt2 = sqrt(2._r8)
real, parameter :: sqrtpi = sqrt(2._r8*3.14_r8)
real, parameter :: basetemp = 300._r8

real, parameter :: maxlen = 2000.0_r8
real, parameter :: maxtke = 5.0_r8

!==============================================================
contains
!==============================================================

subroutine shoc_init( &
             gravit, rair, rh2o, cpair, &
	     latvap)

  implicit none
  
  ! Purpose:  Initialize constants for SHOC
	     
  real(r8), intent(in)  :: gravit
  real(r8), intent(in)  :: rair
  real(r8), intent(in)  :: rh2o
  real(r8), intent(in)  :: cpair
  real(r8), intent(in)  :: latvap
  
  ggr = gravit
  rgas = rair
  rv = rh2o
  cp = cpair
  lcond = latvap
  
end subroutine shoc_init	     				   

!==============================================================
! Main driver for the SHOC scheme

subroutine shoc_main ( &
     shcol, nlev, dtime, &
     host_dx, host_dy, &
     zt_grid,zm_grid,pres,&
     tke, thetal, qw, w_field,&
     wthl_sfc, wqw_sfc, uw_sfc, vw_sfc, &
     u_wind, v_wind,cldliq,qtracers,&
     num_qtracers,wthv_sec,&
     shoc_cldfrac,shoc_ql )
     
  implicit none
  
  ! input variables
  integer, intent(in) :: shcol  ! number of columns [-]
  integer, intent(in) :: nlev	! number of levels [-]
  integer, intent(in) :: num_qtracers  ! number of tracers [-]
  
  real(r8), intent(in) :: dtime	! SHOC timestep [s]
  real(r8), intent(in) :: host_dx(shcol) ! grid spacing of host model in x direction [m]
  real(r8), intent(in) :: host_dy(shcol) ! grid spacing of host model in y direction [m]
  real(r8), intent(in) :: zt_grid(shcol,nlev)	! heights, for thermo grid [m]
  real(r8), intent(in) :: zm_grid(shcol,nlev)   ! heights, for momentum grid [m]
  real(r8), intent(in) :: pres(shcol,nlev) ! pressure levels on thermo grid [hPa]

  real(r8), intent(in) :: cldliq(shcol,nlev) ! cloud liquid mixing ratio [kg/kg]
  real(r8), intent(in) :: w_field(shcol,nlev)
  
  real(r8), intent(in) :: wthl_sfc(shcol)
  real(r8), intent(in) :: wqw_sfc(shcol)
  real(r8), intent(in) :: uw_sfc(shcol)
  real(r8), intent(in) :: vw_sfc(shcol)

  ! input/output variables  
  real(r8), intent(inout) :: tke(shcol,nlev)  ! turbulent kinetic energy [m2/s2]
  real(r8), intent(inout) :: thetal(shcol,nlev) ! liquid water potential 
                                            !   temperature [K]
  real(r8), intent(inout) :: qw(shcol,nlev) ! total water mixing ratio [kg/kg]					    
  real(r8), intent(inout) :: u_wind(shcol,nlev) ! u wind component [m/s]
  real(r8), intent(inout) :: v_wind(shcol,nlev) ! v wind component [m/s]
  real(r8), intent(inout) :: wthv_sec(shcol,nlev) ! buoyancy flux [K m/s]
  real(r8), intent(inout) :: qtracers(shcol,nlev,num_qtracers) ! tracers [varies]
  
  ! output variables
  real(r8), intent(out) :: shoc_cldfrac(shcol,nlev)
  real(r8), intent(out) :: shoc_ql(shcol,nlev)
  
  ! Local variables
  real(r8) :: wtracer_sec(shcol,nlev,num_qtracers)
  
  real(r8) :: shoc_mix(shcol,nlev) ! Turbulent length scale [m]
  real(r8) :: tk(shcol,nlev) ! eddy viscosity
  real(r8) :: tkh(shcol,nlev) 
  real(r8) :: isotropy(shcol,nlev)
  real(r8) :: w_sec(shcol,nlev)
  real(r8) :: thl_sec(shcol,nlev)
  real(r8) :: qw_sec(shcol,nlev)
  real(r8) :: qwthl_sec(shcol,0:nlev)
  real(r8) :: wthl_sec(shcol,0:nlev)
  real(r8) :: wqw_sec(shcol,0:nlev)
  real(r8) :: wtke_sec(shcol,0:nlev)
  real(r8) :: uw_sec(shcol,0:nlev)
  real(r8) :: vw_sec(shcol,0:nlev)
  real(r8) :: w3(shcol,nlev)
  real(r8) :: brunt(shcol,nlev)
  
  real(r8) :: adz_zt(shcol,nlev), adz_zm(shcol,nlev)
  real(r8) :: dz(shcol)

  ! Define vertical grid arrays needed for 
  !   vertical derivatives     
  call shoc_grid(shcol,nlev,zt_grid,zm_grid,& ! Input
         adz_zt,adz_zm,dz)  		      ! Output

  ! Update the turbulent length scale	 
  call shoc_length(shcol,nlev,tke,&	      ! Input
         host_dx, host_dy, &                  ! Input
         cldliq,zt_grid,dz,adz_zt,adz_zm,&    ! Input
	 thetal,wthv_sec,&                    ! Input
	 brunt,shoc_mix)  		      ! Output

  ! Advance the SGS TKE equation	 
  call shoc_tke(shcol,nlev,dtime,&            ! Input
         wthv_sec,shoc_mix,dz,adz_zm,&        ! Input
	 u_wind,v_wind,brunt,&                ! Input
         tke,&				      ! Input/Output
	 tk,tkh,isotropy)                     ! Output
  
  ! Diagnose the second order moments
  call diag_second_shoc_moments(shcol,nlev,&  ! Input
         num_qtracers,thetal,qw,&             ! Input
	 u_wind,v_wind,qtracers,tke,&         ! Input
	 isotropy,tkh,tk,shoc_mix,&           ! Input
	 adz_zt,adz_zm, dz,&                  ! Input
	 zt_grid,zm_grid,&                    ! Input
	 wthl_sfc,wqw_sfc,uw_sfc,vw_sfc,&     ! Input
         w_sec,thl_sec,qw_sec,&               ! Output
	 wthl_sec,wqw_sec,&                   ! Output
	 qwthl_sec,uw_sec,vw_sec,wtke_sec,&   ! Output
	 wtracer_sec)

  ! Diagnose the third moment of vertical velocity	 
  call diag_third_shoc_moments(shcol,nlev,&   ! Input
         w_sec,thl_sec,qw_sec,qwthl_sec,&     ! Input
	 wthl_sec,isotropy,brunt,&            ! Input
	 thetal,tke,&                         ! Input
	 adz_zt,adz_zm,dz,&                   ! Input
	 zt_grid,zm_grid,&                    ! Input
	 w3)                                  ! Output
	 
  ! Update thetal, qw, tracers, and wind components
  !   based on SGS mixing
  call update_prognostics(shcol,nlev,num_qtracers,& ! Input
         dtime,adz_zt,adz_zm,dz,wthl_sec,&    ! Input
         wqw_sec,wtke_sec,uw_sec,&            ! Input
	 vw_sec,wtracer_sec,&                 ! Input
         thetal,qw,qtracers,tke,&             ! Input/Output	
	 u_wind,v_wind)                       ! Input/Output
  
  ! Call the PDF to close on SGS cloud and turbulence
  call shoc_assumed_pdf(shcol,nlev,&          ! Input
         thetal,qw,w_field,thl_sec,qw_sec,&   ! Input
	 wthl_sec,w_sec,&                     ! Input
	 wqw_sec,qwthl_sec,w3,pres,&          ! Input
	 zt_grid,zm_grid,&                    ! Input
	 shoc_cldfrac,shoc_ql,wthv_sec)       ! Output
	
  return
     
end subroutine shoc_main

!==============================================================
! Define grid variables needed for the parameterization

subroutine shoc_grid(shcol,nlev,zt_grid,zm_grid,& ! Input 
              adz_zt,adz_zm,dz)  ! Output

  implicit none

  ! input variables
  integer, intent(in) :: shcol
  integer, intent(in) :: nlev
  real(r8), intent(in) :: zt_grid(shcol,nlev)
  real(r8), intent(in) :: zm_grid(shcol,nlev)
  
  ! output variables 
  real(r8), intent(out) :: adz_zt(shcol,nlev)
  real(r8), intent(out) :: adz_zm(shcol,nlev)
  real(r8), intent(out) :: dz(shcol)
  
  ! local variables
  integer :: i, k
  
  ! difference in lowest two layers
  do i=1,shcol
    dz(i) = 0.5_r8 * (zt_grid(i,3) + zt_grid(i,2))
  enddo
  
  do k=1,nlev-1
    do i=1,shcol
      adz_zt(i,k) = (zm_grid(i,k+1) - zm_grid(i,k))/dz(i)
      
      if (k .eq. 1) then
        adz_zm(i,k) = 1._r8
      else
        adz_zm(i,k) = (zt_grid(i,k) - zt_grid(i,k-1))/dz(i)
      endif
      
    enddo
  enddo
  
  adz_zm(:,nlev) = adz_zm(:,nlev-1)
  adz_zt(:,nlev) = adz_zt(:,nlev-1)
  
  return
  
end subroutine shoc_grid

!==============================================================
! Update T, q, tracers, tke, u, and v based on SGS mixing

subroutine update_prognostics(shcol,nlev,num_tracer,& ! Input
             dtime,adz_zt,adz_zm,dz,wthl_sec,&        ! Input
	     wqw_sec,wtke_sec,uw_sec,&                ! Input
	     vw_sec,wtracer_sec,&                     ! Input
	     thetal,qw,tracer,tke,&                   ! Input/Output
	     u_wind,v_wind)                           ! Input/Output
	     
  implicit none
  
! Input variables
  integer, intent(in) :: shcol ! number of SHOC columns
  integer, intent(in) :: nlev  ! number of vertical levels
  integer, intent(in) :: num_tracer ! number of tracers
  real(r8), intent(in) :: dtime
  real(r8), intent(in) :: adz_zt(shcol,nlev)
  real(r8), intent(in) :: adz_zm(shcol,nlev)
  real(r8), intent(in) :: dz(shcol)
  real(r8), intent(in) :: wthl_sec(shcol,0:nlev)
  real(r8), intent(in) :: wqw_sec(shcol,0:nlev)
  real(r8), intent(in) :: wtracer_sec(shcol,0:nlev,num_tracer)
  real(r8), intent(in) :: uw_sec(shcol,0:nlev)
  real(r8), intent(in) :: vw_sec(shcol,0:nlev)
  real(r8), intent(in) :: wtke_sec(shcol,0:nlev)

! In/out variables  
  real(r8), intent(inout) :: thetal(shcol,nlev)
  real(r8), intent(inout) :: qw(shcol,nlev) ! total water mixing ratio   
  real(r8), intent(inout) :: tracer(shcol,nlev,num_tracer)
  real(r8), intent(inout) :: u_wind(shcol,nlev)
  real(r8), intent(inout) :: v_wind(shcol,nlev)
  real(r8), intent(inout) :: tke(shcol,nlev)
  
! Local variables  
  integer :: kb, k, i, p  
  real(r8) :: thedz, thedz_zm
  
  do i=1,shcol
    do k=1,nlev
      kb=k-1
      thedz=1._r8/(dz(i)*adz_zt(i,k))
      thedz_zm=1._r8/(dz(i)*adz_zm(i,k))

      thetal(i,k)=thetal(i,k)-dtime*(wthl_sec(i,k)-wthl_sec(i,kb))*thedz
      qw(i,k)=qw(i,k)-dtime*(wqw_sec(i,k)-wqw_sec(i,kb))*thedz

      tke(i,k)=tke(i,k)-dtime*(wtke_sec(i,k)-wtke_sec(i,kb))*thedz

      do p=1,num_tracer
        tracer(i,k,p)=tracer(i,k,p)-dtime*(wtracer_sec(i,k,p)-wtracer_sec(i,k,p))*thedz
      enddo
      
      u_wind(i,k)=u_wind(i,k)-dtime*(uw_sec(i,k)-uw_sec(i,kb))*thedz_zm
      v_wind(i,k)=v_wind(i,k)-dtime*(vw_sec(i,k)-vw_sec(i,kb))*thedz_zm
      
    enddo
  enddo
  
  return
  
end subroutine update_prognostics

!==============================================================
! SHOC Diagnose the second order moments

subroutine diag_second_shoc_moments(shcol,nlev, &   ! Input
             num_tracer,thetal,qw, &                ! Input
	     u_wind,v_wind,tracer,tke, &            ! Input
	     isotropy,tkh,tk,shoc_mix,&             ! Input
	     adz,adzw,dz,&                          ! Input   
	     zt_grid,zm_grid,&                      ! Input
	     wthl_sfc, wqw_sfc, uw_sfc, vw_sfc, &   ! Input    
	     w_sec, thl_sec, qw_sec,&               ! Output
	     wthl_sec,wqw_sec,&                     ! Output
	     qwthl_sec, uw_sec, vw_sec, wtke_sec, & ! Output
	     wtracer_sec)                           ! Output

  implicit none

! Input variables
  integer, intent(in) :: shcol ! number of SHOC columns
  integer, intent(in) :: nlev  ! number of vertical levels
  integer, intent(in) :: num_tracer ! number of tracers
  
  real(r8), intent(in) :: thetal(shcol,nlev) ! liquid water potential temp.
  real(r8), intent(in) :: qw(shcol,nlev) ! total water mixing ratio
  real(r8), intent(in) :: u_wind(shcol,nlev) ! u wind component
  real(r8), intent(in) :: v_wind(shcol,nlev) ! v wind component
  real(r8), intent(in) :: tke(shcol,nlev) ! turbulent kinetic energy
  real(r8), intent(in) :: isotropy(shcol,nlev) ! return to isotropy timescale
  real(r8), intent(in) :: tkh(shcol,nlev) ! thermal conductivity
  real(r8), intent(in) :: tk(shcol,nlev) ! momentum conductivity
  real(r8), intent(in) :: shoc_mix(shcol,nlev) ! mixing length
  real(r8), intent(in) :: tracer(shcol,nlev,num_tracer)
  
  real(r8), intent(in) :: zt_grid(shcol,nlev)
  real(r8), intent(in) :: zm_grid(shcol,nlev)
  real(r8), intent(in) :: adzw(shcol,nlev)
  real(r8), intent(in) :: adz(shcol,nlev) ! grid difference ratio
  real(r8), intent(in) :: dz(shcol) ! difference of lowest layers
  
  real(r8), intent(in) :: wthl_sfc(shcol)
  real(r8), intent(in) :: wqw_sfc(shcol)
  real(r8), intent(in) :: uw_sfc(shcol)
  real(r8), intent(in) :: vw_sfc(shcol)  

! Output variables
  real(r8), intent(out) :: w_sec(shcol,nlev)	! second order vertical velocity
  real(r8), intent(out) :: thl_sec(shcol,nlev)  ! second order liquid wat. potential temp.
  real(r8), intent(out) :: qw_sec(shcol,nlev)   ! second order total water mixing rat.
  real(r8), intent(out) :: qwthl_sec(shcol,0:nlev) ! covariance of temp and moisture
  real(r8), intent(out) :: wthl_sec(shcol,0:nlev) ! vertical flux of heat
  real(r8), intent(out) :: wqw_sec(shcol,0:nlev) ! vertical flux of total water
  real(r8), intent(out) :: uw_sec(shcol,0:nlev) ! vertical flux of u
  real(r8), intent(out) :: vw_sec(shcol,0:nlev) ! vertical flux of v
  real(r8), intent(out) :: wtke_sec(shcol,0:nlev) ! vertical flux of tke
  real(r8), intent(out) :: wtracer_sec(shcol,0:nlev,num_tracer) ! vertical flux
                                                              ! of tracer

! Local variables  
  integer :: kb, k, i, p
  real(r8) :: grid_dz2,grid_dz,grid_dzw
  real(r8) :: gr1,gr2,grw1
  real(r8) :: isotropy_zm(shcol,nlev)
  real(r8) :: tkh_zm(shcol,nlev)
  real(r8) :: sm ! Mixing coefficient
  
  call linear_interp(zt_grid,zm_grid,isotropy,isotropy_zm,nlev,nlev,shcol)
  call linear_interp(zt_grid,zm_grid,tkh,tkh_zm,nlev,nlev,shcol)
  
  w_sec = (2._r8/3._r8)*tke
  
  do k=1,nlev
    do i=1,shcol
    
      kb=k-1
      grid_dz=dz(i)*adz(i,k)
      grid_dzw=dz(i)*adzw(i,k)
      if (k .eq. 1) then
        kb=1
        grid_dz=dz(i)*adz(i,2)
	grid_dzw=dz(i)*adzw(i,2)
      endif
      if (k .eq. nlev) then
        kb=nlev-1
        grid_dz=dz(i)*adz(i,k)
	grid_dzw=dz(i)*adzw(i,k)
      endif
    
      gr1=(1._r8/grid_dz)
      gr2=(1._r8/grid_dz)**2
      grw1=(1._r8/grid_dzw)
    
      sm=isotropy_zm(i,k)*tkh_zm(i,k)
      
      thl_sec(i,k)=thl2tune*sm*gr2*(thetal(i,k)-thetal(i,kb))**2
      qw_sec(i,k)=qw2tune*sm*gr2*(qw(i,k)-qw(i,kb))**2
      qwthl_sec(i,k)=qwthl2tune*sm*gr2*((thetal(i,k)-thetal(i,kb))* &
        (qw(i,k)-qw(i,kb)) )
	
      ! diagnose vertical heat flux
      wthl_sec(i,k)=-1._r8*tkh_zm(i,k)*gr1*(thetal(i,k)-thetal(i,kb))
      wqw_sec(i,k)=-1._r8*tkh_zm(i,k)*gr1*(qw(i,k)-qw(i,kb))
      
      wtke_sec(i,k)=-1._r8*tkh_zm(i,k)*gr1*(tke(i,k)-tke(i,kb))
      
      ! diagnose vertical momentum transport
      uw_sec(i,k)=-1._r8*tk(i,k)*grw1*(u_wind(i,k)-u_wind(i,kb))
      vw_sec(i,k)=-1._r8*tk(i,k)*grw1*(v_wind(i,k)-v_wind(i,kb))
      
      do p=1,num_tracer
        wtracer_sec(i,k,p)=tkh_zm(i,k)*gr1*(tracer(i,k,p)-tracer(i,kb,p))
      enddo
    
    enddo ! end i loop (column loop)
  enddo  ! end k loop (vertical loop)
  
  ! apply the surface fluxes (these may need to go into a 0 k dimension)
  do i=1,shcol
    wthl_sec(i,0) = wthl_sfc(i)
    wqw_sec(i,0) = wqw_sfc(i)
    uw_sec(i,0) = uw_sfc(i)
    vw_sec(i,0) = vw_sfc(i)
    wtracer_sec(i,0,:) = 0._r8
    wtke_sec(i,:) = 0._r8
  enddo
  
  return

end subroutine diag_second_shoc_moments

!==============================================================
! SHOC Diagnose the third order moments

subroutine diag_third_shoc_moments(shcol,nlev, &  ! Input 
	      w_sec, thl_sec, qw_sec, qwthl_sec,& ! Input
	      wthl_sec, isotropy, brunt,&         ! Input
	      thetal,tke,&                        ! Input
	      adz, adzw, dz, &                    ! Input
	      zt_grid,zm_grid,&                   ! Input
	      w3)                                 ! Output
	      
  implicit none

! Input variables
  integer, intent(in) :: shcol ! number of SHOC columns
  integer, intent(in) :: nlev  ! number of vertical levels
  
  real(r8), intent(in) :: w_sec(shcol,nlev)	! second order vertical velocity
  real(r8), intent(in) :: thl_sec(shcol,nlev)  ! second order liquid wat. potential temp.
  real(r8), intent(in) :: qw_sec(shcol,nlev)   ! second order total water mixing rat.
  real(r8), intent(in) :: qwthl_sec(shcol,nlev) ! covariance of temp and moisture  
  real(r8), intent(in) :: wthl_sec(shcol,0:nlev) ! vertical flux of heat 
  real(r8), intent(in) :: isotropy(shcol,nlev)
  real(r8), intent(in) :: brunt(shcol,nlev)
  real(r8), intent(in) :: thetal(shcol,nlev)
  real(r8), intent(in) :: tke(shcol,nlev)
  
  real(r8), intent(in) :: adz(shcol,nlev) ! grid difference ratio
  real(r8), intent(in) :: adzw(shcol,nlev) 
  real(r8), intent(in) :: dz(shcol) ! difference of lowest layers	
  
  real(r8), intent(in) :: zt_grid(shcol,nlev)
  real(r8), intent(in) :: zm_grid(shcol,nlev)
  
! Output variables
  real(r8), intent(out) :: w3(shcol,nlev)
  
! Local variables
  integer i, j, k, kb, kc
  real(r8) :: omega0, omega1, omega2
  real(r8) :: X0, Y0, X1, Y1, AA0, AA1
  real(r8) :: zvar, x5var, iso, thedz, thedz2
  real(r8) :: theterm, cond, tsign
  real(r8) :: isosqrt, dthl2, dwthl, dtke, dw2, aw2  
  real(r8) :: a0, a1, a2, a3, a4, a5 
  real(r8) :: buoy_sgs2, c, grd, bet2, bet
  real(r8) :: f0, f1, f2, f3, f4, f5
  
  real(r8) :: w_sec_zm(shcol,nlev)	! second order vertical velocity 
  real(r8) :: isotropy_zm(shcol,nlev)
  real(r8) :: brunt_zm(shcol,nlev)  

  call linear_interp(zt_grid,zm_grid,isotropy,isotropy_zm,nlev,nlev,shcol)
  call linear_interp(zt_grid,zm_grid,brunt,brunt_zm,nlev,nlev,shcol)
  call linear_interp(zt_grid,zm_grid,w_sec,w_sec_zm,nlev,nlev,shcol)
  
  c=7._r8
  a0=(0.52_r8*c**(-2))/(c-2._r8)
  a1=0.87_r8/(c**2)
  a2=0.5_r8/c
  a3=0.6_r8/(c*(c-2._r8))
  a4=2.4_r8/(3._r8*c+5._r8)
  a5=0.6_r8/(c*(3._r8+5._r8*c))    
  
  do k=1,nlev  
    do i=1,shcol
    
      if(k.eq.1) then
        kb=1
        kc=2
        thedz=dz(i)*adz(i,kc)
        thedz2=dz(i)*adz(i,kc)
      elseif(k.eq.nlev) then
        kb=nlev-1
        kc=nlev
        thedz=dz(i)*adz(i,k)
        thedz2=dz(i)*adz(i,k)
      else 
        kb=k-1
        kc=k+1
        thedz=dz(i)*adz(i,k)
        thedz2=dz(i)*(adz(i,kc)+adz(i,k))
     endif 

     grd=dz(i)*adz(i,k) 
    
      thedz=1._r8/thedz
      thedz2=1._r8/thedz2    
      
      iso=isotropy_zm(i,k)
      isosqrt=iso**2
      buoy_sgs2=isosqrt*brunt_zm(i,k)
      bet2=0.5_r8*(ggr/(thetal(i,k))+ ggr/(thetal(i,kb)))
      
      f0=thedz2 * bet2**3 * iso**4 * wthl_sec(i,k) * &
        (thl_sec(i,kc)-thl_sec(i,kb))
	
      f1=thedz2 * bet2**2 * iso**3 * wthl_sec(i,k) * &
        (wthl_sec(i,kc)-wthl_sec(i,kb)) + 0.5 * &
	w_sec_zm(i,k)*(thl_sec(i,kc)-thl_sec(i,kb))
	
      f2=thedz * bet2 * isosqrt * wthl_sec(i,k) * &
        (w_sec(i,k)-w_sec(i,kb))+ 2._r8 * thedz2 * bet2 * &    ! +DPAB check bet2 here
	isosqrt * w_sec_zm(i,k) * (wthl_sec(i,kc) - wthl_sec(i,kb)) 
	
      f3=thedz * bet2 * isosqrt * w_sec_zm(i,k) * &
        (wthl_sec(i,kc) - wthl_sec(i,kb)) + thedz * &
	bet2 * isosqrt * (wthl_sec(i,k) * (tke(i,k) - tke(i,kb)))
	
      f4=thedz * iso * w_sec_zm(i,k) * ((w_sec(i,k) - w_sec(i,kb) + &
        (tke(i,k) - tke(i,kb))))
	
      f5=thedz * iso * w_sec_zm(i,k) * (w_sec(i,k) - w_sec(i,kb))
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Compute the omega terms     
      
      omega0 = a4 / (1._r8 - a5 * buoy_sgs2)
      omega1 = omega0/(2._r8 * c)
      omega2 = omega1 * f3 + (5._r8/4._r8) * omega0 * f4
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Compute the X0, Y0, X1, Y1 terms  
      
      X0 = (a2 * buoy_sgs2 * (1._r8 - a3 * buoy_sgs2)) / &
        (1._r8 - (a1 + a3) * buoy_sgs2)
      Y0 = (2._r8 * a2 * buoy_sgs2 * X0) / (1._r8 - a3 * buoy_sgs2)
      X1 = (a0 * f0 + a1 * f1 + a2 * (1._r8 - a3 * buoy_sgs2) * f2) / &
        (1._r8 - (a1 + a3) * buoy_sgs2)
      Y1 = (2._r8 * a2 * (buoy_sgs2 * X1 + (a0/a1) * f0 + f1) / &
        (1._r8 - a3* buoy_sgs2))     

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
      ! Compute the A0, A1 terms	
      
      AA0 = omega0 * X0 + omega1 * Y0
      AA1 = omega0 * X1 + omega1 * Y1 + omega2      
      
      ! Finally, we have the third moment of w
      w3(i,k)=(AA1-1.2*X1-1.5*f5)/(c-1.2*X0+AA0)      
      
    enddo  ! end i loop
  
  enddo  ! end k loop
  
  ! perform clipping to prevent unrealistically large values from occuring
  
  do k=1,nlev
    do i=1,shcol
    
      tsign = 1._r8
      theterm = w_sec_zm(i,k)
      cond = 1.2_r8 * sqrt(2._r8 * theterm**3)
      if (w3(i,k) .lt. 0) tsign = -1._r8
      if (tsign * w3(i,k) .gt. cond) w3(i,k) = tsign * cond
      
    enddo ! end i loop
  enddo ! end k loop
  
  return
  
end subroutine diag_third_shoc_moments

!==============================================================
! Assumed PDF closure for the SHOC scheme

subroutine shoc_assumed_pdf(shcol,nlev, &       ! Input
             thetal,qw,w_field,thl_sec,qw_sec,& ! Input
	     wthl_sec,w_sec, &                  ! Input
	     wqw_sec,qwthl_sec,w3,pres, &       ! Input
	     zt_grid,zm_grid,&                  ! Input
	     shoc_cldfrac,shoc_ql,wthv_sec)     ! Output

  implicit none

  ! Input variables
  integer, intent(in) :: shcol
  integer, intent(in) :: nlev
  
  real(r8), intent(in) :: thetal(shcol,nlev)
  real(r8), intent(in) :: qw(shcol,nlev)
  real(r8), intent(in) :: thl_sec(shcol,nlev)
  real(r8), intent(in) :: qw_sec(shcol,nlev)
  real(r8), intent(in) :: wthl_sec(shcol,0:nlev)
  real(r8), intent(in) :: w_sec(shcol,nlev)
  real(r8), intent(in) :: wqw_sec(shcol,0:nlev)
  real(r8), intent(in) :: qwthl_sec(shcol,nlev)
  real(r8), intent(in) :: w3(shcol,nlev)
  real(r8), intent(in) :: w_field(shcol,nlev)
  real(r8), intent(in) :: pres(shcol,nlev)
  real(r8), intent(in) :: zt_grid(shcol,nlev)
  real(r8), intent(in) :: zm_grid(shcol,nlev)
  
  ! Output variables
  real(r8), intent(out) :: shoc_cldfrac(shcol,nlev)
  real(r8), intent(out) :: shoc_ql(shcol,nlev)
  real(r8), intent(out) :: wthv_sec(shcol,nlev)

  ! Local variables
  integer i,j,k,dothis,nmicro_fields
  real(r8) skew_w,skew_thl,skew_qw,a
  real(r8) w1_1,w1_2,w2_1,w2_2,w3var
  real(r8) thl1_1,thl1_2,thl2_1,thl2_2
  real(r8) qw1_1,qw1_2,qw2_1,qw2_2
  real(r8) r_qwthl_1,r_wqw_1,r_wthl_1
  real(r8) testvar,s1,s2,std_s1,std_s2,C1,C2
  real(r8) ql1,ql2,flux_ws1,flux_ws2
  real(r8) w_ql1,w_ql2,thl_first,qw_first,w_first
  real(r8) Tl1_1,Tl1_2,betatest,pval
  real(r8) w2thl,w2qw,w2ql,w2ql_1,w2ql_2
  real(r8) s,thec,thlsec,qwsec,qwthlsec,wqwsec,wthlsec
  real(r8) thestd,erf,exp,dum
  real(r8) cqt1,cthl1,cqt2,cthl2
  real(r8) qn1,qn2,qi1,qi2,omn1,omn2, omn
  real(r8) lstarn, test
  real(r8) m, bigm, basetemp2
  real(r8) beta1, beta2, qs1, qs2
  real(r8) esval1_1, esval2_1, esval1_2, esval2_2, om1, om2
  real(r8) lstarn1, lstarn2, sqrtw2, sqrtthl, sqrtqt
  real(r8) sqrtstd1, sqrtstd2, tsign, tvar, sqrtw2t
  real(r8) wqls, wqis, skip, epsterm
  real(r8) sqrtqw2_1, sqrtqw2_2, sqrtthl2_1, sqrtthl2_2
  real(r8) corrtest1, corrtest2, thl_tol, rt_tol, w_tol_sqd, w_thresh
  
  real(r8) :: wthl_sec_zt(shcol,nlev)
  real(r8) :: wqw_sec_zt(shcol,nlev)
  real(r8) :: w3_zt(shcol,nlev)
  real(r8) :: thl_sec_zt(shcol,nlev)
  real(r8) :: qwthl_sec_zt(shcol,nlev)
  real(r8) :: qw_sec_zt(shcol,nlev)
  real(r8) :: w_field_zt(shcol,nlev)
  
  epsterm=rgas/rv 
  
  thl_tol=1.e-2_r8
  rt_tol=1.e-4_r8
  w_tol_sqd=(2.e-2_r8)**2
  w_thresh=0.0_r8
  
  ! Initialize cloud variables to zero  
  shoc_cldfrac(:,:)=0._r8
  shoc_ql(:,1)=0._r8 

  call linear_interp(zm_grid,zt_grid,w3,w3_zt,nlev,nlev,shcol)  
  call linear_interp(zm_grid,zt_grid,thl_sec,thl_sec_zt,nlev,nlev,shcol)
  call linear_interp(zm_grid,zt_grid,wthl_sec(:,1:nlev),wthl_sec_zt,nlev,nlev,shcol)
  call linear_interp(zm_grid,zt_grid,qwthl_sec,qwthl_sec_zt,nlev,nlev,shcol)
  call linear_interp(zm_grid,zt_grid,wqw_sec(:,1:nlev),wqw_sec_zt,nlev,nlev,shcol)
  call linear_interp(zm_grid,zt_grid,qw_sec,qw_sec_zt,nlev,nlev,shcol)  
  call linear_interp(zm_grid,zt_grid,w_field,w_field_zt,nlev,nlev,shcol)
  
! Initialize the buoyancy flux (probably doesn't need to be here)
  wthv_sec(:,1)=wthl_sec_zt(:,1)+((1.-epsterm)/epsterm)&
     *basetemp*wqw_sec(:,1)
    
  do k=1,nlev
    do i=1,shcol

      pval = pres(i,k) 

      ! Get all needed input moments for the PDf
      !  at this particular point
      thl_first = thetal(i,k)
      w_first = w_field_zt(i,k)
      qw_first = qw(i,k)
      
      w3var = w3_zt(i,k)
      thlsec = thl_sec_zt(i,k)
      qwsec = qw_sec_zt(i,k)
      qwthlsec = qwthl_sec_zt(i,k)
      wqwsec = wqw_sec_zt(i,k)
      wthlsec = wthl_sec_zt(i,k)
      
      ! Compute square roots of some variables so we don't
      !  have to compute these again
      sqrtw2 = sqrt(w_sec(i,k))
      sqrtthl = max(thl_tol,sqrt(thlsec))
      sqrtqt = max(rt_tol,sqrt(qwsec))
 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !  FIND PARAMETERS FOR VERTICAL VELOCITY 
      
      Skew_w=w3var/w_sec(i,k)**(3./2.)

      if (w_sec(i,k) .le. w_tol_sqd) then
        Skew_w=0.
        w1_1=w_first
        w1_2=w_first
        w2_1=0.
        w2_2=0.
        a=0.5
      else

        w2_1=0.4
        w2_2=0.4

        a=max(0.01,min(0.5*(1.-Skew_w*sqrt(1./(4.*(1.-w2_1)**3+Skew_w**2))),0.99))

        sqrtw2t=sqrt(1.-w2_1)

        w1_1=sqrt((1.-a)/a)*sqrtw2t
        w1_2=-1.*sqrt(a/(1.-a))*sqrtw2t

        w2_1=w2_1*w_sec(i,k)
        w2_2=w2_2*w_sec(i,k)


      endif  
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !  FIND PARAMETERS FOR THETAL
    
      corrtest1=max(-1.0,min(1.0,wthlsec/(sqrtw2*sqrtthl)))

      if (thlsec .le. thl_tol**2 .or. abs(w1_2-w1_1) .le. w_thresh) then
        thl1_1=thl_first
        thl1_2=thl_first
        thl2_1=0.
        thl2_2=0.
        sqrtthl2_1=0.
        sqrtthl2_2=0.
      else

        thl1_1=(-1.*corrtest1)/w1_2
        thl1_2=(-1.*corrtest1)/w1_1
      
        if (dothetal_skew) then
          tsign=abs(thl1_2-thl1_1)
	
	  if (tsign .gt. 0.4) then
	    Skew_thl=1.2*Skew_w
	  else if (tsign .le. 0.2) then
	    Skew_thl=0.0
	  else
	    Skew_thl=((1.2*Skew_w)/0.2)*(tsign-0.2)
	  endif 
        else
          Skew_thl = 0.0
        endif
	
        thl2_1=min(100.,max(0.,(3.*thl1_2*(1.-a*thl1_1**2-(1.-a)*thl1_2**2) &     
	    -(Skew_thl-a*thl1_1**3-(1.-a)*thl1_2**3))/ &
	    (3.*a*(thl1_2-thl1_1))))*thlsec
      
        thl2_2=min(100.,max(0.,(-3.*thl1_1*(1.-a*thl1_1**2-(1.-a)*thl1_2**2) &
          +(Skew_thl-a*thl1_1**3-(1.-a)*thl1_2**3))/ &
          (3.*(1.-a)*(thl1_2-thl1_1))))*thlsec

        thl1_1=thl1_1*sqrtthl+thl_first
        thl1_2=thl1_2*sqrtthl+thl_first

        sqrtthl2_1=sqrt(thl2_1)
        sqrtthl2_2=sqrt(thl2_2)

      endif	                  
 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  FIND PARAMETERS FOR TOTAL WATER MIXING RATIO

      corrtest2=max(-1.0,min(1.0,wqwsec/(sqrtw2*sqrtqt)))

      if (qwsec .le. rt_tol**2 .or. abs(w1_2-w1_1) .le. w_thresh) then
        qw1_1=qw_first
        qw1_2=qw_first
        qw2_1=0.
        qw2_2=0.
        sqrtqw2_1=0.
        sqrtqw2_2=0.
      else

        qw1_1=(-1.*corrtest2)/w1_2
        qw1_2=(-1.*corrtest2)/w1_1      

        tsign=abs(qw1_2-qw1_1)

        if (tsign .gt. 0.4) then
          Skew_qw=1.2*Skew_w
        else if (tsign .le. 0.2) then
          Skew_qw=0.
        else
          Skew_qw=((1.2*Skew_w)/0.2)*(tsign-0.2)
        endif

      qw2_1=min(100.,max(0.,(3.*qw1_2*(1.-a*qw1_1**2-(1.-a)*qw1_2**2) &
        -(Skew_qw-a*qw1_1**3-(1.-a)*qw1_2**3))/ &
        (3.*a*(qw1_2-qw1_1))))*qwsec

      qw2_2=min(100.,max(0.,(-3.*qw1_1*(1.-a*qw1_1**2-(1.-a)*qw1_2**2) &
      +(Skew_qw-a*qw1_1**3-(1.-a)*qw1_2**3))/ &
      (3.*(1.-a)*(qw1_2-qw1_1))))*qwsec

      qw1_1=qw1_1*sqrtqt+qw_first
      qw1_2=qw1_2*sqrtqt+qw_first

      sqrtqw2_1=sqrt(qw2_1)
      sqrtqw2_2=sqrt(qw2_2)

      endif
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !  CONVERT FROM TILDA VARIABLES TO "REAL" VARIABLES

      w1_1=w1_1*sqrtw2+w_first
      w1_2=w1_2*sqrtw2+w_first
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !  FIND WITHIN-PLUME CORRELATIONS 

      testvar=(a*sqrtqw2_1*sqrtthl2_1+(1.-a)*sqrtqw2_2*sqrtthl2_2)

      if (testvar .eq. 0) then
        r_qwthl_1=0.
      else
        r_qwthl_1=max(-1.0,min(1.0,(qwthlsec-a*(qw1_1-qw_first) &
	  *(thl1_1-thl_first)-(1.-a)*(qw1_2-qw_first) &
	  *(thl1_2-thl_first))/testvar))
      endif    
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !  BEGIN TO COMPUTE CLOUD PROPERTY STATISTICS

      ! THIS DEFINITION WILL NEED TO BE CHANGED! +DPAB
!      Tl1_1=thl1_1-ggr/cp*z(k)+fac_cond*qpl(i,j,k)+fac_sub*(qci(i,j,k)+qpi(i,j,k))
!      Tl1_2=thl1_2-ggr/cp*z(k)+fac_cond*qpl(i,j,k)+fac_sub*(qci(i,j,k)+qpi(i,j,k))

      Tl1_1=thl1_1/((100000._r8/pval)**(rgas/cp))
      Tl1_2=thl1_2/((100000._r8/pval)**(rgas/cp))

      ! Now compute qs

      esval1_1=0.
      esval1_2=0.
      esval2_1=0.
      esval2_2=0.
      om1=1.
      om2=1. 
     
      ! +DPAB, change call to consistent esatw_crm
      esval1_1=esatw_shoc(Tl1_1)*100.
      lstarn1=lcond 
	
      esval2_1=0._r8   ! +DPAB is this right?
      
      qs1=0.622_r8*esval1_1/max(esval1_1,pval-esval1_1)
      beta1=(rgas/rv)*(lstarn1/(rgas*Tl1_1))*(lstarn1/(cp*Tl1_1))
      
      ! Are the two plumes equal?  If so then set qs and beta
      ! in each column to each other to save computation
      lstarn2=lcond
      if (Tl1_1 .eq. Tl1_2) then
        qs2=qs1     
        beta2=beta1
      else
        esval1_2=esatw_shoc(Tl1_2)*100._r8 
	esval2_2 = 0._r8
      endif
      
      qs2=0.622*esval1_2/max(esval1_2,pval-esval1_2)
      beta2=(rgas/rv)*(lstarn2/(rgas*Tl1_2))*(lstarn2/(cp*Tl1_2)) 
      
      !!!!!  Now compute cloud stuff
      !!!!!!  compute s term    
      
      s1=qw1_1-qs1*((1.+beta1*qw1_1)/(1.+beta1*qs1))
      cthl1=((1.+beta1*qw1_1)/(1.+beta1*qs1)**2)*(cp/lcond) &
        *beta1*qs1*(pval/100000.)**(rgas/cp)

      cqt1=1./(1.+beta1*qs1)
      std_s1=sqrt(max(0.,cthl1**2*thl2_1+cqt1**2*qw2_1-2.*cthl1 &
        *sqrtthl2_1*cqt1*sqrtqw2_1*r_qwthl_1))
	
      qn1=0._r8
      C1=0._r8
      
      if (std_s1 .ne. 0) then
        C1=0.5*(1.+erf(s1/(sqrt2*std_s1)))
	IF (C1 .ne. C1) C1 = 0._r8
        IF (C1 .ne. 0._r8) qn1=s1*C1+(std_s1/sqrtpi)*exp(-0.5*(s1/std_s1)**2)
      else
        if (s1 .gt. 0._r8) then
          C1=1.0_r8
          qn1=s1
        endif
      endif
      
      !!!!! now compute non-precipitating cloud condensate 

      ! If two plumes exactly equal, then just set many of these 
      ! variables to themselves to save on computation.
      if (qw1_1 .eq. qw1_2 .and. thl2_1 .eq. thl2_2 .and. qs1 .eq. qs2) then
        s2=s1
        cthl2=cthl1
        cqt2=cqt1
        std_s2=std_s1
        C2=C1
        qn2=qn1
      else
        
        s2=qw1_2-qs2*((1.+beta2*qw1_2)/(1.+beta2*qs2))
        cthl2=((1.+beta2*qw1_2)/(1.+beta2*qs2)**2)*(cp/lcond) &
	  *beta2*qs2*(pval/100000.)**(rgas/cp)
        cqt2=1./(1.+beta2*qs2)
        std_s2=sqrt(max(0.,cthl2**2*thl2_2+cqt2**2*qw2_2-2.*cthl2* &
	  sqrtthl2_2*cqt2*sqrtqw2_2*r_qwthl_1))

        qn2=0._r8
        C2=0._r8
	
        if (std_s2 .ne. 0._r8) then
          C2=0.5_r8*(1.+erf(s2/(sqrt2*std_s2)))
	  if (C2 .ne. C2) C2 = 0._r8
          if (C2 .ne. 0._r8) qn2=s2*C2+(std_s2/sqrtpi)*exp(-0.5_r8*(s2/std_s2)**2)
        else
          if (s2 .gt. 0._r8) then
            C2=1.0_r8
            qn2=s2
          endif
        endif

      endif
      
      shoc_cldfrac(i,k) = min(1.,a*C1+(1._r8-a)*C1)
      
      ql1=min(qn1,qw1_1)
      ql2=min(qn2,qw1_2)
      
      shoc_ql(i,k) = max(0._r8,a*ql1+(1._r8+a)*ql2)
      
      ! update temperature here? +PAB
      
      ! Begin to compute the buoyancy flux +DPAB, CHECK DEFINITION!!!!
      wqls=a*((w1_1-w_first)*ql1)+(1._r8-a)*((w1_2-w_first)*ql2)
      wthv_sec(i,k)=wthlsec+((1._r8-epsterm)/epsterm)*basetemp* &
        wqwsec+((lcond/cp)-(1._r8/epsterm)*basetemp)*wqls
	
	
    enddo  ! end i loop here
  enddo	  ! end k loop here
	                      
  return                              
     
end subroutine shoc_assumed_pdf

!==============================================================
! Compute the turbulent kinetic energy

subroutine shoc_tke(shcol,nlev,dtime,& ! Input
             wthv_sec,shoc_mix,dz,adzw,&        ! Input  
	     u_wind,v_wind,brunt,&              ! Input
             tke, &                             ! Input/Output
             tk,tkh,isotropy)                   ! Output

  ! input variables
  integer, intent(in) :: shcol
  integer, intent(in) :: nlev
  real(r8), intent(in) :: dtime
  real(r8), intent(in) :: wthv_sec(shcol,nlev)
  real(r8), intent(in) :: shoc_mix(shcol,nlev)
  real(r8), intent(in) :: u_wind(shcol,nlev)
  real(r8), intent(in) :: v_wind(shcol,nlev)
  real(r8), intent(in) :: adzw(shcol,nlev)
  real(r8), intent(in) :: dz(shcol)
  real(r8), intent(in) :: brunt(shcol,nlev)
  
  ! output variables
  real(r8), intent(out) :: tkh(shcol,nlev)
  real(r8), intent(out) :: isotropy(shcol,nlev)  

  ! input/output variables
  real(r8), intent(inout) :: tke(shcol,nlev)
  real(r8), intent(inout) :: tk(shcol,nlev)
  
  ! local variables	     
  real(r8) :: shear_prod(shcol,nlev)	
  real(r8) :: grd,betdz,Ck,Ce,Ces,Ce1,Ce2,smix,Pr,Cee,Cs
  real(r8) :: buoy_sgs,ratio,a_prod_sh,a_prod_bu,a_diss
  real(r8) :: lstarn, lstarp, bbb, omn, omp
  real(r8) :: qsatt,dqsat,tk_in, uw_sec, vw_sec
  real(r8) :: tscale1,lambda,buoy_sgs_save,grid_dzw,grw1
  integer i,j,k,kc,kb	 
  
  Cs = 0.15_r8
  Ck=0.1_r8
  Ce=Ck**3/Cs**4
  Ces=Ce/0.7_r8*3.0_r8
  Pr=1.0_r8   
  
  Ce1=Ce/0.7_r8*0.19_r8
  Ce2=Ce/0.7_r8*0.51_r8

  ! Compute shear production term
  do k=1,nlev
    do i=1,shcol
    
      kb=k-1
      grid_dzw=dz(i)*adzw(i,k)
      if (k .eq. 1) then
        kb=1
	grid_dzw=dz(i)*adzw(i,2)
      endif
      if (k .eq. nlev) then
        kb=nlev-1
	grid_dzw=dz(i)*adzw(i,k)
      endif
 
      grw1=(1._r8/grid_dzw)      
   
      tk_in=Ck*shoc_mix(i,k)*sqrt(tke(i,k))
      uw_sec=grw1*(u_wind(i,k)-u_wind(i,kb))
      vw_sec=grw1*(v_wind(i,k)-v_wind(i,kb))  
      shear_prod(i,k)=uw_sec**2+vw_sec**2 
    enddo
  enddo
  
  do k=1,nlev
    do i=1,shcol
    
      tk_in=Ck*shoc_mix(i,k)*sqrt(tke(i,k))
      smix=shoc_mix(i,k)
      a_prod_bu=(ggr/basetemp)*wthv_sec(i,k)
      buoy_sgs=-1._r8*a_prod_bu/(Pr*(tk_in+0.001_r8))
      buoy_sgs_save=brunt(i,k)  
    
      if (buoy_sgs .le. 0._r8) then 
        smix=grid_dzw
      else
        smix=min(grid_dzw,max(0.1_r8*grid_dzw, sqrt(0.76_r8*tk_in/Ck/sqrt(buoy_sgs+1.0e-10_r8))))
      endif
      
      ratio=smix/grid_dzw
      Pr=1. 
      Cee=Ce1+Ce2*ratio  
      
      tke(i,k)=max(0._r8,tke(i,k))
      a_prod_sh=(tk_in+0.001_r8)*shear_prod(i,k)
      a_diss=Cee/shoc_mix(i,k)*tke(i,k)**1.5
      tke(i,k)=max(0._r8,tke(i,k)+dtime*(max(0._r8,a_prod_sh+a_prod_bu)-a_diss))  
      
      tke(i,k)=min(tke(i,k),5.0)
   
      tscale1=(2._r8*tke(i,k))/a_diss
      lambda=0.04_r8
      if (buoy_sgs_save .le. 0) lambda=0._r8
      isotropy(i,k)=min(2000._r8,tscale1/(1._r8+lambda*buoy_sgs_save*tscale1**2))
   
      tk(i,k)=Ck*smix*sqrt(tke(i,k))
      tkh(i,k)=Ck*isotropy(i,k)*tke(i,k)              
   
      tke(i,k) = max(0.0004,tke(i,k))
 
    enddo ! end i loop
  enddo ! end k loop
  
  return
  
end subroutine shoc_tke

!==============================================================
! Compute the turbulent length scale

subroutine shoc_length(shcol,nlev,tke, &   ! Input
             host_dx,host_dy, &        ! Input
             cldin,z_in,dz,adz,adzw,&  ! Input
	     thetal,wthv_sec,&		! Input
	     brunt,shoc_mix)  ! Output

  implicit none

  ! input variables
  integer, intent(in) :: shcol
  integer, intent(in) :: nlev
  
  real(r8), intent(in) :: host_dx(shcol)
  real(r8), intent(in) :: host_dy(shcol)
  real(r8), intent(in) :: tke(shcol,nlev)
  real(r8), intent(in) :: cldin(shcol,nlev)
  real(r8), intent(in) :: z_in(shcol,nlev)
  real(r8), intent(in) :: dz(shcol)
  real(r8), intent(in) :: adz(shcol,nlev)
  real(r8), intent(in) :: adzw(shcol,nlev)
  real(r8), intent(in) :: wthv_sec(shcol,nlev)
  real(r8), intent(in) :: thetal(shcol,nlev)

  ! output variables
  real(r8), intent(out) :: brunt(shcol,nlev)
  real(r8), intent(out) :: shoc_mix(shcol,nlev)

  ! local variables
  integer i, j, k, kk
  integer kl, ku, kb, kc, dothis, kli, kui
  real(r8) :: deep_thresh, deep_thick, cloud_thick, lstarn, thresh
  real(r8) :: cldmix, vonk, thedel, depth
  real(r8) :: omn, betdz, bbb, term, qsatt, dqsat, bet
  real(r8) :: thv_up, thv_dn, thedz, tscale, thefac, thecoef, thegam, norm
  real(r8) :: stabterm, conv_var, tkes, mmax, cldthresh
  logical lf, indexr
  real(r8) :: conv_vel(shcol,nlev)
  
  real(r8) :: numer(shcol)
  real(r8) :: denom(shcol)
  real(r8) :: cldarr(shcol) 
  real(r8) :: l_inf(shcol)
  real(r8) :: brunt2(shcol,nlev)
  
  vonk=0.35_r8
  tscale=400._r8 ! time scale set based on similarity results
  brunt2(:,:) = 0._r8
  
  ! Find length scale outside of clouds
  do k=1,nlev
    do i=1,shcol
    
      if (cldin(i,k) .eq. 0) then
        tkes=sqrt(tke(i,k))
	numer(i)=numer(i)+tkes*z_in(i,k)*dz(i)*adz(i,k)
	denom(i)=denom(i)+tkes*dz(i)*adz(i,k)
      else
        cldarr(i)=1
      endif
    
    enddo
  enddo
  
  do i=1,shcol
    if (denom(i) .gt. 0._r8) then
      l_inf(i)=0.1_r8*(numer(i)/denom(i))
    else
      l_inf(i)=100._r8
    endif
  enddo
  
  do k=1,nlev
    do i=1,shcol
    
      if (k .eq. 1) then
        kb=1
        kc=2
        thedz=adzw(i,kc)*dz(i)
      elseif (k .eq. nlev) then
        kb=nlev-1
        kc=nlev
        thedz=adzw(i,k)*dz(i)
      else
        kb=k-1
        kc=k+1
        thedz=(adzw(i,kc)+adzw(i,k))*dz(i)
      endif
    
      tkes=sqrt(tke(i,k))
      
      ! Compute the local stability
      brunt(i,k) = (ggr/basetemp) * (thetal(i,kc) - thetal(i,kb))/thedz
      
      if (brunt(i,k) .ge. 0) brunt2(i,k) = brunt(i,k)
      
      ! End of computation of local stability
      if (k .eq. 1) then
        term=600._r8 * tkes
	shoc_mix(i,k)=term+(0.4_r8*z_in(i,k)-term)*exp(-z_in(i,k)/100._r8)
      else
        shoc_mix(i,k)=min(maxlen,(2.8284*sqrt(1._r8/((1._r8/(tscale*tkes*vonk*z_in(i,k))) &
	  +(1._r8/(tscale*tkes*l_inf(i)))+0.01_r8*(brunt2(i,k)/tke(i,k)))))/0.3_r8)  
      endif      
      
    enddo  ! end i loop
  enddo ! end k loop
  
  ! Now find length scale in clouds
  cldthresh=0._r8
  
  ! determine the convective velocity scale at
  !   the top of the cloud
  
  conv_vel(:,1)=0._r8
  do k=2,nlev
    do i=1,shcol
      conv_vel(i,k) = conv_vel(i,k-1)+2.5_r8*dz(i)*adzw(i,k)*(ggr/basetemp)*wthv_sec(i,k)
    enddo
  enddo
  
  do i=1,shcol
  
    if (cldarr(i) .eq. 1) then
      
      kl=0
      ku=0
      
      do k=2,nlev-2
        
	! Look for cloud top in this column
	if (cldin(i,k) .gt. cldthresh .and. kl .eq. 0) then
          kl=k
        endif
      
        ! Look for cloud base in this column
	if (cldin(i,k) .gt. cldthresh .and. cldin(i,k+1) .le. cldthresh) then 
	  ku=k
	  conv_var=conv_vel(i,k)**(1._r8/3._r8)
	endif
	
	! Compute the mixing length for the layer just determined
	if (kl .gt. 0 .and. ku .gt. 0 .and. ku-kl .gt. 1) then
	
	  if (conv_var .gt. 0) then
	    
	    depth=(z_in(i,ku) - z_in(i,kl)) + dz(i)*adz(i,kl)
	    mmax=maxlen
	    if (z_in(i,ku) .gt. maxlen) mmax=maxlen
	    
            shoc_mix(i,kl:ku)=min(mmax,sqrt(1._r8/(((conv_var)/ &
              (depth*sqrt(tke(i,kl:ku))))**2+0.01_r8* &
              (brunt2(i,kl:ku)/tke(i,kl:ku))))/0.3_r8)	
	      
	  endif
	    
	  kl=0
	  ku=0
	
	endif 
	
      enddo ! end k loop
	
    endif ! end cldarr conditional  
	
  enddo ! end i loop 
	   
  do k=1,nlev
    do i=1,shcol
      
      shoc_mix(i,k)=min(maxlen,shoc_mix(i,k))
      shoc_mix(i,k)=max(0.1*dz(i)*adz(i,k),shoc_mix(i,k))

      shoc_mix(i,k)=min(sqrt(host_dx(i)*host_dy(i)),shoc_mix(i,k))

    enddo
  enddo 
  
  return
  
end subroutine shoc_length

!==============================================================
! Linear interpolation to get values on various grids

subroutine linear_interp(x1,x2,y1,y2,km1,km2,ncol)
  implicit none

  integer :: ncol
  integer :: km1, km2, k1, k2, i
  real(KIND=8) :: x1(km1), y1(km1,ncol)
  real(KIND=8) :: x2(km2), y2(km2,ncol)

  do i=1,ncol
    do k2=1,km2
      if( x2(k2) <= x1(1) ) then
        y2(k2,i) = y1(1,i) + (y1(2,i)-y1(1,i))*(x2(k2)-x1(1))/(x1(2)-x1(1))
      elseif( x2(k2) >= x1(km1) ) then
        y2(k2,i) = y1(km1,i) + (y1(km1,i)-y1(km1-1,i))*(x2(k2)-x1(km1))/(x1(km1)-x1(km1-1))    
      else
        do k1 = 2,km1
          if( (x2(k2)>=x1(k1-1)).and.(x2(k2)<x1(k1)) ) then
            y2(k2,i) = y1(k1-1,i) + (y1(k1,i)-y1(k1-1,i))*(x2(k2)-x1(k1-1))/(x1(k1)-x1(k1-1))
          endif
        enddo ! end k1 loop
      endif 
    enddo ! end k2 loop
  enddo ! i loop
  
  return

end subroutine linear_interp

!==============================================================
! 
! Saturation vapor pressure and mixing ratio. 
! Based on Flatau et.al, (JAM, 1992:1507) - valid for T > -80C
! sat. vapor over ice below -80C - used Murphy and Koop (2005)
! For water below -80C simply assumed esw/esi = 2.
! des/dT below -80C computed as a finite difference of es

real(r8) function esatw_shoc(t)
implicit none
real(r8) t	! temperature (K)
real(r8) a0,a1,a2,a3,a4,a5,a6,a7,a8 
data a0,a1,a2,a3,a4,a5,a6,a7,a8 /&
	6.105851_r8, 0.4440316_r8, 0.1430341e-1_r8, &
        0.2641412e-3_r8, 0.2995057e-5_r8, 0.2031998e-7_r8, &
        0.6936113e-10_r8, 0.2564861e-13_r8,-0.3704404e-15_r8/
!     	6.11239921, 0.443987641, 0.142986287e-1, &
!       0.264847430e-3, 0.302950461e-5, 0.206739458e-7, &
!       0.640689451e-10, -0.952447341e-13,-0.976195544e-15/
real(r8) dt
 dt = t-273.16_r8
if(dt.gt.-80._r8) then
 esatw_shoc = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt))))))) 
else
 esatw_shoc = 2._r8*0.01_r8*exp(9.550426_r8 - 5723.265_r8/t + 3.53068_r8*Log(t) - 0.00728332_r8*t)
end if
end     

end module shoc
