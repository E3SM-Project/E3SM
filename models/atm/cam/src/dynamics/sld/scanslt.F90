
module scanslt
!-----------------------------------------------------------------------
!
! Module to handle Semi-Lagrangian transport in the context of
! Semi-Lagrangian Dynamics.
!
!-----------------------------------------------------------------------
!
! $Id$
!
!-----------------------------------------------------------------------
   use shr_kind_mod,    only: r8 => shr_kind_r8
   use pmgrid,          only: plon, plev, plevp, plat, beglat, endlat
   use constituents,    only: pcnst
   use abortutils,      only: endrun
   use hycoef,          only: hyai, hybi, hypm, hypi, ps0, hybm, nprlev, hybd
   use perf_mod
   use cam_logfile,  only: iulog
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
   private
!
! Public interfaces
!
   public slt_alloc     ! Initialization
   public slt_initial   ! Initialization
   public slt_run_setup ! setup internal data for run
   public slt_run       ! run method
   public slt_final     ! Finalization
!
! Public data
!
  public plond, platd, i1, j1, nxpt, beglatex, endlatex, nlonex, jintmx

  integer, parameter :: nxpt       = 1     ! no.of pts outside active domain of interpolnt
  integer, parameter :: jintmx     = 2     ! # of extra lats in polar region
  integer, parameter :: plond      = plon + 1 + 2*nxpt ! slt extended domain longitude
  integer, parameter :: platd      = plat + 2*nxpt + 2*jintmx ! slt extended domain lat.
  integer, parameter :: i1         = 1 + nxpt          ! model starting longitude index
  integer, parameter :: j1         = 1 + nxpt + jintmx ! model starting latitude index
  integer, parameter :: numbnd     = nxpt + jintmx ! no.of lats passed N/S of forecast lat
  integer, parameter :: plndlv     = plond*plev    ! Length of multilevel 3-d field slice
  integer :: numlatsex  ! number of ext lats owned by a given proc
  integer :: beglatex   ! extended grid beglat
  integer :: endlatex   ! extended grid endlat
#if ( ! defined SPMD )
  parameter (numlatsex  = platd)
  parameter (beglatex   = 1)
  parameter (endlatex   = platd)
#endif
  integer :: nlonex(platd)           ! Number of longitudes for each latitude on
                                     ! extended grid
  integer, public, parameter :: plondfft = plon + 2

  public numbnd, plndlv
!
! SPMD decomposition information
!

  integer, allocatable :: cutex(:,:)        ! extended partition 
  integer :: neighs                         ! number of south neighbors to comm guardcells
  integer, allocatable :: neighs_proc(:)    ! sorted south process neighbors
  integer :: neighn                         ! number of north neighbors to comm guardcells
  integer, allocatable :: neighn_proc(:)    ! sorted north process neighbors
!
! Public types
!
  public advection_state

  type advection_state

     real(r8), pointer :: u3(:,:,:)          ! u wind
     real(r8), pointer :: v3(:,:,:)          ! v wind
     real(r8), pointer :: t3(:,:,:)          ! Temperature
     real(r8), pointer :: q3(:,:,:,:)        ! Specific humidity
     real(r8), pointer :: u3sld(:,:,:)       ! u-wind (time n) (used for advection)
     real(r8), pointer :: v3sld(:,:,:)       ! v-wind (time n) (used for advection)
     real(r8), pointer :: lnpssld(:,:,:)     ! RHS Ps term for SLD
     real(r8), pointer :: prhssld(:,:,:)     ! RHS Ps term for SLD
     real(r8), pointer :: etadot(:,:,:)      ! vertical motion (on previous time-step)

  end type advection_state

  public phigs, levkntl   ! Can become private as soon as slttraj is part of scanslt
  public engy1lat, hw1lat ! Used in scan2, restart, and linemsdyn
  public epssld
  public qfcst            ! Used in restart

  real(r8) :: hw1lat (pcnst,plat)         ! lat contribution to const. mass integral
  real(r8) :: engy1lat(plat)              ! lat contribution to total energy integral
  real(r8), parameter :: epssld = 0.2_r8  ! "epsilon" for SLD decentering algorithm
  real(r8), allocatable, target :: qfcst(:,:,:,:) ! slt forecast of moisture and constituents

!
! Global integrals
!
  public hw1, hw2, hw3, alpha              ! Needed for restart

  real(r8) :: hw1    (pcnst)               ! Pre-SLT global integral of constitent
  real(r8) :: hw2    (pcnst)               ! Post-SLT global integral of const.
  real(r8) :: hw3    (pcnst)               ! Global integral for denom. of expr. for alpha
  real(r8) :: alpha  (pcnst)               ! alpha(m) = ( hw1(m) - hw2(m) )/hw3(m)
!
! Private module data
!-----------------------------------------------------------------------
!
!------------------------------Parameters-------------------------------
!
   integer, parameter :: pmap = 20000    ! max dimension of evenly spaced vert. grid used
!                                        ! by SLT code to map the departure pts into true
!                                        ! model levels.
!
   integer :: kdpmpf  (pmap)             ! artificial full vert grid indices
   integer :: kdpmph  (pmap)             ! artificial half vert grid indices
   real(r8) :: lam    (plond,platd)      ! longitude coords of extended grid
   real(r8) :: phi    (platd)            ! latitude  coords of extended grid
   real(r8) :: dphi   (platd)            ! latitude intervals (radians)
   real(r8) :: sinlam (plond,platd)      ! sin(lam) model domain only
   real(r8) :: coslam (plond,platd)      ! cos(lam) model domain only
   real(r8) :: lbasdy (4,2,platd)        ! latitude derivative weights
   real(r8) :: lbasdz (4,2,plev)         ! vert (full levels) deriv wghts
   real(r8) :: lbassd (4,2,plevp)        ! vert (half levels) deriv wghts
   real(r8) :: lbasiy (4,2,platd)        ! Lagrange cubic interp wghts (lat.)
   real(r8) :: lbasiz (4,2,plev)         ! Lagrange cubic interp wghts (vert)
   real(r8) :: lbassi (4,2,plevp)        ! Lagrange cubic interp wghts (vert)
   real(r8) :: detai  (plevp)            ! intervals between vert half levs.
   real(r8) :: dlam   (platd)            ! longitudinal grid interval (radians)
   real(r8) :: etaint (plevp)            ! vertical coords at interfaces
   integer :: kstep                      ! Keep track of number of steps ran
   real(r8) :: gamma  (plev,plev)        ! SLD coefficient
   real(r8) :: levknt (plev,plev)        ! Counter for departure point binning statistics
   real(r8) :: levkntl(plev,plev,plat)   ! Counter for departure point binning statistics
   real(r8) :: phigs                     ! Latitude cutoff for using local geodesic algorithm

contains
!
!-----------------------------------------------------------------------
!

#if (defined SPMD)
subroutine slt_spmd()
!-----------------------------------------------------------------------
!
! Purpose:
! Setup SPMD information regarding SLT halos will be used in bndexch.
!
!-----------------------------------------------------------------------
  use spmd_utils,   only: iam, masterproc, npes
  use spmd_dyn, only: cut, dyn_npes, dyn_npes_stride

!
! Local variables
!
  integer :: lat, procid            ! Loop indices
  integer :: procid_s               ! strided process id
  integer :: neighn_minlat(plat)    ! minimum latitude in north neighbor
  integer :: neighs_maxlat(plat)    ! maximum latitude in south neighbor
!
! The extended regions are simply "numbnd" wider at each
! side. The extended region do not go beyond 1 and plat, though
!
  allocate (cutex(2,0:npes-1))
  cutex(1,0:npes-1) = 1
  cutex(2,0:npes-1) = 0
!
!  Initialization for inactive processes
!
  beglatex = 1
  endlatex = 0
  numlatsex = 0
!
  do procid=0,dyn_npes-1
     procid_s = dyn_npes_stride*procid
     cutex(1,procid_s) = cut(1,procid_s) - numbnd
     cutex(2,procid_s) = cut(2,procid_s) + numbnd
     if (iam == procid_s) then
        beglatex = cutex(1,procid_s) + numbnd
        endlatex = cutex(2,procid_s) + numbnd
        numlatsex = endlatex - beglatex + 1
     end if
  end do
!
! Determine neighbor processors needed for boundary communication.
! North first.
!
  neighn = 0
  neighn_minlat(:) = -1
  do procid=0,dyn_npes-1
     procid_s = dyn_npes_stride*procid
     if (procid_s /= iam) then
        if ((cut(1,procid_s) > cut(2,iam)) .and. &
            (cut(1,procid_s) <= cut(2,iam)+numbnd)) then
           neighn_minlat(cut(1,procid_s)) = procid_s
           neighn = neighn + 1
        endif
     endif
  enddo
!
! Sort north processes by increasing latitude
!
  allocate (neighn_proc (neighn))
  neighn = 0
  do lat=1,plat
     if (neighn_minlat(lat) /= -1) then
        neighn = neighn + 1
        neighn_proc(neighn) = neighn_minlat(lat)
     endif
  enddo
!
! South next.
!
  neighs = 0
  neighs_maxlat(:) = -1
  do procid=0,dyn_npes-1
     procid_s = dyn_npes_stride*procid
     if (procid_s /= iam) then
        if ((cut(2,procid_s) < cut(1,iam)) .and. &
            (cut(2,procid_s) >= cut(1,iam)-numbnd)) then
           neighs_maxlat(cut(2,procid_s)) = procid_s
           neighs = neighs + 1
        endif
     endif
  enddo
!
! Sort south processes by decreasing latitude
!
  allocate (neighs_proc (neighs))
  neighs = 0
  do lat=plat,1,-1
     if (neighs_maxlat(lat) /= -1) then
        neighs = neighs + 1
        neighs_proc(neighs) = neighs_maxlat(lat)
     endif
  enddo
!
  if (masterproc) then
     write(iulog,*)'-----------------------------------------'
     write(iulog,*)'Number of lats passed north & south = ',numbnd
     write(iulog,*)'Node  Extended Partition'
     write(iulog,*)'-----------------------------------------'
     do procid=0,npes-1
        write(iulog,200) procid,cutex(1,procid), cutex(2,procid)
200     format(i3,4x,i3,'-',i3)
     end do
  end if
! write(iulog,*)'iam=',iam,'Number of south neighbors needed for bndry exchange = ',neighs
! write(iulog,*)'iam=',iam,'Number of north neighbors needed for bndry exchange = ',neighn

end subroutine slt_spmd
#endif

!
!-----------------------------------------------------------------------
!

subroutine slt_alloc()
!-----------------------------------------------------------------------
!
! Purpose:
! slt initialization
!
!-----------------------------------------------------------------------
  use infnan,       only: nan, assignment(=)

  nlonex(:) = huge(1)
  allocate (qfcst(plon,plev,pcnst,beglat:endlat))
  qfcst (:,:,:,:) = nan
end subroutine slt_alloc

!
!-----------------------------------------------------------------------
!

subroutine slt_initial( detam, gw, cwava, etamid, adv_state )
!-----------------------------------------------------------------------
!
! Purpose:
! slt initialization
!
!-----------------------------------------------------------------------
  use pmgrid,       only: plev, plevp
  use commap,       only: t0, clat
  use rgrid,        only: nlon
  use physconst,    only: cappa, gravit
  use prognostics,  only: u3, v3, n3, div, dpsl, dpsm, ps, ed1
  use time_manager, only: is_first_step
  use sld_control_mod, only: pdela
!
! Arguments
!
   real(r8), intent(in) :: etamid (plev)             ! vertical coords at midpoints
   real(r8), intent(out) :: gw     (plat)            ! Gaussian weights
   real(r8), intent(out) :: detam  (plev)            ! intervals between vert full levs.
   real(r8), intent(out) :: cwava  (plat)            ! weight applied to global integrals
   type(advection_state), intent(out) :: adv_state   ! Advection state
!
! Local variables
!
  real(r8) :: coslat(plon)              ! cosine of latitude
  real(r8) :: rcoslat(plon)             ! 1 over cosine of latitude
  real(r8) :: hyad   (plev)             ! del (A)
  integer  :: i, j, k, kk, l, c, lat    ! Indices
  real(r8) :: tmp1                      ! temp variable
  real(r8) :: pdel(plon,plev)           ! Pressure difference over layer
  real(r8) :: pint(plon,plevp)          ! Pressure at interfaces
  real(r8) :: pmid(plon,plev)           ! Pressure at layer mid-points
!-----------------------------------------------------------------------
#if ( defined SPMD )
  call slt_spmd()
#endif

  call alloc_advstate( adv_state )
!
! Eta coordinates on interfaces: Used for calculation etadot vertical velocity
!
  do k=1,plevp
      etaint(k) = hyai(k) + hybi(k)
   end do
!
! Initialize matrix gamma (used in sld T computation)
!
   gamma(:,:) = 0._r8
   do k = 1,plev
      tmp1 = cappa*t0(k)*hypi(plevp)/hypm(k)
      gamma(k,k) = 0.5_r8*tmp1
      do l=1,k-1
         gamma(l,k) = tmp1
      end do
   end do
!
! Set slt common block variables
!
   call grdini(pmap    ,etamid  ,etaint  ,gravit  ,dlam    , &
               lam     ,phi     ,dphi    ,gw      ,sinlam  , &
               coslam  ,lbasdy  ,lbasdz  ,lbassd  ,lbasiy  , &
               lbasiz  ,lbassi  ,detam   ,detai   ,kdpmpf  , &
               kdpmph  ,cwava   ,phigs   )
!
! Compute pdel from "A" portion of hybrid vertical grid
!
   do k=1,plev
      hyad(k) = hyai(k+1) - hyai(k)
   end do
   do k=1,plev
      do i=1,plon
         pdela(i,k) = hyad(k)*ps0
      end do
   end do
   if (is_first_step()) then
      do lat=beglat,endlat
         do i=1,nlon(lat)
            coslat(i) = cos(clat(lat))
            rcoslat(i) = 1._r8/coslat(i)
         end do
!    
! Set current time pressure arrays for model levels etc.
!
         call plevs0(nlon(lat), plon, plev, ps(1,lat,n3), pint, pmid, pdel)
!
! Compute appropriate (1/ps)etadot(dp/deta)
!
         call etadt0 (lat, nlon(lat), &
                      rcoslat ,div(1,1,lat,n3), u3(1,1,lat,n3), v3(1,1,lat,n3), dpsl(1,lat), &
                      dpsm(1,lat), pdel, ps(1,lat,n3), ed1(1,1,lat))
      end do
   end if
!
!----------------------------------------------------------
! Initialize departure point bin arrays to 0.
!----------------------------------------------------------
!
   kstep = 0
   do lat = beglat,endlat
      do k = 1,plev
         do kk = 1,plev
            levkntl(kk,k,lat) = 0._r8
         end do
      end do
   end do

end subroutine slt_initial
!
!-----------------------------------------------------------------------
!

subroutine slt_run( ztodt   ,detam   ,cwava   ,etamid  ,          &
                    grlps1  ,grlps2  ,grt1    ,grt2    ,grq1    , &
                    grq2    ,grfu1   ,grfu2   ,grfv1   ,grfv2   , &
                    grfu    ,grfv    ,lnpssld ,prhssld ,adv_state )
!-----------------------------------------------------------------------
!
! Purpose:
! slt run method
!
!-----------------------------------------------------------------------
  use pmgrid,       only: plat, plon
  use comspe,       only: maxm
  use constituents, only: pcnst
  use rgrid,        only: nlon
  use prognostics,  only: parrsld, tarrsld, u3, v3, etadot, n3

#if ( defined SPMD )
  use mpishorthand, only: mpicom
#endif
!

!------------------------------Arguments--------------------------------
!
  real(r8), intent(in)   :: ztodt                     ! twice the timestep unless nstep=0
  real(r8), intent(in)   :: detam  (plev)             ! delta eta at levels
  real(r8), intent(in)   :: cwava  (plat)             ! weight for global water vapor int.
  real(r8), intent(in)   :: etamid (plev)             ! eta at levels

  real(r8), intent(out)   :: grlps1(2*maxm,plat/2)      ! ------------------------------
  real(r8), intent(out)   :: grlps2(2*maxm,plat/2)      ! |
  real(r8), intent(out)   :: grt1  (2*maxm,plev,plat/2) ! |
  real(r8), intent(out)   :: grt2  (2*maxm,plev,plat/2) ! |
  real(r8), intent(out)   :: grq1  (2*maxm,plev,plat/2) ! |
  real(r8), intent(out)   :: grq2  (2*maxm,plev,plat/2) ! |- see quad for definitions
  real(r8), intent(out)   :: grfu1 (2*maxm,plev,plat/2) ! |
  real(r8), intent(out)   :: grfu2 (2*maxm,plev,plat/2) ! |
  real(r8), intent(out)   :: grfv1 (2*maxm,plev,plat/2) ! |
  real(r8), intent(out)   :: grfv2 (2*maxm,plev,plat/2) ! ------------------------------
  real(r8), intent(inout) :: grfu  (plon,plev,beglat:endlat)   ! nonlinear term - u momentum eqn
  real(r8), intent(inout) :: grfv  (plon,plev,beglat:endlat)   ! nonlinear term - v momentum eqn
  real(r8), intent(in) :: lnpssld(plon,plev,beglat:endlat)
  real(r8), intent(in) :: prhssld(plon,plev,beglat:endlat)
  type(advection_state), intent(inout) :: adv_state
!---------------------------Local workspace-----------------------------
!
  integer, parameter :: iter=1  ! number of iterations to be used in dep pt calc.

  integer :: c                  ! constituent index  
  integer :: irow               ! latitude pair index
  integer :: lat,latn,lats      ! latitude indices
  integer :: nlon_fft_in        ! FFT work array inner dimension
  integer :: nlon_fft_out       ! FFT work array inner dimension
  real(r8) :: onepeps           ! 1 + epssld
  real(r8) :: dtr               ! 1/dt
!
! FFT buffers
!
  real(r8), allocatable:: fftbuf_in(:,:,:,:)          ! fftbuf_in(nlon_fft_in,plev,5,beglat:endlat) 
  real(r8), allocatable:: fftbuf_out(:,:,:,:)         ! fftbuf_out(nlon_fft_out,plev,5,plat)

  call t_startf ('da_coupling')
  call da_coupling( adv_state, cwava, lnpssld, prhssld )
  call t_stopf ('da_coupling')

#if ( defined SPMD )
!
! Communicate boundary information
!
  call t_barrierf ('sync_bndexch', mpicom)
  call t_startf ('bndexch')
  call bndexch( adv_state )
  call t_stopf ('bndexch')
#endif

  onepeps = 1._r8 + epssld
  dtr     = 2._r8/(ztodt*onepeps)
!
! Fill latitude/longitude extensions of constituents and dynamics terms
!
  if (beglatex .le. endlatex) then
!$OMP PARALLEL DO PRIVATE (C)
     do c = 1, pcnst
        call extys(1    ,plev                      ,adv_state%q3(:,:,beglatex:endlatex,c),    1)
     end do
     call extyv(1       ,plev    ,coslam  ,sinlam  ,adv_state%u3sld  (:,:,beglatex:endlatex     ), &
                                                    adv_state%v3sld  (:,:,beglatex:endlatex     ))
     call extys(1       ,plevp                     ,adv_state%etadot (:,:,beglatex:endlatex), 1)
!
     call extyv(1       ,plev    ,coslam  ,sinlam  ,adv_state%u3     (:,:,beglatex:endlatex), &
                                                    adv_state%v3     (:,:,beglatex:endlatex))
     call extys(1       ,plev                      ,adv_state%t3     (:,:,beglatex:endlatex), 1)
     call extys(1       ,plev                      ,adv_state%lnpssld(:,:,beglatex:endlatex     ), 1)
     call extys(1       ,plev                      ,adv_state%prhssld(:,:,beglatex:endlatex     ), 1)
!
!$OMP PARALLEL DO PRIVATE (C)
     do C = 1, pcnst
        call extx (1    ,plev                      ,adv_state%q3     (:,:,beglatex:endlatex,c),    1)
     end do
     call extx (1       ,plev                      ,adv_state%u3sld  (:,:,beglatex:endlatex     ), 1)
     call extx (1       ,plev                      ,adv_state%v3sld  (:,:,beglatex:endlatex     ), 1)
     call extx (1       ,plevp                     ,adv_state%etadot (:,:,beglatex:endlatex), 1)
!
     call extx (1       ,plev                      ,adv_state%u3     (:,:,beglatex:endlatex), 1)
     call extx (1       ,plev                      ,adv_state%v3     (:,:,beglatex:endlatex), 1)
     call extx (1       ,plev                      ,adv_state%t3     (:,:,beglatex:endlatex), 1)
     call extx (1       ,plev                      ,adv_state%lnpssld(:,:,beglatex:endlatex     ), 1)
     call extx (1       ,plev                      ,adv_state%prhssld(:,:,beglatex:endlatex     ), 1)
  endif
!
! Begin SLT interpolation
!
!
  call t_startf ('scanslt_alloc')
  nlon_fft_in = plondfft
  allocate(fftbuf_in(nlon_fft_in,plev,5,beglat:endlat))

#if ( defined SPMD )
#ifdef NEC_SX
  nlon_fft_out = 2*maxm + 1
#else
  nlon_fft_out = 2*maxm
#endif
  allocate(fftbuf_out(nlon_fft_out,plev,5,plat))
#else
  nlon_fft_out = 1
  allocate(fftbuf_out(1,1,1,1))
#endif
  call t_stopf ('scanslt_alloc')

  call t_startf ('scanslt_bft')
!$OMP PARALLEL DO PRIVATE(LAT)

  do lat = beglat,endlat
     call scanslt_bft(ztodt   ,lat      ,dtr     ,iter       ,       &
                      detam   ,etamid   ,grfu    ,grfv       ,       &
                      tarrsld ,parrsld  ,u3(1,1,lat,n3)      ,       &
                      v3(1,1,lat,n3)    ,etadot(1,1,lat)     ,       &
                      nlon(lat)         ,nlon_fft_in         ,       &
                      fftbuf_in(1,1,1,lat), adv_state )
  end do                    ! end lat-loop
  call t_stopf ('scanslt_bft')
!
  call t_startf ('scanslt_fft')
  call scanslt_fft(nlon_fft_in,nlon_fft_out,fftbuf_in,fftbuf_out)
  call t_stopf ('scanslt_fft')
!
  call t_startf ('scanslt_aft')
!$OMP PARALLEL DO PRIVATE (IROW, LATN, LATS)

  do irow=1,plat/2

      lats = irow
      latn = plat - irow + 1
#if ( defined SPMD )
      call scanslt_aft(irow, nlon_fft_out, &
                       fftbuf_out(1,1,1,lats), fftbuf_out(1,1,1,latn), &
                       grlps1(1,irow)  ,grlps2(1,irow) , &
                       grt1(1,1,irow)  ,grt2(1,1,irow) , &
                       grq1(1,1,irow)  ,grq2(1,1,irow) , &
                       grfu1(1,1,irow) ,grfu2(1,1,irow), &
                       grfv1(1,1,irow) ,grfv2(1,1,irow)  )
#else
      call scanslt_aft(irow, nlon_fft_in, &
                       fftbuf_in(1,1,1,lats), fftbuf_in(1,1,1,latn), &
                       grlps1(1,irow)  ,grlps2(1,irow) , &
                       grt1(1,1,irow)  ,grt2(1,1,irow) , &
                       grq1(1,1,irow)  ,grq2(1,1,irow) , &
                       grfu1(1,1,irow) ,grfu2(1,1,irow), &
                       grfv1(1,1,irow) ,grfv2(1,1,irow)  )
#endif
  end do                    ! end irow-loop
  call t_stopf ('scanslt_aft')

  call t_startf ('scanslt_dealloc')
  deallocate(fftbuf_in)
  deallocate(fftbuf_out)
  call t_stopf ('scanslt_dealloc')

  call t_startf ('ad_coupling')
  call ad_coupling( adv_state )
  call t_stopf ('ad_coupling')

end subroutine slt_run
!
!-----------------------------------------------------------------------
!

subroutine alloc_advstate( adv_state)
!-----------------------------------------------------------------------
!
! Purpose:
! Allocate advection_state structure
!
!-----------------------------------------------------------------------
    use infnan, only: nan, assignment(=)
    type(advection_state), intent(out) :: adv_state

    character(len=*), parameter :: sub = 'alloc_advstate'

!pw    if ( (beglatex == 0) .or. (endlatex == 0) ) then
!pw       call endrun(sub//': extended grid begining and ending latitude have NOT been set yet!')
!pw    end if
    allocate (adv_state%u3     (plond,plev, beglatex:endlatex)      )
    allocate (adv_state%v3     (plond,plev, beglatex:endlatex)      )
    allocate (adv_state%t3     (plond,plev, beglatex:endlatex)      )
    allocate (adv_state%q3     (plond,plev, beglatex:endlatex,pcnst))
    allocate (adv_state%u3sld  (plond,plev, beglatex:endlatex)      )
    allocate (adv_state%v3sld  (plond,plev, beglatex:endlatex)      )
    allocate (adv_state%lnpssld(plond,plev, beglatex:endlatex)      )
    allocate (adv_state%prhssld(plond,plev, beglatex:endlatex)      )
    allocate (adv_state%etadot (plond,plevp,beglatex:endlatex)      )

    if (beglatex .le. endlatex) then
       adv_state%u3(:,:,:)      = nan
       adv_state%v3(:,:,:)      = nan
       adv_state%t3(:,:,:)      = nan
       adv_state%q3(:,:,:,:)    = nan
       adv_state%u3sld(:,:,:)   = nan
       adv_state%v3sld(:,:,:)   = nan
       adv_state%lnpssld(:,:,:) = nan
       adv_state%prhssld(:,:,:) = nan
       adv_state%etadot(:,:,:)  = nan
    end if
end subroutine alloc_advstate

!
!-----------------------------------------------------------------------
!

subroutine dealloc_advstate( adv_state )
!-----------------------------------------------------------------------
!
! Purpose:
! De-allocate advection_state structure
!
!-----------------------------------------------------------------------
    type(advection_state), intent(inout) :: adv_state

    deallocate (adv_state%u3)
    deallocate (adv_state%v3)
    deallocate (adv_state%t3)
    deallocate (adv_state%q3)
    deallocate (adv_state%u3sld)
    deallocate (adv_state%v3sld)
    deallocate (adv_state%lnpssld)
    deallocate (adv_state%prhssld)
    deallocate (adv_state%etadot)

end subroutine dealloc_advstate

!
!-----------------------------------------------------------------------
!

subroutine da_coupling ( adv_state, cwava, lnpssld, prhssld )
!-----------------------------------------------------------------------
!
! Purpose:
! Dynamics to advection coupling
!
!-----------------------------------------------------------------------
   use prognostics,  only: u3, v3, t3, q3, etadot, n3, n3m1, ps
   use time_manager, only: is_first_step
   use rgrid,        only: nlon
   use commap,       only: w
   use qmassa,      only: qmassarun

   type(advection_state), intent(inout) :: adv_state
   real(r8), intent(in) :: cwava  (plat)                    ! weight applied to global integrals
   real(r8), intent(in) :: lnpssld(plon,plev,beglat:endlat)
   real(r8), intent(in) :: prhssld(plon,plev,beglat:endlat)
   real(r8) pmid (plon,plev)          ! pressure at model levels
   real(r8) pint (plon,plevp)         ! pressure at interfaces
   real(r8) pdel (plon,plev)          ! pressure difference between
   integer :: i, j, k, c, jcen, icen  ! Indices
   integer :: irow                    ! latitude pair index

   if(is_first_step()) then
!$OMP PARALLEL DO PRIVATE (I,K,J,JCEN,ICEN)
      do j = beglat, endlat
         jcen = j1 - 1 + j
         do k = 1, plevp
            do i = 1, nlon(j)
               icen = i1 + i - 1
               adv_state%etadot(icen,k,jcen) = etadot(i,k,j)
            end do
         end do
      end do
   end if
!$OMP PARALLEL DO PRIVATE (I,K,J,ICEN,JCEN,IROW,PINT,PMID,PDEL)
   do j = beglat, endlat
      jcen = j1 - 1 + j
      if( j <= plat/2 ) then
         irow = j
      else
         irow = plat + 1 - j
      end if

      do k = 1, plev
         do i = 1, nlon(j)
            icen = i1 + i - 1
            adv_state%u3(icen,k,jcen)      = u3(i,k,j,n3m1)
            adv_state%v3(icen,k,jcen)      = v3(i,k,j,n3m1)
            adv_state%t3(icen,k,jcen)      = t3(i,k,j,n3m1)
            adv_state%lnpssld(icen,k,jcen) = lnpssld(i,k,j)
            adv_state%prhssld(icen,k,jcen) = prhssld(i,k,j)
         end do
      end do
!
! Calculate mass of moisture in field being advected by slt.
!
      call plevs0(nlon(j), plon, plev, ps(1,j,n3), pint    ,pmid    ,pdel)
      call qmassarun(cwava(j)   , w(irow) , q3(1,1,1,j,n3), pdel    ,&
                  hw1lat(1,j), nlon(j),  q3(1,1,1,j,n3m1), j)
!
! The modified etadot @n3m1 will be used later for trajectory calculation in SCANSLT
!
      do k = 1, plevp
         do i = 1, nlon(j)
            icen = i1 + i - 1
            adv_state%etadot(icen,k,jcen) = 2._r8*etadot(i,k,j) - &
                                            adv_state%etadot(icen,k,jcen)
         end do
      end do
   end do
!$OMP PARALLEL DO PRIVATE (I,K,J,C,ICEN,JCEN)
   do c = 1, pcnst
      do j = beglat, endlat
         jcen = j1 - 1 + j
         do k = 1, plev
            do i = 1, nlon(j)
               icen = i1 + i - 1
               adv_state%q3(icen,k,jcen,c) = q3(i,k,c,j,n3)
            end do
         end do
      end do
   end do

end subroutine da_coupling

!
!-----------------------------------------------------------------------
!

subroutine slt_run_setup( adv_state )
!-----------------------------------------------------------------------
!
! Purpose:
! Setup internal advection state data for running scanslt.
!
!-----------------------------------------------------------------------
   use prognostics,  only: u3, v3, n3, n3m1, etadot
   use rgrid,        only: nlon
   use time_manager, only: is_first_restart_step
   type(advection_state), intent(inout) :: adv_state

   integer :: i, j, k, icen, jcen    ! Indices
!
! Compute U/V for time n + 1/2
! This will be used later for trajectory calculation
!

!$OMP PARALLEL DO PRIVATE (I,K,J,JCEN,ICEN)
   do j = beglat, endlat
      jcen = j1 - 1 + j
      do k = 1, plev
         do i = 1, nlon(j)
            icen = i1 + i - 1
            adv_state%u3sld(icen,k,jcen) = 2._r8*u3(i,k,j,n3) - u3(i,k,j,n3m1)
            adv_state%v3sld(icen,k,jcen) = 2._r8*v3(i,k,j,n3) - v3(i,k,j,n3m1)
         end do
      end do
   end do
   if(is_first_restart_step()) then
!$OMP PARALLEL DO PRIVATE (I,K,J,ICEN,JCEN)
      do j = beglat, endlat
         jcen = j1 - 1 + j
         do k = 1, plevp
            do i = 1, nlon(j)
               icen = i1 + i - 1
               adv_state%etadot(icen,k,jcen) = etadot(i,k,j)
            end do
         end do
      end do
   end if
end subroutine slt_run_setup


!
!-----------------------------------------------------------------------
!

subroutine ad_coupling ( adv_state )
!-----------------------------------------------------------------------
!
! Purpose:
! Advection to dynamics coupling
!
!-----------------------------------------------------------------------
   use prognostics,  only: u3, v3, t3, q3, n3, n3m1, etadot
   use rgrid,        only: nlon
   type(advection_state), intent(inout) :: adv_state

   integer :: i, j, k, c, icen, jcen

!$OMP PARALLEL DO PRIVATE (I,K,J,ICEN,JCEN)
   do j = beglat, endlat
      jcen = j1 - 1 + j
      do k = 1, plev
         do i = 1, nlon(j)
            icen = i1 + i - 1
            u3(i,k,j,n3m1) = adv_state%u3(icen,k,jcen)
            v3(i,k,j,n3m1) = adv_state%v3(icen,k,jcen)
            t3(i,k,j,n3m1) = adv_state%t3(icen,k,jcen)
         end do
      end do
      do k = 1, plevp
         do i = 1, nlon(j)
            icen = i1 + i - 1
            adv_state%etadot(icen,k,jcen) = etadot(i,k,j)
         end do
      end do
   end do
!$OMP PARALLEL DO PRIVATE (I,K,J,C,ICEN,JCEN)
   do c = 1, pcnst
      do j = beglat, endlat
         jcen = j1 - 1 + j
         do k = 1, plev
            do i = 1, nlon(j)
               icen = i1 + i - 1
               q3(i,k,c,j,n3) = adv_state%q3(icen,k,jcen,c)
            end do
         end do
      end do
   end do
!
! Increment kstep
!
   kstep = kstep + 1

end subroutine ad_coupling

!
!-----------------------------------------------------------------------
!

subroutine scanslt_bft (ztodt   ,lat     ,dtr     ,iter    ,detam   , &
                        etamid  ,grfu    ,grfv    ,tarrsld ,parrsld , &
                        u3      ,v3      ,etadot  ,nlon    ,nlon_fft, &
                        fftbuf  ,adv_state )
!-----------------------------------------------------------------------
!
! Purpose:
! Interpolate terms for semi-lagrangian transport and SLD dynamics.
! One latitude slice only
!
! Author:  J. Olson
!
!-----------------------------------------------------------------------

  use rgrid,       only: nmmax
  use pspect,      only: pmmax
  use commap,      only: w, clat, t0
  use physconst,   only: cappa, ra

!------------------------------Arguments--------------------------------
!
  real(r8), intent(in)   :: ztodt                 ! twice the time step unless nstep = 0
  integer , intent(in)   :: lat                   ! latitude index
  real(r8), intent(in)   :: dtr                   ! 1/dt
  integer , intent(in)   :: iter                  ! number of iterations for trajectory
  real(r8), intent(in)   :: detam   (plev)        ! delta eta at levels
  real(r8), intent(in)   :: etamid  (plev)        ! eta at levels
  real(r8), intent(in)   :: grfu    (plon,plev,beglat:endlat) ! nonlinear term - u momentum eqn
  real(r8), intent(in)   :: grfv    (plon,plev,beglat:endlat) ! nonlinear term - v momentum eqn
  real(r8), intent(in)   :: tarrsld (plon,plev,beglat:endlat)   ! T  at arr. pt.(SLD)
  real(r8), intent(inout):: parrsld (plon,plev,beglat:endlat)   ! Ps at arr. pt.(SLD)
  real(r8), intent(in)   :: u3      (plon,plev)                 ! u-wind at current time
  real(r8), intent(in)   :: v3      (plon,plev)                 ! v-wind at current time
  real(r8), intent(in)   :: etadot  (plon,plevp)                ! vertical motion at current time

  integer , intent(in)   :: nlon                                ! # of longitudes
  integer, intent(in)    :: nlon_fft              ! first dimension of FFT work array
  real(r8), intent(out)  :: fftbuf(nlon_fft,plev,5)   ! buffer used for in-place FFTs
  type(advection_state), intent(in) :: adv_state
!
!---------------------------Local workspace-----------------------------
!
  integer i                    ! index
  integer k                    ! index
  integer l                    ! index
  integer m                    ! constituent index
  integer irow                 ! N/S latitude pair index
  integer jcen                 ! lat index (extended grid)
!                              ! of forecast
  real(r8) qtmp (plond,plev,beglatex:endlatex)  ! Temporary for q array
  real(r8) fdp  (plon,plev,2)  ! interpolant
  real(r8) lamdp(plon,plev)    ! x-coord of dep pt
  real(r8) phidp(plon,plev)    ! y-coord of dep pt
  real(r8) sigdp(plon,plev)    ! z-coord of dep pt

  integer idp   (plon,plev,4)  ! zonal      dep point index
  integer jdp   (plon,plev)    ! meridional dep point index
  integer kdp   (plon,plev)    ! vertical   dep point index
  integer kkdp  (plon,plev)    ! index of z-coordinate of dep pt (alt)

  real(r8) xl   (plon,plev,4)  ! weight for x-interpolants (left)
  real(r8) xr   (plon,plev,4)  ! weight for x-interpolants (right)
  real(r8) wgt1x(plon,plev,4)  ! weight for x-interpolants (Lag Cubic)
  real(r8) wgt2x(plon,plev,4)  ! weight for x-interpolants (Lag Cubic)
  real(r8) wgt3x(plon,plev,4)  ! weight for x-interpolants (Lag Cubic)
  real(r8) wgt4x(plon,plev,4)  ! weight for x-interpolants (Lag Cubic)
  real(r8) hl   (plon,plev,4)  ! weight for x-interpolants (Hermite)
  real(r8) hr   (plon,plev,4)  ! weight for x-interpolants (Hermite)
  real(r8) dhl  (plon,plev,4)  ! weight for x-interpolants (Hermite)
  real(r8) dhr  (plon,plev,4)  ! weight for x-interpolants (Hermite)

  real(r8) ys   (plon,plev)    ! weight for y-interpolants (south)
  real(r8) yn   (plon,plev)    ! weight for y-interpolants (north)
  real(r8) wgt1y(plon,plev)    ! weight for y-interpolants (Lag Cubic)
  real(r8) wgt2y(plon,plev)    ! weight for y-interpolants (Lag Cubic)
  real(r8) wgt3y(plon,plev)    ! weight for y-interpolants (Lag Cubic)
  real(r8) wgt4y(plon,plev)    ! weight for y-interpolants (Lag Cubic)
  real(r8) hs   (plon,plev)    ! weight for y-interpolants (Hermite)
  real(r8) hn   (plon,plev)    ! weight for y-interpolants (Hermite)
  real(r8) dhs  (plon,plev)    ! weight for y-interpolants (Hermite)
  real(r8) dhn  (plon,plev)    ! weight for y-interpolants (Hermite)
  real(r8) rdphi(plon,plev)    ! reciprocal of y-interval

  real(r8) wgt1z(plon,plev)    ! weight for z-interpolants (Lag Cubic)
  real(r8) wgt2z(plon,plev)    ! weight for z-interpolants (Lag Cubic)
  real(r8) wgt3z(plon,plev)    ! weight for z-interpolants (Lag Cubic)
  real(r8) wgt4z(plon,plev)    ! weight for z-interpolants (Lag Cubic)
  real(r8) hb   (plon,plev)    ! weight for z-interpolants (Hermite)
  real(r8) ht   (plon,plev)    ! weight for z-interpolants (Hermite)
  real(r8) dhb  (plon,plev)    ! weight for z-interpolants (Hermite)
  real(r8) dht  (plon,plev)    ! weight for z-interpolants (Hermite)
  real(r8) rdz  (plon,plev)    ! reciprocal of z-interval
  real(r8) zt   (plon,plev)    ! top vertical interpolation weight 
  real(r8) zb   (plon,plev)    ! bot vertical interpolation weight 

  real(r8) lampr(plon,plev)    ! trajectory increment (x-direction)
  real(r8) phipr(plon,plev)    ! trajectory increment (y-direction)
  real(r8) upr  (plon,plev)    ! interpolated u field (local geodesic)
  real(r8) vpr  (plon,plev)    ! interpolated v field (local geodesic)
!
  real(r8) pd   (plon)             ! RHS term for Ps and (1/ps)etadot(dp/deta) 
  real(r8) pdsum(plon)             ! RHS term for Ps and (1/ps)etadot(dp/deta) 
  real(r8) pd1  (plon)             ! RHS term for Ps and (1/ps)etadot(dp/deta) 
  real(r8) pdsm1(plon)             ! RHS term for Ps and (1/ps)etadot(dp/deta) 
  real(r8) pa   (plon)             ! RHS term for Ps and (1/ps)etadot(dp/deta) 
  real(r8) pasum(plon)             ! RHS term for Ps and (1/ps)etadot(dp/deta) 
  real(r8) coslat                  ! cos(latitude)
  real(r8) tmp1                    ! temp space
!
  logical limdrh                   ! horizontal derivative limiter flag
  logical limdrv                   ! vertical   derivative limiter flag
  logical lhrzint                  ! horizontal interp flag
  logical lvrtint                  ! vertical   interp flag
  logical lhrzwgt                  ! flag to compute horizontal weights
  logical lvrtwgt                  ! flag to compute vertical   weights
!
  real(r8) grfulat(plon,plev)      ! non-linear terms for u-momentum 
  real(r8) grfvlat(plon,plev)      ! non-linear terms for u-momentum 
  real(r8) grtlat (plon,plev)      ! RHS of T-eqn
  real(r8) grqlat (plon,plev)      ! q
  real(r8) grpslat(plon)           ! RHS of Ps-eqn
!
  integer, parameter :: fudex = 1  ! indices into fftbuf 
  integer, parameter :: fvdex = 2
  integer, parameter :: tdex  = 3
  integer, parameter :: qdex  = 4
  integer, parameter :: psdex = 5
!
!-----------------------------------------------------------------------
!
  if(lat.le.plat/2) then
     irow = lat
  else
     irow = plat + 1 - lat
  end if
  jcen = j1 - 1 + lat
  coslat = cos(clat(lat))
!
! Initial guess for trajectory midpoints in spherical coords.
! Use arrival points as initial guess for trajectory midpoints.
!
  do k=1,plev
     do i=1,nlon
        phidp(i,k) = clat(lat)
        sigdp(i,k) = etamid(k)
     end do
  end do
!        
! Offset bottom level departure point first guess by epsilon
!
  do i = 1,nlon
     sigdp(i,plev) = sigdp(i,plev)*(1._r8 - 10._r8*epsilon(sigdp))
  end do
!
! Loop through latitudes producing departure point calculation
!
  call slttraj(pmap             ,jcen    ,lat     ,ztodt   ,ra                , &
               iter             ,lam     ,phi     ,dphi    ,etamid            , &
               etaint           ,detam   ,detai   ,lbasiy  ,lbasiz            , &
               lbassi           ,kdpmpf  ,kdpmph  ,idp     ,jdp               , &
               kdp              ,kkdp    ,xl      ,xr      ,wgt1x             , &
               wgt2x            ,wgt3x   ,wgt4x   ,hl      ,hr                , &
               dhl              ,dhr     ,ys      ,yn      ,wgt1y             , &
               wgt2y            ,wgt3y   ,wgt4y   ,hs      ,hn                , &
               dhs              ,dhn     ,rdphi   ,wgt1z   ,wgt2z             , &
               wgt3z            ,wgt4z   ,hb      ,ht      ,dhb               , &
               dht              ,rdz     ,lampr   ,phipr   ,upr               , &
               vpr              ,lamdp   ,phidp   ,sigdp   ,adv_state%u3      , &
               adv_state%v3     ,adv_state%u3sld  ,adv_state%v3sld   , &
               adv_state%etadot ,u3      ,v3      ,etadot            , &
               dlam    ,nlon    )
!
! Compute constituent forecast
!
  lhrzwgt = .true.
  lvrtwgt = .true.
  lhrzint = .true.
  lvrtint = .true.
  limdrh  = .true.
  limdrv  = .true.
  call bandij (dlam    ,phi     ,lamdp   ,phidp   ,idp     , &
               jdp     ,nlon    )
  call kdpfnd (plev    ,pmap    ,etamid  ,sigdp   ,kdpmpf  , &
               kdp     ,nlon    )
  call sltwgts(limdrh  ,limdrv  ,lhrzwgt ,lvrtwgt ,plev    , &
               idp     ,jdp     ,kdp     ,lam     ,phi     , &
               etamid  ,dphi    ,detam   ,lamdp   ,phidp   , &
               sigdp   ,lbasiy  ,lbasiz  ,kkdp    ,xl      , &
               xr      ,wgt1x   ,wgt2x   ,wgt3x   ,wgt4x   , &
               hl      ,hr      ,dhl     ,dhr     ,ys      , &
               yn      ,wgt1y   ,wgt2y   ,wgt3y   ,wgt4y   , &
               hs      ,hn      ,dhs     ,dhn     ,rdphi   , &
               wgt1z   ,wgt2z   ,wgt3z   ,wgt4z   ,hb      , &
               ht      ,dhb     ,dht     ,rdz     ,zt      , &
               zb      ,nlon    )
  do m = 1,pcnst
     call sltint (jcen    ,adv_state%q3(:,:,beglatex:endlatex,m),lam     ,rdphi   , &
                  rdz     ,lbasdy  ,lbasdz  ,xl      ,xr      , &
                  wgt1x   ,wgt2x   ,wgt3x   ,wgt4x   ,hl      , &
                  hr      ,dhl     ,dhr     ,ys      ,yn      , &
                  wgt1y   ,wgt2y   ,wgt3y   ,wgt4y   ,hs      , &
                  hn      ,dhs     ,dhn     ,wgt1z   ,wgt2z   , &
                  wgt3z   ,wgt4z   ,hb      ,ht      ,dhb     , &
                  dht     ,idp     ,jdp     ,kdp     ,kkdp    , &
                  lhrzint ,lvrtint ,limdrh  ,limdrv  ,qfcst(1,1,m,lat), &
                  nlon  )
     if ( m == 1 ) grqlat(:nlon,:plev) = qfcst(:nlon,:plev,m,lat)
  end do
!
! Accumulate P-interpolants into T equation
!
  do i = 1,nlon
     grpslat(i) = 0._r8
     pasum  (i) = 0._r8
     pa     (i) = 0._r8
  end do
  do k = 1,plev
     do i = 1,nlon
        grtlat(i,k) = tarrsld (i,k,lat)
     end do
  end do
  do k = 1,plev
     do l = 1,k
        do i = 1,nlon
           grtlat(i,k) = grtlat(i,k) + parrsld(i,l,lat)*gamma(l,k)
        end do
     end do
  end do
!
! Accumulate Ps interpolants in 3-D array
!
  do k = 1,plev
     do i = 1,nlon
        grpslat(i) = grpslat(i) + parrsld(i,k,lat)
        pasum  (i) = pasum  (i) + parrsld(i,k,lat)
     end do
  end do
!
! Compute first part of (1/ps)etadot(dp/deta)
!
  do k = 1,plev-1
     do i = 1,nlon
        pa(i) = pa(i) + parrsld(i,k,lat)
        parrsld(i,k,lat) = pa(i)
     end do
!
     if(k.ge.nprlev) then
        do i = 1,nlon
           parrsld(i,k,lat) = parrsld(i,k,lat) - hybi(k+1)*( pasum(i) )
        end do
     end if
  end do
!
! Compute U, V interpolants:  Non-monotonic, 3-D interpolation
!
  limdrh  = .false.
  limdrv  = .true.
  lhrzint = .true.
  lvrtint = .true.
  call sltint (jcen    ,adv_state%u3(:,:,beglatex:endlatex),lam , &
               rdphi   ,rdz     ,lbasdy  ,lbasdz  ,xl      , &
               xr      ,wgt1x   ,wgt2x   ,wgt3x   ,wgt4x   , &
               hl      ,hr      ,dhl     ,dhr     ,ys      , &
               yn      ,wgt1y   ,wgt2y   ,wgt3y   ,wgt4y   , &
               hs      ,hn      ,dhs     ,dhn     ,wgt1z   , &
               wgt2z   ,wgt3z   ,wgt4z   ,hb      ,ht      , &
               dhb     ,dht     ,idp     ,jdp     ,kdp     , &
               kkdp    ,lhrzint ,lvrtint ,limdrh  ,limdrv  , &
               fdp     ,nlon  )
  call sltint (jcen    ,adv_state%v3(:,:,beglatex:endlatex),lam     , &
               rdphi   ,rdz     ,lbasdy  ,lbasdz  ,xl      , &
               xr      ,wgt1x   ,wgt2x   ,wgt3x   ,wgt4x   , &
               hl      ,hr      ,dhl     ,dhr     ,ys      , &
               yn      ,wgt1y   ,wgt2y   ,wgt3y   ,wgt4y   , &
               hs      ,hn      ,dhs     ,dhn     ,wgt1z   , &
               wgt2z   ,wgt3z   ,wgt4z   ,hb      ,ht      , &
               dhb     ,dht     ,idp     ,jdp     ,kdp     , &
               kkdp    ,lhrzint ,lvrtint ,limdrh  ,limdrv  , &
               fdp(1,1,2),nlon  )
!
! Evaluate last half of grfu and grfv (Nu,Nv)
!
  call nunv1(lam(i1,jcen) ,phi(jcen),lamdp        ,phidp              ,fdp(1,1,1), &
             fdp(1,1,2)   ,coslat   ,grfu(1,1,lat),grfv(1,1,lat)      ,grfulat   , &
             grfvlat      ,nlon)
!
! Compute T interpolants:  Non-monotonic, 3-D interpolation
!
  limdrh  = .false.
  limdrv  = .true.
  lhrzint = .true.
  lvrtint = .true.
  call sltint (jcen    ,adv_state%t3(:,:,beglatex:endlatex),lam     , &
               rdphi   ,rdz     ,lbasdy  ,lbasdz  ,xl      , &
               xr      ,wgt1x   ,wgt2x   ,wgt3x   ,wgt4x   , &
               hl      ,hr      ,dhl     ,dhr     ,ys      , &
               yn      ,wgt1y   ,wgt2y   ,wgt3y   ,wgt4y   , &
               hs      ,hn      ,dhs     ,dhn     ,wgt1z   , &
               wgt2z   ,wgt3z   ,wgt4z   ,hb      ,ht      , &
               dhb     ,dht     ,idp     ,jdp     ,kdp     , &
               kkdp    ,lhrzint ,lvrtint ,limdrh  ,limdrv  , &
               fdp(1,1,1),nlon  )
!
! Accumulate T interpolants in 3-D array
!
  do k = 1,plev
     do i = 1,nlon
        grtlat(i,k) = grtlat(i,k) + fdp(i,k,1)
     end do
  end do
!
! Reset "kdp" to arrival indices everywhere so that we can do true
! horizontal interpolation rather than "vertical non-interpolation"
!
  do k = 1,plev
     do i = 1,nlon
        kdp(i,k) = k
     end do
  end do
!
! Compute Ps and remaining T interpolants:
! Non-monotonic, 2-D interpolation
!
  limdrh  = .false.
  limdrv  = .false.
  lhrzint = .true.
  lvrtint = .false.
  call sltint (jcen    ,adv_state%lnpssld(:,:,beglatex:endlatex),lam     ,rdphi   , &
               rdz     ,lbasdy  ,lbasdz  ,xl      ,xr      , &
               wgt1x   ,wgt2x   ,wgt3x   ,wgt4x   ,hl      , &
               hr      ,dhl     ,dhr     ,ys      ,yn      , &
               wgt1y   ,wgt2y   ,wgt3y   ,wgt4y   ,hs      , &
               hn      ,dhs     ,dhn     ,wgt1z   ,wgt2z   , &
               wgt3z   ,wgt4z   ,hb      ,ht      ,dhb     , &
               dht     ,idp     ,jdp     ,kdp     ,kkdp    , &
               lhrzint ,lvrtint ,limdrh  ,limdrv  ,fdp(1,1,1), &
               nlon  )
  call sltint (jcen    ,adv_state%prhssld(:,:,beglatex:endlatex),lam     ,rdphi   , &
               rdz     ,lbasdy  ,lbasdz  ,xl      ,xr      , &
               wgt1x   ,wgt2x   ,wgt3x   ,wgt4x   ,hl      , &
               hr      ,dhl     ,dhr     ,ys      ,yn      , &
               wgt1y   ,wgt2y   ,wgt3y   ,wgt4y   ,hs      , &
               hn      ,dhs     ,dhn     ,wgt1z   ,wgt2z   , &
               wgt3z   ,wgt4z   ,hb      ,ht      ,dhb     , &
               dht     ,idp     ,jdp     ,kdp     ,kkdp    , &
               lhrzint ,lvrtint ,limdrh  ,limdrv  ,fdp(1,1,2), &
               nlon )
  do i = 1,nlon
     pdsum(i) = 0._r8
     pd   (i) = 0._r8
     pdsm1(i) = 0._r8
     pd1  (i) = 0._r8
  end do
!
! Accumulate P-interpolants into T equation
!
  do k = nprlev,plev
     tmp1 = cappa*t0(k)*hypi(plevp)/hypm(k)
     do i = 1,nlon
        grtlat(i,k) = grtlat(i,k) - fdp(i,k,1)*tmp1*hybm(k) 
     end do
  end do
  do k = nprlev,plev
     do l = nprlev,k
        do i = 1,nlon
           grtlat(i,k) = grtlat(i,k) + fdp(i,l,1)*hybd(l)*gamma(l,k)
        end do
     end do
  end do
  do k = 1,plev
     do l = 1,k
        do i = 1,nlon
           grtlat(i,k) = grtlat(i,k) + fdp(i,l,2)*gamma(l,k)
        end do
     end do
  end do
!
! Accumulate Ps interpolants in 3-D array
!
  do k = 1,plev
     do i = 1,nlon
        grpslat(i) = grpslat(i) + fdp(i,k,2)
        pdsum(i) = pdsum(i) + fdp(i,k,2)
     end do
  end do
  do k = nprlev,plev
     do i = 1,nlon
        grpslat(i) = grpslat(i) + fdp(i,k,1)*hybd(k)
        pdsm1(i) = pdsm1(i) + fdp(i,k,1)*hybd(k)
     end do
  end do
!
! Compute remainder of (1/ps)etadot(dp/deta)
!
  do k = 1,plev-1
     do i = 1,nlon
        pd (i) = pd (i) + fdp(i,k,2)
        parrsld(i,k,lat) = parrsld(i,k,lat) + pd (i)
     end do
!
     if(k.ge.nprlev) then
        do i=1,nlon
           pd1(i) = pd1(i) + fdp(i,k,1)*hybd(k)
           parrsld(i,k,lat) = parrsld(i,k,lat) + pd1(i) - &
                hybi(k+1)*(pdsum(i) + pdsm1(i))
        end do
     end if
     do i = 1,nlon
        parrsld(i,k,lat) = parrsld(i,k,lat)*dtr
     end do
  end do
!
! Copy fu,fv,T,q,Ps into FFT buffer
!
  do k = 1,plev
     do i = 1,nlon
        fftbuf(i,k,fudex) = grfulat(i,k)
        fftbuf(i,k,fvdex) = grfvlat(i,k)
        fftbuf(i,k,tdex)  = grtlat(i,k)
        fftbuf(i,k,qdex)  = grqlat(i,k)
     enddo
  enddo
  do i = 1,nlon
     fftbuf(i,1,psdex)    = grpslat(i)
  enddo

  return
end subroutine scanslt_bft

!
!-----------------------------------------------------------------------
!

subroutine scanslt_fft (nlon_fft,nlon_fft2,fftbuf,fftbuf2)
!-----------------------------------------------------------------------
!
! Purpose:
! Compute FFT of non-linear dynamical terms
! in preparation for Fourier -> spectral quadrature.
!
! Author:  J. Olson
! Modified: P. Worley
!
!-----------------------------------------------------------------------

  use pmgrid,      only: plon, plat
  use rgrid,       only: nlon
#if (defined SPMD)
   use mpishorthand, only: mpicom
   use comspe
#endif
   use sld_control_mod, only  : pcray, trig, ifax
!     
! Input arguments
!     
   integer, intent(in) :: nlon_fft         ! first dimension of first FFT work array
   integer, intent(in) :: nlon_fft2        ! first dimension of second FFT work array
!     
! Input/Output arguments
!     
   real(r8), intent(inout) :: fftbuf(nlon_fft,plev,5,beglat:endlat) 
                            ! buffer used for in-place FFTs
!     
! Output arguments
!     
#if (defined SPMD)
   real(r8), intent(out) :: fftbuf2(nlon_fft2,plev,5,plat) 
                            ! buffer for returning reorderd Fourier coefficients
#else
   real(r8), intent(in) :: fftbuf2(1) 
                            ! buffer unused
#endif
!     
!---------------------------Local workspace-----------------------------
!     
! The "work" array has a different size requirement depending upon whether
! the proprietary Cray assembly language version of the FFT library
! routines, or the all-Fortran version, is being used.
!     
#if ( ! defined USEFFTLIB )
   real(r8) work((plon+1)*5*plev)
#else 
   real(r8) work((plon+1)*pcray) ! workspace array for fft991
#endif
   integer lat               ! latitude index
   integer inc               ! increment for fft991
   integer isign             ! flag indicates transform direction
   integer ntr               ! number of transforms to perform
!
   inc = 1
   isign = -1
   ntr = 4*plev + 1
!$OMP PARALLEL DO PRIVATE (LAT,WORK)
   do lat=beglat,endlat
      fftbuf(nlon(lat)+1:nlon_fft,:,:,lat) = 0.0_r8
      call fft991(fftbuf(1,1,1,lat)     ,work    ,trig(1,lat),ifax(1,lat),inc     ,&
                  nlon_fft ,nlon(lat)   ,ntr     ,isign   )
   enddo
!
#if ( defined SPMD )
!
!  reorder Fourier coefficients
!
   call t_barrierf ('sync_realloc4a', mpicom)
   call t_startf('realloc4a')
   call realloc4a(nlon_fft, nlon_fft2, fftbuf, fftbuf2)
   call t_stopf('realloc4a')
#endif

   return
end subroutine scanslt_fft

!
!-----------------------------------------------------------------------
!

subroutine scanslt_aft (irow    ,nlon_fft,fftbufs ,fftbufn , &
                        grlps1  ,grlps2  ,grt1    ,grt2    , &
                        grq1    ,grq2    ,grfu1   ,grfu2   , &
                        grfv1   ,grfv2   )
!-----------------------------------------------------------------------
!
! Purpose:
! Combine terms in preparation for Fourier -> spectral quadrature.
!
! Author:  J. Olson
! Modified: P. Worley
!
!-----------------------------------------------------------------------

  use spmd_utils,      only: iam
  use pspect,      only: pmmax
#if (defined SPMD)
   use comspe, only: numm, maxm
#else
   use comspe, only: maxm
   use rgrid, only: nmmax
#endif


!
! Input arguments
!     
  integer, intent(in)  :: irow                ! latitude pair index
  integer, intent(in)  :: nlon_fft            ! first dimension of FFT work arrays

  real(r8), intent(in) :: fftbufs(nlon_fft,plev,5) ! southern latitude Fourier coefficients
  real(r8), intent(in) :: fftbufn(nlon_fft,plev,5) ! northern latitude Fourier coefficients

  real(r8), intent(out)  :: grlps1(2*maxm) ! ------------------------------
  real(r8), intent(out)  :: grlps2(2*maxm) ! |
  real(r8), intent(out)  :: grt1  (2*maxm,plev) ! |
  real(r8), intent(out)  :: grt2  (2*maxm,plev) ! |
  real(r8), intent(out)  :: grq1  (2*maxm,plev) ! |- see quad for definitions
  real(r8), intent(out)  :: grq2  (2*maxm,plev) ! |
  real(r8), intent(out)  :: grfu1 (2*maxm,plev) ! |
  real(r8), intent(out)  :: grfu2 (2*maxm,plev) ! |
  real(r8), intent(out)  :: grfv1 (2*maxm,plev) ! |
  real(r8), intent(out)  :: grfv2 (2*maxm,plev) ! ------------------------------
!
!---------------------------Local workspace-----------------------------
!
  integer i                    ! index
  integer k                    ! index
  integer mlength, mstrt       ! number of wavenumbers and index offset
!
  integer, parameter :: fudex = 1  ! indices into fftbuf 
  integer, parameter :: fvdex = 2
  integer, parameter :: tdex  = 3
  integer, parameter :: qdex  = 4
  integer, parameter :: psdex = 5
!
!-----------------------------------------------------------------------
!
#if (defined SPMD)
   mlength = numm(iam)
#else
   mlength = nmmax(irow)
#endif
  do k = 1,plev
     do i = 1,2*mlength
        grfu1(i,k) = 0.5_r8*(fftbufn(i,k,fudex)+fftbufs(i,k,fudex))
        grfu2(i,k) = 0.5_r8*(fftbufn(i,k,fudex)-fftbufs(i,k,fudex))

        grfv1(i,k) = 0.5_r8*(fftbufn(i,k,fvdex)+fftbufs(i,k,fvdex))
        grfv2(i,k) = 0.5_r8*(fftbufn(i,k,fvdex)-fftbufs(i,k,fvdex))

        grt1 (i,k) = 0.5_r8*(fftbufn(i,k,tdex)+fftbufs(i,k,tdex))
        grt2 (i,k) = 0.5_r8*(fftbufn(i,k,tdex)-fftbufs(i,k,tdex))

        grq1 (i,k) = 0.5_r8*(fftbufn(i,k,qdex)+fftbufs(i,k,qdex))
        grq2 (i,k) = 0.5_r8*(fftbufn(i,k,qdex)-fftbufs(i,k,qdex))
     end do
  end do

  do i = 1,2*mlength
     grlps1(i) = 0.5_r8*(fftbufn(i,1,psdex)+fftbufs(i,1,psdex))
     grlps2(i) = 0.5_r8*(fftbufn(i,1,psdex)-fftbufs(i,1,psdex))
  end do
!
  return
end subroutine scanslt_aft

!
!-----------------------------------------------------------------------
!

subroutine slt_final( adv_state )
!-----------------------------------------------------------------------
!
! Purpose:
! slt finalization
!
!-----------------------------------------------------------------------
  use rgrid,        only: nlon
  use spmd_utils,   only: masterproc, npes
#if ( defined SPMD )
  use mpishorthand, only: mpicom, mpir8
  use spmd_dyn,     only: compute_gsfactors
#endif

  type(advection_state), intent(inout) :: adv_state
!
! Local variables
!
  integer :: i, lat, k             ! |
  integer :: kk                    ! |
  integer :: kk1                   ! |
  integer :: kk2                   ! | - indices
  integer :: l                     ! |
  integer :: lvsum                 ! counter
  integer :: lsum                  ! counter
  real(r8) :: vsum                 ! accumulator for SLD binning statistics
  logical :: vflag                 ! logical flag to indicate that binning has been done
!                                  ! correctly
#ifdef SPMD
  integer :: numsend               ! number of items to be sent
  integer :: numrecv(0:npes-1)     ! number of items to be received
  integer :: displs(0:npes-1)      ! displacement array
  integer :: numperlat             ! number of items per latitude band
#endif
!-----------------------------------------------------------------------

  call dealloc_advstate( adv_state )
!
! Compute some statistical info
!
#if ( defined SPMD )
  numperlat = plev*plev
  call compute_gsfactors (numperlat, numsend, numrecv, displs)
  call mpiallgatherv ( levkntl(:,:,beglat:endlat), numsend, mpir8, &
                       levkntl(:,:,:), numrecv, displs, mpir8, mpicom)
#endif
  if (masterproc) then
     do k = 1,plev
        do kk = 1,plev
            levknt(kk,k) = 0._r8
            do lat = 1,plat
               levknt(kk,k) = levknt(kk,k) + levkntl(kk,k,lat)
            end do
         end do
      end do
      lsum = 0
      do lat = 1,plat
         lsum = lsum + nlon(lat)
      end do

      vflag = .false.
      do k = 1,plev
         vsum = 0._r8
         do kk = 1,plev
            vsum = vsum + levknt(kk,k)
         end do
         lvsum = vsum + 0.01_r8
         if( lvsum .ne. lsum*kstep ) vflag = .true.
      end do

      do k = 1,plev
         do kk = 1,plev
            levknt(kk,k) = levknt(kk,k)/real(lsum*kstep,r8)
            if (levknt(kk,k) .eq. 0._r8) levknt(kk,k) = 1.e+36_r8
         end do
      end do

      if(vflag) then
         write(iulog,*) '********************************************'
         write(iulog,1000)
         write(iulog,*) '********************************************'
      else
         write(iulog,2000)
         k = 1
         write(iulog,3001) hypm(k)/100._r8,(levknt(kk,k),kk = 1, 6)
         k = 2
         write(iulog,3002) hypm(k)/100._r8,(levknt(kk,k),kk = 1, 7)
         k = 3
         write(iulog,3003) hypm(k)/100._r8,(levknt(kk,k),kk = 1, 8)
         k = 4
         write(iulog,3004) hypm(k)/100._r8,(levknt(kk,k),kk = 1, 9)
         k = 5
         write(iulog,3005) hypm(k)/100._r8,(levknt(kk,k),kk = 1,10)
         do k = 6,plev
            kk1 = k-5
            kk2 = k+5
            if(kk2 .gt. plev) kk2 = plev
            write(iulog,3000) hypm(k)/100._r8,(levknt(kk,k),kk = kk1,kk2)
         end do
      end if
  end if

1000 format(' ERROR in binning departure points; sums were done'// &
             ' incorrectly. not printing results.')
2000  format(/40x,' BINNING STATISTICS FOR VERTICAL DEPARTURE POINTS'// &
                  '    level          -5        -4        -3        -2     ', &
                  '   -1    0 (arr pt)     1         2         3         4 ', &
                  ' 5'/)
3000  format(' ',f9.4, 4x,11f10.5)
3001  format(' ',f9.4,54x,11f10.5)
3002  format(' ',f9.4,44x,11f10.5)
3003  format(' ',f9.4,34x,11f10.5)
3004  format(' ',f9.4,24x,11f10.5)
3005  format(' ',f9.4,14x,11f10.5)

end subroutine slt_final

!
!-----------------------------------------------------------------------
!

subroutine bndexch( adv_state )

!-----------------------------------------------------------------------
!
! Purpose:
! Pack and Exchange initial prognostic information among all the 
! processors
!
!-----------------------------------------------------------------------
#ifdef SPMD
  use spmd_utils,only: iam
  use spmd_dyn,  only: cut, dyn_npes
#endif
!
! Arguments
!
  type(advection_state), intent(inout) :: adv_state
#ifdef SPMD
!
! Local workspace
!
  integer ns, nn
  integer inreg ( 2 )
  integer outreg( 2 )
  integer others,othern   ! Other node
!
! Return if number of processors is less than 2
!
  if (dyn_npes .lt. 2) return
!
! For each partition (south and north) communicate boundaries
! on each side of partition among however many neighbors necessary
!
! send south, receive north
!
  ns = 1
  nn = 1
  do while (ns .le. neighs .or. nn .le. neighn)
     if (ns .le. neighs) then
        others = neighs_proc(ns)
!
! Intersection of my cuts and neighbor processor's extended
! cuts tells if this node needs to send data to neighbor 
!
        call intersct(cut(1,iam),cutex(1,others),outreg  )
        ns = ns + 1
     else
        others = -1
        outreg(1) = 0
        outreg(2) = 0
     end if

     if (nn .le. neighn) then
        othern = neighn_proc(nn)
!
! Intersection of neighbor cuts and this node's extended
! cut tells if this node receives data from neighbor 
!
        call intersct(cut(1,othern),cutex(1,iam),inreg   )
        nn = nn + 1
     else
        othern = -1
        inreg(1) = 0
        inreg(2) = 0
     end if

     call bndexch_mpi(others  ,outreg  ,othern  ,inreg, adv_state   )
  end do

!
! send north, receive south
!
  ns = 1
  nn = 1
  do while (ns .le. neighs .or. nn .le. neighn)
     if (nn .le. neighn) then
        othern = neighn_proc(nn)
!
! Intersection of my cuts and neighbor processor's extended
! cuts tells if this node needs to send data to neighbor 
!
        call intersct(cut(1,iam),cutex(1,othern),outreg  )  
        nn = nn + 1
     else
        othern = -1
        outreg(1) = 0
        outreg(2) = 0
     end if

     if (ns .le. neighs) then
        others = neighs_proc(ns)
!
! Intersection of neighbor cuts and this node's extended
! cut tells if this node receives data from neighbor 
!
        call intersct(cut(1,others),cutex(1,iam),inreg   )  
        ns = ns + 1
     else
        others = -1
        inreg(1) = 0
        inreg(2) = 0
     end if

     call bndexch_mpi(othern  ,outreg  ,others  ,inreg, adv_state   )
  end do
#endif
  return
end subroutine bndexch

!
!-----------------------------------------------------------------------
!

#ifdef SPMD
subroutine bndexch_mpi(othero  ,outreg  ,otheri  ,inreg, adv_state   )
!-----------------------------------------------------------------------
! Send initial prognostic information to my peer process
!-----------------------------------------------------------------------
  use mpishorthand

!---------------------------- Parameters -------------------------------
  integer, parameter :: msgtype = 6000 ! id for message passing
  integer, parameter :: j1m = j1 - 1   ! lat index just before first "real" model lat
!------------------------------Arguments--------------------------------
!
  integer , intent(in)   :: othero
  integer , intent(in)   :: outreg(2)
  integer , intent(in)   :: otheri
  integer , intent(in)   :: inreg(2)
  type(advection_state), intent(inout) :: adv_state
!
!---------------------------Local workspace-----------------------------
!
  integer num
  integer msg

  integer, parameter :: numvars = 8+pcnst
  integer reqs(numvars*(plat+1))
  integer stats(MPI_STATUS_SIZE, numvars*(plat+1))

  integer reqr(numvars*(plat+1))
  integer statr(MPI_STATUS_SIZE, numvars*(plat+1))

  integer i,j,c
  integer reqs_i,reqr_i
!
!-----------------------------------------------------------------------
!
  reqs_i = 0
  if (othero .ne. -1) then
     do i = outreg(1), outreg(2)
        j = numvars*(i-outreg(1))
        msg = msgtype + j
        reqs_i = reqs_i + 1
        call mpiisend (adv_state%u3(1,1,j1m+i),plndlv,mpir8, &
                       othero,msg,mpicom,reqs(reqs_i))

        msg = msgtype + j + 1
        reqs_i = reqs_i + 1
        call mpiisend (adv_state%v3(1,1,j1m+i),plndlv,mpir8, &
                       othero,msg,mpicom,reqs(reqs_i))

        msg = msgtype + j + 2
        reqs_i = reqs_i + 1
        call mpiisend (adv_state%t3(1,1,j1m+i),plndlv,mpir8, &
                       othero,msg,mpicom,reqs(reqs_i))

        num = plndlv
        do c = 1, pcnst
           msg = msgtype + j + 2 + c
           reqs_i = reqs_i + 1
           call mpiisend (adv_state%q3(1,1,j1m+i,c),num,mpir8, &
                          othero,msg,mpicom,reqs(reqs_i))
        end do

        msg = msgtype + j + 3 + pcnst
        reqs_i = reqs_i + 1
        call mpiisend (adv_state%lnpssld(1,1,j1m+i),plndlv,mpir8, &
                       othero,msg,mpicom,reqs(reqs_i))

        msg = msgtype + j + 4 + pcnst
        reqs_i = reqs_i + 1
        call mpiisend (adv_state%prhssld(1,1,j1m+i),plndlv,mpir8, &
                       othero,msg,mpicom,reqs(reqs_i))

        msg = msgtype + j + 5 + pcnst
        reqs_i = reqs_i + 1
        call mpiisend (adv_state%u3sld(1,1,j1m+i),plndlv,mpir8, &
                       othero,msg,mpicom,reqs(reqs_i))

        msg = msgtype + j + 6 + pcnst
        reqs_i = reqs_i + 1
        call mpiisend (adv_state%v3sld(1,1,j1m+i),plndlv,mpir8, &
                       othero,msg,mpicom,reqs(reqs_i))

        msg = msgtype + j + 7 + pcnst
        reqs_i = reqs_i + 1
        call mpiisend (adv_state%etadot(1,1,j1m+i),plond*plevp,mpir8, &
                       othero,msg,mpicom,reqs(reqs_i))

     end do
  end if

  reqr_i = 0
  if (otheri .ne. -1) then
     do i = inreg(1), inreg(2)
        j = numvars*(i-inreg(1))
        msg = msgtype + j
        reqr_i = reqr_i + 1
        call mpiirecv (adv_state%u3(1,1,j1m+i),plndlv,mpir8, &
                       otheri,msg,mpicom,reqr(reqr_i))

        msg = msgtype + j + 1
        reqr_i = reqr_i + 1
        call mpiirecv (adv_state%v3(1,1,j1m+i),plndlv,mpir8, &
                       otheri,msg,mpicom,reqr(reqr_i))

        msg = msgtype + j + 2
        reqr_i = reqr_i + 1
        call mpiirecv (adv_state%t3(1,1,j1m+i),plndlv,mpir8, &
                       otheri,msg,mpicom,reqr(reqr_i))

        num = plndlv
        do c = 1, pcnst
           msg = msgtype + j + 2 + c
           reqr_i = reqr_i + 1
           call mpiirecv (adv_state%q3(1,1,j1m+i,c),num,mpir8, &
                          otheri,msg,mpicom,reqr(reqr_i))
        end do

        msg = msgtype + j + 3 + pcnst
        reqr_i = reqr_i + 1
        call mpiirecv (adv_state%lnpssld(1,1,j1m+i),plndlv,mpir8, &
                       otheri,msg,mpicom,reqr(reqr_i))

        msg = msgtype + j + 4 + pcnst
        reqr_i = reqr_i + 1
        call mpiirecv (adv_state%prhssld(1,1,j1m+i),plndlv,mpir8, &
                       otheri,msg,mpicom,reqr(reqr_i))

        msg = msgtype + j + 5 + pcnst
        reqr_i = reqr_i + 1
        call mpiirecv (adv_state%u3sld(1,1,j1m+i),plndlv,mpir8, &
                       otheri,msg,mpicom,reqr(reqr_i))

        msg = msgtype + j + 6 + pcnst
        reqr_i = reqr_i + 1
        call mpiirecv (adv_state%v3sld(1,1,j1m+i),plndlv,mpir8, &
                       otheri,msg,mpicom,reqr(reqr_i))

        msg = msgtype + j + 7 + pcnst
        reqr_i = reqr_i + 1
        call mpiirecv (adv_state%etadot(1,1,j1m+i),plond*plevp,mpir8, &
                       otheri,msg,mpicom,reqr(reqr_i))
     end do
  end if

  if (reqs_i .ne. 0) then
     call mpiwaitall(reqs_i,reqs,stats)
  end if

  if (reqr_i .ne. 0) then
     call mpiwaitall(reqr_i,reqr,statr)
  end if

  return
end subroutine bndexch_mpi

!
!-----------------------------------------------------------------------
!

subroutine intersct (regiona ,regionb ,regionc )
!
!-----------------------------------------------------------------------
!
! Purpose:
! Given two regions (a,b) output the intersection (common latitudes)  
! of these two sets.  The routine is used in bndexch to determine which
! latitudes need to be communicated to neighboring processors.  Typically
! this routine is invoked as the intersection of the set of resident 
! latitudes on processor A with the set of extended latitudes (needed for 
! the SLT) of processor B.  Any common latitudes will need to be 
! communicated to B to complete SLT processing. 
!
! Original version:  J. Rosinski
!
!------------------------------Arguments--------------------------------
!
  integer , intent(in)   :: regiona( 2 )
  integer , intent(in)   :: regionb( 2 )
  integer , intent(out)  :: regionc( 2 )
!
!-----------------------------------------------------------------------
!
  regionc( 1 ) = max( regiona( 1 ), regionb( 1 ) )
  regionc( 2 ) = min( regiona( 2 ), regionb( 2 ) )

  return
end subroutine intersct
#endif

!
!-----------------------------------------------------------------------
!

subroutine grdini(pmap    ,etamid  ,etaint  ,gravit  ,dlam    , &
                  lam     ,phi     ,dphi    ,gw      ,sinlam  , &
                  coslam  ,lbasdy  ,lbasdz  ,lbassd  ,lbasiy  , &
                  lbasiz  ,lbassi  ,detam   ,detai   ,kdpmpf  , &
                  kdpmph  ,cwava   ,phigs   )
!-----------------------------------------------------------------------
!
! Purpose:
! Initialize model and extended grid parameters
! Initialize weights for Lagrange cubic derivative estimates
! Initialize weights for Lagrange cubic interpolant
!
! Author:  J. Olson
!
!-----------------------------------------------------------------------
  use rgrid,  only: nlon
!------------------------------Arguments--------------------------------
!
  integer , intent(in)   :: pmap                 ! dimension of artificial vert. grid
  real(r8), intent(in)   :: etamid(plev)         ! full-level model vertical grid
  real(r8), intent(in)   :: etaint(plevp)        ! half-level model vertical grid
  real(r8), intent(in)   :: gravit               ! gravitational constant
  real(r8), intent(out)  :: dlam  (platd)        ! longitudinal grid interval (radians)
  real(r8), intent(out)  :: lam   (plond,platd)  ! longitudinal coords of extended grid
  real(r8), intent(out)  :: phi   (platd)        ! latitudinal  coords of extended grid
  real(r8), intent(out)  :: dphi  (platd)        ! latitude intervals (radians)
  real(r8), intent(out)  :: gw    (plat)         ! Gaussian weights
  real(r8), intent(out)  :: sinlam(plond,platd)  ! sin(lam) model domain only
  real(r8), intent(out)  :: coslam(plond,platd)  ! cos(lam) model domain only
  real(r8), intent(out)  :: lbasdy(4,2,platd)    ! latitude derivative weights
  real(r8), intent(out)  :: lbasdz(4,2,plev)     ! vertical (full levels) deriv weights
  real(r8), intent(out)  :: lbassd(4,2,plevp)    ! vertical (half levels) deriv weights
  real(r8), intent(out)  :: lbasiy(4,2,platd)    ! Lagrange cubic interp weights (lat.)
  real(r8), intent(out)  :: lbasiz(4,2,plev)     ! Lagrange cubic interp wghts (full lev)
  real(r8), intent(out)  :: lbassi(4,2,plevp)    ! Lagrange cubic interp wghts (half lev)
  real(r8), intent(out)  :: detam (plev)         ! intervals between vertical full levs.
  real(r8), intent(out)  :: detai (plevp)        ! intervals between vertical half levs.
  integer , intent(out)  :: kdpmpf(pmap)         ! artificial full vertical grid indices
  integer , intent(out)  :: kdpmph(pmap)         ! artificial half vertical grid indices
  real(r8), intent(out)  :: cwava(plat)          ! weight applied to global integrals
  real(r8), intent(out)  :: phigs                ! Cutoff latitude for local geodesic
!
!---------------------------Local variables-----------------------------
!
  integer j                 ! index
  integer k                 ! index
  real(r8) etamln (plev)    ! log(etamid)
  real(r8) etailn (plevp)   ! log(etaint)
  real(r8) detamln(plev)    ! dlog(etamid)
  real(r8) detailn(plevp)   ! dlog(etaint)
!
!-----------------------------------------------------------------------
!
! Initialize extended horizontal grid coordinates.
!
  call grdxy(dlam    ,lam     ,phi     ,gw      ,sinlam  , &
             coslam  )
!
! Basis functions for computing Lagrangian cubic derivatives
! on unequally spaced latitude and vertical grids.
!
  call basdy(phi     ,lbasdy  )
  call basdz(plev    ,etamid  ,lbasdz  )
  call basdz(plevp   ,etaint  ,lbassd  )
!
! Basis functions for computing weights for Lagrangian cubic
! interpolation on unequally spaced latitude and height grids.
!
  call basiy(phi     ,lbasiy  )
  call basiz(plev    ,etamid  ,lbasiz  )
  call basiz(plevp   ,etaint  ,lbassi  )
!
! Compute interval lengths in latitudinal grid
!
  do j = 1,platd-1
     dphi(j) = phi(j+1) - phi(j)
  end do
!
! Compute interval lengths in vertical grids.
!
  do k = 1,plev
     etamln(k) = log(etamid(k))
  end do
  do k = 1,plevp
     etailn(k) = log(etaint(k))
  end do
  do k = 1,plev-1
     detam  (k) = etamid(k+1) - etamid(k)
     detamln(k) = etamln(k+1) - etamln(k)
  end do
  do k = 1,plev
     detai  (k) = etaint(k+1) - etaint(k)
     detailn(k) = etailn(k+1) - etailn(k)
  end do
!
! Build artificial evenly spaced vertical grid for use in determining
! vertical position of departure point.
! Build one grid for full model levels and one for half levels.
!
  call vrtmap(plev    ,pmap    ,etamln  ,detamln ,kdpmpf  )
  call vrtmap(plevp   ,pmap    ,etailn  ,detailn ,kdpmph  )
!
! Compute tracer integral constant
!
  do j=1,plat
     cwava(j) = 1._r8/(nlon(j)*gravit)
  end do
!
! Initialize cutoff latitude (poleward of which local geodesic
! coordinates will be used); radians
!
  phigs = 1.221730_r8
!
  return
end subroutine grdini

!
!-----------------------------------------------------------------------
!

subroutine grdxy(dlam    ,lam     ,phi     ,w       ,sinlam  , &
                 coslam  )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Define the "extended" grid used in the semi-Lagrangian transport
! scheme.  The longitudes are equally spaced and the latitudes are
! Gaussian.  The global grid is extended to include "wraparound" points
! on all sides.
! 
! Method: 
! 
! Author: J. Olson
! 
!-----------------------------------------------------------------------
  use gauaw_mod, only: gauaw
!------------------------------Parameters-------------------------------
  integer, parameter :: istart = nxpt+1         ! index for first model long.
  integer, parameter :: jstart = nxpt+jintmx+1  ! index for first model lat.
  integer, parameter :: jstop  = jstart-1+plat  ! index for last  model lat.
!-----------------------------------------------------------------------

!------------------------------Arguments--------------------------------
  real(r8), intent(out) :: dlam(platd)          ! longitudinal increment
  real(r8), intent(out) :: lam   (plond,platd)  ! long. coords. in extended grid
  real(r8), intent(out) :: phi   (platd)        ! lat.  coords. in extended grid
  real(r8), intent(out) :: w     (plat)         ! Gaussian weights
  real(r8), intent(out) :: sinlam(plond,platd)  ! sin(lam)
  real(r8), intent(out) :: coslam(plond,platd)  ! cos(lam)
!
! dlam    Length of increment in longitude grid.
! lam     Longitude values in the extended grid.
! phi     Latitude values in the extended grid.
! w       Gauss weights for latitudes in the global grid.  (These sum
!         to 2.0 like the ones in CCM1.)
! sinlam  Sine of longitudes in global grid (no extension points).
! coslam  Cosine of longitudes in global grid (no extension points).
!-----------------------------------------------------------------------

!---------------------------Local variables-----------------------------
  integer i,j,ig            ! indices
  integer nlond             ! extended long dim
  real(r8) lam0             ! lamda = 0
  real(r8) pi               ! 3.14...
  real(r8) wrk(platd)       ! work space
!-----------------------------------------------------------------------
!
  lam0 = 0.0_r8
  pi = 4._r8*atan(1._r8)
!
! Interval length in equally spaced longitude grid.
!
  do j=1,platd
     dlam(j) = 2._r8*pi/real(nlonex(j),r8)
!
! Longitude values on extended grid.
!
     nlond = nlonex(j) + 1 + 2*nxpt
     do i = 1,nlond
        lam(i,j) = real(i-istart,r8)*dlam(j) + lam0
     end do
  end do
!
! Compute Gauss latitudes and weights.  On return; phi contains the
! sine of the latitudes starting closest to the north pole and going
! toward the south; w contains the corresponding Gauss weights.
!
  call gauaw(phi     ,w       ,plat    )
!
! Reorder and compute latitude values.
!
  do j = jstart,jstop
     wrk(j) = asin( phi(jstop-j+1) )
  end do
  phi(jstart:jstop) = wrk(jstart:jstop)
!
! North and south poles.
!
  phi(jstart-1) = -pi/2.0_r8
  phi(jstop +1) =  pi/2.0_r8
!
! Extend Gauss latitudes below south pole so that the spacing above
! the pole is symmetric, and phi is decreasing, i.e., phi < -pi/2
!
  if( jstart > 2 )then
     do j = 1,jstart-2
        phi(j) = -pi - phi(2*jstart-2-j)
     end do
  end if
!
! Analogously for Northern Hemisphere
!
  if( platd > jstop+1 )then
     do j = jstop+2,platd
        phi(j) = pi - phi(2*jstop+2-j)
     end do
  end if
!
! Sine and cosine of longitude.
!
  do j=1,platd
     ig = 0
     do i = istart,nlonex(j)+nxpt
        ig = ig + 1
        sinlam(ig,j) = sin( lam(i,j) )
        coslam(ig,j) = cos( lam(i,j) )
     end do
  end do

  return
end subroutine grdxy

!
!-----------------------------------------------------------------------
!

subroutine sltlinint(kdim    ,ikcnst  ,fb      ,xl      ,xr      , &
                     ys      ,yn      ,zt      ,zb      ,idp     , &
                     jdp     ,kdp     ,fdp     ,nlon    )
!
!-----------------------------------------------------------------------
!
! Purpose:
! Interpolate field to departure points using tri-linear interpolation
!
! Author:  J. Olson
!
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
!
  integer , intent(in)   :: kdim             ! vertical dimension
  integer , intent(in)   :: ikcnst           ! constituent index

  real(r8), intent(in)   :: fb (plond,kdim*ikcnst,beglatex:endlatex) ! input field
  real(r8), intent(in)   :: xl (plon,plev,4) ! weight for x-interpolants (left)
  real(r8), intent(in)   :: xr (plon,plev,4) ! weight for x-interpolants (right)
  real(r8), intent(in)   :: ys (plon,plev)   ! weight for y-interpolants (south)
  real(r8), intent(in)   :: yn (plon,plev)   ! weight for y-interpolants (north)
  real(r8), intent(in)   :: zt (plon,plev)   ! top vertical interpolation weight 
  real(r8), intent(in)   :: zb (plon,plev)   ! bot vertical interpolation weight 
  integer , intent(in)   :: idp(plon,plev,4) ! index of x-coordinate of dep pt
  integer , intent(in)   :: jdp(plon,plev)   ! index of y-coordinate of dep pt
  integer , intent(in)   :: kdp(plon,plev)   ! index of z-coordinate of dep pt
  real(r8), intent(out)  :: fdp(plon,plev)   ! interpolant
  integer , intent(in)   :: nlon             ! number of longitudes for this latitude
!
!---------------------------Local workspace-----------------------------
!
  integer i, k, ii2, ii3, jj, kk             ! indices
!
!-----------------------------------------------------------------------
!
  do k=1,plev
     do i = 1,nlon
        ii2 = idp(i,k,2)
        ii3 = idp(i,k,3)
        jj  = jdp(i,k)
        kk  = kdp(i,k)
        fdp(i,k) = (( fb (ii2  ,kk  ,jj  )*xl (i,k,2) &
                 +    fb (ii2+1,kk  ,jj  )*xr (i,k,2) )*ys(i,k) &
                 +  ( fb (ii3  ,kk  ,jj+1)*xl (i,k,3) &
                 +    fb (ii3+1,kk  ,jj+1)*xr (i,k,3) )*yn(i,k) )*zt(i,k) &
                 + (( fb (ii2  ,kk+1,jj  )*xl (i,k,2) &
                 +    fb (ii2+1,kk+1,jj  )*xr (i,k,2) )*ys(i,k) &
                 +  ( fb (ii3  ,kk+1,jj+1)*xl (i,k,3) &
                 +    fb (ii3+1,kk+1,jj+1)*xr (i,k,3) )*yn(i,k) )*zb(i,k)
     end do
  end do

  return
end subroutine sltlinint

!
!-----------------------------------------------------------------------
!

subroutine slttraj(pmap    ,jcen    ,lat     ,ztodt   ,ra      , &
                   iter    ,lam     ,phi     ,dphi    ,etamid  , &
                   etaint  ,detam   ,detai   ,lbasiy  ,lbasiz  , &
                   lbassi  ,kdpmpf  ,kdpmph  ,idp     ,jdp     , &
                   kdp     ,kkdp    ,xl      ,xr      ,wgt1x   , &
                   wgt2x   ,wgt3x   ,wgt4x   ,hl      ,hr      , &
                   dhl     ,dhr     ,ys      ,yn      ,wgt1y   , &
                   wgt2y   ,wgt3y   ,wgt4y   ,hs      ,hn      , &
                   dhs     ,dhn     ,rdphi   ,wgt1z   ,wgt2z   , &
                   wgt3z   ,wgt4z   ,hb      ,ht      ,dhb     , &
                   dht     ,rdz     ,lampr   ,phipr   ,upr     , &
                   vpr     ,lamdp   ,phidp   ,sigdp   ,u3      , &
                   v3      ,u3sld   ,v3sld   ,etadot  ,u3p     , &
                   v3p     ,etadotp ,dlam    ,nlon    )
!
!-----------------------------------------------------------------------
!
! Purpose:
! Determine trajectory departure point for each arrival point.
! Assume that the first guess for the dep point coordinates are the
! arrival point coordinates
!
! Author:  J. Olson
!
!-----------------------------------------------------------------------
#if (!defined UNICOSMP)
  use srchutil
#endif
!------------------------------Arguments--------------------------------
!
  integer , intent(in)   :: pmap                ! artificial vert grid dim.
  integer , intent(in)   :: jcen                ! index of lat slice(extend)
  integer , intent(in)   :: lat                 ! index of lat slice (model)
  real(r8), intent(in)   :: ztodt               ! time step (seconds)
  real(r8), intent(in)   :: ra                  ! 1./(radius of earth)
  integer , intent(in)   :: iter                ! iteration count
  real(r8), intent(in)   :: lam   (plond,platd) ! long. coord of model grid
  real(r8), intent(in)   :: phi   (platd)       ! lat.  coord of model grid
  real(r8), intent(in)   :: dphi  (platd)       ! increment between lats.
  real(r8), intent(in)   :: etamid(plev)        ! vertical full levels
  real(r8), intent(in)   :: etaint(plevp)       ! vertical half levels
  real(r8), intent(in)   :: detam (plev)        ! increment between full levs
  real(r8), intent(in)   :: detai (plevp)       ! increment between half levs
  real(r8), intent(in)   :: lbasiy(4,2,platd)   ! lat interp wts(lagrng)
  real(r8), intent(in)   :: lbasiz(4,2,plev)    ! vert interp wghts (lagrng)
  real(r8), intent(in)   :: lbassi(4,2,plevp)   ! vert interp wghts (lagrng)

  integer , intent(in)   :: kdpmpf(pmap)        ! artificial vert grid indices
  integer , intent(in)   :: kdpmph(pmap)        ! artificial vert grid indices
  integer , intent(out)  :: idp   (plon,plev,4) ! zonal      dep point index
  integer , intent(out)  :: jdp   (plon,plev)   ! meridional dep point index
  integer , intent(out)  :: kdp   (plon,plev)   ! vertical   dep point index
  integer , intent(in)   :: kkdp  (plon,plev)   ! index of z-coordinate of dep pt (alt)

  real(r8), intent(out)  :: xl    (plon,plev,4) ! weight for x-interpolants (left)
  real(r8), intent(out)  :: xr    (plon,plev,4) ! weight for x-interpolants (right)
  real(r8), intent(in)   :: wgt1x (plon,plev,4) ! weight for x-interpolants (Lag Cubic)
  real(r8), intent(in)   :: wgt2x (plon,plev,4) ! weight for x-interpolants (Lag Cubic)
  real(r8), intent(in)   :: wgt3x (plon,plev,4) ! weight for x-interpolants (Lag Cubic)
  real(r8), intent(in)   :: wgt4x (plon,plev,4) ! weight for x-interpolants (Lag Cubic)
  real(r8), intent(in)   :: hl    (plon,plev,4) ! weight for x-interpolants (Hermite)
  real(r8), intent(in)   :: hr    (plon,plev,4) ! weight for x-interpolants (Hermite)
  real(r8), intent(in)   :: dhl   (plon,plev,4) ! weight for x-interpolants (Hermite)
  real(r8), intent(in)   :: dhr   (plon,plev,4) ! weight for x-interpolants (Hermite)

  real(r8), intent(out)  :: ys    (plon,plev)   ! weight for y-interpolants (south)
  real(r8), intent(out)  :: yn    (plon,plev)   ! weight for y-interpolants (north)
  real(r8), intent(in)   :: wgt1y (plon,plev)   ! weight for y-interpolants (Lag Cubic)
  real(r8), intent(in)   :: wgt2y (plon,plev)   ! weight for y-interpolants (Lag Cubic)
  real(r8), intent(in)   :: wgt3y (plon,plev)   ! weight for y-interpolants (Lag Cubic)
  real(r8), intent(in)   :: wgt4y (plon,plev)   ! weight for y-interpolants (Lag Cubic)
  real(r8), intent(in)   :: hs    (plon,plev)   ! weight for y-interpolants (Hermite)
  real(r8), intent(in)   :: hn    (plon,plev)   ! weight for y-interpolants (Hermite)
  real(r8), intent(in)   :: dhs   (plon,plev)   ! weight for y-interpolants (Hermite)
  real(r8), intent(in)   :: dhn   (plon,plev)   ! weight for y-interpolants (Hermite)
  real(r8), intent(in)   :: rdphi (plon,plev)   ! reciprocal of y-interval

  real(r8), intent(in)   :: wgt1z (plon,plev)   ! weight for z-interpolants (Lag Cubic)
  real(r8), intent(in)   :: wgt2z (plon,plev)   ! weight for z-interpolants (Lag Cubic)
  real(r8), intent(in)   :: wgt3z (plon,plev)   ! weight for z-interpolants (Lag Cubic)
  real(r8), intent(in)   :: wgt4z (plon,plev)   ! weight for z-interpolants (Lag Cubic)
  real(r8), intent(in)   :: hb    (plon,plev)   ! weight for z-interpolants (Hermite)
  real(r8), intent(in)   :: ht    (plon,plev)   ! weight for z-interpolants (Hermite)
  real(r8), intent(in)   :: dhb   (plon,plev)   ! weight for z-interpolants (Hermite)
  real(r8), intent(in)   :: dht   (plon,plev)   ! weight for z-interpolants (Hermite)
  real(r8), intent(in)   :: rdz   (plon,plev)   ! reciprocal of z-interval

  real(r8), intent(out)  :: lampr (plon,plev)   ! trajectory increment (x-direction)
  real(r8), intent(out)  :: phipr (plon,plev)   ! trajectory increment (y-direction)
  real(r8), intent(out)  :: upr   (plon,plev)   ! interpolated u field (local geodesic)
  real(r8), intent(out)  :: vpr   (plon,plev)   ! interpolated v field (local geodesic)
  real(r8), intent(out)  :: lamdp (plon,plev)   ! zonal      departure pt. coord.
  real(r8), intent(inout):: phidp (plon,plev)   ! meridional departure pt. coord.
  real(r8), intent(inout):: sigdp (plon,plev)   ! vertical   departure pt. coord.
  real(r8), intent(in)   :: u3    (plond, plev, beglatex:endlatex)   ! u wind component
  real(r8), intent(in)   :: v3    (plond, plev, beglatex:endlatex)   ! v wind component
  real(r8), intent(in)   :: u3sld (plond,plev ,beglatex:endlatex)    ! u3 inpt for SLD int
  real(r8), intent(in)   :: v3sld (plond,plev ,beglatex:endlatex)    ! v3 inpt for SLD int

  real(r8), intent(in)   :: etadot(plond,plevp,beglatex:endlatex)    ! Vertical motion
  real(r8), intent(in)   :: u3p   (plon,plev)                        ! u3 wind at current time
  real(r8), intent(in)   :: v3p   (plon,plev)                        ! v3 wind at current time
  real(r8), intent(in)   :: etadotp(plon,plevp)                      ! Vertical motion at current time
  real(r8), intent(in)   :: dlam  (platd)                            ! d-long for each lat
  integer , intent(in)   :: nlon                                     ! # of longitudes
#if (defined UNICOSMP)
  integer,external :: ismin
#endif
!
!---------------------------Local workspace-----------------------------
!
  real(r8) fac                  ! 1-eps
  real(r8) tmp                  ! temp variable
  real(r8) ump    (plon,plev)   ! interpolated u field
  real(r8) vmp    (plon,plev)   ! interpolated v field
  real(r8) wmp    (plon,plev)   ! interpolated w field
  real(r8) wmpa   (plon,plev)   ! interpolated w field (arrival level)
  real(r8) sigpr  (plon,plev)   ! trajectory increment (vertical)
  real(r8) etamin               ! minimum eta
  real(r8) etamax               ! maximum eta
  real(r8) delsig (plev)        ! dist between dep point and model levs
  real(r8) zt     (plon,plev)   ! top vertical interpolation weight 
  real(r8) zb     (plon,plev)   ! bot vertical interpolation weight 
  real(r8) zttmp                ! top vertical interpolation weight 
  real(r8) zbtmp                ! bot vertical interpolation weight 
  real(r8) rdx    (platd)       ! reciprocal of del-x

  logical locgeo                ! local geodesic flag
  logical lvert                 ! flag to compute vertical trajectory
  logical larrival(plon)        ! flag to indicate whether or not to
!                               ! place the alternative dep point at the
!                               ! arrival point.

  integer kount                 ! index counter
  integer kkmin                 ! index of minimum "delsig"
  integer i                     ! |
  integer ii                    ! |
  integer j                     ! |
  integer jj2                   ! | - indices
  integer jj3                   ! |
  integer k                     ! |
  integer kk                    ! |
  integer n                     ! |
!
!-----------------------------------------------------------------------
!
  fac     = 1._r8 - 10._r8*epsilon (fac)
  locgeo  = .false.
  if(abs(phi(jcen)) .ge. phigs) locgeo = .true.
!
! Interpolate arrival and departure vertical velocities
!
  do k = 1,plev
     do i = 1,nlon
        zttmp = ( etaint(k+1) - sigdp(i,k) )/detai(k)
        zbtmp = 1._r8 - zttmp
        wmpa(i,k) = etadotp(i,k           )*zttmp &
                  + etadotp(i,k+1         )*zbtmp
        wmp (i,k) = etadot(i1+i-1,k  ,jcen)*zttmp &
                  + etadot(i1+i-1,k+1,jcen)*zbtmp
     end do
  end do
!
! Set up computation of trajectory
!
  do k = 1,plev
     do i = 1,nlon
        ii = i + i1 - 1
!
! Place u/v on unit sphere
!
        ump(i,k) = 0.5_r8*(u3sld(ii,k,jcen) + u3p(i,k))*ra
        vmp(i,k) = 0.5_r8*(v3sld(ii,k,jcen) + v3p(i,k))*ra
        wmp(i,k) = 0.5_r8*(wmp(i,k)         + wmpa(i,k)       )
!
! Estimate departure point of parcel trajectory.
!
        lampr(i,k) = -ztodt*ump(i,k)
        phipr(i,k) = -ztodt*vmp(i,k)
        sigpr(i,k) = -ztodt*wmp(i,k)
        if (.not. locgeo) lampr(i,k) = lampr(i,k)/cos( phidp(i,k) )
!
! Initialize winds for use in g.c. calculations.
!
        if (locgeo) then
           upr(i,k) = ump(i,k)
           vpr(i,k) = vmp(i,k)
        end if
     end do
  end do
!
! Estimate initial trajectory departure points
!
  lvert = .true.
  call trajdp(lat     ,jcen    ,ztodt   ,locgeo  ,lvert   , &
              lam(1,jcen),phi  ,etamid  ,upr     ,vpr     , &
              lampr   ,phipr   ,sigpr   ,lamdp   ,phidp   , &
              sigdp   ,nlon    )
!
! Loop over departure point iterates.
!
  do n = 1,iter
!
! Determine departure point indicies and interpolate u,v,w
!
     call bandij (dlam    ,phi     ,lamdp   ,phidp   ,idp     , &
                  jdp     ,nlon    )
     call kdpfnd (plev    ,pmap    ,etamid  ,sigdp   ,kdpmpf  , &
                  kdp     ,nlon    )
!
! Compute weights for x,y,z dimensions
!
     do j = 1,platd
        rdx(j) = 1._r8/(lam(nxpt+2,j) - lam(nxpt+1,j))
     end do
     do k = 1,plev
        do i = 1,nlon
           jj2 = jdp(i,k)
           jj3 = jdp(i,k) + 1
           xl(i,k,2) = (lam(idp(i,k,2)+1,jj2) - lamdp(i,k))*rdx(jj2)
           xl(i,k,3) = (lam(idp(i,k,3)+1,jj3) - lamdp(i,k))*rdx(jj3)
           xr(i,k,2) = 1._r8 - xl(i,k,2)
           xr(i,k,3) = 1._r8 - xl(i,k,3)
           ys(i,k)   = ( phi(jdp(i,k)+1) - phidp(i,k) )/dphi(jdp(i,k))
           yn(i,k)   = 1._r8 - ys(i,k)
           zt(i,k)   = (etamid(kdp(i,k)+1)-sigdp(i,k))/detam(kdp(i,k))
           zb(i,k)   = 1._r8 - zt(i,k)
        end do
     end do
!
     call sltlinint(plev    ,1       ,u3sld(:,:,beglatex:endlatex),xl      ,xr      , &
                    ys      ,yn      ,zt                 ,zb      ,idp     , &
                    jdp     ,kdp     ,ump                ,nlon    )
     call sltlinint(plev    ,1       ,v3sld(:,:,beglatex:endlatex),xl      ,xr      , &
                    ys      ,yn      ,zt                 ,zb      ,idp     , &
                    jdp     ,kdp     ,vmp                ,nlon    )
     call kdpfnd (plevp   ,pmap    ,etaint  ,sigdp   ,kdpmph  , &
                  kdp     ,nlon    )
     do k = 1,plev
        do i = 1,nlon
           zt(i,k)   = (etaint(kdp(i,k)+1)-sigdp(i,k))/detai(kdp(i,k))
           zb(i,k)   = 1._r8 - zt(i,k)
        end do
     end do
     call sltlinint(plevp   ,1       ,etadot(:,:,beglatex:endlatex),xl      ,xr      , &
                    ys      ,yn      ,zt                       ,zb      ,idp     , &
                    jdp     ,kdp     ,wmp                      ,nlon    )
!
! Add arrival point velocities (n) to interpolated velocities
!
     do k = 1,plev
        do i = 1,nlon
           ii = i + i1 - 1
           ump(i,k) = 0.5_r8*(ump(i,k) + u3p(i,k))
           vmp(i,k) = 0.5_r8*(vmp(i,k) + v3p(i,k))
           wmp(i,k) = 0.5_r8*(wmp(i,k) + wmpa(i,k)       )
        end do
     end do
!
! Compute new trajectory departure points
!
     lvert = .true.
     call depinc (jcen    ,ztodt   ,ra      ,locgeo  ,lvert   , &
                  lam(1,jcen),phi  ,ump     ,vmp     ,wmp     , &
                  upr     ,vpr     ,lamdp   ,phidp   ,lampr   , &
                  phipr   ,sigpr   ,nlon    )
     call trajdp (lat     ,jcen    ,ztodt   ,locgeo  ,lvert   , &
                  lam(1,jcen),phi  ,etamid  ,upr     ,vpr     , &
                  lampr   ,phipr   ,sigpr   ,lamdp   ,phidp   , &
                  sigdp   ,nlon    )
  end do
!
! Compile departure point binning statistics
!
  do k = 1,plev
     if(k .eq. 1) then
        etamin = etamid(k  )
        etamax = etamid(k+1)
     elseif(k .eq. plev) then
        etamin = etamid(k-1)
        etamax = etamid(k  )
     else
        etamin = etamid(k-1)
        etamax = etamid(k+1)
     end if
     do i = 1,nlon
        larrival(i) = sigdp(i,k) .ge. etamin .and. sigdp(i,k) .le. etamax
     end do
!
! If:  Departure point is within one grid interval from the arrival
! point, bin it in the ARRIVAL point bin
!
     tmp = 0._r8
     do i = 1,nlon
        if (larrival(i)) then
           tmp           = tmp + 1._r8
        end if
     end do
     levkntl(k,k,lat) = levkntl(k,k,lat) + tmp
!
! Else:  find departure point bin
!
     kount = 0
     do i=1,nlon
        if (larrival(i)) kount = kount + 1
     end do
     if (kount .ne. nlon) then
        do i = 1,nlon
           if (.not. larrival(i)) then
              delsig(1) = 1.e+35_r8
              do kk = 2,plev
                 delsig(kk) = abs( sigdp(i,k) - etaint(kk) )
              end do
              kkmin = ismin( plev,delsig,1 )
              if(kkmin .gt. k) kkmin = kkmin-1
              levkntl(kkmin,k,lat) = levkntl(kkmin,k,lat) + 1._r8
           end if
        end do
     end if
  end do
!
  return
end subroutine slttraj

!
!-----------------------------------------------------------------------
!

subroutine sltwgts(limdrh  ,limdrv  ,lhrzwgt ,lvrtwgt ,kdim    , &
                   idp     ,jdp     ,kdp     ,lam     ,phi     , &
                   z       ,dphi    ,dz      ,lamdp   ,phidp   , &
                   sigdp   ,lbasiy  ,wiz     ,kkdp    ,xl      , &
                   xr      ,wgt1x   ,wgt2x   ,wgt3x   ,wgt4x   , &
                   hl      ,hr      ,dhl     ,dhr     ,ys      , &
                   yn      ,wgt1y   ,wgt2y   ,wgt3y   ,wgt4y   , &
                   hs      ,hn      ,dhs     ,dhn     ,rdphi   , &
                   wgt1z   ,wgt2z   ,wgt3z   ,wgt4z   ,hb      , &
                   ht      ,dhb     ,dht     ,rdz     ,zt      , &
                   zb      ,nlon    )
!-----------------------------------------------------------------------
!
! Purpose: 
! Compute weights for SLT interpolation
!
! Author:  J. Olson
!
!-----------------------------------------------------------------------
#if (!defined UNICOSMP)
  use srchutil, only: whenieq
#endif
!------------------------------Arguments--------------------------------
!
  logical , intent(in)   :: limdrh              ! horizontal derivative limiter flag
  logical , intent(in)   :: limdrv              ! vertical   derivative limiter flag
  logical , intent(in)   :: lhrzwgt             ! flag to compute horizontal weights
  logical , intent(in)   :: lvrtwgt             ! flag to compute vertical   weights

  integer , intent(in)   :: kdim                ! vertical coordinate
  integer , intent(in)   :: idp   (plon,plev,4) ! index of x-coordinate of dep pt
  integer , intent(in)   :: jdp   (plon,plev)   ! index of y-coordinate of dep pt
  integer , intent(in)   :: kdp   (plon,plev)   ! index of z-coordinate of dep pt

  real(r8), intent(in)   :: lam   (plond,platd) ! longitude coordinates of model grid
  real(r8), intent(in)   :: phi   (platd)       ! latitude  coordinates of model grid
  real(r8), intent(in)   :: z     (kdim)        ! vertical  coordinates of model grid
  real(r8), intent(in)   :: dphi  (platd)       ! latitudinal grid increments
  real(r8), intent(in)   :: dz    (kdim)        ! vertical grid increments
  real(r8), intent(in)   :: lamdp (plon,plev)   ! x-coordinates of dep pts.
  real(r8), intent(in)   :: phidp (plon,plev)   ! y-coordinates of dep pt
  real(r8), intent(in)   :: sigdp (plon,plev)   ! z-coordinates of dep pt
  real(r8), intent(in)   :: lbasiy(4,2,platd)   ! y-interpolation weights for Lag. cubic
  real(r8), intent(in)   :: wiz   (4,2,kdim)    ! z-interpolation weights for Lag. cubic
  integer , intent(out)  :: kkdp  (plon,plev)   ! index of z-coordinate of dep pt (alt)

  real(r8), intent(out)  :: xl    (plon,plev,4) ! weight for x-interpolants (left)
  real(r8), intent(out)  :: xr    (plon,plev,4) ! weight for x-interpolants (right)
  real(r8), intent(out)  :: wgt1x (plon,plev,4) ! weight for x-interpolants (Lag Cubic)
  real(r8), intent(out)  :: wgt2x (plon,plev,4) ! weight for x-interpolants (Lag Cubic)
  real(r8), intent(out)  :: wgt3x (plon,plev,4) ! weight for x-interpolants (Lag Cubic)
  real(r8), intent(out)  :: wgt4x (plon,plev,4) ! weight for x-interpolants (Lag Cubic)
  real(r8), intent(out)  :: hl    (plon,plev,4) ! weight for x-interpolants (Hermite)
  real(r8), intent(out)  :: hr    (plon,plev,4) ! weight for x-interpolants (Hermite)
  real(r8), intent(out)  :: dhl   (plon,plev,4) ! weight for x-interpolants (Hermite)
  real(r8), intent(out)  :: dhr   (plon,plev,4) ! weight for x-interpolants (Hermite)

  real(r8), intent(out)  :: ys    (plon,plev)   ! weight for y-interpolants (south)
  real(r8), intent(out)  :: yn    (plon,plev)   ! weight for y-interpolants (north)
  real(r8), intent(out)  :: wgt1y (plon,plev)   ! weight for y-interpolants (Lag Cubic)
  real(r8), intent(out)  :: wgt2y (plon,plev)   ! weight for y-interpolants (Lag Cubic)
  real(r8), intent(out)  :: wgt3y (plon,plev)   ! weight for y-interpolants (Lag Cubic)
  real(r8), intent(out)  :: wgt4y (plon,plev)   ! weight for y-interpolants (Lag Cubic)
  real(r8), intent(out)  :: hs    (plon,plev)   ! weight for y-interpolants (Hermite)
  real(r8), intent(out)  :: hn    (plon,plev)   ! weight for y-interpolants (Hermite)
  real(r8), intent(out)  :: dhs   (plon,plev)   ! weight for y-interpolants (Hermite)
  real(r8), intent(out)  :: dhn   (plon,plev)   ! weight for y-interpolants (Hermite)
  real(r8), intent(out)  :: rdphi (plon,plev)   ! reciprocal of y-interval

  real(r8), intent(out)  :: wgt1z (plon,plev)   ! weight for z-interpolants (Lag Cubic)
  real(r8), intent(out)  :: wgt2z (plon,plev)   ! weight for z-interpolants (Lag Cubic)
  real(r8), intent(out)  :: wgt3z (plon,plev)   ! weight for z-interpolants (Lag Cubic)
  real(r8), intent(out)  :: wgt4z (plon,plev)   ! weight for z-interpolants (Lag Cubic)
  real(r8), intent(out)  :: hb    (plon,plev)   ! weight for z-interpolants (Hermite)
  real(r8), intent(out)  :: ht    (plon,plev)   ! weight for z-interpolants (Hermite)
  real(r8), intent(out)  :: dhb   (plon,plev)   ! weight for z-interpolants (Hermite)
  real(r8), intent(out)  :: dht   (plon,plev)   ! weight for z-interpolants (Hermite)
  real(r8), intent(out)  :: rdz   (plon,plev)   ! reciprocal of z-interval
  real(r8), intent(out)  :: zt    (plon,plev)   ! linear interpolation weight
  real(r8), intent(out)  :: zb    (plon,plev)   ! linear interpolation weight
  integer , intent(in)   :: nlon                ! number of longitudes for this latitude
!
!---------------------------Local workspace-----------------------------
!
  integer i                 ! |
  integer ii,jj(plon,plev,4)! |
  integer k,n,j             ! |
  integer icount            ! |
  integer jdpval            ! |
  integer kdpval            ! | -- indices
  integer jmin              ! |
  integer jmax              ! |
  integer kdimm2            ! |
  integer nval              ! |
  integer indx (plond)      ! |
!
  real(r8) dx  (platd)      ! |
  real(r8) rdx (platd)      ! |
  real(r8) dyj              ! |
  real(r8) tmp1             ! |
  real(r8) tmp2             ! |
  real(r8) tmp3             ! |
  real(r8) tmp4             ! | -- tmp variables
  real(r8) dzk              ! |
  real(r8) denom1           ! |
  real(r8) denom2           ! |
  real(r8) denom3           ! |
  real(r8) denom4           ! |
  real(r8) coef12           ! |
  real(r8) coef34           ! |
!
!-----------------------------------------------------------------------
!
  denom1 = -1._r8/6._r8
  denom2 =  0.5_r8
  denom3 = -0.5_r8
  denom4 =  1._r8/6._r8
!
! HORIZONTAL weights
!
  if (lhrzwgt) then
!
! Determine N/S extent of all departure points in this latitude slice
!
     jmin =  1000000
     jmax = -1000000
     do k=1,plev
        do i=1,nlon
           if(jdp(i,k) .lt. jmin) jmin = jdp(i,k)
           if(jdp(i,k) .gt. jmax) jmax = jdp(i,k)
        end do
     end do
!
! Compute weights for x-direction
!
     do j = 1,platd
        dx (j) = lam(nxpt+2,j) - lam(nxpt+1,j)
        rdx(j) = 1._r8/dx(j)
     end do
!
     do n=1,4
        do k=1,plev
           do i=1,nlon
              jj(i,k,n) = jdp(i,k) - 2 + n
              xl(i,k,n) = (lam(idp(i,k,n)+1,jj(i,k,n)) - lamdp(i,k))*rdx(jj(i,k,n))
              xr(i,k,n) = 1._r8 - xl(i,k,n)
           end do
        end do
     end do
     do n=2,3
           do k=1,plev
              do i=1,nlon
                 hl (i,k,n)   = ( 3.0_r8 - 2.0_r8*xl(i,k,n) )*xl(i,k,n)**2
                 hr (i,k,n)   = ( 3.0_r8 - 2.0_r8*xr(i,k,n) )*xr(i,k,n)**2
                 dhl(i,k,n)   = -dx(jj(i,k,n))*( xl(i,k,n) - 1._r8 )*xl(i,k,n)**2
                 dhr(i,k,n)   =  dx(jj(i,k,n))*( xr(i,k,n) - 1._r8 )*xr(i,k,n)**2
              end do
           end do
           do k=1,plev
              do i=1,nlon
                 tmp1         =  xr(i,k,n) + 1._r8
                 tmp4         =  xr(i,k,n) - 2._r8
                 coef12       = -xl(i,k,n)*tmp4
                 coef34       =  xr(i,k,n)*tmp1
                 wgt1x(i,k,n) =  denom1*coef12*xr(i,k,n)
                 wgt2x(i,k,n) =  denom2*coef12*tmp1
                 wgt3x(i,k,n) =  denom3*coef34*tmp4
                 wgt4x(i,k,n) = -denom4*coef34*xl(i,k,n)
              end do
           end do
     end do
!
! Compute weights for y-direction
!
     icount = 0
     do jdpval=jmin,jmax
        do k=1,plev
           call whenieq(nlon    ,jdp(1,k),1       ,jdpval  ,indx    ,nval    )
           icount = icount + nval
           dyj    = dphi(jdpval)
           do ii = 1,nval
              i = indx(ii)
              ys(i,k) = ( phi(jdpval+1) - phidp(i,k) )/dyj
              yn(i,k) = 1._r8 - ys(i,k)
           end do
              do ii = 1,nval
                 i = indx(ii)
                 rdphi(i,k) = 1._r8/dyj
                 hs   (i,k) = ( 3.0_r8 - 2.0_r8*ys(i,k) )*ys(i,k)**2
                 hn   (i,k) = ( 3.0_r8 - 2.0_r8*yn(i,k) )*yn(i,k)**2
                 dhs  (i,k) = -dyj*( ys(i,k) - 1._r8 )*ys(i,k)**2
                 dhn  (i,k) =  dyj*( yn(i,k) - 1._r8 )*yn(i,k)**2
              end do
              do ii = 1,nval
                 i          = indx(ii)
                 tmp1       = phidp(i,k) - lbasiy(1,1,jdpval)
                 tmp2       = phidp(i,k) - lbasiy(2,1,jdpval)
                 tmp3       = phidp(i,k) - lbasiy(3,1,jdpval)
                 tmp4       = phidp(i,k) - lbasiy(4,1,jdpval)
                 coef12     = tmp3*tmp4   
                 coef34     = tmp1*tmp2   
                 wgt1y(i,k) = coef12*tmp2*lbasiy(1,2,jdpval)
                 wgt2y(i,k) = coef12*tmp1*lbasiy(2,2,jdpval)
                 wgt3y(i,k) = coef34*tmp4*lbasiy(3,2,jdpval)
                 wgt4y(i,k) = coef34*tmp3*lbasiy(4,2,jdpval)
              end do
        end do
     end do
     if (icount.ne.nlon*plev) then
        call endrun ('SLTWGTS:  Did not complete computations for all departure points')
     end if
  end if
!
! VERTICAL weights
!
  if (lvrtwgt) then
!
! Limit kdp to between "2" and "kdim-2" when computing weights and
! derivatives.
!
     kdimm2 = kdim - 2
     do k=1,plev
        do i=1,nlon
           kkdp(i,k) = min0( kdimm2,max0( 2,kdp(i,k) ) )
           dzk       = dz(kdp(i,k))
           rdz(i,k)  = 1._r8/dzk
           zt(i,k)   = ( z  (kdp(i,k)+1) - sigdp(i,k) )/dzk
           zb(i,k)   = 1._r8 - zt(i,k)
           ht (i,k)  = ( 3.0_r8 - 2.0_r8*zt(i,k) )*zt(i,k)**2
           hb (i,k)  = ( 3.0_r8 - 2.0_r8*zb(i,k) )*zb(i,k)**2
           dht(i,k)  = -dzk*( zt(i,k) - 1._r8 )*zt(i,k)**2
           dhb(i,k)  =  dzk*( zb(i,k) - 1._r8 )*zb(i,k)**2
        end do
     end do
!
     if(.not. limdrv) then
        icount = 0
        do kdpval=2,kdimm2
           do k=1,plev
              call whenieq(nlon    ,kkdp(1,k),1       ,kdpval  ,indx    ,nval    )
              icount = icount + nval
              do ii = 1,nval
                 i          = indx(ii)
                 tmp1       = sigdp(i,k) -  wiz(1,1,kdpval)
                 tmp2       = sigdp(i,k) -  wiz(2,1,kdpval)
                 tmp3       = sigdp(i,k) -  wiz(3,1,kdpval)
                 tmp4       = sigdp(i,k) -  wiz(4,1,kdpval)
                 coef12     = tmp3*tmp4   
                 coef34     = tmp1*tmp2   
                 wgt1z(i,k) = coef12*tmp2*wiz(1,2,kdpval)
                 wgt2z(i,k) = coef12*tmp1*wiz(2,2,kdpval)
                 wgt3z(i,k) = coef34*tmp4*wiz(3,2,kdpval)
                 wgt4z(i,k) = coef34*tmp3*wiz(4,2,kdpval)
              end do
           end do
        end do
        if (icount.ne.nlon*plev) then
           call endrun ('SLTWGTS:  Did not complete computations for all departure points')
        end if
     end if
  end if
!
  return
end subroutine sltwgts

!
!-----------------------------------------------------------------------
!

subroutine trajdp(lat     ,jcen    ,dt      ,locgeo  ,lvert   , &
                  lam     ,phi     ,etamid  ,upr     ,vpr     , &
                  lampr   ,phipr   ,sigpr   ,lamdp   ,phidp   , &
                  sigdp   ,nlon    )
!
!-----------------------------------------------------------------------
!
! Purpose:
! Determine trajectory departure point for each arrival point based upon
! departure point increments
!
! Author:  J. Olson
!
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
!
  integer , intent(in)   :: lat               ! index of lat slice (model)
  integer , intent(in)   :: jcen              ! index of lat slice(extend)
  real(r8), intent(in)   :: dt                ! time step (seconds)
  logical , intent(in)   :: locgeo            ! local geodesic flag
  logical , intent(in)   :: lvert             ! flag to compute vertical trajectory
  real(r8), intent(in)   :: lam   (plond)     ! long. coord of model grid
  real(r8), intent(in)   :: phi   (platd)     ! lat.  coord of model grid
  real(r8), intent(in)   :: etamid(plev)      ! vertical full levels
  real(r8), intent(in)   :: upr   (plon,plev) ! interpolated u field (local geodesic)
  real(r8), intent(in)   :: vpr   (plon,plev) ! interpolated v field (local geodesic)
  real(r8), intent(in)   :: lampr (plon,plev) ! trajectory increment (x-direction)
  real(r8), intent(in)   :: phipr (plon,plev) ! trajectory increment (y-direction)
  real(r8), intent(in)   :: sigpr (plon,plev) ! trajectory increment (vertical)
  real(r8), intent(out)  :: lamdp (plon,plev) ! zonal      departure pt. coord.
  real(r8), intent(out)  :: phidp (plon,plev) ! meridional departure pt. coord.
  real(r8), intent(out)  :: sigdp (plon,plev) ! vertical   departure pt. coord.
  integer , intent(in)   :: nlon              ! number of longitudes for this latitude
!
!---------------------------Local workspace-----------------------------
!
  real(r8) fac              ! 1 - eps
  real(r8) botlim           ! bottom limit of the trajectory
  real(r8) phi0             ! Current latitude (radians)
  real(r8) cphi0            ! cos latitude
  real(r8) sphi0            ! sin latitude
  real(r8) dist             ! computed distance**2
  real(r8) distmx           ! max distance**2
  real(r8) pi2              ! pi/2
  real(r8) phipi2           ! pi/2
  real(r8) sgnphi0          ! holds sign of phi0
  real(r8) pi               ! pi
  real(r8) twopi            ! 2*pi
!
  real(r8) clamgc           ! |
  real(r8) cphigc           ! |
  real(r8) dlamsc           ! | 
  real(r8) slamgc           ! | -- temporary variables
  real(r8) slam2            ! |
  real(r8) sphigc           ! |
!
  integer i                 ! |
  integer ii                ! |
  integer k                 ! | -- indices
  integer ibad              ! |
  integer kbad              ! |
!
  logical lstop             ! flag to stop run if departure point is over the pole
!
!-----------------------------------------------------------------------
!
  fac     = 1._r8 - 10._r8*epsilon(fac)
  pi      = 4._r8*atan(1._r8)
  twopi   = 2._r8*pi
  pi2     = pi/2._r8
  phi0    = phi(jcen)
  cphi0   = cos( phi0 )
  sphi0   = sin( phi0 )
  sgnphi0 = sign( 1._r8, phi0 )
  distmx  = (sign(pi2,phi0) - phi0)/(1.1_r8*dt)
  distmx  = distmx*distmx
!
! Compute coordinates of departure point.
!
  lstop = .false.
  do k = 1,plev
!
! if near pole, convert from g.c. to spherical
!
     if (locgeo) then
        do i = 1,nlon
           ii = i + i1 - 1
           sphigc = sin( phipr(i,k) )
           cphigc = cos( phipr(i,k) )
           slamgc = sin( lampr(i,k) )
           clamgc = cos( lampr(i,k) )
           phidp(i,k)  = asin((sphigc*cphi0 + cphigc*sphi0*clamgc)*fac)
           if ( abs(phidp(i,k)) .ge. phi(j1+plat)*fac ) &
                                          phidp(i,k) = sign( phi(j1+plat),phidp(i,k) )*fac
           dlamsc = asin((slamgc*cphigc/cos(phidp(i,k)))*fac)
!
! If traj is over pole, check for proper branch of arcsin
!
           dist   = upr(i,k)**2 + vpr(i,k)**2
           if( dist .gt. distmx) then
              slam2  = slamgc**2
              phipi2 = asin((sqrt((slam2 - 1._r8)/(slam2 - 1._r8/cphi0**2)))*fac)
              if(sgnphi0*phipr(i,k) .gt. phipi2) dlamsc = sign(pi,lampr(i,k)) - dlamsc
           end if
           lamdp(i,k) = lam(ii) + dlamsc
        end do
     else
        do i = 1,nlon
           ii = i + i1 - 1
           lamdp(i,k) = lam(ii) + lampr(i,k)
           phidp(i,k) = phi0    + phipr(i,k)
!
! If traj is over pole, STOP
!
           if (phidp(i,k) >= phi(j1+plat)) then
              lstop = .true.
              ibad  = i
              kbad  = k
           end if
           if (phidp(i,k) < phi(j1-1   )) then
              lstop = .true.
              ibad  = i
              kbad  = k
           end if
        end do
     end if
#if ( defined SPMD )
!
! If traj goes off-processor (out-of-bounds), STOP
!
        do i = 1,nlon
           if (phidp(i,k) >= phi(endlatex-nxpt) ) then
              lstop = .true.
              ibad  = i
              kbad  = k
           end if
           if (phidp(i,k) < phi(beglatex+nxpt) ) then
              lstop = .true.
              ibad  = i
              kbad  = k
           end if
        end do
#endif
!
! Apply appropriate limiters
!
     do i = 1,nlon
        if(lamdp(i,k) .ge. twopi) lamdp(i,k) = lamdp(i,k) - twopi + 10._r8*epsilon(lamdp)
        if(lamdp(i,k) .lt. 0._r8) lamdp(i,k) = lamdp(i,k) + twopi - 10._r8*epsilon(lamdp)
     end do
!
! Compute vertical departure points and apply appropriate limiters
!
     if (lvert) then
        botlim  = etamid(plev)*fac
        do i = 1,nlon
           sigdp(i,k) = etamid(k) + sigpr(i,k)
           if(sigdp(i,k) .lt. etamid(   1)) sigdp(i,k) = etamid(1)
           if(sigdp(i,k) .ge. etamid(plev)) sigdp(i,k) = botlim
        end do
     end if
  end do
!
! Test that the latitudinal extent of trajectory is NOT over the poles
!
  if (lstop) then
     write(iulog,*)'ERROR IN TRAJDP: ****** MODEL IS BLOWING UP *********'
     write(iulog,9000) ibad,kbad,lat
     call endrun
  end if
!
  return
!
! Formats
!
9000 format(//'Departure point out of bounds.  Parcel associated with longitude ' &
              ,i5,', level ',i5, ' and latitude ',i5/' is outside the model domain. ')
end subroutine trajdp

!
!-----------------------------------------------------------------------
!

end module scanslt
