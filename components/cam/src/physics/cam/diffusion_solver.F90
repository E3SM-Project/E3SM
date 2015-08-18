
  module diffusion_solver

  !------------------------------------------------------------------------------------ !
  ! Module to solve vertical diffusion equations using a tri-diagonal solver.           !
  ! The module will also apply countergradient fluxes, and apply molecular              ! 
  ! diffusion for constituents.                                                         !
  !                                                                                     !
  ! Public interfaces :                                                                 ! 
  !    init_vdiff       initializes time independent coefficients                       !
  !    compute_vdiff    solves diffusion equations                                      !
  !    vdiff_selector   type for storing fields selected to be diffused                 !
  !    vdiff_select     selects fields to be diffused                                   !
  !    operator(.not.)  extends .not. to operate on type vdiff_selector                 !
  !    any              provides functionality of intrinsic any for type vdiff_selector !
  !                                                                                     !
  !------------------------------------ Code History ---------------------------------- !
  ! Initial subroutines :  B. Boville and others, 1991-2004                             !
  ! Modularization      :  J. McCaa, September 2004                                     !
  ! Most Recent Code    :  Sungsu Park, Aug. 2006, Dec. 2008, Jan. 2010.                !
  !------------------------------------------------------------------------------------ !
    use module_perturb

  implicit none
  private       
  save

  integer, parameter :: r8 = selected_real_kind(12)      ! 8 byte real

  ! ----------------- !
  ! Public interfaces !
  ! ----------------- !

  public init_vdiff                                      ! Initialization
  public new_fieldlist_vdiff                             ! Returns an empty fieldlist
  public compute_vdiff                                   ! Full routine
  public vdiff_selector                                  ! Type for storing fields selected to be diffused
  public vdiff_select                                    ! Selects fields to be diffused
  public operator(.not.)                                 ! Extends .not. to operate on type vdiff_selector
  public any                                             ! Provides functionality of intrinsic any for type vdiff_selector
 
  ! Below stores logical array of fields to be diffused

  type vdiff_selector 
       private
       logical, allocatable, dimension(:) :: fields
  end type vdiff_selector

  ! Below extends .not. to operate on type vdiff_selector

  interface operator(.not.)
       module procedure not
  end interface

  ! Below provides functionality of intrinsic any for type vdiff_selector

  interface any                           
       module procedure my_any
  end interface

  ! ------------ !
  ! Private data !
  ! ------------ !

  ! Unit number for log output
  integer :: iulog = -1

  real(r8), private   :: cpair                           ! Specific heat of dry air
  real(r8), private   :: gravit                          ! Acceleration due to gravity
  real(r8), private   :: rair                            ! Gas constant for dry air
  real(r8), private   :: zvir                            ! rh2o/rair - 1
  real(r8), private   :: latvap                          ! Latent heat of vaporization
  real(r8), private   :: karman                          ! von Karman constant

  logical,  private   :: do_iss                          ! Use implicit turbulent surface stress computation

  ! Parameters used for Turbulent Mountain Stress

  real(r8), parameter :: z0fac   = 0.025_r8              ! Factor determining z_0 from orographic standard deviation
  real(r8), parameter :: z0max   = 100._r8               ! Max value of z_0 for orography
  real(r8), parameter :: horomin = 10._r8                ! Min value of subgrid orographic height for mountain stress
  real(r8), parameter :: dv2min  = 0.01_r8               ! Minimum shear squared
  real(r8), private   :: oroconst                        ! Converts from standard deviation to height

  contains

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine init_vdiff( kind, iulog_in, rair_in, gravit_in, do_iss_in, &
                         errstring )

    integer,              intent(in)  :: kind            ! Kind used for reals
    integer,              intent(in)  :: iulog_in        ! Unit number for log output.
    real(r8),             intent(in)  :: rair_in         ! Input gas constant for dry air
    real(r8),             intent(in)  :: gravit_in       ! Input gravitational acceleration
    logical,              intent(in)  :: do_iss_in       ! Input ISS flag
    character(128),       intent(out) :: errstring       ! Output status
    
    errstring = ''
    iulog = iulog_in
    if( kind .ne. r8 ) then
        write(iulog,*) 'KIND of reals passed to init_vdiff -- exiting.'
        errstring = 'init_vdiff'
        return
    endif

    rair   = rair_in     
    gravit = gravit_in 
    do_iss = do_iss_in

  end subroutine init_vdiff

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  type(vdiff_selector) pure function new_fieldlist_vdiff(ncnst)

    integer,              intent(in)  :: ncnst           ! Number of constituents

    allocate( new_fieldlist_vdiff%fields( 3 + ncnst ) )
    new_fieldlist_vdiff%fields = .false.

  end function new_fieldlist_vdiff

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine compute_vdiff( lchnk           ,                                                                   &
                            pcols           , pver               , ncnst         , ncol         , pmid        , &
                            pint            , rpdel              , t             , ztodt        , taux        , &
                            tauy            , shflx              , cflx          , ntop         , nbot        , &
                            kvh             , kvm                , kvq           , cgs          , cgh         , &
                            zi              , ksrftms            , qmincg        , fieldlist    , fieldlistm  , &
                            u               , v                  , q             , dse          ,               &
                            tautmsx         , tautmsy            , dtk           , topflx       , errstring   , &
                            tauresx         , tauresy            , itaures       , cpairv       , rairi       , &
                            do_molec_diff  , compute_molec_diff, vd_lu_qdecomp, kvt, flb )

    !-------------------------------------------------------------------------- !
    ! Driver routine to compute vertical diffusion of momentum, moisture, trace !
    ! constituents and dry static energy. The new temperature is computed from  !
    ! the diffused dry static energy.                                           ! 
    ! Turbulent diffusivities and boundary layer nonlocal transport terms are   !
    ! obtained from the turbulence module.                                      !
    !-------------------------------------------------------------------------- !

!    Used for CAM debugging.
!    use phys_debug_util,    only : phys_debug_col
!    use time_manager,       only : is_first_step, get_nstep
    use vdiff_lu_solver,     only : lu_decomp, vd_lu_decomp, vd_lu_solve
  
  ! Modification : Ideally, we should diffuse 'liquid-ice static energy' (sl), not the dry static energy.
  !                Also, vertical diffusion of cloud droplet number concentration and aerosol number
  !                concentration should be done very carefully in the future version.

    ! --------------- !
    ! Input Arguments !
    ! --------------- !
    integer, optional :: flb
    integer,  intent(in)    :: lchnk                   
    integer,  intent(in)    :: pcols
    integer,  intent(in)    :: pver
    integer,  intent(in)    :: ncnst
    integer,  intent(in)    :: ncol                      ! Number of atmospheric columns
    integer,  intent(in)    :: ntop                      ! Top    interface level to which vertical diffusion is applied ( = 1 ).
    integer,  intent(in)    :: nbot                      ! Bottom interface level to which vertical diffusion is applied ( = pver ).
    integer,  intent(in)    :: itaures                   ! Indicator determining whether 'tauresx,tauresy'
                                                         ! is updated (1) or non-updated (0) in this subroutine.

    real(r8), intent(in)    :: pmid(pcols,pver)          ! Mid-point pressures [ Pa ]
    real(r8), intent(in)    :: pint(pcols,pver+1)        ! Interface pressures [ Pa ]
    real(r8), intent(in)    :: rpdel(pcols,pver)         ! 1./pdel
    real(r8), intent(in)    :: t(pcols,pver)             ! Temperature [ K ]
    real(r8), intent(in)    :: ztodt                     ! 2 delta-t [ s ]
    real(r8), intent(in)    :: taux(pcols)               ! Surface zonal      stress.
                                                         ! Input u-momentum per unit time per unit area into the atmosphere [ N/m2 ]
    real(r8), intent(in)    :: tauy(pcols)               ! Surface meridional stress.
                                                         ! Input v-momentum per unit time per unit area into the atmosphere [ N/m2 ]
    real(r8), intent(in)    :: shflx(pcols)              ! Surface sensible heat flux [ W/m2 ]
    real(r8), intent(in)    :: cflx(pcols,ncnst)         ! Surface constituent flux [ kg/m2/s ]
    real(r8), intent(in)    :: zi(pcols,pver+1)          ! Interface heights [ m ]
    real(r8), intent(in)    :: ksrftms(pcols)            ! Surface drag coefficient for turbulent mountain stress. > 0. [ kg/s/m2 ]
    real(r8), intent(in)    :: qmincg(ncnst)             ! Minimum constituent mixing ratios from cg fluxes
    real(r8), intent(in)    :: cpairv(pcols,pver)        ! Specific heat at constant pressure
    real(r8), intent(in)    :: rairi(pcols,pver+1)       ! Gas constant on interface levels
    real(r8), intent(in)    :: kvh(pcols,pver+1)         ! Eddy diffusivity for heat [ m2/s ]

    logical,  intent(in)    :: do_molec_diff             ! Flag indicating multiple constituent diffusivities

    integer,  external, optional :: compute_molec_diff   ! Constituent-independent moleculuar diffusivity routine
    integer,  external, optional :: vd_lu_qdecomp        ! Constituent-dependent moleculuar diffusivity routine
    real(r8), intent(out), optional :: kvt(pcols,pver+1) ! Kinematic molecular conductivity
    type(vdiff_selector), intent(in) :: fieldlist        ! Array of flags selecting which fields to diffuse
    type(vdiff_selector), intent(in) :: fieldlistm       ! Array of flags selecting which fields for molecular diffusion

    ! ---------------------- !
    ! Input-Output Arguments !
    ! ---------------------- !

    real(r8), intent(inout) :: kvm(pcols,pver+1)         ! Eddy viscosity ( Eddy diffusivity for momentum ) [ m2/s ]
    real(r8), intent(inout) :: kvq(pcols,pver+1)         ! Eddy diffusivity for constituents
    real(r8), intent(inout) :: cgs(pcols,pver+1)         ! Counter-gradient star [ cg/flux ]
    real(r8), intent(inout) :: cgh(pcols,pver+1)         ! Counter-gradient term for heat

    real(r8), intent(inout) :: u(pcols,pver)             ! U wind. This input is the 'raw' input wind to
                                                         ! PBL scheme without iterative provisional update. [ m/s ]
    real(r8), intent(inout) :: v(pcols,pver)             ! V wind. This input is the 'raw' input wind to PBL scheme
                                                         ! without iterative provisional update. [ m/s ]
    real(r8), intent(inout) :: q(pcols,pver,ncnst)       ! Moisture and trace constituents [ kg/kg, #/kg ? ]
    real(r8), intent(inout) :: dse(pcols,pver)           ! Dry static energy [ J/kg ]

    real(r8), intent(inout) :: tauresx(pcols)            ! Input  : Reserved surface stress at previous time step
    real(r8), intent(inout) :: tauresy(pcols)            ! Output : Reserved surface stress at current  time step

    ! ---------------- !
    ! Output Arguments !
    ! ---------------- !

    real(r8), intent(out)   :: dtk(pcols,pver)           ! T tendency from KE dissipation
    real(r8), intent(out)   :: tautmsx(pcols)            ! Implicit zonal      turbulent mountain surface stress
                                                         ! [ N/m2 = kg m/s /s/m2 ]
    real(r8), intent(out)   :: tautmsy(pcols)            ! Implicit meridional turbulent mountain surface stress
                                                         ! [ N/m2 = kg m/s /s/m2 ]
    real(r8), intent(out)   :: topflx(pcols)             ! Molecular heat flux at the top interface
    character(128), intent(out) :: errstring             ! Output status

    ! --------------- !
    ! Local Variables ! 
    ! --------------- !

    integer  :: i, k, m, icol                            ! Longitude, level, constituent indices
    integer  :: status                                   ! Status indicator
    integer  :: nbot_molec                               ! Bottom level where molecular diffusivity is applied
    integer  :: ntop_molec                               ! Top level where molecular diffusivity is applied
    logical  :: lqtst(pcols)                             ! Adjust vertical profiles
    logical  :: need_decomp                              ! Whether to compute a new decomposition
    logical  :: cnst_fixed_ubc(ncnst)                    ! Whether upper boundary condition is fixed
    logical  :: cnst_fixed_ubflx(ncnst)                  ! Whether upper boundary flux is a fixed non-zero value

    ! LU decomposition information.
    type(lu_decomp) :: decomp

    real(r8) :: tmp1(pcols)                              ! Temporary storage
    real(r8) :: tmpi1(pcols,pver+1)                      ! Interface KE dissipation
    real(r8) :: tint(pcols,pver+1)                       ! Interface temperature
    real(r8) :: rhoi(pcols,pver+1)                       ! rho at interfaces
    real(r8) :: tmpi2(pcols,pver+1)                      ! dt*(g*rho)**2/dp at interfaces
    real(r8) :: rrho(pcols)                              ! 1./bottom level density 

    real(r8) :: zero(pcols)                              ! Zero array for surface heat exchange coefficients 
    real(r8) :: tautotx(pcols)                           ! Total surface stress ( zonal )
    real(r8) :: tautoty(pcols)                           ! Total surface stress ( meridional )

    real(r8) :: dinp_u(pcols,pver+1)                     ! Vertical difference at interfaces, input u
    real(r8) :: dinp_v(pcols,pver+1)                     ! Vertical difference at interfaces, input v
    real(r8) :: dout_u                                   ! Vertical difference at interfaces, output u
    real(r8) :: dout_v                                   ! Vertical difference at interfaces, output v
    real(r8) :: dse_top(pcols)                           ! dse on top boundary
    real(r8) :: cc_top(pcols)                            ! Lower diagonal at top interface
    real(r8) :: cd_top(pcols)                            ! 
    real(r8) :: rghd(pcols,pver+1)                       ! (1/H_i - 1/H) *(g*rho)^(-1)

    real(r8) :: qtm(pcols,pver)                          ! Temporary copy of q
    real(r8) :: kq_scal(pcols,pver+1)                    ! kq_fac*sqrt(T)*m_d/rho for molecular diffusivity
    real(r8) :: mw_fac(ncnst)                            ! sqrt(1/M_q + 1/M_d) for this constituent
    real(r8) :: cnst_mw(ncnst)                           ! Molecular weight [ kg/kmole ]
    real(r8) :: ubc_mmr(pcols,ncnst)                     ! Upper boundary mixing ratios [ kg/kg ]
    real(r8) :: ubc_flux(ncnst)                          ! Upper boundary flux [ kg/s/m^2 ]
    real(r8) :: ubc_t(pcols)                             ! Upper boundary temperature [ K ]

    real(r8) :: ws(pcols)                                ! Lowest-level wind speed [ m/s ]
    real(r8) :: tau(pcols)                               ! Turbulent surface stress ( not including mountain stress )
    real(r8) :: ksrfturb(pcols)                          ! Surface drag coefficient of 'normal' stress. > 0.
                                                         ! Virtual mass input per unit time per unit area [ kg/s/m2 ]
    real(r8) :: ksrf(pcols)                              ! Surface drag coefficient of 'normal' stress +
                                                         ! Surface drag coefficient of 'tms' stress.  > 0. [ kg/s/m2 ] 
    real(r8) :: usum_in(pcols)                           ! Vertical integral of input u-momentum. Total zonal
                                                         ! momentum per unit area in column  [ sum of u*dp/g = kg m/s m-2 ]
    real(r8) :: vsum_in(pcols)                           ! Vertical integral of input v-momentum. Total meridional
                                                         ! momentum per unit area in column [ sum of v*dp/g = kg m/s m-2 ]
    real(r8) :: usum_mid(pcols)                          ! Vertical integral of u-momentum after adding explicit residual stress
    real(r8) :: vsum_mid(pcols)                          ! Vertical integral of v-momentum after adding explicit residual stress
    real(r8) :: usum_out(pcols)                          ! Vertical integral of u-momentum after doing implicit diffusion
    real(r8) :: vsum_out(pcols)                          ! Vertical integral of v-momentum after doing implicit diffusion
    real(r8) :: tauimpx(pcols)                           ! Actual net stress added at the current step other than mountain stress
    real(r8) :: tauimpy(pcols)                           ! Actual net stress added at the current step other than mountain stress
    real(r8) :: wsmin                                    ! Minimum sfc wind speed for estimating frictional
                                                         ! transfer velocity ksrf. [ m/s ]
    real(r8) :: ksrfmin                                  ! Minimum surface drag coefficient [ kg/s/m^2 ]
    real(r8) :: timeres                                  ! Relaxation time scale of residual stress ( >= dt ) [ s ]
    real(r8) :: ramda                                    ! dt/timeres [ no unit ]
    real(r8) :: psum
    real(r8) :: u_in, u_res, tauresx_in
    real(r8) :: v_in, v_res, tauresy_in  
    
    real(r8) :: mw_fac_loc(pcols,pver+1,ncnst)           ! Local sqrt(1/M_q + 1/M_d) for this constituent

    !--------------------------------
    ! Variables needed for WACCM-X
    !--------------------------------
    logical  :: kvt_returned
    real(r8) :: ttemp(pcols,pver)			 ! temporary temperature array
    real(r8) :: ttemp0(pcols,pver)			 ! temporary temperature array
    
    if(present(flb) .and. icolprnt(lchnk) > 0)write(202,*)'-->diffsol:dse_1:',dse(icolprnt(lchnk),kprnt)

    dtk(:ncol,:pver) = 0.0_r8 !BSINGH(10/15/2014):Initialized to zero;This variable was uninitialized when we diffuse anything except than 'u' and 'v'
    ! ------------------------------------------------ !
    ! Parameters for implicit surface stress treatment !
    ! ------------------------------------------------ !

    wsmin    = 1._r8                                     ! Minimum wind speed for ksrfturb computation        [ m/s ]
    ksrfmin  = 1.e-4_r8                                  ! Minimum surface drag coefficient                   [ kg/s/m^2 ]
    timeres  = 7200._r8                                  ! Relaxation time scale of residual stress ( >= dt ) [ s ]

    ! ----------------------- !
    ! Main Computation Begins !
    ! ----------------------- !

    errstring = ''
    if( ( diffuse(fieldlist,'u') .or. diffuse(fieldlist,'v') ) .and. .not. diffuse(fieldlist,'s') ) then
          errstring = 'diffusion_solver.compute_vdiff: must diffuse s if diffusing u or v'
          return
    end if
    zero(:) = 0._r8

    ! Compute 'rho' and 'dt*(g*rho)^2/dp' at interfaces

    tint(:ncol,1) = t(:ncol,1)
    rhoi(:ncol,1) = pint(:ncol,1) / (rairi(:ncol,1)*tint(:ncol,1))
    do k = 2, pver
       do i = 1, ncol
          tint(i,k)  = 0.5_r8 * ( t(i,k) + t(i,k-1) )
          rhoi(i,k)  = pint(i,k) / (rairi(i,k)*tint(i,k))
          tmpi2(i,k) = ztodt * ( gravit*rhoi(i,k) )**2 / ( pmid(i,k) - pmid(i,k-1) )
       end do
    end do
    tint(:ncol,pver+1) = t(:ncol,pver)
    rhoi(:ncol,pver+1) = pint(:ncol,pver+1) / ( rair*tint(:ncol,pver+1) )

    rrho(:ncol) = rair  * t(:ncol,pver) / pmid(:ncol,pver)
    tmp1(:ncol) = ztodt * gravit * rpdel(:ncol,pver)

    !--------------------------------------- !
    ! Computation of Molecular Diffusivities !
    !--------------------------------------- !

  ! Modification : Why 'kvq' is not changed by molecular diffusion ? 

    if( do_molec_diff ) then

        if( (.not.present(compute_molec_diff)) .or. (.not.present(vd_lu_qdecomp)) .or. (.not.present(kvt)) ) then
              errstring = 'compute_vdiff: do_molec_diff true but compute_molec_diff or vd_lu_qdecomp or kvt missing'
              return
        endif

      ! The next subroutine 'compute_molec_diff' :
      !     Modifies : kvh, kvm, tint, rhoi, and tmpi2
      !     Returns  : kq_scal, ubc_t, ubc_mmr, dse_top, cc_top, cd_top, cnst_mw, 
      !                cnst_fixed_ubc , mw_fac , ntop_molec 

        mw_fac_loc(:,:,:) = 0._r8

        !--------------------------------------------------------------------------------------------------------
        ! In Extended WACCM, kvt is calculated rather than kvh. This is because molecular diffusion operates on 
        ! temperature, while eddy diffusion operates on dse.  Also, pass in constituent dependent "constants"
        !--------------------------------------------------------------------------------------------------------

        status = compute_molec_diff( lchnk         ,                                                                              &
                                     pcols         , pver            , ncnst      , ncol      , t      , pmid   , pint   ,        &
                                     zi            , ztodt           , kvm        , kvt       , tint   , rhoi   , tmpi2  ,        &
                                     kq_scal       , ubc_t           , ubc_mmr    , ubc_flux  , dse_top, cc_top , cd_top ,cnst_mw,&
                                     cnst_fixed_ubc, cnst_fixed_ubflx, mw_fac_loc , ntop_molec, nbot_molec,       kvt_returned )

    else

        kq_scal(:,:) = 0._r8
        cd_top(:)    = 0._r8
        cc_top(:)    = 0._r8
        kvt_returned = .false.

    endif

    !---------------------------- !
    ! Diffuse Horizontal Momentum !
    !---------------------------- !

    if( diffuse(fieldlist,'u') .or. diffuse(fieldlist,'v') ) then

        ! Compute the vertical upward differences of the input u,v for KE dissipation
        ! at each interface.
        ! Velocity = 0 at surface, so difference at the bottom interface is -u,v(pver)
        ! These 'dinp_u, dinp_v' are computed using the non-diffused input wind.

        do i = 1, ncol
           dinp_u(i,1) = 0._r8
           dinp_v(i,1) = 0._r8
           dinp_u(i,pver+1) = -u(i,pver)
           dinp_v(i,pver+1) = -v(i,pver)
        end do
        do k = 2, pver
           do i = 1, ncol
              dinp_u(i,k) = u(i,k) - u(i,k-1)
              dinp_v(i,k) = v(i,k) - v(i,k-1)
           end do
        end do

       ! -------------------------------------------------------------- !
       ! Do 'Implicit Surface Stress' treatment for numerical stability !
       ! in the lowest model layer.                                     !
       ! -------------------------------------------------------------- !

       if( do_iss ) then

         ! Compute surface drag coefficient for implicit diffusion 
         ! including turbulent turbulent mountain stress. 

           do i = 1, ncol
              ws(i)       = max( sqrt( u(i,pver)**2._r8 + v(i,pver)**2._r8 ), wsmin )
              tau(i)      = sqrt( taux(i)**2._r8 + tauy(i)**2._r8 )
              ksrfturb(i) = max( tau(i) / ws(i), ksrfmin )
           end do              
           ksrf(:ncol) = ksrfturb(:ncol) + ksrftms(:ncol)  ! Do all surface stress ( normal + tms ) implicitly

         ! Vertical integration of input momentum. 
         ! This is total horizontal momentum per unit area [ kg*m/s/m2 ] in each column.
         ! Note (u,v) are the raw input to the PBL scheme, not the
         ! provisionally-marched ones within the iteration loop of the PBL scheme.  

           do i = 1, ncol
              usum_in(i) = 0._r8
              vsum_in(i) = 0._r8
              do k = 1, pver
                 usum_in(i) = usum_in(i) + (1._r8/gravit)*u(i,k)/rpdel(i,k)
                 vsum_in(i) = vsum_in(i) + (1._r8/gravit)*v(i,k)/rpdel(i,k)
              end do
           end do              

         ! Add residual stress of previous time step explicitly into the lowest
         ! model layer with a relaxation time scale of 'timeres'.

           ramda         = ztodt / timeres
           u(:ncol,pver) = u(:ncol,pver) + tmp1(:ncol)*tauresx(:ncol)*ramda
           v(:ncol,pver) = v(:ncol,pver) + tmp1(:ncol)*tauresy(:ncol)*ramda

         ! Vertical integration of momentum after adding explicit residual stress
         ! into the lowest model layer.

           do i = 1, ncol
              usum_mid(i) = 0._r8
              vsum_mid(i) = 0._r8
              do k = 1, pver
                 usum_mid(i) = usum_mid(i) + (1._r8/gravit)*u(i,k)/rpdel(i,k)
                 vsum_mid(i) = vsum_mid(i) + (1._r8/gravit)*v(i,k)/rpdel(i,k)
              end do
           end do              

         ! Debug 
         ! icol = phys_debug_col(lchnk) 
         ! if ( icol > 0 .and. get_nstep() .ge. 1 ) then
         !      tauresx_in = tauresx(icol)
         !      tauresy_in = tauresy(icol)
         !      u_in  = u(icol,pver) - tmp1(icol) * tauresx(icol) * ramda
         !      v_in  = v(icol,pver) - tmp1(icol) * tauresy(icol) * ramda
         !      u_res = u(icol,pver)
         !      v_res = v(icol,pver)
         ! endif
         ! Debug

       else

         ! In this case, do 'turbulent mountain stress' implicitly, 
         ! but do 'normal turbulent stress' explicitly.
         ! In this case, there is no 'redisual stress' as long as 'tms' is
         ! treated in a fully implicit wway, which is true.

         ! 1. Do 'tms' implicitly

           ksrf(:ncol) = ksrftms(:ncol) 

         ! 2. Do 'normal stress' explicitly

           u(:ncol,pver) = u(:ncol,pver) + tmp1(:ncol)*taux(:ncol)
           v(:ncol,pver) = v(:ncol,pver) + tmp1(:ncol)*tauy(:ncol)

       end if  ! End of 'do iss' ( implicit surface stress )

       ! --------------------------------------------------------------------------------------- !
       ! Diffuse horizontal momentum implicitly using tri-diagnonal matrix.                      !
       ! The 'u,v' are input-output: the output 'u,v' are implicitly diffused winds.             !
       !    For implicit 'normal' stress : ksrf = ksrftms + ksrfturb,                            !
       !                                   u(pver) : explicitly include 'redisual normal' stress !
       !    For explicit 'normal' stress : ksrf = ksrftms                                        !
       !                                   u(pver) : explicitly include 'normal' stress          !
       ! Note that in all the two cases above, 'tms' is fully implicitly treated.                !
       ! --------------------------------------------------------------------------------------- !

       call vd_lu_decomp( pcols , pver , ncol  ,                         &
                          ksrf  , kvm  , tmpi2 , rpdel , ztodt , gravit, &
                          zero  , ntop , nbot  , decomp)

       call vd_lu_solve(  pcols , pver  , ncol  ,                        &
                          u     , decomp, ntop  , nbot , zero )

       call vd_lu_solve(  pcols , pver  , ncol  ,                        &
                          v     , decomp, ntop  , nbot , zero )

       ! ---------------------------------------------------------------------- !
       ! Calculate 'total' ( tautotx ) and 'tms' ( tautmsx ) stresses that      !
       ! have been actually added into the atmosphere at the current time step. ! 
       ! Also, update residual stress, if required.                             !
       ! ---------------------------------------------------------------------- !

       do i = 1, ncol

          ! Compute the implicit 'tms' using the updated winds.
          ! Below 'tautmsx(i),tautmsy(i)' are pure implicit mountain stresses
          ! that has been actually added into the atmosphere both for explicit
          ! and implicit approach. 

          tautmsx(i) = -ksrftms(i)*u(i,pver)
          tautmsy(i) = -ksrftms(i)*v(i,pver)

          if( do_iss ) then

            ! Compute vertical integration of final horizontal momentum

              usum_out(i) = 0._r8
              vsum_out(i) = 0._r8
              do k = 1, pver
                 usum_out(i) = usum_out(i) + (1._r8/gravit)*u(i,k)/rpdel(i,k)
                 vsum_out(i) = vsum_out(i) + (1._r8/gravit)*v(i,k)/rpdel(i,k)
              end do

            ! Compute net stress added into the atmosphere at the current time step.
            ! Note that the difference between 'usum_in' and 'usum_out' are induced
            ! by 'explicit residual stress + implicit total stress' for implicit case, while
            ! by 'explicit normal   stress + implicit tms   stress' for explicit case. 
            ! Here, 'tautotx(i)' is net stress added into the air at the current time step.

              tauimpx(i) = ( usum_out(i) - usum_in(i) ) / ztodt
              tauimpy(i) = ( vsum_out(i) - vsum_in(i) ) / ztodt

              tautotx(i) = tauimpx(i) 
              tautoty(i) = tauimpy(i) 

            ! Compute redisual stress and update if required.
            ! Note that the total stress we should have added at the current step is
            ! the sum of 'taux(i) - ksrftms(i)*u(i,pver) + tauresx(i)'.

              if( itaures .eq. 1 ) then
                  tauresx(i) = taux(i) + tautmsx(i) + tauresx(i) - tauimpx(i)
                  tauresy(i) = tauy(i) + tautmsy(i) + tauresy(i) - tauimpy(i)
              endif

          else

              tautotx(i) = tautmsx(i) + taux(i)
              tautoty(i) = tautmsy(i) + tauy(i)
              tauresx(i) = 0._r8
              tauresy(i) = 0._r8

          end if  ! End of 'do_iss' if

       end do ! End of 'do i = 1, ncol' loop

     ! Debug 
     ! icol = phys_debug_col(lchnk) 
     ! if ( icol > 0 .and. get_nstep() .ge. 1 ) then
     !      write(iulog,*)
     !      write(iulog,*)  'diffusion_solver debug'  
     !      write(iulog,*)
     !      write(iulog,*)  'u_in, u_res, u_out'
     !      write(iulog,*)   u_in, u_res, u(icol,pver)
     !      write(iulog,*)  'tauresx_in, tautmsx, tauimpx(actual), tauimpx(derived), tauresx_out, taux'
     !      write(iulog,*)   tauresx_in, tautmsx(icol), tauimpx(icol), -ksrf(icol)*u(icol,pver), tauresx(icol), taux(icol)
     !      write(iulog,*)
     !      write(iulog,*)  'v_in, v_res, v_out'
     !      write(iulog,*)   v_in, v_res, v(icol,pver)
     !      write(iulog,*)  'tauresy_in, tautmsy, tauimpy(actual), tauimpy(derived), tauresy_out, tauy'
     !      write(iulog,*)   tauresy_in, tautmsy(icol), tauimpy(icol), -ksrf(icol)*v(icol,pver), tauresy(icol), tauy(icol)
     !      write(iulog,*)
     !      write(iulog,*)  'itaures, ksrf, ksrfturb, ksrftms'
     !      write(iulog,*)   itaures, ksrf(icol), ksrfturb(icol), ksrftms(icol)
     !      write(iulog,*) 
     ! endif
     ! Debug

       ! ------------------------------------ !
       ! Calculate kinetic energy dissipation !
       ! ------------------------------------ !       

     ! Modification : In future, this should be set exactly same as 
     !                the ones in the convection schemes 

       ! 1. Compute dissipation term at interfaces
       !    Note that 'u,v' are already diffused wind, and 'tautotx,tautoty' are 
       !    implicit stress that has been actually added. On the other hand,
       !    'dinp_u, dinp_v' were computed using non-diffused input wind.

     ! Modification : I should check whether non-consistency between 'u' and 'dinp_u'
     !                is correctly intended approach. I think so.

       k = pver + 1
       do i = 1, ncol
          tmpi1(i,1) = 0._r8
          tmpi1(i,k) = 0.5_r8 * ztodt * gravit * &
                       ( (-u(i,k-1) + dinp_u(i,k))*tautotx(i) + (-v(i,k-1) + dinp_v(i,k))*tautoty(i) )
       end do

       do k = 2, pver
          do i = 1, ncol
             dout_u = u(i,k) - u(i,k-1)
             dout_v = v(i,k) - v(i,k-1)
             tmpi1(i,k) = 0.25_r8 * tmpi2(i,k) * kvm(i,k) * &
                          ( dout_u**2 + dout_v**2 + dout_u*dinp_u(i,k) + dout_v*dinp_v(i,k) )
          end do
       end do

       ! 2. Compute dissipation term at midpoints, add to dry static energy

       do k = 1, pver
          do i = 1, ncol
             dtk(i,k) = ( tmpi1(i,k+1) + tmpi1(i,k) ) * rpdel(i,k)
             if(present(flb) .and. icolprnt(lchnk) ==i .and. k==kprnt)then
                dtk(i,k) = 0.0_r8 !this is WRONG
                write(202,*)'-->diffsol:dse_2:',dse(i,k),'WRONG!!!'
             endif
             dse(i,k) = dse(i,k) + dtk(i,k)
             if(present(flb) .and. icolprnt(lchnk) ==i .and. k==kprnt)write(202,*)'-->diffsol:dse_3:',dse(i,k),dtk(i,k)
          end do
       end do

    end if ! End of diffuse horizontal momentum, diffuse(fieldlist,'u') routine

    !-------------------------- !
    ! Diffuse Dry Static Energy !
    !-------------------------- !

  ! Modification : In future, we should diffuse the fully conservative 
  !                moist static energy,not the dry static energy.

    if( diffuse(fieldlist,'s') ) then

      ! Add counter-gradient to input static energy profiles

        do k = 1, pver
           dse(:ncol,k) = dse(:ncol,k) + ztodt * rpdel(:ncol,k) * gravit  *                &
                                       ( rhoi(:ncol,k+1) * kvh(:ncol,k+1) * cgh(:ncol,k+1) &
                                       - rhoi(:ncol,k  ) * kvh(:ncol,k  ) * cgh(:ncol,k  ) )
           if(present(flb) .and. icolprnt(lchnk) >0 .and. k==kprnt)write(202,*)'-->diffsol:dse_4:',dse(icolprnt(lchnk),k)
       end do

     ! Add the explicit surface fluxes to the lowest layer

       dse(:ncol,pver) = dse(:ncol,pver) + tmp1(:ncol) * shflx(:ncol)
       if(present(flb) .and. icolprnt(lchnk) >0 )write(202,*)'-->diffsol:dse_5:',dse(icolprnt(lchnk),kprnt)

     ! Diffuse dry static energy
     !----------------------------------------------------------------------------------------------------
     ! In Extended WACCM, kvt is calculated rather kvh. This is because molecular diffusion operates on 
     ! temperature, while eddy diffusion operates on dse.  Also, pass in constituent dependent "constants"
     !----------------------------------------------------------------------------------------------------
       if ( kvt_returned ) then 
         cc_top(:) = 0._r8      
         cd_top(:) = 0._r8
       endif

       call vd_lu_decomp( pcols , pver , ncol  ,                         &
                          zero  , kvh  , tmpi2 , rpdel , ztodt , gravit, &
                          cc_top, ntop , nbot  , decomp,flb=flb,lchnk=lchnk )

       call vd_lu_solve(  pcols , pver  , ncol  ,                         &
                          dse   , decomp, ntop  , nbot  , cd_top,flb=flb,lchnk=lchnk )
       if(present(flb) .and. icolprnt(lchnk) >0)write(202,*)'-->diffsol:dse_6:',dse(icolprnt(lchnk),kprnt)
     ! Calculate flux at top interface
     
     ! Modification : Why molecular diffusion does not work for dry static energy in all layers ?

       if( do_molec_diff ) then
           topflx(:ncol) =  - kvh(:ncol,ntop_molec) * tmpi2(:ncol,ntop_molec) / (ztodt*gravit) * &
                            ( dse(:ncol,ntop_molec) - dse_top(:ncol) )
       end if

       !---------------------------------------------------
       ! Solve for temperature using thermal conductivity 
       !---------------------------------------------------
       if ( kvt_returned ) then 

         call vd_lu_decomp( pcols , pver , ncol  ,                         &
                            zero  , kvt  , tmpi2 , rpdel , ztodt , gravit, &
                            zero  , ntop , nbot  , decomp, cpairv )

         do k = 1,pver
           ttemp0(:ncol,k) = t(:ncol,k)
           ttemp(:ncol,k) = ttemp0(:ncol,k)
         enddo

         call vd_lu_solve(                                                &
            pcols, pver, ncol     ,                                       &
            ttemp      ,decomp    ,ntop     ,nbot,zero    )   ! upper boundary is zero flux for extended model
    
         !-------------------------------------
         !  Update dry static energy
         !-------------------------------------
         do k = 1,pver
           dse(:ncol,k) = dse(:ncol,k) + cpairv(:ncol,k)*(ttemp(:ncol,k) - ttemp0(:ncol,k))
           if(present(flb) .and. icolprnt(lchnk) >0 .and. k==kprnt)write(202,*)'-->diffsol:dse_7:',dse(icolprnt(lchnk),k)
         enddo

       endif

    endif

    !---------------------------- !
    ! Diffuse Water Vapor Tracers !
    !---------------------------- !

  ! Modification : For aerosols, I need to use separate treatment 
  !                for aerosol mass and aerosol number. 

    ! Loop through constituents

    need_decomp = .true.

    do m = 1, ncnst

       if( diffuse(fieldlist,'q',m) ) then

           ! Add the nonlocal transport terms to constituents in the PBL.
           ! Check for neg q's in each constituent and put the original vertical
           ! profile back if a neg value is found. A neg value implies that the
           ! quasi-equilibrium conditions assumed for the countergradient term are
           ! strongly violated.

           qtm(:ncol,:pver) = q(:ncol,:pver,m)

           do k = 1, pver
              q(:ncol,k,m) = q(:ncol,k,m) + &
                             ztodt * rpdel(:ncol,k) * gravit  * ( cflx(:ncol,m) * rrho(:ncol) ) * &
                           ( rhoi(:ncol,k+1) * kvh(:ncol,k+1) * cgs(:ncol,k+1)                    &
                           - rhoi(:ncol,k  ) * kvh(:ncol,k  ) * cgs(:ncol,k  ) )
           end do
           lqtst(:ncol) = all(q(:ncol,1:pver,m) >= qmincg(m), 2)
           do k = 1, pver
              q(:ncol,k,m) = merge( q(:ncol,k,m), qtm(:ncol,k), lqtst(:ncol) )
           end do

           ! Add the explicit surface fluxes to the lowest layer

           q(:ncol,pver,m) = q(:ncol,pver,m) + tmp1(:ncol) * cflx(:ncol,m)

           ! Diffuse constituents.

           if( need_decomp ) then

               call vd_lu_decomp( pcols , pver , ncol  ,                         &
                                  zero  , kvq  , tmpi2 , rpdel , ztodt , gravit, &
                                  zero  , ntop , nbot  , decomp)

               if( do_molec_diff ) then

                 ! This is for solving molecular diffusion of minor species, thus, for WACCM-X, bypass O and O2 (major species)
                 ! Major species diffusion is calculated separately.  -Hanli Liu
                 !
                   if ( diffuse(fieldlistm,'q',m) ) then                  ! decide if diffuse this constituent

                     status = vd_lu_qdecomp( pcols , pver   , ncol              , cnst_fixed_ubc(m), cnst_mw(m), ubc_mmr(:,m), &
                                             kvq   , kq_scal, mw_fac_loc(:,:,m) ,tmpi2             , rpdel     ,               &
                                             decomp, rhoi      ,               &
                                             tint  , ztodt  , ntop_molec        , nbot_molec       , cd_top    ,               &
                                             lchnk , pmid   , pint              , t                , m         )


                     ! This to calculate the upper boundary flux of H.    -Hanli Liu
                     if ((cnst_fixed_ubflx(m))) &
                            q(:ncol,1,m) = q(:ncol,1,m) - ztodt * rpdel(:ncol,1) * gravit * ubc_flux(m)

                   endif

               else
                   need_decomp =  .false.
               endif
           end if

           call vd_lu_solve(  pcols , pver , ncol  ,                         &
                              q(1,1,m) , decomp    , ntop  , nbot  , cd_top )
       end if
    end do

    ! Probably unnecessary, but just in case there's an OpenMP problem,
    ! deallocate the decomp object.
    call decomp%finalize()

    return
  end subroutine compute_vdiff

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !
  
  character(128) function vdiff_select( fieldlist, name, qindex )
    ! --------------------------------------------------------------------- !
    ! This function sets the field with incoming name as one to be diffused !
    ! --------------------------------------------------------------------- !
    type(vdiff_selector), intent(inout)        :: fieldlist
    character(*),         intent(in)           :: name
    integer,              intent(in), optional :: qindex
    
    vdiff_select = ''
    select case (name)
    case ('u','U')
       fieldlist%fields(1) = .true.
    case ('v','V')
       fieldlist%fields(2) = .true.
    case ('s','S')
       fieldlist%fields(3) = .true.
    case ('q','Q')
       if( present(qindex) ) then
           fieldlist%fields(3 + qindex) = .true.
       else
           fieldlist%fields(4) = .true.
       endif
    case default
       write(vdiff_select,*) 'Bad argument to vdiff_index: ', name
    end select
    return
    
  end function vdiff_select

  type(vdiff_selector) function not(a)
    ! ------------------------------------------------------------- !
    ! This function extends .not. to operate on type vdiff_selector !
    ! ------------------------------------------------------------- !    
    type(vdiff_selector), intent(in)  :: a
    allocate(not%fields(size(a%fields)))
    not%fields = .not. a%fields
  end function not

  logical function my_any(a)
    ! -------------------------------------------------- !
    ! This function extends the intrinsic function 'any' ! 
    ! to operate on type vdiff_selector                  ! 
    ! -------------------------------------------------- !
    type(vdiff_selector), intent(in) :: a
    my_any = any(a%fields)
  end function my_any

  logical function diffuse(fieldlist,name,qindex)
    ! ---------------------------------------------------------------------------- !
    ! This function reports whether the field with incoming name is to be diffused !
    ! ---------------------------------------------------------------------------- !
    type(vdiff_selector), intent(in)           :: fieldlist
    character(*),         intent(in)           :: name
    integer,              intent(in), optional :: qindex
    
    select case (name)
    case ('u','U')
       diffuse = fieldlist%fields(1)
    case ('v','V')
       diffuse = fieldlist%fields(2)
    case ('s','S')
       diffuse = fieldlist%fields(3)
    case ('q','Q')
       if( present(qindex) ) then
           diffuse = fieldlist%fields(3 + qindex)
       else
           diffuse = fieldlist%fields(4)
       endif
    case default
       diffuse = .false.
    end select
    return
  end function diffuse

end module diffusion_solver
