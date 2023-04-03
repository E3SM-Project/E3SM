!===============================================================================
! CLDERA stratospheric aerosol injection tracers
! provides dissipation rates, heatings rates, and injection forcings for diagnostic constituents
!
! Joe Hollowed
! June 2022
! This module written based on the structure of aoa_tracers, and enables the advection and
! evolution of 3 tracers constituents:
! - SO2
! - ASH
! - SULFATE
!===============================================================================

module cldera_sai_tracers

  use shr_kind_mod, only: r8 => shr_kind_r8
  use spmd_utils,   only: masterproc
  use ppgrid,       only: pcols, pver
  use constituents, only: pcnst, cnst_add, cnst_name, cnst_longname
  use cam_logfile,  only: iulog
  use physconst,    only: pi

  implicit none
  private
  save

  ! Public interfaces
  public :: cldera_sai_tracers_register         ! register constituents
  public :: cldera_sai_tracers_implements_cnst  ! true if constituent is implemented by this package
  public :: cldera_sai_tracers_init_cnst        ! initialize constituent field
  public :: cldera_sai_tracers_init             ! initialize history fields, datasets
  public :: cldera_sai_tracers_timestep_tend    ! calculate tendencies
  public :: cldera_sai_tracers_readnl           ! read namelist options

  !----  Private module data

  integer, parameter :: ncnst=3  ! number of constituents implemented by this module

  ! constituent names, indices
  character(len=7), parameter :: c_names(ncnst) = (/'SO2    ', 'ASH    ', 'SULFATE'/)
  integer :: ifirst ! global index of first constituent
  integer :: ixso2  ! global index for SO2 tracer
  integer :: ixash  ! global index for ASH tracer
  integer :: ixsulf ! global index for SULFATE tracer

  ! Data from namelist variables; defaults set in bld/namelist_files/namelist_defaults_eam.xml
  logical  :: cldera_sai_tracers_flag      ! true => activate module, set namelist variable
  logical  :: cldera_sai_read_from_ic_file ! true => tracers initialized from IC file
  real(r8) :: cldera_sai_lat0              ! desired lat of injection (deg)
  real(r8) :: cldera_sai_lon0              ! desired lon of injection (deg)
  real(r8) :: cldera_sai_z0                ! peak of initial injection vertical distribution (km)
  real(r8) :: cldera_sai_MSO2              ! total SO2 mass (Mt)
  real(r8) :: cldera_sai_Mash              ! total ash mass (Mt)
  real(r8) :: cldera_sai_w                 ! SO2->sulfate reaction mass weighting (dimensionless)
  real(r8) :: cldera_sai_duration          ! injection duration (hours)
  real(r8) :: cldera_sai_t0                ! time of injection after simulation start (days)
  real(r8) :: cldera_sai_rkSO2             ! SO2 e-folding (1/day)
  real(r8) :: cldera_sai_rkash             ! ash e-folding (1/day)
  real(r8) :: cldera_sai_rksulf            ! ash e-folding (1/day)
  real(r8) :: cldera_sai_z_cool            ! max height to apply surface cooling by AOD (m)
  real(r8) :: cldera_sai_zeta              ! surface cooling efficiency (dimensionless)
  real(r8) :: cldera_sai_bsw_so2           ! shortwave mass extinction coeff. for SO2 (m2/kg)
  real(r8) :: cldera_sai_bsw_sulf          ! shortwave mass extinction coeff. for sulfate (m2/kg)
  real(r8) :: cldera_sai_bsw_ash           ! shortwave mass extinction coeff. for ash (m2/kg)
  real(r8) :: cldera_sai_blw_so2           ! longwave mass extinction coeff. for SO2 (m2/kg)
  real(r8) :: cldera_sai_blw_sulf          ! longwave mass extinction coeff. for sulfate (m2/kg)
  real(r8) :: cldera_sai_blw_ash           ! longwave mass extinction coeff. for ash (m2/kg)
  logical  :: cldera_sai_formSulfate       ! true => activate sulfate formation
  logical  :: cldera_sai_stratHeating      ! true => activate strat. heating
  logical  :: cldera_sai_surfCooling       ! true => activate surface cooling
    
  ! --- parameters

  ! constants
  real(r8), parameter :: deg2rad = pi/180._r8
  real(r8), parameter :: rad2deg = 180._r8/pi
  real(r8), parameter :: day2s   = 86400._r8
  real(r8), parameter :: s2day   = 1.0_r8/day2s
  real(r8), parameter :: hr2s    = 3600._r8
  real(r8), parameter :: s2hr    = 1.0_r8/hr2s
  real(r8), parameter :: Mt2kg   = 1.0e9_r8
  real(r8), parameter :: km2m    = 1000._r8
  real(r8), parameter :: hPa2Pa  = 100._r8
  real(r8), parameter :: P0      = 1000 * hPa2Pa    ! reference presssure (Pa)
  
  ! for injection source localization
  real(r8) :: lat0                      ! lat of injection (rad)
  real(r8) :: lon0                      ! lon of injection (rad)           
  real(r8) :: dur                       ! injection duration (s)
  real(r8) :: t0                        ! initial time of injection (s)
  real(r8) :: tf                        ! final time of injection (s)
  real(r8) :: ttol  = 60.0_r8           ! time tolerance (s)
  real(r8) :: z0                        ! peak of injection vertical distribution (m) 
  real(r8) :: sigma = 1.5_r8 * km2m      ! vertical width parameter (m)
  
  ! for injected species 
  real(r8) :: M_so2                     ! total SO2 mass (kg)
  real(r8) :: M_ash                     ! total ash mass (kg)
  real(r8) :: k_so2                     ! SO2 e-folding (1/s)
  real(r8) :: k_ash                     ! ash e-folding (1/s)
  real(r8) :: k_sulf                    ! sulfate e-folding (1/s)
  real(r8) :: w                         ! SO2->sulfate reaction mass weighting
  
  ! for mass estimate
  real(r8) :: col_mass_check            ! total column mass (kg)
  
  ! for stratospheric heating
  real(r8) :: blw_so2                   ! longwave mass extinction coefficient for SO2 (m2/kg)
  real(r8) :: blw_sulf                  ! longwave mass extinction coefficient for sulfate (m2/kg)
  real(r8) :: blw_ash                   ! longwave mass extinction coefficient for ash (m2/kg)
  
  ! for surface cooling
  real(r8) :: zeta                      ! surface cooling efficiency (dimensionless)
  real(r8) :: z_cool                    ! max height to apply surf. cooling (m)
  real(r8) :: bsw_so2                   ! shortwave mass extinction coefficient for SO2 (m2/kg)
  real(r8) :: bsw_sulf                  ! shortwave mass extinction coefficient for sulfate (m2/kg)
  real(r8) :: bsw_ash                   ! shortwave mass extinction coefficient for ash (m2/kg)
  

!===============================================================================
!===============================================================================

contains

!================================================================================
!================================================================================


  subroutine cldera_sai_tracers_readnl(nlfile)

    use namelist_utils,     only: find_group_name
    use units,              only: getunit, freeunit
    use mpishorthand
    use cam_abortutils,     only: endrun

    implicit none

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'cldera_sai_tracers_readnl'

    ! call READ with cldera_sai_tracers_nl as the 'namelist group'; names in the namelist 
    ! will be matched to values of those same names in the file connected to unitn
    ! see "Input Actions", "Data Syntax" at the following page:
    ! https://docs.oracle.com/cd/E19957-01/805-4939/6j4m0vnc6/index.html
    ! The default values must be constants, and cannot contain arithemtic operators
    namelist /cldera_sai_tracers_nl/ cldera_sai_tracers_flag, cldera_sai_read_from_ic_file, &
                                     cldera_sai_lat0, cldera_sai_lon0, cldera_sai_z0, &
                                     cldera_sai_MSO2, cldera_sai_Mash, cldera_sai_w, &
                                     cldera_sai_duration, cldera_sai_t0, &
                                     cldera_sai_rkSO2, cldera_sai_rkash, cldera_sai_rksulf, &
                                     cldera_sai_z_cool, cldera_sai_zeta, &
                                     cldera_sai_bsw_so2, cldera_sai_blw_so2, & 
                                     cldera_sai_bsw_sulf, cldera_sai_blw_sulf, & 
                                     cldera_sai_bsw_ash, cldera_sai_blw_ash, & 
                                     cldera_sai_formSulfate, cldera_sai_stratHeating, &
                                     cldera_sai_surfCooling

    !-----------------------------------------------------------------------------

    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'cldera_sai_tracers_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, cldera_sai_tracers_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    
       ! set local values, scale to SI units
       lat0   = cldera_sai_lat0     * deg2rad     
       lon0   = cldera_sai_lon0     * deg2rad     
       z0     = cldera_sai_z0       * km2m    
       M_so2  = cldera_sai_MSO2     * Mt2kg      
       M_ash  = cldera_sai_Mash     * Mt2kg      
       w      = cldera_sai_w                 
       dur    = cldera_sai_duration * hr2s
       t0     = cldera_sai_t0       * day2s          
       tf     = t0 + dur
       
       k_so2  = 1.0_r8 / (cldera_sai_rkSO2  * day2s)
       k_ash  = 1.0_r8 / (cldera_sai_rkash  * day2s)     
       k_sulf = 1.0_r8 / (cldera_sai_rksulf * day2s)     
       
       z_cool   = cldera_sai_z_cool
       zeta     = cldera_sai_zeta
       bsw_so2  = cldera_sai_bsw_so2
       bsw_sulf = cldera_sai_bsw_sulf
       bsw_ash  = cldera_sai_bsw_ash
       blw_so2  = cldera_sai_blw_so2
       blw_sulf = cldera_sai_blw_sulf
       blw_ash  = cldera_sai_blw_ash
       
    end if

#ifdef SPMD
    ! Broadcast namelist variables
    call mpibcast(cldera_sai_tracers_flag,      1, mpilog, 0, mpicom)
    call mpibcast(cldera_sai_read_from_ic_file, 1, mpilog, 0, mpicom)
    call mpibcast(cldera_sai_formSulfate,       1, mpilog, 0, mpicom)
    call mpibcast(cldera_sai_stratHeating,      1, mpilog, 0, mpicom)
    call mpibcast(cldera_sai_surfCooling,       1, mpilog, 0, mpicom)
    call mpibcast(lat0,          1, mpir8, 0, mpicom)
    call mpibcast(lon0,          1, mpir8, 0, mpicom)
    call mpibcast(z0,            1, mpir8, 0, mpicom)
    call mpibcast(M_so2,         1, mpir8, 0, mpicom)
    call mpibcast(M_ash,         1, mpir8, 0, mpicom)
    call mpibcast(w,             1, mpir8, 0, mpicom)
    call mpibcast(dur,           1, mpir8, 0, mpicom)
    call mpibcast(t0,            1, mpir8, 0, mpicom)
    call mpibcast(tf,            1, mpir8, 0, mpicom)
    call mpibcast(k_so2,         1, mpir8, 0, mpicom)
    call mpibcast(k_ash,         1, mpir8, 0, mpicom)
    call mpibcast(k_sulf,        1, mpir8, 0, mpicom)
    call mpibcast(z_cool,        1, mpir8, 0, mpicom)
    call mpibcast(zeta,          1, mpir8, 0, mpicom)
    call mpibcast(bsw_so2,       1, mpir8, 0, mpicom)
    call mpibcast(bsw_sulf,      1, mpir8, 0, mpicom)
    call mpibcast(bsw_ash,       1, mpir8, 0, mpicom)
    call mpibcast(blw_so2,       1, mpir8, 0, mpicom)
    call mpibcast(blw_sulf,      1, mpir8, 0, mpicom)
    call mpibcast(blw_ash,       1, mpir8, 0, mpicom)
#endif

  endsubroutine cldera_sai_tracers_readnl


!================================================================================
!================================================================================


  subroutine cldera_sai_tracers_register
    !-----------------------------------------------------------------------
    !
    ! Purpose: register advected constituents
    !
    !-----------------------------------------------------------------------
    use physconst,  only: cpair, mwdry
    !-----------------------------------------------------------------------

    if (masterproc) then
      write(iulog,*) 'CLDERA SAI: REGISTERING CONSTITUENTS'
    endif
    if (.not. cldera_sai_tracers_flag) return

    call cnst_add(c_names(1), mwdry, cpair, 0._r8, ixso2, readiv=cldera_sai_read_from_ic_file, &
                  mixtype='dry', longname='Stratospheric aerosol injection plume SO2')
    ifirst = ixso2
    call cnst_add(c_names(2), mwdry, cpair, 0._r8, ixash, readiv=cldera_sai_read_from_ic_file, &
                  mixtype='dry', longname='Stratospheric aerosol injection plume ash')
    call cnst_add(c_names(3), mwdry, cpair, 0._r8, ixsulf, readiv=cldera_sai_read_from_ic_file, &
                  mixtype='dry', longname='Stratospheric aerosol injection plume sulfate')

  end subroutine cldera_sai_tracers_register


!===============================================================================
!===============================================================================


  function cldera_sai_tracers_implements_cnst(name)
    !-----------------------------------------------------------------------
    !
    ! Purpose: return true if specified constituent is implemented by this package
    !
    !-----------------------------------------------------------------------

    character(len=*), intent(in) :: name   ! constituent name
    logical :: cldera_sai_tracers_implements_cnst        ! return value

    !---------------------------Local workspace-----------------------------
    integer :: m
    !-----------------------------------------------------------------------

    cldera_sai_tracers_implements_cnst = .false.

    if (.not. cldera_sai_tracers_flag) return

    do m = 1, ncnst
       if (name == c_names(m)) then
          cldera_sai_tracers_implements_cnst = .true.
          return
       end if
    end do

  end function cldera_sai_tracers_implements_cnst


!===============================================================================
!===============================================================================


  subroutine cldera_sai_tracers_init_cnst(name, q, gcid)

    !-----------------------------------------------------------------------
    !
    ! Purpose: initialize test tracers mixing ratio fields
    ! This subroutine is called at the beginning of an initial run ONLY
    !
    !-----------------------------------------------------------------------

    character(len=*), intent(in)  :: name
    real(r8),         intent(out) :: q(:,:)   ! kg tracer/kg dry air (gcol, plev)
    integer,          intent(in)  :: gcid(:)  ! global column id

    integer :: m
    !-----------------------------------------------------------------------

    if (.not. cldera_sai_tracers_flag) return

    do m = 1, ncnst
       if (name ==  c_names(m))  then
          ! pass global constituent index
          call init_cnst_3d(ifirst+m-1, q, gcid)
       endif
    end do

  end subroutine cldera_sai_tracers_init_cnst


!===============================================================================
!===============================================================================


  subroutine cldera_sai_tracers_init

    !-----------------------------------------------------------------------
    !
    ! Purpose: initialize SAI constituents
    !          (declare history variables)
    !-----------------------------------------------------------------------

    use cam_history,    only: addfld, add_default, horiz_only

    integer :: m, mm
    !-----------------------------------------------------------------------

    if (.not. cldera_sai_tracers_flag) return

    ! Set names of tendencies and declare them as history variables
    do m = 1, ncnst
       mm = ifirst+m-1
       call addfld(cnst_name(mm), (/ 'lev' /), 'A', 'kg/kg', cnst_longname(mm))   
       call add_default (cnst_name(mm), 1, ' ')
    end do
    
    ! add other output fields 
    call addfld('AREA',  horiz_only, 'A', 'm2', 'area of grid cell' )
    call add_default('AREA', 1, ' ')
    call addfld('AIR_MASS',  (/ 'lev' /), 'A', 'kg', 'mass of grid box' )
    call add_default('AIR_MASS', 1, ' ')
    call addfld('COL_MASS',  horiz_only, 'A', 'kg', 'total air mass of column' )
    call add_default('COL_MASS', 1, ' ')
    call addfld('SAI_HEAT',  (/ 'lev' /), 'A', 'K/day', 'radiative heating rate' )
    call add_default('SAI_HEAT', 1, ' ')
    ! strat. heating quantities
    call addfld('I_LW',  horiz_only, 'A', 'W/m2', 'incident longwave flux density' )
    call add_default('I_LW', 1, ' ')
    call addfld('ATTEN_LW',  (/ 'lev' /), 'A', 'dimensionless', 'longwave extinction' )
    call add_default('ATTEN_LW', 1, ' ')
    ! surface cooling quantities
    call addfld('I_SW',  horiz_only, 'A', 'W/m2', 'incident shortwave flux density' )
    call add_default('I_SW', 1, ' ')
    call addfld('AOD',  horiz_only, 'A', 'dimensionless', 'total AOD from all constituents' )
    call add_default('AOD', 1, ' ')
    call addfld('COL_SO2',  horiz_only, 'A', 'kg', 'total SO2 mass in column' )
    call add_default('COL_SO2', 1, ' ')
    call addfld('COL_SULF',  horiz_only, 'A', 'kg', 'total sulfate mass in column' )
    call add_default('COL_SULF', 1, ' ')
    call addfld('COL_ASH',  horiz_only, 'A', 'kg', 'total ash mass in column' )
    call add_default('COL_ASH', 1, ' ')
    call addfld('COL_SAI',  horiz_only, 'A', 'kg', 'total of all tracer masses in column' )
    call add_default('COL_SAI', 1, ' ')
    call addfld('ATTEN_SW',  horiz_only, 'A', 'dimensionless', 'shortwave extinction' )
    call add_default('ATTEN_SW', 1, ' ')


  end subroutine cldera_sai_tracers_init


!===============================================================================
!===============================================================================


  subroutine cldera_sai_tracers_timestep_tend(state, ptend, dt, ncol)

    use physics_types,   only: physics_state, physics_ptend, physics_ptend_init
    use phys_grid,       only: get_area_p, get_area_all_p, phys_grid_find_col  
    use physconst,       only: rearth, cpair, rair, rgrav=>rga, rearth, sb=>stebol
    use time_manager,    only: get_curr_time
    use time_manager,    only: get_nstep
    use cam_history,     only: outfld
    use cam_abortutils,  only: endrun 
    use ref_pres,        only: ptop_ref, psurf_ref
    use mpishorthand
#if defined(CLDERA_PROFILING)
    use ppgrid,         only: begchunk
    use cldera_interface_mod, only: cldera_set_field_part_data
#endif

    ! Arguments
    type(physics_state), intent(inout) :: state              ! state variables
    type(physics_ptend), intent(out)   :: ptend              ! package tendencies
    real(r8),            intent(in)    :: dt                 ! timestep
    integer,             intent(in)    :: ncol               ! no. of column in chunk

    !----------------- Local workspace -------------------------------

    integer  :: i, k
    integer  :: lchnk             ! chunk identifier
    integer  :: nstep             ! current timestep number
    logical  :: lq(pcnst)
    
    ! state variables
    integer  :: day,sec
    real(r8) :: t                 ! time (s)           
    real(r8) :: lat               ! latitude (rad)  
    real(r8) :: lon               ! longitude (rad)  
    real(r8) :: zz                ! midpoint geopotential height above surface (m)
    real(r8) :: pp                ! midpoint pressure (Pa)
    
    ! for injection source localization
    real(r8) :: V(pver)                  ! vertical profile (1/m)
    real(r8) :: sVk = 0.0                ! vertical profile sum (1/m)
    logical  :: inject                   ! true => this rank contains the injection column
    integer  :: inj_owner                ! MPI rank owning injection column
    integer  :: inji                     ! injection column index within chunk
    integer  :: injc                     ! injection local chunk index
    integer  :: myrank                   ! rank of this MPI process
    integer  :: ier                      ! return error status
    real(r8) :: source_on                ! switch for source term

    ! for tracer mass normalization
    real(r8) :: A_so2 = 0.0              ! SO2 normalization (kg m / s)
    real(r8) :: A_ash = 0.0              ! ash normalization (kg m / s)
    
    ! for mass estimates
    real(r8) :: area(ncol)               ! column horizontal area (m^2)
    real(r8) :: grid_mass(ncol,pver)     ! grid cell air mass (kg)
    real(r8) :: col_mass(ncol)           ! total column air mass (kg)
    real(r8) :: col_mass_so2(ncol)       ! total column SO2 mass (kg)
    real(r8) :: col_mass_sulf(ncol)      ! total column sulfate mass (kg)
    real(r8) :: col_mass_ash(ncol)       ! total column ash mass (kg)
    
    ! for stratospheric heating
    real(r8) :: Ilw(ncol)                ! longwave incident flux density (W/m^2)
    real(r8) :: cell_od = 0.0            ! optical depth of a single cell (dimensionless)
    real(r8) :: cell_od_sum = 0.0        ! total optical depth below a cell (dimensionless)
    real(r8) :: dI_lw(ncol,pver)         ! longwave attenuation (W/m2)
    real(r8) :: local_heating(ncol,pver) ! stratospheric heating (J/kg/s)
    
    ! for surface cooling
    real(r8) :: Isw(ncol)               ! shortwave incident flux density (W/m^2)
    real(r8) :: aod(ncol)               ! total AOD (dimensionless)
    real(r8) :: aod_so2(ncol)           ! SO2 AOD (dimensionless)
    real(r8) :: aod_ash(ncol)           ! ash AOD (dimensionless)
    real(r8) :: aod_sulf(ncol)          ! sulfate AOD (dimensionless)
    real(r8) :: mass_cool(ncol)         ! air mass within z_cool os surface (kg)
    real(r8) :: dI_sw(ncol)             ! shortwave attenuation (W/m2)
    real(r8) :: surf_cooling(ncol,pver) ! surface cooling (J/kg/s)
    
    V = 0.0
    dI_sw = 0.0
    dI_lw = 0.0
    local_heating = 0.0
    surf_cooling = 0.0
    mass_cool = 0


    !------------------------------------------------------------------

    if (.not. cldera_sai_tracers_flag) then
       !Initialize an empty ptend for use with physics_update
       call physics_ptend_init(ptend,state%psetcols,'none')
       return
    end if

    lq(:)       = .FALSE.
    lq(ixso2)   = .TRUE.
    lq(ixash)   = .TRUE.
    lq(ixsulf)  = .TRUE.
    call physics_ptend_init(ptend, state%psetcols, 'cldera_sai_tracers', lq=lq, ls=.true.)

    nstep = get_nstep()
    lchnk = state%lchnk
    
    call get_curr_time(day,sec)
    t = (day * 86400.0) + sec   ! current time in seconds
    
    !------------------------------------------------------------------
   

    ! =============== LOCATE INJECTION COLUMN ===============
    ! uses phys_grid_find_col from cam/physics/phys_grid.F90
    ! this locates the process containing the column nearest the requested
    ! injection site. 
    ! For processes where myrank /= inj_owner, indexing with the returned injc, inji
    ! may cause a segfault, don't do it
    call phys_grid_find_col(lat0 * rad2deg, lon0 * rad2deg, inj_owner, injc, inji)
    call mpi_comm_rank (mpicom, myrank, ier)
    ! turn on 'inject' bool if this process contains the injection column, 
    ! and if the injection time has not elapsed
    inject = myrank == inj_owner .and. lchnk == injc .and. &
             t > (t0-ttol) .and. t < (tf-ttol)
   
    ! =============== COMPUTE VERTICAL PROFILE & NORMALIZATION ===============
    ! Compute V(k) with input z(inji, k) for myrank == inj_owner
    ! For myrank /= inj_owner, V(:), A_so2, A_ash, sVk = 0
    V(:) = 0.0
    if(inject) then
        do k = 1,pver
            zz  = state%zm(inji, k)
            V(k) = EXP(-(1.0/2.0) * (zz - z0)**2/sigma**2)
        enddo
        sVk = SUM(V(:))
        A_so2 = M_so2 / (dur * sVk) 
        A_ash = M_ash / (dur * sVk)
    endif 

    ! =============== COMPUTE GRID AND COLUMN MASS, AOD, LW AND SW PROFILES  =============== 
    call get_area_all_p(lchnk, ncol, area)
    area = area * rearth**2  ! rad^2 to m^2
    
    do i = 1, ncol
        
        lat = state%lat(i)
        ! shortwave and longwave profiles
        Isw(i) = sb*(315 - 60*SIN(lat)*SIN(lat))**4
        Ilw(i) = 558.54*COS(lat)

        do k = 1,pver
            ! mass in grid cell via hydrostatic approximation
            ! pdel = Pa, area = m^2, rgrav = s**2/m ===> mass = kg
            grid_mass(i,k) = state%pdel(i,k) * area(i) * rgrav
            
            ! sum column mass above the surface for air, tracers
            if(k == 1) then ! model top
                col_mass(i)      = grid_mass(i, k)
                col_mass_so2(i)  = grid_mass(i, k) * state%q(i,k,ixso2)
                col_mass_sulf(i) = grid_mass(i, k) * state%q(i,k,ixsulf)
                col_mass_ash(i)  = grid_mass(i, k) * state%q(i,k,ixash)
            else
                col_mass(i)      = col_mass(i) + grid_mass(i, k)
                col_mass_so2(i)  = col_mass_so2(i)  + grid_mass(i, k) * state%q(i,k,ixso2)
                col_mass_sulf(i) = col_mass_sulf(i) + grid_mass(i, k) * state%q(i,k,ixsulf)
                col_mass_ash(i)  = col_mass_ash(i)  + grid_mass(i, k) * state%q(i,k,ixash)
            endif
            
            ! ---- get total air mass below z_cool in column
            if(state%zm(i,k) <= z_cool) then
                mass_cool(i) = mass_cool(i) + grid_mass(i,k)
            end if

        enddo

        ! compute AOD; mass extinction coefficients times cumulative sum of tracer mass, 
        ! normalized by column area 
        aod_so2(i)  = bsw_so2  * col_mass_so2(i)  / area(i)
        aod_sulf(i) = bsw_sulf * col_mass_sulf(i) / area(i)
        aod_ash(i)  = bsw_ash  * col_mass_ash(i)  / area(i)
        aod(i) = aod_so2(i) + aod_sulf(i) + aod_ash(i)

        ! check consistency of cumulative mass calculation against reference pressures at
        ! model top, surface; raise warning if disagrees more than 5%
        !col_mass_check = (psurf_ref - ptop_ref) * area(i) * rgrav
        !if ( ABS(col_mass(i) - col_mass_check) / col_mass_check  > 0.05_r8) then
        !     write(iulog,*) "CLDERA_SAI_TRACERS cumulative mass error at col ", i
        !     write(iulog,*) "   expected:               ", col_mass_check
        !     write(iulog,*) "   measured cumulative:    ", col_mass(i)
        !     write(iulog,*) "   relative err:           ",   &
        !                    100.0 * ABS(col_mass(i) - col_mass_check)/col_mass_check, "%"
        !     ! uncomment to end model run if this check fails
        !     !call endrun()
        !end if
    enddo

    ! record air mass on grid to history files
    call outfld('AOD',      aod(:), ncol, lchnk)
    call outfld('I_SW',     Isw(:), ncol, lchnk)
    call outfld('I_LW',     Ilw(:), ncol, lchnk)
    call outfld('AREA',     area(:), ncol, lchnk)
    call outfld('AIR_MASS', grid_mass(:,:), ncol, lchnk)
    call outfld('COL_MASS', col_mass(:), ncol, lchnk)
    call outfld('COL_SO2',  col_mass_so2(:), ncol, lchnk)
    call outfld('COL_SULF', col_mass_sulf(:), ncol, lchnk)
    call outfld('COL_ASH',  col_mass_ash(:), ncol, lchnk)
    call outfld('COL_SAI',  col_mass_ash(:)+col_mass_so2(:)+col_mass_sulf(:), ncol, lchnk)

#if defined(CLDERA_PROFILING)
    call cldera_set_field_part_data("aod_so2" ,lchnk-begchunk+1,aod_so2)
    call cldera_set_field_part_data("aod_ash" ,lchnk-begchunk+1,aod_ash)
    call cldera_set_field_part_data("aod_sulf",lchnk-begchunk+1,aod_sulf)
    call cldera_set_field_part_data("aod"     ,lchnk-begchunk+1,aod)
#endif

    ! =============== COMPUTE TENDENCIES ===============
    do i = 1, ncol
        
        ! check that this column of this process is the injection column
        source_on = 0.0
        if(inject .and. i==inji) source_on = 1.0
        
        do k = pver, 1, -1
          ! index k decrements from pver to 1, 
          ! so that cell_od_sum accumulates optical depth starting at the surface
          
          lat = state%lat(i)
          lon = state%lon(i)

          ! =============== SO2 ===============
          ! ---- source + decay
          ptend%q(i,k,ixso2) = 1/grid_mass(i, k) * &
                               (-k_so2 * state%q(i, k, ixso2) * grid_mass(i, k) + &
                                A_so2 * V(k) * source_on)

          ! =============== ASH ===============
          ! ---- source + decay
          ptend%q(i,k,ixash) = 1/grid_mass(i, k) * &
                               (-k_ash * state%q(i, k, ixash) * grid_mass(i, k) + &
                                A_ash * V(k) * source_on)

          ! =============== SULFATE ===============
          ptend%q(i,k,ixsulf) = 0.0
          if(cldera_sai_formSulfate) then 
              ! ---- source + decay
              ptend%q(i,k,ixsulf) = -k_sulf * state%q(i, k, ixsulf) + &
                                    w * k_so2 * state%q(i, k, ixso2)
          endif
          
          ! =============== LOCAL HEATING ===============
          local_heating(i, k) = 0.0
          if(cldera_sai_stratHeating) then
              
              ! compute optical depth of this grid cell
              cell_od = (grid_mass(i,k)/area(i)) * (blw_so2 * state%q(i,k,ixso2) + &
                                                    blw_sulf * state%q(i,k,ixsulf) + &
                                                    blw_ash * state%q(i,k,ixash)) 
              ! update cumulative optical depth of column beneath grid cell; 
              ! reset to at 0 for each new column
              if(k == pver) then 
                  cell_od_sum = 0.0
              else 
                  cell_od_sum = cell_od_sum + (grid_mass(i,k+1)/area(i)) * &
                                               (blw_so2 * state%q(i,k+1,ixso2) + &
                                                blw_sulf * state%q(i,k+1,ixsulf) + &
                                                blw_ash * state%q(i,k+1,ixash))
              endif
    
              dI_lw(i,k) = Ilw(i) * EXP(-cell_od_sum) * (1.0 - EXP(-cell_od))
              local_heating(i, k) = (area(i) / grid_mass(i,k)) * dI_lw(i,k)
       
          endif
    
          ! =============== SURFACE COOLING ===============
          surf_cooling(i, k) = 0.0
          if(cldera_sai_surfCooling .and. state%zm(i,k) <= z_cool) then
              
              dI_sw(i) = Isw(i) * (EXP(-aod(i)) - 1.0)
              surf_cooling(i, k) = (zeta * area(i) / mass_cool(i)) * dI_sw(i)
          
          endif

          ! --- total heating in J/kg/s
          ptend%s(i, k) = local_heating(i, k) + surf_cooling(i, k)

       end do
    end do
    
    ! record heating rates on grid to history files,
    ! convert from heating rate in J/kg/s to temperature rate of change in K/day
    call outfld('SAI_HEAT', (local_heating(:,:)+surf_cooling(:,:)) / cpair / s2day, ncol, lchnk)
    call outfld('ATTEN_LW', dI_lw(:,:), ncol, lchnk)
    call outfld('ATTEN_SW', dI_sw(:), ncol, lchnk)
    
  end subroutine cldera_sai_tracers_timestep_tend


!===========================================================================
!===========================================================================


  subroutine init_cnst_3d(m, q, gcid)
    
    use dyn_grid, only : get_horiz_grid_d, get_horiz_grid_dim_d
    use dycore,   only : dycore_is

    integer,  intent(in)  :: m       ! global constituent index
    real(r8), intent(out) :: q(:,:)  ! kg tracer/kg dry air (gcol,plev)
    integer,  intent(in)  :: gcid(:) ! global column id

    real(r8), allocatable :: lat(:)
    integer :: plon, plat, ngcols
    integer :: j, k, gsize
    !-----------------------------------------------------------------------

    if (masterproc) then
      write(iulog,*) 'CLDERA SAI CONSTITUENTS: INITIALIZING ',cnst_name(m),m
    end if

    ! initialize everything to zero 
    if (m == ixso2) then
       q(:,:) = 0.0_r8
    else if (m == ixash) then
       q(:,:) = 0.0_r8
    else if (m == ixsulf) then
       q(:,:) = 0.0_r8
    end if
    
  end subroutine init_cnst_3d


!=====================================================================
!=====================================================================


end module cldera_sai_tracers
