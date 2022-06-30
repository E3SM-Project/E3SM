!===============================================================================
! CLDERA stratospheric aerosol injection tracers
! provides dissipation rates, heatings rates, and injection forcings for diagnostic constituents
!
! Joe Hollowed
! June 2022
! This module written based on the structure of aoa_tracers, and enables the advection and
! evolution of 5 tracers constituents:
! - SO2
! - ASH
! - SULFATE
! - SAI_PT  : potential temperature
! _ SAI_PV  : potential vorticity
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
  public :: cldera_sai_tracers_implements_cnst  ! true if named constituent is implemented by this package
  public :: cldera_sai_tracers_init_cnst        ! initialize constituent field
  public :: cldera_sai_tracers_init             ! initialize history fields, datasets
  public :: cldera_sai_tracers_timestep_init    ! place to perform per timestep initialization
  public :: cldera_sai_tracers_timestep_tend    ! calculate tendencies
  public :: cldera_sai_tracers_readnl           ! read namelist options

  !----  Private module data

  integer, parameter :: ncnst=5  ! number of constituents implemented by this module

  ! constituent names, indices
  character(len=8), parameter :: c_names(ncnst) = (/'SO2', 'ASH', 'SULFATE','SAI_PT', 'SAI_PV'/)
  integer :: ifirst ! global index of first constituent
  integer :: ixso2  ! global index for SO2 tracer
  integer :: ixash  ! global index for ASH tracer
  integer :: ixsulf ! global index for SULFATE tracer
  integer :: ixpt   ! global index for SAI_PT tracer
  integer :: ixpv   ! global index for SAI_PV tracer

  ! Data from namelist variables; defaults set in bld/namelist_files/namelist_defaults_eam.xml
  logical  :: cldera_sai_tracers_flag      ! true => activate module, set namelist variable
  logical  :: cldera_sai_read_from_ic_file ! true => tracers initialized from IC file
  real(r8) :: cldera_sai_lat0              ! desired lat of injection (deg)
  real(r8) :: cldera_sai_lon0              ! desired lon of injection (deg)
  real(r8) :: cldera_sai_MSO2              ! total SO2 mass (Mt)
  real(r8) :: cldera_sai_Mash              ! total ash mass (Mt)
  real(r8) :: cldera_sai_tf                ! injection duration (hours)
  real(r8) :: cldera_sai_rkSO2             ! SO2 e-folding (1/day)
  real(r8) :: cldera_sai_rkash             ! ash e-folding (1/day)
  real(r8) :: cldera_sai_rksulf            ! ash e-folding (1/day)
  real(r8) :: cldera_sai_w                 ! SO2->sulfate reaction mass weighting
  real(r8) :: cldera_sai_dTstrat           ! stratospheric heating rate (K/day)
  real(r8) :: cldera_sai_dTsurf            ! surface cooling rate (K/day)
  real(r8) :: cldera_sai_qstar             ! mixing ratio normalization for strat. heating
  real(r8) :: cldera_sai_taustar           ! AOD normalizaxtion for surface cooling
  real(r8) :: cldera_sai_zAOD              ! max height to apply surface cooling by AOD (km)
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
  real(r8), parameter :: Mt2kg   = 1.0e9_r8
  real(r8), parameter :: km2m    = 1000._r8
  real(r8), parameter :: hPa2Pa  = 100._r8
  real(r8), parameter :: P0      = 1000 * hPa2Pa    ! PT reference presssure (Pa)
  
  ! for injection source localization
  real(r8) :: lat0                      ! lat of injection (rad)
  real(r8) :: lon0                      ! lon of injection (rad)           
  real(r8) :: tf                        ! injection duration (s)
  real(r8) :: ttol  = 60.0_r8           ! time tolerance (s)
  real(r8) :: alpha = -2.0_r8           ! vertical skewness of plume 
  real(r8) :: mu    = 22.59_r8 * km2m   ! peak injection altitude (m)
  real(r8) :: sigma = 4.0_r8 * km2m     ! vertical width parameter (m)
  real(r8) :: zmin  = 17.0_r8 * km2m    ! lower plume truncation (m)
  real(r8) :: zmax  = 30.0_r8 * km2m    ! upper plume truncation (m)
  real(r8) :: rr_tmp                    ! great circle distance (m)
  real(r8) :: source_cutoff            ! masking
  
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
  real(r8) :: dTstrat                   ! stratospheric heating rate (K/s)
  real(r8) :: qstar                     ! mixing ratio norm. for strat. heating
  
  ! for surface cooling
  real(r8) :: dTsurf                    ! surface cooling rate (K/s) 
  real(r8) :: taustar                   ! AOD norm. for surface cooling 
  real(r8) :: zAOD                      ! max height to apply surf. cooling (m)
  real(r8) :: b_so2  = 1.0_r8           ! SO2 mass extinction coeff. (1/kg)
  real(r8) :: b_ash  = 1.0_r8           ! ash mass extinction coeff. (1/kg)
  real(r8) :: b_sulf = 1.0_r8           ! sulfate mass extinction coeff. (1/kg)
  real(r8) :: AOD_cutoff               ! masking
  

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
                                     cldera_sai_lat0, cldera_sai_lon0, cldera_sai_tf, &
                                     cldera_sai_MSO2, cldera_sai_Mash, &
                                     cldera_sai_rkSO2, cldera_sai_rkash, cldera_sai_rksulf, &
                                     cldera_sai_dTsurf, cldera_sai_dTstrat, &
                                     cldera_sai_w, cldera_sai_taustar, cldera_sai_zAOD, &
                                     cldera_sai_qstar, & 
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
       lat0    = cldera_sai_lat0    * deg2rad     
       lon0    = cldera_sai_lon0    * deg2rad     
       tf      = cldera_sai_tf      * hr2s          
       M_so2   = cldera_sai_MSO2    * Mt2kg      
       M_ash   = cldera_sai_Mash    * Mt2kg      
       dTstrat = cldera_sai_dTstrat / day2s  
       dTsurf  = cldera_sai_dTsurf  / day2s   
       zAOD    = cldera_sai_zAOD    * km2m
       k_so2   = 1.0_r8 / (cldera_sai_rkSO2  * day2s)
       k_ash   = 1.0_r8 / (cldera_sai_rkash  * day2s)     
       k_sulf  = 1.0_r8 / (cldera_sai_rksulf * day2s)     
       qstar   = cldera_sai_qstar            
       w       = cldera_sai_w                 
       taustar = cldera_sai_taustar          

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
    call mpibcast(M_so2,         1, mpir8, 0, mpicom)
    call mpibcast(M_ash,         1, mpir8, 0, mpicom)
    call mpibcast(k_so2,         1, mpir8, 0, mpicom)
    call mpibcast(k_ash,         1, mpir8, 0, mpicom)
    call mpibcast(k_sulf,        1, mpir8, 0, mpicom)
    call mpibcast(tf,            1, mpir8, 0, mpicom)
    call mpibcast(w,             1, mpir8, 0, mpicom)
    call mpibcast(dTstrat,       1, mpir8, 0, mpicom)
    call mpibcast(dTsurf,        1, mpir8, 0, mpicom)
    call mpibcast(qstar,         1, mpir8, 0, mpicom)
    call mpibcast(taustar,       1, mpir8, 0, mpicom)
    call mpibcast(zAOD,          1, mpir8, 0, mpicom)
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
                  longname='Stratospheric aerosol injection plume SO2')
    ifirst = ixso2
    call cnst_add(c_names(2), mwdry, cpair, 0._r8, ixash, readiv=cldera_sai_read_from_ic_file, &
                  longname='Stratospheric aerosol injection plume ash')
    call cnst_add(c_names(3), mwdry, cpair, 0._r8, ixsulf, readiv=cldera_sai_read_from_ic_file, &
                  longname='Stratospheric aerosol injection plume sulfate')
    call cnst_add(c_names(4), mwdry, cpair, 0._r8, ixpt,   readiv=cldera_sai_read_from_ic_file, &
                  longname='potential temperature at initial time of stratospheric aerosol injection')
    call cnst_add(c_names(5), mwdry, cpair, 0._r8, ixpv,   readiv=cldera_sai_read_from_ic_file, &
                  longname='potential vorticity at initial time of stratospheric aerosol injection')
    !call cnst_add(c_names(6), mwdry, cpair, 0._r8, ixaoa,   readiv=cldera_sai_read_from_ic_file, &
    !              longname='stratospheric aeosol injection plume clock tracer')

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
    ! Purpose: initialize age of air constituents
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
    
    call addfld('AREA',  horiz_only, 'A', 'kg', 'area of grid cell' )
    call add_default('AREA', 1, ' ')
    call addfld('AIR_MASS',  (/ 'lev' /), 'A', 'kg', 'mass of grid box' )
    call add_default('AIR_MASS', 1, ' ')
    call addfld('SAI_HEAT',  (/ 'lev' /), 'A', 'K/day', 'stratospheric heating by sulfur' )
    call add_default('SAI_HEAT', 1, ' ')
    call addfld('SAI_COOL',  (/ 'lev' /), 'A', 'K/day', 'surface cooling by AOD' )
    call add_default('SAI_COOL', 1, ' ')
    call addfld('COL_MASS',  (/ 'lev' /), 'A', 'kg', 'total mass in column above position' )
    call add_default('COL_MASS', 1, ' ')
    call addfld('SAI_AOD',  (/ 'lev' /), 'A', 'dimensionless', 'total AOD from all constituents' )
    call add_default('SAI_AOD', 1, ' ')

  end subroutine cldera_sai_tracers_init


!===============================================================================
!===============================================================================


  subroutine cldera_sai_tracers_timestep_init( phys_state )
    !-----------------------------------------------------------------------
    ! Provides a place to reinitialize diagnostic constituents at each timestep
    ! Currently does nothing; this is a template
    !-----------------------------------------------------------------------

    use time_manager,   only: get_curr_date
    use ppgrid,         only: begchunk, endchunk
    use physics_types,  only: physics_state

    type(physics_state), intent(inout), dimension(begchunk:endchunk), optional :: phys_state


    integer c, i, k, ncol
    integer yr, mon, day, tod
    !--------------------------------------------------------------------------

    if (.not. cldera_sai_tracers_flag) return

    ! currently does nothing, return
    return

    call get_curr_date (yr,mon,day,tod)

    if ( day == 1 .and. tod == 0) then
       if (masterproc) then
         write(iulog,*) 'CLDERA SAI: RE-INITIALIZING CONSTITUENTS'
       endif

       do c = begchunk, endchunk
          ncol = phys_state(c)%ncol
          do k = 1, pver
             do i = 1, ncol
                ! do something?
                continue
             end do
          end do
       end do

    end if

  end subroutine cldera_sai_tracers_timestep_init


!===============================================================================
!===============================================================================


  subroutine cldera_sai_tracers_timestep_tend(state, ptend, dt, ncol)

    use physics_types,   only: physics_state, physics_ptend, physics_ptend_init
    use phys_grid,       only: get_area_p, get_area_all_p, phys_grid_find_col  
    use physconst,       only: rearth, cpair, rair, rgrav=>rga, rearth
    use time_manager,    only: get_curr_time
    use time_manager,    only: get_nstep
    use cam_history,     only: outfld
    use cam_abortutils,  only: endrun 
    use ref_pres,        only: ptop_ref, psurf_ref
    use mpishorthand

    ! Arguments
    type(physics_state), intent(inout) :: state              ! state variables
    type(physics_ptend), intent(out)   :: ptend              ! package tendencies
    real(r8),            intent(in)    :: dt                 ! timestep
    integer,             intent(in)    :: ncol               ! no. of column in chunk

    !----------------- Local workspace-------------------------------

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
    real(r8) :: V(ncol, pver)                 ! vertical profile (1/m)
    real(r8) :: sVk   = 0.0_r8                ! vertical profile sum (1/m)
    real(r8) :: rr    = 1.0e20                ! great circle distance (m)
    integer  :: inj_owner                     ! MPI rank owning injection column
    integer  :: inji                          ! injection column index within chunk
    integer  :: injc                          ! injection local chunk index
    real(r8) :: inj_dist                      ! distance of matched col from requested location (m)
    logical  :: inject                        ! true => this rank contains the injection column
    integer  :: myrank                        ! rank of this MPI process
    integer  :: ier                           ! return error status

    ! for mass normalization
    real(r8) :: A_so2 = 0.0                   ! SO2 normalization (kg m / s)
    real(r8) :: A_ash = 0.0                   ! ash normalization (kg m / s)
    
    ! for mass estimate
    real(r8) :: area(ncol)                    ! column horizontal area (m^2)
    real(r8) :: grid_mass(ncol,pver)          ! grid cell mass (kg)
    real(r8) :: col_mass(ncol,pver)           ! cumulative column mass (kg)
    
    ! for stratospheric heating
    real(r8) :: strat_heating(ncol,pver)      ! stratospheric heating (J/kg/s)
    
    ! for surface cooling
    real(r8) :: surf_cooling(ncol,pver)       ! surface cooling (J/kg/s)
    real(r8) :: aod_so2(ncol,pver)       ! SO2 AOD (kg)
    real(r8) :: aod_ash(ncol,pver)       ! ash AOD (kg)
    real(r8) :: aod_sulf(ncol,pver)      ! sulfate AOD (kg)
    
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
    lq(ixpt)    = .TRUE.
    lq(ixpv)    = .TRUE.
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
    inject = myrank == inj_owner .and. lchnk == injc

    ! great circle distance checks how far column is from requested injection site
    inj_dist = rearth * ACOS( SIN(state%lat(inji))*SIN(lat0) + COS(state%lat(inji))*COS(lat0) * &
                                                               COS(ABS(state%lon(inji)-lon0))) 
    if(inject) then
        ! likely not masterproc, so output goes to e3sm.log rather than atm.log
        write(iulog,*) "CLDERA_SAI_TRACERS injection owner", inj_owner, ", ", myrank
        write(iulog,*) "CLDERA_SAI_TRACERS injection owner col", inji
        write(iulog,*) "CLDERA_SAI_TRACERS injection owner chunk", injc
        write(iulog,*) "CLDERA_SAI_TRACERS injection owner lat0", state%lat(inji) * rad2deg 
        write(iulog,*) "CLDERA_SAI_TRACERS injection owner lon0", state%lon(inji) * rad2deg
        write(iulog,*) "CLDERA_SAI_TRACERS injection owner dist (km)", inj_dist / 1000.0
    endif
   

    ! =============== COMPUTE VERTICAL PROFILE, NORMALIZATION ===============
    ! all entires in V are initialized to zero, then update values at levels k
    ! only for the injection column at inji of rank inj_owner. 
    ! For myrank == inj_owner , V(:,k) will effectively serve as a 2d "mask" 
    ! for the injection at level k on the local chunk of columns, where one
    ! column at index inji will have nonzero values 
    ! For myrank /= inj_owner, V(:,:) = 0, and also the normalization 
    ! A_so2, A_ash is left as it's initialized value of 0.0
    V(:,:) = 0.0_r8
    if(inject) then
        do k = 1,pver
            zz  = state%zm(inji, k)
            if(zz >= zmin .and. zz <= zmax) then
                V(inji, k) = 1.0_r8/(SQRT(2.0_r8*pi)*sigma) * &
                             EXP(-(zz - mu)**2.0_r8/(2.0_r8*sigma**2.0_r8)) * &
                             (1.0_r8- ERF(alpha * (mu - zz)/(SQRT(2.0_r8)*sigma)))
            endif
        enddo
        sVk = SUM(V(inji, :))
        A_so2 = M_so2 / (tf * sVk) 
        A_ash = M_ash / (tf * sVk)
        write(iulog,*) "CLDERA_SAI_TRACERS normalization A (kg m / s)", A_so2
        write(iulog,*) "CLDERA_SAI_TRACERS normalization sVk (1/m)", sVk
    endif
         
        
    ! =============== COMPUTE GRID, COLUMN MASS, AOD =============== 
    call get_area_all_p(lchnk, ncol, area)
    area = area * rearth**2
    
    do i = 1, ncol
        do k = 1,pver
            ! mass in grid cell via hydrostatic approximation
            ! pdel = Pa, area = rad^2, rearth = m, rgrav = s**2/m ===> mass = kg
            grid_mass(i,k) = state%pdel(i,k) * area(i) * rgrav
            
            ! take cumulative sum for the column mass above position k
            if(k == 1) then ! model top
                col_mass(i,k) = grid_mass(i, k)
            else
                col_mass(i,k) = col_mass(i, k-1) + grid_mass(i, k)
            endif

            ! mass extinction coefficients times cumulative sum of tracer mass 
            ! above position k for AOD
            if(k == 1) then
                aod_so2(i, k)  = b_so2  * grid_mass(i,k) * state%q(i,k,ixso2)
                aod_sulf(i, k) = b_sulf * grid_mass(i,k) * state%q(i,k,ixsulf)
                aod_ash(i, k)  = b_ash  * grid_mass(i,k) * state%q(i,k,ixash)
            else
                aod_so2(i, k)  = aod_so2(i,k-1)  + b_so2  * grid_mass(i,k) * state%q(i,k,ixso2)
                aod_sulf(i, k) = aod_sulf(i,k-1) + b_sulf * grid_mass(i,k) * state%q(i,k,ixsulf)
                aod_ash(i, k)  = aod_ash(i,k-1)  + b_ash  * grid_mass(i,k) * state%q(i,k,ixash)
            endif
        enddo

        ! check consistency of cumulative mass calculation against reference pressures at
        ! model top, surface; raise warning if disagrees more than 5%
        col_mass_check = (psurf_ref - ptop_ref) * area(i) * rgrav
        if ( ABS(col_mass(i,pver) -col_mass_check) / col_mass_check  > 0.05_r8) then
             write(iulog,*) "CLDERA_SAI_TRACERS cumulative mass error at col ", i
             write(iulog,*) "   expected:               ", col_mass_check
             write(iulog,*) "   measured cumulative:    ", col_mass(i, pver)
             write(iulog,*) "   relative err:           ",   &
                            100.0 * ABS(col_mass(i, pver) - col_mass_check)/col_mass_check, "%"
             ! uncomment to end model run if this check fails
             !call endrun()
        end if
    enddo
    
    ! record air mass on grid to history files
    call outfld('AREA',     area(:), ncol, lchnk)
    call outfld('AIR_MASS', grid_mass(:,:), ncol, lchnk)
    call outfld('COL_MASS', col_mass(:,:), ncol, lchnk)
    call outfld('SAI_AOD',  aod_so2(:,:) + aod_sulf(:,:) + aod_ash(:,:), ncol, lchnk)


    ! =============== COMPUTE TENDENCIES ===============
    do i = 1, ncol
       do k = 1, pver
          
          lat = state%lat(i)
          lon = state%lon(i)
          
          ! ---- constrain injection to duration tf, bound AOD to zmax
          ! tf expected to be divisible by timestep, so compare with tolerance ttol << dt
          if(ABS(t - tf) > ttol) then
              source_cutoff = 0.0
          else
              source_cutoff = 1.0
          end if
          if(state%zm(i,k) > zAOD) then
              AOD_cutoff = 0.0
          else
              AOD_cutoff = 1.0
          end if


          ! =============== SO2 ===============
          ! ---- source + decay
          ptend%q(i,k,ixso2) = 1/grid_mass(i, k) * &
                               (-k_so2 * state%q(i, k, ixso2) * grid_mass(i, k) + &
                                A_so2 * V(i, k) * source_cutoff)
          ! for debugging
          !if(ptend%q(i, k, ixso2) * grid_mass(i,k) > 1) then
          !    write(iulog,*) "CLDERA_SAI_TRACERS PTENDSUM (kg/s) ", inject, myrank, lchnk, i, &
          !                                         ptend%q(i, k, ixso2) * grid_mass(i,k), t, tf
          !endif
          

          ! =============== ASH ===============
          ! ---- source + decay
          ptend%q(i,k,ixash) = 1/grid_mass(i, k) * &
                               (-k_ash * state%q(i, k, ixash) * grid_mass(i, k) + &
                                A_ash * V(i, k) * source_cutoff)


          ! =============== SULFATE ===============
          if(cldera_sai_formSulfate) then 
              ! ---- source + decay
              ptend%q(i,k,ixsulf) = -k_sulf * state%q(i, k, ixsulf) + &
                                    w * k_so2 * state%q(i, k, ixso2)
          else
              ptend%q(i,k,ixsulf) = 0.0_r8
          endif

          
          ! =============== POTENTIAL TEMP ===============
          ! initialize within the first minute of the run
          if (t < 60.0_r8) then
              state%q(i,k,ixpt) = state%t(i,k) * (P0 / state%pmid(i, k))**(rair/cpair)
          end if
          ptend%q(i,k,ixpt) = 0.0_r8
                   

          ! =============== POTENTIAL VORT ===============
          ! currently does nothing...
          ptend%q(i,k,ixpv) = 0.0_r8


          ! =============== DIABATIC HEATING ===============
          ! ---- stratospheric heating
          if(cldera_sai_stratHeating) then 
              strat_heating(i, k) = (state%q(i,k,ixso2) + state%q(i,k,ixsulf)) * &
                                    (1/qstar) * cpair * dTstrat
          else
              strat_heating(i, k) = 0.0_r8
          endif
    
          ! ---- surface cooling
          if(cldera_sai_surfCooling) then 
              surf_cooling(i, k) = (aod_so2(i,k) + aod_sulf(i,k) + aod_ash(i,k)) * &
                                   (1/taustar) * cpair * dTsurf * AOD_cutoff
          else
              surf_cooling(i, k) = 0.0_r8
          endif

          ! --- total heating
          ptend%s(i, k) = strat_heating(i, k) + surf_cooling(i, k)


       end do
    end do
    
    ! record heating rates on grid to history files,
    ! convert from heating rate in J/kg/s to temperature rate of change in K/day
    call outfld('SAI_HEAT', strat_heating(:,:) / cpair / s2day, ncol, lchnk)
    call outfld('SAI_COOL', surf_cooling(:,:) / cpair / s2day, ncol, lchnk)
    
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
    else if (m == ixpt) then
       q(:,:) = 0.0_r8
    else if (m == ixpv) then
       q(:,:) = 0.0_r8
    end if
    
  end subroutine init_cnst_3d


!=====================================================================
!=====================================================================


end module cldera_sai_tracers
