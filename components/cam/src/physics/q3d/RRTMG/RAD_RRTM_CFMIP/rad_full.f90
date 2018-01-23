      SUBROUTINE rad_full ()
 
      use rad ! note qrad, pres_input, tabs_slice, insolation_TOA, lwUp, etc from this module
      use rad_driver, only: &
                      p_factor_xy, p_coszrs_xy, &
                      rad_driver_rrtm, initialize_radiation, &
                      isInitialized_RadDriver, tracesini
      use parkind, only: &
                   kind_rb ! RRTM expects reals with this kind parameter (8 byte reals) 

!===========================================================================
!  Modified for use with the Lorenz grid physics model.
!  Thomas Cram and Celal Konor, CSU, October 2009.
!===========================================================================

!===========================================================================
!================== BEGIN CHANGES REQUIRED HERE ============================
!===========================================================================
! PORTABILITY NOTE:  Here, all of the stuff needed to call the radiation is
!    drawn from various modules in SAM as an example.  You will need to bring
!    all of the associated variables into this routine for your model 
!    (either by using the appropriate modules, passing them in as arguments 
!    or defining them here).

! the following are logicals
  use rrtm_grid, only: &
       dolongwave, &       ! do longwave radiation
       doshortwave, &      ! do shortwave radiation
       doperpetual, &      ! use diurnally-averaged insolation
       dosolarconstant, &  ! specify mean insolation, zenith angle for perpetual insolation
       doseasons, &        ! allow diurnally-varying insolation to vary with time of year
       ocean, &            ! if true, run is over ocean.
       masterproc, &       ! true if MPI rank==0.
       dostatisrad         ! accumulate radiation statistics at this time step

! the following are characters
  use rrtm_grid, only: &
       iopfile, &    ! name of SCAM IOP forcing file, e.g. 'ctl_s11.nc'
       case          ! used to construct path of SCAM IOP file in SAM

! the following are integers
  use rrtm_grid, only: &
       nx, ny, nzm, &     ! number of grid points (x,y,z)
       nstep, icycle, &   ! model step and substep number
       nrad, &            ! call radiation every nrad timesteps
       iyear              ! current year

! the following are reals or real arrays
  use rrtm_params, only: &
       cp, &                   ! specific heat of dry air at constant pressure, J/kg/K
       ggr, &                  ! gravitational acceleration, m/s2
       secday, &               ! seconds in one day
       coszrs, &               ! cosine of solar zenith angle
       latitude, &            ! latitude in degrees
       longitude              ! longitude in degrees

  use rrtm_grid, only: &
       day, day0, &       ! model day (current and at start of run) day=0. for 00Z, Jan 1.
       dz, adz, &         ! vertical grid spacing is dz*adz(k) in this model
       solar_constant, &  ! modified solar constant if doperpetual==dosolarconstant==.true.
       zenith_angle       ! zenith angle if doperpetual==dosolarconstant==.true.

! NOTE:  when dosolarconstat==.true, insolation = solar_constant*cos(zenith_angle)

  use rrtm_vars, only: &
!       t, &       ! model thermodynamic variable
       tabs, &    ! absolute temperature (K)
       qv, &      ! vapor mixing ratio
       qcl, &     ! cloud mixing ratio
       qci, &     ! ice mixing ratio
       sstxy, &   ! sea surface temperature
       pres, &    ! model layer pressure (mb)
       presi, &   ! model interface pressure (mb)
!       rho, &     ! density profile.  In this anelastic model, rho=rho(z).
       radswup, &  ! SW upward flux statistic, summed in x and y
       radswdn, &  ! SW downward flux statistic, summed in x and y
       radqrsw, &  ! SW heating rate statistic, summed in x and y
       radqrcsw, & ! SW clear sky heating rate statistic, summed in x and y.
       radlwup, &  ! LW upward flux statistic, summed in x and y
       radlwdn, &  ! LW downward flux statistic, summed in x and y
       radqrlw, &  ! LW heating rate statistic, summed in x and y
       radqrclw, & ! LW clear sky heating rate statistic, summed in x and y.
       s_flnt, &   ! Net LW upward flux at TOA statistic
       s_fsnt, &   ! Net SW downward flux at TOA statistic
       s_flntc, &  ! Net LW clearsky upward flux at TOA statistic
       s_fsntc, &  ! Net SW clearsky downward flux at TOA statistic
       s_solin, &  ! Insolation at TOA statistic 
       s_flns, &   ! Net LW upward flux at surface statistic
       s_fsns, &   ! Net SW downward flux at surface statistic
       s_flnsc, &  ! Net LW clearsky upward flux at surface statistic
       s_fsnsc, &  ! Net SW clearsky downward flux at surface statistic
       s_fsds, &   ! Downwelling SW flux at surface statistic
       s_flds, &   ! Downwelling LW flux at surface statistic
       lwnt_xy, &  ! Time-accumulated x-y field of net LW flux at TOA
       swnt_xy, &  ! Time-accumulated x-y field of net SW flux at TOA
       lwntc_xy, & ! Time-accumulated x-y field of net clearsky LW flux at TOA
       swntc_xy, & ! Time-accumulated x-y field of net clearsky SW flux at TOA
       solin_xy, & ! Time-accumulated x-y field of insolation at TOA
       lwns_xy, &  ! Time-accumulated x-y field of net LW upward flux at surface
       swns_xy, &  ! Time-accumulated x-y field of net SW downward flux at surface
       lwnsc_xy, & ! Time-accumulated x-y field of net clearsky LW upward flux at surface
       swnsc_xy    ! Time-accumulated x-y field of net clearsky SW downward flux at surface

!-----------------------------------------------------------------------------
! VVM modules

      USE trace_gases, only: &
                       o3, co2, n2o, ch4, o2, cfc11, cfc12, cfc22, ccl4
      
  !===========================================================================
  !================== END CHANGES REQUIRED HERE ==============================
  !===========================================================================

  implicit none

! parameters
!  integer, parameter :: iyear = 2003

! local variables

  logical :: updateRadiation = .true.
  logical :: update_TraceGases_PatchedSounding = .true.
  
  logical :: isRestart ! true if run is restarted from previous simulation.
  character(LEN=250) :: RestartFileName, SoundingFileName

!  real(kind=kind_rb), dimension(nx,nzm) :: &
!       swHeatingRate, & ! units: K/s
!       lwHeatingRate, & !  units: K/s
!       swHeatingRateClearSky, & !  units: K/s
!       lwHeatingRateClearSky !  units: K/s
  
  integer :: ierr, i, j, k

!=====================================================================
! PORTABILITY NOTE: all processors should have same pressure soundings
!    (pres and presi) in hPa.  This is important for the consistency
!     of trace soundings across the different processors.
!=====================================================================
  
! only call radiation roughly every minute.
!   Each model will likely have their own way of doing this.
!   NOTE: In SAM, nstep is the timestep number,
!          and icycle is the substep number within a timestep.

  if(((mod(nstep,nrad).eq.1).AND.(icycle.eq.1)) &
       .OR.(.NOT.isInitialized_RadDriver)) then
    updateRadiation = .true.
  else
    updateRadiation = .false.
  end if
  if(masterproc)print *,'updateradiation ',updateradiation

  if(updateRadiation) then
    
!========== CALL INITIALIZATION ROUTINE FOR RADIATION ============

!------------------------------------------------------
! Initialize if necessary 

    if(.not. isInitialized_RadDriver) then

! First allocate perpetual factor and perpetual cosine solar zenith angle
      ALLOCATE(p_factor_xy(nx,ny), p_coszrs_xy(nx,ny), STAT = ierr)
      IF (ierr .NE. 0) THEN
        WRITE(*,*) 'Could not allocate fixed size arrays p_factor_xy and p_coszrs_xy in rad_full'
        CALL rad_error()
      ENDIF

      CALL initialize_radiation(nx, ny, cp, iyear, day0, latitude, longitude, doperpetual)
      update_TraceGases_PatchedSounding = .true.
            
!------------------------------------------------------
! allocate arrays that are fixed in size for the whole simulation.
      ALLOCATE( &
           qrad(nx,ny,nzm), &
           lwUp_3d(nx,ny,nzm), &
           lwDown_3d(nx,ny,nzm), &
           swUp_3d(nx,ny,nzm), &
           swDown_3d(nx,ny,nzm), &
           lwHeatingRate_3d(nx,ny,nzm), &
           swHeatingRate_3d(nx,ny,nzm), &
           lwp_3d(nx,ny,nzm), &
           iwp_3d(nx,ny,nzm), &
           reliq_3d(nx,ny,nzm), &
           reice_3d(nx,ny,nzm), &
           NetlwUpSurface(nx,ny), &
           NetlwUpSurfaceClearSky(nx,ny), &
           NetlwUpToa(nx,ny), &
           NetlwUpToaClearSky(nx,ny), &
           NetswDownSurface(nx,ny), &
           NetswDownSurfaceClearSky(nx,ny), &
           NetswDownToa(nx,ny), &
           NetswDownToaClearSky(nx,ny), &
           NetswUpToa(nx,ny), &
           insolation_TOA(nx,ny), &
           swDownSurface(nx,ny), &
           lwDownSurface(nx,ny), &
           swnsxy(nx,ny), &
           lwnsxy(nx,ny), &
           STAT=ierr)
        if(ierr.ne.0) then
          write(*,*) 'Could not allocate fixed size arrays in rad_full'
          call rad_error()
        end if

      isInitialized_RadDriver = .true.
    end if ! if(.NOT.isInitialized_RadDriver)
    
! update trace gas sounding and patch to background sounding hourly

!tcram: trace gas profiles are held constant throughout simulation
!    if(day-day_when_patch_tracegases_last_updated.gt.3599.999/86400.) then
!      update_TraceGases_PatchedSounding = .true.
!    end if

!==============================================================================
    if(.NOT.isAllocated_RadInputsOutputs.OR.update_TraceGases_PatchedSounding) then

! Initialize or update two things:
!    - patch in a sounding above the model top if necessary.
!    - interpolate the trace gas concentrations to grid used to compute radiation.

! Read in sounding (Minghua) and define set up patching.
!   This will be stored in radiation restart file for restarts.
! PORTABILIY: CHANGE SOUNDING FILENAME

!tcram: trace gas soundings are acquired via the rad_variables_tendencies module

!      SoundingFileName = './'//trim(case)//'/'//trim(iopfile) ! e.g., './RUNDATA/ctl_s11.nc' OR './RUNDATA/p2k_s11.nc'
!      call read_patch_background_sounding(SoundingFileName,presi(nzm+1), &
!           npatch_start,npatch_end,nzsnd,psnd,tsnd,qsnd,masterproc)
!      if(npatch_end.ne.npatch_start) then
!        nzpatch = npatch_end - npatch_start + 1
!      else
        nzpatch = 0
!      end if

      nzrad = nzm + nzpatch

      if(isAllocated_RadInputsOutputs.AND.nzrad.ne.nzrad_old) then

! deallocate old arrays

!tcram: Added the following allocatable variables: o3_slice, co2_slice,
!       ch4_slic, n2o_slice, o2_slice, cfc11_slice, cfc12_slice, cfc22_slice,
!       ccl4_slice, latitude_slice, longitude_slice, p_factor_slice, p_coszrs_slice

        DEALLOCATE( &
             tabs_slice, &
             qv_slice, &
             qcl_slice, &
             qci_slice, &
             tg_slice, &
             o3_slice, &
             co2_slice, &
             ch4_slice, &
             n2o_slice, &
             o2_slice, &
             cfc11_slice, &
             cfc12_slice, &
             cfc22_slice, &
             ccl4_slice, &
             latitude_slice, &
             longitude_slice, &
             p_factor_slice, &
             p_coszrs_slice, &
             pres_input, &
             presi_input, &
             lwUp, &
             lwDown, &
             lwUpClearSky, &
             lwDownClearSky, &
             swUp, &
             swDown, &
             swUpClearSky, &
             swDownClearSky, &
             swHeatingRate, &
             swHeatingRateClearSky, &
             lwHeatingRate, &
             lwHeatingRateClearSky, &
             LWP, IWP, &
             liquidRe, iceRe, &
             STAT=ierr)
        if(ierr.ne.0) then
          write(*,*) 'Could not deallocate input/output arrays in rad_full'
          call rad_error()
        else
          isAllocated_RadInputsOutputs = .false.
        end if
      end if

      if(.NOT.isAllocated_RadInputsOutputs) then
! allocate arrays.
        ALLOCATE( &
             tabs_slice(nx,nzrad), &
             qv_slice(nx,nzrad), &
             qcl_slice(nx,nzrad), &
             qci_slice(nx,nzrad), &
             tg_slice(nx), &
             o3_slice(nzrad+1), &
             co2_slice(nzrad+1), &
             ch4_slice(nzrad+1), &
             n2o_slice(nzrad+1), &
             o2_slice(nzrad+1), &
             cfc11_slice(nzrad+1), &
             cfc12_slice(nzrad+1), &
             cfc22_slice(nzrad+1), &
             ccl4_slice(nzrad+1), &
             latitude_slice(nx), &
             longitude_slice(nx), &
             p_factor_slice(nx), &
             p_coszrs_slice(nx), &
             pres_input(nzrad), &
             presi_input(nzrad+1), &
             lwUp(nx,nzrad+2), &
             lwDown(nx,nzrad+2), &
             lwUpClearSky(nx,nzrad+2), &
             lwDownClearSky(nx,nzrad+2), &
             swUp(nx,nzrad+2), &
             swDown(nx,nzrad+2), &
             swUpClearSky(nx,nzrad+2), &
             swDownClearSky(nx,nzrad+2), &
             swHeatingRate(nx,nzrad+1), &
             swHeatingRateClearSky(nx,nzrad+1), &
             lwHeatingRate(nx,nzrad+1), &
             lwHeatingRateClearSky(nx,nzrad+1), &
             LWP(nx, nzm+1), IWP(nx, nzm+1), &
             liquidRe(nx,nzm+1), iceRe(nx, nzm+1), &
             STAT=ierr)
        if(ierr.ne.0) then
          write(*,*) 'Could not allocate input/output arrays in rad_full'
          call rad_error()
        else
          isAllocated_RadInputsOutputs = .true.
          nzrad_old = nzrad
        end if
      end if

! set up pressure inputs to radiation -- needed for initialize_radiation
      pres_input(1:nzm) = pres(1:nzm)
      presi_input(1:nzm+1) = presi(1:nzm+1)
      
      if(nzpatch.gt.0) then
        pres_input(nzm+1:nzrad) = psnd(npatch_start:npatch_end)
        presi_input(nzm+2:nzrad) = &
             0.5*(psnd(npatch_start:npatch_end-1) &
             + psnd(npatch_start+1:npatch_end))
        presi_input(nzrad+1) = MAX(0.5*psnd(npatch_end), &
             1.5*psnd(npatch_end) - 0.5*psnd(npatch_end-1))
      end if

! interpolates standard sounding of trace gas concentrations to grid for radiation.
!tcram: trace gases acquired from main model via the rad_variables_tendencies module
!      call tracesini(nzrad,pres_input,presi_input,ggr,masterproc)
 
      update_TraceGases_PatchedSounding = .false.
      day_when_patch_tracegases_last_updated = day

    end if ! if(update_TraceGases_PatchedSounding.OR..NOT.isAllocated_RadInputsOutputs)
!==============================================================================

! zero out radiation statistics that are summed in x- and y-directions.

    radlwup(:) = 0.
    radlwdn(:) = 0.
    radqrlw(:) = 0.
    radqrclw(:) = 0.
    radswup(:) = 0.
    radswdn(:) = 0.
    radqrsw(:) = 0.
    radqrcsw(:) = 0.

! set up pressure inputs to radiation

    pres_input(1:nzm) = pres(1:nzm)
    presi_input(1:nzm+1) = presi(1:nzm+1)

    if(nzpatch.gt.0) then
      pres_input(nzm+1:nzrad) = psnd(npatch_start:npatch_end) ! layer pressures
      presi_input(nzm+2:nzrad) = & ! interface pressures.
           0.5*(psnd(npatch_start:npatch_end-1) &
           + psnd(npatch_start+1:npatch_end))
      presi_input(nzrad+1) = MAX(0.5*psnd(npatch_end), &
                                 1.5*psnd(npatch_end) - 0.5*psnd(npatch_end-1))
    end if

!==============================================================================
!tcram: The RRTMG code takes a 1-D vector of columns, so loop over the y-index
!==============================================================================
!$omp parallel do default(shared) &
!$omp      private(j, tabs_slice, qv_slice, qcl_slice, qci_slice, tg_slice,    &
!$omp              o3_slice, co2_slice, ch4_slice, n2o_slice, o2_slice,        &
!$omp              cfc11_slice, cfc12_slice, cfc22_slice, ccl4_slice,          &
!$omp              latitude_slice, longitude_slice, p_factor_slice, p_coszrs_slice,   &
!$omp              lwUp, lwDown, lwUpClearSky, lwDownClearSky,                 &
!$omp              swUp, swDown, swUpClearSky, swDownClearSky,                 &
!$omp              swHeatingRate, swHeatingRateClearSky, lwHeatingRate, lwHeatingRateClearSky, &
!$omp              LWP, IWP, liquidRe, iceRe)
    do 1000 j = 1,ny

      ! extract a slice from the three-dimensional domain on this processor.
      !   We need absolute temperature (K), mass mixing ratios (kg/kg) of
      !   water vapor, cloud liquid water and cloud ice, along with SST (K).

      tabs_slice(1:nx,1:nzm) = tabs(1:nx,j,1:nzm)
      qv_slice(1:nx,1:nzm) = qv(1:nx,j,1:nzm)
      qcl_slice(1:nx,1:nzm) = qcl(1:nx,j,1:nzm)
      qci_slice(1:nx,1:nzm) = qci(1:nx,j,1:nzm)
      tg_slice(1:nx) = sstxy(1:nx,j)
      o3_slice(1:nzm) = o3(1,j,1:nzm)
      co2_slice(1:nzm) = co2(1,j,1:nzm)
      ch4_slice(1:nzm) = ch4(1,j,1:nzm)
      n2o_slice(1:nzm) = n2o(1,j,1:nzm)
      o2_slice(1:nzm) = o2(1,j,1:nzm)
      cfc11_slice(1:nzm) = cfc11(1,j,1:nzm)
      cfc12_slice(1:nzm) = cfc12(1,j,1:nzm)
      cfc22_slice(1:nzm) = cfc22(1,j,1:nzm)
      ccl4_slice(1:nzm) = ccl4(1,j,1:nzm)
      
      latitude_slice(1:nx) = latitude(1:nx,j)
      longitude_slice(1:nx) = longitude(1:nx,j)
      
      p_factor_slice(1:nx) = p_factor_xy(1:nx,j)
      p_coszrs_slice(1:nx) = p_coszrs_xy(1:nx,j)

! patch sounding on top of model sounding for more complete radiation calculation.
      if(nzpatch.gt.0) then
        tabs_slice(1:nx,nzm+1:nzrad) = spread( tsnd(npatch_start:npatch_end), dim=1, ncopies=nx )
        qv_slice(1:nx,nzm+1:nzrad) = spread( qsnd(npatch_start:npatch_end), dim=1, ncopies=nx )
        qcl_slice(1:nx,nzm+1:nzrad) = 0.
        qci_slice(1:nx,nzm+1:nzrad) = 0.
      end if

!------------------------------------------------------------------------------
! Make call to wrapper routine for RRTMG (v.4.8 for LW, v.3.8 for SW)

      call rad_driver_rrtm(nx,nzrad,j,pres_input,presi_input, &
           tabs_slice,qv_slice,qcl_slice,qci_slice,tg_slice, &
           o3_slice,co2_slice,ch4_slice,n2o_slice,o2_slice, &
           cfc11_slice,cfc12_slice,cfc22_slice,ccl4_slice, &
           dolongwave,doshortwave,doperpetual,doseasons, &
           dosolarconstant,solar_constant,zenith_angle, &
           day,day0,latitude_slice,longitude_slice, &
           p_factor_slice, p_coszrs_slice, &
           ocean,ggr,cp, masterproc, &
           lwUp,lwDown,lwUpClearSky,lwDownClearSky, &
           swUp,swDown,swUpClearSky,swDownClearSky, &
           swHeatingRate, &
           swHeatingRateClearSky, &
           lwHeatingRate, &
           lwHeatingRateClearSky, &
           coszrs, &
           LWP, IWP, liquidRe, iceRe )

!------------------------------------------------------------------------------

! Compute heating rates from fluxes using local density in model.
! Results in heating rates in K/s.
! PORTABILIY NOTE: CHANGE THERMAL MASS TO THAT USED IN YOUR MODEL.
!      Units below are cp*rho*deltaz ~ J/kg/K * kg/m3 * m
!      where delta z = dz*adz(k) in SAM.

!tcram: Use heating rates returned from the RRTMG calculation

!      do k = 1,nzm ! loop over model levels
!        swHeatingRate(1:nx,k) = &
!             (swDown(:,k+1) - swDown(:,k) + swUp(:,k) - swUp(:,k+1)) &
!             /(cp*rho(k)*dz*adz(k))
!        lwHeatingRate(1:nx,k) = &
!             (lwDown(:,k+1) - lwDown(:,k) + lwUp(:,k) - lwUp(:,k+1)) &
!             /(cp*rho(k)*dz*adz(k))
!        swHeatingRateClearSky(1:nx,k) = &
!             (swDownClearSky(:,k+1) - swDownClearSky(:,k) &
!             + swUpClearSky(:,k) - swUpClearSky(:,k+1)) &
!             /(cp*rho(k)*dz*adz(k))
!        lwHeatingRateClearSky(1:nx,k) = &
!             (lwDownClearSky(:,k+1) - lwDownClearSky(:,k) &
!             + lwUpClearSky(:,k) - lwUpClearSky(:,k+1)) &
!             /(cp*rho(k)*dz*adz(k))
!      end do

!------------------------------------------------------------------------------
! update total radiative heating rate of model

!tcram: convert units to K/s

      qrad(1:nx,j,1:nzm) = (swHeatingRate(1:nx,1:nzm) + lwHeatingRate(1:nx,1:nzm)) / secday

!---------------------------------------------------------------------
! Load shortwave and longwave heating rates into proper VVM 3-D arrays

      swHeatingRate_3d(1:nx, j, 1:nzm) = swHeatingRate(1:nx,1:nzm) / secday
      lwHeatingRate_3d(1:nx, j, 1:nzm) = lwHeatingRate(1:nx,1:nzm) / secday
      lwUp_3d(1:nx, j, 1:nzm)          = lwUp(1:nx,1:nzm)
      lwDown_3d(1:nx, j, 1:nzm)        = lwDown(1:nx,1:nzm)
      swUp_3d(1:nx, j, 1:nzm)          = swUp(1:nx,1:nzm)
      swDown_3d(1:nx, j, 1:nzm)        = swDown(1:nx,1:nzm)
      lwp_3d(1:nx, j, 1:nzm)           = LWP(1:nx,1:nzm)
      iwp_3d(1:nx, j, 1:nzm)           = IWP(1:nx,1:nzm)
      reliq_3d(1:nx, j, 1:nzm)         = liquidRe(1:nx,1:nzm)
      reice_3d(1:nx, j, 1:nzm)         = iceRe(1:nx,1:nzm)
       
!------------------------------------------------------------------------------
! accumulate heating rates and fluxes for horizontally-averaged statistics
!$omp critical
      radlwup(:) = radlwup(:) + sum(lwUp(:, 1:nzm+1),   dim = 1)
      radlwdn(:) = radlwdn(:) + sum(lwDown(:, 1:nzm+1), dim = 1)
      radqrlw(1:nzm) = radqrlw(1:nzm) + sum(lwHeatingRate(:, 1:nzm), dim = 1)
      radqrclw(1:nzm) = radqrclw(1:nzm) + sum(lwHeatingRateClearSky(:, 1:nzm), dim = 1)
      radswup(:) = radswup(:) + sum(swUp(:, 1:nzm+1),   dim = 1)
      radswdn(:) = radswdn(:) + sum(swDown(:, 1:nzm+1), dim = 1)
      radqrsw(1:nzm) = radqrsw(1:nzm) + sum(swHeatingRate(:, 1:nzm), dim = 1)
      radqrcsw(1:nzm) = radqrcsw(1:nzm) + sum(swHeatingRateClearSky(:, 1:nzm), dim = 1)
!$omp end critical

      ! shortwave fluxes at top-of-atmosphere (TOA) and surface -- NOTE POSITIVE DOWNWARDS
      insolation_TOA(1:nx,j) = swDown(:,nzrad+2) ! shortwave down at TOA
      swDownSurface(1:nx,j) = swDown(:,1) ! shortwave down at surface

      NetswUpToa(1:nx,j) = swUp(:,nzrad+2) - swDown(:,nzrad+2) ! net shortwave up at TOA

      NetswDownToa(1:nx,j) = swDown(:,nzrad+2) - swUp(:,nzrad+2) ! net shortwave down at TOA
      NetswDownToaClearSky(1:nx,j) = swDownClearSky(:,nzrad+2) - swUpClearSky(:,nzrad+2) ! net clearsky shortwave down at TOA

      NetswDownSurface(1:nx,j) = swDown(:,1) - swUp(:,1) ! net shortwave down at surface
      NetswDownSurfaceClearSky(1:nx,j) = swDownClearSky(:,1) - swUpClearSky(:,1) ! net clearsky shortwave down at surface

      ! longwave fluxes at top-of-atmosphere (TOA) and surface -- NOTE POSITIVE UPWARDS
      lwDownSurface(1:nx,j) = lwDown(:,1) ! longwave down at surface

      NetlwUpToa(1:nx,j) = lwUp(:,nzrad+2) - lwDown(:,nzrad+2) ! net longwave up at TOA
      NetlwUpToaClearSky(1:nx,j) = lwUpClearSky(:,nzrad+2) - lwDownClearSky(:,nzrad+2) ! net clearsky longwave up at TOA

      NetlwUpSurface(1:nx,j) = lwUp(:,1) - lwDown(:,1) ! net longwave up at surface
      NetlwUpSurfaceClearSky(1:nx,j) = lwUpClearSky(:,1) - lwDownClearSky(:,1) ! net clearsky longwave up at surface

 1000 continue ! j = 1,ny
!$omp end parallel do
!==============================================================================

  end if ! if(updateRadiation)

!------------------------------------------------------------------------------

  if(icycle.eq.1) then
    ! Update 2d diagnostic fields 
    !
    ! Net surface and toa fluxes
    !
    ! First two for ocean evolution

    lwnsxy(:, :) = NetlwUpSurface(:, :) ! instantaneous
    swnsxy(:, :) = NetswDownSurface(:, :)

    ! net full sky radiative fluxes (varying in x and y, time-accumulated)
    lwns_xy(:, :) = lwns_xy(:, :) + NetlwUpSurface(:, :)
    swns_xy(:, :) = swns_xy(:, :) + NetswDownSurface(:, :)
    lwnt_xy(:, :) = lwnt_xy(:, :) + NetlwUpToa(:, :) 
    swnt_xy(:, :) = swnt_xy(:, :) + NetswDownToa(:, :)

    ! net clear sky radiative fluxes (varying in x and y, time-accumulated)
    lwnsc_xy(:, :) = lwnsc_xy(:, :) + NetlwUpSurfaceClearSky(:, :)
    swnsc_xy(:, :) = swnsc_xy(:, :) + NetswDownSurfaceClearSky(:, :)
    lwntc_xy(:, :) = lwntc_xy(:, :) + NetlwUpToaClearSky(:, :) 
    swntc_xy(:, :) = swntc_xy(:, :) + NetswDownToaClearSky(:, :)

    ! TOA Insolation
    solin_xy(:, :) = solin_xy(:, :) + insolation_Toa(:, :) 
  end if

!------------------------------------------------------------------------------

! Update 1D diagnostics

  if(dostatisrad) then
    s_flns = s_flns + sum(NetlwUpSurface(:, :))   ! lwnsxy
    s_fsns = s_fsns + sum(NetswDownSurface(:, :)) ! swnsxy
    s_flnt = s_flnt + sum(NetlwUpToa(:, :))       ! lwntxy
    s_fsnt = s_fsnt + sum(NetswDownToa(:, :))     ! swntxy
    s_flnsc = s_flnsc + sum(NetlwUpSurfaceClearSky(:, :))   ! lwnscxy
    s_fsnsc = s_fsnsc + sum(NetswDownSurfaceClearSky(:, :)) ! swnscxy 
    s_flntc = s_flntc + sum(NetlwUpToaClearSky(:, :))       ! lwntcxy
    s_fsntc = s_fsntc + sum(NetswDownToaClearSky(:, :))     ! swntcxy
    s_solin = s_solin + sum(insolation_TOA(:, :))           ! solinxy 
    ! 
    s_fsds = s_fsds + sum(swDownSurface(:, :)) ! downwelling SW at surface
    s_flds = s_flds + sum(lwDownSurface(:, :)) ! downwelling LW at surface
  end if ! if(dostatis)

end subroutine rad_full

