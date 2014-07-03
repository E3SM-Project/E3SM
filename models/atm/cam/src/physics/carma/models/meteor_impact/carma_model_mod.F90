!! This module is used to define a particular CARMA microphysical model. For 
!! simple cases, this may be the only code that needs to be modified. This module
!! defines several constants and has three methods:
!!
!!   - CARMA_DefineModel()
!!   - CARMA_EmitParticle()
!!   - CARMA_InitializeParticle()
!!
!! These methods define the microphysical model, the particle emissions and
!! the initial conditions of the particles. Each realization of CARMA
!! microphysics has its own version of this file.
!!
!! This file is used to model a meteor impact upon the land. This model is
!! preliminary. Please talk to Chuck Bardeen (bardeenc@ucar.edu) if you are
!! interested in this model.
!!
!! @version Oct-2012 
!! @author  Chuck Bardeen 
module carma_model_mod

  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carmaelement_mod
  use carmagas_mod
  use carmagroup_mod
  use carmasolute_mod
  use carmastate_mod
  use carma_mod
  use carma_flags_mod
  use carma_model_flags_mod
  
  use shr_kind_mod,   only: r8 => shr_kind_r8
  use abortutils,     only: endrun
  use physics_types,  only: physics_state, physics_ptend
  use ppgrid,         only: pcols, pver
  use physics_buffer, only: physics_buffer_desc

  implicit none

  private

  ! Declare the public methods.
  public CARMA_DefineModel
  public CARMA_Detrain
  public CARMA_DiagnoseBins
  public CARMA_DiagnoseBulk
  public CARMA_EmitParticle
  public CARMA_InitializeModel
  public CARMA_InitializeParticle
  public CARMA_WetDeposition
  
  ! Declare public constants
  integer, public, parameter      :: NGROUP   = 2               !! Number of particle groups
  integer, public, parameter      :: NELEM    = 2               !! Number of particle elements
  integer, public, parameter      :: NBIN     = 21              !! Number of particle bins
  integer, public, parameter      :: NSOLUTE  = 0               !! Number of particle solutes
  integer, public, parameter      :: NGAS     = 0               !! Number of gases


  !! Relative humidities for mie and radiation calculations. The RRTMG radiation code will interpolate
  !! based upon the current relative humidity from a table built using the specified relative
  !! humidities.
  integer, public, parameter      :: NMIE_RH  = 8               !! Number of relative humidities for mie calculations
  real(kind=f), public            :: mie_rh(NMIE_RH) = (/ 0._f, 0.5_f, 0.7_f, 0.8_f, 0.9_f, 0.95_f, 0.98_f, 0.99_f /)
  
  ! Defines whether the groups should undergo deep convection in phase 1 or phase 2.
  ! Water vapor and cloud particles are convected in phase 1, while all other constituents
  ! are done in phase 2.
  logical, public                 :: is_convtran1(NGROUP) = .false.  !! Should the group be transported in the first phase?

  ! Define any particle compositions that are used. Each composition type
  ! should have a unique number.
  integer, public, parameter      :: I_DUST         = 1         !! dust composition
  integer, public, parameter      :: I_SOOT         = 2         !! soot composition

  ! Define group, element, solute and gas indexes.
  integer, public, parameter      :: I_GRP_DUST     = 1         !! dust aerosol group
  integer, public, parameter      :: I_GRP_SOOT     = 2         !! soot aerosol group

  integer, public, parameter      :: I_ELEM_DUST    = 1         !! dust aerosol element
  integer, public, parameter      :: I_ELEM_SOOT    = 2         !! soot aerosol element
  
  
  integer                         :: carma_dustmap(NBIN)        !! mapping of the CARMA dust bins to the surface dust bins.
  real(kind=f)                    :: carma_dustbinfactor(NBIN)  !! bin weighting factor for dust emissions
  real(kind=f)                    :: carma_sootbinfactor(NBIN)  !! bin weighting factor for soot emissions
  real(kind=f)                    :: carma_emis_area            !! surface area where emissions are happening (m2)
  real(kind=f)                    :: carma_emis_dtime           !! duration of the event (s)
contains


  !! Defines all the CARMA components (groups, elements, solutes and gases) and process
  !! (coagulation, growth, nucleation) that will be part of the microphysical model.
  !!
  !!  @version May-2009 
  !!  @author  Chuck Bardeen 
  subroutine CARMA_DefineModel(carma, rc)
    type(carma_type), intent(inout)    :: carma     !! the carma object
    integer, intent(out)               :: rc        !! return code, negative indicates failure
    
    ! Local variables
    real(kind=f), parameter            :: RHO_DUST = 2.0_f     ! density of dust particles (g/cm)
    real(kind=f)                       :: RHO_SOOT             ! density of soot particles (g/cm)
    real(kind=f), parameter            :: dust_rmin = 20.e-7_f ! dust minimum radius (cm)
    real(kind=f), parameter            :: dust_vmrat = 2.49_f  ! dust volume ratio
    real(kind=f), parameter            :: soot_rmin = 20.e-7_f ! dust minimum radius (cm)
    real(kind=f), parameter            :: soot_vmrat = 2.49_f  ! dust volume ratio
    complex(kind=f)                    :: refidx(NWAVE)        ! refractice indices

    integer                            :: LUNOPRT               ! logical unit number for output
    logical                            :: do_print              ! do print output?
    real(kind=f)                       :: soot_rmon = 30.e-7_f  ! soot monomer radius (cm)
    real(kind=f)                       :: soot_df(NBIN) = 2.2_f ! soot fractal dimension
    real(kind=f)                       :: soot_falpha = 1._f    ! soot fractal packing coefficient

    ! Default return code.
    rc = RC_OK

    ! Adjust longitudes to be 0 to 360 rather than +- 180.
    if (carma_emis_minlon < 0._f) carma_emis_minlon = 360._f + carma_emis_minlon
    if (carma_emis_maxlon < 0._f) carma_emis_maxlon = 360._f + carma_emis_maxlon
    
    if (carma_emis_minlat > carma_emis_maxlat) then
      if (do_print) write(LUNOPRT,*) 'CARMA_DefineModel::ERROR - carma_emis_minlat greater than carma_emis_maxlat' 
    end if
    
    ! Report model specific namelist configuration parameters.
    if (masterproc) then
      call CARMA_Get(carma, rc, do_print=do_print, LUNOPRT=LUNOPRT)
      if (rc < 0) call endrun("CARMA_InitializeModel: CARMA_Get failed.")
    
      if (do_print) write(LUNOPRT,*) ''
      if (do_print) write(LUNOPRT,*) 'CARMA ', trim(carma_model), ' specific settings :'
      if (do_print) write(LUNOPRT,*) '  carma_emis_dust       = ', carma_emis_dust, ' (kg)'
      if (do_print) write(LUNOPRT,*) '  carma_emis_soot       = ', carma_emis_soot, ' (kg)'
      if (do_print) write(LUNOPRT,*) '  carma_emis_startdate  = ', carma_emis_startdate
      if (do_print) write(LUNOPRT,*) '  carma_emis_starttime  = ', carma_emis_starttime
      if (do_print) write(LUNOPRT,*) '  carma_emis_stopdate   = ', carma_emis_stopdate
      if (do_print) write(LUNOPRT,*) '  carma_emis_stoptime   = ', carma_emis_stoptime
      if (do_print) write(LUNOPRT,*) '  carma_emis_minlat     = ', carma_emis_minlat
      if (do_print) write(LUNOPRT,*) '  carma_emis_maxlat     = ', carma_emis_maxlat
      if (do_print) write(LUNOPRT,*) '  carma_emis_minlon     = ', carma_emis_minlon
      if (do_print) write(LUNOPRT,*) '  carma_emis_maxlon     = ', carma_emis_maxlon
      if (do_print) write(LUNOPRT,*) '  carma_fractal_soot    = ', carma_fractal_soot
    end if

    ! Define the Groups
    !
    ! NOTE: If NWAVE > 0 then the group should have refractive indices defined.
    !
    ! NOTE: For CAM, the optional do_wetdep and do_drydep flags should be
    ! defined. If wetdep is defined, then the optional solubility factor
    ! should also be defined.

    ! Use the same refractive index at all wavelengths. This value is typical of soot and
    ! is recommended by Toon et al. 2012. TBD Wagner et al. 2011 shows variability in the
    ! real part (0.003 (IR) to 0.05 (UV)).
    refidx(:) = (1.53_f, 0.008_f)
        
    call CARMAGROUP_Create(carma, I_GRP_DUST, "Dust", dust_rmin, dust_vmrat, I_SPHERE, 1._f, .false., &
                           rc, do_wetdep=.true., do_drydep=.true., solfac=0.3_f, &
                           scavcoef=0.1_f, shortname="CRDUST", refidx=refidx, do_mie=.true.)
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddGroup failed.')

    ! Use the same refractive index at all wavelengths. This value is typical of soot and
    ! is recommended by Toon et al. 2012.
    refidx(:) = (1.8_f, 0.67_f)
    
    if (carma_fractal_soot) then
      RHO_SOOT = 1.8_f

      ! This matches the df profile used by Wolf and Toon [2010].
      soot_df(:) = (/ 3.0000_f, 3.0000_f, 1.5033_f, 1.5082_f, 1.5494_f, 1.6168_f, 1.7589_f, &
                      1.9957_f, 2.2519_f, 2.3840_f, 2.4000_f, 2.4000_f, 2.4000_f, 2.4000_f, &
                      2.4000_f, 2.4000_f, 2.4000_f, 2.4000_f, 2.4000_f, 2.4000_f, 2.4000_f /)

      call CARMAGROUP_Create(carma, I_GRP_SOOT, "Soot", soot_rmin, soot_vmrat, I_SPHERE, 1._f, .false., &
                             rc, do_wetdep=.true., do_drydep=.true., solfac=0.1_f, &
                             scavcoef=0.1_f, shortname="CRSOOT", refidx=refidx, do_mie=.true., &
                             is_fractal=.true., rmon=soot_rmon, df=soot_df, falpha=soot_falpha, &
                             imiertn=I_MIERTN_BOTET1997)
    else
      RHO_SOOT = 1.0_f
      call CARMAGROUP_Create(carma, I_GRP_SOOT, "Soot", soot_rmin, soot_vmrat, I_SPHERE, 1._f, .false., &
                             rc, do_wetdep=.true., do_drydep=.true., solfac=0.1_f, &
                             scavcoef=0.1_f, shortname="CRSOOT", refidx=refidx, do_mie=.true.)
    end if
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddGroup failed.')
   
    
    ! Define the Elements
    !
    ! NOTE: For CAM, the optional shortname needs to be provided for the group. These names
    ! should be 6 characters or less and without spaces.
    call CARMAELEMENT_Create(carma, I_ELEM_DUST, I_GRP_DUST, "Dust", RHO_DUST, I_INVOLATILE, I_DUST, rc, shortname="CRDUST")
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddElement failed.')

    call CARMAELEMENT_Create(carma, I_ELEM_SOOT, I_GRP_SOOT, "Soot", RHO_SOOT, I_INVOLATILE, I_SOOT, rc, shortname="CRSOOT")
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddElement failed.')

    
    ! Define the Solutes

    
    ! Define the Gases

    
    ! Define the Processes
    call CARMA_AddCoagulation(carma, I_GRP_DUST, I_GRP_DUST, I_GRP_DUST, I_COLLEC_DATA, rc)
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddCoagulation failed.')

    call CARMA_AddCoagulation(carma, I_GRP_SOOT, I_GRP_SOOT, I_GRP_SOOT, I_COLLEC_DATA, rc)
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddCoagulation failed.')

    return
  end subroutine CARMA_DefineModel


  !! Defines all the CARMA components (groups, elements, solutes and gases) and process
  !! (coagulation, growth, nucleation) that will be part of the microphysical model.
  !!
  !!  @version May-2009 
  !!  @author  Chuck Bardeen 
  !!
  !!  @see CARMASTATE_SetDetrain
  subroutine CARMA_Detrain(carma, cstate, cam_in, dlf, state, icol, dt, rc, rliq, prec_str, snow_str, &
      tnd_qsnow, tnd_nsnow)
    use camsrfexch,         only: cam_in_t
    use physconst,          only: latice, latvap, cpair

    implicit none
    
    type(carma_type), intent(in)         :: carma            !! the carma object
    type(carmastate_type), intent(inout) :: cstate           !! the carma state object
    type(cam_in_t),  intent(in)          :: cam_in           !! surface input
    real(r8), intent(in)                 :: dlf(pcols, pver) !! Detraining cld H20 from convection (kg/kg/s)
    type(physics_state), intent(in)      :: state            !! physics state variables
    integer, intent(in)                  :: icol             !! column index
    real(r8), intent(in)                 :: dt               !! time step (s)
    integer, intent(out)                 :: rc               !! return code, negative indicates failure
    real(r8), intent(inout), optional    :: rliq(pcols)      !! vertical integral of liquid not yet in q(ixcldliq)
    real(r8), intent(inout), optional    :: prec_str(pcols)  !! [Total] sfc flux of precip from stratiform (m/s) 
    real(r8), intent(inout), optional    :: snow_str(pcols)  !! [Total] sfc flux of snow from stratiform (m/s)
    real(r8), intent(out), optional      :: tnd_qsnow(pcols,pver) !! snow mass tendency (kg/kg/s)
    real(r8), intent(out), optional      :: tnd_nsnow(pcols,pver) !! snow number tendency (#/kg/s)
    
    ! Default return code.
    rc = RC_OK
        
    return
  end subroutine CARMA_Detrain


  !! For diagnostic groups, sets up up the CARMA bins based upon the CAM state.
  !!
  !!  @version July-2009 
  !!  @author  Chuck Bardeen 
  subroutine CARMA_DiagnoseBins(carma, cstate, state, pbuf, icol, dt, rc, rliq, prec_str, snow_str)
    use time_manager,     only: is_first_step

    implicit none
    
    type(carma_type), intent(in)          :: carma        !! the carma object
    type(carmastate_type), intent(inout)  :: cstate       !! the carma state object
    type(physics_state), intent(in)       :: state        !! physics state variables
    type(physics_buffer_desc), pointer     :: pbuf(:)      !! physics buffer
    integer, intent(in)                   :: icol         !! column index
    real(r8), intent(in)                  :: dt           !! time step
    integer, intent(out)                  :: rc           !! return code, negative indicates failure
    real(r8), intent(in), optional        :: rliq(pcols)      !! vertical integral of liquid not yet in q(ixcldliq)
    real(r8), intent(inout), optional     :: prec_str(pcols)  !! [Total] sfc flux of precip from stratiform (m/s) 
    real(r8), intent(inout), optional     :: snow_str(pcols)  !! [Total] sfc flux of snow from stratiform (m/s)
    
    real(r8)                             :: mmr(pver) !! elements mass mixing ratio
    integer                              :: ibin      !! bin index
    
    ! Default return code.
    rc = RC_OK
    
    ! By default, do nothing. If diagnosed groups exist, this needs to be replaced by
    ! code to determine the mass in each bin from the CAM state.
    
    return
  end subroutine CARMA_DiagnoseBins
  
  
  !! For diagnostic groups, determines the tendencies on the CAM state from the CARMA bins.
  !!
  !!  @version July-2009 
  !!  @author  Chuck Bardeen 
  subroutine CARMA_DiagnoseBulk(carma, cstate, cam_out, state, pbuf, ptend, icol, dt, rc, rliq, prec_str, snow_str, &
    prec_sed, snow_sed, tnd_qsnow, tnd_nsnow, re_ice)
    use camsrfexch,       only: cam_out_t

    implicit none
    
    type(carma_type), intent(in)         :: carma     !! the carma object
    type(carmastate_type), intent(inout) :: cstate    !! the carma state object
    type(cam_out_t),      intent(inout)  :: cam_out   !! cam output to surface models
    type(physics_state), intent(in)      :: state     !! physics state variables
    type(physics_buffer_desc), pointer   :: pbuf(:)   !! physics buffer
    type(physics_ptend), intent(inout)   :: ptend     !! constituent tendencies
    integer, intent(in)                  :: icol      !! column index
    real(r8), intent(in)                 :: dt        !! time step
    integer, intent(out)                 :: rc        !! return code, negative indicates failure
    real(r8), intent(inout), optional    :: rliq(pcols)      !! vertical integral of liquid not yet in q(ixcldliq)
    real(r8), intent(inout), optional    :: prec_str(pcols)  !! [Total] sfc flux of precip from stratiform (m/s) 
    real(r8), intent(inout), optional    :: snow_str(pcols)  !! [Total] sfc flux of snow from stratiform (m/s)
    real(r8), intent(inout), optional    :: prec_sed(pcols)       !! total precip from cloud sedimentation (m/s)
    real(r8), intent(inout), optional    :: snow_sed(pcols)       !! snow from cloud ice sedimentation (m/s)
    real(r8), intent(inout), optional    :: tnd_qsnow(pcols,pver) !! snow mass tendency (kg/kg/s)
    real(r8), intent(inout), optional    :: tnd_nsnow(pcols,pver) !! snow number tendency (#/kg/s)
    real(r8), intent(out), optional      :: re_ice(pcols,pver)    !! ice effective radius (m)
    
    integer                              :: ielem     ! element index
    integer                              :: ibin      ! bin index
    real(r8)                             :: mmr(pver) ! mass mixing ration (kg/kg)
    real(r8)                             :: sflx      ! surface flux (kg/m2/s)

    ! Default return code.
    rc = RC_OK

    ! Add the sedimentation and dry deposition fluxes to the hydrophilic black carbon.
    !
    ! NOTE: Don't give the surface model negative values for the surface fluxes.
    ielem = I_ELEM_SOOT
    do ibin = 1, NBIN
    
      call CARMASTATE_GetBin(cstate, ielem, ibin, mmr, rc, sedimentationFlux=sflx)
      if (rc < 0) call endrun('CARMA_DiagnoseBulk::CARMA_GetBin failed.')
      
      cam_out%bcphidry(icol) = cam_out%bcphidry(icol) + max(sflx, 0._r8)
    end do

    ielem = I_ELEM_DUST
    do ibin = 1, NBIN
    
      call CARMASTATE_GetBin(cstate, ielem, ibin, mmr, rc, sedimentationFlux=sflx)
      if (rc < 0) call endrun('CARMA_DiagnoseBulk::CARMA_GetBin failed.')
      
      if (carma_dustmap(ibin) == 1) then
        cam_out%dstdry1(icol) = cam_out%dstdry1(icol) + max(sflx, 0._r8)
      else if (carma_dustmap(ibin) == 2) then
        cam_out%dstdry2(icol) = cam_out%dstdry2(icol) + max(sflx, 0._r8)
      else if (carma_dustmap(ibin) == 3) then
        cam_out%dstdry3(icol) = cam_out%dstdry3(icol) + max(sflx, 0._r8)
      else if (carma_dustmap(ibin) == 4) then
        cam_out%dstdry4(icol) = cam_out%dstdry4(icol) + max(sflx, 0._r8)
      end if
    end do
    
    return
  end subroutine CARMA_DiagnoseBulk


  !! Calculates the emissions for CARMA aerosol particles. By default, there is no
  !! emission, but this routine can be overridden for models that wish to have
  !! an aerosol emission.
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  subroutine CARMA_EmitParticle(carma, ielem, ibin, icnst, dt, state, cam_in, tendency, surfaceFlux, rc)
    use shr_kind_mod,  only: r8 => shr_kind_r8
    use ppgrid,        only: pcols, pver
    use physics_types, only: physics_state
    use time_manager,  only: get_curr_date, get_perp_date, get_curr_calday, &
                             is_perpetual, is_first_step
    use camsrfexch,    only: cam_in_t
    use tropopause,    only: tropopause_find
    use physconst,     only: gravit
    
    implicit none
    
    type(carma_type), intent(in)       :: carma                 !! the carma object
    integer, intent(in)                :: ielem                 !! element index
    integer, intent(in)                :: ibin                  !! bin index
    integer, intent(in)                :: icnst                 !! consituent index
    real(r8), intent(in)               :: dt                    !! time step (s)
    type(physics_state), intent(in)    :: state                 !! physics state
    type(cam_in_t), intent(in)         :: cam_in                !! surface inputs
    real(r8), intent(out)              :: tendency(pcols, pver) !! constituent tendency (kg/kg/s)
    real(r8), intent(out)              :: surfaceFlux(pcols)    !! constituent surface flux (kg/m^2/s)
    integer, intent(out)               :: rc                    !! return code, negative indicates failure
    
    real(r8), parameter                :: mu_dust_gnd  = 1._r8  ! width parameter, dust, ground (km)
    real(r8), parameter                :: mu_dust_trop = 3._r8  ! width parameter, dust, tropopause (km)
    real(r8), parameter                :: mu_soot_gnd  = 1._r8  ! width parameter, soot, ground (km)
    real(r8), parameter                :: mu_soot_trop = 3._r8  ! width parameter, soot, tropopause (km)

    integer       :: tropLev(pcols)           ! tropopause level index   
    real(r8)      :: tropP(pcols)             ! tropopause pressure (Pa)  
    real(r8)      :: tropT(pcols)             ! tropopause temperature (K) 
    real(r8)      :: tropZ(pcols)             ! tropopause height (m) 

    real(r8)     :: lon(state%ncol)         ! longitude
    real(r8)     :: lat(state%ncol)         ! latitude
    integer      :: igroup                  ! group index
    integer      :: ncol                    ! number of columns in chunk
    integer      :: icol                    ! column index
    integer      :: k                       ! vertical index
    real(r8)     :: calday                  ! current calendar day
    integer      :: currentDate             ! current date (yyyydoy)
    integer      :: yr                      ! year
    integer      :: mon                     ! month
    integer      :: day                     ! day of month
    integer      :: ncsec                   ! time of day (seconds)
    integer      :: doy                     ! day of year
    real(r8)     :: startyear               ! start year
    real(r8)     :: stopyear                ! stop year
    real(r8)     :: startdoy                ! start year
    real(r8)     :: stopdoy                 ! stop year
    integer      :: emis_time               ! length of time for emission
    real(r8)     :: vfunc(pver)             ! scaling factor to preserve total emission
    character(len=32) :: shortname          ! the shortname of the group
    real(r8)     :: zmid                    ! layer midpoint altitude (km)
    real(r8)     :: ztrop                   ! tropopause altitude (km)
    real(r8)     :: rate                    ! emission rate (kg/s/m)
    real(r8)     :: massflux                ! mass flux (kg/m3/s)
    real(r8)     :: thickness               ! layer thickness (m)

    ! Default return code.
    rc = RC_OK

    ! Determine the day of year.
    calday = get_curr_calday()
    if ( is_perpetual() ) then
      call get_perp_date(yr, mon, day, ncsec)
    else
      call get_curr_date(yr, mon, day, ncsec)
    end if
    doy = floor(calday)

    ! Determine the latitude and longitude of each column.
    ncol = state%ncol
    
    lat  = state%lat(:ncol) * RAD2DEG
    lon  = state%lon(:ncol) * RAD2DEG
    
    
    ! Add any surface flux here.
    surfaceFlux(:ncol) = 0.0_r8
    
    ! For emissions into the atmosphere, put the emission here.
    !
    ! Use Toon et al. [2012] as the source function for soot and dust
    ! from a 1 km meteor impact.
    !
    ! For soot, it is assumed that the soot is emitted in one column
    ! containing the impact and that there are two gaussian
    ! distributions: one centered at the surface and one centered at
    ! the tropopause. The emission rate of soot is given as g/s/km and
    ! we assume that the total mass is delivered in one time step.
    !
    ! NOTE: Perhaps some of these fields should end up in the CARMA
    ! model namelist, so different experiments can be run more easily.
    tendency(:ncol, :pver) = 0.0_r8
    
    ! Determine the start and stop year and day of year from the namelist
    ! variables.
    currentDate = yr * 1000 + doy
    startyear = carma_emis_startdate / 1000
    stopyear  = carma_emis_stopdate  / 1000
    
    startdoy  = mod(carma_emis_startdate, 1000)
    stopdoy   = mod(carma_emis_stopdate, 1000)

    ! Only emit particles during the specified time interval.
    if (((currentDate > carma_emis_startdate) .or. &
         ((currentDate == carma_emis_startdate) .and. (ncsec >= carma_emis_starttime))) .and. &
        ((currentDate < carma_emis_stopdate) .or. &
         ((currentDate == carma_emis_stopdate) .and. (ncsec < carma_emis_stoptime)))) then
    
      ! Make sure to emit for at least one timestep and in multiples of the time
      ! step length.
      ! TBD - This has a leap year problem, but works otherwise ...
      carma_emis_dtime = INT((((stopyear - startyear) * 365._f + (stopdoy - startdoy)) * 24._f * 3600._f + &
                   (carma_emis_stoptime - carma_emis_starttime)) / dt) * dt
  
      ! For simplicity, calculate the emission function at the cell midpoint and
      ! assume that rate is used throughout the cell.
      call CARMAELEMENT_GET(carma, ielem, rc, igroup=igroup)
      if (RC < RC_ERROR) return
      
      call CARMAGROUP_GET(carma, igroup, rc, shortname=shortname)
      if (RC < RC_ERROR) return
      
      if ((shortname == "CRDUST") .or. (shortname == "CRSOOT")) then

        ! Find the tropopause using the default algorithm backed by the climatology.
        call tropopause_find(state, tropLev, tropZ=tropZ)
  
        ! Loop over all of the columns.
        do icol = 1, ncol
  
          ! Is the column one of the ones over which there should be emissions>
          if ((lat(icol) > carma_emis_minlat) .and. (lat(icol) < carma_emis_maxlat) .and. &
              (((carma_emis_minlon <= carma_emis_maxlon) .and. (lon(icol) >= carma_emis_minlon) .and. &
              (lon(icol) <= carma_emis_maxlon)) .or. &
              ((carma_emis_minlon > carma_emis_maxlon)  .and. &
              ((lon(icol) >= carma_emis_minlon) .or. (lon(icol) <= carma_emis_maxlon))))) then
  
            ! Set tendencies for any sources or sinks in the atmosphere.
            do k = 1, pver
            
              ! Get the cell midpoint and height
              zmid  = state%zm(icol, k) / 1000._f
    
              ! Get the tropopause height.
              ztrop = tropZ(icol) / 1000._f
  
              ! Use the dust emission from Toon et al. 2012.      
              if (shortname == "CRDUST") then
            
                ! Determine the total emission rate for this grid box using equation 2
                ! from Toon et al. [2012] and also adjust for the fraction of the
                ! mass that goes into the specified bin based on the assumed size
                ! distribution also from Toon et al. [2012].
                vfunc(k) = 1._f / (2._f * sqrt(2._f * PI)) * &
                      (1._f / mu_dust_gnd  * exp(-0.5_f * ((zmid           / mu_dust_gnd)**2)) + &
                       1._f / (2._f * mu_dust_trop) * exp(-0.5_f * (((zmid - ztrop) / mu_dust_trop)**2))) * &
                       (state%zi(icol, k) - state%zi(icol, k+1))
                       
                rate = carma_emis_dust * carma_dustbinfactor(ibin)
              end if
    
              ! Use the soot emissions from Toon et al. 2012.
              if (shortname == "CRSOOT") then
            
                ! Determine the total emission rate for this grid box using equation 2
                ! from Toon et al. [2012] and also adjust for the fraction of the
                ! mass that goes into the specified bin based on the assumed size
                ! distribution also from Toon et al. [2012].
                vfunc(k) = 1._f / (2._f * sqrt(2._f * PI)) * &
                      (1._f / mu_soot_gnd  * exp(-0.5_f * ((zmid           / mu_soot_gnd)**2)) + &
                       1._f / (2._f * mu_soot_trop) * exp(-0.5_f * (((zmid - ztrop) / mu_soot_trop)**2))) * &
                       (state%zi(icol, k) - state%zi(icol, k+1))
                
                
                rate = carma_emis_soot * carma_sootbinfactor(ibin)
              end if
                
              ! Calculate a rate by dividing by total emission time.
              rate = rate  * vfunc(k) / carma_emis_dtime
                
              ! Scale for the fraction of the total surface area that is emitting and
              ! convert to kg/m2/s
              massflux = rate / carma_emis_area
              
              ! Convert the mass flux to a tendency on the mass mixing ratio.
              tendency(icol, k) = massflux / (state%pdel(icol, k) / gravit)
            end do
            
            ! Now normalize in the vertical to preserve the total mass.
            tendency(icol, :) = tendency(icol, :) / sum(vfunc(:))
          end if
        end do
      end if
    end if
    
    return
  end subroutine CARMA_EmitParticle


  !! Allows the model to perform its own initialization in addition to what is done
  !! by default in CARMA_init.
  !!
  !! NOTE: If CARMA constituents appear in the initial condition file, then those
  !! values will override anything set here.
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  subroutine CARMA_InitializeModel(carma, lq_carma, rc)
    use constituents, only: pcnst
    use dyn_grid, only: get_horiz_grid_dim_d, get_horiz_grid_d

    implicit none
    
    type(carma_type), intent(in)       :: carma                 !! the carma object
    logical, intent(inout)             :: lq_carma(pcnst)       !! flags to indicate whether the constituent
                                                                !! could have a CARMA tendency
    integer, intent(out)               :: rc                    !! return code, negative indicates failure
    
    ! NOTE: The dust distribution has not been specified yet, but it should be different
    ! from the soot.
    real(kind=f), parameter            :: rm_dust    = 0.11     ! dust mean radius (um)
    real(kind=f), parameter            :: sigma_dust = 1.6      ! dust variance
    real(kind=f), parameter            :: rm_soot    = 0.11     ! soot mean radius (um)
    real(kind=f), parameter            :: sigma_soot = 1.6      ! soot variance

    integer                            :: i
    integer                            :: hdim1_d
    integer                            :: hdim2_d
    integer                            :: ngcols
    real(kind=f)                       :: r(NBIN)
    real(kind=f)                       :: dr(NBIN)
    real(kind=f)                       :: rmass(NBIN)
    real(kind=f)                       :: dM(NBIN)
    real(kind=f), allocatable          :: lat(:)
    real(kind=f), allocatable          :: lon(:)
    real(kind=f), allocatable          :: colarea(:)
    character(len=32)                  :: shortname             ! the shortname of the group
    
    integer                            :: LUNOPRT               ! logical unit number for output
    logical                            :: do_print              ! do print output?

  1 format(i3,5x,i3,4x,e10.3,4x,e10.3) 

    ! Default return code.
    rc = RC_OK

    ! Create a mapping of the CARMA dust bins to the dust sizes assumed at the
    ! surface. The sizes of the dust bins at the surface are from Mahowald et al.
    ! [2006].
    !
    !   1 :  0.1 - 1.0 um
    !   2 :  1.0 - 2.5 um
    !   3 :  2.5 - 5.0 um
    !   4 :  5.0 - 10.0 um
    call CARMAGROUP_GET(carma, I_GRP_DUST, rc, r=r)
    if (RC < RC_ERROR) return
    
    do i = 1, NBIN
      if (r(i) .le. 1e-4_f) then
        carma_dustmap(i)  = 1
      else if (r(i) .le. 2.5e-4_f) then
        carma_dustmap(i) = 2
      else if (r(i) .le. 5e-4_f) then
        carma_dustmap(i) = 3
      else
        carma_dustmap(i) = 4
      end if
    end do
    
    ! Determine the weight of mass in each bin based upon the size distribution specified
    ! in Toon et al. [2012], for soot and dust. They are lognormal for the smaller sizes
    ! and dust is lognormal for larger sizes.
    
    call CARMAGROUP_GET(carma, I_GRP_DUST, rc, shortname=shortname, r=r, dr=dr, rmass=rmass)
    if (RC < RC_ERROR) return
    
    dM(:)        = rmass(:) * &
         exp(-(log(r(:) * 1e4_f / rm_dust) ** 2) / (2._f * (log(sigma_dust) ** 2))) / &
         log(sigma_dust) * (dr(:) / r(:))
    carma_dustbinfactor(:)  = dM / sum(dM) 

    call CARMAGROUP_GET(carma, I_GRP_SOOT, rc, shortname=shortname, r=r, dr=dr, rmass=rmass)
    if (RC < RC_ERROR) return
    
    dM(:)        = rmass(:) * &
         exp(-(log(r(:) * 1e4_f / rm_soot) ** 2) / (2._f * (log(sigma_soot) ** 2))) / &
         log(sigma_soot) * (dr(:) / r(:))
    carma_sootbinfactor(:)  = dM / sum(dM)

    
    ! Determine the total area in which debris will be emitted. This is used to scale
    ! the emission per column, based upon the fraction of surface area. This assumes a
    ! regular physics grid.
    call get_horiz_grid_dim_d(hdim1_d, hdim2_d)
    
    ngcols = hdim1_d*hdim2_d
    
    allocate(lat(ngcols))
    allocate(lon(ngcols))
    allocate(colarea(ngcols))
    
    call get_horiz_grid_d(ngcols, clat_d_out=lat, clon_d_out=lon, area_d_out=colarea)

    lat  = lat * RAD2DEG
    lon  = lon * RAD2DEG

    ! rad2 -> m2
    colarea = colarea * REARTH * REARTH / 1e4

    ! Integrate surface area with same checks as in the emission routine to determine
    ! the area where the emissions come from (m2). Assume that the grid box is either
    ! all in or all out based upon the center lat/lon. Don't include fractions of a 
    ! grid box.
    carma_emis_area = 0._f
    
    do i = 1, ngcols
      if ((lat(i) >= carma_emis_minlat) .and. (lat(i) <= carma_emis_maxlat) .and. &
          (((carma_emis_minlon <= carma_emis_maxlon) .and. (lon(i) >= carma_emis_minlon) .and. &
          (lon(i) <= carma_emis_maxlon)) .or. &
          ((carma_emis_minlon > carma_emis_maxlon)  .and. &
          ((lon(i) >= carma_emis_minlon) .or. (lon(i) <= carma_emis_maxlon))))) then
        carma_emis_area = carma_emis_area + colarea(i)
      end if
    end do
    
    carma_emis_area = carma_emis_area 
    
    deallocate(lat)
    deallocate(lon)
    deallocate(colarea)

    ! Report model specific namelist configuration parameters.
    if (masterproc) then
      call CARMA_Get(carma, rc, do_print=do_print, LUNOPRT=LUNOPRT)
      if (rc < 0) call endrun("CARMA_InitializeModel: CARMA_Get failed.")
    
      
      if (do_print) then
        write(LUNOPRT,*) ''
        write(LUNOPRT,*) 'CARMA Initialization ...'
        
        write(LUNOPRT,*) ''
        write(LUNOPRT,*) 'ibin  dustmap  dustfactor    sootfactor'

        do i = 1, NBIN
          write(LUNOPRT,1) i, carma_dustmap(i), carma_dustbinfactor(i), carma_sootbinfactor(i)
        end do

        write(LUNOPRT,*) ''
        write(LUNOPRT,*) '  Emission area     :  ', carma_emis_area / 1e6_f, ' (km^2)'
        write(LUNOPRT,*) ''

      end if
    end if

    return
  end subroutine CARMA_InitializeModel


  !! Sets the initial condition for CARMA aerosol particles. By default, there are no
  !! particles, but this routine can be overridden for models that wish to have an
  !! initial value.
  !!
  !! NOTE: If CARMA constituents appear in the initial condition file, then those
  !! values will override anything set here.
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  subroutine CARMA_InitializeParticle(carma, ielem, ibin, q, gcid, rc)
    use shr_kind_mod,   only: r8 => shr_kind_r8
    use pmgrid,         only: plat, plev, plon

    implicit none
    
    type(carma_type), intent(in)       :: carma                 !! the carma object
    integer, intent(in)                :: ielem                 !! element index
    integer, intent(in)                :: ibin                  !! bin index
    real(r8), intent(inout)            :: q(:,:)   ! kg tracer/kg dry air (gcol, plev)
    integer, intent(in)                :: gcid(:)  ! global column id
    integer, intent(out)               :: rc                    !! return code, negative indicates failure

    ! Default return code.
    rc = RC_OK

    ! Add initial condition here.
    
    return
  end subroutine CARMA_InitializeParticle

    
  !!  Called after wet deposition has been performed. Allows the specific model to add
  !!  wet deposition of CARMA aerosols to the aerosols being communicated to the surface.
  !!
  !!  @version July-2011 
  !!  @author  Chuck Bardeen 
  subroutine CARMA_WetDeposition(carma, ielem, ibin, sflx, cam_out, state, rc)
    use camsrfexch,       only: cam_out_t

    implicit none
    
    type(carma_type), intent(in)         :: carma       !! the carma object
    integer, intent(in)                  :: ielem       !! element index
    integer, intent(in)                  :: ibin        !! bin index
    real(r8), intent(in)                 :: sflx(pcols) !! surface flux (kg/m2/s)
    type(cam_out_t), intent(inout)       :: cam_out     !! cam output to surface models
    type(physics_state), intent(in)      :: state       !! physics state variables
    integer, intent(out)                 :: rc          !! return code, negative indicates failure
    
    integer    :: icol
 
    ! Default return code.
    rc = RC_OK
    
    ! Add the wet deposition fluxes to the hydrophilic black carbon.
    !
    ! NOTE: Don't give the surface model negative values for the surface fluxes.
    if (ielem == I_ELEM_SOOT) then
      do icol = 1, state%ncol
        cam_out%bcphiwet(icol) = cam_out%bcphiwet(icol) + max(sflx(icol), 0._r8)
      end do
    end if

    if (ielem == I_ELEM_DUST) then
      do icol = 1, state%ncol
        if (carma_dustmap(ibin) == 1) then
          cam_out%dstwet1(icol) = cam_out%dstwet1(icol) + max(sflx(icol), 0._r8)
        else if (carma_dustmap(ibin) == 2) then
          cam_out%dstwet2(icol) = cam_out%dstwet2(icol) + max(sflx(icol), 0._r8)
        else if (carma_dustmap(ibin) == 3) then
          cam_out%dstwet3(icol) = cam_out%dstwet3(icol) + max(sflx(icol), 0._r8)
        else if (carma_dustmap(ibin) == 4) then
          cam_out%dstwet4(icol) = cam_out%dstwet4(icol) + max(sflx(icol), 0._r8)
        end if
      end do
    end if
    
    return
  end subroutine CARMA_WetDeposition 
  
end module
