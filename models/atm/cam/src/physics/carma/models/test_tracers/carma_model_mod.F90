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
!! This file is a test case that uses CARMA groups and bins to implement a
!! tracer trajectory test for the Asian Monsoon region. This is the reverse of 
!! back trajectory calculations being done by John Bergman. In this model each
!! group is a region of the model and each bin represents a day. Emissions
!! start on the carma_launch_doy and then continue for NBINS days.
!!
!! NOTE: This test can use a lot of advected constituents. If you want to reduce
!! the number of regions or days tracked, you also need to reduce the number of
!! advected constituents added in configure.
!!
!! @version April-2011 
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
  use radconstants,   only: nswbands, nlwbands
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
  integer, public, parameter      :: NGROUP   = 6               !! Number of particle groups
  integer, public, parameter      :: NELEM    = 6               !! Number of particle elements
  integer, public, parameter      :: NBIN     = 62              !! Number of particle bins
  integer, public, parameter      :: NSOLUTE  = 0               !! Number of particle solutes
  integer, public, parameter      :: NGAS     = 0               !! Number of gases
  
  ! These need to be defined, but are only used when the particles are radiatively active.
  integer, public, parameter      :: NMIE_RH  = 8               !! Number of relative humidities for mie calculations
  real(kind=f), public            :: mie_rh(NMIE_RH)
  
  ! Defines whether the groups should undergo deep convection in phase 1 or phase 2.
  ! Water vapor and cloud particles are convected in phase 1, while all other constituents
  ! are done in phase 2.
  logical, public                 :: is_convtran1(NGROUP) = .False.  !! Should the group be transported in the first phase?

  ! Define any particle compositions that are used. Each composition type
  ! should have a unique number.
  integer, public, parameter      :: I_INERT   = 1               !! tracer composition
    
  real(kind=f), public            :: rgn_minlat(NELEM-1)    = (/  0._f,  0._f,  0._f,  0._f,  0._f /)
  real(kind=f), public            :: rgn_maxlat(NELEM-1)    = (/ 40._f, 40._f, 40._f, 40._f, 40._f /)
  
  real(kind=f), public            :: rgn_minlon(NELEM-1)    = (/  60._f,  60._f, 105._f,  60._f, 105._f /)
  real(kind=f), public            :: rgn_maxlon(NELEM-1)    = (/ 105._f, 105._f, 140._f, 105._f, 140._f /)

  real(kind=f), public            :: rgn_ps(NELEM)          = (/ -75000._f, 75000._f, 0._f, 0._f, 0._f, 0._f /)

  logical, public                 :: rgn_doLand(NELEM)      = (/ .True.,  .True.,  .True.,  .False., .False., .True. /)
  logical, public                 :: rgn_doOcean(NELEM)     = (/ .False., .False., .False., .True.,  .True.,  .True. /)
  logical, public                 :: rgn_doSeaIce(NELEM)    = (/ .False., .False., .False., .True.,  .True.,  .True. /)

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
    real(kind=f), parameter            :: rmin     = 2.5e-4_f  ! minimum radius (cm)
    real(kind=f), parameter            :: vmrat    = 1.00001_f ! volume ratio
    integer                            :: LUNOPRT
    logical                            :: do_print
    
    ! Default return code.
    rc = RC_OK
    
    call CARMA_Get(carma, rc, do_print=do_print, LUNOPRT=LUNOPRT)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_Get failed.')
    
    ! Report model specific configuration parameters.
    if (masterproc) then
      if (do_print) then
        write(LUNOPRT,*) ''
        write(LUNOPRT,*) 'CARMA ', trim(carma_model), ' specific settings :'
        write(LUNOPRT,*) '  carma_launch_doy     = ', carma_launch_doy
        write(LUNOPRT,*) '  carma_emission_rate  = ', carma_emission_rate
      end if
    end if

    ! Define the Groups
    !
    ! NOTE: If NWAVE > 0 then the group should have refractive indices defined.
    !
    ! NOTE: For CAM, the optional do_wetdep and do_drydep flags should be
    ! defined. If wetdep is defined, then the optional solubility factor
    ! should also be defined.
    call CARMAGROUP_Create(carma, 1, "Region 1", rmin, vmrat, I_SPHERE, 1._f, .True., rc, shortname="CRRG1")
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAGROUP_Create failed.')
    
    call CARMAGROUP_Create(carma, 2, "Region 2", rmin, vmrat, I_SPHERE, 1._f, .True., rc, shortname="CRRG2")
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAGROUP_Create failed.')
   
    call CARMAGROUP_Create(carma, 3, "Region 3", rmin, vmrat, I_SPHERE, 1._f, .True., rc, shortname="CRRG3")
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAGROUP_Create failed.')

    call CARMAGROUP_Create(carma, 4, "Region 4", rmin, vmrat, I_SPHERE, 1._f, .True., rc, shortname="CRRG4")
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAGROUP_Create failed.')

    call CARMAGROUP_Create(carma, 5, "Region 5", rmin, vmrat, I_SPHERE, 1._f, .True., rc, shortname="CRRG5")
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAGROUP_Create failed.')

    call CARMAGROUP_Create(carma, 6, "Rest of World", rmin, vmrat, I_SPHERE, 1._f, .True., rc, shortname="CRRG6")
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAGROUP_Create failed.')
   
    
    ! Define the Elements
    !
    ! NOTE: For CAM, the optional shortname needs to be provided for the group. These names
    ! should be 6 characters or less and without spaces.
    call CARMAELEMENT_Create(carma, 1, 1, "Region 1", WTMOL_AIR, I_INVOLATILE, I_INERT, rc, shortname="CRRG1")
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAElement_Create failed.')
    
    call CARMAELEMENT_Create(carma, 2, 2, "Region 2", WTMOL_AIR, I_INVOLATILE, I_INERT, rc, shortname="CRRG2")
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAElement_Create failed.')
    
    call CARMAELEMENT_Create(carma, 3, 3, "Region 3", WTMOL_AIR, I_INVOLATILE, I_INERT, rc, shortname="CRRG3")
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAElement_Create failed.')

    call CARMAELEMENT_Create(carma, 4, 4, "Region 4", WTMOL_AIR, I_INVOLATILE, I_INERT, rc, shortname="CRRG4")
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAElement_Create failed.')

    call CARMAELEMENT_Create(carma, 5, 5, "Region 5", WTMOL_AIR, I_INVOLATILE, I_INERT, rc, shortname="CRRG5")
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAElement_Create failed.')

    call CARMAELEMENT_Create(carma, 6, 6, "Rest of World", WTMOL_AIR, I_INVOLATILE, I_INERT, rc, shortname="CRRG6")
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAElement_Create failed.')

    
    ! Define the Solutes

    
    ! Define the Gases

    
    ! Define the Processes


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
    use physconst,          only: latice, latvap, cpair, cappa
    use time_manager,       only: get_curr_date, get_perp_date, get_curr_calday, is_perpetual

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
    type(physics_buffer_desc), pointer    :: pbuf(:)      !! physics buffer
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
  !! When the tracer hits at the surface at a time other than on its launch day,
  !! it will be removed from the model.
  !!
  !!  @version July-2009 
  !!  @author  Chuck Bardeen 
  subroutine CARMA_DiagnoseBulk(carma, cstate, cam_out, state, pbuf, ptend, icol, dt, rc, rliq, prec_str, snow_str, &
    prec_sed, snow_sed, tnd_qsnow, tnd_nsnow, re_ice)
    use shr_kind_mod,  only: r8 => shr_kind_r8
    use ppgrid,        only: pcols, pver
    use physics_types, only: physics_state
    use time_manager,  only: get_curr_date, get_perp_date, get_curr_calday, &
                             is_perpetual
    use camsrfexch,    only: cam_out_t

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
    
    real(r8)     :: calday                  ! current calendar day
    integer      :: yr                      ! year
    integer      :: mon                     ! month
    integer      :: day                     ! day of month
    integer      :: ncsec                   ! time of day (seconds)
    integer      :: doy                     ! day of year
    integer      :: elapsed                 ! days since launch
    
 
    ! Default return code.
    rc = RC_OK

    ! By default, do nothing. If diagnosed groups exist, this needs to be replaced by
    ! code to determine the bulk mass from the CARMA state.
    
    if (present(re_ice)) re_ice(:,:) = 0.0_f
    
    ! Determine the day of year.
    calday = get_curr_calday()
    if ( is_perpetual() ) then
      call get_perp_date(yr, mon, day, ncsec)
    else
      call get_curr_date(yr, mon, day, ncsec)
    end if
    doy = floor(calday)
    
    ! Any material that has made it to the surface from a previous day should be removed.
    elapsed = doy - carma_launch_doy
    
    if (elapsed > 1) then
      cstate%f_pc(pver, 1:min(NBIN,elapsed-1), :NELEM) = 0._f
    end if
    
    return
  end subroutine CARMA_DiagnoseBulk
 
 
  !! Calculates the emissions for CARMA aerosol particles.
  !!
  !! Emit particles after the specified launch day, with each bin being used
  !! for only one day. Each element is a different regions. Regions can be defined by:
  !!
  !!   latitude and longitude range
  !!   surface pressure
  !!   surface type (land, ocean, sea ice)
  !!
  !! The tracer is emitted as a constant column mass (kg/m2/s), so regions are weighted
  !! by their surface area. Mixed surface types have are scaled by the fraction of the
  !! surface type included in the region.
  !!
  !! A negative surface pressure means that the surface pressure must be less than or
  !! equal to the number specified, While a positive number means it must be greater
  !! than the specified value.
  !!
  !! One extra region is defined that is all of the areas excluded (via lat/lon) from
  !! all of the other regions (i.e. the rest of the world).
  !!
  !! NOTE: Launch days that wrap around are not currently supported
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  subroutine CARMA_EmitParticle(carma, ielem, ibin, icnst, dt, state, cam_in, tendency, surfaceFlux, rc)
    use shr_kind_mod,  only: r8 => shr_kind_r8
    use ppgrid,        only: pcols, pver
    use physics_types, only: physics_state
    use time_manager,  only: get_curr_date, get_perp_date, get_curr_calday, &
                             is_perpetual
    use camsrfexch,    only: cam_in_t
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

    real(r8)     :: lat(state%ncol)         ! latitude (degrees)
    real(r8)     :: lon(state%ncol)         ! longitude (degrees)
    integer      :: ncol                    ! number of columns in chunk
    integer      :: icol                    ! column index
    real(r8)     :: calday                  ! current calendar day
    integer      :: yr                      ! year
    integer      :: mon                     ! month
    integer      :: day                     ! day of month
    integer      :: ncsec                   ! time of day (seconds)
    integer      :: doy                     ! day of year
    integer      :: elapsed                 ! days since launch
    logical      :: doPS                    ! is pressure in correct range?
    logical      :: doRegion                ! are lat/lon in correct range?
    integer      :: i
    real(r8)     :: frac                    ! scaling fraction from land type

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

    elapsed = doy - carma_launch_doy

    ! Determine the latitude and longitude of each column.
    ncol = state%ncol

    surfaceFlux(:ncol)      = 0.0_f
    tendency(:ncol, :pver)  = 0.0_f

    ! Is this a day to launch more material?
    if ((elapsed + 1) == ibin) then

      lat = state%lat(:ncol) / DEG2RAD
      lon = state%lon(:ncol) / DEG2RAD

      do icol = 1, ncol

        ! Determine the region based upon latitude and longitude. The last region is
        ! defined to be rest of the world (i.e. all regions not in another region).
        doRegion = .False.
        
        if (ielem == NELEM) then
          doRegion = .True.
          
          do i = 1, NELEM-1
            if ((rgn_minlat(i) < lat(icol)) .and. (lat(icol) <= rgn_maxlat(i)) .and. &
                (rgn_minlon(i) < lon(icol)) .and. (lon(icol) <= rgn_maxlon(i))) then
              doRegion = .False.
            end if
          end do
        else
          if ((rgn_minlat(ielem) < lat(icol)) .and. (lat(icol) <= rgn_maxlat(ielem)) .and. &
              (rgn_minlon(ielem) < lon(icol)) .and. (lon(icol) <= rgn_maxlon(ielem))) then
            doRegion = .True.
          end if
        end if

        ! Check the surface pressure.
        doPS = .False.
        if (rgn_ps(ielem) == 0._f) then
          doPS = .True.
        else 
          if (rgn_ps(ielem) > 0._f) then
            if (state%ps(icol) > rgn_ps(ielem)) then
              doPS = .True.
            end if
          else
            if (state%ps(icol) <= abs(rgn_ps(ielem))) then
              doPS = .True.
            end if            
          end if
        end if
        
        ! Calculate the emission rate as a constant mass.
        if (doRegion .and. doPS) then

          ! A negative emission rate means to treat it as an mmr (kg/kg/s) and a
          ! postive value means to treat it as a column mass (kg/m2/s).
          if (carma_emission_rate > 0._f) then
            tendency(icol, pver) = carma_emission_rate / state%pdel(icol, pver) / gravit
          else
            ! For mmr, calculate a tendecy to keep the surface at that emitted value,
            ! rather than having a constant emission rate.
!            tendency(icol, pver) = -carma_emission_rate  
            tendency(icol, pver) = ((-carma_emission_rate * dt) - state%q(icol, pver, icnst)) / dt  
          end if
        end if

        ! Scale with the land/ocean fraction.
        frac = 0._f
        
        if (rgn_doLand(ielem)) then
          frac = frac + cam_in%landfrac(icol)
        end if

        if (rgn_doOcean(ielem)) then
          frac = frac + cam_in%ocnfrac(icol)
        end if

        if (rgn_doSeaIce(ielem)) then
          frac = frac + cam_in%icefrac(icol)
        end if

        tendency(icol, pver) = tendency(icol, pver) * frac
      end do
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
    use constituents, only : pcnst
    implicit none
    
    type(carma_type), intent(inout)    :: carma                 !! the carma object
    logical, intent(inout)             :: lq_carma(pcnst)       !! flags to indicate whether the constituent
                                                                !! could have a CARMA tendency
    integer, intent(out)               :: rc                    !! return code, negative indicates failure

    ! Default return code.
    rc = RC_OK
    
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

    ! Add initial condition here (default is 0.)
    q(:,:) = 0._f
    
    return
  end subroutine CARMA_InitializeParticle
  
  
  !!  Called after wet deposition has been performed. Allows the specific model to add
  !!  wet deposition of CARMA aerosols to the aerosols being communicated to the surface.
  !!
  !!  @version July-2011 
  !!  @author  Chuck Bardeen 
  subroutine CARMA_WetDeposition(carma, ielem, ibin, sflx, cam_out, state, rc)
    use camsrfexch, only: cam_out_t

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
    
    return
  end subroutine CARMA_WetDeposition 

end module
