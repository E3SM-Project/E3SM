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
!! This file is a simple test case involving one group of dust particles and
!! 8 size bins. Optical properties are calculated, assuming a constant refractive
!! index of (1.55, 4e-3). The particles are not subject to particle swelling, but
!! do coagulate.
!!
!! @version May-2009 
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
  integer, public, parameter      :: NGROUP   = 3               !! Number of particle groups
  integer, public, parameter      :: NELEM    = 3               !! Number of particle elements
  integer, public, parameter      :: NBIN     = 22              !! Number of particle bins
  integer, public, parameter      :: NSOLUTE  = 0               !! Number of particle solutes
  integer, public, parameter      :: NGAS     = 0               !! Number of gases
  
  ! These need to be defined, but are only used when the particles are radiatively active.
  integer, public, parameter      :: NMIE_RH  = 8               !! Number of relative humidities for mie calculations
  real(kind=f), public            :: mie_rh(NMIE_RH)
  
  ! Defines whether the groups should undergo deep convection in phase 1 or phase 2.
  ! Water vapor and cloud particles are convected in phase 1, while all other constituents
  ! are done in phase 2.
  logical, public                 :: is_convtran1(NGROUP) = .false.  !! Should the group be transported in the first phase?

  ! Define any particle compositions that are used. Each composition type
  ! should have a unique number.
  integer, public, parameter      :: I_ICE   = 1               !! dust composition

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

    ! Default return code.
    rc = RC_OK

    ! Define the Groups
    !
    ! NOTE: If NWAVE > 0 then the group should have refractive indices defined.
    !
    ! NOTE: For CAM, the optional do_wetdep and do_drydep flags should be
    ! defined. If wetdep is defined, then the optional solubility factor
    ! should also be defined.
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddGroup failed.')
    call CARMAGROUP_Create(carma, 1, "Mass", rmin, 1.00001_f, I_SPHERE, 1._f, .true., &
                           rc, shortname="CRMASS")
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAGROUP_Create failed.')
    
    call CARMAGROUP_Create(carma, 2, "Number", rmin, 1.00001_f, I_SPHERE, 1._f, .true., &
                           rc, shortname="CRNUMB")
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAGROUP_Create failed.')
   
    call CARMAGROUP_Create(carma, 3, "Always", rmin, 1.000001_f, I_SPHERE, 1._f, .true., &
                           rc, shortname="CRALL")
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAGROUP_Create failed.')
   
    
    ! Define the Elements
    !
    ! NOTE: For CAM, the optional shortname needs to be provided for the group. These names
    ! should be 6 characters or less and without spaces.
    call CARMAELEMENT_Create(carma, 1, 1, "Mass", RHO_I, I_VOLATILE, I_ICE, rc, shortname="CRMASS")
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAElement_Create failed.')
    
    call CARMAELEMENT_Create(carma, 2, 2, "Number", RHO_I, I_VOLATILE, I_ICE, rc, shortname="CRNUMB")
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAElement_Create failed.')
    
    call CARMAELEMENT_Create(carma, 3, 3, "Always", RHO_I, I_VOLATILE, I_ICE, rc, shortname="CRALL")
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
    
    integer, parameter                   :: ielemMass  = 1
    integer, parameter                   :: ielemNumber  = 2
    integer, parameter                   :: ielemAlways  = 3
    real(r8), parameter                  :: theta_min = 320._r8
    real(r8), parameter                  :: theta_max = 400._r8
    real(r8), parameter                  :: dtheta = (theta_max - theta_min) / (NBIN - 2)
    real(r8), parameter                  :: startDay = 20._r8
    real(r8), parameter                  :: stopDay = 30._r8
!    real(r8), parameter                  :: startDay = 335._r8
!    real(r8), parameter                  :: stopDay = 345._r8
    real(r8), parameter                  :: mmrThreshold = 1e-8_r8
    integer                              :: k
    integer                              :: ibin
    real(r8)                             :: p(pver)
    real(r8)                             :: t(pver)
    real(r8)                             :: theta(pver)
    real(r8)                             :: mmr(pver)
    real(r8)                             :: cday

    ! Default return code.
    rc = RC_OK
    
    ! Only detrain condensate during certain days, so we get a pulse that we can then follow around.
    cday = get_curr_calday()
    
    if ((cday >= startDay) .and. (cday <= stopDay)) then
    
      call CARMASTATE_GetState(cstate, rc, t=t, p=p)
      if (rc < RC_OK) call endrun('CARMA_Detrain::CARMASTATE_GetState failed.')
      
      ! Calculate the potential temperature.
      theta(:) = t(:) * (1e5_f / p(:)) ** (cappa)       

      ! Put all of the detraining cloud water from convection into the large scale cloud.
      ! put detraining cloud water into liq and ice based on temperature partition
      do k = 1, pver
      
        ! Put the detrained mass into a bin sorted by potential temperature.
        ibin = ((theta(k) - theta_min) / dtheta) + 2
        
        if (theta(k) < theta_min) then
          ibin = 1
        else if (theta(k) > theta_max) then
          ibin = NBIN
        endif

        
        ! Mass based tracer reflects detrained mass.
        call CARMASTATE_GetBin(cstate, ielemMass, ibin, mmr, rc)
        if (rc < RC_OK) call endrun('CARMA_Detrain::CARMASTATE_GetBin failed.')
        
        mmr(k) = mmr(k) + dlf(icol, k) * dt
        
        call CARMASTATE_SetBin(cstate, ielemMass, ibin, mmr, rc)
        if (rc < RC_OK) call endrun('CARMA_Detrain::CARMASTATE_SetBin failed.')

        
        ! Event based tracer reflects the number of detrainment events, but apply a
        ! threshold so that not everything counts.
        if ((dlf(icol, k) * dt) >= mmrThreshold) then

          call CARMASTATE_GetBin(cstate, ielemNumber, ibin, mmr, rc)
          if (rc < RC_OK) call endrun('CARMA_Detrain::CARMASTATE_GetBin failed.')
        
          mmr(k) = mmr(k) + 1e-6_r8
        
          call CARMASTATE_SetBin(cstate, ielemNumber, ibin, mmr, rc)
          if (rc < RC_OK) call endrun('CARMA_Detrain::CARMASTATE_SetBin failed.')
        end if
        

        ! This group puts in tracer at all potential temperatures, just to see what happens
        ! independent of whether detrainment actually puts mass there.
        call CARMASTATE_GetBin(cstate, ielemAlways, ibin, mmr, rc)
        if (rc < RC_OK) call endrun('CARMA_Detrain::CARMASTATE_GetBin failed.')
      
        mmr(k) = mmr(k) + 1e-6_r8
      
        call CARMASTATE_SetBin(cstate, ielemAlways, ibin, mmr, rc)
        if (rc < RC_OK) call endrun('CARMA_Detrain::CARMASTATE_SetBin failed.')
 
      end do
    end if

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
 
    ! Default return code.
    rc = RC_OK

    ! By default, do nothing. If diagnosed groups exist, this needs to be replaced by
    ! code to determine the bulk mass from the CARMA state.
    
    if (present(re_ice)) re_ice(:,:) = 0.0_f
    
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
                             is_perpetual
    use camsrfexch,       only: cam_in_t

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

    integer      :: ncol                    ! number of columns in chunk
    integer      :: icol                    ! column index
    real(r8)     :: calday                  ! current calendar day
    integer      :: yr                      ! year
    integer      :: mon                     ! month
    integer      :: day                     ! day of month
    integer      :: ncsec                   ! time of day (seconds)
    integer      :: doy                     ! day of year

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

    ncol = state%ncol
    
    ! Add any surface flux here.
    surfaceFlux(:ncol) = 0.0_r8
    
    ! For emissions into the atmosphere, put the emission here.
    tendency(:ncol, :pver) = 0.0_r8
    
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
    
    return
  end subroutine CARMA_WetDeposition 

end module
