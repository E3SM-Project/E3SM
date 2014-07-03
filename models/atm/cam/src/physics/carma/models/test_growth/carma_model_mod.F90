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
  integer, public, parameter      :: NGROUP   = 2               !! Number of particle groups
  integer, public, parameter      :: NELEM    = 3               !! Number of particle elements
  integer, public, parameter      :: NBIN     = 16              !! Number of particle bins
  integer, public, parameter      :: NSOLUTE  = 1               !! Number of particle solutes
  integer, public, parameter      :: NGAS     = 1               !! Number of gases
  
  ! These need to be defined, but are only used when the particles are radiatively active.
  integer, public, parameter      :: NMIE_RH  = 8               !! Number of relative humidities for mie calculations
  real(kind=f), public            :: mie_rh(NMIE_RH)
  
  ! Defines whether the groups should undergo deep convection in phase 1 or phase 2.
  ! Water vapor and cloud particles are convected in phase 1, while all other constituents
  ! are done in phase 2.
  logical, public                 :: is_convtran1(NGROUP) = .false.  !! Should the group be transported in the first phase?

  ! Define any particle compositions that are used. Each composition type
  ! should have a unique number.
  integer, public, parameter      :: I_H2SO4   = 1              !! sulfate aerosol composition
  integer, public, parameter      :: I_ICE     = 2              !! ice

  ! Define group, element, solute and gas indexes.
  integer, public, parameter      :: I_GRP_CRCN     = 1         !! sulfate aerosol
  integer, public, parameter      :: I_GRP_CRICE    = 2         !! ice

  integer, public, parameter      :: I_ELEM_CRCN    = 1         !! sulfate aerosol
  integer, public, parameter      :: I_ELEM_CRICE   = 2         !! ice
  integer, public, parameter      :: I_ELEM_CRCORE  = 3         !! sulfate core

  integer, public, parameter      :: I_SOL_CRH2SO4  = 1         !! H2SO4

  integer, public, parameter      :: I_GAS_H2O      = 1         !! water vapor

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
    real(kind=f), parameter            :: RHO_CN   = 2.65_f    ! dry density of sea salt particles (g/cm)
    real(kind=f), parameter            :: rmin_ice = 5.e-5_f   ! min radius for ice bins (cm)
    real(kind=f), parameter            :: rmin_cn  = 1.e-7_f   ! min radius for sulfate bins (cm)

    ! Default return code.
    rc = RC_OK
    
    ! Define the Groups
    !
    ! NOTE: If NWAVE > 0 then the group should have refractive indices defined.
    !
    ! NOTE: For CAM, the optional do_wetdep and do_drydep flags should be
    ! defined. If wetdep is defined, then the optional solubility factor
    ! should also be defined.
    
    ! Since these sulfates are prescribed, don't sediment them. This will save some
    ! processing time.
    call CARMAGROUP_Create(carma, I_GRP_CRCN, "Sulfate CN", rmin_cn, 4.0_f, I_SPHERE, 1._f, .false., &
                           rc, shortname="CRCN", cnsttype=I_CNSTTYPE_DIAGNOSTIC, &
                           do_vtran=.false.)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAGROUP_Create failed.')

    call CARMAGROUP_Create(carma, I_GRP_CRICE, "Ice", rmin_ice, 2.8_f, I_HEXAGON, 1._f / 6._f, .true., &
                           rc, shortname="CRICE", ifallrtn=I_FALLRTN_STD_SHAPE)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAGROUP_Create failed.')

   
    
    ! Define the Elements
    !
    ! NOTE: For CAM, the optional shortname needs to be provided for the group. These names
    ! should be 6 characters or less and without spaces.
    call CARMAELEMENT_Create(carma, I_ELEM_CRCN, I_GRP_CRCN, "Sulfate CN", RHO_CN, &
         I_INVOLATILE, I_H2SO4, rc, shortname="CRCN", isolute=I_SOL_CRH2SO4)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAElement_Create failed.')

    call CARMAELEMENT_Create(carma, I_ELEM_CRICE, I_GRP_CRICE, "Ice", RHO_I, &
         I_VOLATILE, I_ICE, rc, shortname="CRICE")
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAElement_Create failed.')
    
    call CARMAELEMENT_Create(carma, I_ELEM_CRCORE, I_GRP_CRICE, "Core Mass", RHO_CN, &
         I_COREMASS, I_H2SO4, rc, shortname="CRCORE", isolute=1)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAElement_Create failed.')

    
    ! Define the Solutes
    call CARMASOLUTE_Create(carma, I_SOL_CRH2SO4, "Sulfuric Acid", 2, 98._f, 1.38_f, rc, shortname="CRH2SO4")
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMASOLUTE_Create failed.')

    
    ! Define the Gases
    call CARMAGAS_Create(carma, I_GAS_H2O, "Water Vapor", WTMOL_H2O, I_VAPRTN_H2O_MURPHY2005, I_GCOMP_H2O, rc, shortname="Q")
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAGAS_Create failed.')

    
    ! Define the Processes
    call CARMA_AddGrowth(carma, I_ELEM_CRICE, I_GAS_H2O, rc)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_AddGrowth failed.')

    ! NOTE: For now, assume the latent heat for nucleation is the latent of of fusion of
    ! water, using the CAM constant (scaled from J/kg to erg/g).
    !
    ! NOTE: SInce the sulfates are not seen as part of the water/energy budget in CAM, don't
    ! include any latent heat from the freezing of the sulfate liquid. The latent heat of
    ! the gas associated with nucleation is accounted for.
    call CARMA_AddNucleation(carma, I_ELEM_CRCN, I_ELEM_CRCORE, &
         I_AERFREEZE + I_AF_KOOP_2000, 0._f, rc, igas=I_GAS_H2O, ievp2elem=I_ELEM_CRCN)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_AddNucleation failed.')

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
    type(physics_buffer_desc), pointer    :: pbuf(:)      !! physics buffer
    integer, intent(in)                   :: icol         !! column index
    real(r8), intent(in)                  :: dt           !! time step
    integer, intent(out)                  :: rc           !! return code, negative indicates failure
    real(r8), intent(in), optional        :: rliq(pcols)      !! vertical integral of liquid not yet in q(ixcldliq)
    real(r8), intent(inout), optional     :: prec_str(pcols)  !! [Total] sfc flux of precip from stratiform (m/s) 
    real(r8), intent(inout), optional     :: snow_str(pcols)  !! [Total] sfc flux of snow from stratiform (m/s)
    
    integer                               :: igroup             ! group index
    integer                               :: ielem              ! element index
    integer                               :: ibin               ! bin index
    
    ! Sulfate size distribution parameters
    real(r8), parameter                   :: n    = 100._r8     ! concentration (cm-3) 
    real(r8), parameter                   :: r0   = 2.5e-6_r8   ! mean radius (cm)
    real(r8), parameter                   :: rsig = 1.5_r8      ! distribution width
    
    real(r8)                              :: arg1(NBIN)
    real(r8)                              :: arg2(NBIN)
    real(r8)                              :: rhop(NBIN)         ! particle mass density (kg/m3)
    real(kind=f)                          :: rhoa_wet(pver)     ! air density (g/cm3)
    real(kind=f)                          :: r(NBIN)            ! bin mean radius
    real(kind=f)                          :: dr(NBIN)           ! bin radius width
    real(kind=f)                          :: rmass(NBIN)        ! bin mass
    real(r8)                              :: mmr(NBIN,pver)     ! elements mass mixing ratio

    ! Default return code.
    rc = RC_OK
    
    ! Get the air density.
    call CARMASTATE_GetState(cstate, rc, rhoa_wet=rhoa_wet)
    if (rc < RC_OK) call endrun('CARMA_DiagnoseBins::CARMASTATE_GetState failed.')

    ! Use a fixed sulfate size distribution. By doing this as a diagnostic group,
    ! the constituents for the sulfate bins do not need to be advected, which
    ! improves the speed of the model.
    igroup = 1
    ielem  = 1
    
    call CARMAGROUP_Get(carma, igroup, rc, r=r, dr=dr, rmass=rmass)
    if (rc < RC_OK) call endrun('CARMA_DiagnoseBins::CARMAGROUP_Get failed.')
    
    arg1(:) = n * dr(:) / (sqrt(2._f*PI) * r(:) * log(rsig))
    arg2(:) = -((log(r(:)) - log(r0))**2) / (2._f*(log(rsig))**2)

    rhop(:) = arg1(:) * exp(arg2(:)) * rmass(:) * 1e6_f / 1e3_f
      
    do ibin = 1, NBIN
      mmr(ibin, :) = rhop(ibin) / rhoa_wet(:)
      call CARMASTATE_SetBin(cstate, ielem, ibin, mmr(ibin, :), rc)
      if (rc < RC_OK) call endrun('CARMA_DiagnoseBins::CARMAGROUP_SetBin failed.')
    end do
    
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
    
    type(carma_type), intent(in)       :: carma                 !! the carma object
    logical, intent(inout)             :: lq_carma(pcnst)       !! flags to indicate whether the constituent
                                                                !! could have a CARMA tendency
    integer, intent(out)               :: rc                    !! return code, negative indicates failure

    ! Default return code.
    rc = RC_OK

    ! Add initialization here.

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

    ! Put a horizontally uniform layer of the smallest bin size
    ! in the model.
    if (ibin == 1) then
!      q(:, 1)        = 100e-9_r8    ! top
      q(:, plev/4)   = 100e-9_r8    ! 1/4
!      q(:, plev/2)   = 100e-9_r8    ! middle
!      q(:, 3*plev/4) = 100e-9_r8    ! 3/4
!      q(:, plev-1)   = 100e-9_r8    ! bottom
    end if 
    
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
