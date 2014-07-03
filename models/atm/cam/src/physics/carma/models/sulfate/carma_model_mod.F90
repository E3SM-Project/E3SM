!! This CARMA model is for sulfate aerosols and is based upon work done by Mike Mills
!! and Jason English, which is described in English et al. 2011.
!!
!! These aerosols are not currently radiatively active and do not replace the sulfate
!! aerosols in CAM; however, this is something that could be done in the future.
!!
!! This module defines several constants needed by CARMA, extends a couple of CARMA
!! interface methods:
!!
!!   - CARMA_DefineModel()
!!   - CARMA_EmitParticle()
!!
!! @version Dec-2010
!! @author  Tianyi Fan, Chuck Bardeen 
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

  use spmd_utils,     only: masterproc
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
  integer, public, parameter      :: NGROUP   = 1               !! Number of particle groups
  integer, public, parameter      :: NELEM    = 1               !! Number of particle elements
  integer, public, parameter      :: NBIN     = 30              !! Number of particle bins
  integer, public, parameter      :: NSOLUTE  = 0               !! Number of particle solutes
  integer, public, parameter      :: NGAS     = 2               !! Number of gases

  ! These need to be defined, but are only used when the particles are radiatively active.
  integer, public, parameter      :: NMIE_RH  = 8               !! Number of relative humidities for mie calculations
  real(kind=f), public            :: mie_rh(NMIE_RH)

  ! Defines whether the groups should undergo deep convection in phase 1 or phase 2.
  ! Water vapor and cloud particles are convected in phase 1, while all other constituents
  ! are done in phase 2.
  logical, public                 :: is_convtran1(NGROUP) = .false.  !! Should the group be transported in the first phase?

  ! Define any particle compositions that are used. Each composition type
  ! should have a unique number.
  integer, public, parameter      :: I_H2SO4   = 1               !! sulfate aerosol composition
  integer, public, parameter      :: I_WATER   = 2               !! water
    
  ! Define group, element, solute and gas indexes.
  integer, public, parameter      :: I_GRP_SULFATE  = 1              !! sulfate aerosol

  integer, public, parameter      :: I_ELEM_SULFATE = 1              !! sulfate aerosol

  integer, public, parameter      :: I_GAS_H2O      = 1              !! water vapor
  integer, public, parameter      :: I_GAS_H2SO4    = 2              !! sulphuric acid

  real(kind=f), public, parameter :: WTMOL_H2SO4    = 98.078479_f    !! molecular weight of sulphuric acid

  ! Physics buffer index for sulfate surface area density
  integer                         :: ipbuf4sad

contains


  !! Defines all the CARMA components (groups, elements, solutes and gases) and process
  !! (coagulation, growth, nucleation) that will be part of the microphysical model.
  !!
  !!  @version May-2009 
  !!  @author  Chuck Bardeen 
  subroutine CARMA_DefineModel(carma, rc)
    use physics_buffer, only: pbuf_add_field, dtype_r8

    type(carma_type), intent(inout)    :: carma     !! the carma object
    integer, intent(out)               :: rc        !! return code, negative indicates failure
    
    ! Local variables
    real(kind=f), parameter            :: RHO_SULFATE = 1.923_f    ! dry density of sulfate particles (g/cm3)
!  Set radius of smallest bin such that mass is that of 2 molecules of H2SO4:    
    real(kind=f), parameter            :: rmin        = 3.43230298e-8_f  ! minimum radius (cm)
    real(kind=f), parameter            :: vmrat       = 2.4_f    ! volume ratio
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
        write(LUNOPRT,*) '  carma_do_hetchem      = ', carma_do_hetchem
      end if
    end if

    ! Define the Groups
    !
    ! NOTE: For CAM, the optional do_wetdep and do_drydep flags should be
    ! defined. If wetdep is defined, then the optional solubility factor
    ! should also be defined.
    ! solfac was formerly set to 0.3, changed to 1.0 because it seems physical.
    ! This change needs to be validated -MJM 12/1/2011
    call CARMAGROUP_Create(carma, I_GRP_SULFATE, "sulfate", rmin, vmrat, I_SPHERE, 1._f, .false., &
                           rc, irhswell=I_WTPCT_H2SO4, do_wetdep=.true., do_drydep=.true., solfac=1.0_f, &
                           scavcoef=0.1_f, is_sulfate=.true., shortname="PURSUL")
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddGroup failed.')

    
    ! Define the Elements
    !
    ! NOTE: For CAM, the optional shortname needs to be provided for the group. These names
    ! should be 6 characters or less and without spaces.
    call CARMAELEMENT_Create(carma, I_ELEM_SULFATE, I_GRP_SULFATE, "Sulfate", RHO_SULFATE, &
         I_VOLATILE, I_H2SO4, rc, shortname="PURSUL")
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddElement failed.')
        
    ! Define the Solutes
    
    
    ! Define the Gases
    call CARMAGAS_Create(carma, I_GAS_H2O, "Water Vapor", WTMOL_H2O, I_VAPRTN_H2O_MURPHY2005, I_GCOMP_H2O, &
                         rc, shortname = "Q", ds_threshold=-0.2_f)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAGAS_Create failed.')
    
    call CARMAGAS_Create(carma, I_GAS_H2SO4, "Sulfuric Acid", WTMOL_H2SO4, I_VAPRTN_H2SO4_AYERS1980, &
                         I_GCOMP_H2SO4, rc, shortname = "H2SO4", ds_threshold=-0.2_f)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAGAS_Create failed.')
    
    ! Define the Processes

    ! Set H2SO4 to be the condensing gas, water vapor is assumed to be in equilibrium
    ! and will be used to define the wet particle radius.
    call CARMA_AddGrowth(carma, I_ELEM_SULFATE, I_GAS_H2SO4, rc)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_AddGrowth failed.')

    call CARMA_AddNucleation(carma, I_ELEM_SULFATE, I_ELEM_SULFATE, I_HOMNUC, 0._f, rc, igas=I_GAS_H2SO4)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_AddNucleation failed.')

    call CARMA_AddCoagulation(carma, I_GRP_SULFATE, I_GRP_SULFATE, I_GRP_SULFATE, I_COLLEC_FUCHS, rc)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_AddCoagulation failed.')

    call pbuf_add_field('SADSULF', 'global', dtype_r8, (/pcols, pver/), ipbuf4sad)

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
    use camsrfexch,    only: cam_out_t
    use physics_buffer, only: pbuf_get_field

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

    ! Local variables
    real(r8)                             :: newstate(cstate%f_NZ)
    real(r8)                             :: numberDensity(cstate%f_NZ)
    real(r8)                             :: ad(cstate%f_NZ)       ! stratospheric aerosol surface area density (cm2/cm3)
    real(r8)                             :: r_wet(cstate%f_NZ)    ! Sulfate aerosol bin wet radius
    real(r8), pointer, dimension(:,:)    :: sadsulf_ptr           ! Sulfate surface area density pointer
    integer                              :: ibin, igroup

    ! Default return code.
    rc = RC_OK

    call CARMAELEMENT_Get(carma, I_ELEM_SULFATE, rc, igroup=igroup)
    if (rc < 0) call endrun('CARMA_DiagnoseBulk::CARMAELEMENT_Get failed.')

    ad(:)  = 0.0_r8

    do ibin = 1, NBIN
      call CARMASTATE_GetBin(cstate, I_ELEM_SULFATE, ibin, newstate(:), rc, &
                             numberDensity=numberDensity, r_wet=r_wet)
      if (rc < 0) call endrun('CARMA_DiagnoseBulk::CARMASTATE_GetBin failed.')

      ! Calculate the total densities.
      !
      ! NOTE: Calculate AD in cm2/cm3.
      if (numberDensity(1) /= CAM_FILL) then
        ad(:)  = ad(:)  + numberDensity(:) * 4.0_r8 * PI * (r_wet(:)**2)
      end if
    end do
    
    call pbuf_get_field(pbuf, ipbuf4sad, sadsulf_ptr)
    sadsulf_ptr(icol, :cstate%f_NZ) = ad(:cstate%f_NZ)


  end subroutine CARMA_DiagnoseBulk


  !! Calculates the emissions for CARMA aerosol particles. By default, there is no
  !! emission, but this routine can be overridden for models that wish to have
  !! an aerosol emission.
  !!
  !! @author  Tianyi Fan, Chuck Bardeen
  !! @version Dec-2010
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
    integer,  intent(out)              :: rc                    !! return code, negative indicates failure
        
    ! Default return code.
    rc = RC_OK   
             
    ! Add any surface flux here.
    surfaceFlux = 0._r8
    
    ! For emissions into the atmosphere, put the emission here.
    tendency = 0._r8

    return
  end subroutine CARMA_EmitParticle


  !! Allows the model to perform its own initialization in addition to what is done
  !! by default in CARMA_init.
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
    !
    ! NOTE: Initialized to 0. by the caller, so nothing needs to be done.

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
