!! This CARMA model is for polar mesospheric clouds and meteor smoke aerosols and
!! is based upon Bardeen et al., JGR, 2010.
!!
!! This module defines several constants needed by CARMA, extends a couple of CARMA
!! interface methods:
!!
!!   - CARMA_DefineModel()
!!   - CARMA_EmitParticle()
!!   - CARMA_InitializeModel()
!!
!! and adds some local functions used to do sea salt emission:
!!
!!
!! NOTE: This model is still under development and is not intended to be released as
!! part of the standard CAM distribution. Please contact Charles Bardeen at
!! bardeenc@ucar.edu if you are interested in using or deriving off of this
!! work.
!!
!! @version Jan-2011
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

  use spmd_utils,     only: masterproc
  use shr_kind_mod,   only: r8 => shr_kind_r8
  use radconstants,   only: nswbands, nlwbands
  use abortutils,     only: endrun
  use physics_types,  only: physics_state, physics_ptend
  use ppgrid,         only: pcols, pver
  use physics_buffer, only: physics_buffer_desc

#if ( defined SPMD )
  use mpishorthand
#endif  

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
  integer, public, parameter      :: NELEM    = 5               !! Number of particle elements
  integer, public, parameter      :: NBIN     = 28              !! Number of particle bins
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
  integer, public, parameter      :: I_METEOR_SMOKE   = 1       !! meteor smoke
  integer, public, parameter      :: I_ICE            = 2       !! ice
  integer, public, parameter      :: I_H2SO4          = 3       !! sulfuric acid

  ! Define group, element, solute and gas indexes.
  integer, public, parameter      :: I_GRP_DUST     = 1         !! meteor smoke
  integer, public, parameter      :: I_GRP_CRICE    = 2         !! ice
  integer, public, parameter      :: I_GRP_SULFATE  = 3         !! sulfate aerosol

  integer, public, parameter      :: I_ELEM_DUST    = 1         !! meteor smoke
  integer, public, parameter      :: I_ELEM_CRICE   = 2         !! ice
  integer, public, parameter      :: I_ELEM_CRCORE  = 3         !! meteor smoke core in ice
  integer, public, parameter      :: I_ELEM_SULFATE = 4         !! sulfate aerosol
  integer, public, parameter      :: I_ELEM_SULCORE = 5         !! meteor smoke core in sulfate

  integer, public, parameter      :: I_SOL_CRH2SO4  = 1         !! sulfuric acid

  integer, public, parameter      :: I_GAS_H2O      = 1         !! water vapor
  integer, public, parameter      :: I_GAS_H2SO4    = 2         !! sulphuric acid

  real(kind=f), public, parameter :: WTMOL_H2SO4    = 98.078479_f    !! molecular weight of sulphuric acid

  ! These variables are all set during initialization and are used to calculate
  ! emission tendencies.
  integer                             :: carma_emis_nLevs        ! number of emission levels
  real(r8), allocatable, dimension(:) :: carma_emis_lev          ! emission levels (Pa)
  real(r8), allocatable, dimension(:) :: carma_emis_rate         ! emission rate lookup table (# cm-3 s-1)
  integer                             :: carma_emis_ilev_min     ! index of minimum level in table 
  integer                             :: carma_emis_ilev_max     ! index of maximum level in table 
  integer                             :: carma_emis_ilev_incr    ! index increment to increase level 
  real(r8)                            :: carma_emis_expected     ! Expected emission rate per column (kg/m2/s)

  integer                             :: carma_escale_nLats      ! number of emission scale latitudes
  integer                             :: carma_escale_nTimes     ! number of emission scale times
  integer                             :: carma_escale_nLTimes    ! number of emission scale local times
  real(r8), allocatable, dimension(:,:) :: carma_escale_grf      ! global relative flux
  real(r8), allocatable, dimension(:) :: carma_escale_lat        ! global relative flux latitudes
  real(r8), allocatable, dimension(:) :: carma_escale_lrf        ! locat time realtive flux
  real(r8), allocatable, dimension(:) :: carma_escale_ltime      ! local time relative flux times

  integer                             :: warren_nwave            ! number of wavelengths in file
  real(r8), allocatable, dimension(:) :: warren_wave             ! Warren & Brandt 2008, wavelengths
  real(r8), allocatable, dimension(:) :: warren_real             ! Warren & Brandt 2008, real part of m
  real(r8), allocatable, dimension(:) :: warren_imag             ! Warren & Brandt 2008, imag part of m

contains


  !! Defines all the CARMA components (groups, elements, solutes and gases) and process
  !! (coagulation, growth, nucleation) that will be part of the microphysical model.
  !!
  !!  @version May-2009 
  !!  @author  Chuck Bardeen 
  subroutine CARMA_DefineModel(carma, rc)
    use ioFileMod,    only: getfil
    use wrap_nf

    type(carma_type), intent(inout)    :: carma     !! the carma object
    integer, intent(out)               :: rc        !! return code, negative indicates failure
    
    ! Local variables
    real(kind=f), parameter            :: RHO_METEOR_SMOKE = 2.0_f      ! density of meteor smoke particles (g/cm)
    real(kind=f), parameter            :: rmin             = 2e-8_f     ! minimum radius (cm)
    real(kind=f), parameter            :: RHO_SULFATE      = 1.923_f    ! dry density of sulfate particles (g/cm3)
!  Set radius of smallest bin such that mass is that of 2 molecules of H2SO4:    
    real(kind=f), parameter            :: rmin_sulfate     = 3.43230298e-8_f  ! minimum radius (cm)
    real(kind=f), parameter            :: vmrat_sulfate    = 2.56_f     ! volume ratio (adjusted for 28 bins from 2.4 for 30 bins)


    integer                            :: i
    integer                            :: j
    real(kind=f)                       :: wave(NWAVE)               ! CAM band wavelength centers (cm)
    integer                            :: fid
    integer                            :: wave_did
    integer                            :: wave_vid
    integer                            :: real_vid
    integer                            :: imag_vid
    character(len=256)                 :: efile                     ! refractive index file name
    real(kind=f)                       :: interp
    complex(kind=f)                    :: refidx_ice(NWAVE)         ! the refractive index at each CAM wavelength    
    integer                            :: LUNOPRT
    logical                            :: do_print
    
    ! Default return code.
    rc = RC_OK
    
    call CARMA_Get(carma, rc, do_print=do_print, LUNOPRT=LUNOPRT, wave=wave)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_Get failed.')
    
    ! Report model specific configuration parameters.
    if (masterproc) then
      if (do_print) then
        write(LUNOPRT,*) ''
        write(LUNOPRT,*) 'CARMA ', trim(carma_model), ' specific settings :'
        write(LUNOPRT,*) '  carma_do_escale     = ', carma_do_escale
        write(LUNOPRT,*) '  carma_neutral_h2s04 = ', carma_neutral_h2so4
        write(LUNOPRT,*) '  carma_emis_file     = ', trim(carma_emis_file)
        write(LUNOPRT,*) '  carma_escale_file   = ', trim(carma_escale_file)
        write(LUNOPRT,*) '  carma_mice_file     = ', trim(carma_mice_file)
      end if
    end if
    
    
    ! Define the Groups
    !
    ! NOTE: For CAM, the optional do_wetdep and do_drydep flags should be
    ! defined. If wetdep is defined, then the optional solubility factor
    ! should also be defined.
    call CARMAGROUP_Create(carma, I_GRP_DUST, "meteor smoke", rmin, 2.0_f, I_SPHERE, 1._f, .false., &
                          rc, do_wetdep=.true., do_drydep=.true., solfac=0.3_f, &
                           scavcoef=0.1_f, shortname="DUST")
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddGroup failed.')                           

    ! Get the refractive index for ice as a function of wavelength for particle heating
    ! calculations.
    !
    ! NOTE: These values probably should be a band average, but for now just do band centers.
    
    ! Read the values in from Warren et al. 2008.
    if (carma_do_pheat) then
      if (masterproc) then 
      
        ! Open the netcdf file (read only)
        call getfil(carma_mice_file, efile, fid)
        if (do_print) write(LUNOPRT,*) 'carma_init(): Reading ice refractive indexes from ', efile
  
        call wrap_open(efile, 0, fid)
  
        ! Alocate the table arrays
        call wrap_inq_dimid(fid, "wavelength", wave_did)
        call wrap_inq_dimlen(fid, wave_did, warren_nwave)
      endif
      
#if ( defined SPMD )
        call mpibcast(warren_nwave, 1, mpiint, 0, mpicom)
#endif
  
        allocate(warren_wave(warren_nwave))
        allocate(warren_real(warren_nwave))
        allocate(warren_imag(warren_nwave))
  
        if (masterproc) then
          
          ! Read in the tables.
          call wrap_inq_varid(fid, 'wavelength', wave_vid)
          call wrap_get_var_realx(fid, wave_vid, warren_wave)
          warren_wave = warren_wave * 1e-4          ! um -> cm
  
          call wrap_inq_varid(fid, 'm_real', real_vid)
          call wrap_get_var_realx(fid, real_vid, warren_real)
  
          call wrap_inq_varid(fid, 'm_imag', imag_vid)
          call wrap_get_var_realx(fid, imag_vid, warren_imag)
  
          ! Close the file.
          call wrap_close(fid)
        end if
  
#if ( defined SPMD )
        call mpibcast(warren_wave,  warren_nwave, mpir8, 0, mpicom)
        call mpibcast(warren_real,  warren_nwave, mpir8, 0, mpicom)
        call mpibcast(warren_imag,  warren_nwave, mpir8, 0, mpicom)
#endif
      
      ! Interpolate the values.
      do i = 1, NWAVE
        do j = 1, warren_nwave
          if (wave(i) > warren_wave(j)) then
            if (j > 1) then
              interp = (wave(i) - warren_wave(j-1)) / (warren_wave(j) - warren_wave(j-1))
              refidx_ice(i) = cmplx(warren_real(j-1) + interp*(warren_real(j) - warren_real(j-1)), &
                   warren_imag(j-1) + interp*(warren_imag(j) - warren_imag(j-1)))
            else
              refidx_ice(i) = cmplx(warren_real(j), warren_imag(j))
            endif
            
            exit
          end if
        end do
      end do
    end if

    
    call CARMAGROUP_Create(carma, I_GRP_CRICE, "ice crystal", rmin, 2.2_f, I_SPHERE, 1._f, .true., &
                          rc, do_mie=carma_do_pheat, refidx=refidx_ice, shortname="CRICE")
    if (rc < 0) call endrun('CARMA_DefineModel::CARMAGROUP_Create failed.')

    ! NOTE: For CAM, the optional do_wetdep and do_drydep flags should be
    ! defined. If wetdep is defined, then the optional solubility factor
    ! should also be defined.
    ! solfac was formerly set to 0.3, changed to 1.0 because it seems physical.
    ! This change needs to be validated -MJM 12/1/2011
    !
    ! NOTE: Add flag to enable neutralization for sulfates?
    call CARMAGROUP_Create(carma, I_GRP_SULFATE, "sulfate", rmin_sulfate, vmrat_sulfate, I_SPHERE, 1._f, .false., &
                           rc, irhswell=I_WTPCT_H2SO4, do_wetdep=.true., do_drydep=.true., solfac=1.0_f, &
                           scavcoef=0.1_f, is_sulfate=.true., shortname="SULF")
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddGroup failed.')

    
    ! Define the Elements
    !
    ! NOTE: For CAM, the optional shortname needs to be provided for the group. These names
    ! should be 6 characters or less and without spaces.
    call CARMAELEMENT_Create(carma, I_ELEM_DUST, I_GRP_DUST, "meteor smoke", RHO_METEOR_SMOKE, &
         I_INVOLATILE, I_METEOR_SMOKE, rc, shortname="DUST")
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddElement failed.')

    call CARMAELEMENT_Create(carma, I_ELEM_CRICE, I_GRP_CRICE, "ice crystal", RHO_I, &
         I_VOLATILE, I_ICE, rc, shortname="CRICE")
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddElement failed.')

    call CARMAELEMENT_Create(carma, I_ELEM_CRCORE, I_GRP_CRICE, "ice core", RHO_METEOR_SMOKE, &
         I_COREMASS, I_METEOR_SMOKE, rc, shortname="CRCORE")
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddElement failed.')
    
    call CARMAELEMENT_Create(carma, I_ELEM_SULFATE, I_GRP_SULFATE, "sulfate", RHO_SULFATE, &
         I_VOLATILE, I_H2SO4, rc, shortname="SULF")
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddElement failed.')
    
    call CARMAELEMENT_Create(carma, I_ELEM_SULCORE, I_GRP_SULFATE, "sulfate core", RHO_METEOR_SMOKE, &
         I_COREMASS, I_METEOR_SMOKE, rc, shortname="SFCORE")
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddElement failed.')
    
    ! Define the Solutes
    !
    ! Should this be a RHO_SULFATE_WET of 1.38?
    
    
    ! Define the Gases
    call CARMAGAS_Create(carma, I_GAS_H2O, "Water Vapor", WTMOL_H2O, &
         I_VAPRTN_H2O_MURPHY2005, I_GCOMP_H2O, rc, shortname="Q", ds_threshold=0.2_f)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAGAS_Create failed.')
    
    call CARMAGAS_Create(carma, I_GAS_H2SO4, "Sulfuric Acid", WTMOL_H2SO4, I_VAPRTN_H2SO4_AYERS1980, &
                         I_GCOMP_H2SO4, rc, shortname = "H2SO4", ds_threshold=-0.2_f)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAGAS_Create failed.')

    
    ! Define the Processes
    call CARMA_AddCoagulation(carma, I_GRP_DUST, I_GRP_DUST, I_GRP_DUST, I_COLLEC_DATA, rc)
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddCoagulation failed.')

    call CARMA_AddNucleation(carma, I_ELEM_DUST, I_ELEM_CRCORE, I_HETNUC, 0._f, rc, igas=I_GAS_H2O, ievp2elem=I_ELEM_DUST)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_AddNucleation failed.')

    call CARMA_AddGrowth(carma, I_ELEM_CRICE, I_GAS_H2O, rc)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_AddGrowth failed.')

    call CARMA_AddCoagulation(carma, I_GRP_DUST, I_GRP_CRICE, I_GRP_CRICE, I_COLLEC_DATA, rc)
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddCoagulation failed.')

    ! Set H2SO4 to be the condensing gas, water vapor is assumed to be in equilibrium
    ! and will be used to define the wet particle radius.
    call CARMA_AddGrowth(carma, I_ELEM_SULFATE, I_GAS_H2SO4, rc)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_AddGrowth failed.')

    call CARMA_AddNucleation(carma, I_ELEM_SULFATE, I_ELEM_SULFATE, I_HOMNUC, 0._f, rc, igas=I_GAS_H2SO4)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_AddNucleation failed.')
    
    ! Also need nucleation with meteor smoke.
    call CARMA_AddNucleation(carma, I_ELEM_DUST, I_ELEM_SULCORE, I_HETNUCSULF, 0._f, rc, igas=I_GAS_H2SO4, ievp2elem=I_ELEM_DUST)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_AddNucleation failed.')

    call CARMA_AddCoagulation(carma, I_GRP_SULFATE, I_GRP_SULFATE, I_GRP_SULFATE, I_COLLEC_FUCHS, rc)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_AddCoagulation failed.')
    
    ! Dust-Sulfate Coagulation?
    call CARMA_AddCoagulation(carma, I_GRP_DUST, I_GRP_SULFATE, I_GRP_SULFATE, I_COLLEC_FUCHS, rc)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_AddCoagulation failed.')
    
    
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
    
    return
  end subroutine CARMA_DiagnoseBulk


  subroutine CARMA_EmitParticle(carma, ielem, ibin, icnst, dt, state, cam_in, tendency, surfaceFlux, rc)
    use shr_kind_mod,  only: r8 => shr_kind_r8
    use ppgrid,        only: pcols, pver
    use physics_types, only: physics_state
    use camsrfexch,       only: cam_in_t
    use time_manager,  only: get_curr_calday, is_perpetual, get_perp_date, get_curr_date
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
    
    integer                            :: ilat                  ! latitude index 
    integer                            :: iltime                ! local time index
    integer                            :: ncol                  ! number of columns in chunk
    integer                            :: icol                  ! column index
    integer                            :: igroup                ! the index of the carma aerosol group
    integer                            :: k                     ! vertical index
    integer                            :: ilev                  ! level index in emissions data
    character(len=32)                  :: shortname             ! the shortname of the group
    real(r8)                           :: r(NBIN)               ! bin center
    real(r8)                           :: dr(NBIN)              ! bin width
    real(r8)                           :: rmass(NBIN)           ! bin mass
    real(r8)                           :: pressure              ! pressure (Pa)
    real(r8)                           :: thickness             ! layer thickness (m)
    real(r8)                           :: rate                  ! emission rate (#/cm-3/s)
    real(r8)                           :: massflux              ! emission mass flux (kg/m2/s)
    real(r8)                           :: columnMass            ! mass of the total column (kg/m2/s)
    real(r8)                           :: scale                 ! scaling factor to conserve the expected mass
    real(r8)                           :: rfScale(pcols)        ! scaling factor from global and local relative flux

    real(r8)                           :: calday                ! current calendar day
    integer                            :: yr, mon, day, ncsec, doy
    integer                            :: ncdate
    real(r8)                           :: ltime                 ! local time

    
    ! Default return code.
    rc = RC_OK

    ! Get the current date and time.
    calday = get_curr_calday()
    if ( is_perpetual() ) then
      call get_perp_date(yr, mon, day, ncsec)
    else
      call get_curr_date(yr, mon, day, ncsec)
    end if
    doy = floor(calday)

    ! NOTE: The global relative flux file is based upon a noleap calendar, so don't
    ! let the doy get too big.
    doy = min(365, doy)

    ncol = state%ncol

    ! Add any surface flux here.
    surfaceFlux(:ncol) = 0.0_r8
    
    ! For emissions into the atmosphere, put the emission here.
    !
    ! NOTE: Do not set tendency to be the surface flux. Surface source is put in to
    ! the bottom layer by vertical diffusion. See vertical_solver module, line 355.            
    tendency(:ncol, :pver) = 0.0_r8


    ! Only do emission for the first bin of the meteor smoke group.
    call CARMAELEMENT_GET(carma, ielem, rc, igroup=igroup)
    if (RC < RC_ERROR) return
    
    call CARMAGROUP_GET(carma, igroup, rc, shortname=shortname, r=r, dr=dr, rmass=rmass)
    if (RC < RC_ERROR) return
    
    ! For meteoritic dust, the source from the smoke only goes into the
    ! smallest bin (~1.3 nm). The depth that the micrometeorite penetrates
    ! is proportional to the pressure, so the emission is a function of
    ! pressure. 
    if ((shortname .eq. "DUST") .and. (ibin .eq. 1)) then

      ! Change this to a scaling if appropriate.
      rfScale = 1.0_r8

      ! Set tendencies for any sources or sinks in the atmosphere.
      do k = 1, pver
        do icol = 1, ncol
      
          pressure = state%pmid(icol, k)
          
          ! This is roughly a log-normal approximation to the production
          ! rate, but only applies from about 70 to 110 km.
          !
          ! NOTE: Based upon US Standard Atmosphere 1976.
          if ((pressure >= carma_emis_lev(carma_emis_ilev_min)) .and. &
              (pressure <= carma_emis_lev(carma_emis_ilev_max))) then

            ! The rates are in terms of # cm-3 s-1, but were really derived
            ! from the mass flux of meteoritic dust. Since we are using a
            ! size different that 1.3 nm for the smallest bin, scale the
            ! number appropriately.
            !
            ! The values are in a lookup table, so find the two numbers
            ! surrounding the pressure and do a linear interpolation on the
            ! rate. This linear search is kind of expensive, particularly if
            ! there are a lot of points.
            ! 
            ! NOTE: The tendency is on a mass mixing ratio (kg/kg/s)
            do ilev = carma_emis_ilev_min, (carma_emis_ilev_max - carma_emis_ilev_incr), carma_emis_ilev_incr
              if ((pressure >= carma_emis_lev(ilev)) .and. (pressure <= carma_emis_lev(ilev+carma_emis_ilev_incr))) then
                rate = carma_emis_rate(ilev)
                
                if (pressure > carma_emis_lev(ilev)) then
                  rate = rate + &
                    ((carma_emis_rate(ilev+carma_emis_ilev_incr) - carma_emis_rate(ilev)) / &
                    (carma_emis_lev(ilev+carma_emis_ilev_incr) - carma_emis_lev(ilev))) * &
                    (pressure - carma_emis_lev(ilev))
                end if
                
                rate = rate * (((1.3e-7_r8)**3) / (r(ibin)**3))
                exit
              end if
            end do
            
            ! Calculate the mass flux in terms of kg/m3/s
            massflux = (rate * rmass(ibin) * 1.0e-3_r8 * 1.0e6_r8)
             
            if (carma_do_escale) then
            
              ! Global Scaling
              !
              ! Interpolate the global scale by latitude.
              !
              ! NOTE: It would be better to interpolate the  table once in init to the
              ! latitude structure and just look up by index.
              !
              ! NOTE: The latitudes have some small significant digits at the end, which makes
              ! exact comparisons to latitude values fail.
              do ilat = 1, carma_escale_nLats
                if ((state%lat(icol) / DEG2RAD) <= carma_escale_lat(ilat)) then
                  if (abs((state%lat(icol) / DEG2RAD) - carma_escale_lat(ilat)) <= 0.00001_r8) then
                    rfScale(icol) = carma_escale_grf(ilat, doy)
                  else
                    rfScale(icol) = carma_escale_grf(ilat-1, doy) + &
                    (((state%lat(icol) / DEG2RAD) - carma_escale_lat(ilat-1)) / &
                    (carma_escale_lat(ilat) - carma_escale_lat(ilat-1))) * &
                    (carma_escale_grf(ilat, doy) - carma_escale_grf(ilat-1, doy))
                  endif
                  exit
                end if
              end do

              if (abs((state%lat(icol) / DEG2RAD) - 90.0) <= 0.00001_r8) then
                 rfScale(icol) = carma_escale_grf(carma_escale_nLats, doy)
              end if
              
              ! Local Time Scaling
              !
              ! Interpolate the local scale by local time.
              ltime = abs((ncsec / 3600._r8) + (24._r8 * (state%lon(icol) / DEG2RAD) / 360._r8))
              if (ltime > 24._r8) then
                ltime = ltime - 24._r8
              end if

              do iltime = 1, carma_escale_nLTimes
                if (ltime <= carma_escale_ltime(iltime)) then
                  if (abs(ltime - carma_escale_ltime(iltime)) <= 0.00001_r8) then
                    rfScale(icol) = rfScale(icol) * carma_escale_lrf(iltime)
                  else
                    rfScale(icol) = rfScale(icol) * (carma_escale_lrf(iltime-1) + &
                    ((iltime - carma_escale_ltime(iltime-1)) / (carma_escale_ltime(iltime) - carma_escale_ltime(iltime-1))) * &
                    (carma_escale_lrf(iltime) - carma_escale_lrf(iltime-1)))
                  endif
                  exit
                end if
              end do
            endif
            
            ! Convert the mass flux to a tendency on the mass mixing ratio.
            thickness = state%zi(icol, k) - state%zi(icol, k+1)
            tendency(icol, k) = (massflux * thickness) / (state%pdel(icol, k) / gravit)        
          end if
        enddo
      enddo

      ! Scale the columns to keep the total mass influx in the column a
      ! constant.
      do icol = 1, ncol
        columnMass = sum(tendency(icol, :) * (state%pdel(icol, :) / gravit))

        ! Protect against divide-by-zero (but not overflow).
        if (columnMass /= 0._r8) then
           scale = carma_emis_expected / columnMass
        else
           scale = 0._r8
        end if

        ! Also apply the relative flux scaling. This needs to be done after
        ! the normalization
        tendency(icol, :) = tendency(icol, :) * scale * rfScale(icol)
      end do
    end if
    
    return
  end subroutine CARMA_EmitParticle


  !! Allows the model to perform its own initialization in addition to what is done
  !! by default in CARMA_init.
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  subroutine CARMA_InitializeModel(carma, lq_carma, rc)
    use constituents, only: pcnst
    use ioFileMod,    only: getfil
    use wrap_nf

    implicit none

    type(carma_type), intent(in)       :: carma                 !! the carma object
    logical, intent(inout)             :: lq_carma(pcnst)       !! flags to indicate whether the constituent
                                                                !! could have a CARMA tendency
    integer, intent(out)               :: rc                    !! return code, negative indicates failure

    integer                            :: ilev                  ! level index
    integer                            :: fid                   ! file id
    integer                            :: lev_did               ! level dimension id
    integer                            :: lev_vid               ! level variable id
    integer                            :: rate_vid              ! rate variable
    integer                            :: tmp
    integer                            :: lat_did               ! latitude dimension id
    integer                            :: ltime_did             ! local time dimension id
    integer                            :: time_did              ! time
    integer                            :: lat_vid               ! latitude variable id
    integer                            :: lrf_vid               ! local relative flux variable id
    integer                            :: grf_vid               ! global relative flux variable id
    integer                            :: ltime_vid             ! local time variable id
    character(len=256)                 :: efile                 ! emission file name

    integer                            :: LUNOPRT               ! logical unit number for output
    logical                            :: do_print              ! do print output?

    ! Default return code.
    rc = RC_OK

    ! Add initialization here.
    call CARMA_Get(carma, rc, do_print=do_print, LUNOPRT=LUNOPRT)
    if (rc < 0) call endrun("CARMA_InitializeModel: CARMA_Get failed.")
    
    ! Initialize the emissions rate table.
    if (carma_do_emission) then 
      if (masterproc) then

        ! Open the netcdf file (read only)
        call getfil(carma_emis_file, efile, fid)
        if (do_print) write(LUNOPRT,*) 'carma_init(): Reading particle emission rates from ', efile

        call wrap_open(efile, 0, fid)

        ! Alocate the table arrays
        call wrap_inq_dimid(fid, "lev", lev_did)
        call wrap_inq_dimlen(fid, lev_did, carma_emis_nLevs)
      endif
    
#if ( defined SPMD )
      call mpibcast(carma_emis_nLevs, 1, mpiint, 0, mpicom)
#endif

      allocate(carma_emis_lev(carma_emis_nLevs))
      allocate(carma_emis_rate(carma_emis_nLevs))

      if (masterproc) then
        ! Read in the tables.
        call wrap_inq_varid(fid, 'MSMOKE', rate_vid)
        call wrap_get_var_realx(fid, rate_vid, carma_emis_rate)

        call wrap_inq_varid(fid, 'lev', lev_vid)
        call wrap_get_var_realx(fid, lev_vid, carma_emis_lev)

        ! Close the file.
        call wrap_close(fid)

        ! Find out where the bounds of the table are and in what order
        ! the pressures levels are in.
        carma_emis_ilev_min = 1
        carma_emis_ilev_max = carma_emis_nLevs

        do ilev = 1, carma_emis_nLevs
          if (carma_emis_rate(ilev) <= 0.0) then
            carma_emis_ilev_min  = ilev + 1
          else
            exit  
          endif
        end do

        do ilev = carma_emis_nLevs, 1, -1
          if (carma_emis_rate(ilev) <= 0.0) then
            carma_emis_ilev_max  = ilev - 1
          else
            exit  
          endif
        end do

        if (carma_emis_lev(carma_emis_ilev_min) < carma_emis_lev(carma_emis_ilev_max)) then
          carma_emis_ilev_incr = 1
        else
          carma_emis_ilev_incr = -1
          tmp = carma_emis_ilev_min
          carma_emis_ilev_min = carma_emis_ilev_max
          carma_emis_iLev_max = tmp 
        endif

        if (do_print) write(LUNOPRT,*) ''
        if (do_print) write(LUNOPRT,*) 'carma_init(): carma_emis_nLevs     = ', carma_emis_nLevs
        if (do_print) write(LUNOPRT,*) 'carma_init(): carma_emis_ilev_min  = ', carma_emis_ilev_min 
        if (do_print) write(LUNOPRT,*) 'carma_init(): carma_emis_ilev_max  = ', carma_emis_ilev_max 
        if (do_print) write(LUNOPRT,*) 'carma_init(): carma_emis_ilev_incr = ', carma_emis_ilev_incr 
        if (do_print) write(LUNOPRT,*) ''
        
        if (do_print) write(LUNOPRT,*) 'level, pressure (Pa), emission rate (# cm-3 sec-1)'
        do ilev = carma_emis_ilev_min, carma_emis_ilev_max, carma_emis_ilev_incr
          if (do_print) write(LUNOPRT,*) ilev, carma_emis_lev(ilev), carma_emis_rate(ilev)
        enddo
        
        if (do_print) write(LUNOPRT, *) 'carma_init(): Total Emission = ', carma_emis_total, ' (kt/yr)'
        carma_emis_expected = ((carma_emis_total * 1e6_r8) / (3600.0_r8 * 24.0_r8 * 365.0_r8)) / &
             (4.0_r8 * PI * ((REARTH / 100._r8) ** 2))
        if (do_print) write(LUNOPRT,*) 'carma_init(): Done with emission table.'

      endif

#if ( defined SPMD )
      call mpibcast(carma_emis_lev,  carma_emis_nLevs, mpir8, 0, mpicom)
      call mpibcast(carma_emis_rate, carma_emis_nLevs, mpir8, 0, mpicom)

      call mpibcast(carma_emis_expected,  1, mpir8,  0, mpicom)

      call mpibcast(carma_emis_ilev_min,  1, mpiint, 0, mpicom)
      call mpibcast(carma_emis_ilev_max,  1, mpiint, 0, mpicom)
      call mpibcast(carma_emis_ilev_incr, 1, mpiint, 0, mpicom)
#endif

    endif
    

    ! Initialize the emissions scaling table.
    if (carma_do_escale) then 
      if (masterproc) then

        ! Open the netcdf file (read only)
        call getfil(carma_escale_file, efile, fid)
        if (do_print) write(LUNOPRT,*) 'carma_init(): Reading particle emission scaling from ', efile

        call wrap_open(efile, 0, fid)

        ! Alocate the table arrays
        call wrap_inq_dimid(fid, "lat", lat_did)
        call wrap_inq_dimlen(fid, lat_did, carma_escale_nLats)

        call wrap_inq_dimid(fid, "time", time_did)
        call wrap_inq_dimlen(fid, time_did, carma_escale_nTimes)
        
        ! There should be one time for each day of the year, so
        ! quit if it isn't correct.
        if (carma_escale_nTimes .ne. 365) then
          call endrun("CARMA_InitializeModel: Emission scaling file should have entries for 365 days, but doesn't.")
        endif
        
        call wrap_inq_dimid(fid, "ltime", ltime_did)
        call wrap_inq_dimlen(fid, ltime_did, carma_escale_nLTimes)
     endif
    
#if ( defined SPMD )
      call mpibcast(carma_escale_nLats,   1, mpiint, 0, mpicom)
      call mpibcast(carma_escale_nTimes,  1, mpiint, 0, mpicom)
      call mpibcast(carma_escale_nLTimes, 1, mpiint, 0, mpicom)
#endif

      allocate(carma_escale_lat(carma_escale_nLats))
      allocate(carma_escale_grf(carma_escale_nLats, carma_escale_nTimes))
      allocate(carma_escale_ltime(carma_escale_nLTimes))
      allocate(carma_escale_lrf(carma_escale_nLTimes))

      if (masterproc) then
        ! Read in the tables.
        call wrap_inq_varid(fid, 'SGRF', grf_vid)
        tmp = nf90_get_var (fid, grf_vid, carma_escale_grf)
        if (tmp/=NF90_NOERR) then
           write(iulog,*) 'CARMA_InitializeModel: error reading varid =', grf_vid
           call handle_error (tmp)
        end if

        call wrap_inq_varid(fid, 'lat', lat_vid)
        call wrap_get_var_realx(fid, lat_vid, carma_escale_lat)

        call wrap_inq_varid(fid, 'SLRF', lrf_vid)
        call wrap_get_var_realx(fid, lrf_vid, carma_escale_lrf)

        call wrap_inq_varid(fid, 'ltime', ltime_vid)
        call wrap_get_var_realx(fid, ltime_vid, carma_escale_ltime)
        
        ! Close the file.
        call wrap_close(fid)

        if (do_print) write(LUNOPRT,*) ''
        if (do_print) write(LUNOPRT,*) 'carma_init(): carma_escale_nLats   = ', carma_escale_nLats
        if (do_print) write(LUNOPRT,*) 'carma_init(): carma_escale_nTimes  = ', carma_escale_nTimes
        if (do_print) write(LUNOPRT,*) 'carma_init(): carma_escale_nLTimes = ', carma_escale_nLTimes
        if (do_print) write(LUNOPRT,*) ''
        
        if (do_print) write(LUNOPRT,*) 'carma_init(): Done with emission scaling tables.'

      endif

#if ( defined SPMD )
      call mpibcast(carma_escale_lat,   carma_escale_nLats, mpir8, 0, mpicom)
      call mpibcast(carma_escale_grf,   carma_escale_nLats*carma_escale_nTimes, mpir8, 0, mpicom)
      call mpibcast(carma_escale_ltime, carma_escale_nLTimes, mpir8, 0, mpicom)
      call mpibcast(carma_escale_lrf,   carma_escale_nLTimes, mpir8, 0, mpicom)
#endif

    endif
    
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
