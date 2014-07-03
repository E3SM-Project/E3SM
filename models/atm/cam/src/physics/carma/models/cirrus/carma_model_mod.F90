!! This module is used to define a particular CARMA microphysical model. For 
!! simple cases, this may be the only code that needs to be modified. This module
!! defines several constants and has the following methods:
!!
!!   - CARMA_DiagnoseBins()
!!   - CARMA_DiagnoseBulk()
!!   - CARMA_DefineModel()
!!   - CARMA_Detrain()
!!   - CARMA_EmitParticle()
!!   - CARMA_InitializeModel()
!!   - CARMA_InitializeParticle()
!!
!! These methods define the microphysical model, the particle emissions and
!! the initial conditions of the particles. For diagnostic groups, there are
!! also routines that diagnose the mass in the bins of that group from the
!! parent model's state inforamtion and that calculate the tendency on the
!! parent model's state based upon changes in the bins.
!!
!! This cirrus cloud model allows CARMA bin microphysics to do the ice microphysics
!! while MG does the liquid microphysics. The MG microphysics here should not update
!! CLDICE or NUMICE, since those values will not be reflected in the CARMA ice
!! bins, which are the true state variables for ice. In this situation, CLDICE and
!! NUMICE are merely diagnostic variables available as input to the rest of CAM.
!!
!! The CARMA microphysics will run before MG and will handle:
!!   - Detrainment (liquid and ice)
!!   - Homogeneous ice nucleation (currently with prescribed sulfates)
!!   - Heterogeneous ice nucleation (future)
!!   - Bergeron process
!!   - Melting of detrained ice
!!   - Freezing of cloud drops
!!   - Autoconversion (ice -> snow)
!!   - Variable ice density (function of particle size)
!!   - In-cloud values (dividing by cloud fraction)
!!
!! Some potential issues that are not currently handled by CARMA:
!!   - collection of ice by snow
!!   - aggregation of ice
!!   - sub-grid vertical velocity for CARMA
!!   - Goff & Gratch vs. Murphy & Koop vapor pressures
!!   - Radiation using CARMA size distribution (each bin as tracer)
!!   - Hallet-Mossop Process
!!
!! The following variables will have been set by CARMA:
!!   - (S) CLDICE, (S) NUMICE
!!   - (S) CLDLIQ, (S) NUMLIQ
!!   - (S) T
!!   - (P) TNDQSNOW,  (P) TNDNSNOW
!!   - (P) REICE
!!
!! Varaibles with an S will be in the physics_state and variables with a P are
!! parameters passed into the MG microphysics.
!!
!! The module carma_intr defines a few flags that indicate what portion of the
!! cloud microphysics is handled by CARMA:
!!
!!   - carma_do_cldice  - CARMA does ice clouds
!!   - carma_do_cldliq  - CARMA does liquid clouds
!!
!!---------------------------------------------------------------------------------


!! Each realization of CARMA microphysics has its own version of this file.
!!
!! This model replaces the ice microphysics from the MG two-moment scheme with
!! a CARMA bin microphysics representation of the ice. The purpose of this
!! model is to provide a more detail description of the thin cirrus clouds that
!! form in the TTL and to investigate the impact of these clouds on radiative
!! forcing, troposphere-to-stratosphere transport, and control of water vapor
!! in the UT/LS.
!!
!! @version July-2009 
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
  use physics_buffer, only: physics_buffer_desc, pbuf_old_tim_idx, pbuf_get_field, pbuf_get_index
  use physconst,      only: gravit
  
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
  integer, public, parameter      :: NGROUP   = 4               !! Number of particle groups
  integer, public, parameter      :: NELEM    = 5               !! Number of particle elements
  integer, public, parameter      :: NBIN     = 28              !! Number of particle bins
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
  integer, public, parameter      :: I_H2SO4   = 1               !! sulfate aerosol composition
  integer, public, parameter      :: I_ICE     = 2               !! ice
  integer, public, parameter      :: I_WATER   = 3               !! water
    
  ! Define group, element, solute and gas indexes.
  integer, public, parameter      :: I_GRP_CRCN     = 1             !! sulfate aerosol
  integer, public, parameter      :: I_GRP_CRDICE   = 2             !! detrained ice
  integer, public, parameter      :: I_GRP_CRSICE   = 3             !! in-situ ice
  integer, public, parameter      :: I_GRP_CRLIQ    = 4             !! liquid drop

  integer, public, parameter      :: I_ELEM_CRCN    = 1             !! sulfate
  integer, public, parameter      :: I_ELEM_CRDICE  = 2             !! detrained ice
  integer, public, parameter      :: I_ELEM_CRSICE  = 3             !! in-situ ice
  integer, public, parameter      :: I_ELEM_CRCORE  = 4             !! sulfate core
  integer, public, parameter      :: I_ELEM_CRLIQ   = 5             !! water vapor

  integer, public, parameter      :: I_SOL_CRH2SO4  = 1             !! sulfuric acid

  integer, public, parameter      :: I_GAS_H2O      = 1             !! water vapor


  ! From Morrison & Gettelman [2008] and micro_mg.F90 (formerly cldwat2m_micro.F90)
  !
  ! NOTE: In the bin model, the bin boundaries are also important for determining the threshold,
  ! since the whole bin is autoconverted if the threshold is less than the bin midpoint radius.
  real(kind=f), public, parameter :: CAM_RHOCI        = 0.5_f    !! (g/cm3) MG bulk density for cloud ice
  real(kind=f), public, parameter :: CAM_RHOSN        = 0.1_f    !! (g/cm3) MG bulk density for snow

  
  ! Parameters and variabls that control the detrainment process.
  integer, parameter              :: NINTS_BINS           = 10        !! number of steps to integrate bin fractions
  integer, parameter              :: NINTS_SNOW           = 100       !! number of steps to integrate snow fractions
  
  
  real(kind=f), parameter         :: r_dliq_lnd           = 8e-4_f   !! detrained liquid radius (cm)
!  real(kind=f), parameter         :: r_dliq_lnd           = 18e-4_f   !! detrained liquid radius (cm)
  real(kind=f), parameter         :: r_dliq_ocn           = 8e-4_f  !! detrained liquid radius (cm)
!  real(kind=f), parameter         :: r_dliq_ocn           = 14e-4_f  !! detrained liquid radius (cm)
!  real(kind=f), parameter         :: r_dliq_ocn           = 18e-4_f  !! detrained liquid radius (cm)

  real(kind=f), parameter         :: snow_max_d           = 10000._f   !! maximum diameter for snow integration (um)

!  integer, parameter              :: MIN_DTEMP             = -60            !! Miniumum detrainment temperature (C)
  integer, parameter              :: MIN_DTEMP             = -90            !! Miniumum detrainment temperature (C)
  integer, parameter              :: NDTEMP                = -MIN_DTEMP + 1 !! Number of detrainment temperature bins
  
!  character(len=12), parameter    :: carma_dice_method    = "mono"
  real(kind=f), parameter         :: dice_snow_reff_mono  = 348e-4_f  !! Effective Radius of snow for monodisperse detrainment (cm)
  real(kind=f), parameter         :: r_dice_mono          = 25e-4_f   !! detrained ice radius, monodisperse (cm)
  real(kind=f), parameter         :: dice_loss            = 0.0_f     !! detrained fraction lost to precipitation, mondisperse
!  real(kind=f), parameter         :: dice_loss            = 0.004_f   !! detrained fraction lost to precipitation, mondisperse


  ! This distribution varies the size disribution as a function of temperature, with the
  ! distribution biased toards larer particles at warm temperature and small particles at
  ! cold temperatures. This fit is from eq. 7 of Heymsfield and Schmitt [2010]. The Jensen
  ! fit used above is similar to the cold end of this range.
  character(len=12), parameter    :: carma_dice_method    = "dist_hym2010"
!  real(kind=f), parameter         :: dist_hym2010_alpha   = 14.26_f       !! alpha (stratiform) from eq 7 in Heymsefield & Schmitt [2010] (cm -1)
!  real(kind=f), parameter         :: dist_hym2010_beta    = -0.0538_f     !! beta (stratiform)  from eq 7 in Heymsefield & Schmitt [2010] (cm -1)
  real(kind=f), parameter         :: dist_hym2010_alpha   = 2.425_f       !! alpha (convective) from eq 7 in Heymsefield & Schmitt [2010] (cm -1)
  real(kind=f), parameter         :: dist_hym2010_beta    = -0.088_f      !! beta (convective)  from eq 7 in Heymsefield & Schmitt [2010] (cm -1)

  real(kind=f)                    :: dice_snow_rmass(NDTEMP)              !! snow particle mass (kg)
  real(kind=f)                    :: dice_snow_fraction(NDTEMP)           !! detrained mass fraction, snow
  real(kind=f)                    :: dice_bin_fraction(NBIN, NDTEMP)      !! detrained mass fraction, ice bin

  logical, public, parameter      :: carma_do_mass_check  = .false.  ! If .true. then CARMA will check for mass loss by CARMA
  logical, public, parameter      :: carma_do_mass_check2 = .false.  ! If .true. then CARMA will check for mass loss (internal steps, e.g. detrain, diagnoseBIns, ...)
  logical, public, parameter      :: carma_do_mass_check3 = .false.  ! If .true. then CARMA will check for incoming mass loss (CAM -> CARMA)
  logical, public, parameter      :: carma_do_mass_fix    = .true.   ! If .true. then CARMA will fix for mass loss between cldice and ice bins
  logical, public, parameter      :: carma_do_print_fix   = .false.  ! If .true. then CARMA will print the value of the mass fix

  logical, public, parameter      :: carma_do_initice    = .true.   ! If .true. then CARMA carma prognositic bins are set from the bulk ice on the first timestep
  logical, public, parameter      :: carma_do_bulk_tend  = .true.   ! If .true. then update CAM bulk tendencies
  logical, public, parameter      :: carma_do_autosnow    = .false.  ! If .true. then the largest ice bin is autoconverted to snow at the end of the timestep.

  integer                         :: ixcldice
  integer                         :: ixnumice
  integer                         :: ixcldliq
  integer                         :: ixnumliq

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
    use physconst,          only: latice, latvap
    use ioFileMod,          only: getfil
    use wrap_nf

    implicit none
    
    type(carma_type), intent(inout)    :: carma     !! the carma object
    integer, intent(out)               :: rc        !! return code, negative indicates failure
        
    ! Local variables
    real(kind=f), parameter            :: rmin_ice   = 5.e-5_f  ! min radius for ice bins (cm)
    real(kind=f), parameter            :: rmin_cn    = 1.e-7_f  ! min radius for sulfate bins (cm)
    real(kind=f), parameter            :: RHO_CN     = 1.78_f   ! density of sulfate particles (g/cm)
    real(kind=f)                       :: rmassmin              ! mass of the first radius bin (g)
    real(kind=f)                       :: vmrat                 ! volume ratio between adjacent bin
    real(kind=f)                       :: rhoelem(NBIN)         ! element density per bin (g/cm3)
    real(kind=f)                       :: arat(NBIN)            ! projected area ratio
    integer                            :: maxbin                ! the bin number of the largest prognostic ice bin
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
        write(LUNOPRT,*) '  carma_mice_file       = ', trim(carma_mice_file)
        write(LUNOPRT,*) '  carma_sulfate_method  = ', trim(carma_sulfate_method)
      end if
    end if

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
          if (wave(i) <= warren_wave(j)) then
            if ((j > 1) .and. (wave(i) /= warren_wave(j))) then
              interp = (wave(i) - warren_wave(j-1)) / (warren_wave(j) - warren_wave(j-1))
              refidx_ice(i) = cmplx(warren_real(j-1) + interp*(warren_real(j) - warren_real(j-1)), warren_imag(j-1) + interp*(warren_imag(j) - warren_imag(j-1)))
            else
              refidx_ice(i) = cmplx(warren_real(j), warren_imag(j))
            endif
            
            exit
          end if
        end do
      end do
    end if


    ! Define the Groups
    !
    ! NOTE: If NWAVE > 0 then the group should have refractive indices defined.
    !
    ! NOTE: For CAM, the optional do_wetdep and do_drydep flags should be
    ! defined. If wetdep is defined, then the optional solubility factor
    ! should also be defined.
    rmassmin = (4._f / 3._f) * PI * (rmin_cn ** 3) * RHO_CN
!    vmrat = 4.0_f     ! For 16 bins
!    vmrat = 2.8_f     ! For 21 bins
    vmrat = 2.16_f     ! For 28 bins
!    vmrat = 2.0_f     ! For 32 bins
    
    ! Since these sulfates are prescribed, don't sediment them. This will save some
    ! processing time.
    call CARMAGROUP_Create(carma, I_GRP_CRCN, "Sulfate CN", rmin_cn, vmrat, I_SPHERE, 1._f, .false., &
                           rc, shortname="CRCN", rmassmin=rmassmin, do_mie=.false., &
                           cnsttype=I_CNSTTYPE_DIAGNOSTIC, do_vtran=.false.)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAGROUP_Create failed.')

    ! NOTE: For freezing and melting, the ice and water bins need to have the same mass.
    rmassmin = (4._f / 3._f) * PI * (rmin_ice ** 3) * RHO_I
    vmrat = 2.055_f     ! For 28 bins, Heysmfield Ice Density, cold

    ! If doing autoconversion of ice to snow, then the last bin will always be zero and
    ! there is no point making it an advected constituent.
    if (carma_do_autosnow) then
      maxbin = NBIN-1
    else
      maxbin = NBIN
    end if

    ! Make the aged detrained ice have a variable density to represent the complex set of
    ! possible shapes that we can't represent. This is based upon Heymsfield and 
    ! Westfield [2010] and Heysfield and Schmitt [2010].
    call CARMAGROUP_Create(carma, I_GRP_CRDICE, "Detrained Ice, Aged", rmin_ice, vmrat, I_SPHERE, 1._f, .true., &
                           rc, shortname="CRDICE", rmassmin=rmassmin, do_mie=carma_do_pheat, refidx=refidx_ice, &
                           ifallrtn=I_FALLRTN_HEYMSFIELD2010, imiertn=I_MIERTN_BOHREN1983, is_cloud=.true., maxbin=maxbin)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAGROUP_Create failed.')
    is_convtran1(2) = .true.

    ! Make the in-situ ice a plate, AR=6. This is based upon observations from Lawson
    ! et al. [2008]. AR=6 is for larger particles, so AR=3 is a compromise that is
    ! part way between that and more spheroidal particles that are likely at smaller sizes.
    !
    ! NOTE: All cloud particles should be convectively transported in the first phase of
    ! convection.
    !
    ! NOTE: All ice particles have the last bin as the one that gets autoconverted to
    ! snow at the end of the timestep and thus it does not need to be a prognostic bin.
!    call CARMAGROUP_Create(carma, I_GRP_CRSICE, "In-situ Ice", rmin_ice, vmrat, I_SPHERE, 1._f, .true., &
!    call CARMAGROUP_Create(carma, I_GRP_CRSICE, "In-situ Ice", rmin_ice, vmrat, I_HEXAGON, 1._f / 6._f, .true., &
    call CARMAGROUP_Create(carma, I_GRP_CRSICE, "In-situ Ice", rmin_ice, vmrat, I_HEXAGON, 1._f / 3._f, .true., &
                           rc, shortname="CRSICE", rmassmin=rmassmin, do_mie=carma_do_pheat, refidx=refidx_ice, &
                           ifallrtn=I_FALLRTN_HEYMSFIELD2010, imiertn=I_MIERTN_BOHREN1983, is_cloud=(.not. carma_do_clearsky), maxbin=maxbin)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAGROUP_Create failed.')
    is_convtran1(3) = .true.

    ! Water drops are spherical.
    call CARMAGROUP_Create(carma, I_GRP_CRLIQ, "Water Drop", rmin_ice, vmrat, I_SPHERE, 1._f, .false., &
                           rc, shortname="CRLIQ", rmassmin=rmassmin, do_mie=.false., &
                           cnsttype=I_CNSTTYPE_DIAGNOSTIC, is_cloud=.true., do_vtran=.false.)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAGROUP_Create failed.')
    is_convtran1(4) = .true.


    ! Define the Elements
    !
    ! NOTE: For CAM, the optional shortname needs to be provided for the group. These names
    ! should be 6 characters or less and without spaces.
    call CARMAELEMENT_Create(carma, I_ELEM_CRCN, I_GRP_CRCN, "Sulfate CN", RHO_CN, I_INVOLATILE, I_H2SO4, rc, shortname="CRCN", isolute=I_SOL_CRH2SO4)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAElement_Create failed.')

    ! The density of ice is changed based on the maximum dimensions of ice particles
    ! as a function of mass from Heymsfield and Schmitt [2010].
!    call rhoice_heymsfield2010(carma, RHO_I, I_GRP_CRDICE, "conv", rhoelem, arat, rc)
    call rhoice_heymsfield2010(carma, RHO_I, I_GRP_CRDICE, "warm", rhoelem, arat, rc)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::rhoice_heymsfield2010 failed.')
    
    call CARMAELEMENT_Create(carma, I_ELEM_CRDICE, I_GRP_CRDICE, "Detrained Ice", RHO_I, I_VOLATILE, I_ICE, rc, shortname="CRDICE", rhobin=rhoelem, arat=arat)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAElement_Create failed.')
    

    call CARMAELEMENT_Create(carma, I_ELEM_CRSICE, I_GRP_CRSICE, "In-situ Ice", RHO_I, I_VOLATILE, I_ICE, rc, shortname="CRSICE")
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAElement_Create failed.')


    call CARMAELEMENT_Create(carma, I_ELEM_CRCORE, I_GRP_CRSICE, "Core Mass", RHO_CN, I_COREMASS, I_H2SO4, rc, shortname="CRCORE", isolute=I_SOL_CRH2SO4)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAElement_Create failed.')


    call CARMAELEMENT_Create(carma, I_ELEM_CRLIQ, I_GRP_CRLIQ, "Water Drop", RHO_W, I_VOLATILE, I_WATER, rc, shortname="CRLIQ")
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAElement_Create failed.')


    ! Define the Solutes
    call CARMASOLUTE_Create(carma, I_SOL_CRH2SO4, "Sulfuric Acid", 2, 98._f, 1.38_f, rc, shortname="CRH2SO4")
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMASOLUTE_Create failed.')

    
    ! Define the Gases
    call CARMAGAS_Create(carma, I_GAS_H2O, "Water Vapor", WTMOL_H2O, I_VAPRTN_H2O_MURPHY2005, I_GCOMP_H2O, rc, shortname="Q", ds_threshold=-0.2_f)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAGAS_Create failed.')
 
    
    ! Define the Processes
    
    ! Detrained Ice, Aged
    call CARMA_AddGrowth(carma, I_ELEM_CRDICE, I_GAS_H2O, rc)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_AddGrowth failed.')

    call CARMA_AddNucleation(carma, I_ELEM_CRDICE, I_ELEM_CRLIQ, I_ICEMELT, -latice*1e4_f, rc)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_AddNucleation failed.')

    call CARMA_AddCoagulation(carma, I_GRP_CRDICE, I_GRP_CRDICE, I_GRP_CRDICE, I_COLLEC_DATA, rc)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_AddCoagulation failed.')


    ! In-Situ Ice
    call CARMA_AddGrowth(carma, I_ELEM_CRSICE, I_GAS_H2O, rc)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_AddGrowth failed.')

    ! NOTE: For now, assume the latent heat for nucleation is the latent of of fusion of
    ! water, using the CAM constant (scaled from J/kg to erg/g).
    !
    ! NOTE: Since the sulfates are not seen as part of the water/energy budget in CAM, don't
    ! include any latent heat from the freezing of the sulfate liquid. The latent heat of
    ! the gas associated with nucleation is accounted for.
    call CARMA_AddNucleation(carma, 1, 4, I_AERFREEZE + I_AF_KOOP_2000, 0._f, rc, igas=1, ievp2elem=1)
!    call CARMA_AddNucleation(carma, I_ELEM_CRCN, I_ELEM_CRCORE, I_AERFREEZE + I_AF_KOOP_2000 + I_AF_MURRAY_2010, 0._f, rc, igas=I_GAS_H2O, ievp2elem=I_ELEM_CRCN)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_AddNucleation failed.')

    call CARMA_AddNucleation(carma, I_ELEM_CRSICE, I_ELEM_CRLIQ, I_ICEMELT, -latice*1e4_f, rc)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_AddNucleation failed.')

    call CARMA_AddCoagulation(carma, I_GRP_CRSICE, I_GRP_CRSICE, I_GRP_CRSICE, I_COLLEC_DATA, rc)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_AddCoagulation failed.')


    ! Water Drop
    call CARMA_AddGrowth(carma, I_ELEM_CRLIQ, I_GAS_H2O, rc)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_AddGrowth failed.')

    call CARMA_AddNucleation(carma, I_ELEM_CRLIQ, I_ELEM_CRDICE, I_DROPFREEZE, latice*1e4_f, rc)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_AddNucleation failed.')

    return
  end subroutine CARMA_DefineModel
  

  !! Defines all the CARMA components (groups, elements, solutes and gases) and process
  !! (coagulation, growth, nucleation) that will be part of the microphysical model.
  !!
  !!  @version May-2009 
  !!  @author  Chuck Bardeen 
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
    
    real(kind=f)                         :: t(pver)             ! temperature (K)
    real(kind=f)                         :: mmr_ice(NBIN, pver) ! ice mass mixing ratio (kg/kg)
    real(kind=f)                         :: mmr_liq(NBIN, pver) ! liquid mass mixing ratio (kg/kg)
    real(kind=f)                         :: r_ice(NBIN)         ! ice radius bins (cm)
    real(kind=f)                         :: r_liq(NBIN)         ! liquid radius bins (cm)
    
    real(kind=f)                         :: ice_fraction        ! fraction of detrained condensate that is ice
    real(kind=f)                         :: mass_liq            ! detrainment rate of liquid (kg/kg/s)
    real(kind=f)                         :: mass_ice            ! detrainment rate of ice (kg/kg/s)
    real(kind=f)                         :: mass_snow           ! detrainment rate of snow (kg/kg/s)
    real(kind=f)                         :: mass_dlf            ! detrained mass (m/s)
    integer                              :: k                   ! vertical index
    integer                              :: ibin                ! bin index
    integer                              :: itemp               ! termperature index

    real(r8)                             :: iceMass(pver)       ! ice mass mixing ratio (kg/kg)
    real(r8)                             :: iceNumber(pver)     ! ice number mixing ratio (#/kg)
    real(r8)                             :: snowMass(pver)      ! snow mass mixing ratio (kg/kg)
    real(r8)                             :: snowNumber(pver)    ! snow number (#/kg)
    real(r8)                             :: snowSurface         ! snow on surface (kg/m2)
    real(r8)                             :: waterMass(pver)     ! ice mass mixing ratio (kg/kg)
    real(r8)                             :: waterNumber(pver)   ! ice number mixing ratio (#/kg)
    real(r8)                             :: rainSurface         ! rain on surface (kg/m2)
    real(r8)                             :: newSnow             ! snow mass (kg)
    real(r8)                             :: newRain             ! rain mass mass (kg)
    
    logical                              :: do_thermo           ! do thermodynamics?


    ! Default return code.
    rc = RC_OK
    
    call CARMA_Get(carma, rc, do_thermo=do_thermo)
    if (rc < RC_OK) call endrun('CARMA_Detrain::CARMA_Get failed.')
    
    ! Put all of the detraining cloud water from convection into the large scale cloud.
    ! put detraining cloud water into liq and ice based on temperature partition
    call CARMAGROUP_Get(carma, I_GRP_CRDICE, rc, r=r_ice(:))
    if (rc < RC_OK) call endrun('CARMA_Detrain::CARMAGROUP_Get failed.')

    call CARMAGROUP_Get(carma, I_GRP_CRLIQ, rc, r=r_liq(:))
    if (rc < RC_OK) call endrun('CARMA_Detrain::CARMAGROUP_Get failed.')
    
    ! Account for the reserved ice that is being detrained in the precipitation.
!    prec_str(icol) = prec_str(icol) - rliq(icol)
!    rliq(icol)     = 0._f
    
    call CARMASTATE_GetState(cstate, rc, t=t)
    if (rc < RC_OK) call endrun('CARMA_Detrain::CARMAGROUP_GetState failed.')
    
    ! Determine the amount of detrainment that could be used to saturate the
    ! atmosphere with respect to liquid. For GCM scales, assume that three things
    ! happen to detrained condensate:
    !
    !   1) large particles will fallout as snow or rain
    !   2) will be converted to vapor
    !   3) will remain as ice
    !
    ! Because of the large scales of the GCM and because this is a stratiform
    ! parameterization, a lot of the condensate that hasn't fallen out will
    ! increase the humidity (i.e. detrained anvil evaporates or falls out entirely
    ! with 100 km of the convection).
    mmr_ice(:, :)  = 0._f
    mmr_liq(:, :)  = 0._f
    
    do k = 1,pver
    
      ! Remove amount being detrained from rliq and prec_str.
      mass_dlf = dlf(icol, k) * (state%pdel(icol, k) / gravit) / 1000._f
      prec_str(icol) = prec_str(icol) - mass_dlf 
      rliq(icol)     = rliq(icol)     - mass_dlf
      
    
      if (t(k) > 268.15_f) then
        ice_fraction = 0.0_f
      else if (t(k) < 238.15_f) then
        ice_fraction = 1.0_f
      else
        ice_fraction = (268.15_f - t(k)) / 30._f
      end if
      
      itemp = max(-max(MIN_DTEMP, nint(t(k) - T0)), 0) + 1
      
      mass_liq  = dlf(icol, k) * (1._f - ice_fraction)
      mass_ice  = dlf(icol, k) * ice_fraction * (1._f - dice_snow_fraction(itemp))
      mass_snow = dlf(icol, k) * ice_fraction * dice_snow_fraction(itemp)

      ! Calculate the detrainment of ice and liquid into the appropriate CARMA
      ! bins.
      !      
      ! Scale the size based on whether the surface is land or ocean. This
      ! assumes that there are more aerosols over land, reducing the detrainment
      ! size. This is similar to the c0_lnd and c0_ocn parameter split done in
      ! the convective parameterization.
      !
      ! NOTE: This should really be tied to aerosol amount, not land fraction.
      do ibin = 1, NBIN
      
        ! Assume detrained cloud water is monodisperse.
        if (r_liq(ibin) >= r_dliq_ocn) then
          mmr_liq(ibin, k) = mmr_liq(ibin, k) + (mass_liq * dt) * (1._f - cam_in%landfrac(icol))
          exit
        end if
      end do
      
      do ibin = 1, NBIN
      
        ! Assume detrained cloud water is monodisperse.
        if (r_liq(ibin) >= r_dliq_lnd) then
          
          mmr_liq(ibin, k) = mmr_liq(ibin, k) + (mass_liq * dt) * cam_in%landfrac(icol)
          exit
        end if
      end do

      ! Detrain cloud ice into the bins according to the predefined distribution.
      do ibin = 1, NBIN
   
        ! Detrain using a size distribution (log-normal in mass). The table has
        ! already bin setup during initialization indicating the fraction of the mass
        ! that goes into each bin.
        !
        ! NOTE: Since snow has already been removed, but was part of the fractions
        ! in the bins, scale the bin fractions so that it sums to 1.
        mmr_ice(ibin, k) = mmr_ice(ibin, k) + dice_bin_fraction(ibin, itemp) * (mass_ice * dt)
      end do

      ! The large portion of the distribution can go directly to snow,
      ! since it is too big to be represented in the bin strucutre.
      tnd_qsnow(icol, k) = mass_snow
      tnd_nsnow(icol, k) = mass_snow / dice_snow_rmass(itemp)
      
      rliq(icol) = rliq(icol) + mass_snow * (state%pdel(icol, k) / gravit) / 1000._f
           
      ! Account for latent heat release during freezing. By default the detrained
      ! condensate is assumed to be liquid for energy balance.
      t(k) = t(k) + ((mass_ice + mass_snow) * latice * dt / cpair)
    end do


    do ibin = 1, NBIN
      call CARMASTATE_SetDetrain(cstate, I_ELEM_CRLIQ, ibin, mmr_liq(ibin, :), rc)
      if (rc < RC_OK) call endrun('CARMA_Detrain::CARMAState_SetBin failed.')
      
  
      call CARMASTATE_SetDetrain(cstate, I_ELEM_CRDICE, ibin, mmr_ice(ibin, :), rc)
      if (rc < RC_OK) call endrun('CARMA_Detrain::CARMAState_SetBin failed.')
    end do
    
        
    if (do_thermo) then
      call CARMASTATE_SetState(cstate, rc, t(:))
    end if
    
    ! Check for total water conservation by CARMA.
    if (carma_do_mass_check2) then
      call CARMA_GetTotalWaterAndRain(carma, cstate, waterMass, waterNumber, rainSurface, rc)
      call CARMA_GetTotalIceAndSnow(carma, cstate, .false., iceMass, iceNumber, snowMass, snowNumber, snowSurface, rc)

      call CARMA_CheckMassAndEnergy(carma, cstate, .false., "CARMA_Detrain", state, &
        icol, dt, rliq, prec_str, snow_str, waterMass, iceMass, snowMass, rc)    
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
    
    real(r8)                              :: mu(pver)       ! spectral width parameter of droplet size distr
    real(r8)                              :: lambda(pver)   ! slope of cloud liquid size distr
    real(r8)                              :: mmr(NBIN,pver) ! elements mass mixing ratio

    real(kind=f)                          :: r(NBIN)      ! bin mean radius
    real(kind=f)                          :: dr(NBIN)     ! bin radius width
    real(kind=f)                          :: rmass(NBIN)  ! bin mass
    
    integer                               :: igroup       ! group index
    integer                               :: ielem        ! element index
    integer                               :: ibin         ! bin index
    integer                               :: k            ! vertical index
    
    real(r8)                              :: iceMass(pver)      ! ice mass mixing ratio (kg/kg)
    real(r8)                              :: iceNumber(pver)    ! ice number mixing ratio (#/kg)
    real(r8)                              :: snowMass(pver)     ! snow mass mixing ratio (kg/kg)
    real(r8)                              :: snowNumber(pver)   ! snow number (#/kg)
    real(r8)                              :: snowSurface        ! snow on surface (kg/m2)
    real(r8)                              :: carma_ice          ! total cldice from CARMA bins (kg/kg)
    real(r8)                              :: waterMass(pver)    ! ice mass mixing ratio (kg/kg)
    real(r8)                              :: waterNumber(pver)  ! ice number mixing ratio (#/kg)
    real(r8)                              :: rainSurface        ! rain on surface (kg/m2)
    real(r8)                              :: carma_water        ! total cldliq from CARMA bins (kg/kg)
    real(r8)                              :: diff

    ! Aerosol size distribution
    real(r8), parameter                   :: n    = 100._r8     ! concentration (cm-3) 
    real(r8), parameter                   :: r0   = 2.5e-6_r8   ! mean radius (cm)
    real(r8), parameter                   :: rsig = 1.5_r8      ! distribution width
    
    real(r8)                              :: arg1(NBIN)
    real(r8)                              :: arg2(NBIN)
    real(r8)                              :: rhop(NBIN)         ! particle mass density (kg/m3)
    real(r8)                              :: totalrhop          ! total particle mass density (kg/m3)
    real(kind=f)                          :: rhoa_wet(pver)     ! air density (g/cm3)

    real(r8)                              :: rliq_new(pcols)    ! vertical integral of liquid not yet in q(ixcldliq)

    integer                               :: LUNOPRT
    logical                               :: do_print
    real(r8)                              :: lat
    real(r8)                              :: lon

    real(r8), pointer, dimension(:, :)    :: sulf               ! last saturation wrt ice
    integer                               :: lchnk              ! chunk identifier
    integer                               :: itim
    
    character(len=8)                      :: c_name             ! constituent name
    

  1 format(/,'CARMA_DiagnoseBins::ERROR - CAM ice mass conservation error, icol=',i4,', iz=',i4,',lat=',&
              f7.2,',lon=',f7.2,',cam=',e16.10,',carma=',e16.10,',rer=',e9.3)
  2 format(/,'CARMA_DiagnoseBins::ERROR - CAM liquid mass conservation error, icol=',i4,', iz=',i4,',lat=',&
              f7.2,',lon=',f7.2,',cam=',e16.10,',carma=',e16.10,',rer=',e9.3)

    ! Default return code.
    rc = RC_OK
    
    call CARMA_Get(carma, rc, do_print=do_print, LUNOPRT=LUNOPRT)
    
    ! Get the air density.
    call CARMASTATE_GetState(cstate, rc, rhoa_wet=rhoa_wet)
    if (rc < RC_OK) call endrun('CARMA_DiagnoseBins::CARMASTATE_GetState failed.')


    ! Aerosols
    !
    igroup = 1
    ielem  = 1
    
    ! Use a fixed aerosols distribution.
    if ((carma_sulfate_method == "fixed") .or. (carma_sulfate_method == "bulk")) then
    
      call CARMAGROUP_Get(carma, igroup, rc, r=r, dr=dr, rmass=rmass)
      if (rc < RC_OK) call endrun('CARMA_DiagnoseBins::CARMAGROUP_Get failed.')
    
      arg1(:) = n * dr(:) / (sqrt(2._f*PI) * r(:) * log(rsig))
      arg2(:) = -((log(r(:)) - log(r0))**2) / (2._f*(log(rsig))**2)

      ! kg/m3
      rhop(:)   = arg1(:) * exp(arg2(:)) * rmass(:) * 1e6_f / 1e3_f

    
      if (carma_sulfate_method == "bulk") then
        totalrhop = sum(rhop(:))

        ! Get the index for the prescribed sulfates. This gives the mmr that should be
        ! present at this location. Use this to scale the size distribution that CARMA
        ! will generate.
        lchnk = state%lchnk
        itim  = pbuf_old_tim_idx()

        call pbuf_get_field(pbuf, pbuf_get_index('sulf'), sulf, (/1,1,itim/),(/pcols,pver,1/))
      end if
    end if
    
    do ibin = 1, NBIN
    
      ! Use a fixed mixing ration.
      if (carma_sulfate_method == "fixed") then
        mmr(ibin, :) = rhop(ibin) / rhoa_wet(:)
      end if
      
      
      ! Since bulk aerosols don't have a size distribution, use the fixed
      ! distribution for the shape of the distribution, but scale the total
      ! mass to the prescribed value.
      if (carma_sulfate_method == "bulk") then
        mmr(ibin, :) = rhop(ibin) / totalrhop * sulf(icol, :)
      end if

      ! Use the CRCNxx fields from a special prescribed aerosol file that has
      ! results from a CARMA simulation of sulfates. This will set the magnitude
      ! and the size distribution.
      if (carma_sulfate_method == "carma") then
        ! Get the index for the prescribed sulfates.
        lchnk = state%lchnk
        itim  = pbuf_old_tim_idx()
        write(c_name, '(A, I2.2)') "CRCN", ibin

        call pbuf_get_field(pbuf, pbuf_get_index(c_name), sulf, (/1,1,itim/),(/pcols,pver,1/))
        mmr(ibin, :) = sulf(icol, :)
      end if
      
      
      call CARMASTATE_SetBin(cstate, ielem, ibin, mmr(ibin, :), rc)
      if (rc < RC_OK) call endrun('CARMA_DiagnoseBins::CARMAGROUP_SetBin failed.')
    end do
    

    ! Cloud Ice & Snow
    !
    ! Cloud ice is maintained in advected species in CARMA, and we are only
    ! concerned with snow production by CARMA.
    !
    ! NOTE: To allow this code to be tested when not doing the cloud ice, but
    ! either doing nothing or doing detrainment, use the ice properties to convert
    ! from the 2 moment values to a size distribution.
    !
    ! NOTE: To keep mass and energy conservation happy, on the first step when
    ! camra_do_cldice is true, we take the bulk values and convert them
    ! into bins; however, this might cause issues with the CARMA growth code.
    if ((.not. carma_do_cldice) .or. (is_first_step() .and. carma_do_initice)) then
      igroup = I_GRP_CRDICE
      ielem  = I_ELEM_CRDICE
      
      call CARMAGROUP_Get(carma, igroup, rc, r=r, dr=dr, rmass=rmass)
      if (rc < RC_OK) call endrun('CARMA_DiagnoseBins::CARMAGROUP_Get failed.')
  
      ! Need to determine the shape parameters for the size distribution. It
      ! would be nice if these routines came from the MG microphysics module,
      ! but until then their code with be duplicated here.
      call CARMA_GetGammaParmsForIce(carma, state%q(icol, :, ixcldice), state%q(icol, :, ixnumice), rhoa_wet(:), mu(:), lambda(:), rc)
      if (rc < RC_OK) call endrun('CARMA_DiagnoseBins::CARMA_GetGammaParmsForIce failed.')
  
      call CARMA_GetMmrFromGamma(carma, r(:), dr(:), rmass(:), state%q(icol, :, ixcldice), state%q(icol, :, ixnumice), mu(:), lambda(:), mmr(:, :), rc)
      if (rc < RC_OK) call endrun('CARMA_DiagnoseBins::CARMA_GetMmrFromGamma failed.')
    
      do ibin = 1, NBIN
        call CARMASTATE_SetBin(cstate, ielem, ibin, mmr(ibin, :), rc)
        if (rc < RC_OK) call endrun('CARMA_DiagnoseBulk::CARMASTATE_SetBin failed.')
      end do
    else
    
      ! If CARMA is keeping track of ice, then total up the detrained and in-situ ice to
      ! make sure that no one else has meesed with the ice fields. If changes were made,
      ! then adjustments need to be made to the totals to prevent mass and energy conservation
      ! errors within CAM. The difference could be accounted for in snow_str and prec_str.
      !
      ! NOTE: Advection, diffusion, ... may have affected the tracer values since the
      ! previous time step, however, the tracer correlations need to remain intact for
      ! CARMA to work properly. Also, no special processing can occur on the cldice fields
      ! outside of CARMA, since they are merely diagnostic fields of the CARMA state.
      if (carma_do_mass_check3 .or. carma_do_mass_fix) then
      
        call CARMA_GetTotalIceAndSnow(carma, cstate, .false., iceMass, iceNumber, snowMass, snowNumber, snowSurface, rc)
        if (rc < RC_OK) call endrun('CARMA_DiagnoseBins::CARMA_GetTotalIceAndSnow failed.')
      
        do k = 1, pver
          
          ! NOTE: CAM resets cloud ice less than 1e-36 to 0 in physics_update, so ignore values smaller
          ! than that.
          carma_ice = iceMass(k) + snowMass(k)
          if (carma_ice < 1e-36_r8) then
            carma_ice = 0._r8
          end if
          
          if (carma_ice /= state%q(icol, k, ixcldice)) then

            if (carma_do_mass_check3) then
              if (abs(carma_ice - state%q(icol, k, ixcldice)) / max(abs(carma_ice), &
                  abs(state%q(icol, k, ixcldice))) >= 1e-10_r8)  then
                if (do_print) then
                  call CARMASTATE_Get(cstate, rc, lat=lat, lon=lon)
                  if (rc < RC_OK) call endrun('CARMA_DiagnoseBins::CARMASTATE_Get failed.')
                  
                  write(LUNOPRT,1) icol, k, lat, lon, state%q(icol, k, ixcldice), &
                    carma_ice, (carma_ice - state%q(icol, k, ixcldice)) / max(abs(carma_ice), &
                    abs(state%q(icol, k, ixcldice)))
        
                  write(LUNOPRT,*) "  CAM cldice   :  ", state%q(icol, k, ixcldice)
                  write(LUNOPRT,*) ""
                  write(LUNOPRT,*) "  CARMA cldice :  ", iceMass(k) + snowMass(k)
                  write(LUNOPRT,*) "  CARMA ice    :  ", iceMass(k)
                  write(LUNOPRT,*) "  CARMA snow   :  ", snowMass(k)
                end if
              end if
            end if
          
            if (carma_do_mass_fix) then
         
              diff = (state%q(icol, k, ixcldice) - (iceMass(k) + snowMass(k))) * (state%pdel(icol, k) / gravit) / dt / 1000._r8
  
              snow_str(icol) = snow_str(icol) + diff
              prec_str(icol) = prec_str(icol) + diff
              
              if (carma_do_print_fix) then
                if (do_print) write(LUNOPRT,*) "  CARMA_DiagnoseBins::WARNING - Adjusting prec_str for ice mass difference", icol, k, (state%q(icol, k, ixcldice) - (iceMass(k) + snowMass(k)))
              end if
            end if
          end if
        end do
      end if
    end if
    

    ! Water Drops
    !
    ! Use the CAM mass and number (CLDLIQ and NUMLIQ) to determine an initial
    ! size distribution.
    igroup = I_GRP_CRLIQ
    ielem  = I_ELEM_CRLIQ
    
    call CARMAGROUP_Get(carma, igroup, rc, r=r, dr=dr, rmass=rmass)
    if (rc < RC_OK) call endrun('CARMA_DiagnoseBins::CARMAGROUP_Get failed.')

    ! Need to determine the shape parameters for the size distribution. It
    ! would be nice if these routines came from the MG microphysics module,
    ! but until then their code with be duplicated here.
    call CARMA_GetGammaParmsForLiq(carma, state%q(icol, :, ixcldliq), state%q(icol, :, ixnumliq), &
      rhoa_wet(:), mu(:), lambda(:), rc)
    if (rc < RC_OK) call endrun('CARMA_DiagnoseBins::CARMA_GetGammaParmsForLiq failed.')

    call CARMA_GetMmrFromGamma(carma, r(:), dr(:), rmass(:), state%q(icol, :, ixcldliq), &
      state%q(icol, :, ixnumliq), mu(:), lambda(:), mmr(:, :), rc)
    if (rc < RC_OK) call endrun('CARMA_DiagnoseBins::CARMA_GetMmrFromGamma failed.')
  
    do ibin = 1, NBIN
      call CARMASTATE_SetBin(cstate, ielem, ibin, mmr(ibin, :), rc)
      if (rc < RC_OK) call endrun('CARMA_DiagnoseBins::CARMASTATE_SetBin failed.')
    end do
    
    
    if (carma_do_mass_check2 .or. carma_do_mass_check3) then
    
      ! Check to see of the mass that we get back adds up.
      call CARMA_GetTotalWaterAndRain(carma, cstate, waterMass, waterNumber, rainSurface, rc)
      if (rc < RC_OK) call endrun('CARMA_DiagnoseBins::CARMA_GetTotalWaterAndRain failed.')
      
      
      if (carma_do_mass_check3) then
        do k = 1, pver
                
          carma_water = waterMass(k)
          if (carma_water < 1e-38_r8) then
            carma_water = 0._r8
          end if
    
          ! The routine that provides the modal properties for water has a miniumum of 1e-18.
          ! This causes problems in comparisons, since smaller qc values are seen in the data,
          ! but CARMA's bins won't have values that small.
          if (carma_water /= state%q(icol, k, ixcldliq)) then
            if (abs(carma_water - state%q(icol, k, ixcldliq)) / max(abs(carma_water), &
                abs(state%q(icol, k, ixcldliq))) >= 1e-10_r8)  then
              if (do_print) then
                call CARMASTATE_Get(cstate, rc, lat=lat, lon=lon)
                if (rc < RC_OK) call endrun('CARMA_DiagnoseBins::CARMASTATE_Get failed.')

                write(LUNOPRT,2) icol, k, lat, lon, state%q(icol, k, ixcldliq), &
                  carma_water, (carma_water - state%q(icol, k, ixcldliq)) / max(abs(carma_water), &
                  abs(state%q(icol, k, ixcldliq)))
      
                write(LUNOPRT,*) "  CAM cldliq   :  ", state%q(icol, k, ixcldliq)
                write(LUNOPRT,*) ""
                write(LUNOPRT,*) "  CARMA cldliq :  ", waterMass(k)
              end if
            end if
          end if
        end do
      end if
      
  
      ! Check for total water conservation by CARMA.
      if (carma_do_mass_check2) then
        call CARMA_GetTotalIceAndSnow(carma, cstate, .false., iceMass, iceNumber, snowMass, snowNumber, snowSurface, rc)
        if (rc < RC_OK) call endrun('CARMA_DiagnoseBins::CARMA_GetTotalIceAndSnow failed.')

        ! The detrained ice is not include yet, so ignore rliq.
        rliq_new(:) = 0._f
        call CARMA_CheckMassAndEnergy(carma, cstate, .false., "CARMA_DiagnoseBins", state, &
          icol, dt, rliq_new, prec_str, snow_str, waterMass, iceMass, snowMass, rc)  
      if (rc < RC_OK) call endrun('CARMA_DiagnoseBins::CARMA_CheckMassAndEnergy failed.')
      end if
    end if
    
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
 
    ! These values are chosen to match up with how small cloud ice values are handled in
    ! micro_mg.
    real(r8), parameter                  :: qsmall     = 1.e-18_r8     ! min mixing ratio 
    real(r8), parameter                  :: omsm       = 0.99999_r8    ! Prevents roundoff errors


    integer                              :: igroup    ! group index
    integer                              :: ielem     ! element index
    integer                              :: ibin      ! bin index
    integer                              :: icore     ! core index
    integer                              :: icorelem(NELEM) ! core indexes for group
    integer                              :: ncore     ! number of core elements
    integer                              :: itim
    
    real(kind=f)                         :: iceMass(pver)      ! ice mass mixing ratio (kg/kg)
    real(kind=f)                         :: iceNumber(pver)    ! ice number mixing ratio (#/kg)
    real(kind=f)                         :: snowMass(pver)     ! snow mass mixing ratio (kg/kg)
    real(kind=f)                         :: snowNumber(pver)   ! snow number (#/kg)
    real(kind=f)                         :: snowSurface        ! snow on surface (kg/m2)
    real(kind=f)                         :: waterMass(pver)    ! water mass mixing ratio (kg/kg)
    real(kind=f)                         :: waterNumber(pver)  ! water number mixing ratio (#/kg)
    real(kind=f)                         :: rainSurface        ! rain on surface (kg/m2)
    real(kind=f)                         :: iceRe(pver)        ! ice effective radius (m)

    real(r8)                             :: newRain            ! [Total] sfc flux of rain from stratiform (m/s)
    real(r8)                             :: newSnow            ! [Total] sfc flux of snow from stratiform (m/s)

    real(kind=f)                         :: mmr(pver)          ! mass mixing ratio (#/kg)
    real(kind=f)                         :: mmrcore(pver)      ! core mass mixing ratio (#/kg)
    real(kind=f)                         :: nmr(pver)          ! number mixing ratio (#/kg)
    real(kind=f)                         :: r(NBIN)            ! radius (cm)
    real(kind=f)                         :: sfc                ! surface mass (kg/m2)
    real(kind=f)                         :: sfccore            ! core surface mass (kg/m2)


    ! Default return code.
    rc = RC_OK
    
    
    ! Aerosols
    !
    ! Currently, we are just using a fixed aerosol size distribution, but in the
    ! future this could be linked to the model aerosols.
    

    ! Cloud Ice & Snow
    !
    ! Determine the changes to cloud ice (mass and number) and snow (mass and number)
    ! by looking at the totals if the detrained and in-situ ice.
    !
    ! Ice particles in the largest bin are treated as snow rather than ice.
    
    ! Get the total ice.
    call CARMA_GetTotalIceAndSnow(carma, cstate, .true., iceMass, iceNumber, snowMass, snowNumber, snowSurface, rc, iceRe=iceRe)    
    
    ! Calculate the tendencies on CLDICE, NUMICE, QSNOW and NSNOW
    if (carma_do_bulk_tend) then
    
      ptend%q(icol, :, ixcldice) = (iceMass(:) - state%q(icol, :, ixcldice)) / dt
      ptend%q(icol, :, ixnumice) = ((iceNumber(:) - state%q(icol, :, ixnumice)) / dt)

      where(iceMass(:) < qsmall)
        ptend%q(icol, :, ixcldice) = (-state%q(icol, :, ixcldice)) / dt
        ptend%q(icol, :, ixnumice) = (-state%q(icol, :, ixnumice)) / dt
      end where

      ! Snow is not a constituent, so write this information into the physics buffer.      
      tnd_qsnow(icol, :) = tnd_qsnow(icol, :) + snowMass(:) / dt
      tnd_nsnow(icol, :) = tnd_nsnow(icol, :) + snowNumber(:) / dt
      
      ! Now we need to change the reserve liquid. This was indicating the amount of
      ! water than was not in the atmosphere because it was in convection (dlf). Now
      ! we have included that water, but we have removed water representing snow in
      ! the atmosphere. This needs to be communicated to the CAM microphysics which
      ! will take care of actually precipitating or evaporating the snow.
      rliq(icol) = rliq(icol) + sum(snowMass(:) * (state%pdel(icol, :) / gravit)) / dt / 1000._r8
      
      ! The ice effective radius is used by the radiation code; however, it uses a mass
      ! weighted effective diameter in um.
      re_ice(icol, :) = iceRe
    end if

          
    ! Water Drops
    !
    ! Calcualte the total mass and total number of the water drops, and then
    ! determine the appropriate tendencies.
    call CARMA_GetTotalWaterAndRain(carma, cstate, waterMass, waterNumber, rainSurface, rc)    

    ! Calculate the tendencies on CLDLIQ and NUMLIQ
    if (carma_do_bulk_tend) then
    
      ! In CAM in cldwat2m, a couple of things are done:
      !
      !   1) If cldliq < qsmall, then the number desnity is set to 0.
      !   2) to keep from overshooting into negative values, they don't try to drive
      !      the value all the way to 0.
      ptend%q(icol, :, ixcldliq) = (waterMass(:) - state%q(icol, :, ixcldliq)) / dt
      
      ptend%q(icol, :, ixnumliq) = (waterNumber(:) - state%q(icol, :, ixnumliq)) / dt

      where(waterMass(:) < qsmall)
        ptend%q(icol, :, ixnumliq) = (-state%q(icol, :, ixnumliq)) / dt
      end where
    end if
        
        
    ! For mass balance, we also need to supply the total precipation and snow. Not
    ! all of the snow may make the ground, but that will be determined later in the
    ! MG microphysics. For now, we need to account for all condensate that is not
    ! in CLDICE or CLDLIQ.
    !
    ! Need the 1000. to convert from kg/m2/s to m/s
    newSnow = snowSurface
    newRain = rainSurface

    snow_sed(icol) = snow_sed(icol) + newSnow / dt / 1000._r8
    prec_sed(icol) = prec_sed(icol) + (newRain + newSnow) / dt / 1000._r8
    
    snow_str(icol) = snow_str(icol) + newSnow / dt / 1000._r8
    prec_str(icol) = prec_str(icol) + (newRain + newSnow) / dt / 1000._r8

    ! Check for total water conservation by CARMA.
    if (carma_do_mass_check) then
      
      ! The CAM state has not been updated yet, so compare the original CAM state
      ! with the new CARMA state.
      call CARMA_CheckMassAndEnergy(carma, cstate, .true., "CARMA_DiagnoseBulk", state, &
        icol, dt, rliq, prec_str, snow_str, waterMass, iceMass, snowMass, rc)
    end if

    return
  end subroutine CARMA_DiagnoseBulk


  !! Allows the model to perform its own initialization in addition to what is done
  !! by default in CARMA_init.
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  subroutine CARMA_InitializeModel(carma, lq_carma, rc)
    use constituents,     only: cnst_get_ind, pcnst

    implicit none
    
    type(carma_type), intent(in)       :: carma                 !! the carma object
    logical, intent(inout)             :: lq_carma(pcnst)       !! flags to indicate whether the constituent could have a CARMA tendency
    integer, intent(out)               :: rc                    !! return code, negative indicates failure

    integer                            :: ibin                  ! bin index
    integer                            :: i
    integer                            :: itemp                 ! temperature index
    integer                            :: LUNOPRT
    
    logical                            :: do_print_init
    logical                            :: do_grow
    logical                            :: do_detrain
    logical                            :: do_thermo

    real(kind=f)                       :: r(NBIN)               ! bin center radius (cm)
    real(kind=f)                       :: dr(NBIN)              ! bin width (cm)
    real(kind=f)                       :: rmass(NBIN)           ! bin mass (g)
    real(kind=f)                       :: sub_d                 ! integration substep diameter (um)
    real(kind=f)                       :: sub_dd                ! integration substep width (um)
    real(kind=f)                       :: snow_d                ! starting snow diameter (um)
    real(kind=f)                       :: nsnow                 ! number of snow particles (um)
    real(kind=f)                       :: snow_r3               ! snow N*r^3 (um^3)
    real(kind=f)                       :: snow_r2               ! snow N*r^3 (um^2)
    real(kind=f)                       :: snow_reff(NDTEMP)     ! snow effective radius (cm)
    real(kind=f)                       :: eshape                ! particle aspect ratio (> 0 is prolate)
    real(kind=f)                       :: shapeFactor           ! shape factor for maximum radius
    real(kind=f)                       :: remainder
    real(kind=f)                       :: lambda                ! fit factor for H&S 2010 size distribution
    real(kind=f)                       :: temp                  ! temperature (C)
     
    ! Default return code.
    rc = 0
    
    call CARMA_Get(carma, rc, do_print_init=do_print_init, LUNOPRT=LUNOPRT, do_grow=do_grow, do_detrain=do_detrain, do_thermo=do_thermo)
    if (rc < RC_OK) call endrun('CARMA_CheckMassAndEnergy::CARMA_Get failed.') 
    
    ! Lookup indices to other constituents that are needed.
    call cnst_get_ind('CLDICE', ixcldice)
    call cnst_get_ind('NUMICE', ixnumice)
    call cnst_get_ind('CLDLIQ', ixcldliq)
    call cnst_get_ind('NUMLIQ', ixnumliq)
    
    ! Add the CAM ice and liquid fields as some that could be modified by CARMA.
    lq_carma(ixcldice) = .true.
    lq_carma(ixnumice) = .true.
    lq_carma(ixcldliq) = .true.
    lq_carma(ixnumliq) = .true.      

    if (do_print_init) then
      write(LUNOPRT,*) ""
      write(LUNOPRT,*) "Initializing CARMA Detrainment"
      write(LUNOPRT,*) ""
      write(LUNOPRT,*) "  Using ice method = ", carma_dice_method
    end if
    
    ! For detrainment of ice, setup the fractions of ice that go into each bin and
    ! into snow. This can be done different ways:
    !
    !   - monodisperse
    !   - temperature dependent size distribution
    !
    ! In any of these, a fraction of the ice can go directly to snow, rather than
    ! going into bins first.
    !
    ! Puts all of the detraining cloud water from convection into the large scale cloud,
    ! and puts detraining cloud water into liquid and ice based on temperature partition
    call CARMAGROUP_Get(carma, I_GRP_CRDICE, rc, r=r(:), dr=dr(:), eshape=eshape, rmass=rmass(:))
    if (rc < RC_OK) call endrun('CARMA_InitializeModel::CARMAGROUP_Get failed.')

    ! This size distribution in based upon the maximum diameter, so the ice particles have
    ! a shape then pass the largest dimension to the size distribution.
    !
    ! NOTE: This is assuming the shape is a spheroid. Should consider passing shape
    ! parameters out of setupvfall, so that f1 is available for this.
    if (eshape >= 1._f) then
      shapeFactor = eshape**(1._f / 3._f)
    else
      shapeFactor = eshape**(- 1._f / 3._f)
    end if
        
    dice_bin_fraction(:, :) = 0._f
    dice_snow_fraction(:)   = 0._f
    dice_snow_rmass(:)      = 0._f
    
    
    ! Heymsfield & Schmitt [2010] tmeperature dependent distribution
    if (carma_dice_method == "dist_hym2010") then
    
      ! Integrate over the defined temperaure range.
      do itemp = 1, NDTEMP
      
        temp   = 1._f - itemp       

        ! Determine the exponentianal factor of the number distribution from H&S eq. 7.
        lambda = dist_hym2010_alpha * exp(temp * dist_hym2010_beta)

        ! Determine a mass distribution using from a size distribution using this
        ! lambda. The number distribution is N = N0 * exp(-lambda * D), with D in
        ! cm from H&S eq. 1. Since this is just to generate a PDF, just use an N0
        ! of 1.
        !
        ! NOTE: This mass distribution (dMdD) is based on the diameter in cm.
        do ibin = 1, NBIN
      
          ! Determine the fraction in each bin.
          !
          ! NOTE: The bins are wide realtive to this function, so sum over an interval
          sub_dd = 2._f * dr(ibin) * shapeFactor
          sub_d =  2._f * r(ibin)  * shapeFactor

          dice_bin_fraction(ibin, itemp) = dice_bin_fraction(ibin, itemp) + &
            rmass(ibin) / lambda * &
            (exp(-lambda * (sub_d - (sub_dd / 2._f))) - exp(-lambda * (sub_d + (sub_dd / 2._f))))
        end do
 
        ! Integrate to determine how much mass exits outside of the bins.
        ! Now integrate the snow distribution. We know the snow amount, but need an effective radius
        ! to determine the snow number.
        sub_d  = 2._f * (r(NBIN) + (dr(NBIN) / 2._f)) * shapeFactor
        sub_dd = (snow_max_d * 1e-4 - sub_d) / NINTS_SNOW
        sub_d  = sub_d + sub_dd / 2._f
        
        remainder = 0._f
        
        do i = 1, NINTS_SNOW
        
          ! Determine the number.
          !
          ! NOTE: Use the unscaled diameter and assume a sphere to get the volume of the particle.
          nsnow = exp(-lambda * (sub_d - (sub_dd / 2._f))) - exp(-lambda * (sub_d + (sub_dd / 2._f))) * sub_dd
  
          ! Assume density from Heymsfield & Schmitt [2010]. This assumes that:
          !
          !   m = aD^2.1
          !
          ! NOTE: This needs to match the density assumption made in the detrained ice bins.
          remainder = remainder + nsnow / lambda * 4.22e-3_f * (sub_d**2.1)
  
          sub_d = sub_d + sub_dd
        end do
        
        ! The sum of the integral may not be exactly 1, so scale the total so as not to skew
        ! the amount going straight to snow.
        dice_bin_fraction(:, itemp) = dice_bin_fraction(:, itemp) / (sum(dice_bin_fraction(:, itemp)) + remainder)
      
        
        ! Now integrate the snow distribution. We know the snow amount, but need an effective radius
        ! to determine the snow number.
        snow_d = 2._f * ((r(NBIN) + dr(NBIN) / 2._f)) 
        sub_dd = (snow_max_d * 1e-4 - snow_d) / NINTS_SNOW
        sub_d  = snow_d + (sub_dd / 2._f)
       
        snow_r3 = 0._f
        snow_r2 = 0._f
        
        do i = 1, NINTS_SNOW
        
          ! Determine the number.
          !
          ! NOTE: Use the unscaled diameter and assume a sphere to get the volume of the particle.
          nsnow = exp(-lambda * (sub_d - (sub_dd / 2._f))) - exp(-lambda * (sub_d + (sub_dd / 2._f))) * sub_dd
            
          snow_r3 = snow_r3 + nsnow / lambda * (sub_d / 2._f)**3
          snow_r2 = snow_r2 + nsnow / lambda * (sub_d / 2._f)**2
  
          sub_d = sub_d + sub_dd
        end do
        
        if (snow_r2 <= 0._f) then
          snow_reff(itemp) = 0.1_f
        else
          snow_reff(itemp) = snow_r3 / snow_r2
        end if
        
        ! If autoconversion is on then, detrain directly to snow. Otherwise, add the extra
        ! mass to the largest bin.
        if (carma_do_autosnow) then
          dice_snow_fraction(itemp) = 1._f - sum(dice_bin_fraction(:, itemp))
        else
          dice_snow_fraction(itemp) = 0._f
          dice_bin_fraction(NBIN, itemp) = dice_bin_fraction(NBIN, itemp) + 1._f - sum(dice_bin_fraction(:, itemp))
        end if

        ! The sum of the integral may not be exactly 1, so scale the total so as not to skew
        ! the amount going straight to snow. 
        dice_bin_fraction(:, itemp) = dice_bin_fraction(:, itemp) / sum(dice_bin_fraction(:, itemp))          
      end do
    
    ! Default to monodisperse
    else
    
      do ibin = 1, NBIN
        if (r(ibin) >= r_dice_mono) then
          dice_bin_fraction(ibin, :) = 1._f - dice_loss
          
          exit
        end if
      end do
       
      dice_snow_fraction(:) = 1._f - sum(dice_bin_fraction(:, 1))
      snow_reff(:) = dice_snow_reff_mono
    end if
    
    
    ! Determine the amount that goes into snow.
    dice_snow_rmass(:)    = 4._f / 3._f * PI * (snow_reff(:)**3) * CAM_RHOSN / 1e3_f
    
    if (do_print_init) then
      do itemp = 1, NDTEMP, 10
      
        if ((itemp == 1) .or. (carma_dice_method == "dist_hym2010")) then
        
          if (carma_dice_method == "dist_hym2010") then
            write(LUNOPRT,*) ""
            write(LUNOPRT,*) "  Temperature = ", 1 - itemp, " C"
            write(LUNOPRT,*) ""
          end if
        
          write(LUNOPRT,*) ""
          write(LUNOPRT,*) "        ibin       r (um)                   fraction"

          do ibin = 1, NBIN
            write(LUNOPRT,*) ibin, r(ibin)*1e4_f, dice_bin_fraction(ibin, itemp)
          end do
      
          write(LUNOPRT,*) ""
          write(LUNOPRT,*) "  Total fractions"
          write(LUNOPRT,*) "    ice  = ", 1._f - dice_snow_fraction(itemp)
          write(LUNOPRT,*) "    snow = ", dice_snow_fraction(itemp)

          write(LUNOPRT,*) ""
          write(LUNOPRT,*) "  Snow"
          write(LUNOPRT,*) "    min_r (um) = ", snow_d / 2._f
          write(LUNOPRT,*) "    rmass (kg) = ", dice_snow_rmass(itemp)
          write(LUNOPRT,*) "    reff  (um) = ", snow_reff(itemp)*1e4_f
          write(LUNOPRT,*) ""
        end if
      end do
    end if
    
    ! Log a warning message if doing growth or detrainment and not doing
    ! thermodynamics. This will cause an energy error to be reported by CAM.
    if ((do_grow .or. do_detrain) .and. .not. do_thermo) then
      if (do_print_init) write(LUNOPRT,*) "CARMA_InitializeModel:WARNING - do_grow and/or do_detrain are selected without do_thermo which may result in energy conservation errors."
    end if

    return
  end subroutine CARMA_InitializeModel


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
    use phys_grid,     only: get_lon_all_p, get_lat_all_p, get_rlat_all_p
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

    integer      :: lat(pcols)              ! latitude index 
    integer      :: lon(pcols)              ! longitude index
    real(r8)     :: clat(pcols)             ! latitude 
    integer      :: lchnk                   ! chunk identifier
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

    ! Determine the latitude and longitude of each column.
    ncol = state%ncol    
    
    ! Add any surface flux here.
    surfaceFlux(:ncol) = 0.0_r8
    
    ! For emissions into the atmosphere, put the emission here.
    tendency(:ncol, :pver) = 0.0_r8
    
    return
  end subroutine CARMA_EmitParticle


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
    
    type(carma_type), intent(in)       :: carma              !! the carma object
    integer, intent(in)                :: ielem              !! element index
    integer, intent(in)                :: ibin               !! bin index
    real(r8), intent(inout)            :: q(:, :)            !! mass mixing ratio (gcol, lev)
    integer, intent(in)                :: gcid(:)            !! global column id
    integer, intent(out)               :: rc                 !! return code, negative indicates failure
        

    ! Default return code.
    rc = RC_OK

    ! Add initial condition here.
    
    return
  end subroutine CARMA_InitializeParticle
  
  
  !! This routine is used to determine the shape parameters (pgam and lamc) for
  !! cloud ice.
  !!
  !! This code is taken from cldwat2m.F90, and ideally, there would be a routine
  !! in the cldwat2m module available for this purpose rather than duplicating the
  !! code and the parameters here.
  subroutine CARMA_GetGammaParmsForIce(carma, qiic, niic, rho, pgam, lami, rc)
    use shr_spfn_mod, only           : gamma => shr_spfn_gamma_nonintrinsic

    implicit none
    
    type(carma_type), intent(in)       :: carma         !! the carma object
    real(r8), intent(in)               :: qiic(pver)    !! in-cloud cloud liquid mixing ratio
    real(r8), intent(in)               :: niic(pver)    !! in-cloud droplet number conc
    real(r8), intent(in)               :: rho(pver)     !! air density (kg m-3)
    real(r8), intent(out)              :: pgam(pver)    !! spectral width parameter of droplet size distr
    real(r8), intent(out)              :: lami(pver)    !! slope of cloud liquid size distr
    integer, intent(out)               :: rc            !! return code, negative indicates failure


    real(r8), parameter                :: qsmall     = 1.e-36_r8     ! min mixing ratio 
!    real(r8), parameter                :: qsmall     = 1.e-18_r8     ! min mixing ratio 
    real(r8), parameter                :: pi         = 3.1415927_r8
    real(r8), parameter                :: dcs        = 250.e-6_r8    ! autoconversion size threshold for cloud ice to snow (m)
    real(r8), parameter                :: rhoi       = 500._r8       ! bulk density ice

    ! cloud ice mass-diameter relationship
    real(r8), parameter                :: ci = rhoi*pi/6._r8  
    real(r8), parameter                :: di = 3._r8
    
    integer                            :: k
    real(r8)                           :: n0i(pver)     ! intercept of cloud ice size distr
    real(r8)                           :: lammax        ! maximum allowed slope of size distr
    real(r8)                           :: lammin        ! minimum allowed slope of size distr
    real(r8)                           :: nc(pver)      ! in-cloud droplet number conc


    ! Default return code.
    rc = RC_OK
    

   !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ! get size distribution parameters based on in-cloud cloud water/ice 
    ! these calculations also ensure consistency between number and mixing ratio
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    
    ! NOTE: For ice, pgam is assumed to be 0.
    pgam(:) = 0._r8

    ! check for slope
    lammax = 1._r8/10.e-6_r8
    lammin = 1._r8/(2._r8*dcs)

    do k = 1, pver
  
      if (qiic(k).ge.qsmall) then

        ! add upper limit to in-cloud number concentration to prevent numerical error
        nc(k)=min(niic(k),qiic(k)*1.e20_r8)
  
        lami(k) = (gamma(1._r8+di)*ci*nc(k)/qiic(k))**(1._r8/di)
        n0i(k) = nc(k)*lami(k)

        ! adjust vars
        if (lami(k).lt.lammin) then
  
          lami(k) = lammin
          n0i(k) = lami(k)**(di+1._r8)*qiic(k)/(ci*gamma(1._r8+di))
          nc(k) = n0i(k)/lami(k)
        else if (lami(k).gt.lammax) then
          lami(k) = lammax
          n0i(k) = lami(k)**(di+1._r8)*qiic(k)/(ci*gamma(1._r8+di))
          nc(k) = n0i(k)/lami(k)
        end if
      else
        lami(k) = 0._r8
        n0i(k)  = 0._r8
      end if
    end do

    return
  end subroutine CARMA_GetGammaParmsForIce
  
    
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

  
  !! This routine is used to determine the shape parameters (pgam and lamc) for
  !! cloud water.
  !!
  !! This code is taken from cldwat2m.F90, and ideally, there would be a routine
  !! in the cldwat2m module available for this purpose rather than duplicating the
  !! code and the parameters here.
  subroutine CARMA_GetGammaParmsForLiq(carma, qcic, ncic, rho, pgam, lamc, rc)
    use shr_spfn_mod, only           : gamma => shr_spfn_gamma_nonintrinsic

    implicit none
    
    type(carma_type), intent(in)       :: carma         !! the carma object
    real(r8), intent(in)               :: qcic(pver)    !! in-cloud cloud liquid mixing ratio
    real(r8), intent(in)               :: ncic(pver)    !! in-cloud droplet number conc
    real(r8), intent(in)               :: rho(pver)     !! air density (kg m-3)
    real(r8), intent(out)              :: pgam(pver)    !! spectral width parameter of droplet size distr
    real(r8), intent(out)              :: lamc(pver)    !! slope of cloud liquid size distr
    integer, intent(out)               :: rc            !! return code, negative indicates failure


    real(r8), parameter                :: rhow       = 1000._r8      ! bulk density liquid (kg/m3)
!    real(r8), parameter                :: qsmall     = 1.e-18_r8     ! min mixing ratio 
    real(r8), parameter                :: qsmall     = 1.e-36_r8     ! min mixing ratio 
    real(r8), parameter                :: pi         = 3.1415927_r8
    real(r8), parameter                :: cdnl       = 0.e6_r8    ! cloud droplet number limiter
    
    integer                            :: k
    real(r8)                           :: n0c(pver)     ! intercept of cloud liquid size distr
    real(r8)                           :: lams(pver)    ! slope of snow size distr
    real(r8)                           :: n0s(pver)     ! intercept of snow size distr
    real(r8)                           :: lamr(pver)    ! slope of rain size distr
    real(r8)                           :: n0r(pver)     ! intercept of rain size distr
    real(r8)                           :: lammax        ! maximum allowed slope of size distr
    real(r8)                           :: lammin        ! minimum allowed slope of size distr
    real(r8)                           :: cdist1(pver)  ! size distr parameter to calculate droplet freezing
    real(r8)                           :: nc(pver)      ! in-cloud droplet number conc


    ! Default return code.
    rc = RC_OK
    

   !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ! get size distribution parameters based on in-cloud cloud water/ice 
    ! these calculations also ensure consistency between number and mixing ratio
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    do k = 1, pver
  
      if (qcic(k).ge.qsmall) then
  
        ! add upper limit to in-cloud number concentration to prevent numerical error
        nc(k) = min(ncic(k),qcic(k)*1.e20_r8)
        nc(k)=max(nc(k),cdnl/rho(k)) ! sghan minimum in #/cm  
    
        ! get pgam from fit to observations of martin et al. 1994    
        pgam(k) = 0.0005714_r8*(nc(k)/1.e6_r8*rho(k))+0.2714_r8
        pgam(k) = 1._r8/(pgam(k)**2)-1._r8
        pgam(k) = max(pgam(k),2._r8)
        pgam(k) = min(pgam(k),15._r8)
    
        ! calculate lamc
        lamc(k) = (pi/6._r8*rhow*nc(k)*gamma(pgam(k)+4._r8) / &
                  (qcic(k)*gamma(pgam(k)+1._r8)))**(1._r8/3._r8)
    
        ! lammin, 50 micron diameter max mean size
        lammin = (pgam(k)+1._r8)/50.e-6_r8
        lammax = (pgam(k)+1._r8)/2.e-6_r8
    
        if (lamc(k).lt.lammin) then
          lamc(k) = lammin
          nc(k) = 6._r8*lamc(k)**3*qcic(k)* &
                    gamma(pgam(k)+1._r8)/ &
                    (pi*rhow*gamma(pgam(k)+4._r8))
        else if (lamc(k).gt.lammax) then
          lamc(k) = lammax
          nc(k) = 6._r8*lamc(k)**3*qcic(k)* &
                    gamma(pgam(k)+1._r8)/ &
                    (pi*rhow*gamma(pgam(k)+4._r8))
        end if    
    
        ! parameter to calculate droplet freezing
        cdist1(k) = nc(k)/gamma(pgam(k)+1._r8) 
        
      else
        pgam(k)   = 0._r8
        lamc(k)   = 0._r8
        cdist1(k) = 0._r8
      end if
    end do
    
    return
  end subroutine CARMA_GetGammaParmsForLiq
  
  
  ! Using the specified parameters for the gamma distribution, determine the mass mixing ratio of particles
  subroutine CARMA_GetMmrFromGamma(carma, r, dr, rmass, qic, nic, mu, lambda, mmr, rc)
    use shr_spfn_mod, only           : gamma => shr_spfn_gamma_nonintrinsic

    implicit none
    
    type(carma_type), intent(in)       :: carma           !! the carma object
    real(kind=f), intent(in)           :: r(NBIN)         !! bin mean radius
    real(kind=f), intent(in)           :: dr(NBIN)        !! bin radius width
    real(kind=f), intent(in)           :: rmass(NBIN)     !! bin mass
    real(r8), intent(in)               :: qic(pver)       !! in-cloud cloud liquid mixing ratio
    real(r8), intent(in)               :: nic(pver)       !! in-cloud droplet number conc
    real(r8), intent(in)               :: mu(pver)        !! spectral width parameter of droplet size distr
    real(r8), intent(in)               :: lambda(pver)    !! slope of cloud liquid size distr
    real(r8), intent(out)              :: mmr(NBIN,pver)  !! elements mass mixing ratio
    integer, intent(out)               :: rc              !! return code, negative indicates failure
    
    integer                            :: k               ! z index
    integer                            :: ibin            ! bin index
    real(kind=f)                       :: totalMass       ! mmr of all particles (kg/kg)
    real(kind=f)                       :: n               ! number of particles (#/kg)
    real(kind=f)                       :: n0              ! number parameter for gamma distribution
    real(kind=f)                       :: d(NBIN)         ! bin diameter (m)
    real(kind=f)                       :: dd(NBIN)        ! diameter width of bin (m)

!    real(r8), parameter                :: qsmall     = 1.e-18_r8     ! min mixing ratio 
    real(r8), parameter                :: qsmall     = 1.e-36_r8     ! min mixing ratio 

  
    ! Default return code.
    rc = RC_OK
    
    ! Their equations are in terms of diameter (in m)
    d(:)  = 2._r8 * r(:) * 1e-2_r8
    dd(:) = 2._r8 * dr(:) * 1e-2_r8

    do k = 1, pver
    
      ! From Morisson & Gettelman [2008] and cldwat2m
      !
      ! If there is a small mass, then ther are no particles.
      if (qic(k) < qsmall) then
        mmr(:, k) = 0._r8
      else
        n0 = (nic(k) * (lambda(k) ** (mu(k) + 1._r8)) / (gamma(mu(k) + 1._r8)))    
  
      
        ! Iterate over the bins.
        !
        ! NOTE: Just the functional fit can go negative for some bins with larger diameter, but this is not physical.
        do ibin = 1, NBIN
          n = n0 * (d(ibin)**mu(k)) * exp(-lambda(k) * d(ibin)) * dd(ibin)
          mmr(ibin, k) = n * rmass(ibin) * 1e-3_r8
        end do
      
        ! Adjust the number density so that we don't create mass. This will adjust for
        ! problems fitting the size distribution and for differences in the assumptions
        ! of the bulk density of the particles.
        totalMass = sum(mmr(:, k))
        if (totalMass /= 0._r8) then
          mmr(:, k) = mmr(:, k) * (qic(k) / totalMass)
        else
          mmr(:, k) = 0._r8
        end if
      end if
    end do

    return
  end subroutine CARMA_GetMmrFromGamma


  !! Detemrine the total cloud ice concentration and number stored in the bins that represent
  !! water within the CARMA model.
  !!
  !! For snow, it is assumed that the largest ice bin in the in situ and detrained ice are
  !! snow. The mass of these bins is the same, but the dimensions are different since there
  !! are different shape assumptions for the different types.
  !!
  !!  @version Nov-2009 
  !!  @author  Chuck Bardeen 
  subroutine CARMA_GetTotalIceAndSnow(carma, cstate, makeSnow, iceMass, iceNumber, snowMass, snowNumber, snowSurface, rc, iceRe)
    implicit none
    
    type(carma_type), intent(in)         :: carma     !! the carma object
    type(carmastate_type), intent(inout) :: cstate    !! the carma state object
    logical, intent(in)                  :: makeSnow  !! should bins be changed because of snow?
    real(kind=f), intent(out)            :: iceMass(pver)      !! ice mass mixing ratio (kg/kg)
    real(kind=f), intent(out)            :: iceNumber(pver)    !! ice number mixing ratio (#/kg)
    real(kind=f), intent(out)            :: snowMass(pver)     !! snow mass mixing ratio (kg/kg)
    real(kind=f), intent(out)            :: snowNumber(pver)   !! snow number (#/kg)
    real(kind=f), intent(out)            :: snowSurface        !! snow on surface (kg/m2)
    integer, intent(out)                 :: rc        !! return code, negative indicates failure
    real(kind=f), intent(out), optional  :: iceRe(pver)   !! ice effective radius (m)
 
    integer                              :: LUNOPRT              ! logical unit number for output
    logical                              :: do_print             ! do print output?

    integer                              :: igroup    ! group index
    integer                              :: ielem     ! element index
    integer                              :: ibin      ! bin index
    integer                              :: iz        ! vertical index
    integer                              :: icore     ! core index
    integer                              :: icorelem(NELEM) ! core indexes for group
    integer                              :: ncore     ! number of core elements
    integer                              :: maxbin    ! maximum prognostic bin
    
    real(kind=f)                         :: coreMass(pver)     ! core mass mixing ratio (kg/kg)
    real(kind=f)                         :: coreSurface        ! core on surface (kg/kg)

    real(kind=f)                         :: newSnow            ! [Total] sfc flux of snow from stratiform (m/s)

    real(kind=f)                         :: mmr(pver)          ! mass mixing ratio (#/kg)
    real(kind=f)                         :: mmrcore(pver)      ! core mass mixing ratio (#/kg)
    real(kind=f)                         :: nmr(pver)          ! number mixing ratio (#/kg)
    real(kind=f)                         :: r(NBIN)            ! radius (cm)
    real(kind=f)                         :: rmass(NBIN)        ! mass (g)
    real(kind=f)                         :: rrat(NBIN)         ! particle maximum radius ratio ()
    real(kind=f)                         :: arat(NBIN)         ! particle area ration ()
    real(kind=f)                         :: sfc                ! surface mass (kg/m2)
    real(kind=f)                         :: sfccore            ! core surface mass (kg/m2)
    real(kind=f)                         :: nd(pver)           ! number density (#/cm3)
    real(kind=f)                         :: pa(pver)           ! projected area (cm2)
    real(kind=f)                         :: md(pver)           ! mass density (g/cm3)


    ! Default return code.
    rc = RC_OK
    
    call CARMA_Get(carma, rc, do_print=do_print, LUNOPRT=LUNOPRT)
    if (rc < RC_OK) call endrun('CARMA_CheckMassAndEnergy::CARMA_Get failed.') 
    
    iceMass(:)     = 0._f
    iceNumber(:)   = 0._f
    snowMass(:)    = 0._f
    snowNumber(:)  = 0._f
    snowSurface    = 0._f
    pa(:)          = 0._f
    md(:)          = 0._f

    if (present(iceRe)) iceRe(:) = 0._f
      
      
    ! Detrained Ice, Aged
    igroup = I_GRP_CRDICE
    ielem  = I_ELEM_CRDICE
    
    call CARMAGROUP_Get(carma, igroup, rc, r=r, rmass=rmass, arat=arat, rrat=rrat, maxbin=maxbin)
    if (rc < RC_OK) call endrun('GetTotalIceAndSnow::CARMAGROUP_Get failed.')
    
    do ibin = 1, NBIN
      call CARMASTATE_GetBin(cstate, ielem, ibin, mmr, rc, nmr=nmr, surface=sfc, numberDensity=nd)
      if (rc < RC_OK) call endrun('GetTotalIceAndSnow::CARMASTATE_GetBin failed.')

      ! Only calculate snow if CARMA is responsible for the cloud ice.
      if (carma_do_cldice .and. carma_do_autosnow .and. (ibin > maxbin)) then
        snowMass(:)    = snowMass(:)   + mmr(:)
        snowNumber(:)  = snowNumber(:) + nmr(:)
      
        if (makeSnow) then
          ! This ice is now snow, so zero it out in the ice bins.
          mmr(:) = 0._f
 
          call CARMASTATE_SetBin(cstate, ielem, ibin, mmr, rc)
          if (rc < RC_OK) call endrun('GetTotalIceAndSnow::CARMASTATE_SetBin failed.')
        end if
      else
        iceMass(:)   = iceMass(:)   + mmr(:)
        iceNumber(:) = iceNumber(:) + nmr(:)

        where (nd(:) > SMALL_PC)
        
          ! NOTE: This is following the definition of Dave Mitchell for effective diameter,
          ! Mitchell [2002], which indicates it needs to be scaled based on the effective
          ! ice density.
          pa(:) = pa(:) + nd(:) * PI * ((r(ibin) * rrat(ibin))**2) * arat(ibin)
          md(:)  = md(:)  + nd(:) * rmass(ibin)
        end where
      end if
    
      ! The particles that sedimented out of the bottom layer need to be included
      ! in the mass of snow.
      snowSurface = snowSurface + sfc
    end do


    ! Detrained Ice, Fresh
  
    do ibin = 1, NBIN
      call CARMASTATE_GetDetrain(cstate, ielem, ibin, mmr, rc, nmr=nmr, numberDensity=nd)
      if (rc < RC_OK) call endrun('GetTotalIceAndSnow::CARMASTATE_GetBin failed.')

      ! Only calculate snow if CARMA is responsible for the cloud ice.
      if (carma_do_cldice .and. carma_do_autosnow .and. (ibin > maxbin)) then
        snowMass(:)    = snowMass(:)   + mmr(:)
        snowNumber(:)  = snowNumber(:) + nmr(:)
      
        if (makeSnow) then
          ! This ice is now snow, so zero it out in the ice bins.
          mmr(:) = 0._f
 
          call CARMASTATE_SetDetrain(cstate, ielem, ibin, mmr, rc)
          if (rc < RC_OK) call endrun('GetTotalIceAndSnow::CARMASTATE_SetBin failed.')
        end if
      else
        iceMass(:)   = iceMass(:)   + mmr(:)
        iceNumber(:) = iceNumber(:) + nmr(:)

        where (nd(:) > SMALL_PC)
          pa(:) = pa(:) + nd(:) * PI * ((r(ibin) * rrat(ibin))**2) * arat(ibin)
          md(:)  = md(:)  + nd(:) * rmass(ibin)
        end where
      end if
    end do
    
    
    ! In-situ ice.
    igroup = I_GRP_CRSICE
    ielem  = I_ELEM_CRSICE
    
    call CARMAGROUP_Get(carma, igroup, rc, r=r, ncore=ncore, icorelem=icorelem, arat=arat, rrat=rrat, maxbin=maxbin)
    if (rc < RC_OK) call endrun('GetTotalIceAndSnow::CARMAGROUP_Get failed.')
    
    do ibin = 1, NBIN
      call CARMASTATE_GetBin(cstate, ielem, ibin, mmr, rc, nmr=nmr, surface=sfc, numberDensity=nd)
      if (rc < RC_OK) call endrun('CARMA_DiagnoseBulk::CARMASTATE_GetBin failed.')
      
      ! Determine how much of the mmr is related to core mass. This needs to
      ! be subtracted to get the amount of water in the ice.
      coreMass(:) = 0.0_f
      coreSurface = 0.0_f
      
      do icore = 1, ncore
        call CARMASTATE_GetBin(cstate, icorelem(icore), ibin, mmrcore, rc, surface=sfccore)
        if (rc < RC_OK) call endrun('GetTotalIceAndSnow::CARMASTATE_GetBin failed.')
      
        coreMass(:) = coreMass(:) + mmrcore(:)
        coreSurface = coreSurface + sfccore
      end do
      
      ! The core mass can't be more than the particle mass. If so, this indicates
      ! that are problem happened, perhaps during advection and the particle masses
      ! should be ignored. This should never happen from CARMA itself.
      if (carma_do_mass_fix) then
        do iz = 1, pver
        
          if (coreMass(iz) > mmr(iz)) then
            if (carma_do_mass_fix) then
            
              if (carma_do_print_fix .and. do_print) write(LUNOPRT,*) &
                 "  GetTotalIceAndSnow::WARNING - Adjusting particle for core mass error", &
                 iz, ielem, ibin, mmr(iz), coreMass(iz)

              ! It is hard to know what the right fix should be. You could reset
              ! the particle mass to the coremass, but this will create lots of
              ! small particles. It may be safer just to zero out both the particle
              ! count and all of the core masses, assuming that this is a particle
              ! that was created by diffusion in the transport and shouldn't really exist.
              mmr(iz) = coreMass(iz)
            end if
          end if
        end do
            
        call CARMASTATE_SetBin(cstate, ielem, ibin, mmr, rc)
        if (rc < RC_OK) call endrun('GetTotalIceAndSnow::CARMASTATE_SetBin failed.')
      end if
            
      ! Only calculate snow if CARMA is responsible for the cloud ice.
      if (carma_do_cldice .and. carma_do_autosnow .and. (ibin > maxbin)) then
        snowMass(:)    = snowMass(:)   + mmr(:) - coreMass(:)
        snowNumber(:)  = snowNumber(:) + nmr(:)
      
        if (makeSnow) then
        
          ! This ice is now snow, so zero it out in the ice bins.
          mmr(:) = 0._f
        
          call CARMASTATE_SetBin(cstate, ielem, ibin, mmr, rc)
          if (rc < RC_OK) call endrun('GetTotalIceAndSnow::CARMASTATE_SetBin failed.')
        
          ! Also zero out the core mass.
          !
          ! In the future, you could try to keep track of the mass of the cores
          ! in the snow and communicate that to CAM.
          do icore = 1, ncore
            call CARMASTATE_SetBin(cstate, icorelem(icore), ibin, mmr, rc)
            if (rc < RC_OK) call endrun('GetTotalIceAndSnow::CARMASTATE_SetBin failed.')
          end do
        end if
      else
        iceMass(:)   = iceMass(:)   + mmr(:) - coreMass(:)
        iceNumber(:) = iceNumber(:) + nmr(:)

        where (nd(:) > SMALL_PC)
          pa(:) = pa(:) + nd(:) * PI * ((r(ibin) * rrat(ibin))**2) * arat(ibin)
          md(:)  = md(:)  + nd(:) * rmass(ibin)
        end where
      end if
    
    
      ! The particles that sedimented out of the bottom layer need to be included
      ! in the mass of snow.
      snowSurface = snowSurface + sfc - sfccore
      
      
      ! Calculate the effective radius (total volume / total area).
      ! NOTE: cm -> m.
      if (present(iceRe)) then
        where (pa(:) > 0.0_r8)
          iceRe(:) = (3._f / 4._f) * (md(:)  / (0.917_f * pa(:))) * 1e-2_f
        end where
        
      end if
    end do

    return
  end subroutine CARMA_GetTotalIceAndSnow
    
    
  !! Detemrine the total cloud water concentration and number stored in the bins that represent
  !! water within the CARMA model.
  !!
  !!  @version Nov-2009 
  !!  @author  Chuck Bardeen 
  subroutine CARMA_GetTotalWaterAndRain(carma, cstate, waterMass, waterNumber, rainSurface, rc)    
    implicit none
    
    type(carma_type), intent(in)         :: carma              !! the carma object
    type(carmastate_type), intent(inout) :: cstate             !! the carma state object
    real(kind=f), intent(out)            :: waterMass(pver)    !! water mass mixing ratio (kg/kg)
    real(kind=f), intent(out)            :: waterNumber(pver)  !! water number mixing ratio (#/kg)
    real(kind=f), intent(out)            :: rainSurface        !! rain on surface (kg/m2)
    integer, intent(out)                 :: rc                 !! return code, negative indicates failure
 
    integer                              :: igroup             ! group index
    integer                              :: ielem              ! element index
    integer                              :: ibin               ! bin index
    
    real(kind=f)                         :: mmr(pver)          ! mass mixing ratio (#/kg)
    real(kind=f)                         :: nmr(pver)          ! number mixing ratio (#/kg)
    real(kind=f)                         :: sfc                ! surface mass (kg/m2)

    rc = RC_OK

    waterMass(:)    = 0._f
    waterNumber(:)  = 0._f
    rainSurface     = 0._f
    
    igroup = I_GRP_CRLIQ
    ielem  = I_ELEM_CRLIQ
    
    do ibin = 1, NBIN
      call CARMASTATE_GetBin(cstate, ielem, ibin, mmr, rc, nmr=nmr, surface=sfc)
      if (rc < RC_OK) call endrun('CARMA_GetTotalWaterAndRain::CARMASTATE_GetBin failed.')
      
      waterMass(:)   = waterMass(:)   + mmr(:)
      waterNumber(:) = waterNumber(:) + nmr(:)
      
      ! The particles that sedimented out of the bottom layer need to be included
      ! in the mass of rain.
      rainSurface = rainSurface + sfc

      ! Include the detrained liquid that hasn't been added to the particle bins yet.
      call CARMASTATE_GetDetrain(cstate, ielem, ibin, mmr, rc, nmr=nmr)
      if (rc < RC_OK) call endrun('CARMA_GetTotalWaterAndRain::CARMASTATE_GetDetrain failed.')
      
      waterMass(:)   = waterMass(:)   + mmr(:)
      waterNumber(:) = waterNumber(:) + nmr(:)
    end do

    return
  end subroutine CARMA_GetTotalWaterAndRain
  


  subroutine CARMA_CheckMassAndEnergy(carma, cstate, madeSnow, name, state, icol, dt, rliq, prec_str, snow_str, waterMass, iceMass, snowMass, rc)
    implicit none
    
    type(carma_type), intent(in)         :: carma            !! the carma object
    type(carmastate_type), intent(inout) :: cstate           !! the carma state object
    logical, intent(in)                  :: madeSnow         !! should bins be changed because of snow?
    character*(*),intent(in)             :: name             !! test name
    type(physics_state), intent(in)      :: state            !! physics state variables
    integer, intent(in)                  :: icol             !! column index
    real(kind=f), intent(in)             :: dt               !! time step
    real(kind=f), intent(in)             :: rliq(pcols)      !! vertical integral of liquid not yet in q(ixcldliq)
    real(kind=f), intent(in)             :: prec_str(pcols)  !! [Total] sfc flux of precip from stratiform (m/s) 
    real(kind=f), intent(in)             :: snow_str(pcols)  !! [Total] sfc flux of snow from stratiform (m/s)
    real(kind=f), intent(in)             :: waterMass(pver)  !! water mass mixing ratio (kg/kg)
    real(kind=f), intent(in)             :: iceMass(pver)    !! ice mass mixing ratio (kg/kg)
    real(kind=f), intent(in)             :: snowMass(pver)   !! snow mass mixing ratio (kg/kg)
    integer, intent(out)                 :: rc               !! return code, negative indicates failure
 
 
    integer                              :: LUNOPRT              ! logical unit number for output
    logical                              :: do_print             ! do print output?
    logical                              :: do_detrain           ! do convective detrainment?

    real(kind=f)                         :: mmr(pver)          ! mass mixing ratio (#/kg)
    real(kind=f)                         :: totalMass
    real(kind=f)                         :: totalMass2

    real(r8)                             :: lat
    real(r8)                             :: lon
    

  1 format(/,'CARMA_CheckMassAndEnergy::ERROR - CARMA mass conservation error, ',a,',icol=',i4,',lat=',&
              f7.2,',lon=',f7.2,',cam=',e16.10,',carma=',e16.10,',diff=',e16.10,',rer=',e9.3)

    ! Default return code.
    rc = RC_OK
    
    call CARMA_Get(carma, rc, do_print=do_print, LUNOPRT=LUNOPRT, do_detrain=do_detrain)
    if (rc < RC_OK) call endrun('CARMA_CheckMassAndEnergy::CARMA_Get failed.') 
    
    totalMass = sum(state%q(icol, :, ixcldliq) * (state%pdel(icol, :) / gravit))
    totalMass = totalMass + sum(state%q(icol, :, ixcldice) * (state%pdel(icol, :) / gravit))
    totalMass = totalMass + sum(state%q(icol, :, 1) * (state%pdel(icol, :) / gravit))
        
    if (abs((totalMass - state%tw_cur(icol))) / state%tw_cur(icol) > 1e14_f) then
      if (do_print) write(LUNOPRT,*) "CARMA_CheckMassAndEnergy::WARNING Total water not conserved, ", totalMass, state%tw_cur, (totalMass - state%tw_cur(icol)), (totalMass - state%tw_cur(icol)) / state%tw_cur(icol)
    end if
    
    
    ! Get the total water coming out of CARMA
    call CARMASTATE_GetGas(cstate, I_GAS_H2O, mmr(:), rc)
    if (rc < RC_OK) call endrun('CARMA_CheckMassAndEnergy::CARMASTATE_GetGas failed.')

    totalMass2 = sum(waterMass(:) * (state%pdel(icol, :) / gravit))
    totalMass2 = totalMass2 + sum((iceMass(:)) * (state%pdel(icol, :) / gravit))
    
    ! If snow has been made, that means it has been removed from the cloud ice that is
    ! in the atmosphere and is now accounted for by the prec_str and snow_str fields.
    ! Prior to that, it is still considered as part of the atmospheric ice total.
    if (.not. madeSnow) then
      totalMass2 = totalMass2 + sum(snowMass(:) * (state%pdel(icol, :) / gravit))
    end if
    
    totalMass2 = totalMass2 + sum(mmr(:) * (state%pdel(icol, :) / gravit))
    totalMass2 = totalMass2 + prec_str(icol) * dt * 1000._f
    
    if (do_detrain) totalMass2 = totalMass2 + rliq(icol) * dt * 1000._f

    if (totalMass /= totalMass2) then

      if (totalMass /= 0._f) then

        if (abs((totalMass - totalMass2) / totalMass) > 1e-10_f)  then
          if (do_print) then
            call CARMASTATE_Get(cstate, rc, lat=lat, lon=lon)
            if (rc < RC_OK) call endrun('CARMA_DiagnoseBins::CARMASTATE_Get failed.')

            write(LUNOPRT,1) name, icol, lat, lon, totalMass, totalMass2, totalMass2-TotalMass, (totalMass - totalMass2) / totalMass

            write(LUNOPRT,*) "  state tw :  ", state%tw_cur(icol)
            write(LUNOPRT,*) ""
            write(LUNOPRT,*) "  old vap  :  ", sum(state%q(icol, :, 1) * (state%pdel(icol, :) / gravit))
            write(LUNOPRT,*) "  old liq  :  ", sum(state%q(icol, :, ixcldliq) * (state%pdel(icol, :) / gravit))
            write(LUNOPRT,*) "  old ice  :  ", sum(state%q(icol, :, ixcldice) * (state%pdel(icol, :) / gravit))
            write(LUNOPRT,*) ""
            write(LUNOPRT,*) "  new vap  :  ", sum(mmr(:) * (state%pdel(icol, :) / gravit))
            write(LUNOPRT,*) "  new liq  :  ", sum(waterMass(:) * (state%pdel(icol, :) / gravit))
            write(LUNOPRT,*) "  new ice  :  ", sum(iceMass(:) * (state%pdel(icol, :) / gravit))
            write(LUNOPRT,*) "  new snow :  ", sum(snowMass(:) * (state%pdel(icol, :) / gravit))
            write(LUNOPRT,*) "  rliq     :  ", rliq(icol) * dt * 1000._f
            write(LUNOPRT,*) "  prec_str :  ", prec_str(icol) * dt * 1000._f
          end if
        end if
      end if
    end if

    return
  end subroutine CARMA_CheckMassAndEnergy

end module
