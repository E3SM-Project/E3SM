
!---------------------------------------------------------------------------------
! this file is based on structure of gw_drag from WACCMX 3.5.48
!---------------------------------------------------------------------------------

module ionosphere

!---------------------------------------------------------------------------------
! Purpose:
!
! Module to compute the relevant physics and electrodynamics needed to improve and
! add a realistic simulation of the ionosphere.  Main modules are initialization,
! ion/electron temperature calculation, ambipolar diffusion, dynamo, .... 
!
! Authors: Joe McInerney/Hanli Liu/Art Richmond
!
!---------------------------------------------------------------------------------
  use shr_kind_mod,   only : r8 => shr_kind_r8            ! Real kind to declare variables
  use ppgrid,         only : pcols, pver, pverp           ! Dimensions and chunk bounds
  use cam_history,    only : outfld                       ! Routine to output fields to history files
  use physics_types,  only : physics_state, physics_ptend !Structures containing physics state and tendency variables
  use physics_buffer, only : pbuf_add_field, pbuf_get_index,dtype_r8, &
       physics_buffer_desc, pbuf_get_field ! Needed to get variables from physics buffer
  use cam_logfile,    only : iulog                        ! Output unit for run.out file
  use mo_jeuv,        only : nIonRates                    ! Number of ionization rates in mo_photo
  use shr_const_mod,  only : kboltz => shr_const_boltz, pi => shr_const_pi ! Boltzmann constant and pi

  implicit none

  save
  
  private   ! Make default type private to the module
!
! PUBLIC: interfaces
!
  public ionos_init     ! Initialization
  public ionos_register ! Registration of ionosphere variables in pbuf physics buffer
  public ionos_intr     ! interface to actual ionosphere simulation routines
  public tridag         ! Generic tridiagonal solver routine
!
! PUBLIC variables
!
!
!  Define public variables
!
!
! PRIVATE: Rest of the data and interfaces are private to this module
!   
  real(r8), parameter               :: kboltz_ev = 8.617E-5_r8 ! Boltzmann constant (eV/K)
  real(r8), parameter               :: temax = 7.0E3_r8        ! maximum electron temperature.
                
  integer                           :: ionBot                  ! bottom of ionosphere calculations
  integer                           :: ionBotP                 ! bottom of ionosphere calculations

  real(r8)                          :: calDay                  ! current calendar day
  real(r8)                          :: rads2Degs               ! radians to degrees
  real(r8)                          :: degs2Rads               ! degrees to radians

  real(r8)                          :: f107                    ! 10.7 cm solar flux

  type ionos_state

    real(r8), dimension(pcols)      :: cosZenAngR              ! cosine of zenith angle (radians)
    real(r8), dimension(pcols)      :: zenAngD                 ! zenith angle (degrees)

    real(r8), dimension(pcols,pver) :: bNorth3d  ! northward component of magnetic field units?
    real(r8), dimension(pcols,pver) :: bEast3d   ! eastward component of magnetic field
    real(r8), dimension(pcols,pver) :: bDown3d   ! downward component of magnetic field

    real(r8), dimension(pcols,pver,nIonRates) :: ionPRates    ! ionization rates temporary array (s-1 cm-3)
    real(r8), dimension(pcols,pver)           :: sumIonPRates ! Sum of ionization rates for O+,O2+,N+,N2+,NO+ (s-2 cm-3)

    real(r8), dimension(pcols,pver)  :: ue2d  ! horizontal x(eastward) component of ExB drift with added vertical dimension
    real(r8), dimension(pcols,pver)  :: ve2d  ! horizontal y(northward) component of ExB drift with added vertical dimension
    real(r8), dimension(pcols,pver)  :: we2d  ! vertical z component of ExB drift with added vertical dimension - midpoints

    real(r8), dimension(pcols,pverp) :: wei2d    ! vertical z component of ExB drift with added vertical dimension - interfaces

    real(r8), dimension(pcols,pverp) :: omegai   ! vertical velocity on interface levels (Pa/s)
 
    real(r8), dimension(pcols,pver)  :: dipMag   ! dip angle for each column (radians)
    real(r8), dimension(pcols,pver)  :: dipMagD  ! dip angle for each column (degrees)

    real(r8), dimension(pcols,pverp) :: tNInt    ! Interface Temperature (K)

    real(r8), dimension(pcols,pver)  :: ndensN2  ! N2 number density (cm-3)
    real(r8), dimension(pcols,pver)  :: ndensO2  ! O2 number density (cm-3)
    real(r8), dimension(pcols,pver)  :: ndensO1  ! O number density (cm-3)
    real(r8), dimension(pcols,pver)  :: ndensNO  ! NO number density (cm-3)
    real(r8), dimension(pcols,pver)  :: ndensN1  ! N number density  (cm-3)
    real(r8), dimension(pcols,pver)  :: ndensE   ! E electron number density (cm-3)
    real(r8), dimension(pcols,pver)  :: ndensOp  ! O plus number density (cm-3)
    real(r8), dimension(pcols,pver)  :: ndensO2p ! O2 plus ion number density (cm-3)
    real(r8), dimension(pcols,pver)  :: ndensNOp ! NO plus ion number density  (cm-3)

    real(r8), dimension(pcols,pver)  :: sourceg4 ! g4 source term for electron/ion temperature update

    real(r8), dimension(pcols,pverp) :: rairvi   ! Constituent dependent gas constant on interface levels

  end type ionos_state

contains

!==============================================================================

  subroutine ionos_init()
  
!-----------------------------------------------------------------------
! Time independent initialization for ionosphere simulation.
!-----------------------------------------------------------------------
    use cam_history,      only : addfld, add_default ! Routines and variables for adding fields to history output
    use dycore,           only : get_resolution
    use interpolate_data, only : lininterp
     
    implicit none

!-------------------------------------------------------------------------------
!  Add history variables for ionosphere 
!-------------------------------------------------------------------------------
    call addfld ('IONOS_NDENSN2',(/ 'lev' /), 'I', '1/m3','N2 Number Density-Ionos')
    call addfld ('IONOS_NDENSO2',(/ 'lev' /), 'I', '1/m3','O2 Number Density-Ionos')
    call addfld ('IONOS_NDENSO1',(/ 'lev' /), 'I', '1/m3','O1 Number Density-Ionos')
    call addfld ('IONOS_NDENSNO',(/ 'lev' /), 'I', '1/m3','NO Number Density-Ionos')
    call addfld ('IONOS_NDENSN1',(/ 'lev' /), 'I', '1/m3','NO Number Density-Ionos')
    call addfld ('IONOS_NDENSE' ,(/ 'lev' /), 'I', '1/m3','E Number Density-Ionos')
    call addfld ('IONOS_NDENSOP' ,(/ 'lev' /), 'I','1/m3','OP Number Density-Ionos')
    call addfld ('IONOS_NDENSO2P',(/ 'lev' /), 'I','1/m3','O2P Number Density-Ionos')
    call addfld ('IONOS_NDENSNOP',(/ 'lev' /), 'I','1/m3','NOP Number Density-Ionos')

    call addfld ('TE'           ,(/ 'lev' /), 'I', 'K','Electron Temperature')
    call addfld ('TI'           ,(/ 'lev' /), 'I', 'K','Ion Temperature')
    call addfld ('QIN'          ,(/ 'lev' /), 'I', 'J/kg/s','JOULE+IN Heating')
    call addfld ('LOSS_g3'     ,(/ 'lev' /)         , 'I', ' ','Loss Term g3')
    call addfld ('LOSS_EI'      ,(/ 'lev' /)         , 'I', ' ','Loss Term EI')
    call addfld ('LOSS_IN'      ,(/ 'lev' /)         , 'I', ' ','Loss Term IN')
    call addfld ('IONI_RATE'    ,(/ 'lev' /)         , 'I', ' ','Total ionization')
    call addfld ('IONI_Eff'     ,(/ 'lev' /)         , 'I', ' ','ionization efficiency')
    call addfld ('SOURCE_g4'    ,(/ 'lev' /)         , 'I', ' ','SOURCE g4')
    call addfld ('AUR_IRATESUM' ,(/ 'lev' /)         , 'I', ' ','Auroral ionization')

!-------------------------------------------------------------------------------
!  Set default values for ionosphere history variables
!-------------------------------------------------------------------------------
    call add_default ('IONOS_NDENSN2' , 1, ' ')
    call add_default ('IONOS_NDENSO2' , 1, ' ')
    call add_default ('IONOS_NDENSO1' , 1, ' ')
    call add_default ('IONOS_NDENSNO' , 1, ' ')
    call add_default ('IONOS_NDENSN1' , 1, ' ')
    call add_default ('IONOS_NDENSE'  , 1, ' ')
    call add_default ('IONOS_NDENSOP' , 1, ' ')
    call add_default ('IONOS_NDENSO2P', 1, ' ')
    call add_default ('IONOS_NDENSNOP', 1, ' ')
    call add_default ('TE'    , 1, ' ')
    call add_default ('TI'    , 1, ' ')
    call add_default ('QIN'   , 1, ' ')
    call add_default ('LOSS_g3'      , 1, ' ')
    call add_default ('LOSS_EI'       , 1, ' ')
    call add_default ('LOSS_IN'       , 1, ' ')
    call add_default ('IONI_RATE'     , 1, ' ')
    call add_default ('IONI_Eff'      , 1, ' ')
    call add_default ('SOURCE_g4'     , 1, ' ')
    call add_default ('AUR_IRATESUM'  , 1, ' ')
     
   return

  end subroutine ionos_init

!==============================================================================     

  subroutine ionos_register

!-----------------------------------------------------------------------
! Register ionosphere variables with physics buffer:
!
! Ion production rates pcols,pver,nIonRates,
!   so firstdim = 1 middledim = pver lastdim = nIonRates.
! 
! pcols dimension and lchnk assumed here
!
!-----------------------------------------------------------------------
     
    implicit none

    integer :: idx

    ! Electron temperature (global so can write to history files) 
    call pbuf_add_field('ElecTemp','global',dtype_r8,(/pcols,pver/),idx)

    ! Ion temperature (global so can write to history files)
    call pbuf_add_field('IonTemp', 'global',dtype_r8,(/pcols,pver/),idx)

  end subroutine ionos_register
!
!==============================================================================

  subroutine ionos_intr(state, ptend, pbuf, ztodt)

!-----------------------------------------------------------------------
! Interface for improved ionosphere simulation 
!-----------------------------------------------------------------------
!
!------------------------------Arguments--------------------------------

    use time_manager, only : get_nstep,get_step_size
     
    implicit none

    type(physics_state), intent(in)    :: state               ! physics state structure
    type(physics_ptend), intent(inout) :: ptend               ! parameterization tendency structure
    type(physics_buffer_desc),pointer  :: pbuf(:)             ! physics buffer

    type(ionos_state)                  :: istate              ! ionosphere state structure

    real(r8),            intent(in)    :: ztodt               ! Two times model timestep (2 delta-t)

!---------------------------Local storage-------------------------------
    integer :: lchnk                              ! chunk identifier
    integer :: ncol                               ! number of atmospheric columns
                
    integer :: ionBot			          ! bottom of ionosphere calculations

!----------------------------------------------------------------
!  Get current chunk number and number of columns in this chunk
!----------------------------------------------------------------
    lchnk = state%lchnk
    ncol  = state%ncol
    
!------------------------------------------------------------
!  Initialize data needed in the ionosphere calculations
!------------------------------------------------------------
    call ionos_datinit(state, lchnk, ncol, pbuf, istate, ionBot)

!-------------------------------
!  Get electron temperature
!-------------------------------
    call ionos_tempei(state, ptend, lchnk, ncol, ztodt, pbuf, istate, ionBot)

!---------------------------------------
!  Update heating of neutral atmosphere
!--------------------------------------

    return

  end subroutine ionos_intr

!===============================================================================

  subroutine ionos_datinit(state, lchnk, ncol, pbuf, istate, ionBot)
  
!------------------------------------------------------------------------------
! Time independent initialization for ionosphere simulation called in phys_init
! of physpkg module which is called in cam_comp module
!------------------------------------------------------------------------------
    use cam_history,      only : addfld, add_default ! Routines/variables needed for outputting fields
    use dycore,           only : get_resolution
    use interpolate_data, only : lininterp
    use constituents,     only : cnst_get_ind, cnst_mw            ! Routines to get molecular weights for constituents
    use hycoef,           only : hypm                             ! Model vertical pressure grid in Pascals
    use mo_apex,          only : bnorth, beast, bdown             ! Magnetic field components
    use time_manager,     only : get_curr_calday                  ! Routine to get current calendar day
    use mo_solar_parms,   only : get_solar_parms                  ! Routine to get solar parameters, i.e. f107
    use cam_control_mod,  only : nsrest        ! Variable to determine if this is an initial run or a restart/branch
    use time_manager,     only : get_nstep                        ! Routine to get current time step
    use physconst,        only : rairv, mbarv                     ! Constituent dependent rair and mbar
    use orbit,            only: zenith
     
    implicit none

    type(physics_buffer_desc), pointer  :: pbuf(:)             ! physics buffer
    type(physics_state), intent(in)     :: state               ! physics state structure
    type(ionos_state),   intent(inout)  :: istate              ! ionosphere state structure

    integer, intent(in)  :: lchnk   ! Chunk number 
    integer, intent(in)  :: ncol    ! Number of columns in current chunk 

    integer, intent(out) :: ionBot  ! bottom of ionosphere calculations 

!---------------------------Local storage-------------------------------

    integer :: indxUe     ! pbuf index for eastward ExB drift
    integer :: indxVe     ! pbuf index for northward ExB drift
    integer :: indxWe     ! pbuf index for upward ExB drift

    integer :: indxIR     ! pbuf index for ionization rates
    integer :: indxAIPRS  ! pbuf index for aurora ion production rate sum
    integer :: indxTe     ! pbuf index for electron temperature
    integer :: indxTi     ! pbuf index for ion temperature

    integer :: indxN2     ! cnst index for N2
    integer :: indxO1     ! cnst index for O
    integer :: indxO2     ! cnst index for O2
    integer :: indxNO     ! cnst index for NO
    integer :: indxH      ! cnst index for H
    integer :: indxN1     ! cnst index for N
    integer :: indxE      ! cnst index for electrons
    integer :: indxOp     ! cnst index for Op
    integer :: indxO2p    ! cnst index for O2p
    integer :: indxNOp    ! cnst index for NOp

    integer :: iVer       ! Counter for vertical loops
    integer :: iCol       ! Counter for column loops
    integer :: iIonR      ! Counter for ionization rates loops
   
    integer :: indxSP     ! pbuf index for Pedersen Conductivity
    integer :: indxSH     ! pbuf index for Hall Conductivity

    integer :: nstep      ! current timestep number

    real(r8), dimension(:), pointer     :: ue  ! pointer to eastward ExB drift in pbuf (from module iondrag)
    real(r8), dimension(:), pointer     :: ve  ! pointer to northward ExB drift in pbuf
    real(r8), dimension(:), pointer     :: we  ! pointer to upward ExB drift in pbuf


    real(r8), dimension(pcols,pver) :: tE                ! Pointer to electron temperature in pbuf (K) 
    real(r8), dimension(pcols,pver) :: ti                ! Pointer to ion temperature in pbuf (K) 

    real(r8), dimension(:,:,:), pointer :: ionRates     ! Pointer to ionization rates for O+,O2+,N+,N2+,NO+ in pbuf (s-1 from modules mo_jeuv and mo_jshort)


    real(r8), dimension(:,:), pointer   :: sigma_ped    ! Pointer to Pedersen Conductivity in pbuf (siemens/m) from module iondrag
    real(r8), dimension(:,:), pointer   :: sigma_hall   ! Pointer to Hall Conductivity in pbuf (siemens/m)


    real(r8), dimension(ncol)        :: geoLatR  ! Latitude (radians)  Make ncol because zenith aurora are ncol
    real(r8), dimension(ncol)        :: geoLonR  ! Longitude (radians)

    real(r8), dimension(pcols,pver)  :: pMid     ! Midpoint pressure (Pa)
    real(r8), dimension(pcols,pver)  :: tN       ! Neutral temperature (K)

    real(r8), dimension(pcols,pverp) :: pInt    ! Interface pressure (Pa)
    real(r8), dimension(pcols,pverp) :: tNInt   ! Interface Temperture (K)

    real(r8), dimension(ncol)       :: cosZenAngR            ! cosine of zenith angle (radians)
    real(r8), dimension(ncol)       :: zenAngD               ! zenith angle (degrees)

    real(r8), dimension(pcols,pver)  :: bNorth3d  ! northward component of magnetic field units?
    real(r8), dimension(pcols,pver)  :: bEast3d   ! eastward component of magnetic field
    real(r8), dimension(pcols,pver)  :: bDown3d          ! downward component of magnetic field

    real(r8), dimension(pcols,pver)  :: ue2d  ! horizontal x(eastward) component of ExB drift with added vertical dimension
    real(r8), dimension(pcols,pver)  :: ve2d  ! horizontal y(northward) component of ExB drift with added vertical dimension
    real(r8), dimension(pcols,pver)  :: we2d  ! vertical z component of ExB drift with added vertical dimension - midpoints
 
    real(r8), dimension(pcols,pver)  :: omega         ! vertical velocity at midpoint levels (Pa/s)
    real(r8), dimension(pcols,pver)  :: sourceR       ! R term of source g4 calculation
    real(r8), dimension(pcols,pver)  :: sourceEff     ! Efficiency term of source g4 calculation

    real(r8), dimension(pcols,pverp) :: wei2d   ! vertical z component of ExB drift with added vertical dimension - interfaces

    real(r8), dimension(pcols,pverp) :: omegai  ! vertical velocity on interface levels (Pa/s)
    
    real(r8), dimension(pcols,pverp) :: rairvi        ! Constituent dependent gas constant
 
    real(r8), dimension(pcols,pver)  :: dipMag  ! dip angle for each column (radians)
    real(r8), dimension(pcols,pver)  :: dipMagD ! dip angle for each column (degrees)

    real(r8) :: rmassN2    ! N2 molecular weight kg/kmol
    real(r8) :: rmassO2    ! O2 molecular weight kg/kmol
    real(r8) :: rmassO1    ! O atomic weight kg/kmol
    real(r8) :: rmassNO    ! NO molecular weight kg/kmol
    real(r8) :: rmassH     ! H atomic weight kg/kmol
    real(r8) :: rmassN1    ! N atomic weight kg/kmol
    real(r8) :: rmassE     ! Electron mass kg/kmol
    real(r8) :: rmassOp    ! O plus ion molecular weight kg/kmol
    real(r8) :: rmassO2p   ! O2 plus ion molecular weight kg/kmol
    real(r8) :: rmassNOp   ! NO plus ion molecular weight kg/kmol
 
    real(r8), dimension(pcols,pver) :: mmrN2        ! N2 mass mixing ratio kg/kg
    real(r8), dimension(pcols,pver) :: mmrO2    ! O2 mass mixing ratio kg/kg
    real(r8), dimension(pcols,pver) :: mmrO1    ! O mass mixing ratio kg/kg
    real(r8), dimension(pcols,pver) :: mmrNO    ! NO mass mixing ratio kg/kg
    real(r8), dimension(pcols,pver) :: mmrN1    ! N mass mixing ratio kg/kg
    real(r8), dimension(pcols,pver) :: mmrE     ! E electron mass mixing ratio kg/kg
    real(r8), dimension(pcols,pver) :: mmrOp    ! O plus ion mass mixing ratio kg/kg
    real(r8), dimension(pcols,pver) :: mmrO2p   ! O2 plus ion mass mixing ratio kg/kg
    real(r8), dimension(pcols,pver) :: mmrNOp        ! NO plus ion mass mixing ratio kg/kg

    real(r8), dimension(pcols,pver)  :: ndensN2         ! N2 number density (cm-3)
    real(r8), dimension(pcols,pver)  :: ndensO2  ! O2 number density (cm-3)
    real(r8), dimension(pcols,pver)  :: ndensO1  ! O number density (cm-3)
    real(r8), dimension(pcols,pver)  :: ndensNO  ! NO number density (cm-3)
    real(r8), dimension(pcols,pver)  :: ndensN1  ! N number density  (cm-3)
    real(r8), dimension(pcols,pver)  :: ndensE   ! E electron number density (cm-3)
    real(r8), dimension(pcols,pver)  :: ndensOp  ! O plus number density (cm-3)
    real(r8), dimension(pcols,pver)  :: ndensO2p ! O2 plus ion number density (cm-3)
    real(r8), dimension(pcols,pver)  :: ndensNOp ! NO plus ion number density  (cm-3)

    real(r8), dimension(pcols,pver,nIonRates) :: ionPRates    ! ionization rates temporary array (s-1 cm-3)
    real(r8), dimension(pcols,pver)           :: sumIonPRates ! Sum of ionization rates for O+,O2+,N+,N2+,NO+ (s-2 cm-3)

    real(r8), dimension(:,:), pointer         :: aurIPRateSum ! Auroral ion production sum for O2+,O+,N2+ (s-1 cm-3 from module mo_aurora)

    real(r8), dimension(pcols,pver)           :: sourceg4     ! g4 source term for electron/ion temperature update

    logical, dimension(ncol)                  :: do_aurora    ! Logical to test if aurora ion production rates calculated in mo_aurora (set to zero)

!--------------------------------------------------------------------------------
     
    tNInt(:,:)        = 0._r8
    zenAngD(:)        = 0._r8
    bNorth3d(:,:)     = 0._r8
    bEast3d(:,:)      = 0._r8
    bDown3d(:,:)      = 0._r8
    ue2d(:,:)         = 0._r8
    ve2d(:,:)         = 0._r8
    we2d(:,:)         = 0._r8
    wei2d(:,:)        = 0._r8
    omegai(:,:)       = 0._r8
    rairvi(:,:)       = 0._r8
    dipMag(:,:)       = 0._r8
    dipMagD(:,:)      = 0._r8
    ionPRates(:,:,:)  = 0._r8
    sumIonPRates(:,:) = 0._r8

    ndensN2(:,:)      = 0._r8
    ndensO2(:,:)      = 0._r8
    ndensO1(:,:)      = 0._r8
    ndensE(:,:)       = 0._r8
    ndensOp(:,:)      = 0._r8
    ndensO2p(:,:)     = 0._r8
    ndensNOp(:,:)     = 0._r8

    sourceR(:,:)      = 0._r8
    sourceEff(:,:)    = 0._r8
    sourceg4(:,:)     = 0._r8

!------------------------------------------------------------------------------------------------------
!  Set the bottom of the ionosphere calculations at around 0.01 hectopascals(millibars).  hypm is in 
!  Pascals.  Since the model vertical range goes from top down this can be used as the counter of the
!  number of levels in the range from the top down to the ionBot level. Also, set the interface 
!  variable ionBotP one more than ionBot
!------------------------------------------------------------------------------------------------------
    do iVer = 1, pver
   
       if (hypm(iVer) < 50._r8) ionBot = iVer
   
    enddo    
    
!--------------------------------------------------------------
!  Set radians to degrees variable and get zenith angle
!--------------------------------------------------------------
    rads2Degs   = 180._r8/pi
    degs2Rads   = pi/180._r8

!----------------------------------------------------------------
!  Get latitude and longitude of each column in this chunk
!----------------------------------------------------------------    
    geoLatR(1:ncol) = state%lat(1:ncol)
    geoLonR(1:ncol) = state%lon(1:ncol)

!-------------------------------------------------------------------------------------------------------
!  Need to get midpoint and interface pressure and neutral temperature from state structure (pcols,pver)
!-------------------------------------------------------------------------------------------------------
    pMid(1:ncol,1:pver) = state%pmid(1:ncol,1:pver)
    tN(1:ncol,1:pver)   = state%t(1:ncol,1:pver)
    
!-------------------------------------------------------------------------------------
!  Calculate neutral temperature on interface levels.  tN vertical dimension is pver
!-------------------------------------------------------------------------------------   
    do iVer = 2, pver
 
      do iCol = 1, ncol

        tNInt(iCol,iVer) = 0.5_r8 * tN(iCol,iVer) + 0.5_r8 * tN(iCol,iVer-1)

      enddo
    enddo

    do iCol = 1, ncol
        tNInt(iCol,1) = 1.5_r8 * tNInt(iCol,2) - 0.5_r8 * tNInt(iCol,3) 
    enddo
    do iCol = 1, ncol
        tNInt(iCol,pverp) = 1.5_r8 * tNInt(iCol,pver) - 0.5_r8 * tNInt(iCol,pver-1) 
    enddo
    
    istate%tNInt(1:ncol,1:pverp) = tNInt(1:ncol,1:pverp)

!--------------------------------------------------------------
!  Get zenith angle
!-------------------------------------------------------------- 
    calDay = get_curr_calday()    
    call zenith(calDay,geoLatR,geoLonR,cosZenAngR,ncol)

    do iCol = 1, ncol

      zenAngD(iCol) = ACOS(cosZenAngR(iCol)) * rads2Degs
    
    enddo
 
    istate%cosZenAngR(1:ncol) = cosZenAngR(1:ncol)
    istate%zenAngD(1:ncol) = zenAngD(1:ncol)

!--------------------------------------------------------------
!  Get F10.7 solar flux
!--------------------------------------------------------------
    call get_solar_parms( f107_s = f107 )

!---------------------------------------------------------------------------------------
!  Expand magnetic field components in vertical to make 3D, pcols,pver,begchunk:endchunk
!---------------------------------------------------------------------------------------
    do iVer = 1, pver
    
      do iCol = 1, ncol

        bNorth3d(iCol,iVer) = bnorth(iCol,lchnk)
        bEast3d(iCol,iVer) = beast(iCol,lchnk)
        bDown3d(iCol,iVer) = bdown(iCol,lchnk)

      enddo
    
    enddo

    istate%bNorth3d(1:ncol,1:pver) = bNorth3d(1:ncol,1:pver)
    istate%bEast3d(1:ncol,1:pver)  = bEast3d(1:ncol,1:pver)
    istate%bDown3d(1:ncol,1:pver)  = bDown3d(1:ncol,1:pver)

!??? why do we need 3d mag field???  Needed later for ambipolar diffusion?

!---------------------------------------------------------------------------------
! Get ion ExB drift from physics buffer (they were defined by exbdrift module)
! Ion drifts are 2d output arrays, i.e., no vertical dimension but 1d here (pcols)
!---------------------------------------------------------------------------------
    indxUe = pbuf_get_index( 'UE' )
    indxVe = pbuf_get_index( 'VE' )
    indxWe = pbuf_get_index( 'VE' )
    call pbuf_get_field(pbuf, indxUe, ue)
    call pbuf_get_field(pbuf, indxVe, ve)
    call pbuf_get_field(pbuf, indxWe, we)

!-------------------------------------------------------------------------------
!  Form 2D drifts needed for this module.  Vertical level goes to ionBotP since
!  needed below to get vertical drift on interface levels
!-------------------------------------------------------------------------------
    do iVer = 1, pver
    
      do iCol = 1, ncol

        ue2d(iCol,iVer) = ue(iCol)
        ve2d(iCol,iVer) = ve(iCol)
        we2d(iCol,iVer) = we(iCol)

      enddo
    
    enddo
 
    istate%ue2d(1:ncol,1:pver) = ue2d(1:ncol,1:pver)
    istate%ve2d(1:ncol,1:pver) = ve2d(1:ncol,1:pver)
    istate%we2d(1:ncol,1:pver) = we2d(1:ncol,1:pver)
   
!-------------------------------------------------------------------------------
!  Need vertical ExB drift on interface levels
!-------------------------------------------------------------------------------
    do iVer = 2, pver
      do iCol = 1, ncol
        wei2d(iCol,iVer) = 0.5_r8 * (we2d(iCol,iVer-1) + we2d(iCol,iVer))
      enddo    
    enddo
    do iCol = 1, ncol
        wei2d(iCol,1) = 1.5_r8 * wei2d(iCol,2) - 0.5_r8 * wei2d(iCol,3)
    enddo
    do iCol = 1, ncol
        wei2d(iCol,pverp) = 1.5_r8 * wei2d(iCol,pver) - 0.5_r8 * wei2d(iCol,pver-1)
    enddo
    
    istate%wei2d(1:ncol,1:pverp) = wei2d(1:ncol,1:pverp)

!--------------------------------------------------------------------------------------
!  Need to get vertical velocity on interface levels handling top level specifically
!--------------------------------------------------------------------------------------
    omega(1:ncol,1:pver) = state%omega(1:ncol,1:pver)
    
    do iVer = 2, pver
      do iCol = 1, ncol
        omegai(iCol,iVer) = 0.5_r8 * (omega(iCol,iVer-1) + omega(iCol,iVer))
      enddo
    enddo
    do iCol = 1, ncol
        omegai(iCol,1) = 1.5_r8 * omegai(iCol,2) - 0.5_r8 * omegai(iCol,3)
    enddo
    do iCol = 1, ncol
        omegai(iCol,pverp) = 1.5_r8 * omegai(iCol,pver) - 0.5_r8 * omegai(iCol,pver-1)
    enddo

    istate%omegai(1:ncol,1:pverp) = omegai(1:ncol,1:pverp)

!------------------------------------------------------------------------
!  Get constituent dependent gas constant and derive on interface levels
!------------------------------------------------------------------------        
    do iVer = 2, pver
      do iCol = 1, ncol
        rairvi(iCol,iVer) = 0.5_r8 * rairv(iCol,iVer-1,lchnk) + 0.5_r8 * rairv(iCol,iVer,lchnk)
      enddo
    enddo

    do iCol = 1, ncol
       rairvi(iCol,1) = 1.5_r8 * rairvi(iCol,2) - 0.5_r8 * rairvi(iCol,3)
    enddo
    do iCol = 1, ncol
       rairvi(iCol,pverp) = 1.5_r8 * rairvi(iCol,pver) - 0.5_r8 * rairvi(iCol,pver-1)
    enddo

    istate%rairvi(1:ncol,1:pverp) = rairvi(1:ncol,1:pverp)

!-------------------------------------------------------------------------------
!  Need to get dip angle from magnetic field components
!-------------------------------------------------------------------------------     
    do iVer = 1, pver
      do iCol = 1, ncol
        dipMag(iCol,iVer) = ATAN(bDown3d(iCol,iVer) / SQRT(bNorth3d(iCol,iVer)**2 + bEast3d(iCol,iVer)**2))
        if (dipMag(iCol,iVer) < 0.17_r8 .and. dipMag(iCol,iVer) > 0._r8 ) dipMag(iCol,iVer) = 0.17_r8
        if (dipMag(iCol,iVer) > -0.17_r8 .and. dipMag(iCol,iVer) < 0._r8 ) dipMag(iCol,iVer) = 0.17_r8
        dipMagD(iCol,iVer) = dipMag(iCol,iVer) * rads2Degs
      enddo
    enddo
 
    istate%dipMag(1:ncol,1:pver)  = dipMag(1:ncol,1:pver)
    istate%dipMagD(1:ncol,1:pver) = dipMagD(1:ncol,1:pver)

!----------------------------------------------------------------------------------------
!  Get molecular weights using constituent indices (kg/kmole) except N2 which is hard-wired
!----------------------------------------------------------------------------------------
    call cnst_get_ind( 'O',   indxO1 )
    call cnst_get_ind( 'O2',  indxO2 )
    call cnst_get_ind( 'NO',  indxNO )
    call cnst_get_ind( 'H',   indxH )
    call cnst_get_ind( 'N',   indxN1 )
    call cnst_get_ind( 'e',   indxE )
    call cnst_get_ind( 'Op',  indxOp )
    call cnst_get_ind( 'O2p', indxO2p )
    call cnst_get_ind( 'NOp', indxNOp )

    rmassO2 = cnst_mw(indxO2)
    rmassO1 = cnst_mw(indxO1)
    rmassNO = cnst_mw(indxNO)
    rmassH  = cnst_mw(indxH) 
    rmassN1  = cnst_mw(indxN1) 
    rmassE  = cnst_mw(indxE)        !This is hard-wired in 'iondrag.F90' module 
    rmassOP = cnst_mw(indxOp)
    rmassO2P = cnst_mw(indxO2p)
    rmassNOP = cnst_mw(indxNOp)

    rmassN2 = 28._r8

!-----------------------------------------------------------------------------------------------------
!  Get mass mixing ratios from state structure q array (kg/kg) except N2 which is calculated from O,O2
!-----------------------------------------------------------------------------------------------------

    mmrO2(1:ncol,1:pver)  = state%q(1:ncol,1:pver,indxO2)
    mmrO1(1:ncol,1:pver)  = state%q(1:ncol,1:pver,indxO1)
    mmrNO(1:ncol,1:pver)  = state%q(1:ncol,1:pver,indxNO)
    mmrN1(1:ncol,1:pver)  = state%q(1:ncol,1:pver,indxN1)
    mmrE(1:ncol,1:pver)   = state%q(1:ncol,1:pver,indxE)                  
    mmrOp(1:ncol,1:pver)  = state%q(1:ncol,1:pver,indxOp)
    mmrO2p(1:ncol,1:pver) = state%q(1:ncol,1:pver,indxO2p)
    mmrNOp(1:ncol,1:pver) = state%q(1:ncol,1:pver,indxNOp)    

    mmrN2(1:ncol,1:pver) = 1._r8 - (mmrO2(1:ncol,1:pver) + mmrO1(1:ncol,1:pver)) 
    mmrN2(1:ncol,1:pver) = MAX(1.e-20_r8,mmrN2(1:ncol,1:pver))

!-------------------------------------------------------------------------------------------------
!  Need to get number density from mass mixing ratio.  state%mbarv is g/mole, same as rmass units
!  kg/kg * (g/mole)/(g/mole) * (Pa or N/m*m)/(Joules or N*m/molecule*K) * (K) = m-3 * 1E-06 = cm-3
!--------------------------------------------------------------------------------------------------
    do iVer = 1, pver
      do iCol = 1, ncol
            ndensO2(iCol,iVer)  = mmrO2(iCol,iVer) * mbarv(iCol,iVer,lchnk) / rmassO2 * &
                                           pMid(iCol,iVer) / (kboltz * tN(iCol,iVer)) * 1.E-06_r8
            ndensO1(iCol,iVer)  = mmrO1(iCol,iVer) * mbarv(iCol,iVer,lchnk) / rmassO1 * &
                                           pMid(iCol,iVer) / (kboltz * tN(iCol,iVer)) * 1.E-06_r8
            ndensNO(iCol,iVer)  = mmrNO(iCol,iVer) * mbarv(iCol,iVer,lchnk) / rmassNO * &
                                           pMid(iCol,iVer) / (kboltz * tN(iCol,iVer)) * 1.E-06_r8
            ndensN1(iCol,iVer)  = mmrN1(iCol,iVer) * mbarv(iCol,iVer,lchnk) / rmassN1 * &
                                        pMid(iCol,iVer) / (kboltz * tN(iCol,iVer)) * 1.E-06_r8
            ndensE(iCol,iVer)   = mmrE(iCol,iVer) * mbarv(iCol,iVer,lchnk) / rmassE * &
                                        pMid(iCol,iVer) / (kboltz * tN(iCol,iVer)) * 1.E-06_r8
            ndensOp(iCol,iVer)  = mmrOp(iCol,iVer) * mbarv(iCol,iVer,lchnk) / rmassOp * &
                                        pMid(iCol,iVer) / (kboltz * tN(iCol,iVer)) * 1.E-06_r8
            ndensO2p(iCol,iVer) = mmrO2p(iCol,iVer) * mbarv(iCol,iVer,lchnk) / rmassO2p * &
                                      pMid(iCol,iVer) / (kboltz * tN(iCol,iVer)) * 1.E-06_r8
            ndensNOp(iCol,iVer) = mmrNOp(iCol,iVer) * mbarv(iCol,iVer,lchnk) / rmassNOp * &
                                        pMid(iCol,iVer) / (kboltz * tN(iCol,iVer)) * 1.E-06_r8
            ndensN2(iCol,iVer)  = mmrN2(iCol,iVer) * mbarv(iCol,iVer,lchnk) / rmassN2 * &
                                    pMid(iCol,iVer) / (kboltz * tN(iCol,iVer)) * 1.E-06_r8             
      enddo
    enddo
 
    istate%ndensO2(1:ncol,1:pver)  = ndensO2(1:ncol,1:pver)
    istate%ndensO1(1:ncol,1:pver)  = ndensO1(1:ncol,1:pver)
    istate%ndensNO(1:ncol,1:pver)  = ndensNO(1:ncol,1:pver)
    istate%ndensN1(1:ncol,1:pver)  = ndensN1(1:ncol,1:pver)
    istate%ndensE(1:ncol,1:pver)   = ndensE(1:ncol,1:pver)
    istate%ndensOp(1:ncol,1:pver)  = ndensOp(1:ncol,1:pver)
    istate%ndensO2p(1:ncol,1:pver) = ndensO2p(1:ncol,1:pver)
    istate%ndensNOp(1:ncol,1:pver) = ndensNOp(1:ncol,1:pver)
    istate%ndensN2(1:ncol,1:pver)  = ndensN2(1:ncol,1:pver)

!-------------------------------------------------------------------------------
!  Write output fields to history file
!-------------------------------------------------------------------------------
    call outfld ('IONOS_NDENSN2' , ndensN2  , pcols, lchnk)
    call outfld ('IONOS_NDENSO1' , ndensO1  , pcols, lchnk)
    call outfld ('IONOS_NDENSO2' , ndensO2  , pcols, lchnk)
    call outfld ('IONOS_NDENSNO' , ndensNO  , pcols, lchnk)
    call outfld ('IONOS_NDENSN1' , ndensN1  , pcols, lchnk)
    call outfld ('IONOS_NDENSE'  , ndensE   , pcols, lchnk)
    call outfld ('IONOS_NDENSOP' , ndensOp  , pcols, lchnk)
    call outfld ('IONOS_NDENSO2P', ndensO2p , pcols, lchnk)
    call outfld ('IONOS_NDENSNOP', ndensNOp , pcols, lchnk)

!------------------------------------------------------------------------------------
! Get ionization rates from physics buffer which were calculated in mo_jeuv and 
! mo_jshort modules.  Rates array dimensions are pcols, pver, nIonRates.  Units s-1
!------------------------------------------------------------------------------------
    indxIR = pbuf_get_index( 'IonRates' )
    call pbuf_get_field(pbuf, indxIR, ionRates, start=(/1,1,1/), kount=(/pcols,pver,nIonRates/))

!----------------------------------------------------------------------------------------------
!  Need to convert these ionization rates to ion production rates by multiplying number density 
!  of neutral species appropriate from reactions in mo_jeuv(jeuv) and mo_jshort(jshort)(for NO)  
!----------------------------------------------------------------------------------------------         
    do iVer = 1, pver
      do iCol = 1, ncol
        do iIonR = 1, nIonRates
          IF (iIonR <= 3) ionPRates(iCol,iVer,iIonR) = ionRates(iCol,iVer,iIonR) * ndensO1(iCol,iVer)
          IF (iIonR == 4) ionPRates(iCol,iVer,iIonR) = ionRates(iCol,iVer,iIonR) * ndensN1(iCol,iVer)
          IF ((iIonR == 5) .OR. (iIonR >= 7 .AND. iIonR <= 9)) &
                                    ionPRates(iCol,iVer,iIonR) = ionRates(iCol,iVer,iIonR) * ndensO2(iCol,iVer)
          IF (iIonR == 6 .OR. iIonR == 10 .OR. iIonR == 11) &
                                    ionPRates(iCol,iVer,iIonR) = ionRates(iCol,iVer,iIonR) * ndensN2(iCol,iVer)
        enddo
                                    
!----------------------------------------------
!  Sum ion production rates all reactions
!----------------------------------------------   
        sumIonPRates(iCol,iVer) = SUM(ionPRates(iCol,iVer,1:11))
      enddo
    enddo

    istate%ionPRates(1:ncol,1:pver,1:nIonRates) = ionPRates(1:ncol,1:pver,1:nIonRates)
    istate%sumIonPRates(1:ncol,1:pver) = sumIonPRates(1:ncol,1:pver)
        
!-------------------------------------------------------------------------------------------
! Get aurora ion production rate sum from physics buffer which were calculated in mo_aurora 
! module.  Rate array dimensions are pcols, pver.  Units s-1 cm-3
!-------------------------------------------------------------------------------------------
    indxAIPRS = pbuf_get_index( 'AurIPRateSum' )
    call pbuf_get_field(pbuf, indxAIPRS, aurIPRateSum)

!-------------------------------------------------------------------------------------
! Check latitudes, and set aurora ion production rates for all columns and levels 
! to zero if all below 32.5 deg since not calculated in mo_aurora for those latitudes 
! (same as criteria in mo_aurora)
!-------------------------------------------------------------------------------------
    do_aurora(:) = abs( geoLatR(:) ) > pi/6._r8
    if( all( .not. do_aurora(:) ) ) then
       aurIPRateSum(:,:) = 0._r8
    end if
 
!-------------------------------------------------------------------------------------------------
!  Calculate electron heating rate which is source in electron/ion temperature derivation
!-------------------------------------------------------------------------------------------------
    do iVer = 1, ionBot
      do iCol = 1, ncol
        sourceR(iCol,iVer) = LOG( ndensE(iCol,iVer) / (ndensO2(iCol,iVer) + ndensN2(iCol,iVer) + &
                                                                                  0.1_r8 * ndensO1(iCol,iVer)) )
        sourceEff(iCol,iVer) = EXP( -(12.75_r8 + 6.941_r8 * sourceR(iCol,iVer) + 1.166_r8 * sourceR(iCol,iVer)**2 + &
                                        0.08043_r8 * sourceR(iCol,iVer)**3 + 0.001996_r8 * sourceR(iCol,iVer)**4) )

!-------------------------------------------------------------------------------
!  Calculate g4 source term for electron temperature update
!-------------------------------------------------------------------------------         
        sourceg4(iCol,iVer) = (sumIonPRates(iCol,iVer) + aurIPRateSum(iCol,iVer)) * sourceEff(iCol,iVer)

      enddo

    enddo

    istate%sourceg4(1:ncol,1:pver) = sourceg4(1:ncol,1:pver)
      call outfld ('IONI_RATE'      , sumIonPRates   , pcols, lchnk)
      call outfld ('IONI_Eff'      , sourceEff   , pcols, lchnk)
      call outfld ('SOURCE_g4'      , sourceg4       , pcols, lchnk)
      call outfld ('AUR_IRATESUM'      , aurIPRateSum       , pcols, lchnk)

!----------------------------------------------------------------------------------------------
! Get Pedersen and Hall Conductivities from physics buffer which were calculated in iondrag 
! module.  Conductivity array dimensions are pcols, pver
!-------------------------------------------------------------------------------
    indxSP = pbuf_get_index( 'PedConduct' )
    indxSH = pbuf_get_index( 'HallConduct' )
    call pbuf_get_field(pbuf, indxSP, sigma_ped)
    call pbuf_get_field(pbuf, indxSH, sigma_hall)

    return

  end subroutine ionos_datinit
!
!===============================================================================

  subroutine ionos_tempei(state, ptend, lchnk, ncol, ztodt, pbuf, istate, ionBot)

!-----------------------------------------------------------------------
! Routine to compute the electron and ion temperature needed for
! ambipolar diffusion calculations.
!-----------------------------------------------------------------------

    use cam_history,     only : addfld, add_default
    use phys_grid,       only : get_lat_p, get_lon_p, get_rlat_p,get_rlon_p
    use physconst,       only : gravit ! Gravity (m/s2)
    use cam_control_mod, only : nsrest        ! Variable to determine if this is an initial run or a restart/branch
    use time_manager,    only : get_nstep                        ! Routine to get current time step
    use physconst,       only : rairv, mbarv                     ! Constituent dependent rair and mbar
     
    implicit none

!------------------------------Arguments--------------------------------

    type(physics_buffer_desc), pointer   :: pbuf(:)             ! physics buffer
    type(physics_state),   intent(in)    :: state               ! physics state structure
    type(physics_ptend),   intent(inout) :: ptend               ! parameterization tendency structure
    type(ionos_state),     intent(inout) :: istate              ! ionosphere state structure

    real(r8), intent(in) :: ztodt                            ! Two times model timestep (2 delta-t)


    integer, intent(in) :: lchnk   ! chunk identifier
    integer, intent(in) :: ncol    ! number of atmospheric columns
    integer, intent(in) :: ionBot  ! bottom of ionosphere calculations 

!---------------------------Local storage-------------------------------
    integer :: ionBotP                                  ! bottom of ionosphere calculations plus one more level 

    integer :: i,k                                      ! loop indexes
    integer :: iVer                                     ! Counter for vertical loops
    integer :: iCol                                     ! Counter for column loops
    integer :: iter                                     ! Counter for iteration loop

    integer :: nstep                                    ! current timestep number

    integer :: indxTe                                   ! pbuf index for electron temperature
    integer :: indxTi                                   ! pbuf index for ion temperature

    integer, parameter  :: maxIter  = 6                 ! maximum number of iterations to solve for electron/ion temperature

    real(r8), parameter :: Kec1   = 7.5E5_r8            ! c1 constant for calculation of electron conductivity(Ke)
    real(r8), parameter :: Kec2   = 3.22E4_r8           ! c2 constant for calculation of electron conductivity(Ke)
    real(r8), parameter :: stepweight  = 1.0_r8         ! weight of previous and current times step for diagonals

    real(r8) :: wrk1                                    ! 2/3/kboltz_ev
    real(r8) :: sqrt2

    real(r8) :: lossc5                                  ! c5 prime of Lc(eN2) component of loss term
    real(r8) :: lossc7                                  ! c7 of Lc(eO2) component of loss term equation

    real(r8) :: lossc9                                  ! c9 of Lc(eO) component of loss term equation
    real(r8) :: lossc13                                 ! c13 of Lc(eO)f component of loss term equation
    real(r8) :: FeDB                                    ! B term of electron heat flux of UB
    real(r8) :: FeD                                     ! Day time flux
    real(r8) :: FeN                                     ! Night time flux
    real(r8) :: f1Ted1                                  ! d1 of f1(Te) calculation used to get electron conductivity
    real(r8) :: f1Ted2                                  ! d2 of f1(Te) calculation used to get electron conductivity
    real(r8) :: f1Ted3                                  ! d3 of f1(Te) calculation used to get electron conductivity
    real(r8) :: f1Te

    real(r8), dimension(:,:), pointer   :: tE           ! Pointer to electron temperature in pbuf (K) 
    real(r8), dimension(:,:), pointer   :: ti           ! Pointer to ion temperature in pbuf (K) 

    real(r8), dimension(pcols,pver)  	:: pMid 	! Midpoint pressure (Pa)
    real(r8), dimension(pcols,pver)  	:: tN		! Neutral temperature (K)
    real(r8), dimension(pcols,pver)  	:: tE0  	! Electron temperature from last time step (K)
    real(r8), dimension(pcols,pver)  	:: tEPrevI	! Electron temperature from previous iteration (K)

    real(r8), dimension(pcols,pverp) 	:: pInt 	! Interface pressure (Pa)
    real(r8), dimension(pcols,pverp) 	:: tNInt	! Interface Temperture (K)
    real(r8), dimension(pcols,pverp) 	:: rairvi	! Constituent dependent gas constant on interface levels

    real(r8), dimension(pcols,pver)  	:: ndensN2	! N2 number density (cm-3)
    real(r8), dimension(pcols,pver)  	:: ndensO2	! O2 number density (cm-3)
    real(r8), dimension(pcols,pver)  	:: ndensO1	! O number density (cm-3)
    real(r8), dimension(pcols,pver)  	:: ndensNO	! NO number density (cm-3)
    real(r8), dimension(pcols,pver)  	:: ndensN1	! N number density  (cm-3)
    real(r8), dimension(pcols,pver)  	:: ndensE	! E electron number density (cm-3)
    real(r8), dimension(pcols,pver)  	:: ndensOp	! O plus number density (cm-3)
    real(r8), dimension(pcols,pver)  	:: ndensO2p	! O2 plus ion number density (cm-3)
    real(r8), dimension(pcols,pver)  	:: ndensNOp	! NO plus ion number density  (cm-3)
 
    real(r8), dimension(pcols,pver)  	:: dipMag	! dip angle for each column (radians)
    real(r8), dimension(pcols,pver)  	:: dipMagD	! dip angle for each column (degrees)

    real(r8), dimension(ncol)       	:: zenAngD	! zenith angle (degrees)

    real(r8), dimension(pcols)       	:: Feub 	! electron heat flux at upper boundary
 
    real(r8), dimension(pcols,pver)  	:: sqrtTE	! Square root of electron temperature

    real(r8), dimension(pcols,pver)  	:: sqrtTN	! Square root of electron temperature
 
    real(r8), dimension(pver)        	:: Ke		! electron conductivity

    real(r8), dimension(pverp)       	:: Kei  	! electron conductivity interface levels

    real(r8), dimension(pcols,pver)  	:: lossc4p	! c4 prime of Lc(eN2) component of loss term
    real(r8), dimension(pcols,pver)  	:: lossceN2	! Lc(eN2) component of loss term equation
 
    real(r8), dimension(pcols,pver)  	:: lossc6p	! c6 prime of Lc(eO2) component of loss term equation
    real(r8), dimension(pcols,pver)  	:: lossceO2	! Lc(eO2) component of loss term equation
 
    real(r8), dimension(pcols,pver)  	:: lossc8p	! c8 prime of Lc(eO) component of loss term equation
    real(r8), dimension(pcols,pver)  	:: lossceO1	! Lc(eO) component of loss term equation

    real(r8), dimension(pcols,pver)  	:: lossc10p	! c10 prime of Lc(eN2) component of loss term equation
    real(r8), dimension(pcols,pver)  	:: losscA	! A of Lc(eN2)v component of loss term equation
    real(r8), dimension(pcols,pver)  	:: tENDiff	! Difference between electron and neutral temperatures
    real(r8), dimension(pcols,pver)  	:: lossceN2v	! Lc(eN2)v component of loss term equation

    real(r8), dimension(pcols,pver)  	:: lossc11p	! c11 prime of Lc(eO2)v component of loss term equation
    real(r8), dimension(pcols,pver)  	:: lossceO2v	! Lc(eO2)v component of loss term equation
 
    real(r8), dimension(pcols,pver)  	:: lossc12p	! c12 prime of Lc(eO)f component of loss term equation
    real(r8), dimension(pcols,pver)  	:: lossceOf	! Lc(eO)f component of loss term equation
 
    real(r8), dimension(pcols,pver)  	:: lossc14p	! c14 prime of Lc(eO)1D component of loss term equation
    real(r8), dimension(pcols,pver)  	:: losscf2d	! d of f2 of Lc(eO)1D component of loss term equation
    real(r8), dimension(pcols,pver)  	:: losscf2	! f2 of Lc(eO)1D component of loss term equation
    real(r8), dimension(pcols,pver)  	:: losscf3	! f3 of Lc(eO)1D component of loss term equation
    real(r8), dimension(pcols,pver)  	:: lossceO1D	! Lc(eO)1D component of loss term equation
 
    real(r8), dimension(pcols,pver)  	:: lossc15p	! c15 prime of Lc(eN2)Rot component of loss term equation
    real(r8), dimension(pcols,pver)  	:: lossceN2Rot  ! Lc(eN2)Rot component of loss term equation
 
    real(r8), dimension(pcols,pver)  	:: lossc16p	! c16 prime of Lc(eO2)Rot component of loss term equation
    real(r8), dimension(pcols,pver)  	:: lossceO2Rot  ! Lc(eO2)Rot component of loss term equation
 
    real(r8), dimension(pcols,pver)  	:: lossc3p	! c3 prime of Lc(ei) component of loss term equation
    real(r8), dimension(pcols,pver)  	:: losscei	! Lc(ei) component of loss term equation
    real(r8), dimension(pcols,pver)  	:: losscin	! ion-neutral heating coeff.

    real(r8), dimension(pcols,pver)  	:: lossg3	! g3 loss term for Te tendency

    real(r8), dimension(pcols,pver)  	:: delTEN	! Difference between electron and neutral temperatures from production/loss

    real(r8), dimension(pcols,pverp) 	:: delZi	! Delta z: interfaces
    real(r8), dimension(pcols,pver)  	:: delZ 	! Delta z: midpoints

    real(r8), dimension(pcols,pver)  	:: sourceg4	! g4 source term for electron/ion temperature update

 
    real(r8), dimension(pcols,pver)  	:: qjoule	! joule heating
    real(r8), dimension(pcols,pver)  	:: qen  	! electron-neutral heating
    real(r8), dimension(pcols,pver)  	:: qei  	! electron-ion Coulomb heating
    real(r8), dimension(pcols,pver)  	:: rho  	! mass density

    real(r8), dimension(pcols,pver)  	:: wrk2

    real(r8), dimension(ionBot)      	:: subdiag	! subdiagonal values for Te tendency solving
    real(r8), dimension(ionBot)      	:: superdiag	! superdiagonal values for Te tendency solving
    real(r8), dimension(ionBot)      	:: diag 	! diagonal values for Te tendency solving
    real(r8), dimension(ionBot)      	:: rHSInit	! initial RHS of electron temperature update
    real(r8), dimension(ionBot)      	:: rHSH 	! h for RHS of electron temperature update
    real(r8), dimension(ionBot)      	:: rHS  	! RHS of electron temperature update
    real(r8), dimension(ionBot)      	:: tETemp	! temporary electron temperature array for input to tridag

    logical, dimension(pcols)        	:: colConv	! flag for column converging = 1 if converged otherwise = 0

    logical                          	:: converged	!Flag for convergence in electron temperature calculation iteration loop
         
!-----------------------------------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------------------------- 
!  Initialize arrays to zero
!---------------------------------------------------------------------------------------------------------
    pMid(:,:)           = 0._r8
    tN(:,:)             = 0._r8
    pInt(:,:)           = 0._r8
    tNInt(:,:)          = 0._r8
    rairvi(:,:)         = 0._r8
    ndensO2(:,:)        = 0._r8
    ndensO1(:,:)        = 0._r8
    ndensNO(:,:)        = 0._r8
    ndensN1(:,:)        = 0._r8
    ndensE(:,:)         = 0._r8
    ndensOp(:,:)        = 0._r8
    ndensO2p(:,:)       = 0._r8
    ndensNOp(:,:)       = 0._r8
    ndensN2(:,:)        = 0._r8
    zenAngD(:)          = 0._r8
    sqrtTE(:,:)         = 0._r8
    sqrtTN(:,:)         = 0._r8
    Ke(:)               = 0._r8
    Kei(:)              = 0._r8
    lossc4p(:,:)        = 0._r8
    lossceN2(:,:)       = 0._r8
    lossc6p(:,:)        = 0._r8
    lossc6p(:,:)        = 0._r8
    lossceO2(:,:)       = 0._r8
    lossc8p(:,:)        = 0._r8
    lossceO1(:,:)       = 0._r8
    lossc10p(:,:)       = 0._r8
    losscA(:,:)         = 0._r8
    tENDiff(:,:)        = 0._r8
    lossceN2v(:,:)      = 0._r8
    lossc11p(:,:)       = 0._r8
    lossceO2v(:,:)      = 0._r8
    lossc12p(:,:)       = 0._r8
    lossceOf(:,:)       = 0._r8
    lossc14p(:,:)       = 0._r8
    losscf2d(:,:)       = 0._r8
    losscf2(:,:)        = 0._r8
    losscf3(:,:)        = 0._r8
    lossceO1D(:,:)      = 0._r8
    lossc15p(:,:)       = 0._r8
    lossceN2Rot(:,:)    = 0._r8
    lossc16p(:,:)       = 0._r8
    lossceO2Rot(:,:)    = 0._r8
    lossc3p(:,:)        = 0._r8
    losscei(:,:)        = 0._r8
    losscin(:,:)        = 0._r8
    lossg3(:,:)         = 0._r8
    delTEN(:,:)         = 0._r8 
    delZi(:,:)          = 0._r8
    delZ(:,:)           = 0._r8 
    subDiag(:)          = 0._r8         
    superDiag(:)        = 0._r8        
    diag(:)             = 0._r8 
    sourceg4(:,:)       = 0._r8 
    dipMag(:,:)         = 0._r8
    dipMagD(:,:)        = 0._r8
    rHSInit(:)          = 0._r8 
    rHSH(:)             = 0._r8
    rHS(:)              = 0._r8 
    teTemp(:)           = 0._r8 
    qjoule(:,:)         = 0._r8
    qei(:,:)            = 0._r8
    qen(:,:)            = 0._r8
    rho(:,:)            = 0._r8
    colConv(:)          = .false.
    
!-------------------------------------------
!  Calculate some commonly used variables
!-------------------------------------------
    wrk1 = 2._r8 / 3._r8/ kboltz_ev
    sqrt2 = SQRT(2._r8)
    ionBotP = ionBot + 1

!-------------------------------------------------------------------------------------------------------
!  Need to get midpoint and interface pressure and neutral temperature from state structure (ncol,ionBot)
!-------------------------------------------------------------------------------------------------------
    pMid(1:ncol,1:pver)  = state%pmid(1:ncol,1:pver)
    tN(1:ncol,1:pver)    = state%t(1:ncol,1:pver)
    rho(1:ncol,1:pver) = pMid(1:ncol,1:pver)/rairv(1:ncol,1:pver,lchnk)/tN(1:ncol,1:pver) * 1.E-3_r8     ! convert to g/cm3

    qjoule(1:ncol,1:ionBot) = ptend%s(1:ncol,1:ionBot) * 6.24E15_r8     ! convert from J/kg/s to ev/g/s

    pInt(1:ncol,1:pverp)  = state%pint(1:ncol,1:pverp)
    tNInt(1:ncol,1:pverp)   = istate%tNInt(1:ncol,1:pverp)        
    rairvi(1:ncol,1:pverp)  = istate%rairvi(1:ncol,1:pverp)

    sqrtTN(1:ncol,1:pver) = SQRT(tN(1:ncol,1:pver))

!----------------------------------------------------------------
!  Get variables needed from the ionosphere state structure
!----------------------------------------------------------------
    ndensO2(1:ncol,1:pver)  = istate%ndensO2(1:ncol,1:pver) 
    ndensO1(1:ncol,1:pver)  = istate%ndensO1(1:ncol,1:pver)
    ndensNO(1:ncol,1:pver)  = istate%ndensNO(1:ncol,1:pver) 
    ndensN1(1:ncol,1:pver)  = istate%ndensN1(1:ncol,1:pver) 
    ndensE(1:ncol,1:pver)   = istate%ndensE(1:ncol,1:pver)  
    ndensOp(1:ncol,1:pver)  = istate%ndensOp(1:ncol,1:pver) 
    ndensO2p(1:ncol,1:pver) = istate%ndensO2p(1:ncol,1:pver)
    ndensNOp(1:ncol,1:pver) = istate%ndensNOp(1:ncol,1:pver)
    ndensN2(1:ncol,1:pver)  = istate%ndensN2(1:ncol,1:pver) 
 
    sourceg4(1:ncol,1:pver)  = istate%sourceg4(1:ncol,1:pver)

    dipMag(1:ncol,1:pver)   = istate%dipMag(1:ncol,1:pver)
    dipMagD(1:ncol,1:pver)  = istate%dipMagD(1:ncol,1:pver)

    zenAngD(1:ncol)  = istate%zenAngD(1:ncol) 
      
!-------------------------------------------------------------------------------------------------------------------
!  Get electron temperature from physics buffer and if this is the first time calculated then initialize to neutral
!  temperature.  nsrest=0 means this is an initial run and nstep=0 means this is the first time step. 
!-------------------------------------------------------------------------------------------------------------------
    indxTe = pbuf_get_index( 'ElecTemp' )
    call pbuf_get_field(pbuf, indxTe, tE)

    indxTi = pbuf_get_index( 'IonTemp' )
    call pbuf_get_field(pbuf, indxTi, ti)

    nstep = get_nstep()

    if(nsrest == 0 .and. nstep == 0) then     
      
      ti(1:ncol,1:pver) = tN(1:ncol,1:pver)
      tE(1:ncol,1:pver) = tN(1:ncol,1:pver)
    
    else

      tE(1:ncol,1:pver) = MAX(tN(1:ncol,1:pver),tE(1:ncol,1:pver))
      tE(1:ncol,1:pver) = MIN(temax,tE(1:ncol,1:pver))

      ti(1:ncol,1:pver) = MAX(tN(1:ncol,1:pver),ti(1:ncol,1:pver))
      ti(1:ncol,1:pver) = MIN(ti(1:ncol,1:pver),tE(1:ncol,1:pver))

      tE(1:ncol,ionBotP:pver) = tN(1:ncol,ionBotP:pver)
      ti(1:ncol,ionBotP:pver) = tN(1:ncol,ionBotP:pver)

    endif

    tE0(1:ncol,1:pver) = tE(1:ncol,1:pver)

    wrk2(1:ncol,1:ionBot) =  ndensE(1:ncol,1:ionBot)/wrk1/(SIN(dipMag(1:ncol,1:ionBot)))**2._r8
    
!-----------------------------------------------------------------------------
!  Get constant terms needed for loss term g3 for electron temperature update 
!-----------------------------------------------------------------------------
    lossc5  = 1.21E-4_r8
    lossc7  = 3.6E-2_r8
    lossc9  = 5.7E-4_r8
    lossc13 = 7.E-5_r8
    
!-----------------------------------------------------------------------------
!  Get terms needed for loss term g3 for electron temperature update which do 
!  not need to be updated in iteration loop.  
!-----------------------------------------------------------------------------
    do iCol = 1, ncol

      if (.not. colConv(iCol)) then

        do iVer = 1, ionBot

          lossc4p(iCol,iVer)  = 1.77E-19_r8 * ndensN2(iCol,iVer) * ndensE(iCol,iVer)
          lossc6p(iCol,iVer)  = 1.21E-18_r8 * ndensO2(iCol,iVer) * ndensE(iCol,iVer)
          lossc8p(iCol,iVer)  = 7.9E-19_r8 * ndensO1(iCol,iVer) * ndensE(iCol,iVer)
          lossc10p(iCol,iVer) = 1.3E-4_r8 * ndensN2(iCol,iVer) * ndensE(iCol,iVer)
          lossc11p(iCol,iVer) = 3.125E-21_r8 * ndensO2(iCol,iVer) * ndensE(iCol,iVer)
          lossc12p(iCol,iVer) = 3.4E-12_r8 * ndensO1(iCol,iVer) * ndensE(iCol,iVer)
          lossc14p(iCol,iVer) = 1.57E-12_r8 * ndensO1(iCol,iVer) * ndensE(iCol,iVer)
          lossc15p(iCol,iVer) = 2.9E-14_r8 * ndensN2(iCol,iVer) * ndensE(iCol,iVer)
          lossc16p(iCol,iVer) = 6.9E-14_r8 * ndensO2(iCol,iVer) * ndensE(iCol,iVer)
          lossc3p(iCol,iVer)  = 3.2E-8_r8 * 15._r8 * (ndensOP(iCol,iVer) + &
                              0.5_r8 * ndensO2P(iCol,iVer) + 0.53_r8 * ndensNOP(iCol,iVer)) * ndensE(iCol,iVer)

          losscin(iCol,iVer) = (6.6e-14_r8*ndensN2(iCol,iVer) + 5.8e-14_r8*ndensO2(iCol,iVer)                    &
                             + 0.21e-14_r8*ndensO1(iCol,iVer)*sqrt2*sqrtTN(iCol,iVer))*ndensOP(iCol,iVer)    &
                             +(5.9e-14_r8*ndensN2(iCol,iVer) + 5.45e-14_r8*ndensO2(iCol,iVer)                   &
                             + 4.5e-14_r8*ndensO1(iCol,iVer))*ndensOP(iCol,iVer)                             &
                             +(5.8e-14_r8*ndensN2(iCol,iVer) + 0.14e-14_r8*ndensO2(iCol,iVer)*sqrtTN(iCol,iVer) &
                             + 4.4e-14_r8*ndensO1(iCol,iVer)) * ndensO2P(iCol,iVer)

        enddo !iVer loop

!----------------------------------------------------------------------------------
!  Calculate upper boundary heat flux 
!----------------------------------------------------------------------------------
        if (ABS(dipMagD(iCol,1)) < 60.0_r8) FeDB = 0.5_r8 * &
                                        (1._r8 + SIN(pi * (ABS(dipMagD(iCol,1)) - 30.0_r8) /60.0_r8))

        if (ABS(dipMagD(iCol,1)) >= 60.0_r8) FeDB = 1._r8 

        FeD = -4.5E7_r8 * f107 * 0.5_r8 * FeDB
        FeN = .2_r8 * FeD

!---------------------------------------------------
!  Set upper boundary condition for right hand side
!---------------------------------------------------
        if (zenAngD(iCol) <= 80.0_r8) Feub(iCol) = FeD
        if (zenAngD(iCol) > 80.0_r8 .AND. zenAngD(iCol) < 100.0_r8) Feub(iCol) = 0.5_r8 * (FeD + FeN) &
                                                                         + 0.5_r8 * (FeD - FeN) * &
                                                                 COS(pi * ((zenAngD(iCol) - 80.0_r8) / 20.0_r8))
        if (zenAngD(iCol) >= 100.0_r8) Feub(iCol) = FeN

!------------------------------------------------------------------------------------------
!  Calculate thickness terms for vertical derivative
!------------------------------------------------------------------------------------------
        do iVer = 1, ionBot

          delZ(iCol,iVer) = (pInt(iCol,iVer+1) - &
               pInt(iCol,iVer)) * rairv(iCol,iVer,lchnk) * tN(iCol,iVer) / pMid(iCol,iVer) / gravit

        enddo

        do iVer = 2, ionBotP		! Assuming ionBotP < pverp
          delZi(iCol,iVer) = (pMid(iCol,iVer) - &
               pMid(iCol,iVer-1)) * rairvi(iCol,iVer) * tNInt(iCol,iVer) / pInt(iCol,iVer) / gravit
        enddo
        delZi(iCol,1) = 1.5_r8*delZi(iCol,2) - .5_r8*delZi(iCol,3)

!----------------------------------------------------------
!  Convert delZ variables from meters to centimeters
!----------------------------------------------------------
        delZi(iCol,1:ionBotP) = delZi(iCol,1:ionBotP)*100._r8
        delZ(iCol,1:ionBot) = delZ(iCol,1:ionBot)*100._r8
  
      endif ! Column not converged

    enddo !iCol loop

!-------------------------------------------------------------------------------------------------------
!  Iterate to calculate new electron temperature. 
!  Time splitting is used: first solve the heating/cooling equation, then solve the diffusion equations.
!  Also, set convergence flag to false and iterate until true or 6 iterations, whichever comes first 
!-------------------------------------------------------------------------------------------------------
    converged = .false. 
    iter = 0
    do while (.not. converged .and. iter < maxIter)
    
!--------------------------------------------------------------------------------------------------------
!  Increment iteration loop counter and save electron temperature from previous iteration for convergence 
!  test at end of interation loop.  Also, take square root of electron temperature to be used later
!--------------------------------------------------------------------------------------------------------        
      iter = iter + 1
      
      tEPrevI(1:ncol,1:ionBot) = tE(1:ncol,1:ionBot)

      sqrtTE(1:ncol,1:ionBot) = SQRT(tE(1:ncol,1:ionBot))

!--------------------------------------------------------------------------------------------------------
!  Loop over columns then vertical levels and call tridiagonal solver for each column to get electron 
!  temperature
!--------------------------------------------------------------------------------------------------------
      do iCol = 1, ncol

        if (.not. colConv(iCol)) then

          do iVer = 1, ionBot
 
!-----------------------------------------------------------------------------
!  Get loss term g3 for electron temperature update.  Need to calculate 
!  constituent dependent loss terms which make up g3
!-----------------------------------------------------------------------------
            lossceN2(iCol,iVer) = lossc4p(iCol,iVer) * (1._r8 - lossc5 * tE(iCol,iVer)) * tE(iCol,iVer)
            lossceO2(iCol,iVer) = lossc6p(iCol,iVer) * (1._r8 + lossc7 * sqrtTE(iCol,iVer)) * sqrtTE(iCol,iVer)
            lossceO1(iCol,iVer) = lossc8p(iCol,iVer) * (1._r8 + lossc9 * tE(iCol,iVer)) * sqrtTE(iCol,iVer)

            if (tE(iCol,iVer) < 1000.0_r8) losscA(iCol,iVer) = 5.71E-8_r8 * EXP(-3352.6_r8 / tE(iCol,iVer)) 
            if (tE(iCol,iVer) >= 1000.0_r8 .AND. tE(iCol,iVer) <= 2000.0_r8) &
                                          losscA(iCol,iVer) = 2.0E-7_r8 * EXP(-4605.2_r8 / tE(iCol,iVer)) 
            if (tE(iCol,iVer) > 2000.0_r8) losscA(iCol,iVer) = 2.53E-6_r8 * sqrtTE(iCol,iVer) * &
                                                                          EXP(-17620._r8 / tE(iCol,iVer))

            tENDiff(iCol,iVer) = tE(iCol,iVer) -  tN(iCol,iVer)

            if (ABS(tENDiff(iCol,iVer)) < 0.1_r8) tENDiff(iCol,iVer) = 0.1_r8

            lossceN2v(iCol,iVer) = lossc10p(iCol,iVer) * losscA(iCol,iVer) * &
               (1._r8 - EXP(3200._r8 * (1._r8 / tE(iCol,iVer) - 1._r8 / tN(iCol,iVer)))) / tENDiff(iCol,iVer)
            lossceO2v(iCol,iVer) = lossc11p(iCol,iVer) * tE(iCol,iVer)**2
            lossceOf(iCol,iVer) = lossc12p(iCol,iVer) * (1._r8 - lossc13 * tE(iCol,iVer)) * &
                                               (0.4_r8 + 150._r8 / tE(iCol,iVer)) / tN(iCol,iVer)
            losscf2d(iCol,iVer) = 2.4E+4_r8 + 0.3_r8 * (tE(iCol,iVer) - 1500._r8) + &
                          1.947E-5_r8 * (tE(iCol,iVer) - 1500._r8) * (tE(iCol,iVer) - 4000._r8)
            losscf2(iCol,iVer) = losscf2d(iCol,iVer) * (1._r8 / 3000._r8 - 1._r8 / tE(iCol,iVer))
            losscf3(iCol,iVer) = -22713._r8 * (1._r8 / tN(iCol,iVer) - 1._r8 / tE(iCol,iVer))
            lossceO1D(iCol,iVer) = lossc14p(iCol,iVer) * EXP(losscf2(iCol,iVer)) * & 
                                                   (1._r8 - EXP(losscf3(iCol,iVer))) / tENDiff(iCol,iVer)
            lossceN2Rot(iCol,iVer) = lossc15p(iCol,iVer) / sqrtTE(iCol,iVer)
            lossceO2Rot(iCol,iVer) = lossc16p(iCol,iVer) / sqrtTE(iCol,iVer)
            losscei(iCol,iVer) = lossc3p(iCol,iVer) / tE(iCol,iVer)**1.5_r8
            lossg3(iCol,iVer) =  lossceN2(iCol,iVer) + lossceO2(iCol,iVer) + lossceO1(iCol,iVer) + lossceN2v(iCol,iVer)   &
                             + lossceO2v(iCol,iVer) + lossceOf(iCol,iVer) + lossceO1D(iCol,iVer)                        &
                             + lossceN2Rot(iCol,iVer) + lossceO2Rot(iCol,iVer)

            tE(iCol,iVer) = (wrk2(iCol,iVer)/ztodt * tE0(iCol,iVer) + (lossg3(iCol,iVer) * tN(iCol,iVer)                  &
                             + losscei(iCol,iVer) * ti(iCol,iVer) + sourceg4(iCol,iVer))                                &
                             /(SIN(dipMag(iCol,iVer)))**2)   /                 &
                          (wrk2(iCol,iVer)/ztodt + (lossg3(iCol,iVer) + losscei(iCol,iVer))/(SIN(dipMag(iCol,iVer)))**2._r8)

!-------------------------------------------------------------------------------------------------------------------
!  Limit maximum value of electron temperature to be a locally set maximum (temax) and minimum value of electron 
!  temperature to be at least the neutral temperature for a given column and level
!-------------------------------------------------------------------------------------------------------------------      
            tE(iCol,iVer) = min(temax,tE(iCol,iVer))
            tE(iCol,iVer) = max(tN(iCol,iVer),tE(iCol,iVer))

          enddo !iVer loop
	  
	endif ! Column not converged
	  
      enddo ! End of column loop

! fvitt -- make sqrtTE consistent with new tE
      sqrtTE(1:ncol,1:ionBot) = SQRT(tE(1:ncol,1:ionBot))

!-----------------------------------------------------
!  Calculate thermal conductivity of electron gas   
!-----------------------------------------------------
      do iCol = 1, ncol
      
        if (.not. colConv(iCol)) then
      
          do iVer = 1, ionBot

            f1Ted1 = 2.82E-17_r8 * sqrtTE(iCol,iVer) - 3.41E-21_r8 * tE(iCol,iVer)**1.5_r8
            f1Ted2 = 2.2E-16_r8 + 7.92E-18_r8 * sqrtTE(iCol,iVer)
            f1Ted3 = 1.1E-16_r8 * (1._r8 + 5.7E-4_r8 * tE(iCol,iVer))

            f1Te = ndensN2(iCol,iVer) / ndensE(iCol,iVer) * f1Ted1 + ndensO2(iCol,iVer) / &
                        ndensE(iCol,iVer) * f1Ted2 + ndensO1(iCol,iVer) / ndensE(iCol,iVer) * f1Ted3

!-----------------------------------------------------------------------------
!  Calculate electron conductivity using parameters set in module and f1(Te)
!-----------------------------------------------------------------------------
            Ke(iVer) = Kec1 * tE(iCol,iVer)**2.5_r8 / (1._r8 + Kec2 * tE(iCol,iVer)**2._r8 * f1Te)

          enddo !iVer loop

!----------------------------------------------------------------------
!  Get electron conductivity at interface levels to be used later
!----------------------------------------------------------------------
          do iVer = 2,ionBot
            Kei(iVer) = SQRT(Ke(iVer-1)*Ke(iVer))
          enddo
          Kei(1) = 1.5_r8*Kei(2)-.5_r8*Kei(3)
          Kei(ionBotP) = 1.5_r8*Kei(ionBot)-.5_r8*Kei(ionBot-1)

!------------------------------------------------------------------------------------------------------
!  Derive subdiagonal, superdiagonal, and diagonal as input to solver for electron temperature tendency
!------------------------------------------------------------------------------------------------------
          do iVer = 2, ionBot-1
              subDiag(iVer) = -Kei(iVer) / delZi(iCol,iVer) / delZ(iCol,iVer)
              superDiag(iVer) = -Kei(iVer+1) / delZi(iCol,iVer+1) / delZ(iCol,iVer)
              diag(iVer) = wrk2(iCol,iVer)/ztodt-subDiag(iVer)-superDiag(iVer)
              rHS(iVer) = tE(iCol,iVer) * wrk2(iCol,iVer)/ztodt
          enddo !iVer loop

!-------------------------------------------------------------------------------------
!  Calculate diagonal, superdiagonal, and right hand side upper boundary values
!-------------------------------------------------------------------------------------
          superDiag(1)  = -Kei(2) / delZi(iCol,2) / delZ(iCol,1)
          diag(1) = wrk2(iCol,1)/ztodt - superDiag(1)
          rHS(1) = tE(iCol,1) * wrk2(iCol,1)/ztodt-Feub(iCol) / delZ(iCol,1)

!---------------------------------------------------------------------------------------------
!  Calculate subdiagonal, diagonal, superdiagonal, and right hand side lower boundary values
!---------------------------------------------------------------------------------------------
          subDiag(ionBot) = -Kei(ionBot) / delZi(iCol,ionBot) / delZ(iCol,ionBot)
          superDiag(ionBot) = -Kei(ionBotP) / delZi(iCol,ionBotP) / delZ(iCol,ionBot)
          diag(ionBot) = wrk2(iCol,ionBot)/ztodt-subDiag(ionBot)-superDiag(ionBot)
          rHS(ionBot) = tE(iCol,ionBot) * wrk2(iCol,ionBot)/ztodt -superDiag(ionBot)*tN(iCol,ionBotP)

!-------------------------------------------------
! Call solver to get electron temperature update
!-------------------------------------------------
          call tridag(subDiag,diag,superDiag,rHS,tETemp,ionBot)

          tE(iCol,1:ionBot) = tETemp(1:ionBot)
          do iVer = 1,ionBot
             tE(iCol,iVer) = min(temax,tE(iCol,iVer))
             tE(iCol,iVer) = max(tN(iCol,iVer),tE(iCol,iVer))
          enddo
!---------------------------------------------------------------------------------------------------------
!  Calculate ion temperature from electron temperature, ion-neutral and electron-ion loss terms, neutral 
!  temperature, mass density and joule heating.  Set minimum value to neutral temperature and maximum 
!  value to electron temperature for each column and vertical level
!---------------------------------------------------------------------------------------------------------
          do iVer = 1,ionBot
              ti(iCol,iVer) = (losscei(iCol,iVer) * tE(iCol,iVer) + losscin(iCol,iVer) * tN(iCol,iVer) +   &
                               rho(iCol,iVer) * qjoule(iCol,iVer))/(losscei(iCol,iVer) + losscin(iCol,iVer))
              ti(iCol,iVer) = max(tN(iCol,iVer),ti(iCol,iVer))
              ti(iCol,iVer) = min(tE(iCol,iVer),ti(iCol,iVer))
          enddo
      
!--------------------------------------------------------------------------------------------------------
! Check for convergence which is a change of electron temperature ratio to previous loop for all levels
! and columns of less than 0.05K.  Had to modify this to do convergence check on each column since 
! checking all columns in a chunk gives different answers depending on number of tasks and tasks per node.
!--------------------------------------------------------------------------------------------------------
          if (ALL(ABS(tE(iCol,1:ionBot) / tEPrevI(iCol,1:ionBot) - 1._r8) < 0.05_r8)) then 
	
	    colConv(iCol) = .true.

          endif

        endif ! Column not converged

      enddo ! iCol loop
!--------------------------------------------------------------
!  Check to see if all columns have converged and set flag
!--------------------------------------------------------------
      if (ALL(colConv(1:ncol))) converged = .true.

    enddo ! End of iteration loop
      
    write(iulog,*)'lchnk, Number of tempei iterations, converged flag ', lchnk, iter, converged

!------------------------------------------------------------------
!  Store calculated electron and ion temperatures in pbuf structure
!------------------------------------------------------------------
    call pbuf_get_field(pbuf, indxTe, tE, start=(/1,1/), kount=(/ncol,pver/))
    call pbuf_get_field(pbuf, indxTi, ti, start=(/1,1/), kount=(/ncol,pver/))

!--------------------------------------------------------------------------------------------------------
! Calculate electron-neutral heating and electron-ion Coulomb heating.  Then update dry static energy.
!--------------------------------------------------------------------------------------------------------
    do iVer = 1, ionBot
       do iCol = 1, ncol
          sqrtTE(iCol,iVer) = SQRT(tE(iCol,iVer))
          lossceN2(iCol,iVer) = lossc4p(iCol,iVer) * (1._r8 - lossc5 * tE(iCol,iVer)) * tE(iCol,iVer)
          lossceO2(iCol,iVer) = lossc6p(iCol,iVer) * (1._r8 + lossc7 * sqrtTE(iCol,iVer)) * sqrtTE(iCol,iVer)
          lossceO1(iCol,iVer) = lossc8p(iCol,iVer) * (1._r8 + lossc9 * tE(iCol,iVer)) * sqrtTE(iCol,iVer)
       enddo
    enddo
    qen(1:ncol,1:ionBot) = (lossceN2(1:ncol,1:ionBot)+lossceO2(1:ncol,1:ionBot)+lossceO1(1:ncol,1:ionBot)) *  &
                           (tE(1:ncol,1:ionBot)-tN(1:ncol,1:ionBot)) / rho(1:ncol,1:ionBot)
    qei(1:ncol,1:ionBot) = losscei(1:ncol,1:ionBot) * (tE(1:ncol,1:ionBot)-ti(1:ncol,1:ionBot)) / rho(1:ncol,1:ionBot)
    ptend%s(1:ncol,1:ionBot) = ptend%s(1:ncol,1:ionBot) + (qei(1:ncol,1:ionBot)+qen(1:ncol,1:ionBot))/6.24E15_r8    

    call outfld ('TE'            , tE           , pcols, lchnk)
    call outfld ('TI'            , ti           , pcols, lchnk)
    call outfld ('LOSS_g3'      , lossg3       , pcols, lchnk)
    call outfld ('LOSS_EI'      , losscei       , pcols, lchnk)
    call outfld ('LOSS_IN'      , losscin       , pcols, lchnk)
    call outfld ('QIN'            , qei           , pcols, lchnk)

    return

  end subroutine ionos_tempei


!-----------------------------------------------------------------------
! Simple tridiagonal solver routine
!-----------------------------------------------------------------------
  SUBROUTINE tridag(a,b,c,r,u,n)

    use cam_abortutils,   only: endrun

    INTEGER n
    REAL(r8) :: a(n),b(n),c(n),r(n),u(n)
    INTEGER j
    REAL(r8) :: bet,gam(n)

    if(b(1).eq.0._r8) call endrun('ionosphere: bt(1)=0 in tridag')
    bet=b(1)
    u(1)=r(1)/bet
    do j=2,n
      gam(j)=c(j-1)/bet
      bet=b(j)-a(j)*gam(j)
      if(bet.eq.0._r8) call endrun('ionosphere: bet=0 in tridag')
      u(j)=(r(j)-a(j)*u(j-1))/bet
    end do

    do j=n-1,1,-1
      u(j)=u(j)-gam(j+1)*u(j+1)
    end do

    return

  END SUBROUTINE tridag

end module ionosphere
