module clm_bgc_interfaceMod
!#define FLEXIBLE_POOLS
!----------------------------------------------------------------------------------------------
! CLM-PFLOTRAN soil bgc coupling interface
! authors:
!
!          Climate Change Science Institute & Environmental Science Division
!          Oak Ridge National Laboratory
!
! date: 8/26/2015

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  !BOP
  !
  ! !MODULE: clm_bgc_interfaceMod
  !
  ! !DESCRIPTION:
  ! Performs
  !
  ! !USES:

  ! clm g/l/c/p constants
  use shr_log_mod         , only : errMsg => shr_log_errMsg
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use GridcellType        , only : grc
  use LandunitType        , only : lun
  use ColumnType          , only : col
  use PatchType           , only : pft

  use decompMod           , only : bounds_type

  ! (dummy) variable definitions
  use atm2lndType         , only : atm2lnd_type
  use SoilStateType       , only : soilstate_type
  use WaterStateType      , only : waterstate_type
  use WaterFluxType       , only : waterflux_type
  use SoilHydrologyType   , only : soilhydrology_type
  use TemperatureType     , only : temperature_type
  use EnergyFluxType      , only : energyflux_type
  use SoilWaterRetentionCurveMod, only : soil_water_retention_curve_type

  use CNStateType         , only : cnstate_type
  use CNCarbonFluxType    , only : carbonflux_type
  use CNCarbonStateType   , only : carbonstate_type
  use CNNitrogenFluxType  , only : nitrogenflux_type
  use CNNitrogenStateType , only : nitrogenstate_type

  use ch4Mod              , only : ch4_type
  ! bgc interface: use_clm_bgc
  use PhotosynthesisType     , only : photosyns_type
  use cropType               , only : crop_type
  use CanopyStateType        , only : canopystate_type
  use PhosphorusStateType    , only : phosphorusstate_type
  use PhosphorusFluxType     , only : phosphorusflux_type

  use clm_bgc_interface_data , only : clm_bgc_interface_data_type


!  ! PFLOTRAN thc module switchs for coupling with CLM45-CN
!  use clm_varctl          , only : use_pflotran, pf_tmode, pf_hmode, pf_cmode
!  use clm_varctl          , only : pf_surfaceflow, pf_frzmode
!  use clm_varctl          , only : initth_pf2clm

  ! most used constants in this module
  use clm_varpar          , only : nlevsoi, nlevgrnd, nlevdecomp, nlevdecomp_full
  use clm_varpar          , only : ndecomp_pools
  use clm_varpar          , only : max_patch_per_col
  use clm_varcon          , only : denh2o, denice, tfrz, dzsoi_decomp
  use landunit_varcon     , only : istsoil, istcrop

  ! misc.
  use abortutils          , only : endrun


!#ifdef CLM_PFLOTRAN
!  use clm_pflotran_interface_data
!  use pflotran_clm_main_module
!  use pflotran_clm_setmapping_module
!#endif

  ! !PUBLIC TYPES:
  implicit none

  save

  private    ! By default everything is private

!#ifdef CLM_PFLOTRAN
!  type(pflotran_model_type), pointer, public :: pflotran_m
!#endif
!  !
!  character(len=256), private:: pflotran_prefix = ''
!  character(len=32), private :: restart_stamp = ''

  real(r8), parameter :: rgas = 8.3144621d0                 ! m3 Pa K-1 mol-1

  ! !PUBLIC MEMBER FUNCTIONS:
!  public :: clm_pf_readnl

  ! wrappers around '#ifdef CLM_PFLOTRAN .... #endif' block statements to maintain sane runtime behavior
  ! when pflotran is not available.
!  public :: clm_pf_interface_init
!  public :: clm_pf_set_restart_stamp
!  public :: clm_pf_run
!  public :: clm_pf_write_restart
!  public :: clm_pf_finalize


  !! run clm bgc (CN or BGC) through interface
  public    :: clm_bgc_run
  !! update clm variables by pflotran
  public    :: update_bgc_data_pf2clm

  !! pass clm variables to clm_bgc_data
  public    :: get_clm_bgc_data

  !! pass clm variables to clm_bgc_data, called by get_clm_bgc_data
  private   :: get_clm_soil_property
  private   :: get_clm_soil_thermohydro
  private   :: get_clm_bgc_state
  private   :: get_clm_bgc_flux

  !! update clm variables by clm_bgc_data,
  !! specific bgc-module may require certain combination of these subroutines
  private   :: update_bgc_state_decomp
  private   :: update_bgc_state_smin
  private   :: update_bgc_flux_decomp_sourcesink
  private   :: update_bgc_flux_decomp_cascade
  private   :: update_bgc_flux_smin
  private   :: update_bgc_flux_gas_pf
  private   :: update_soil_moisture
  private   :: update_soil_temperature


contains

!-----------------------------------------------------------------------
!
! public interface functions allowing runtime behavior regardless of
! whether pflotran is compiled in.
!
!-----------------------------------------------------------------------
  subroutine get_clm_bgc_data(clm_bgc_data,bounds,       &
           num_soilc, filter_soilc,                               &
           num_soilp, filter_soilp,                               &
           atm2lnd_vars,                                          &
           waterstate_vars, waterflux_vars,                       &
           soilstate_vars,  temperature_vars, energyflux_vars,    &
           soilhydrology_vars, soil_water_retention_curve,        &
           cnstate_vars, carbonflux_vars, carbonstate_vars,       &
           nitrogenflux_vars, nitrogenstate_vars,                 &
           phosphorusflux_vars, phosphorusstate_vars,             &
           ch4_vars                                               &
           )

    implicit none

    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:)   ! filter for soil columns
    integer                  , intent(in)    :: num_soilp         ! number of soil patches in filter
    integer                  , intent(in)    :: filter_soilp(:)   ! filter for soil patches
    type(atm2lnd_type)       , intent(in)    :: atm2lnd_vars
    type(waterstate_type)    , intent(inout) :: waterstate_vars
    type(waterflux_type)     , intent(inout) :: waterflux_vars
    type(soilstate_type)     , intent(inout) :: soilstate_vars
    type(temperature_type)   , intent(inout) :: temperature_vars
    type(soilhydrology_type) , intent(inout) :: soilhydrology_vars
    type(energyflux_type)    , intent(inout) :: energyflux_vars
    class(soil_water_retention_curve_type), intent(in) :: soil_water_retention_curve

    type(cnstate_type)          , intent(inout) :: cnstate_vars
    type(carbonflux_type)       , intent(inout) :: carbonflux_vars
    type(carbonstate_type)      , intent(inout) :: carbonstate_vars
    type(nitrogenflux_type)     , intent(inout) :: nitrogenflux_vars
    type(nitrogenstate_type)    , intent(inout) :: nitrogenstate_vars
    type(phosphorusflux_type)   , intent(inout) :: phosphorusflux_vars
    type(phosphorusstate_type)  , intent(inout) :: phosphorusstate_vars
    type(ch4_type)              , intent(inout) :: ch4_vars
    type(clm_bgc_interface_data_type), intent(inout) :: clm_bgc_data

    !-----------------------------------------------------------------------

    character(len=256) :: subname = "get_clm_bgc_data"

    call get_clm_soil_property(clm_bgc_data,                &
                    bounds, num_soilc, filter_soilc,        &
                    soilstate_vars)

    call get_clm_soil_thermohydro(clm_bgc_data,             &
                   bounds, num_soilc, filter_soilc,         &
                   atm2lnd_vars, soilstate_vars,            &
                   waterstate_vars, waterflux_vars,         &
                   temperature_vars, energyflux_vars,       &
                   soil_water_retention_curve)

    call get_clm_bgc_state(clm_bgc_data,                    &
                    bounds, num_soilc, filter_soilc,        &
                    carbonstate_vars, nitrogenstate_vars,   &
                    phosphorusstate_vars)

    call get_clm_bgc_flux(clm_bgc_data,                     &
                    bounds, num_soilc, filter_soilc,        &
                    cnstate_vars, carbonflux_vars,          &
                    nitrogenflux_vars, phosphorusflux_vars)

  end subroutine get_clm_bgc_data
  !-----------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: clm_pf_readnl
  !
  ! !INTERFACE:
!  subroutine clm_pf_readnl( NLFilename )
!  !
!  ! !DESCRIPTION:
!  ! Read namelist for clm-pflotran interface
!  !
!  ! !USES:
!    use clm_varctl    , only : iulog
!    use spmdMod       , only : masterproc, mpicom
!    use fileutils     , only : getavu, relavu, opnfil
!    use clm_nlUtilsMod, only : find_nlgroup_name
!    use shr_nl_mod    , only : shr_nl_find_group_name
!    use shr_mpi_mod   , only : shr_mpi_bcast
!
!    implicit none
!
!  ! !ARGUMENTS:
!    character(len=*), intent(IN) :: NLFilename ! Namelist filename
!  ! !LOCAL VARIABLES:
!    integer :: ierr                 ! error code
!    integer :: unitn                ! unit for namelist file
!    character(len=32) :: subname = 'clm_pf_readnl'  ! subroutine name
!  !EOP
!  !-----------------------------------------------------------------------
!    namelist / clm_pflotran_inparm / pflotran_prefix
!
!    ! ----------------------------------------------------------------------
!    ! Read namelist from standard namelist file.
!    ! ----------------------------------------------------------------------
!
!    if ( masterproc )then
!
!       unitn = getavu()
!       write(iulog,*) 'Read in clm-pflotran namelist'
!       call opnfil (NLFilename, unitn, 'F')
!       call shr_nl_find_group_name(unitn, 'clm_pflotran_inparm', status=ierr)
!       if (ierr == 0) then
!          read(unitn, clm_pflotran_inparm, iostat=ierr)
!          if (ierr /= 0) then
!             call endrun(msg=subname //':: ERROR: reading clm_pflotran_inparm namelist.'//&
!                         errMsg(__FILE__, __LINE__))
!          end if
!       end if
!       call relavu( unitn )
!       write(iulog, '(/, A)') " clm-pflotran namelist:"
!       write(iulog, '(A, " : ", A,/)') "   pflotran_prefix", trim(pflotran_prefix)
!    end if
!
!    ! Broadcast namelist variables read in
!    call shr_mpi_bcast(pflotran_prefix, mpicom)
!
!  end subroutine clm_pf_readnl
!
!  !-----------------------------------------------------------------------
!  !BOP
!  !
!  ! !IROUTINE: clm_pf_set_restart_stamp
!  !
!  ! !INTERFACE:
!  subroutine clm_pf_set_restart_stamp(clm_restart_filename)
!  !
!  ! !DESCRIPTION: Set the pflotran restart date stamp. Note we do NOT
!  ! restart here, that gets handled by pflotran's internal
!  ! initialization during interface_init_clm_pf()
!  !
!  ! !USES:
!  ! !ARGUMENTS:
!    character(len=256), intent(in) :: clm_restart_filename
!  ! !LOCAL VARIABLES:
!    integer :: name_length, start_pos, end_pos
!    character(len=32) :: clm_stamp
!  !EOP
!  !-----------------------------------------------------------------------
!
!    ! clm restart file name is of the form:
!    !     ${CASE_NAME}.clm2.r.YYYY-MM-DD-SSSSS.nc
!    ! we need to extract the: YYYY-MM-DD-SSSSS
!    write(*, '("clm-pf : clm restart file name : ", A/)') trim(clm_restart_filename)
!    name_length = len(trim(clm_restart_filename))
!    start_pos = name_length - 18
!    end_pos = name_length - 3
!    clm_stamp = clm_restart_filename(start_pos : end_pos)
!    write(*, '("clm-pf : clm date stamp : ", A/)') trim(clm_stamp)
!    restart_stamp = clm_stamp
!
!  end subroutine clm_pf_set_restart_stamp
!
!
!  !-----------------------------------------------------------------------------
!  !BOP
!  !
!  ! !IROUTINE: pflotran_not_available
!  !
!  ! !INTERFACE:
!  subroutine pflotran_not_available(subname)
!  !
!  ! !DESCRIPTION:
!  ! Print an error message and abort.
!  !
!  ! !USES:
!
!  ! !ARGUMENTS:
!    implicit none
!    character(len=*), intent(in) :: subname
!  ! !LOCAL VARIABLES:
!  !EOP
!  !-----------------------------------------------------------------------
!    call endrun(trim(subname) // ": ERROR: CLM-PFLOTRAN interface has not been compiled " // &
!         "into this version of CLM.")
!  end subroutine pflotran_not_available
!
!
!!******************************************************************************************!
!!
!! public interface function wrappers
!!
!!------------------------------------------------------------------------------------------!
!
!  !-----------------------------------------------------------------------------
!  subroutine clm_pf_interface_init(bounds)
!
!    implicit none
!
!    type(bounds_type), intent(in) :: bounds  ! bounds
!
!    character(len=256) :: subname = "clm_pf_interface_init()"
!
!#ifdef CLM_PFLOTRAN
!    call interface_init(bounds)
!#else
!    call pflotran_not_available(subname)
!#endif
!  end subroutine clm_pf_interface_init
!
!  !--------------------------------------------------------------------------------------------
!
!  subroutine clm_pf_run(bounds,       &
!          ! pflotran only works for 'soilc', i.e. (natural/cropped soil columns)
!           ! at this coding stage
!           num_soilc, filter_soilc,                               &
!           num_soilp, filter_soilp,                               &
!           ! soil thermal-hydrology (TODO: will update when testing with th coupling)
!           atm2lnd_vars,                                          &
!           waterstate_vars, waterflux_vars,                       &
!           soilstate_vars,  temperature_vars, energyflux_vars,    &
!           soilhydrology_vars, soil_water_retention_curve,        &
!           ! soil bgc
!           cnstate_vars, carbonflux_vars, carbonstate_vars,       &
!           nitrogenflux_vars, nitrogenstate_vars,                 &
!           ch4_vars                                               &
!           )
!
!    implicit none
!
!    ! !ARGUMENTS:
!    type(bounds_type)        , intent(in)    :: bounds
!    integer                  , intent(in)    :: num_soilc         ! number of soil columns in filter
!    integer                  , intent(in)    :: filter_soilc(:)   ! filter for soil columns
!    integer                  , intent(in)    :: num_soilp         ! number of soil patches in filter
!    integer                  , intent(in)    :: filter_soilp(:)   ! filter for soil patches
!    type(atm2lnd_type)       , intent(in)    :: atm2lnd_vars
!    type(waterstate_type)    , intent(inout) :: waterstate_vars
!    type(waterflux_type)     , intent(inout) :: waterflux_vars
!    type(soilstate_type)     , intent(inout) :: soilstate_vars
!    type(temperature_type)   , intent(inout) :: temperature_vars
!    type(soilhydrology_type) , intent(inout) :: soilhydrology_vars
!    type(energyflux_type)    , intent(inout) :: energyflux_vars
!    class(soil_water_retention_curve_type), intent(in) :: soil_water_retention_curve
!
!    type(cnstate_type)       , intent(inout) :: cnstate_vars
!    type(carbonflux_type)    , intent(inout) :: carbonflux_vars
!    type(carbonstate_type)   , intent(inout) :: carbonstate_vars
!    type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars
!    type(nitrogenstate_type) , intent(inout) :: nitrogenstate_vars
!    type(ch4_type)           , intent(inout) :: ch4_vars
!
!    !-----------------------------------------------------------------------
!
!    character(len=256) :: subname = "clm_pf_run"
!
!#ifdef CLM_PFLOTRAN
!    call pflotran_run_onestep(bounds,                             &
!           num_soilc, filter_soilc,                               &
!           num_soilp, filter_soilp,                               &
!           ! soil thermal-hydrology (TODO: will update when testing with th coupling)
!           atm2lnd_vars,                                          &
!           waterstate_vars, waterflux_vars,                       &
!           soilstate_vars,  temperature_vars, energyflux_vars,    &
!           soilhydrology_vars, soil_water_retention_curve,        &
!           ! soil bgc
!           cnstate_vars, carbonflux_vars, carbonstate_vars,       &
!           nitrogenflux_vars, nitrogenstate_vars,                 &
!           ch4_vars                                               &
!           )
!#else
!    call pflotran_not_available(subname)
!#endif
!  end subroutine clm_pf_run
!
!  !-----------------------------------------------------------------------------
!  subroutine clm_pf_write_restart(date_stamp)
!
!    implicit none
!    character(len=*), intent(in) :: date_stamp
!
!    character(len=32) :: subname = "clm_pf_write_restart"
!
!#ifdef CLM_PFLOTRAN
!    call pflotran_write_checkpoint(date_stamp)
!#else
!    call pflotran_not_available(subname)
!#endif
!  end subroutine clm_pf_write_restart
!
!
!  !-----------------------------------------------------------------------------
!  !BOP
!  !
!  ! !ROUTINE: clm_pf_finalize
!  !
!  ! !INTERFACE:
!  subroutine clm_pf_finalize()
!
!    implicit none
!    character(len=256) :: subname = "clm_pf_finalize"
!
!#ifdef CLM_PFLOTRAN
!    call pflotran_finalize()
!#else
!    call pflotran_not_available(subname)
!#endif
!  end subroutine clm_pf_finalize




!#ifdef CLM_PFLOTRAN

!************************************************************************************!
! The following is private and requires explicit coupling between CLM and PFLOTRAN
!
!------------------------------------------------------------------------------------!
!

  !-----------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: interface_init
  !
  ! !INTERFACE:
!  subroutine interface_init(bounds)
!    !
!    ! !DESCRIPTION:
!    ! initialize the pflotran iterface
!    !
!    ! !USES:
!    use clm_varctl      , only : iulog
!    use decompMod       , only : get_proc_global, ldecomp
!    use spmdMod         , only : mpicom, masterproc
!    use domainMod       , only : ldomain
!
!    use CNDecompCascadeConType , only : decomp_cascade_con
!
!    use abortutils      , only : endrun
!
!    use clm_time_manager, only : nsstep, nestep
!
!    ! pflotran
!    use Option_module!, only : printErrMsg
!    use Simulation_Base_class, only : simulation_base_type
!    use Simulation_Subsurface_class, only : subsurface_simulation_type
!    use Realization_class, only : realization_type
!    use PFLOTRAN_Constants_module
!    !
!    ! !ARGUMENTS:
!
!    implicit none
!
!#include "finclude/petscsys.h"
!#include "finclude/petscvec.h"
!#include "finclude/petscvec.h90"
!#include "finclude/petscviewer.h"
!
!    !
!    ! !REVISION HISTORY:
!    ! Created by Gautam Bisht
!    ! Revised by Fengming Yuan, CCSI-ORNL
!    !
!    !EOP
!    !
!    ! LOCAL VARAIBLES:
!
!    type(bounds_type), intent(in) :: bounds  ! bounds
!
!    integer  :: numg             ! total number of gridcells across all processors
!    integer  :: numl             ! total number of landunits across all processors
!    integer  :: numc             ! total number of columns across all processors
!    integer  :: nump             ! total number of pfts across all processors
!    integer  :: g,l,c,j          ! indices
!    integer  :: gcount, cellcount
!
!    integer, pointer :: clm_cell_ids_nindex(:)
!    integer, pointer :: clm_2dtop_cell_ids_nindex(:)
!    integer, pointer :: clm_2dbot_cell_ids_nindex(:)
!    integer :: clm_cell_npts
!    integer :: clm_2dtop_npts
!    integer :: clm_2dbot_npts
!
!    !type(pflotran_model_type), pointer:: pflotran_m
!    class(realization_type), pointer  :: realization
!    type(option_type), pointer :: option
!    PetscErrorCode :: ierr
!
!    character(len= 32) :: subname = 'interface_init' ! subroutine name
!
!    associate( &
!         ! Assign local pointers to derived subtypes components (gridcell-level)
!         latdeg     =>  grc%latdeg     , & !  [real(r8) (:)]  latitude (degree)
!         londeg     =>  grc%londeg     , & !  [real(r8) (:)]  longitude (degree)
!         area       =>  grc%area       , & !  [real(r8) (:)]  total land area per gridcell (km^2)
!         gindex     =>  grc%gindex     , & !  [real(r8) (:)]  longitude (degree)
!         ! Assign local pointers to derived subtypes components (landunit-level)
!         ltype      =>  lun%itype      , & !  [integer (:)]  landunit type index
!         ! Assign local pointer to derived subtypes components (column-level)
!         cgridcell  =>  col%gridcell   , & !  [integer (:)]  gridcell index of column
!         clandunit  =>  col%landunit   , & !  [integer (:)]  landunit index of column
!         cwtgcell   =>  col%wtgcell    , & !  [real(r8) (:)]  weight (relative to gridcell
!         ctype      =>  col%itype      , & !  [integer (:)]  column type index
!         topo       =>  col%glc_topo   , &  ! surface elevation (m)
!         micro_sigma=>  col%micro_sigma, &  ! microtopography pdf sigma (m)
!         slope      =>  col%topo_slope , &  ! gridcell topographic slope
!         topo_std   =>  col%topo_std     &  ! gridcell elevation standard deviation
!         )
!
!    !------------------------------------------------------------------------
!
!    do c = bounds%begc, bounds%endc
!      l = col%landunit(c)
!      if (.not.(ltype(l)==istsoil .or. ltype(l)==istcrop) .and. col%active(c)) then
!        write (iulog,*) 'WARNING: Land Unit type of active Columns of non-soil/crop type found within the domain'
!        write (iulog,*) 'CLM-CN -- PFLOTRAN does not support this land unit presently'
!        write (iulog,*) 'So, DEactive the column: ', c
!
!        col%active(c) = .false.
!
!      endif
!
!    enddo ! do c = bounds%begc, bounds%endc
!
!    if (masterproc) then
!      write(iulog,*) '%%-----------------------------------------------------%%'
!      write(iulog,*) '%%               clm_pf_interface_init                 %%'
!      write(iulog,*) '%%-----------------------------------------------------%%'
!      write(iulog,*) ' '
!    endif
!
!
!    ! Determine necessary indices
!    call get_proc_global(numg, numl, numc, nump)
!
!    !------------------------------------------------------------------------
!    allocate(pflotran_m)
!
!    ! Create PFLOTRAN model
!    call pflotranModelCreate(mpicom, pflotran_prefix, pflotran_m)
!    option => pflotran_m%option
!
!    call pflotranModelSetupMappingfiles(pflotran_m)
!
!    ! initialize pflotran model
!    call pflotranModelStepperRunInit(pflotran_m)
!
!    ! PFLOTRAN ck file NOT works well when coupling with CLM.
!    ! So it's off now and restart PF from CLM.
!    !restart_stamp = ""
!    !call pflotranModelSetupRestart(pflotran_m, restart_stamp)
!    initth_pf2clm = .false.
!    !if (restart_stamp .ne. "") then
!    !   initth_pf2clm = .true.
!    !endif
!
!    ! Initialize PETSc vector for data transfer between CLM and PFLOTRAN
!    call CLMPFLOTRANIDataInit()
!
!    clm_pf_idata%nzclm_mapped  = nlevsoi        ! the soil layer no. mapped btw CLM and PF for data-passing
!
!    select type (simulation => pflotran_m%simulation)
!      class is (subsurface_simulation_type)
!         realization => simulation%realization
!      class default
!         pflotran_m%option%io_buffer = "This version of clm-pflotran only works with subsurface simulations."
!         write(*, '(/A/)') pflotran_m%option%io_buffer
!         call printErrMsg(pflotran_m%option)
!    end select
!
!    if (pflotran_m%option%iflowmode == TH_MODE .or. pflotran_m%option%iflowmode == Richards_MODE) then
!       pflotran_m%option%io_buffer = "This version of clm-pflotran Richards_mode or TH_mode has NOT yet well tested."
!       write(*, '(/A/)') pflotran_m%option%io_buffer
!       call printErrMsg(pflotran_m%option)
!    endif
!
!    if(pflotran_m%option%nsurfflowdof > 0) then
!       pflotran_m%option%io_buffer = "This version of clm-pflotran DOES NOT work with PF Surface simulation."
!       write(*, '(/A/)') pflotran_m%option%io_buffer
!       call printErrMsg(pflotran_m%option)
!    endif
!    pf_surfaceflow = .false.
!
!    !------------------------------------------------
!
!    ! Compute number of cells in CLM domain.
!    clm_cell_npts  = (bounds%endg - bounds%begg + 1)*clm_pf_idata%nzclm_mapped
!    clm_2dtop_npts = (bounds%endg - bounds%begg + 1)
!    clm_2dbot_npts = (bounds%endg - bounds%begg + 1)
!    allocate(clm_cell_ids_nindex( 1:clm_cell_npts))
!    allocate(clm_2dtop_cell_ids_nindex(1:clm_2dtop_npts))
!    allocate(clm_2dbot_cell_ids_nindex(1:clm_2dbot_npts))
!
!    ! Save cell IDs of CLM grid
!    cellcount = 0
!    do g = bounds%begg, bounds%endg
!       gcount = g - bounds%begg + 1
!
!       do j = 1, clm_pf_idata%nzclm_mapped
!          cellcount = cellcount + 1
!          !clm_cell_ids_nindex(cellcount) = (ldecomp%gdc2glo(g)-1)*nlevsoi + j - 1
!          ! maintain CLM each processor's grids to match up with PF (this appears saving a lot of data communications)
!          clm_cell_ids_nindex(cellcount) = (g-1)*clm_pf_idata%nzclm_mapped + j - 1
!       enddo
!
!       clm_2dtop_cell_ids_nindex(gcount) = (g-1)*clm_pf_idata%nzclm_mapped
!
!       clm_2dbot_cell_ids_nindex(gcount) = g*clm_pf_idata%nzclm_mapped-1
!    enddo
!
!    ! CLM: 3-D Subsurface domain (local and ghosted cells)
!    clm_pf_idata%nlclm_sub = clm_cell_npts
!    clm_pf_idata%ngclm_sub = clm_cell_npts
!
!    ! CLM: Surface/Bottom cells of subsurface domain (local and ghosted cells)
!    clm_pf_idata%nlclm_2dtop = (bounds%endg - bounds%begg + 1)
!    clm_pf_idata%ngclm_2dtop = (bounds%endg - bounds%begg + 1)
!
!    ! CLM: bottom face of subsurface domain
!    clm_pf_idata%nlclm_2dbot = (bounds%endg - bounds%begg + 1)
!    clm_pf_idata%ngclm_2dbot = (bounds%endg - bounds%begg + 1)
!
!    ! PFLOTRAN: 3-D Subsurface domain (local and ghosted cells)
!    clm_pf_idata%nlpf_sub = realization%patch%grid%nlmax
!    clm_pf_idata%ngpf_sub = realization%patch%grid%ngmax
!
!    ! PFLOTRAN: surface of subsurface domain (local and ghosted cells)
!    !if (pflotran_m%option%iflowmode == RICHARDS_MODE) then
!    !  clm_pf_idata%nlpf_2dtop = pflotranModelNSurfCells3DDomain(pflotran_m)   ! this will be overwritten in Mapping below
!    !  clm_pf_idata%ngpf_2dtop = pflotranModelNSurfCells3DDomain(pflotran_m)
!    !endif
!
!    ! For CLM/PF: ground surface NOT defined, so need to set the following to zero.
!    clm_pf_idata%nlclm_srf = 0
!    clm_pf_idata%ngclm_srf = 0
!    clm_pf_idata%nlpf_srf = 0
!    clm_pf_idata%ngpf_srf = 0
!
!    ! Initialize maps for transferring data between CLM and PFLOTRAN.
!    call pflotranModelInitMapping(pflotran_m, clm_cell_ids_nindex, &
!                                  clm_cell_npts, CLM_SUB_TO_PF_SUB)
!    call pflotranModelInitMapping(pflotran_m, clm_cell_ids_nindex, &
!                                  clm_cell_npts, PF_SUB_TO_CLM_SUB)
!
!    !
!    if (pflotran_m%option%iflowmode == RICHARDS_MODE .or. &
!        pflotran_m%option%iflowmode == TH_MODE) then
!       call pflotranModelInitMapping(pflotran_m, clm_2dtop_cell_ids_nindex,   &
!                                      clm_2dtop_npts, CLM_2DTOP_TO_PF_2DTOP)
!       call pflotranModelInitMapping(pflotran_m, clm_2dtop_cell_ids_nindex,   &
!                                      clm_2dtop_npts, PF_2DTOP_TO_CLM_2DTOP)
!
!       call pflotranModelInitMapping(pflotran_m, clm_2dbot_cell_ids_nindex, &
!                                     clm_2dbot_npts, CLM_2DBOT_TO_PF_2DBOT)
!       call pflotranModelInitMapping(pflotran_m, clm_2dbot_cell_ids_nindex, &
!                                     clm_2dbot_npts, PF_2DBOT_TO_CLM_2DBOT)
!    endif
!
!    ! the CLM-CN/BGC decomposing pool size
!    clm_pf_idata%ndecomp_pools     = ndecomp_pools
!
!    ! Allocate vectors for data transfer between CLM and PFLOTRAN.
!    call CLMPFLOTRANIDataCreateVec(MPI_COMM_WORLD)
!
!    ! the CLM-CN/BGC decomposing pool if varying cn ratios, and names.
!    clm_pf_idata%floating_cn_ratio = decomp_cascade_con%floating_cn_ratio_decomp_pools(1:ndecomp_pools)
!    clm_pf_idata%decomp_pool_name  = decomp_cascade_con%decomp_pool_name_history(1:ndecomp_pools)
!
!write(*,'(A40,10L10)')"DEBUG | floating_cn_ratio=",clm_pf_idata%floating_cn_ratio
!write(*,'(A40,10A10)')"DEBUG | decomp_pool_name=",clm_pf_idata%decomp_pool_name
!    ! if BGC is on
!    if(pflotran_m%option%ntrandof > 0) then
!      call pflotranModelGetRTspecies(pflotran_m)
!    endif
!
!    ! Get pflotran top surface area
!    call pflotranModelGetTopFaceArea(pflotran_m)
!
!    deallocate(clm_cell_ids_nindex)
!    deallocate(clm_2dtop_cell_ids_nindex)
!    deallocate(clm_2dbot_cell_ids_nindex)
!
!!-------------------------------------------------------------------------------------
!    ! coupled module controls betweeen PFLOTRAN and CLM45 (F.-M. Yuan, Aug. 2013)
!    if(pflotran_m%option%iflowmode==RICHARDS_MODE) then
!      pf_hmode = .true.
!      pf_tmode = .false.
!      pf_frzmode = .false.
!
!    elseif(pflotran_m%option%iflowmode==TH_MODE) then
!      pf_hmode = .true.
!      pf_tmode = .true.
!      if (pflotran_m%option%use_th_freezing) then
!        pf_frzmode = .true.
!      else
!        pf_frzmode = .false.
!      endif
!
!    endif
!
!    if(pflotran_m%option%ntrandof.gt.0) then
!      pf_cmode = .true.        ! initialized as '.false.' in clm initialization
!    endif
!
!!-------------------------------------------------------------------------------------
!
!    end associate
!  end subroutine interface_init
!
!  !-----------------------------------------------------------------------------
!  !
!  ! !SUBROUTINE: pflotran_run_onestep
!  !
!  ! !INTERFACE:
!  subroutine pflotran_run_onestep(bounds,                         &
!           num_soilc, filter_soilc,                               &
!           num_soilp, filter_soilp,                               &
!           ! soil thermal-hydrology (TODO: will update when testing with th coupling)
!           atm2lnd_vars,                                          &
!           waterstate_vars, waterflux_vars,                       &
!           soilstate_vars,  temperature_vars, energyflux_vars,    &
!           soilhydrology_vars, soil_water_retention_curve,        &
!           ! soil bgc
!           cnstate_vars, carbonflux_vars, carbonstate_vars,       &
!           nitrogenflux_vars, nitrogenstate_vars,                 &
!           ch4_vars)
!  !
!  ! !DESCRIPTION:
!  !
!  !  F.-M. YUAN: based on Gautam's 'step_th_clm_pf',
!  !              'chemistry' (PF_CMODE) added (Sept. 6, 2013)
!  !
!  ! !USES:
!
!
!    use clm_pflotran_interface_data
!    use PFLOTRAN_Constants_module
!
!    use clm_time_manager  , only : get_step_size, get_nstep, nsstep, nestep,       &
!                                   is_first_step, is_first_restart_step
!
!    !use clm_varctl           , only : use_century_decomp
!    !use CNDecompCascadeCNMod , only : decomp_rate_constants_cn
!    !use CNDecompCascadeBGCMod, only : decomp_rate_constants_bgc
!
!    implicit none
!
!    type(bounds_type) , intent(in)  :: bounds
!    integer, intent(in) :: num_soilc                  ! number of soil columns in filter
!    integer, intent(in) :: filter_soilc(:)            ! filter for soil columns
!    integer, intent(in) :: num_soilp                  ! number of soil pfts in filter
!    integer, intent(in) :: filter_soilp(:)            ! filter for soil pfts
!
!    type(atm2lnd_type)      , intent(in) :: atm2lnd_vars
!
!    type(soilstate_type)    , intent(in) :: soilstate_vars
!    type(soilhydrology_type), intent(in) :: soilhydrology_vars
!    class(soil_water_retention_curve_type), intent(in) :: soil_water_retention_curve
!
!    type(waterstate_type)   , intent(inout) :: waterstate_vars
!    type(waterflux_type)    , intent(inout) :: waterflux_vars
!    type(temperature_type)  , intent(inout) :: temperature_vars
!    type(energyflux_type)   , intent(inout) :: energyflux_vars
!
!    type(cnstate_type)      , intent(inout) :: cnstate_vars
!    type(carbonstate_type)  , intent(inout) :: carbonstate_vars
!    type(carbonflux_type)   , intent(inout) :: carbonflux_vars
!    type(nitrogenstate_type), intent(inout) :: nitrogenstate_vars
!    type(nitrogenflux_type) , intent(inout) :: nitrogenflux_vars
!
!    type(ch4_type)          , intent(inout) :: ch4_vars
!
!    !LOCAL VARIABLES:
!    real(r8) :: dtime                         ! land model time step (sec)
!    integer  :: nstep                         ! time step number
!    integer  :: total_clmstep                 ! total clm time step number
!    logical  :: ispfprint= .TRUE.             ! let PF printout or not
!    logical  :: isinitpf = .FALSE.            ! (re-)initialize PF from CLM or not
!
!#ifdef CLM_PF_DEBUG
!    real(r8) :: t0, t1, t2, t3
!#endif
!
!  !-----------------------------------------------------------------------
!
!    nstep = get_nstep() - nsstep
!    dtime = get_step_size()
!
!    if (is_first_step() .or. is_first_restart_step()) then
!       isinitpf = .TRUE.
!    else
!       isinitpf = .FALSE.
!    endif
!
!#ifdef CLM_PF_DEBUG
!if(nstep>=48*210 .and. nstep<=48*211) then
!  call cpu_time(t0)
!  if(pflotran_m%option%myrank .eq. pflotran_m%option%io_rank) then
!      write(pflotran_m%option%myrank+200,*) '------------------------------------------------------- '
!      write(pflotran_m%option%myrank+200,*) '------- checking CLM-PFLOTRAN timing - nstep = ', nstep
!      write(pflotran_m%option%myrank+200,*) 'CPU_time @check-point 0: ', t0
!  endif
!endif
!#endif
!
!write(*,'(/,A,50(1h-))')">>>DEBUG | pflotran_run_onestep:"
!write(*,'(10A20)')"pf_cmode", "pf_tmode", "pf_hmode", "pf_frzmode","isinitpf", "initth_pf2clm"
!write(*,'(10L20)')pf_cmode,pf_tmode, pf_hmode, pf_frzmode,isinitpf, initth_pf2clm
!
!    ! (0)
!    if (isinitpf) then
!       total_clmstep = nestep - nsstep
!       ispfprint = .true.               ! turn-on or shut-off PF's *.h5 output
!       call pflotranModelUpdateFinalWaypoint(pflotran_m, total_clmstep*dtime, ispfprint)
!
!       ! Set CLM soil properties onto PFLOTRAN grid
!       call get_clm_soil_properties(bounds,  &
!              num_soilc, filter_soilc,       &
!              soilstate_vars)
!       call pflotranModelSetSoilProp(pflotran_m)
!
!       ! wgs: beg
!       ! if initializing soil 'TH' states from CLM to pflotran
!       if (.not.initth_pf2clm) then
!          call get_clm_soil_th(initth_pf2clm, initth_pf2clm, &
!                    bounds, num_soilc, filter_soilc,         &
!                    atm2lnd_vars, soilstate_vars,            &
!                    waterstate_vars, temperature_vars,       &
!                    soil_water_retention_curve)
!
!          if (pf_hmode) then
!          ! directly pass TH to internal PF vec (field%, work%)
!              call pflotranModelSetInternalTHStatesfromCLM(pflotran_m)
!          else
!          ! pass TH to global_auxvar
!              call pflotranModelUpdateTHfromCLM(pflotran_m, pf_hmode, pf_tmode)
!          end if
!
!       ! if initializaing CLM's H states from pflotran (only H mode now)
!       else
!          call pflotranModelGetSaturationFromPF(pflotran_m)   ! hydrological states
!          call update_soil_moisture_pf2clm(                 &
!                    bounds, num_soilc, filter_soilc,        &
!                     soilstate_vars, waterstate_vars)
!       end if
!       ! wgs: end
!
!       ! the following is for some specific PF's hydrological parameters useful to constrain H source/sink or BC
!       if (pf_hmode) then
!          call pflotranModelGetSoilPropFromPF(pflotran_m)
!       end if
!
!    endif !!if (isinitpf)
!
!    ! (1) passing TH states from CLM to PF, if not H/TH mode NOT on in PF, every CLM time-step
!
!    ! if PF T/H mode not available, have to pass those from CLM to global variable in PF to drive BGC/H
!    if ((.not.pf_tmode .or. .not.pf_hmode) .and. (.not.isinitpf)) then
!        call get_clm_soil_th(pf_tmode, pf_hmode,  &
!           bounds, num_soilc, filter_soilc,         &
!           atm2lnd_vars, soilstate_vars,            &
!           waterstate_vars, temperature_vars,       &
!           soil_water_retention_curve)
!
!        call pflotranModelUpdateTHfromCLM(pflotran_m, pf_hmode, pf_tmode)
!    endif
!
!    ! ice-len adjusted porostiy
!    if (.not.pf_frzmode) then
!        call get_clm_iceadj_porosity(               &
!           bounds, num_soilc, filter_soilc,         &
!           soilstate_vars, waterstate_vars)
!
!        call pflotranModelResetSoilPorosityFromCLM(pflotran_m)
!!write(*,*)">>>DEBUG | pflotranModelResetSoilPorosityFromCLM"
!    endif
!
!    ! (2) pass CLM water fluxes to CLM-PFLOTRAN interface
!    if (pf_hmode) then      !if coupled 'H' mode between CLM45 and PFLOTRAN
!        call get_clm_bcwflx(                        &
!            bounds, num_soilc, filter_soilc,        &
!            atm2lnd_vars, soilstate_vars,           &
!            temperature_vars, energyflux_vars,      &
!            waterstate_vars, waterflux_vars)
!
!        ! pass flux 'vecs' from CLM to pflotran
!        call pflotranModelUpdateHSourceSink(pflotran_m)    ! H SrcSink
!        call pflotranModelSetSoilHbcsFromCLM(pflotran_m)   ! H BC
!    end if
!
!    ! (3) CLM thermal BC to PFLOTRAN
!    if (pf_tmode) then
!        call get_clm_bceflx(                        &
!            bounds, num_soilc, filter_soilc,        &
!            atm2lnd_vars, waterstate_vars,          &
!            temperature_vars, energyflux_vars)
!
!        call pflotranModelUpdateSubsurfTCond(pflotran_m)   ! SrcSink and T bc
!    end if
!
!    ! (4)
!    if (pf_cmode) then
!      ! (4a) for checking CLM's T/H response functions (TODO - also for passing decomposition rate to PF if needed)
!      !if (use_century_decomp) then
!      !    call decomp_rate_constants_bgc(bounds, num_soilc, filter_soilc)
!      !else
!      !    call decomp_rate_constants_cn(bounds, num_soilc, filter_soilc)
!      !end if
!
!      ! (4b) reset PFLOTRAN bgc state variables from CLM-CN, every CLM time-step
!!      if (isinitpf) then     ! NOTE: if only initialize ONCE, uncomment this 'if...endif' block
!        call get_clm_bgc_conc(bounds,    &
!           num_soilc, filter_soilc,      &
!           carbonstate_vars,             &
!           nitrogenstate_vars,           &
!           ch4_vars)
!
!        call pflotranModelSetBgcConcFromCLM(pflotran_m)
!!write(*,*)">>>DEBUG | pflotranModelSetBgcConcFromCLM"
!        if ( (.not.pf_hmode .or. .not.pf_frzmode)) then
!        ! this is needed, because at step 0, PF's interface data is empty
!        !which causes Aq. conc. adjustment balance issue
!           call pflotranModelGetSaturationFromPF(pflotran_m)
!!write(*,*)">>>DEBUG | pflotranModelGetSaturationFromPF"
!        endif
!!      endif     ! NOTE: if only initialize ONCE, uncomment this 'if...endif' block
!
!      ! MUST reset PFLOTRAN soil aq. bgc state variables due to liq. water volume change
!      ! when NOT coupled with PF Hydrology or NOT in freezing-mode (porosity will be forced to vary from CLM)
!        if (.not.pf_hmode .or. .not.pf_frzmode) then
!           call pflotranModelUpdateAqConcFromCLM(pflotran_m)
!!write(*,*)">>>DEBUG | pflotranModelUpdateAqConcFromCLM"
!        endif
!
!      ! (4c) bgc rate (source/sink) from CLM to PFLOTRAN
!        call get_clm_bgc_rate(bounds,   &
!            num_soilc, filter_soilc,    &
!            cnstate_vars,               &
!            carbonflux_vars,            &
!            nitrogenflux_vars)
!
!        call pflotranModelSetBgcRatesFromCLM(pflotran_m)
!!write(*,*)">>>DEBUG | pflotranModelSetBgcRatesFromCLM"
!    endif
!
!#ifdef CLM_PF_DEBUG
!if(nstep>=48*210 .and. nstep<=48*211) then
!  call cpu_time(t1)
!  if(pflotran_m%option%myrank .eq. pflotran_m%option%io_rank) then
!      write(pflotran_m%option%myrank+200,*) 'CPU_time elapsed @check-point 1 - 0: ', t1-t0
!  endif
!endif
!#endif
!
!write(*,*)">>>DEBUG | pflotranModelStepperRunTillPauseTime: BEG"
!    ! (5) the main callings of PFLOTRAN
!    call pflotranModelStepperRunTillPauseTime( pflotran_m, (nstep+1.0d0)*dtime, dtime, .false. )
!write(*,*)">>>DEBUG | pflotranModelStepperRunTillPauseTime: END"
!
!#ifdef CLM_PF_DEBUG
!if(nstep>=48*210 .and. nstep<=48*211) then
!  call cpu_time(t2)
!  if(pflotran_m%option%myrank .eq. pflotran_m%option%io_rank) then
!      write(pflotran_m%option%myrank+200,*) 'CPU_time elapsed @check-point 2 - 1: ', t2-t1
!  endif
!endif
!#endif
!
!    ! (6) update CLM variables from PFLOTRAN
!    if (pf_hmode) then
!        call pflotranModelGetSaturationFromPF(pflotran_m)   ! hydrological states
!
!        call update_soil_moisture_pf2clm(        &
!           bounds, num_soilc, filter_soilc,      &
!           soilstate_vars, waterstate_vars)
!    endif
!
!    if (pf_tmode) then
!        call pflotranModelGetTemperatureFromPF(pflotran_m)  ! thermal states
!
!        call update_soil_temperature_pf2clm(     &
!           bounds, num_soilc, filter_soilc,      &
!           soilstate_vars, temperature_vars)
!    endif
!
!    ! bgc variables
!    if (pf_cmode) then
!        call pflotranModelGetBgcVariablesFromPF(pflotran_m)
!!write(*,*)">>>DEBUG | pflotranModelGetBgcVariablesFromPF"
!        call update_soil_bgc_pf2clm(bounds,      &
!           num_soilc, filter_soilc,              &
!           atm2lnd_vars, waterstate_vars,        &
!           soilstate_vars, cnstate_vars,         &
!           carbonstate_vars, carbonflux_vars,    &
!           nitrogenstate_vars,nitrogenflux_vars, &
!           ch4_vars                              &
!           )
!!write(*,*)">>>DEBUG | update_soil_bgc_pf2clm"
!        ! need to save the current time-step PF porosity/liq. saturation for bgc species mass conservation
!        ! if CLM forced changing them into PF at NEXT timestep
!        if (.not.pf_hmode .or. .not.pf_frzmode) then
!           call pflotranModelGetSaturationFromPF(pflotran_m)
!        endif
!!write(*,*)">>>DEBUG | pflotranModelGetSaturationFromPF"
!    endif
!
!    ! the actual infiltration/runoff/drainage and solute flux with BC, if defined,
!    ! are retrieving from PFLOTRAN using 'update_bcflow_pf2clm' subroutine
!    if (pf_hmode) then
!        call pflotranModelGetBCMassBalanceDeltaFromPF(pflotran_m )
!
!        call update_bcflow_pf2clm(              &
!           bounds, num_soilc, filter_soilc,     &
!           atm2lnd_vars,                        &
!           temperature_vars, energyflux_vars,   &
!           waterstate_vars, waterflux_vars)
!
!        !(TODO) bgc BC flows (e.g. runoff/leaching)
!    endif
!
!
!#ifdef CLM_PF_DEBUG
!if(nstep>=48*210 .and. nstep<=48*211) then
!  call cpu_time(t3)
!  if(pflotran_m%option%myrank .eq. pflotran_m%option%io_rank) then
!      write(pflotran_m%option%myrank+200,*) 'CPU_time elapsed @check-point 3 - 2: ', t3-t2
!      write(pflotran_m%option%myrank+200,*) 'CPU_time elapsed @check-point 3 - 0: ', t3-t0
!      write(pflotran_m%option%myrank+200,*) '------------------------------------------------------- '
!  endif
!endif
!#endif
!
!  end subroutine pflotran_run_onestep
!
!  !-----------------------------------------------------------------------
!  !BOP
!  !
!  ! !ROUTINE: write_checkpoint
!  !
!  ! !INTERFACE:
!  subroutine pflotran_write_checkpoint(date_stamp)
!  !
!  ! !DESCRIPTION:
!  ! Trigger a pflotran checkpoint file to be written
!  !
!  ! !USES:
!  ! !ARGUMENTS:
!    character(len=32), intent(in) :: date_stamp ! file name date stamp
!
!  ! !LOCAL VARIABLES:
!
!  !EOP
!  !-----------------------------------------------------------------------
!
!    ! temporarily OFF - it's not working well for BGC
!    ! So, now must initializing PF variables from CLM each start/restart.
!
!    !call pflotranModelStepperCheckpoint(pflotran_m, date_stamp)
!
!  end subroutine pflotran_write_checkpoint
!
!  !-----------------------------------------------------------------------------
!  !
!  ! !IROUTINE: pflotran_finalize
!  !
!  ! !INTERFACE:
!  subroutine pflotran_finalize()
!  !
!  ! !DESCRIPTION:
!  !
!  ! finalizing pflotran runs and destroying objects
!  !
!  ! !USES:
!
!    implicit none
!
!  !-----------------------------------------------------------------------
!
!    if (use_pflotran) then
!       call pflotranModelDestroy(pflotran_m)
!    endif
!
!  end subroutine pflotran_finalize


  ! ============================= GET CLM initial/src-sink/BC to PFLOTRAN ==================================

  !-----------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: get_clm_soil_properties
  !
  ! !INTERFACE:
  subroutine get_clm_soil_property(clm_bgc_data,            &
                        bounds, num_soilc, filter_soilc,    &
                        soilstate_vars)
    !
    ! !DESCRIPTION:
    ! get soil column physical properties to PFLOTRAN
    !
    ! !USES:
    use CNDecompCascadeConType, only : decomp_cascade_con
    ! pflotran
    !
    ! !ARGUMENTS:

    implicit none

    !
    ! !REVISION HISTORY:
    ! Created by Gautam Bisht
    ! Revised by Fengming Yuan, CCSI-ORNL
    !
    !EOP
    !
    ! LOCAL VARAIBLES:

    type(bounds_type)        , intent(in) :: bounds           ! bounds
    integer                  , intent(in) :: num_soilc        ! number of column soil points in column filter
    integer                  , intent(in) :: filter_soilc(:)  ! column filter for soil points
    type(soilstate_type)     , intent(in) :: soilstate_vars
    type(clm_bgc_interface_data_type), intent(inout) :: clm_bgc_data

    integer  :: fc, g, l, c, j      ! indices
    integer  :: gcount, cellcount

    character(len= 32) :: subname = 'get_clm_soil_property' ! subroutine name

    associate( &
         ! Assign local pointer to derived subtypes components (column-level)
         clandunit  =>  col%landunit   , & !  [integer (:)]  landunit index of column
         cgridcell  =>  col%gridcell   , & !  [integer (:)]  gridcell index of column
         wtgcell    =>  col%wtgcell    , & !  [real(r8) (:)]  weight (relative to gridcell
         cactive    =>  col%active     , & !  [logical (:)]  column active or not
         z          =>  col%z          , & !  [real(r8) (:,:)]  layer depth (m)
         dz         =>  col%dz         , & !  [real(r8) (:,:)]  layer thickness depth (m)
         zi         =>  col%zi         , & !  [real(r8) (:,:)]  interface level below a "z" level (m)
         !
         bd         =>  soilstate_vars%bd_col           , & !
         bsw        =>  soilstate_vars%bsw_col          , & !  [real(r8) (:,:)]  Clapp and Hornberger "b" (nlevgrnd)
         hksat      =>  soilstate_vars%hksat_col        , & !  [real(r8) (:,:)]  hydraulic conductivity at saturation (mm H2O /s) (nlevgrnd)
         sucsat     =>  soilstate_vars%sucsat_col       , & !  [real(r8) (:,:)]  minimum soil suction (mm) (nlevgrnd)
         watsat     =>  soilstate_vars%watsat_col       , & !  [real(r8) (:,:)]  volumetric soil water at saturation (porosity) (nlevgrnd)
         watfc      =>  soilstate_vars%watfc_col        , & !  [real(r8) (:,:)]  volumetric soil water at saturation (porosity) (nlevgrnd)

         porosity   =>  soilstate_vars%porosity_col     , &
         eff_porosity=> soilstate_vars%eff_porosity_col , &

         initial_cn_ratio => decomp_cascade_con%initial_cn_ratio , &
         initial_cp_ratio => decomp_cascade_con%initial_cp_ratio   &
         )

!-------------------------------------------------------------------------------------

!
    do fc = 1, num_soilc
        c = filter_soilc(fc)
        do j = 1,nlevsoi
            clm_bgc_data%z(c,j)             = z(c,j)
            clm_bgc_data%dz(c,j)            = dz(c,j)
            clm_bgc_data%bd_col(c,j)        = bd(c,j)
            clm_bgc_data%bsw_col(c,j)       = bsw(c,j)
            clm_bgc_data%hksat_col(c,j)     = hksat(c,j)
            clm_bgc_data%sucsat_col(c,j)    = sucsat(c,j)
            clm_bgc_data%watsat_col(c,j)    = watsat(c,j)
            clm_bgc_data%watfc_col(c,j)     = watfc(c,j)

            clm_bgc_data%porosity_col(c,j)      = porosity(c,j)
            clm_bgc_data%eff_porosity_col(c,j)  = eff_porosity(c,j)
        end do
    end do

    clm_bgc_data%initial_cn_ratio(:)= initial_cn_ratio(:)
    clm_bgc_data%initial_cp_ratio(:)= initial_cp_ratio(:)

!!write(*,'(A40,10E14.6)')">>>DEBUG | hksat_x=",(hksat_x_clm_loc(1:10))
!!write(*,'(A40,10E14.6)')">>>DEBUG | hksat_y=",(hksat_y_clm_loc(1:10))
!!write(*,'(A40,10E14.6)')">>>DEBUG | hksat_z=",(hksat_z_clm_loc(1:10))
!!write(*,'(A40,10E14.6)')">>>DEBUG | sucsat=",(sucsat_clm_loc(1:10))
!!write(*,'(A40,10E14.6)')">>>DEBUG | watsat=",(watsat_clm_loc(1:10))
!!write(*,'(A40,10E14.6)')">>>DEBUG | watfc=",(watfc_clm_loc(1:10))
!!write(*,'(A40,10E14.6)')">>>DEBUG | bulkdensity=",(bulkdensity_dry_clm_loc(1:10))

    end associate
  end subroutine get_clm_soil_property


  !-----------------------------------------------------------------------------
  !BOP
  !
  ! !ROUTINE: get_clm_soil_th
  !
  ! !INTERFACE:
  subroutine get_clm_soil_thermohydro(clm_bgc_data,           &
                       bounds, num_soilc, filter_soilc,     &
                       atm2lnd_vars, soilstate_vars,        &
                       waterstate_vars, waterflux_vars,     &
                       temperature_vars, energyflux_vars,   &
                       soil_water_retention_curve)
  !
  ! !DESCRIPTION:
  !  update soil temperature/saturation from CLM to PFLOTRAN for driving PF's BGC
  !  if either NOT available inside PFLOTRAN
  !
  ! !USES:
    use clm_time_manager    , only : get_nstep
    use shr_const_mod       , only : SHR_CONST_G

!    use PFLOTRAN_Constants_module

  ! !ARGUMENTS:
    implicit none

!    logical                  , intent(in) :: pftmode, pfhmode
    type(bounds_type)        , intent(in) :: bounds           ! bounds
    integer                  , intent(in) :: num_soilc        ! number of column soil points in column filter
    integer                  , intent(in) :: filter_soilc(:)  ! column filter for soil points
    type(atm2lnd_type)       , intent(in) :: atm2lnd_vars
    type(soilstate_type)     , intent(in) :: soilstate_vars
    type(waterstate_type)    , intent(in) :: waterstate_vars
    type(waterflux_type)     , intent(in) :: waterflux_vars
    type(temperature_type)   , intent(in) :: temperature_vars
    type(energyflux_type)    , intent(in) :: energyflux_vars
    class(soil_water_retention_curve_type), intent(in) :: soil_water_retention_curve
    type(clm_bgc_interface_data_type), intent(inout) :: clm_bgc_data

  ! !LOCAL VARIABLES:
    integer  :: fc, c, g, j, gcount, cellcount      ! indices
    integer  :: pftindex, p
    real(r8) :: sattmp, psitmp, itheta
    real(r8) :: watmin(num_soilc, nlevsoi)
    real(r8) :: sucmin(num_soilc, nlevsoi)

  !EOP
  !-----------------------------------------------------------------------
    associate ( &
      gridcell        => col%gridcell     , & ! column's gridcell
      wtgcell         => col%wtgcell      , & ! column's weight relative to gridcell
      cactive         => col%active       , & ! [logical (:)]  column active or not
      dz              => col%dz           , & ! layer thickness depth (m)
      zi              => col%zi           , & ! interface depth (m)
    !
!      sucsat          => soilstate_vars%sucsat_col    , & ! minimum soil suction (mm) (nlevgrnd)
!      bsw             => soilstate_vars%bsw_col       , & ! Clapp and Hornberger "b"
!      watsat          => soilstate_vars%watsat_col    , & ! volumetric soil water at saturation (porosity) (nlevgrnd)
      soilpsi         => soilstate_vars%soilpsi_col   , & ! soil water matric potential in each soil layer (MPa)
      rootfr          => soilstate_vars%rootfr_col             , & ! pft-level effective fraction of roots in each soil layer
    !
      h2osoi_liq            => waterstate_vars%h2osoi_liq_col   , & ! liquid water (kg/m2)
      h2osoi_ice            => waterstate_vars%h2osoi_ice_col   , & ! ice lens (kg/m2)
      frac_sno              => waterstate_vars%frac_sno_eff_col            , & ! Input: fraction of ground covered by snow (0 to 1)
      frac_h2osfc           => waterstate_vars%frac_h2osfc_col             , & ! Input: fraction of ground covered by surface water (0 to 1)
    !
      t_soisno              => temperature_vars%t_soisno_col    , & ! snow-soil temperature (Kelvin)
      t_grnd                => temperature_vars%t_grnd_col                  , & ! Input:  [real(r8) (:)]  ground surface temperature [K]
    !
      forc_pbot             => atm2lnd_vars%forc_pbot_not_downscaled_grc  , &  ! atmospheric pressure (Pa)
      forc_pco2             => atm2lnd_vars%forc_pco2_grc                  , & ! partial pressure co2 (Pa)
      forc_pch4             => atm2lnd_vars%forc_pch4_grc                  , & ! partial pressure ch4 (Pa)

      htvp                  => energyflux_vars%htvp_col                    , & ! Input:  [real(r8) (:)]  latent heat of vapor of water (or sublimation) [j/kg]
      eflx_bot              => energyflux_vars%eflx_bot_col                   , &! heat flux from beneath column (W/m**2) [+ = upward]
      eflx_gnet_patch       => energyflux_vars%eflx_gnet_patch                , &! net ground heat flux into the surface (W/m**2) per patch
      eflx_soil_grnd_patch  => energyflux_vars%eflx_soil_grnd_patch           , &! soil heat flux (W/m**2) [+ = into soil]

      qflx_top_soil         => waterflux_vars%qflx_top_soil_col           , & ! Input: net water input into soil from top (mm/s)
      qflx_ev_h2osfc        => waterflux_vars%qflx_ev_h2osfc_col          , & ! Input: column-level evaporation flux from h2osfc (W/m2) [+ to atm] : checking unit
      qflx_evap_soi         => waterflux_vars%qflx_evap_soi_col           , & ! Input: column-level soil evaporation (mm H2O/s) (+ = to atm)
      qflx_sub_snow         => waterflux_vars%qflx_sub_snow_col           , & ! Input: column-level evaporation flux from snow (mm H2O/s) [+ to atm]
      qflx_tran_veg         => waterflux_vars%qflx_tran_veg_col         &   ! Input: pft-level vegetation transpiration (mm H2O/s) (+ = to atm)

    )

    !--------------------------------------------------------------------------------------
!
!    watmin(:,:) = 0.01_r8
!    sucmin(:,:) = 1.e8_r8

    !! grid:
    clm_bgc_data%forc_pbot_not_downscaled_grc(:)    = forc_pbot(:)
    clm_bgc_data%forc_pco2_grc(:)                   = forc_pco2(:)
    clm_bgc_data%forc_pch4_grc(:)                   = forc_pch4(:)


    do fc = 1,num_soilc
        c = filter_soilc(fc)

        clm_bgc_data%frac_sno_eff_col(c)  = frac_sno(c)
        clm_bgc_data%frac_h2osfc_col(c)   = frac_h2osfc(c)

        clm_bgc_data%t_grnd_col(c)        = t_grnd(c)

        clm_bgc_data%qflx_top_soil_col(c) = qflx_top_soil(c)
        clm_bgc_data%qflx_ev_h2osfc_col(c)= qflx_ev_h2osfc(c)
        clm_bgc_data%qflx_evap_soi_col(c) = qflx_evap_soi(c)
        clm_bgc_data%qflx_sub_snow_col(c) = qflx_sub_snow(c)
        clm_bgc_data%qflx_tran_veg_col(c) = qflx_tran_veg(c)

        clm_bgc_data%htvp_col(c)    = htvp(c)
        clm_bgc_data%eflx_bot_col(c)    = eflx_bot(c)

        do j = 1, nlevsoi
            clm_bgc_data%soilpsi_col(c,j)       = soilpsi(c,j)
            clm_bgc_data%rootfr_col(c,j)        = rootfr(c,j)

            clm_bgc_data%h2osoi_liq_col(c,j)    = h2osoi_liq(c,j)
            clm_bgc_data%h2osoi_ice_col(c,j)    = h2osoi_ice(c,j)

            clm_bgc_data%t_soisno_col(c,j)      = t_soisno(c,j)

        end do
    end do

    ! CLM appears NO column-level ground-heat-flux variable, instead by 'patch'
    do fc = 1, num_soilc
        c = filter_soilc(fc)
        clm_bgc_data%eflx_soil_grnd_col(c)  = 0._r8
        clm_bgc_data%eflx_gnet_col(c)       = 0._r8
        do pftindex = 1, max_patch_per_col
            if (pftindex <= col%npfts(c)) then
                p = col%pfti(c) + pftindex - 1
                clm_bgc_data%eflx_soil_grnd_col(c)  = clm_bgc_data%eflx_soil_grnd_col(c) &
                                                    + eflx_soil_grnd_patch(p) * pft%wtcol(p)           ! W/m2
                clm_bgc_data%eflx_gnet_col(c)       = clm_bgc_data%eflx_gnet_col(c) &
                                                    + eflx_gnet_patch(p) * pft%wtcol(p)
            end if
       end do
    end do
!
!write(*,'(A30,12E14.6)')">>>DEBUG | soillsat=", soillsat_clmp_loc(1:10)
!write(*,'(A30,12E14.6)')">>>DEBUG | gsoilpsi[Pa]=", soilpsi_clmp_loc(1:10)
!write(*,'(A30,12E14.6)')">>>DEBUG | soilt[oC]=", soilt_clmp_loc(1:10)


   end associate
  end subroutine get_clm_soil_thermohydro

  !-----------------------------------------------------------------------------
  !BOP
  !
  ! !ROUTINE: get_clm_iceadj_porosity
  !
  ! !INTERFACE:
!  subroutine get_clm_iceadj_porosity(bounds, &
!           num_soilc, filter_soilc,          &
!           soilstate_vars, waterstate_vars)
!  !
!  ! !DESCRIPTION:
!  !  update soil effective porosity from CLM to PFLOTRAN if PF freezing mode is off
!  !
!  ! !USES:
!
!    use PFLOTRAN_Constants_module
!    use clm_varctl          , only : pf_frzmode
!
!  ! !ARGUMENTS:
!    implicit none
!
!#include "finclude/petscsys.h"
!#include "finclude/petscvec.h"
!#include "finclude/petscvec.h90"
!
!    type(bounds_type)        , intent(in) :: bounds           ! bounds
!    integer                  , intent(in) :: num_soilc        ! number of column soil points in column filter
!    integer                  , intent(in) :: filter_soilc(:)  ! column filter for soil points
!    type(soilstate_type)     , intent(in) :: soilstate_vars
!    type(waterstate_type)    , intent(in) :: waterstate_vars
!
!  ! !LOCAL VARIABLES:
!    integer  :: fc, c, g, j, gcount, cellcount       ! indices
!    real(r8) :: itheta
!
!    PetscScalar, pointer :: adjporosity_clmp_loc(:)  !
!    PetscScalar, pointer :: soilisat_clmp_loc(:)  !
!    PetscErrorCode :: ierr
!
!  !EOP
!  !-----------------------------------------------------------------------
!    associate ( &
!    gridcell        => col%gridcell     , & ! column's gridcell
!    wtgcell         => col%wtgcell      , & ! column's weight relative to gridcell
!    cactive         => col%active       , & ! column's active or not
!    dz              => col%dz           , & ! layer thickness depth (m)
!    watsat          => soilstate_vars%watsat_col       , & ! volumetric soil water at saturation (porosity) (nlevgrnd)
!    h2osoi_ice      => waterstate_vars%h2osoi_ice_col    & ! ice lens (kg/m2)
!    )
!
!    ! if 'pf_tmode' is NOT using freezing option, the phase-change of soil water done in 'SoilTemperatureMod.F90' in 'bgp2'
!    ! must be included to adjust porosity (effective porosity) in pflotran
!    ! This is doing prior to the real liquid water source/sink, because 'h2osoi_liq' will be updated during those calls after 'bgp2'.
!    if (.not. pf_frzmode) then
!
!        ! re-calculate the effective porosity (CLM ice-len adjusted), which should be pass to pflotran
!        call VecGetArrayF90(clm_pf_idata%effporosity_clmp, adjporosity_clmp_loc,  ierr)
!        CHKERRQ(ierr)
!        call VecGetArrayF90(clm_pf_idata%soilisat_clmp, soilisat_clmp_loc,  ierr)
!        CHKERRQ(ierr)
!
!        adjporosity_clmp_loc(:) = 0._r8
!        soilisat_clmp_loc(:) = 0._r8
!        do fc = 1,num_soilc
!           c = filter_soilc(fc)
!           if ( wtgcell(c) <= 0._r8 .or. (.not.cactive(c)) ) cycle     ! don't assign data from PF for inactive cell
!
!           g = gridcell(c)
!           gcount = g - bounds%begg
!           do j = 1, nlevsoi
!             cellcount = gcount*clm_pf_idata%nzclm_mapped + j
!
!             if (j<=clm_pf_idata%nzclm_mapped) then
!               itheta = h2osoi_ice(c,j) / (dz(c,j) * denice)
!               itheta = min(itheta, 0.99_r8*watsat(c,j))
!               adjporosity_clmp_loc(cellcount) = adjporosity_clmp_loc(cellcount) &
!                 + (watsat(c,j) - itheta) * wtgcell(c)
!               soilisat_clmp_loc(cellcount)    = soilisat_clmp_loc(cellcount) &
!                 + itheta * wtgcell(c)
!             endif
!           end do
!        end do
!        call VecRestoreArrayF90(clm_pf_idata%effporosity_clmp,  adjporosity_clmp_loc,  ierr)
!        CHKERRQ(ierr)
!        call VecRestoreArrayF90(clm_pf_idata%soilisat_clmp, soilisat_clmp_loc,  ierr)
!        CHKERRQ(ierr)
!
!    end if
!
!    end associate
!  end subroutine get_clm_iceadj_porosity


  !-----------------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: get_clm_bcwflx
  !
  ! !INTERFACE:
!  subroutine get_clm_bcwflx(bounds,          &
!       num_soilc, filter_soilc,              &
!       clm_a2l, soils_vars,                  &
!       ces_vars, cef_vars, cws_vars, cwf_vars)
!  !
!  ! !DESCRIPTION:
!  !
!  !  F.-M. YUAN: the water fluxes in CLM4.5 are separately calculated in a few subroutines
!  !        in 'SoilHydrologyMod.F90'. When coupled with pflotran, it's hard to get those together
!  !        like GB does in 'step_th_clm_pf' subroutine. So, this subroutine is a collective call of
!  !        that and others in 'Hydrology2Mod.F90' so that pflotran can be called out of 'hydrology2'.
!  !
!  ! !USES:
!
!    use shr_const_mod   , only : SHR_CONST_G
!    use clm_time_manager, only : get_step_size, get_nstep
!
!  ! !ARGUMENTS:
!    implicit none
!
!#include "finclude/petscsys.h"
!#include "finclude/petscvec.h"
!#include "finclude/petscvec.h90"
!
!    type(bounds_type) , intent(in)  :: bounds
!    integer, intent(in) :: num_soilc                 ! number of column non-lake points in column filter
!    integer, intent(in) :: filter_soilc(:)           ! column filter for non-lake points
!
!    type(atm2lnd_type)       , intent(in) :: clm_a2l
!    type(soilstate_type)     , intent(in) :: soils_vars
!    type(temperature_type)   , intent(in) :: ces_vars
!    type(energyflux_type)    , intent(in) :: cef_vars
!    type(waterstate_type)    , intent(in) :: cws_vars
!    type(waterflux_type)     , intent(in) :: cwf_vars
!
!  ! !LOCAL VARIABLES:
!    integer  :: fc, g, c, j, p             ! do loop indices
!    integer  :: gcount                     ! gridcell index (0-based)
!    integer  :: pftindex                   ! pft index
!    real(r8) :: dtime                      ! land model time step (sec)
!    integer  :: nstep                      ! time step number
!    real(r8) :: area
!    integer  :: cellcount                  ! 3-D cell index (1-based)
!
!    real(r8) :: rootfr_col(bounds%begc:bounds%endc, nlevsoi)
!    real(r8) :: rootfr_sum(bounds%begc:bounds%endc)           ! accumulator for rootr weighting
!    real(r8) :: qflx_evap_col(bounds%begc:bounds%endc)        ! soil surface evaporation (mmH2O/s)
!    real(r8) :: qflx_tran_col(bounds%begc:bounds%endc)        ! veg. transpiration (mmH2O/s)
!
!    real(r8) :: qflx, qflx_sink, qflx_source, soilvwc
!    real(r8) :: dsoilliq1 = 0._r8, dsoilliq2 = 0._r8, dsoilliq3 = 0._r8
!
!    real(r8) :: qflx_ground, kbot
!    real(r8) :: reference_pressure, ponding_pressure      ! Pa
!    real(r8) :: pondmax(bounds%begc:bounds%endc)          ! mm H2O: max. ponding depth for column
!    real(r8) :: sr      = 0.10_r8
!    real(r8) :: tempreal
!
!    ! for PF --> CLM (seq.)
!    PetscScalar, pointer :: press_clms_loc(:)       !
!    PetscScalar, pointer :: soillsat_clms_loc(:)    !
!    PetscScalar, pointer :: porosity_clms_loc(:)    !
!    PetscScalar, pointer :: sr_pcwmax_clms_loc(:)   !
!
!    PetscScalar, pointer :: area_clms_loc(:)         !
!
!    ! for CLM (mpi) --> PF
!    PetscScalar, pointer :: zsoi_clmp_loc(:)         !
!    PetscScalar, pointer :: qflx_clmp_loc(:)         !   source/sink term for plant Transpiration: unit in mass rate (kgH2O/sec)
!    PetscScalar, pointer :: press_top_clmp_loc(:)    !   BC in pressure type: unit in Pa
!    PetscScalar, pointer :: press_base_clmp_loc(:)   !
!    PetscScalar, pointer :: qflux_top_clmp_loc(:)    !   BC in neumann flux type: unit in m/s
!    PetscScalar, pointer :: qflux_base_clmp_loc(:)   !
!    PetscScalar, pointer :: press_maxponding_clmp_loc(:)   !
!    PetscErrorCode :: ierr
!
!  !EOP
!  !-----------------------------------------------------------------------
!    associate ( &
!    ltype             => lun%itype             , & ! landunit type
!    cgridcell         => col%gridcell          , & ! column's gridcell
!    clandunit         => col%landunit          , & ! column's landunit
!    zi                => col%zi                , & ! Input: (:,:) soil layer interface depth (m)
!    dz                => col%dz                , & ! Input: (:,:) soil layer thickness (m)
!    pfti              => col%pfti                                , &! beginning pft index for each column
!    pwtgcell          => pft%wtgcell                             , &! weight relative to gridcell for each pft
!    pwtcol            => pft%wtcol                               , &! weight relative to column for each pft
!    !
!    bsw               => soils_vars%bsw_col                  , &! Clapp and Hornberger "b" (nlevgrnd)
!    hksat             => soils_vars%hksat_col                , &! hydraulic conductivity at saturation (mm H2O /s) (nlevgrnd)
!    watsat            => soils_vars%watsat_col               , &! volumetric soil water at saturation (porosity) (nlevgrnd)
!    sucsat            => soils_vars%sucsat_col               , &! minimum soil suction (mm) (nlevgrnd)
!    rootfr_pft        => soils_vars%rootfr_patch             , & ! pft-level effective fraction of roots in each soil layer
!
!    forc_pbot         => clm_a2l%forc_pbot_not_downscaled_grc , & ! Input:  [real(r8) (:)]  atmospheric pressure (Pa)
!    t_grnd            => ces_vars%t_grnd_col                  , & ! Input:  [real(r8) (:)]  ground surface temperature [K]
!    htvp              => cef_vars%htvp_col                    , & ! Input:  [real(r8) (:)]  latent heat of vapor of water (or sublimation) [j/kg]
!    !
!    frac_sno          => cws_vars%frac_sno_eff_col            , & ! Input: fraction of ground covered by snow (0 to 1)
!    frac_h2osfc       => cws_vars%frac_h2osfc_col             , & ! Input: fraction of ground covered by surface water (0 to 1)
!    h2osoi_liq        => cws_vars%h2osoi_liq_col              , & ! Input: liquid water (kg/m2)
!    h2osoi_ice        => cws_vars%h2osoi_ice_col              , & ! Input: ice lens (kg/m2)
!    !
!    qflx_top_soil     => cwf_vars%qflx_top_soil_col           , & ! Input: net water input into soil from top (mm/s)
!    qflx_ev_h2osfc    => cwf_vars%qflx_ev_h2osfc_col          , & ! Input: column-level evaporation flux from h2osfc (W/m2) [+ to atm] : checking unit
!    qflx_evap_soil    => cwf_vars%qflx_evap_soi_col           , & ! Input: column-level soil evaporation (mm H2O/s) (+ = to atm)
!    qflx_subl_snow    => cwf_vars%qflx_sub_snow_col           , & ! Input: column-level evaporation flux from snow (mm H2O/s) [+ to atm]
!    qflx_tran_veg_pft => cwf_vars%qflx_tran_veg_patch         &   ! Input: pft-level vegetation transpiration (mm H2O/s) (+ = to atm)
!
!    )
!
!!----------------------------------------------------------------------------
!    nstep = get_nstep()
!    dtime = get_step_size()
!
!    ! (1) soil surface evaporation: needs further checking here? -
!    do fc = 1, num_soilc
!       c = filter_soilc(fc)
!       if (ltype(clandunit(c)) == istsoil .or. ltype(clandunit(c))==istcrop) then
!           ! not sure if using 'qflx_evap_soi' as a whole ground-surface better than individual surfaces, i.e. 'qflx_ev_snow/soi/h2osfc'
!           ! all of those 4 variables are calculated in 'BareGroundFluxesMod', 'CanopyFluxesMod', and then adjusted in 'Biogeophysics2' after soil temperature call
!           ! note that: all 4 variables could be negative (i.e., dew formation on ground)
!           qflx_evap_col(c)=(1.0_r8 - frac_sno(c) - frac_h2osfc(c))*qflx_evap_soil(c) + &
!                          frac_h2osfc(c)*qflx_ev_h2osfc(c)/htvp(c)
!                          !frac_sno(c)*qflx_ev_snow(c)   ! snow-covered area should be excluded (see SoilHydrologyMod:: infiltration)
!       else
!          ! for other types of landunits
!          qflx_evap_col(c) = (1.0_r8 - frac_sno(c))*qflx_evap_soil(c)
!       end if
!
!       if (t_grnd(c) <= tfrz) qflx_evap_col(c) = max(0._r8, qflx_evap_col(c))   ! frozen ground, no dew contribution to subsurface infiltration
!    end do
!
!    ! (3) Compute the vegetation Transpiration (originally those are in 'SoilHydrologyMod.F90')
!    do j = 1, nlevsoi
!        do fc = 1, num_soilc
!            c = filter_soilc(fc)
!            rootfr_col(c,j) = 0._r8
!        end do
!    end do
!    rootfr_sum(:)    = 0._r8
!    qflx_tran_col(:) = 0._r8
!
!    do pftindex = 1, max_patch_per_col
!       do fc = 1, num_soilc
!          c = filter_soilc(fc)
!          if (pftindex <= col%npfts(c)) then
!             p = pfti(c) + pftindex - 1
!
!             if (pwtgcell(p)>0._r8) then
!                do j = 1,nlevsoi
!                       rootfr_col(c,j) = rootfr_col(c,j) + &
!                            rootfr_pft(p,j) * qflx_tran_veg_pft(p) * pwtcol(p)
!                       rootfr_sum(c) = rootfr_sum(c) + qflx_tran_veg_pft(p) * pwtcol(p)
!                end do
!             end if
!
!             qflx_tran_col(c) = qflx_tran_col(c) + qflx_tran_veg_pft(p) * pwtcol(p)
!         end if
!
!       end do
!    end do
!    ! Compute the Transpiration sink vertical distribution
!    do j = 1, nlevsoi
!        do fc = 1, num_soilc
!           c = filter_soilc(fc)
!           if (rootfr_sum(c) /= 0._r8) then
!                rootfr_col(c,j) = rootfr_col(c,j)/rootfr_sum(c)
!           end if
!        end do
!    end do
!
!!----------------------------------------------------------------------------------------------------------
!    ! (4) pass the clm_qflx to the vecs
!    ! NOTE the following unit conversions:
!    ! qflx_soil_top and qflx_tran_veg are in [mm/sec] from CLM;
!    ! qflx_clm_loc is in [kgH2O/sec] as mass rate for pflotran (as input)
!
!    ! previous time-step soil water pressure and saturation for adjusting qflx
!    ! note that this is a temporary workaround - waiting for PF's solution
!    call VecGetArrayF90(clm_pf_idata%press_clms, press_clms_loc, ierr)
!    CHKERRQ(ierr)
!    call VecGetArrayF90(clm_pf_idata%soillsat_clms, soillsat_clms_loc, ierr)
!    CHKERRQ(ierr)
!    call VecGetArrayF90(clm_pf_idata%effporosity_clms, porosity_clms_loc, ierr)
!    CHKERRQ(ierr)
!    call VecGetArrayF90(clm_pf_idata%sr_pcwmax_clms, sr_pcwmax_clms_loc, ierr)
!    CHKERRQ(ierr)
!
!    call VecGetArrayF90(clm_pf_idata%area_top_face_clms, area_clms_loc, ierr)
!    CHKERRQ(ierr)
!
!    call VecGetArrayF90(clm_pf_idata%zsoi_clmp, zsoi_clmp_loc, ierr)
!    CHKERRQ(ierr)
!
!    call VecGetArrayF90(clm_pf_idata%qflux_clmp, qflx_clmp_loc, ierr)
!    CHKERRQ(ierr)
!    call VecGetArrayF90(clm_pf_idata%press_subsurf_clmp, press_top_clmp_loc, ierr)
!    CHKERRQ(ierr)
!    call VecGetArrayF90(clm_pf_idata%press_subbase_clmp, press_base_clmp_loc, ierr)
!    CHKERRQ(ierr)
!    call VecGetArrayF90(clm_pf_idata%qflux_subsurf_clmp, qflux_top_clmp_loc, ierr)
!    CHKERRQ(ierr)
!    call VecGetArrayF90(clm_pf_idata%qflux_subbase_clmp, qflux_base_clmp_loc, ierr)
!    CHKERRQ(ierr)
!    call VecGetArrayF90(clm_pf_idata%press_maxponding_clmp, press_maxponding_clmp_loc, ierr)
!    CHKERRQ(ierr)
!
!    ! Initialize flux variables and bcs
!    do fc = 1, num_soilc
!
!      c = filter_soilc(fc)
!      g = cgridcell(c)
!      gcount = g - bounds%begg
!
!      do j = 1, nlevsoi
!        cellcount = gcount*clm_pf_idata%nzclm_mapped + j
!
!        if(j<=clm_pf_idata%nzclm_mapped) then
!          qflx_clmp_loc(cellcount) = 0.0_r8
!
!          if (j .eq. 1) then
!             press_top_clmp_loc(gcount+1)  = press_clms_loc(cellcount)   ! same as the first top layer
!          end if
!
!          if (j .eq. clm_pf_idata%nzclm_mapped) then
!             press_base_clmp_loc(gcount+1) = press_clms_loc(cellcount)   ! same as the bottom layer
!          end if
!
!        endif
!      end do
!
!    end do
!
!    pondmax(:) = 0._r8   ! this is temporarily set (not yet figure out how CLM get this value)
!    do fc = 1, num_soilc
!
!       c = filter_soilc(fc)
!       g = cgridcell(c)
!       gcount = g - bounds%begg
!       area = area_clms_loc(gcount*clm_pf_idata%nzclm_mapped+1)
!       reference_pressure = clm_pf_idata%pressure_reference
!       ponding_pressure = pondmax(c)*SHR_CONST_G              ! max. ponding water depth (mm) ==> pressure (Pa)
!
!       press_maxponding_clmp_loc(gcount+1) = reference_pressure+ponding_pressure*col%wtgcell(c)
!
!       do j = 1, nlevsoi
!         cellcount = gcount*clm_pf_idata%nzclm_mapped + j
!
!         if(j<=clm_pf_idata%nzclm_mapped) then
!           qflx = 0._r8
!           ! top BC
!           if (j .eq. 1) then
!
!             ! net liq water input/output to soil column
!             qflx_ground = qflx_top_soil(c) - qflx_evap_col(c)  ! unit: mm/sec
!
!             ! if net input potential, it's forming TOP BC of pressure type (water ponding potetial)
!             ! both waterhead and flux calcuated here, but not applied in PFLOTRAN in the same time (upon BC type picked-up by PF)
!             if ( qflx_ground .gt. 0._r8) then
!                ! Newly ADDED mmH2O ==> pressure (Pa) as top BC (dirichlet) by forming a layer of surface water column
!                ! AND, the actual infiltration/runoff are retrieving from PFLOTRAN using 'update_surflow_pf2clm' subroutine
!                if (soillsat_clms_loc(cellcount) >= 1._r8) then
!                   ! water-head formed on saturated below-ground soil layer
!                   press_top_clmp_loc(gcount+1) = press_clms_loc(cellcount) + &
!                          qflx_ground*col%wtgcell(c)*dtime*SHR_CONST_G
!                else
!                   ! ground-water-head discontinued from below-ground (atm. pressure applied at both ends)
!                   press_top_clmp_loc(gcount+1) = press_top_clmp_loc(gcount+1) +  &
!                          qflx_ground*col%wtgcell(c)*dtime*SHR_CONST_G
!                endif
!
!                ! mmH2O/sec ==> mH2O/sec of potential infiltration (flux) rate as top BC (neumann)
!                ! must be used together with 'seepage' as top BC in the meantime to removing upwarding water
!                ! AND, the actual infiltration/runoff are retrieving from PFLOTRAN using 'update_bcflow_pf2clm' subroutine
!                qflux_top_clmp_loc(gcount+1) = qflx_ground*1.e-3
!
!             ! if net loss potential, it's as source/sink term of soil column
!             else
!                qflx = qflx + min(0._r8, qflx_ground)              ! unit here: mmH2O/sec
!                qflux_top_clmp_loc(gcount+1) = 0._r8
!             end if
!          end if
!
!          ! adding plant root extraction of water (transpiration)
!          qflx = qflx - qflx_tran_col(c)*rootfr_col(c,j)  ! by this point: unit: mmH2O/sec for CLM column
!
!          qflx = qflx * col%wtgcell(c)                    ! from now on: per PF 3-D cells
!          qflx = qflx * area * 1.e-3 *denh2o              ! unit: mmH2O/sec ==> kgH2O/sec
!
!          ! previous time-step soil water saturation for adjusting qflx to avoid too wet or too dry to cause PF math issue
!          ! (this is a temporary workaround - waiting for PF's solution)
!          soilvwc = soillsat_clms_loc(cellcount) *  &
!                    porosity_clms_loc(cellcount)      ! PF saturation ==> real vwc (using adjusted porosity???)
!
!          ! checking if over-filled when sinking (excluding infiltration)
!          qflx_sink = max(0._r8, qflx)             ! sink (+) only (kgH2O/sec)
!          dsoilliq1 = (0.99_r8*porosity_clms_loc(cellcount)-soilvwc) &
!                  * zsoi_clmp_loc(cellcount) &
!                  * area * denh2o /dtime    ! mH2O ==> kgH2O/sec to be filled at most (1% for error)
!          qflx_sink = min(qflx_sink, max(0._r8,dsoilliq1))
!
!          ! checking if too dry to be ETed (or other sourced): lower than 'sr_pcwmax'
!          qflx_source = min(0._r8, qflx)           ! source (-) only (kgH2O/sec)
!          sr = 1.001_r8*sr_pcwmax_clms_loc(cellcount) * &              ! '1.001' will give 0.1% for holding error in the calculation
!                    watsat(c,j)                                                 ! PF saturation ==> 'real' vwc
!          dsoilliq2 = (sr-soilvwc)*zsoi_clmp_loc(cellcount)*area*denh2o/dtime                     ! mH2O ==> kgH2O/sec to be extracted at most (-)
!          qflx_source = max(qflx_source, min(0._r8,dsoilliq2))
!
!          qflx_clmp_loc(cellcount) = qflx_clmp_loc(cellcount)+ (qflx_sink+qflx_source)             ! source/sink unit: kg/sec
!
!          ! bottom BC (neumman type): m/sec
!          if (j .eq. clm_pf_idata%nzclm_mapped) then
!             ! available water flux-out rate (-) adjusted by source(-)/sink(+) term
!             dsoilliq3 = min(0._r8, dsoilliq2 - qflx_clmp_loc(cellcount)) &
!                             /area/denh2o                                      ! kgH2O/sec ==> mH2O/sec
!
!             ! free drainage at bottom
!             tempreal = soilvwc/watsat(c,j)                                    ! using 'real' saturation
!             kbot = hksat(c,j)*(tempreal**(2._r8*bsw(c,j)+3._r8))*1.e-3        ! mmH2O/sec ==> mH2O/sec
!             qflux_base_clmp_loc(gcount+1) =  qflux_base_clmp_loc(gcount+1) &
!                    + max(dsoilliq3, -kbot * col%wtgcell(c))             ! mH2O/sec
!
!          end if
!
!         endif !if(j<=clm_pf_idata%nzclm_mapped)
!
!       end do  ! do j=1,nlevsoi
!
!    end do ! do c=1, num_soilc
!
!    call VecRestoreArrayF90(clm_pf_idata%press_clms, press_clms_loc, ierr)
!    CHKERRQ(ierr)
!    call VecRestoreArrayF90(clm_pf_idata%soillsat_clms, soillsat_clms_loc, ierr)
!    CHKERRQ(ierr)
!    call VecRestoreArrayF90(clm_pf_idata%sr_pcwmax_clms, sr_pcwmax_clms_loc, ierr)
!    CHKERRQ(ierr)
!    call VecRestoreArrayF90(clm_pf_idata%effporosity_clms, porosity_clms_loc, ierr)
!    CHKERRQ(ierr)
!
!    call VecRestoreArrayF90(clm_pf_idata%area_top_face_clms, area_clms_loc, ierr)
!    CHKERRQ(ierr)
!
!    call VecRestoreArrayF90(clm_pf_idata%zsoi_clmp, zsoi_clmp_loc, ierr)
!    CHKERRQ(ierr)
!
!    call VecRestoreArrayF90(clm_pf_idata%qflux_clmp, qflx_clmp_loc, ierr)
!    CHKERRQ(ierr)
!    call VecRestoreArrayF90(clm_pf_idata%press_subsurf_clmp, press_top_clmp_loc, ierr)
!    CHKERRQ(ierr)
!    call VecRestoreArrayF90(clm_pf_idata%press_subbase_clmp, press_base_clmp_loc, ierr)
!    CHKERRQ(ierr)
!    call VecRestoreArrayF90(clm_pf_idata%qflux_subsurf_clmp, qflux_top_clmp_loc, ierr)
!    CHKERRQ(ierr)
!    call VecRestoreArrayF90(clm_pf_idata%qflux_subbase_clmp, qflux_base_clmp_loc, ierr)
!    CHKERRQ(ierr)
!    call VecRestoreArrayF90(clm_pf_idata%press_maxponding_clmp, press_maxponding_clmp_loc, ierr)
!    CHKERRQ(ierr)
!
!!----------------------------------------------------------------------------------------------------------
!
!  end associate
!  end subroutine get_clm_bcwflx
!
!  !-----------------------------------------------------------------------------
!  !BOP
!  !
!  ! !IROUTINE: get_clm_bceflx
!  !
!  ! !INTERFACE:
!  subroutine get_clm_bceflx(bounds,          &
!       num_soilc, filter_soilc,              &
!       clm_a2l, cws_vars, ces_vars, cef_vars)
!  !
!  ! !DESCRIPTION:
!  !
!  !  F.-M. YUAN: the boundary heat fluxes in CLM4.5 are extracted to drive pflotran TH mode.
!  !        GB only defined ground-heaf-flux.
!  !        So, this subroutine is a collective setting on either heat-flux (neumann type)
!  !        or interface thermal state (temperature) (dirichlet type)
!  !        at both ground and bottom interface (BC).
!  !
!  ! !USES:
!
!  ! !ARGUMENTS:
!    implicit none
!
!#include "finclude/petscsys.h"
!#include "finclude/petscvec.h"
!#include "finclude/petscvec.h90"
!
!    type(bounds_type) , intent(in)  :: bounds
!    integer, intent(in) :: num_soilc                 ! number of column non-lake points in column filter
!    integer, intent(in) :: filter_soilc(:)           ! column filter for non-lake points
!
!    type(atm2lnd_type)       , intent(in) :: clm_a2l
!    type(waterstate_type)    , intent(in) :: cws_vars
!    type(temperature_type)   , intent(in) :: ces_vars
!    type(energyflux_type)    , intent(in) :: cef_vars
!
!  ! !LOCAL VARIABLES:
!    integer  :: fc, c, g, p, gcount            ! do loop indices
!    integer  :: pftindex
!
!    ! for CLM (mpi) --> PF
!    PetscScalar, pointer :: gflux_subsurf_clmp_loc(:)    !   BC in neumman type: unit W/m2
!    PetscScalar, pointer :: gtemp_subsurf_clmp_loc(:)    !   BC in dirichlet type: unit in degC
!    PetscScalar, pointer :: gflux_subbase_clmp_loc(:)    !   BC in neumman type: unit W/m2
!    PetscScalar, pointer :: gtemp_subbase_clmp_loc(:)    !   BC in dirichlet type: unit in degC
!    PetscErrorCode :: ierr
!
!  !EOP
!  !-----------------------------------------------------------------------
!    associate ( &
!    cgridcell         => col%gridcell                            , &! column's gridcell
!    clandunit         => col%landunit                            , &! column's landunit
!    dz                => col%dz                                  , &! layer thickness depth (m)
!    pfti              => col%pfti                                , &! beginning pft index for each column
!    pwtgcell          => pft%wtgcell                             , &! weight relative to gridcell for each pft
!    pwtcol            => pft%wtcol                               , &! weight relative to column for each pft
!    !
!    frac_sno          => cws_vars%frac_sno_eff_col               , & ! Input: fraction of ground covered by snow (0 to 1)
!    frac_h2osfc       => cws_vars%frac_h2osfc_col                , & ! Input: fraction of ground covered by surface water (0 to 1)
!    !
!    eflx_bot          => cef_vars%eflx_bot_col                   , &! heat flux from beneath column (W/m**2) [+ = upward]
!    eflx_gnet         => cef_vars%eflx_gnet_patch                , &! net ground heat flux into the surface (W/m**2) per patch
!    eflx_soil_grnd    => cef_vars%eflx_soil_grnd_patch           , &! soil heat flux (W/m**2) [+ = into soil]
!    t_grnd            => ces_vars%t_grnd_col                   & ! ground surface temperature [K]
!    )
!
! !----------------------------------------------------------------------------------------------------------
!    ! (1) pass the clm_gflux/gtemp to the vec
!
!    call VecGetArrayF90(clm_pf_idata%gflux_subsurf_clmp, gflux_subsurf_clmp_loc, ierr)
!    CHKERRQ(ierr)
!    call VecGetArrayF90(clm_pf_idata%gflux_subbase_clmp, gflux_subbase_clmp_loc, ierr)
!    CHKERRQ(ierr)
!    call VecGetArrayF90(clm_pf_idata%gtemp_subsurf_clmp, gtemp_subsurf_clmp_loc, ierr)
!    CHKERRQ(ierr)
!    call VecGetArrayF90(clm_pf_idata%gtemp_subbase_clmp, gtemp_subbase_clmp_loc, ierr)
!    CHKERRQ(ierr)
!
!    gflux_subsurf_clmp_loc(:) = 0._r8
!    gflux_subbase_clmp_loc(:) = 0._r8
!    gtemp_subsurf_clmp_loc(:) = 0._r8
!    gtemp_subbase_clmp_loc(:) = 0._r8
!
!    ! CLM appears NO column-level ground-heat-flux variable, instead by 'patch'
!    do fc = 1, num_soilc
!       c = filter_soilc(fc)
!       if ( col%wtgcell(c) <= 0._r8 .or. (.not.col%active(c)) ) cycle     ! don't assign data from PF for inactive cell
!       g = cgridcell(c)
!       gcount = g - bounds%begg
!
!       do pftindex = 1, max_patch_per_col
!          if (pftindex <= col%npfts(c)) then
!             p = pfti(c) + pftindex - 1
!             if (pwtgcell(p)>0._r8) then
!               gflux_subsurf_clmp_loc(gcount+1)  = gflux_subsurf_clmp_loc(gcount+1) &
!                  + eflx_soil_grnd(p)*1.d-3 * pwtcol(p)           ! (TODO - checking) from W/m2 --> kJ/m2/s
!            end if
!         end if
!       end do
!    end do
!
!    ! CLM column-level variables available to PFLOTRAN
!    do fc = 1,num_soilc
!       c = filter_soilc(fc)
!       if ( col%wtgcell(c) <= 0._r8 .or. (.not.col%active(c)) ) cycle     ! don't assign data from PF for inactive cell
!
!       g = cgridcell(c)
!       gcount = g - bounds%begg
!
!       gflux_subbase_clmp_loc(gcount+1)  =  gflux_subbase_clmp_loc(gcount+1) &
!          + eflx_bot(c)*1.d-3 * col%wtgcell(c)
!       gtemp_subsurf_clmp_loc(gcount+1)  = gflux_subbase_clmp_loc(gcount+1) &
!          + (t_grnd(c) - tfrz) * col%wtgcell(c)
!
!       gtemp_subbase_clmp_loc(gcount+1)  = 0._r8                 ! not yet get it from CLM (i.e.,dirichlet type not available)
!
!    end do
!
!    call VecRestoreArrayF90(clm_pf_idata%gflux_subsurf_clmp, gflux_subsurf_clmp_loc, ierr)
!    CHKERRQ(ierr)
!    call VecRestoreArrayF90(clm_pf_idata%gflux_subbase_clmp, gflux_subbase_clmp_loc, ierr)
!    CHKERRQ(ierr)
!    call VecRestoreArrayF90(clm_pf_idata%gtemp_subsurf_clmp, gtemp_subsurf_clmp_loc, ierr)
!    CHKERRQ(ierr)
!    call VecRestoreArrayF90(clm_pf_idata%gtemp_subbase_clmp, gtemp_subbase_clmp_loc, ierr)
!    CHKERRQ(ierr)
!
!  end associate
!  end subroutine get_clm_bceflx


  !-----------------------------------------------------------------------------
  !BOP
  !
  ! !ROUTINE: get_clm_bgc_conc(bounds)
  !
  ! !INTERFACE:

  subroutine get_clm_bgc_state(clm_bgc_data,                    &
                        bounds, num_soilc, filter_soilc,        &
                        carbonstate_vars, nitrogenstate_vars,   &
                        phosphorusstate_vars)


    implicit none

    type(bounds_type)           , intent(in) :: bounds
    integer                     , intent(in) :: num_soilc         ! number of soil columns in filter
    integer                     , intent(in) :: filter_soilc(:)   ! filter for soil columns

    type(carbonstate_type)      , intent(in) :: carbonstate_vars
    type(nitrogenstate_type)    , intent(in) :: nitrogenstate_vars
    type(phosphorusstate_type)  , intent(in) :: phosphorusstate_vars
!    type(ch4_type)           , intent(in) :: ch4_vars
    type(clm_bgc_interface_data_type), intent(inout) :: clm_bgc_data

    character(len=256) :: subname = "get_clm_bgc_state"

    ! Local variables
    integer  :: fc, c, j, k
!    integer  :: gcount, cellcount
!    real(r8) :: wtgcell, realc_gcell, realn_gcell


    !------------------------------------------------------------------------------------------
    !
    associate ( &
       decomp_cpools_vr=> carbonstate_vars%decomp_cpools_vr_col     , &      ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
       decomp_npools_vr=> nitrogenstate_vars%decomp_npools_vr_col   , &      ! (gN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
       decomp_ppools_vr=> phosphorusstate_vars%decomp_ppools_vr_col , & ! [real(r8) (:,:,:) ! col (gP/m3) vertically-resolved decomposing (litter, cwd, soil) P pools

       smin_no3_vr     => nitrogenstate_vars%smin_no3_vr_col        , &      ! (gN/m3) vertically-resolved soil mineral NO3
       smin_nh4_vr     => nitrogenstate_vars%smin_nh4_vr_col        , &      ! (gN/m3) vertically-resolved soil mineral NH4
       smin_nh4sorb_vr => nitrogenstate_vars%smin_nh4sorb_vr_col    , &      ! (gN/m3) vertically-resolved soil mineral NH4 absorbed

       solutionp_vr    => phosphorusstate_vars%solutionp_vr_col     , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil solution P
       labilep_vr      => phosphorusstate_vars%labilep_vr_col       , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil labile mineral P
       secondp_vr      => phosphorusstate_vars%secondp_vr_col       , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil secondary mineralP
       sminp_vr        => phosphorusstate_vars%sminp_vr_col         , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil mineral P = solutionp + labilep + secondp
       occlp_vr        => phosphorusstate_vars%occlp_vr_col         , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil occluded mineral P
       primp_vr        => phosphorusstate_vars%primp_vr_col           & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil primary mineral P
    )
!
    do fc = 1, num_soilc
        c = filter_soilc(fc)
        do j = 1, nlevdecomp
            do k = 1, ndecomp_pools
                clm_bgc_data%decomp_cpools_vr_col(c,j,k)    = decomp_cpools_vr(c,j,k)
                clm_bgc_data%decomp_npools_vr_col(c,j,k)    = decomp_npools_vr(c,j,k)
                clm_bgc_data%decomp_ppools_vr_col(c,j,k)    = decomp_ppools_vr(c,j,k)
            end do

            clm_bgc_data%smin_no3_vr_col(c,j)           = smin_no3_vr(c,j)
            clm_bgc_data%smin_nh4_vr_col(c,j)           = smin_nh4_vr(c,j)
            clm_bgc_data%smin_nh4sorb_vr_col(c,j)       = smin_nh4sorb_vr(c,j)

            clm_bgc_data%solutionp_vr_col(c,j)          = solutionp_vr(c,j)
            clm_bgc_data%labilep_vr_col(c,j)            = labilep_vr(c,j)
            clm_bgc_data%secondp_vr_col(c,j)            = secondp_vr(c,j)
            clm_bgc_data%sminp_vr_col(c,j)              = solutionp_vr(c,j) + labilep_vr(c,j) + secondp_vr(c,j)
            clm_bgc_data%occlp_vr_col(c,j)              = occlp_vr(c,j)
            clm_bgc_data%primp_vr_col(c,j)              = primp_vr(c,j)
        end do
    end do
!

!write(*,'(A,10E14.6)')">>>DEBUG | dzsoi_decomp(nlevsoi)=",dzsoi_decomp(1:nlevsoi)
!write(*,'(A,50(1h-))')">>>DEBUG | get_clm_bgc_conc,lev=1 for C & N"
!write(*,'(12A14)')"lit1","lit2","lit3","cwd","som1","som2","som3","som4","no3","nh4","nh4sorb"
!write(*,'(12E14.6)')decomp_cpools_vr(1,1,1:8)
!write(*,'(12E14.6)')decomp_npools_vr(1,1,1:8),smin_no3_vr(1,1),smin_nh4_vr(1,1),smin_nh4sorb_vr(1,1)
!write(*,'(A30,12I14)')"DEBUG | nlevdecomp=",(j,j=1,nlevdecomp)
!write(*,'(A30,12E14.6)')"DEBUG | cpool_lit1=",decomp_cpools_vr_lit1_clm_loc(1:nlevdecomp)*clm_pf_idata%C_molecular_weight
!write(*,'(A30,12E14.6)')"DEBUG | cpool_lit2=",decomp_cpools_vr_lit2_clm_loc(1:nlevdecomp)*clm_pf_idata%C_molecular_weight
!write(*,'(A30,12E14.6)')"DEBUG | cpool_lit3=",decomp_cpools_vr_lit3_clm_loc(1:nlevdecomp)*clm_pf_idata%C_molecular_weight
!write(*,'(A30,12E14.6)')"DEBUG | cpool_cwd=",decomp_cpools_vr_cwd_clm_loc(1:nlevdecomp)*clm_pf_idata%C_molecular_weight
!write(*,'(A30,12E14.6)')"DEBUG | cpool_som1=",decomp_cpools_vr_som1_clm_loc(1:nlevdecomp)*clm_pf_idata%C_molecular_weight
!write(*,'(A30,12E14.6)')"DEBUG | cpool_som2=",decomp_cpools_vr_som2_clm_loc(1:nlevdecomp)*clm_pf_idata%C_molecular_weight
!write(*,'(A30,12E14.6)')"DEBUG | cpool_som3=",decomp_cpools_vr_som3_clm_loc(1:nlevdecomp)*clm_pf_idata%C_molecular_weight
!write(*,'(A30,12E14.6)')"DEBUG | cpool_som4=",decomp_cpools_vr_som4_clm_loc(1:nlevdecomp)*clm_pf_idata%C_molecular_weight
!
!write(*,'(A30,12E14.6)')"DEBUG | npool_lit1=",decomp_npools_vr_lit1_clm_loc(1:nlevdecomp)*clm_pf_idata%N_molecular_weight
!write(*,'(A30,12E14.6)')"DEBUG | npool_lit2=",decomp_npools_vr_lit2_clm_loc(1:nlevdecomp)*clm_pf_idata%N_molecular_weight
!write(*,'(A30,12E14.6)')"DEBUG | npool_lit3=",decomp_npools_vr_lit3_clm_loc(1:nlevdecomp)*clm_pf_idata%N_molecular_weight
!write(*,'(A30,12E14.6)')"DEBUG | npool_cwd=",decomp_npools_vr_cwd_clm_loc(1:nlevdecomp)*clm_pf_idata%N_molecular_weight
!write(*,'(A30,12E14.6)')"DEBUG | npool_som1=",decomp_npools_vr_som1_clm_loc(1:nlevdecomp)*clm_pf_idata%N_molecular_weight
!write(*,'(A30,12E14.6)')"DEBUG | npool_som2=",decomp_npools_vr_som2_clm_loc(1:nlevdecomp)*clm_pf_idata%N_molecular_weight
!write(*,'(A30,12E14.6)')"DEBUG | npool_som3=",decomp_npools_vr_som3_clm_loc(1:nlevdecomp)*clm_pf_idata%N_molecular_weight
!write(*,'(A30,12E14.6)')"DEBUG | npool_som4=",decomp_npools_vr_som4_clm_loc(1:nlevdecomp)*clm_pf_idata%N_molecular_weight
!
!write(*,'(A30,12E14.6)')"DEBUG | no3=",smin_no3_vr_clm_loc(1:nlevdecomp)*clm_pf_idata%N_molecular_weight
!write(*,'(A30,12E14.6)')"DEBUG | nh4=",smin_nh4_vr_clm_loc(1:nlevdecomp)*clm_pf_idata%N_molecular_weight
!write(*,'(A30,12E14.6)')"DEBUG | nh4sorb=",smin_nh4sorb_vr_clm_loc(1:nlevdecomp)*clm_pf_idata%N_molecular_weight

  end associate
  end subroutine get_clm_bgc_state

  !-----------------------------------------------------------------------------
  !
  ! !IROUTINE: get_clm_bgc_rate()
  !
  ! !INTERFACE:
  subroutine get_clm_bgc_flux(clm_bgc_data,                                 &
                        bounds, num_soilc, filter_soilc,                    &
                        cnstate_vars, carbonflux_vars, nitrogenflux_vars,   &
                        phosphorusflux_vars)

  !
  ! !DESCRIPTION:
  !
  !
  ! !USES:

    use clm_time_manager, only : get_step_size, get_nstep
#ifndef FLEXIBLE_POOLS
    use clm_varpar,       only : i_met_lit, i_cel_lit, i_lig_lit, i_cwd
#endif

  ! !ARGUMENTS:
    implicit none

    type(bounds_type) , intent(in) :: bounds
    integer           , intent(in) :: num_soilc       ! number of soil columns in filter
    integer           , intent(in) :: filter_soilc(:) ! filter for soil columns

    type(cnstate_type)                  , intent(in) :: cnstate_vars
    type(carbonflux_type)               , intent(in) :: carbonflux_vars
    type(nitrogenflux_type)             , intent(in) :: nitrogenflux_vars
    type(phosphorusflux_type)           , intent(inout) :: phosphorusflux_vars
    type(clm_bgc_interface_data_type)   , intent(inout) :: clm_bgc_data

    character(len=256) :: subname = "get_clm_bgc_flux"


 ! !LOCAL VARIABLES:
    integer  :: fc, c, g, j, k                         ! do loop indices
!    integer  :: gcount, cellcount
!    real(r8) :: wtgcell, realc_gcell, realn_gcell

    real(r8) :: dtime                               ! land model time step (sec)

    ! ratios of NH4:NO3 in N deposition and fertilization (temporarily set here, will be as inputs)
    real(r8) :: r_nh4_no3_dep(bounds%begc:bounds%endc)
    real(r8) :: r_nh4_no3_fert(bounds%begc:bounds%endc)
    real(r8) :: fnh4_dep, fnh4_fert

    ! C/N source/sink rates as inputs for soil bgc

    !
    !---------------------------------------------------------------------------
    !
    associate ( &
      ! plant litering and removal + SOM/LIT vertical transport
      externalc_to_decomp_cpools_vr    => carbonflux_vars%externalc_to_decomp_cpools_col   , &
      externaln_to_decomp_npools_vr    => nitrogenflux_vars%externaln_to_decomp_npools_col , &
      externalp_to_decomp_ppools_vr    => phosphorusflux_vars%externalp_to_decomp_ppools_col , &
      ! inorg. nitrogen source
      ndep_to_sminn                  => nitrogenflux_vars%ndep_to_sminn_col                 , &
      nfix_to_sminn                  => nitrogenflux_vars%nfix_to_sminn_col                 , &
      fert_to_sminn                  => nitrogenflux_vars%fert_to_sminn_col                 , &
      soyfixn_to_sminn               => nitrogenflux_vars%soyfixn_to_sminn_col              , &
      supplement_to_sminn_vr         => nitrogenflux_vars%supplement_to_sminn_vr_col        , &
      !
      nfixation_prof                 => cnstate_vars%nfixation_prof_col                     , &
      ndep_prof                      => cnstate_vars%ndep_prof_col                          , &
!      activeroot_prof                => cnstate_vars%activeroot_prof_col                    , &
      ! inorg. nitrogen sink (if not going to be done in PF)
      no3_net_transport_vr           => nitrogenflux_vars%no3_net_transport_vr_col          , &
      ! inorg. nitrogen sink potential
      col_plant_ndemand_vr           => nitrogenflux_vars%plant_ndemand_vr_col              , &

      pdep_to_sminp                   => phosphorusflux_vars%pdep_to_sminp_col             , &
      ! assume pdep_prof = ndep_prof
      fert_p_to_sminp                 => phosphorusflux_vars%fert_p_to_sminp_col           , &
      supplement_to_sminp_vr          => phosphorusflux_vars%supplement_to_sminp_vr_col    , &

      sminp_net_transport_vr         => phosphorusflux_vars%sminp_net_transport_vr_col      , &
      col_plant_pdemand_vr           => phosphorusflux_vars%plant_pdemand_vr_col              &
    )

    dtime = get_step_size()


!
    r_nh4_no3_dep(:)  = 1.0_r8      ! temporarily assuming half of N dep is in NH4 and another half in NO3
    r_nh4_no3_fert(:) = 1.0_r8      ! temporarily assiming half of N fertilization is in NH4 and another half in NO3
!
    do fc = 1,num_soilc
        c = filter_soilc(fc)
        fnh4_dep  = max(0._r8, min(1.0_r8, 1._r8/(r_nh4_no3_dep(c)+1._r8)))
        fnh4_fert = max(0._r8, min(1.0_r8, 1._r8/(r_nh4_no3_fert(c)+1._r8)))

        do j = 1, nlevdecomp
            do k = 1, ndecomp_pools
                clm_bgc_data%externalc_to_decomp_cpools_col(c,j,k)  = externalc_to_decomp_cpools_vr(c,j,k)
                clm_bgc_data%externaln_to_decomp_npools_col(c,j,k)  = externaln_to_decomp_npools_vr(c,j,k)
                clm_bgc_data%externalp_to_decomp_ppools_col(c,j,k)  = externalp_to_decomp_ppools_vr(c,j,k)
            end do

            clm_bgc_data%externaln_to_nh4_col(c,j)          =   fnh4_dep*ndep_to_sminn(c) * ndep_prof(c,j) +  &
                                                                fnh4_fert*fert_to_sminn(c) * ndep_prof(c,j) + &
                                                                fnh4_fert*supplement_to_sminn_vr(c,j) +       &
                                                                nfix_to_sminn(c) * nfixation_prof(c,j) +      &
                                                                soyfixn_to_sminn(c) * nfixation_prof(c,j)

            clm_bgc_data%externaln_to_no3_col(c,j)          =   (1._r8-fnh4_dep)*ndep_to_sminn(c) * ndep_prof(c, j) +  &
                                                                (1._r8-fnh4_fert)*fert_to_sminn(c) * ndep_prof(c, j) + &
                                                                (1._r8-fnh4_fert)*supplement_to_sminn_vr(c,j)

            clm_bgc_data%externalp_to_primp_col(c,j)        =   pdep_to_sminp(c)*ndep_prof(c, j)
            clm_bgc_data%externalp_to_labilep_col(c,j)      =   fert_p_to_sminp(c)*ndep_prof(c, j)
            clm_bgc_data%externalp_to_solutionp_col(c,j)    =   supplement_to_sminp_vr(c,j)

            !! net flux to no3 = externaln_to_no3_col(c,j) - no3_net_transport_vr_col(c,j)
            clm_bgc_data%no3_net_transport_vr_col(c,j)      = no3_net_transport_vr(c,j)
            clm_bgc_data%sminp_net_transport_vr_col(c,j)    = sminp_net_transport_vr(c,j)  !!from solutionp

            clm_bgc_data%plant_ndemand_vr_col(c,j)          = col_plant_ndemand_vr(c,j)
            clm_bgc_data%plant_pdemand_vr_col(c,j)          = col_plant_pdemand_vr(c,j)

        end do
    end do

!write(*,'(A,50(1h-))')">>>DEBUG | get_clm_bgc_rate,lev=1 for C & N"
!write(*,'(12A14)')"lit1","lit2","lit3","cwd","som1","som2","som3","som4","no3","nh4","plantNdemand"
!write(*,'(12E14.6)')col_net_to_decomp_cpools_vr(1,1,1:ndecomp_pools)
!write(*,'(12E14.6)')col_net_to_decomp_npools_vr(1,1,1:ndecomp_pools),&
!                (rate_smin_no3_clm_loc(1))*clm_pf_idata%N_molecular_weight, &
!                (rate_smin_nh4_clm_loc(1))*clm_pf_idata%N_molecular_weight, &
!                (rate_plantndemand_clm_loc(1))*clm_pf_idata%N_molecular_weight

!write(*,'(A30,10E14.6)')">>>DEBUG | rate_lit1c=",rate_lit1c_clm_loc(1:10)*clm_pf_idata%C_molecular_weight
!write(*,'(A30,10E14.6)')">>>DEBUG | rate_lit2c=",rate_lit2c_clm_loc(1:10)*clm_pf_idata%C_molecular_weight
!write(*,'(A30,10E14.6)')">>>DEBUG | rate_lit3c=",rate_lit3c_clm_loc(1:10)*clm_pf_idata%C_molecular_weight
!write(*,'(A30,10E14.6)')">>>DEBUG | rate_cwdc=", rate_cwdc_clm_loc(1:10)*clm_pf_idata%C_molecular_weight
!write(*,'(A30,10E14.6)')">>>DEBUG | rate_som1c=",rate_som1c_clm_loc(1:10)*clm_pf_idata%C_molecular_weight
!write(*,'(A30,10E14.6)')">>>DEBUG | rate_som2c=",rate_som2c_clm_loc(1:10)*clm_pf_idata%C_molecular_weight
!write(*,'(A30,10E14.6)')">>>DEBUG | rate_som3c=",rate_som3c_clm_loc(1:10)*clm_pf_idata%C_molecular_weight
!write(*,'(A30,10E14.6)')">>>DEBUG | rate_som4c=",rate_som4c_clm_loc(1:10)*clm_pf_idata%C_molecular_weight
!!
!write(*,'(A30,10E14.6)')">>>DEBUG | rate_lit1n=",rate_lit1n_clm_loc(1:10)*clm_pf_idata%N_molecular_weight
!write(*,'(A30,10E14.6)')">>>DEBUG | rate_lit2n=",rate_lit2n_clm_loc(1:10)*clm_pf_idata%N_molecular_weight
!write(*,'(A30,10E14.6)')">>>DEBUG | rate_lit3n=",rate_lit3n_clm_loc(1:10)*clm_pf_idata%N_molecular_weight
!write(*,'(A30,10E14.6)')">>>DEBUG | rate_cwdn=", rate_cwdn_clm_loc(1:10)*clm_pf_idata%N_molecular_weight
!write(*,'(A30,10E14.6)')">>>DEBUG | rate_som1n=",rate_som1n_clm_loc(1:10)*clm_pf_idata%N_molecular_weight
!write(*,'(A30,10E14.6)')">>>DEBUG | rate_som2n=",rate_som2n_clm_loc(1:10)*clm_pf_idata%N_molecular_weight
!write(*,'(A30,10E14.6)')">>>DEBUG | rate_som3n=",rate_som3n_clm_loc(1:10)*clm_pf_idata%N_molecular_weight
!write(*,'(A30,10E14.6)')">>>DEBUG | rate_som4n=",rate_som4n_clm_loc(1:10)*clm_pf_idata%N_molecular_weight
!
!write(*,'(A30,10E14.6)')">>>DEBUG | rate_nh4=",(rate_smin_nh4_clm_loc(1:10))*clm_pf_idata%N_molecular_weight
!write(*,'(A30,10E14.6)')">>>DEBUG | rate_no3=",(rate_smin_no3_clm_loc(1:10))*clm_pf_idata%N_molecular_weight
!write(*,'(A30,10E14.6)')">>>DEBUG | rate_plantndemand=",(rate_plantndemand_clm_loc(1:10))*clm_pf_idata%N_molecular_weight

    end associate
  end subroutine get_clm_bgc_flux


  ! =============================UPDATE evolving variables to CLM ==========================================


  !-----------------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: update_soil_moisture
  !
  ! !INTERFACE:
  subroutine update_soil_moisture(clm_bgc_data,     &
           bounds, num_soilc, filter_soilc,   &
           waterstate_vars)

  !
  ! !DESCRIPTION:
  !
  !
  ! !USES:

  ! !ARGUMENTS:
    implicit none

    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_soilc        ! number of column soil points in column filter
    integer, intent(in) :: filter_soilc(:)  ! column filter for soil points
!    type(soilstate_type) , intent(in)    :: soilstate_vars
    type(waterstate_type), intent(inout) :: waterstate_vars
    type(clm_bgc_interface_data_type), intent(in) :: clm_bgc_data

  ! !LOCAL VARIABLES:
    integer  :: fc, c, j, g, gcount      ! indices

  !EOP
  !-----------------------------------------------------------------------
    associate ( &
!         gridcell   =>  col%gridcell   , & !  [integer (:)]  gridcell index of column
!         wtgcell    =>  col%wtgcell    , & !  [real(r8) (:)]  weight (relative to gridcell
!         dz         =>  col%dz         , & !  [real(r8) (:,:)]  layer thickness depth (m)
         !
         h2osoi_liq_col =>  waterstate_vars%h2osoi_liq_col      , &
         h2osoi_ice_col =>  waterstate_vars%h2osoi_ice_col      , &
         h2osoi_vol_col =>  waterstate_vars%h2osoi_vol_col        &
!         watsat_col     =>  soilstate_vars%watsat_col      &
    )

    do fc = 1,num_soilc
        c = filter_soilc(fc)
        do j = 1, nlevsoi
            h2osoi_liq_col(c,j) =  clm_bgc_data%h2osoi_liq_col(c,j)
            h2osoi_ice_col(c,j) =  clm_bgc_data%h2osoi_ice_col(c,j)
            h2osoi_vol_col(c,j) =  clm_bgc_data%h2osoi_vol_col(c,j)
        end do
    end do

    end associate
  end subroutine update_soil_moisture
!
!  !-----------------------------------------------------------------------------
!  !BOP
!  !
!  ! !IROUTINE: update_soil_temperature_pf2clm
!  !
!  ! !INTERFACE:
  subroutine update_soil_temperature(clm_bgc_data,     &
           bounds, num_soilc, filter_soilc,   &
           temperature_vars)

  !
  ! !DESCRIPTION:
  !
  !
  ! !USES:

  ! !ARGUMENTS:
    implicit none

    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_soilc        ! number of column soil points in column filter
    integer, intent(in) :: filter_soilc(:)  ! column filter for soil points
!    type(soilstate_type) , intent(in)    :: soilstate_vars
    type(temperature_type), intent(inout):: temperature_vars
    type(clm_bgc_interface_data_type), intent(in) :: clm_bgc_data

  ! !LOCAL VARIABLES:
    integer  :: fc, c, j, g, gcount      ! indices

  !EOP
  !-----------------------------------------------------------------------
    associate ( &
         gridcell   =>  col%gridcell   , & !  [integer (:)]  gridcell index of column
         wtgcell    =>  col%wtgcell    , & !  [real(r8) (:)]  weight (relative to gridcell
         dz         =>  col%dz         , & !  [real(r8) (:,:)]  layer thickness depth (m)
         !
         t_soisno   => temperature_vars%t_soisno_col   & ! snow-soil temperature (Kelvin)
    )
    do fc = 1,num_soilc
        c = filter_soilc(fc)
        do j = 1, nlevsoi
            t_soisno(c,j)   = temperature_vars%t_soisno_col(c,j)
        end do
    end do

    end associate
  end subroutine update_soil_temperature

  !-----------------------------------------------------------------------------
  !
  ! !ROUTINE: update_soil_bgc_pf2clm()
  !
  ! !INTERFACE:
  !
  ! ! calculating BGC state variable changes over one time-step (rates)
  !   NOTE: Don't update the organic C/N state variables, which will be updated in those 'update' subroutines
  !           and the 'CNSoilLittVertTranspMod.F90' after 'update1'.
  !
  subroutine update_bgc_state_decomp(clm_bgc_data,    &
           bounds, num_soilc, filter_soilc,         &
           carbonstate_vars, nitrogenstate_vars,    &
           phosphorusstate_vars                     &
           )

    use CNDecompCascadeConType, only : decomp_cascade_con
    use clm_time_manager, only : get_step_size
#ifndef FLEXIBLE_POOLS
    use clm_varpar, only : i_met_lit, i_cel_lit, i_lig_lit, i_cwd
#endif

    implicit none

    type(bounds_type)        , intent(in) :: bounds
    integer                  , intent(in) :: num_soilc         ! number of soil columns in filter
    integer                  , intent(in) :: filter_soilc(:)   ! filter for soil columns

!    type(atm2lnd_type)       , intent(in) :: clm_a2l
!    type(waterstate_type)    , intent(in) :: waterstate_vars
!    type(soilstate_type)     , intent(in) :: soilstate_vars
!    type(cnstate_type)       , intent(in) :: cnstate_vars
    type(carbonstate_type)      , intent(inout) :: carbonstate_vars
    type(nitrogenstate_type)    , intent(inout) :: nitrogenstate_vars
    type(phosphorusstate_type)  , intent(inout) :: phosphorusstate_vars
!    type(carbonflux_type)    , intent(inout) :: carbonflux_vars
!    type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars
!    type(ch4_type)           , intent(inout) :: ch4_vars

    type(clm_bgc_interface_data_type), intent(in) :: clm_bgc_data

    character(len=256) :: subname = "update_soil_bgc_state"

    integer  :: fc,c,j,k
!    integer  :: gcount, cellcount
!    real(r8) :: wtgcell

!    real(r8) :: dtime            ! land model time step (sec)

!------------------------------------------------------------------------------------
     !
     associate ( &
!     initial_cn_ratio             => decomp_cascade_con%initial_cn_ratio             , &
     decomp_cpools_vr             => carbonstate_vars%decomp_cpools_vr_col           , &
     decomp_npools_vr             => nitrogenstate_vars%decomp_npools_vr_col         , &
     decomp_ppools_vr             => phosphorusstate_vars%decomp_ppools_vr_col         &

!     sminn_vr                     => nitrogenstate_vars%sminn_vr_col                 , &
!     smin_no3_vr                  => nitrogenstate_vars%smin_no3_vr_col              , &
!     smin_nh4_vr                  => nitrogenstate_vars%smin_nh4_vr_col              , &
!     smin_nh4sorb_vr              => nitrogenstate_vars%smin_nh4sorb_vr_col          , &
!
!     solutionp_vr    => phosphorusstate_vars%solutionp_vr_col     , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil solution P
!     labilep_vr      => phosphorusstate_vars%labilep_vr_col       , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil labile mineral P
!     secondp_vr      => phosphorusstate_vars%secondp_vr_col       , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil secondary mineralP
!     sminp_vr        => phosphorusstate_vars%sminp_vr_col         , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil mineral P = solutionp + labilep + second
!     occlp_vr        => phosphorusstate_vars%occlp_vr_col         , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil occluded mineral P
!     primp_vr        => phosphorusstate_vars%primp_vr_col           & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil primary mineral P

!     decomp_cpools_delta_vr       => carbonflux_vars%decomp_cpools_sourcesink_col    , &
!     decomp_npools_delta_vr       => nitrogenflux_vars%decomp_npools_sourcesink_col  , &
!     sminn_to_plant_vr            => nitrogenflux_vars%sminn_to_plant_vr_col         , &
!     smin_no3_to_plant_vr         => nitrogenflux_vars%smin_no3_to_plant_vr_col      , &
!     smin_nh4_to_plant_vr         => nitrogenflux_vars%smin_nh4_to_plant_vr_col      , &
!     potential_immob_vr           => nitrogenflux_vars%potential_immob_vr_col        , &
!     actual_immob_vr              => nitrogenflux_vars%actual_immob_vr_col           , &
!     gross_nmin_vr                => nitrogenflux_vars%gross_nmin_vr_col               &
     )
! ------------------------------------------------------------------------
!     dtime = get_step_size()

     ! soil C/N pool increments set to the previous timestep (i.e., not yet updated)
!     decomp_cpools_delta_vr  = 0._r8-decomp_cpools_vr
!     decomp_npools_delta_vr  = 0._r8-decomp_npools_vr
!
    do fc = 1, num_soilc
        c = filter_soilc(fc)
        do j = 1, nlevdecomp
            do k = 1, ndecomp_pools
                decomp_cpools_vr(c,j,k) = clm_bgc_data%decomp_cpools_vr_col(c,j,k)
                decomp_npools_vr(c,j,k) = clm_bgc_data%decomp_npools_vr_col(c,j,k)
                decomp_ppools_vr(c,j,k) = clm_bgc_data%decomp_ppools_vr_col(c,j,k)
            end do
        end do
    end do

    end associate
  end subroutine update_bgc_state_decomp

    !-----------------------------------------------------------------------------
    subroutine update_bgc_state_smin(clm_bgc_data,      &
           bounds, num_soilc, filter_soilc,             &
           nitrogenstate_vars, phosphorusstate_vars)

    use CNDecompCascadeConType, only : decomp_cascade_con
    use clm_time_manager, only : get_step_size
!#ifndef FLEXIBLE_POOLS
!    use clm_varpar, only : i_met_lit, i_cel_lit, i_lig_lit, i_cwd
!#endif

    implicit none

    type(bounds_type)        , intent(in) :: bounds
    integer                  , intent(in) :: num_soilc         ! number of soil columns in filter
    integer                  , intent(in) :: filter_soilc(:)   ! filter for soil columns

!    type(atm2lnd_type)       , intent(in) :: clm_a2l
!    type(waterstate_type)    , intent(in) :: waterstate_vars
!    type(soilstate_type)     , intent(in) :: soilstate_vars
!    type(cnstate_type)       , intent(in) :: cnstate_vars
!    type(carbonstate_type)      , intent(inout) :: carbonstate_vars
    type(nitrogenstate_type)    , intent(inout) :: nitrogenstate_vars
    type(phosphorusstate_type)  , intent(inout) :: phosphorusstate_vars
!    type(carbonflux_type)    , intent(inout) :: carbonflux_vars
!    type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars
!    type(ch4_type)           , intent(inout) :: ch4_vars

    type(clm_bgc_interface_data_type), intent(in) :: clm_bgc_data

    character(len=256) :: subname = "update_soil_bgc_state"

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

    integer  :: fc,c,j
!    integer  :: gcount, cellcount
!    real(r8) :: wtgcell
!
!    real(r8) :: dtime            ! land model time step (sec)

!------------------------------------------------------------------------------------
     !
     associate ( &
!     initial_cn_ratio             => decomp_cascade_con%initial_cn_ratio             , &
!     decomp_cpools_vr             => carbonstate_vars%decomp_cpools_vr_col           , &
!     decomp_npools_vr             => nitrogenstate_vars%decomp_npools_vr_col         , &
!     decomp_ppools_vr             => phosphorusstate_vars%decomp_ppools_vr_col       , &

     sminn_vr           => nitrogenstate_vars%sminn_vr_col                 , &
     smin_no3_vr        => nitrogenstate_vars%smin_no3_vr_col              , &
     smin_nh4_vr        => nitrogenstate_vars%smin_nh4_vr_col              , &
     smin_nh4sorb_vr    => nitrogenstate_vars%smin_nh4sorb_vr_col          , &

     solutionp_vr       => phosphorusstate_vars%solutionp_vr_col     , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil solution P
     labilep_vr         => phosphorusstate_vars%labilep_vr_col       , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil labile mineral P
     secondp_vr         => phosphorusstate_vars%secondp_vr_col       , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil secondary mineralP
     sminp_vr           => phosphorusstate_vars%sminp_vr_col         , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil mineral P = solutionp + labilep + second
     occlp_vr           => phosphorusstate_vars%occlp_vr_col         , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil occluded mineral P
     primp_vr           => phosphorusstate_vars%primp_vr_col           & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil primary mineral P

     )
! ------------------------------------------------------------------------
!     dtime = get_step_size()

     ! soil C/N pool increments set to the previous timestep (i.e., not yet updated)
!     decomp_cpools_delta_vr  = 0._r8-decomp_cpools_vr
!     decomp_npools_delta_vr  = 0._r8-decomp_npools_vr
!
    do fc = 1, num_soilc
        c = filter_soilc(fc)
        do j = 1, nlevdecomp
            smin_no3_vr(c,j)        = clm_bgc_data%smin_no3_vr_col(c,j)
            smin_nh4_vr(c,j)        = clm_bgc_data%smin_nh4_vr_col(c,j)
            smin_nh4sorb_vr(c,j)    = clm_bgc_data%smin_nh4sorb_vr_col(c,j)
            sminn_vr                = clm_bgc_data%sminn_vr_col(c,j)

            solutionp_vr(c,j)       = clm_bgc_data%solutionp_vr_col(c,j)
            labilep_vr(c,j)         = clm_bgc_data%labilep_vr_col(c,j)
            secondp_vr(c,j)         = clm_bgc_data%secondp_vr_col(c,j)
            sminp_vr(c,j)           = clm_bgc_data%sminp_vr_col(c,j)
            occlp_vr(c,j)           = clm_bgc_data%occlp_vr_col(c,j)
            primp_vr(c,j)           = clm_bgc_data%primp_vr_col(c,j)
        end do
    end do
write(*,'(A30,12E14.6)')"DEBUG | clm UPDATE no3=",smin_no3_vr(1,1:nlevdecomp)
write(*,'(A30,12E14.6)')"DEBUG | clm UPDATE nh4=",smin_nh4_vr(1,1:nlevdecomp)
    end associate
  end subroutine update_bgc_state_smin

    !-----------------------------------------------------------------------------
    subroutine update_bgc_flux_decomp_sourcesink(clm_bgc_data,       &
           bounds, num_soilc, filter_soilc,              &
           carbonflux_vars, nitrogenflux_vars, &
           phosphorusflux_vars)

    use CNDecompCascadeConType, only : decomp_cascade_con
    use clm_time_manager, only : get_step_size
!#ifndef FLEXIBLE_POOLS
!    use clm_varpar, only : i_met_lit, i_cel_lit, i_lig_lit, i_cwd
!#endif

    implicit none

    type(bounds_type)        , intent(in) :: bounds
    integer                  , intent(in) :: num_soilc         ! number of soil columns in filter
    integer                  , intent(in) :: filter_soilc(:)   ! filter for soil columns

!    type(atm2lnd_type)       , intent(in) :: clm_a2l
!    type(waterstate_type)    , intent(in) :: waterstate_vars
!    type(soilstate_type)     , intent(in) :: soilstate_vars
!    type(cnstate_type)       , intent(in) :: cnstate_vars
!    type(carbonstate_type)   , intent(inout) :: carbonstate_vars
    type(carbonflux_type)    , intent(inout) :: carbonflux_vars
!    type(nitrogenstate_type) , intent(inout) :: nitrogenstate_vars
    type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars
    type(phosphorusflux_type)  , intent(inout) :: phosphorusflux_vars
!    type(ch4_type)           , intent(inout) :: ch4_vars
    type(clm_bgc_interface_data_type), intent(in) :: clm_bgc_data

    integer :: fc, c, j, k
    character(len=256) :: subname = "update_soil_bgc_pf2clm"

    associate ( &
!     initial_cn_ratio             => decomp_cascade_con%initial_cn_ratio             , &
!     decomp_cpools_vr             => carbonstate_vars%decomp_cpools_vr_col           , &
!     decomp_npools_vr             => nitrogenstate_vars%decomp_npools_vr_col         , &
!     decomp_ppools_vr             => phosphorusstate_vars%decomp_ppools_vr_col       , &
!
!     sminn_vr                     => nitrogenstate_vars%sminn_vr_col                 , &
!     smin_no3_vr                  => nitrogenstate_vars%smin_no3_vr_col              , &
!     smin_nh4_vr                  => nitrogenstate_vars%smin_nh4_vr_col              , &
!     smin_nh4sorb_vr              => nitrogenstate_vars%smin_nh4sorb_vr_col          , &
!
!     solutionp_vr    => phosphorusstate_vars%solutionp_vr_col     , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil solution P
!     labilep_vr      => phosphorusstate_vars%labilep_vr_col       , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil labile mineral P
!     secondp_vr      => phosphorusstate_vars%secondp_vr_col       , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil secondary mineralP
!     sminp_vr        => phosphorusstate_vars%sminp_vr_col         , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil mineral P = solutionp + labilep + second
!     occlp_vr        => phosphorusstate_vars%occlp_vr_col         , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil occluded mineral P
!     primp_vr        => phosphorusstate_vars%primp_vr_col           & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil primary mineral P

     decomp_cpools_sourcesink_vr  => carbonflux_vars%decomp_cpools_sourcesink_col    , &
     decomp_npools_sourcesink_vr  => nitrogenflux_vars%decomp_npools_sourcesink_col  , &
     decomp_ppools_sourcesink_vr  => phosphorusflux_vars%decomp_ppools_sourcesink_col  &

!     sminn_to_plant_vr            => nitrogenflux_vars%sminn_to_plant_vr_col         , &
!     smin_no3_to_plant_vr         => nitrogenflux_vars%smin_no3_to_plant_vr_col      , &
!     smin_nh4_to_plant_vr         => nitrogenflux_vars%smin_nh4_to_plant_vr_col      , &
!     potential_immob_vr           => nitrogenflux_vars%potential_immob_vr_col        , &
!     actual_immob_vr              => nitrogenflux_vars%actual_immob_vr_col           , &
!     gross_nmin_vr                => nitrogenflux_vars%gross_nmin_vr_col               &
     )

    do fc = 1, num_soilc
        c = filter_soilc(fc)
        do j = 1, nlevdecomp
            do k = 1, ndecomp_pools
                decomp_cpools_sourcesink_vr(c,j,k) = clm_bgc_data%decomp_cpools_sourcesink_col(c,j,k)
                decomp_npools_sourcesink_vr(c,j,k) = clm_bgc_data%decomp_npools_sourcesink_col(c,j,k)
                decomp_ppools_sourcesink_vr(c,j,k) = clm_bgc_data%decomp_ppools_sourcesink_col(c,j,k)
            end do
        end do
    end do
    end associate
    end subroutine update_bgc_flux_decomp_sourcesink

    !-----------------------------------------------------------------------------
    subroutine update_bgc_flux_decomp_cascade(clm_bgc_data, &
           bounds, num_soilc, filter_soilc,                 &
           carbonflux_vars, nitrogenflux_vars,              &
           phosphorusflux_vars)

    use CNDecompCascadeConType, only : decomp_cascade_con
    use clm_time_manager, only : get_step_size
!#ifndef FLEXIBLE_POOLS
!    use clm_varpar, only : i_met_lit, i_cel_lit, i_lig_lit, i_cwd
!#endif

    implicit none

    type(bounds_type)        , intent(in) :: bounds
    integer                  , intent(in) :: num_soilc         ! number of soil columns in filter
    integer                  , intent(in) :: filter_soilc(:)   ! filter for soil columns

!    type(atm2lnd_type)       , intent(in) :: clm_a2l
!    type(waterstate_type)    , intent(in) :: waterstate_vars
!    type(soilstate_type)     , intent(in) :: soilstate_vars
!    type(cnstate_type)       , intent(in) :: cnstate_vars
!    type(carbonstate_type)   , intent(inout) :: carbonstate_vars
    type(carbonflux_type)    , intent(inout) :: carbonflux_vars
!    type(nitrogenstate_type) , intent(inout) :: nitrogenstate_vars
    type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars
    type(phosphorusflux_type)  , intent(inout) :: phosphorusflux_vars
!    type(ch4_type)           , intent(inout) :: ch4_vars
    type(clm_bgc_interface_data_type), intent(in) :: clm_bgc_data

    integer :: fc, c, j, k
    character(len=256) :: subname = "update_soil_bgc_pf2clm"

    associate ( &
!     initial_cn_ratio             => decomp_cascade_con%initial_cn_ratio             , &
!     decomp_cpools_vr             => carbonstate_vars%decomp_cpools_vr_col           , &
!     decomp_npools_vr             => nitrogenstate_vars%decomp_npools_vr_col         , &
!     decomp_ppools_vr             => phosphorusstate_vars%decomp_ppools_vr_col       , &
!
!     sminn_vr                     => nitrogenstate_vars%sminn_vr_col                 , &
!     smin_no3_vr                  => nitrogenstate_vars%smin_no3_vr_col              , &
!     smin_nh4_vr                  => nitrogenstate_vars%smin_nh4_vr_col              , &
!     smin_nh4sorb_vr              => nitrogenstate_vars%smin_nh4sorb_vr_col          , &
!
!     solutionp_vr    => phosphorusstate_vars%solutionp_vr_col     , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil solution P
!     labilep_vr      => phosphorusstate_vars%labilep_vr_col       , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil labile mineral P
!     secondp_vr      => phosphorusstate_vars%secondp_vr_col       , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil secondary mineralP
!     sminp_vr        => phosphorusstate_vars%sminp_vr_col         , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil mineral P = solutionp + labilep + second
!     occlp_vr        => phosphorusstate_vars%occlp_vr_col         , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil occluded mineral P
!     primp_vr        => phosphorusstate_vars%primp_vr_col           & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil primary mineral P

     decomp_cascade_hr_vr_col         => carbonflux_vars%decomp_cascade_hr_vr_col           , &

     decomp_cascade_ctransfer_vr_col  => carbonflux_vars%decomp_cascade_ctransfer_vr_col    , &
     decomp_cascade_ntransfer_vr_col  => nitrogenflux_vars%decomp_cascade_ntransfer_vr_col  , &
     decomp_cascade_ptransfer_vr_col  => phosphorusflux_vars%decomp_cascade_ptransfer_vr_col, &

     decomp_cascade_sminn_flux_vr_col => nitrogenflux_vars%decomp_cascade_sminn_flux_vr_col     , & ! Output: [real(r8) (:,:,:) ]  vert-res mineral N flux for transition along decomposition cascade (gN/m3/s)
     decomp_cascade_sminp_flux_vr_col => phosphorusflux_vars%decomp_cascade_sminp_flux_vr_col     , & ! Output: [real(r8) (:,:,:) ]  vert-res mineral P flux for transition along decomposition cascade (gP/m3/s)

     sminn_to_denit_decomp_cascade_vr_col => nitrogenflux_vars%sminn_to_denit_decomp_cascade_vr_col & ! Output: [real(r8) (:,:,:) ]

!     sminn_to_plant_vr            => nitrogenflux_vars%sminn_to_plant_vr_col         , &
!     smin_no3_to_plant_vr         => nitrogenflux_vars%smin_no3_to_plant_vr_col      , &
!     smin_nh4_to_plant_vr         => nitrogenflux_vars%smin_nh4_to_plant_vr_col      , &
!     potential_immob_vr           => nitrogenflux_vars%potential_immob_vr_col        , &
!     actual_immob_vr              => nitrogenflux_vars%actual_immob_vr_col           , &
!     gross_nmin_vr                => nitrogenflux_vars%gross_nmin_vr_col               &
     )

    do fc = 1, num_soilc
        c = filter_soilc(fc)
        do j = 1, nlevdecomp
            do k = 1, ndecomp_pools
                decomp_cascade_hr_vr_col(c,j,k)         = clm_bgc_data%decomp_cascade_hr_vr_col(c,j,k)
                decomp_cascade_ctransfer_vr_col(c,j,k)  = clm_bgc_data%decomp_cascade_ctransfer_vr_col(c,j,k)
                decomp_cascade_ntransfer_vr_col(c,j,k)  = clm_bgc_data%decomp_cascade_ntransfer_vr_col(c,j,k)
                decomp_cascade_ptransfer_vr_col(c,j,k)  = clm_bgc_data%decomp_cascade_ptransfer_vr_col(c,j,k)

                decomp_cascade_sminn_flux_vr_col(c,j,k)     = clm_bgc_data%decomp_cascade_sminn_flux_vr_col(c,j,k)
                decomp_cascade_sminp_flux_vr_col(c,j,k)     = clm_bgc_data%decomp_cascade_sminp_flux_vr_col(c,j,k)
                sminn_to_denit_decomp_cascade_vr_col(c,j,k) = clm_bgc_data%sminn_to_denit_decomp_cascade_vr_col(c,j,k)
            end do
        end do
    end do
    end associate
    end subroutine update_bgc_flux_decomp_cascade

    !-----------------------------------------------------------------------------
    subroutine update_bgc_flux_smin(clm_bgc_data,   &
           bounds, num_soilc, filter_soilc,         &
           nitrogenflux_vars, phosphorusflux_vars)

    use CNDecompCascadeConType, only : decomp_cascade_con
    use clm_time_manager, only : get_step_size
#ifndef FLEXIBLE_POOLS
    use clm_varpar, only : i_met_lit, i_cel_lit, i_lig_lit, i_cwd
#endif

    implicit none

    type(bounds_type)        , intent(in) :: bounds
    integer                  , intent(in) :: num_soilc         ! number of soil columns in filter
    integer                  , intent(in) :: filter_soilc(:)   ! filter for soil columns

!    type(atm2lnd_type)       , intent(in) :: clm_a2l
!    type(waterstate_type)    , intent(in) :: waterstate_vars
!    type(soilstate_type)     , intent(in) :: soilstate_vars
!    type(cnstate_type)       , intent(in) :: cnstate_vars
!    type(carbonstate_type)   , intent(inout) :: carbonstate_vars
!    type(carbonflux_type)    , intent(inout) :: carbonflux_vars
!    type(nitrogenstate_type) , intent(inout) :: nitrogenstate_vars
    type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars
    type(phosphorusflux_type)  , intent(inout) :: phosphorusflux_vars
!    type(ch4_type)           , intent(inout) :: ch4_vars
    type(clm_bgc_interface_data_type), intent(in) :: clm_bgc_data

    integer :: fc, c, j
    character(len=256) :: subname = "update_soil_bgc_pf2clm"

    associate ( &
!     initial_cn_ratio             => decomp_cascade_con%initial_cn_ratio             , &
!     decomp_cpools_vr             => carbonstate_vars%decomp_cpools_vr_col           , &
!     decomp_npools_vr             => nitrogenstate_vars%decomp_npools_vr_col         , &
!     decomp_ppools_vr             => phosphorusstate_vars%decomp_ppools_vr_col       , &
!
!     sminn_vr                     => nitrogenstate_vars%sminn_vr_col                 , &
!     smin_no3_vr                  => nitrogenstate_vars%smin_no3_vr_col              , &
!     smin_nh4_vr                  => nitrogenstate_vars%smin_nh4_vr_col              , &
!     smin_nh4sorb_vr              => nitrogenstate_vars%smin_nh4sorb_vr_col          , &
!
!     solutionp_vr    => phosphorusstate_vars%solutionp_vr_col     , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil solution P
!     labilep_vr      => phosphorusstate_vars%labilep_vr_col       , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil labile mineral P
!     secondp_vr      => phosphorusstate_vars%secondp_vr_col       , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil secondary mineralP
!     sminp_vr        => phosphorusstate_vars%sminp_vr_col         , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil mineral P = solutionp + labilep + second
!     occlp_vr        => phosphorusstate_vars%occlp_vr_col         , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil occluded mineral P
!     primp_vr        => phosphorusstate_vars%primp_vr_col           & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil primary mineral P

!     decomp_cpools_sourcesink_vr  => carbonflux_vars%decomp_cpools_sourcesink_col    , &
!     decomp_npools_sourcesink_vr  => nitrogenflux_vars%decomp_npools_sourcesink_col  , &
!     decomp_ppools_sourcesink_vr  => phosphorusflux_vars%decomp_ppools_sourcesink_col, &

     sminn_to_plant_vr            => nitrogenflux_vars%sminn_to_plant_vr_col         , &
     smin_no3_to_plant_vr         => nitrogenflux_vars%smin_no3_to_plant_vr_col      , &
     smin_nh4_to_plant_vr         => nitrogenflux_vars%smin_nh4_to_plant_vr_col      , &
     potential_immob_vr           => nitrogenflux_vars%potential_immob_vr_col        , &
     actual_immob_vr              => nitrogenflux_vars%actual_immob_vr_col           , &
     gross_nmin_vr                => nitrogenflux_vars%gross_nmin_vr_col             , &
     net_nmin_vr                  => nitrogenflux_vars%net_nmin_vr_col               , & ! Output: [real(r8) (:,:)   ]

     sminp_to_plant_vr            => phosphorusflux_vars%sminp_to_plant_vr_col       , &
     potential_immob_p_vr         => phosphorusflux_vars%potential_immob_p_vr_col    , & ! Output: [real(r8) (:,:)   ]
     actual_immob_p_vr            => phosphorusflux_vars%actual_immob_p_vr_col       , &
     gross_pmin_vr                => phosphorusflux_vars%gross_pmin_vr_col           , & ! Output: [real(r8) (:,:)   ]
     net_pmin_vr                  => phosphorusflux_vars%net_pmin_vr_col               & ! Output: [real(r8) (:,:)   ]
     )

    do fc = 1, num_soilc
        c = filter_soilc(fc)
        do j = 1, nlevdecomp

            sminn_to_plant_vr(c,j)      = clm_bgc_data%sminn_to_plant_vr_col(c,j)
            smin_no3_to_plant_vr(c,j)   = clm_bgc_data%smin_no3_to_plant_vr_col(c,j)
            smin_nh4_to_plant_vr(c,j)   = clm_bgc_data%smin_nh4_to_plant_vr_col(c,j)

            potential_immob_vr(c,j)     = clm_bgc_data%potential_immob_vr_col(c,j)
            actual_immob_vr(c,j)        = clm_bgc_data%actual_immob_vr_col(c,j)
            gross_nmin_vr(c,j)          = clm_bgc_data%gross_nmin_vr_col(c,j)
            net_nmin_vr(c,j)            = clm_bgc_data%net_nmin_vr_col(c,j)     !!NOT available in PF

            sminp_to_plant_vr(c,j)      = clm_bgc_data%sminp_to_plant_vr_col(c,j)
            potential_immob_p_vr(c,j)   = clm_bgc_data%potential_immob_p_vr_col(c,j)
            actual_immob_p_vr(c,j)      = clm_bgc_data%actual_immob_p_vr_col(c,j)
            gross_pmin_vr(c,j)          = clm_bgc_data%gross_pmin_vr_col(c,j)
            net_pmin_vr(c,j)            = clm_bgc_data%net_pmin_vr_col(c,j)     !!NOT available in PF

        end do
    end do
    end associate
    end subroutine update_bgc_flux_smin
  !-----------------------------------------------------------------------------
  !
  ! !ROUTINE: update_bgc_gaslosses_pf2clm()
  !
  ! !INTERFACE:
  !
  ! This is a temporary solution to estimate pflotran bgc gaseous emission and transport loss
  ! from their aq. phase states
  ! (due to not yet available in pflotran bgc)
  !
  subroutine update_bgc_flux_gas_pf(clm_bgc_data,  &
     bounds, num_soilc, filter_soilc,           &
     carbonflux_vars, nitrogenflux_vars)

     use clm_time_manager, only : get_step_size, get_nstep

     !
     implicit none

     type(bounds_type) , intent(in)  :: bounds
     integer, intent(in) :: num_soilc       ! number of soil columns in filter
     integer, intent(in) :: filter_soilc(:) ! filter for soil columns

!     type(atm2lnd_type)       , intent(in) :: clm_a2l
!     type(waterstate_type)    , intent(in) :: waterstate_vars
     type(carbonflux_type)    , intent(inout) :: carbonflux_vars
     type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars
     type(clm_bgc_interface_data_type), intent(in) :: clm_bgc_data

     !character(len=256) :: subname = "get_pf_bgc_gaslosses"

     integer  :: fc, c, g, j
!     integer  :: gcount, cellcount
!     real(r8) :: dtime            ! land model time step (sec)
!     integer  :: nstep


!------------------------------------------------------------------------------------
    associate ( &
!     forc_pco2                    => clm_a2l%forc_pco2_grc                  , & ! partial pressure co2 (Pa)
!     forc_pch4                    => clm_a2l%forc_pch4_grc                  , & ! partial pressure ch4 (Pa)
!     forc_pbot                    => clm_a2l%forc_pbot_not_downscaled_grc   , & ! atmospheric pressure (Pa)
!     frac_sno                     => waterstate_vars%frac_sno_eff_col       , & ! fraction of ground covered by snow (0 to 1)
!     frac_h2osfc                  => waterstate_vars%frac_h2osfc_col        , & ! fraction of ground covered by surface water (0 to 1)
!     dz                           => col%dz                                 , & ! soil layer thickness depth (m)
     hr_vr                        => carbonflux_vars%hr_vr_col              , &
     f_co2_soil_vr                => carbonflux_vars%f_co2_soil_vr_col      , &
     f_n2o_soil_vr                => nitrogenflux_vars%f_n2o_soil_vr_col    , &
     f_n2_soil_vr                 => nitrogenflux_vars%f_n2_soil_vr_col     , &
     f_ngas_decomp_vr             => nitrogenflux_vars%f_ngas_decomp_vr_col , &
     f_ngas_nitri_vr              => nitrogenflux_vars%f_ngas_nitri_vr_col  , &
     f_ngas_denit_vr              => nitrogenflux_vars%f_ngas_denit_vr_col    &
     )
! ------------------------------------------------------------------------
    do fc = 1,num_soilc
        c = filter_soilc(fc)
        do j = 1, nlevdecomp
!
              f_co2_soil_vr(c,j)         = clm_bgc_data%f_co2_soil_vr_col(c,j)
              f_n2_soil_vr(c,j)          = clm_bgc_data%f_n2_soil_vr_col(c,j)
              f_n2o_soil_vr(c,j)         = clm_bgc_data%f_n2o_soil_vr_col(c,j)

              hr_vr(c,j)                 = clm_bgc_data%hr_vr_col(c,j)
              f_ngas_decomp_vr(c,j)      = clm_bgc_data%f_ngas_decomp_vr_col(c,j)
              f_ngas_nitri_vr(c,j)       = clm_bgc_data%f_ngas_nitri_vr_col(c,j)
              f_ngas_denit_vr(c,j)       = clm_bgc_data%f_ngas_denit_vr_col(c,j)
       enddo
     enddo ! do c = begc, endc
!
!

    end associate
  end subroutine update_bgc_flux_gas_pf

!-----------------------------------------------------------------------
  subroutine update_bgc_data_pf2clm(clm_bgc_data,bounds,       &
           num_soilc, filter_soilc,                               &
           num_soilp, filter_soilp,                               &
           atm2lnd_vars,                                          &
           waterstate_vars, waterflux_vars,                       &
           soilstate_vars,  temperature_vars, energyflux_vars,    &
           soilhydrology_vars, soil_water_retention_curve,        &
           cnstate_vars, carbonflux_vars, carbonstate_vars,       &
           nitrogenflux_vars, nitrogenstate_vars,                 &
           phosphorusflux_vars, phosphorusstate_vars,             &
           ch4_vars)
    !! USES
    use clm_varctl          , only : use_pflotran, pf_tmode, pf_hmode, pf_cmode

    implicit none

    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:)   ! filter for soil columns
    integer                  , intent(in)    :: num_soilp         ! number of soil patches in filter
    integer                  , intent(in)    :: filter_soilp(:)   ! filter for soil patches
    type(atm2lnd_type)       , intent(in)    :: atm2lnd_vars
    type(waterstate_type)    , intent(inout) :: waterstate_vars
    type(waterflux_type)     , intent(inout) :: waterflux_vars
    type(soilstate_type)     , intent(inout) :: soilstate_vars
    type(temperature_type)   , intent(inout) :: temperature_vars
    type(soilhydrology_type) , intent(inout) :: soilhydrology_vars
    type(energyflux_type)    , intent(inout) :: energyflux_vars
    class(soil_water_retention_curve_type), intent(in) :: soil_water_retention_curve

    type(cnstate_type)          , intent(inout) :: cnstate_vars
    type(carbonflux_type)       , intent(inout) :: carbonflux_vars
    type(carbonstate_type)      , intent(inout) :: carbonstate_vars
    type(nitrogenflux_type)     , intent(inout) :: nitrogenflux_vars
    type(nitrogenstate_type)    , intent(inout) :: nitrogenstate_vars
    type(phosphorusflux_type)   , intent(inout) :: phosphorusflux_vars
    type(phosphorusstate_type)  , intent(inout) :: phosphorusstate_vars
    type(ch4_type)              , intent(inout) :: ch4_vars
    type(clm_bgc_interface_data_type), intent(in) :: clm_bgc_data

    !-----------------------------------------------------------------------

    character(len=256) :: subname = "get_clm_bgc_data"


    if (pf_cmode) then
        !! bgc_state_decomp is updated in CLM
        !! by passing bgc_flux_decomp_sourcesink into CNSoilLittVertTransp
        call update_bgc_flux_decomp_sourcesink(clm_bgc_data,    &
                    bounds, num_soilc, filter_soilc,            &
                    carbonflux_vars, nitrogenflux_vars,         &
                    phosphorusflux_vars)

        call update_bgc_state_smin(clm_bgc_data,                &
                    bounds, num_soilc, filter_soilc,            &
                    nitrogenstate_vars, phosphorusstate_vars)

        call update_bgc_flux_smin(clm_bgc_data,                 &
                    bounds, num_soilc, filter_soilc,            &
                    nitrogenflux_vars, phosphorusflux_vars)

        call update_bgc_flux_gas_pf(clm_bgc_data,               &
                    bounds, num_soilc, filter_soilc,            &
                    carbonflux_vars, nitrogenflux_vars)

    end if

    if (pf_tmode) then
        call update_soil_temperature(clm_bgc_data,      &
                   bounds, num_soilc, filter_soilc,     &
                   temperature_vars)
    end if

    if (pf_hmode) then
        call update_soil_moisture(clm_bgc_data,         &
                   bounds, num_soilc, filter_soilc,     &
                   waterstate_vars)
    end if

  end subroutine update_bgc_data_pf2clm


  !-----------------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: update_bcflow_pf2clm
  !
  ! !INTERFACE:
!  subroutine update_bcflow_pf2clm(         &
!       bounds, num_soilc, filter_soilc,    &
!       clm_a2l, ces_vars, cef_vars, cws_vars, cwf_vars)
!  !
!  ! !DESCRIPTION:
!  ! update qflx_surf, qflx_infl from PF's
!  ! 'mass_balance' retrieving from PFLOTRAN
!  !
!  ! !USES:
!    use clm_time_manager, only : get_step_size
!
!    type(bounds_type), intent(in) :: bounds
!    integer, intent(in) :: num_soilc        ! number of column soil points in column filter
!    integer, intent(in) :: filter_soilc(:)  ! column filter for soil points
!
!    type(atm2lnd_type)       , intent(in) :: clm_a2l
!    type(temperature_type)   , intent(in) :: ces_vars
!    type(energyflux_type)    , intent(in) :: cef_vars
!    type(waterstate_type)    , intent(in) :: cws_vars
!    type(waterflux_type)     , intent(inout) :: cwf_vars
!
!  ! !LOCAL VARIABLES:
!#include "finclude/petscsys.h"
!#include "finclude/petscvec.h"
!#include "finclude/petscvec.h90"
!
!    integer  :: fc, c, g, gcount           ! indices
!    real(r8) :: area                       ! top face area
!    real(r8) :: dtime                      ! land model time step (sec)
!    real(r8) :: dew
!    real(r8) :: qflx_evap(bounds%begc:bounds%endc)            ! soil surface evaporation (mmH2O/s)
!
!    PetscScalar, pointer :: area_clm_loc(:)
!    PetscScalar, pointer :: qinfl_subsurf_clm_loc(:)        ! kgH2O/time-step
!    PetscScalar, pointer :: qsurf_subsurf_clm_loc(:)        ! kgH2O/time-step
!    PetscScalar, pointer :: qflux_subbase_clm_loc(:)        ! kgH2O/time-step
!    PetscErrorCode :: ierr
!
!    character(len=32) :: subname = 'update_bcflow_pf2clm'  ! subroutine name
!
!    !-----------------------------------------------------------------------
!
!    associate(&
!    clandunit         =>    col%landunit          , & ! column's landunit
!    ltype             =>    lun%itype             , & ! landunit type
!    zi                =>    col%zi                , & ! Input: (:,:) soil layer interface depth (m)
!    dz                =>    col%dz                , & ! Input: (:,:) soil layer thickness (m)
!    !
!    forc_pbot         =>    clm_a2l%forc_pbot_not_downscaled_grc , & ! Input:  [real(r8) (:)]  atmospheric pressure (Pa)
!    t_grnd            =>    ces_vars%t_grnd_col                  , & ! Input:  [real(r8) (:)]  ground surface temperature [K]
!    htvp              =>    cef_vars%htvp_col                    , & ! Input:  [real(r8) (:)]  latent heat of vapor of water (or sublimation) [j/kg]
!    !
!    frac_sno          =>    cws_vars%frac_sno_eff_col            , & ! Input: fraction of ground covered by snow (0 to 1)
!    frac_h2osfc       =>    cws_vars%frac_h2osfc_col             , & ! Input: fraction of ground covered by surface water (0 to 1)
!    !
!    qflx_top_soil     =>    cwf_vars%qflx_top_soil_col           , & ! Input:  [real(r8) (:)]  net water input into soil from top (mm/s)
!    qflx_ev_h2osfc    =>    cwf_vars%qflx_ev_h2osfc_col          , & ! Input: column-level evaporation flux from h2osfc (W/m2) [+ to atm] : checking unit
!    qflx_evap_soil    =>    cwf_vars%qflx_evap_soi_col           , & ! Input: column-level soil evaporation (mm H2O/s) (+ = to atm)
!    qflx_subl_snow    =>    cwf_vars%qflx_sub_snow_col           , & ! Input: column-level evaporation flux from snow (mm H2O/s) [+ to atm]
!    qflx_surf         =>    cwf_vars%qflx_surf_col               , & ! Output: [real(r8) (:)]  surface runoff (mm H2O /s)
!    qflx_infl         =>    cwf_vars%qflx_infl_col               , & ! Output: [real(r8) (:)]  soil infiltration (mm H2O /s)
!    qflx_drain        =>    cwf_vars%qflx_drain_col                & ! Output: [real(r8) (:)]  sub-surface runoff (drainage) (mm H2O /s)
!    )
!
!    dtime = get_step_size()
!
!    ! the following was actually duplicated from 'get_clm_bcwflx' to calculate total water evap from 'qflx_topsoil'
!    ! in order to get potential infiltration from CLM (not yet run-off)
!    do fc = 1, num_soilc
!       c = filter_soilc(fc)
!       if (ltype(clandunit(c)) == istsoil .or. ltype(clandunit(c))==istcrop) then
!           ! not sure if using 'qflx_evap_soi_col' as a whole ground-surface better than individual surfaces, i.e. 'qflx_ev_snow/soi/h2osfc'
!           ! all of those 4 variables are calculated in 'BareGroundFluxesMod', 'CanopyFluxesMod', and then adjusted in 'Biogeophysics2' after soil temperature call
!           ! note that: all 4 variables could be negative (i.e., dew formation on ground)
!           qflx_evap(c)=(1.0_r8 - frac_sno(c) - frac_h2osfc(c))*qflx_evap_soil(c) + &
!                          frac_h2osfc(c)*qflx_ev_h2osfc(c)/htvp(c)
!                          !frac_sno(c)*qflx_evap_snow(c)   ! snow-covered area should be excluded (see SoilHydrologyMod:: infiltration)
!       else
!          ! for other types of landunits
!          qflx_evap(c) = (1.0_r8 - frac_sno(c))*qflx_evap_soil(c)
!       end if
!
!       if (t_grnd(c) <= tfrz) qflx_evap(c) = max(0._r8, qflx_evap(c))   ! frozen ground, no dew contribution to subsurface infiltration
!    end do
!
!
!    ! from PF==>CLM
!    call VecGetArrayF90(clm_pf_idata%area_top_face_clms, area_clm_loc, ierr)
!    CHKERRQ(ierr)
!
!    call VecGetArrayF90(clm_pf_idata%qinfl_subsurf_clms,qinfl_subsurf_clm_loc,ierr)
!    CHKERRQ(ierr)
!    call VecGetArrayF90(clm_pf_idata%qsurf_subsurf_clms,qsurf_subsurf_clm_loc,ierr)
!    CHKERRQ(ierr)
!    call VecGetArrayF90(clm_pf_idata%qflux_subbase_clms,qflux_subbase_clm_loc,ierr)
!    CHKERRQ(ierr)
!
!    do fc = 1, num_soilc
!
!      c = filter_soilc(fc)
!      g = col%gridcell(c)
!      gcount = g - bounds%begg
!
!      dew = 0._r8
!      if (t_grnd(c) > tfrz) then   ! frozen ground, no dew contribution to subsurface infiltration
!        if (ltype(clandunit(c)) == istsoil .or. ltype(clandunit(c))==istcrop) then
!           dew = -min(0._r8, (1.0_r8 - frac_sno(c) - frac_h2osfc(c))*qflx_evap_soil(c) + &
!                             frac_h2osfc(c)*qflx_ev_h2osfc(c)/htvp(c) )
!        else
!           dew = -min(0._r8, qflx_evap_soil(c))
!        end if
!      endif
!
!      !'from PF: qinfl_subsurf_clm_loc: positive - in, negative - out
!      area = area_clm_loc(gcount*clm_pf_idata%nzclm_mapped+1)
!      qflx_infl(c) = qinfl_subsurf_clm_loc(gcount+1)  &
!                       /dtime/(area*denh2o*1.e-3)     ! kgH2O/time-step ==> mmH2O/sec
!
!      qflx_surf(c) = dew + max(0._r8, qflx_top_soil(c)-qflx_evap(c)) - qflx_infl(c)
!      qflx_surf(c) = max(0._r8, qflx_surf(c))
!
!      !'from PF: qflux_subbase_clm_loc: positive - in, negative - out)
!      area = area_clm_loc((gcount+1)*clm_pf_idata%nzclm_mapped)                              ! note: this 'area_clm_loc' is in 3-D for all subsurface domain
!      qflx_drain(c) = -qflux_subbase_clm_loc(gcount+1)  &
!                       /dtime/(area*denh2o*1.e-3)     ! kgH2O/time-step ==> mmH2O/sec (+ drainage, - upward-in)
!
!    end do
!    call VecRestoreArrayF90(clm_pf_idata%area_top_face_clms, area_clm_loc, ierr)
!    CHKERRQ(ierr)
!    call VecRestoreArrayF90(clm_pf_idata%qinfl_subsurf_clms,qinfl_subsurf_clm_loc,ierr)
!    CHKERRQ(ierr)
!    call VecRestoreArrayF90(clm_pf_idata%qsurf_subsurf_clms,qsurf_subsurf_clm_loc,ierr)
!    CHKERRQ(ierr)
!    call VecRestoreArrayF90(clm_pf_idata%qflux_subbase_clms,qflux_subbase_clm_loc,ierr)
!    CHKERRQ(ierr)
!
!    end associate
!  end subroutine update_bcflow_pf2clm
!
!#endif
!----------------------------------------------------------------------------------------------
! END of CLM-PFLOTRAN bgc coupling interface
!----------------------------------------------------------------------------------------------



!----------------------------------------------------------------------------------------------
! BEG of CLM-bgc through interface
!----------------------------------------------------------------------------------------------
  ! !INTERFACE:
  subroutine clm_bgc_run(bounds,                                        &
            num_soilc, filter_soilc, num_soilp, filter_soilp,           &
            photosyns_vars, canopystate_vars,                           &
            soilstate_vars, temperature_vars, waterstate_vars,          &
            cnstate_vars, ch4_vars,                                     &
            carbonstate_vars, carbonflux_vars,                          &
            c13_carbonflux_vars, c14_carbonflux_vars,                   &
            nitrogenstate_vars, nitrogenflux_vars, crop_vars,           &
            phosphorusstate_vars,phosphorusflux_vars)

    !! USES:
    use CNDecompMod          , only: CNDecompAlloc1

    !! ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_soilc          ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:)    ! filter for soil columns
    integer                  , intent(in)    :: num_soilp          ! number of soil patches in filter
    integer                  , intent(in)    :: filter_soilp(:)    ! filter for soil patches
    type(photosyns_type)     , intent(in)    :: photosyns_vars
    type(canopystate_type)   , intent(in)    :: canopystate_vars
    type(soilstate_type)     , intent(in)    :: soilstate_vars
    type(temperature_type)   , intent(in)    :: temperature_vars
    type(waterstate_type)    , intent(in)    :: waterstate_vars
    type(cnstate_type)       , intent(inout) :: cnstate_vars
    type(ch4_type)           , intent(in)    :: ch4_vars
    type(carbonstate_type)   , intent(inout) :: carbonstate_vars
    type(carbonflux_type)    , intent(inout) :: carbonflux_vars
    type(carbonflux_type)    , intent(inout) :: c13_carbonflux_vars
    type(carbonflux_type)    , intent(inout) :: c14_carbonflux_vars
    type(nitrogenstate_type) , intent(inout) :: nitrogenstate_vars
    type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars
    type(crop_type)          , intent(in)    :: crop_vars
    type(phosphorusstate_type) , intent(inout) :: phosphorusstate_vars
    type(phosphorusflux_type)  , intent(inout) :: phosphorusflux_vars


    !! wgs: TODO: data pass from ALM to INTERFACE to ALM
    call CNDecompAlloc1(bounds,                                         &
            num_soilc, filter_soilc, num_soilp, filter_soilp,           &
            photosyns_vars, canopystate_vars,                           &
            soilstate_vars, temperature_vars, waterstate_vars,          &
            cnstate_vars, ch4_vars,                                     &
            carbonstate_vars, carbonflux_vars,                          &
            c13_carbonflux_vars, c14_carbonflux_vars,                   &
            nitrogenstate_vars, nitrogenflux_vars, crop_vars,           &
            phosphorusstate_vars,phosphorusflux_vars)

  end subroutine clm_bgc_run
!----------------------------------------------------------------------------------------------
! END of CLM-bgc through interface
!----------------------------------------------------------------------------------------------

end module clm_bgc_interfaceMod

