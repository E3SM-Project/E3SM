module clm_pflotran_interfaceMod
!#define FLEXIBLE_POOLS
!----------------------------------------------------------------------------------------------
! CLM-PFLOTRAN soil bgc coupling interface
! authors: Fengming YUAN1, Gautam Bishit1,2, and Guoping Tang1
!
!          1.Climate Change Science Institute & Environmental Science Division
!          Oak Ridge National Laboratory
!
!          2.Lawrence Berkley National Laboratory
!
! date: 2012 - 2015
! modified (based on clm_bgc_interface): 8/28/2015, wgs
!----------------------------------------------------------------------------------------------

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  !BOP
  !
  ! !MODULE: clm_pflotran_interfaceMod
  !
  ! !DESCRIPTION:
  ! Performs
  !
  ! !USES:

  ! clm g/l/c/p constants
  use shr_log_mod         , only : errMsg => shr_log_errMsg
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use GridcellType        , only : grc_pp
  use LandunitType        , only : lun_pp
  use ColumnType          , only : col_pp 
  use PatchType           , only : pft_pp

  use decompMod           , only : bounds_type

  ! (dummy) variable definitions
  !! (wgs,8/28/2015): variables below are replaced by clm_bgc_data

!  use atm2lndType         , only : atm2lnd_type
!  use SoilStateType       , only : soilstate_type
!  use WaterStateType      , only : waterstate_type
!  use WaterFluxType       , only : waterflux_type
!  use SoilHydrologyType   , only : soilhydrology_type
!  use TemperatureType     , only : temperature_type
!  use EnergyFluxType      , only : energyflux_type
!  use SoilWaterRetentionCurveMod, only : soil_water_retention_curve_type
!
!  use CNStateType         , only : cnstate_type
!  use CNCarbonFluxType    , only : carbonflux_type
!  use CNCarbonStateType   , only : carbonstate_type
!  use CNNitrogenFluxType  , only : nitrogenflux_type
!  use CNNitrogenStateType , only : nitrogenstate_type
!
!  use ch4Mod              , only : ch4_type
!
!  use PhotosynthesisType     , only : photosyns_type
!  use cropType               , only : crop_type
!  use CanopyStateType        , only : canopystate_type
!  use PhosphorusStateType    , only : phosphorusstate_type
!  use PhosphorusFluxType     , only : phosphorusflux_type


  ! PFLOTRAN thc module switchs for coupling with CLM45-CN
  use clm_varctl          , only : use_pflotran, pf_tmode, pf_hmode, pf_cmode
  use clm_varctl          , only : pf_surfaceflow, pf_frzmode
  use clm_varctl          , only : initth_pf2clm

  ! most used constants in this module
  use clm_varpar          , only : nlevsoi, nlevgrnd, nlevdecomp, nlevdecomp_full
  use clm_varpar          , only : ndecomp_pools
  use clm_varpar          , only : max_patch_per_col
  use clm_varcon          , only : denh2o, denice, tfrz, dzsoi_decomp
  use landunit_varcon     , only : istsoil, istcrop

  ! misc.
  use abortutils          , only : endrun

!  use clm_bgc_interfaceMod
  use clm_bgc_interface_data, only : clm_bgc_interface_data_type


#ifdef CLM_PFLOTRAN
  use clm_pflotran_interface_data
  use pflotran_clm_main_module
  use pflotran_clm_setmapping_module
#endif

  ! !PUBLIC TYPES:
  implicit none

  save

  private    ! By default everything is private

#ifdef CLM_PFLOTRAN
  type(pflotran_model_type), pointer, public :: pflotran_m
#endif
  !
  character(len=256), private:: pflotran_prefix = ''
  character(len=32), private :: restart_stamp = ''

  real(r8), parameter :: rgas = 8.3144621d0                 ! m3 Pa K-1 mol-1

  ! !PUBLIC MEMBER FUNCTIONS:
  public :: clm_pf_readnl

  ! wrappers around '#ifdef CLM_PFLOTRAN .... #endif' block statements to maintain sane runtime behavior
  ! when pflotran is not available.
  public :: clm_pf_interface_init
  public :: clm_pf_set_restart_stamp
  public :: clm_pf_run
  public :: clm_pf_write_restart
  public :: clm_pf_finalize


#ifdef CLM_PFLOTRAN
  ! private work functions that truely require '#ifdef CLM_PFLOTRAN .... #endif'
  private :: interface_init
  private :: pflotran_run_onestep
  private :: pflotran_write_checkpoint
  private :: pflotran_finalize
  private :: get_clm_soil_properties
  private :: get_clm_soil_th
  private :: get_clm_iceadj_porosity
  private :: get_clm_bcwflx
  private :: get_clm_bceflx
  private :: get_clm_bgc_conc
  private :: get_clm_bgc_rate
  private :: update_soil_temperature_pf2clm
  private :: update_soil_moisture_pf2clm
  private :: update_soil_bgc_pf2clm
  private :: update_bgc_gaslosses_pf2clm
! private :: update_bcflow_pf2clm !!comment out currently for bgc-only

#endif

contains

!-----------------------------------------------------------------------
!
! public interface functions allowing runtime behavior regardless of
! whether pflotran is compiled in.
!
!-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: clm_pf_readnl
  !
  ! !INTERFACE:
  subroutine clm_pf_readnl( NLFilename )
  !
  ! !DESCRIPTION:
  ! Read namelist for clm-pflotran interface
  !
  ! !USES:
    use clm_varctl    , only : iulog
    use spmdMod       , only : masterproc, mpicom
    use fileutils     , only : getavu, relavu, opnfil
    use clm_nlUtilsMod, only : find_nlgroup_name
    use shr_nl_mod    , only : shr_nl_find_group_name
    use shr_mpi_mod   , only : shr_mpi_bcast

    implicit none

  ! !ARGUMENTS:
    character(len=*), intent(IN) :: NLFilename ! Namelist filename
  ! !LOCAL VARIABLES:
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file
    character(len=32) :: subname = 'clm_pf_readnl'  ! subroutine name
  !EOP
  !-----------------------------------------------------------------------
    namelist / clm_pflotran_inparm / pflotran_prefix

    ! ----------------------------------------------------------------------
    ! Read namelist from standard namelist file.
    ! ----------------------------------------------------------------------

    if ( masterproc )then

       unitn = getavu()
       write(iulog,*) 'Read in clm-pflotran namelist'
       call opnfil (NLFilename, unitn, 'F')
       call shr_nl_find_group_name(unitn, 'clm_pflotran_inparm', status=ierr)
       if (ierr == 0) then
          read(unitn, clm_pflotran_inparm, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg=subname //':: ERROR: reading clm_pflotran_inparm namelist.'//&
                         errMsg(__FILE__, __LINE__))
          end if
       end if
       call relavu( unitn )
       write(iulog, '(/, A)') " clm-pflotran namelist:"
       write(iulog, '(A, " : ", A,/)') "   pflotran_prefix", trim(pflotran_prefix)
    end if

    ! Broadcast namelist variables read in
    call shr_mpi_bcast(pflotran_prefix, mpicom)

  end subroutine clm_pf_readnl

  !-----------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: clm_pf_set_restart_stamp
  !
  ! !INTERFACE:
  subroutine clm_pf_set_restart_stamp(clm_restart_filename)
  !
  ! !DESCRIPTION: Set the pflotran restart date stamp. Note we do NOT
  ! restart here, that gets handled by pflotran's internal
  ! initialization during interface_init_clm_pf()
  !
  ! !USES:
  ! !ARGUMENTS:
    character(len=256), intent(in) :: clm_restart_filename
  ! !LOCAL VARIABLES:
    integer :: name_length, start_pos, end_pos
    character(len=32) :: clm_stamp
  !EOP
  !-----------------------------------------------------------------------

    ! clm restart file name is of the form:
    !     ${CASE_NAME}.clm2.r.YYYY-MM-DD-SSSSS.nc
    ! we need to extract the: YYYY-MM-DD-SSSSS
    write(*, '("clm-pf : clm restart file name : ", A/)') trim(clm_restart_filename)
    name_length = len(trim(clm_restart_filename))
    start_pos = name_length - 18
    end_pos = name_length - 3
    clm_stamp = clm_restart_filename(start_pos : end_pos)
    write(*, '("clm-pf : clm date stamp : ", A/)') trim(clm_stamp)
    restart_stamp = clm_stamp

  end subroutine clm_pf_set_restart_stamp


  !-----------------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: pflotran_not_available
  !
  ! !INTERFACE:
  subroutine pflotran_not_available(subname)
  !
  ! !DESCRIPTION:
  ! Print an error message and abort.
  !
  ! !USES:

  ! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: subname
  ! !LOCAL VARIABLES:
  !EOP
  !-----------------------------------------------------------------------
    call endrun(trim(subname) // ": ERROR: CLM-PFLOTRAN interface has not been compiled " // &
         "into this version of CLM.")
  end subroutine pflotran_not_available


!******************************************************************************************!
!
! public interface function wrappers
!
!------------------------------------------------------------------------------------------!

  !-----------------------------------------------------------------------------
  subroutine clm_pf_interface_init(bounds)

    implicit none

    type(bounds_type), intent(in) :: bounds  ! bounds

    character(len=256) :: subname = "clm_pf_interface_init()"

#ifdef CLM_PFLOTRAN
    call interface_init(bounds)
#else
    call pflotran_not_available(subname)
#endif
  end subroutine clm_pf_interface_init

  !--------------------------------------------------------------------------------------------

    subroutine clm_pf_run(clm_bgc_data,bounds,       &
           num_soilc, filter_soilc)

    implicit none

    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:)   ! filter for soil columns
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
    type(clm_bgc_interface_data_type), intent(inout) :: clm_bgc_data

    !-----------------------------------------------------------------------

    character(len=256) :: subname = "clm_pf_run"

#ifdef CLM_PFLOTRAN
    call pflotran_run_onestep(clm_bgc_data,bounds,  &
           num_soilc, filter_soilc)
#else
    call pflotran_not_available(subname)
#endif
  end subroutine clm_pf_run

  !-----------------------------------------------------------------------------
  subroutine clm_pf_write_restart(date_stamp)

    implicit none
    character(len=*), intent(in) :: date_stamp

    character(len=32) :: subname = "clm_pf_write_restart"

#ifdef CLM_PFLOTRAN
    call pflotran_write_checkpoint(date_stamp)
#else
    call pflotran_not_available(subname)
#endif
  end subroutine clm_pf_write_restart


  !-----------------------------------------------------------------------------
  !BOP
  !
  ! !ROUTINE: clm_pf_finalize
  !
  ! !INTERFACE:
  subroutine clm_pf_finalize()

    implicit none
    character(len=256) :: subname = "clm_pf_finalize"

#ifdef CLM_PFLOTRAN
    call pflotran_finalize()
#else
    call pflotran_not_available(subname)
#endif
  end subroutine clm_pf_finalize




#ifdef CLM_PFLOTRAN

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
  subroutine interface_init(bounds)
    !
    ! !DESCRIPTION:
    ! initialize the pflotran iterface
    !
    ! !USES:
    use clm_varctl      , only : iulog
    use decompMod       , only : get_proc_global, ldecomp
    use spmdMod         , only : mpicom, masterproc
    use domainMod       , only : ldomain

    use CNDecompCascadeConType , only : decomp_cascade_con

    use abortutils      , only : endrun

    use clm_time_manager, only : nsstep, nestep

    ! pflotran
    use Option_module!, only : printErrMsg
    use Simulation_Base_class, only : simulation_base_type
    use Simulation_Subsurface_class, only : subsurface_simulation_type
    use Realization_class, only : realization_type
    use PFLOTRAN_Constants_module
    !
    ! !ARGUMENTS:

    implicit none

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscviewer.h"

    !
    ! !REVISION HISTORY:
    ! Created by Gautam Bisht
    ! Revised by Fengming Yuan, CCSI-ORNL
    !
    !EOP
    !
    ! LOCAL VARAIBLES:

    type(bounds_type), intent(in) :: bounds  ! bounds

    integer  :: numg             ! total number of gridcells across all processors
    integer  :: numl             ! total number of landunits across all processors
    integer  :: numc             ! total number of columns across all processors
    integer  :: nump             ! total number of pfts across all processors
    integer  :: g,l,c,j          ! indices
    integer  :: gcount, cellcount

    integer, pointer :: clm_cell_ids_nindex(:)
    integer, pointer :: clm_2dtop_cell_ids_nindex(:)
    integer, pointer :: clm_2dbot_cell_ids_nindex(:)
    integer :: clm_cell_npts
    integer :: clm_2dtop_npts
    integer :: clm_2dbot_npts

    !type(pflotran_model_type), pointer:: pflotran_m
    class(realization_type), pointer  :: realization
    type(option_type), pointer :: option
    PetscErrorCode :: ierr

    character(len= 32) :: subname = 'interface_init' ! subroutine name

    associate( &
         ! Assign local pointers to derived subtypes components (gridcell-level)
         latdeg     =>  grc_pp%latdeg     , & !  [real(r8) (:)]  latitude (degree)
         londeg     =>  grc_pp%londeg     , & !  [real(r8) (:)]  longitude (degree)
         area       =>  grc_pp%area       , & !  [real(r8) (:)]  total land area per gridcell (km^2)
         gindex     =>  grc_pp%gindex     , & !  [real(r8) (:)]  longitude (degree)
         ! Assign local pointers to derived subtypes components (landunit-level)
         ltype      =>  lun_pp%itype      , & !  [integer (:)]  landunit type index
         ! Assign local pointer to derived subtypes components (column-level)
         cgridcell  =>  col_pp%gridcell   , & !  [integer (:)]  gridcell index of column
         clandunit  =>  col_pp%landunit   , & !  [integer (:)]  landunit index of column
         cwtgcell   =>  col_pp%wtgcell    , & !  [real(r8) (:)]  weight (relative to gridcell
         ctype      =>  col_pp%itype      , & !  [integer (:)]  column type index
         topo       =>  col_pp%glc_topo   , &  ! surface elevation (m)
         micro_sigma=>  col_pp%micro_sigma, &  ! microtopography pdf sigma (m)
         slope      =>  col_pp%topo_slope , &  ! gridcell topographic slope
         topo_std   =>  col_pp%topo_std     &  ! gridcell elevation standard deviation
         )

    !------------------------------------------------------------------------

    do c = bounds%begc, bounds%endc
      l = col_pp%landunit(c)
      if (.not.(ltype(l)==istsoil .or. ltype(l)==istcrop) .and. col_pp%active(c)) then
        write (iulog,*) 'WARNING: Land Unit type of active Columns of non-soil/crop type found within the domain'
        write (iulog,*) 'CLM-CN -- PFLOTRAN does not support this land unit presently'
        write (iulog,*) 'So, DEactive the column: ', c

        col_pp%active(c) = .false.

      endif

    enddo ! do c = bounds%begc, bounds%endc

    if (masterproc) then
      write(iulog,*) '%%-----------------------------------------------------%%'
      write(iulog,*) '%%               clm_pf_interface_init                 %%'
      write(iulog,*) '%%-----------------------------------------------------%%'
      write(iulog,*) ' '
    endif


    ! Determine necessary indices
    call get_proc_global(numg, numl, numc, nump)

    !------------------------------------------------------------------------
    allocate(pflotran_m)

    ! Create PFLOTRAN model
    call pflotranModelCreate(mpicom, pflotran_prefix, pflotran_m)
    option => pflotran_m%option

    call pflotranModelSetupMappingfiles(pflotran_m)

    ! initialize pflotran model
    call pflotranModelStepperRunInit(pflotran_m)

    ! PFLOTRAN ck file NOT works well when coupling with CLM.
    ! So it's off now and restart PF from CLM.
    !restart_stamp = ""
    !call pflotranModelSetupRestart(pflotran_m, restart_stamp)
    initth_pf2clm = .false.
    !if (restart_stamp .ne. "") then
    !   initth_pf2clm = .true.
    !endif

    ! Initialize PETSc vector for data transfer between CLM and PFLOTRAN
    call CLMPFLOTRANIDataInit()

    clm_pf_idata%nzclm_mapped  = nlevsoi        ! the soil layer no. mapped btw CLM and PF for data-passing

    select type (simulation => pflotran_m%simulation)
      class is (subsurface_simulation_type)
         realization => simulation%realization
      class default
         pflotran_m%option%io_buffer = "This version of clm-pflotran only works with subsurface simulations."
         write(*, '(/A/)') pflotran_m%option%io_buffer
         call printErrMsg(pflotran_m%option)
    end select

    if (pflotran_m%option%iflowmode == TH_MODE .or. pflotran_m%option%iflowmode == Richards_MODE) then
       pflotran_m%option%io_buffer = "This version of clm-pflotran Richards_mode or TH_mode has NOT yet well tested."
       write(*, '(/A/)') pflotran_m%option%io_buffer
       call printErrMsg(pflotran_m%option)
    endif

    if(pflotran_m%option%nsurfflowdof > 0) then
       pflotran_m%option%io_buffer = "This version of clm-pflotran DOES NOT work with PF Surface simulation."
       write(*, '(/A/)') pflotran_m%option%io_buffer
       call printErrMsg(pflotran_m%option)
    endif
    pf_surfaceflow = .false.

    !------------------------------------------------

    ! Compute number of cells in CLM domain.
    clm_cell_npts  = (bounds%endg - bounds%begg + 1)*clm_pf_idata%nzclm_mapped
    clm_2dtop_npts = (bounds%endg - bounds%begg + 1)
    clm_2dbot_npts = (bounds%endg - bounds%begg + 1)
    allocate(clm_cell_ids_nindex( 1:clm_cell_npts))
    allocate(clm_2dtop_cell_ids_nindex(1:clm_2dtop_npts))
    allocate(clm_2dbot_cell_ids_nindex(1:clm_2dbot_npts))

    ! Save cell IDs of CLM grid
    cellcount = 0
    do g = bounds%begg, bounds%endg
       gcount = g - bounds%begg + 1

       do j = 1, clm_pf_idata%nzclm_mapped
          cellcount = cellcount + 1
          !clm_cell_ids_nindex(cellcount) = (ldecomp%gdc2glo(g)-1)*nlevsoi + j - 1
          ! maintain CLM each processor's grids to match up with PF (this appears saving a lot of data communications)
          clm_cell_ids_nindex(cellcount) = (g-1)*clm_pf_idata%nzclm_mapped + j - 1
       enddo

       clm_2dtop_cell_ids_nindex(gcount) = (g-1)*clm_pf_idata%nzclm_mapped

       clm_2dbot_cell_ids_nindex(gcount) = g*clm_pf_idata%nzclm_mapped-1
    enddo

    ! CLM: 3-D Subsurface domain (local and ghosted cells)
    clm_pf_idata%nlclm_sub = clm_cell_npts
    clm_pf_idata%ngclm_sub = clm_cell_npts

    ! CLM: Surface/Bottom cells of subsurface domain (local and ghosted cells)
    clm_pf_idata%nlclm_2dtop = (bounds%endg - bounds%begg + 1)
    clm_pf_idata%ngclm_2dtop = (bounds%endg - bounds%begg + 1)

    ! CLM: bottom face of subsurface domain
    clm_pf_idata%nlclm_2dbot = (bounds%endg - bounds%begg + 1)
    clm_pf_idata%ngclm_2dbot = (bounds%endg - bounds%begg + 1)

    ! PFLOTRAN: 3-D Subsurface domain (local and ghosted cells)
    clm_pf_idata%nlpf_sub = realization%patch%grid%nlmax
    clm_pf_idata%ngpf_sub = realization%patch%grid%ngmax

    ! PFLOTRAN: surface of subsurface domain (local and ghosted cells)
    !if (pflotran_m%option%iflowmode == RICHARDS_MODE) then
    !  clm_pf_idata%nlpf_2dtop = pflotranModelNSurfCells3DDomain(pflotran_m)   ! this will be overwritten in Mapping below
    !  clm_pf_idata%ngpf_2dtop = pflotranModelNSurfCells3DDomain(pflotran_m)
    !endif

    ! For CLM/PF: ground surface NOT defined, so need to set the following to zero.
    clm_pf_idata%nlclm_srf = 0
    clm_pf_idata%ngclm_srf = 0
    clm_pf_idata%nlpf_srf = 0
    clm_pf_idata%ngpf_srf = 0

    ! Initialize maps for transferring data between CLM and PFLOTRAN.
    call pflotranModelInitMapping(pflotran_m, clm_cell_ids_nindex, &
                                  clm_cell_npts, CLM_SUB_TO_PF_SUB)
    call pflotranModelInitMapping(pflotran_m, clm_cell_ids_nindex, &
                                  clm_cell_npts, PF_SUB_TO_CLM_SUB)

    !
    if (pflotran_m%option%iflowmode == RICHARDS_MODE .or. &
        pflotran_m%option%iflowmode == TH_MODE) then
       call pflotranModelInitMapping(pflotran_m, clm_2dtop_cell_ids_nindex,   &
                                      clm_2dtop_npts, CLM_2DTOP_TO_PF_2DTOP)
       call pflotranModelInitMapping(pflotran_m, clm_2dtop_cell_ids_nindex,   &
                                      clm_2dtop_npts, PF_2DTOP_TO_CLM_2DTOP)

       call pflotranModelInitMapping(pflotran_m, clm_2dbot_cell_ids_nindex, &
                                     clm_2dbot_npts, CLM_2DBOT_TO_PF_2DBOT)
       call pflotranModelInitMapping(pflotran_m, clm_2dbot_cell_ids_nindex, &
                                     clm_2dbot_npts, PF_2DBOT_TO_CLM_2DBOT)
    endif

    ! the CLM-CN/BGC decomposing pool size
    clm_pf_idata%ndecomp_pools     = ndecomp_pools

    ! Allocate vectors for data transfer between CLM and PFLOTRAN.
    call CLMPFLOTRANIDataCreateVec(MPI_COMM_WORLD)

    ! the CLM-CN/BGC decomposing pool if varying cn ratios, and names.
    clm_pf_idata%floating_cn_ratio = decomp_cascade_con%floating_cn_ratio_decomp_pools(1:ndecomp_pools)
    clm_pf_idata%decomp_pool_name  = decomp_cascade_con%decomp_pool_name_history(1:ndecomp_pools)

    ! if BGC is on
    if(pflotran_m%option%ntrandof > 0) then
      call pflotranModelGetRTspecies(pflotran_m)
    endif

    ! Get pflotran top surface area
    call pflotranModelGetTopFaceArea(pflotran_m)

    deallocate(clm_cell_ids_nindex)
    deallocate(clm_2dtop_cell_ids_nindex)
    deallocate(clm_2dbot_cell_ids_nindex)

!-------------------------------------------------------------------------------------
    ! coupled module controls betweeen PFLOTRAN and CLM45 (F.-M. Yuan, Aug. 2013)
    if(pflotran_m%option%iflowmode==RICHARDS_MODE) then
      pf_hmode = .true.
      pf_tmode = .false.
      pf_frzmode = .false.

    elseif(pflotran_m%option%iflowmode==TH_MODE) then
      pf_hmode = .true.
      pf_tmode = .true.
      if (pflotran_m%option%use_th_freezing) then
        pf_frzmode = .true.
      else
        pf_frzmode = .false.
      endif

    endif

    if(pflotran_m%option%ntrandof.gt.0) then
      pf_cmode = .true.        ! initialized as '.false.' in clm initialization
    endif

!-------------------------------------------------------------------------------------

    end associate
  end subroutine interface_init

  !-----------------------------------------------------------------------------
  !
  ! !SUBROUTINE: pflotran_run_onestep
  !
  ! !INTERFACE:

  subroutine pflotran_run_onestep(clm_bgc_data,bounds,            &
           num_soilc, filter_soilc)
  !
  ! !DESCRIPTION:
  !
  !  F.-M. YUAN: based on Gautam's 'step_th_clm_pf',
  !              'chemistry' (PF_CMODE) added (Sept. 6, 2013)
  !
  ! !USES:


!    use clm_pflotran_interface_data
    use PFLOTRAN_Constants_module

    use clm_time_manager  , only : get_step_size, get_nstep, nsstep, nestep,       &
                                   is_first_step, is_first_restart_step

    !use clm_varctl           , only : use_century_decomp
    !use CNDecompCascadeCNMod , only : decomp_rate_constants_cn
    !use CNDecompCascadeBGCMod, only : decomp_rate_constants_bgc

    implicit none

    type(bounds_type) , intent(in)  :: bounds
    integer, intent(in) :: num_soilc                  ! number of soil columns in filter
    integer, intent(in) :: filter_soilc(:)            ! filter for soil columns
!    integer, intent(in) :: num_soilp                  ! number of soil pfts in filter
!    integer, intent(in) :: filter_soilp(:)            ! filter for soil pfts

    type(clm_bgc_interface_data_type), intent(inout) :: clm_bgc_data

    !LOCAL VARIABLES:
    real(r8) :: dtime                         ! land model time step (sec)
    integer  :: nstep                         ! time step number
    integer  :: total_clmstep                 ! total clm time step number
    logical  :: ispfprint= .TRUE.             ! let PF printout or not
    logical  :: isinitpf = .FALSE.            ! (re-)initialize PF from CLM or not

#ifdef CLM_PF_DEBUG
    real(r8) :: t0, t1, t2, t3
#endif

  !-----------------------------------------------------------------------

    nstep = get_nstep() - nsstep
    dtime = get_step_size()

    if (is_first_step() .or. is_first_restart_step()) then
       isinitpf = .TRUE.
    else
       isinitpf = .FALSE.
    endif

#ifdef CLM_PF_DEBUG
if(nstep>=48*210 .and. nstep<=48*211) then
  call cpu_time(t0)
  if(pflotran_m%option%myrank .eq. pflotran_m%option%io_rank) then
      write(pflotran_m%option%myrank+200,*) '------------------------------------------------------- '
      write(pflotran_m%option%myrank+200,*) '------- checking CLM-PFLOTRAN timing - nstep = ', nstep
      write(pflotran_m%option%myrank+200,*) 'CPU_time @check-point 0: ', t0
  endif
endif
#endif


    ! (0)
    if (isinitpf) then
       total_clmstep = nestep - nsstep
       ispfprint = .true.               ! turn-on or shut-off PF's *.h5 output
       call pflotranModelUpdateFinalWaypoint(pflotran_m, total_clmstep*dtime, ispfprint)

       ! Set CLM soil properties onto PFLOTRAN grid

       call get_clm_soil_properties(clm_bgc_data, &
                    bounds, num_soilc, filter_soilc)
       call pflotranModelSetSoilProp(pflotran_m)

       ! if initializing soil 'TH' states from CLM to pflotran
       if (.not.initth_pf2clm) then
          call get_clm_soil_th(clm_bgc_data,initth_pf2clm, initth_pf2clm, &
                                bounds, num_soilc, filter_soilc)

          if (pf_hmode) then
          ! directly pass TH to internal PF vec (field%, work%)
              call pflotranModelSetInternalTHStatesfromCLM(pflotran_m)
          else
          ! pass TH to global_auxvar
              call pflotranModelUpdateTHfromCLM(pflotran_m, pf_hmode, pf_tmode)
          end if

       ! if initializaing CLM's H states from pflotran (only H mode now)
       else
          call pflotranModelGetSaturationFromPF(pflotran_m)   ! hydrological states
          call update_soil_moisture_pf2clm(clm_bgc_data,         &
                                bounds, num_soilc, filter_soilc)
       end if

       ! the following is for some specific PF's hydrological parameters useful to constrain H source/sink or BC
       if (pf_hmode) then
          call pflotranModelGetSoilPropFromPF(pflotran_m)
       end if

    endif !!if (isinitpf)

    ! (1) passing TH states from CLM to PF, if not H/TH mode NOT on in PF, every CLM time-step

    ! if PF T/H mode not available, have to pass those from CLM to global variable in PF to drive BGC/H
    if ((.not.pf_tmode .or. .not.pf_hmode) .and. (.not.isinitpf)) then
        call get_clm_soil_th(clm_bgc_data,pf_tmode, pf_hmode,  &
                            bounds, num_soilc, filter_soilc)

        call pflotranModelUpdateTHfromCLM(pflotran_m, pf_hmode, pf_tmode)
    endif

    ! ice-len adjusted porostiy
    if (.not.pf_frzmode) then
        call get_clm_iceadj_porosity(clm_bgc_data, &
                    bounds, num_soilc, filter_soilc)
        call pflotranModelResetSoilPorosityFromCLM(pflotran_m)
    endif

    ! (2) pass CLM water fluxes to CLM-PFLOTRAN interface
    if (pf_hmode) then      !if coupled 'H' mode between CLM45 and PFLOTRAN
        call get_clm_bcwflx(clm_bgc_data,           &
                    bounds, num_soilc, filter_soilc)

        ! pass flux 'vecs' from CLM to pflotran
        call pflotranModelUpdateHSourceSink(pflotran_m)    ! H SrcSink
        call pflotranModelSetSoilHbcsFromCLM(pflotran_m)   ! H BC
    end if

    ! (3) CLM thermal BC to PFLOTRAN
    if (pf_tmode) then
        call get_clm_bceflx(clm_bgc_data,           &
                    bounds, num_soilc, filter_soilc)

        call pflotranModelUpdateSubsurfTCond(pflotran_m)   ! SrcSink and T bc
    end if

    ! (4)
    if (pf_cmode) then
      ! (4a) for checking CLM's T/H response functions (TODO - also for passing decomposition rate to PF if needed)
      !if (use_century_decomp) then
      !    call decomp_rate_constants_bgc(bounds, num_soilc, filter_soilc)
      !else
      !    call decomp_rate_constants_cn(bounds, num_soilc, filter_soilc)
      !end if

      ! (4b) reset PFLOTRAN bgc state variables from CLM-CN, every CLM time-step
!      if (isinitpf) then     ! NOTE: if only initialize ONCE, uncomment this 'if...endif' block

        call get_clm_bgc_conc(clm_bgc_data,     &
                    bounds, num_soilc, filter_soilc)

        call pflotranModelSetBgcConcFromCLM(pflotran_m)

        if ( (.not.pf_hmode .or. .not.pf_frzmode)) then
        ! this is needed, because at step 0, PF's interface data is empty
        !which causes Aq. conc. adjustment balance issue
           call pflotranModelGetSaturationFromPF(pflotran_m)
        endif
!      endif     ! NOTE: if only initialize ONCE, uncomment this 'if...endif' block

      ! MUST reset PFLOTRAN soil aq. bgc state variables due to liq. water volume change
      ! when NOT coupled with PF Hydrology or NOT in freezing-mode (porosity will be forced to vary from CLM)
        if (.not.pf_hmode .or. .not.pf_frzmode) then
           call pflotranModelUpdateAqConcFromCLM(pflotran_m)
        endif

      ! (4c) bgc rate (source/sink) from CLM to PFLOTRAN
        call get_clm_bgc_rate(clm_bgc_data,  &
                    bounds, num_soilc, filter_soilc)

        call pflotranModelSetBgcRatesFromCLM(pflotran_m)
!write(*,*)">>>DEBUG | pflotranModelSetBgcRatesFromCLM"
    endif

#ifdef CLM_PF_DEBUG
if(nstep>=48*210 .and. nstep<=48*211) then
  call cpu_time(t1)
  if(pflotran_m%option%myrank .eq. pflotran_m%option%io_rank) then
      write(pflotran_m%option%myrank+200,*) 'CPU_time elapsed @check-point 1 - 0: ', t1-t0
  endif
endif
#endif

write(*,*)">>>DEBUG | pflotranModelStepperRunTillPauseTime: BEG...PFLOTRAN"
    ! (5) the main callings of PFLOTRAN
    call pflotranModelStepperRunTillPauseTime( pflotran_m, (nstep+1.0d0)*dtime, dtime, .false. )

write(*,*)">>>DEBUG | pflotranModelStepperRunTillPauseTime: END...CONTINUE..."

#ifdef CLM_PF_DEBUG
if(nstep>=48*210 .and. nstep<=48*211) then
  call cpu_time(t2)
  if(pflotran_m%option%myrank .eq. pflotran_m%option%io_rank) then
      write(pflotran_m%option%myrank+200,*) 'CPU_time elapsed @check-point 2 - 1: ', t2-t1
  endif
endif
#endif

    ! (6) update CLM variables from PFLOTRAN
    if (pf_hmode) then
        call pflotranModelGetSaturationFromPF(pflotran_m)   ! hydrological states

        call update_soil_moisture_pf2clm(clm_bgc_data,        &
                        bounds, num_soilc, filter_soilc)
    endif

    if (pf_tmode) then
        call pflotranModelGetTemperatureFromPF(pflotran_m)  ! thermal states

        call update_soil_temperature_pf2clm(clm_bgc_data,     &
                        bounds, num_soilc, filter_soilc)
    endif

    ! bgc variables
    if (pf_cmode) then
        call pflotranModelGetBgcVariablesFromPF(pflotran_m)

        call update_soil_bgc_pf2clm(clm_bgc_data,       &
                        bounds, num_soilc, filter_soilc)

        ! need to save the current time-step PF porosity/liq. saturation for bgc species mass conservation
        ! if CLM forced changing them into PF at NEXT timestep
        if (.not.pf_hmode .or. .not.pf_frzmode) then
           call pflotranModelGetSaturationFromPF(pflotran_m)
        endif
    endif

    ! the actual infiltration/runoff/drainage and solute flux with BC, if defined,
    ! are retrieving from PFLOTRAN using 'update_bcflow_pf2clm' subroutine
    ! TODO: comment out 'update_bcflow_pf2clm' currently
!     if (pf_hmode) then
!         call pflotranModelGetBCMassBalanceDeltaFromPF(pflotran_m )

!         call update_bcflow_pf2clm(              &
!            bounds, num_soilc, filter_soilc,     &
!            atm2lnd_vars,                        &
!            temperature_vars, energyflux_vars,   &
!            waterstate_vars, waterflux_vars)

!         !(TODO) bgc BC flows (e.g. runoff/leaching)
!     endif


#ifdef CLM_PF_DEBUG
if(nstep>=48*210 .and. nstep<=48*211) then
  call cpu_time(t3)
  if(pflotran_m%option%myrank .eq. pflotran_m%option%io_rank) then
      write(pflotran_m%option%myrank+200,*) 'CPU_time elapsed @check-point 3 - 2: ', t3-t2
      write(pflotran_m%option%myrank+200,*) 'CPU_time elapsed @check-point 3 - 0: ', t3-t0
      write(pflotran_m%option%myrank+200,*) '------------------------------------------------------- '
  endif
endif
#endif

  end subroutine pflotran_run_onestep

  !-----------------------------------------------------------------------
  !BOP
  !
  ! !ROUTINE: write_checkpoint
  !
  ! !INTERFACE:
  subroutine pflotran_write_checkpoint(date_stamp)
  !
  ! !DESCRIPTION:
  ! Trigger a pflotran checkpoint file to be written
  !
  ! !USES:
  ! !ARGUMENTS:
    character(len=32), intent(in) :: date_stamp ! file name date stamp

  ! !LOCAL VARIABLES:

  !EOP
  !-----------------------------------------------------------------------

    ! temporarily OFF - it's not working well for BGC
    ! So, now must initializing PF variables from CLM each start/restart.

    !call pflotranModelStepperCheckpoint(pflotran_m, date_stamp)

  end subroutine pflotran_write_checkpoint

  !-----------------------------------------------------------------------------
  !
  ! !IROUTINE: pflotran_finalize
  !
  ! !INTERFACE:
  subroutine pflotran_finalize()
  !
  ! !DESCRIPTION:
  !
  ! finalizing pflotran runs and destroying objects
  !
  ! !USES:

    implicit none

  !-----------------------------------------------------------------------

    if (use_pflotran) then
       call pflotranModelDestroy(pflotran_m)
    endif

  end subroutine pflotran_finalize


  ! ============================= GET CLM initial/src-sink/BC to PFLOTRAN ==================================

  !-----------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: get_clm_soil_properties
  !
  ! !INTERFACE:
  subroutine get_clm_soil_properties(clm_bgc_data, &
                    bounds, num_soilc, filter_soilc)
    !
    ! !DESCRIPTION:
    ! get soil column physical properties to PFLOTRAN
    !
    ! !USES:

    ! pflotran
    !
    ! !ARGUMENTS:

    implicit none

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscviewer.h"

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
!    type(soilstate_type)     , intent(in) :: soilstate_vars
    type(clm_bgc_interface_data_type), intent(in) :: clm_bgc_data

    integer  :: fc, g, l, c, j      ! indices
    integer  :: gcount, cellcount

    character(len= 32) :: subname = 'get_clm_soil_properties' ! subroutine name


    PetscScalar, pointer :: hksat_x_clm_loc(:) ! hydraulic conductivity in x-dir at saturation (mm H2O /s)
    PetscScalar, pointer :: hksat_y_clm_loc(:) ! hydraulic conductivity in y-dir at saturation (mm H2O /s)
    PetscScalar, pointer :: hksat_z_clm_loc(:) ! hydraulic conductivity in z-dir at saturation (mm H2O /s)
    PetscScalar, pointer :: watsat_clm_loc(:)  ! minimum soil suction (mm)
    PetscScalar, pointer :: sucsat_clm_loc(:)  ! volumetric soil water at saturation (porosity)
    PetscScalar, pointer :: bsw_clm_loc(:)     ! Clapp and Hornberger "b"
    PetscScalar, pointer :: watfc_clm_loc(:)
    PetscScalar, pointer :: bulkdensity_dry_clm_loc(:)
    PetscScalar, pointer :: zsoi_clm_loc(:)

    PetscErrorCode :: ierr

    associate( &
         ! Assign local pointer to derived subtypes components (column-level)
         clandunit  =>  col_pp%landunit   , & !  [integer (:)]  landunit index of column
         cgridcell  =>  col_pp%gridcell   , & !  [integer (:)]  gridcell index of column
         wtgcell    =>  col_pp%wtgcell    , & !  [real(r8) (:)]  weight (relative to gridcell
         cactive    =>  col_pp%active     , & !  [logical (:)]  column active or not
         z          =>  clm_bgc_data%z          , & !  [real(r8) (:,:)]  layer depth (m)
         dz         =>  clm_bgc_data%dz         , & !  [real(r8) (:,:)]  layer thickness depth (m)
         zi         =>  col_pp%zi         , & !  [real(r8) (:,:)]  interface level below a "z" level (m)
         !
         bd         =>  clm_bgc_data%bd_col         , & !
         bsw        =>  clm_bgc_data%bsw_col        , & !  [real(r8) (:,:)]  Clapp and Hornberger "b" (nlevgrnd)
         hksat      =>  clm_bgc_data%hksat_col      , & !  [real(r8) (:,:)]  hydraulic conductivity at saturation (mm H2O /s) (nlevgrnd)
         sucsat     =>  clm_bgc_data%sucsat_col     , & !  [real(r8) (:,:)]  minimum soil suction (mm) (nlevgrnd)
         watsat     =>  clm_bgc_data%watsat_col     , & !  [real(r8) (:,:)]  volumetric soil water at saturation (porosity) (nlevgrnd)
         watfc      =>  clm_bgc_data%watfc_col        & !  [real(r8) (:,:)]  volumetric soil water at saturation (porosity) (nlevgrnd)
         )

!-------------------------------------------------------------------------------------

    call VecGetArrayF90(clm_pf_idata%hksat_x_clmp, hksat_x_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%hksat_y_clmp, hksat_y_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%hksat_z_clmp, hksat_z_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%sucsat_clmp,  sucsat_clm_loc,  ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%watsat_clmp,  watsat_clm_loc,  ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%bsw_clmp,     bsw_clm_loc,     ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%watfc_clmp,  watfc_clm_loc,  ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%bulkdensity_dry_clmp,  bulkdensity_dry_clm_loc,  ierr)
    CHKERRQ(ierr)

    call VecGetArrayF90(clm_pf_idata%zsoi_clmp,  zsoi_clm_loc,  ierr)
    CHKERRQ(ierr)

    do fc = 1, num_soilc
      c = filter_soilc(fc)
      if ( wtgcell(c) <= 0._r8 .or. (.not.cactive(c)) ) cycle     ! don't assign data to PF for inactive cell

      ! Set gridcell and landunit indices
      g = cgridcell(c)
      l = clandunit(c)
      gcount = g - bounds%begg

      do j = 1,nlevsoi
          cellcount = gcount*clm_pf_idata%nzclm_mapped + j

          if (j <= clm_pf_idata%nzclm_mapped) then
            hksat_x_clm_loc(cellcount) = hksat_x_clm_loc(cellcount) &
                + hksat(c,j)*wtgcell(c)
            hksat_y_clm_loc(cellcount) = hksat_y_clm_loc(cellcount) &
                + hksat(c,j)*wtgcell(c)
            hksat_z_clm_loc(cellcount) = hksat_z_clm_loc(cellcount) &
                + hksat(c,j)*wtgcell(c)
            sucsat_clm_loc(cellcount) = sucsat_clm_loc(cellcount)   &
                + sucsat(c,j)*wtgcell(c)
            watsat_clm_loc(cellcount) = watsat_clm_loc(cellcount)   &
                + watsat(c,j)*wtgcell(c)
            bsw_clm_loc(cellcount)    = bsw_clm_loc(cellcount)      &
                + bsw(c,j)*wtgcell(c)
            watfc_clm_loc(cellcount)  = watfc_clm_loc(cellcount)    &
                + watfc(c,j)*wtgcell(c)
            bulkdensity_dry_clm_loc(cellcount) = bulkdensity_dry_clm_loc(cellcount) &
                + bd(c,j)*wtgcell(c)

            zsoi_clm_loc(cellcount) = z(c, j)     ! make sure this is right for multiple columns' situation

          else
            ! may need to further checking here
          endif

       enddo

    enddo ! do c = 1, numsoilc

!write(*,'(A40,10E14.6)')">>>DEBUG | hksat_x=",(hksat_x_clm_loc(1:10))
!write(*,'(A40,10E14.6)')">>>DEBUG | hksat_y=",(hksat_y_clm_loc(1:10))
!write(*,'(A40,10E14.6)')">>>DEBUG | hksat_z=",(hksat_z_clm_loc(1:10))
!write(*,'(A40,10E14.6)')">>>DEBUG | sucsat=",(sucsat_clm_loc(1:10))
!write(*,'(A40,10E14.6)')">>>DEBUG | watsat=",(watsat_clm_loc(1:10))
!write(*,'(A40,10E14.6)')">>>DEBUG | watfc=",(watfc_clm_loc(1:10))
!write(*,'(A40,10E14.6)')">>>DEBUG | bulkdensity=",(bulkdensity_dry_clm_loc(1:10))

    call VecRestoreArrayF90(clm_pf_idata%hksat_x_clmp, hksat_x_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%hksat_y_clmp, hksat_y_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%hksat_z_clmp, hksat_z_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%sucsat_clmp,  sucsat_clm_loc,  ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%watsat_clmp,  watsat_clm_loc,  ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%bsw_clmp,     bsw_clm_loc,     ierr)
    CHKERRQ(ierr)

    call VecRestoreArrayF90(clm_pf_idata%watfc_clmp,  watfc_clm_loc,  ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%bulkdensity_dry_clmp,  bulkdensity_dry_clm_loc,  ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%zsoi_clmp,  zsoi_clm_loc,  ierr)
    CHKERRQ(ierr)

    end associate
  end subroutine get_clm_soil_properties


  !-----------------------------------------------------------------------------
  !BOP
  !
  ! !ROUTINE: get_clm_soil_th
  !
  ! !INTERFACE:
  subroutine get_clm_soil_th(clm_bgc_data,pftmode, pfhmode,        &
           bounds, num_soilc, filter_soilc)

  !
  ! !DESCRIPTION:
  !  update soil temperature/saturation from CLM to PFLOTRAN for driving PF's BGC
  !  if either NOT available inside PFLOTRAN
  !
  ! !USES:
    use clm_time_manager    , only : get_nstep
    use shr_const_mod       , only : SHR_CONST_G

    use PFLOTRAN_Constants_module

  ! !ARGUMENTS:
    implicit none

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
    logical                  , intent(in) :: pftmode, pfhmode
    type(bounds_type)        , intent(in) :: bounds           ! bounds
    integer                  , intent(in) :: num_soilc        ! number of column soil points in column filter
    integer                  , intent(in) :: filter_soilc(:)  ! column filter for soil points
!    type(atm2lnd_type)       , intent(in) :: atm2lnd_vars
!    type(soilstate_type)     , intent(in) :: soilstate_vars
!    type(waterstate_type)    , intent(in) :: waterstate_vars
!    type(temperature_type)   , intent(in) :: temperature_vars
!    class(soil_water_retention_curve_type), intent(in) :: soil_water_retention_curve
    type(clm_bgc_interface_data_type), intent(in) :: clm_bgc_data

  ! !LOCAL VARIABLES:
    integer  :: fc, c, g, j, gcount, cellcount      ! indices
    real(r8) :: sattmp, psitmp, itheta
    real(r8) :: watmin(num_soilc, nlevsoi)
    real(r8) :: sucmin(num_soilc, nlevsoi)

    PetscScalar, pointer :: soilpress_clmp_loc(:)
    PetscScalar, pointer :: soilpsi_clmp_loc(:)
    PetscScalar, pointer :: soillsat_clmp_loc(:)
    PetscScalar, pointer :: soilisat_clmp_loc(:)
    PetscScalar, pointer :: soilt_clmp_loc(:)
    PetscErrorCode :: ierr

  !EOP
  !-----------------------------------------------------------------------
    associate ( &
      gridcell        => col_pp%gridcell     , & ! column's gridcell
      wtgcell         => col_pp%wtgcell      , & ! column's weight relative to gridcell
      cactive         => col_pp%active       , & ! [logical (:)]  column active or not
      dz              => clm_bgc_data%dz           , & ! layer thickness depth (m)
!      zi              => clm_bgc_data%zi           , & ! interface depth (m)
    !
      sucsat          => clm_bgc_data%sucsat_col    , & ! minimum soil suction (mm) (nlevgrnd)
      bsw             => clm_bgc_data%bsw_col       , & ! Clapp and Hornberger "b"
      watsat          => clm_bgc_data%watsat_col    , & ! volumetric soil water at saturation (porosity) (nlevgrnd)
      soilpsi         => clm_bgc_data%soilpsi_col   , & ! soil water matric potential in each soil layer (MPa)
    !
      h2osoi_liq      => clm_bgc_data%h2osoi_liq_col   , & ! liquid water (kg/m2)
      h2osoi_ice      => clm_bgc_data%h2osoi_ice_col   , & ! ice lens (kg/m2)
    !
      t_soisno        => clm_bgc_data%t_soisno_col    & ! snow-soil temperature (Kelvin)
!    !
!      forc_pbot       => atm2lnd_vars%forc_pbot_not_downscaled_grc  &  ! atmospheric pressure (Pa)
    )

    !--------------------------------------------------------------------------------------

    call VecGetArrayF90(clm_pf_idata%press_clmp, soilpress_clmp_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%soilpsi_clmp, soilpsi_clmp_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%soillsat_clmp, soillsat_clmp_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%soilisat_clmp, soilisat_clmp_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%soilt_clmp, soilt_clmp_loc, ierr)
    CHKERRQ(ierr)

    watmin(:,:) = 0.01_r8
    sucmin(:,:) = 1.e8_r8

    soilisat_clmp_loc(:) = 0._r8
    soillsat_clmp_loc(:) = 0._r8
    soilpsi_clmp_loc(:)  = 0._r8
    soilpress_clmp_loc(:)= clm_pf_idata%pressure_reference
    soilt_clmp_loc(:)    = 0._r8

    do fc = 1,num_soilc
      c = filter_soilc(fc)
      if ( wtgcell(c) <= 0._r8 .or. (.not.cactive(c)) ) cycle     ! don't assign data from PF for inactive cell

      g = gridcell(c)
      gcount = g - bounds%begg
      do j = 1, nlevsoi
        cellcount = gcount*clm_pf_idata%nzclm_mapped + j

        if (j<=clm_pf_idata%nzclm_mapped) then
          ! this adjusting should be done first, if PF-freezing-mode off, so that the following calculation can be done correctly
          itheta = 0._r8
          if(.not.pf_frzmode) then
             itheta = h2osoi_ice(c,j) / (dz(c,j) * denice)
             itheta = min(itheta, watsat(c,j)-watmin(c,j)/(dz(c,j)*denh2o))
             soilisat_clmp_loc(cellcount) = soilisat_clmp_loc(cellcount) + itheta/watsat(c,j)*wtgcell(c)
          endif

          if (.not.pfhmode) then
             ! porosity will be ice-adjusted for PF, if PF freezing-mode is off,
             ! so need to adjust 'psi' so that 'saturation' in PF is correct
             sattmp = h2osoi_liq(c,j) / ((watsat(c,j)-itheta)*dz(c,j)*denh2o)
             sattmp = min(max(sattmp, watmin(c,j)/watsat(c,j)),1._r8)
             soillsat_clmp_loc(cellcount) = soillsat_clmp_loc(cellcount) + sattmp*wtgcell(c)

             ! soil matric potential by Clapp-Hornburger method (this is the default used by CLM)
             ! but, this value IS different from what CLM used (not ice-content adjusted)
             ! So that in PF, if not ice-adjusted, the PSI is very small (negative) which implies possible water movement
             psitmp = sucsat(c,j) * (-SHR_CONST_G) * (sattmp**(-bsw(c,j)))  ! -Pa
             psitmp = min(max(psitmp,-sucmin(c,j)/SHR_CONST_G),0._r8)
             soilpsi_clmp_loc(cellcount)   = soilpsi_clmp_loc(cellcount) + psitmp*wtgcell(c)
             soilpress_clmp_loc(cellcount) = soilpress_clmp_loc(cellcount) + psitmp*wtgcell(c)

          endif

          if (.not.pftmode) then
             soilt_clmp_loc(cellcount) = soilt_clmp_loc(cellcount) + (t_soisno(c,j)-tfrz)*wtgcell(c)
          endif

        endif

      enddo
    enddo

!-----------------------------------------------------------------------------
!write(*,'(A30,12E14.6)')">>>DEBUG | soillsat=", soillsat_clmp_loc(1:10)
!write(*,'(A30,12E14.6)')">>>DEBUG | gsoilpsi[Pa]=", soilpsi_clmp_loc(1:10)
!write(*,'(A30,12E14.6)')">>>DEBUG | soilt[oC]=", soilt_clmp_loc(1:10)
!-----------------------------------------------------------------------------

    call VecRestoreArrayF90(clm_pf_idata%press_clmp, soilpress_clmp_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%soilpsi_clmp, soilpsi_clmp_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%soillsat_clmp, soillsat_clmp_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%soilisat_clmp, soilisat_clmp_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%soilt_clmp, soilt_clmp_loc, ierr)
    CHKERRQ(ierr)

   end associate
  end subroutine get_clm_soil_th


  !-----------------------------------------------------------------------------
  !BOP
  !
  ! !ROUTINE: get_clm_iceadj_porosity
  !
  ! !INTERFACE:
  subroutine get_clm_iceadj_porosity(clm_bgc_data, &
                    bounds, num_soilc, filter_soilc)
  !
  ! !DESCRIPTION:
  !  update soil effective porosity from CLM to PFLOTRAN if PF freezing mode is off
  !
  ! !USES:

    use PFLOTRAN_Constants_module
    use clm_varctl          , only : pf_frzmode

  ! !ARGUMENTS:
    implicit none

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

    type(bounds_type)        , intent(in) :: bounds           ! bounds
    integer                  , intent(in) :: num_soilc        ! number of column soil points in column filter
    integer                  , intent(in) :: filter_soilc(:)  ! column filter for soil points
!    type(soilstate_type)     , intent(in) :: soilstate_vars
!    type(waterstate_type)    , intent(in) :: waterstate_vars
    type(clm_bgc_interface_data_type), intent(in) :: clm_bgc_data

  ! !LOCAL VARIABLES:
    integer  :: fc, c, g, j, gcount, cellcount       ! indices
    real(r8) :: itheta

    PetscScalar, pointer :: adjporosity_clmp_loc(:)  !
    PetscScalar, pointer :: soilisat_clmp_loc(:)  !
    PetscErrorCode :: ierr

  !EOP
  !-----------------------------------------------------------------------
    associate ( &
    gridcell        => col_pp%gridcell     , & ! column's gridcell
    wtgcell         => col_pp%wtgcell      , & ! column's weight relative to gridcell
    cactive         => col_pp%active       , & ! column's active or not
    dz              => col_pp%dz           , & ! layer thickness depth (m)
    watsat          => clm_bgc_data%watsat_col       , & ! volumetric soil water at saturation (porosity) (nlevgrnd)
    h2osoi_ice      => clm_bgc_data%h2osoi_ice_col    & ! ice lens (kg/m2)
    )

    ! if 'pf_tmode' is NOT using freezing option, the phase-change of soil water done in 'SoilTemperatureMod.F90' in 'bgp2'
    ! must be included to adjust porosity (effective porosity) in pflotran
    ! This is doing prior to the real liquid water source/sink, because 'h2osoi_liq' will be updated during those calls after 'bgp2'.
    if (.not. pf_frzmode) then

        ! re-calculate the effective porosity (CLM ice-len adjusted), which should be pass to pflotran
        call VecGetArrayF90(clm_pf_idata%effporosity_clmp, adjporosity_clmp_loc,  ierr)
        CHKERRQ(ierr)
        call VecGetArrayF90(clm_pf_idata%soilisat_clmp, soilisat_clmp_loc,  ierr)
        CHKERRQ(ierr)

        adjporosity_clmp_loc(:) = 0._r8
        soilisat_clmp_loc(:) = 0._r8
        do fc = 1,num_soilc
           c = filter_soilc(fc)
           if ( wtgcell(c) <= 0._r8 .or. (.not.cactive(c)) ) cycle     ! don't assign data from PF for inactive cell

           g = gridcell(c)
           gcount = g - bounds%begg
           do j = 1, nlevsoi
             cellcount = gcount*clm_pf_idata%nzclm_mapped + j

             if (j<=clm_pf_idata%nzclm_mapped) then
               itheta = h2osoi_ice(c,j) / (dz(c,j) * denice)
               itheta = min(itheta, 0.99_r8*watsat(c,j))
               adjporosity_clmp_loc(cellcount) = adjporosity_clmp_loc(cellcount) &
                 + (watsat(c,j) - itheta) * wtgcell(c)
               soilisat_clmp_loc(cellcount)    = soilisat_clmp_loc(cellcount) &
                 + itheta * wtgcell(c)
             endif
           end do
        end do
        call VecRestoreArrayF90(clm_pf_idata%effporosity_clmp,  adjporosity_clmp_loc,  ierr)
        CHKERRQ(ierr)
        call VecRestoreArrayF90(clm_pf_idata%soilisat_clmp, soilisat_clmp_loc,  ierr)
        CHKERRQ(ierr)

    end if

    end associate
  end subroutine get_clm_iceadj_porosity


  !-----------------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: get_clm_bcwflx
  !
  ! !INTERFACE:
  subroutine get_clm_bcwflx(clm_bgc_data,           &
       bounds, num_soilc, filter_soilc)
  !
  ! !DESCRIPTION:
  !
  !  F.-M. YUAN: the water fluxes in CLM4.5 are separately calculated in a few subroutines
  !        in 'SoilHydrologyMod.F90'. When coupled with pflotran, it's hard to get those together
  !        like GB does in 'step_th_clm_pf' subroutine. So, this subroutine is a collective call of
  !        that and others in 'Hydrology2Mod.F90' so that pflotran can be called out of 'hydrology2'.
  !
  ! !USES:

    use shr_const_mod   , only : SHR_CONST_G
    use clm_time_manager, only : get_step_size, get_nstep

  ! !ARGUMENTS:
    implicit none

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

    type(bounds_type) , intent(in)  :: bounds
    integer, intent(in) :: num_soilc                 ! number of column non-lake points in column filter
    integer, intent(in) :: filter_soilc(:)           ! column filter for non-lake points

!    type(atm2lnd_type)       , intent(in) :: clm_a2l
!    type(soilstate_type)     , intent(in) :: soils_vars
!    type(temperature_type)   , intent(in) :: ces_vars
!    type(energyflux_type)    , intent(in) :: cef_vars
!    type(waterstate_type)    , intent(in) :: cws_vars
!    type(waterflux_type)     , intent(in) :: cwf_vars
    type(clm_bgc_interface_data_type), intent(in) :: clm_bgc_data

  ! !LOCAL VARIABLES:
    integer  :: fc, g, c, j, p             ! do loop indices
    integer  :: gcount                     ! gridcell index (0-based)
    integer  :: pftindex                   ! pft index
    real(r8) :: dtime                      ! land model time step (sec)
    integer  :: nstep                      ! time step number
    real(r8) :: area
    integer  :: cellcount                  ! 3-D cell index (1-based)

    real(r8) :: rootfr_col(bounds%begc:bounds%endc, nlevsoi)
    real(r8) :: rootfr_sum(bounds%begc:bounds%endc)           ! accumulator for rootr weighting
    real(r8) :: qflx_evap_col(bounds%begc:bounds%endc)        ! soil surface evaporation (mmH2O/s)
    real(r8) :: qflx_tran_col(bounds%begc:bounds%endc)        ! veg. transpiration (mmH2O/s)

    real(r8) :: qflx, qflx_sink, qflx_source, soilvwc
    real(r8) :: dsoilliq1 = 0._r8, dsoilliq2 = 0._r8, dsoilliq3 = 0._r8

    real(r8) :: qflx_ground, kbot
    real(r8) :: reference_pressure, ponding_pressure      ! Pa
    real(r8) :: pondmax(bounds%begc:bounds%endc)          ! mm H2O: max. ponding depth for column
    real(r8) :: sr      = 0.10_r8
    real(r8) :: tempreal

    ! for PF --> CLM (seq.)
    PetscScalar, pointer :: press_clms_loc(:)       !
    PetscScalar, pointer :: soillsat_clms_loc(:)    !
    PetscScalar, pointer :: porosity_clms_loc(:)    !
    PetscScalar, pointer :: sr_pcwmax_clms_loc(:)   !

    PetscScalar, pointer :: area_clms_loc(:)         !

    ! for CLM (mpi) --> PF
    PetscScalar, pointer :: zsoi_clmp_loc(:)         !
    PetscScalar, pointer :: qflx_clmp_loc(:)         !   source/sink term for plant Transpiration: unit in mass rate (kgH2O/sec)
    PetscScalar, pointer :: press_top_clmp_loc(:)    !   BC in pressure type: unit in Pa
    PetscScalar, pointer :: press_base_clmp_loc(:)   !
    PetscScalar, pointer :: qflux_top_clmp_loc(:)    !   BC in neumann flux type: unit in m/s
    PetscScalar, pointer :: qflux_base_clmp_loc(:)   !
    PetscScalar, pointer :: press_maxponding_clmp_loc(:)   !
    PetscErrorCode :: ierr

  !EOP
  !-----------------------------------------------------------------------
    associate ( &
    ltype             => lun_pp%itype             , & ! landunit type
    cgridcell         => col_pp%gridcell          , & ! column's gridcell
    clandunit         => col_pp%landunit          , & ! column's landunit
    zi                => col_pp%zi                , & ! Input: (:,:) soil layer interface depth (m)
    dz                => col_pp%dz                , & ! Input: (:,:) soil layer thickness (m)
    pfti              => col_pp%pfti                                , &! beginning pft index for each column
    pwtgcell          => pft_pp%wtgcell                             , &! weight relative to gridcell for each pft
    pwtcol            => pft_pp%wtcol                               , &! weight relative to column for each pft
    !
    bsw               => clm_bgc_data%bsw_col                  , &! Clapp and Hornberger "b" (nlevgrnd)
    hksat             => clm_bgc_data%hksat_col                , &! hydraulic conductivity at saturation (mm H2O /s) (nlevgrnd)
    watsat            => clm_bgc_data%watsat_col               , &! volumetric soil water at saturation (porosity) (nlevgrnd)
    sucsat            => clm_bgc_data%sucsat_col               , &! minimum soil suction (mm) (nlevgrnd)
    rootfr            => clm_bgc_data%rootfr_col               , &
!    rootfr_pft        => soils_vars%rootfr_patch             , & ! pft-level effective fraction of roots in each soil layer

    forc_pbot         => clm_bgc_data%forc_pbot_not_downscaled_grc , & ! Input:  [real(r8) (:)]  atmospheric pressure (Pa)
    t_grnd            => clm_bgc_data%t_grnd_col                  , & ! Input:  [real(r8) (:)]  ground surface temperature [K]
    htvp              => clm_bgc_data%htvp_col                    , & ! Input:  [real(r8) (:)]  latent heat of vapor of water (or sublimation) [j/kg]
    !
    frac_sno          => clm_bgc_data%frac_sno_eff_col            , & ! Input: fraction of ground covered by snow (0 to 1)
    frac_h2osfc       => clm_bgc_data%frac_h2osfc_col             , & ! Input: fraction of ground covered by surface water (0 to 1)
    h2osoi_liq        => clm_bgc_data%h2osoi_liq_col              , & ! Input: liquid water (kg/m2)
    h2osoi_ice        => clm_bgc_data%h2osoi_ice_col              , & ! Input: ice lens (kg/m2)
    !
    qflx_top_soil     => clm_bgc_data%qflx_top_soil_col           , & ! Input: net water input into soil from top (mm/s)
    qflx_ev_h2osfc    => clm_bgc_data%qflx_ev_h2osfc_col          , & ! Input: column-level evaporation flux from h2osfc (W/m2) [+ to atm] : checking unit
    qflx_evap_soil    => clm_bgc_data%qflx_evap_soi_col           , & ! Input: column-level soil evaporation (mm H2O/s) (+ = to atm)
    qflx_subl_snow    => clm_bgc_data%qflx_sub_snow_col           , & ! Input: column-level evaporation flux from snow (mm H2O/s) [+ to atm]
!    qflx_tran_veg_pft => cwf_vars%qflx_tran_veg_patch         &   ! Input: pft-level vegetation transpiration (mm H2O/s) (+ = to atm)
    qflx_tran_veg     => clm_bgc_data%qflx_tran_veg_col             &

    )

!----------------------------------------------------------------------------
    nstep = get_nstep()
    dtime = get_step_size()

    ! (1) soil surface evaporation: needs further checking here? -
    do fc = 1, num_soilc
       c = filter_soilc(fc)
       if (ltype(clandunit(c)) == istsoil .or. ltype(clandunit(c))==istcrop) then
           ! not sure if using 'qflx_evap_soi' as a whole ground-surface better than individual surfaces, i.e. 'qflx_ev_snow/soi/h2osfc'
           ! all of those 4 variables are calculated in 'BareGroundFluxesMod', 'CanopyFluxesMod', and then adjusted in 'Biogeophysics2' after soil temperature call
           ! note that: all 4 variables could be negative (i.e., dew formation on ground)
           qflx_evap_col(c)=(1.0_r8 - frac_sno(c) - frac_h2osfc(c))*qflx_evap_soil(c) + &
                          frac_h2osfc(c)*qflx_ev_h2osfc(c)/htvp(c)
                          !frac_sno(c)*qflx_ev_snow(c)   ! snow-covered area should be excluded (see SoilHydrologyMod:: infiltration)
       else
          ! for other types of landunits
          qflx_evap_col(c) = (1.0_r8 - frac_sno(c))*qflx_evap_soil(c)
       end if

       if (t_grnd(c) <= tfrz) qflx_evap_col(c) = max(0._r8, qflx_evap_col(c))   ! frozen ground, no dew contribution to subsurface infiltration
    end do

    ! (3) Compute the vegetation Transpiration (originally those are in 'SoilHydrologyMod.F90')
    do j = 1, nlevsoi
        do fc = 1, num_soilc
            c = filter_soilc(fc)
            rootfr_col(c,j) = 0._r8
        end do
    end do
    rootfr_sum(:)    = 0._r8
    qflx_tran_col(:) = 0._r8

!    do pftindex = 1, max_patch_per_col
!       do fc = 1, num_soilc
!          c = filter_soilc(fc)
!          if (pftindex <= col_pp%npfts(c)) then
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

    do fc = 1, num_soilc
        c = filter_soilc(fc)
            qflx_tran_col(c)  = qflx_tran_veg(c)
        do j = 1, nlevsoi
            rootfr_col(c,j) = rootfr(c,j)
        end do

    end do

!----------------------------------------------------------------------------------------------------------
    ! (4) pass the clm_qflx to the vecs
    ! NOTE the following unit conversions:
    ! qflx_soil_top and qflx_tran_veg are in [mm/sec] from CLM;
    ! qflx_clm_loc is in [kgH2O/sec] as mass rate for pflotran (as input)

    ! previous time-step soil water pressure and saturation for adjusting qflx
    ! note that this is a temporary workaround - waiting for PF's solution
    call VecGetArrayF90(clm_pf_idata%press_clms, press_clms_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%soillsat_clms, soillsat_clms_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%effporosity_clms, porosity_clms_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%sr_pcwmax_clms, sr_pcwmax_clms_loc, ierr)
    CHKERRQ(ierr)

    call VecGetArrayF90(clm_pf_idata%area_top_face_clms, area_clms_loc, ierr)
    CHKERRQ(ierr)

    call VecGetArrayF90(clm_pf_idata%zsoi_clmp, zsoi_clmp_loc, ierr)
    CHKERRQ(ierr)

    call VecGetArrayF90(clm_pf_idata%qflux_clmp, qflx_clmp_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%press_subsurf_clmp, press_top_clmp_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%press_subbase_clmp, press_base_clmp_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%qflux_subsurf_clmp, qflux_top_clmp_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%qflux_subbase_clmp, qflux_base_clmp_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%press_maxponding_clmp, press_maxponding_clmp_loc, ierr)
    CHKERRQ(ierr)

    ! Initialize flux variables and bcs
    do fc = 1, num_soilc

      c = filter_soilc(fc)
      g = cgridcell(c)
      gcount = g - bounds%begg

      do j = 1, nlevsoi
        cellcount = gcount*clm_pf_idata%nzclm_mapped + j

        if(j<=clm_pf_idata%nzclm_mapped) then
          qflx_clmp_loc(cellcount) = 0.0_r8

          if (j .eq. 1) then
             press_top_clmp_loc(gcount+1)  = press_clms_loc(cellcount)   ! same as the first top layer
          end if

          if (j .eq. clm_pf_idata%nzclm_mapped) then
             press_base_clmp_loc(gcount+1) = press_clms_loc(cellcount)   ! same as the bottom layer
          end if

        endif
      end do

    end do

    pondmax(:) = 0._r8   ! this is temporarily set (not yet figure out how CLM get this value)
    do fc = 1, num_soilc

       c = filter_soilc(fc)
       g = cgridcell(c)
       gcount = g - bounds%begg
       area = area_clms_loc(gcount*clm_pf_idata%nzclm_mapped+1)
       reference_pressure = clm_pf_idata%pressure_reference
       ponding_pressure = pondmax(c)*SHR_CONST_G              ! max. ponding water depth (mm) ==> pressure (Pa)

       press_maxponding_clmp_loc(gcount+1) = reference_pressure+ponding_pressure*col_pp%wtgcell(c)

       do j = 1, nlevsoi
         cellcount = gcount*clm_pf_idata%nzclm_mapped + j

         if(j<=clm_pf_idata%nzclm_mapped) then
           qflx = 0._r8
           ! top BC
           if (j .eq. 1) then

             ! net liq water input/output to soil column
             qflx_ground = qflx_top_soil(c) - qflx_evap_col(c)  ! unit: mm/sec

             ! if net input potential, it's forming TOP BC of pressure type (water ponding potetial)
             ! both waterhead and flux calcuated here, but not applied in PFLOTRAN in the same time (upon BC type picked-up by PF)
             if ( qflx_ground .gt. 0._r8) then
                ! Newly ADDED mmH2O ==> pressure (Pa) as top BC (dirichlet) by forming a layer of surface water column
                ! AND, the actual infiltration/runoff are retrieving from PFLOTRAN using 'update_surflow_pf2clm' subroutine
                if (soillsat_clms_loc(cellcount) >= 1._r8) then
                   ! water-head formed on saturated below-ground soil layer
                   press_top_clmp_loc(gcount+1) = press_clms_loc(cellcount) + &
                          qflx_ground*col_pp%wtgcell(c)*dtime*SHR_CONST_G
                else
                   ! ground-water-head discontinued from below-ground (atm. pressure applied at both ends)
                   press_top_clmp_loc(gcount+1) = press_top_clmp_loc(gcount+1) +  &
                          qflx_ground*col_pp%wtgcell(c)*dtime*SHR_CONST_G
                endif

                ! mmH2O/sec ==> mH2O/sec of potential infiltration (flux) rate as top BC (neumann)
                ! must be used together with 'seepage' as top BC in the meantime to removing upwarding water
                ! AND, the actual infiltration/runoff are retrieving from PFLOTRAN using 'update_bcflow_pf2clm' subroutine
                qflux_top_clmp_loc(gcount+1) = qflx_ground*1.e-3

             ! if net loss potential, it's as source/sink term of soil column
             else
                qflx = qflx + min(0._r8, qflx_ground)              ! unit here: mmH2O/sec
                qflux_top_clmp_loc(gcount+1) = 0._r8
             end if
          end if

          ! adding plant root extraction of water (transpiration)
          qflx = qflx - qflx_tran_col(c)*rootfr_col(c,j)  ! by this point: unit: mmH2O/sec for CLM column

          qflx = qflx * col_pp%wtgcell(c)                    ! from now on: per PF 3-D cells
          qflx = qflx * area * 1.e-3 *denh2o              ! unit: mmH2O/sec ==> kgH2O/sec

          ! previous time-step soil water saturation for adjusting qflx to avoid too wet or too dry to cause PF math issue
          ! (this is a temporary workaround - waiting for PF's solution)
          soilvwc = soillsat_clms_loc(cellcount) *  &
                    porosity_clms_loc(cellcount)      ! PF saturation ==> real vwc (using adjusted porosity???)

          ! checking if over-filled when sinking (excluding infiltration)
          qflx_sink = max(0._r8, qflx)             ! sink (+) only (kgH2O/sec)
          dsoilliq1 = (0.99_r8*porosity_clms_loc(cellcount)-soilvwc) &
                  * zsoi_clmp_loc(cellcount) &
                  * area * denh2o /dtime    ! mH2O ==> kgH2O/sec to be filled at most (1% for error)
          qflx_sink = min(qflx_sink, max(0._r8,dsoilliq1))

          ! checking if too dry to be ETed (or other sourced): lower than 'sr_pcwmax'
          qflx_source = min(0._r8, qflx)           ! source (-) only (kgH2O/sec)
          sr = 1.001_r8*sr_pcwmax_clms_loc(cellcount) * &              ! '1.001' will give 0.1% for holding error in the calculation
                    watsat(c,j)                                                 ! PF saturation ==> 'real' vwc
          dsoilliq2 = (sr-soilvwc)*zsoi_clmp_loc(cellcount)*area*denh2o/dtime                     ! mH2O ==> kgH2O/sec to be extracted at most (-)
          qflx_source = max(qflx_source, min(0._r8,dsoilliq2))

          qflx_clmp_loc(cellcount) = qflx_clmp_loc(cellcount)+ (qflx_sink+qflx_source)             ! source/sink unit: kg/sec

          ! bottom BC (neumman type): m/sec
          if (j .eq. clm_pf_idata%nzclm_mapped) then
             ! available water flux-out rate (-) adjusted by source(-)/sink(+) term
             dsoilliq3 = min(0._r8, dsoilliq2 - qflx_clmp_loc(cellcount)) &
                             /area/denh2o                                      ! kgH2O/sec ==> mH2O/sec

             ! free drainage at bottom
             tempreal = soilvwc/watsat(c,j)                                    ! using 'real' saturation
             kbot = hksat(c,j)*(tempreal**(2._r8*bsw(c,j)+3._r8))*1.e-3        ! mmH2O/sec ==> mH2O/sec
             qflux_base_clmp_loc(gcount+1) =  qflux_base_clmp_loc(gcount+1) &
                    + max(dsoilliq3, -kbot * col_pp%wtgcell(c))             ! mH2O/sec

          end if

         endif !if(j<=clm_pf_idata%nzclm_mapped)

       end do  ! do j=1,nlevsoi

    end do ! do c=1, num_soilc

    call VecRestoreArrayF90(clm_pf_idata%press_clms, press_clms_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%soillsat_clms, soillsat_clms_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%sr_pcwmax_clms, sr_pcwmax_clms_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%effporosity_clms, porosity_clms_loc, ierr)
    CHKERRQ(ierr)

    call VecRestoreArrayF90(clm_pf_idata%area_top_face_clms, area_clms_loc, ierr)
    CHKERRQ(ierr)

    call VecRestoreArrayF90(clm_pf_idata%zsoi_clmp, zsoi_clmp_loc, ierr)
    CHKERRQ(ierr)

    call VecRestoreArrayF90(clm_pf_idata%qflux_clmp, qflx_clmp_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%press_subsurf_clmp, press_top_clmp_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%press_subbase_clmp, press_base_clmp_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%qflux_subsurf_clmp, qflux_top_clmp_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%qflux_subbase_clmp, qflux_base_clmp_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%press_maxponding_clmp, press_maxponding_clmp_loc, ierr)
    CHKERRQ(ierr)

!----------------------------------------------------------------------------------------------------------

  end associate
  end subroutine get_clm_bcwflx

  !-----------------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: get_clm_bceflx
  !
  ! !INTERFACE:
  subroutine get_clm_bceflx(clm_bgc_data,           &
       bounds, num_soilc, filter_soilc)
  !
  ! !DESCRIPTION:
  !
  !  F.-M. YUAN: the boundary heat fluxes in CLM4.5 are extracted to drive pflotran TH mode.
  !        GB only defined ground-heaf-flux.
  !        So, this subroutine is a collective setting on either heat-flux (neumann type)
  !        or interface thermal state (temperature) (dirichlet type)
  !        at both ground and bottom interface (BC).
  !
  ! !USES:

  ! !ARGUMENTS:
    implicit none

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

    type(bounds_type) , intent(in)  :: bounds
    integer, intent(in) :: num_soilc                 ! number of column non-lake points in column filter
    integer, intent(in) :: filter_soilc(:)           ! column filter for non-lake points

!     type(atm2lnd_type)       , intent(in) :: clm_a2l
!     type(waterstate_type)    , intent(in) :: cws_vars
!     type(temperature_type)   , intent(in) :: ces_vars
!     type(energyflux_type)    , intent(in) :: cef_vars
    type(clm_bgc_interface_data_type), intent(in) :: clm_bgc_data

  ! !LOCAL VARIABLES:
    integer  :: fc, c, g, p, gcount            ! do loop indices
    integer  :: pftindex

    ! for CLM (mpi) --> PF
    PetscScalar, pointer :: gflux_subsurf_clmp_loc(:)    !   BC in neumman type: unit W/m2
    PetscScalar, pointer :: gtemp_subsurf_clmp_loc(:)    !   BC in dirichlet type: unit in degC
    PetscScalar, pointer :: gflux_subbase_clmp_loc(:)    !   BC in neumman type: unit W/m2
    PetscScalar, pointer :: gtemp_subbase_clmp_loc(:)    !   BC in dirichlet type: unit in degC
    PetscErrorCode :: ierr

  !EOP
  !-----------------------------------------------------------------------
    associate ( &
    cgridcell         => col_pp%gridcell                            , &! column's gridcell
    clandunit         => col_pp%landunit                            , &! column's landunit
    dz                => col_pp%dz                                  , &! layer thickness depth (m)
    pfti              => col_pp%pfti                                , &! beginning pft index for each column
    pwtgcell          => pft_pp%wtgcell                             , &! weight relative to gridcell for each pft
    pwtcol            => pft_pp%wtcol                               , &! weight relative to column for each pft
    !
    frac_sno          => clm_bgc_data%frac_sno_eff_col               , & ! Input: fraction of ground covered by snow (0 to 1)
    frac_h2osfc       => clm_bgc_data%frac_h2osfc_col                , & ! Input: fraction of ground covered by surface water (0 to 1)
    !
    eflx_bot          => clm_bgc_data%eflx_bot_col                   , &! heat flux from beneath column (W/m**2) [+ = upward]
!     eflx_gnet         => cef_vars%eflx_gnet_patch                , &! net ground heat flux into the surface (W/m**2) per patch
!     eflx_soil_grnd    => cef_vars%eflx_soil_grnd_patch           , &! soil heat flux (W/m**2) [+ = into soil]
    eflx_gnet         => clm_bgc_data%eflx_gnet_col               , &
    eflx_soil_grnd    => clm_bgc_data%eflx_soil_grnd_col          , &
    t_grnd            => clm_bgc_data%t_grnd_col                   & ! ground surface temperature [K]
    )

 !----------------------------------------------------------------------------------------------------------
    ! (1) pass the clm_gflux/gtemp to the vec

    call VecGetArrayF90(clm_pf_idata%gflux_subsurf_clmp, gflux_subsurf_clmp_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%gflux_subbase_clmp, gflux_subbase_clmp_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%gtemp_subsurf_clmp, gtemp_subsurf_clmp_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%gtemp_subbase_clmp, gtemp_subbase_clmp_loc, ierr)
    CHKERRQ(ierr)

    gflux_subsurf_clmp_loc(:) = 0._r8
    gflux_subbase_clmp_loc(:) = 0._r8
    gtemp_subsurf_clmp_loc(:) = 0._r8
    gtemp_subbase_clmp_loc(:) = 0._r8

    ! CLM appears NO column-level ground-heat-flux variable, instead by 'patch'
    do fc = 1, num_soilc
       c = filter_soilc(fc)
       if ( col_pp%wtgcell(c) <= 0._r8 .or. (.not.col_pp%active(c)) ) cycle     ! don't assign data from PF for inactive cell
       g = cgridcell(c)
       gcount = g - bounds%begg

       gflux_subsurf_clmp_loc(gcount+1)  = gflux_subsurf_clmp_loc(gcount+1) &
                                         + eflx_soil_grnd(c)*1.d-3 * col_pp%wtgcell(c)     ! 1.d-3: from W/m2 --> kJ/m2/s
!        do pftindex = 1, max_patch_per_col
!           if (pftindex <= col_pp%npfts(c)) then
!              p = pfti(c) + pftindex - 1
!              if (pwtgcell(p)>0._r8) then
!                gflux_subsurf_clmp_loc(gcount+1)  = gflux_subsurf_clmp_loc(gcount+1) &
!                   + eflx_soil_grnd(p)*1.d-3 * pwtcol(p)           ! (TODO - checking) from W/m2 --> kJ/m2/s
!             end if
!          end if
!        end do
    end do

    ! CLM column-level variables available to PFLOTRAN
    do fc = 1,num_soilc
       c = filter_soilc(fc)
       if ( col_pp%wtgcell(c) <= 0._r8 .or. (.not.col_pp%active(c)) ) cycle     ! don't assign data from PF for inactive cell

       g = cgridcell(c)
       gcount = g - bounds%begg

       gflux_subbase_clmp_loc(gcount+1)  =  gflux_subbase_clmp_loc(gcount+1) &
          + eflx_bot(c)*1.d-3 * col_pp%wtgcell(c)
       gtemp_subsurf_clmp_loc(gcount+1)  = gflux_subbase_clmp_loc(gcount+1) &
          + (t_grnd(c) - tfrz) * col_pp%wtgcell(c)

       gtemp_subbase_clmp_loc(gcount+1)  = 0._r8                 ! not yet get it from CLM (i.e.,dirichlet type not available)

    end do

    call VecRestoreArrayF90(clm_pf_idata%gflux_subsurf_clmp, gflux_subsurf_clmp_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%gflux_subbase_clmp, gflux_subbase_clmp_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%gtemp_subsurf_clmp, gtemp_subsurf_clmp_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%gtemp_subbase_clmp, gtemp_subbase_clmp_loc, ierr)
    CHKERRQ(ierr)

  end associate
  end subroutine get_clm_bceflx


  !-----------------------------------------------------------------------------
  !BOP
  !
  ! !ROUTINE: get_clm_bgc_conc(bounds)
  !
  ! !INTERFACE:

  subroutine get_clm_bgc_conc(clm_bgc_data, &
           bounds, num_soilc, filter_soilc)
!! TODO: add phosphorus vars
#ifndef FLEXIBLE_POOLS
    use clm_varpar, only : i_met_lit, i_cel_lit, i_lig_lit, i_cwd
#endif

    implicit none

    type(bounds_type)        , intent(in) :: bounds
    integer                  , intent(in) :: num_soilc         ! number of soil columns in filter
    integer                  , intent(in) :: filter_soilc(:)   ! filter for soil columns

!    type(carbonstate_type)      , intent(in) :: carbonstate_vars
!    type(nitrogenstate_type)    , intent(in) :: nitrogenstate_vars
!    type(phosphorusstate_type)  , intent(in) :: phosphorusstate_vars
!    type(ch4_type)              , intent(in) :: ch4_vars
    type(clm_bgc_interface_data_type), intent(in) :: clm_bgc_data

    character(len=256) :: subname = "get_clm_bgc_concentration"

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

    ! Local variables
    integer  :: g, fc, c, j, k
    integer  :: gcount, cellcount
    real(r8) :: wtgcell, realc_gcell, realn_gcell

#ifdef FLEXIBLE_POOLS
    integer  :: vec_offset
    PetscScalar, pointer :: decomp_cpools_vr_clm_loc(:)      ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
    PetscScalar, pointer :: decomp_npools_vr_clm_loc(:)      ! (gN/m3) vertically-resolved decomposing (litter, cwd, soil) N pools
#else
    integer  :: isom
    PetscScalar, pointer :: decomp_cpools_vr_lit1_clm_loc(:) ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
    PetscScalar, pointer :: decomp_cpools_vr_lit2_clm_loc(:) ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
    PetscScalar, pointer :: decomp_cpools_vr_lit3_clm_loc(:) ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
    PetscScalar, pointer :: decomp_cpools_vr_cwd_clm_loc(:)  ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
    PetscScalar, pointer :: decomp_cpools_vr_som1_clm_loc(:) ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
    PetscScalar, pointer :: decomp_cpools_vr_som2_clm_loc(:) ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
    PetscScalar, pointer :: decomp_cpools_vr_som3_clm_loc(:) ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
    PetscScalar, pointer :: decomp_cpools_vr_som4_clm_loc(:) ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
    PetscScalar, pointer :: decomp_npools_vr_lit1_clm_loc(:) ! (gN/m3) vertically-resolved decomposing (litter, cwd, soil) N pools
    PetscScalar, pointer :: decomp_npools_vr_lit2_clm_loc(:) ! (gN/m3) vertically-resolved decomposing (litter, cwd, soil) N pools
    PetscScalar, pointer :: decomp_npools_vr_lit3_clm_loc(:) ! (gN/m3) vertically-resolved decomposing (litter, cwd, soil) N pools
    PetscScalar, pointer :: decomp_npools_vr_cwd_clm_loc(:)  ! (gN/m3) vertically-resolved decomposing (litter, cwd, soil) N pools
    PetscScalar, pointer :: decomp_npools_vr_som1_clm_loc(:) ! (gN/m3) vertically-resolved decomposing (litter, cwd, soil) n pools
    PetscScalar, pointer :: decomp_npools_vr_som2_clm_loc(:) ! (gN/m3) vertically-resolved decomposing (litter, cwd, soil) n pools
    PetscScalar, pointer :: decomp_npools_vr_som3_clm_loc(:) ! (gN/m3) vertically-resolved decomposing (litter, cwd, soil) n pools
    PetscScalar, pointer :: decomp_npools_vr_som4_clm_loc(:) ! (gN/m3) vertically-resolved decomposing (litter, cwd, soil) n pools
#endif

    PetscScalar, pointer :: smin_no3_vr_clm_loc(:)           ! (gN/m3) vertically-resolved soil mineral NO3
    PetscScalar, pointer :: smin_nh4_vr_clm_loc(:)           ! (gN/m3) vertically-resolved soil mineral NH4
    PetscScalar, pointer :: smin_nh4sorb_vr_clm_loc(:)       ! (gN/m3) vertically-resolved soil mineral NH4 absorbed

    PetscErrorCode :: ierr
    !
    !------------------------------------------------------------------------------------------
    !
    associate ( &
       decomp_cpools_vr=> clm_bgc_data%decomp_cpools_vr_col     , &      ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
       decomp_npools_vr=> clm_bgc_data%decomp_npools_vr_col     , &      ! (gN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
       smin_no3_vr     => clm_bgc_data%smin_no3_vr_col          , &      ! (gN/m3) vertically-resolved soil mineral NO3
       smin_nh4_vr     => clm_bgc_data%smin_nh4_vr_col          , &      ! (gN/m3) vertically-resolved soil mineral NH4
       smin_nh4sorb_vr => clm_bgc_data%smin_nh4sorb_vr_col      , &       ! (gN/m3) vertically-resolved soil mineral NH4 absorbed

       decomp_ppools_vr=> clm_bgc_data%decomp_ppools_vr_col     , & ! [real(r8) (:,:,:) ! col (gP/m3) vertically-resolved decomposing (litter, cwd, soil) P pools
       solutionp_vr    => clm_bgc_data%solutionp_vr_col         , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil solution P
       labilep_vr      => clm_bgc_data%labilep_vr_col           , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil labile mineral P
       secondp_vr      => clm_bgc_data%secondp_vr_col           , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil secondary mineralP
       occlp_vr        => clm_bgc_data%occlp_vr_col             , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil occluded mineral P
       primp_vr        => clm_bgc_data%primp_vr_col             , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil primary mineral P
       sminp_vr        => clm_bgc_data%sminp_vr_col               & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil mineral P = solutionp + labilep + secondp
    )

#ifdef FLEXIBLE_POOLS
    call VecGetArrayF90(clm_pf_idata%decomp_cpools_vr_clmp, decomp_cpools_vr_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%decomp_npools_vr_clmp, decomp_npools_vr_clm_loc, ierr)
    CHKERRQ(ierr)
#else
    call VecGetArrayF90(clm_pf_idata%decomp_cpools_vr_lit1_clmp, decomp_cpools_vr_lit1_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%decomp_cpools_vr_lit2_clmp, decomp_cpools_vr_lit2_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%decomp_cpools_vr_lit3_clmp, decomp_cpools_vr_lit3_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%decomp_cpools_vr_cwd_clmp,  decomp_cpools_vr_cwd_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%decomp_cpools_vr_som1_clmp, decomp_cpools_vr_som1_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%decomp_cpools_vr_som2_clmp, decomp_cpools_vr_som2_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%decomp_cpools_vr_som3_clmp, decomp_cpools_vr_som3_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%decomp_cpools_vr_som4_clmp, decomp_cpools_vr_som4_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%decomp_npools_vr_lit1_clmp, decomp_npools_vr_lit1_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%decomp_npools_vr_lit2_clmp, decomp_npools_vr_lit2_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%decomp_npools_vr_lit3_clmp, decomp_npools_vr_lit3_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%decomp_npools_vr_cwd_clmp,  decomp_npools_vr_cwd_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%decomp_npools_vr_som1_clmp, decomp_npools_vr_som1_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%decomp_npools_vr_som2_clmp, decomp_npools_vr_som2_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%decomp_npools_vr_som3_clmp, decomp_npools_vr_som3_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%decomp_npools_vr_som4_clmp, decomp_npools_vr_som4_clm_loc, ierr)
    CHKERRQ(ierr)
#endif

    call VecGetArrayF90(clm_pf_idata%smin_no3_vr_clmp, smin_no3_vr_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%smin_nh4_vr_clmp, smin_nh4_vr_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%smin_nh4sorb_vr_clmp, smin_nh4sorb_vr_clm_loc, ierr)
    CHKERRQ(ierr)

#ifdef FLEXIBLE_POOLS
    decomp_cpools_vr_clm_loc(:) = 0._r8
    decomp_npools_vr_clm_loc(:) = 0._r8
#else
    decomp_cpools_vr_lit1_clm_loc(:) = 0._r8
    decomp_cpools_vr_lit2_clm_loc(:) = 0._r8
    decomp_cpools_vr_lit3_clm_loc(:) = 0._r8
    decomp_cpools_vr_cwd_clm_loc(:)  = 0._r8
    decomp_cpools_vr_som1_clm_loc(:) = 0._r8
    decomp_cpools_vr_som2_clm_loc(:) = 0._r8
    decomp_cpools_vr_som3_clm_loc(:) = 0._r8
    decomp_cpools_vr_som4_clm_loc(:) = 0._r8
    decomp_npools_vr_lit1_clm_loc(:) = 0._r8
    decomp_npools_vr_lit2_clm_loc(:) = 0._r8
    decomp_npools_vr_lit3_clm_loc(:) = 0._r8
    decomp_npools_vr_cwd_clm_loc(:)  = 0._r8
    decomp_npools_vr_som1_clm_loc(:) = 0._r8
    decomp_npools_vr_som2_clm_loc(:) = 0._r8
    decomp_npools_vr_som3_clm_loc(:) = 0._r8
    decomp_npools_vr_som4_clm_loc(:) = 0._r8
#endif

    smin_no3_vr_clm_loc(:)      = 0._r8
    smin_nh4_vr_clm_loc(:)      = 0._r8
    smin_nh4sorb_vr_clm_loc(:)  = 0._r8

    do fc = 1, num_soilc  ! will need to extend to multiple columns?
       c = filter_soilc(fc)

       if ( col_pp%wtgcell(c) <= 0._r8 .or. (.not.col_pp%active(c)) ) cycle     ! don't assign data to PF for inactive cell

       g       = col_pp%gridcell(c)
       wtgcell = col_pp%wtgcell(c)

       gcount = g - bounds%begg
       do j = 1, nlevdecomp

          ! note: all clm-pf soil layers are 'nzclm_mapped' for both TH/BGC,
          !       but in CLM, TH is within nlevsoi, bgc within 'nlevdecomp'

          cellcount = gcount*clm_pf_idata%nzclm_mapped+j

          if(j <= clm_pf_idata%nzclm_mapped) then

             do k = 1, ndecomp_pools
                 realc_gcell = decomp_cpools_vr(c,j,k) &
                              /clm_pf_idata%C_molecular_weight * wtgcell
                 realn_gcell = decomp_npools_vr(c,j,k) &
                              /clm_pf_idata%N_molecular_weight * wtgcell

#ifdef FLEXIBLE_POOLS
                 vec_offset = (k-1)*clm_pf_idata%ngclm_sub       ! decomp_pool vec: 'cell' first, then 'species'

                 decomp_cpools_vr_clm_loc(vec_offset+cellcount) = realc_gcell &
                         + decomp_cpools_vr_clm_loc(vec_offset+cellcount)
                 decomp_npools_vr_clm_loc(vec_offset+cellcount) = realn_gcell &
                         + decomp_npools_vr_clm_loc(vec_offset+cellcount)

#else
                 if (k==i_met_lit) then
                    decomp_cpools_vr_lit1_clm_loc(cellcount) = realc_gcell &
                         + decomp_cpools_vr_lit1_clm_loc(cellcount)
                    decomp_npools_vr_lit1_clm_loc(cellcount) = realn_gcell &
                         + decomp_npools_vr_lit1_clm_loc(cellcount)

                 elseif (k==i_cel_lit) then
                    decomp_cpools_vr_lit2_clm_loc(cellcount) = realc_gcell &
                         + decomp_cpools_vr_lit2_clm_loc(cellcount)
                    decomp_npools_vr_lit2_clm_loc(cellcount) = realn_gcell &
                         + decomp_npools_vr_lit2_clm_loc(cellcount)

                 elseif (k==i_lig_lit) then
                    decomp_cpools_vr_lit3_clm_loc(cellcount) = realc_gcell &
                         + decomp_cpools_vr_lit3_clm_loc(cellcount)
                    decomp_npools_vr_lit3_clm_loc(cellcount) = realn_gcell &
                         + decomp_npools_vr_lit3_clm_loc(cellcount)

                 elseif (k==i_cwd) then
                    decomp_cpools_vr_cwd_clm_loc(cellcount) = realc_gcell &
                         + decomp_cpools_vr_cwd_clm_loc(cellcount)
                    decomp_npools_vr_cwd_clm_loc(cellcount) = realn_gcell &
                         + decomp_npools_vr_cwd_clm_loc(cellcount)

                 else
                    isom = k-i_cwd
                    if (isom==1 .and. isom<=ndecomp_pools) then
                       decomp_cpools_vr_som1_clm_loc(cellcount) = realc_gcell &
                         + decomp_cpools_vr_som1_clm_loc(cellcount)
                       decomp_npools_vr_som1_clm_loc(cellcount) = realn_gcell &
                         + decomp_npools_vr_som1_clm_loc(cellcount)

                    elseif (isom==2 .and. isom<=ndecomp_pools) then
                       decomp_cpools_vr_som2_clm_loc(cellcount) = realc_gcell &
                         + decomp_cpools_vr_som2_clm_loc(cellcount)
                       decomp_npools_vr_som2_clm_loc(cellcount) = realn_gcell &
                         + decomp_npools_vr_som2_clm_loc(cellcount)

                    elseif (isom==3 .and. isom<=ndecomp_pools) then
                       decomp_cpools_vr_som3_clm_loc(cellcount) = realc_gcell &
                         + decomp_cpools_vr_som3_clm_loc(cellcount)
                       decomp_npools_vr_som3_clm_loc(cellcount) = realn_gcell &
                         + decomp_npools_vr_som3_clm_loc(cellcount)

                    elseif (isom==4 .and. isom<=ndecomp_pools) then       ! if using 'century' type, will end here
                       decomp_cpools_vr_som4_clm_loc(cellcount) = realc_gcell &
                         + decomp_cpools_vr_som4_clm_loc(cellcount)
                       decomp_npools_vr_som4_clm_loc(cellcount) = realn_gcell &
                         + decomp_npools_vr_som4_clm_loc(cellcount)
                    end if

                end if
#endif

             enddo ! do k=1, ndecomp_pools

             realn_gcell = smin_no3_vr(c,j)/clm_pf_idata%N_molecular_weight * wtgcell
             smin_no3_vr_clm_loc(cellcount) = realn_gcell &
                         + smin_no3_vr_clm_loc(cellcount)

             realn_gcell = smin_nh4_vr(c,j)/clm_pf_idata%N_molecular_weight * wtgcell
             smin_nh4_vr_clm_loc(cellcount) = realn_gcell &
                         + smin_nh4_vr_clm_loc(cellcount)

             realn_gcell = smin_nh4sorb_vr(c,j)/clm_pf_idata%N_molecular_weight * wtgcell
             smin_nh4sorb_vr_clm_loc(cellcount) = realn_gcell &
                         + smin_nh4sorb_vr_clm_loc(cellcount)

           endif

       enddo ! do j = 1, nlevdecomp

    enddo ! do c = begc, endc

#ifdef FLEXIBLE_POOLS
    call VecRestoreArrayF90(clm_pf_idata%decomp_cpools_vr_clmp, decomp_cpools_vr_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%decomp_npools_vr_clmp, decomp_npools_vr_clm_loc, ierr)
    CHKERRQ(ierr)
#else
!-----------------------------------------------------------------------------
!write(*,'(A,50(1h-))')">>>DEBUG | get_clm_bgc_conc,lev=1 for C & N"
!write(*,'(12A14)')"lit1","lit2","lit3","cwd","som1","som2","som3","som4","no3","nh4","nh4sorb"
!write(*,'(12E14.6)')decomp_cpools_vr(1,1,1:8)
!write(*,'(12E14.6)')decomp_npools_vr(1,1,1:8),smin_no3_vr(1,1),smin_nh4_vr(1,1),smin_nh4sorb_vr(1,1)
!-----------------------------------------------------------------------------

    call VecRestoreArrayF90(clm_pf_idata%decomp_cpools_vr_lit1_clmp, decomp_cpools_vr_lit1_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%decomp_cpools_vr_lit2_clmp, decomp_cpools_vr_lit2_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%decomp_cpools_vr_lit3_clmp, decomp_cpools_vr_lit3_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%decomp_cpools_vr_cwd_clmp,  decomp_cpools_vr_cwd_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%decomp_cpools_vr_som1_clmp, decomp_cpools_vr_som1_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%decomp_cpools_vr_som2_clmp, decomp_cpools_vr_som2_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%decomp_cpools_vr_som3_clmp, decomp_cpools_vr_som3_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%decomp_cpools_vr_som4_clmp, decomp_cpools_vr_som4_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%decomp_npools_vr_lit1_clmp, decomp_npools_vr_lit1_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%decomp_npools_vr_lit2_clmp, decomp_npools_vr_lit2_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%decomp_npools_vr_lit3_clmp, decomp_npools_vr_lit3_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%decomp_npools_vr_cwd_clmp,  decomp_npools_vr_cwd_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%decomp_npools_vr_som1_clmp, decomp_npools_vr_som1_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%decomp_npools_vr_som2_clmp, decomp_npools_vr_som2_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%decomp_npools_vr_som3_clmp, decomp_npools_vr_som3_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%decomp_npools_vr_som4_clmp, decomp_npools_vr_som4_clm_loc, ierr)
    CHKERRQ(ierr)
#endif
    call VecRestoreArrayF90(clm_pf_idata%smin_no3_vr_clmp, smin_no3_vr_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%smin_nh4_vr_clmp, smin_nh4_vr_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%smin_nh4sorb_vr_clmp, smin_nh4sorb_vr_clm_loc, ierr)
    CHKERRQ(ierr)

  end associate
  end subroutine get_clm_bgc_conc

  !-----------------------------------------------------------------------------
  !
  ! !IROUTINE: get_clm_bgc_rate()
  !
  ! !INTERFACE:
  subroutine get_clm_bgc_rate(clm_bgc_data,  &
          bounds, num_soilc, filter_soilc)
!! TODO: add phosphorus vars
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

    type(clm_bgc_interface_data_type), intent(in) :: clm_bgc_data

!    type(cnstate_type)       , intent(in) :: cnstate_vars
!    type(carbonflux_type)    , intent(in) :: carbonflux_vars
!    type(nitrogenflux_type)  , intent(in) :: nitrogenflux_vars

    character(len=256) :: subname = "get_clm_bgc_rate"

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

 ! !LOCAL VARIABLES:
    integer  :: fc, c, g, j, k                         ! do loop indices
    integer  :: gcount, cellcount
    real(r8) :: wtgcell, realc_gcell, realn_gcell

    real(r8) :: dtime                               ! land model time step (sec)

    ! ratios of NH4:NO3 in N deposition and fertilization (temporarily set here, will be as inputs)
!    real(r8) :: r_nh4_no3_dep(bounds%begc:bounds%endc)
!    real(r8) :: r_nh4_no3_fert(bounds%begc:bounds%endc)
!    real(r8) :: fnh4_dep, fnh4_fert

    ! C/N source/sink rates as inputs for pflotran: Units - moles/m3/s (note: do unit conversion here for input rates)
#ifdef FLEXIBLE_POOLS
    integer  :: vec_offset
    PetscScalar, pointer :: rate_decomp_c_clm_loc(:)       !
    PetscScalar, pointer :: rate_decomp_n_clm_loc(:)       !
#else
    integer  :: isom
    PetscScalar, pointer :: rate_lit1c_clm_loc(:)   !
    PetscScalar, pointer :: rate_lit2c_clm_loc(:)   !
    PetscScalar, pointer :: rate_lit3c_clm_loc(:)   !
    PetscScalar, pointer :: rate_cwdc_clm_loc(:)    !
    PetscScalar, pointer :: rate_som1c_clm_loc(:)   !
    PetscScalar, pointer :: rate_som2c_clm_loc(:)   !
    PetscScalar, pointer :: rate_som3c_clm_loc(:)   !
    PetscScalar, pointer :: rate_som4c_clm_loc(:)   !
    PetscScalar, pointer :: rate_lit1n_clm_loc(:)   !
    PetscScalar, pointer :: rate_lit2n_clm_loc(:)   !
    PetscScalar, pointer :: rate_lit3n_clm_loc(:)   !
    PetscScalar, pointer :: rate_cwdn_clm_loc(:)    !
    PetscScalar, pointer :: rate_som1n_clm_loc(:)   !
    PetscScalar, pointer :: rate_som2n_clm_loc(:)   !
    PetscScalar, pointer :: rate_som3n_clm_loc(:)   !
    PetscScalar, pointer :: rate_som4n_clm_loc(:)   !
#endif

    PetscScalar, pointer :: rate_plantndemand_clm_loc(:)   !
    PetscScalar, pointer :: rate_smin_no3_clm_loc(:)   !
    PetscScalar, pointer :: rate_smin_nh4_clm_loc(:)   !

    PetscErrorCode :: ierr

    !
    !---------------------------------------------------------------------------
    !
    associate ( &
      ! plant litering and removal + SOM/LIT vertical transport
      col_net_to_decomp_cpools_vr       => clm_bgc_data%externalc_to_decomp_cpools_col  , &
      col_net_to_decomp_npools_vr       => clm_bgc_data%externaln_to_decomp_npools_col  , &
      ! inorg. nitrogen source
!      ndep_to_sminn                  => nitrogenflux_vars%ndep_to_sminn_col                 , &
!      nfix_to_sminn                  => nitrogenflux_vars%nfix_to_sminn_col                 , &
!      fert_to_sminn                  => nitrogenflux_vars%fert_to_sminn_col                 , &
!      soyfixn_to_sminn               => nitrogenflux_vars%soyfixn_to_sminn_col              , &
!      supplement_to_sminn_vr         => nitrogenflux_vars%supplement_to_sminn_vr_col        , &
!      !
!      nfixation_prof                 => cnstate_vars%nfixation_prof_col                     , &
!      ndep_prof                      => cnstate_vars%ndep_prof_col                          , &
!      activeroot_prof                => cnstate_vars%activeroot_prof_col                    , &
      ! inorg. nitrogen sink (if not going to be done in PF)
      no3_net_transport_vr              => clm_bgc_data%no3_net_transport_vr_col        , &
      ! inorg. nitrogen sink potential
      col_plant_ndemand_vr              => clm_bgc_data%plant_ndemand_vr_col            , &

      externaln_to_nh4_vr               => clm_bgc_data%externaln_to_nh4_col            , &
      externaln_to_no3_vr               => clm_bgc_data%externaln_to_no3_col            , &

      col_net_to_decomp_ppools_vr       => clm_bgc_data%externalp_to_decomp_ppools_col  , &
      externalp_to_primp_vr             => clm_bgc_data%externalp_to_primp_col          , &
      externalp_to_labilep_vr           => clm_bgc_data%externalp_to_labilep_col        , &
      externalp_to_solutionp            => clm_bgc_data%externalp_to_solutionp_col      , &
      sminp_net_transport_vr            => clm_bgc_data%sminp_net_transport_vr_col      , &
      col_plant_pdemand_vr              => clm_bgc_data%plant_pdemand_vr_col              &
    )

    dtime = get_step_size()

#ifdef FLEXIBLE_POOLS
    call VecGetArrayF90(clm_pf_idata%rate_decomp_c_clmp, rate_decomp_c_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%rate_decomp_n_clmp, rate_decomp_n_clm_loc, ierr)
    CHKERRQ(ierr)
#else
    call VecGetArrayF90(clm_pf_idata%rate_lit1c_clmp, rate_lit1c_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%rate_lit2c_clmp, rate_lit2c_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%rate_lit3c_clmp, rate_lit3c_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%rate_cwdc_clmp, rate_cwdc_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%rate_som1c_clmp, rate_som1c_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%rate_som2c_clmp, rate_som2c_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%rate_som3c_clmp, rate_som3c_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%rate_som4c_clmp, rate_som4c_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%rate_lit1n_clmp, rate_lit1n_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%rate_lit2n_clmp, rate_lit2n_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%rate_lit3n_clmp, rate_lit3n_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%rate_cwdn_clmp, rate_cwdn_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%rate_som1n_clmp, rate_som1n_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%rate_som2n_clmp, rate_som2n_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%rate_som3n_clmp, rate_som3n_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%rate_som4n_clmp, rate_som4n_clm_loc, ierr)
#endif

    call VecGetArrayF90(clm_pf_idata%rate_plantndemand_clmp, rate_plantndemand_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%rate_smin_no3_clmp, rate_smin_no3_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%rate_smin_nh4_clmp, rate_smin_nh4_clm_loc, ierr)
    CHKERRQ(ierr)

    ! Initialize to ZERO
#ifdef FLEXIBLE_POOLS
    rate_decomp_c_clm_loc(:) = 0.0_r8
    rate_decomp_n_clm_loc(:) = 0.0_r8
#else
    rate_lit1c_clm_loc(:) = 0.0_r8
    rate_lit2c_clm_loc(:) = 0.0_r8
    rate_lit3c_clm_loc(:) = 0.0_r8
    rate_cwdc_clm_loc (:) = 0.0_r8
    rate_som1c_clm_loc(:) = 0.0_r8
    rate_som2c_clm_loc(:) = 0.0_r8
    rate_som3c_clm_loc(:) = 0.0_r8
    rate_som4c_clm_loc(:) = 0.0_r8
    rate_lit1n_clm_loc(:) = 0.0_r8
    rate_lit2n_clm_loc(:) = 0.0_r8
    rate_lit3n_clm_loc(:) = 0.0_r8
    rate_cwdn_clm_loc (:) = 0.0_r8
    rate_som1n_clm_loc(:) = 0.0_r8
    rate_som2n_clm_loc(:) = 0.0_r8
    rate_som3n_clm_loc(:) = 0.0_r8
    rate_som4n_clm_loc(:) = 0.0_r8
#endif
    rate_smin_no3_clm_loc(:) = 0.0_r8
    rate_smin_nh4_clm_loc(:) = 0.0_r8
    rate_plantndemand_clm_loc(:) = 0.0_r8

!    r_nh4_no3_dep(:)  = 1.0_r8      ! temporarily assuming half of N dep is in NH4 and another half in NO3
!    r_nh4_no3_fert(:) = 1.0_r8      ! temporarily assiming half of N fertilization is in NH4 and another half in NO3

    do fc = 1,num_soilc
       c = filter_soilc(fc)
       if ( col_pp%wtgcell(c) <= 0._r8 .or. (.not.col_pp%active(c)) ) cycle     ! don't assign data from PF for inactive cell

       g = col_pp%gridcell(c)
       gcount = g - bounds%begg
       wtgcell = col_pp%wtgcell(c)

       do j = 1, nlevdecomp
          cellcount = gcount*clm_pf_idata%nzclm_mapped+j

          ! note: all clm-pf soil layers are 'clm_pf_idata%nzclm_mapped' for both TH/BGC,
          ! but in CLM, TH is within clm_pf_idata%nzclm_mapped, bgc within 'nlevdecomp'
          if(j <= clm_pf_idata%nzclm_mapped) then

              do k = 1, ndecomp_pools
                ! need more checking here: how to weight column data onto grid (F.-M. Yuan)
                 realc_gcell = col_net_to_decomp_cpools_vr(c,j,k) &
                              /clm_pf_idata%C_molecular_weight * wtgcell
                 realn_gcell = col_net_to_decomp_npools_vr(c,j,k) &
                              /clm_pf_idata%N_molecular_weight * wtgcell

#ifdef FLEXIBLE_POOLS
                 vec_offset = (k-1)*clm_pf_idata%ngclm_sub       ! decomp_pool vec: 'cell' first, then 'species'

                 rate_decomp_c_clm_loc(vec_offset+cellcount) = realc_gcell &
                         + rate_decomp_c_clm_loc(vec_offset+cellcount)
                 rate_decomp_n_clm_loc(vec_offset+cellcount) = realn_gcell &
                         + rate_decomp_n_clm_loc(vec_offset+cellcount)
#else
                 if (k==i_met_lit) then
                    rate_lit1c_clm_loc(cellcount) = realc_gcell &
                         + rate_lit1c_clm_loc(cellcount)
                    rate_lit1n_clm_loc(cellcount) = realn_gcell &
                         + rate_lit1n_clm_loc(cellcount)
                 else if(k==i_cel_lit) then
                    rate_lit2c_clm_loc(cellcount) = realc_gcell &
                         + rate_lit2c_clm_loc(cellcount)
                    rate_lit2n_clm_loc(cellcount) = realn_gcell &
                         + rate_lit2n_clm_loc(cellcount)
                 else if(k==i_lig_lit) then
                    rate_lit3c_clm_loc(cellcount) = realc_gcell &
                         + rate_lit3c_clm_loc(cellcount)
                    rate_lit3n_clm_loc(cellcount) = realn_gcell &
                         + rate_lit3n_clm_loc(cellcount)
                 else if(k==i_cwd) then
                    rate_cwdc_clm_loc(cellcount) = realc_gcell &
                         + rate_cwdc_clm_loc(cellcount)
                    rate_cwdn_clm_loc(cellcount) = realn_gcell &
                         + rate_cwdn_clm_loc(cellcount)
                 !
                 else
                    isom = k-i_cwd
                    if(isom==1 .and. isom<=ndecomp_pools) then
                       rate_som1c_clm_loc(cellcount) = realc_gcell &
                         + rate_som1c_clm_loc(cellcount)
                       rate_som1n_clm_loc(cellcount) = realn_gcell &
                         + rate_som1n_clm_loc(cellcount)
                    else if(isom==2 .and. isom<=ndecomp_pools) then
                       rate_som2c_clm_loc(cellcount) = realc_gcell &
                         + rate_som2c_clm_loc(cellcount)
                       rate_som2n_clm_loc(cellcount) = realn_gcell &
                         + rate_som2n_clm_loc(cellcount)
                    else if(isom==3 .and. isom<=ndecomp_pools) then
                       rate_som3c_clm_loc(cellcount) = realc_gcell &
                         + rate_som3c_clm_loc(cellcount)
                       rate_som3n_clm_loc(cellcount) = realn_gcell &
                         + rate_som3n_clm_loc(cellcount)
                    else if(isom==4 .and. isom<=ndecomp_pools) then       ! if using 'century' type, will end here
                       rate_som4c_clm_loc(cellcount) = realc_gcell &
                         + rate_som4c_clm_loc(cellcount)
                       rate_som4n_clm_loc(cellcount) = realn_gcell &
                         + rate_som4n_clm_loc(cellcount)
                    end if

                 endif
#endif
              enddo ! do k=1, ndecomp_pools

!              fnh4_dep  = max(0._r8, min(1.0_r8, 1._r8/(r_nh4_no3_dep(c)+1._r8)))
!              fnh4_fert = max(0._r8, min(1.0_r8, 1._r8/(r_nh4_no3_fert(c)+1._r8)))

!              realn_gcell = &
!                        ( fnh4_dep*ndep_to_sminn(c) * ndep_prof(c, j) +  &
!                          fnh4_fert*fert_to_sminn(c) * ndep_prof(c, j) + &
!                          fnh4_fert*supplement_to_sminn_vr(c,j) +        &
!                          nfix_to_sminn(c) * nfixation_prof(c, j) +      &
!                          soyfixn_to_sminn(c) * nfixation_prof(c, j)    &
!                         )/ clm_pf_idata%N_molecular_weight * wtgcell

              realn_gcell = externaln_to_nh4_vr(c,j)/ clm_pf_idata%N_molecular_weight * wtgcell
              rate_smin_nh4_clm_loc(cellcount) = realn_gcell + rate_smin_nh4_clm_loc(cellcount)

!              realn_gcell = &
!                         ( (1._r8-fnh4_dep)*ndep_to_sminn(c) * ndep_prof(c, j) +  &
!                           (1._r8-fnh4_fert)*fert_to_sminn(c) * ndep_prof(c, j) + &
!                           (1._r8-fnh4_fert)*supplement_to_sminn_vr(c,j) &
!                         )/ clm_pf_idata%N_molecular_weight * wtgcell

              realn_gcell = externaln_to_no3_vr(c,j)/ clm_pf_idata%N_molecular_weight * wtgcell
              ! PF hydrological mode is OFF, then NO3 transport NOT to calculate in PF
              ! then it's done in CLM, so need to pass those to PF as source/sink term (RT mass transfer)
              if(.not.pf_hmode) then
                realn_gcell = realn_gcell - (no3_net_transport_vr(c,j) &
                         )/ clm_pf_idata%N_molecular_weight * wtgcell
              endif
              rate_smin_no3_clm_loc(cellcount) = realn_gcell + rate_smin_no3_clm_loc(cellcount)

              ! plant N uptake rate here IS the N demand (potential uptake)
              realn_gcell = col_plant_ndemand_vr(c,j) &
                           /clm_pf_idata%N_molecular_weight * wtgcell
              rate_plantndemand_clm_loc(cellcount) = rate_plantndemand_clm_loc(cellcount) + realn_gcell

#ifdef CLM_PF_DEBUG
      write(pflotran_m%option%myrank+200,*) 'checking bgc-mass-rate - clm: ', &
        'rank=',pflotran_m%option%myrank, 'column=',c, 'layer_id=',j, &
        'rate_nh4_clm(layer_id)=',rate_smin_nh4_clm_loc(cellcount)
#endif

          endif ! if (j<=clm_pf_idata%nzclm_mapped)
       enddo ! do j=1, nlevdecomp
    enddo ! do fc=1,numsoic

#ifdef FLEXIBLE_POOLS
    call VecRestoreArrayF90(clm_pf_idata%rate_decomp_c_clmp, rate_decomp_c_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%rate_decomp_n_clmp, rate_decomp_n_clm_loc, ierr)
    CHKERRQ(ierr)
#else

!-----------------------------------------------------------------------------
!write(*,'(A,50(1h-))')">>>DEBUG | get_clm_bgc_rate,lev=1 for C & N"
!write(*,'(12A14)')"lit1","lit2","lit3","cwd","som1","som2","som3","som4","no3","nh4","plantNdemand"
!write(*,'(12E14.6)')col_net_to_decomp_cpools_vr(1,1,1:ndecomp_pools)
!write(*,'(12E14.6)')col_net_to_decomp_npools_vr(1,1,1:ndecomp_pools),&
!                (rate_smin_no3_clm_loc(1))*clm_pf_idata%N_molecular_weight, &
!                (rate_smin_nh4_clm_loc(1))*clm_pf_idata%N_molecular_weight, &
!                (rate_plantndemand_clm_loc(1))*clm_pf_idata%N_molecular_weight
!-----------------------------------------------------------------------------

    call VecRestoreArrayF90(clm_pf_idata%rate_lit1c_clmp, rate_lit1c_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%rate_lit2c_clmp, rate_lit2c_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%rate_lit3c_clmp, rate_lit3c_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%rate_cwdc_clmp, rate_cwdc_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%rate_som1c_clmp, rate_som1c_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%rate_som2c_clmp, rate_som2c_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%rate_som3c_clmp, rate_som3c_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%rate_som4c_clmp, rate_som4c_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%rate_lit1n_clmp, rate_lit1n_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%rate_lit2n_clmp, rate_lit2n_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%rate_lit3n_clmp, rate_lit3n_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%rate_cwdn_clmp, rate_cwdn_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%rate_som1n_clmp, rate_som1n_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%rate_som2n_clmp, rate_som2n_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%rate_som3n_clmp, rate_som3n_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%rate_som4n_clmp, rate_som4n_clm_loc, ierr)
#endif

    call VecRestoreArrayF90(clm_pf_idata%rate_plantndemand_clmp, rate_plantndemand_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%rate_smin_no3_clmp, rate_smin_no3_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%rate_smin_nh4_clmp, rate_smin_nh4_clm_loc, ierr)
    CHKERRQ(ierr)

    end associate
  end subroutine get_clm_bgc_rate


  ! =============================UPDATE PFLOTRAN evolving variables to CLM ==========================================


  !-----------------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: update_soil_moisture_pf2clm
  !
  ! !INTERFACE:
  subroutine update_soil_moisture_pf2clm(clm_bgc_data,     &
           bounds, num_soilc, filter_soilc)

  !
  ! !DESCRIPTION:
  !
  !
  ! !USES:

  ! !ARGUMENTS:
    implicit none

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_soilc        ! number of column soil points in column filter
    integer, intent(in) :: filter_soilc(:)  ! column filter for soil points
!     type(soilstate_type) , intent(in)    :: soilstate_vars
!     type(waterstate_type), intent(inout) :: waterstate_vars
    type(clm_bgc_interface_data_type), intent(inout) :: clm_bgc_data

  ! !LOCAL VARIABLES:
    integer  :: fc, c, j, g, gcount      ! indices

    PetscScalar, pointer :: sat_ice_clm_loc(:)
    PetscScalar, pointer :: sat_clm_loc(:)
    PetscScalar, pointer :: watsat_clm_loc(:)
    PetscErrorCode :: ierr

  !EOP
  !-----------------------------------------------------------------------
    associate ( &
         gridcell   =>  col_pp%gridcell   , & !  [integer (:)]  gridcell index of column
         wtgcell    =>  col_pp%wtgcell    , & !  [real(r8) (:)]  weight (relative to gridcell
         dz         =>  col_pp%dz         , & !  [real(r8) (:,:)]  layer thickness depth (m)
         !
         h2osoi_liq_col =>  clm_bgc_data%h2osoi_liq_col      , &
         h2osoi_ice_col =>  clm_bgc_data%h2osoi_ice_col      , &
         h2osoi_vol_col =>  clm_bgc_data%h2osoi_vol_col      , &
         watsat_col     =>  clm_bgc_data%watsat_col      &
    )

    call VecGetArrayF90(clm_pf_idata%soillsat_clms, sat_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%effporosity_clms, watsat_clm_loc, ierr)
    CHKERRQ(ierr)
    if (pf_frzmode) &
       call VecGetArrayF90(clm_pf_idata%soilisat_clms, sat_ice_clm_loc, ierr)
    CHKERRQ(ierr)

    do fc = 1,num_soilc
      c = filter_soilc(fc)
      if ( col_pp%wtgcell(c) <= 0._r8 .or. (.not.col_pp%active(c)) ) cycle     ! don't assign data from PF for inactive cell

      g = col_pp%gridcell(c)
      gcount = g - bounds%begg
      do j = 1, nlevsoi

        if (j<=clm_pf_idata%nzclm_mapped) then
          h2osoi_liq_col(c,j) = sat_clm_loc(gcount*nlevsoi + j) * &
                              watsat_clm_loc(gcount*nlevsoi + j) * dz(c,j) * denh2o    ! 'watsat_clm_loc' may be effective porosity
          if (pf_frzmode) then
             h2osoi_ice_col(c,j) = sat_ice_clm_loc(gcount*nlevsoi + j) * &
                              watsat_clm_loc(gcount*nlevsoi + j) * dz(c,j) * denice
          end if

        else
           h2osoi_liq_col(c,j) = h2osoi_liq_col(c,clm_pf_idata%nzclm_mapped)
          if (pf_frzmode) then
             h2osoi_ice_col(c,j) = h2osoi_ice_col(c,clm_pf_idata%nzclm_mapped)
          end if

        end if

        h2osoi_vol_col(c,j) = h2osoi_liq_col(c,j) / dz(c,j) / denh2o + &
                              h2osoi_ice_col(c,j) / dz(c,j) / denice
        h2osoi_vol_col(c,j) = min(h2osoi_vol_col(c,j), watsat_col(c,j))

      enddo

    enddo

    call VecRestoreArrayF90(clm_pf_idata%soillsat_clms, sat_clm_loc, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%effporosity_clms, watsat_clm_loc, ierr)
    CHKERRQ(ierr)
    if (pf_frzmode) &
       call VecRestoreArrayF90(clm_pf_idata%soilisat_clms, sat_ice_clm_loc, ierr)
    CHKERRQ(ierr)

    end associate
  end subroutine update_soil_moisture_pf2clm

  !-----------------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: update_soil_temperature_pf2clm
  !
  ! !INTERFACE:
  subroutine update_soil_temperature_pf2clm(clm_bgc_data,     &
           bounds, num_soilc, filter_soilc)

  !
  ! !DESCRIPTION:
  !
  !
  ! !USES:

  ! !ARGUMENTS:
    implicit none

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_soilc        ! number of column soil points in column filter
    integer, intent(in) :: filter_soilc(:)  ! column filter for soil points
!     type(soilstate_type) , intent(in)    :: soilstate_vars
!     type(temperature_type), intent(inout):: temperature_vars
    type(clm_bgc_interface_data_type), intent(inout) :: clm_bgc_data

  ! !LOCAL VARIABLES:
    integer  :: fc, c, j, g, gcount      ! indices

    PetscScalar, pointer :: soilt_clms_loc(:)
    PetscErrorCode :: ierr

  !EOP
  !-----------------------------------------------------------------------
    associate ( &
         gridcell   =>  col_pp%gridcell   , & !  [integer (:)]  gridcell index of column
         wtgcell    =>  col_pp%wtgcell    , & !  [real(r8) (:)]  weight (relative to gridcell
         dz         =>  col_pp%dz         , & !  [real(r8) (:,:)]  layer thickness depth (m)
         !
         t_soisno   => clm_bgc_data%t_soisno_col   & ! snow-soil temperature (Kelvin)
    )

    call VecGetArrayReadF90(clm_pf_idata%soilt_clms, soilt_clms_loc, ierr)
    CHKERRQ(ierr)

    do fc = 1,num_soilc
      c = filter_soilc(fc)
      if ( col_pp%wtgcell(c) <= 0._r8 .or. (.not.col_pp%active(c)) ) cycle     ! don't assign data from PF for inactive cell

      g = col_pp%gridcell(c)
      gcount = g - bounds%begg
      do j = 1, nlevsoi

        if (j<=clm_pf_idata%nzclm_mapped) then
             t_soisno(c,j) = soilt_clms_loc(gcount*clm_pf_idata%nzclm_mapped+j) + tfrz

        else
             t_soisno(c,j) = t_soisno(c, clm_pf_idata%nzclm_mapped)

        end if

      enddo

    enddo

    call VecRestoreArrayReadF90(clm_pf_idata%soilt_clms, soilt_clms_loc, ierr)
    CHKERRQ(ierr)

    end associate
  end subroutine update_soil_temperature_pf2clm

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
  subroutine update_soil_bgc_pf2clm(clm_bgc_data,   &
           bounds, num_soilc, filter_soilc)
!! TODO: add phosphorus vars
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
!    type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars
!    type(ch4_type)           , intent(inout) :: ch4_vars
    type(clm_bgc_interface_data_type), intent(inout) :: clm_bgc_data

    character(len=256) :: subname = "update_soil_bgc_pf2clm"

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

    integer  :: fc,c,g,j,k
    integer  :: gcount, cellcount
    real(r8) :: wtgcell

    real(r8) :: dtime            ! land model time step (sec)

#ifdef FLEXIBLE_POOLS
    integer  :: vec_offset
    PetscScalar, pointer :: decomp_cpools_vr_clm_loc(:)      ! (moleC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
    PetscScalar, pointer :: decomp_npools_vr_clm_loc(:)      ! (moleN/m3) vertically-resolved decomposing (litter, cwd, soil) N pools
#else
    integer  :: isom
    PetscScalar, pointer :: decomp_cpools_vr_lit1_clm_loc(:) ! (moleC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
    PetscScalar, pointer :: decomp_cpools_vr_lit2_clm_loc(:) ! (moleC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
    PetscScalar, pointer :: decomp_cpools_vr_lit3_clm_loc(:) ! (moleC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
    PetscScalar, pointer :: decomp_cpools_vr_cwd_clm_loc(:)  ! (moleC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
    PetscScalar, pointer :: decomp_cpools_vr_som1_clm_loc(:) ! (moleC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
    PetscScalar, pointer :: decomp_cpools_vr_som2_clm_loc(:) ! (moleC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
    PetscScalar, pointer :: decomp_cpools_vr_som3_clm_loc(:) ! (moleC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
    PetscScalar, pointer :: decomp_cpools_vr_som4_clm_loc(:) ! (moleC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
    PetscScalar, pointer :: decomp_npools_vr_lit1_clm_loc(:) ! (moleN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
    PetscScalar, pointer :: decomp_npools_vr_lit2_clm_loc(:) ! (moleN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
    PetscScalar, pointer :: decomp_npools_vr_lit3_clm_loc(:) ! (moleN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
    PetscScalar, pointer :: decomp_npools_vr_cwd_clm_loc(:)  ! (moleN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
#endif

    PetscScalar, pointer :: smin_no3_vr_clm_loc(:)           ! (moleN/m3) vertically-resolved soil mineral NO3
    PetscScalar, pointer :: smin_nh4_vr_clm_loc(:)           ! (moleN/m3) vertically-resolved total soil mineral NH4
    PetscScalar, pointer :: smin_nh4sorb_vr_clm_loc(:)       ! (moleN/m3) vertically-resolved absorbed soil mineral NH4

    ! 'accextrn_vr' - accumulative (root) extracted N, i.e., actual plant N uptake from each soil layer, within a CLM timestep
    PetscScalar, pointer :: accextrnh4_vr_clm_loc(:)         ! (moleN/m3/timestep) vertically-resolved soil mineral N root-extraction (accumulated)
    PetscScalar, pointer :: accextrno3_vr_clm_loc(:)         ! (moleN/m3/timestep) vertically-resolved soil mineral N root-extraction (accumulated)

    ! 'accnmin_vr' - accumulative gross N mineralization within a CLM timestep
    PetscScalar, pointer :: accnmin_vr_clm_loc(:)               ! (moleN/m3/timestep) vertically-resolved soil N mineralization (accumulated)

    ! 'accnimm_vr' - accumulative N immobilization within a CLM timestep
    PetscScalar, pointer :: accnimmp_vr_clm_loc(:)              ! (moleN/m3/timestep) vertically-resolved soil N potential immoblilization (accumulated)
    PetscScalar, pointer :: accnimm_vr_clm_loc(:)               ! (moleN/m3/timestep) vertically-resolved soil N immoblilization (accumulated)

    PetscErrorCode :: ierr

!------------------------------------------------------------------------------------
     !
     associate ( &
     initial_cn_ratio             => clm_bgc_data%initial_cn_ratio             , &
     decomp_cpools_vr             => clm_bgc_data%decomp_cpools_vr_col         , &
     decomp_npools_vr             => clm_bgc_data%decomp_npools_vr_col         , &
     sminn_vr                     => clm_bgc_data%sminn_vr_col                 , &
     smin_no3_vr                  => clm_bgc_data%smin_no3_vr_col              , &
     smin_nh4_vr                  => clm_bgc_data%smin_nh4_vr_col              , &
     smin_nh4sorb_vr              => clm_bgc_data%smin_nh4sorb_vr_col          , &
     decomp_cpools_delta_vr       => clm_bgc_data%decomp_cpools_sourcesink_col , &
     decomp_npools_delta_vr       => clm_bgc_data%decomp_npools_sourcesink_col  , &
     
     sminn_to_plant_vr            => clm_bgc_data%sminn_to_plant_vr_col         , &
     smin_no3_to_plant_vr         => clm_bgc_data%smin_no3_to_plant_vr_col      , &
     smin_nh4_to_plant_vr         => clm_bgc_data%smin_nh4_to_plant_vr_col      , &
     potential_immob_vr           => clm_bgc_data%potential_immob_vr_col        , &
     actual_immob_vr              => clm_bgc_data%actual_immob_vr_col           , &
     gross_nmin_vr                => clm_bgc_data%gross_nmin_vr_col               &
     )
! ------------------------------------------------------------------------
     dtime = get_step_size()

     ! soil C/N pool increments set to the previous timestep (i.e., not yet updated)
     decomp_cpools_delta_vr  = 0._r8-decomp_cpools_vr
     decomp_npools_delta_vr  = 0._r8-decomp_npools_vr

     ! clm-pf interface data updated
#ifdef FLEXIBLE_POOLS
     call VecGetArrayReadF90(clm_pf_idata%decomp_cpools_vr_clms, decomp_cpools_vr_clm_loc, ierr)
     CHKERRQ(ierr)
     call VecGetArrayReadF90(clm_pf_idata%decomp_npools_vr_clms, decomp_npools_vr_clm_loc, ierr)
     CHKERRQ(ierr)
#else
     call VecGetArrayReadF90(clm_pf_idata%decomp_cpools_vr_lit1_clms, decomp_cpools_vr_lit1_clm_loc, ierr)
     CHKERRQ(ierr)
     call VecGetArrayReadF90(clm_pf_idata%decomp_cpools_vr_lit2_clms, decomp_cpools_vr_lit2_clm_loc, ierr)
     CHKERRQ(ierr)
     call VecGetArrayReadF90(clm_pf_idata%decomp_cpools_vr_lit3_clms, decomp_cpools_vr_lit3_clm_loc, ierr)
     CHKERRQ(ierr)
     call VecGetArrayReadF90(clm_pf_idata%decomp_cpools_vr_cwd_clms,  decomp_cpools_vr_cwd_clm_loc, ierr)
     CHKERRQ(ierr)
     call VecGetArrayReadF90(clm_pf_idata%decomp_cpools_vr_som1_clms, decomp_cpools_vr_som1_clm_loc, ierr)
     CHKERRQ(ierr)
     call VecGetArrayReadF90(clm_pf_idata%decomp_cpools_vr_som2_clms, decomp_cpools_vr_som2_clm_loc, ierr)
     CHKERRQ(ierr)
     call VecGetArrayReadF90(clm_pf_idata%decomp_cpools_vr_som3_clms, decomp_cpools_vr_som3_clm_loc, ierr)
     CHKERRQ(ierr)
     call VecGetArrayReadF90(clm_pf_idata%decomp_cpools_vr_som4_clms, decomp_cpools_vr_som4_clm_loc, ierr)
     CHKERRQ(ierr)
     call VecGetArrayReadF90(clm_pf_idata%decomp_npools_vr_lit1_clms, decomp_npools_vr_lit1_clm_loc, ierr)
     CHKERRQ(ierr)
     call VecGetArrayReadF90(clm_pf_idata%decomp_npools_vr_lit2_clms, decomp_npools_vr_lit2_clm_loc, ierr)
     CHKERRQ(ierr)
     call VecGetArrayReadF90(clm_pf_idata%decomp_npools_vr_lit3_clms, decomp_npools_vr_lit3_clm_loc, ierr)
     CHKERRQ(ierr)
     call VecGetArrayReadF90(clm_pf_idata%decomp_npools_vr_cwd_clms,  decomp_npools_vr_cwd_clm_loc, ierr)
     CHKERRQ(ierr)
#endif

     call VecGetArrayReadF90(clm_pf_idata%smin_no3_vr_clms, smin_no3_vr_clm_loc, ierr)
     CHKERRQ(ierr)
     call VecGetArrayReadF90(clm_pf_idata%smin_nh4_vr_clms, smin_nh4_vr_clm_loc, ierr)
     CHKERRQ(ierr)
     call VecGetArrayReadF90(clm_pf_idata%smin_nh4sorb_vr_clms, smin_nh4sorb_vr_clm_loc, ierr)
     CHKERRQ(ierr)

     call VecGetArrayReadF90(clm_pf_idata%accextrnh4_vr_clms, accextrnh4_vr_clm_loc, ierr)
     CHKERRQ(ierr)
     call VecGetArrayReadF90(clm_pf_idata%accextrno3_vr_clms, accextrno3_vr_clm_loc, ierr)
     CHKERRQ(ierr)
     call VecGetArrayReadF90(clm_pf_idata%accnmin_vr_clms, accnmin_vr_clm_loc, ierr)
     CHKERRQ(ierr)
     call VecGetArrayReadF90(clm_pf_idata%accnimmp_vr_clms, accnimmp_vr_clm_loc, ierr)
     CHKERRQ(ierr)
     call VecGetArrayReadF90(clm_pf_idata%accnimm_vr_clms, accnimm_vr_clm_loc, ierr)
     CHKERRQ(ierr)

    ! soil C/N pool increments, and actual plant N uptake, gross N mineralization and immobilization
    do fc = 1,num_soilc
       ! only operating on soil column, which then back to CLM-CN
       ! (TODO) NOT YET do columna-level down-scaling from PF's grid-cell variables

       c = filter_soilc(fc)
       g = col_pp%gridcell(c)
       wtgcell = col_pp%wtgcell(c)
       if ( col_pp%wtgcell(c) <= 0._r8 .or. (.not.col_pp%active(c)) ) cycle     ! don't assign data from PF for inactive cell

       gcount = g - bounds%begg
       do j = 1, nlevdecomp

          cellcount = gcount*clm_pf_idata%nzclm_mapped+j
          if(j <= clm_pf_idata%nzclm_mapped) then

              ! updates the 'decomp_pools' src/sink terms,
              ! which then used in 'CNSoilLittVertTranspMod.F90' to update the pools
              do k=1, ndecomp_pools
#ifdef FLEXIBLE_POOLS
                 vec_offset = (k-1)*clm_pf_idata%nlclm_sub       ! decomp_pool vec: 'cell' first, then 'species'

                 decomp_cpools_delta_vr(c,j,k) = decomp_cpools_delta_vr(c,j,k)    &
                                 + decomp_cpools_vr_clm_loc(vec_offset+cellcount) &
                                 * clm_pf_idata%C_molecular_weight

                 if (clm_pf_idata%floating_cn_ratio(k)) then
                    decomp_npools_delta_vr(c,j,k) = decomp_npools_delta_vr(c,j,k) &
                                 + decomp_npools_vr_clm_loc(vec_offset+cellcount) &
                                 * clm_pf_idata%N_molecular_weight

                 else
                    decomp_npools_delta_vr(c,j,k) = decomp_cpools_delta_vr(c,j,k) &
                                                        /initial_cn_ratio(k)
                 endif

#else
                 if (k==i_met_lit) then
                    decomp_cpools_delta_vr(c,j,k) = decomp_cpools_delta_vr(c,j,k) &
                                       + decomp_cpools_vr_lit1_clm_loc(cellcount) &
                                       * clm_pf_idata%C_molecular_weight
                    decomp_npools_delta_vr(c,j,k) = decomp_npools_delta_vr(c,j,k) &
                                       + decomp_npools_vr_lit1_clm_loc(cellcount) &
                                       * clm_pf_idata%N_molecular_weight
                 else if (k==i_cel_lit) then
                    decomp_cpools_delta_vr(c,j,k) = decomp_cpools_delta_vr(c,j,k) &
                                       + decomp_cpools_vr_lit2_clm_loc(cellcount) &
                                       * clm_pf_idata%C_molecular_weight
                    decomp_npools_delta_vr(c,j,k) = decomp_npools_delta_vr(c,j,k) &
                                       + decomp_npools_vr_lit2_clm_loc(cellcount) &
                                       * clm_pf_idata%N_molecular_weight
                 else if (k==i_lig_lit) then
                    decomp_cpools_delta_vr(c,j,k) = decomp_cpools_delta_vr(c,j,k) &
                                       + decomp_cpools_vr_lit3_clm_loc(cellcount) &
                                       * clm_pf_idata%C_molecular_weight
                    decomp_npools_delta_vr(c,j,k) = decomp_npools_delta_vr(c,j,k) &
                                       + decomp_npools_vr_lit3_clm_loc(cellcount) &
                                       * clm_pf_idata%N_molecular_weight
                 else if (k==i_cwd) then
                    decomp_cpools_delta_vr(c,j,k) = decomp_cpools_delta_vr(c,j,k) &
                                       + decomp_cpools_vr_cwd_clm_loc(cellcount)  &
                                       * clm_pf_idata%C_molecular_weight
                    decomp_npools_delta_vr(c,j,k) = decomp_npools_delta_vr(c,j,k) &
                                       + decomp_npools_vr_cwd_clm_loc(cellcount)  &
                                       * clm_pf_idata%N_molecular_weight
                !
                else

                    isom = k-i_cwd

                    if (isom==1 .and. isom<=ndecomp_pools) then
                       decomp_cpools_delta_vr(c,j,k) = decomp_cpools_delta_vr(c,j,k) &
                                       + decomp_cpools_vr_som1_clm_loc(cellcount)    &
                                       * clm_pf_idata%C_molecular_weight
                    elseif (isom==2 .and. isom<=ndecomp_pools) then
                       decomp_cpools_delta_vr(c,j,k) = decomp_cpools_delta_vr(c,j,k) &
                                       + decomp_cpools_vr_som2_clm_loc(cellcount)    &
                                       * clm_pf_idata%C_molecular_weight
                    elseif (isom==3 .and. isom<=ndecomp_pools) then
                       decomp_cpools_delta_vr(c,j,k) = decomp_cpools_delta_vr(c,j,k) &
                                       + decomp_cpools_vr_som3_clm_loc(cellcount)    &
                                       * clm_pf_idata%C_molecular_weight
                    elseif (isom==4 .and. isom<=ndecomp_pools) then
                       decomp_cpools_delta_vr(c,j,k) = decomp_cpools_delta_vr(c,j,k) &
                                       + decomp_cpools_vr_som4_clm_loc(cellcount)    &
                                       * clm_pf_idata%C_molecular_weight
                    end if

                   if (isom>0) then
                       decomp_npools_delta_vr(c,j,k) = decomp_cpools_delta_vr(c,j,k) &
                                                        /initial_cn_ratio(k)
                   end if

                endif

#endif
                 if (abs(decomp_cpools_delta_vr(c,j,k))<=1.d-20) decomp_cpools_delta_vr(c,j,k)=0._r8
                 if (abs(decomp_npools_delta_vr(c,j,k))<=1.d-21) decomp_npools_delta_vr(c,j,k)=0._r8

              enddo

              ! directly update the 'smin' N pools (SO, must bypass the 'CNNStateUpdate1,2,3' relevant to soil N)
              smin_no3_vr(c,j) = &
                           smin_no3_vr_clm_loc(cellcount)*clm_pf_idata%N_molecular_weight

              smin_nh4_vr(c,j) = &
                           smin_nh4_vr_clm_loc(cellcount)*clm_pf_idata%N_molecular_weight

              smin_nh4sorb_vr(c,j) = &
                           smin_nh4sorb_vr_clm_loc(cellcount)*clm_pf_idata%N_molecular_weight

              sminn_vr(c,j) = smin_no3_vr(c,j) + smin_nh4_vr(c,j) + smin_nh4sorb_vr(c,j)

              ! flows or changes (unit: g/m3/s)
              smin_nh4_to_plant_vr(c,j) = (accextrnh4_vr_clm_loc(cellcount) &
                                         * clm_pf_idata%N_molecular_weight)/dtime
              smin_no3_to_plant_vr(c,j) = (accextrno3_vr_clm_loc(cellcount) &
                                         * clm_pf_idata%N_molecular_weight)/dtime
              sminn_to_plant_vr(c,j) = smin_nh4_to_plant_vr(c,j) + smin_no3_to_plant_vr(c,j)

              gross_nmin_vr(c,j) = (accnmin_vr_clm_loc(cellcount) &
                                         * clm_pf_idata%N_molecular_weight)/dtime

              potential_immob_vr(c,j) = (accnimmp_vr_clm_loc(cellcount) &
                                         * clm_pf_idata%N_molecular_weight)/dtime
              actual_immob_vr(c,j)    = (accnimm_vr_clm_loc(cellcount)  &
                                         * clm_pf_idata%N_molecular_weight)/dtime

          else    ! just in case 'clm_pf_idata%nzclm_mapped<nlevdcomp'

              do k=1, ndecomp_pools
                 decomp_cpools_delta_vr(c,j,k) = 0._r8
                 decomp_npools_delta_vr(c,j,k) = 0._r8
              enddo

              sminn_vr(c,j)        = 0._r8
              smin_no3_vr(c,j)     = 0._r8
              smin_nh4_vr(c,j)     = 0._r8
              smin_nh4sorb_vr(c,j) = 0._r8

              smin_nh4_to_plant_vr(c,j)  = 0._r8
              smin_no3_to_plant_vr(c,j)  = 0._r8
              sminn_to_plant_vr(c,j)     = 0._r8
              gross_nmin_vr(c,j)         = 0._r8
              potential_immob_vr(c,j)    = 0._r8
              actual_immob_vr(c,j)       = 0._r8

          endif

       enddo
     enddo ! do c = 1, numsoilc
!-----------------------------------------------------------------------------
!write(*,'(A30,12E14.6)')"DEBUG | pf UPDATE no3=",smin_no3_vr(1,1:nlevdecomp)
!write(*,'(A30,12E14.6)')"DEBUG | pf UPDATE nh4=",smin_nh4_vr(1,1:nlevdecomp)
!write(*,'(A30,12E14.6)')"DEBUG | pf UPDATE nh4sorb=",smin_nh4sorb_vr(1,1:nlevdecomp)
!-----------------------------------------------------------------------------

#ifdef FLEXIBLE_POOLS
     call VecRestoreArrayReadF90(clm_pf_idata%decomp_cpools_vr_clms, decomp_cpools_vr_clm_loc, ierr)
     CHKERRQ(ierr)
     call VecRestoreArrayReadF90(clm_pf_idata%decomp_npools_vr_clms, decomp_npools_vr_clm_loc, ierr)
     CHKERRQ(ierr)
#else
     call VecRestoreArrayReadF90(clm_pf_idata%decomp_cpools_vr_lit1_clms, decomp_cpools_vr_lit1_clm_loc, ierr)
     CHKERRQ(ierr)
     call VecRestoreArrayReadF90(clm_pf_idata%decomp_cpools_vr_lit2_clms, decomp_cpools_vr_lit2_clm_loc, ierr)
     CHKERRQ(ierr)
     call VecRestoreArrayReadF90(clm_pf_idata%decomp_cpools_vr_lit3_clms, decomp_cpools_vr_lit3_clm_loc, ierr)
     CHKERRQ(ierr)
     call VecRestoreArrayReadF90(clm_pf_idata%decomp_cpools_vr_cwd_clms,  decomp_cpools_vr_cwd_clm_loc, ierr)
     CHKERRQ(ierr)
     call VecRestoreArrayReadF90(clm_pf_idata%decomp_cpools_vr_som1_clms, decomp_cpools_vr_som1_clm_loc, ierr)
     CHKERRQ(ierr)
     call VecRestoreArrayReadF90(clm_pf_idata%decomp_cpools_vr_som2_clms, decomp_cpools_vr_som2_clm_loc, ierr)
     CHKERRQ(ierr)
     call VecRestoreArrayReadF90(clm_pf_idata%decomp_cpools_vr_som3_clms, decomp_cpools_vr_som3_clm_loc, ierr)
     CHKERRQ(ierr)
     call VecRestoreArrayReadF90(clm_pf_idata%decomp_cpools_vr_som4_clms, decomp_cpools_vr_som4_clm_loc, ierr)
     CHKERRQ(ierr)
     call VecRestoreArrayReadF90(clm_pf_idata%decomp_npools_vr_lit1_clms, decomp_npools_vr_lit1_clm_loc, ierr)
     CHKERRQ(ierr)
     call VecRestoreArrayReadF90(clm_pf_idata%decomp_npools_vr_lit2_clms, decomp_npools_vr_lit2_clm_loc, ierr)
     CHKERRQ(ierr)
     call VecRestoreArrayReadF90(clm_pf_idata%decomp_npools_vr_lit3_clms, decomp_npools_vr_lit3_clm_loc, ierr)
     CHKERRQ(ierr)
     call VecRestoreArrayReadF90(clm_pf_idata%decomp_npools_vr_cwd_clms,  decomp_npools_vr_cwd_clm_loc, ierr)
     CHKERRQ(ierr)
#endif

     call VecRestoreArrayF90(clm_pf_idata%smin_no3_vr_clms, smin_no3_vr_clm_loc, ierr)
     CHKERRQ(ierr)
     call VecRestoreArrayF90(clm_pf_idata%smin_nh4_vr_clms, smin_nh4_vr_clm_loc, ierr)
     CHKERRQ(ierr)
     call VecRestoreArrayF90(clm_pf_idata%smin_nh4sorb_vr_clms, smin_nh4sorb_vr_clm_loc, ierr)
     CHKERRQ(ierr)

     call VecRestoreArrayF90(clm_pf_idata%accextrnh4_vr_clms, accextrnh4_vr_clm_loc, ierr)
     CHKERRQ(ierr)
     call VecRestoreArrayF90(clm_pf_idata%accextrno3_vr_clms, accextrno3_vr_clm_loc, ierr)
     CHKERRQ(ierr)
     call VecRestoreArrayF90(clm_pf_idata%accnmin_vr_clms, accnmin_vr_clm_loc, ierr)
     CHKERRQ(ierr)
     call VecRestoreArrayF90(clm_pf_idata%accnimmp_vr_clms, accnimmp_vr_clm_loc, ierr)
     CHKERRQ(ierr)
     call VecRestoreArrayF90(clm_pf_idata%accnimm_vr_clms, accnimm_vr_clm_loc, ierr)
     CHKERRQ(ierr)

     ! update bgc gas losses
     call update_bgc_gaslosses_pf2clm(clm_bgc_data, &
        bounds, num_soilc, filter_soilc)

    end associate
  end subroutine update_soil_bgc_pf2clm

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
  subroutine update_bgc_gaslosses_pf2clm(clm_bgc_data,  &
     bounds, num_soilc, filter_soilc)

     use clm_time_manager, only : get_step_size, get_nstep

     !
     implicit none

     type(bounds_type) , intent(in)  :: bounds
     integer, intent(in) :: num_soilc       ! number of soil columns in filter
     integer, intent(in) :: filter_soilc(:) ! filter for soil columns

!     type(atm2lnd_type)       , intent(in) :: clm_a2l
!     type(waterstate_type)    , intent(in) :: waterstate_vars
!     type(carbonflux_type)    , intent(inout) :: carbonflux_vars
!     type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars
     type(clm_bgc_interface_data_type), intent(inout) :: clm_bgc_data

     !character(len=256) :: subname = "get_pf_bgc_gaslosses"

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

     integer  :: fc, c, g, j
     integer  :: gcount, cellcount
     real(r8) :: dtime            ! land model time step (sec)
     integer  :: nstep

     ! for gas species
     real(r8) :: tc, tk, total_p       ! temperature (oC, K), total air pressure (Pa)
     real(r8) :: co2_p, n2_p, n2o_p    ! partial pressure (Pa) of CO2, N2, N2O
     real(r8) :: cgas, cgas_p          ! mole-C(N)/m3(bulk soil)
     real(r8) :: air_vol, air_molar, wfps
     integer :: lair_barrier(1:num_soilc)      ! toppest soil layer that little air space for air flow into deep soil (-1: no, 0: ground, >0: soil layer)

     ! gases from PFLOTRAN are timely accumulated, so gas fluxes are calculated here if over atm. partial pressure (no explicit transport available from PF now)
     PetscScalar, pointer :: gco2_vr_clms_loc(:)               ! (M: molC/m3 bulk soil) vertically-resolved soil gas CO2 from PF's evolution
     PetscScalar, pointer :: gn2_vr_clms_loc(:)                ! (M: molN/m3 bulk soil) vertically-resolved soil gas N2 from PF's evolution
     PetscScalar, pointer :: gn2o_vr_clms_loc(:)               ! (M: molN/m3 bulk soil) vertically-resolved soil gas N2O from PF's evolution
     PetscScalar, pointer :: gco2_vr_clmp_loc(:)               ! (M: molC/m3 bulk soil) vertically-resolved soil gas CO2 to reset PF's CO2g
     PetscScalar, pointer :: gn2_vr_clmp_loc(:)                ! (M: molN/m3 bulk soil) vertically-resolved soil gas N2 to reset PF's N2g
     PetscScalar, pointer :: gn2o_vr_clmp_loc(:)               ! (M: molN/m3 bulk soil) vertically-resolved soil gas N2O to reset PF's N2Og

     ! 'acchr_vr' - accumulative CO2 proudction from decompositon (for tracking HR, not involving mass-balance)
     PetscScalar, pointer :: acchr_vr_clm_loc(:)               ! (moleC/m3/timestep) vertically-resolved soil CO2 production (accumulated)
     ! 'accngasmin_vr' - accumulative N gas proudction from mineralization (for tracking, not involving mass-balance)
     PetscScalar, pointer :: accngasmin_vr_clm_loc(:)          ! (moleN/m3/timestep) vertically-resolved soil gaseous N  production (accumulated)
     ! 'accngasnitr_vr' - accumulative N gas proudction from nitrification (for tracking)
     PetscScalar, pointer :: accngasnitr_vr_clm_loc(:)         ! (moleN/m3/timestep) vertically-resolved soil gaseous N  production (accumulated)
     ! 'accngasdeni_vr' - accumulative N gas proudction from denitrification (for tracking)
     PetscScalar, pointer :: accngasdeni_vr_clm_loc(:)         ! (moleN/m3/timestep) vertically-resolved soil gaseous N  production (accumulated)

     !
     PetscScalar, pointer :: soillsat_clm_loc(:)
     PetscScalar, pointer :: soilisat_clm_loc(:)
     PetscScalar, pointer :: soilpor_clm_loc(:)
     PetscScalar, pointer :: soilt_clm_loc(:)
     PetscScalar, pointer :: soilpress_clm_loc(:)
     PetscErrorCode :: ierr

     real(r8), parameter :: rgas = 8.3144621       ! m3 Pa K-1 mol-1

!------------------------------------------------------------------------------------
    associate ( &
     forc_pco2                    => clm_bgc_data%forc_pco2_grc                  , & ! partial pressure co2 (Pa)
     forc_pch4                    => clm_bgc_data%forc_pch4_grc                  , & ! partial pressure ch4 (Pa)
     forc_pbot                    => clm_bgc_data%forc_pbot_not_downscaled_grc   , & ! atmospheric pressure (Pa)
     frac_sno                     => clm_bgc_data%frac_sno_eff_col       , & ! fraction of ground covered by snow (0 to 1)
     frac_h2osfc                  => clm_bgc_data%frac_h2osfc_col        , & ! fraction of ground covered by surface water (0 to 1)
     dz                           => clm_bgc_data%dz                                 , & ! soil layer thickness depth (m)
     hr_vr                        => clm_bgc_data%hr_vr_col              , &
     f_co2_soil_vr                => clm_bgc_data%f_co2_soil_vr_col      , &
     f_n2o_soil_vr                => clm_bgc_data%f_n2o_soil_vr_col    , &
     f_n2_soil_vr                 => clm_bgc_data%f_n2_soil_vr_col     , &
     f_ngas_decomp_vr             => clm_bgc_data%f_ngas_decomp_vr_col , &
     f_ngas_nitri_vr              => clm_bgc_data%f_ngas_nitri_vr_col  , &
     f_ngas_denit_vr              => clm_bgc_data%f_ngas_denit_vr_col    &
     )
! ------------------------------------------------------------------------
     dtime = get_step_size()
     nstep = get_nstep()

     ! get the current time-step state variables of aq. phase of interested species
     call VecGetArrayF90(clm_pf_idata%gco2_vr_clms, gco2_vr_clms_loc, ierr)
     CHKERRQ(ierr)
     call VecGetArrayF90(clm_pf_idata%gn2_vr_clms, gn2_vr_clms_loc, ierr)
     CHKERRQ(ierr)
     call VecGetArrayF90(clm_pf_idata%gn2o_vr_clms, gn2o_vr_clms_loc, ierr)
     CHKERRQ(ierr)

     call VecGetArrayF90(clm_pf_idata%gco2_vr_clmp, gco2_vr_clmp_loc, ierr)
     CHKERRQ(ierr)
     call VecGetArrayF90(clm_pf_idata%gn2_vr_clmp, gn2_vr_clmp_loc, ierr)
     CHKERRQ(ierr)
     call VecGetArrayF90(clm_pf_idata%gn2o_vr_clmp, gn2o_vr_clmp_loc, ierr)
     CHKERRQ(ierr)

     call VecGetArrayF90(clm_pf_idata%acchr_vr_clms, acchr_vr_clm_loc, ierr)
     CHKERRQ(ierr)
     call VecGetArrayF90(clm_pf_idata%accngasmin_vr_clms, accngasmin_vr_clm_loc, ierr)
     CHKERRQ(ierr)
     call VecGetArrayF90(clm_pf_idata%accngasnitr_vr_clms, accngasnitr_vr_clm_loc, ierr)
     CHKERRQ(ierr)
     call VecGetArrayF90(clm_pf_idata%accngasdeni_vr_clms, accngasdeni_vr_clm_loc, ierr)
     CHKERRQ(ierr)

     ! env. variables to properties of gases
     if (pf_tmode) then
         call VecGetArrayF90(clm_pf_idata%soilt_clms, soilt_clm_loc, ierr)
         CHKERRQ(ierr)     ! PF evolved 'soilt'
     else
         call VecGetArrayF90(clm_pf_idata%soilt_clmp, soilt_clm_loc, ierr)
         CHKERRQ(ierr)   ! CLM evolved 'soilt' - for CLM, MPI vecs and Seq. vecs should be same
     end if

     if (pf_frzmode) then
         call VecGetArrayF90(clm_pf_idata%soilisat_clms, soilisat_clm_loc, ierr)
         CHKERRQ(ierr)     ! PF evolved 'soil ice saturation'
     end if

     if (pf_hmode) then
         call VecGetArrayF90(clm_pf_idata%soillsat_clms, soillsat_clm_loc, ierr)
         CHKERRQ(ierr)     ! PF evolved 'soil liq. saturation'
         call VecGetArrayF90(clm_pf_idata%press_clms, soilpress_clm_loc, ierr)
         CHKERRQ(ierr)     ! PF evolved 'soil liq. saturation'
     else
         call VecGetArrayF90(clm_pf_idata%soillsat_clmp, soillsat_clm_loc, ierr)
         CHKERRQ(ierr)! CLM evolved 'soilt liq. saturation'
         call VecGetArrayF90(clm_pf_idata%press_clmp, soilpress_clm_loc, ierr)
         CHKERRQ(ierr)! CLM evolved 'soilt liq. saturation'
     endif
     call VecGetArrayF90(clm_pf_idata%effporosity_clms, soilpor_clm_loc, ierr)  ! PF evolved 'soil porosity'
     CHKERRQ(ierr)

     ! find the toppest air barrier layer
     lair_barrier(:) = -1            ! (-1: no barrier, 0: ground snow/ice/water-layer barrier, >=1: barrier in soil column)
     do fc = 1,num_soilc
       c = filter_soilc(fc)
       if ( col_pp%wtgcell(c) <= 0._r8 .or. (.not.col_pp%active(c)) ) cycle     ! don't assign data from PF for inactive cell

       g = col_pp%gridcell(c)
       gcount = g - bounds%begg

       if((frac_sno(c)+ frac_h2osfc(c))>=0.90_r8) then
         lair_barrier(c) = 0
       endif

       do j = 1, clm_pf_idata%nzclm_mapped

          cellcount = gcount*clm_pf_idata%nzclm_mapped+j

          wfps = 0._r8
          if (lair_barrier(c) >= 0) exit

          if (pf_frzmode) then
             wfps =soillsat_clm_loc(cellcount)  &
                    + soilisat_clm_loc(cellcount)
          else
             wfps =soillsat_clm_loc(cellcount)   ! note: 'lsat' from PF has been adjusted by 'isat' reduced porosity
          endif
          if (wfps > 0.95_r8) then             ! 95% total saturation as a critical-point for air-flow into deep soil
              lair_barrier(c) = j
          endif
       enddo
    enddo

    !
    do fc = 1,num_soilc          ! operating on soil column, which then back to CLM-CN
       c = filter_soilc(fc)
       if ( col_pp%wtgcell(c) <= 0._r8 .or. (.not.col_pp%active(c)) ) cycle     ! don't assign data from PF for inactive cell

       g = col_pp%gridcell(c)
       gcount = g - bounds%begg

       total_p = forc_pbot(g)

       do j = 1, nlevdecomp
          cellcount = gcount*clm_pf_idata%nzclm_mapped+j
          if(j <= clm_pf_idata%nzclm_mapped) then

              tc = soilt_clm_loc(cellcount)    ! soil layer tc (oC)
              tk = tc+tfrz

              ! total_p is soil air pressure (soil water pressure if not less than atm. pressure)
              total_p = max(total_p, soilpress_clm_loc(cellcount))

              ! the following is for adjusting air space in soil (seems not right ?? -- off)
              !air_vol = (1.0_r8 - soillsat_clm_loc(cellcount)) * &
              !                    soilpor_clm_loc(cellcount)                 ! m3 air/m3 soil
              !air_vol = max(air_vol, 0.01d0)            ! min. 0.01 to avoid math. issue

              air_vol = 1.0_r8                          ! atm used if not commented out
              air_molar = total_p*air_vol/rgas/tk       ! moles of air in a cell

              ! gas fluxes from  immobile PFLOTRAN evolving CO2imm, N2Oimm and N2imm, which are cumulative
              ! CO2 -
              cgas = gco2_vr_clms_loc(cellcount)       ! mol/m3 soil (evolving in PF, but not yet transport)
              co2_p = forc_pco2(g)                              ! assuming atm. pco2 (pa) as directly equilibrated with soil CO2(g)
              cgas_p = co2_p/forc_pbot(g) * air_molar

              f_co2_soil_vr(c,j) = cgas-cgas_p
              cgas = cgas - f_co2_soil_vr(c,j)
              if (j <= lair_barrier(c) .or. lair_barrier(c) < 0) then     ! above barrier OR no-barrier(-1)
                 gco2_vr_clmp_loc(cellcount) = cgas_p  ! this refreshed-air will pass back to PF
              else
                 gco2_vr_clmp_loc(cellcount) = cgas    ! currently don't have air transport (TODO)
              endif

              f_co2_soil_vr(c,j) = f_co2_soil_vr(c,j)*clm_pf_idata%C_molecular_weight                  ! moleCO2/m3 --> gC/m3 soil
              f_co2_soil_vr(c,j) = f_co2_soil_vr(c,j)/dtime                               ! gC/m3/s

              ! N2
              cgas = gn2_vr_clms_loc(cellcount)       ! mol/m3 soil (evolving in PF, but not yet transport)
              n2_p = 0.78084_r8                                ! assuming atm. pn2 as directly equilibrated with soil n2(g)
              cgas_p = n2_p * air_molar           ! moleN2/m3

              f_n2_soil_vr(c,j) = cgas-cgas_p
              cgas = cgas - f_n2_soil_vr(c,j)
              if (j <= lair_barrier(c) .or. lair_barrier(c) < 0) then     ! above barrier OR no-barrier(-1)
                 gn2_vr_clmp_loc(cellcount) = cgas_p  ! this refreshed-air will pass back to PF
              else
                 gn2_vr_clmp_loc(cellcount) = cgas    ! currently don't have air transport (TODO)
              endif

              f_n2_soil_vr(c,j) = f_n2_soil_vr(c,j)*clm_pf_idata%N_molecular_weight*2._r8    ! mole-N2/m3 --> g-N/m3 soil
              f_n2_soil_vr(c,j) = f_n2_soil_vr(c,j)/dtime                 ! gN/m3/s

              ! N2O
              cgas = gn2o_vr_clms_loc(cellcount)
              n2o_p = 310e-9_r8                   ! assuming general atm. pN2O (310ppbv in 1990) as directly equilibrated with soil N2(aq)
              cgas_p = n2o_p * air_molar          ! moleN2O/m3

              f_n2o_soil_vr(c,j) = cgas-cgas_p
              cgas = cgas - f_n2o_soil_vr(c,j)
              if (j <= lair_barrier(c) .or. lair_barrier(c) < 0) then     ! above barrier OR no-barrier(-1)
                 gn2o_vr_clmp_loc(cellcount) = cgas_p  ! this refreshed-air will pass back to PF
              else
                 gn2o_vr_clmp_loc(cellcount) = cgas    ! currently don't have air transport (TODO)
              endif

              f_n2o_soil_vr(c,j) = f_n2o_soil_vr(c,j)*clm_pf_idata%N_molecular_weight*2._r8   ! mole-N2O/m3 --> g-N/m3 soil
              f_n2o_soil_vr(c,j) = f_n2o_soil_vr(c,j)/dtime                ! gN/m3/s

              ! tracking HR from SOM-C reaction network
              hr_vr(c,j)           = (acchr_vr_clm_loc(cellcount) &
                                    * clm_pf_idata%C_molecular_weight)/dtime

              ! tracking gaseous N production from N reaction network
              f_ngas_decomp_vr(c,j)= (accngasmin_vr_clm_loc (cellcount) &
                                    * clm_pf_idata%N_molecular_weight)/dtime

              f_ngas_nitri_vr(c,j) = (accngasnitr_vr_clm_loc(cellcount) &
                                    * clm_pf_idata%N_molecular_weight)/dtime

              f_ngas_denit_vr(c,j) = (accngasdeni_vr_clm_loc(cellcount) &
                                    * clm_pf_idata%N_molecular_weight)/dtime

          else    ! just in case 'clm_pf_idata%nzclm_mapped<nlevdecomp'

              f_co2_soil_vr(c,j)         = f_co2_soil_vr(c,clm_pf_idata%nzclm_mapped)
              f_n2_soil_vr(c,j)          = f_n2_soil_vr(c,clm_pf_idata%nzclm_mapped)
              f_n2o_soil_vr(c,j)         = f_n2o_soil_vr(c,clm_pf_idata%nzclm_mapped)

              hr_vr(c,j)                 = 0._r8
              f_ngas_decomp_vr(c,j)      = 0._r8
              f_ngas_nitri_vr(c,j)       = 0._r8
              f_ngas_denit_vr(c,j)       = 0._r8

          endif
       enddo
     enddo ! do c = begc, endc

     call VecRestoreArrayF90(clm_pf_idata%gco2_vr_clms, gco2_vr_clms_loc, ierr)
     CHKERRQ(ierr)
     call VecRestoreArrayF90(clm_pf_idata%gn2_vr_clms, gn2_vr_clms_loc, ierr)
     CHKERRQ(ierr)
     call VecRestoreArrayF90(clm_pf_idata%gn2o_vr_clms, gn2o_vr_clms_loc, ierr)
     CHKERRQ(ierr)
     call VecRestoreArrayF90(clm_pf_idata%gco2_vr_clmp, gco2_vr_clmp_loc, ierr)
     CHKERRQ(ierr)
     call VecRestoreArrayF90(clm_pf_idata%gn2_vr_clmp, gn2_vr_clmp_loc, ierr)
     CHKERRQ(ierr)
     call VecRestoreArrayF90(clm_pf_idata%gn2o_vr_clmp, gn2o_vr_clmp_loc, ierr)
     CHKERRQ(ierr)

     call VecRestoreArrayF90(clm_pf_idata%acchr_vr_clms, acchr_vr_clm_loc, ierr)
     CHKERRQ(ierr)
     call VecRestoreArrayF90(clm_pf_idata%accngasmin_vr_clms, accngasmin_vr_clm_loc, ierr)
     CHKERRQ(ierr)
     call VecRestoreArrayF90(clm_pf_idata%accngasnitr_vr_clms, accngasnitr_vr_clm_loc, ierr)
     CHKERRQ(ierr)
     call VecRestoreArrayF90(clm_pf_idata%accngasdeni_vr_clms, accngasdeni_vr_clm_loc, ierr)
     CHKERRQ(ierr)

     if (pf_tmode) then
         call VecRestoreArrayF90(clm_pf_idata%soilt_clms, soilt_clm_loc, ierr)
         CHKERRQ(ierr)
     else
         call VecRestoreArrayF90(clm_pf_idata%soilt_clmp, soilt_clm_loc, ierr)
         CHKERRQ(ierr)  ! for CLM, MPI vecs and Seq. vecs should be same
     end if
     if (pf_frzmode) then
         call VecRestoreArrayF90(clm_pf_idata%soilisat_clms, soilisat_clm_loc, ierr)
         CHKERRQ(ierr)
     end if

     if (pf_hmode) then
         call VecRestoreArrayF90(clm_pf_idata%soillsat_clms, soillsat_clm_loc, ierr)
         CHKERRQ(ierr)     ! PF evolved 'soil liq. saturation'
         call VecRestoreArrayF90(clm_pf_idata%press_clms, soilpress_clm_loc, ierr)
         CHKERRQ(ierr)     ! PF evolved 'soil liq. saturation'
     else
         call VecRestoreArrayF90(clm_pf_idata%soillsat_clmp, soillsat_clm_loc, ierr)
         CHKERRQ(ierr)! CLM evolved 'soilt liq. saturation'
         call VecRestoreArrayF90(clm_pf_idata%press_clmp, soilpress_clm_loc, ierr)
         CHKERRQ(ierr)! CLM evolved 'soilt liq. saturation'
     endif
     call VecRestoreArrayF90(clm_pf_idata%effporosity_clms, soilpor_clm_loc, ierr)
     CHKERRQ(ierr)

     ! need to reset the PF's internal gas concentration (CLM ==> PF)
     call pflotranModelUpdateAqGasesfromCLM(pflotran_m)

    end associate
  end subroutine update_bgc_gaslosses_pf2clm


  !-----------------------------------------------------------------------------
  !BOP
  ! comment out this subroutine currently for bgc-only
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
!    clandunit         =>    col_pp%landunit          , & ! column's landunit
!    ltype             =>    lun_pp%itype             , & ! landunit type
!    zi                =>    col_pp%zi                , & ! Input: (:,:) soil layer interface depth (m)
!    dz                =>    col_pp%dz                , & ! Input: (:,:) soil layer thickness (m)
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
!      g = col_pp%gridcell(c)
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

#endif
!----------------------------------------------------------------------------------------------
! END of CLM-PFLOTRAN bgc coupling interface
!----------------------------------------------------------------------------------------------

end module clm_pflotran_interfaceMod

