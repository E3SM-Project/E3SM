module FatesRestartInterfaceMod


  use FatesConstantsMod , only : r8 => fates_r8
  use FatesConstantsMod , only : fates_avg_flag_length
  use FatesConstantsMod , only : fates_short_string_length
  use FatesConstantsMod , only : fates_long_string_length
  use FatesGlobals      , only : fates_log
  use FatesGlobals      , only : endrun => fates_endrun
  use FatesIODimensionsMod, only : fates_io_dimension_type
  use FatesIOVariableKindMod, only : fates_io_variable_kind_type
  use FatesRestartVariableMod, only : fates_restart_variable_type
  use FatesInterfaceMod, only : bc_in_type

  ! CIME GLOBALS
  use shr_log_mod       , only : errMsg => shr_log_errMsg


  implicit none

  ! ------------------------------------------------------------
  ! A note on variable naming conventions.
  ! Many variables in this restart IO portion of the code will
  ! follow the conventions:
  !
  ! <use_case>_<description>_<dimension>
  !
  ! For instance we use an index for restart variable "ir_"
  ! to point the object that contains the number of patches per
  ! site "npatch" and this value is relevant to all sites "si"
  ! thus:   ir_npatch_si
  !
  ! We also use associations to the data arrays of restart IO
  ! variables "rio", for example the leaf litter "leaf_litter"
  ! is retrieved for every patch and every functional type "paft"
  ! thus: rio_leaf_litter_paft
  !
  ! si: site dimension
  ! pa: patch dimension
  ! co: cohort dimension
  ! ft: functional type dimension
  ! cl: canopy layer dimension (upper, lower, etc)
  ! ls: layer sublayer dimension (fine discretization of upper,lower)
  ! wm: the number of memory slots for water (currently 10)
  ! -------------------------------------------------------------
  
  
  ! Indices to the restart variable object
  integer, private :: ir_npatch_si 
  integer, private :: ir_oldstock_si
  integer, private :: ir_cd_status_si
  integer, private :: ir_dd_status_si
  integer, private :: ir_nchill_days_si
  integer, private :: ir_leafondate_si
  integer, private :: ir_leafoffdate_si
  integer, private :: ir_dleafondate_si
  integer, private :: ir_dleafoffdate_si
  integer, private :: ir_acc_ni_si
  integer, private :: ir_gdd_si
  integer, private :: ir_nep_timeintegrated_si
  integer, private :: ir_npp_timeintegrated_si
  integer, private :: ir_hr_timeintegrated_si
  integer, private :: ir_cbal_error_fates_si
  integer, private :: ir_cbal_error_bgc_si
  integer, private :: ir_cbal_error_total_si
  integer, private :: ir_totecosysc_old_si
  integer, private :: ir_totfatesc_old_si
  integer, private :: ir_totbgcc_old_si
  integer, private :: ir_fates_to_bgc_this_ts_si
  integer, private :: ir_fates_to_bgc_last_ts_si
  integer, private :: ir_seedrainflux_si
  integer, private :: ir_trunk_product_si
  integer, private :: ir_ncohort_pa
  integer, private :: ir_balive_co
  integer, private :: ir_bdead_co
  integer, private :: ir_bleaf_co
  integer, private :: ir_broot_co
  integer, private :: ir_bstore_co
  integer, private :: ir_canopy_layer_co
  integer, private :: ir_canopy_layer_yesterday_co
  integer, private :: ir_canopy_trim_co
  integer, private :: ir_dbh_co
  integer, private :: ir_height_co
  integer, private :: ir_laimemory_co
  integer, private :: ir_leaf_md_co
  integer, private :: ir_root_md_co
  integer, private :: ir_nplant_co
  integer, private :: ir_gpp_acc_co
  integer, private :: ir_npp_acc_co
  integer, private :: ir_gpp_acc_hold_co
  integer, private :: ir_npp_acc_hold_co
  integer, private :: ir_npp_leaf_co
  integer, private :: ir_npp_froot_co
  integer, private :: ir_npp_sw_co
  integer, private :: ir_npp_dead_co
  integer, private :: ir_npp_seed_co
  integer, private :: ir_npp_store_co
  integer, private :: ir_bmort_co
  integer, private :: ir_hmort_co
  integer, private :: ir_cmort_co
  integer, private :: ir_imort_co
  integer, private :: ir_fmort_co

   !Logging
  integer, private :: ir_lmort_logging_co
  integer, private :: ir_lmort_collateral_co
  integer, private :: ir_lmort_infra_co


  integer, private :: ir_ddbhdt_co
  integer, private :: ir_dbalivedt_co
  integer, private :: ir_dbdeaddt_co
  integer, private :: ir_dbstoredt_co
  integer, private :: ir_resp_tstep_co
  integer, private :: ir_pft_co
  integer, private :: ir_status_co
  integer, private :: ir_isnew_co
  integer, private :: ir_cwd_ag_pacw
  integer, private :: ir_cwd_bg_pacw
  integer, private :: ir_leaf_litter_paft
  integer, private :: ir_root_litter_paft
  integer, private :: ir_leaf_litter_in_paft
  integer, private :: ir_root_litter_in_paft
  integer, private :: ir_seed_bank_sift
  integer, private :: ir_spread_si
  integer, private :: ir_livegrass_pa
  integer, private :: ir_age_pa
  integer, private :: ir_area_pa
  integer, private :: ir_fsun_paclftls
  integer, private :: ir_fabd_sun_paclftls
  integer, private :: ir_fabi_sun_paclftls
  integer, private :: ir_fabd_sha_paclftls
  integer, private :: ir_fabi_sha_paclftls
  integer, private :: ir_watermem_siwm

  ! The number of variable dim/kind types we have defined (static)
  integer, parameter :: fates_restart_num_dimensions = 2   !(cohort,column)
  integer, parameter :: fates_restart_num_dim_kinds = 4    !(cohort-int,cohort-r8,site-int,site-r8)

  ! integer constants for storing logical data
  integer, parameter :: old_cohort = 0
  integer, parameter :: new_cohort = 1  

  ! Local debug flag
  logical, parameter :: DEBUG=.false.

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

  ! This structure is allocated by thread, and must be calculated after the FATES
  ! sites are allocated, and their mapping to the HLM is identified.  This structure
  ! is not combined with iovar_bounds, because that one is multi-instanced.  This
  ! structure is used more during the update phase, wherease _bounds is used
  ! more for things like flushing
  type restart_map_type
     integer, allocatable :: site_index(:)   ! maps site indexes to the HIO site position
     integer, allocatable :: cohort1_index(:) ! maps site index to the HIO cohort 1st position
  end type restart_map_type



  type, public :: fates_restart_interface_type

     type(fates_restart_variable_type),allocatable :: rvars(:)
     integer,private :: num_restart_vars_

     ! Instanteate one registry of the different dimension/kinds (dk)
     ! All output variables will have a pointer to one of these dk's
     type(fates_io_variable_kind_type) :: dim_kinds(fates_restart_num_dim_kinds)
     
     ! This is a structure that explains where FATES patch boundaries
     ! on each thread point to in the host IO array, this structure is
     ! allocated by number of threads. This could be dynamically
     ! allocated, but is unlikely to change...?
     ! Note: history io also instanteates fates_io_dimension_type
     type(fates_io_dimension_type) :: dim_bounds(fates_restart_num_dimensions)
     
     type(restart_map_type), pointer :: restart_map(:)

     integer, private :: cohort_index_, column_index_

   contains
     
     procedure, public :: Init
     procedure, public :: SetThreadBoundsEach
     procedure, public :: assemble_restart_output_types
     procedure, public :: initialize_restart_vars
     procedure, public :: num_restart_vars
     procedure, public :: column_index
     procedure, public :: cohort_index
     procedure, public :: set_restart_vectors
     procedure, public :: create_patchcohort_structure
     procedure, public :: get_restart_vectors
     ! private work functions
     procedure, private :: init_dim_kinds_maps
     procedure, private :: set_dim_indices
     procedure, private :: set_cohort_index
     procedure, private :: set_column_index
     procedure, private :: flush_rvars
     procedure, private :: define_restart_vars
     procedure, private :: set_restart_var

  end type fates_restart_interface_type

  


contains

  ! =====================================================================================
  
  subroutine Init(this, num_threads, fates_bounds)
    
    use FatesIODimensionsMod, only : fates_bounds_type, column, cohort

    implicit none

    class(fates_restart_interface_type), intent(inout) :: this
    integer, intent(in) :: num_threads
    type(fates_bounds_type), intent(in) :: fates_bounds

    integer :: dim_count = 0

    dim_count = dim_count + 1
    call this%set_cohort_index(dim_count)
    call this%dim_bounds(dim_count)%Init(cohort, num_threads, &
         fates_bounds%cohort_begin, fates_bounds%cohort_end)

    dim_count = dim_count + 1
    call this%set_column_index(dim_count)
    call this%dim_bounds(dim_count)%Init(column, num_threads, &
         fates_bounds%column_begin, fates_bounds%column_end)

    ! FIXME(bja, 2016-10) assert(dim_count == FatesIOdimensionsmod::num_dimension_types)

    ! Allocate the mapping between FATES indices and the IO indices
    allocate(this%restart_map(num_threads))
    
  end subroutine Init  

  ! ======================================================================

  subroutine SetThreadBoundsEach(this, thread_index, thread_bounds)
    
    use FatesIODimensionsMod, only : fates_bounds_type

    implicit none

    class(fates_restart_interface_type), intent(inout) :: this

    integer, intent(in) :: thread_index
    type(fates_bounds_type), intent(in) :: thread_bounds

    integer :: index
    
    index = this%cohort_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%cohort_begin, thread_bounds%cohort_end)
    
    index = this%column_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%column_begin, thread_bounds%column_end)
    
  end subroutine SetThreadBoundsEach

  ! ===================================================================================

  subroutine assemble_restart_output_types(this)
    
    use FatesIOVariableKindMod, only : site_r8, site_int, cohort_r8, cohort_int

    implicit none
   
    class(fates_restart_interface_type), intent(inout) :: this

    call this%init_dim_kinds_maps()

    call this%set_dim_indices(cohort_r8, 1, this%cohort_index())
    call this%set_dim_indices(cohort_int, 1, this%cohort_index())

    call this%set_dim_indices(site_r8, 1, this%column_index())
    call this%set_dim_indices(site_int, 1, this%column_index())

  end subroutine assemble_restart_output_types

 ! ===================================================================================
  
  subroutine set_dim_indices(this, dk_name, idim, dim_index)

    use FatesIOVariableKindMod , only : iotype_index

    implicit none

    ! arguments
    class(fates_restart_interface_type), intent(inout) :: this
    character(len=*), intent(in)     :: dk_name
    integer, intent(in)              :: idim  ! dimension index
    integer, intent(in) :: dim_index


    ! local
    integer :: ityp

    ityp = iotype_index(trim(dk_name), fates_restart_num_dim_kinds, this%dim_kinds)

    ! First check to see if the dimension is allocated
    if (this%dim_kinds(ityp)%ndims < idim) then
       write(fates_log(), *) 'Trying to define dimension size to a dim-type structure'
       write(fates_log(), *) 'but the dimension index does not exist'
       write(fates_log(), *) 'type: ',dk_name,' ndims: ',this%dim_kinds(ityp)%ndims,' input dim:',idim
       stop
       !end_run
    end if

    if (idim == 1) then
       this%dim_kinds(ityp)%dim1_index = dim_index
    else if (idim == 2) then
       this%dim_kinds(ityp)%dim2_index = dim_index
    end if

    ! With the map, we can set the dimension size
    this%dim_kinds(ityp)%dimsize(idim) = this%dim_bounds(dim_index)%upper_bound - &
         this%dim_bounds(dim_index)%lower_bound + 1

 end subroutine set_dim_indices


  ! =======================================================================

  subroutine set_cohort_index(this, index)
    implicit none
    class(fates_restart_interface_type), intent(inout) :: this
    integer, intent(in) :: index
    this%cohort_index_ = index
  end subroutine set_cohort_index
  
  integer function cohort_index(this)
    implicit none
    class(fates_restart_interface_type), intent(in) :: this
    cohort_index = this%cohort_index_
  end function cohort_index
  
  ! =======================================================================

  subroutine set_column_index(this, index)
    implicit none
    class(fates_restart_interface_type), intent(inout) :: this
    integer, intent(in) :: index
    this%column_index_ = index
  end subroutine set_column_index
  
  integer function column_index(this)
    implicit none
    class(fates_restart_interface_type), intent(in) :: this
    column_index = this%column_index_
 end function column_index
 
 ! =======================================================================

 subroutine init_dim_kinds_maps(this)
    
    ! ----------------------------------------------------------------------------------
    ! This subroutine simply initializes the structures that define the different
    ! array and type formats for different IO variables
    !
    ! CO_R8   : 1D cohort scale 8-byte reals
    ! SI_R8   : 1D site scale 8-byte reals
    ! CO_INT  : 1D cohort scale integers
    ! SI_INT  : 1D site scale integers
    !
    ! The allocation on the structures is not dynamic and should only add up to the
    ! number of entries listed here.
    !
    ! ----------------------------------------------------------------------------------
    use FatesIOVariableKindMod, only : site_r8, site_int, cohort_r8, cohort_int
    
    implicit none
    
    ! Arguments
    class(fates_restart_interface_type), intent(inout) :: this

    integer :: index

    ! 1d cohort r8
    index = 1
    call this%dim_kinds(index)%Init(cohort_r8, 1)

    ! 1d Site r8
    index = index + 1
    call this%dim_kinds(index)%Init(site_r8, 1)

    ! cohort int
    index = index + 1
    call this%dim_kinds(index)%Init(cohort_int, 1)

    ! site int
    index = index + 1
    call this%dim_kinds(index)%Init(site_int, 1)

    ! FIXME(bja, 2016-10) assert(index == fates_num_dim_kinds)
  end subroutine init_dim_kinds_maps


  ! ====================================================================================

  integer function num_restart_vars(this)
    
    implicit none

    class(fates_restart_interface_type), intent(in) :: this

    num_restart_vars = this%num_restart_vars_
    
  end function num_restart_vars
  
  ! ====================================================================================
  
  subroutine initialize_restart_vars(this)

    implicit none

    class(fates_restart_interface_type), intent(inout) :: this

   ! Determine how many of the restart IO variables registered in FATES
   ! are going to be allocated
   call this%define_restart_vars(initialize_variables=.false.)

   ! Allocate the list of restart output variable objects
   allocate(this%rvars(this%num_restart_vars()))
   
   ! construct the object that defines all of the IO variables
   call this%define_restart_vars(initialize_variables=.true.)
   
 end subroutine initialize_restart_vars

  ! ======================================================================================

 subroutine flush_rvars(this,nc)
 
   class(fates_restart_interface_type)        :: this
   integer,intent(in)                         :: nc

   integer                                   :: ivar
   type(fates_restart_variable_type),pointer :: rvar
   integer                      :: lb1,ub1,lb2,ub2

   do ivar=1,ubound(this%rvars,1)
      associate( rvar => this%rvars(ivar) )
        call rvar%Flush(nc, this%dim_bounds, this%dim_kinds)
      end associate
   end do
   
 end subroutine flush_rvars

 

 ! ====================================================================================
 
 subroutine define_restart_vars(this, initialize_variables)
    
    ! ---------------------------------------------------------------------------------
    ! 
    !                    REGISTRY OF RESTART OUTPUT VARIABLES
    !
    ! Please add any restart variables to this registry. This registry will handle
    ! all variables that can make use of 1D column dimensioned or 1D cohort dimensioned
    ! variables.  Note that restarts are only using 1D vectors in ALM and CLM.  If you
    ! have a multi-dimensional variable that is below the cohort scale, then pack
    ! that variable into a cohort-sized output array by giving it a vtype "cohort_r8"
    ! or "cohort_int".  
    !
    ! Unlike history variables, restarts flush to zero.
    ! ---------------------------------------------------------------------------------
   
    use FatesIOVariableKindMod, only : site_r8, site_int, cohort_int, cohort_r8
    implicit none
    
    class(fates_restart_interface_type), intent(inout) :: this
    logical, intent(in) :: initialize_variables  ! are we 'count'ing or 'initializ'ing?
    integer :: ivar
    real(r8), parameter :: flushinvalid = -9999.0
    real(r8), parameter :: flushzero = 0.0
    real(r8), parameter :: flushone  = 1.0

    
    ivar=0

    ! -----------------------------------------------------------------------------------
    ! Site level variables
    ! -----------------------------------------------------------------------------------

    call this%set_restart_var(vname='fates_PatchesPerSite', vtype=site_int, &
         long_name='Total number of FATES patches per column', units='none', flushval = flushinvalid, &
          hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_npatch_si )

    call this%set_restart_var(vname='fates_old_stock',  vtype=site_r8, &
         long_name='biomass stock in each site (previous step)', units='kgC/site', &
         flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_oldstock_si )

    call this%set_restart_var(vname='fates_cold_dec_status', vtype=site_r8, &
         long_name='status flag for cold deciduous plants', units='unitless', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_cd_status_si )

    call this%set_restart_var(vname='fates_drought_dec_status', vtype=site_r8, &
         long_name='status flag for drought deciduous plants', units='unitless', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_dd_status_si )

    call this%set_restart_var(vname='fates_chilling_days', vtype=site_r8, &
         long_name='chilling day counter', units='unitless', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_nchill_days_si )

    call this%set_restart_var(vname='fates_leafondate', vtype=site_r8, &
         long_name='the day of year for leaf on', units='day of year', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_leafondate_si )

    call this%set_restart_var(vname='fates_leafoffdate', vtype=site_r8, &
         long_name='the day of year for leaf off', units='day of year', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_leafoffdate_si )

    call this%set_restart_var(vname='fates_drought_leafondate', vtype=site_r8, &
         long_name='the day of year for drought based leaf-on', units='day of year', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_dleafondate_si )

    call this%set_restart_var(vname='fates_drought_leafoffdate', vtype=site_r8, &
         long_name='the day of year for drought based leaf-off', units='day of year', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_dleafoffdate_si )

    call this%set_restart_var(vname='fates_acc_nesterov_id', vtype=site_r8, &
         long_name='a nesterov index accumulator', units='unitless', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_acc_ni_si )
    
    call this%set_restart_var(vname='fates_gdd_site', vtype=site_r8, &
         long_name='growing degree days at each site', units='degC days', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_gdd_si )

    call this%set_restart_var(vname='fates_nep_timeintegrated_site', vtype=site_r8, &
         long_name='NEP integrated over model time-steps', units='gc/m2', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_nep_timeintegrated_si )

    call this%set_restart_var(vname='fates_npp_timeintegrated_site', vtype=site_r8, &
         long_name='NPP integrated over model time-steps', units='gc/m2', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_npp_timeintegrated_si )

    call this%set_restart_var(vname='fates_hr_timeintegrated_site', vtype=site_r8, &
         long_name='heterotrophic respiration integrated over model time-steps', &
         units='gc/m2', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_hr_timeintegrated_si )

    call this%set_restart_var(vname='fates_cbal_err_fatesite', vtype=site_r8, &
         long_name='the carbon accounting error for FATES processes', &
         units='gC/m2/s', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_cbal_error_fates_si )

    call this%set_restart_var(vname='fates_cbal_err_bgcsite', vtype=site_r8, &
         long_name='the carbon accounting error for (fates relevant) BGC processes', &
         units='gC/m2/s', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_cbal_error_bgc_si )

    call this%set_restart_var(vname='fates_cbal_err_totsite', vtype=site_r8, &
         long_name='the carbon accounting error for fates and bgc processes', &
         units='gC/m2/s', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_cbal_error_total_si )

    call this%set_restart_var(vname='fates_totecosysc_old_site', vtype=site_r8, &
         long_name='total ecosystem carbon above and below ground (previous time-step)', &
         units='gC/m2', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_totecosysc_old_si )

    call this%set_restart_var(vname='fates_totfatesc_old_site', vtype=site_r8, &
         long_name='total carbon tracked in FATES, (previous time-step)', &
         units='gc/m2', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_totfatesc_old_si )

    call this%set_restart_var(vname='fates_totbgcc_old_site', vtype=site_r8, &
         long_name='total carbon tracked in the BGC module', &
         units='gc/m2', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_totbgcc_old_si )

    call this%set_restart_var(vname='fates_to_bgc_this_edts_col', vtype=site_r8, &
         long_name='total flux of carbon from FATES to BGC models on current timestep', &
         units='gC/m2/s', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_fates_to_bgc_this_ts_si )

    call this%set_restart_var(vname='fates_to_bgc_last_edts_col', vtype=site_r8, &
         long_name='total flux of carbon from FATES to BGC models on previous timestep', &
         units='gC/m2/s', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_fates_to_bgc_last_ts_si )

    call this%set_restart_var(vname='fates_seed_rain_flux_site', vtype=site_r8, &
         long_name='flux of seeds from exterior', &
         units='kgC/m2/year', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_seedrainflux_si )

    call this%set_restart_var(vname='fates_trunk_product_site', vtype=site_r8, &
         long_name='Accumulate trunk product flux at site', &
         units='kgC/m2', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_trunk_product_si )


    ! -----------------------------------------------------------------------------------
    ! Variables stored within cohort vectors
    ! Note: Some of these are multi-dimensional variables in the patch/site dimension
    ! that are collapsed into the cohort vectors for storage and transfer
    ! -----------------------------------------------------------------------------------

    ! This variable may be confusing, because it is a patch level variables
    ! but it is using the cohort IO vector to hold data
    call this%set_restart_var(vname='fates_CohortsPerPatch', vtype=cohort_int, &
         long_name='the number of cohorts per patch', units='unitless', flushval = flushinvalid, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_ncohort_pa )

    ! 1D cohort Variables
    ! -----------------------------------------------------------------------------------

    call this%set_restart_var(vname='fates_balive', vtype=cohort_r8, &
         long_name='ed cohort alive biomass', units='kgC/indiv', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_balive_co )

    call this%set_restart_var(vname='fates_bdead', vtype=cohort_r8, &
         long_name='ed cohort - dead (structural) biomass in living plants', &
         units='kgC/indiv', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_bdead_co )

    call this%set_restart_var(vname='fates_bl', vtype=cohort_r8, &
         long_name='ed cohort - leaf biomass', units='kgC/indiv', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_bleaf_co )

    call this%set_restart_var(vname='fates_br', vtype=cohort_r8, &
         long_name='ed cohort - fine root biomass', units='kgC/indiv', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_broot_co )

    call this%set_restart_var(vname='fates_bstore', vtype=cohort_r8, &
         long_name='ed cohort - storage biomass', units='kgC/indiv', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_bstore_co )

    call this%set_restart_var(vname='fates_canopy_layer', vtype=cohort_r8, &
         long_name='ed cohort - canopy_layer', units='unitless', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_canopy_layer_co )

    call this%set_restart_var(vname='fates_canopy_layer_yesterday', vtype=cohort_r8, &
         long_name='ed cohort - canopy_layer_yesterday', units='unitless', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_canopy_layer_yesterday_co )

    call this%set_restart_var(vname='fates_canopy_trim', vtype=cohort_r8, &
         long_name='ed cohort - canopy_trim', units='fraction', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_canopy_trim_co )

    call this%set_restart_var(vname='fates_dbh', vtype=cohort_r8, &
         long_name='ed cohort - diameter at breast height', units='cm', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_dbh_co )

    call this%set_restart_var(vname='fates_height', vtype=cohort_r8, &
         long_name='ed cohort - plant height', units='m', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_height_co )

    call this%set_restart_var(vname='fates_laimemory', vtype=cohort_r8, &
         long_name='ed cohort - target leaf biomass set from prev year', &
         units='kgC/indiv', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_laimemory_co )

    call this%set_restart_var(vname='fates_leaf_maint_dmnd', vtype=cohort_r8, &
         long_name='ed cohort - leaf maintenance demand', &
         units='kgC/indiv/year', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_leaf_md_co )

    call this%set_restart_var(vname='fates_root_maint_dmnd', vtype=cohort_r8, &
         long_name='ed cohort - fine root maintenance demand', &
         units='kgC/indiv', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_root_md_co )

    call this%set_restart_var(vname='fates_nplant', vtype=cohort_r8, &
         long_name='ed cohort - number of plants in the cohort', &
         units='/patch', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_nplant_co )

    call this%set_restart_var(vname='fates_gpp_acc', vtype=cohort_r8, &
         long_name='ed cohort - accumulated gpp over dynamics step', &
         units='kgC/indiv', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_gpp_acc_co )

    call this%set_restart_var(vname='fates_npp_acc', vtype=cohort_r8, &
         long_name='ed cohort - accumulated npp over dynamics step', &
         units='kgC/indiv', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_npp_acc_co )

    call this%set_restart_var(vname='fates_gpp_acc_hold', vtype=cohort_r8, &
         long_name='ed cohort - current step gpp', &
         units='kgC/indiv/year', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_gpp_acc_hold_co )

    call this%set_restart_var(vname='fates_npp_acc_hold', vtype=cohort_r8, &
         long_name='ed cohort - current step npp', &
         units='kgC/indiv/year', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_npp_acc_hold_co )

    call this%set_restart_var(vname='fates_npp_leaf', vtype=cohort_r8, &
         long_name='ed cohort - npp sent to leaves', &
         units='kgC/indiv/day', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_npp_leaf_co )

    call this%set_restart_var(vname='fates_npp_froot', vtype=cohort_r8, &
         long_name='ed cohort - npp sent to fine roots', &
         units='kgC/indiv/day', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_npp_froot_co )

    call this%set_restart_var(vname='fates_npp_sapwood', vtype=cohort_r8, &
         long_name='ed cohort - npp sent to sapwood', &
         units='kgC/indiv/day', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_npp_sw_co )

    call this%set_restart_var(vname='fates_npp_bdead', vtype=cohort_r8, &
         long_name='ed cohort - npp sent to dead (structure) biomass in live plants', &
         units='kgC/indiv/day', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_npp_dead_co )

    call this%set_restart_var(vname='fates_npp_seed', vtype=cohort_r8, &
         long_name='ed cohort - npp sent to seed biomass', &
         units='kgC/indiv/day', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_npp_seed_co )

    call this%set_restart_var(vname='fates_npp_store', vtype=cohort_r8, &
         long_name='ed cohort - npp sent to storage biomass', &
         units='kgC/indiv/day', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_npp_store_co )

    call this%set_restart_var(vname='fates_bmort', vtype=cohort_r8, &
         long_name='ed cohort - background mortality rate', &
         units='/year', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_bmort_co )

    call this%set_restart_var(vname='fates_hmort', vtype=cohort_r8, &
         long_name='ed cohort - hydraulic mortality rate', &
         units='/year', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_hmort_co )

    call this%set_restart_var(vname='fates_cmort', vtype=cohort_r8, &
         long_name='ed cohort - carbon starvation mortality rate', &
         units='/year', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_cmort_co )

    call this%set_restart_var(vname='fates_imort', vtype=cohort_r8, &
         long_name='ed cohort - impact mortality rate', &
         units='/year', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_imort_co )

    call this%set_restart_var(vname='fates_fmort', vtype=cohort_r8, &
         long_name='ed cohort - frost mortality rate', &
         units='/year', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_fmort_co )


    call this%set_restart_var(vname='fates_lmort_logging', vtype=cohort_r8, &
         long_name='ed cohort - directly logging mortality rate', &
         units='%/event', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_lmort_logging_co )

    call this%set_restart_var(vname='fates_lmort_collateral', vtype=cohort_r8, &
         long_name='ed cohort - collateral mortality rate', &
         units='%/event', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_lmort_collateral_co ) 
  
    call this%set_restart_var(vname='fates_lmort_in', vtype=cohort_r8, &
         long_name='ed cohort - mechanical mortality rate', &
         units='%/event', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_lmort_infra_co ) 


    call this%set_restart_var(vname='fates_ddbhdt', vtype=cohort_r8, &
         long_name='ed cohort - differential: ddbh/dt', &
         units='cm/year', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_ddbhdt_co )

    call this%set_restart_var(vname='fates_dbalivedt', vtype=cohort_r8, &
         long_name='ed cohort - differential: ddbh/dt', &
         units='cm/year', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_dbalivedt_co )

    call this%set_restart_var(vname='fates_dbdeaddt', vtype=cohort_r8, &
         long_name='ed cohort - differential: ddbh/dt', &
         units='cm/year', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_dbdeaddt_co )

    call this%set_restart_var(vname='fates_dbstoredt', vtype=cohort_r8, &
         long_name='ed cohort - differential: ddbh/dt', &
         units='cm/year', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_dbstoredt_co )

    call this%set_restart_var(vname='fates_resp_tstep', vtype=cohort_r8, &
         long_name='ed cohort - autotrophic respiration over timestep', &
         units='kgC/indiv/timestep', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_resp_tstep_co )

    call this%set_restart_var(vname='fates_pft', vtype=cohort_int, &
         long_name='ed cohort - plant functional type', units='index', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_pft_co )

    call this%set_restart_var(vname='fates_status_coh', vtype=cohort_int, &
         long_name='ed cohort - plant phenology status', units='unitless', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_status_co )

    call this%set_restart_var(vname='fates_isnew', vtype=cohort_int, &
         long_name='ed cohort - binary flag specifying if a plant has experienced a full day cycle', &
         units='0/1', flushval = flushone, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_isnew_co )


    ! Mixed dimension variables using the cohort vector
    ! -----------------------------------------------------------------------------------

    call this%set_restart_var(vname='fates_cwd_ag', vtype=cohort_r8, &
         long_name='coarse woody debris above ground (non-respiring), by patch x cw class', &
         units='kgC/m2', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_cwd_ag_pacw )

    call this%set_restart_var(vname='fates_cwd_bg', vtype=cohort_r8, &
         long_name='coarse woody debris below ground (non-respiring), by patch x cw class', &
         units='kgC/m2', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_cwd_bg_pacw )

    call this%set_restart_var(vname='fates_leaf_litter', vtype=cohort_r8, &
         long_name='leaf litter, by patch x pft (non-respiring)', &
         units='kgC/m2', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_leaf_litter_paft )

    call this%set_restart_var(vname='fates_root_litter', vtype=cohort_r8, &
         long_name='root litter, by patch x pft (non-respiring)', &
         units='kgC/m2', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_root_litter_paft )

    call this%set_restart_var(vname='fates_leaf_litter_in', vtype=cohort_r8, &
         long_name='leaf litter flux from turnover and mort, by patch x pft', &
         units='kgC/m2', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_leaf_litter_in_paft )

    call this%set_restart_var(vname='fates_root_litter_in', vtype=cohort_r8, &
         long_name='root litter flux from turnover and mort, by patch x pft', &
         units='kgC/m2', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_root_litter_in_paft )

    call this%set_restart_var(vname='fates_seed_bank', vtype=cohort_r8, &
         long_name='seed pool for each functional type, by site x pft', &
         units='kgC/m2/year', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_seed_bank_sift )

    call this%set_restart_var(vname='fates_spread', vtype=site_r8, &
         long_name='dynamic ratio of dbh to canopy area, by patch x canopy-layer', &
         units='cm/m2', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_spread_si )

    call this%set_restart_var(vname='fates_livegrass', vtype=cohort_r8, &
         long_name='total AGB from grass, by patch', &
         units='kgC/m2', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_livegrass_pa )

    call this%set_restart_var(vname='fates_age', vtype=cohort_r8, &
         long_name='age of the ED patch', units='yr', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_age_pa )

    call this%set_restart_var(vname='fates_area', vtype=cohort_r8, &
         long_name='are of the ED patch', units='m2', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_area_pa )

    ! These dimensions are pa "patch" cl "canopy layer" ft "functional type" ls "layer sublevel"
    call this%set_restart_var(vname='fates_f_sun', vtype=cohort_r8, &
         long_name='fraction of sunlit leaves, by patch x can-layer x pft x sublayer', &
         units='fraction', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_fsun_paclftls )

    call this%set_restart_var(vname='fates_fabd_sun_z', vtype=cohort_r8, &
         long_name='sun fraction of direct light absorbed, by patch x can-layer x pft x sublayer', &
         units='fraction', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_fabd_sun_paclftls )

    call this%set_restart_var(vname='fates_fabi_sun_z', vtype=cohort_r8, &
         long_name='sun fraction of indirect light absorbed, by patch x can-layer x pft x sublayer', &
         units='fraction', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_fabi_sun_paclftls )

    call this%set_restart_var(vname='fates_fabd_sha_z', vtype=cohort_r8, &
         long_name='shade fraction of direct light absorbed, by patch x can-layer x pft x sublayer', &
         units='fraction', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_fabd_sha_paclftls )

    call this%set_restart_var(vname='fates_fabi_sha_z', vtype=cohort_r8, &
         long_name='shade fraction of indirect light absorbed, by patch x can-layer x pft x sublayer', &
         units='fraction', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_fabi_sha_paclftls )

    !
    ! site x time level vars
    !

    call this%set_restart_var(vname='fates_water_memory', vtype=cohort_r8, &
         long_name='last 10 days of volumetric soil water, by site x day-index', &
         units='m3/m3', flushval = flushzero, &
         hlms='CLM:ALM', initialize=initialize_variables, ivar=ivar, index = ir_watermem_siwm )
         
    
    ! Must be last thing before return
    this%num_restart_vars_ = ivar
    
  end subroutine define_restart_vars
    

  ! =====================================================================================
   
  subroutine set_restart_var(this,vname,vtype,long_name,units,flushval, &
        hlms,initialize,ivar,index)

    use FatesUtilsMod, only : check_hlm_list
    use FatesInterfaceMod, only : hlm_name

    ! arguments
    class(fates_restart_interface_type) :: this
    character(len=*),intent(in)  :: vname
    character(len=*),intent(in)  :: vtype
    character(len=*),intent(in)  :: units 
    real(r8), intent(in)         :: flushval
    character(len=*),intent(in)  :: long_name
    character(len=*),intent(in)  :: hlms
    logical, intent(in)          :: initialize
    integer, intent(inout)       :: ivar
    integer, intent(inout)       :: index  ! This is the index for the variable of
                                           ! interest that is associated with an
                                           ! explict name (for fast reference during update)
                                           ! A zero is passed back when the variable is
                                           ! not used

    
    type(fates_restart_variable_type),pointer :: rvar
    integer :: ub1,lb1,ub2,lb2    ! Bounds for allocating the var
    integer :: ityp
    
    logical :: use_var
    
    use_var = check_hlm_list(trim(hlms), trim(hlm_name))


    if( use_var ) then
       
       ivar  = ivar+1
       index = ivar    
       
       if( initialize )then
          
          call this%rvars(ivar)%Init(vname, units, long_name, vtype, flushval, &
               fates_restart_num_dim_kinds, this%dim_kinds, this%dim_bounds)

       end if
    else
       
       index = 0
    end if
    
    return
 end subroutine set_restart_var

 ! =====================================================================================

 subroutine set_restart_vectors(this,nc,nsites,sites)

   use EDTypesMod, only : nclmax
   use EDTypesMod, only : nlevleaf
   use FatesInterfaceMod, only : fates_maxElementsPerPatch
   use FatesInterfaceMod, only : numpft
   use EDTypesMod, only : ed_site_type
   use EDTypesMod, only : ed_cohort_type
   use EDTypesMod, only : ed_patch_type
   use EDTypesMod, only : ncwd
   use EDTypesMod, only : numWaterMem

    ! Arguments
    class(fates_restart_interface_type)             :: this
    integer                 , intent(in)            :: nc   ! clump index
    integer                 , intent(in)            :: nsites
    type(ed_site_type)      , intent(inout), target :: sites(nsites)

    ! Locals
    integer  :: s        ! The local site index

    ! ----------------------------------------------------------------------------------
    ! The following group of integers indicate the positional index (idx)
    ! of variables at different scales inside the I/O arrays (io)
    ! Keep in mind that many of these variables have a composite dimension
    ! at the patch scale.  To hold this memory, we borrow the cohort
    ! vector.  Thus the head of each array points to the first cohort
    ! of each patch. "io_idx_co_1st"
    ! ----------------------------------------------------------------------------------
    integer  :: io_idx_si      ! site index
    integer  :: io_idx_co_1st  ! 1st cohort of each patch
    integer  :: io_idx_co      ! cohort index
    integer  :: io_idx_pa_pft  ! each pft within each patch (pa_pft)
    integer  :: io_idx_pa_cwd  ! each cwd class within each patch (pa_cwd)
    integer  :: io_idx_pa_sunz ! index for the combined dimensions for radiation
    integer  :: io_idx_si_wmem ! each water memory class within each site
    
    ! Some counters (for checking mostly)
    integer  :: totalcohorts   ! total cohort count on this thread (diagnostic)
    integer  :: patchespersite   ! number of patches per site
    integer  :: cohortsperpatch  ! number of cohorts per patch 

    integer  :: ft               ! functional type index
    integer  :: k,j,i            ! indices to the radiation matrix
    
    type(fates_restart_variable_type) :: rvar
    type(ed_patch_type),pointer  :: cpatch
    type(ed_cohort_type),pointer :: ccohort


    associate( rio_npatch_si           => this%rvars(ir_npatch_si)%int1d, &
           rio_old_stock_si             => this%rvars(ir_oldstock_si)%r81d, &
           rio_cd_status_si            => this%rvars(ir_cd_status_si)%r81d, &
           rio_dd_status_si            => this%rvars(ir_dd_status_si)%r81d, &
           rio_nchill_days_si          => this%rvars(ir_nchill_days_si)%r81d, &
           rio_leafondate_si           => this%rvars(ir_leafondate_si)%r81d, &
           rio_leafoffdate_si          => this%rvars(ir_leafoffdate_si)%r81d, &
           rio_dleafondate_si          => this%rvars(ir_dleafondate_si)%r81d, &
           rio_dleafoffdate_si         => this%rvars(ir_dleafoffdate_si)%r81d, &
           rio_acc_ni_si               => this%rvars(ir_acc_ni_si)%r81d, &
           rio_gdd_si                  => this%rvars(ir_gdd_si)%r81d, &
           rio_nep_timeintegrated_si   => this%rvars(ir_nep_timeintegrated_si)%r81d, &
           rio_npp_timeintegrated_si   => this%rvars(ir_npp_timeintegrated_si)%r81d, &
           rio_hr_timeintegrated_si    => this%rvars(ir_hr_timeintegrated_si)%r81d, &
           rio_cbal_err_fates_si       => this%rvars(ir_cbal_error_fates_si)%r81d, &
           rio_cbal_err_bgc_si         => this%rvars(ir_cbal_error_bgc_si)%r81d, &
           rio_cbal_err_tot_si         => this%rvars(ir_cbal_error_total_si)%r81d, &
           rio_totecosysc_old_si       => this%rvars(ir_totecosysc_old_si)%r81d, &
           rio_totfatesc_old_si        => this%rvars(ir_totfatesc_old_si)%r81d, &
           rio_totbgcc_old_si          => this%rvars(ir_totbgcc_old_si)%r81d, &
           rio_fates_to_bgc_this_ts_si => this%rvars(ir_fates_to_bgc_this_ts_si)%r81d, &
           rio_fates_to_bgc_last_ts_si => this%rvars(ir_fates_to_bgc_last_ts_si)%r81d, &
           rio_seedrainflux_si         => this%rvars(ir_seedrainflux_si)%r81d, &
           rio_trunk_product_si        => this%rvars(ir_trunk_product_si)%r81d, &
           rio_ncohort_pa              => this%rvars(ir_ncohort_pa)%int1d, &
           rio_balive_co               => this%rvars(ir_balive_co)%r81d, &
           rio_bdead_co                => this%rvars(ir_bdead_co)%r81d, &
           rio_bleaf_co                => this%rvars(ir_bleaf_co)%r81d, &
           rio_broot_co                => this%rvars(ir_broot_co)%r81d, &
           rio_bstore_co               => this%rvars(ir_bstore_co)%r81d, &
           rio_canopy_layer_co         => this%rvars(ir_canopy_layer_co)%r81d, &
           rio_canopy_layer_yesterday_co    => this%rvars(ir_canopy_layer_yesterday_co)%r81d, &
           rio_canopy_trim_co          => this%rvars(ir_canopy_trim_co)%r81d, &
           rio_dbh_co                  => this%rvars(ir_dbh_co)%r81d, &
           rio_height_co               => this%rvars(ir_height_co)%r81d, &
           rio_laimemory_co            => this%rvars(ir_laimemory_co)%r81d, &
           rio_leaf_md_co              => this%rvars(ir_leaf_md_co)%r81d, &
           rio_root_md_co              => this%rvars(ir_root_md_co)%r81d, &
           rio_nplant_co               => this%rvars(ir_nplant_co)%r81d, &
           rio_gpp_acc_co              => this%rvars(ir_gpp_acc_co)%r81d, &
           rio_npp_acc_co              => this%rvars(ir_npp_acc_co)%r81d, &
           rio_gpp_acc_hold_co         => this%rvars(ir_gpp_acc_hold_co)%r81d, &
           rio_npp_acc_hold_co         => this%rvars(ir_npp_acc_hold_co)%r81d, &
           rio_npp_leaf_co             => this%rvars(ir_npp_leaf_co)%r81d, &
           rio_npp_froot_co            => this%rvars(ir_npp_froot_co)%r81d, &
           rio_npp_sw_co               => this%rvars(ir_npp_sw_co)%r81d, &
           rio_npp_dead_co             => this%rvars(ir_npp_dead_co)%r81d, &
           rio_npp_seed_co             => this%rvars(ir_npp_seed_co)%r81d, &
           rio_npp_store_co            => this%rvars(ir_npp_store_co)%r81d, &
           rio_bmort_co                => this%rvars(ir_bmort_co)%r81d, &
           rio_hmort_co                => this%rvars(ir_hmort_co)%r81d, &
           rio_cmort_co                => this%rvars(ir_cmort_co)%r81d, &
           rio_imort_co                => this%rvars(ir_imort_co)%r81d, &
           rio_fmort_co                => this%rvars(ir_fmort_co)%r81d, &



	   rio_lmort_logging_co                => this%rvars(ir_lmort_logging_co)%r81d, &
	   rio_lmort_collateral_co              => this%rvars(ir_lmort_collateral_co)%r81d, &
	   rio_lmort_infra_co                => this%rvars(ir_lmort_infra_co)%r81d, &


           rio_ddbhdt_co               => this%rvars(ir_ddbhdt_co)%r81d, &
           rio_dbalivedt_co            => this%rvars(ir_dbalivedt_co)%r81d, &
           rio_dbdeaddt_co             => this%rvars(ir_dbdeaddt_co)%r81d, &
           rio_dbstoredt_co            => this%rvars(ir_dbstoredt_co)%r81d, &
           rio_resp_tstep_co           => this%rvars(ir_resp_tstep_co)%r81d, &
           rio_pft_co                  => this%rvars(ir_pft_co)%int1d, &
           rio_status_co               => this%rvars(ir_status_co)%int1d, &
           rio_isnew_co                => this%rvars(ir_isnew_co)%int1d, &
           rio_cwd_ag_pacw             => this%rvars(ir_cwd_ag_pacw)%r81d, &
           rio_cwd_bg_pacw             => this%rvars(ir_cwd_bg_pacw)%r81d, &
           rio_leaf_litter_paft        => this%rvars(ir_leaf_litter_paft)%r81d, &
           rio_root_litter_paft        => this%rvars(ir_root_litter_paft)%r81d, &
           rio_leaf_litter_in_paft     => this%rvars(ir_leaf_litter_in_paft)%r81d, &
           rio_root_litter_in_paft     => this%rvars(ir_root_litter_in_paft)%r81d, &
           rio_seed_bank_sift          => this%rvars(ir_seed_bank_sift)%r81d, &
           rio_spread_si               => this%rvars(ir_spread_si)%r81d, &
           rio_livegrass_pa            => this%rvars(ir_livegrass_pa)%r81d, &
           rio_age_pa                  => this%rvars(ir_age_pa)%r81d, &
           rio_area_pa                 => this%rvars(ir_area_pa)%r81d, &
           rio_fsun_paclftls           => this%rvars(ir_fsun_paclftls)%r81d, &
           rio_fabd_sun_z_paclftls     => this%rvars(ir_fabd_sun_paclftls)%r81d, &
           rio_fabi_sun_z_paclftls     => this%rvars(ir_fabi_sun_paclftls)%r81d, &
           rio_fabd_sha_z_paclftls     => this%rvars(ir_fabd_sha_paclftls)%r81d, &
           rio_fabi_sha_z_paclftls     => this%rvars(ir_fabi_sha_paclftls)%r81d, &
           rio_watermem_siwm           => this%rvars(ir_watermem_siwm)%r81d )

       totalCohorts = 0
       
       ! ---------------------------------------------------------------------------------
       ! Flush arrays to values defined by %flushval (see registry entry in
       ! subroutine define_history_vars()
       ! ---------------------------------------------------------------------------------
       call this%flush_rvars(nc)
       
       do s = 1,nsites
          
          ! Calculate the offsets
          ! fcolumn is the global column index of the current site.
          ! For the first site, if that site aligns with the first column index
          ! in the clump, than the offset should be be equal to begCohort
          
          io_idx_si      = this%restart_map(nc)%site_index(s)
          io_idx_co_1st  = this%restart_map(nc)%cohort1_index(s)

          io_idx_co      = io_idx_co_1st
          io_idx_pa_pft  = io_idx_co_1st
          io_idx_pa_cwd  = io_idx_co_1st
          io_idx_si_wmem = io_idx_co_1st
          io_idx_pa_sunz = io_idx_co_1st
          
          ! write seed_bank info(site-level, but PFT-resolved)
          do i = 1,numpft
             rio_seed_bank_sift(io_idx_co_1st+i-1) = sites(s)%seed_bank(i)
          end do

          ! canopy spread term
          rio_spread_si(io_idx_si)   = sites(s)%spread
          
          cpatch => sites(s)%oldest_patch
          
          ! new column, reset num patches
          patchespersite = 0
          
          do while(associated(cpatch))
             
             ! found patch, increment
             patchespersite = patchespersite + 1
             
             ccohort => cpatch%shortest
             
             ! new patch, reset num cohorts
             cohortsperpatch = 0
             
             do while(associated(ccohort))
                
                ! found cohort, increment
                cohortsperpatch = cohortsperpatch + 1
                totalCohorts    = totalCohorts + 1
             
                if ( DEBUG ) then
                   write(fates_log(),*) 'CLTV io_idx_co ', io_idx_co
                   write(fates_log(),*) 'CLTV lowerbound ', lbound(rio_npp_acc_co,1) 
                   write(fates_log(),*) 'CLTV upperbound  ', ubound(rio_npp_acc_co,1)
                endif
             
                rio_balive_co(io_idx_co)       = ccohort%balive
                rio_bdead_co(io_idx_co)        = ccohort%bdead
                rio_bleaf_co(io_idx_co)        = ccohort%bl
                rio_broot_co(io_idx_co)        = ccohort%br
                rio_bstore_co(io_idx_co)       = ccohort%bstore
                rio_canopy_layer_co(io_idx_co) = ccohort%canopy_layer
                rio_canopy_layer_yesterday_co(io_idx_co) = ccohort%canopy_layer_yesterday
                rio_canopy_trim_co(io_idx_co)  = ccohort%canopy_trim
                rio_dbh_co(io_idx_co)          = ccohort%dbh
                rio_height_co(io_idx_co)       = ccohort%hite
                rio_laimemory_co(io_idx_co)    = ccohort%laimemory
                rio_leaf_md_co(io_idx_co)      = ccohort%leaf_md
                rio_root_md_co(io_idx_co)      = ccohort%root_md
                rio_nplant_co(io_idx_co)       = ccohort%n
                rio_gpp_acc_co(io_idx_co)      = ccohort%gpp_acc
                rio_npp_acc_co(io_idx_co)      = ccohort%npp_acc
                rio_gpp_acc_hold_co(io_idx_co) = ccohort%gpp_acc_hold
                rio_npp_acc_hold_co(io_idx_co) = ccohort%npp_acc_hold
                rio_npp_leaf_co(io_idx_co)     = ccohort%npp_leaf
                rio_npp_froot_co(io_idx_co)    = ccohort%npp_froot
                rio_npp_sw_co(io_idx_co)       = ccohort%npp_bsw
                rio_npp_dead_co(io_idx_co)     = ccohort%npp_bdead
                rio_npp_seed_co(io_idx_co)     = ccohort%npp_bseed
                rio_npp_store_co(io_idx_co)    = ccohort%npp_store
                rio_bmort_co(io_idx_co)        = ccohort%bmort
                rio_hmort_co(io_idx_co)        = ccohort%hmort
                rio_cmort_co(io_idx_co)        = ccohort%cmort
                rio_imort_co(io_idx_co)        = ccohort%imort
                rio_fmort_co(io_idx_co)        = ccohort%fmort

                !Logging
	        rio_lmort_logging_co(io_idx_co)        = ccohort%lmort_logging
	        rio_lmort_collateral_co(io_idx_co)     = ccohort%lmort_collateral
	        rio_lmort_infra_co(io_idx_co)        = ccohort%lmort_infra



                rio_ddbhdt_co(io_idx_co)       = ccohort%ddbhdt
                rio_dbalivedt_co(io_idx_co)    = ccohort%dbalivedt
                rio_dbdeaddt_co(io_idx_co)     = ccohort%dbdeaddt
                rio_dbstoredt_co(io_idx_co)    = ccohort%dbstoredt
                rio_resp_tstep_co(io_idx_co)   = ccohort%resp_tstep
                rio_pft_co(io_idx_co)          = ccohort%pft
                rio_status_co(io_idx_co)       = ccohort%status_coh
                if ( ccohort%isnew ) then
                   rio_isnew_co(io_idx_co)     = new_cohort
                else
                   rio_isnew_co(io_idx_co)     = old_cohort
                endif
                
                if ( DEBUG ) then
                   write(fates_log(),*) 'CLTV offsetNumCohorts II ',io_idx_co, &
                         cohortsperpatch
                endif
             
                io_idx_co = io_idx_co + 1
                
                ccohort => ccohort%taller
                
             enddo ! ccohort do while
             
             !
             ! deal with patch level fields here
             !
             rio_livegrass_pa(io_idx_co_1st)   = cpatch%livegrass
             rio_age_pa(io_idx_co_1st)         = cpatch%age
             rio_area_pa(io_idx_co_1st)        = cpatch%area
             
             ! set cohorts per patch for IO
             rio_ncohort_pa( io_idx_co_1st )   = cohortsperpatch
             
             if ( DEBUG ) then
                write(fates_log(),*) 'offsetNumCohorts III ' &
                      ,io_idx_co,cohortsperpatch
             endif
             !
             ! deal with patch level fields of arrays here
             !
             ! these are arrays of length numpft, each patch contains one
             ! vector so we increment 
             do i = 1,numpft
                rio_leaf_litter_paft(io_idx_pa_pft)    = cpatch%leaf_litter(i)
                rio_root_litter_paft(io_idx_pa_pft)    = cpatch%root_litter(i)
                rio_leaf_litter_in_paft(io_idx_pa_pft) = cpatch%leaf_litter_in(i)
                rio_root_litter_in_paft(io_idx_pa_pft) = cpatch%root_litter_in(i)
                io_idx_pa_pft = io_idx_pa_pft + 1
             end do
             
             do i = 1,ncwd ! ncwd currently 4
                rio_cwd_ag_pacw(io_idx_pa_cwd) = cpatch%cwd_ag(i)
                rio_cwd_bg_pacw(io_idx_pa_cwd) = cpatch%cwd_bg(i)
                io_idx_pa_cwd = io_idx_pa_cwd + 1
             end do
             
             if ( DEBUG ) write(fates_log(),*) 'CLTV io_idx_pa_sunz 1 ',io_idx_pa_sunz
             
             if ( DEBUG ) write(fates_log(),*) 'CLTV 1186 ',nlevleaf,numpft,nclmax
             
             do k = 1,nlevleaf     ! nlevleaf currently 40
                do j = 1,numpft    ! dependent on parameter file
                   do i = 1,nclmax ! nclmax currently 2
                      rio_fsun_paclftls(io_idx_pa_sunz)        = cpatch%f_sun(i,j,k)
                      rio_fabd_sun_z_paclftls(io_idx_pa_sunz)  = cpatch%fabd_sun_z(i,j,k)
                      rio_fabi_sun_z_paclftls(io_idx_pa_sunz)  = cpatch%fabi_sun_z(i,j,k)
                      rio_fabd_sha_z_paclftls(io_idx_pa_sunz)  = cpatch%fabd_sha_z(i,j,k)
                      rio_fabi_sha_z_paclftls(io_idx_pa_sunz)  = cpatch%fabi_sha_z(i,j,k)
                      io_idx_pa_sunz = io_idx_pa_sunz + 1
                   end do
                end do
             end do
             
             if ( DEBUG ) write(fates_log(),*) 'CLTV io_idx_pa_sunz 2 ',io_idx_pa_sunz


             ! Set the first cohort index to the start of the next patch, increment
             ! by the maximum number of cohorts per patch
             io_idx_co_1st = io_idx_co_1st + fates_maxElementsPerPatch
             
             ! reset counters so that they are all advanced evenly.
             io_idx_pa_pft  = io_idx_co_1st
             io_idx_pa_cwd  = io_idx_co_1st
             io_idx_co      = io_idx_co_1st
             io_idx_pa_sunz = io_idx_co_1st
             
             if ( DEBUG ) then
                write(fates_log(),*) 'CLTV io_idx_co_1st ', io_idx_co_1st
                write(fates_log(),*) 'CLTV numCohort ', cohortsperpatch
                write(fates_log(),*) 'CLTV totalCohorts ', totalCohorts
             end if
             
             cpatch => cpatch%younger
             
          enddo ! cpatch do while
          
          rio_old_stock_si(io_idx_si)    = sites(s)%old_stock
          rio_cd_status_si(io_idx_si)    = sites(s)%status
          rio_dd_status_si(io_idx_si)    = sites(s)%dstatus
          rio_nchill_days_si(io_idx_si)  = sites(s)%ncd 
          rio_leafondate_si(io_idx_si)   = sites(s)%leafondate
          rio_leafoffdate_si(io_idx_si)  = sites(s)%leafoffdate
          rio_dleafondate_si(io_idx_si)  = sites(s)%dleafondate
          rio_dleafoffdate_si(io_idx_si) = sites(s)%dleafoffdate
          rio_acc_ni_si(io_idx_si)       = sites(s)%acc_NI
          rio_gdd_si(io_idx_si)          = sites(s)%ED_GDD_site
          
          ! Carbon Balance and Checks
          rio_nep_timeintegrated_si(io_idx_si) = sites(s)%nep_timeintegrated 
          rio_npp_timeintegrated_si(io_idx_si) = sites(s)%npp_timeintegrated
          rio_hr_timeintegrated_si(io_idx_si)  = sites(s)%hr_timeintegrated 
          rio_totecosysc_old_si(io_idx_si)     = sites(s)%totecosysc_old
          rio_totfatesc_old_si(io_idx_si)      = sites(s)%totfatesc_old
          rio_totbgcc_old_si(io_idx_si)        = sites(s)%totbgcc_old
          rio_cbal_err_fates_si(io_idx_si)     = sites(s)%cbal_err_fates
          rio_cbal_err_bgc_si(io_idx_si)       = sites(s)%cbal_err_bgc
          rio_cbal_err_tot_si(io_idx_si)       = sites(s)%cbal_err_tot
          rio_fates_to_bgc_this_ts_si(io_idx_si) = sites(s)%fates_to_bgc_this_ts
          rio_fates_to_bgc_last_ts_si(io_idx_si) = sites(s)%fates_to_bgc_last_ts
          rio_seedrainflux_si(io_idx_si)         = sites(s)%tot_seed_rain_flux

          ! Accumulated trunk product
          rio_trunk_product_si(io_idx_si) = sites(s)%resources_management%trunk_product_site
          ! set numpatches for this column

          rio_npatch_si(io_idx_si)  = patchespersite
          
          do i = 1,numWaterMem ! numWaterMem currently 10
             rio_watermem_siwm( io_idx_si_wmem ) = sites(s)%water_memory(i)
             io_idx_si_wmem = io_idx_si_wmem + 1
          end do
          
       enddo
       
       if ( DEBUG ) then
          write(fates_log(),*) 'CLTV total cohorts ',totalCohorts
       end if
       
       return
     end associate
   end subroutine set_restart_vectors

   ! ====================================================================================

   subroutine create_patchcohort_structure(this, nc, nsites, sites, bc_in) 

     ! ----------------------------------------------------------------------------------
     ! This subroutine takes a peak at the restart file to determine how to allocate
     ! memory for the state structure, and then makes those allocations. This
     ! subroutine is called prior to the transfer of the restart vectors into the
     ! linked-list state structure.
     ! ---------------------------------------------------------------------------------
     use EDTypesMod,           only : ed_site_type
     use EDTypesMod,           only : ed_cohort_type
     use EDTypesMod,           only : ed_patch_type
     use EDTypesMod,           only : ncwd
     use EDTypesMod,           only : nlevleaf
     use EDTypesMod,           only : nclmax
     use FatesInterfaceMod,    only : fates_maxElementsPerPatch
     use EDTypesMod,           only : maxpft
     use EDTypesMod,           only : area
     use EDPatchDynamicsMod,   only : zero_patch
     use EDGrowthFunctionsMod, only : Dbh
     use EDCohortDynamicsMod,  only : create_cohort
     use EDInitMod,            only : zero_site
     use EDInitMod,            only : init_site_vars
     use EDPatchDynamicsMod,   only : create_patch
     use EDPftvarcon,            only : EDPftvarcon_inst

     ! !ARGUMENTS:
     class(fates_restart_interface_type) , intent(inout) :: this
     integer                     , intent(in)            :: nc
     integer                     , intent(in)            :: nsites
     type(ed_site_type)          , intent(inout), target :: sites(nsites)
     type(bc_in_type)            , intent(in)            :: bc_in(nsites)

     ! local variables
     
     type(ed_patch_type) , pointer     :: newp
     type(ed_cohort_type), allocatable :: temp_cohort
     real(r8)                          :: cwd_ag_local(ncwd)
     real(r8)                          :: cwd_bg_local(ncwd)
     real(r8)                          :: leaf_litter_local(maxpft)
     real(r8)                          :: root_litter_local(maxpft)
     real(r8)                          :: patch_age
     integer                           :: cohortstatus
     integer                           :: s        ! site index
     integer                           :: idx_pa        ! local patch index
     integer                           :: io_idx_si     ! global site index in IO vector
     integer                           :: io_idx_co_1st ! global cohort index in IO vector

     integer                           :: fto
     integer                           :: ft

     ! Dummy arguments used for calling create patch, these will be overwritten before
     ! run-time.  Just used now for allocation.
     cwd_ag_local(:)      = 0.0_r8
     cwd_bg_local(:)      = 0.0_r8
     leaf_litter_local(:) = 0.0_r8
     root_litter_local(:) = 0.0_r8
     patch_age            = 0.0_r8

     ! ----------------------------------------------------------------------------------
     ! We really only need the counts for the number of patches per site
     ! and the number of cohorts per patch. These values tell us how much
     ! space to allocate.
     ! ----------------------------------------------------------------------------------
     
     associate( rio_npatch_si  => this%rvars(ir_npatch_si)%int1d , &
                rio_ncohort_pa => this%rvars(ir_ncohort_pa)%int1d )
            
       do s = 1,nsites
          
          io_idx_si  = this%restart_map(nc)%site_index(s)
          io_idx_co_1st  = this%restart_map(nc)%cohort1_index(s)
          
          call init_site_vars( sites(s) )
          call zero_site( sites(s) )
          
          ! 
          ! set a few items that are necessary on restart for ED but not on the 
          ! restart file
          !
          
          sites(s)%ncd = 0.0_r8

          if ( rio_npatch_si(io_idx_si)<0 .or. rio_npatch_si(io_idx_si) > 10000 ) then
             write(fates_log(),*) 'a column was expected to contain a valid number of patches'
             write(fates_log(),*) '0 is a valid number, but this column seems uninitialized',rio_npatch_si(io_idx_si)
             call endrun(msg=errMsg(sourcefile, __LINE__))
          end if
       
          ! Initialize the site pointers to null
          sites(s)%youngest_patch         => null()                 
          sites(s)%oldest_patch           => null()

          do idx_pa = 1,rio_npatch_si(io_idx_si)

             if ( DEBUG ) then
                write(fates_log(),*) 'create patch ',idx_pa
                write(fates_log(),*) 'idx_pa 1-cohortsperpatch : ', rio_ncohort_pa( io_idx_co_1st )
             end if
             
             ! create patch
             allocate(newp)    
             
             ! make new patch
             call create_patch(sites(s), newp, patch_age, area, &
                  cwd_ag_local, cwd_bg_local,  &
                  leaf_litter_local, root_litter_local) 
             
             newp%siteptr => sites(s)

             ! give this patch a unique patch number
             newp%patchno = idx_pa

             do fto = 1, rio_ncohort_pa( io_idx_co_1st )

                allocate(temp_cohort)
                
                temp_cohort%n = 700.0_r8
                temp_cohort%balive = 0.0_r8
                temp_cohort%bdead = 0.0_r8
                temp_cohort%bstore = 0.0_r8
                temp_cohort%laimemory = 0.0_r8
                temp_cohort%canopy_trim = 0.0_r8
                temp_cohort%canopy_layer = 1.0_r8
                temp_cohort%canopy_layer_yesterday = 1.0_r8

                ! set the pft (only 2 used in ed) based on odd/even cohort
                ! number
                ft=2
                if ((mod(fto, 2)  ==  0 )) then
                   ft=1
                endif
                temp_cohort%pft = ft

                cohortstatus = newp%siteptr%status

                if(EDPftvarcon_inst%stress_decid(ft) == 1)then !drought decidous, override status. 
                   cohortstatus = newp%siteptr%dstatus
                endif

                temp_cohort%hite = 1.25_r8
                ! the dbh function should only take as an argument, the one
                ! item it needs, not the entire cohort...refactor
                temp_cohort%dbh = Dbh(temp_cohort) + 0.0001_r8*ft
                
                if (DEBUG) then
                   write(fates_log(),*) 'EDRestVectorMod.F90::createPatchCohortStructure call create_cohort '
                end if
                
                call create_cohort(newp, ft, temp_cohort%n, temp_cohort%hite, temp_cohort%dbh, &
                     temp_cohort%balive, temp_cohort%bdead, temp_cohort%bstore,  &
                     temp_cohort%laimemory, cohortstatus, temp_cohort%canopy_trim, newp%NCL_p, &
                     bc_in(s))
                
                deallocate(temp_cohort)
                
             enddo ! ends loop over fto
             
             !
             ! insert this patch with cohorts into the site pointer.  At this
             ! point just insert the new patch in the youngest position
             !
             if (idx_pa == 1) then ! nothing associated yet. first patch is pointed to by youngest and oldest
                
                if ( DEBUG ) write(fates_log(),*) 'idx_pa = 1 ',idx_pa
                
                sites(s)%youngest_patch         => newp                   
                sites(s)%oldest_patch           => newp                        
                sites(s)%youngest_patch%younger => null()
                sites(s)%youngest_patch%older   => null()
                sites(s)%oldest_patch%younger   => null()
                sites(s)%oldest_patch%older     => null()
                
             else if (idx_pa == 2) then ! add second patch to list
                
                if ( DEBUG ) write(fates_log(),*) 'idx_pa = 2 ',idx_pa
                
                sites(s)%youngest_patch         => newp
                sites(s)%youngest_patch%younger => null()
                sites(s)%youngest_patch%older   => sites(s)%oldest_patch
                sites(s)%oldest_patch%younger   => sites(s)%youngest_patch
                sites(s)%oldest_patch%older     => null()

             else ! more than 2 patches, insert patch into youngest slot
                
                if ( DEBUG ) write(fates_log(),*) 'idx_pa > 2 ',idx_pa
                
                newp%older                      => sites(s)%youngest_patch
                sites(s)%youngest_patch%younger => newp
                newp%younger                    => null()
                sites(s)%youngest_patch         => newp
                
             endif
             
             io_idx_co_1st = io_idx_co_1st + fates_maxElementsPerPatch

          enddo ! ends loop over idx_pa

       enddo ! ends loop over s
       
     end associate
   end subroutine create_patchcohort_structure
   
   ! ====================================================================================

   subroutine get_restart_vectors(this, nc, nsites, sites)

     use EDTypesMod, only : ed_site_type
     use EDTypesMod, only : ed_cohort_type
     use EDTypesMod, only : ed_patch_type
     use EDTypesMod, only : ncwd
     use EDTypesMod, only : nlevleaf
     use EDTypesMod, only : nclmax
     use FatesInterfaceMod, only : numpft
     use FatesInterfaceMod, only : fates_maxElementsPerPatch
     use EDTypesMod, only : numWaterMem

     ! !ARGUMENTS:
     class(fates_restart_interface_type) , intent(inout) :: this
     integer                     , intent(in)            :: nc
     integer                     , intent(in)            :: nsites
     type(ed_site_type)          , intent(inout), target :: sites(nsites)


     ! locals
     ! ----------------------------------------------------------------------------------
     ! LL pointers
     type(ed_patch_type),pointer  :: cpatch      ! current patch
     type(ed_cohort_type),pointer :: ccohort     ! current cohort

     ! loop indices
     integer :: s, i, j, k

     ! ----------------------------------------------------------------------------------
     ! The following group of integers indicate the positional index (idx)
     ! of variables at different scales inside the I/O arrays (io)
     ! Keep in mind that many of these variables have a composite dimension
     ! at the patch scale.  To hold this memory, we borrow the cohort
     ! vector.  Thus the head of each array points to the first cohort
     ! of each patch. "io_idx_co_1st"
     ! ----------------------------------------------------------------------------------
     integer  :: io_idx_si      ! site index
     integer  :: io_idx_co_1st  ! 1st cohort of each patch
     integer  :: io_idx_co      ! cohort index
     integer  :: io_idx_pa_pft  ! each pft within each patch (pa_pft)
     integer  :: io_idx_pa_cwd  ! each cwd class within each patch (pa_cwd)
     integer  :: io_idx_pa_sunz ! index for the combined dimensions for radiation
     integer  :: io_idx_si_wmem ! each water memory class within each site

     ! Some counters (for checking mostly)
     integer  :: totalcohorts   ! total cohort count on this thread (diagnostic)
     integer  :: patchespersite   ! number of patches per site
     integer  :: cohortsperpatch  ! number of cohorts per patch 
     


     associate( rio_npatch_si         => this%rvars(ir_npatch_si)%int1d, &
          rio_old_stock_si            => this%rvars(ir_oldstock_si)%r81d, &
          rio_cd_status_si            => this%rvars(ir_cd_status_si)%r81d, &
          rio_dd_status_si            => this%rvars(ir_dd_status_si)%r81d, &
          rio_nchill_days_si          => this%rvars(ir_nchill_days_si)%r81d, &
          rio_leafondate_si           => this%rvars(ir_leafondate_si)%r81d, &
          rio_leafoffdate_si          => this%rvars(ir_leafoffdate_si)%r81d, &
          rio_dleafondate_si          => this%rvars(ir_dleafondate_si)%r81d, &
          rio_dleafoffdate_si         => this%rvars(ir_dleafoffdate_si)%r81d, &
          rio_acc_ni_si               => this%rvars(ir_acc_ni_si)%r81d, &
          rio_gdd_si                  => this%rvars(ir_gdd_si)%r81d, &
          rio_nep_timeintegrated_si   => this%rvars(ir_nep_timeintegrated_si)%r81d, &
          rio_npp_timeintegrated_si   => this%rvars(ir_npp_timeintegrated_si)%r81d, &
          rio_hr_timeintegrated_si    => this%rvars(ir_hr_timeintegrated_si)%r81d, &
          rio_cbal_err_fates_si       => this%rvars(ir_cbal_error_fates_si)%r81d, &
          rio_cbal_err_bgc_si         => this%rvars(ir_cbal_error_bgc_si)%r81d, &
          rio_cbal_err_tot_si         => this%rvars(ir_cbal_error_total_si)%r81d, &
          rio_totecosysc_old_si       => this%rvars(ir_totecosysc_old_si)%r81d, &
          rio_totfatesc_old_si        => this%rvars(ir_totfatesc_old_si)%r81d, &
          rio_totbgcc_old_si          => this%rvars(ir_totbgcc_old_si)%r81d, &
          rio_fates_to_bgc_this_ts_si => this%rvars(ir_fates_to_bgc_this_ts_si)%r81d, &
          rio_fates_to_bgc_last_ts_si => this%rvars(ir_fates_to_bgc_last_ts_si)%r81d, &
          rio_seedrainflux_si         => this%rvars(ir_seedrainflux_si)%r81d, &
          rio_trunk_product_si        => this%rvars(ir_trunk_product_si)%r81d, &
          rio_ncohort_pa              => this%rvars(ir_ncohort_pa)%int1d, &
          rio_balive_co               => this%rvars(ir_balive_co)%r81d, &
          rio_bdead_co                => this%rvars(ir_bdead_co)%r81d, &
          rio_bleaf_co                => this%rvars(ir_bleaf_co)%r81d, &
          rio_broot_co                => this%rvars(ir_broot_co)%r81d, &
          rio_bstore_co               => this%rvars(ir_bstore_co)%r81d, &
          rio_canopy_layer_co         => this%rvars(ir_canopy_layer_co)%r81d, &
          rio_canopy_layer_yesterday_co         => this%rvars(ir_canopy_layer_yesterday_co)%r81d, &
          rio_canopy_trim_co          => this%rvars(ir_canopy_trim_co)%r81d, &
          rio_dbh_co                  => this%rvars(ir_dbh_co)%r81d, &
          rio_height_co               => this%rvars(ir_height_co)%r81d, &
          rio_laimemory_co            => this%rvars(ir_laimemory_co)%r81d, &
          rio_leaf_md_co              => this%rvars(ir_leaf_md_co)%r81d, &
          rio_root_md_co              => this%rvars(ir_root_md_co)%r81d, &
          rio_nplant_co               => this%rvars(ir_nplant_co)%r81d, &
          rio_gpp_acc_co              => this%rvars(ir_gpp_acc_co)%r81d, &
          rio_npp_acc_co              => this%rvars(ir_npp_acc_co)%r81d, &
          rio_gpp_acc_hold_co         => this%rvars(ir_gpp_acc_hold_co)%r81d, &
          rio_npp_acc_hold_co         => this%rvars(ir_npp_acc_hold_co)%r81d, &
          rio_npp_leaf_co             => this%rvars(ir_npp_leaf_co)%r81d, &
          rio_npp_froot_co            => this%rvars(ir_npp_froot_co)%r81d, &
          rio_npp_sw_co               => this%rvars(ir_npp_sw_co)%r81d, &
          rio_npp_dead_co             => this%rvars(ir_npp_dead_co)%r81d, &
          rio_npp_seed_co             => this%rvars(ir_npp_seed_co)%r81d, &
          rio_npp_store_co            => this%rvars(ir_npp_store_co)%r81d, &
          rio_bmort_co                => this%rvars(ir_bmort_co)%r81d, &
          rio_hmort_co                => this%rvars(ir_hmort_co)%r81d, &
          rio_cmort_co                => this%rvars(ir_cmort_co)%r81d, &
          rio_imort_co                => this%rvars(ir_imort_co)%r81d, &
          rio_fmort_co                => this%rvars(ir_fmort_co)%r81d, &

	  rio_lmort_logging_co                => this%rvars(ir_lmort_logging_co)%r81d, &
	  rio_lmort_collateral_co                => this%rvars(ir_lmort_collateral_co)%r81d, &
	  rio_lmort_infra_co                => this%rvars(ir_lmort_infra_co)%r81d, &



          rio_ddbhdt_co               => this%rvars(ir_ddbhdt_co)%r81d, &
          rio_dbalivedt_co            => this%rvars(ir_dbalivedt_co)%r81d, &
          rio_dbdeaddt_co             => this%rvars(ir_dbdeaddt_co)%r81d, &
          rio_dbstoredt_co            => this%rvars(ir_dbstoredt_co)%r81d, &
          rio_resp_tstep_co           => this%rvars(ir_resp_tstep_co)%r81d, &
          rio_pft_co                  => this%rvars(ir_pft_co)%int1d, &
          rio_status_co               => this%rvars(ir_status_co)%int1d, &
          rio_isnew_co                => this%rvars(ir_isnew_co)%int1d, &
          rio_cwd_ag_pacw             => this%rvars(ir_cwd_ag_pacw)%r81d, &
          rio_cwd_bg_pacw             => this%rvars(ir_cwd_bg_pacw)%r81d, &
          rio_leaf_litter_paft        => this%rvars(ir_leaf_litter_paft)%r81d, &
          rio_root_litter_paft        => this%rvars(ir_root_litter_paft)%r81d, &
          rio_leaf_litter_in_paft     => this%rvars(ir_leaf_litter_in_paft)%r81d, &
          rio_root_litter_in_paft     => this%rvars(ir_root_litter_in_paft)%r81d, &
          rio_seed_bank_sift          => this%rvars(ir_seed_bank_sift)%r81d, &
          rio_spread_si               => this%rvars(ir_spread_si)%r81d, &
          rio_livegrass_pa            => this%rvars(ir_livegrass_pa)%r81d, &
          rio_age_pa                  => this%rvars(ir_age_pa)%r81d, &
          rio_area_pa                 => this%rvars(ir_area_pa)%r81d, &
          rio_fsun_paclftls           => this%rvars(ir_fsun_paclftls)%r81d, &
          rio_fabd_sun_z_paclftls     => this%rvars(ir_fabd_sun_paclftls)%r81d, &
          rio_fabi_sun_z_paclftls     => this%rvars(ir_fabi_sun_paclftls)%r81d, &
          rio_fabd_sha_z_paclftls     => this%rvars(ir_fabd_sha_paclftls)%r81d, &
          rio_fabi_sha_z_paclftls     => this%rvars(ir_fabi_sha_paclftls)%r81d, &
          rio_watermem_siwm           => this%rvars(ir_watermem_siwm)%r81d )
     
       totalcohorts = 0
     
       do s = 1,nsites
          
          io_idx_si      = this%restart_map(nc)%site_index(s)
          io_idx_co_1st  = this%restart_map(nc)%cohort1_index(s)
          
          io_idx_co      = io_idx_co_1st
          io_idx_pa_pft  = io_idx_co_1st
          io_idx_pa_cwd  = io_idx_co_1st
          io_idx_pa_sunz = io_idx_co_1st
          io_idx_si_wmem = io_idx_co_1st
          
          ! read seed_bank info(site-level, but PFT-resolved)
          do i = 1,numpft 
             sites(s)%seed_bank(i) = rio_seed_bank_sift(io_idx_co_1st+i-1)
          enddo

          sites(s)%spread = rio_spread_si(io_idx_si) 
          
          ! Perform a check on the number of patches per site
          patchespersite = 0
          
          cpatch => sites(s)%oldest_patch
          do while(associated(cpatch))
             
             patchespersite = patchespersite + 1
             
             ccohort => cpatch%shortest
             
             ! new patch, reset num cohorts
             cohortsperpatch = 0
             
             do while(associated(ccohort))        
                
                ! found cohort, increment
                cohortsperpatch  = cohortsperpatch    + 1
                totalcohorts     = totalcohorts + 1
                
                if ( DEBUG ) then
                   write(fates_log(),*) 'CVTL io_idx_co ',io_idx_co
                endif

                ccohort%balive       = rio_balive_co(io_idx_co)
                ccohort%bdead        = rio_bdead_co(io_idx_co)
                ccohort%bl           = rio_bleaf_co(io_idx_co)
                ccohort%br           = rio_broot_co(io_idx_co)
                ccohort%bstore       = rio_bstore_co(io_idx_co)
                ccohort%canopy_layer = rio_canopy_layer_co(io_idx_co)
                ccohort%canopy_layer_yesterday = rio_canopy_layer_yesterday_co(io_idx_co)
                ccohort%canopy_trim  = rio_canopy_trim_co(io_idx_co)
                ccohort%dbh          = rio_dbh_co(io_idx_co)
                ccohort%hite         = rio_height_co(io_idx_co)
                ccohort%laimemory    = rio_laimemory_co(io_idx_co)
                ccohort%leaf_md      = rio_leaf_md_co(io_idx_co)
                ccohort%root_md      = rio_root_md_co(io_idx_co)
                ccohort%n            = rio_nplant_co(io_idx_co)
                ccohort%gpp_acc      = rio_gpp_acc_co(io_idx_co)
                ccohort%npp_acc      = rio_npp_acc_co(io_idx_co)
                ccohort%gpp_acc_hold = rio_gpp_acc_hold_co(io_idx_co)
                ccohort%npp_acc_hold = rio_npp_acc_hold_co(io_idx_co)
                ccohort%npp_leaf     = rio_npp_leaf_co(io_idx_co)
                ccohort%npp_froot    = rio_npp_froot_co(io_idx_co)
                ccohort%npp_bsw      = rio_npp_sw_co(io_idx_co)
                ccohort%npp_bdead    = rio_npp_dead_co(io_idx_co)
                ccohort%npp_bseed    = rio_npp_seed_co(io_idx_co)
                ccohort%npp_store    = rio_npp_store_co(io_idx_co)
                ccohort%bmort        = rio_bmort_co(io_idx_co)
                ccohort%hmort        = rio_hmort_co(io_idx_co)
                ccohort%cmort        = rio_cmort_co(io_idx_co)
                ccohort%imort        = rio_imort_co(io_idx_co)
                ccohort%fmort        = rio_fmort_co(io_idx_co)

		!Logging
		ccohort%lmort_logging        = rio_lmort_logging_co(io_idx_co)
		ccohort%lmort_collateral     = rio_lmort_collateral_co(io_idx_co)
		ccohort%lmort_infra        = rio_lmort_infra_co(io_idx_co)



                ccohort%ddbhdt       = rio_ddbhdt_co(io_idx_co)
                ccohort%dbalivedt    = rio_dbalivedt_co(io_idx_co)
                ccohort%dbdeaddt     = rio_dbdeaddt_co(io_idx_co)
                ccohort%dbstoredt    = rio_dbstoredt_co(io_idx_co)
                ccohort%resp_tstep   = rio_resp_tstep_co(io_idx_co)
                ccohort%pft          = rio_pft_co(io_idx_co)
                ccohort%status_coh   = rio_status_co(io_idx_co)
                ccohort%isnew        = ( rio_isnew_co(io_idx_co) .eq. new_cohort )
                
                io_idx_co = io_idx_co + 1
             
                ccohort => ccohort%taller
                
             enddo ! current cohort do while

             if(cohortsperpatch .ne. rio_ncohort_pa(io_idx_co_1st)) then
                write(fates_log(),*) 'Number of cohorts per patch during retrieval'
                write(fates_log(),*) 'does not match allocation'
                call endrun(msg=errMsg(sourcefile, __LINE__))
             end if

          
             ! FIX(SPM,032414) move to init if you can...or make a new init function
             cpatch%leaf_litter(:)    = 0.0_r8
             cpatch%root_litter(:)    = 0.0_r8
             cpatch%leaf_litter_in(:) = 0.0_r8
             cpatch%root_litter_in(:) = 0.0_r8
             
             !
             ! deal with patch level fields here
             !
             cpatch%livegrass  = rio_livegrass_pa(io_idx_co_1st)
             cpatch%age        = rio_age_pa(io_idx_co_1st) 
             cpatch%area       = rio_area_pa(io_idx_co_1st) 
             
             ! set cohorts per patch for IO
             
             if ( DEBUG ) then
                write(fates_log(),*) 'CVTL III ' &
                     ,io_idx_co,cohortsperpatch
             endif
             !
             ! deal with patch level fields of arrays here
             !
             ! these are arrays of length numpft, each patch contains one
             ! vector so we increment 

             do i = 1,numpft
                cpatch%leaf_litter(i)    = rio_leaf_litter_paft(io_idx_pa_pft)    
                cpatch%root_litter(i)    = rio_root_litter_paft(io_idx_pa_pft)    
                cpatch%leaf_litter_in(i) = rio_leaf_litter_in_paft(io_idx_pa_pft) 
                cpatch%root_litter_in(i) = rio_root_litter_in_paft(io_idx_pa_pft) 
                io_idx_pa_pft = io_idx_pa_pft + 1
             enddo
          
             do i = 1,ncwd ! ncwd currently 4
                cpatch%cwd_ag(i) = rio_cwd_ag_pacw(io_idx_pa_cwd)
                cpatch%cwd_bg(i) = rio_cwd_bg_pacw(io_idx_pa_cwd)
                io_idx_pa_cwd = io_idx_pa_cwd + 1
             enddo
             
             if ( DEBUG ) write(fates_log(),*) 'CVTL io_idx_pa_sunz 1 ',io_idx_pa_sunz
             
             do k = 1,nlevleaf ! nlevleaf currently 40
                do j = 1,numpft
                   do i = 1,nclmax ! nclmax currently 2
                      cpatch%f_sun(i,j,k)      = rio_fsun_paclftls(io_idx_pa_sunz) 
                      cpatch%fabd_sun_z(i,j,k) = rio_fabd_sun_z_paclftls(io_idx_pa_sunz)
                      cpatch%fabi_sun_z(i,j,k) = rio_fabi_sun_z_paclftls(io_idx_pa_sunz)
                      cpatch%fabd_sha_z(i,j,k) = rio_fabd_sha_z_paclftls(io_idx_pa_sunz)
                      cpatch%fabi_sha_z(i,j,k) = rio_fabi_sha_z_paclftls(io_idx_pa_sunz)
                      io_idx_pa_sunz = io_idx_pa_sunz + 1
                   end do
                end do
             end do
             
             if ( DEBUG ) write(fates_log(),*) 'CVTL io_idx_pa_sunz 2 ',io_idx_pa_sunz
             
             ! Now increment the position of the first cohort to that of the next
             ! patch
             
             io_idx_co_1st = io_idx_co_1st + fates_maxElementsPerPatch
             
             ! and max the number of allowed cohorts per patch
             io_idx_pa_pft  = io_idx_co_1st
             io_idx_pa_cwd  = io_idx_co_1st
             io_idx_co      = io_idx_co_1st
             io_idx_pa_sunz = io_idx_co_1st
             
             if ( DEBUG ) then
                write(fates_log(),*) 'CVTL io_idx_co_1st ', io_idx_co_1st
                write(fates_log(),*) 'CVTL cohortsperpatch ', cohortsperpatch
                write(fates_log(),*) 'CVTL totalCohorts ', totalCohorts
             end if
             
             cpatch => cpatch%younger
             
          enddo ! patch do while
          
          if(patchespersite .ne. rio_npatch_si(io_idx_si)) then
             write(fates_log(),*) 'Number of patches per site during retrieval does not match allocation'
             call endrun(msg=errMsg(sourcefile, __LINE__))
          end if
          
          do i = 1,numWaterMem
             sites(s)%water_memory(i) = rio_watermem_siwm( io_idx_si_wmem )
             io_idx_si_wmem = io_idx_si_wmem + 1
          end do
          
          sites(s)%old_stock      = rio_old_stock_si(io_idx_si)
          sites(s)%status         = rio_cd_status_si(io_idx_si)
          sites(s)%dstatus        = rio_dd_status_si(io_idx_si)
          sites(s)%ncd            = rio_nchill_days_si(io_idx_si)
          sites(s)%leafondate     = rio_leafondate_si(io_idx_si)
          sites(s)%leafoffdate    = rio_leafoffdate_si(io_idx_si)
          sites(s)%dleafondate    = rio_dleafondate_si(io_idx_si)
          sites(s)%dleafoffdate   = rio_dleafoffdate_si(io_idx_si)
          sites(s)%acc_NI         = rio_acc_ni_si(io_idx_si)
          sites(s)%ED_GDD_site    = rio_gdd_si(io_idx_si)

          ! Carbon Balance and Checks
          sites(s)%nep_timeintegrated   = rio_nep_timeintegrated_si(io_idx_si)
          sites(s)%npp_timeintegrated   = rio_npp_timeintegrated_si(io_idx_si)
          sites(s)%hr_timeintegrated    = rio_hr_timeintegrated_si(io_idx_si)
          sites(s)%totecosysc_old       = rio_totecosysc_old_si(io_idx_si)
          sites(s)%totfatesc_old        = rio_totfatesc_old_si(io_idx_si)
          sites(s)%totbgcc_old          = rio_totbgcc_old_si(io_idx_si)
          sites(s)%cbal_err_fates       = rio_cbal_err_fates_si(io_idx_si)
          sites(s)%cbal_err_bgc         = rio_cbal_err_bgc_si(io_idx_si)
          sites(s)%cbal_err_tot         = rio_cbal_err_tot_si(io_idx_si)
          sites(s)%fates_to_bgc_this_ts = rio_fates_to_bgc_this_ts_si(io_idx_si)
          sites(s)%fates_to_bgc_last_ts = rio_fates_to_bgc_last_ts_si(io_idx_si)
          sites(s)%tot_seed_rain_flux   = rio_seedrainflux_si(io_idx_si)
          sites(s)%resources_management%trunk_product_site = rio_trunk_product_si(io_idx_si)

       end do
       
       if ( DEBUG ) then
          write(fates_log(),*) 'CVTL total cohorts ',totalCohorts
       end if
       
     end associate
   end subroutine get_restart_vectors
   
 end module FatesRestartInterfaceMod
