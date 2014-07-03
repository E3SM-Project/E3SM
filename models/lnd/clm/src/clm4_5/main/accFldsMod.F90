module accFldsMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: accFldsMod
!
! !DESCRIPTION:
! This module contains subroutines that initialize, update and extract
! the user-specified fields over user-defined intervals. Each interval
! and accumulation type is unique to each field processed.
! Subroutine [initAccumFlds] defines the fields to be processed
! and the type of accumulation. Subroutine [updateAccumFlds] does
! the actual accumulation for a given field. Fields are accumulated
! by calls to subroutine [update_accum_field]. To accumulate a field,
! it must first be defined in subroutine [initAccumFlds] and then
! accumulated by calls to [updateAccumFlds].
! Four types of accumulations are possible:
!   o average over time interval
!   o running mean over time interval
!   o running accumulation over time interval
! Time average fields are only valid at the end of the averaging interval.
! Running means are valid once the length of the simulation exceeds the
! averaging interval. Accumulated fields are continuously accumulated.
! The trigger value "-99999." resets the accumulation to zero.
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use abortutils,   only: endrun
  use clm_varctl,   only: iulog
  use surfrdMod,    only: crop_prog
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: initAccFlds     ! Initialization accumulator fields
  public :: initAccClmtype  ! Initialize clmtype variables obtained from accum fields
  public :: updateAccFlds   ! Update accumulator fields
!
! !REVISION HISTORY:
! Created by M. Vertenstein 03/2003
! F. Li and S. Levis (11/06/12) for prec10, prec60
!
!EOP

contains

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: initAccFlds()
!
! !INTERFACE:
  subroutine initAccFlds()
!
! !DESCRIPTION:
! Initializes accumulator and sets up array of accumulated fields
!
! !USES:
    use accumulMod   , only : init_accum_field, print_accum_fields
    use clm_time_manager , only : get_step_size
    use shr_const_mod, only : SHR_CONST_CDAY, SHR_CONST_TKFRZ
    use nanMod       , only : bigint
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY::
! Created by M. Vertenstein 03/2003
!
!
! !LOCAL VARIABLES:
!EOP
!
    integer :: dtime                     !time step size
    integer, parameter :: not_used = bigint
!------------------------------------------------------------------------

    ! Hourly average of 2m temperature.

    dtime = get_step_size()
    call init_accum_field(name='TREFAV', units='K', &
         desc='average over an hour of 2-m temperature', &
         accum_type='timeavg', accum_period=nint(3600._r8/dtime), &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    ! Hourly average of Urban 2m temperature.

    call init_accum_field(name='TREFAV_U', units='K', &
         desc='average over an hour of urban 2-m temperature', &
         accum_type='timeavg', accum_period=nint(3600._r8/dtime), &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    ! Hourly average of Rural 2m temperature.

    call init_accum_field(name='TREFAV_R', units='K', &
         desc='average over an hour of rural 2-m temperature', &
         accum_type='timeavg', accum_period=nint(3600._r8/dtime), &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    ! 24hr average of vegetation temperature (heald, 04/06)
    call init_accum_field (name='T_VEG24', units='K', &
         desc='24hr average of vegetation temperature', &
         accum_type='runmean', accum_period=-1, &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    ! 240hr average of vegetation temperature (heald, 04/06)
    call init_accum_field (name='T_VEG240', units='K', &
         desc='240hr average of vegetation temperature', &
         accum_type='runmean', accum_period=-10, &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    ! 24hr average of direct solar radiation (heald, 04/06)
    call init_accum_field (name='FSD24', units='W/m2', &
         desc='24hr average of direct solar radiation', &
         accum_type='runmean', accum_period=-1, &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    ! 240hr average of direct solar radiation (heald, 04/06)
    call init_accum_field (name='FSD240', units='W/m2', &
         desc='240hr average of direct solar radiation', &
         accum_type='runmean', accum_period=-10, &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    ! 24hr average of diffuse solar radiation (heald, 04/06)
    call init_accum_field (name='FSI24', units='W/m2', &
         desc='24hr average of diffuse solar radiation', &
         accum_type='runmean', accum_period=-1, &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    ! 240hr average of diffuse solar radiation (heald, 04/06)
    call init_accum_field (name='FSI240', units='W/m2', &
         desc='240hr average of diffuse solar radiation', &
         accum_type='runmean', accum_period=-10, &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    ! 24hr average of fraction of canopy that is sunlit (heald, 04/06)
    call init_accum_field (name='FSUN24', units='fraction', &
         desc='24hr average of diffuse solar radiation', &
         accum_type='runmean', accum_period=-1, &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    ! 240hr average of fraction of canopy that is sunlit (heald, 04/06)
    call init_accum_field (name='FSUN240', units='fraction', &
         desc='240hr average of diffuse solar radiation', &
         accum_type='runmean', accum_period=-10, &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    ! Average of LAI from previous and current timestep (heald, 04/06)
    call init_accum_field (name='LAIP', units='m2/m2', &
         desc='leaf area index average over timestep', &
         accum_type='runmean', accum_period=1, &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    ! The following is a running mean.
    ! The accumulation period is set to -10 for a 10-day running mean.

    call init_accum_field (name='T10', units='K', &
         desc='10-day running mean of 2-m temperature', &
         accum_type='runmean', accum_period=-10, &
         subgrid_type='pft', numlev=1,init_value=SHR_CONST_TKFRZ+20._r8)

#if (defined CNDV)
    ! 30-day average of 2m temperature.

    call init_accum_field (name='TDA', units='K', &
         desc='30-day average of 2-m temperature', &
         accum_type='timeavg', accum_period=-30, &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    ! The following are running means.
    ! The accumulation period is set to -365 for a 365-day running mean.

    call init_accum_field (name='PREC365', units='MM H2O/S', &
         desc='365-day running mean of total precipitation', &
         accum_type='runmean', accum_period=-365, &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    ! The following are accumulated fields.
    ! These types of fields are accumulated until a trigger value resets
    ! the accumulation to zero (see subroutine update_accum_field).
    ! Hence, [accper] is not valid.

    call init_accum_field (name='AGDDTW', units='K', &
         desc='growing degree-days base twmax', &
         accum_type='runaccum', accum_period=not_used, &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    call init_accum_field (name='AGDD', units='K', &
         desc='growing degree-days base 5C', &
         accum_type='runaccum', accum_period=not_used,  &
         subgrid_type='pft', numlev=1, init_value=0._r8)
#endif

      call init_accum_field (name='PREC60', units='MM H2O/S', &
         desc='60-day running mean of total precipitation', &
         accum_type='runmean', accum_period=-60, &
         subgrid_type='pft', numlev=1, init_value=0._r8)

         call init_accum_field (name='PREC10', units='MM H2O/S', &
         desc='10-day running mean of total precipitation', &
         accum_type='runmean', accum_period=-10, &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    if ( crop_prog )then
       ! 10-day average of min 2m temperature.

       call init_accum_field (name='TDM10', units='K', &
            desc='10-day running mean of min 2-m temperature', &
            accum_type='runmean', accum_period=-10, &
            subgrid_type='pft', numlev=1, init_value=SHR_CONST_TKFRZ)

       ! 5-day average of min 2m temperature.

       call init_accum_field (name='TDM5', units='K', &
            desc='5-day running mean of min 2-m temperature', &
            accum_type='runmean', accum_period=-5, &
            subgrid_type='pft', numlev=1, init_value=SHR_CONST_TKFRZ)

       ! All GDD summations are relative to the planting date
       ! (Kucharik & Brye 2003)

       call init_accum_field (name='GDD0', units='K', &
            desc='growing degree-days base 0C from planting', &
            accum_type='runaccum', accum_period=not_used, &
            subgrid_type='pft', numlev=1, init_value=0._r8)

       call init_accum_field (name='GDD8', units='K', &
            desc='growing degree-days base 8C from planting', &
            accum_type='runaccum', accum_period=not_used, &
            subgrid_type='pft', numlev=1, init_value=0._r8)

       call init_accum_field (name='GDD10', units='K', &
            desc='growing degree-days base 10C from planting', &
            accum_type='runaccum', accum_period=not_used,  &
            subgrid_type='pft', numlev=1, init_value=0._r8)

       call init_accum_field (name='GDDPLANT', units='K', &
            desc='growing degree-days from planting', &
            accum_type='runaccum', accum_period=not_used,  &
            subgrid_type='pft', numlev=1, init_value=0._r8)

       call init_accum_field (name='GDDTSOI', units='K', &
            desc='growing degree-days from planting (top two soil layers)', &
            accum_type='runaccum', accum_period=not_used,  &
            subgrid_type='pft', numlev=1, init_value=0._r8)
    end if

    ! Print output of accumulated fields

    call print_accum_fields()

  end subroutine initAccFlds

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: updateAccFlds
!
! !INTERFACE:
  subroutine updateAccFlds()
!
! !DESCRIPTION:
! Update and/or extract accumulated fields
!
! !USES:
    use clmtype
    use clm_atmlnd   , only : clm_a2l
    use decompMod    , only : get_proc_bounds
    use clm_varcon   , only : spval
    use shr_const_mod, only : SHR_CONST_CDAY, SHR_CONST_TKFRZ
    use clm_time_manager , only : get_step_size, get_nstep, is_end_curr_day, get_curr_date
    use accumulMod   , only : update_accum_field, extract_accum_field
    use pftvarcon    , only : nwcereal, nwcerealirrig, mxtmp, baset
    use clm_time_manager , only : get_start_date
    use pftvarcon    , only : ndllf_dcd_brl_tree
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by M. Vertenstein 03/2003
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
    integer , pointer :: itype(:)            ! pft vegetation
    integer , pointer :: pgridcell(:)        ! index into gridcell level quantities
    real(r8), pointer :: forc_t(:)           ! atmospheric temperature (Kelvin)
    real(r8), pointer :: forc_rain(:)        ! rain rate [mm/s]
    real(r8), pointer :: forc_snow(:)        ! snow rate [mm/s]
    real(r8), pointer :: t_ref2m(:)          ! 2 m height surface air temperature (Kelvin)
    real(r8), pointer :: t_ref2m_u(:)        ! Urban 2 m height surface air temperature (Kelvin)
    real(r8), pointer :: t_ref2m_r(:)        ! Rural 2 m height surface air temperature (Kelvin)
    logical , pointer :: urbpoi(:)           ! true => landunit is an urban point
    logical , pointer :: ifspecial(:)        ! true => landunit is not vegetated
    integer , pointer :: plandunit(:)        ! landunit index associated with each pft
    real(r8), pointer :: vf(:)               ! vernalization factor
    real(r8), pointer :: t_soisno(:,:)       ! soil temperature (K)
    real(r8), pointer :: h2osoi_liq(:,:)     ! liquid water (kg/m2)
    real(r8), pointer :: watsat(:,:)         ! volumetric soil water at saturation (porosity) (nlevgrnd)
    real(r8), pointer :: dz(:,:)             ! layer thickness depth (m)
    real(r8), pointer :: latdeg(:)           ! latitude (radians)
    logical , pointer :: croplive(:)         ! Flag, true if planted, not harvested
    integer , pointer :: pcolumn(:)          ! index into column level quantities
!
! local pointers to implicit out arguments
!
    ! heald (04/06): variables to be accumulated for VOC emissions
    real(r8), pointer :: t_veg(:)            ! pft vegetation temperature (Kelvin) 
    real(r8), pointer :: forc_solad(:,:)     ! direct beam radiation (visible only)
    real(r8), pointer :: forc_solai(:,:)     ! diffuse radiation     (visible only)
    real(r8), pointer :: fsun(:)             ! sunlit fraction of canopy 
    real(r8), pointer :: elai(:)             ! one-sided leaf area index with burying by snow 
    ! heald (04/06): accumulated variables for VOC emissions
    real(r8), pointer :: t_veg24(:)          ! 24hr average vegetation temperature (K)
    real(r8), pointer :: t_veg240(:)         ! 240hr average vegetation temperature (Kelvin)
    real(r8), pointer :: fsd24(:)            ! 24hr average of direct beam radiation 
    real(r8), pointer :: fsd240(:)           ! 240hr average of direct beam radiation 
    real(r8), pointer :: fsi24(:)            ! 24hr average of diffuse beam radiation 
    real(r8), pointer :: fsi240(:)           ! 240hr average of diffuse beam radiation 
    real(r8), pointer :: fsun24(:)           ! 24hr average of sunlit fraction of canopy 
    real(r8), pointer :: fsun240(:)          ! 240hr average of sunlit fraction of canopy
    real(r8), pointer :: elai_p(:)           ! leaf area index average over timestep 

    real(r8), pointer :: t_ref2m_min(:)      ! daily minimum of average 2 m height surface air temperature (K)
    real(r8), pointer :: t_ref2m_max(:)      ! daily maximum of average 2 m height surface air temperature (K)
    real(r8), pointer :: t_ref2m_min_inst(:) ! instantaneous daily min of average 2 m height surface air temp (K)
    real(r8), pointer :: t_ref2m_max_inst(:) ! instantaneous daily max of average 2 m height surface air temp (K)
    real(r8), pointer :: t_ref2m_min_u(:)    ! Urban daily minimum of average 2 m height surface air temperature (K)
    real(r8), pointer :: t_ref2m_min_r(:)    ! Rural daily minimum of average 2 m height surface air temperature (K)
    real(r8), pointer :: t_ref2m_max_u(:)    ! Urban daily maximum of average 2 m height surface air temperature (K)
    real(r8), pointer :: t_ref2m_max_r(:)    ! Rural daily maximum of average 2 m height surface air temperature (K)
    real(r8), pointer :: t_ref2m_min_inst_u(:) ! Urban instantaneous daily min of average 2 m height surface air temp (K)
    real(r8), pointer :: t_ref2m_min_inst_r(:) ! Rural instantaneous daily min of average 2 m height surface air temp (K)
    real(r8), pointer :: t_ref2m_max_inst_u(:) ! Urban instantaneous daily max of average 2 m height surface air temp (K)
    real(r8), pointer :: t_ref2m_max_inst_r(:) ! Rural instantaneous daily max of average 2 m height surface air temp (K)
    real(r8), pointer :: t10(:)              ! 10-day running mean of the 2 m temperature (K)
#if (defined CNDV)
    real(r8), pointer :: t_mo(:)             ! 30-day average temperature (Kelvin)
    real(r8), pointer :: t_mo_min(:)         ! annual min of t_mo (Kelvin)
    real(r8), pointer :: prec365(:)          ! 365-day running mean of tot. precipitation
    real(r8), pointer :: agddtw(:)           ! accumulated growing degree days above twmax
    real(r8), pointer :: agdd(:)             ! accumulated growing degree days above 5
    real(r8), pointer :: twmax(:)            ! upper limit of temperature of the warmest month
#endif
     real(r8), pointer :: prec10(:)          ! 10-day running mean of tot. precipitation
    real(r8), pointer :: prec60(:)          ! 60-day running mean of tot. precipitation
    real(r8), pointer :: gdd0(:)             ! growing degree-days base 0C'
    real(r8), pointer :: gdd8(:)             ! growing degree-days base 8C from planting
    real(r8), pointer :: gdd10(:)            ! growing degree-days base 10C from planting
    real(r8), pointer :: gddplant(:)         ! growing degree-days from planting
    real(r8), pointer :: gddtsoi(:)          ! growing degree-days from planting (top two soil layers)
    real(r8), pointer :: a10tmin(:)          ! 10-day running mean of min 2-m temperature
    real(r8), pointer :: a5tmin(:)           ! 5-day running mean of min 2-m temperature
!
!
! !OTHER LOCAL VARIABLES:
!EOP
    integer :: g,l,c,p                   ! indices
    integer :: itypveg                   ! vegetation type
    integer :: dtime                     ! timestep size [seconds]
    integer :: nstep                     ! timestep number
    integer :: year                      ! year (0, ...) for nstep
    integer :: month                     ! month (1, ..., 12) for nstep
    integer :: day                       ! day of month (1, ..., 31) for nstep
    integer :: secs                      ! seconds into current date for nstep
    logical :: end_cd                    ! temporary for is_end_curr_day() value
    integer :: ier                       ! error status
    integer :: begp, endp                !  per-proc beginning and ending pft indices
    integer :: begc, endc                !  per-proc beginning and ending column indices
    integer :: begl, endl                !  per-proc beginning and ending landunit indices
    integer :: begg, endg                !  per-proc gridcell ending gridcell indices
    real(r8), pointer :: rbufslp(:)      ! temporary single level - pft level
!------------------------------------------------------------------------

    ! Determine necessary indices

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    ! Assign local pointers to derived subtypes components (gridcell-level)

    forc_t     => clm_a2l%forc_t
    forc_rain  => clm_a2l%forc_rain
    forc_snow  => clm_a2l%forc_snow
    forc_solad => clm_a2l%forc_solad 	    ! (heald 04/06)
    forc_solai => clm_a2l%forc_solai 	    ! (heald 04/06)

    ! Assign local pointers to derived subtypes components (landunit-level)
    ifspecial  => clm3%g%l%ifspecial
    urbpoi     => clm3%g%l%urbpoi

    ! Assign local pointers to derived subtypes components (pft-level)

    itype            => clm3%g%l%c%p%itype
    pgridcell        => clm3%g%l%c%p%gridcell
    t_ref2m          => clm3%g%l%c%p%pes%t_ref2m
    t_ref2m_max_inst => clm3%g%l%c%p%pes%t_ref2m_max_inst
    t_ref2m_min_inst => clm3%g%l%c%p%pes%t_ref2m_min_inst
    t_ref2m_max      => clm3%g%l%c%p%pes%t_ref2m_max
    t_ref2m_min      => clm3%g%l%c%p%pes%t_ref2m_min
    t_ref2m_u        => clm3%g%l%c%p%pes%t_ref2m_u
    t_ref2m_r        => clm3%g%l%c%p%pes%t_ref2m_r
    t_ref2m_max_u    => clm3%g%l%c%p%pes%t_ref2m_max_u
    t_ref2m_max_r    => clm3%g%l%c%p%pes%t_ref2m_max_r
    t_ref2m_min_u    => clm3%g%l%c%p%pes%t_ref2m_min_u
    t_ref2m_min_r    => clm3%g%l%c%p%pes%t_ref2m_min_r
    t_ref2m_max_inst_u => clm3%g%l%c%p%pes%t_ref2m_max_inst_u
    t_ref2m_max_inst_r => clm3%g%l%c%p%pes%t_ref2m_max_inst_r
    t_ref2m_min_inst_u => clm3%g%l%c%p%pes%t_ref2m_min_inst_u
    t_ref2m_min_inst_r => clm3%g%l%c%p%pes%t_ref2m_min_inst_r
    plandunit        => clm3%g%l%c%p%landunit
    t10              => clm3%g%l%c%p%pes%t10
    a10tmin          => clm3%g%l%c%p%pes%a10tmin
    a5tmin           => clm3%g%l%c%p%pes%a5tmin
#if (defined CNDV)
    t_mo             => clm3%g%l%c%p%pdgvs%t_mo
    t_mo_min         => clm3%g%l%c%p%pdgvs%t_mo_min
    prec365          => clm3%g%l%c%p%pdgvs%prec365
    agddtw           => clm3%g%l%c%p%pdgvs%agddtw
    agdd             => clm3%g%l%c%p%pdgvs%agdd
    twmax            => dgv_pftcon%twmax
#endif
    prec60           => clm3%g%l%c%p%pps%prec60
    prec10           => clm3%g%l%c%p%pps%prec10
    gdd0             => clm3%g%l%c%p%pps%gdd0
    gdd8             => clm3%g%l%c%p%pps%gdd8
    gdd10            => clm3%g%l%c%p%pps%gdd10
    gddplant         => clm3%g%l%c%p%pps%gddplant
    gddtsoi          => clm3%g%l%c%p%pps%gddtsoi
    vf               => clm3%g%l%c%p%pps%vf
    t_soisno         => clm3%g%l%c%ces%t_soisno
    h2osoi_liq       => clm3%g%l%c%cws%h2osoi_liq
    watsat           => clm3%g%l%c%cps%watsat
    dz               => clm3%g%l%c%cps%dz
    latdeg           => clm3%g%latdeg
    croplive         => clm3%g%l%c%p%pps%croplive
    pcolumn          => clm3%g%l%c%p%column
    t_veg24          => clm3%g%l%c%p%pvs%t_veg24           ! (heald 04/06)
    t_veg240         => clm3%g%l%c%p%pvs%t_veg240          ! (heald 04/06)
    fsd24            => clm3%g%l%c%p%pvs%fsd24             ! (heald 04/06)
    fsd240           => clm3%g%l%c%p%pvs%fsd240            ! (heald 04/06)
    fsi24            => clm3%g%l%c%p%pvs%fsi24             ! (heald 04/06)
    fsi240           => clm3%g%l%c%p%pvs%fsi240            ! (heald 04/06)
    fsun24           => clm3%g%l%c%p%pvs%fsun24            ! (heald 04/06)
    fsun240          => clm3%g%l%c%p%pvs%fsun240           ! (heald 04/06)
    elai_p           => clm3%g%l%c%p%pvs%elai_p            ! (heald 04/06)
    t_veg            => clm3%g%l%c%p%pes%t_veg 	           ! (heald 04/06)
    fsun             => clm3%g%l%c%p%pps%fsun 	           ! (heald 04/06)
    elai             => clm3%g%l%c%p%pps%elai 	           ! (heald 04/06)

    ! Determine calendar information

    dtime = get_step_size()
    nstep = get_nstep()
    call get_curr_date (year, month, day, secs)

    ! Don't do any accumulation if nstep is zero
    ! (only applies to coupled or cam mode)

    if (nstep == 0) return

    ! NOTE: currently only single level pft fields are used below
    ! Variables are declared above that should make it easy to incorporate
    ! multi-level or single-level fields of any subgrid type

    ! Allocate needed dynamic memory for single level pft field

    allocate(rbufslp(begp:endp), stat=ier)
    if (ier/=0) then
       write(iulog,*)'update_accum_hist allocation error for rbuf1dp'
       call endrun
    endif

    ! Accumulate and extract TREFAV - hourly average 2m air temperature
    ! Used to compute maximum and minimum of hourly averaged 2m reference
    ! temperature over a day. Note that "spval" is returned by the call to
    ! accext if the time step does not correspond to the end of an
    ! accumulation interval. First, initialize the necessary values for
    ! an initial run at the first time step the accumulator is called

    call update_accum_field  ('TREFAV', t_ref2m, nstep)
    call extract_accum_field ('TREFAV', rbufslp, nstep)
    end_cd = is_end_curr_day()
    do p = begp,endp
       if (rbufslp(p) /= spval) then
          t_ref2m_max_inst(p) = max(rbufslp(p), t_ref2m_max_inst(p))
          t_ref2m_min_inst(p) = min(rbufslp(p), t_ref2m_min_inst(p))
       endif
       if (end_cd) then
          t_ref2m_max(p) = t_ref2m_max_inst(p)
          t_ref2m_min(p) = t_ref2m_min_inst(p)
          t_ref2m_max_inst(p) = -spval
          t_ref2m_min_inst(p) =  spval
       else if (secs == int(dtime)) then
          t_ref2m_max(p) = spval
          t_ref2m_min(p) = spval
       endif
    end do

    ! Accumulate and extract TREFAV_U - hourly average urban 2m air temperature
    ! Used to compute maximum and minimum of hourly averaged 2m reference
    ! temperature over a day. Note that "spval" is returned by the call to
    ! accext if the time step does not correspond to the end of an
    ! accumulation interval. First, initialize the necessary values for
    ! an initial run at the first time step the accumulator is called

    call update_accum_field  ('TREFAV_U', t_ref2m_u, nstep)
    call extract_accum_field ('TREFAV_U', rbufslp, nstep)
    do p = begp,endp
       l = plandunit(p)
       if (rbufslp(p) /= spval) then
          t_ref2m_max_inst_u(p) = max(rbufslp(p), t_ref2m_max_inst_u(p))
          t_ref2m_min_inst_u(p) = min(rbufslp(p), t_ref2m_min_inst_u(p))
       endif
       if (end_cd) then
         if (urbpoi(l)) then
          t_ref2m_max_u(p) = t_ref2m_max_inst_u(p)
          t_ref2m_min_u(p) = t_ref2m_min_inst_u(p)
          t_ref2m_max_inst_u(p) = -spval
          t_ref2m_min_inst_u(p) =  spval
         end if
       else if (secs == int(dtime)) then
          t_ref2m_max_u(p) = spval
          t_ref2m_min_u(p) = spval
       endif
    end do

    ! Accumulate and extract TREFAV_R - hourly average rural 2m air temperature
    ! Used to compute maximum and minimum of hourly averaged 2m reference
    ! temperature over a day. Note that "spval" is returned by the call to
    ! accext if the time step does not correspond to the end of an
    ! accumulation interval. First, initialize the necessary values for
    ! an initial run at the first time step the accumulator is called

    call update_accum_field  ('TREFAV_R', t_ref2m_r, nstep)
    call extract_accum_field ('TREFAV_R', rbufslp, nstep)
    do p = begp,endp
       l = plandunit(p)
       if (rbufslp(p) /= spval) then
          t_ref2m_max_inst_r(p) = max(rbufslp(p), t_ref2m_max_inst_r(p))
          t_ref2m_min_inst_r(p) = min(rbufslp(p), t_ref2m_min_inst_r(p))
       endif
       if (end_cd) then
         if (.not.(ifspecial(l))) then
          t_ref2m_max_r(p) = t_ref2m_max_inst_r(p)
          t_ref2m_min_r(p) = t_ref2m_min_inst_r(p)
          t_ref2m_max_inst_r(p) = -spval
          t_ref2m_min_inst_r(p) =  spval
         end if
       else if (secs == int(dtime)) then
          t_ref2m_max_r(p) = spval
          t_ref2m_min_r(p) = spval
       endif
    end do

    ! Accumulate and extract T_VEG24 & T_VEG240 (heald 04/06)
    do p = begp,endp
       rbufslp(p) = t_veg(p)
    end do
    call update_accum_field  ('T_VEG24', rbufslp, nstep)
    call extract_accum_field ('T_VEG24', t_veg24, nstep)
    call update_accum_field  ('T_VEG240', rbufslp, nstep)
    call extract_accum_field ('T_VEG240', t_veg240, nstep)

    ! Accumulate and extract forc_solad24 & forc_solad240 (heald 04/06)
    do p = begp,endp
       g = pgridcell(p)
       rbufslp(p) = forc_solad(g,1)
    end do
    call update_accum_field  ('FSD240', rbufslp, nstep)
    call extract_accum_field ('FSD240', fsd240, nstep)
    call update_accum_field  ('FSD24', rbufslp, nstep)
    call extract_accum_field ('FSD24', fsd24, nstep)

    ! Accumulate and extract forc_solai24 & forc_solai240 (heald 04/06)
    do p = begp,endp
       g = pgridcell(p)
       rbufslp(p) = forc_solai(g,1)
    end do
    call update_accum_field  ('FSI24', rbufslp, nstep)
    call extract_accum_field ('FSI24', fsi24, nstep)
    call update_accum_field  ('FSI240', rbufslp, nstep)
    call extract_accum_field ('FSI240', fsi240, nstep)

    ! Accumulate and extract fsun24 & fsun240 (heald 04/06)
    do p = begp,endp
       rbufslp(p) = fsun(p)
    end do
    call update_accum_field  ('FSUN24', rbufslp, nstep)
    call extract_accum_field ('FSUN24', fsun24, nstep)
    call update_accum_field  ('FSUN240', rbufslp, nstep)
    call extract_accum_field ('FSUN240', fsun240, nstep)

    ! Accumulate and extract elai_p (heald 04/06)
    do p = begp,endp
       rbufslp(p) = elai(p)
    end do
    call update_accum_field  ('LAIP', rbufslp, nstep)
    call extract_accum_field ('LAIP', elai_p, nstep)

    ! Accumulate and extract T10
    !(acumulates TSA as 10-day running mean)

    call update_accum_field  ('T10', t_ref2m, nstep)
    call extract_accum_field ('T10', t10, nstep)

#if (defined CNDV)
    ! Accumulate and extract TDA
    ! (accumulates TBOT as 30-day average)
    ! Also determine t_mo_min

    do p = begp,endp
       g = pgridcell(p)
       rbufslp(p) = forc_t(g)
    end do
    call update_accum_field  ('TDA', rbufslp, nstep)
    call extract_accum_field ('TDA', rbufslp, nstep)
    do p = begp,endp
       t_mo(p) = rbufslp(p)
       t_mo_min(p) = min(t_mo_min(p), rbufslp(p))
    end do

    ! Accumulate and extract PREC365
    ! (accumulates total precipitation as 365-day running mean)

    do p = begp,endp
       g = pgridcell(p)
       rbufslp(p) = forc_rain(g) + forc_snow(g)
    end do
    call update_accum_field  ('PREC365', rbufslp, nstep)
    call extract_accum_field ('PREC365', prec365, nstep)

    ! Accumulate growing degree days based on 10-day running mean temperature.
    ! The trigger to reset the accumulated values to zero is -99999.

    ! Accumulate and extract AGDDTW (gdd base twmax, which is 23 deg C
    ! for boreal woody pfts)

    do p = begp,endp
       rbufslp(p) = max(0._r8, (t10(p) - SHR_CONST_TKFRZ - twmax(ndllf_dcd_brl_tree)) &
                    * dtime/SHR_CONST_CDAY)
       if (month==1 .and. day==1 .and. secs==int(dtime)) rbufslp(p) = -99999._r8
    end do
    call update_accum_field  ('AGDDTW', rbufslp, nstep)
    call extract_accum_field ('AGDDTW', agddtw, nstep)

    ! Accumulate and extract AGDD

    do p = begp,endp
       rbufslp(p) = max(0.0_r8, (t_ref2m(p) - (SHR_CONST_TKFRZ + 5.0_r8)) &
            * dtime/SHR_CONST_CDAY)
    end do
    call update_accum_field  ('AGDD', rbufslp, nstep)
    call extract_accum_field ('AGDD', agdd, nstep)
#endif
    
     do p = begp,endp
       g = pgridcell(p)
       rbufslp(p) = forc_rain(g) + forc_snow(g)
    end do
    call update_accum_field  ('PREC60', rbufslp, nstep)
    call extract_accum_field ('PREC60', prec60, nstep)

    ! Accumulate and extract PREC10
    ! (accumulates total precipitation as 10-day running mean)
     do p = begp,endp
       g = pgridcell(p)
       rbufslp(p) = forc_rain(g) + forc_snow(g)
    end do
    call update_accum_field  ('PREC10', rbufslp, nstep)
    call extract_accum_field ('PREC10', prec10, nstep)
 
    if ( crop_prog )then
       ! Accumulate and extract TDM10

       do p = begp,endp
          rbufslp(p) = min(t_ref2m_min(p),t_ref2m_min_inst(p)) !slevis: ok choice?
          if (rbufslp(p) > 1.e30_r8) rbufslp(p) = SHR_CONST_TKFRZ !and were 'min'&
       end do                                         !'min_inst' not initialized?
       call update_accum_field  ('TDM10', rbufslp, nstep)
       call extract_accum_field ('TDM10', a10tmin, nstep)

       ! Accumulate and extract TDM5

       do p = begp,endp
          rbufslp(p) = min(t_ref2m_min(p),t_ref2m_min_inst(p)) !slevis: ok choice?
          if (rbufslp(p) > 1.e30_r8) rbufslp(p) = SHR_CONST_TKFRZ !and were 'min'&
       end do                                         !'min_inst' not initialized?
       call update_accum_field  ('TDM5', rbufslp, nstep)
       call extract_accum_field ('TDM5', a5tmin, nstep)

       ! Accumulate and extract GDD0

       do p = begp,endp
          itypveg = itype(p)
          g = pgridcell(p)
          if (month==1 .and. day==1 .and. secs==int(dtime)) then
             rbufslp(p) = -99999._r8 ! reset gdd
          else if (( month > 3 .and. month < 10 .and. latdeg(g) >= 0._r8) .or. &
                   ((month > 9 .or.  month < 4) .and. latdeg(g) <  0._r8)     ) then
             rbufslp(p) = max(0._r8, min(26._r8, t_ref2m(p)-SHR_CONST_TKFRZ)) &
                          * dtime/SHR_CONST_CDAY
          else
             rbufslp(p) = 0._r8      ! keeps gdd unchanged at other times (eg, through Dec in NH)
          end if
       end do
       call update_accum_field  ('GDD0', rbufslp, nstep)
       call extract_accum_field ('GDD0', gdd0, nstep)

       ! Accumulate and extract GDD8

       do p = begp,endp
          itypveg = itype(p)
          g = pgridcell(p)
          if (month==1 .and. day==1 .and. secs==int(dtime)) then
             rbufslp(p) = -99999._r8 ! reset gdd
          else if (( month > 3 .and. month < 10 .and. latdeg(g) >= 0._r8) .or. &
                   ((month > 9 .or.  month < 4) .and. latdeg(g) <  0._r8)     ) then
             rbufslp(p) = max(0._r8, min(30._r8, &
                                         t_ref2m(p)-(SHR_CONST_TKFRZ + 8._r8))) &
                          * dtime/SHR_CONST_CDAY
          else
             rbufslp(p) = 0._r8      ! keeps gdd unchanged at other times (eg, through Dec in NH)
          end if
       end do
       call update_accum_field  ('GDD8', rbufslp, nstep)
       call extract_accum_field ('GDD8', gdd8, nstep)

       ! Accumulate and extract GDD10

       do p = begp,endp
          itypveg = itype(p)
          g = pgridcell(p)
          if (month==1 .and. day==1 .and. secs==int(dtime)) then
             rbufslp(p) = -99999._r8 ! reset gdd
          else if (( month > 3 .and. month < 10 .and. latdeg(g) >= 0._r8) .or. &
                   ((month > 9 .or.  month < 4) .and. latdeg(g) <  0._r8)     ) then
             rbufslp(p) = max(0._r8, min(30._r8, &
                                         t_ref2m(p)-(SHR_CONST_TKFRZ + 10._r8))) &
                          * dtime/SHR_CONST_CDAY
          else
             rbufslp(p) = 0._r8      ! keeps gdd unchanged at other times (eg, through Dec in NH)
          end if
       end do
       call update_accum_field  ('GDD10', rbufslp, nstep)
       call extract_accum_field ('GDD10', gdd10, nstep)

       ! Accumulate and extract GDDPLANT

       do p = begp,endp
          if (croplive(p)) then ! relative to planting date
             itypveg = itype(p)
             rbufslp(p) = max(0._r8, min(mxtmp(itypveg), &
                                         t_ref2m(p)-(SHR_CONST_TKFRZ + baset(itypveg)))) &
                          * dtime/SHR_CONST_CDAY
             if (itypveg == nwcereal .or. itypveg == nwcerealirrig) rbufslp(p) = rbufslp(p)*vf(p)
          else
             rbufslp(p) = -99999._r8
          end if
       end do
       call update_accum_field  ('GDDPLANT', rbufslp, nstep)
       call extract_accum_field ('GDDPLANT', gddplant, nstep)

       ! Accumulate and extract GDDTSOI
       ! In agroibis this variable is calculated
       ! to 0.05 m, so here we use the top two soil layers
   
       do p = begp,endp
          if (croplive(p)) then ! relative to planting date
             itypveg = itype(p)
             c = pcolumn(p)
             rbufslp(p) = max(0._r8, min(mxtmp(itypveg), &
              ((t_soisno(c,1)*dz(c,1)+t_soisno(c,2)*dz(c,2))/(dz(c,1)+dz(c,2))) - &
              (SHR_CONST_TKFRZ + baset(itypveg)))) * dtime/SHR_CONST_CDAY
             if (itypveg == nwcereal .or. itypveg == nwcerealirrig) rbufslp(p) = rbufslp(p)*vf(p)
          else
             rbufslp(p) = -99999._r8
          end if
       end do
       call update_accum_field  ('GDDTSOI', rbufslp, nstep)
       call extract_accum_field ('GDDTSOI', gddtsoi, nstep)

    end if

    ! Deallocate dynamic memory

    deallocate(rbufslp)

  end subroutine updateAccFlds

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: initAccClmtype
!
! !INTERFACE:
  subroutine initAccClmtype
!
! !DESCRIPTION:
! Initialize clmtype variables that are associated with
! time accumulated fields. This routine is called in an initial run
! at nstep=0 for cam and csm mode.
! This routine is also always called for a restart run and
! therefore must be called after the restart file is read in
! and the accumulated fields are obtained.
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clmtype
    use decompMod   , only : get_proc_bounds, get_proc_global
    use accumulMod  , only : extract_accum_field
    use clm_time_manager, only : get_nstep
    use clm_varctl  , only : nsrest, nsrStartup
    use clm_varcon  , only : spval
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
! !LOCAL VARIABLES:
!
! local pointers to implicit out arguments
!
    real(r8), pointer :: t_ref2m_min(:)      ! daily minimum of average 2 m height surface air temperature (K)
    real(r8), pointer :: t_ref2m_max(:)      ! daily maximum of average 2 m height surface air temperature (K)
    real(r8), pointer :: t_ref2m_min_inst(:) ! instantaneous daily min of average 2 m height surface air temp (K)
    real(r8), pointer :: t_ref2m_max_inst(:) ! instantaneous daily max of average 2 m height surface air temp (K)
    real(r8), pointer :: t_ref2m_min_u(:)    ! Urban daily minimum of average 2 m height surface air temperature (K)
    real(r8), pointer :: t_ref2m_min_r(:)    ! Rural daily minimum of average 2 m height surface air temperature (K)
    real(r8), pointer :: t_ref2m_max_u(:)    ! Urban daily maximum of average 2 m height surface air temperature (K)
    real(r8), pointer :: t_ref2m_max_r(:)    ! Rural daily maximum of average 2 m height surface air temperature (K)
    real(r8), pointer :: t_ref2m_min_inst_u(:) ! Urban instantaneous daily min of average 2 m height surface air temp (K)
    real(r8), pointer :: t_ref2m_min_inst_r(:) ! Rural instantaneous daily min of average 2 m height surface air temp (K)
    real(r8), pointer :: t_ref2m_max_inst_u(:) ! Urban instantaneous daily max of average 2 m height surface air temp (K)
    real(r8), pointer :: t_ref2m_max_inst_r(:) ! Rural instantaneous daily max of average 2 m height surface air temp (K)
    real(r8), pointer :: t10(:)              ! 10-day running mean of the 2 m temperature (K)
#ifdef CNDV
    real(r8), pointer :: t_mo(:)             ! 30-day average temperature (Kelvin)
    real(r8), pointer :: prec365(:)          ! 365-day running mean of tot. precipitation
    real(r8), pointer :: agddtw(:)           ! accumulated growing degree days above twmax
    real(r8), pointer :: agdd(:)             ! accumulated growing degree days above 5
#endif
    real(r8), pointer :: prec60(:)            ! 60-day running mean of tot. precipitation
    real(r8), pointer :: prec10(:)            ! 10-day running mean of tot. precipitation
    real(r8), pointer :: gdd0(:)             ! growing degree-days base 0C'
    real(r8), pointer :: gdd8(:)             ! growing degree-days base 8C from planting
    real(r8), pointer :: gdd10(:)            ! growing degree-days base 10C from planting
    real(r8), pointer :: gddplant(:)         ! growing degree-days from planting
    real(r8), pointer :: gddtsoi(:)          ! growing degree-days from planting (top two soil layers)
    real(r8), pointer :: a10tmin(:)          ! 10-day running mean of min 2-m temperature
    real(r8), pointer :: a5tmin(:)           ! 5-day running mean of min 2-m temperature
    ! heald (04/06): accumulated variables for VOC emissions
    real(r8), pointer :: t_veg24(:)          ! 24hr average vegetation temperature (K)
    real(r8), pointer :: t_veg240(:)         ! 240hr average vegetation temperature (Kelvin)
    real(r8), pointer :: fsd24(:)            ! 24hr average of direct beam radiation 
    real(r8), pointer :: fsd240(:)           ! 240hr average of direct beam radiation 
    real(r8), pointer :: fsi24(:)            ! 24hr average of diffuse beam radiation 
    real(r8), pointer :: fsi240(:)           ! 240hr average of diffuse beam radiation 
    real(r8), pointer :: fsun24(:)           ! 24hr average of sunlit fraction of canopy 
    real(r8), pointer :: fsun240(:)          ! 240hr average of sunlit fraction of canopy
    real(r8), pointer :: elai_p(:)           ! leaf area index average over timestep 
!
! !LOCAL VARIABLES:
!
!
! !OTHER LOCAL VARIABLES:
!EOP
    integer :: p            ! indices
    integer :: nstep        ! time step
    integer :: ier          ! error status
    integer :: begp, endp   ! per-proc beginning and ending pft indices
    integer :: begc, endc   ! per-proc beginning and ending column indices
    integer :: begl, endl   ! per-proc beginning and ending landunit indices
    integer :: begg, endg   ! per-proc gridcell ending gridcell indices
    real(r8), pointer :: rbufslp(:)  ! temporary
    character(len=32) :: subname = 'initAccClmtype'  ! subroutine name
!-----------------------------------------------------------------------

    ! Assign local pointers to derived subtypes components (pft-level)

    t_ref2m_max_inst => clm3%g%l%c%p%pes%t_ref2m_max_inst
    t_ref2m_min_inst => clm3%g%l%c%p%pes%t_ref2m_min_inst
    t_ref2m_max      => clm3%g%l%c%p%pes%t_ref2m_max
    t_ref2m_min      => clm3%g%l%c%p%pes%t_ref2m_min
    t_ref2m_max_inst_u => clm3%g%l%c%p%pes%t_ref2m_max_inst_u
    t_ref2m_max_inst_r => clm3%g%l%c%p%pes%t_ref2m_max_inst_r
    t_ref2m_min_inst_u => clm3%g%l%c%p%pes%t_ref2m_min_inst_u
    t_ref2m_min_inst_r => clm3%g%l%c%p%pes%t_ref2m_min_inst_r
    t_ref2m_max_u      => clm3%g%l%c%p%pes%t_ref2m_max_u
    t_ref2m_max_r      => clm3%g%l%c%p%pes%t_ref2m_max_r
    t_ref2m_min_u      => clm3%g%l%c%p%pes%t_ref2m_min_u
    t_ref2m_min_r      => clm3%g%l%c%p%pes%t_ref2m_min_r
    t10              => clm3%g%l%c%p%pes%t10
    a10tmin          => clm3%g%l%c%p%pes%a10tmin
    a5tmin           => clm3%g%l%c%p%pes%a5tmin
#if (defined CNDV)
    t_mo             => clm3%g%l%c%p%pdgvs%t_mo
    prec365          => clm3%g%l%c%p%pdgvs%prec365
    agddtw           => clm3%g%l%c%p%pdgvs%agddtw
    agdd             => clm3%g%l%c%p%pdgvs%agdd
#endif
    prec60           => clm3%g%l%c%p%pps%prec60
    prec10           => clm3%g%l%c%p%pps%prec10
    gdd0             => clm3%g%l%c%p%pps%gdd0
    gdd8             => clm3%g%l%c%p%pps%gdd8
    gdd10            => clm3%g%l%c%p%pps%gdd10
    gddplant         => clm3%g%l%c%p%pps%gddplant
    gddtsoi          => clm3%g%l%c%p%pps%gddtsoi
    ! heald (04/06): accumulated variables for VOC emissions
    t_veg24          => clm3%g%l%c%p%pvs%t_veg24
    t_veg240         => clm3%g%l%c%p%pvs%t_veg240
    fsd24            => clm3%g%l%c%p%pvs%fsd24
    fsd240           => clm3%g%l%c%p%pvs%fsd240
    fsi24            => clm3%g%l%c%p%pvs%fsi24
    fsi240           => clm3%g%l%c%p%pvs%fsi240
    fsun24           => clm3%g%l%c%p%pvs%fsun24
    fsun240          => clm3%g%l%c%p%pvs%fsun240
    elai_p           => clm3%g%l%c%p%pvs%elai_p

    ! Determine necessary indices

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    ! Determine time step

    nstep = get_nstep()

    ! Initialize 2m ref temperature max and min values

    if (nsrest == nsrStartup) then ! Why not restart&branch? These vars are not in clmr.
       do p = begp,endp
          t_ref2m_max(p) = spval
          t_ref2m_min(p) = spval
          t_ref2m_max_inst(p) = -spval
          t_ref2m_min_inst(p) =  spval
          t_ref2m_max_u(p) = spval
          t_ref2m_max_r(p) = spval
          t_ref2m_min_u(p) = spval
          t_ref2m_min_r(p) = spval
          t_ref2m_max_inst_u(p) = -spval
          t_ref2m_max_inst_r(p) = -spval
          t_ref2m_min_inst_u(p) =  spval
          t_ref2m_min_inst_r(p) =  spval
       end do
    end if

    ! Allocate needed dynamic memory for single level pft field

    allocate(rbufslp(begp:endp), stat=ier)
    if (ier/=0) then
       write(iulog,*)'extract_accum_hist allocation error for rbufslp in '//subname
       call endrun
    endif

    ! Initialize clmtype variables that are to be time accumulated

    call extract_accum_field ('T_VEG24', rbufslp, nstep)
    do p = begp,endp
       t_veg24(p) = rbufslp(p)
    end do

    call extract_accum_field ('T_VEG240', rbufslp, nstep)
    do p = begp,endp
       t_veg240(p) = rbufslp(p)
    end do

    call extract_accum_field ('FSD24', rbufslp, nstep)
    do p = begp,endp
       fsd24(p) = rbufslp(p)
    end do

    call extract_accum_field ('FSD240', rbufslp, nstep)
    do p = begp,endp
       fsd240(p) = rbufslp(p)
    end do

    call extract_accum_field ('FSI24', rbufslp, nstep)
    do p = begp,endp
       fsi24(p) = rbufslp(p)
    end do

    call extract_accum_field ('FSI240', rbufslp, nstep)
    do p = begp,endp
       fsi240(p) = rbufslp(p)
    end do

    call extract_accum_field ('FSUN24', rbufslp, nstep)
    do p = begp,endp
       fsun24(p) = rbufslp(p)
    end do

    call extract_accum_field ('FSUN240', rbufslp, nstep)
    do p = begp,endp
       fsun240(p) = rbufslp(p)
    end do

    call extract_accum_field ('LAIP', rbufslp, nstep)
    do p = begp,endp
       elai_p(p) = rbufslp(p)
    end do

    if ( crop_prog )then

       call extract_accum_field ('GDD0', rbufslp, nstep)
       do p = begp,endp
          gdd0(p) = rbufslp(p)
       end do
   
       call extract_accum_field ('GDD8', rbufslp, nstep)
       do p = begp,endp
          gdd8(p) = rbufslp(p)
       end do

       call extract_accum_field ('GDD10', rbufslp, nstep)
       do p = begp,endp
          gdd10(p) = rbufslp(p)
       end do

       call extract_accum_field ('GDDPLANT', rbufslp, nstep)
       do p = begp,endp
          gddplant(p) = rbufslp(p)
       end do

       call extract_accum_field ('GDDTSOI', rbufslp, nstep)
       do p = begp,endp
          gddtsoi(p) = rbufslp(p)
       end do

       call extract_accum_field ('TDM10', rbufslp, nstep)
       do p = begp,endp
          a10tmin(p) = rbufslp(p)
       end do

       call extract_accum_field ('TDM5', rbufslp, nstep)
       do p = begp,endp
          a5tmin(p) = rbufslp(p)
       end do

    end if

    call extract_accum_field ('T10', rbufslp, nstep)
    do p = begp,endp
       t10(p) = rbufslp(p)
    end do

#if (defined CNDV)

    call extract_accum_field ('TDA', rbufslp, nstep)
    do p = begp,endp
       t_mo(p) = rbufslp(p)
    end do

    call extract_accum_field ('PREC365', rbufslp, nstep)
    do p = begp,endp
       prec365(p) = rbufslp(p)
    end do

    call extract_accum_field ('AGDDTW', rbufslp, nstep)
    do p = begp,endp
       agddtw(p) = rbufslp(p)
    end do

    call extract_accum_field ('AGDD', rbufslp, nstep)
    do p = begp,endp
       agdd(p) = rbufslp(p)
    end do

#endif
    call extract_accum_field ('PREC60', rbufslp, nstep)
    do p = begp,endp
       prec60(p) = rbufslp(p)
    end do

    call extract_accum_field ('PREC10', rbufslp, nstep)
    do p = begp,endp
       prec10(p) = rbufslp(p)
    end do

    deallocate(rbufslp)

  end subroutine initAccClmtype

end module accFldsMod
