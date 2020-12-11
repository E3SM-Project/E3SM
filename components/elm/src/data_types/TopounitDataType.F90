module TopounitDataType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Topounit data type allocation and initialization
  ! -------------------------------------------------------- 
  !
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use abortutils     , only : endrun
  use elm_varcon     , only : spval, ispval
  use elm_varctl     , only : iulog, use_cn, use_fates, use_lch4
  use elm_varpar     , only : numrad
  use histFileMod    , only : hist_addfld1d, hist_addfld2d
  use ncdio_pio      , only : file_desc_t, ncd_double
  use decompMod      , only : bounds_type
  use restUtilMod
  use TopounitType   , only : top_pp
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private

  !-----------------------------------------------------------------------
  ! Define the data structure where land model receives atmospheric state information.
  type, public :: topounit_atmospheric_state
    real(r8), pointer :: tbot       (:) => null() ! temperature of air at atmospheric forcing height (K)
    real(r8), pointer :: thbot      (:) => null() ! potential temperature of air at atmospheric forcing height (K)
    real(r8), pointer :: pbot       (:) => null() ! air pressure at atmospheric forcing height (Pa)
    real(r8), pointer :: rhobot     (:) => null() ! air density at atmospheric forcing height (kg/m**3)
    real(r8), pointer :: qbot       (:) => null() ! specific humidity at atmospheric forcing height (kg H2O/kg moist air)
    real(r8), pointer :: rhbot      (:) => null() ! relative humidity at atmospheric forcing height (%)
    real(r8), pointer :: ubot       (:) => null() ! wind speed in U (east) direction at atmospheric forcing height (m/s)
    real(r8), pointer :: vbot       (:) => null() ! wind speed in V (north) direction at atmospheric forcing height (m/s)
    real(r8), pointer :: windbot    (:) => null() ! horizontal component of wind at atmospheric forcing height (m/s)
    real(r8), pointer :: zbot       (:) => null() ! atmospheric forcing height (m)
    real(r8), pointer :: po2bot     (:) => null() ! partial pressure of O2 at atmospheric forcing height (Pa)
    real(r8), pointer :: pco2bot    (:) => null() ! partial pressure of CO2 at atmospheric forcing height (Pa) 
    real(r8), pointer :: pc13o2bot  (:) => null() ! partial pressure of C13O2 at atmospheric forcing height (Pa) 
    real(r8), pointer :: pch4bot    (:) => null() ! partial pressure of CH4 at atmospheric forcing height (Pa)
    ! Accumulated fields
    real(r8), pointer :: rh24h      (:) => null() ! 24-hour running mean of relative humidity at atmospheric forcing height (%)
    real(r8), pointer :: wind24h    (:) => null() ! 24-hour running mean of horizontal wind at atmospheric forcing height (m/s)
  contains
    procedure, public :: Init  => init_top_as
    procedure, public :: Clean => clean_top_as
    procedure, public :: InitAccBuffer => init_acc_buffer_top_as
    procedure, public :: InitAccVars   => init_acc_vars_top_as
    procedure, public :: UpdateAccVars => update_acc_vars_top_as
  end type topounit_atmospheric_state

  !-----------------------------------------------------------------------
  ! Define the data structure that where land model receives atmospheric flux information.
  type, public :: topounit_atmospheric_flux
    real(r8), pointer :: rain      (:)   => null() ! rain rate (kg H2O/m**2/s, equivalent to mm liquid H2O/s)
    real(r8), pointer :: snow      (:)   => null() ! snow rate (kg H2O/m**2/s, equivalent to mm liquid H2O/s)
    real(r8), pointer :: solad     (:,:) => null() ! direct beam radiation (numrad) (vis=forc_sols , nir=forc_soll) (W/m**2)
    real(r8), pointer :: solai     (:,:) => null() ! diffuse radiation (numrad) (vis=forc_solsd, nir=forc_solld) (W/m**2)
    real(r8), pointer :: solar     (:)   => null() ! incident solar radiation (W/m**2)
    real(r8), pointer :: lwrad     (:)   => null() ! atm downwrd IR longwave radiation (W/m**2) 
    ! Accumulated fields           
    real(r8), pointer :: prec24h   (:)   => null() ! 24-hour mean precip rate (kg H2O/m**2/s, equivalent to mm liquid H2O/s)
    real(r8), pointer :: prec10d   (:)   => null() ! 10-day mean precip rate (kg H2O/m**2/s, equivalent to mm liquid H2O/s)
    real(r8), pointer :: prec60d   (:)   => null() ! 60-day mean precip rate (kg H2O/m**2/s, equivalent to mm liquid H2O/s)
    real(r8), pointer :: fsd24h    (:)   => null() ! 24hr average of direct beam radiation (W/m**2)
    real(r8), pointer :: fsd240h   (:)   => null() ! 240hr average of direct beam radiation (W/m**2) 
    real(r8), pointer :: fsi24h    (:)   => null() ! 24hr average of diffuse beam radiation (W/m**2) 
    real(r8), pointer :: fsi240h   (:)   => null() ! 240hr average of diffuse beam radiation (W/m**2) 
    
  contains
    procedure, public :: Init  => init_top_af
    procedure, public :: Clean => clean_top_af
    procedure, public :: InitAccBuffer => init_acc_buffer_top_af
    procedure, public :: InitAccVars   => init_acc_vars_top_af
    procedure, public :: UpdateAccVars => update_acc_vars_top_af
  end type topounit_atmospheric_flux
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds energy state information for land at the level of topographic unit.
  type, public :: topounit_energy_state
    real(r8), pointer :: t_rad      (:) => null() ! mean radiative temperature of land surface (K)
  contains
    procedure, public :: Init  => init_top_es
    procedure, public :: Clean => clean_top_es
  end type topounit_energy_state

  !-----------------------------------------------------------------------
  ! declare the public instances of topounit data types
  type(topounit_atmospheric_state),    public, target :: top_as
  type(topounit_atmospheric_flux),     public, target :: top_af
  type(topounit_energy_state),         public, target :: top_es

  contains

  !-----------------------------------------------------------------------
  subroutine init_top_as(this, begt, endt)
    class(topounit_atmospheric_state) :: this
    integer, intent(in) :: begt   ! beginning topographic unit index
    integer, intent(in) :: endt   ! ending topographic unit index
   
    ! Allocate for atmospheric state forcing variables, initialize to special value
    allocate(this%tbot     (begt:endt)) ; this%tbot      (:) = spval
    allocate(this%thbot    (begt:endt)) ; this%thbot     (:) = spval
    allocate(this%pbot     (begt:endt)) ; this%pbot      (:) = spval
    allocate(this%rhobot   (begt:endt)) ; this%rhobot    (:) = spval
    allocate(this%qbot     (begt:endt)) ; this%qbot      (:) = spval
    allocate(this%rhbot    (begt:endt)) ; this%rhbot     (:) = spval
    allocate(this%ubot     (begt:endt)) ; this%ubot      (:) = spval
    allocate(this%vbot     (begt:endt)) ; this%vbot      (:) = spval
    allocate(this%windbot  (begt:endt)) ; this%windbot   (:) = spval
    allocate(this%zbot     (begt:endt)) ; this%zbot      (:) = spval
    allocate(this%po2bot   (begt:endt)) ; this%po2bot    (:) = spval
    allocate(this%pco2bot  (begt:endt)) ; this%pco2bot   (:) = spval
    allocate(this%pc13o2bot(begt:endt)) ; this%pc13o2bot (:) = spval
    allocate(this%pch4bot  (begt:endt)) ; this%pch4bot   (:) = spval
    if (use_fates) then
      allocate(this%rh24h  (begt:endt)) ; this%rh24h     (:) = spval
      allocate(this%wind24h(begt:endt)) ; this%wind24h   (:) = spval
    end if
    
    !-----------------------------------------------------------------------
    ! initialize history fields for select members of top_as
    !-----------------------------------------------------------------------
    this%tbot(begt:endt) = spval
    call hist_addfld1d (fname='TBOT', units='K',  &
         avgflag='A', long_name='atmospheric air temperature', &
         ptr_lnd=this%tbot)

    this%thbot(begt:endt) = spval
    call hist_addfld1d (fname='THBOT', units='K',  &
         avgflag='A', long_name='atmospheric air potential temperature', &
         ptr_lnd=this%thbot)

    this%pbot(begt:endt) = spval
    call hist_addfld1d (fname='PBOT', units='Pa',  &
         avgflag='A', long_name='atmospheric pressure', &
         ptr_lnd=this%pbot)

    this%qbot(begt:endt) = spval
    call hist_addfld1d (fname='QBOT', units='kg/kg',  &
         avgflag='A', long_name='atmospheric specific humidity', &
         ptr_lnd=this%qbot)

    this%rhbot(begt:endt) = spval
    call hist_addfld1d (fname='RH', units='%',  &
         avgflag='A', long_name='atmospheric relative humidity', &
         ptr_gcell=this%rhbot, default='inactive')

    this%windbot(begt:endt) = spval
    call hist_addfld1d (fname='WIND', units='m/s',  &
         avgflag='A', long_name='atmospheric wind velocity magnitude', &
         ptr_lnd=this%windbot)

    this%zbot(begt:endt) = spval
    call hist_addfld1d (fname='ZBOT', units='m',  &
         avgflag='A', long_name='atmospheric reference height', &
         ptr_lnd=this%zbot)

    this%pco2bot(begt:endt) = spval
    call hist_addfld1d (fname='PCO2', units='Pa',  &
         avgflag='A', long_name='atmospheric partial pressure of CO2', &
         ptr_lnd=this%pco2bot)

    if (use_lch4) then
       this%pch4bot(begt:endt) = spval
       call hist_addfld1d (fname='PCH4', units='Pa',  &
            avgflag='A', long_name='atmospheric partial pressure of CH4', &
            ptr_lnd=this%pch4bot)
    end if
  end subroutine init_top_as

  !-----------------------------------------------------------------------
  subroutine clean_top_as(this, begt, endt)
    class(topounit_atmospheric_state) :: this
    integer, intent(in) :: begt   ! beginning topographic unit index
    integer, intent(in) :: endt   ! ending topographic unit index

    deallocate(this%tbot)
    deallocate(this%thbot)
    deallocate(this%pbot)
    deallocate(this%rhobot)
    deallocate(this%qbot)
    deallocate(this%rhbot)
    deallocate(this%ubot)
    deallocate(this%vbot)
    deallocate(this%windbot)
    deallocate(this%zbot)
    deallocate(this%po2bot)
    deallocate(this%pco2bot)
    deallocate(this%pc13o2bot)
    deallocate(this%pch4bot)
    if (use_fates) then
      deallocate(this%rh24h)
      deallocate(this%wind24h)
    end if
  end subroutine clean_top_as
  
  !-----------------------------------------------------------------------
  subroutine init_acc_buffer_top_as (this, bounds)
    ! !DESCRIPTION:
    ! Initialize accumulation buffer for accumulated fields for atmospheric state
    ! This routine set defaults values that are then overwritten by the
    ! restart file for restart or branch runs
    !
    ! !USES 
    use elm_varcon  , only : spval
    use accumulMod  , only : init_accum_field
    !
    ! !ARGUMENTS:
    class(topounit_atmospheric_state) :: this
    type(bounds_type), intent(in) :: bounds  
    !---------------------------------------------------------------------
    ! Accumulator state variables used by FATES
    if (use_fates) then
       call init_accum_field (name='RH24H', units='%', &
            desc='24hr running mean of relative humidity', accum_type='runmean', accum_period=-1, &
            subgrid_type='topounit', numlev=1, init_value=0._r8)
       call init_accum_field (name='WIND24H', units='m/s', &
            desc='24hr running mean of wind', accum_type='runmean', accum_period=-1, &
            subgrid_type='topounit', numlev=1, init_value=0._r8)
    end if
  end subroutine init_acc_buffer_top_as

  !-----------------------------------------------------------------------
  subroutine init_acc_vars_top_as(this, bounds)
    ! !DESCRIPTION:
    ! Initialize variables associated with atmospheric state
    ! time accumulated fields. This routine is called for both an initial run
    ! and a restart run (and must therefore must be called after the restart file 
    ! is read in and the accumulation buffer is obtained)
    !
    ! !USES 
    use accumulMod       , only : extract_accum_field
    use clm_time_manager , only : get_nstep
    !
    ! !ARGUMENTS:
    class(topounit_atmospheric_state) :: this
    type(bounds_type), intent(in)    :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer  :: begt, endt
    integer  :: nstep
    integer  :: ier
    real(r8), pointer :: rbufslt(:)  ! temporary
    !---------------------------------------------------------------------

    begt = bounds%begt; endt = bounds%endt

    ! Allocate needed dynamic memory for single level topounit field
    allocate(rbufslt(begt:endt), stat=ier)
    if (ier/=0) then
       write(iulog,*)' in '
       call endrun(msg="extract_accum_hist allocation error for rbufslt"//&
            errMsg(__FILE__, __LINE__))
    endif

    ! Determine time step
    nstep = get_nstep()

    if (use_fates) then
       call extract_accum_field ('RH24H', rbufslt, nstep)
       this%rh24h(begt:endt) = rbufslt(begt:endt)

       call extract_accum_field ('WIND24H', rbufslt, nstep)
       this%wind24h(begt:endt) = rbufslt(begt:endt)
    end if

    deallocate(rbufslt)
  end subroutine init_acc_vars_top_as

  !-----------------------------------------------------------------------
  subroutine update_acc_vars_top_as (this, bounds)
    !
    ! USES
    use clm_time_manager, only : get_nstep
    use accumulMod      , only : update_accum_field, extract_accum_field
    !
    ! !ARGUMENTS:
    class(topounit_atmospheric_state)    :: this
    type(bounds_type)      , intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: g,t,c,p                   ! indices
    integer :: dtime                     ! timestep size [seconds]
    integer :: nstep                     ! timestep number
    integer :: ier                       ! error status
    integer :: begt, endt
    real(r8), pointer :: rbufslt(:)      ! temporary single level - topounit level
    !---------------------------------------------------------------------

    begt = bounds%begt; endt = bounds%endt

    nstep = get_nstep()

    ! Allocate needed dynamic memory for single level topounit field

    allocate(rbufslt(begt:endt), stat=ier)
    if (ier/=0) then
       write(iulog,*)'update_accum_hist allocation error for rbufslt'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    ! Accumulate 24-hour running mean of relative humidity
    do t = begt,endt
       rbufslt(t) = this%rhbot(t)
    end do
    if (use_fates) then
       call update_accum_field  ('RH24H', rbufslt, nstep)
       call extract_accum_field ('RH24H', this%rh24h, nstep)
    end if
    
    ! Accumulate 24-hour running mean of wind speed
    do t = begt,endt
       rbufslt(t) = this%windbot(t)
    end do
    if (use_fates) then
       call update_accum_field  ('WIND24H', rbufslt, nstep)
       call extract_accum_field ('WIND24H', this%wind24h, nstep)
    end if

    deallocate(rbufslt)
  end subroutine update_acc_vars_top_as

  !-----------------------------------------------------------------------
  subroutine init_top_af(this, begt, endt)
    class(topounit_atmospheric_flux) :: this
    integer, intent(in) :: begt   ! beginning topographic unit index
    integer, intent(in) :: endt   ! ending topographic unit index
   
    ! Allocate for atmospheric flux forcing variables, initialize to special value
    allocate(this%rain     (begt:endt))          ; this%rain      (:) = spval
    allocate(this%snow     (begt:endt))          ; this%snow      (:) = spval
    allocate(this%solad    (begt:endt, numrad))  ; this%solad     (:,:) = spval
    allocate(this%solai    (begt:endt, numrad))  ; this%solai     (:,:) = spval
    allocate(this%solar    (begt:endt))          ; this%solar     (:) = spval
    allocate(this%lwrad    (begt:endt))          ; this%lwrad     (:) = spval
    if (use_fates) then
      allocate(this%prec24h  (begt:endt)) ; this%prec24h   (:) = spval
    end if
    if (use_cn) then
      allocate(this%prec10d  (begt:endt)) ; this%prec10d   (:) = spval
      allocate(this%prec60d  (begt:endt)) ; this%prec60d   (:) = spval
    end if
    allocate(this%fsd24h   (begt:endt))          ; this%fsd24h    (:) = spval
    allocate(this%fsd240h  (begt:endt))          ; this%fsd240h   (:) = spval
    allocate(this%fsi24h   (begt:endt))          ; this%fsi24h    (:) = spval
    allocate(this%fsi240h  (begt:endt))          ; this%fsi240h   (:) = spval
    
    
    !-----------------------------------------------------------------------
    ! initialize history fields for select members of top_af
    !-----------------------------------------------------------------------
    this%rain(begt:endt) = spval
    call hist_addfld1d (fname='RAIN', units='mm/s',  &
         avgflag='A', long_name='atmospheric rain', &
         ptr_lnd=this%rain)

    this%snow(begt:endt) = spval
    call hist_addfld1d (fname='SNOW', units='mm/s',  &
         avgflag='A', long_name='atmospheric snow', &
         ptr_lnd=this%snow)

    this%solar(begt:endt) = spval
    call hist_addfld1d (fname='FSDS', units='W/m^2',  &
         avgflag='A', long_name='atmospheric incident solar radiation', &
         ptr_lnd=this%solar)

    this%lwrad(begt:endt) = spval
    call hist_addfld1d (fname='FLDS', units='W/m^2',  &
         avgflag='A', long_name='atmospheric longwave radiation', &
         ptr_lnd=this%lwrad)

    end subroutine init_top_af

  !-----------------------------------------------------------------------
  subroutine clean_top_af(this, begt, endt)
    class(topounit_atmospheric_flux) :: this
    integer, intent(in) :: begt   ! beginning topographic unit index
    integer, intent(in) :: endt   ! ending topographic unit index

    deallocate(this%rain)
    deallocate(this%snow)
    deallocate(this%solad)
    deallocate(this%solai)
    deallocate(this%solar)
    deallocate(this%lwrad)
    if (use_fates) then
      deallocate(this%prec24h)
    end if
    if (use_cn) then
      deallocate(this%prec10d)
      deallocate(this%prec60d)
    end if
    deallocate(this%fsd24h)
    deallocate(this%fsd240h)
    deallocate(this%fsi24h)
    deallocate(this%fsi240h)
  end subroutine clean_top_af
  
  !-----------------------------------------------------------------------
  subroutine init_acc_buffer_top_af (this, bounds)
    ! !DESCRIPTION:
    ! Initialize accumulation buffer accumulated fields for atmospheric flux
    ! This routine set defaults values that are then overwritten by the
    ! restart file for restart or branch runs
    !
    ! !USES 
    use elm_varcon  , only : spval
    use accumulMod  , only : init_accum_field
    !
    ! !ARGUMENTS:
    class(topounit_atmospheric_flux) :: this
    type(bounds_type), intent(in) :: bounds  
    !---------------------------------------------------------------------
    ! Accumulator flux variables used by CN
    if (use_cn) then
      call init_accum_field (name='PREC10D', units='MM H2O/S', &
         desc='10-day running mean of total precipitation', accum_type='runmean', accum_period=-10, &
         subgrid_type='topounit', numlev=1, init_value=0._r8)
      call init_accum_field (name='PREC60D', units='MM H2O/S', &
         desc='60-day running mean of total precipitation', accum_type='runmean', accum_period=-60, &
         subgrid_type='topounit', numlev=1, init_value=0._r8)
    end if
    ! Accumulator flux variables used by FATES
    if (use_fates) then
       call init_accum_field (name='PREC24H', units='MM H2O/S', &
            desc='24hr running mean of total precipitation', accum_type='runmean', accum_period=-1, &
            subgrid_type='topounit', numlev=1, init_value=0._r8)
    end if
    ! Accumulator variables for radiation fluxes
    !call init_accum_field (name='FSD24H', units='W/m2',                                             &
    !     desc='24hr average of direct solar radiation',  accum_type='runmean', accum_period=-1,    &
    !     subgrid_type='topounit', numlev=1, init_value=0._r8)

    !call init_accum_field (name='FSD240H', units='W/m2',                                            &
    !     desc='240hr average of direct solar radiation',  accum_type='runmean', accum_period=-10,  &
    !     subgrid_type='topounit', numlev=1, init_value=0._r8)

    !call init_accum_field (name='FSI24H', units='W/m2',                                             &
    !     desc='24hr average of diffuse solar radiation',  accum_type='runmean', accum_period=-1,   &
    !     subgrid_type='topounit', numlev=1, init_value=0._r8)

    !call init_accum_field (name='FSI240H', units='W/m2',                                            &
    !     desc='240hr average of diffuse solar radiation',  accum_type='runmean', accum_period=-10, &
    !     subgrid_type='topounit', numlev=1, init_value=0._r8)
  end subroutine init_acc_buffer_top_af
  
  !-----------------------------------------------------------------------
  subroutine init_acc_vars_top_af(this, bounds)
    ! !DESCRIPTION:
    ! Initialize variables associated with atmospheric flux
    ! time accumulated fields. This routine is called for both an initial run
    ! and a restart run (and must therefore must be called after the restart file 
    ! is read in and the accumulation buffer is obtained)
    !
    ! !USES 
    use accumulMod       , only : extract_accum_field
    use clm_time_manager , only : get_nstep
    !
    ! !ARGUMENTS:
    class(topounit_atmospheric_flux) :: this
    type(bounds_type), intent(in)    :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer  :: begt, endt
    integer  :: nstep
    integer  :: ier
    real(r8), pointer :: rbufslt(:)  ! temporary
    !---------------------------------------------------------------------

    begt = bounds%begt; endt = bounds%endt

    ! Allocate needed dynamic memory for single level topounit field
    allocate(rbufslt(begt:endt), stat=ier)
    if (ier/=0) then
       write(iulog,*)' in '
       call endrun(msg="extract_accum_hist allocation error for rbufslt"//&
            errMsg(__FILE__, __LINE__))
    endif

    ! Determine time step
    nstep = get_nstep()

    !call extract_accum_field ('FSD24H', rbufslt, nstep)
    !this%fsd24h(begt:endt) = rbufslt(begt:endt)

    !call extract_accum_field ('FSD240H', rbufslt, nstep)
    !this%fsd240h(begt:endt) = rbufslt(begt:endt)

    !call extract_accum_field ('FSI24H', rbufslt, nstep)
    !this%fsi24h(begt:endt) = rbufslt(begt:endt)

    !call extract_accum_field ('FSI240H', rbufslt, nstep)
    !this%fsi240h(begt:endt) = rbufslt(begt:endt)

    if (use_cn) then
       call extract_accum_field ('PREC10D', rbufslt, nstep)
       this%prec10d(begt:endt) = rbufslt(begt:endt)

       call extract_accum_field ('PREC60D', rbufslt, nstep)
       this%prec60d(begt:endt) = rbufslt(begt:endt)
    end if

    if (use_fates) then
       call extract_accum_field ('PREC24H', rbufslt, nstep)
       this%prec24h(begt:endt) = rbufslt(begt:endt)
    end if

    deallocate(rbufslt)
  end subroutine init_acc_vars_top_af

  !-----------------------------------------------------------------------
  subroutine update_acc_vars_top_af (this, bounds)
    !
    ! USES
    use clm_time_manager, only : get_nstep
    use accumulMod      , only : update_accum_field, extract_accum_field
    !
    ! !ARGUMENTS:
    class(topounit_atmospheric_flux)    :: this
    type(bounds_type)      , intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: g,t,c,p                   ! indices
    integer :: dtime                     ! timestep size [seconds]
    integer :: nstep                     ! timestep number
    integer :: ier                       ! error status
    integer :: begt, endt
    real(r8), pointer :: rbufslt(:)      ! temporary single level - topounit level
    !---------------------------------------------------------------------

    begt = bounds%begt; endt = bounds%endt

    nstep = get_nstep()

    ! Allocate needed dynamic memory for single level topounit field

    allocate(rbufslt(begt:endt), stat=ier)
    if (ier/=0) then
       write(iulog,*)'update_accum_hist allocation error for rbufslt'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    ! Accumulate and extract forc_solad24 & forc_solad240 
    do t = begt,endt
       rbufslt(t) = this%solad(t,1)
    end do
    !call update_accum_field  ('FSD240H', rbufslt              , nstep)
    !call extract_accum_field ('FSD240H', this%fsd240h         , nstep)
    !call update_accum_field  ('FSD24H' , rbufslt              , nstep)
    !call extract_accum_field ('FSD24H' , this%fsd24h          , nstep)

    ! Accumulate and extract forc_solai24 & forc_solai240 
    do t = begt,endt
       rbufslt(t) = this%solai(t,1)
    end do
    !call update_accum_field  ('FSI24H' , rbufslt              , nstep)
    !call extract_accum_field ('FSI24H' , this%fsi24h          , nstep)
    !call update_accum_field  ('FSI240H', rbufslt              , nstep)
    !call extract_accum_field ('FSI240H', this%fsi240h         , nstep)

    ! Accumulate and extract total precip
    do t = begt,endt
       rbufslt(t) = this%rain(t) + this%snow(t)
    end do
    if (use_cn) then
       ! Accumulate and extract PREC60D (accumulates total precipitation as 60-day running mean)
       call update_accum_field  ('PREC60D', rbufslt, nstep)
       call extract_accum_field ('PREC60D', this%prec60d, nstep)

       ! Accumulate and extract PREC10D (accumulates total precipitation as 10-day running mean)
       call update_accum_field  ('PREC10D', rbufslt, nstep)
       call extract_accum_field ('PREC10D', this%prec10d, nstep)
    end if

    if (use_fates) then
       call update_accum_field  ('PREC24H', rbufslt, nstep)
       call extract_accum_field ('PREC24H', this%prec24h, nstep)
    end if

    deallocate(rbufslt)

  end subroutine update_acc_vars_top_af

  !-----------------------------------------------------------------------
  subroutine init_top_es(this, begt, endt)
    class(topounit_energy_state) :: this
    integer, intent(in) :: begt   ! beginning topographic unit index
    integer, intent(in) :: endt   ! ending topographic unit index

    allocate(this%t_rad   (begt:endt)) ; this%t_rad   (:) = nan
  end subroutine init_top_es
  
  !-----------------------------------------------------------------------
  subroutine clean_top_es(this, begt, endt)
    class(topounit_energy_state) :: this
    integer, intent(in) :: begt   ! beginning topographic unit index
    integer, intent(in) :: endt   ! ending topographic unit index
    
    deallocate(this%t_rad    )
  end subroutine clean_top_es
  

end module TopounitDataType
