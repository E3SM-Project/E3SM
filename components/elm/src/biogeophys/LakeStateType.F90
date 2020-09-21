module LakeStateType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Lake data types and associated procesures
  !
  ! !USES:
  use shr_kind_mod , only : r8 => shr_kind_r8
  use shr_log_mod  , only : errMsg => shr_log_errMsg
  use elm_varcon   , only : spval, grlnd
  use decompMod    , only : bounds_type
  use spmdMod      , only : masterproc
  use abortUtils   , only : endrun
  use LandunitType , only : lun_pp                
  use ColumnType   , only : col_pp                
  !
  implicit none
  save
  private
  !
  ! !PUBLIC TYPES:
  type, public :: lakestate_type
     ! Time constant variables
     real(r8), pointer :: lakefetch_col     (:)   ! col lake fetch from surface data (m)                    
     real(r8), pointer :: etal_col          (:)   ! col lake extinction coefficient from surface data (1/m) 

     ! Time varying variables
     real(r8), pointer :: lake_raw_col      (:)   ! col aerodynamic resistance for moisture (s/m)
     real(r8), pointer :: ks_col            (:)   ! col coefficient for calculation of decay of eddy diffusivity with depth
     real(r8), pointer :: ws_col            (:)   ! col surface friction velocity (m/s)
     real(r8), pointer :: ust_lake_col      (:)   ! col friction velocity (m/s)          
     real(r8), pointer :: betaprime_col     (:)   ! col effective beta: sabg_lyr(p,jtop) for snow layers, beta otherwise
     real(r8), pointer :: savedtke1_col     (:)   ! col top level eddy conductivity from previous timestep (W/mK)
     real(r8), pointer :: lake_icefrac_col  (:,:) ! col mass fraction of lake layer that is frozen
     real(r8), pointer :: lake_icethick_col (:)   ! col ice thickness (m) (integrated if lakepuddling)
     real(r8), pointer :: lakeresist_col    (:)   ! col [s/m] (Needed for calc. of grnd_ch4_cond)
     real(r8), pointer :: ram1_lake_patch   (:)   ! patch aerodynamical resistance (s/m)

   contains

     procedure, public  :: Init         
     procedure, public  :: Restart      
     procedure, private :: InitAllocate 
     procedure, private :: InitHistory  
     procedure, private :: InitCold     

  end type lakestate_type
  !-----------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(lakestate_type) :: this
    type(bounds_type), intent(in) :: bounds  

    call this%InitAllocate ( bounds )
    call this%InitHistory ( bounds )
    call this%InitCold ( bounds ) 

  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Allocate module variables and data structures
    !
    ! !USES:
    use shr_infnan_mod, only: nan => shr_infnan_nan, assignment(=)
    use clm_varpar    , only: nlevlak, nlevsno
    !
    ! !ARGUMENTS:
    class(lakestate_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    !---------------------------------------------------------------------

    ! Initialize savedtke1 to spval so that c->g averaging will be done correctly
    ! TODO: can this be now be set to nan???
    ! Initialize ust_lake to spval to detect input from restart file if not arbinit 
    ! TODO: can this be removed now???

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc = bounds%endc

    allocate(this%etal_col           (begc:endc))           ; this%etal_col           (:)   = nan
    allocate(this%lakefetch_col      (begc:endc))           ; this%lakefetch_col      (:)   = nan
    allocate(this%lakeresist_col     (begc:endc))           ; this%lakeresist_col     (:)   = nan
    allocate(this%savedtke1_col      (begc:endc))           ; this%savedtke1_col      (:)   = spval  
    allocate(this%lake_icefrac_col   (begc:endc,1:nlevlak)) ; this%lake_icefrac_col   (:,:) = nan
    allocate(this%lake_icethick_col  (begc:endc))           ; this%lake_icethick_col  (:)   = nan
    allocate(this%ust_lake_col       (begc:endc))           ; this%ust_lake_col       (:)   = spval   
    allocate(this%ram1_lake_patch      (begp:endp))           ; this%ram1_lake_patch      (:)   = nan
    allocate(this%lake_raw_col       (begc:endc))           ; this%lake_raw_col       (:)   = nan
    allocate(this%ks_col             (begc:endc))           ; this%ks_col             (:)   = nan
    allocate(this%ws_col             (begc:endc))           ; this%ws_col             (:)   = nan
    allocate(this%betaprime_col      (begc:endc))           ; this%betaprime_col      (:)   = nan

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! History fields initialization
    !
    ! !USES:
    use shr_infnan_mod, only: nan => shr_infnan_nan, assignment(=)
    use histFileMod   , only: hist_addfld1d, hist_addfld2d
    !
    ! !ARGUMENTS:
    class(lakestate_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    !---------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc

    this%lake_icefrac_col(begc:endc,:) = spval
    call hist_addfld2d (fname='LAKEICEFRAC',  units='unitless', type2d='levlak', &
         avgflag='A', long_name='lake layer ice mass fraction', &
         ptr_col=this%lake_icefrac_col)
    
    this%lake_icethick_col(begc:endc) = spval ! This will be more useful than LAKEICEFRAC for many users.
    call hist_addfld1d (fname='LAKEICETHICK', units='m', &
         avgflag='A', long_name='thickness of lake ice (including physical expansion on freezing)', &
         ptr_col=this%lake_icethick_col, set_nolake=spval)

    this%savedtke1_col(begc:endc) = spval
    call hist_addfld1d (fname='TKE1',  units='W/(mK)', &
         avgflag='A', long_name='top lake level eddy thermal conductivity', &
         ptr_col=this%savedtke1_col)

    this%ram1_lake_patch(begp:endp) = spval
    call hist_addfld1d (fname='RAM_LAKE', units='s/m', &
         avgflag='A', long_name='aerodynamic resistance for momentum (lakes only)', &
         ptr_patch=this%ram1_lake_patch, set_nolake=spval, default='inactive')

    this%ust_lake_col(begc:endc) = spval
    call hist_addfld1d (fname='UST_LAKE', units='m/s', &
         avgflag='A', long_name='friction velocity (lakes only)', &
         ptr_col=this%ust_lake_col, set_nolake=spval, default='inactive')

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize time constant and time varying module variables
    !
    ! !USES:
    use clm_varctl , only : fsurdat
    use clm_varctl , only : iulog
    use clm_varpar , only : nlevlak
    use elm_varcon , only : tkwat 
    use fileutils  , only : getfil
    use ncdio_pio  , only : file_desc_t, ncd_defvar, ncd_io, ncd_double, ncd_int, ncd_inqvdlen
    use ncdio_pio  , only : ncd_pio_openfile, ncd_inqfdims, ncd_pio_closefile, ncd_inqdid, ncd_inqdlen
    !
    ! !ARGUMENTS:
    class(lakestate_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer             :: c,g,i,j,l,lev
    logical             :: readvar 
    type(file_desc_t)   :: ncid                   ! netcdf id
    character(len=256)  :: locfn                  ! local filename
    real(r8)            :: depthratio             ! ratio of lake depth to standard deep lake depth
    real(r8) ,pointer   :: lakefetch_in (:)       ! read in - lakefetch 
    real(r8) ,pointer   :: etal_in (:)            ! read in - etal 
    !-----------------------------------------------------------------------

    !-------------------------------------------------
    ! Initialize time constant variables
    !-------------------------------------------------

    call getfil (fsurdat, locfn, 0)
    call ncd_pio_openfile (ncid, locfn, 0)

    ! Read lake eta
    allocate(etal_in(bounds%begg:bounds%endg))
    call ncd_io(ncid=ncid, varname='ETALAKE', flag='read', data=etal_in, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       if (masterproc) then
          write(iulog,*) 'WARNING:: ETALAKE not found on surface data set. All lake columns will have eta', &
               ' set equal to default value as a function of depth.'
       end if
       etal_in(:) = -1._r8
    end if
    do c = bounds%begc, bounds%endc
       g = col_pp%gridcell(c)
       this%etal_col(c) = etal_in(g)
    end do
    deallocate(etal_in)

    ! Read lake fetch
    allocate(lakefetch_in(bounds%begg:bounds%endg))
    call ncd_io(ncid=ncid, varname='LAKEFETCH', flag='read', data=lakefetch_in, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       if (masterproc) then
          write(iulog,*) 'WARNING:: LAKEFETCH not found on surface data set. All lake columns will have fetch', &
               ' set equal to default value as a function of depth.'
       end if
       lakefetch_in(:) = -1._r8
    end if
    do c = bounds%begc, bounds%endc
       g = col_pp%gridcell(c)
       this%lakefetch_col(c) = lakefetch_in(g)
    end do
    deallocate(lakefetch_in)

    call ncd_pio_closefile(ncid)

    !-------------------------------------------------
    ! Initialize time varying variables
    !-------------------------------------------------
         
    do c = bounds%begc, bounds%endc
       l = col_pp%landunit(c)
       if (lun_pp%lakpoi(l)) then

          ! Set lake ice fraction and top eddy conductivity from previous timestep
          ! Always initialize with no ice to prevent excessive ice sheets from forming when
          ! starting with old lake model that has unrealistically cold lake conseratures.
          ! Keep lake temperature as is, and the energy deficit below freezing (which is no smaller
          ! than it would have been with prognostic ice, as the temperature would then have been higher
          ! and more heat would have flowed out of the lake) will be converted to ice in the first timestep.
          this%lake_icefrac_col(c,1:nlevlak) = 0._r8

          ! Set lake top eddy conductivity from previous timestep
          this%savedtke1_col(c) = tkwat

          ! Set column friction vlocity 
          this%ust_lake_col(c)  = 0.1_r8
       end if
    end do

  end subroutine InitCold

  !------------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    ! 
    ! !DESCRIPTION:
    ! Read/Write module information to/from restart file.
    !
    ! !USES:
    use ncdio_pio  , only : file_desc_t, ncd_defvar, ncd_io, ncd_double, ncd_int, ncd_inqvdlen
    use restUtilMod
    !
    ! !ARGUMENTS:
    class(lakestate_type) :: this
    type(bounds_type), intent(in)    :: bounds  
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    character(len=*) , intent(in)    :: flag   ! 'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    integer :: j,c ! indices
    logical :: readvar      ! determine if variable is on initial file
    !-----------------------------------------------------------------------

    call restartvar(ncid=ncid, flag=flag, varname='LAKE_ICEFRAC', xtype=ncd_double,  &
         dim1name='column', dim2name='levlak', switchdim=.true., &
         long_name='lake layer ice fraction', units='kg/kg', &
         interpinic_flag='interp', readvar=readvar, data=this%lake_icefrac_col)

    call restartvar(ncid=ncid, flag=flag, varname='SAVEDTKE1', xtype=ncd_double,  &
         dim1name='column', &
         long_name='top lake layer eddy conductivity', units='W/(m K)', &
         interpinic_flag='interp', readvar=readvar, data=this%savedtke1_col)

    call restartvar(ncid=ncid, flag=flag, varname='USTLAKE', xtype=ncd_double,  &
         dim1name='column', &
         long_name='friction velocity for lakes', units='m/s', &
         interpinic_flag='interp', readvar=readvar, data=this%ust_lake_col)

  end subroutine Restart

end module LakeStateType

