module LandunitDataType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Landunit data type allocation and initialization
  ! -------------------------------------------------------- 
  !
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use elm_varpar     , only : nlevsno, nlevgrnd, nlevlak, nlevurb
  use elm_varcon     , only : spval, ispval
  use elm_varctl     , only : use_vancouver, use_mexicocity
  use histFileMod    , only : hist_addfld1d
  use ncdio_pio      , only : file_desc_t, ncd_double
  use decompMod      , only : bounds_type
  use restUtilMod
  use LandunitType   , only : lun_pp
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds energy state information at the landunit level.
  !-----------------------------------------------------------------------
  type, public :: landunit_energy_state
    ! temperature variables
    real(r8), pointer :: t_building       (:)   ! internal building temperature (K)
    real(r8), pointer :: taf              (:)   ! urban canopy air temperature (K)

  contains
    procedure, public :: Init    => lun_es_init
    procedure, public :: Restart => lun_es_restart
    procedure, public :: Clean   => lun_es_clean
  end type landunit_energy_state
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds energy flux information at the landunit level.
  !-----------------------------------------------------------------------
  type, public :: landunit_energy_flux
    ! temperature variables
    real(r8), pointer :: eflx_traffic      (:)   ! traffic sensible heat flux (W/m**2)
    real(r8), pointer :: eflx_wasteheat    (:)   ! sensible heat flux from domestic heating/cooling sources of waste heat (W/m**2)
    real(r8), pointer :: eflx_heat_from_ac (:)   ! sensible heat flux to be put back into canyon due to removal by AC (W/m**2)

  contains
    procedure, public :: Init    => lun_ef_init
    procedure, public :: Clean   => lun_ef_clean
  end type landunit_energy_flux
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds water state information at the landunit level.
  !-----------------------------------------------------------------------
  type, public :: landunit_water_state
    ! temperature variables
    real(r8), pointer :: qaf               (:)   ! urban canopy air specific humidity (kg H2O/kg moist air)

  contains
    procedure, public :: Init    => lun_ws_init
    procedure, public :: Restart => lun_ws_restart
    procedure, public :: Clean   => lun_ws_clean
  end type landunit_water_state
  
  !-----------------------------------------------------------------------
  ! declare the public instances of landunit-level data types
  !-----------------------------------------------------------------------
  type(landunit_energy_state)          , public, target :: lun_es    ! landunit energy state
  type(landunit_energy_flux )          , public, target :: lun_ef    ! landunit energy flux
  type(landunit_water_state )          , public, target :: lun_ws    ! landunit water state

  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean landunit energy state data structure
  !------------------------------------------------------------------------
  subroutine lun_es_init(this, begl, endl)
    !
    ! !ARGUMENTS:
    class(landunit_energy_state) :: this
    integer, intent(in) :: begl,endl
    !------------------------------------------------------------------------
    ! !LOCAL VARIABLES:
    integer :: l                        ! indices

    !-----------------------------------------------------------------------
    ! allocate for each member of lun_es
    !-----------------------------------------------------------------------
    allocate(this%t_building            (begl:endl))                      ; this%t_building            (:)   = nan
    allocate(this%taf                   (begl:endl))                      ; this%taf                   (:)   = nan

    !-----------------------------------------------------------------------
    ! initialize history fields for select members of lun_es
    !-----------------------------------------------------------------------
    this%t_building(begl:endl) = spval
    call hist_addfld1d(fname='TBUILD', units='K',  &
         avgflag='A', long_name='internal urban building temperature', &
         ptr_lunit=this%t_building, set_nourb=spval, l2g_scale_type='unity')
    
    !-----------------------------------------------------------------------
    ! cold-start initial conditions for lun_es
    !-----------------------------------------------------------------------
    do l = begl, endl 
       if (lun_pp%urbpoi(l)) then
          if (use_vancouver) then
             this%taf(l) = 297.56_r8
          else if (use_mexicocity) then
             this%taf(l) = 289.46_r8
          else
             this%taf(l) = 283._r8
          end if
       end if
    end do
    

  end subroutine lun_es_init

  !------------------------------------------------------------------------
  subroutine lun_es_restart(this, bounds, ncid, flag)
    ! 
    ! !DESCRIPTION:
    ! Read/Write landunit energy state information to/from restart file.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(landunit_energy_state) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   
    character(len=*) , intent(in)    :: flag   
    !
    ! !LOCAL VARIABLES:
    logical :: readvar   ! determine if variable is on initial file
    !-----------------------------------------------------------------------

    call restartvar(ncid=ncid, flag=flag, varname='taf', xtype=ncd_double, dim1name='landunit',                       &
         long_name='urban canopy air temperature', units='K',                                                         &
         interpinic_flag='interp', readvar=readvar, data=this%taf)
  end subroutine lun_es_restart

  !------------------------------------------------------------------------
  subroutine lun_es_clean(this)
    !
    ! !ARGUMENTS:
    class(landunit_energy_state) :: this
    !------------------------------------------------------------------------
    deallocate(this%t_building)
    deallocate(this%taf)
    
  end subroutine lun_es_clean

  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean landunit energy flux data structure
  !------------------------------------------------------------------------
  subroutine lun_ef_init(this, begl, endl)
    !
    ! !ARGUMENTS:
    class(landunit_energy_flux) :: this
    integer, intent(in) :: begl,endl
    !------------------------------------------------------------------------
    ! !LOCAL VARIABLES:
    integer :: l                        ! indices

    !-----------------------------------------------------------------------
    ! allocate for each member of lun_ef
    !-----------------------------------------------------------------------
    allocate( this%eflx_heat_from_ac   (begl:endl))             ; this%eflx_heat_from_ac   (:)   = nan
    allocate( this%eflx_traffic        (begl:endl))             ; this%eflx_traffic        (:)   = nan
    allocate( this%eflx_wasteheat      (begl:endl))             ; this%eflx_wasteheat      (:)   = nan

    !-----------------------------------------------------------------------
    ! cold-start initial conditions for lun_ef
    !-----------------------------------------------------------------------
    do l = begl, endl 
       if (lun_pp%urbpoi(l)) then
          this%eflx_traffic(l)   = spval
          this%eflx_wasteheat(l) = spval
       end if
    end do

  end subroutine lun_ef_init

  !------------------------------------------------------------------------
  subroutine lun_ef_clean(this)
    !
    ! !ARGUMENTS:
    class(landunit_energy_flux) :: this
    !------------------------------------------------------------------------
    deallocate(this%eflx_heat_from_ac)
    deallocate(this%eflx_traffic)
    deallocate(this%eflx_wasteheat)
    
  end subroutine lun_ef_clean
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean landunit water state data structure
  !------------------------------------------------------------------------
  subroutine lun_ws_init(this, begl, endl)
    !
    ! !ARGUMENTS:
    class(landunit_water_state) :: this
    integer, intent(in) :: begl,endl
    !------------------------------------------------------------------------
    ! !LOCAL VARIABLES:
    integer :: l                        ! indices

    !-----------------------------------------------------------------------
    ! allocate for each member of lun_ws
    !-----------------------------------------------------------------------
    allocate(this%qaf          (begl:endl))               ; this%qaf         (:)   = nan

    !-----------------------------------------------------------------------
    ! cold-start initial conditions for lun_ws
    !-----------------------------------------------------------------------
    do l = begl, endl 
       if (lun_pp%urbpoi(l)) then
          if (use_vancouver) then
             this%qaf(l) = 0.0111_r8
          else if (use_mexicocity) then
             this%qaf(l) = 0.00248_r8
          else
             this%qaf(l) = 1.e-4_r8 ! Arbitrary set since forc_q is not yet available
          end if
       end if
    end do

  end subroutine lun_ws_init

  !------------------------------------------------------------------------
  subroutine lun_ws_restart(this, bounds, ncid, flag)
    ! 
    ! !DESCRIPTION:
    ! Read/Write landunit water state information to/from restart file.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(landunit_water_state)      :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   
    character(len=*) , intent(in)    :: flag   
    !
    ! !LOCAL VARIABLES:
    logical :: readvar   ! determine if variable is on initial file
    !-----------------------------------------------------------------------

    call restartvar(ncid=ncid, flag=flag, varname='qaf', xtype=ncd_double, dim1name='landunit',                       &
         long_name='urban canopy specific humidity', units='kg/kg',                                                   &
         interpinic_flag='interp', readvar=readvar, data=this%qaf)

  end subroutine lun_ws_restart

  !------------------------------------------------------------------------
  subroutine lun_ws_clean(this)
    !
    ! !ARGUMENTS:
    class(landunit_water_state) :: this
    !------------------------------------------------------------------------
    deallocate(this%qaf)
    
  end subroutine lun_ws_clean


  

end module LandunitDataType

