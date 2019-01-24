module LandunitDataType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Landunit data type allocation and initialization
  ! -------------------------------------------------------- 
  !
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use clm_varpar     , only : nlevsno, nlevgrnd, nlevlak, nlevurb
  use clm_varcon     , only : spval, ispval
  use clm_varctl     , only : use_vancouver, use_mexicocity
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
    real(r8), pointer :: t_building           (:)   ! internal building temperature (K)
    real(r8), pointer :: taf                  (:)   ! urban canopy air temperature (K)

  contains
    procedure, public :: Init    => lun_es_init
    procedure, public :: Restart => lun_es_restart
    procedure, public :: Clean   => lun_es_clean
  end type landunit_energy_state
  
  !-----------------------------------------------------------------------
  ! declare the public instances of landunit-level data types
  !-----------------------------------------------------------------------
  type(landunit_energy_state)          , public, target :: lun_es    ! landunit energy state

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
  

end module LandunitDataType

