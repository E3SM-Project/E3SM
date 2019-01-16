module VegetationDataType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Vegetation data type allocation and initialization
  ! -------------------------------------------------------- 
  !
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use clm_varpar     , only : nlevsno, nlevgrnd, nlevlak, nlevurb
  use clm_varcon     , only : spval, ispval
  use histFileMod    , only : hist_addfld1d, hist_addfld2d, no_snow_normal
  use ncdio_pio      , only : file_desc_t, ncd_double
  use decompMod      , only : bounds_type
  use restUtilMod
  use VegetationType , only : veg_pp
  use LandunitType   , only : lun_pp
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds energy state information at the vegetation level.
  !-----------------------------------------------------------------------
  type, public :: vegetation_energy_state
    ! temperature variables
    real(r8), pointer :: t_veg              (:)   ! vegetation temperature (K)
    real(r8), pointer :: t_ref2m            (:)   ! 2 m height surface air temperature (K)
    real(r8), pointer :: t_ref2m_r          (:)   ! rural 2 m height surface air temperature (K)
    real(r8), pointer :: t_ref2m_u          (:)   ! urban 2 m height surface air temperature (K)
    
  contains
    procedure, public :: Init    => veg_es_init
    procedure, public :: Restart => veg_es_restart
    procedure, public :: Clean   => veg_es_clean
  end type vegetation_energy_state
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds water state information at the vegetation level.
  !-----------------------------------------------------------------------
  type, public :: vegetation_water_state
    real(r8), pointer :: xxx      (:) => null() ! xxx (xxx)
  contains
    procedure, public :: Init  => init_veg_ws
    procedure, public :: Clean => clean_veg_ws
  end type vegetation_water_state
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds carbon state information at the vegetation level.
  !-----------------------------------------------------------------------
  type, public :: vegetation_carbon_state
    real(r8), pointer :: xxx      (:) => null() ! xxx (xxx)
  contains
    procedure, public :: Init  => init_veg_cs
    procedure, public :: Clean => clean_veg_cs
  end type vegetation_carbon_state
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds nitrogen state information at the vegetation level.
  !-----------------------------------------------------------------------
  type, public :: vegetation_nitrogen_state
    real(r8), pointer :: xxx      (:) => null() ! xxx (xxx)
  contains
    procedure, public :: Init  => init_veg_ns
    procedure, public :: Clean => clean_veg_ns
  end type vegetation_nitrogen_state
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds phosphorus state information at the vegetation level.
  !-----------------------------------------------------------------------
  type, public :: vegetation_phosphorus_state
    real(r8), pointer :: xxx      (:) => null() ! xxx (xxx)
  contains
    procedure, public :: Init  => init_veg_ps
    procedure, public :: Clean => clean_veg_ps
  end type vegetation_phosphorus_state

  !-----------------------------------------------------------------------
  ! Define the data structure that holds energy flux information at the vegetation level.
  !-----------------------------------------------------------------------
  type, public :: vegetation_energy_flux
    real(r8), pointer :: xxx      (:) => null() ! xxx (xxx)
  contains
    procedure, public :: Init  => init_veg_ef
    procedure, public :: Clean => clean_veg_ef
  end type vegetation_energy_flux
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds water flux information at the vegetation level.
  !-----------------------------------------------------------------------
  type, public :: vegetation_water_flux
    real(r8), pointer :: xxx      (:) => null() ! xxx (xxx)
  contains
    procedure, public :: Init  => init_veg_wf
    procedure, public :: Clean => clean_veg_wf
  end type vegetation_water_flux
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds carbon flux information at the vegetation level.
  !-----------------------------------------------------------------------
  type, public :: vegetation_carbon_flux
    real(r8), pointer :: xxx      (:) => null() ! xxx (xxx)
  contains
    procedure, public :: Init  => init_veg_cf
    procedure, public :: Clean => clean_veg_cf
  end type vegetation_carbon_flux
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds nitrogen flux information at the vegetation level.
  !-----------------------------------------------------------------------
  type, public :: vegetation_nitrogen_flux
    real(r8), pointer :: xxx      (:) => null() ! xxx (xxx)
  contains
    procedure, public :: Init  => init_veg_nf
    procedure, public :: Clean => clean_veg_nf
  end type vegetation_nitrogen_flux
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds phosphorus flux information at the vegetation level.
  !-----------------------------------------------------------------------
  type, public :: vegetation_phosphorus_flux
    real(r8), pointer :: xxx      (:) => null() ! xxx (xxx)
  contains
    procedure, public :: Init  => init_veg_pf
    procedure, public :: Clean => clean_veg_pf
  end type vegetation_phosphorus_flux

   
  !-----------------------------------------------------------------------
  ! declare the public instances of vegetation-level data types
  !-----------------------------------------------------------------------
  type(vegetation_energy_state)          , public, target :: veg_es    ! vegetation energy state
  type(vegetation_water_state)           , public, target :: veg_ws    ! vegetation water state
  type(vegetation_carbon_state)          , public, target :: veg_cs    ! vegetation carbon state
  type(vegetation_nitrogen_state)        , public, target :: veg_ns    ! vegetation nitrogen state
  type(vegetation_phosphorus_state)      , public, target :: veg_ps    ! vegetation phosphorus state
  type(vegetation_energy_flux)           , public, target :: veg_ef    ! vegetation energy flux
  type(vegetation_water_flux)            , public, target :: veg_wf    ! vegetation water flux
  type(vegetation_carbon_flux)           , public, target :: veg_cf    ! vegetation carbon flux
  type(vegetation_nitrogen_flux)         , public, target :: veg_nf    ! vegetation nitrogen flux
  type(vegetation_phosphorus_flux)       , public, target :: veg_pf    ! vegetation phosphorus flux

  !------------------------------------------------------------------------
  
  contains

  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean vegetation energy state data structure
  !------------------------------------------------------------------------
  subroutine veg_es_init(this, begp, endp)
    !
    ! !USES:
    use clm_varctl     , only : iulog, use_vancouver, use_mexicocity
    !
    ! !ARGUMENTS:
    class(vegetation_energy_state) :: this
    integer, intent(in) :: begp,endp
    !------------------------------------------------------------------------
    !
    ! !LOCAL VARIABLES:
    integer           :: p,c,l,j                        ! indices

    !-----------------------------------------------------------------------
    ! allocate for each member of veg_es
    !-----------------------------------------------------------------------
    allocate(this%t_veg            (begp:endp))                     ; this%t_veg              (:)   = nan
    allocate(this%t_ref2m          (begp:endp))                     ; this%t_ref2m            (:)   = nan
    allocate(this%t_ref2m_r        (begp:endp))                     ; this%t_ref2m_r          (:)   = nan
    allocate(this%t_ref2m_u        (begp:endp))                     ; this%t_ref2m_u          (:)   = nan

    !-----------------------------------------------------------------------
    ! initialize history fields for select members of veg_es
    !-----------------------------------------------------------------------
    this%t_veg(begp:endp) = spval
    call hist_addfld1d (fname='TV', units='K',  &
         avgflag='A', long_name='vegetation temperature', &
         ptr_patch=this%t_veg)

    this%t_ref2m(begp:endp) = spval
    call hist_addfld1d (fname='TSA', units='K',  &
         avgflag='A', long_name='2m air temperature', &
         ptr_patch=this%t_ref2m)

    this%t_ref2m_r(begp:endp) = spval
    call hist_addfld1d (fname='TSA_R', units='K',  &
         avgflag='A', long_name='Rural 2m air temperature', &
         ptr_patch=this%t_ref2m_r, set_spec=spval)

    this%t_ref2m_u(begp:endp) = spval
    call hist_addfld1d (fname='TSA_U', units='K',  &
         avgflag='A', long_name='Urban 2m air temperature', &
         ptr_patch=this%t_ref2m_u, set_nourb=spval)

    !-----------------------------------------------------------------------
    ! set cold-start initial values for select members of veg_es
    !-----------------------------------------------------------------------
    
    ! Set t_veg, t_ref2m, t_ref2m_u and tref2m_r 

    do p = begp, endp
       c = veg_pp%column(p)
       l = veg_pp%landunit(p)

       if (use_vancouver) then
          this%t_veg(p)   = 297.56
       else if (use_mexicocity) then
          this%t_veg(p)   = 289.46
       else
          this%t_veg(p)   = 283._r8
       end if

       if (use_vancouver) then
          this%t_ref2m(p) = 297.56
       else if (use_mexicocity) then
          this%t_ref2m(p) = 289.46
       else
          this%t_ref2m(p) = 283._r8
       end if

       if (lun_pp%urbpoi(l)) then
          if (use_vancouver) then
             this%t_ref2m_u(p) = 297.56
          else if (use_mexicocity) then
             this%t_ref2m_u(p) = 289.46
          else
             this%t_ref2m_u(p) = 283._r8
          end if
       else 
          if (.not. lun_pp%ifspecial(l)) then 
             if (use_vancouver) then
                this%t_ref2m_r(p) = 297.56
             else if (use_mexicocity) then
                this%t_ref2m_r(p) = 289.46
             else
                this%t_ref2m_r(p) = 283._r8
             end if
          else 
             this%t_ref2m_r(p) = spval
          end if
       end if

    end do ! veg loop

  end subroutine veg_es_init
    
  !------------------------------------------------------------------------
  subroutine veg_es_restart(this, bounds, ncid, flag)
    ! 
    ! !DESCRIPTION:
    ! Read/Write vegetation energy state information to/from restart file.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(vegetation_energy_state) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   
    character(len=*) , intent(in)    :: flag   
    !
    ! !LOCAL VARIABLES:
    logical :: readvar   ! determine if variable is on initial file
    !-----------------------------------------------------------------------

    call restartvar(ncid=ncid, flag=flag, varname='T_VEG', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='vegetation temperature', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_veg)

    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='2m height surface air temperature', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_ref2m)

    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_R', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='Rural 2m height surface air temperature', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_ref2m_r)

    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_U', xtype=ncd_double, dim1name='pft',                      &
         long_name='Urban 2m height surface air temperature', units='K',                                              &
         interpinic_flag='interp', readvar=readvar, data=this%t_ref2m_u)

  end subroutine veg_es_restart

  !------------------------------------------------------------------------
  subroutine veg_es_clean(this)
    !
    ! !ARGUMENTS:
    class(vegetation_energy_state) :: this
    !------------------------------------------------------------------------
  end subroutine veg_es_clean
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean vegetation water state data structure
  !------------------------------------------------------------------------
  subroutine init_veg_ws(this, begp, endp)
    !
    ! !ARGUMENTS:
    class(vegetation_water_state) :: this
    integer, intent(in) :: begp,endp
    !------------------------------------------------------------------------
    
  end subroutine init_veg_ws
    
  !------------------------------------------------------------------------
  subroutine clean_veg_ws(this)
    !
    ! !ARGUMENTS:
    class(vegetation_water_state) :: this
    !------------------------------------------------------------------------
    
  end subroutine clean_veg_ws
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean vegetation carbon state data structure
  !------------------------------------------------------------------------
  subroutine init_veg_cs(this, begp, endp)
    !
    ! !ARGUMENTS:
    class(vegetation_carbon_state) :: this
    integer, intent(in) :: begp,endp
    !------------------------------------------------------------------------
    
  end subroutine init_veg_cs
    
  !------------------------------------------------------------------------
  subroutine clean_veg_cs(this)
    !
    ! !ARGUMENTS:
    class(vegetation_carbon_state) :: this
    !------------------------------------------------------------------------
    
  end subroutine clean_veg_cs
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean vegetation nitrogen state data structure
  !------------------------------------------------------------------------
  subroutine init_veg_ns(this, begp, endp)
    !
    ! !ARGUMENTS:
    class(vegetation_nitrogen_state) :: this
    integer, intent(in) :: begp,endp
    !------------------------------------------------------------------------
    
  end subroutine init_veg_ns
    
  !------------------------------------------------------------------------
  subroutine clean_veg_ns(this)
    !
    ! !ARGUMENTS:
    class(vegetation_nitrogen_state) :: this
    !------------------------------------------------------------------------
    
  end subroutine clean_veg_ns
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean vegetation phosphorus state data structure
  !------------------------------------------------------------------------
  subroutine init_veg_ps(this, begp, endp)
    !
    ! !ARGUMENTS:
    class(vegetation_phosphorus_state) :: this
    integer, intent(in) :: begp,endp
    !------------------------------------------------------------------------
    
  end subroutine init_veg_ps
    
  !------------------------------------------------------------------------
  subroutine clean_veg_ps(this)
    !
    ! !ARGUMENTS:
    class(vegetation_phosphorus_state) :: this
    !------------------------------------------------------------------------
    
  end subroutine clean_veg_ps

  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean vegetation energy flux data structure
  !------------------------------------------------------------------------
  subroutine init_veg_ef(this, begp, endp)
    !
    ! !ARGUMENTS:
    class(vegetation_energy_flux) :: this
    integer, intent(in) :: begp,endp
    !------------------------------------------------------------------------
    
  end subroutine init_veg_ef
    
  !------------------------------------------------------------------------
  subroutine clean_veg_ef(this)
    !
    ! !ARGUMENTS:
    class(vegetation_energy_flux) :: this
    !------------------------------------------------------------------------
    
  end subroutine clean_veg_ef
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean vegetation water flux data structure
  !------------------------------------------------------------------------
  subroutine init_veg_wf(this, begp, endp)
    !
    ! !ARGUMENTS:
    class(vegetation_water_flux) :: this
    integer, intent(in) :: begp,endp
    !------------------------------------------------------------------------
    
  end subroutine init_veg_wf
    
  !------------------------------------------------------------------------
  subroutine clean_veg_wf(this)
    !
    ! !ARGUMENTS:
    class(vegetation_water_flux) :: this
    !------------------------------------------------------------------------
    
  end subroutine clean_veg_wf
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean vegetation carbon flux data structure
  !------------------------------------------------------------------------
  subroutine init_veg_cf(this, begp, endp)
    !
    ! !ARGUMENTS:
    class(vegetation_carbon_flux) :: this
    integer, intent(in) :: begp,endp
    !------------------------------------------------------------------------
    
  end subroutine init_veg_cf
    
  !------------------------------------------------------------------------
  subroutine clean_veg_cf(this)
    !
    ! !ARGUMENTS:
    class(vegetation_carbon_flux) :: this
    !------------------------------------------------------------------------
    
  end subroutine clean_veg_cf
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean vegetation nitrogen flux data structure
  !------------------------------------------------------------------------
  subroutine init_veg_nf(this, begp, endp)
    !
    ! !ARGUMENTS:
    class(vegetation_nitrogen_flux) :: this
    integer, intent(in) :: begp,endp
    !------------------------------------------------------------------------
    
  end subroutine init_veg_nf
    
  !------------------------------------------------------------------------
  subroutine clean_veg_nf(this)
    !
    ! !ARGUMENTS:
    class(vegetation_nitrogen_flux) :: this
    !------------------------------------------------------------------------
    
  end subroutine clean_veg_nf
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean vegetation phosphorus flux data structure
  !------------------------------------------------------------------------
  subroutine init_veg_pf(this, begp, endp)
    !
    ! !ARGUMENTS:
    class(vegetation_phosphorus_flux) :: this
    integer, intent(in) :: begp,endp
    !------------------------------------------------------------------------
    
  end subroutine init_veg_pf
    
  !------------------------------------------------------------------------
  subroutine clean_veg_pf(this)
    !
    ! !ARGUMENTS:
    class(vegetation_phosphorus_flux) :: this
    !------------------------------------------------------------------------
    
  end subroutine clean_veg_pf
  
    !------------------------------------------------------------------------
    
end module vegetationDataType
