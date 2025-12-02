module h3DType

  !
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use elm_varpar     , only : nlevsno, nlevgrnd, nlevlak ,nh3dc_per_lunit
  use elm_varcon     , only : spval, ispval
  use column_varcon  , only : is_hydrologically_active
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  type, public :: h3d_type


     real(r8), pointer :: qflx_evap_tot_lun(:,:)  ! 
     real(r8), pointer :: qflx_tran_veg_lun(:,:)  ! 
     real(r8), pointer :: qflx_rsub_sat_lun(:,:)  !soil saturation excess (mm)
     real(r8), pointer :: qflx_drain_lun(:,:)     !sub-surface drainage (mm/s)
     real(r8), pointer :: qflx_surf_lun(:,:)      !surface runoff (mm/s)
     real(r8), pointer :: qflx_charge_lun(:,:)    !aquifer recharge rate (mm/s)
     real(r8), pointer :: h2osfc_lun       (:,:)  !surface water (mm)
     real(r8), pointer :: h2osoi_liq_lun   (:,:)  !soil liquid water (kg/m2)
     real(r8), pointer :: eflx_lh_tot_lun  (:,:)  !total latent heat flux (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_sh_tot_lun  (:,:)  !total sensible heat flux (W/m**2) [+ to atm]

   contains

     procedure, public :: Init => h3d_vars_init
     procedure, public :: Clean => h3d_vars_clean

  end type h3d_type

  type(h3d_type), public, target :: h3d_vars !h3d data structure 
  !------------------------------------------------------------------------

contains
  
  !------------------------------------------------------------------------
  subroutine h3d_vars_init(this, begl, endl, begc, endc)
    use histFileMod    , only : hist_addfld1d,hist_addfld2d
    !
    ! !ARGUMENTS:
    class(h3d_type)  :: this
    integer, intent(in) :: begl,endl,begc,endc
    !------------------------------------------------------------------------

    allocate(this%qflx_evap_tot_lun (begl:endl,1:nh3dc_per_lunit))       ; this%qflx_evap_tot_lun (:,:)   = spval
    allocate(this%qflx_tran_veg_lun (begl:endl,1:nh3dc_per_lunit))       ; this%qflx_tran_veg_lun (:,:)   = spval
    allocate(this%qflx_rsub_sat_lun (begl:endl,1:nh3dc_per_lunit))       ; this%qflx_rsub_sat_lun (:,:)   = spval
    allocate(this%qflx_drain_lun    (begl:endl,1:nh3dc_per_lunit))       ; this%qflx_drain_lun    (:,:)   = spval
    allocate(this%qflx_surf_lun     (begl:endl,1:nh3dc_per_lunit))       ; this%qflx_surf_lun     (:,:)   = spval
    allocate(this%qflx_charge_lun   (begl:endl,1:nh3dc_per_lunit))       ; this%qflx_charge_lun   (:,:)   = spval
    allocate(this%h2osfc_lun        (begl:endl,1:nh3dc_per_lunit))       ; this%h2osfc_lun        (:,:)   = spval
    allocate(this%h2osoi_liq_lun    (begl:endl,1:nh3dc_per_lunit))       ; this%h2osoi_liq_lun    (:,:)   = spval
    allocate(this%eflx_lh_tot_lun   (begl:endl,1:nh3dc_per_lunit))       ; this%eflx_lh_tot_lun   (:,:)   = spval
    allocate(this%eflx_sh_tot_lun   (begl:endl,1:nh3dc_per_lunit))       ; this%eflx_sh_tot_lun   (:,:)   = spval

    call hist_addfld2d (fname='h3d_evap_tot',  units='mm SH2O/s', type2d='h3dc', &
         avgflag='A', long_name='h3d total evapotranspiration (vegetated landunits only)', &
         ptr_lunit=this%qflx_evap_tot_lun, l2g_scale_type='veg')
 
    call hist_addfld2d (fname='h3d_tran_veg',  units='mm SH2O/s', type2d='h3dc', &
         avgflag='A', long_name='h3d vegetation transpiration (vegetated landunits only)', &
         ptr_lunit=this%qflx_tran_veg_lun, l2g_scale_type='veg')

    call hist_addfld2d (fname='h3d_rsub_sat',  units='mm SH2O/s', type2d='h3dc', &
         avgflag='A', long_name='soil saturation excess (vegetated landunits only)', &
         ptr_lunit=this%qflx_rsub_sat_lun, l2g_scale_type='veg')
    
    call hist_addfld2d (fname='h3d_h2osfc',  units='mm SH2O', type2d='h3dc', &
         avgflag='A', long_name='soil saturation excess (vegetated landunits only)', &
         ptr_lunit=this%h2osfc_lun, l2g_scale_type='veg')

    call hist_addfld2d (fname='h3d_h2osoi_liq',  units='mm SH2O', type2d='h3dc', &
         avgflag='A', long_name='soil liquid water (vegetated landunits only)', &
         ptr_lunit=this%h2osoi_liq_lun, l2g_scale_type='veg')

    call hist_addfld2d (fname='h3d_lh_tot',  units='W/m2', type2d='h3dc', &
         avgflag='A', long_name='total latent heat flux  (vegetated landunits only)', &
         ptr_lunit=this%eflx_lh_tot_lun, l2g_scale_type='veg')

    call hist_addfld2d (fname='h3d_sh_tot',  units='W/m2', type2d='h3dc', &
         avgflag='A', long_name='total sensible heat flux  (vegetated landunits only)', &
         ptr_lunit=this%eflx_sh_tot_lun, l2g_scale_type='veg')

    call hist_addfld2d (fname='h3d_qdrain',  units='mm/s', type2d='h3dc', &
         avgflag='A', long_name='sub-surface drainage', &
         ptr_lunit=this%qflx_drain_lun, l2g_scale_type='veg')

    call hist_addfld2d (fname='h3d_qsurf',  units='mm/s', type2d='h3dc', &
         avgflag='A', long_name='surface runoff', &
         ptr_lunit=this%qflx_surf_lun, l2g_scale_type='veg')

    call hist_addfld2d (fname='h3d_qcharge',  units='mm/s', type2d='h3dc', &
         avgflag='A', long_name='aquifer recharge rate (vegetated landunits only)', &
         ptr_lunit=this%qflx_charge_lun, l2g_scale_type='veg')

  end subroutine h3d_vars_init

  !------------------------------------------------------------------------
  subroutine h3d_vars_clean(this)
    !
    ! !ARGUMENTS:
    class(h3d_type) :: this
    !------------------------------------------------------------------------

    deallocate(this%qflx_evap_tot_lun  )
    deallocate(this%qflx_tran_veg_lun  )
    deallocate(this%qflx_rsub_sat_lun  )
    deallocate(this%h2osfc_lun         )
    deallocate(this%h2osoi_liq_lun     )
    deallocate(this%eflx_lh_tot_lun    )
    deallocate(this%eflx_sh_tot_lun    )

  end subroutine h3d_vars_clean


end module h3DType
