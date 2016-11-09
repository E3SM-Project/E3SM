
module upper_bc

!---------------------------------------------------------------------------------
! Module to compute the upper boundary condition for temperature (dry static energy)
! and trace gases. Uses the MSIS model, and SNOE and TIME GCM data.
!
! original code by Stacy Walters
! adapted by B. A. Boville
!---------------------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid,       only: pcols, pver, pverp
  use constituents, only: pcnst
  use physconst,    only: rair, mwdry
  use cam_logfile,  only: iulog
  use spmd_utils,   only: masterproc

  implicit none
  private
  save
!
! Public interfaces
!
  public :: ubc_defaultopts    ! set default values of namelist variables
  public :: ubc_setopts        ! get namelist input
  public :: ubc_init           ! global initialization
  public :: ubc_timestep_init  ! time step initialization
  public :: ubc_get_vals       ! get ubc values for this step

! Namelist variables
  character(len=256) :: snoe_ubc_file = ' '
  real(r8)           :: t_pert_ubc  = 0._r8
  real(r8)           :: no_xfac_ubc = 1._r8

  character(len=256) :: tgcm_ubc_file = ' '
  integer            :: tgcm_ubc_cycle_yr = 0
  integer            :: tgcm_ubc_fixed_ymd = 0
  integer            :: tgcm_ubc_fixed_tod = 0
  character(len=32)  :: tgcm_ubc_data_type = 'CYCLICAL'

#if ( defined WACCM_MOZART || defined WACCM_GHG )
  logical, parameter :: waccm = .true.
#else
  logical, parameter :: waccm = .false.
#endif

!================================================================================================
contains
!================================================================================================

subroutine ubc_defaultopts(tgcm_ubc_file_out, tgcm_ubc_data_type_out, tgcm_ubc_cycle_yr_out, tgcm_ubc_fixed_ymd_out, &
     tgcm_ubc_fixed_tod_out, snoe_ubc_file_out, t_pert_ubc_out, no_xfac_ubc_out)
!----------------------------------------------------------------------- 
! Purpose: Return default runtime options
!-----------------------------------------------------------------------

   real(r8), intent(out), optional         :: t_pert_ubc_out
   real(r8), intent(out), optional         :: no_xfac_ubc_out
   character(len=*), intent(out), optional :: tgcm_ubc_file_out
   character(len=*), intent(out), optional :: snoe_ubc_file_out
   integer         , intent(out), optional :: tgcm_ubc_cycle_yr_out
   integer         , intent(out), optional :: tgcm_ubc_fixed_ymd_out
   integer         , intent(out), optional :: tgcm_ubc_fixed_tod_out
   character(len=*), intent(out), optional :: tgcm_ubc_data_type_out

!-----------------------------------------------------------------------

   if ( present(tgcm_ubc_file_out) ) then
      tgcm_ubc_file_out = tgcm_ubc_file
   endif
   if ( present(tgcm_ubc_data_type_out) ) then
      tgcm_ubc_data_type_out = tgcm_ubc_data_type
   endif
   if ( present(tgcm_ubc_cycle_yr_out) ) then
      tgcm_ubc_cycle_yr_out = tgcm_ubc_cycle_yr
   endif
   if ( present(tgcm_ubc_fixed_ymd_out) ) then
      tgcm_ubc_fixed_ymd_out = tgcm_ubc_fixed_ymd
   endif
   if ( present(tgcm_ubc_fixed_tod_out) ) then
      tgcm_ubc_fixed_tod_out = tgcm_ubc_fixed_tod
   endif
   if ( present(snoe_ubc_file_out) ) then
      snoe_ubc_file_out = snoe_ubc_file
   endif
   if ( present(t_pert_ubc_out) ) then
      t_pert_ubc_out = t_pert_ubc
   endif
   if ( present(no_xfac_ubc_out) ) then
      no_xfac_ubc_out = no_xfac_ubc
   endif

end subroutine ubc_defaultopts

!================================================================================================

subroutine ubc_setopts(tgcm_ubc_file_in, tgcm_ubc_data_type_in, tgcm_ubc_cycle_yr_in, tgcm_ubc_fixed_ymd_in, &
     tgcm_ubc_fixed_tod_in, snoe_ubc_file_in, t_pert_ubc_in, no_xfac_ubc_in)
!----------------------------------------------------------------------- 
! Purpose: Set runtime options
!-----------------------------------------------------------------------

   use cam_abortutils, only : endrun

   real(r8), intent(in), optional         :: t_pert_ubc_in
   real(r8), intent(in), optional         :: no_xfac_ubc_in
   character(len=*), intent(in), optional :: tgcm_ubc_file_in
   character(len=*), intent(in), optional :: snoe_ubc_file_in
   integer         , intent(in), optional :: tgcm_ubc_cycle_yr_in
   integer         , intent(in), optional :: tgcm_ubc_fixed_ymd_in
   integer         , intent(in), optional :: tgcm_ubc_fixed_tod_in
   character(len=*), intent(in), optional :: tgcm_ubc_data_type_in

!-----------------------------------------------------------------------

   if ( present(tgcm_ubc_file_in) ) then
      tgcm_ubc_file = tgcm_ubc_file_in
   endif
   if ( present(tgcm_ubc_data_type_in) ) then
      tgcm_ubc_data_type = tgcm_ubc_data_type_in
   endif
   if ( present(tgcm_ubc_cycle_yr_in) ) then
      tgcm_ubc_cycle_yr = tgcm_ubc_cycle_yr_in
   endif
   if ( present(tgcm_ubc_fixed_ymd_in) ) then
      tgcm_ubc_fixed_ymd = tgcm_ubc_fixed_ymd_in
   endif
   if ( present(tgcm_ubc_fixed_tod_in) ) then
      tgcm_ubc_fixed_tod = tgcm_ubc_fixed_tod_in
   endif
   if ( present(snoe_ubc_file_in) ) then
      snoe_ubc_file = snoe_ubc_file_in
   endif
   if ( present(t_pert_ubc_in) ) then
      t_pert_ubc = t_pert_ubc_in
   endif
   if ( present(no_xfac_ubc_in) ) then
      no_xfac_ubc = no_xfac_ubc_in
      if( no_xfac_ubc < 0._r8 ) then
         write(iulog,*) 'ubc_setopts: no_xfac_ubc = ',no_xfac_ubc,' must be >= 0'
         call endrun
      end if
   endif

end subroutine ubc_setopts

!===============================================================================

  subroutine ubc_init()
!-----------------------------------------------------------------------
! Initialization of time independent fields for the upper boundary condition
! Calls initialization routine for MSIS, TGCM and SNOE
!-----------------------------------------------------------------------
    use mo_tgcm_ubc, only: tgcm_ubc_inti
    use mo_snoe,     only: snoe_inti
    use mo_msis_ubc, only: msis_ubc_inti
    use time_manager,only: get_step_size
    use constituents,only: cnst_get_ind, cnst_mw !Needed for ubc_flux
    use physics_buffer, only : physics_buffer_desc

!---------------------------Local workspace-----------------------------
    integer :: steps_per_day
    logical :: zonal_avg
!-----------------------------------------------------------------------

    if (.not.waccm) return

    steps_per_day = 86400/get_step_size()
    zonal_avg     = .false.

!-----------------------------------------------------------------------
!	... initialize the tgcm upper boundary module
!-----------------------------------------------------------------------
    call tgcm_ubc_inti( tgcm_ubc_file, tgcm_ubc_data_type, tgcm_ubc_cycle_yr, tgcm_ubc_fixed_ymd, tgcm_ubc_fixed_tod)
    if (masterproc) write(iulog,*) 'ubc_init: after tgcm_ubc_inti'

!-----------------------------------------------------------------------
!	... initialize the snoe module
!-----------------------------------------------------------------------
    call snoe_inti(snoe_ubc_file)
    if (masterproc) write(iulog,*) 'ubc_init: after snoe_inti'

!-----------------------------------------------------------------------
!	... initialize the msis module
!-----------------------------------------------------------------------
    call msis_ubc_inti( zonal_avg, steps_per_day )
    if (masterproc) write(iulog,*) 'ubc_init: after msis_ubc_inti'

  end subroutine ubc_init

!===============================================================================

  subroutine ubc_timestep_init(pbuf2d, state)
!-----------------------------------------------------------------------
! timestep dependent setting
!-----------------------------------------------------------------------

    use mo_solar_parms, only: solar_parms_get
    use mo_msis_ubc,    only: msis_timestep_init
    use mo_tgcm_ubc,    only: tgcm_timestep_init
    use mo_snoe,        only: snoe_timestep_init
    use physics_types,  only: physics_state
    use ppgrid,         only: begchunk, endchunk
    use physics_buffer, only : physics_buffer_desc

    type(physics_state), intent(in) :: state(begchunk:endchunk)                 
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    real(r8) :: ap, kp, f107, f107a

    if (.not.waccm) return

    call solar_parms_get( kp_s = kp, f107_s = f107, f107a_s = f107a, ap_s = ap )
    call msis_timestep_init( ap, f107, f107a )
    call tgcm_timestep_init( pbuf2d, state )
    call snoe_timestep_init( kp, f107 )

  end subroutine ubc_timestep_init

!===============================================================================

  subroutine ubc_get_vals (lchnk, ncol, ntop_molec, pint, zi, msis_temp, ubc_mmr, ubc_flux)
!-----------------------------------------------------------------------
! interface routine for vertical diffusion and pbl scheme
!-----------------------------------------------------------------------
    use mo_msis_ubc,  only: get_msis_ubc
    use mo_snoe,      only: set_no_ubc, ndx_no
    use mo_tgcm_ubc,  only: set_tgcm_ubc
    use cam_abortutils,   only: endrun
    use physconst,    only: rairv, mbarv
    use phys_control, only: waccmx_is
    use physconst,    only: avogad
    use constituents, only: cnst_get_ind, cnst_mw  ! Needed for ubc_flux
  
!------------------------------Arguments--------------------------------
    integer,  intent(in)  :: lchnk                 ! chunk identifier
    integer,  intent(in)  :: ncol                  ! number of atmospheric columns
    integer,  intent(in)  :: ntop_molec            ! top of molecular diffusion region (=1)
    real(r8), intent(in)  :: pint(pcols,pverp)     ! interface pressures
    real(r8), intent(in)  :: zi(pcols,pverp)       ! interface geoptl height above sfc

    real(r8), intent(out) :: ubc_mmr(pcols,pcnst)  ! upper bndy mixing ratios (kg/kg)
    real(r8), intent(out) :: msis_temp(pcols)      ! upper bndy temperature (K)
    real(r8), intent(out) :: ubc_flux(pcnst)       ! upper bndy flux (kg/s/m^2)

!---------------------------Local storage-------------------------------
    integer  :: ierr                               ! error flag for allocates
    integer  :: indx_H                             ! cnst index for H

    real(r8), parameter :: m2km = 1.e-3_r8         ! meter to km
    real(r8) :: rho_top(pcols)                     ! density at top interface
    real(r8) :: z_top(pcols)                       ! height of top interface (km)

!-----------------------------------------------------------------------

    if (.not.waccm) return

    ubc_mmr(:,:) = 0._r8
    ubc_flux(:) = 0._r8

    call get_msis_ubc( lchnk, ncol, msis_temp, ubc_mmr )
    if( t_pert_ubc /= 0._r8 ) then
       msis_temp(:ncol) = msis_temp(:ncol) + t_pert_ubc
       if( any( msis_temp(:ncol) < 0._r8 ) ) then
          write(iulog,*) 'ubc_get_vals: msis temp < 0 after applying offset = ',t_pert_ubc
          call endrun
       end if
    end if

    rho_top(:ncol) = pint(:ncol,ntop_molec) / (rairv(:ncol,ntop_molec,lchnk)*msis_temp(:ncol))
    z_top(:ncol)   = m2km * zi(:ncol,ntop_molec)
    
    if ( any( z_top(:ncol) < 180.0_r8 ) ) then  
       call set_no_ubc  ( lchnk, ncol, z_top, ubc_mmr, rho_top )
       if( ndx_no > 0 .and. no_xfac_ubc /= 1._r8 ) then
          ubc_mmr(:ncol,ndx_no) = no_xfac_ubc * ubc_mmr(:ncol,ndx_no)
       end if
    else
       ! waccmx
       ubc_mmr(:ncol,ndx_no) = 0.0_r8
    endif
    
    !--------------------------------------------------------------------------------------------
    ! For WACCM-X, calculate upper boundary H flux
    !--------------------------------------------------------------------------------------------
    if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then 
    
      call cnst_get_ind('H',  indx_H)

      ubc_flux(indx_H) = cnst_mw(indx_H)/ avogad * 0.9e12_r8      ! upward mass flux of H, kg/s/m^2
   
    endif     

    call set_tgcm_ubc( lchnk, ncol, ubc_mmr, mbarv(:,ntop_molec,lchnk))

  end subroutine ubc_get_vals

end module upper_bc
