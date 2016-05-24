module cam_diagnostics

!---------------------------------------------------------------------------------
! Module to compute a variety of diagnostics quantities for history files
!---------------------------------------------------------------------------------

use shr_kind_mod,  only: r8 => shr_kind_r8
use camsrfexch,    only: cam_in_t, cam_out_t
use physics_types, only: physics_state, physics_tend
use ppgrid,        only: pcols, pver, pverp, begchunk, endchunk
use physics_buffer, only: physics_buffer_desc, pbuf_add_field, dtype_r8, dyn_time_lvls, &
                          pbuf_get_field, pbuf_get_index, pbuf_old_tim_idx



use cam_history,   only: outfld, write_inithist, hist_fld_active
use constituents,  only: pcnst, cnst_name, cnst_longname, cnst_cam_outfld, ptendnam, dmetendnam, apcnst, bpcnst, &
                         cnst_get_ind
use chemistry,     only: chem_is
use dycore,        only: dycore_is
use phys_control,  only: phys_getopts
use wv_saturation, only: qsat, qsat_water, svp_ice
use time_manager,  only: is_first_step

use scamMod,       only: single_column, wfld
use cam_abortutils,    only: endrun

implicit none
private
save

! Public interfaces

public :: &
   diag_register,      &! register pbuf space
   diag_init,          &! initialization
   diag_allocate,      &! allocate memory for module variables
   diag_deallocate,    &! deallocate memory for module variables
   diag_conv_tend_ini, &! initialize convective tendency calcs
   diag_phys_writeout, &! output diagnostics of the dynamics
   diag_phys_tend_writeout, & ! output physics tendencies
   diag_state_b4_phys_write,& ! output state before physics execution
   diag_conv,          &! output diagnostics of convective processes
   diag_surf,          &! output diagnostics of the surface
   diag_export,        &! output export state
   diag_physvar_ic,    &
   diag_readnl          ! read namelist options

logical, public :: inithist_all = .false. ! Flag to indicate set of fields to be 
                                          ! included on IC file
                                          !  .false.  include only required fields
                                          !  .true.   include required *and* optional fields

! Private data

integer :: dqcond_num                     ! number of constituents to compute convective
character(len=16) :: dcconnam(pcnst)      ! names of convection tendencies
                                          ! tendencies for
real(r8), allocatable :: dtcond(:,:,:)    ! temperature tendency due to convection
type dqcond_t
   real(r8), allocatable :: cnst(:,:,:)   ! constituent tendency due to convection
end type dqcond_t
type(dqcond_t), allocatable :: dqcond(:)

character(len=8) :: diag_cnst_conv_tend = 'q_only' ! output constituent tendencies due to convection
                                                   ! 'none', 'q_only' or 'all'

logical          :: history_amwg                   ! output the variables used by the AMWG diag package
logical          :: history_vdiag                  ! output the variables used by the AMWG variability diag package
logical          :: history_eddy                   ! output the eddy variables
logical          :: history_budget                 ! output tendencies and state variables for CAM4
                                                   ! temperature, water vapor, cloud ice and cloud
                                                   ! liquid budgets.
integer          :: history_budget_histfile_num    ! output history file number for budget fields

!Physics buffer indices
integer  ::      qcwat_idx  = 0 
integer  ::      tcwat_idx  = 0 
integer  ::      lcwat_idx  = 0 
integer  ::      cld_idx    = 0 
integer  ::      concld_idx = 0 
integer  ::      tke_idx    = 0 
integer  ::      kvm_idx    = 0 
integer  ::      kvh_idx    = 0 
integer  ::      cush_idx   = 0 
integer  ::      t_ttend_idx = 0

integer  ::      prec_dp_idx  = 0
integer  ::      snow_dp_idx  = 0
integer  ::      prec_sh_idx  = 0
integer  ::      snow_sh_idx  = 0
integer  ::      prec_sed_idx = 0
integer  ::      snow_sed_idx = 0
integer  ::      prec_pcw_idx = 0
integer  ::      snow_pcw_idx = 0


integer :: tpert_idx=-1, qpert_idx=-1, pblh_idx=-1
logical :: prog_modal_aero 
contains

! ===============================================================================

subroutine diag_register
    
   ! Request physics buffer space for fields that persist across timesteps.
   call pbuf_add_field('T_TTEND', 'global', dtype_r8, (/pcols,pver,dyn_time_lvls/), t_ttend_idx)

end subroutine diag_register

!===============================================================================

subroutine diag_readnl(nlfile)
  use namelist_utils,  only: find_group_name
  use units,           only: getunit, freeunit
  use mpishorthand
  use spmd_utils,      only: masterproc

  character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

  ! Local variables
  integer :: unitn, ierr
  character(len=*), parameter :: subname = 'diag_readnl'

  namelist /cam_diag_opts/ diag_cnst_conv_tend
  !-----------------------------------------------------------------------------

  if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'cam_diag_opts', status=ierr)
      if (ierr == 0) then
         read(unitn, cam_diag_opts, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(diag_cnst_conv_tend, len(diag_cnst_conv_tend), mpichar,  0, mpicom)
#endif
   
end subroutine diag_readnl

!================================================================================================

subroutine diag_init()

  ! Declare the history fields for which this module contains outfld calls.

   use cam_history,        only: addfld, horiz_only, add_default
   use constituent_burden, only: constituent_burden_init
   use cam_control_mod,    only: moist_physics, ideal_phys
   use tidal_diag,         only: tidal_diag_init 

   integer :: k, m
   ! Note - this is a duplication of information in ice_constants 
   ! Cannot put in a use statement if want to swap ice models to cice4
   integer, parameter :: plevmx = 4       ! number of subsurface levels
   character(len=8), parameter :: tsnam(plevmx) = (/ 'TS1', 'TS2', 'TS3', 'TS4' /)
   integer :: ixcldice, ixcldliq ! constituent indices for cloud liquid and ice water.
   integer :: ierr

   call phys_getopts(prog_modal_aero_out = prog_modal_aero )

   ! outfld calls in diag_phys_writeout

   call addfld ('NSTEP',horiz_only,    'A','timestep','Model timestep')
   call addfld ('PHIS',horiz_only,    'I','m2/s2','Surface geopotential')

   call addfld ('PS',horiz_only,    'A','Pa','Surface pressure')
   call addfld ('T',(/ 'lev' /), 'A','K','Temperature')
   call addfld ('U',(/ 'lev' /), 'A','m/s','Zonal wind')
   call addfld ('V',(/ 'lev' /), 'A','m/s','Meridional wind')
   call addfld (cnst_name(1),(/ 'lev' /), 'A','kg/kg',cnst_longname(1))

   ! State before physics
   call addfld ('TBP',(/ 'lev' /), 'A','K','Temperature (before physics)'       )
   call addfld (bpcnst(1) ,(/ 'lev' /), 'A','kg/kg',cnst_longname(1)//' (before physics)')
   ! State after physics
   call addfld ('TAP',(/ 'lev' /), 'A','K','Temperature (after physics)'       )
   call addfld ('UAP',(/ 'lev' /), 'A','m/s','Zonal wind (after physics)'        )
   call addfld ('VAP',(/ 'lev' /), 'A','m/s','Meridional wind (after physics)'   )
   call addfld (apcnst(1) ,(/ 'lev' /), 'A','kg/kg',cnst_longname(1)//' (after physics)')
   if ( dycore_is('LR') .or. dycore_is('SE') ) then
      call addfld ('TFIX',horiz_only,    'A'     ,'K/s','T fixer (T equivalent of Energy correction)')
      call addfld ('PTTEND_RESID',(/ 'lev' /), 'A'     ,'K/s',&
                   'T-tendency due to BAB kluge at end of tphysac (diagnostic not part of T-budget)' )
   end if
   call addfld ('TTEND_TOT',(/ 'lev' /), 'A','K/s' ,'Total temperature tendency'   )
  
   ! column burdens for all constituents except water vapor
   call constituent_burden_init

   call addfld ('Z3',(/ 'lev' /), 'A','m','Geopotential Height (above sea level)')
   call addfld ('Z1000',horiz_only,    'A','m','Geopotential Z at 700 mbar pressure surface')
   call addfld ('Z700',horiz_only,    'A','m','Geopotential Z at 700 mbar pressure surface')
   call addfld ('Z500',horiz_only,    'A','m','Geopotential Z at 500 mbar pressure surface')
   call addfld ('Z300',horiz_only,    'A','m','Geopotential Z at 300 mbar pressure surface')
   call addfld ('Z200',horiz_only,    'A','m','Geopotential Z at 200 mbar pressure surface')
   call addfld ('Z100',horiz_only,    'A','m','Geopotential Z at 100 mbar pressure surface')
   call addfld ('Z050',horiz_only,    'A','m','Geopotential Z at 50 mbar pressure surface')

   call addfld ('ZZ',(/ 'lev' /), 'A','m2','Eddy height variance' )
   call addfld ('VZ',(/ 'lev' /), 'A','m2/s','Meridional transport of geopotential energy')
   call addfld ('VT',(/ 'lev' /), 'A','K m/s   ','Meridional heat transport')
   call addfld ('VU',(/ 'lev' /), 'A','m2/s2','Meridional flux of zonal momentum' )
   call addfld ('VV',(/ 'lev' /), 'A','m2/s2','Meridional velocity squared' )

   if(prog_modal_aero) then !Only for prognostic aerosols
      call addfld ('bc_a1_2',(/ 'lev' /), 'A','kg2/kg2','bc_a1 squared')
      call addfld ('dst_a1_2',(/ 'lev' /), 'A','kg2/kg2','dst_a1 squared')
      call addfld ('dst_a3_2',(/ 'lev' /), 'A','kg2/kg2','dst_a3 squared')
      call addfld ('ncl_a1_2',(/ 'lev' /), 'A','kg2/kg2','ncl_a1 squared')
      call addfld ('ncl_a2_2',(/ 'lev' /), 'A','kg2/kg2','ncl_a2 squared')
      call addfld ('ncl_a3_2',(/ 'lev' /), 'A','kg2/kg2','ncl_a3 squared')
      call addfld ('so4_a1_2',(/ 'lev' /), 'A','kg2/kg2','so4_a1 squared')
      call addfld ('so4_a2_2',(/ 'lev' /), 'A','kg2/kg2','so4_a2 squared')
      call addfld ('so4_a3_2',(/ 'lev' /), 'A','kg2/kg2','so4_a3 squared')
      call addfld ('soa_a1_2',(/ 'lev' /), 'A','kg2/kg2','soa_a1 squared')
      call addfld ('soa_a2_2',(/ 'lev' /), 'A','kg2/kg2','soa_a2 squared')
      call addfld ('pom_a1_2',(/ 'lev' /), 'A','kg2/kg2','pom_a1 squared')
      call addfld ('Vbc_a1',(/ 'lev' /), 'A','m/skg/kg','Meridional bc_a1 transport')
      call addfld ('Vdst_a1',(/ 'lev' /), 'A','m/skg/kg','Meridional dst_a1 transport')
      call addfld ('Vdst_a3',(/ 'lev' /), 'A','m/skg/kg','Meridional dst_a3 transport')
      call addfld ('Vncl_a1',(/ 'lev' /), 'A','m/skg/kg','Meridional ncl_a1 transport')
      call addfld ('Vncl_a2',(/ 'lev' /), 'A','m/skg/kg','Meridional ncl_a2 transport')
      call addfld ('Vncl_a3',(/ 'lev' /), 'A','m/skg/kg','Meridional ncl_a3 transport')
      call addfld ('Vso4_a1',(/ 'lev' /), 'A','m/skg/kg','Meridional so4_a1 transport')
      call addfld ('Vso4_a2',(/ 'lev' /), 'A','m/skg/kg','Meridional so4_a2 transport')
      call addfld ('Vso4_a3',(/ 'lev' /), 'A','m/skg/kg','Meridional so4_a3 transport')
      call addfld ('Vsoa_a1',(/ 'lev' /), 'A','m/skg/kg','Meridional soa_a1 transport')
      call addfld ('Vsoa_a2',(/ 'lev' /), 'A','m/skg/kg','Meridional soa_a2 transport')
      call addfld ('Vpom_a1',(/ 'lev' /), 'A','m/skg/kg','Meridional pom_a1 transport')
   endif
   call addfld ('VQ',(/ 'lev' /), 'A','m/skg/kg','Meridional water transport')
   call addfld ('QQ',(/ 'lev' /), 'A','kg2/kg2','Eddy moisture variance')
   call addfld ('OMEGAV',(/ 'lev' /) ,'A','m Pa/s2 ','Vertical flux of meridional momentum' )
   call addfld ('OMGAOMGA',(/ 'lev' /) ,'A','Pa2/s2','Vertical flux of vertical momentum' )
   call addfld ('OMEGAQ',(/ 'lev' /) ,'A','kgPa/kgs','Vertical water transport' )

   call addfld ('UU',(/ 'lev' /), 'A','m2/s2','Zonal velocity squared' )
   call addfld ('WSPEED',(/ 'lev' /), 'X','m/s','Horizontal total wind speed maximum' )
   call addfld ('WSPDSRFMX',horiz_only,    'X','m/s','Horizontal total wind speed maximum at the surface' )
   call addfld ('WSPDSRFAV',horiz_only,    'A','m/s','Horizontal total wind speed average at the surface' )

   call addfld ('OMEGA',(/ 'lev' /), 'A','Pa/s','Vertical velocity (pressure)')
   call addfld ('OMEGAT',(/ 'lev' /), 'A','K Pa/s  ','Vertical heat flux' )
   call addfld ('OMEGAU',(/ 'lev' /), 'A','m Pa/s2 ','Vertical flux of zonal momentum' )
   call addfld ('OMEGA850',horiz_only,    'A','Pa/s','Vertical velocity at 850 mbar pressure surface')
   call addfld ('OMEGA500',horiz_only,    'A','Pa/s','Vertical velocity at 500 mbar pressure surface')

   call addfld ('MQ',(/ 'lev' /), 'A','kg/m2','Water vapor mass in layer')
   call addfld ('TMQ',horiz_only,    'A','kg/m2','Total (vertically integrated) precipitable water')
   call addfld ('TUQ',horiz_only,    'A','kg/m/s','Total (vertically integrated) zonal water flux')
   call addfld ('TVQ',horiz_only,    'A','kg/m/s','Total (vertically integrated) meridional water flux')
   call addfld ('RELHUM',(/ 'lev' /), 'A','percent','Relative humidity')
   call addfld ('RHW',(/ 'lev' /), 'A','percent'   ,'Relative humidity with respect to liquid')
   call addfld ('RHI',(/ 'lev' /), 'A','percent'   ,'Relative humidity with respect to ice')
   call addfld ('RHCFMIP',(/ 'lev' /), 'A','percent' ,'Relative humidity with respect to water above 273 K, ice below 273 K')
   call addfld ('PSL',horiz_only,    'A','Pa','Sea level pressure')

   call addfld ('T850',horiz_only,    'A','K','Temperature at 850 mbar pressure surface')
   call addfld ('T500',horiz_only,    'A','K','Temperature at 500 mbar pressure surface')
   call addfld ('T300',horiz_only,    'A','K','Temperature at 300 mbar pressure surface')
   call addfld ('T200',horiz_only,    'A','K','Temperature at 200 mbar pressure surface')
   call addfld ('Q850',horiz_only,    'A','kg/kg','Specific Humidity at 850 mbar pressure surface')
   call addfld ('Q200',horiz_only,    'A','kg/kg','Specific Humidity at 700 mbar pressure surface')
   call addfld ('U850',horiz_only,    'A','m/s','Zonal wind at 850 mbar pressure surface')
   call addfld ('U250',horiz_only,    'A','m/s','Zonal wind at 250 mbar pressure surface')
   call addfld ('U200',horiz_only,    'A','m/s','Zonal wind at 200 mbar pressure surface')
   call addfld ('U010',horiz_only,    'A','m/s','Zonal wind at  10 mbar pressure surface')
   call addfld ('V850',horiz_only,    'A','m/s','Meridional wind at 850 mbar pressure surface')
   call addfld ('V200',horiz_only,    'A','m/s','Meridional wind at 200 mbar pressure surface')
   call addfld ('V250',horiz_only,    'A','m/s','Meridional wind at 250 mbar pressure surface')

   call addfld ('TT',(/ 'lev' /), 'A','K2','Eddy temperature variance' )

   call addfld ('UBOT',horiz_only,    'A','m/s','Lowest model level zonal wind')
   call addfld ('VBOT',horiz_only,    'A','m/s','Lowest model level meridional wind')
   call addfld ('QBOT',horiz_only,    'A','kg/kg','Lowest model level water vapor mixing ratio')
   call addfld ('ZBOT',horiz_only,    'A','m','Lowest model level height')

   call addfld ('ATMEINT',horiz_only, 'A','J/m2','Vertically integrated total atmospheric energy ')

   call addfld ('T1000',horiz_only,   'A','K','Temperature at 1000 mbar pressure surface')
   call addfld ('T925',horiz_only,   'A','K','Temperature at 925 mbar pressure surface')   
   call addfld ('T700',horiz_only,   'A','K','Temperature at 700 mbar pressure surface')
   call addfld ('T010',horiz_only,   'A','K','Temperature at 10 mbar pressure surface')
   call addfld ('Q1000',horiz_only,   'A','kg/kg','Specific Humidity at 1000 mbar pressure surface')   
   call addfld ('Q925',horiz_only,   'A','kg/kg','Specific Humidity at 925 mbar pressure surface')

   call addfld ('T7001000',horiz_only,   'A','K','Temperature difference 700 mb - 1000 mb')
   call addfld ('TH7001000',horiz_only,   'A','K','Theta difference 700 mb - 1000 mb')
   call addfld ('THE7001000',horiz_only,   'A','K','ThetaE difference 700 mb - 1000 mb')

   call addfld ('T8501000',horiz_only,   'A','K','Temperature difference 850 mb - 1000 mb')
   call addfld ('TH8501000',horiz_only,   'A','K','Theta difference 850 mb - 1000 mb')   
   call addfld ('THE8501000',horiz_only,   'A','K','ThetaE difference 850 mb - 1000 mb')
   call addfld ('T9251000',horiz_only,   'A','K','Temperature difference 925 mb - 1000 mb') 
   call addfld ('TH9251000',horiz_only,   'A','K','Theta difference 925 mb - 1000 mb')   
   call addfld ('THE9251000',horiz_only,   'A','K','ThetaE difference 925 mb - 1000 mb') 

   ! This field is added by radiation when full physics is used
   if ( ideal_phys )then
      call addfld('QRS', (/ 'lev' /), 'A', 'K/s', 'Solar heating rate')
   end if
 
   ! ----------------------------
   ! determine default variables
   ! ----------------------------
   call phys_getopts(history_amwg_out   = history_amwg    , &
                     history_vdiag_out  = history_vdiag   , &
                     history_eddy_out   = history_eddy    , &
                     history_budget_out = history_budget  , &
                     history_budget_histfile_num_out = history_budget_histfile_num)

   if (history_amwg) then
      call add_default ('PHIS    '  , 1, ' ')
      call add_default ('PS      '  , 1, ' ')
      call add_default ('T       '  , 1, ' ')
      call add_default ('U       '  , 1, ' ')
      call add_default ('V       '  , 1, ' ')
      call add_default (cnst_name(1), 1, ' ')
      call add_default ('Z3      '  , 1, ' ')
      call add_default ('OMEGA   '  , 1, ' ')
      call add_default ('VT      ', 1, ' ')
      call add_default ('VU      ', 1, ' ')
      call add_default ('VV      ', 1, ' ')
      call add_default ('VQ      ', 1, ' ')

      if(prog_modal_aero) then !Only for prognostic aerosols
         call add_default ('Vbc_a1  ', 1, ' ')
         call add_default ('Vdst_a1 ', 1, ' ')
         call add_default ('Vdst_a3 ', 1, ' ')
         call add_default ('Vncl_a1 ', 1, ' ')
         call add_default ('Vncl_a2 ', 1, ' ')
         call add_default ('Vncl_a3 ', 1, ' ')
         call add_default ('Vso4_a1 ', 1, ' ')
         call add_default ('Vso4_a2 ', 1, ' ')
         call add_default ('Vso4_a3 ', 1, ' ')
         call add_default ('Vsoa_a1 ', 1, ' ')
         call add_default ('Vsoa_a2 ', 1, ' ')
         call add_default ('Vpom_a1 ', 1, ' ')
      endif
      call add_default ('UU      ', 1, ' ')
      call add_default ('OMEGAT  ', 1, ' ')
      if(prog_modal_aero) then !Only for prognostic aerosols
         call add_default ('bc_a1_2 ', 1, ' ')
         call add_default ('dst_a1_2', 1, ' ')
         call add_default ('dst_a3_2', 1, ' ')
         call add_default ('ncl_a1_2', 1, ' ')
         call add_default ('ncl_a2_2', 1, ' ')
         call add_default ('ncl_a3_2', 1, ' ')
         call add_default ('so4_a1_2', 1, ' ')
         call add_default ('so4_a2_2', 1, ' ')
         call add_default ('so4_a3_2', 1, ' ')
         call add_default ('soa_a1_2', 1, ' ')
         call add_default ('soa_a2_2', 1, ' ')
         call add_default ('pom_a1_2', 1, ' ')
      endif
      call add_default ('TMQ     ', 1, ' ')
      call add_default ('PSL     ', 1, ' ')
      if (moist_physics) then
         call add_default ('RELHUM  ', 1, ' ')
      end if
   end if
   
   if (history_vdiag) then
     call add_default ('U200', 2, ' ')
     call add_default ('V200', 2, ' ')
     call add_default ('U850', 2, ' ')
     call add_default ('U200', 3, ' ')
     call add_default ('U850', 3, ' ')
     call add_default ('OMEGA500', 3, ' ')
   end if

   if (history_eddy) then
      call add_default ('VT      ', 1, ' ')
      call add_default ('VU      ', 1, ' ')
      call add_default ('VV      ', 1, ' ')
      call add_default ('VQ      ', 1, ' ')
      call add_default ('UU      ', 1, ' ')
      call add_default ('OMEGAT  ', 1, ' ')
      call add_default ('OMEGAQ  ', 1, ' ')
      call add_default ('OMEGAU  ', 1, ' ')
      call add_default ('OMEGAV  ', 1, ' ')
   endif

   if ( history_budget ) then
      call add_default ('PHIS    '  , history_budget_histfile_num, ' ')
      call add_default ('PS      '  , history_budget_histfile_num, ' ')
      call add_default ('T       '  , history_budget_histfile_num, ' ')
      call add_default ('U       '  , history_budget_histfile_num, ' ')
      call add_default ('V       '  , history_budget_histfile_num, ' ')
      call add_default (cnst_name(1), history_budget_histfile_num, ' ')
      call add_default ('TTEND_TOT' , history_budget_histfile_num, ' ')

      ! State before physics (FV)
      call add_default ('TBP     '  , history_budget_histfile_num, ' ')
      call add_default (bpcnst(1)   , history_budget_histfile_num, ' ')
      ! State after physics (FV)
      call add_default ('TAP     '  , history_budget_histfile_num, ' ')
      call add_default ('UAP     '  , history_budget_histfile_num, ' ')
      call add_default ('VAP     '  , history_budget_histfile_num, ' ')
      call add_default (apcnst(1)   , history_budget_histfile_num, ' ')
      if ( dycore_is('LR') .or. dycore_is('SE') ) then
         call add_default ('TFIX    '    , history_budget_histfile_num, ' ')
         call add_default ('PTTEND_RESID', history_budget_histfile_num, ' ')
      end if
   end if

   ! This field is added by radiation when full physics is used
   if ( ideal_phys )then
      call add_default('QRS     ', 1, ' ')
   end if

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Exit here for adiabatic/ideal physics cases !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if (.not. moist_physics) return


   call addfld ('PDELDRY',(/ 'lev' /), 'A','Pa','Dry pressure difference between levels')
   call addfld ('PSDRY',horiz_only,    'A','Pa','Surface pressure')

   if (chem_is('waccm_ghg') .or. chem_is('waccm_mozart') .or. chem_is('waccm_mozart_mam3')) then
      call add_default ('PS      ', 2, ' ')
      call add_default ('T       ', 2, ' ')
   end if

   ! outfld calls in diag_conv

   call cnst_get_ind('CLDLIQ', ixcldliq)
   call cnst_get_ind('CLDICE', ixcldice)
   call addfld ('DTCOND',(/ 'lev' /), 'A','K/s','T tendency - moist processes')
   call addfld ('DTCOND_24_COS',(/ 'lev' /), 'A','K/s','T tendency - moist processes 24hr. cos coeff.')
   call addfld ('DTCOND_24_SIN',(/ 'lev' /), 'A','K/s','T tendency - moist processes 24hr. sin coeff.')
   call addfld ('DTCOND_12_COS',(/ 'lev' /), 'A','K/s','T tendency - moist processes 12hr. cos coeff.')
   call addfld ('DTCOND_12_SIN',(/ 'lev' /), 'A','K/s','T tendency - moist processes 12hr. sin coeff.')

   ! determine number of constituents for which convective tendencies must be computed
   if (history_budget) then
      dqcond_num = pcnst
   else
      if (diag_cnst_conv_tend == 'none')   dqcond_num = 0
      if (diag_cnst_conv_tend == 'q_only') dqcond_num = 1
      if (diag_cnst_conv_tend == 'all')    dqcond_num = pcnst
   end if

   do m = 1, dqcond_num
      dcconnam(m) = 'DC'//cnst_name(m)
   end do

   if (diag_cnst_conv_tend == 'q_only' .or. diag_cnst_conv_tend == 'all' .or. history_budget) then
      call addfld (dcconnam(1),(/ 'lev' /),'A', 'kg/kg/s',trim(cnst_name(1))//' tendency due to moist processes')
      if ( diag_cnst_conv_tend == 'q_only' .or. diag_cnst_conv_tend == 'all' ) then
         call add_default (dcconnam(1),                           1, ' ')
      end if
      if( history_budget ) then
         call add_default (dcconnam(1), history_budget_histfile_num, ' ')
      end if
      if (diag_cnst_conv_tend == 'all' .or. history_budget) then
         do m = 2, pcnst
            call addfld (dcconnam(m),(/ 'lev' /),'A', 'kg/kg/s',trim(cnst_name(m))//' tendency due to moist processes')
            if( diag_cnst_conv_tend == 'all' ) then
               call add_default (dcconnam(m),                           1, ' ')
            end if
            if( history_budget .and. (m == ixcldliq .or. m == ixcldice) ) then
               call add_default (dcconnam(m), history_budget_histfile_num, ' ')
            end if
         end do
      end if
   end if

   call addfld ('PRECL',horiz_only,    'A','m/s','Large-scale (stable) precipitation rate (liq + ice)'                )
   call addfld ('PRECC',horiz_only,    'A','m/s','Convective precipitation rate (liq + ice)'                          )
   call addfld ('PRECT',horiz_only,    'A','m/s','Total (convective and large-scale) precipitation rate (liq + ice)'  )
   call addfld ('PREC_PCW',horiz_only,    'A','m/s','LS_pcw precipitation rate')
   call addfld ('PREC_zmc',horiz_only,    'A','m/s','CV_zmc precipitation rate')
   call addfld ('PRECTMX',horiz_only,    'X','m/s','Maximum (convective and large-scale) precipitation rate (liq+ice)'  )
   call addfld ('PRECSL',horiz_only,    'A','m/s','Large-scale (stable) snow rate (water equivalent)'                  )
   call addfld ('PRECSC',horiz_only,    'A','m/s','Convective snow rate (water equivalent)'                            )
   call addfld ('PRECCav',horiz_only,    'A','m/s','Average large-scale precipitation (liq + ice)'                      )
   call addfld ('PRECLav',horiz_only,    'A','m/s','Average convective precipitation  (liq + ice)'                      )

   ! outfld calls in diag_surf

   call addfld ('SHFLX',horiz_only,    'A','W/m2','Surface sensible heat flux')
   call addfld ('LHFLX',horiz_only,    'A','W/m2','Surface latent heat flux')
   call addfld ('QFLX',horiz_only,    'A','kg/m2/s','Surface water flux')

   call addfld ('TAUX',horiz_only,    'A','N/m2','Zonal surface stress')
   call addfld ('TAUY',horiz_only,    'A','N/m2','Meridional surface stress')
   call addfld ('TREFHT',horiz_only,    'A','K','Reference height temperature')
   call addfld ('TREFHTMN',horiz_only,    'M','K','Minimum reference height temperature over output period')
   call addfld ('TREFHTMX',horiz_only,    'X','K','Maximum reference height temperature over output period')
   call addfld ('QREFHT',horiz_only,    'A','kg/kg','Reference height humidity')
   call addfld ('U10',horiz_only,    'A','m/s','10m wind speed')
   call addfld ('RHREFHT',horiz_only,    'A','fraction','Reference height relative humidity')

   call addfld ('LANDFRAC',horiz_only,    'A','fraction','Fraction of sfc area covered by land')
   call addfld ('ICEFRAC',horiz_only,    'A','fraction','Fraction of sfc area covered by sea-ice')
   call addfld ('OCNFRAC',horiz_only,    'A','fraction','Fraction of sfc area covered by ocean')

   call addfld ('TREFMNAV',horiz_only,    'A','K','Average of TREFHT daily minimum')
   call addfld ('TREFMXAV',horiz_only,    'A','K','Average of TREFHT daily maximum')

   call addfld ('TS',horiz_only,    'A','K','Surface temperature (radiative)')
   call addfld ('TSMN',horiz_only,    'M','K','Minimum surface temperature over output period')
   call addfld ('TSMX',horiz_only,    'X','K','Maximum surface temperature over output period')
   call addfld ('SNOWHLND',horiz_only,    'A','m','Water equivalent snow depth')
   call addfld ('SNOWHICE',horiz_only,    'A','m','Snow depth over ice', fill_value = 1.e30_r8)
   call addfld ('TBOT',horiz_only,    'A','K','Lowest model level temperature')

   call addfld ('ASDIR',       horiz_only,    'A',   '1','albedo: shortwave, direct')
   call addfld ('ASDIF',       horiz_only,    'A',   '1','albedo: shortwave, diffuse')
   call addfld ('ALDIR',       horiz_only,    'A',   '1','albedo: longwave, direct')
   call addfld ('ALDIF',       horiz_only,    'A',   '1','albedo: longwave, diffuse')
   call addfld ('SST',       horiz_only,    'A',     'K','sea surface temperature')

   ! defaults
   if (history_amwg) then
       call add_default ('DTCOND  ', 1, ' ')
       call add_default ('PRECL   ', 1, ' ')
       call add_default ('PRECC   ', 1, ' ')
       call add_default ('PRECSL  ', 1, ' ')
       call add_default ('PRECSC  ', 1, ' ')
       call add_default ('SHFLX   ', 1, ' ')
       call add_default ('LHFLX   ', 1, ' ')
       call add_default ('QFLX    ', 1, ' ')
       call add_default ('TAUX    ', 1, ' ')
       call add_default ('TAUY    ', 1, ' ')
       call add_default ('TREFHT  ', 1, ' ')
       call add_default ('LANDFRAC', 1, ' ')
       call add_default ('OCNFRAC ', 1, ' ')
       call add_default ('QREFHT  ', 1, ' ')
       call add_default ('U10     ', 1, ' ')
       call add_default ('ICEFRAC ', 1, ' ')
       call add_default ('TS      ', 1, ' ')
       call add_default ('TSMN    ', 1, ' ')
       call add_default ('TSMX    ', 1, ' ')
       call add_default ('SNOWHLND', 1, ' ')
       call add_default ('SNOWHICE', 1, ' ')
    endif

    if (history_vdiag) then
        call add_default ('PRECT   ', 2, ' ')
        call add_default ('PRECT   ', 3, ' ')
        call add_default ('PRECT   ', 4, ' ')
    end if

   ! outfld calls in diag_phys_tend_writeout

   call addfld ('PTTEND'   ,(/ 'lev' /), 'A','K/s','T total physics tendency'                             )
   call addfld (ptendnam(       1),(/ 'lev' /), 'A',  'kg/kg/s',trim(cnst_name(       1))//' total physics tendency '      )
   call addfld (ptendnam(ixcldliq),(/ 'lev' /), 'A',  'kg/kg/s',trim(cnst_name(ixcldliq))//' total physics tendency '      )
   call addfld (ptendnam(ixcldice),(/ 'lev' /), 'A',  'kg/kg/s',trim(cnst_name(ixcldice))//' total physics tendency '      )
   if ( dycore_is('LR') )then
      call addfld (dmetendnam(       1),(/ 'lev' /), 'A','kg/kg/s', &
           trim(cnst_name(       1))//' dme adjustment tendency (FV) ')
      call addfld (dmetendnam(ixcldliq),(/ 'lev' /), 'A','kg/kg/s', &
           trim(cnst_name(ixcldliq))//' dme adjustment tendency (FV) ')
      call addfld (dmetendnam(ixcldice),(/ 'lev' /), 'A','kg/kg/s', &
           trim(cnst_name(ixcldice))//' dme adjustment tendency (FV) ')
   end if

   if ( history_budget ) then
      call add_default ('PTTEND'          , history_budget_histfile_num, ' ')
      call add_default (ptendnam(       1), history_budget_histfile_num, ' ')
      call add_default (ptendnam(ixcldliq), history_budget_histfile_num, ' ')
      call add_default (ptendnam(ixcldice), history_budget_histfile_num, ' ')
      if ( dycore_is('LR') )then
         call add_default(dmetendnam(1)       , history_budget_histfile_num, ' ')
         call add_default(dmetendnam(ixcldliq), history_budget_histfile_num, ' ')
         call add_default(dmetendnam(ixcldice), history_budget_histfile_num, ' ')
      end if
      if( history_budget_histfile_num > 1 ) then
         call add_default ('DTCOND  '         , history_budget_histfile_num, ' ')
      end if
   end if

   ! outfld calls in diag_physvar_ic

   call addfld ('QCWAT&IC',(/ 'lev' /), 'I','kg/kg','q associated with cloud water'                   )
   call addfld ('TCWAT&IC',(/ 'lev' /), 'I','kg/kg','T associated with cloud water'                   )
   call addfld ('LCWAT&IC',(/ 'lev' /), 'I','kg/kg','Cloud water (ice + liq'                          )
   call addfld ('CLOUD&IC',(/ 'lev' /), 'I','fraction','Cloud fraction'                                  )
   call addfld ('CONCLD&IC',(/ 'lev' /), 'I','fraction','Convective cloud fraction'                      )
   call addfld ('TKE&IC',(/ 'ilev' /),'I','m2/s2','Turbulent Kinetic Energy'                        )
   call addfld ('CUSH&IC',horiz_only,    'I','m','Convective Scale Height'                         )
   call addfld ('KVH&IC',(/ 'ilev' /),'I','m2/s','Vertical diffusion diffusivities (heat/moisture)')
   call addfld ('KVM&IC',(/ 'ilev' /),'I','m2/s','Vertical diffusion diffusivities (momentum)'     )
   call addfld ('PBLH&IC',horiz_only,    'I','m','PBL height'                                      )
   call addfld ('TPERT&IC',horiz_only,    'I','K','Perturbation temperature (eddies in PBL)'        )
   call addfld ('QPERT&IC',horiz_only,    'I','kg/kg','Perturbation specific humidity (eddies in PBL)'  )
   call addfld ('TBOT&IC',horiz_only,    'I','K','Lowest model level temperature'                  )


   ! Initial file - Optional fields

   if (inithist_all) then
      call add_default ('CONCLD&IC  ',0, 'I')
      call add_default ('QCWAT&IC   ',0, 'I')
      call add_default ('TCWAT&IC   ',0, 'I')
      call add_default ('LCWAT&IC   ',0, 'I')
      call add_default ('PBLH&IC    ',0, 'I')
      call add_default ('TPERT&IC   ',0, 'I')
      call add_default ('QPERT&IC   ',0, 'I')
      call add_default ('CLOUD&IC   ',0, 'I')
      call add_default ('TKE&IC     ',0, 'I')
      call add_default ('CUSH&IC    ',0, 'I')
      call add_default ('KVH&IC     ',0, 'I')
      call add_default ('KVM&IC     ',0, 'I')
      call add_default ('TBOT&IC    ',0, 'I')
   end if

   ! CAM export state 
   call addfld('a2x_BCPHIWET', horiz_only, 'A', 'kg/m2/s', 'wetdep of hydrophilic black carbon')
   call addfld('a2x_BCPHIDRY', horiz_only, 'A', 'kg/m2/s', 'drydep of hydrophilic black carbon')
   call addfld('a2x_BCPHODRY', horiz_only, 'A', 'kg/m2/s', 'drydep of hydrophobic black carbon')
   call addfld('a2x_OCPHIWET', horiz_only, 'A', 'kg/m2/s', 'wetdep of hydrophilic organic carbon')
   call addfld('a2x_OCPHIDRY', horiz_only, 'A', 'kg/m2/s', 'drydep of hydrophilic organic carbon')
   call addfld('a2x_OCPHODRY', horiz_only, 'A', 'kg/m2/s', 'drydep of hydrophobic organic carbon')
   call addfld('a2x_DSTWET1', horiz_only, 'A',  'kg/m2/s', 'wetdep of dust (bin1)')
   call addfld('a2x_DSTDRY1', horiz_only, 'A',  'kg/m2/s', 'drydep of dust (bin1)')
   call addfld('a2x_DSTWET2', horiz_only, 'A',  'kg/m2/s', 'wetdep of dust (bin2)')
   call addfld('a2x_DSTDRY2', horiz_only, 'A',  'kg/m2/s', 'drydep of dust (bin2)')
   call addfld('a2x_DSTWET3', horiz_only, 'A',  'kg/m2/s', 'wetdep of dust (bin3)')
   call addfld('a2x_DSTDRY3', horiz_only, 'A',  'kg/m2/s', 'drydep of dust (bin3)')
   call addfld('a2x_DSTWET4', horiz_only, 'A',  'kg/m2/s', 'wetdep of dust (bin4)')
   call addfld('a2x_DSTDRY4', horiz_only, 'A',  'kg/m2/s', 'drydep of dust (bin4)')

   !---------------------------------------------------------
   ! CAM history fields for CAM-DOM/CAM-CSIM 
   !---------------------------------------------------------

   ! CAM-DOM history fields
#ifdef COUP_DOM
   call addfld ('TSOCN&IC',horiz_only,    'I','m','Ocean tempertare')
   call add_default ('TSOCN&IC   ',0, 'I')
#endif

  ! CAM-CSIM history fields

  do k=1,plevmx
     call addfld (tsnam(k),horiz_only,'A','K',tsnam(k)//' subsoil temperature')
  end do
  call addfld ('SICTHK'   ,horiz_only,'A','m','Sea ice thickness')
  call addfld ('TSICE'   ,horiz_only,'A','K','Ice temperature')
  do k = 1,plevmx
     call addfld (trim(tsnam(k))//'&IC',horiz_only,'I','K',tsnam(k)//' subsoil temperature')
  end do
  call addfld ('SICTHK&IC',horiz_only,'I','m','Sea ice thickness'                      )
  call addfld ('TSICE&IC',horiz_only,'I','K','Ice temperature'                        )
  call addfld ('SNOWHICE&IC',horiz_only,'I','m','Water equivalent snow depth'            )
  call addfld ('ICEFRAC&IC',horiz_only,'I','fraction','Fraction of sfc area covered by sea-ice')
  call addfld ('TSICERAD&IC',horiz_only,'I','K','Radiatively equivalent ice temperature' )
  do k = 1,plevmx
     call add_default(trim(tsnam(k))//'&IC',0, 'I')
  end do
  call add_default ('SICTHK&IC  ',0, 'I')
  call add_default ('TSICE&IC   ',0, 'I')
  call add_default ('SNOWHICE&IC',0, 'I')
  call add_default ('ICEFRAC&IC ',0, 'I')
  if (inithist_all) then
     call add_default ('TSICERAD&IC',0, 'I')
  end if

  !---------------------------------------------------------
  ! WACCM diagnostic history fields 
  !---------------------------------------------------------

  if (chem_is('waccm_ghg') .or. chem_is('waccm_mozart') .or. chem_is('waccm_mozart_mam3')) then

    ! create history variables for fourier coefficients of the diurnal 
    ! and semidiurnal tide in T, U, V, and Z3

    call tidal_diag_init() 

  endif

  qcwat_idx  = pbuf_get_index('QCWAT',ierr)
  tcwat_idx  = pbuf_get_index('TCWAT',ierr)
  lcwat_idx  = pbuf_get_index('LCWAT',ierr)
  cld_idx    = pbuf_get_index('CLD')
  concld_idx = pbuf_get_index('CONCLD')
  
  tke_idx  = pbuf_get_index('tke')
  kvm_idx  = pbuf_get_index('kvm')
  kvh_idx  = pbuf_get_index('kvh')
  cush_idx = pbuf_get_index('cush')
  
  pblh_idx  = pbuf_get_index('pblh')
  tpert_idx = pbuf_get_index('tpert')
  qpert_idx = pbuf_get_index('qpert',ierr)

  prec_dp_idx  = pbuf_get_index('PREC_DP') 
  snow_dp_idx  = pbuf_get_index('SNOW_DP') 
  prec_sh_idx  = pbuf_get_index('PREC_SH')
  snow_sh_idx  = pbuf_get_index('SNOW_SH')
  prec_sed_idx = pbuf_get_index('PREC_SED')
  snow_sed_idx = pbuf_get_index('SNOW_SED')
  prec_pcw_idx = pbuf_get_index('PREC_PCW')
  snow_pcw_idx = pbuf_get_index('SNOW_PCW')

end subroutine diag_init

!===============================================================================

subroutine diag_allocate()
   use infnan, only: nan, assignment(=)

   ! Allocate memory for module variables.
   ! Done at the begining of a physics step at same point as the pbuf allocate for
   ! variables with "physpkg" scope.

   ! Local variables
   character(len=*), parameter :: sub = 'diag_allocate'
   integer :: i, istat

   allocate(dtcond(pcols,pver,begchunk:endchunk), stat=istat)
   if ( istat /= 0 ) call endrun (sub//': ERROR: allocate failed')
   dtcond = nan

   if (dqcond_num > 0) then
      allocate(dqcond(dqcond_num))
      do i = 1, dqcond_num
         allocate(dqcond(i)%cnst(pcols,pver,begchunk:endchunk), stat=istat)
         if ( istat /= 0 ) call endrun (sub//': ERROR: allocate failed')
         dqcond(i)%cnst = nan
      end do
   end if

end subroutine diag_allocate

!===============================================================================

subroutine diag_deallocate()

! Deallocate memory for module variables.
! Done at the end of a physics step at same point as the pbuf deallocate for
! variables with "physpkg" scope.

! Local variables
   character(len=*), parameter :: sub = 'diag_deallocate'
   integer :: i, istat

   deallocate(dtcond, stat=istat)
   if ( istat /= 0 ) call endrun (sub//': ERROR: deallocate failed')

   if (dqcond_num > 0) then
      do i = 1, dqcond_num
         deallocate(dqcond(i)%cnst, stat=istat)
         if ( istat /= 0 ) call endrun (sub//': ERROR: deallocate failed')
      end do
      deallocate(dqcond, stat=istat)
      if ( istat /= 0 ) call endrun (sub//': ERROR: deallocate failed')
   end if

end subroutine diag_deallocate
!===============================================================================

subroutine diag_conv_tend_ini(state,pbuf)

! Initialize convective tendency calcs.

   
! Argument:

   type(physics_state), intent(in) :: state
   
   type(physics_buffer_desc), pointer :: pbuf(:)

! Local variables:

   integer :: i, k, m, lchnk, ncol
   real(r8), pointer, dimension(:,:) :: t_ttend

   lchnk = state%lchnk
   ncol  = state%ncol

   do k = 1, pver
      do i = 1, ncol
         dtcond(i,k,lchnk) = state%s(i,k)
      end do
   end do

   do m = 1, dqcond_num
      do k = 1, pver
         do i = 1, ncol
            dqcond(m)%cnst(i,k,lchnk) = state%q(i,k,m)
         end do
      end do
   end do

   !! initialize to pbuf T_TTEND to temperature at first timestep
   if (is_first_step()) then
      do m = 1, dyn_time_lvls
         call pbuf_get_field(pbuf, t_ttend_idx, t_ttend, start=(/1,1,m/), kount=(/pcols,pver,1/))
         t_ttend(:ncol,:) = state%t(:ncol,:)
      end do
   end if

end subroutine diag_conv_tend_ini
!===============================================================================

  subroutine diag_phys_writeout(state, psl)

!----------------------------------------------------------------------- 
! 
! Purpose: record dynamics variables on physics grid
!
!-----------------------------------------------------------------------
    use physconst,          only: gravit, rga, rair, cpair, latvap, rearth, pi, cappa
    use time_manager,       only: get_nstep
    use interpolate_data,   only: vertinterp
    use constituent_burden, only: constituent_burden_comp
    use cam_control_mod,    only: moist_physics
    use co2_cycle,          only: c_i, co2_transport

    use tidal_diag,         only: tidal_diag_write
!-----------------------------------------------------------------------
!
! Arguments
!
   type(physics_state), intent(inout) :: state
   real(r8), optional , intent(out)   :: psl(pcols) 
!
!---------------------------Local workspace-----------------------------
!
    real(r8) ftem(pcols,pver) ! temporary workspace
    real(r8) ftem1(pcols,pver) ! another temporary workspace
    real(r8) ftem2(pcols,pver) ! another temporary workspace
    real(r8) psl_tmp(pcols)   ! Sea Level Pressure
    real(r8) z3(pcols,pver)   ! geo-potential height
    real(r8) p_surf(pcols)    ! data interpolated to a pressure surface
    real(r8) p_surf_t1(pcols)    ! data interpolated to a pressure surface
    real(r8) p_surf_t2(pcols)    ! data interpolated to a pressure surface
    real(r8) p_surf_q1(pcols)    ! data interpolated to a pressure surface    
    real(r8) p_surf_q2(pcols)    ! data interpolated to a pressure surface        
    real(r8) tem2(pcols,pver) ! temporary workspace
    real(r8) timestep(pcols)  ! used for outfld call
    real(r8) esl(pcols,pver)   ! saturation vapor pressures 
    real(r8) esi(pcols,pver)   ! 
    real(r8) dlon(pcols)      ! width of grid cell (meters)
    integer  plon             ! number of longitudes

    integer i, k, m, lchnk, ncol, nstep
!
!-----------------------------------------------------------------------
!
    lchnk = state%lchnk
    ncol  = state%ncol

    ! Output NSTEP for debugging
    nstep = get_nstep()
    timestep(:ncol) = nstep
    call outfld ('NSTEP   ',timestep, pcols, lchnk)

    call outfld('T       ',state%t , pcols   ,lchnk   )
    call outfld('PS      ',state%ps, pcols   ,lchnk   )
    call outfld('U       ',state%u , pcols   ,lchnk   )
    call outfld('V       ',state%v , pcols   ,lchnk   )
    do m=1,pcnst
       if ( cnst_cam_outfld(m) ) then
          call outfld(cnst_name(m),state%q(1,1,m),pcols ,lchnk )
       end if
    end do

    if (co2_transport()) then
       do m = 1,4
          call outfld(trim(cnst_name(c_i(m)))//'_BOT', state%q(1,pver,c_i(m)), pcols, lchnk)
       end do
    end if

    ! column burdens of all constituents except water vapor
    call constituent_burden_comp(state)

    if ( moist_physics) then
       call outfld('PDELDRY ',state%pdeldry, pcols, lchnk)
       call outfld('PSDRY',   state%psdry,   pcols, lchnk) 
    end if

    call outfld('PHIS    ',state%phis,    pcols,   lchnk     )



#if (defined BFB_CAM_SCAM_IOP )
    call outfld('phis    ',state%phis,    pcols,   lchnk     )
#endif

!
! Add height of surface to midpoint height above surface 
!
    do k = 1, pver
       z3(:ncol,k) = state%zm(:ncol,k) + state%phis(:ncol)*rga
    end do
    call outfld('Z3      ',z3,pcols,lchnk)
!           
! Output Z3 on pressure surfaces
!
    if (hist_fld_active('Z1000')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 100000._r8, z3, p_surf)
       call outfld('Z1000    ', p_surf, pcols, lchnk)
    end if
    if (hist_fld_active('Z700')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 70000._r8, z3, p_surf)
       call outfld('Z700    ', p_surf, pcols, lchnk)
    end if
    if (hist_fld_active('Z500')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 50000._r8, z3, p_surf)
       call outfld('Z500    ', p_surf, pcols, lchnk)
    end if
    if (hist_fld_active('Z300')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 30000._r8, z3, p_surf)
       call outfld('Z300    ', p_surf, pcols, lchnk)
    end if
    if (hist_fld_active('Z200')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 20000._r8, z3, p_surf)
       call outfld('Z200    ', p_surf, pcols, lchnk)
    end if
    if (hist_fld_active('Z100')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 10000._r8, z3, p_surf)
       call outfld('Z100    ', p_surf, pcols, lchnk)
    end if
    if (hist_fld_active('Z050')) then
       call vertinterp(ncol, pcols, pver, state%pmid,  5000._r8, z3, p_surf)
       call outfld('Z050    ', p_surf, pcols, lchnk)
    end if
!
! Quadratic height fiels Z3*Z3
!
    ftem(:ncol,:) = z3(:ncol,:)*z3(:ncol,:)
    call outfld('ZZ      ',ftem,pcols,lchnk)

    ftem(:ncol,:) = z3(:ncol,:)*state%v(:ncol,:)*gravit
    call outfld('VZ      ',ftem,  pcols,lchnk)
!
! Meridional advection fields
!
    ftem(:ncol,:) = state%v(:ncol,:)*state%t(:ncol,:)
    call outfld ('VT      ',ftem    ,pcols   ,lchnk     )

    ftem(:ncol,:) = state%v(:ncol,:)*state%q(:ncol,:,1)
    call outfld ('VQ      ',ftem    ,pcols   ,lchnk     )
    if(prog_modal_aero) then !Only for prognostic aerosols
       ftem(:ncol,:) = state%v(:ncol,:)*state%q(:ncol,:,14)
       call outfld ('Vbc_a1  ',ftem    ,pcols   ,lchnk     )
       
       ftem(:ncol,:) = state%v(:ncol,:)*state%q(:ncol,:,15)
       call outfld ('Vdst_a1 ',ftem    ,pcols   ,lchnk     )
       
       ftem(:ncol,:) = state%v(:ncol,:)*state%q(:ncol,:,22)
       call outfld ('Vdst_a3 ',ftem    ,pcols   ,lchnk     )
       
       ftem(:ncol,:) = state%v(:ncol,:)*state%q(:ncol,:,16)
       call outfld ('Vncl_a1 ',ftem    ,pcols   ,lchnk     )
       
       ftem(:ncol,:) = state%v(:ncol,:)*state%q(:ncol,:,20)
       call outfld ('Vncl_a2 ',ftem    ,pcols   ,lchnk     )
       
       ftem(:ncol,:) = state%v(:ncol,:)*state%q(:ncol,:,23)
       call outfld ('Vncl_a3 ',ftem    ,pcols   ,lchnk     )
       
       ftem(:ncol,:) = state%v(:ncol,:)*state%q(:ncol,:,11)
       call outfld ('Vso4_a1 ',ftem    ,pcols   ,lchnk     )
       
       ftem(:ncol,:) = state%v(:ncol,:)*state%q(:ncol,:,18)
       call outfld ('Vso4_a2 ',ftem    ,pcols   ,lchnk     )
       
       ftem(:ncol,:) = state%v(:ncol,:)*state%q(:ncol,:,24)
       call outfld ('Vso4_a3 ',ftem    ,pcols   ,lchnk     )
       
       ftem(:ncol,:) = state%v(:ncol,:)*state%q(:ncol,:,13)
       call outfld ('Vsoa_a1 ',ftem    ,pcols   ,lchnk     )
       
       ftem(:ncol,:) = state%v(:ncol,:)*state%q(:ncol,:,19)
       call outfld ('Vsoa_a2 ',ftem    ,pcols   ,lchnk     )
       
       ftem(:ncol,:) = state%v(:ncol,:)*state%q(:ncol,:,12)
       call outfld ('Vpom_a1 ',ftem    ,pcols   ,lchnk     )
    endif

    ftem(:ncol,:) = state%q(:ncol,:,1)*state%q(:ncol,:,1)
    call outfld ('QQ      ',ftem    ,pcols   ,lchnk     )

    if(prog_modal_aero) then !Only for prognostic aerosols
       ftem(:ncol,:) = state%q(:ncol,:,14)*state%q(:ncol,:,14)
       call outfld ('bc_a1_2 ',ftem    ,pcols   ,lchnk     )
       
       ftem(:ncol,:) = state%q(:ncol,:,15)*state%q(:ncol,:,15)
       call outfld ('dst_a1_2',ftem    ,pcols   ,lchnk     )
       
       ftem(:ncol,:) = state%q(:ncol,:,22)*state%q(:ncol,:,22)
       call outfld ('dst_a3_2',ftem    ,pcols   ,lchnk     )
       
       ftem(:ncol,:) = state%q(:ncol,:,16)*state%q(:ncol,:,16)
       call outfld ('ncl_a1_2',ftem    ,pcols   ,lchnk     )
       
       ftem(:ncol,:) = state%q(:ncol,:,20)*state%q(:ncol,:,20)
       call outfld ('ncl_a2_2',ftem    ,pcols   ,lchnk     )
       
       ftem(:ncol,:) = state%q(:ncol,:,23)*state%q(:ncol,:,23)
       call outfld ('ncl_a3_2',ftem    ,pcols   ,lchnk     )
       
       ftem(:ncol,:) = state%q(:ncol,:,11)*state%q(:ncol,:,11)
       call outfld ('so4_a1_2',ftem    ,pcols   ,lchnk     )
       
       ftem(:ncol,:) = state%q(:ncol,:,18)*state%q(:ncol,:,18)
       call outfld ('so4_a2_2',ftem    ,pcols   ,lchnk     )
       
       ftem(:ncol,:) = state%q(:ncol,:,24)*state%q(:ncol,:,24)
       call outfld ('so4_a3_2',ftem    ,pcols   ,lchnk     )
       
       ftem(:ncol,:) = state%q(:ncol,:,13)*state%q(:ncol,:,13)
       call outfld ('soa_a1_2',ftem    ,pcols   ,lchnk     )
       
       ftem(:ncol,:) = state%q(:ncol,:,19)*state%q(:ncol,:,19)
       call outfld ('soa_a2_2',ftem    ,pcols   ,lchnk     )
       
       ftem(:ncol,:) = state%q(:ncol,:,12)*state%q(:ncol,:,12)
       call outfld ('pom_a1_2',ftem    ,pcols   ,lchnk     )
    endif

    ftem(:ncol,:) = state%v(:ncol,:)**2
    call outfld ('VV      ',ftem    ,pcols   ,lchnk     )

    ftem(:ncol,:) = state%v(:ncol,:) * state%u(:ncol,:)
    call outfld ('VU      ',ftem    ,pcols   ,lchnk     )

! zonal advection

    ftem(:ncol,:) = state%u(:ncol,:)**2
    call outfld ('UU      ',ftem    ,pcols   ,lchnk     )

! Wind speed
    ftem(:ncol,:) = sqrt( state%u(:ncol,:)**2 + state%v(:ncol,:)**2)
    call outfld ('WSPEED  ',ftem    ,pcols   ,lchnk     )
    call outfld ('WSPDSRFMX',ftem(:,pver)   ,pcols   ,lchnk     )
    call outfld ('WSPDSRFAV',ftem(:,pver)   ,pcols   ,lchnk     ) 

! Vertical velocity and advection

    if (single_column) then
       call outfld('OMEGA   ',wfld,    pcols,   lchnk     )
    else
       call outfld('OMEGA   ',state%omega,    pcols,   lchnk     )
    endif

#if (defined BFB_CAM_SCAM_IOP )
    call outfld('omega   ',state%omega,    pcols,   lchnk     )
#endif

    ftem(:ncol,:) = state%omega(:ncol,:)*state%t(:ncol,:)
    call outfld('OMEGAT  ',ftem,    pcols,   lchnk     )
    ftem(:ncol,:) = state%omega(:ncol,:)*state%u(:ncol,:)
    call outfld('OMEGAU  ',ftem,    pcols,   lchnk     )
    ftem(:ncol,:) = state%omega(:ncol,:)*state%v(:ncol,:)
    call outfld('OMEGAV  ',ftem,    pcols,   lchnk     )
    ftem(:ncol,:) = state%omega(:ncol,:)*state%q(:ncol,:,1)
    call outfld('OMEGAQ  ',ftem,    pcols,   lchnk     )
    ftem(:ncol,:) = state%omega(:ncol,:)*state%omega(:ncol,:)
    call outfld('OMGAOMGA',ftem,    pcols,   lchnk     )
!
! Output omega at 850 and 500 mb pressure levels
!
    if (hist_fld_active('OMEGA850')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 85000._r8, state%omega, p_surf)
       call outfld('OMEGA850', p_surf, pcols, lchnk)
    end if
    if (hist_fld_active('OMEGA500')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 50000._r8, state%omega, p_surf)
       call outfld('OMEGA500', p_surf, pcols, lchnk)
    end if
!     
! Mass of q, by layer and vertically integrated
!
    ftem(:ncol,:) = state%q(:ncol,:,1) * state%pdel(:ncol,:) * rga
    call outfld ('MQ      ',ftem    ,pcols   ,lchnk     )

    do k=2,pver
       ftem(:ncol,1) = ftem(:ncol,1) + ftem(:ncol,k)
    end do
    call outfld ('TMQ     ',ftem, pcols   ,lchnk     )
!
! Mass of vertically integrated q flux
!
    ftem(:ncol,:) = state%u(:ncol,:)*state%q(:ncol,:,1)*state%pdel(:ncol,:)*rga
    do k=2,pver
       ftem(:ncol,1) = ftem(:ncol,1) + ftem(:ncol,k)
    end do
    call outfld ('TUQ     ',ftem, pcols   ,lchnk     )

    ftem(:ncol,:) = state%v(:ncol,:)*state%q(:ncol,:,1)*state%pdel(:ncol,:)*rga
    do k=2,pver
       ftem(:ncol,1) = ftem(:ncol,1) + ftem(:ncol,k)
    end do
    call outfld ('TVQ     ',ftem, pcols   ,lchnk     )

    if (moist_physics) then

       ! Relative humidity
       if (hist_fld_active('RELHUM')) then
          call qsat(state%t(:ncol,:), state%pmid(:ncol,:), &
               tem2(:ncol,:), ftem(:ncol,:))
          ftem(:ncol,:) = state%q(:ncol,:,1)/ftem(:ncol,:)*100._r8
          call outfld ('RELHUM  ',ftem    ,pcols   ,lchnk     )
       end if

       if (hist_fld_active('RHW') .or. hist_fld_active('RHI') .or. hist_fld_active('RHCFMIP') ) then
	  
          ! RH w.r.t liquid (water)
          call qsat_water (state%t(:ncol,:), state%pmid(:ncol,:), &
               esl(:ncol,:), ftem(:ncol,:))
          ftem(:ncol,:) = state%q(:ncol,:,1)/ftem(:ncol,:)*100._r8
          call outfld ('RHW  ',ftem    ,pcols   ,lchnk     )

          ! Convert to RHI (ice)
          do i=1,ncol
             do k=1,pver
                esi(i,k)=svp_ice(state%t(i,k))
                ftem1(i,k)=ftem(i,k)*esl(i,k)/esi(i,k)
             end do
          end do
          call outfld ('RHI  ',ftem1    ,pcols   ,lchnk     )

	  ! use temperature to decide if you populate with ftem (liquid, above 0 C) or ftem1 (ice, below 0 C)

	  ftem2(:ncol,:)=ftem(:ncol,:)

          do i=1,ncol
             do k=1,pver
		if (state%t(i,k) .gt. 273) then
                   ftem2(i,k)=ftem(i,k)  !!wrt water
 		else
                   ftem2(i,k)=ftem1(i,k) !!wrt ice
		end if
             end do
          end do
          
          call outfld ('RHCFMIP  ',ftem2    ,pcols   ,lchnk     )

       end if

    end if
!
! Sea level pressure
!
    if (present(psl) .or. hist_fld_active('PSL')) then
       call cpslec (ncol, state%pmid, state%phis, state%ps, state%t,psl_tmp, gravit, rair) 
       call outfld ('PSL     ',psl_tmp  ,pcols, lchnk     )
       if (present(psl)) then	
          psl(:ncol) = psl_tmp(:ncol)
       end if
    end if
!
! Output T,q,u,v fields on pressure surfaces
!
    if (hist_fld_active('T850')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 85000._r8, state%t, p_surf)
       call outfld('T850    ', p_surf, pcols, lchnk )
    end if
    if (hist_fld_active('T500')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 50000._r8, state%t, p_surf)
       call outfld('T500    ', p_surf, pcols, lchnk )
    end if
    if (hist_fld_active('T300')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 30000._r8, state%t, p_surf)
       call outfld('T300    ', p_surf, pcols, lchnk )
    end if
    if (hist_fld_active('T200')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 20000._r8, state%t, p_surf)
       call outfld('T200    ', p_surf, pcols, lchnk )
    end if
    if (hist_fld_active('Q850')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 85000._r8, state%q(1,1,1), p_surf)
       call outfld('Q850    ', p_surf, pcols, lchnk )
    end if
    if (hist_fld_active('Q200')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 20000._r8, state%q(1,1,1), p_surf)
       call outfld('Q200    ', p_surf, pcols, lchnk )
    end if
    if (hist_fld_active('U850')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 85000._r8, state%u, p_surf)
       call outfld('U850    ', p_surf, pcols, lchnk )
    end if
    if (hist_fld_active('U250')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 25000._r8, state%u, p_surf)
       call outfld('U250    ', p_surf, pcols, lchnk )
    end if
    if (hist_fld_active('U200')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 20000._r8, state%u, p_surf)
       call outfld('U200    ', p_surf, pcols, lchnk )
    end if
    if (hist_fld_active('U010')) then
       call vertinterp(ncol, pcols, pver, state%pmid,  1000._r8, state%u, p_surf)
       call outfld('U010    ', p_surf, pcols, lchnk )
    end if
    if (hist_fld_active('V850')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 85000._r8, state%v, p_surf)
       call outfld('V850    ', p_surf, pcols, lchnk )
    end if
    if (hist_fld_active('V250')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 25000._r8, state%v, p_surf)
       call outfld('V250    ', p_surf, pcols, lchnk )
    end if
    if (hist_fld_active('V200')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 20000._r8, state%v, p_surf)
       call outfld('V200    ', p_surf, pcols, lchnk )
    end if

    ftem(:ncol,:) = state%t(:ncol,:)*state%t(:ncol,:)
    call outfld('TT      ',ftem    ,pcols   ,lchnk   )
!
! Output U, V, T, Q, P and Z at bottom level
!
    call outfld ('UBOT    ', state%u(1,pver)  ,  pcols, lchnk)
    call outfld ('VBOT    ', state%v(1,pver)  ,  pcols, lchnk)
    call outfld ('QBOT    ', state%q(1,pver,1),  pcols, lchnk)
    call outfld ('ZBOT    ', state%zm(1,pver) , pcols, lchnk)

! Total energy of the atmospheric column for atmospheric heat storage calculations

    !! temporary variable to get surface geopotential in dimensions of (ncol,pver)
    do k=1,pver
      ftem1(:ncol,k)=state%phis(:ncol)  !! surface geopotential in units (m2/s2)
    end do

    !! calculate sum of sensible, kinetic, latent, and surface geopotential energy
    !! E=CpT+PHIS+Lv*q+(0.5)*(u^2+v^2)
    ftem(:ncol,:) = (cpair*state%t(:ncol,:) +  ftem1(:ncol,:) + latvap*state%q(:ncol,:,1) + &
         0.5_r8*(state%u(:ncol,:)**2+state%v(:ncol,:)**2))*(state%pdel(:ncol,:)/gravit)
    !! vertically integrate
    do k=2,pver       
	ftem(:ncol,1) = ftem(:ncol,1) + ftem(:ncol,k)
    end do
    call outfld ('ATMEINT   ',ftem(:ncol,1)  ,pcols   ,lchnk     )

!! Boundary layer atmospheric stability, temperature, water vapor diagnostics

    if (hist_fld_active('T1000')      .or. &
        hist_fld_active('T9251000')   .or. & 
        hist_fld_active('TH9251000')  .or. &
        hist_fld_active('THE9251000') .or. &
        hist_fld_active('T8501000')   .or. &
        hist_fld_active('TH8501000')  .or. &
        hist_fld_active('THE8501000') .or. &
        hist_fld_active('T7001000')   .or. &
        hist_fld_active('TH7001000')  .or. &
        hist_fld_active('THE7001000')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 100000._r8, state%t, p_surf_t1)
    end if

    if (hist_fld_active('T925')       .or. &
        hist_fld_active('T9251000')   .or. & 
        hist_fld_active('TH9251000')  .or. &
        hist_fld_active('THE9251000')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 92500._r8, state%t, p_surf_t2)
    end if

    if (hist_fld_active('Q1000')      .or. &
        hist_fld_active('THE9251000') .or. &
        hist_fld_active('THE8501000') .or. &
        hist_fld_active('THE7001000')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 100000._r8, state%q(1,1,1), p_surf_q1)
    end if

    if (hist_fld_active('Q925')       .or. &
        hist_fld_active('THE9251000')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 92500._r8, state%q(1,1,1), p_surf_q2)
    end if

    !!! at 1000 mb and 925 mb
    if (hist_fld_active('T1000')) then
       call outfld('T1000    ', p_surf_t1, pcols, lchnk )
    end if

    if (hist_fld_active('T925')) then
       call outfld('T925    ', p_surf_t2, pcols, lchnk )
    end if

    if (hist_fld_active('Q1000')) then
       call outfld('Q1000    ', p_surf_q1, pcols, lchnk ) 
    end if

    if (hist_fld_active('Q925')) then
       call outfld('Q925    ', p_surf_q2, pcols, lchnk )  
    end if

    if (hist_fld_active('T9251000')) then
       p_surf = p_surf_t2-p_surf_t1  
       call outfld('T9251000    ', p_surf, pcols, lchnk ) 
    end if

    if (hist_fld_active('TH9251000')) then
       p_surf = (p_surf_t2*(1000.0_r8/925.0_r8)**cappa)-(p_surf_t1*(1.0_r8)**cappa)   
       call outfld('TH9251000    ', p_surf, pcols, lchnk )    
    end if

    if (hist_fld_active('THE9251000')) then
       p_surf = (p_surf_t2*(1000.0_r8/925.0_r8)**cappa)*exp((2500000.0_r8*p_surf_q2)/(1004.0_r8*p_surf_t2))- &
            (p_surf_t1*(1.0_r8)**cappa)*exp((2500000.0_r8*p_surf_q1)/(1004.0_r8*p_surf_t1))  
       call outfld('THE9251000    ', p_surf, pcols, lchnk )
    end if

    if (hist_fld_active('T8501000')  .or. &
        hist_fld_active('TH8501000') .or. &
        hist_fld_active('THE8501000')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 85000._r8, state%t, p_surf_t2)
    end if

    !!! at 1000 mb and 850 mb
    if (hist_fld_active('T8501000')) then
       p_surf = p_surf_t2-p_surf_t1  
       call outfld('T8501000    ', p_surf, pcols, lchnk ) 
    end if

    if (hist_fld_active('TH8501000')) then
       p_surf = (p_surf_t2*(1000.0_r8/850.0_r8)**cappa)-(p_surf_t1*(1.0_r8)**cappa)   
       call outfld('TH8501000    ', p_surf, pcols, lchnk )   
    end if

    if (hist_fld_active('THE8501000')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 85000._r8, state%q(1,1,1), p_surf_q2)
       p_surf = (p_surf_t2*(1000.0_r8/850.0_r8)**cappa)*exp((2500000.0_r8*p_surf_q2)/(1004.0_r8*p_surf_t2))- &
            (p_surf_t1*(1.0_r8)**cappa)*exp((2500000.0_r8*p_surf_q1)/(1004.0_r8*p_surf_t1))  
       call outfld('THE8501000    ', p_surf, pcols, lchnk ) 
    end if

    if (hist_fld_active('T7001000')  .or. &
        hist_fld_active('TH7001000') .or. &
        hist_fld_active('T700') .or. &
        hist_fld_active('THE7001000')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 70000._r8, state%t, p_surf_t2)
    end if

   !!! at 700 mb
    if (hist_fld_active('T700')) then
       call outfld('T700    ', p_surf_t2, pcols, lchnk )
    end if

    !!! at 1000 mb and 700 mb
    if (hist_fld_active('T7001000')) then
       p_surf = p_surf_t2-p_surf_t1
       call outfld('T7001000    ', p_surf, pcols, lchnk )
    end if

    if (hist_fld_active('TH7001000')) then
       p_surf = (p_surf_t2*(1000.0_r8/700.0_r8)**cappa)-(p_surf_t1*(1.0_r8)**cappa)
       call outfld('TH7001000    ', p_surf, pcols, lchnk )
    end if

    if (hist_fld_active('THE7001000')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 70000._r8, state%q(1,1,1), p_surf_q2)
       p_surf = (p_surf_t2*(1000.0_r8/700.0_r8)**cappa)*exp((2500000.0_r8*p_surf_q2)/(1004.0_r8*p_surf_t2))- &
            (p_surf_t1*(1.0_r8)**cappa)*exp((2500000.0_r8*p_surf_q1)/(1004.0_r8*p_surf_t1))
       call outfld('THE7001000    ', p_surf, pcols, lchnk )
    end if

    if (hist_fld_active('T010')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 1000._r8, state%t, p_surf)
       call outfld('T010           ', p_surf, pcols, lchnk )
    end if


  !---------------------------------------------------------
  ! WACCM tidal diagnostics 
  !---------------------------------------------------------

  if (chem_is('waccm_ghg') .or. chem_is('waccm_mozart') .or. chem_is('waccm_mozart_mam3')) then

    call tidal_diag_write(state) 

  endif

    return
  end subroutine diag_phys_writeout
!===============================================================================

subroutine diag_conv(state, ztodt, pbuf)

!----------------------------------------------------------------------- 
! 
! Output diagnostics associated with all convective processes.
!
!-----------------------------------------------------------------------
   use physconst,     only: cpair
   use tidal_diag,    only: get_tidal_coeffs

! Arguments:

   real(r8),            intent(in) :: ztodt   ! timestep for computing physics tendencies
   type(physics_state), intent(in) :: state
   type(physics_buffer_desc), pointer :: pbuf(:)

! convective precipitation variables
   real(r8), pointer :: prec_dp(:)                 ! total precipitation   from ZM convection
   real(r8), pointer :: snow_dp(:)                 ! snow from ZM   convection
   real(r8), pointer :: prec_sh(:)                 ! total precipitation   from Hack convection
   real(r8), pointer :: snow_sh(:)                 ! snow from   Hack   convection
   real(r8), pointer :: prec_sed(:)                ! total precipitation   from ZM convection
   real(r8), pointer :: snow_sed(:)                ! snow from ZM   convection
   real(r8), pointer :: prec_pcw(:)                ! total precipitation   from Hack convection
   real(r8), pointer :: snow_pcw(:)                ! snow from Hack   convection

! Local variables:
   
   integer :: i, k, m, lchnk, ncol

   real(r8) :: rtdt

   real(r8):: precc(pcols)                ! convective precip rate
   real(r8):: precl(pcols)                ! stratiform precip rate
   real(r8):: snowc(pcols)                ! convective snow rate
   real(r8):: snowl(pcols)                ! stratiform snow rate
   real(r8):: prect(pcols)                ! total (conv+large scale) precip rate
   real(r8) :: dcoef(4)                   ! for tidal component of T tend

   lchnk = state%lchnk
   ncol  = state%ncol

   rtdt = 1._r8/ztodt

   call pbuf_get_field(pbuf, prec_dp_idx, prec_dp)
   call pbuf_get_field(pbuf, snow_dp_idx, snow_dp)
   call pbuf_get_field(pbuf, prec_sh_idx, prec_sh)
   call pbuf_get_field(pbuf, snow_sh_idx, snow_sh)
   call pbuf_get_field(pbuf, prec_sed_idx, prec_sed)
   call pbuf_get_field(pbuf, snow_sed_idx, snow_sed)
   call pbuf_get_field(pbuf, prec_pcw_idx, prec_pcw)
   call pbuf_get_field(pbuf, snow_pcw_idx, snow_pcw)

! Precipitation rates (multi-process)
   precc(:ncol) = prec_dp(:ncol)  + prec_sh(:ncol)
   precl(:ncol) = prec_sed(:ncol) + prec_pcw(:ncol)
   snowc(:ncol) = snow_dp(:ncol)  + snow_sh(:ncol)
   snowl(:ncol) = snow_sed(:ncol) + snow_pcw(:ncol)
   prect(:ncol) = precc(:ncol)    + precl(:ncol)

   call outfld('PRECC   ', precc, pcols, lchnk )
   call outfld('PRECL   ', precl, pcols, lchnk )
   call outfld('PREC_PCW', prec_pcw,pcols   ,lchnk )
   call outfld('PREC_zmc', prec_dp ,pcols   ,lchnk )
   call outfld('PRECSC  ', snowc, pcols, lchnk )
   call outfld('PRECSL  ', snowl, pcols, lchnk )
   call outfld('PRECT   ', prect, pcols, lchnk )
   call outfld('PRECTMX ', prect, pcols, lchnk )

   call outfld('PRECLav ', precl, pcols, lchnk )
   call outfld('PRECCav ', precc, pcols, lchnk )

#if ( defined BFB_CAM_SCAM_IOP )
   call outfld('Prec   ' , prect, pcols, lchnk )
#endif

   ! Total convection tendencies.

   do k = 1, pver
      do i = 1, ncol
         dtcond(i,k,lchnk) = (state%s(i,k) - dtcond(i,k,lchnk))*rtdt / cpair
      end do
   end do
   call outfld('DTCOND  ', dtcond(:,:,lchnk), pcols, lchnk)

   ! output tidal coefficients
   call get_tidal_coeffs( dcoef )
   call outfld( 'DTCOND_24_SIN', dtcond(:ncol,:,lchnk)*dcoef(1), ncol, lchnk )
   call outfld( 'DTCOND_24_COS', dtcond(:ncol,:,lchnk)*dcoef(2), ncol, lchnk )
   call outfld( 'DTCOND_12_SIN', dtcond(:ncol,:,lchnk)*dcoef(3), ncol, lchnk )
   call outfld( 'DTCOND_12_COS', dtcond(:ncol,:,lchnk)*dcoef(4), ncol, lchnk )

   do m = 1, dqcond_num
      if ( cnst_cam_outfld(m) ) then
         do k = 1, pver
            do i = 1, ncol
               dqcond(m)%cnst(i,k,lchnk) = (state%q(i,k,m) - dqcond(m)%cnst(i,k,lchnk))*rtdt
            end do
         end do
         call outfld(dcconnam(m), dqcond(m)%cnst(:,:,lchnk), pcols, lchnk)
      end if
   end do

end subroutine diag_conv

!===============================================================================

subroutine diag_surf (cam_in, cam_out, ps, trefmxav, trefmnav )

!----------------------------------------------------------------------- 
! 
! Purpose: record surface diagnostics
!
!-----------------------------------------------------------------------

   use time_manager,     only: is_end_curr_day
   use co2_cycle,        only: c_i, co2_transport
   use constituents,     only: sflxnam

!-----------------------------------------------------------------------
!
! Input arguments
!
    type(cam_in_t),  intent(in) :: cam_in
    type(cam_out_t), intent(in) :: cam_out

    real(r8), intent(inout) :: trefmnav(pcols) ! daily minimum tref  
    real(r8), intent(inout) :: trefmxav(pcols) ! daily maximum tref

    real(r8), intent(in)    :: ps(pcols)       ! Surface pressure.
!
!---------------------------Local workspace-----------------------------
!
    integer :: i, k, m      ! indexes
    integer :: lchnk        ! chunk identifier
    integer :: ncol         ! longitude dimension
    real(r8) tem2(pcols)    ! temporary workspace
    real(r8) ftem(pcols)    ! temporary workspace
!
!-----------------------------------------------------------------------
!
    lchnk = cam_in%lchnk
    ncol  = cam_in%ncol

    call outfld('SHFLX',    cam_in%shf,       pcols, lchnk)
    call outfld('LHFLX',    cam_in%lhf,       pcols, lchnk)
    call outfld('QFLX',     cam_in%cflx(1,1), pcols, lchnk)

    call outfld('TAUX',     cam_in%wsx,       pcols, lchnk)
    call outfld('TAUY',     cam_in%wsy,       pcols, lchnk)
    call outfld('TREFHT  ', cam_in%tref,      pcols, lchnk)
    call outfld('TREFHTMX', cam_in%tref,      pcols, lchnk)
    call outfld('TREFHTMN', cam_in%tref,      pcols, lchnk)
    call outfld('QREFHT',   cam_in%qref,      pcols, lchnk)
    call outfld('U10',      cam_in%u10,       pcols, lchnk)
! 
! Calculate and output reference height RH (RHREFHT)

   call qsat(cam_in%tref(:ncol), ps(:ncol), tem2(:ncol), ftem(:ncol))
       ftem(:ncol) = cam_in%qref(:ncol)/ftem(:ncol)*100._r8

      
    call outfld('RHREFHT',   ftem,      pcols, lchnk)


#if (defined BFB_CAM_SCAM_IOP )
    call outfld('shflx   ',cam_in%shf,   pcols,   lchnk)
    call outfld('lhflx   ',cam_in%lhf,   pcols,   lchnk)
    call outfld('trefht  ',cam_in%tref,  pcols,   lchnk)
#endif
!
! Ouput ocn and ice fractions
!
    call outfld('LANDFRAC', cam_in%landfrac, pcols, lchnk)
    call outfld('ICEFRAC',  cam_in%icefrac,  pcols, lchnk)
    call outfld('OCNFRAC',  cam_in%ocnfrac,  pcols, lchnk)
!
! Compute daily minimum and maximum of TREF
!
    do i = 1,ncol
       trefmxav(i) = max(cam_in%tref(i),trefmxav(i))
       trefmnav(i) = min(cam_in%tref(i),trefmnav(i))
    end do
    if (is_end_curr_day()) then
       call outfld('TREFMXAV', trefmxav,pcols,   lchnk     )
       call outfld('TREFMNAV', trefmnav,pcols,   lchnk     )
       trefmxav(:ncol) = -1.0e36_r8
       trefmnav(:ncol) =  1.0e36_r8
    endif

    call outfld('TBOT',     cam_out%tbot,     pcols, lchnk)
    call outfld('TS',       cam_in%ts,        pcols, lchnk)
    call outfld('TSMN',     cam_in%ts,        pcols, lchnk)
    call outfld('TSMX',     cam_in%ts,        pcols, lchnk)
    call outfld('SNOWHLND', cam_in%snowhland, pcols, lchnk)
    call outfld('SNOWHICE', cam_in%snowhice,  pcols, lchnk)
    call outfld('ASDIR',    cam_in%asdir,     pcols, lchnk)
    call outfld('ASDIF',    cam_in%asdif,     pcols, lchnk)
    call outfld('ALDIR',    cam_in%aldir,     pcols, lchnk)
    call outfld('ALDIF',    cam_in%aldif,     pcols, lchnk)
    call outfld('SST',      cam_in%sst,       pcols, lchnk)

    if (co2_transport()) then
       do m = 1,4
          call outfld(sflxnam(c_i(m)), cam_in%cflx(:,c_i(m)), pcols, lchnk)
       end do
    end if

end subroutine diag_surf

!===============================================================================

subroutine diag_export(cam_out)

!----------------------------------------------------------------------- 
! 
! Purpose: Write export state to history file
!
!-----------------------------------------------------------------------

   ! arguments
   type(cam_out_t), intent(inout) :: cam_out

   ! Local variables:
   integer :: lchnk        ! chunk identifier
   logical :: atm_dep_flux ! true ==> sending deposition fluxes to coupler.
                           ! Otherwise, set them to zero.
   !-----------------------------------------------------------------------

   lchnk = cam_out%lchnk

   call phys_getopts(atm_dep_flux_out=atm_dep_flux)

   if (.not. atm_dep_flux) then
      ! set the fluxes to zero before outfld and sending them to the
      ! coupler
      cam_out%bcphiwet = 0.0_r8
      cam_out%bcphidry = 0.0_r8
      cam_out%bcphodry = 0.0_r8
      cam_out%ocphiwet = 0.0_r8
      cam_out%ocphidry = 0.0_r8
      cam_out%ocphodry = 0.0_r8
      cam_out%dstwet1  = 0.0_r8
      cam_out%dstdry1  = 0.0_r8
      cam_out%dstwet2  = 0.0_r8
      cam_out%dstdry2  = 0.0_r8
      cam_out%dstwet3  = 0.0_r8
      cam_out%dstdry3  = 0.0_r8
      cam_out%dstwet4  = 0.0_r8
      cam_out%dstdry4  = 0.0_r8
   end if

   call outfld('a2x_BCPHIWET', cam_out%bcphiwet, pcols, lchnk)
   call outfld('a2x_BCPHIDRY', cam_out%bcphidry, pcols, lchnk)
   call outfld('a2x_BCPHODRY', cam_out%bcphodry, pcols, lchnk)
   call outfld('a2x_OCPHIWET', cam_out%ocphiwet, pcols, lchnk)
   call outfld('a2x_OCPHIDRY', cam_out%ocphidry, pcols, lchnk)
   call outfld('a2x_OCPHODRY', cam_out%ocphodry, pcols, lchnk)
   call outfld('a2x_DSTWET1',  cam_out%dstwet1,  pcols, lchnk)
   call outfld('a2x_DSTDRY1',  cam_out%dstdry1,  pcols, lchnk)
   call outfld('a2x_DSTWET2',  cam_out%dstwet2,  pcols, lchnk)
   call outfld('a2x_DSTDRY2',  cam_out%dstdry2,  pcols, lchnk)
   call outfld('a2x_DSTWET3',  cam_out%dstwet3,  pcols, lchnk)
   call outfld('a2x_DSTDRY3',  cam_out%dstdry3,  pcols, lchnk)
   call outfld('a2x_DSTWET4',  cam_out%dstwet4,  pcols, lchnk)
   call outfld('a2x_DSTDRY4',  cam_out%dstdry4,  pcols, lchnk)

end subroutine diag_export

!#######################################################################

   subroutine diag_physvar_ic (lchnk,  pbuf, cam_out, cam_in)
!
!---------------------------------------------
!
! Purpose: record physics variables on IC file
!
!---------------------------------------------
!

!
! Arguments
!
   integer       , intent(in) :: lchnk  ! chunk identifier
   type(physics_buffer_desc), pointer :: pbuf(:)
   
   type(cam_out_t), intent(inout) :: cam_out
   type(cam_in_t),  intent(inout) :: cam_in 
!
!---------------------------Local workspace-----------------------------
!
   integer  :: k                 ! indices
   integer  :: itim_old          ! indices

   real(r8), pointer, dimension(:,:) :: cwat_var
   real(r8), pointer, dimension(:,:) :: conv_var_3d
   real(r8), pointer, dimension(:  ) :: conv_var_2d
   real(r8), pointer :: tpert(:), pblh(:), qpert(:)
!
!-----------------------------------------------------------------------
!
   if( write_inithist() ) then

      !
      ! Associate pointers with physics buffer fields
      !
      itim_old = pbuf_old_tim_idx()

      if (qcwat_idx > 0) then
         call pbuf_get_field(pbuf, qcwat_idx, cwat_var, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
         call outfld('QCWAT&IC   ',cwat_var, pcols,lchnk)
      end if

      if (tcwat_idx > 0) then
         call pbuf_get_field(pbuf, tcwat_idx,  cwat_var, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
         call outfld('TCWAT&IC   ',cwat_var, pcols,lchnk)
      end if

      if (lcwat_idx > 0) then
         call pbuf_get_field(pbuf, lcwat_idx,  cwat_var, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
         call outfld('LCWAT&IC   ',cwat_var, pcols,lchnk)
      end if

      call pbuf_get_field(pbuf, cld_idx,    cwat_var, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
      call outfld('CLOUD&IC   ',cwat_var, pcols,lchnk)

      call pbuf_get_field(pbuf, concld_idx, cwat_var, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
      call outfld('CONCLD&IC   ',cwat_var, pcols,lchnk)

      call pbuf_get_field(pbuf, tke_idx, conv_var_3d)
      call outfld('TKE&IC    ',conv_var_3d, pcols,lchnk)

      call pbuf_get_field(pbuf, kvm_idx,  conv_var_3d)
      call outfld('KVM&IC    ',conv_var_3d, pcols,lchnk)

      call pbuf_get_field(pbuf, kvh_idx,  conv_var_3d)
      call outfld('KVH&IC    ',conv_var_3d, pcols,lchnk)
 
      call pbuf_get_field(pbuf, cush_idx, conv_var_2d ,(/1,itim_old/),  (/pcols,1/))
      call outfld('CUSH&IC   ',conv_var_2d, pcols,lchnk)

      if (qpert_idx > 0) then 
         call pbuf_get_field(pbuf, qpert_idx, qpert)
         call outfld('QPERT&IC   ', qpert, pcols, lchnk)
      end if

      call pbuf_get_field(pbuf, pblh_idx,  pblh)
      call outfld('PBLH&IC    ', pblh,  pcols, lchnk)

      call pbuf_get_field(pbuf, tpert_idx, tpert)
      call outfld('TPERT&IC   ', tpert, pcols, lchnk)

      ! The following is only needed for cam-csim
      call outfld('TBOT&IC    ', cam_out%tbot, pcols, lchnk)
   end if

   end subroutine diag_physvar_ic


!#######################################################################

subroutine diag_phys_tend_writeout(state, pbuf,  tend, ztodt, tmp_q, tmp_cldliq, tmp_cldice, &
                                   tmp_t, qini, cldliqini, cldiceini)

   !---------------------------------------------------------------
   !
   ! Purpose:  Dump physics tendencies for moisture and temperature
   !
   !---------------------------------------------------------------

   use check_energy,    only: check_energy_get_integrals
   use physconst,       only: cpair
   
   ! Arguments

   type(physics_state), intent(in   ) :: state 
   
   type(physics_buffer_desc), pointer :: pbuf(:)
   type(physics_tend ), intent(in   ) :: tend
   real(r8)           , intent(in   ) :: ztodt                  ! physics timestep
   real(r8)           , intent(inout) :: tmp_q     (pcols,pver) ! As input, holds pre-adjusted tracers (FV)
   real(r8)           , intent(inout) :: tmp_cldliq(pcols,pver) ! As input, holds pre-adjusted tracers (FV)
   real(r8)           , intent(inout) :: tmp_cldice(pcols,pver) ! As input, holds pre-adjusted tracers (FV)
   real(r8)           , intent(inout) :: tmp_t     (pcols,pver) ! holds last physics_updated T (FV)
   real(r8)           , intent(in   ) :: qini      (pcols,pver) ! tracer fields at beginning of physics
   real(r8)           , intent(in   ) :: cldliqini (pcols,pver) ! tracer fields at beginning of physics
   real(r8)           , intent(in   ) :: cldiceini (pcols,pver) ! tracer fields at beginning of physics

   !---------------------------Local workspace-----------------------------

   integer  :: m      ! constituent index
   integer  :: lchnk  ! chunk index
   integer  :: ncol   ! number of columns in chunk
   real(r8) :: ftem2(pcols     ) ! Temporary workspace for outfld variables
   real(r8) :: ftem3(pcols,pver) ! Temporary workspace for outfld variables
   real(r8) :: rtdt
   real(r8) :: heat_glob         ! global energy integral (FV only)
   integer  :: ixcldice, ixcldliq! constituent indices for cloud liquid and ice water.
   ! CAM pointers to get variables from the physics buffer
   real(r8), pointer, dimension(:,:) :: t_ttend  
   integer  :: itim_old

   !-----------------------------------------------------------------------

   lchnk = state%lchnk
   ncol  = state%ncol
   rtdt  = 1._r8/ztodt
   call cnst_get_ind('CLDLIQ', ixcldliq)
   call cnst_get_ind('CLDICE', ixcldice)

   ! Dump out post-physics state (FV only)

   if (dycore_is('LR') .or. dycore_is('SE')) then
      tmp_t(:ncol,:pver) = (tmp_t(:ncol,:pver) - state%t(:ncol,:pver))/ztodt
      call outfld('PTTEND_RESID', tmp_t, pcols, lchnk   )
   end if
   call outfld('TAP', state%t, pcols, lchnk   )
   call outfld('UAP', state%u, pcols, lchnk   )
   call outfld('VAP', state%v, pcols, lchnk   )

   if ( cnst_cam_outfld(       1) ) call outfld (apcnst(       1), state%q(1,1,       1), pcols, lchnk)
   if ( cnst_cam_outfld(ixcldliq) ) call outfld (apcnst(ixcldliq), state%q(1,1,ixcldliq), pcols, lchnk)
   if ( cnst_cam_outfld(ixcldice) ) call outfld (apcnst(ixcldice), state%q(1,1,ixcldice), pcols, lchnk)

   ! T-tendency due to FV Energy fixer (remove from total physics tendency diagnostic)

   if (dycore_is('LR') .or. dycore_is('SE')) then
      call check_energy_get_integrals( heat_glob_out=heat_glob )
      ftem2(:ncol)  = heat_glob/cpair
      call outfld('TFIX', ftem2, pcols, lchnk   )
      ftem3(:ncol,:pver)  = tend%dtdt(:ncol,:pver) - heat_glob/cpair
   else
      ftem3(:ncol,:pver)  = tend%dtdt(:ncol,:pver)
   end if

   ! Total physics tendency for Temperature

   call outfld('PTTEND',ftem3, pcols, lchnk )

   ! Tendency for dry mass adjustment of q (valid for FV only)

   if (dycore_is('LR')) then
      tmp_q     (:ncol,:pver) = (state%q(:ncol,:pver,       1) - tmp_q     (:ncol,:pver))*rtdt
      tmp_cldliq(:ncol,:pver) = (state%q(:ncol,:pver,ixcldliq) - tmp_cldliq(:ncol,:pver))*rtdt
      tmp_cldice(:ncol,:pver) = (state%q(:ncol,:pver,ixcldice) - tmp_cldice(:ncol,:pver))*rtdt
      if ( cnst_cam_outfld(       1) ) call outfld (dmetendnam(       1), tmp_q     , pcols, lchnk)
      if ( cnst_cam_outfld(ixcldliq) ) call outfld (dmetendnam(ixcldliq), tmp_cldliq, pcols, lchnk)
      if ( cnst_cam_outfld(ixcldice) ) call outfld (dmetendnam(ixcldice), tmp_cldice, pcols, lchnk)
   end if

   ! Total physics tendency for moisture and other tracers

   if ( cnst_cam_outfld(       1) ) then
      ftem3(:ncol,:pver) = (state%q(:ncol,:pver,       1) - qini     (:ncol,:pver) )*rtdt
      call outfld (ptendnam(       1), ftem3, pcols, lchnk)
   end if
   if ( cnst_cam_outfld(ixcldliq) ) then
      ftem3(:ncol,:pver) = (state%q(:ncol,:pver,ixcldliq) - cldliqini(:ncol,:pver) )*rtdt
      call outfld (ptendnam(ixcldliq), ftem3, pcols, lchnk)
   end if
   if ( cnst_cam_outfld(ixcldice) ) then
      ftem3(:ncol,:pver) = (state%q(:ncol,:pver,ixcldice) - cldiceini(:ncol,:pver) )*rtdt
      call outfld (ptendnam(ixcldice), ftem3, pcols, lchnk)
   end if

   ! Total (physics+dynamics, everything!) tendency for Temperature

   !! get temperature stored in physics buffer
   itim_old = pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, t_ttend_idx, t_ttend, start=(/1,1,itim_old/), kount=(/pcols,pver,1/))

   !! calculate and outfld the total temperature tendency
   ftem3(:ncol,:) = (state%t(:ncol,:) - t_ttend(:ncol,:))/ztodt
   call outfld('TTEND_TOT', ftem3, pcols, lchnk)

   !! update physics buffer with this time-step's temperature
   t_ttend(:ncol,:) = state%t(:ncol,:)

end subroutine diag_phys_tend_writeout

!#######################################################################

   subroutine diag_state_b4_phys_write (state)
!
!---------------------------------------------------------------
!
! Purpose:  Dump state just prior to executing physics
!
!---------------------------------------------------------------
!
! Arguments
!
   type(physics_state), intent(in) :: state 
!
!---------------------------Local workspace-----------------------------
!
   integer :: ixcldice, ixcldliq ! constituent indices for cloud liquid and ice water.
   integer :: lchnk              ! chunk index
!
!-----------------------------------------------------------------------
!
   lchnk = state%lchnk

   call cnst_get_ind('CLDLIQ', ixcldliq)
   call cnst_get_ind('CLDICE', ixcldice)
   call outfld('TBP', state%t, pcols, lchnk   )
   if ( cnst_cam_outfld(       1) ) call outfld (bpcnst(       1), state%q(1,1,       1), pcols, lchnk)
   if ( cnst_cam_outfld(ixcldliq) ) call outfld (bpcnst(ixcldliq), state%q(1,1,ixcldliq), pcols, lchnk)
   if ( cnst_cam_outfld(ixcldice) ) call outfld (bpcnst(ixcldice), state%q(1,1,ixcldice), pcols, lchnk)

   end subroutine diag_state_b4_phys_write

end module cam_diagnostics
