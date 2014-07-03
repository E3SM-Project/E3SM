module aerosol_intr

!---------------------------------------------------------------------------------
! Module to interface the aerosol parameterizations with CAM
! Original version: Phil Rasch, Jan 2003
!---------------------------------------------------------------------------------

use shr_kind_mod,      only: r8 => shr_kind_r8, cl => shr_kind_cl
use spmd_utils,        only: masterproc
use ppgrid,            only: pcols, pver, pverp
use camsrfexch,        only: cam_in_t, cam_out_t, hub2atm_setopts
use physconst,         only: mwdry, mwh2o, gravit, rair
use phys_control,      only: cam_physpkg_is, phys_getopts
use phys_grid,         only: get_lat_all_p, get_rlat_all_p
use physics_types,     only: physics_state, physics_ptend, physics_ptend_init
use physics_buffer,    only: physics_buffer_desc, pbuf_get_field, pbuf_old_tim_idx, pbuf_get_index

use constituents,      only: pcnst, cnst_name
use rad_constituents,  only: rad_cnst_get_info
use dust_intr,         only: dust_initialize, dust_wet_intr, dust_drydep_intr
use progseasalts_intr, only: progseasalts_initialize, progseasalts_wet_intr, progseasalts_drydep_intr
use drydep_mod,        only: inidrydep
use wetdep,            only: clddiag

use time_manager,      only: is_first_step, get_curr_date, get_perp_date, is_perpetual, &
                             get_curr_calday, get_nstep

use modal_aero_calcsize,    only: modal_aero_calcsize_init, modal_aero_calcsize_sub, modal_aero_calcsize_diag
use modal_aero_wateruptake, only: modal_aero_wateruptake_init, modal_aero_wateruptake_dr
use mz_aerosols_intr,       only: mz_aero_wet_intr

use cam_logfile,       only: iulog
use perf_mod,          only: t_startf, t_stopf
use abortutils,        only: endrun

#if ( defined MODAL_AERO )
use mz_aerosols_intr,  only: mz_aero_dry_intr
#endif

implicit none
private
save

! Public interfaces

public ::&
   aerosol_readnl,                &  ! read aerosol_nl namelist group
   aerosol_register_cnst,         &  ! register consituents
   aerosol_init,                  &  ! call aerosol init routines
   aerosol_implements_cnst,       &  ! returns true if consituent is implemented by this package
   aerosol_init_cnst,             &  ! initialize mixing ratios if not read from initial file
   aerosol_drydep_intr,           &  ! interface to dry deposition
   aerosol_wet_intr,              &  ! interface to wet deposition
   aerosol_emis_intr                 ! interface to surface emissions

! prognostic bulk sea salt
logical :: progsslt = .false.

! prognostic bulk dust
logical :: dust = .false.

! modal aerosols
logical :: prog_modal_aero ! if true then prognostic modal aerosols are present 

! modal aerosols
logical :: clim_modal_aero ! if true then modal aerosols are used in the climate calculation.

! is_any_aerosol is set to .true. to indicate that prognostic dust, sea salt, or
! modal aerosols are present.  It is only used to initialize dry deposition module
logical :: is_any_aerosol = .false.

! Physics buffer indices
integer :: dgnum_idx           = 0
integer :: dgnumwet_idx        = 0
integer :: wetdens_ap_idx      = 0
integer :: qaerwat_idx         = 0
integer :: fracis_idx          = 0

integer :: cld_idx             = 0
integer :: qme_idx             = 0 
integer :: prain_idx           = 0 
integer :: nevapr_idx          = 0 
integer :: icwmrdp_idx         = 0 
integer :: rprddp_idx          = 0 
integer :: icwmrsh_idx         = 0 
integer :: rprdsh_idx          = 0 
integer :: sh_frac_idx         = 0 
integer :: dp_frac_idx         = 0 
integer :: nevapr_shcu_idx     = 0 
integer :: nevapr_dpcu_idx     = 0 
integer :: pblh_idx            = 0

integer :: rate1_cw2pr_st_idx  = 0  


! Namelist variables
real(r8)      :: dust_emis_fact = -1.e36_r8   ! tuning parameter for dust emissions
character(cl) :: soil_erod = 'soil_erod'   ! full pathname for soil erodibility dataset

!===============================================================================
contains
!===============================================================================

subroutine aerosol_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'aerosol_readnl'

   namelist /aerosol_nl/ dust_emis_fact, soil_erod
   !-----------------------------------------------------------------------------

   ! Read namelist
   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'aerosol_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, aerosol_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(dust_emis_fact,               1, mpir8,   0, mpicom)
   call mpibcast(soil_erod,       len(soil_erod), mpichar, 0, mpicom)
#endif

end subroutine aerosol_readnl

!===============================================================================

subroutine aerosol_register_cnst

   ! Some misc register actions.  Mostly moved to chemistry code.

   use dust_intr,        only: ndst=>dust_number, dust_names
   use progseasalts_intr,only: nsst, progseasalts_names
   use carma_flags_mod,  only: carma_do_drydep
   use physics_buffer,   only: pbuf_add_field, dtype_r8
   use mo_chem_utls,     only: get_spc_ndx

   !---------------------------Local workspace-----------------------------
   integer :: m
   integer :: nmodes
   integer :: dst_idx(ndst), slt_idx(nsst)
   !-----------------------------------------------------------------------

   do m=1,ndst
      dst_idx(m) = get_spc_ndx( dust_names(m) )
   enddo
   do m=1,nsst
      slt_idx(m) = get_spc_ndx( progseasalts_names(m) )
   enddo

   dust      = all( dst_idx > 0 )
   progsslt  = all( slt_idx > 0 )

   ! prog_modal_aero determines whether prognostic modal aerosols are present in the run.
   call phys_getopts(prog_modal_aero_out=prog_modal_aero)

   ! ***N.B.*** The wet radius and water uptake calculations are currently only being done
   !            for the aerosol modes in the climate calc.  This could be extended to compute
   !            these quantities for the modes used in diagnostic radiative calcs as well.
   !
   ! clim_modal_aero determines whether modal aerosols are used in the climate calculation.
   ! The modal aerosols can be either prognostic or prescribed.
   call rad_cnst_get_info(0, nmodes=nmodes)
   clim_modal_aero = (nmodes > 0)

   is_any_aerosol = (dust .or. progsslt .or. prog_modal_aero)

   if (is_any_aerosol .or. carma_do_drydep) then

      ! tell camsrfexch to allocate fv & ram1 -- needed by prodsslts and dust
      call hub2atm_setopts(aero_dust_in=.true.)

   end if

   ! The following fields are diagnosed for either prognostic or prescribed
   ! modal aerosols.  Currently only modes that affect the climate are
   ! accounted for.
   if (clim_modal_aero) then
      call pbuf_add_field('DGNUM',      'global',  dtype_r8, (/pcols, pver, nmodes/), dgnum_idx)    
      call pbuf_add_field('DGNUMWET',   'global',  dtype_r8, (/pcols, pver, nmodes/), dgnumwet_idx)
      call pbuf_add_field('WETDENS_AP', 'physpkg', dtype_r8, (/pcols, pver, nmodes/), wetdens_ap_idx)
      call pbuf_add_field('QAERWAT',    'physpkg', dtype_r8, (/pcols, pver, nmodes/), qaerwat_idx)  ! 1st order rate for direct conversion of
                                                                                                    ! strat. cloud water to precip (1/s)
   end if

   call pbuf_add_field('FRACIS','physpkg',dtype_r8,(/pcols,pver,pcnst/),fracis_idx)

end subroutine aerosol_register_cnst

!=======================================================================

subroutine aerosol_init(pbuf2d)

   ! initialize aerosol parameterizations -- mostly moved into chemistry code

   use physics_buffer,only   : pbuf_set_field

   type(physics_buffer_desc), pointer :: pbuf2d(:,:)

   !----------------------------------------------------------------------- 
    
   if (dust .or. prog_modal_aero) then
      call dust_initialize(dust_emis_fact, soil_erod)
   end if
   
   if (progsslt .or. prog_modal_aero) then
      call progseasalts_initialize
   end if

   if (clim_modal_aero) then

      if (.not. prog_modal_aero) then
         ! If climate calculations are affected by prescribed modal aerosols, the
         ! the initialization routine for the dry mode radius calculation is called
         ! here.  For prognostic MAM the initialization is called from
         ! modal_aero_initialize
         call modal_aero_calcsize_init()
      endif

      call modal_aero_wateruptake_init()

   end if

   ! Dry deposition needs to be initialized if any of the aerosols
   ! are running.
   if (is_any_aerosol) then
      call inidrydep(rair, gravit)
   endif

   cld_idx             = pbuf_get_index('CLD')    
   qme_idx             = pbuf_get_index('QME')    
   prain_idx           = pbuf_get_index('PRAIN')  
   nevapr_idx          = pbuf_get_index('NEVAPR') 
   icwmrdp_idx         = pbuf_get_index('ICWMRDP') 
   rprddp_idx          = pbuf_get_index('RPRDDP')  
   icwmrsh_idx         = pbuf_get_index('ICWMRSH') 
   rprdsh_idx          = pbuf_get_index('RPRDSH')  
   sh_frac_idx         = pbuf_get_index('SH_FRAC' )
   dp_frac_idx         = pbuf_get_index('DP_FRAC') 
   nevapr_shcu_idx     = pbuf_get_index('NEVAPR_SHCU') 
   nevapr_dpcu_idx     = pbuf_get_index('NEVAPR_DPCU') 
   pblh_idx            = pbuf_get_index('pblh')

   if (prog_modal_aero) then
      rate1_cw2pr_st_idx  = pbuf_get_index('RATE1_CW2PR_ST') 
   end if

   if (clim_modal_aero .and. is_first_step()) then
      ! initialize fields in physics buffer
      call pbuf_set_field(pbuf2d, dgnum_idx,      0.0_r8)
      call pbuf_set_field(pbuf2d, dgnumwet_idx,   0.0_r8)
      call pbuf_set_field(pbuf2d, wetdens_ap_idx, 0.0_r8)
      call pbuf_set_field(pbuf2d, qaerwat_idx, 0.0_r8)
   end if

end subroutine aerosol_init

!===============================================================================

function aerosol_implements_cnst(name)

   ! return true if specified constituent is implemented by this package

   use dust_intr,         only: dust_implements_cnst
   use progseasalts_intr, only: progseasalts_implements_cnst

   !-----------------------------Arguments---------------------------------
   character(len=*), intent(in) :: name   ! constituent name
   logical :: aerosol_implements_cnst        ! return value

   !---------------------------Local workspace-----------------------------
   integer :: m
   !-----------------------------------------------------------------------

   aerosol_implements_cnst = .false.

   if ( dust ) then
      aerosol_implements_cnst = &
         (aerosol_implements_cnst.OR.dust_implements_cnst (name))
   end if

   if ( progsslt ) then
      aerosol_implements_cnst = &
         (aerosol_implements_cnst.OR.progseasalts_implements_cnst (name))
   end if

end function aerosol_implements_cnst

!=======================================================================

subroutine aerosol_init_cnst(name, q, gcid)

   ! Set initial mass mixing ratios.  

   use dust_intr,        only: dust_implements_cnst, dust_init_cnst
   use progseasalts_intr,only: progseasalts_implements_cnst, progseasalts_init_cnst

   !-----------------------------Arguments---------------------------------
   character(len=*), intent(in)  :: name      ! constituent name
   integer,          intent(in)  :: gcid(:)   ! global column id
   real(r8),         intent(out) :: q(:,:)    !  mass mixing ratio
   !-----------------------------------------------------------------------

   if ( dust ) then
      if (dust_implements_cnst(name)) then
         call dust_init_cnst(name, q, gcid)
      end if
   end if

   if ( progsslt ) then
      if (progseasalts_implements_cnst(name)) then
         call progseasalts_init_cnst(name, q, gcid)
      end if
   end if

end subroutine aerosol_init_cnst

!===============================================================================
  subroutine aerosol_wet_intr (state, ptend, dt, pbuf,  cam_out, dlf)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Interface to wet processing of all aerosols
! 
! Method: 
!  use a modified version of the scavenging parameterization described in
!     Barth et al, 2000, JGR (sulfur cycle paper)
!     Rasch et al, 2001, JGR (INDOEX paper)
! 
! Author: Phil Rasch
! 
!-----------------------------------------------------------------------

   ! Arguments:

   type(physics_state), intent(in )   :: state       ! Physics state variables
   type(physics_ptend), intent(out)   :: ptend       ! indivdual parameterization tendencies
   real(r8),            intent(in)    :: dt          ! time step
   type(physics_buffer_desc), pointer :: pbuf(:)
   type(cam_out_t),     intent(inout) :: cam_out     ! export state

   real(r8), intent(in) :: dlf(pcols,pver)           ! shallow+deep convective detrainment [kg/kg/s]

   ! local vars

   integer  :: ncol, lchnk, nstep
   integer  :: nmodes
   real(r8) :: calday        ! current calendar day

   real(r8) :: rainmr(pcols,pver)   ! mixing ratio of rain within cloud volume
   real(r8) :: cldv(pcols,pver)     ! cloudy volume undergoing wet chem and scavenging
   real(r8) :: cldvcu(pcols,pver)   ! Convective precipitation area at the top interface of current layer
   real(r8) :: cldvst(pcols,pver)   ! Stratiform precipitation area at the top interface of current layer 
   integer  :: lat(pcols)                  ! latitude indices
   real(r8) :: clat(pcols)                    ! latitudes
   real(r8) :: conicw(pcols,pver)          ! convective in-cloud water
   real(r8) :: cmfdqr(pcols,pver)          ! convective production of rain
   real(r8) :: cldc(pcols,pver)            ! convective cloud fraction, currently empty
   real(r8) :: clds(pcols,pver)            ! Stratiform cloud fraction
   real(r8) :: evapc(pcols,pver)           ! Evaporation rate of convective precipitation

   ! physics buffer 
   integer itim
   real(r8), pointer :: cldn(:,:)       ! cloud fraction
   real(r8), pointer :: cme(:,:)
   real(r8), pointer :: prain(:,:)
   real(r8), pointer :: evapr(:,:)
   real(r8), pointer :: icwmrdp(:,:)    ! in cloud water mixing ratio, deep convection
   real(r8), pointer :: rprddp(:,:)     ! rain production, deep convection
   real(r8), pointer :: icwmrsh(:,:)    ! in cloud water mixing ratio, deep convection
   real(r8), pointer :: rprdsh(:,:)     ! rain production, deep convection
   real(r8), pointer :: fracis(:,:,:)   ! fraction of transported species that are insoluble
   real(r8), pointer :: sh_frac(:,:)    ! Shallow convective cloud fraction
   real(r8), pointer :: dp_frac(:,:)    ! Deep convective cloud fraction
   real(r8), pointer :: evapcsh(:,:)    ! Evaporation rate of shallow convective precipitation >=0.
   real(r8), pointer :: evapcdp(:,:)    ! Evaporation rate of deep    convective precipitation >=0.
   real(r8), pointer :: dgnumwet_pbuf(:,:,:)
   real(r8), pointer :: qaerwat(:,:,:)  ! aerosol water
   real(r8), pointer :: rate1ord_cw2pr_st(:,:)

   logical           :: lq(pcnst)

   !-----------------------------------------------------------------------

   nstep = get_nstep()

   calday = get_curr_calday()

   ncol = state%ncol
   lchnk = state%lchnk

   ! Allow the program to determine which values of lq are to be set later.  Temporarily turn on
   ! lq so the %q fields are allocated and then set them all back to false.
   lq(:)    = .TRUE.    ! Set to allocate q in ptend
   call physics_ptend_init(ptend, state%psetcols, 'wetdep', lq=lq)
   ptend%lq(:)    = .FALSE.   ! reset to default

   call get_lat_all_p(lchnk, ncol, lat)
   call get_rlat_all_p(lchnk, ncol, clat)

! Associate pointers with physics buffer fields
   itim = pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, cld_idx,         cldn,     start=(/1,1,itim/), kount=(/pcols,pver,1/) )

   call pbuf_get_field(pbuf, qme_idx,         cme     )
   call pbuf_get_field(pbuf, prain_idx,       prain   )
   call pbuf_get_field(pbuf, nevapr_idx,      evapr   )
   call pbuf_get_field(pbuf, fracis_idx,      fracis,   start=(/1,1,1/), kount=(/pcols, pver, pcnst/) )
   call pbuf_get_field(pbuf, icwmrdp_idx,     icwmrdp )
   call pbuf_get_field(pbuf, rprddp_idx,      rprddp  )
   call pbuf_get_field(pbuf, icwmrsh_idx,     icwmrsh )
   call pbuf_get_field(pbuf, rprdsh_idx,      rprdsh  )
   call pbuf_get_field(pbuf, sh_frac_idx,     sh_frac )
   call pbuf_get_field(pbuf, dp_frac_idx,     dp_frac )
   call pbuf_get_field(pbuf, nevapr_shcu_idx, evapcsh )
   call pbuf_get_field(pbuf, nevapr_dpcu_idx, evapcdp )

   cldc(:ncol,:)  = dp_frac(:ncol,:) + sh_frac(:ncol,:)
   evapc(:ncol,:) = evapcsh(:ncol,:) + evapcdp(:ncol,:)
   clds(:ncol,:)  = cldn(:ncol,:) - cldc(:ncol,:)       ! Stratiform cloud fraction

   ! sum deep and shallow convection contributions
   if (cam_physpkg_is('cam5')) then
      ! Dec.29.2009. Sungsu
      conicw(:ncol,:) = (icwmrdp(:ncol,:)*dp_frac(:ncol,:) + icwmrsh(:ncol,:)*sh_frac(:ncol,:))/ &
                        max(0.01_r8, sh_frac(:ncol,:) + dp_frac(:ncol,:))
   else
      conicw(:ncol,:) = icwmrdp(:ncol,:) + icwmrsh(:ncol,:)
   end if

   cmfdqr(:ncol,:) = rprddp(:ncol,:)  + rprdsh(:ncol,:)


   !   fields needed for wet scavenging
   call clddiag( &
      state%t, state%pmid, state%pdel, cmfdqr, evapc, &
      cldn, cldc, clds, cme, evapr,                   &
      prain, cldv, cldvcu, cldvst, rainmr,            &
      ncol)

   if (clim_modal_aero .or. prog_modal_aero) then

      ! Do calculations of mode radius and water uptake if:
      ! 1) modal aerosols are affecting the climate, or
      ! 2) prognostic modal aerosols are enabled

      call t_startf('calcsize')

      if (prog_modal_aero) then
         ! for prognostic modal aerosols the transfer of mass between aitken and accumulation
         ! modes is done in conjunction with the dry radius calculation
         call modal_aero_calcsize_sub(state, ptend, dt, pbuf)
      else
         ! for prescribed aerosols just do a diagnostic dry radius calc
         call modal_aero_calcsize_diag(state, pbuf)
      end if

      call t_stopf('calcsize')

      call t_startf('wateruptake')

      call modal_aero_wateruptake_dr(state, pbuf)

      call t_stopf('wateruptake')

   end if

#if ( defined MODAL_AERO )

   call rad_cnst_get_info(0, nmodes=nmodes)
   call pbuf_get_field(pbuf, dgnumwet_idx,       dgnumwet_pbuf,    start=(/1,1,1/), kount=(/pcols,pver,nmodes/) )
   call pbuf_get_field(pbuf, qaerwat_idx,        qaerwat,          start=(/1,1,1/), kount=(/pcols,pver,nmodes/) )
   call pbuf_get_field(pbuf, rate1_cw2pr_st_idx, rate1ord_cw2pr_st)
   call mz_aero_wet_intr (state, ptend, nstep, dt, cme, prain, &
            evapr, cldv, cldvcu, cldvst, cldc, cldn, fracis, calday, cmfdqr, evapc, conicw, rainmr, &
            rate1ord_cw2pr_st, &   ! rce 2010/05/01
            dgnumwet_pbuf, qaerwat, cam_out, dlf, pbuf)

#else
   if ( dust ) then
      !   wet scavenging for dust
      call dust_wet_intr (state, ptend, nstep, dt, lat, clat, cme, prain, &
         evapr, cldv, cldc, cldn, fracis, calday, cmfdqr, conicw, rainmr, cam_out)
   endif

   if ( progsslt ) then
      !   wet scavenging for prognostic sea salts
      call progseasalts_wet_intr (state, ptend, nstep, dt, lat, clat, cme, prain, &
            evapr, cldv, cldc, cldn, fracis, calday, cmfdqr, conicw, rainmr)
   endif

!   wet scavenging for mozart aerosols
! can not be done under trop_mozart chem_timestep_tend --
! this need to be done before deep convection so that fracis is set correctly
! for the mozart aerosols -- fracis is used in deep convection routine
   call mz_aero_wet_intr (state, ptend, nstep, dt, cme, prain, &
            evapr, cldv, cldvcu, cldvst, cldc, cldn, fracis, calday, cmfdqr, evapc, conicw, rainmr, cam_out, dlf, pbuf)
#endif

end subroutine aerosol_wet_intr

!===============================================================================

subroutine aerosol_drydep_intr (state, ptend, cam_in, cam_out, dt, &
                                fsds, obklen, ustar, prect, pbuf )

   !----------------------------------------------------------------------- 
   ! 
   ! Purpose: 
   ! Interface to dry deposition parameterizions and sedimentation of all aerosols
   ! 
   ! Method: 
   ! Use prescribed dry deposition velocities for sulfur and carbon
   ! Use calculated dry dep velocities from CLM for dust and prognostic seasalt
   ! 
   ! Author: Phil Rasch
   ! 
   !-----------------------------------------------------------------------

   ! Arguments:
   type(physics_state),    intent(in )   :: state   ! Physics state variables
   type(physics_ptend),    intent(out)   :: ptend   ! indivdual parameterization tendencies
   type(cam_in_t), target, intent(in)    :: cam_in  ! import state
   type(cam_out_t),        intent(inout) :: cam_out ! export state
   real(r8),               intent(in)    :: dt      ! time step
   real(r8), intent(in) :: fsds(pcols)                 ! longwave down at sfc
   real(r8), intent(in) :: obklen(pcols)                 ! longwave down at sfc
   real(r8), intent(in) :: ustar(pcols)                  ! sfc fric vel
   real(r8), intent(in) :: prect(pcols)                 ! prect

    type(physics_buffer_desc), pointer, dimension(:) :: pbuf

   ! Local variables:

   integer  :: lchnk
   integer  :: ncol
   integer  :: nmodes
   integer  :: lat(pcols)     ! latitude index for S->N storage
   real(r8) :: clat(pcols)    ! latitude 
   integer :: yr, mon, day, ncsec

   real(r8), pointer, dimension(:)  :: pblh        ! pbl height
   real(r8), pointer, dimension(:)  :: wvflx	   ! (pcols) water vapor flux
   real(r8), pointer, dimension(:)  :: ts	   ! (pcols) sfc temp
   real(r8), pointer, dimension(:)  :: snowh	   ! (pcols) snow depth
   real(r8), pointer, dimension(:)  :: hflx	   ! (pcols) sensible heat flux
   real(r8), pointer, dimension(:)  :: landfrac    ! (pcols) land fraction
   real(r8), pointer, dimension(:)  :: icefrac	   ! (pcols) ice fraction
   real(r8), pointer, dimension(:)  :: ocnfrac	   ! (pcols) ocn fraction

   real(r8), pointer :: dgnumwet_pbuf(:,:,:)
   real(r8), pointer :: wetdens_pbuf(:,:,:)
   real(r8), pointer :: qaerwat(:,:,:)

   logical           :: lq(pcnst)
 
   !-----------------------------------------------------------------------

   lchnk = state%lchnk
   ncol = state%ncol
    
   wvflx    => cam_in%cflx(:,1)
   ts       => cam_in%tref(:)
   snowh    => cam_in%snowhland(:)
   hflx     => cam_in%shf(:)
   landfrac => cam_in%landfrac(:)
   icefrac  => cam_in%icefrac(:)
   ocnfrac  => cam_in%ocnfrac(:)

   if ( is_perpetual() ) then
      call get_perp_date(yr, mon, day, ncsec)
   else
      call get_curr_date(yr, mon, day, ncsec)
   end if

   call get_lat_all_p(lchnk, ncol, lat)
   call get_rlat_all_p(lchnk, ncol, clat)

   ! note that tendencies are not only in sfc layer (because of sedimentation)
   ! and that ptend is updated within each subroutine for different species

   ! Allow the program to determine which values of lq are to be set later.  Temporarily turn on
   ! lq so the %q fields are allocated and then set them all back to false.
   lq(:)= .TRUE. ! Set to allocate q in ptend
   call physics_ptend_init(ptend, state%psetcols, 'drydep', lq=lq)
   ptend%lq(:)= .FALSE. ! reset to default

   call pbuf_get_field(pbuf, pblh_idx, pblh)

#if ( defined MODAL_AERO )
    call rad_cnst_get_info(0, nmodes=nmodes)
    call pbuf_get_field(pbuf, dgnumwet_idx,   dgnumwet_pbuf, start=(/1,1,1/), kount=(/pcols,pver,nmodes/) ) 
    call pbuf_get_field(pbuf, wetdens_ap_idx, wetdens_pbuf,  start=(/1,1,1/), kount=(/pcols,pver,nmodes/) ) 
    call pbuf_get_field(pbuf, qaerwat_idx,    qaerwat,       start=(/1,1,1/), kount=(/pcols,pver,nmodes/) ) 

   call mz_aero_dry_intr(state, ptend, wvflx, dt, lat, clat, &
            fsds, obklen, ts, ustar, prect, snowh, pblh, hflx, mon, &
            landfrac, icefrac, ocnfrac, cam_in%fv, cam_in%ram1, &
            dgnumwet_pbuf, wetdens_pbuf, qaerwat, cam_out, pbuf)
#endif

   if (dust .and. .not. prog_modal_aero) then
      call dust_drydep_intr (state, ptend, wvflx, dt, lat, clat, &
            fsds, obklen, ts, ustar, prect, snowh, pblh, hflx, mon, &
            landfrac, icefrac, ocnfrac, cam_in%fv,  cam_in%ram1, cam_out)
   endif

   if (progsslt .and. .not. prog_modal_aero) then
      call progseasalts_drydep_intr (state, ptend, wvflx, dt, lat, clat, &
            fsds, obklen, ts, ustar, prect, snowh, pblh, hflx, mon, landfrac, &
            icefrac, ocnfrac, cam_in%fv,  cam_in%ram1)
   endif



end subroutine aerosol_drydep_intr

!===============================================================================

subroutine aerosol_emis_intr (state, cam_in)

   ! return surface fluxes of aerosol species and tendencies in surface layer
   ! due to surface emissions

   use dust_intr,          only: dust_emis_intr
   use progseasalts_intr,  only: progseasalts_emis_intr

   ! Arguments:
   type(physics_state), intent(in )   :: state   ! Physics state variables
   type(cam_in_t),      intent(inout) :: cam_in  ! import state

   !-----------------------------------------------------------------------

   if (dust .or. prog_modal_aero) then
      call dust_emis_intr(state, cam_in)
   end if

   if (progsslt .or. prog_modal_aero) then
      call progseasalts_emis_intr(state, cam_in)
   end if

end subroutine aerosol_emis_intr

!===============================================================================

end module aerosol_intr
