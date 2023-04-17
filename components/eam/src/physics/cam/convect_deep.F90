

module convect_deep
!---------------------------------------------------------------------------------
! Purpose:
!
! CAM interface to several deep convection interfaces. Currently includes:
!    Zhang-McFarlane (default)
!    Kerry Emanuel 
!
!
! Author: D.B. Coleman, Sep 2004
!
!---------------------------------------------------------------------------------
   use shr_kind_mod, only: r8=>shr_kind_r8
   use ppgrid,       only: pver, pcols, pverp
   use cam_logfile,  only: iulog
   use perf_mod,     only: t_startf, t_stopf

   implicit none

   save
   private                         ! Make default type private to the module

! Public methods

   public ::&
      convect_deep_register,           &! register fields in physics buffer
      convect_deep_init,               &! initialize donner_deep module
      convect_deep_tend,               &! return tendencies
      convect_deep_tend_2,             &! return tendencies
      deep_scheme_does_scav_trans             ! = .t. if scheme does scavenging and conv. transport
   
! Private module data
   character(len=16) :: deep_scheme    ! default set in phys_control.F90, use namelist to change
! Physics buffer indices 
   integer     ::  icwmrdp_idx      = 0 
   integer     ::  rprddp_idx       = 0 
   integer     ::  nevapr_dpcu_idx  = 0 
   integer     ::  cldtop_idx       = 0 
   integer     ::  cldbot_idx       = 0 
   integer     ::  cld_idx          = 0 
   integer     ::  fracis_idx       = 0 

   integer     ::  pblh_idx        = 0 
   integer     ::  tpert_idx       = 0 
   integer     ::  prec_dp_idx     = 0
   integer     ::  snow_dp_idx     = 0
   
   integer     ::  dp_cldliq_idx   = 0
   integer     ::  dp_cldice_idx   = 0
   integer     ::  dp_flxprc_idx   = 0
   integer     ::  dp_flxsnw_idx   = 0
   
   integer     ::  dp_frac_idx     = 0

   integer     ::  ttend_dp_idx        = 0

!=========================================================================================
  contains 

!=========================================================================================
function deep_scheme_does_scav_trans()
!
! Function called by tphysbc to determine if it needs to do scavenging and convective transport
! or if those have been done by the deep convection scheme. Each scheme could have its own
! identical query function for a less-knowledgable interface but for now, we know that KE 
! does scavenging & transport, and ZM doesn't
!

  logical deep_scheme_does_scav_trans

  deep_scheme_does_scav_trans = .false.

  if ( deep_scheme .eq. 'KE' ) deep_scheme_does_scav_trans = .true.

  return

end function deep_scheme_does_scav_trans

!=========================================================================================
subroutine convect_deep_register

!----------------------------------------
! Purpose: register fields with the physics buffer
!----------------------------------------

  
  use physics_buffer, only : pbuf_add_field, dtype_r8
  use zm_conv_intr, only: zm_conv_register
  use phys_control, only: phys_getopts, use_gw_convect

  implicit none

  integer idx

  ! get deep_scheme setting from phys_control
  call phys_getopts(deep_scheme_out = deep_scheme)

  select case ( deep_scheme )
  case('ZM') !    Zhang-McFarlane (default)
     call zm_conv_register
  end select 

! Add PBUF variables that are related to deep convection that 
!  are expected by the E3SM code, no matter which deep convection scheme
!  is used.

! Flux of precipitation from deep convection (kg/m2/s)
   call pbuf_add_field('DP_FLXPRC','global',dtype_r8,(/pcols,pverp/),dp_flxprc_idx) 

! Flux of snow from deep convection (kg/m2/s) 
   call pbuf_add_field('DP_FLXSNW','global',dtype_r8,(/pcols,pverp/),dp_flxsnw_idx) 

! deep gbm cloud liquid water (kg/kg)
   call pbuf_add_field('DP_CLDLIQ','global',dtype_r8,(/pcols,pver/), dp_cldliq_idx)  

! deep gbm cloud liquid water (kg/kg)    
   call pbuf_add_field('DP_CLDICE','global',dtype_r8,(/pcols,pver/), dp_cldice_idx) 
  
  call pbuf_add_field('ICWMRDP',    'physpkg',dtype_r8,(/pcols,pver/),icwmrdp_idx)
  call pbuf_add_field('RPRDDP',     'physpkg',dtype_r8,(/pcols,pver/),rprddp_idx)
  call pbuf_add_field('NEVAPR_DPCU','physpkg',dtype_r8,(/pcols,pver/),nevapr_dpcu_idx)
  call pbuf_add_field('PREC_DP',    'physpkg',dtype_r8,(/pcols/),     prec_dp_idx)
  call pbuf_add_field('SNOW_DP',   'physpkg',dtype_r8,(/pcols/),      snow_dp_idx)

  ! If WACCM gravity waves are on, output this field.
  if (use_gw_convect) then
     call pbuf_add_field('TTEND_DP','physpkg',dtype_r8,(/pcols,pver/),ttend_dp_idx)
  end if

end subroutine convect_deep_register

!=========================================================================================



subroutine convect_deep_init(pref_edge,pbuf2d)

!----------------------------------------
! Purpose:  declare output fields, initialize variables needed by convection
!----------------------------------------

  use cam_history,   only:  addfld                          
  use pmgrid,        only: plevp
  use spmd_utils,    only: masterproc
  use zm_conv_intr,  only: zm_conv_init
  use cam_abortutils,    only: endrun
  
  use physics_buffer, only: physics_buffer_desc, pbuf_get_index, pbuf_set_field

  implicit none

  real(r8),intent(in) :: pref_edge(plevp)        ! reference pressures at interfaces
  type(physics_buffer_desc), pointer    :: pbuf2d(:,:)

  dp_frac_idx = pbuf_get_index('DP_FRAC')
  icwmrdp_idx = pbuf_get_index('ICWMRDP')

  select case ( deep_scheme )
  case('off') !     ==> no deep convection
     if (masterproc) write(iulog,*)'convect_deep: no deep convection selected'
     call pbuf_set_field(pbuf2d, dp_cldliq_idx, 0._r8)
     call pbuf_set_field(pbuf2d, dp_cldice_idx, 0._r8)
     call pbuf_set_field(pbuf2d, dp_flxprc_idx, 0._r8)
     call pbuf_set_field(pbuf2d, dp_flxsnw_idx, 0._r8)
     call pbuf_set_field(pbuf2d, dp_frac_idx, 0._r8)
     call pbuf_set_field(pbuf2d, icwmrdp_idx, 0._r8)
     call pbuf_set_field(pbuf2d, rprddp_idx, 0._r8)
     call pbuf_set_field(pbuf2d, nevapr_dpcu_idx, 0._r8)
  case('CLUBB_SGS')
     if (masterproc) write(iulog,*)'convect_deep: CLUBB_SGS selected'
  case('ZM') !    1 ==> Zhang-McFarlane (default)
     if (masterproc) write(iulog,*)'convect_deep initializing Zhang-McFarlane convection'
     call zm_conv_init(pref_edge)
  case default
     if (masterproc) write(iulog,*)'WARNING: convect_deep: no deep convection scheme. May fail.'
  end select

  cldtop_idx = pbuf_get_index('CLDTOP')
  cldbot_idx = pbuf_get_index('CLDBOT')
  cld_idx    = pbuf_get_index('CLD')
  fracis_idx = pbuf_get_index('FRACIS')

  pblh_idx   = pbuf_get_index('pblh')
  tpert_idx  = pbuf_get_index('tpert')

  call addfld ('ICWMRDP', (/ 'lev' /), 'A', 'kg/kg', 'Deep Convection in-cloud water mixing ratio '            )

end subroutine convect_deep_init
!=========================================================================================
!subroutine convect_deep_tend(state, ptend, tdt, pbuf)

subroutine convect_deep_tend( &
     mcon    ,cme     ,          &
     dlf     ,pflx    ,zdu      , &
     rliq    , &
     ztodt   , &
     state   ,ptend   ,landfrac ,pbuf, mu, eu, &
     du, md, ed, dp, dsubcld, jt, maxg, ideep,lengath ) 



   use physics_types, only: physics_state, physics_ptend, physics_tend, physics_ptend_init
   
   use cam_history,    only: outfld
   use constituents,   only: pcnst
   use zm_conv_intr,   only: zm_conv_tend
   use cam_history,    only: outfld
   use physconst,      only: cpair
   use physics_buffer, only: physics_buffer_desc, pbuf_get_field

! Arguments
   type(physics_state), intent(in ) :: state   ! Physics state variables
   type(physics_ptend), intent(out) :: ptend   ! individual parameterization tendencies
   

   type(physics_buffer_desc), pointer :: pbuf(:)
   real(r8), intent(in) :: ztodt               ! 2 delta t (model time increment)
   real(r8), intent(in) :: landfrac(pcols)     ! Land fraction
      

   real(r8), intent(out) :: mcon(pcols,pverp)  ! Convective mass flux--m sub c
   real(r8), intent(out) :: dlf(pcols,pver)    ! scattrd version of the detraining cld h2o tend
   real(r8), intent(out) :: pflx(pcols,pverp)  ! scattered precip flux at each level
   real(r8), intent(out) :: cme(pcols,pver)    ! cmf condensation - evaporation
   real(r8), intent(out) :: zdu(pcols,pver)    ! detraining mass flux

   real(r8), intent(out) :: rliq(pcols) ! reserved liquid (not yet in cldliq) for energy integrals
   real(r8), intent(out):: mu(pcols,pver)
   real(r8), intent(out):: eu(pcols,pver) 
   real(r8), intent(out):: du(pcols,pver) 
   real(r8), intent(out):: md(pcols,pver) 
   real(r8), intent(out):: ed(pcols,pver) 
   real(r8), intent(out):: dp(pcols,pver) 
   
   ! wg layer thickness in mbs (between upper/lower interface).
   real(r8), intent(out):: dsubcld(pcols) 
   
   ! wg layer thickness in mbs between lcl and maxi.    
   integer, intent(out) :: jt(pcols)   
   
   ! wg top  level index of deep cumulus convection.
   integer, intent(out) :: maxg(pcols) 
   
   ! wg gathered values of maxi.
   integer, intent(out) :: ideep(pcols)
   
   ! w holds position of gathered points vs longitude index   
   integer, intent(out) :: lengath

   real(r8), pointer :: prec(:)   ! total precipitation
   real(r8), pointer :: snow(:)   ! snow from ZM convection 

   real(r8), pointer, dimension(:) :: jctop
   real(r8), pointer, dimension(:) :: jcbot
   real(r8), pointer, dimension(:,:,:) :: cld        
   real(r8), pointer, dimension(:,:) :: ql        ! wg grid slice of cloud liquid water.
   real(r8), pointer, dimension(:,:) :: rprd      ! rain production rate
   real(r8), pointer, dimension(:,:,:) :: fracis  ! fraction of transported species that are insoluble

   real(r8), pointer, dimension(:,:) :: evapcdp   ! Evaporation of deep convective precipitation

   real(r8), pointer :: pblh(:)                ! Planetary boundary layer height
   real(r8), pointer :: tpert(:)               ! Thermal temperature excess 

   ! Temperature tendency from deep convection (pbuf pointer).
   real(r8), pointer, dimension(:,:) :: ttend_dp

   real(r8) zero(pcols, pver)

   integer i, k

   call pbuf_get_field(pbuf, cldtop_idx,  jctop )
   call pbuf_get_field(pbuf, cldbot_idx,  jcbot )
   call pbuf_get_field(pbuf, icwmrdp_idx, ql    )

  select case ( deep_scheme )
  case('off', 'CLUBB_SGS' ) !    0 ==> no deep convection
    zero = 0     
    mcon = 0
    dlf = 0
    pflx = 0
    cme = 0
    zdu = 0
    rliq = 0

    call physics_ptend_init(ptend, state%psetcols, 'convect_deep')

!
! Associate pointers with physics buffer fields
!

    call pbuf_get_field(pbuf, cld_idx,         cld,    start=(/1,1/),   kount=(/pcols,pver/) ) 
    call pbuf_get_field(pbuf, rprddp_idx,      rprd )
    call pbuf_get_field(pbuf, fracis_idx,      fracis, start=(/1,1,1/), kount=(/pcols, pver, pcnst/) )
    call pbuf_get_field(pbuf, nevapr_dpcu_idx, evapcdp )
    call pbuf_get_field(pbuf, prec_dp_idx,     prec )
    call pbuf_get_field(pbuf, snow_dp_idx,     snow )

   !cld = 0  !!HuiWan 2014-05. Bugfix. cld is an input variable to zm_conv_tend.
    ql = 0
    rprd = 0
    fracis = 0
    evapcdp = 0
    prec=0
    snow=0

    jctop = pver
    jcbot = 1._r8

  case('ZM') !    1 ==> Zhang-McFarlane (default)
     call pbuf_get_field(pbuf, pblh_idx,  pblh)
     call pbuf_get_field(pbuf, tpert_idx, tpert)

     call t_startf('zm_conv_tend')
     call zm_conv_tend( pblh    ,mcon    ,cme     , &
          tpert   ,dlf     ,pflx    ,zdu      , &
          rliq    , &
          ztodt   , &
          jctop, jcbot , &
          state   ,ptend   ,landfrac, pbuf, mu, eu, &
          du, md, ed, dp, dsubcld, jt, maxg, ideep, lengath)
     call t_stopf('zm_conv_tend')


  end select

  ! If we added this, set it.
  if (ttend_dp_idx > 0) then
     call pbuf_get_field(pbuf, ttend_dp_idx, ttend_dp)
     if ( allocated(ptend%s) ) then
        ttend_dp(:state%ncol,:pver) = ptend%s(:state%ncol,:pver)/cpair
     else
        ttend_dp(:state%ncol,:pver) = 0.0_r8
     endif
  end if

  call outfld( 'ICWMRDP ', ql  , pcols, state%lchnk )


end subroutine convect_deep_tend
!=========================================================================================


subroutine convect_deep_tend_2( state,  ptend,  ztodt, pbuf, mu, eu, &
     du, md, ed, dp, dsubcld, jt, maxg, ideep, lengath, species_class)


   use physics_types, only: physics_state, physics_ptend, physics_ptend_init
   
   use physics_buffer,  only: physics_buffer_desc
   use constituents, only: pcnst
   use zm_conv_intr, only: zm_conv_tend_2

! Arguments
   type(physics_state), intent(in ) :: state          ! Physics state variables
   type(physics_ptend), intent(out) :: ptend          ! indivdual parameterization tendencies
   
   type(physics_buffer_desc), pointer :: pbuf(:)

   real(r8), intent(in) :: ztodt                          ! 2 delta t (model time increment)
   real(r8), intent(in):: mu(pcols,pver) 
   real(r8), intent(in):: eu(pcols,pver) 
   real(r8), intent(in):: du(pcols,pver) 
   real(r8), intent(in):: md(pcols,pver) 
   real(r8), intent(in):: ed(pcols,pver) 
   real(r8), intent(in):: dp(pcols,pver) 
   
   ! wg layer thickness in mbs (between upper/lower interface).
   real(r8), intent(in):: dsubcld(pcols) 
   
   ! wg layer thickness in mbs between lcl and maxi.    
   integer, intent(in) :: jt(pcols)   
   
   ! wg top  level index of deep cumulus convection.
   integer, intent(in) :: maxg(pcols) 
   
   ! wg gathered values of maxi.
   integer, intent(in) :: ideep(pcols)
   
   ! w holds position of gathered points vs longitude index   
   integer, intent(in) :: lengath

   integer, intent(in) :: species_class(:)

   if ( deep_scheme .eq. 'ZM' ) then  !    1 ==> Zhang-McFarlane (default)
      call zm_conv_tend_2( state,   ptend,  ztodt,  pbuf,mu, eu, &
     du, md, ed, dp, dsubcld, jt, maxg, ideep, lengath, species_class) 
   else
      call physics_ptend_init(ptend, state%psetcols, 'convect_deep2')
   end if


end subroutine convect_deep_tend_2


end module convect_deep
