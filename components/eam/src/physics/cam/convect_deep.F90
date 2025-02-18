module convect_deep
   !----------------------------------------------------------------------------
   ! Deep convection parameterization interface
   ! Zhang-McFarlane is the only supported scheme
   !
   ! Author: D.B. Coleman, Sep 2004
   !----------------------------------------------------------------------------
   use shr_kind_mod, only: r8=>shr_kind_r8
   use ppgrid,       only: pver, pcols, pverp
   use cam_logfile,  only: iulog
   use spmd_utils,   only: masterproc
   use perf_mod,     only: t_startf, t_stopf

   implicit none
   save
   private

   ! public methods
   public :: convect_deep_register        ! register fields in physics buffer
   public :: convect_deep_init            ! initialize donner_deep module
   public :: convect_deep_tend            ! call the deep convection scheme
   public :: convect_deep_tend_2          ! convective transport of constituents
   public :: deep_scheme_does_scav_trans  ! indicate if deep scheme does scavenging and transport
   
   ! private module data
   character(len=16) :: deep_scheme ! name of deep scheme from namelist

   ! physics buffer field indices
   integer :: icwmrdp_idx      = 0  ! deep convective in cloud water mixing ratio            [kg/m2]
   integer :: icimrdp_idx      = 0  ! deep convective in-cloud ice water content             [kg/m2]
   integer :: rprddp_idx       = 0  ! deep convective rain production rate                   [?]
   integer :: nevapr_dpcu_idx  = 0  ! deep convective evaporation of precip                  [?]
   integer :: cldtop_idx       = 0  ! level index for top of deep convection
   integer :: cldbot_idx       = 0  ! level index for bottom of deep convection
   integer :: fracis_idx       = 0  ! fraction of transported species that are insoluble
   integer :: pblh_idx         = 0  ! planetary boundary layer height
   integer :: tpert_idx        = 0  ! thermal temperature excess                             [K]
   integer :: prec_dp_idx      = 0  ! total surface precip from deep convection              [m/s]
   integer :: snow_dp_idx      = 0  ! frozen surface precip from deep convection             [m/s]
   integer :: dp_cldliq_idx    = 0  ! deep convective cloud liquid water                     [kg/kg]
   integer :: dp_cldice_idx    = 0  ! deep convective cloud liquid water                     [kg/kg]
   integer :: dp_flxprc_idx    = 0  ! deep convective flux of precipitation                  [kg/m2/s]
   integer :: dp_flxsnw_idx    = 0  ! deep convective flux of snow                           [kg/m2/s]
   integer :: dp_frac_idx      = 0  ! deep convection cloud fraction
   integer :: lambdadpcu_idx   = 0  ! droplet size distribution shape parameter for radiation
   integer :: mudpcu_idx       = 0  ! droplet size distribution shape parameter for radiation
   integer :: ttend_dp_idx     = 0  ! convective heating for convective gravity wave scheme  [K/s]

contains 

!===================================================================================================

function deep_scheme_does_scav_trans()
   !----------------------------------------------------------------------------
   ! Function called by tphysbc to determine if it needs to do scavenging and 
   ! convective transport or if those have been done by the deep convection
   !----------------------------------------------------------------------------
   logical :: deep_scheme_does_scav_trans
   deep_scheme_does_scav_trans = .false.
   return
end function deep_scheme_does_scav_trans

!===================================================================================================

subroutine convect_deep_register()
   !----------------------------------------------------------------------------
   ! Purpose: register fields with the physics buffer
   !----------------------------------------------------------------------------
   use phys_control,   only: phys_getopts, use_gw_convect
   use physics_buffer, only: pbuf_add_field, dtype_r8
   use zm_conv_intr,   only: zm_conv_register
   !----------------------------------------------------------------------------
   integer idx

   ! get deep_scheme setting from phys_control
   call phys_getopts(deep_scheme_out = deep_scheme)

   select case ( deep_scheme )
   case('ZM') ! Zhang-McFarlane
      call zm_conv_register
   end select 

   ! Add PBUF variables related to deep convection
   call pbuf_add_field('DP_FLXPRC',  'global',  dtype_r8, (/pcols,pverp/),dp_flxprc_idx)
   call pbuf_add_field('DP_FLXSNW',  'global',  dtype_r8, (/pcols,pverp/),dp_flxsnw_idx)
   call pbuf_add_field('DP_CLDLIQ',  'global',  dtype_r8, (/pcols,pver/), dp_cldliq_idx)
   call pbuf_add_field('DP_CLDICE',  'global',  dtype_r8, (/pcols,pver/), dp_cldice_idx) 
   call pbuf_add_field('ICWMRDP',    'physpkg', dtype_r8, (/pcols,pver/), icwmrdp_idx)
   call pbuf_add_field('ICIMRDP',    'physpkg', dtype_r8, (/pcols,pver/), icimrdp_idx)
   call pbuf_add_field('RPRDDP',     'physpkg', dtype_r8, (/pcols,pver/), rprddp_idx)
   call pbuf_add_field('NEVAPR_DPCU','physpkg', dtype_r8, (/pcols,pver/), nevapr_dpcu_idx)
   call pbuf_add_field('PREC_DP',    'physpkg', dtype_r8, (/pcols/),      prec_dp_idx)
   call pbuf_add_field('SNOW_DP',    'physpkg', dtype_r8, (/pcols/),      snow_dp_idx)
   call pbuf_add_field('LAMBDADPCU', 'physpkg', dtype_r8, (/pcols,pver/), lambdadpcu_idx)
   call pbuf_add_field('MUDPCU',     'physpkg', dtype_r8, (/pcols,pver/), mudpcu_idx)

   if (use_gw_convect) then
      call pbuf_add_field('TTEND_DP','physpkg', dtype_r8, (/pcols,pver/), ttend_dp_idx)
   end if

end subroutine convect_deep_register

!===================================================================================================

subroutine convect_deep_init(pref_edge,pbuf2d)
   !----------------------------------------------------------------------------
   ! Purpose: declare output fields, initialize variables needed by convection
   !----------------------------------------------------------------------------
   use cam_history,     only: addfld
   use pmgrid,          only: plevp
   use zm_conv_intr,    only: zm_conv_init
   use cam_abortutils,  only: endrun
   use physics_buffer,  only: physics_buffer_desc, pbuf_get_index, pbuf_set_field
   !----------------------------------------------------------------------------
   ! Arguments
   real(r8), dimension(plevp), intent(in) :: pref_edge   ! reference pressures at interfaces
   type(physics_buffer_desc),  pointer    :: pbuf2d(:,:) ! physics buffer
   !----------------------------------------------------------------------------
   dp_frac_idx = pbuf_get_index('DP_FRAC')
   icwmrdp_idx = pbuf_get_index('ICWMRDP')
   icimrdp_idx = pbuf_get_index('ICIMRDP')

   select case ( deep_scheme )
   case('off') ! no deep convection
      if (masterproc) write(iulog,*)'convect_deep: no deep convection selected'
      call pbuf_set_field(pbuf2d, dp_cldliq_idx,   0._r8)
      call pbuf_set_field(pbuf2d, dp_cldice_idx,   0._r8)
      call pbuf_set_field(pbuf2d, dp_flxprc_idx,   0._r8)
      call pbuf_set_field(pbuf2d, dp_flxsnw_idx,   0._r8)
      call pbuf_set_field(pbuf2d, dp_frac_idx,     0._r8)
      call pbuf_set_field(pbuf2d, icwmrdp_idx,     0._r8)
      call pbuf_set_field(pbuf2d, icimrdp_idx,     0._r8)
      call pbuf_set_field(pbuf2d, rprddp_idx,      0._r8)
      call pbuf_set_field(pbuf2d, nevapr_dpcu_idx, 0._r8)
   case('CLUBB_SGS')
      if (masterproc) write(iulog,*)'convect_deep: CLUBB_SGS selected'
   case('ZM') ! Zhang-McFarlane
      if (masterproc) write(iulog,*)'convect_deep initializing Zhang-McFarlane convection'
      call zm_conv_init(pref_edge)
   case default
      if (masterproc) write(iulog,*)'WARNING: convect_deep: no deep convection scheme'
   end select

   cldtop_idx = pbuf_get_index('CLDTOP')
   cldbot_idx = pbuf_get_index('CLDBOT')
   fracis_idx = pbuf_get_index('FRACIS')
   pblh_idx   = pbuf_get_index('pblh')
   tpert_idx  = pbuf_get_index('tpert')

   call addfld('ICWMRDP', (/'lev'/), 'A', 'kg/kg', 'Deep Convection in-cloud water mixing ratio ')

end subroutine convect_deep_init

!===================================================================================================
subroutine convect_deep_tend( mcon, cme, dlf, pflx, zdu, rliq, rice, ztodt, state, ptend, &
                              landfrac, pbuf, mu, eu, du, md, ed, dp, dsubcld, jt, maxg, &
                              ideep,lengath ) 
   !----------------------------------------------------------------------------
   ! Purpose: call the deep convection scheme (or zero output fields if no deep scheme)
   !----------------------------------------------------------------------------
   use physics_types,  only: physics_state, physics_ptend, physics_tend, physics_ptend_init
   use physics_buffer, only: physics_buffer_desc, pbuf_get_field
   use zm_conv_intr,   only: zm_conv_tend
   use cam_history,    only: outfld
   use constituents,   only: pcnst
   use physconst,      only: cpair
   !----------------------------------------------------------------------------
   ! Arguments
   type(physics_state),             intent(in ) :: state    ! Physics state variables
   type(physics_ptend),             intent(out) :: ptend    ! individual parameterization tendencies
   type(physics_buffer_desc),       pointer     :: pbuf(:)  ! physics buffer
   real(r8),                        intent(in)  :: ztodt    ! 2 delta t (model time increment)
   real(r8), dimension(pcols),      intent(in)  :: landfrac ! Land fraction
   real(r8), dimension(pcols,pverp),intent(out) :: mcon     ! Convective mass flux--m sub c
   real(r8), dimension(pcols,pver), intent(out) :: dlf      ! scattered detraining cld h2o tend
   real(r8), dimension(pcols,pverp),intent(out) :: pflx     ! scattered precip flux at each level
   real(r8), dimension(pcols,pver), intent(out) :: cme      ! mass tendency due to condensation - evaporation
   real(r8), dimension(pcols,pver), intent(out) :: zdu      ! detraining mass flux
   real(r8), dimension(pcols),      intent(out) :: rliq     ! reserved liq (not yet in cldliq) for energy integrals
   real(r8), dimension(pcols),      intent(out) :: rice     ! reserved ice (not yet in cldice) for energy integrals
   real(r8), dimension(pcols,pver), intent(out) :: mu       ! upward cloud mass flux               [?]
   real(r8), dimension(pcols,pver), intent(out) :: eu       ! entrainment in updraft               [?]
   real(r8), dimension(pcols,pver), intent(out) :: du       ! detrainment in updraft               [?]
   real(r8), dimension(pcols,pver), intent(out) :: md       ! downward cloud mass flux             [?]
   real(r8), dimension(pcols,pver), intent(out) :: ed       ! entrainment in downdraft             [?]
   real(r8), dimension(pcols,pver), intent(out) :: dp       ! layer thickness between interfaces   [mb]
   real(r8), dimension(pcols),      intent(out) :: dsubcld  ! layer thickness between lcl and maxi [mb]
   integer,  dimension(pcols),      intent(out) :: jt       ! top level index of deep convection
   integer,  dimension(pcols),      intent(out) :: maxg     ! gathered values of maxi
   integer,  dimension(pcols),      intent(out) :: ideep    ! flag to indicate active column for gathering
   integer,                         intent(out) :: lengath  ! number of columns in gathered arrays
   !----------------------------------------------------------------------------
   ! Local variables
   integer,           dimension(pcols) :: jctop       ! integer version of cloud top level index
   integer,           dimension(pcols) :: jcbot       ! integer version of cloud bot level index
   real(r8), pointer, dimension(:)     :: jctop_r8    ! real version of cloud top level index for pbuf
   real(r8), pointer, dimension(:)     :: jcbot_r8    ! real version of cloud bot level index for pbuf
   real(r8), pointer, dimension(:)     :: prec        ! total precip at surface
   real(r8), pointer, dimension(:)     :: snow        ! snow at surface
   real(r8), pointer, dimension(:,:)   :: ql          ! grid slice of cloud liquid water
   real(r8), pointer, dimension(:,:)   :: qi          ! wg grid slice of cloud ice
   real(r8), pointer, dimension(:,:)   :: rprd        ! rain production rate
   real(r8), pointer, dimension(:,:)   :: evapcdp     ! evaporation of deep convective precip
   real(r8), pointer, dimension(:)     :: pblh        ! planetary boundary layer height
   real(r8), pointer, dimension(:)     :: tpert       ! thermal temperature excess 
   real(r8), pointer, dimension(:,:,:) :: fracis      ! fraction of transported species that are insoluble
   real(r8), pointer, dimension(:,:)   :: mudpcu      ! Droplet size distribution shape parameter for radiation
   real(r8), pointer, dimension(:,:)   :: lambdadpcu  ! Droplet size distribution shape parameter for radiation
   real(r8), pointer, dimension(:,:)   :: ttend_dp    ! temperature tendency from deep convection
   integer :: i, k ! loop iterators
   !----------------------------------------------------------------------------

   ! Associate pointers with physics buffer fields
   call pbuf_get_field(pbuf, cldtop_idx,  jctop_r8 )
   call pbuf_get_field(pbuf, cldbot_idx,  jcbot_r8 )
   call pbuf_get_field(pbuf, icwmrdp_idx, ql    )

  select case ( deep_scheme )
     case('off', 'CLUBB_SGS' ) ! no deep convection  
         ! initialize physics tendencies
         call physics_ptend_init(ptend, state%psetcols, 'convect_deep')

         ! Associate pointers with physics buffer fields
         call pbuf_get_field(pbuf, fracis_idx,      fracis,    start=(/1,1,1/), kount=(/pcols,pver,pcnst/) )
         call pbuf_get_field(pbuf, icimrdp_idx,     qi         )
         call pbuf_get_field(pbuf, rprddp_idx,      rprd       )
         call pbuf_get_field(pbuf, nevapr_dpcu_idx, evapcdp    )
         call pbuf_get_field(pbuf, prec_dp_idx,     prec       )
         call pbuf_get_field(pbuf, snow_dp_idx,     snow       )
         call pbuf_get_field(pbuf, mudpcu_idx,      mudpcu     )
         call pbuf_get_field(pbuf, lambdadpcu_idx,  lambdadpcu )

         ! initialize pbuf variables to zero
         fracis     = 0
         ql         = 0
         qi         = 0
         rprd       = 0
         evapcdp    = 0
         prec       = 0
         snow       = 0
         mudpcu     = 0
         lambdadpcu = 0
         
         ! initialize output variables to zero
         mcon       = 0
         dlf        = 0
         pflx       = 0
         cme        = 0
         zdu        = 0
         rliq       = 0
         rice       = 0

         ! set cloud bot/top level indices to encompass all levels
         jctop = pver
         jcbot = 1

      case('ZM') ! Zhang-McFarlane
         ! Associate pointers with physics buffer fields
         call pbuf_get_field(pbuf, pblh_idx,  pblh)
         call pbuf_get_field(pbuf, tpert_idx, tpert)
         ! run the deep convection scheme
         call t_startf('zm_conv_tend')
         call zm_conv_tend( pblh, mcon, cme, tpert, dlf, pflx, zdu, rliq, rice, &
                            ztodt, jctop, jcbot, state, ptend, landfrac, pbuf, &
                            mu, eu, du, md, ed, dp, dsubcld, jt, maxg, ideep, lengath )
         call t_stopf('zm_conv_tend')
   end select

   ! set the real version of cloud bot/top level indices
   do i = 1,pcols
      jctop_r8(i) = real(jctop(i), r8)
      jcbot_r8(i) = real(jcbot(i), r8)
   end do

   ! set the deep convective tendency if it was added to pbuf
   if (ttend_dp_idx > 0) then
      call pbuf_get_field(pbuf, ttend_dp_idx, ttend_dp)
      if ( allocated(ptend%s) ) then
         ttend_dp(1:state%ncol,1:pver) = ptend%s(:state%ncol,:pver)/cpair
      else
         ttend_dp(1:state%ncol,1:pver) = 0.0_r8
      endif
   end if

   ! history file output
   call outfld( 'ICWMRDP ', ql, pcols, state%lchnk )

end subroutine convect_deep_tend

!===================================================================================================

subroutine convect_deep_tend_2( state,  ptend,  ztodt, pbuf, mu, eu, du, md, ed, dp, &
                                dsubcld, jt, maxg, ideep, lengath, species_class)
   !----------------------------------------------------------------------------
   ! Purpose: convective transport of constituents not handled by convect_deep_tend
   !----------------------------------------------------------------------------
   use physics_types,   only: physics_state, physics_ptend, physics_ptend_init
   use physics_buffer,  only: physics_buffer_desc
   use constituents,    only: pcnst
   use zm_conv_intr,    only: zm_conv_tend_2
   !----------------------------------------------------------------------------
   ! Arguments
   type(physics_state),             intent(in ) :: state          ! physics state
   type(physics_ptend),             intent(out) :: ptend          ! output tendencies
   real(r8),                        intent(in ) :: ztodt          ! 2 delta t (model time increment)
   type(physics_buffer_desc),       pointer     :: pbuf(:)        ! physics buffer
   real(r8), dimension(pcols,pver), intent(in ) :: mu             ! upward cloud mass flux               [?]
   real(r8), dimension(pcols,pver), intent(in ) :: eu             ! entrainment in updraft               [?]
   real(r8), dimension(pcols,pver), intent(in ) :: du             ! detrainment in updraft               [?]
   real(r8), dimension(pcols,pver), intent(in ) :: md             ! downward cloud mass flux             [?]
   real(r8), dimension(pcols,pver), intent(in ) :: ed             ! entrainment in downdraft             [?]
   real(r8), dimension(pcols,pver), intent(in ) :: dp             ! layer thickness between interfaces   [mb]
   real(r8), dimension(pcols),      intent(in ) :: dsubcld        ! layer thickness between lcl and maxi [mb]
   integer,  dimension(pcols),      intent(in ) :: jt             ! top level index of deep convection
   integer,  dimension(pcols),      intent(in ) :: maxg           ! gathered values of maxi
   integer,  dimension(pcols),      intent(in ) :: ideep          ! flag to indicate active column for gathering
   integer,                         intent(in ) :: lengath        ! number of columns in gathered arrays
   integer,  dimension(:),          intent(in ) :: species_class  ! constituent tracer type
   !----------------------------------------------------------------------------
   if ( deep_scheme .eq. 'ZM' ) then ! Zhang-McFarlane
      call zm_conv_tend_2( state, ptend, ztodt,  pbuf, mu, eu, du, md, ed, dp, &
                           dsubcld, jt, maxg, ideep, lengath, species_class)
   else ! no deep convection  
      call physics_ptend_init(ptend, state%psetcols, 'convect_deep2')
   end if

end subroutine convect_deep_tend_2

!===================================================================================================

end module convect_deep
