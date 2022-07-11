
module cloud_diagnostics
!-------------------------------------------------------------------------------
! Purpose: Put cloud physical specifications on the history tape
!-------------------------------------------------------------------------------
  use shr_kind_mod,  only: r8=>shr_kind_r8
  use ppgrid,        only: begchunk, endchunk, pcols, pver,pverp
  use physconst,     only: gravit
  use cam_history,   only: outfld, addfld, horiz_only, add_default

  implicit none
  private
  save

  public :: cloud_diagnostics_register
  public :: cloud_diagnostics_init
  public :: cloud_diagnostics_calc
   
  ! Local variables
  integer :: ixcldice, ixcldliq             ! constituent indices
  integer :: iciwp_idx, iclwp_idx, cld_idx  ! cloud water path and fraction
  integer :: dei_idx, mu_idx, lambda_idx    ! condensate size/shape parameters
  integer :: rei_idx, rel_idx               ! effective radius of ice/liquid

  logical :: do_cld_diag
  integer :: conv_water_in_rad

  ! Index fields for precipitation efficiency.
  integer :: acpr_idx, acgcme_idx, acnum_idx

contains
!===============================================================================
!===============================================================================

subroutine cloud_diagnostics_register
  use physics_buffer,only: pbuf_add_field, dtype_r8, dtype_i4

  call pbuf_add_field('ICIWP', 'global', dtype_r8,(/pcols,pver/), iciwp_idx) ! In cloud ice water path for radiation
  call pbuf_add_field('ICLWP', 'global', dtype_r8,(/pcols,pver/), iclwp_idx) ! In cloud liquid water path for radiation

  return
end subroutine cloud_diagnostics_register

!===============================================================================
!===============================================================================
subroutine cloud_diagnostics_init()
  use physics_buffer,    only: pbuf_get_index, pbuf_get_field, physics_buffer_desc, pbuf_get_chunk
  use phys_control,      only: phys_getopts
  use physics_types,     only: physics_state
  use constituents,      only: cnst_get_ind
  use cloud_cover_diags, only: cloud_cover_diags_init
  implicit none
  !-----------------------------------------------------------------------------
  ! Local Variables
  type(physics_buffer_desc), pointer :: pbuf_chunk(:)
  character(len=16) :: wpunits, sampling_seq
  logical  :: history_verbose       ! produce verbose history output
  logical  :: use_MMF
  integer  :: i,k,lchnk
  integer  :: ncol 
  !-----------------------------------------------------------------------------
  call phys_getopts(use_MMF_out=use_MMF)
  cld_idx = pbuf_get_index('CLD')
  
  call phys_getopts( history_verbose_out = history_verbose )

  dei_idx    = pbuf_get_index('DEI')
  mu_idx     = pbuf_get_index('MU')
  lambda_idx = pbuf_get_index('LAMBDAC')
  rei_idx    = pbuf_get_index('REI')
  rel_idx    = pbuf_get_index('REL')

  call cnst_get_ind('CLDICE', ixcldice)
  call cnst_get_ind('CLDLIQ', ixcldliq)
  
  call phys_getopts(conv_water_in_rad_out=conv_water_in_rad)

  wpunits = 'kg/m2'
  sampling_seq=''
  
  call addfld('ICLDIWP', (/'lev'/),  'A',wpunits,'In-cloud ice water path'               , sampling_seq=sampling_seq)
  call addfld('ICLDTWP', (/'lev'/),  'A',wpunits,'In-cloud cloud total water path (liquid and ice)', sampling_seq=sampling_seq)
  call addfld('GCLDLWP',(/'lev'/),   'A',wpunits,'Grid-box cloud water path'             , sampling_seq=sampling_seq)
  call addfld('TGCLDCWP',horiz_only, 'A',wpunits,'Total grid-box cloud water path (liquid and ice)', sampling_seq=sampling_seq)
  call addfld('TGCLDLWP',horiz_only, 'A',wpunits,'Total grid-box cloud liquid water path', sampling_seq=sampling_seq)
  call addfld('TGCLDIWP',horiz_only, 'A',wpunits,'Total grid-box cloud ice water path'   , sampling_seq=sampling_seq)
  
  call addfld('lambda_cloud',(/ 'lev' /),'I','1/meter','lambda in cloud')
  call addfld('mu_cloud',(/ 'lev' /),'I','1','mu in cloud')
  call addfld('dei_cloud',(/ 'lev' /),'I','micrometers','ice radiative effective diameter in cloud')
  call addfld('rel_cloud',(/ 'lev' /),'I','1/meter','effective radius of liq in cloud', sampling_seq=sampling_seq)
  call addfld('rei_cloud',(/ 'lev' /),'I','1','effective radius of ice in cloud', sampling_seq=sampling_seq)
  
  call addfld('SETLWP',  (/'lev'/), 'A','gram/m2','Prescribed liquid water path' , sampling_seq=sampling_seq)
  call addfld('LWSH',    horiz_only,'A','m','Liquid water scale height'          , sampling_seq=sampling_seq)

  call addfld ('CLOUD',(/ 'lev' /), 'A','fraction','Cloud fraction' , sampling_seq=sampling_seq)

  call cloud_cover_diags_init(sampling_seq)

  call add_default ('CLOUD'   , 1, ' ')
  call add_default ('TGCLDLWP', 1, ' ')
  call add_default ('TGCLDIWP', 1, ' ')
  call add_default ('TGCLDCWP', 1, ' ')

  return
end subroutine cloud_diagnostics_init

!===============================================================================
!===============================================================================

subroutine cloud_diagnostics_calc(state,  pbuf)
  !-----------------------------------------------------------------------------
  ! Compute (liquid+ice) water path and cloud water/ice diagnostics
  !-----------------------------------------------------------------------------
  use physics_types,     only: physics_state    
  use physics_buffer,    only: physics_buffer_desc, pbuf_get_field, pbuf_old_tim_idx
  use phys_control,      only: phys_getopts
  use pkg_cldoptics,     only: cldovrlap, cldclw,  cldems
  use conv_water,        only: conv_water_4rad
  use radiation,         only: radiation_do
  use cloud_cover_diags, only: cloud_cover_diags_out
  use ref_pres,          only: top_lev=>trop_cloud_top_lev
  implicit none
  !-----------------------------------------------------------------------------
  ! Arguments
  type(physics_state), intent(in)    :: state        ! state variables
  type(physics_buffer_desc), pointer :: pbuf(:)
  !-----------------------------------------------------------------------------
  ! Local variables
  real(r8), pointer :: cld(:,:)       ! cloud fraction
  real(r8), pointer :: iciwp(:,:)     ! in-cloud cloud ice water path
  real(r8), pointer :: iclwp(:,:)     ! in-cloud cloud liquid water path
  real(r8), pointer :: dei(:,:)       ! effective radiative diameter of ice
  real(r8), pointer :: mu(:,:)        ! gamma distribution for liq clouds
  real(r8), pointer :: lambda(:,:)    ! gamma distribution for liq clouds
  real(r8), pointer :: rei(:,:)       ! effective radiative radius of ice
  real(r8), pointer :: rel(:,:)       ! effective radiative radius of liq

  real(r8), pointer :: cldemis(:,:)   ! cloud emissivity
  real(r8), pointer :: cldtau(:,:)    ! cloud optical depth

  integer,  pointer :: nmxrgn(:)      ! Number of maximally overlapped regions
  real(r8), pointer :: pmxrgn(:,:)    ! Maximum values of pressure for each

  integer  :: itim_old

  real(r8) :: ficemr (pcols,pver)     ! Ice fraction from ice and liquid mixing ratios

  real(r8) :: icimr(pcols,pver)       ! In cloud ice mixing ratio
  real(r8) :: icwmr(pcols,pver)       ! In cloud water mixing ratio
  real(r8) :: iwc(pcols,pver)         ! Grid box average ice water content
  real(r8) :: lwc(pcols,pver)         ! Grid box average liquid water content

  ! old data
  real(r8) :: tpw    (pcols)          ! total precipitable water
  real(r8) :: clwpold(pcols,pver)     ! Presribed cloud liq. h2o path
  real(r8) :: hl     (pcols)          ! Liquid water scale height

  integer  :: i,k                      ! loop indexes
  integer  :: ncol, lchnk
  real(r8) :: rgrav

  real(r8) :: allcld_ice (pcols,pver) ! Convective cloud ice
  real(r8) :: allcld_liq (pcols,pver) ! Convective cloud liquid

  logical  :: dosw,dolw
  logical  :: use_MMF

  !-----------------------------------------------------------------------------
  if (.not.do_cld_diag) return

  call phys_getopts(use_MMF_out=use_MMF)

  dosw     = .true.
  dolw     = .true.

  if (.not.(dosw .or. dolw)) return

  ncol  = state%ncol
  lchnk = state%lchnk

  itim_old = pbuf_old_tim_idx()
  call pbuf_get_field(pbuf, cld_idx, cld, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

  call pbuf_get_field(pbuf, iclwp_idx, iclwp )
  call pbuf_get_field(pbuf, iciwp_idx, iciwp )
  call pbuf_get_field(pbuf, dei_idx, dei )
  call pbuf_get_field(pbuf, mu_idx, mu )
  call pbuf_get_field(pbuf, lambda_idx, lambda )
  call pbuf_get_field(pbuf, rei_idx, rei )
  call pbuf_get_field(pbuf, rel_idx, rel )

  call outfld('dei_cloud',dei(:,:),pcols,lchnk)
  call outfld('mu_cloud',mu(:,:),pcols,lchnk)
  call outfld('lambda_cloud',lambda(:,:),pcols,lchnk)
  call outfld('rel_cloud', rel, pcols, lchnk)
  call outfld('rei_cloud', rei, pcols, lchnk)

  ! Compute liquid and ice water paths

  ! ----------------------------------------------------------- !
  ! Adjust in-cloud water values to take account of convective  !
  ! in-cloud water. It is used to calculate the values of       !
  ! iclwp and iciwp to pass to the radiation.                   !
  ! ----------------------------------------------------------- !
  if( conv_water_in_rad /= 0 ) then
    allcld_ice(:ncol,:) = 0._r8 ! Grid-avg all cloud liquid
    allcld_liq(:ncol,:) = 0._r8 ! Grid-avg all cloud ice

    call conv_water_4rad( state, pbuf, conv_water_in_rad, allcld_liq, allcld_ice )
  else
    allcld_liq(:ncol,top_lev:pver) = state%q(:ncol,top_lev:pver,ixcldliq)  ! Grid-ave all cloud liquid
    allcld_ice(:ncol,top_lev:pver) = state%q(:ncol,top_lev:pver,ixcldice)  ! Grid-ave all cloud ice
  end if

  ! ------------------------------------------------------------ !
  ! Compute in cloud ice and liquid mixing ratios                !
  ! Note that 'iclwp, iciwp' are used for radiation computation. !
  ! ------------------------------------------------------------ !
  do k = 1,pver
    do i = 1,ncol
      iciwp(i,k)  = 0._r8
      iclwp(i,k)  = 0._r8
      icimr(i,k)  = 0._r8
      icwmr(i,k)  = 0._r8
      iwc(i,k)    = 0._r8
      lwc(i,k)    = 0._r8
    end do
  end do

  do k = top_lev, pver
    do i = 1, ncol
      ! Limits for in-cloud mixing ratios consistent with MG microphysics
      ! in-cloud mixing ratio maximum limit of 0.005 kg/kg
      icimr(i,k)     = min( allcld_ice(i,k) / max(0.0001_r8,cld(i,k)),0.005_r8 )
      icwmr(i,k)     = min( allcld_liq(i,k) / max(0.0001_r8,cld(i,k)),0.005_r8 )
      iwc(i,k)       = allcld_ice(i,k) * state%pmid(i,k) / (287.15_r8*state%t(i,k))
      lwc(i,k)       = allcld_liq(i,k) * state%pmid(i,k) / (287.15_r8*state%t(i,k))
      ! Calculate total cloud water paths in each layer
      iciwp(i,k)     = icimr(i,k) * state%pdel(i,k) / gravit
      iclwp(i,k)     = icwmr(i,k) * state%pdel(i,k) / gravit
    end do
  end do

  ! Determine parameters for maximum/random overlap
  call cldovrlap(lchnk, ncol, state%pint, cld, nmxrgn, pmxrgn)

  ! Cloud cover diagnostics
  call cloud_cover_diags_out(lchnk, ncol, cld, state%pmid, nmxrgn, pmxrgn )

  ! --------------------------------------------- !
  ! General outfield calls for microphysics       !
  ! --------------------------------------------- !
  call outfld( 'IWC'  , iwc,         pcols, lchnk )
  call outfld( 'LWC'  , lwc,         pcols, lchnk )
  call outfld( 'ICIMR', icimr,       pcols, lchnk )
  call outfld( 'ICWMR', icwmr,       pcols, lchnk )

  ! Compute total preciptable water in column (in mm)
  tpw(:ncol) = 0.0_r8
  rgrav = 1.0_r8/gravit
  do k=1,pver
    do i=1,ncol
      tpw(i) = tpw(i) + state%pdel(i,k)*state%q(i,k,1)*rgrav
    end do
  end do

  ! Diagnostic liquid water path (old specified form)
  call cldclw(lchnk, ncol, state%zi, clwpold, tpw, hl)
  call outfld('SETLWP'  ,clwpold, pcols,lchnk)
  call outfld('LWSH'    ,hl     , pcols,lchnk)

  return
end subroutine cloud_diagnostics_calc

!===============================================================================
!===============================================================================
end module cloud_diagnostics
