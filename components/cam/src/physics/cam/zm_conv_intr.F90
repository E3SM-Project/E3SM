


module zm_conv_intr
!---------------------------------------------------------------------------------
! Purpose:
!
! CAM interface to the Zhang-McFarlane deep convection scheme
!
! Author: D.B. Coleman
! January 2010 modified by J. Kay to add COSP simulator fields to physics buffer
! July 2015 B. Singh Added code for unified convective trasport
!---------------------------------------------------------------------------------
   use shr_kind_mod, only: r8=>shr_kind_r8
   use physconst,    only: cpair                              
   use physconst,    only: latvap, gravit   !songxl 2014-05-20
   use ppgrid,       only: pver, pcols, pverp, begchunk, endchunk
   use zm_conv,      only: zm_conv_evap, zm_convr, convtran, momtran, trigmem
   use cam_history,  only: outfld, addfld, horiz_only, add_default
   use perf_mod
   use cam_logfile,  only: iulog
   
   implicit none
   private
   save

   ! Public methods

   public ::&
      zm_conv_register,           &! register fields in physics buffer
      zm_conv_init,               &! initialize donner_deep module
      zm_conv_tend,               &! return tendencies
      zm_conv_tend_2               ! return tendencies

   ! Private module data
   integer ::& ! indices for fields in the physics buffer
      dp_flxprc_idx, &
      dp_flxsnw_idx, &
      dp_cldliq_idx, &
      dp_cldice_idx, &
      prec_dp_idx,   &
      snow_dp_idx

!<songxl 2014-05-20------------------
   integer :: hu_nm1_idx    !hu_nm1 index in physics buffer
   integer :: cnv_nm1_idx   !cnv_nm1 index in physics buffer
   integer :: tm1_idx,     &!tm1 index in physics buffer
              qm1_idx       !qm1 index in physics buffer
!>songxl 2014-05-20------------------

!  indices for fields in the physics buffer
   integer  ::    cld_idx          = 0    
   integer  ::    icwmrdp_idx      = 0     
   integer  ::    rprddp_idx       = 0    
   integer  ::    fracis_idx       = 0   
   integer  ::    nevapr_dpcu_idx  = 0    
   logical  ::    convproc_do_aer 
   logical  ::    convproc_do_gas 
   logical  ::    clim_modal_aero

!=========================================================================================
contains
!=========================================================================================

subroutine zm_conv_register

!----------------------------------------
! Purpose: register fields with the physics buffer
!----------------------------------------

  use physics_buffer, only : pbuf_add_field, dtype_r8

  implicit none

  integer idx

! Flux of precipitation from deep convection (kg/m2/s)
   call pbuf_add_field('DP_FLXPRC','global',dtype_r8,(/pcols,pverp/),dp_flxprc_idx) 

! Flux of snow from deep convection (kg/m2/s) 
   call pbuf_add_field('DP_FLXSNW','global',dtype_r8,(/pcols,pverp/),dp_flxsnw_idx) 

! deep gbm cloud liquid water (kg/kg)
   call pbuf_add_field('DP_CLDLIQ','global',dtype_r8,(/pcols,pver/), dp_cldliq_idx)  

! deep gbm cloud liquid water (kg/kg)    
   call pbuf_add_field('DP_CLDICE','global',dtype_r8,(/pcols,pver/), dp_cldice_idx)  

!<songxl 2014-05-20-------------
!  if(trigmem)then
! moist static energy at n-1 time step (J/kg)
    call pbuf_add_field('HU_NM1','global',dtype_r8,(/pcols,pver/), hu_nm1_idx)
! moist convection index at n-1 time step(1=yes, 0=no)
    call pbuf_add_field('CNV_NM1','global',dtype_r8,(/pcols,pver/), cnv_nm1_idx)
! temperature at n-1 time step
    call pbuf_add_field('TM1', 'global', dtype_r8,(/pcols,pver/), tm1_idx)
! specific humidity at n-1 time step
    call pbuf_add_field('QM1', 'global', dtype_r8,(/pcols,pver/), qm1_idx) 
!  endif
!>songxl 2014-05-20-------------

end subroutine zm_conv_register

!=========================================================================================

subroutine zm_conv_init(pref_edge)

!----------------------------------------
! Purpose:  declare output fields, initialize variables needed by convection
!----------------------------------------

  use cam_history,    only: outfld, addfld, horiz_only, add_default
  use ppgrid,         only: pcols, pver
  use zm_conv,        only: zm_convi
  use pmgrid,         only: plev,plevp
  use spmd_utils,     only: masterproc
  use error_messages, only: alloc_err	
  use phys_control,   only: phys_deepconv_pbl, phys_getopts, cam_physpkg_is
  use physics_buffer, only: pbuf_get_index
  use rad_constituents, only: rad_cnst_get_info 

  implicit none

  real(r8),intent(in) :: pref_edge(plevp)        ! reference pressures at interfaces


  logical :: no_deep_pbl    ! if true, no deep convection in PBL
  integer  limcnv           ! top interface level limit for convection
  integer k, istat
  logical :: history_budget ! output tendencies and state variables for CAM4
                            ! temperature, water vapor, cloud ice and cloud
                            ! liquid budgets.
  integer :: history_budget_histfile_num ! output history file number for budget fields
  integer :: nmodes 

! 
! Register fields with the output buffer
!


    call addfld ('PRECZ',horiz_only,    'A','m/s','total precipitation from ZM convection')
    call addfld ('ZMDT',(/ 'lev' /), 'A','K/s','T tendency - Zhang-McFarlane moist convection')
    call addfld ('ZMDQ',(/ 'lev' /), 'A','kg/kg/s','Q tendency - Zhang-McFarlane moist convection')
    call addfld ('ZMDICE',(/ 'lev' /), 'A','kg/kg/s','Cloud ice tendency - Zhang-McFarlane convection')
    call addfld ('ZMDLIQ',(/ 'lev' /), 'A','kg/kg/s','Cloud liq tendency - Zhang-McFarlane convection')
    call addfld ('EVAPTZM',(/ 'lev' /), 'A','K/s','T tendency - Evaporation/snow prod from Zhang convection')
    call addfld ('FZSNTZM',(/ 'lev' /), 'A','K/s','T tendency - Rain to snow conversion from Zhang convection')
    call addfld ('EVSNTZM',(/ 'lev' /), 'A','K/s','T tendency - Snow to rain prod from Zhang convection')
    call addfld ('EVAPQZM',(/ 'lev' /), 'A','kg/kg/s','Q tendency - Evaporation from Zhang-McFarlane moist convection')
    
    call addfld ('ZMFLXPRC',(/ 'ilev' /), 'A','kg/m2/s','Flux of precipitation from ZM convection'       )
    call addfld ('ZMFLXSNW',(/ 'ilev' /), 'A','kg/m2/s','Flux of snow from ZM convection'                )
    call addfld ('ZMNTPRPD',(/ 'lev' /) , 'A','kg/kg/s','Net precipitation production from ZM convection')
    call addfld ('ZMNTSNPD',(/ 'lev' /) , 'A','kg/kg/s','Net snow production from ZM convection'         )
    call addfld ('ZMEIHEAT',(/ 'lev' /) , 'A','W/kg'    ,'Heating by ice and evaporation in ZM convection')
    
    call addfld ('CMFMCDZM',(/ 'ilev' /),'A','kg/m2/s','Convection mass flux from ZM deep ')
    call addfld ('PRECCDZM',horiz_only,    'A','m/s','Convective precipitation rate from ZM deep')

    call addfld ('PCONVB',horiz_only , 'A','Pa'    ,'convection base pressure')
    call addfld ('PCONVT',horiz_only , 'A','Pa'    ,'convection top  pressure')

    call addfld ('CAPE',       horiz_only, 'A',   'J/kg', 'Convectively available potential energy')
    call addfld ('FREQZM',horiz_only  ,'A','fraction', 'Fractional occurance of ZM convection') 
    ! To keep track of when ZM is suppressed due to Drizzle issue
    call addfld ('FREQZMSUP',horiz_only  ,'A','fraction', 'Fractional occurrence of suppressed ZM convection')
    ! To keep track of how high deep convection shoots up to
    call addfld ('MAXZM_IDX',horiz_only  ,'A','index', 'Index of highest grid box ZM convection reaches')

    call addfld ('ZMMTT',     (/ 'lev' /), 'A', 'K/s', 'T tendency - ZM convective momentum transport')
    call addfld ('ZMMTU',    (/ 'lev' /), 'A',  'm/s2', 'U tendency - ZM convective momentum transport')
    call addfld ('ZMMTV',    (/ 'lev' /), 'A',  'm/s2', 'V tendency - ZM convective momentum transport')

    call addfld ('ZMMU', (/ 'lev' /), 'A',   'kg/m2/s', 'ZM convection updraft mass flux')
    call addfld ('ZMMD', (/ 'lev' /), 'A',   'kg/m2/s', 'ZM convection downdraft mass flux')

    call addfld ('ZMUPGU',    (/ 'lev' /), 'A', 'm/s2', 'zonal force from ZM updraft pressure gradient term')
    call addfld ('ZMUPGD',    (/ 'lev' /), 'A', 'm/s2', 'zonal force from ZM downdraft pressure gradient term')
    call addfld ('ZMVPGU',    (/ 'lev' /), 'A', 'm/s2', 'meridional force from ZM updraft pressure gradient term')
    call addfld ('ZMVPGD',    (/ 'lev' /), 'A', 'm/s2', 'merdional force from ZM downdraft pressure gradient term')

    call addfld ('ZMICUU',     (/ 'lev' /), 'A', 'm/s', 'ZM in-cloud U updrafts')
    call addfld ('ZMICUD',     (/ 'lev' /), 'A', 'm/s', 'ZM in-cloud U downdrafts')
    call addfld ('ZMICVU',     (/ 'lev' /), 'A', 'm/s', 'ZM in-cloud V updrafts')
    call addfld ('ZMICVD',     (/ 'lev' /), 'A', 'm/s', 'ZM in-cloud V downdrafts')
    
    call phys_getopts( history_budget_out = history_budget, &
                       history_budget_histfile_num_out = history_budget_histfile_num, &
                       convproc_do_aer_out = convproc_do_aer, & 
                       convproc_do_gas_out = convproc_do_gas)   
    ! Determine whether its a 'modal' aerosol simulation  or not
    call rad_cnst_get_info(0, nmodes=nmodes)
    clim_modal_aero = (nmodes > 0)

    if ( history_budget ) then
       call add_default('EVAPTZM  ', history_budget_histfile_num, ' ')
       call add_default('EVAPQZM  ', history_budget_histfile_num, ' ')
       call add_default('ZMDT     ', history_budget_histfile_num, ' ')
       call add_default('ZMDQ     ', history_budget_histfile_num, ' ')
       call add_default('ZMDLIQ   ', history_budget_histfile_num, ' ')
       call add_default('ZMDICE   ', history_budget_histfile_num, ' ')

       if( cam_physpkg_is('cam4') .or. cam_physpkg_is('cam5') ) then
          call add_default('ZMMTT    ', history_budget_histfile_num, ' ')
       end if

    end if
!
! Limit deep convection to regions below 40 mb
! Note this calculation is repeated in the shallow convection interface
!
    limcnv = 0   ! null value to check against below
    if (pref_edge(1) >= 4.e3_r8) then
       limcnv = 1
    else
       do k=1,plev
          if (pref_edge(k) < 4.e3_r8 .and. pref_edge(k+1) >= 4.e3_r8) then
             limcnv = k
             exit
          end if
       end do
       if ( limcnv == 0 ) limcnv = plevp
    end if
    
    if (masterproc) then
       write(iulog,*)'ZM_CONV_INIT: Deep convection will be capped at intfc ',limcnv, &
            ' which is ',pref_edge(limcnv),' pascals'
    end if
        
    no_deep_pbl = phys_deepconv_pbl()
    call zm_convi(limcnv,no_deep_pbl_in = no_deep_pbl)

    cld_idx         = pbuf_get_index('CLD')
    icwmrdp_idx     = pbuf_get_index('ICWMRDP')
    rprddp_idx      = pbuf_get_index('RPRDDP')
    fracis_idx      = pbuf_get_index('FRACIS')
    nevapr_dpcu_idx = pbuf_get_index('NEVAPR_DPCU')
    prec_dp_idx     = pbuf_get_index('PREC_DP')
    snow_dp_idx     = pbuf_get_index('SNOW_DP')

end subroutine zm_conv_init
!=========================================================================================
!subroutine zm_conv_tend(state, ptend, tdt)

subroutine zm_conv_tend(pblh    ,mcon    ,cme     , &
     tpert   ,dlf     ,pflx    ,zdu      , &
     rliq    , &
     ztodt   , &
     jctop   ,jcbot , &
     state   ,ptend_all   ,landfrac,  pbuf, mu, eu, &
     du, md, ed, dp, dsubcld, jt, maxg, ideep, lengath) 

   use cam_history,   only: outfld
   use physics_types, only: physics_state, physics_ptend
   use physics_types, only: physics_ptend_init, physics_update
   use physics_types, only: physics_state_copy, physics_state_dealloc
   use physics_types, only: physics_ptend_sum, physics_ptend_dealloc

   use phys_grid,     only: get_lat_p, get_lon_p
   use time_manager,  only: get_nstep, is_first_step
   use time_manager,  only: is_first_restart_step             !songxl 2011-09-20
   use physics_buffer, only : pbuf_get_field, physics_buffer_desc, pbuf_old_tim_idx
   use constituents,  only: pcnst, cnst_get_ind, cnst_is_convtran1
   use physconst,     only: gravit
   use phys_control,  only: cam_physpkg_is

   ! Arguments

   type(physics_state), intent(in )   :: state          ! Physics state variables
   type(physics_ptend), intent(out)   :: ptend_all      ! individual parameterization tendencies
   type(physics_buffer_desc), pointer :: pbuf(:)

   real(r8), intent(in) :: ztodt                       ! 2 delta t (model time increment)
   real(r8), intent(in) :: pblh(pcols)                 ! Planetary boundary layer height
   real(r8), intent(in) :: tpert(pcols)                ! Thermal temperature excess
   real(r8), intent(in) :: landfrac(pcols)             ! RBN - Landfrac 

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
   integer, intent(out)  :: lengath

   ! Local variables

   integer :: i,k,m
   integer :: ilon                      ! global longitude index of a column
   integer :: ilat                      ! global latitude index of a column
   integer :: nstep
   integer :: ixcldice, ixcldliq      ! constituent indices for cloud liquid and ice water.
   integer :: lchnk                   ! chunk identifier
   integer :: ncol                    ! number of atmospheric columns
   integer :: itim_old                ! for physics buffer fields

   real(r8) :: ftem(pcols,pver)              ! Temporary workspace for outfld variables
   real(r8) :: ntprprd(pcols,pver)    ! evap outfld: net precip production in layer
   real(r8) :: ntsnprd(pcols,pver)    ! evap outfld: net snow production in layer
   real(r8) :: tend_s_snwprd  (pcols,pver) ! Heating rate of snow production
   real(r8) :: tend_s_snwevmlt(pcols,pver) ! Heating rate of evap/melting of snow
   real(r8) :: fake_dpdry(pcols,pver) ! used in convtran call
   real(r8) :: cnv_nm1_b4_zm_convr(pcols,pver)  !used to store convection parameter before zm_convr
   real(r8) :: hu_nm1_b4_zm_convr(pcols,pver)  !used to store moist static energy before zm_convr
   real(r8) :: ql_b4_zm_convr(pcols,pver)  !used to store moist static energy before zm_convr

   ! physics types
   type(physics_state) :: state1        ! locally modify for evaporation to use, not returned
   type(physics_ptend) :: ptend_loc     ! package tendencies

   ! physics buffer fields
   real(r8), pointer, dimension(:)   :: prec         ! total precipitation
   real(r8), pointer, dimension(:)   :: snow         ! snow from ZM convection 
   real(r8), pointer, dimension(:,:) :: cld
   real(r8), pointer, dimension(:,:) :: ql           ! wg grid slice of cloud liquid water.
   real(r8), pointer, dimension(:,:) :: rprd         ! rain production rate
   real(r8), pointer, dimension(:,:,:) :: fracis  ! fraction of transported species that are insoluble
   real(r8), pointer, dimension(:,:) :: evapcdp      ! Evaporation of deep convective precipitation
   real(r8), pointer, dimension(:,:) :: flxprec      ! Convective-scale flux of precip at interfaces (kg/m2/s)
   real(r8), pointer, dimension(:,:) :: flxsnow      ! Convective-scale flux of snow   at interfaces (kg/m2/s)
   real(r8), pointer, dimension(:,:) :: dp_cldliq
   real(r8), pointer, dimension(:,:) :: dp_cldice
!<songxl 2014-05-20----------
   real(r8), pointer, dimension(:,:) :: hu_nm1
   real(r8), pointer, dimension(:,:) :: cnv_nm1
   real(r8), pointer, dimension(:,:) :: tm1   ! intermediate T between n and n-1 time step
   real(r8), pointer, dimension(:,:) :: qm1   ! intermediate q between n and n-1 time step
!>songxl 2014-05-20---------

   real(r8) :: jctop(pcols)  ! o row of top-of-deep-convection indices passed out.
   real(r8) :: jcbot(pcols)  ! o row of base of cloud indices passed out.
   ! Added freqzmsup as a local variable
   real(r8) :: pcont(pcols), pconb(pcols), freqzm(pcols), freqzmsup(pcols) 

   ! history output fields
   real(r8) :: cape(pcols)        ! w  convective available potential energy.
   real(r8) :: mu_out(pcols,pver)
   real(r8) :: md_out(pcols,pver)

   ! used in momentum transport calculation
   real(r8) :: winds(pcols, pver, 2)
   real(r8) :: wind_tends(pcols, pver, 2)
   real(r8) :: pguall(pcols, pver, 2)
   real(r8) :: pgdall(pcols, pver, 2)
   real(r8) :: icwu(pcols,pver, 2)
   real(r8) :: icwd(pcols,pver, 2)
   real(r8) :: seten(pcols, pver)
   logical  :: l_windt(2)
   real(r8) :: tfinal1, tfinal2
   integer  :: ii

   logical  :: lq(pcnst)

   !----------------------------------------------------------------------

   ! initialize
   lchnk = state%lchnk
   ncol  = state%ncol
   nstep = get_nstep()

   ftem = 0._r8   
   mu_out(:,:) = 0._r8
   md_out(:,:) = 0._r8
   wind_tends(:ncol,:pver,:) = 0.0_r8

   call physics_state_copy(state,state1)             ! copy state to local state1.

   lq(:) = .FALSE.
   lq(1) = .TRUE.
   call physics_ptend_init(ptend_loc, state%psetcols, 'zm_convr', ls=.true., lq=lq)! initialize local ptend type

!
! Associate pointers with physics buffer fields
!
   itim_old = pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, cld_idx,         cld,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

   call pbuf_get_field(pbuf, icwmrdp_idx,     ql )
   call pbuf_get_field(pbuf, rprddp_idx,      rprd )
   call pbuf_get_field(pbuf, fracis_idx,      fracis, start=(/1,1,1/),    kount=(/pcols, pver, pcnst/) )
   call pbuf_get_field(pbuf, nevapr_dpcu_idx, evapcdp )
   call pbuf_get_field(pbuf, prec_dp_idx,     prec )
   call pbuf_get_field(pbuf, snow_dp_idx,     snow )

!<songxl 2014-05-20-----------------
!   if(trigmem)then
     call pbuf_get_field(pbuf, cnv_nm1_idx,     cnv_nm1)
     call pbuf_get_field(pbuf, hu_nm1_idx,      hu_nm1 )
     call pbuf_get_field(pbuf, tm1_idx,         tm1 )
     call pbuf_get_field(pbuf, qm1_idx,         qm1 )
   if(trigmem)then
     if ( is_first_step() .or. is_first_restart_step() ) then
       hu_nm1(:ncol,:pver) = cpair*state%t(:ncol,:pver) + gravit*state%zm(:ncol,:pver)   &
                               + latvap*state%q(:ncol,:pver,1)
       cnv_nm1(:ncol,:pver) = 0._r8
       qm1(:ncol,:pver) =  state%q(:ncol,:pver,1)
       tm1(:ncol,:pver) =  state%t(:ncol,:pver)
     end if
   end if
!<songxl 2014-05-20-----------------

!
! Begin with Zhang-McFarlane (1996) convection parameterization
!
   cnv_nm1_b4_zm_convr(:,:)=cnv_nm1(:,:)
   hu_nm1_b4_zm_convr(:,:)=hu_nm1(:,:)
   ql_b4_zm_convr(:,:)=ql(:,:)
   call t_startf ('zm_convr')
   call zm_convr(   lchnk   ,ncol    , &
                    state%t       ,state%q(:,:,1)     ,prec    ,jctop   ,jcbot   , &
                    pblh    ,state%zm      ,state%phis    ,state%zi      ,ptend_loc%q(:,:,1)    , &
                    ptend_loc%s    ,state%pmid     ,state%pint    ,state%pdel     , &
                    .5_r8*ztodt    ,mcon    ,cme     , cape,      &
                    tpert   ,dlf     ,pflx    ,zdu     ,rprd    , &
                    mu,md,du,eu,ed      , &
                    dp ,dsubcld ,jt,maxg,ideep   , &
                    lengath ,ql      ,rliq  ,landfrac, hu_nm1, cnv_nm1, tm1, qm1 )  !songxl 2014-05-20   
   call t_stopf ('zm_convr')

   call outfld('CAPE', cape, pcols, lchnk)        ! RBN - CAPE output
!
! Output fractional occurance of ZM convection
!
   freqzm(:) = 0._r8
   do i = 1,lengath
      freqzm(ideep(i)) = 1.0_r8
   end do
   call outfld('FREQZM  ',freqzm          ,pcols   ,lchnk   )

  ! ---------------------------------------------------------
  ! Huge chunk moved to after zm_conv_evap
  ! ---------------------------------------------------------

  call physics_ptend_init(ptend_all, state%psetcols, 'zm_conv_tend')

  ! add tendency from this process to tendencies from other processes
  call physics_ptend_sum(ptend_loc,ptend_all, ncol)

  ! update physics state type state1 with ptend_loc 
  call physics_update(state1, ptend_loc, ztodt)

  ! initialize ptend for next process
  lq(:) = .FALSE.
  lq(1) = .TRUE.
  call physics_ptend_init(ptend_loc, state1%psetcols, 'zm_conv_evap', ls=.true., lq=lq)

!
! Determine the phase of the precipitation produced and add latent heat of fusion
! Evaporate some of the precip directly into the environment (Sundqvist)
! Allow this to use the updated state1 and the fresh ptend_loc type
! heating and specific humidity tendencies produced
!

    call pbuf_get_field(pbuf, dp_flxprc_idx, flxprec    )
    call pbuf_get_field(pbuf, dp_flxsnw_idx, flxsnow    )
    call pbuf_get_field(pbuf, dp_cldliq_idx, dp_cldliq  )
    call pbuf_get_field(pbuf, dp_cldice_idx, dp_cldice  )
    dp_cldliq(:ncol,:) = 0._r8
    dp_cldice(:ncol,:) = 0._r8

    call t_startf ('zm_conv_evap')
    call zm_conv_evap(state1%ncol,state1%lchnk, &
         state1%t,state1%pmid,state1%pdel,state1%q(:pcols,:pver,1), &
         ptend_loc%s, tend_s_snwprd, tend_s_snwevmlt, ptend_loc%q(:pcols,:pver,1), &
         rprd, cld, ztodt, &
         prec, snow, ntprprd, ntsnprd , flxprec, flxsnow)
    call t_stopf ('zm_conv_evap')
    evapcdp(:ncol,:pver) = ptend_loc%q(:ncol,:pver,1)
    ! Loop that zeroes tendencies if precipitation rate is too low
    freqzmsup(:) = 0._r8
    do i =1, lengath
       ii = ideep(i)
       ! precipitaiton threshold of 5 mm/d
       if (prec(ii).lt.0.0000000578_r8) then
          freqzmsup(ii) = 1.0_r8
          ! Change all tendencies to zero
          ptend_loc%s(ii,:pver)=0._r8
          ptend_loc%q(ii,:pver,1)=0._r8
          tend_s_snwprd(ii,:pver)=0._r8
          tend_s_snwevmlt(ii,:pver)=0._r8
          prec(ii)=0._r8
          snow(ii)=0._r8
          mcon(ii,:)=0._r8
          ntsnprd(ii,:)=0._r8
          ntprprd(ii,:)=0._r8
          flxprec(ii,:)=0._r8
          flxsnow(ii,:)=0._r8
          ptend_all%s(ii,:pver)=0._r8
          ptend_all%q(ii,:pver,1)=0._r8
          rprd(ii,:)=0._r8
          mu(ii,:)=0._r8
          md(ii,:)=0._r8
          eu(ii,:)=0._r8
          du(ii,:)=0._r8
          zdu(ii,:)=0._r8
          ed(ii,:)=0._r8
          evapcdp(ii,:pver)=0._r8
          pflx(ii,:)=0._r8
          dlf(ii,:)=0._r8
          rliq(ii)=0._r8
          cme(ii,:)=0._r8
          maxg(ii)=1
          jt(ii)=pver
          dsubcld(ii)=0._r8
          hu_nm1(ii,:)=hu_nm1_b4_zm_convr(ii,:)
          cnv_nm1(ii,:)=cnv_nm1_b4_zm_convr(ii,:)
          ql(ii,:)=ql_b4_zm_convr(ii,:)
       endif
    end do
    ! End of loop
    

   
   
   !---- Moved chunk from before zm_conv_evap  --------
   ! Convert mass flux from reported mb/s to kg/m^2/s
   !
   mcon(:ncol,:pver) = mcon(:ncol,:pver) * 100._r8/gravit

   ! Store upward and downward mass fluxes in un-gathered arrays
   ! + convert from mb/s to kg/m^2/s
   do i=1,lengath
      do k=1,pver
         ii = ideep(i)
         mu_out(ii,k) = mu(i,k) * 100._r8/gravit
         md_out(ii,k) = md(i,k) * 100._r8/gravit
      end do
   end do


   if(convproc_do_aer .or. convproc_do_gas) then 
      call outfld('ZMMU', mu_out,      pcols, lchnk)
      call outfld('ZMMD', md_out,      pcols, lchnk)
   else
      call outfld('ZMMU', mu_out(1,1), pcols, lchnk)
      call outfld('ZMMD', md_out(1,1), pcols, lchnk)
   endif

   ftem(:ncol,:pver) = ptend_loc%s(:ncol,:pver)/cpair
   call outfld('ZMDT    ',ftem           ,pcols   ,lchnk   )
   call outfld('ZMDQ    ',ptend_loc%q(1,1,1) ,pcols   ,lchnk   )

!    do i = 1,pcols
!    do i = 1,nco
   pcont(:ncol) = state%ps(:ncol)
   pconb(:ncol) = state%ps(:ncol)
   do i = 1,lengath
       if (maxg(i).gt.jt(i)) then
          pcont(ideep(i)) = state%pmid(ideep(i),jt(i))  ! gathered array (or jctop ungathered)
          pconb(ideep(i)) = state%pmid(ideep(i),maxg(i))! gathered array
       endif
       !     write(iulog,*) ' pcont, pconb ', pcont(i), pconb(i), cnt(i), cnb(i)
    end do
    call outfld('PCONVT  ',pcont          ,pcols   ,lchnk   )
    call outfld('PCONVB  ',pconb          ,pcols   ,lchnk   )
    ! ---------- end of chunk that was moved ----------------

!
! Write out variables from zm_conv_evap
!
   ftem(:ncol,:pver) = ptend_loc%s(:ncol,:pver)/cpair
   call outfld('EVAPTZM ',ftem           ,pcols   ,lchnk   )
   ftem(:ncol,:pver) = tend_s_snwprd  (:ncol,:pver)/cpair
   call outfld('FZSNTZM ',ftem           ,pcols   ,lchnk   )
   ftem(:ncol,:pver) = tend_s_snwevmlt(:ncol,:pver)/cpair
   call outfld('EVSNTZM ',ftem           ,pcols   ,lchnk   )
   call outfld('EVAPQZM ',ptend_loc%q(1,1,1) ,pcols   ,lchnk   )
   call outfld('ZMFLXPRC', flxprec, pcols, lchnk)
   call outfld('ZMFLXSNW', flxsnow, pcols, lchnk)
   call outfld('ZMNTPRPD', ntprprd, pcols, lchnk)
   call outfld('ZMNTSNPD', ntsnprd, pcols, lchnk)
   call outfld('ZMEIHEAT', ptend_loc%s, pcols, lchnk)
   call outfld('CMFMCDZM   ',mcon ,  pcols   ,lchnk   )
   call outfld('PRECCDZM   ',prec,  pcols   ,lchnk   )
   call outfld('FREQZMSUP  ',freqzmsup,  pcols, lchnk)
   call outfld('MAXZM_IDX  ', jt,  pcols, lchnk)


   call outfld('PRECZ   ', prec   , pcols, lchnk)

  ! add tendency from this process to tend from other processes here
  call physics_ptend_sum(ptend_loc,ptend_all, ncol)

  ! update physics state type state1 with ptend_loc 
  call physics_update(state1, ptend_loc, ztodt)


  ! Momentum Transport (non-cam3 physics)

  if ( .not. cam_physpkg_is('cam3')) then

     call physics_ptend_init(ptend_loc, state1%psetcols, 'momtran', ls=.true., lu=.true., lv=.true.)

     winds(:ncol,:pver,1) = state1%u(:ncol,:pver)
     winds(:ncol,:pver,2) = state1%v(:ncol,:pver)
   
     l_windt(1) = .true.
     l_windt(2) = .true.

     call t_startf ('momtran')
     call momtran (lchnk, ncol,                                        &
                   l_windt,winds, 2,  mu(1,1), md(1,1),   &
                   du(1,1), eu(1,1), ed(1,1), dp(1,1), dsubcld(1),  &
                   jt(1),maxg(1), ideep(1), 1, lengath,  &
                   nstep,  wind_tends, pguall, pgdall, icwu, icwd, ztodt, seten )  
     call t_stopf ('momtran')

     ptend_loc%u(:ncol,:pver) = wind_tends(:ncol,:pver,1)
     ptend_loc%v(:ncol,:pver) = wind_tends(:ncol,:pver,2)
     ptend_loc%s(:ncol,:pver) = seten(:ncol,:pver)  

     call physics_ptend_sum(ptend_loc,ptend_all, ncol)

     ! update physics state type state1 with ptend_loc 
     call physics_update(state1, ptend_loc, ztodt)

     ftem(:ncol,:pver) = seten(:ncol,:pver)/cpair
     call outfld('ZMMTT', ftem             , pcols, lchnk)
     call outfld('ZMMTU', wind_tends(1,1,1), pcols, lchnk)
     call outfld('ZMMTV', wind_tends(1,1,2), pcols, lchnk)
   
     ! Output apparent force from  pressure gradient
     call outfld('ZMUPGU', pguall(1,1,1), pcols, lchnk)
     call outfld('ZMUPGD', pgdall(1,1,1), pcols, lchnk)
     call outfld('ZMVPGU', pguall(1,1,2), pcols, lchnk)
     call outfld('ZMVPGD', pgdall(1,1,2), pcols, lchnk)

     ! Output in-cloud winds
     call outfld('ZMICUU', icwu(1,1,1), pcols, lchnk)
     call outfld('ZMICUD', icwd(1,1,1), pcols, lchnk)
     call outfld('ZMICVU', icwu(1,1,2), pcols, lchnk)
     call outfld('ZMICVD', icwd(1,1,2), pcols, lchnk)

   end if

   ! Transport cloud water and ice only
   call cnst_get_ind('CLDLIQ', ixcldliq)
   call cnst_get_ind('CLDICE', ixcldice)

   lq(:)  = .FALSE.
   lq(2:) = cnst_is_convtran1(2:)
   call physics_ptend_init(ptend_loc, state1%psetcols, 'convtran1', lq=lq)


   ! dpdry is not used in this call to convtran since the cloud liquid and ice mixing
   ! ratios are moist
   fake_dpdry(:,:) = 0._r8

   call t_startf ('convtran1')
   call convtran (lchnk,                                        &
                  ptend_loc%lq,state1%q, pcnst,  mu, md,   &
                  du, eu, ed, dp, dsubcld,  &
                  jt,maxg, ideep, 1, lengath,  &
                  nstep,   fracis,  ptend_loc%q, fake_dpdry)
   call t_stopf ('convtran1')

   call outfld('ZMDICE ',ptend_loc%q(1,1,ixcldice) ,pcols   ,lchnk   )
   call outfld('ZMDLIQ ',ptend_loc%q(1,1,ixcldliq) ,pcols   ,lchnk   )

   ! add tendency from this process to tend from other processes here
   call physics_ptend_sum(ptend_loc,ptend_all, ncol)

   call physics_state_dealloc(state1)
   call physics_ptend_dealloc(ptend_loc)

end subroutine zm_conv_tend
!=========================================================================================


subroutine zm_conv_tend_2( state,  ptend,  ztodt, pbuf,mu, eu, &
     du, md, ed, dp, dsubcld, jt, maxg, ideep, lengath, species_class) 

   use physics_types, only: physics_state, physics_ptend, physics_ptend_init
   use time_manager,  only: get_nstep
   use physics_buffer, only: pbuf_get_index, pbuf_get_field, physics_buffer_desc
   use constituents,  only: pcnst, cnst_get_ind, cnst_is_convtran1
   use error_messages, only: alloc_err
   use physconst,      only: spec_class_aerosol, spec_class_gas 
 
! Arguments
   type(physics_state), intent(in )   :: state          ! Physics state variables
   type(physics_ptend), intent(out)   :: ptend          ! indivdual parameterization tendencies
   
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
   integer, intent(in)  :: lengath

   integer, intent(in) :: species_class(:)

! Local variables
   integer :: i, lchnk, istat, m 
   integer :: nstep
   real(r8), dimension(pcols,pver) :: dpdry

! physics buffer fields 
   integer ifld
   real(r8), pointer, dimension(:,:,:) :: fracis  ! fraction of transported species that are insoluble
   logical   :: lq(pcnst)

!
! Initialize
!
  lq(:) = .FALSE.
  lq(:) = .not. cnst_is_convtran1(:)
  call physics_ptend_init(ptend, state%psetcols, 'convtran2', lq=lq )

!
! Associate pointers with physics buffer fields
!
   ifld = pbuf_get_index('FRACIS')
   call pbuf_get_field(pbuf, fracis_idx, fracis, start=(/1,1,1/), kount=(/pcols, pver, pcnst/) )

!
! Transport all constituents except cloud water and ice
!

  lchnk = state%lchnk

   nstep = get_nstep()
   if((convproc_do_aer .or. convproc_do_gas) .and. clim_modal_aero) then
      do m = 1, pcnst
         if ( (species_class(m) == spec_class_aerosol .and. convproc_do_aer) .or. &
              (species_class(m) == spec_class_gas     .and. convproc_do_gas) ) then
            ptend%lq(m) = .false.
         end if
      enddo
   endif


   if (any(ptend%lq(:))) then
      ! initialize dpdry for call to convtran
      ! it is used for tracers of dry mixing ratio type
      dpdry = 0._r8
      do i = 1,lengath
         dpdry(i,:) = state%pdeldry(ideep(i),:)/100._r8
      end do

      call t_startf ('convtran2')
      call convtran (lchnk,                                        &
                     ptend%lq,state%q, pcnst,  mu, md,   &
                     du, eu, ed, dp, dsubcld,  &
                     jt,maxg,ideep, 1, lengath,  &
                     nstep,   fracis,  ptend%q, dpdry)
      call t_stopf ('convtran2')
   end if

   if((convproc_do_aer .or. convproc_do_gas) .and. clim_modal_aero) then
      do m = 1, pcnst
         if ( (species_class(m) == spec_class_aerosol .and. convproc_do_aer) .or. &
              (species_class(m) == spec_class_gas     .and. convproc_do_gas) ) then
            ptend%lq(m) = .false.
            ptend%q(:,:,m) = 0.0_r8
         end if
      enddo
   endif

end subroutine zm_conv_tend_2

!=========================================================================================



end module zm_conv_intr
