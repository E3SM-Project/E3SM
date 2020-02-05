
module cloud_diagnostics

!---------------------------------------------------------------------------------
! Purpose:
!
! Put cloud physical specifications on the history tape
!  Modified from code that computed cloud optics
!
! Author: Byron Boville  Sept 06, 2002
!  Modified Oct 15, 2008
!    
!
!---------------------------------------------------------------------------------

   use shr_kind_mod,  only: r8=>shr_kind_r8
   use ppgrid,        only: pcols, pver,pverp
   use physconst,     only: gravit
   use cam_history,   only: outfld
   use cam_history,   only: addfld, horiz_only, add_default

   implicit none
   private
   save

   public :: cloud_diagnostics_init
   public :: cloud_diagnostics_calc
   public :: cloud_diagnostics_register

! Local variables
   integer :: dei_idx, mu_idx, lambda_idx, iciwp_idx, iclwp_idx, cld_idx  ! index into pbuf for cloud fields
   integer :: ixcldice, ixcldliq, rei_idx, rel_idx

   logical :: do_cld_diag, mg_clouds, rk_clouds, camrt_rad
   integer :: conv_water_in_rad
   
   integer :: cicewp_idx = -1
   integer :: cliqwp_idx = -1
   integer :: cldemis_idx = -1
   integer :: cldtau_idx = -1
   integer :: nmxrgn_idx = -1
   integer :: pmxrgn_idx = -1

   ! Index fields for precipitation efficiency.
   integer :: acpr_idx, acgcme_idx, acnum_idx


contains

!===============================================================================
  subroutine cloud_diagnostics_register

    use phys_control,  only: phys_getopts
    use physics_buffer,only: pbuf_add_field, dtype_r8, dtype_i4

    character(len=16) :: rad_pkg, microp_pgk

    call phys_getopts(radiation_scheme_out=rad_pkg,microp_scheme_out=microp_pgk)
    camrt_rad = rad_pkg .eq. 'camrt'
    rk_clouds = microp_pgk == 'RK'
    mg_clouds = microp_pgk == 'MG'

    if (rk_clouds) then
       call pbuf_add_field('CLDEMIS','physpkg', dtype_r8,(/pcols,pver/), cldemis_idx)
       call pbuf_add_field('CLDTAU', 'physpkg', dtype_r8,(/pcols,pver/), cldtau_idx)

       call pbuf_add_field('CICEWP', 'physpkg', dtype_r8,(/pcols,pver/), cicewp_idx)
       call pbuf_add_field('CLIQWP', 'physpkg', dtype_r8,(/pcols,pver/), cliqwp_idx)

       call pbuf_add_field('PMXRGN', 'physpkg', dtype_r8,(/pcols,pverp/), pmxrgn_idx)
       call pbuf_add_field('NMXRGN', 'physpkg', dtype_i4,(/pcols /),      nmxrgn_idx)
    else if (mg_clouds) then
       ! In cloud ice water path for radiation
       call pbuf_add_field('ICIWP',      'global', dtype_r8,(/pcols,pver/), iciwp_idx)
       ! In cloud liquid water path for radiation
       call pbuf_add_field('ICLWP',      'global', dtype_r8,(/pcols,pver/), iclwp_idx)
    endif
  end subroutine cloud_diagnostics_register

!===============================================================================
  subroutine cloud_diagnostics_init()
!-----------------------------------------------------------------------
    use physics_buffer,only: pbuf_get_index
    use phys_control,  only: phys_getopts
    use constituents,  only: cnst_get_ind
    use cloud_cover_diags, only: cloud_cover_diags_init

    implicit none

!-----------------------------------------------------------------------

    character(len=16) :: wpunits, sampling_seq
    logical           :: history_amwg          ! output the variables used by the AMWG diag package
    logical           :: history_verbose       ! produce verbose history output

    logical :: use_MMF
    call phys_getopts(use_MMF_out=use_MMF)

    !-----------------------------------------------------------------------

    cld_idx    = pbuf_get_index('CLD')

    ! ----------------------------
    ! determine default variables
    ! ----------------------------
    call phys_getopts( history_amwg_out = history_amwg, &
                       history_verbose_out = history_verbose )

    if (mg_clouds) then

       call addfld ('ICWMR', (/ 'lev' /), 'A', 'kg/kg', 'Prognostic in-cloud water mixing ratio'                  )
       call addfld ('ICIMR', (/ 'lev' /), 'A', 'kg/kg', 'Prognostic in-cloud ice mixing ratio'                    )
       call addfld ('IWC', (/ 'lev' /), 'A', 'kg/m3', 'Grid box average ice water content'                      )
       call addfld ('LWC', (/ 'lev' /), 'A', 'kg/m3', 'Grid box average liquid water content'                   )


       if (history_amwg) then
          call add_default ('ICWMR', 1, ' ')
          call add_default ('ICIMR', 1, ' ')
          call add_default ('IWC      ', 1, ' ')
       end if

       dei_idx    = pbuf_get_index('DEI')
       mu_idx     = pbuf_get_index('MU')
       lambda_idx = pbuf_get_index('LAMBDAC')

    elseif (rk_clouds) then

       rei_idx    = pbuf_get_index('REI')
       rel_idx    = pbuf_get_index('REL')

    endif

    call cnst_get_ind('CLDICE', ixcldice)
    call cnst_get_ind('CLDLIQ', ixcldliq)

    do_cld_diag = rk_clouds .or. mg_clouds

    if (.not.do_cld_diag) return
    
    call phys_getopts(conv_water_in_rad_out=conv_water_in_rad)

    if (use_MMF) then
          wpunits = 'kg/m2'
          sampling_seq=''
    else 
      if (rk_clouds) then 
         wpunits = 'gram/m2'
         sampling_seq='rad_lwsw'
      else if (mg_clouds) then 
         wpunits = 'kg/m2'
         sampling_seq=''
      endif
    endif
    
    call addfld ('ICLDIWP', (/ 'lev' /), 'A', wpunits,'In-cloud ice water path'               , sampling_seq=sampling_seq)
    call addfld ('ICLDTWP', (/ 'lev' /), 'A',wpunits,'In-cloud cloud total water path (liquid and ice)', &
         sampling_seq=sampling_seq)

    call addfld ('GCLDLWP',(/ 'lev' /), 'A',wpunits,'Grid-box cloud water path'             , &
         sampling_seq=sampling_seq)
    call addfld ('TGCLDCWP',horiz_only,    'A',wpunits,'Total grid-box cloud water path (liquid and ice)', &
         sampling_seq=sampling_seq)
    call addfld ('TGCLDLWP',horiz_only,    'A',wpunits,'Total grid-box cloud liquid water path', &
         sampling_seq=sampling_seq)
    call addfld ('TGCLDIWP',horiz_only,    'A',wpunits,'Total grid-box cloud ice water path'   , &
         sampling_seq=sampling_seq)
    
    if(mg_clouds) then
       call addfld ('lambda_cloud',(/ 'lev' /),'I','1/meter','lambda in cloud')
       call addfld ('mu_cloud',(/ 'lev' /),'I','1','mu in cloud')
       call addfld ('dei_cloud',(/ 'lev' /),'I','micrometers','ice radiative effective diameter in cloud')
    endif

    if(rk_clouds) then
       call addfld ('rel_cloud',(/ 'lev' /),'I','1/meter','effective radius of liq in cloud', sampling_seq=sampling_seq)
       call addfld ('rei_cloud',(/ 'lev' /),'I','1','effective radius of ice in cloud', sampling_seq=sampling_seq)
    endif

    call addfld ('SETLWP',(/ 'lev' /), 'A','gram/m2','Prescribed liquid water path'          , sampling_seq=sampling_seq)
    call addfld ('LWSH',horiz_only,    'A','m','Liquid water scale height'             , sampling_seq=sampling_seq)

    call addfld ('EFFCLD',(/ 'lev' /), 'A','fraction','Effective cloud fraction'              , sampling_seq=sampling_seq)

    if (camrt_rad) then
       call addfld ('EMIS', (/ 'lev' /), 'A', '1','cloud emissivity'                      , sampling_seq=sampling_seq)
    else
       call addfld ('EMISCLD', (/ 'lev' /), 'A', '1','cloud emissivity'                      , sampling_seq=sampling_seq)
    endif

    call cloud_cover_diags_init(sampling_seq)


    if (history_amwg) then
       call add_default ('TGCLDLWP', 1, ' ')
       call add_default ('TGCLDIWP', 1, ' ')
       call add_default ('TGCLDCWP', 1, ' ')
       if (camrt_rad) then
           call add_default ('EMIS', 1, ' ')
       else
           if (history_verbose) call add_default ('EMISCLD', 1, ' ')
       endif
    endif

    return
  end subroutine cloud_diagnostics_init

subroutine cloud_diagnostics_calc(state,  pbuf)
!===============================================================================
!
! Compute (liquid+ice) water path and cloud water/ice diagnostics
! *** soon this code will compute liquid and ice paths from input liquid and ice mixing ratios
! 
! **** mixes interface and physics code temporarily
!-----------------------------------------------------------------------
    use physics_types, only: physics_state    
    use physics_buffer,only: physics_buffer_desc, pbuf_get_field, pbuf_old_tim_idx
    use phys_control,  only: phys_getopts
    use pkg_cldoptics, only: cldovrlap, cldclw,  cldems
    use conv_water,    only: conv_water_4rad
    use radiation,     only: radiation_do
    use cloud_cover_diags, only: cloud_cover_diags_out

    use ref_pres,       only: top_lev=>trop_cloud_top_lev

    implicit none

! Arguments
    type(physics_state), intent(in)    :: state        ! state variables
    type(physics_buffer_desc), pointer :: pbuf(:)

! Local variables

    real(r8), pointer :: cld(:,:)       ! cloud fraction
    real(r8), pointer :: iciwp(:,:)   ! in-cloud cloud ice water path
    real(r8), pointer :: iclwp(:,:)   ! in-cloud cloud liquid water path
    real(r8), pointer :: dei(:,:)       ! effective radiative diameter of ice
    real(r8), pointer :: mu(:,:)        ! gamma distribution for liq clouds
    real(r8), pointer :: lambda(:,:)    ! gamma distribution for liq clouds
    real(r8), pointer :: rei(:,:)       ! effective radiative radius of ice
    real(r8), pointer :: rel(:,:)       ! effective radiative radius of liq

    real(r8), pointer :: cldemis(:,:)   ! cloud emissivity
    real(r8), pointer :: cldtau(:,:)    ! cloud optical depth
    real(r8), pointer :: cicewp(:,:)    ! in-cloud cloud ice water path
    real(r8), pointer :: cliqwp(:,:)    ! in-cloud cloud liquid water path

    integer,  pointer :: nmxrgn(:)      ! Number of maximally overlapped regions
    real(r8), pointer :: pmxrgn(:,:)    ! Maximum values of pressure for each

    integer :: itim_old

    real(r8) :: cwp   (pcols,pver)      ! in-cloud cloud (total) water path
    real(r8) :: gicewp(pcols,pver)      ! grid-box cloud ice water path
    real(r8) :: gliqwp(pcols,pver)      ! grid-box cloud liquid water path
    real(r8) :: gwp   (pcols,pver)      ! grid-box cloud (total) water path
    real(r8) :: tgicewp(pcols)          ! Vertically integrated ice water path
    real(r8) :: tgliqwp(pcols)          ! Vertically integrated liquid water path
    real(r8) :: tgwp   (pcols)          ! Vertically integrated (total) cloud water path

    real(r8) :: ficemr (pcols,pver)     ! Ice fraction from ice and liquid mixing ratios

    real(r8) :: icimr(pcols,pver)       ! In cloud ice mixing ratio
    real(r8) :: icwmr(pcols,pver)       ! In cloud water mixing ratio
    real(r8) :: iwc(pcols,pver)         ! Grid box average ice water content
    real(r8) :: lwc(pcols,pver)         ! Grid box average liquid water content

! old data
    real(r8) :: tpw    (pcols)          ! total precipitable water
    real(r8) :: clwpold(pcols,pver)     ! Presribed cloud liq. h2o path
    real(r8) :: hl     (pcols)          ! Liquid water scale height

    integer :: i,k                      ! loop indexes
    integer :: ncol, lchnk
    real(r8) :: rgrav

    real(r8) :: allcld_ice (pcols,pver) ! Convective cloud ice
    real(r8) :: allcld_liq (pcols,pver) ! Convective cloud liquid

    real(r8) :: effcld(pcols,pver)      ! effective cloud=cld*emis

    logical :: dosw,dolw
    logical :: use_MMF
  
!-----------------------------------------------------------------------
    if (.not.do_cld_diag) return

    call phys_getopts(use_MMF_out=use_MMF)

    if(rk_clouds) then
       dosw     = radiation_do('sw')      ! do shortwave heating calc this timestep?
       dolw     = radiation_do('lw')      ! do longwave heating calc this timestep?
    else
       dosw     = .true.
       dolw     = .true.
    endif

    if (.not.(dosw .or. dolw)) return

    ncol  = state%ncol
    lchnk = state%lchnk

    itim_old = pbuf_old_tim_idx()
    call pbuf_get_field(pbuf, cld_idx, cld, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

    if(mg_clouds)then

       call pbuf_get_field(pbuf, iclwp_idx, iclwp )
       call pbuf_get_field(pbuf, iciwp_idx, iciwp )
       call pbuf_get_field(pbuf, dei_idx, dei )
       call pbuf_get_field(pbuf, mu_idx, mu )
       call pbuf_get_field(pbuf, lambda_idx, lambda )

       call outfld('dei_cloud',dei(:,:),pcols,lchnk)
       call outfld('mu_cloud',mu(:,:),pcols,lchnk)
       call outfld('lambda_cloud',lambda(:,:),pcols,lchnk)

    elseif(rk_clouds) then

       call pbuf_get_field(pbuf, rei_idx, rei )
       call pbuf_get_field(pbuf, rel_idx, rel )

       call outfld('rel_cloud', rel, pcols, lchnk)
       call outfld('rei_cloud', rei, pcols, lchnk)

       if (cldemis_idx>0) then
          call pbuf_get_field(pbuf, cldemis_idx, cldemis )
       else
          allocate(cldemis(pcols,pver))
       endif
       if (cldtau_idx>0) then
          call pbuf_get_field(pbuf, cldtau_idx, cldtau )
       else
          allocate(cldtau(pcols,pver))
       endif

    endif

    if (cicewp_idx>0) then
       call pbuf_get_field(pbuf, cicewp_idx, cicewp )
    else
       allocate(cicewp(pcols,pver))
    endif
    if (cliqwp_idx>0) then
       call pbuf_get_field(pbuf, cliqwp_idx, cliqwp )
    else
       allocate(cliqwp(pcols,pver))
    endif

    if (nmxrgn_idx>0) then
       call pbuf_get_field(pbuf, nmxrgn_idx, nmxrgn )
    else
       allocate(nmxrgn(pcols))
    endif

    if (pmxrgn_idx>0) then
       call pbuf_get_field(pbuf, pmxrgn_idx, pmxrgn )
    else
       allocate(pmxrgn(pcols,pverp))
    endif

! Compute liquid and ice water paths
    if(mg_clouds) then

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
          allcld_ice(:ncol,top_lev:pver) = state%q(:ncol,top_lev:pver,ixcldice)  !           "        ice
       end if

       ! ------------------------------------------------------------ !
       ! Compute in cloud ice and liquid mixing ratios                !
       ! Note that 'iclwp, iciwp' are used for radiation computation. !
       ! ------------------------------------------------------------ !


       iciwp = 0._r8
       iclwp = 0._r8
       icimr = 0._r8
       icwmr = 0._r8
       iwc = 0._r8
       lwc = 0._r8

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

       do k=1,pver
          do i = 1,ncol
             gicewp(i,k) = iciwp(i,k)*cld(i,k)
             gliqwp(i,k) = iclwp(i,k)*cld(i,k)
             cicewp(i,k) = iciwp(i,k)
             cliqwp(i,k) = iclwp(i,k)
          end do
       end do

    elseif(rk_clouds) then

       if (conv_water_in_rad /= 0) then
          call conv_water_4rad(state,pbuf,conv_water_in_rad,allcld_liq,allcld_ice)
       else
          allcld_liq = state%q(:,:,ixcldliq)
          allcld_ice = state%q(:,:,ixcldice)
       end if
    
       do k=1,pver
          do i = 1,ncol
             gicewp(i,k) = allcld_ice(i,k)*state%pdel(i,k)/gravit*1000.0_r8  ! Grid box ice water path.
             gliqwp(i,k) = allcld_liq(i,k)*state%pdel(i,k)/gravit*1000.0_r8  ! Grid box liquid water path.
             cicewp(i,k) = gicewp(i,k) / max(0.01_r8,cld(i,k))               ! In-cloud ice water path.
             cliqwp(i,k) = gliqwp(i,k) / max(0.01_r8,cld(i,k))               ! In-cloud liquid water path.
             ficemr(i,k) = allcld_ice(i,k) / max(1.e-10_r8,(allcld_ice(i,k) + allcld_liq(i,k)))
          end do
       end do
    endif

! Determine parameters for maximum/random overlap
    call cldovrlap(lchnk, ncol, state%pint, cld, nmxrgn, pmxrgn)

! Cloud cover diagnostics (done in radiation_tend for camrt)
    if (.not.camrt_rad) then
       call cloud_cover_diags_out(lchnk, ncol, cld, state%pmid, nmxrgn, pmxrgn )
    endif
    
    tgicewp(:ncol) = 0._r8
    tgliqwp(:ncol) = 0._r8

    do k=1,pver
       tgicewp(:ncol)  = tgicewp(:ncol) + gicewp(:ncol,k)
       tgliqwp(:ncol)  = tgliqwp(:ncol) + gliqwp(:ncol,k)
    end do

    tgwp(:ncol) = tgicewp(:ncol) + tgliqwp(:ncol)
    gwp(:ncol,:pver) = gicewp(:ncol,:pver) + gliqwp(:ncol,:pver)
    cwp(:ncol,:pver) = cicewp(:ncol,:pver) + cliqwp(:ncol,:pver)

    if(rk_clouds) then

       ! Cloud emissivity.
       call cldems(lchnk, ncol, cwp, ficemr, rei, cldemis, cldtau)
       
       ! Effective cloud cover
       do k=1,pver
          do i=1,ncol
             effcld(i,k) = cld(i,k)*cldemis(i,k)
          end do
       end do
       
       call outfld('EFFCLD'  ,effcld , pcols,lchnk)
       if (camrt_rad) then
          call outfld('EMIS' ,cldemis, pcols,lchnk)
       else
          call outfld('EMISCLD' ,cldemis, pcols,lchnk)
       endif

    else if (mg_clouds) then

       ! --------------------------------------------- !
       ! General outfield calls for microphysics       !
       ! --------------------------------------------- !

       call outfld( 'IWC'      , iwc,         pcols, lchnk )
       call outfld( 'LWC'      , lwc,         pcols, lchnk )
       call outfld( 'ICIMR'    , icimr,       pcols, lchnk )
       call outfld( 'ICWMR'    , icwmr,       pcols, lchnk )

    endif

    if (.not. use_MMF) then 
       ! for MMF, these are diagnostics in crm_physics.F90
       call outfld('GCLDLWP' ,gwp    , pcols,lchnk)
       call outfld('TGCLDCWP',tgwp   , pcols,lchnk)
       call outfld('TGCLDLWP',tgliqwp, pcols,lchnk)
       call outfld('TGCLDIWP',tgicewp, pcols,lchnk)
       call outfld('ICLDTWP' ,cwp    , pcols,lchnk)
       call outfld('ICLDIWP' ,cicewp , pcols,lchnk)
    endif

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
    
    if(rk_clouds) then
       if (cldemis_idx<0) deallocate(cldemis)
       if (cldtau_idx<0) deallocate(cldtau)
    endif
    if (cicewp_idx<0) deallocate(cicewp)
    if (cliqwp_idx<0) deallocate(cliqwp)
    if (pmxrgn_idx<0) deallocate(pmxrgn)
    if (nmxrgn_idx<0) deallocate(nmxrgn)

    return
end subroutine cloud_diagnostics_calc

end module cloud_diagnostics
