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
   use cam_history,   only: addfld, add_default, phys_decomp

   implicit none
   private
   save

   public :: cloud_diagnostics_init
   public :: cloud_diagnostics_calc
   public :: cloud_diagnostics_register

! Local variables
   integer :: i_dei, i_mu, i_lambda, i_iciwp, i_iclwp, i_cld  ! index into pbuf for cloud fields
   integer :: i_ciwp, i_clwp, i_rei, i_rel

   logical :: do_cld_diag, mg_clouds, rk_clouds, camrt_rad
   integer :: conv_water_in_rad
   
   integer :: cicewp_idx = -1
   integer :: cliqwp_idx = -1
   integer :: cldemis_idx = -1
   integer :: cldtau_idx = -1
   integer :: nmxrgn_idx = -1
   integer :: pmxrgn_idx = -1

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
    logical           :: history_amwg                  ! output the variables used by the AMWG diag package

    i_cld    = pbuf_get_index('CLD')

    if (mg_clouds) then

       i_iciwp  = pbuf_get_index('ICIWP')
       i_iclwp  = pbuf_get_index('ICLWP')

       i_dei    = pbuf_get_index('DEI')
       i_mu     = pbuf_get_index('MU')
       i_lambda = pbuf_get_index('LAMBDAC')

    elseif (rk_clouds) then

       i_rei    = pbuf_get_index('REI')
       i_rel    = pbuf_get_index('REL')

       call cnst_get_ind('CLDICE', i_ciwp)
       call cnst_get_ind('CLDLIQ', i_clwp)

    endif

    do_cld_diag = rk_clouds .or. mg_clouds

    if (.not.do_cld_diag) return
    
    if (rk_clouds) then 
       wpunits = 'gram/m2'
       sampling_seq='rad_lwsw'
    else if (mg_clouds) then 
       wpunits = 'kg/m2'
       sampling_seq=''
    endif

    call addfld ('ICLDIWP', wpunits, pver, 'A','In-cloud ice water path'               ,phys_decomp, sampling_seq=sampling_seq)
    call addfld ('ICLDTWP ',wpunits, pver, 'A','In-cloud cloud total water path (liquid and ice)',phys_decomp, &
         sampling_seq=sampling_seq)

    call addfld ('GCLDLWP ',wpunits,pver, 'A','Grid-box cloud water path'             ,phys_decomp, &
         sampling_seq=sampling_seq)
    call addfld ('TGCLDCWP',wpunits,1,    'A','Total grid-box cloud water path (liquid and ice)',phys_decomp, &
         sampling_seq=sampling_seq)
    call addfld ('TGCLDLWP',wpunits,1,    'A','Total grid-box cloud liquid water path',phys_decomp, &
         sampling_seq=sampling_seq)
    call addfld ('TGCLDIWP',wpunits,1,    'A','Total grid-box cloud ice water path'   ,phys_decomp, &
         sampling_seq=sampling_seq)
    
    if(mg_clouds) then
       call addfld ('lambda_cloud','1/meter',pver,'I','lambda in cloud', phys_decomp)
       call addfld ('mu_cloud','1',pver,'I','mu in cloud', phys_decomp)
       call addfld ('dei_cloud','micrometers',pver,'I','ice radiative effective diameter in cloud', phys_decomp)
    endif

    if(rk_clouds) then
       call addfld ('rel_cloud','1/meter',pver,'I','effective radius of liq in cloud', phys_decomp, sampling_seq=sampling_seq)
       call addfld ('rei_cloud','1',pver,'I','effective radius of ice in cloud', phys_decomp, sampling_seq=sampling_seq)
       call phys_getopts(conv_water_in_rad_out=conv_water_in_rad)
    endif

    call addfld ('SETLWP  ','gram/m2 ',pver, 'A','Prescribed liquid water path'          ,phys_decomp, sampling_seq=sampling_seq)
    call addfld ('LWSH    ','m       ',1,    'A','Liquid water scale height'             ,phys_decomp, sampling_seq=sampling_seq)

    call addfld ('EFFCLD  ','fraction',pver, 'A','Effective cloud fraction'              ,phys_decomp, sampling_seq=sampling_seq)

    if (camrt_rad) then
       call addfld ('EMIS', '1', pver, 'A','cloud emissivity'                      ,phys_decomp, sampling_seq=sampling_seq)
    else
       call addfld ('EMISCLD', '1', pver, 'A','cloud emissivity'                      ,phys_decomp, sampling_seq=sampling_seq)
    endif

    call cloud_cover_diags_init(sampling_seq)

    ! ----------------------------
    ! determine default variables
    ! ----------------------------
    call phys_getopts( history_amwg_out = history_amwg)

    if (history_amwg) then
       call add_default ('TGCLDLWP', 1, ' ')
       call add_default ('TGCLDIWP', 1, ' ')
       call add_default ('TGCLDCWP', 1, ' ')
       if (camrt_rad) then
           call add_default ('EMIS', 1, ' ')
       else
           call add_default ('EMISCLD', 1, ' ')
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
    use pkg_cldoptics, only: cldovrlap, cldclw,  cldems
    use conv_water,    only: conv_water_4rad
    use radiation,     only: radiation_do
    use cloud_cover_diags, only: cloud_cover_diags_out

    implicit none

! Arguments
    type(physics_state), intent(in)    :: state        ! state variables
    type(physics_buffer_desc), pointer :: pbuf(:)

! Local variables

    real(r8), pointer :: cld(:,:)       ! cloud fraction
    real(r8), pointer :: iciwpth(:,:)   ! in-cloud cloud ice water path
    real(r8), pointer :: iclwpth(:,:)   ! in-cloud cloud liquid water path
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

    integer :: itim

    real(r8) :: cwp   (pcols,pver)      ! in-cloud cloud (total) water path
    real(r8) :: gicewp(pcols,pver)      ! grid-box cloud ice water path
    real(r8) :: gliqwp(pcols,pver)      ! grid-box cloud liquid water path
    real(r8) :: gwp   (pcols,pver)      ! grid-box cloud (total) water path
    real(r8) :: tgicewp(pcols)          ! Vertically integrated ice water path
    real(r8) :: tgliqwp(pcols)          ! Vertically integrated liquid water path
    real(r8) :: tgwp   (pcols)          ! Vertically integrated (total) cloud water path

    real(r8) :: ficemr (pcols,pver)     ! Ice fraction from ice and liquid mixing ratios

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
  
!-----------------------------------------------------------------------
    if (.not.do_cld_diag) return

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

    itim = pbuf_old_tim_idx()
    call pbuf_get_field(pbuf, i_cld, cld, start=(/1,1,itim/), kount=(/pcols,pver,1/) )

    if(mg_clouds)then

       call pbuf_get_field(pbuf, i_iclwp, iclwpth )
       call pbuf_get_field(pbuf, i_iciwp, iciwpth )
       call pbuf_get_field(pbuf, i_dei, dei )
       call pbuf_get_field(pbuf, i_mu, mu )
       call pbuf_get_field(pbuf, i_lambda, lambda )

       call outfld('dei_cloud',dei(:,:),pcols,lchnk)
       call outfld('mu_cloud',mu(:,:),pcols,lchnk)
       call outfld('lambda_cloud',lambda(:,:),pcols,lchnk)

    elseif(rk_clouds) then

       call pbuf_get_field(pbuf, i_rei, rei )
       call pbuf_get_field(pbuf, i_rel, rel )

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

       do k=1,pver
          do i = 1,ncol
             gicewp(i,k) = iciwpth(i,k)*cld(i,k)
             gliqwp(i,k) = iclwpth(i,k)*cld(i,k)
             cicewp(i,k) = iciwpth(i,k)
             cliqwp(i,k) = iclwpth(i,k)
          end do
       end do

    elseif(rk_clouds) then

       if (conv_water_in_rad /= 0) then
          call conv_water_4rad(lchnk,ncol,pbuf,conv_water_in_rad,rei,state%pdel, &
                               state%q(:,:,i_clwp),state%q(:,:,i_ciwp),allcld_liq,allcld_ice)
       else
          allcld_liq = state%q(:,:,i_clwp)
          allcld_ice = state%q(:,:,i_ciwp)
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

    endif

    call outfld('GCLDLWP' ,gwp    , pcols,lchnk)
    call outfld('TGCLDCWP',tgwp   , pcols,lchnk)
    call outfld('TGCLDLWP',tgliqwp, pcols,lchnk)
    call outfld('TGCLDIWP',tgicewp, pcols,lchnk)
    call outfld('ICLDTWP' ,cwp    , pcols,lchnk)
    call outfld('ICLDIWP' ,cicewp , pcols,lchnk)

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
