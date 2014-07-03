!================================================================================================
! output data necessary to drive radiation offline
! Francis Vitt -- Created 15 Dec 2009
!================================================================================================
module radiation_data

  use shr_kind_mod,     only: r8=>shr_kind_r8
  use ppgrid,           only: pcols, pver, pverp
  use cam_history,      only: addfld, add_default, phys_decomp, outfld
  use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_gas, rad_cnst_get_aer_mmr
  use radconstants,     only: nradgas, gaslist
  use cam_history_support, only: fieldname_len, fillvalue
  use spmd_utils,       only: masterproc
  use abortutils,       only: endrun

  implicit none
  private

  public :: output_rad_data
  public :: init_rad_data
  public :: rad_data_readnl

  integer :: cld_ifld,concld_ifld,rel_ifld,rei_ifld
  integer :: dei_ifld,mu_ifld,lambdac_ifld,iciwp_ifld,iclwp_ifld,rel_fn_ifld
  integer :: des_ifld,icswp_ifld,cldfsnow_ifld

  character(len=fieldname_len), public, parameter :: &
       lndfrc_fldn    = 'rad_lndfrc      ' , &
       icefrc_fldn    = 'rad_icefrc      ' , &
       snowh_fldn     = 'rad_snowh       ' , &
       landm_fldn     = 'rad_landm       ' , &
       asdir_fldn     = 'rad_asdir       ' , &
       asdif_fldn     = 'rad_asdif       ' , &
       aldir_fldn     = 'rad_aldir       ' , &
       aldif_fldn     = 'rad_aldif       ' , &
       coszen_fldn    = 'rad_coszen      ' , &
       asdir_pos_fldn = 'rad_asdir_pos   ' , &
       asdif_pos_fldn = 'rad_asdif_pos   ' , &
       aldir_pos_fldn = 'rad_aldir_pos   ' , &
       aldif_pos_fldn = 'rad_aldif_pos   ' , &
       lwup_fldn      = 'rad_lwup        ' , &
       ts_fldn        = 'rad_ts          ' , &
       temp_fldn      = 'rad_temp        ' , &
       pdel_fldn      = 'rad_pdel        ' , &
       pdeldry_fldn   = 'rad_pdeldry     ' , &
       pmid_fldn      = 'rad_pmid        ' , &
       watice_fldn    = 'rad_watice      ' , &
       watliq_fldn    = 'rad_watliq      ' , &
       watvap_fldn    = 'rad_watvap      ' , &
       zint_fldn      = 'rad_zint        ' , &
       pint_fldn      = 'rad_pint        ' , &
       cld_fldn       = 'rad_cld         ' , &
       cldfsnow_fldn  = 'rad_cldfsnow    ' , &
       concld_fldn    = 'rad_concld      ' , &
       rel_fldn       = 'rad_rel         ' , &
       rei_fldn       = 'rad_rei         ' , &
       dei_fldn       = 'rad_dei         ' , &
       des_fldn       = 'rad_des         ' , &
       mu_fldn        = 'rad_mu          ' , &
       lambdac_fldn   = 'rad_lambdac     ' , &
       iciwp_fldn     = 'rad_iciwp       ' , &
       iclwp_fldn     = 'rad_iclwp       ' , &
       icswp_fldn     = 'rad_icswp       '

  ! rad constituents mixing ratios
  integer :: ngas, naer
  character(len=64), allocatable :: gasnames(:)
  character(len=64), allocatable :: aernames(:)
 
  ! control options  
  logical          :: rad_data_output = .false.
  integer          :: rad_data_histfile_num = 2
  character(len=1) :: rad_data_avgflag = 'A'

  ! MG microphys check
  logical, public :: mg_microphys

contains

!================================================================================================
!================================================================================================
  subroutine rad_data_readnl(nlfile)

    ! Read rad_data_nl namelist group.  Parse input.

    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use mpishorthand

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr, i
    character(len=*), parameter :: subname = 'rad_data_readnl'

    namelist /rad_data_nl/ rad_data_output, rad_data_histfile_num, rad_data_avgflag

    !-----------------------------------------------------------------------------

    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'rad_data_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, rad_data_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

#ifdef SPMD
    ! Broadcast namelist variables
    call mpibcast (rad_data_output,       1,   mpilog ,  0, mpicom)
    call mpibcast (rad_data_histfile_num, 1,   mpiint ,  0, mpicom)
    call mpibcast (rad_data_avgflag,      1,   mpichar , 0, mpicom)
#endif
    
  end subroutine rad_data_readnl

  !================================================================================================
  !================================================================================================
  subroutine init_rad_data
    use phys_control,     only: phys_getopts
    use physics_buffer, only: pbuf_get_index
    implicit none
    
    integer :: i
    character(len=64) :: name
    character(len=128):: long_name
    character(len=64) :: long_name_description
    character(len=16)  :: microp_scheme  ! microphysics scheme

    if (.not.rad_data_output) return
   
    call phys_getopts(microp_scheme_out=microp_scheme)
    mg_microphys =  (trim(microp_scheme) == 'MG')

    cld_ifld    = pbuf_get_index('CLD')
    concld_ifld = pbuf_get_index('CONCLD')
    rel_ifld    = pbuf_get_index('REL')
    rei_ifld    = pbuf_get_index('REI')
    if (mg_microphys) then
       dei_ifld      = pbuf_get_index('DEI')
       des_ifld      = pbuf_get_index('DES')
       mu_ifld       = pbuf_get_index('MU')
       lambdac_ifld  = pbuf_get_index('LAMBDAC')
       iciwp_ifld    = pbuf_get_index('ICIWP')
       iclwp_ifld    = pbuf_get_index('ICLWP')
       icswp_ifld    = pbuf_get_index('ICSWP')
       cldfsnow_ifld = pbuf_get_index('CLDFSNOW')
    endif

    call addfld (lndfrc_fldn, 'fraction', 1,    rad_data_avgflag,&
         'radiation input: land fraction',phys_decomp)
    call addfld (icefrc_fldn, 'fraction', 1,    rad_data_avgflag,&
         'radiation input: ice fraction',phys_decomp)
    call addfld (snowh_fldn,  'm',        1,    rad_data_avgflag,&
         'radiation input: water equivalent snow depth',phys_decomp)
    call addfld (landm_fldn,  'none',     1,    rad_data_avgflag,&
         'radiation input: land mask: ocean(0), continent(1), transition(0-1)',phys_decomp)

    call addfld (asdir_fldn,  '1',        1,    rad_data_avgflag,&
         'radiation input: short wave direct albedo',phys_decomp, flag_xyfill=.true.)
    call addfld (asdif_fldn,  '1',        1,    rad_data_avgflag,&
         'radiation input: short wave difuse albedo',phys_decomp, flag_xyfill=.true.)
    call addfld (aldir_fldn,  '1',        1,    rad_data_avgflag,&
         'radiation input: long wave direct albedo', phys_decomp, flag_xyfill=.true.)
    call addfld (aldif_fldn,  '1',        1,    rad_data_avgflag,&
         'radiation input: long wave difuse albedo', phys_decomp, flag_xyfill=.true.)

    call addfld (coszen_fldn,     '1', 1,    rad_data_avgflag,&
         'radiation input: cosine solar zenith when positive', phys_decomp, flag_xyfill=.true.)
    call addfld (asdir_pos_fldn,  '1', 1,    rad_data_avgflag,&
         'radiation input: short wave direct albedo weighted by coszen', phys_decomp, flag_xyfill=.true.)
    call addfld (asdif_pos_fldn,  '1', 1,    rad_data_avgflag,&
         'radiation input: short wave difuse albedo weighted by coszen', phys_decomp, flag_xyfill=.true.)
    call addfld (aldir_pos_fldn,  '1', 1,    rad_data_avgflag,&
         'radiation input: long wave direct albedo weighted by coszen', phys_decomp, flag_xyfill=.true.)
    call addfld (aldif_pos_fldn,  '1', 1,    rad_data_avgflag,&
         'radiation input: long wave difuse albedo weighted by coszen', phys_decomp, flag_xyfill=.true.)
    
    call addfld (lwup_fldn,   'W/m2',     1,    rad_data_avgflag,&
         'radiation input: long wave up radiation flux ',phys_decomp)
    call addfld (ts_fldn,     'K',        1,    rad_data_avgflag,&
         'radiation input: surface temperature',phys_decomp)

    call addfld (temp_fldn,   'K',        pver, rad_data_avgflag,&
         'radiation input: midpoint temperature',phys_decomp)
    call addfld (pdel_fldn,   'Pa',       pver, rad_data_avgflag,&
         'radiation input: pressure layer thickness',phys_decomp)
    call addfld (pdeldry_fldn,'Pa',       pver, rad_data_avgflag,&
         'radiation input: dry pressure layer thickness',phys_decomp)
    call addfld (pmid_fldn,   'Pa',       pver, rad_data_avgflag,&
         'radiation input: midpoint pressure',phys_decomp)
    call addfld (watice_fldn, 'kg/kg',    pver, rad_data_avgflag,&
         'radiation input: cloud ice',phys_decomp)
    call addfld (watliq_fldn, 'kg/kg',    pver, rad_data_avgflag,&
         'radiation input: cloud liquid water',phys_decomp)
    call addfld (watvap_fldn, 'kg/kg',    pver, rad_data_avgflag,&
         'radiation input: water vapor',phys_decomp)

    call addfld (zint_fldn,   'km',       pverp,rad_data_avgflag,&
         'radiation input: interface height',phys_decomp)
    call addfld (pint_fldn,   'Pa',       pverp,rad_data_avgflag,&
         'radiation input: interface pressure',phys_decomp)

    call addfld (cld_fldn,    'fraction', pver, rad_data_avgflag,&
         'radiation input: cloud fraction',phys_decomp)
    call addfld (concld_fldn, 'fraction', pver, rad_data_avgflag,&
         'radiation input: convective cloud fraction',phys_decomp)
    call addfld (rel_fldn,    'micron',   pver, rad_data_avgflag,&
         'radiation input: effective liquid drop radius',phys_decomp)
    call addfld (rei_fldn,    'micron',   pver, rad_data_avgflag,&
         'radiation input: effective ice partical radius',phys_decomp)
    
    if (mg_microphys) then
       call addfld (dei_fldn,    'micron',   pver, rad_data_avgflag,&
            'radiation input: effective ice partical diameter',phys_decomp)
       call addfld (des_fldn,    'micron',   pver, rad_data_avgflag,&
            'radiation input: effective snow partical diameter',phys_decomp)
       call addfld (mu_fldn,     ' ',        pver, rad_data_avgflag,&
            'radiation input: ice gamma parameter for optics (radiation)',phys_decomp)
       call addfld (lambdac_fldn,' ',        pver, rad_data_avgflag,&
            'radiation input: slope of droplet distribution for optics (radiation)',phys_decomp)
       call addfld (iciwp_fldn,  'kg/m2',    pver, rad_data_avgflag,&
            'radiation input: In-cloud ice water path',phys_decomp)
       call addfld (iclwp_fldn,  'kg/m2',    pver, rad_data_avgflag,&
            'radiation input: In-cloud liquid water path',phys_decomp)
       call addfld (icswp_fldn,  'kg/m2',    pver, rad_data_avgflag,&
            'radiation input: In-cloud snow water path',phys_decomp)
       call addfld (cldfsnow_fldn, 'fraction', pver, rad_data_avgflag,&
            'radiation input: cloud liquid drops + snow',phys_decomp)
    endif

    call add_default (lndfrc_fldn,    rad_data_histfile_num, ' ')
    call add_default (icefrc_fldn,    rad_data_histfile_num, ' ')
    call add_default (snowh_fldn,     rad_data_histfile_num, ' ')
    call add_default (landm_fldn,     rad_data_histfile_num, ' ')
    call add_default (asdir_fldn,     rad_data_histfile_num, ' ')
    call add_default (asdif_fldn,     rad_data_histfile_num, ' ')
    call add_default (aldir_fldn,     rad_data_histfile_num, ' ')
    call add_default (aldif_fldn,     rad_data_histfile_num, ' ')

    call add_default (coszen_fldn,    rad_data_histfile_num, ' ')
    call add_default (asdir_pos_fldn, rad_data_histfile_num, ' ')
    call add_default (asdif_pos_fldn, rad_data_histfile_num, ' ')
    call add_default (aldir_pos_fldn, rad_data_histfile_num, ' ')
    call add_default (aldif_pos_fldn, rad_data_histfile_num, ' ')

    call add_default (lwup_fldn,      rad_data_histfile_num, ' ')
    call add_default (ts_fldn,        rad_data_histfile_num, ' ')
    call add_default (temp_fldn,      rad_data_histfile_num, ' ')
    call add_default (pdel_fldn,      rad_data_histfile_num, ' ')
    call add_default (pdeldry_fldn,   rad_data_histfile_num, ' ')
    call add_default (pmid_fldn,      rad_data_histfile_num, ' ')
    call add_default (watice_fldn,    rad_data_histfile_num, ' ')
    call add_default (watliq_fldn,    rad_data_histfile_num, ' ')
    call add_default (watvap_fldn,    rad_data_histfile_num, ' ')
    call add_default (zint_fldn,      rad_data_histfile_num, ' ')
    call add_default (pint_fldn,      rad_data_histfile_num, ' ')

    call add_default (cld_fldn,       rad_data_histfile_num, ' ')
    call add_default (concld_fldn,    rad_data_histfile_num, ' ')
    call add_default (rel_fldn,       rad_data_histfile_num, ' ')
    call add_default (rei_fldn,       rad_data_histfile_num, ' ')
    
    if (mg_microphys) then
       call add_default (dei_fldn,       rad_data_histfile_num, ' ')
       call add_default (des_fldn,       rad_data_histfile_num, ' ')
       call add_default (mu_fldn,        rad_data_histfile_num, ' ')
       call add_default (lambdac_fldn,   rad_data_histfile_num, ' ')
       call add_default (iciwp_fldn,     rad_data_histfile_num, ' ')
       call add_default (iclwp_fldn,     rad_data_histfile_num, ' ')
       call add_default (icswp_fldn,     rad_data_histfile_num, ' ')
       call add_default (cldfsnow_fldn,  rad_data_histfile_num, ' ')
    endif

    ! rad constituents

    call rad_cnst_get_info(0, ngas=ngas, naero=naer)
    long_name_description = ' mass mixing ratio used in rad climate calculation'

    ! The code to output the gases assumes that the rad_constituents module has
    ! ordered them in the same way that they are ordered in the "gaslist" array
    ! in module radconstants, and that there are nradgas of them.  This ordering 
    ! is performed in the internal init_lists routine in rad_constituents.
    if (ngas /= nradgas) then
       call endrun('init_rad_data: ERROR: ngas /= nradgas')
    end if

    allocate( gasnames(ngas) )
    call rad_cnst_get_info(0, gasnames=gasnames)

    do i = 1, ngas
       long_name = trim(gasnames(i))//trim(long_name_description)
       name = 'rad_'//gasnames(i)
       call addfld(trim(name), 'kg/kg', pver, rad_data_avgflag, trim(long_name), phys_decomp)
       call add_default (trim(name), rad_data_histfile_num, ' ')
    end do

    if (naer > 0) then
       allocate( aernames(naer) )
       call rad_cnst_get_info(0, aernames=aernames)

       do i = 1, naer
          long_name = trim(aernames(i))//trim(long_name_description)
          name = 'rad_'//aernames(i)
          call addfld(trim(name), 'kg/kg', pver, rad_data_avgflag, trim(long_name), phys_decomp)
          call add_default (trim(name), rad_data_histfile_num, ' ')
       end do
    end if

  end subroutine init_rad_data

  !================================================================================================
  !================================================================================================
  subroutine output_rad_data(  pbuf, state, cam_in, landm, coszen )

    use physics_types,    only: physics_state
    use camsrfexch,       only: cam_in_t     
    
    use constituents,     only: cnst_get_ind
    use physics_buffer, only : physics_buffer_desc, pbuf_get_field, pbuf_old_tim_idx
    implicit none
    type(physics_buffer_desc), pointer :: pbuf(:)
    
    type(physics_state), intent(in), target :: state
    type(cam_in_t),      intent(in) :: cam_in
    real(r8),            intent(in) :: landm(pcols)
    real(r8),            intent(in) :: coszen(pcols)

    ! Local variables
    integer :: i
    character(len=32) :: name
    real(r8), pointer :: mmr(:,:)

    integer :: lchnk, itim_old, ifld
    integer :: ixcldice              ! cloud ice water index
    integer :: ixcldliq              ! cloud liquid water index
    integer :: icol
    integer :: ncol

    ! surface albedoes weighted by (positive cosine zenith angle)
    real(r8):: coszrs_pos(pcols)    ! = max(coszrs,0)
    real(r8):: asdir_pos (pcols)    !
    real(r8):: asdif_pos (pcols)    !
    real(r8):: aldir_pos (pcols)    !
    real(r8):: aldif_pos (pcols)    !

    real(r8), pointer, dimension(:,:)  :: ptr

    if (.not.rad_data_output) return

    ! get index of (liquid+ice) cloud water
    call cnst_get_ind('CLDICE', ixcldice)
    call cnst_get_ind('CLDLIQ', ixcldliq)

    lchnk = state%lchnk
    ncol = state%ncol

    do icol = 1, ncol
       coszrs_pos(icol)  = max(coszen(icol),0._r8)
    enddo
    asdir_pos(:ncol)  = cam_in%asdir(:ncol) * coszrs_pos(:ncol)
    asdif_pos(:ncol)  = cam_in%asdif(:ncol) * coszrs_pos(:ncol)
    aldir_pos(:ncol)  = cam_in%aldir(:ncol) * coszrs_pos(:ncol)
    aldif_pos(:ncol)  = cam_in%aldif(:ncol) * coszrs_pos(:ncol)

    call outfld(lndfrc_fldn, cam_in%landfrac,  pcols, lchnk)
    call outfld(icefrc_fldn, cam_in%icefrac,   pcols, lchnk)
    call outfld(snowh_fldn,  cam_in%snowhland, pcols, lchnk)
    call outfld(landm_fldn,  landm,            pcols, lchnk)
    call outfld(temp_fldn,   state%t,               pcols, lchnk   )
    call outfld(pdel_fldn,   state%pdel,            pcols, lchnk   )
    call outfld(pdeldry_fldn,state%pdeldry,         pcols, lchnk   )
    call outfld(watice_fldn, state%q(:,:,ixcldice), pcols, lchnk   )
    call outfld(watliq_fldn, state%q(:,:,ixcldliq), pcols, lchnk   )
    call outfld(watvap_fldn, state%q(:,:,1),        pcols, lchnk   )
    call outfld(zint_fldn,   state%zi,              pcols, lchnk   )
    call outfld(pint_fldn,   state%pint,            pcols, lchnk   )
    call outfld(pmid_fldn,   state%pmid,            pcols, lchnk   )

    call outfld(asdir_fldn, cam_in%asdir, pcols, lchnk   )
    call outfld(asdif_fldn, cam_in%asdif, pcols, lchnk   )
    call outfld(aldir_fldn, cam_in%aldir, pcols, lchnk   )
    call outfld(aldif_fldn, cam_in%aldif, pcols, lchnk   )

    call outfld(coszen_fldn, coszrs_pos, pcols, lchnk   )
    call outfld(asdir_pos_fldn, asdir_pos, pcols, lchnk   )
    call outfld(asdif_pos_fldn, asdif_pos, pcols, lchnk   )
    call outfld(aldir_pos_fldn, aldir_pos, pcols, lchnk   )
    call outfld(aldif_pos_fldn, aldif_pos, pcols, lchnk   )

    call outfld(lwup_fldn,  cam_in%lwup,  pcols, lchnk   )
    call outfld(ts_fldn,    cam_in%ts,    pcols, lchnk   )

    itim_old = pbuf_old_tim_idx()

    call pbuf_get_field(pbuf, cld_ifld,    ptr, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
    call outfld(cld_fldn,    ptr, pcols, lchnk )

    call pbuf_get_field(pbuf, concld_ifld, ptr, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
    call outfld(concld_fldn, ptr, pcols, lchnk )

    call pbuf_get_field(pbuf, rel_ifld,    ptr )
    call outfld(rel_fldn,    ptr, pcols, lchnk )

    call pbuf_get_field(pbuf, rei_ifld,    ptr )
    call outfld(rei_fldn,    ptr, pcols, lchnk )

    if (mg_microphys) then

       call pbuf_get_field(pbuf,  dei_ifld, ptr     )
       call outfld(dei_fldn,      ptr, pcols, lchnk   )       

       call pbuf_get_field(pbuf,  des_ifld, ptr     )
       call outfld(des_fldn,      ptr, pcols, lchnk   )       

       call pbuf_get_field(pbuf,  mu_ifld, ptr      )
       call outfld(mu_fldn,       ptr, pcols, lchnk   ) 

       call pbuf_get_field(pbuf,  lambdac_ifld, ptr )
       call outfld(lambdac_fldn,  ptr, pcols, lchnk   )       

       call pbuf_get_field(pbuf,  iciwp_ifld, ptr   )
       call outfld(iciwp_fldn,    ptr, pcols, lchnk   )       

       call pbuf_get_field(pbuf,  iclwp_ifld, ptr   )
       call outfld(iclwp_fldn,    ptr, pcols, lchnk   )       

       call pbuf_get_field(pbuf,  icswp_ifld, ptr   )
       call outfld(icswp_fldn,    ptr, pcols, lchnk   )       

       call pbuf_get_field(pbuf,  cldfsnow_ifld, ptr, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
       call outfld(cldfsnow_fldn, ptr, pcols, lchnk   )

    endif

    ! output mixing ratio of rad constituents 

    do i = 1, ngas
       name = 'rad_'//gasnames(i)
       call rad_cnst_get_gas(0, gaslist(i), state, pbuf, mmr)
       call outfld(trim(name), mmr, pcols, lchnk)
    end do

    do i = 1, naer
       name = 'rad_'//aernames(i)
       call rad_cnst_get_aer_mmr(0, i, state, pbuf, mmr)
       call outfld(trim(name), mmr, pcols, lchnk)
    end do

  end subroutine output_rad_data


end module radiation_data
