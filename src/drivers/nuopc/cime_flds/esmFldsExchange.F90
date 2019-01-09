module esmFldsExchange_mod

  !---------------------------------------------------------------------
  ! This is a mediator specific routine that determines ALL possible
  ! fields exchanged between components and their associated routing,
  ! mapping and merging
  !---------------------------------------------------------------------

  use ESMF
  use NUOPC
  use med_constants_mod     , only : CX, CS, CL
  use shr_nuopc_fldList_mod , only : shr_nuopc_fldList_type
  use shr_nuopc_fldList_mod , only : mapbilnr, mapconsf, mapconsd, mappatch, mapfcopy, mapunset, mapfiler
  use shr_nuopc_scalars_mod , only : flds_scalar_name
  use shr_nuopc_fldList_mod , only : addfld => shr_nuopc_fldList_AddFld
  use shr_nuopc_fldList_mod , only : addmap => shr_nuopc_fldList_AddMap
  use shr_nuopc_fldList_mod , only : addmrg => shr_nuopc_fldList_AddMrg
  use shr_nuopc_methods_mod , only : chkerr => shr_nuopc_methods_chkerr
  use shr_nuopc_methods_mod , only : fldchk => shr_nuopc_methods_FB_FldChk
  use med_internalstate_mod , only : InternalState
  use glc_elevclass_mod     , only : glc_elevclass_as_string
  use esmflds  

  implicit none
  public

  public :: esmFldsExchange

  character(*), parameter :: u_FILE_u = &
       __FILE__

!================================================================================
contains
!================================================================================

  subroutine esmFldsExchange(gcomp, phase, rc)

    ! input/output parameters:
    type(ESMF_GridComp)              :: gcomp
    character(len=*) , intent(in)    :: phase
    integer          , intent(inout) :: rc

    ! local variables:
    type(InternalState) :: is_local
    integer             :: ice_ncat                   ! number of sea ice thickness categories
    integer             :: glc_nec                    ! number of land-ice elevation classes
    integer             :: max_megan
    integer             :: max_ddep
    integer             :: max_fire
    logical             :: flds_i2o_per_cat
    integer             :: dbrc
    integer             :: num, i, n
    integer             :: n1, n2, n3, n4
    character(len=4)    :: iso(4)
    character(len=3)    :: cnum
    character(len=CL)   :: cvalue
    character(len=CS)   :: name, fldname
    character(len=CX)   :: atm2ice_fmap, atm2ice_smap, atm2ice_vmap
    character(len=CX)   :: atm2ocn_fmap, atm2ocn_smap, atm2ocn_vmap
    character(len=CX)   :: atm2lnd_fmap, atm2lnd_smap
    character(len=CX)   :: glc2lnd_smap, glc2lnd_fmap
    character(len=CX)   :: glc2ice_rmap
    character(len=CX)   :: glc2ocn_liq_rmap, glc2ocn_ice_rmap
    character(len=CX)   :: ice2atm_fmap, ice2atm_smap
    character(len=CX)   :: ocn2atm_fmap, ocn2atm_smap
    character(len=CX)   :: lnd2atm_fmap, lnd2atm_smap
    character(len=CX)   :: lnd2glc_fmap, lnd2glc_smap
    character(len=CX)   :: lnd2rof_fmap
    character(len=CX)   :: rof2lnd_fmap
    character(len=CX)   :: rof2ocn_fmap, rof2ocn_ice_rmap, rof2ocn_liq_rmap
    character(len=CX)   :: atm2wav_smap, ice2wav_smap, ocn2wav_smap
    character(len=CX)   :: wav2ocn_smap
    logical             :: flds_co2a  ! use case
    logical             :: flds_co2b  ! use case
    logical             :: flds_co2c  ! use case
    character(len=*), parameter :: subname='(esmFlds_Init)'
    !--------------------------------------

    !---------------------------------------
    ! Get the internal state
    !---------------------------------------

    if (phase /= 'advertise') then
       nullify(is_local%wrap)
       call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    !---------------------------
    ! For now hardwire these - these must be less than or equal
    ! to the values in the esmDict.F90 that generates fd.yaml
    !---------------------------

    max_megan = 20
    max_ddep  = 80
    max_fire  = 10
    glc_nec   = 10
    ice_ncat  =  5
    flds_i2o_per_cat = .true.

    iso(1) = ''
    iso(2) = '_16O'
    iso(3) = '_18O'
    iso(4) = '_HDO'

    rc = ESMF_SUCCESS

    call NUOPC_CompAttributeGet(gcomp, name='flds_co2a', value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) flds_co2a
    call ESMF_LogWrite('flds_co2a = '// trim(cvalue), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='flds_co2b', value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) flds_co2b
    call ESMF_LogWrite('flds_co2b = '// trim(cvalue), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='flds_co2c', value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) flds_co2c
    call ESMF_LogWrite('flds_co2c = '// trim(cvalue), ESMF_LOGMSG_INFO, rc=dbrc)

    !----------------------------------------------------------
    ! Initialize mapping file names
    !----------------------------------------------------------

    ! to atm

    call NUOPC_CompAttributeGet(gcomp, name='ice2atm_fmap', value=ice2atm_fmap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('ice2atm_fmap = '// trim(ice2atm_fmap), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='ice2atm_smap', value=ice2atm_smap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('ice2atm_smap = '// trim(ice2atm_smap), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='lnd2atm_fmap', value=lnd2atm_fmap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('lnd2atm_fmap = '// trim(lnd2atm_fmap), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='ocn2atm_smap', value=ocn2atm_smap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('ocn2atm_smap = '// trim(ocn2atm_smap), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='ocn2atm_fmap', value=ocn2atm_fmap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('ocn2atm_fmap = '// trim(ocn2atm_fmap), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='lnd2atm_smap', value=lnd2atm_smap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('lnd2atm_smap = '// trim(lnd2atm_smap), ESMF_LOGMSG_INFO, rc=dbrc)

    ! to lnd

    call NUOPC_CompAttributeGet(gcomp, name='atm2lnd_fmap', value=atm2lnd_fmap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('atm2lnd_fmap = '// trim(atm2lnd_fmap), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='atm2lnd_smap', value=atm2lnd_smap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('atm2lnd_smap = '// trim(atm2lnd_smap), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='rof2lnd_fmap', value=rof2lnd_fmap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('rof2lnd_fmap = '// trim(rof2lnd_fmap), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='glc2lnd_fmap', value=glc2lnd_fmap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('glc2lnd_smap = '// trim(glc2lnd_fmap), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='glc2lnd_smap', value=glc2lnd_smap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('glc2lnd_smap = '// trim(glc2lnd_smap), ESMF_LOGMSG_INFO, rc=dbrc)

    ! to ice

    call NUOPC_CompAttributeGet(gcomp, name='atm2ice_fmap', value=atm2ice_fmap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('atm2ice_fmap = '// trim(atm2ice_fmap), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='atm2ice_smap', value=atm2ice_smap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('atm2ice_smap = '// trim(atm2ice_smap), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='atm2ice_vmap', value=atm2ice_vmap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('atm2ice_vmap = '// trim(atm2ice_vmap), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='glc2ice_rmap', value=glc2ice_rmap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('glc2ice_rmap = '// trim(glc2ice_rmap), ESMF_LOGMSG_INFO, rc=dbrc)

    ! to ocn

    call NUOPC_CompAttributeGet(gcomp, name='atm2ocn_fmap', value=atm2ocn_fmap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('atm2ocn_fmap = '// trim(atm2ocn_fmap), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='atm2ocn_smap', value=atm2ocn_smap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('atm2ocn_smap = '// trim(atm2ocn_smap), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='atm2ocn_vmap', value=atm2ocn_vmap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('atm2ocn_vmap = '// trim(atm2ocn_vmap), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='glc2ocn_liq_rmap', value=glc2ocn_liq_rmap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('glc2ocn_liq_rmap = '// trim(glc2ocn_liq_rmap), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='glc2ocn_ice_rmap', value=glc2ocn_ice_rmap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('glc2ocn_ice_rmap = '// trim(glc2ocn_ice_rmap), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='wav2ocn_smap', value=wav2ocn_smap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('wav2ocn_smap = '// trim(wav2ocn_smap), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='rof2ocn_fmap', value=rof2ocn_fmap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('rof2ocn_fmap = '// trim(rof2ocn_fmap), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='rof2ocn_liq_rmap', value=rof2ocn_liq_rmap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('rof2ocn_liq_rmap = '// trim(rof2ocn_liq_rmap), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='rof2ocn_ice_rmap', value=rof2ocn_ice_rmap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('rof2ocn_ice_rmap = '// trim(rof2ocn_ice_rmap), ESMF_LOGMSG_INFO, rc=dbrc)

    ! to rof

    call NUOPC_CompAttributeGet(gcomp, name='lnd2rof_fmap', value=lnd2rof_fmap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('lnd2rof_fmap = '// trim(lnd2rof_fmap), ESMF_LOGMSG_INFO, rc=dbrc)

    ! to glc

    call NUOPC_CompAttributeGet(gcomp, name='lnd2glc_fmap', value=lnd2glc_fmap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('lnd2glc_fmap = '// trim(lnd2glc_fmap), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='lnd2glc_smap', value=lnd2glc_smap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('lnd2glc_smap = '// trim(lnd2glc_smap), ESMF_LOGMSG_INFO, rc=dbrc)

    ! to wav

    call NUOPC_CompAttributeGet(gcomp, name='atm2wav_smap', value=atm2wav_smap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('atm2wav_smap = '// trim(atm2wav_smap), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='ice2wav_smap', value=ice2wav_smap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('ice2wav_smap = '// trim(ice2wav_smap), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='ocn2wav_smap', value=ocn2wav_smap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('ocn2wav_smap = '// trim(ocn2wav_smap), ESMF_LOGMSG_INFO, rc=dbrc)

    !----------------------------------------------------------
    ! scalar information
    !----------------------------------------------------------

    if (phase == 'advertise') then
       do n = 1,ncomps
          call addfld(fldListFr(n)%flds, trim(flds_scalar_name))
          call addfld(fldListTo(n)%flds, trim(flds_scalar_name))
       end do
    end if

    !----------------------------------------------------------
    ! Masks from components
    !----------------------------------------------------------

    if (phase == 'advertise') then
       call addfld(fldListFr(complnd)%flds, 'Sl_lfrin')
       call addfld(fldListFr(compocn)%flds, 'So_omask')
       call addfld(fldListFr(compice)%flds, 'Si_imask')
    else
       call addmap(fldListFr(compocn)%flds, 'So_omask', compice,  mapfcopy, 'unset', 'unset')
    end if

    !----------------------------------------------------------
    ! to atm: Fractions
    !----------------------------------------------------------

    if (phase == 'advertise') then
       call addfld(fldListTo(compatm)%flds, 'Sl_lfrac')
       call addfld(fldListTo(compatm)%flds, 'Si_ifrac')
       call addfld(fldListTo(compatm)%flds, 'So_ofrac')
    end if

    !----------------------------------------------------------
    ! From ice: Fractional ice coverage wrt ocean
    !----------------------------------------------------------

    if (phase == 'advertise') then
       call addfld(fldListFr(compice)%flds, 'Si_ifrac')
       call addfld(fldListTo(compocn)%flds, 'Si_ifrac')
       call addfld(fldListTo(compwav)%flds, 'Si_ifrac')
    else
       call addmap(fldListFr(compice)%flds, 'Si_ifrac', compocn,  mapfcopy, 'unset', 'unset')

       call addmrg(fldListTo(compocn)%flds, 'Si_ifrac', mrg_from1=compice, mrg_fld1='Si_ifrac', mrg_type1='copy')
       call addmrg(fldListTo(compwav)%flds, 'Si_ifrac', mrg_from1=compice, mrg_fld1='Si_ifrac', mrg_type1='copy')
    end if

    ! ---------------------------------------------------------------------
    ! From atm: Height at the lowest model level
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Sa_z')
       call addfld(fldListTo(complnd)%flds, 'Sa_z')
       call addfld(fldListTo(compice)%flds, 'Sa_z')
    else
       call addmap(fldListFr(compatm)%flds, 'Sa_z', complnd, mapbilnr, 'one', atm2lnd_smap)
       call addmap(fldListFr(compatm)%flds, 'Sa_z', compice, mapbilnr, 'one', atm2ice_smap)
       call addmap(fldListFr(compatm)%flds, 'Sa_z', compocn, mapbilnr, 'one', atm2ocn_smap)

       call addmrg(fldListTo(complnd)%flds, 'Sa_z', mrg_from1=compatm, mrg_fld1='Sa_z', mrg_type1='copy')
       call addmrg(fldListTo(compice)%flds, 'Sa_z', mrg_from1=compatm, mrg_fld1='Sa_z', mrg_type1='copy')
    end if

    ! ---------------------------------------------------------------------
    ! From atm: Surface height
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Sa_topo')
       call addfld(fldListTo(complnd)%flds, 'Sa_topo')
    else
       call addmap(fldListFr(compatm)%flds, 'Sa_topo', complnd, mapbilnr, 'one', atm2lnd_smap)
       call addmrg(fldListTo(complnd)%flds, 'Sa_topo', mrg_from1=compatm, mrg_fld1='Sa_topo', mrg_type1='copy')
    end if

    ! ---------------------------------------------------------------------
    ! From atm: zonal wind at the lowest model level
    ! ---------------------------------------------------------------------
    ! (cesm only) Sa_u and Sa_v are mapped to the ocean grid in the
    ! mediator - BUT are not sent to the ocean - They are only used in
    ! the atm/ocn flux calculation - a special mapping will be done in
    ! the mediator for these fields -

    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Sa_u')
       call addfld(fldListTo(complnd)%flds, 'Sa_u')
       call addfld(fldListTo(compice)%flds, 'Sa_u')
       call addfld(fldListTo(compwav)%flds, 'Sa_u')
    else
       call addmap(fldListFr(compatm)%flds, 'Sa_u', complnd, mapbilnr, 'one', atm2lnd_smap)
       call addmap(fldListFr(compatm)%flds, 'Sa_u', compice, mappatch, 'one', atm2ice_vmap)
       call addmap(fldListFr(compatm)%flds, 'Sa_u', compocn, mappatch, 'one', atm2ocn_vmap)
       call addmap(fldListFr(compatm)%flds, 'Sa_u', compwav, mapbilnr, 'one', atm2wav_smap)

       call addmrg(fldListTo(complnd)%flds, 'Sa_u', mrg_from1=compatm, mrg_fld1='Sa_u', mrg_type1='copy')
       call addmrg(fldListTo(compice)%flds, 'Sa_u', mrg_from1=compatm, mrg_fld1='Sa_u', mrg_type1='copy')
       call addmrg(fldListTo(compwav)%flds, 'Sa_u', mrg_from1=compatm, mrg_fld1='Sa_u', mrg_type1='copy')
    end if

    ! ---------------------------------------------------------------------
    ! From atm: Meridional wind at the lowest model level
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Sa_v')
       call addfld(fldListTo(complnd)%flds, 'Sa_v')
       call addfld(fldListTo(compice)%flds, 'Sa_v')
       call addfld(fldListTo(compwav)%flds, 'Sa_v')
    else
       call addmap(fldListFr(compatm)%flds, 'Sa_v', complnd, mapbilnr, 'one', atm2lnd_smap)
       call addmap(fldListFr(compatm)%flds, 'Sa_v', compice, mappatch, 'one', atm2ice_vmap)
       call addmap(fldListFr(compatm)%flds, 'Sa_v', compocn, mappatch, 'one', atm2ocn_vmap)
       call addmap(fldListFr(compatm)%flds, 'Sa_v', compwav, mapbilnr, 'one', atm2wav_smap)

       call addmrg(fldListTo(complnd)%flds, 'Sa_v', mrg_from1=compatm, mrg_fld1='Sa_v', mrg_type1='copy')
       call addmrg(fldListTo(compice)%flds, 'Sa_v', mrg_from1=compatm, mrg_fld1='Sa_v', mrg_type1='copy')
       call addmrg(fldListTo(compwav)%flds, 'Sa_v', mrg_from1=compatm, mrg_fld1='Sa_v', mrg_type1='copy')
    end if

    ! ---------------------------------------------------------------------
    ! From atm: Temperature at the lowest model level
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Sa_tbot')
       call addfld(fldListTo(complnd)%flds, 'Sa_tbot')
       call addfld(fldListTo(compice)%flds, 'Sa_tbot')
       call addfld(fldListTo(compwav)%flds, 'Sa_tbot')
    else
       call addmap(fldListFr(compatm)%flds, 'Sa_tbot', complnd, mapbilnr, 'one', atm2lnd_smap)
       call addmap(fldListFr(compatm)%flds, 'Sa_tbot', compice, mapbilnr, 'one', atm2ice_smap)
       call addmap(fldListFr(compatm)%flds, 'Sa_tbot', compocn, mapbilnr, 'one', atm2ocn_smap)
       call addmap(fldListFr(compatm)%flds, 'Sa_tbot', compwav, mapbilnr, 'one', atm2wav_smap)

       call addmrg(fldListTo(complnd)%flds, 'Sa_tbot', mrg_from1=compatm, mrg_fld1='Sa_tbot', mrg_type1='copy')
       call addmrg(fldListTo(compice)%flds, 'Sa_tbot', mrg_from1=compatm, mrg_fld1='Sa_tbot', mrg_type1='copy')
       call addmrg(fldListTo(compwav)%flds, 'Sa_tbot', mrg_from1=compatm, mrg_fld1='Sa_tbot', mrg_type1='copy')
    end if

    ! ---------------------------------------------------------------------
    ! to lnd and to ice: potential temperature at the lowest model level
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds , 'Sa_ptem')
       call addfld(fldListTo(complnd)%flds , 'Sa_ptem')
       call addfld(fldListTo(compice)%flds , 'Sa_ptem')
       call addfld(fldListMed_aoflux_a%flds, 'Sa_ptem')
       call addfld(fldListMed_aoflux_o%flds, 'Sa_ptem')
    else
       if (fldchk(is_local%wrap%FBImp(compatm,compatm), 'Sa_ptem', rc=rc) .and. &
           fldchk(is_local%wrap%FBMed_aoflux_o        , 'Sa_ptem', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Sa_ptem', compocn, mapbilnr, 'one', atm2ocn_smap)
       end if
       if (fldchk(is_local%wrap%FBImp(compatm,compatm), 'Sa_ptem', rc=rc) .and. &
           fldchk(is_local%wrap%FBExp(complnd)        , 'Sa_ptem', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Sa_ptem', complnd, mapbilnr, 'one', atm2lnd_smap)
          call addmrg(fldListTo(complnd)%flds, 'Sa_ptem', mrg_from1=compatm, mrg_fld1='Sa_ptem', mrg_type1='copy')
       end if
       if (fldchk(is_local%wrap%FBImp(compatm,compatm), 'Sa_ptem', rc=rc) .and. &
           fldchk(is_local%wrap%FBExp(compice)        , 'Sa_ptem', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Sa_ptem', compice, mapbilnr, 'one', atm2ice_smap)
          call addmrg(fldListTo(compice)%flds, 'Sa_ptem', mrg_from1=compatm, mrg_fld1='Sa_ptem', mrg_type1='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to lnd and ice: Specific humidity at the lowest model level
    ! ---------------------------------------------------------------------
    do n = 1,size(iso)
       if (phase == 'advertise') then
          call addfld(fldListFr(compatm)%flds , 'Sa_shum'//iso(n))
          call addfld(fldListTo(complnd)%flds , 'Sa_shum'//iso(n))
          call addfld(fldListTo(compice)%flds , 'Sa_shum'//iso(n))
          call addfld(fldListMed_aoflux_a%flds, 'Sa_shum'//iso(n))
          call addfld(fldListMed_aoflux_o%flds, 'Sa_shum'//iso(n))
       else
          if (fldchk(is_local%wrap%FBImp(compatm,compatm), 'Sa_shum'//iso(n), rc=rc) .and. &
              fldchk(is_local%wrap%FBMed_aoflux_o        , 'Sa_shum'//iso(n), rc=rc)) then
             call addmap(fldListFr(compatm)%flds, 'Sa_shum', compocn, mapbilnr, 'one', atm2ocn_smap)
          end if
          if (fldchk(is_local%wrap%FBImp(compatm,compatm), 'Sa_shum'//iso(n), rc=rc) .and. &
              fldchk(is_local%wrap%FBExp(complnd)        , 'Sa_shum'//iso(n), rc=rc)) then
             call addmap(fldListFr(compatm)%flds, 'Sa_shum'//iso(n), complnd, mapbilnr, 'one', atm2lnd_smap)
             call addmrg(fldListTo(complnd)%flds, 'Sa_shum'//iso(n), &
                  mrg_from1=compatm, mrg_fld1='Sa_shum'//iso(n), mrg_type1='copy')
          end if
          if (fldchk(is_local%wrap%FBImp(compatm,compatm), 'Sa_shum'//iso(n), rc=rc) .and. &
              fldchk(is_local%wrap%FBExp(compice)        , 'Sa_shum'//iso(n), rc=rc)) then
             call addmap(fldListFr(compatm)%flds, 'Sa_shum'//iso(n), compice, mapbilnr, 'one', atm2ice_smap)
             call addmrg(fldListTo(compice)%flds, 'Sa_shum'//iso(n), &
                  mrg_from1=compatm, mrg_fld1='Sa_shum'//iso(n), mrg_type1='copy')
          end if
       end if
    end do

    ! ---------------------------------------------------------------------
    ! From atm: Pressure at the lowest model level
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Sa_pbot')
       call addfld(fldListTo(complnd)%flds, 'Sa_pbot')
       call addfld(fldListTo(compice)%flds, 'Sa_pbot')
    else
       call addmap(fldListFr(compatm)%flds, 'Sa_pbot', complnd, mapbilnr, 'one', atm2lnd_smap)
       call addmap(fldListFr(compatm)%flds, 'Sa_pbot', compice, mapbilnr, 'one', atm2ice_smap)
       call addmap(fldListFr(compatm)%flds, 'Sa_pbot', compocn, mapbilnr, 'one', atm2ocn_smap)

       call addmrg(fldListTo(complnd)%flds, 'Sa_pbot', mrg_from1=compatm, mrg_fld1='Sa_pbot', mrg_type1='copy')
       call addmrg(fldListTo(compice)%flds, 'Sa_pbot', mrg_from1=compatm, mrg_fld1='Sa_pbot', mrg_type1='copy')
    end if

    ! ---------------------------------------------------------------------
    ! From atm: Density at the lowest model level
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Sa_dens')
       call addfld(fldListTo(compice)%flds, 'Sa_dens')
    else
       call addmap(fldListFr(compatm)%flds, 'Sa_dens', compice, mapbilnr, 'one', atm2ice_smap)
       call addmap(fldListFr(compatm)%flds,' Sa_dens', compocn, mapbilnr, 'one', atm2ocn_smap)

       call addmrg(fldListTo(compice)%flds, 'Sa_dens', mrg_from1=compatm, mrg_fld1='Sa_dens', mrg_type1='copy')
    end if

    ! ---------------------------------------------------------------------
    !  From atm: Convective and large scale precipitation rate water equivalent
    ! ---------------------------------------------------------------------
    do n = 1,size(iso)
       if (phase == 'advertise') then
          call addfld(fldListFr(compatm)%flds, 'Faxa_rainc'//iso(n))
          call addfld(fldListFr(compatm)%flds, 'Faxa_rainl'//iso(n))
          call addfld(fldListFr(compatm)%flds, 'Faxa_rain' //iso(n))
          call addfld(fldListTo(complnd)%flds, 'Faxa_rainc'//iso(n))
          call addfld(fldListTo(complnd)%flds, 'Faxa_rainl'//iso(n))
          call addfld(fldListTo(compocn)%flds, 'Faxa_rain' //iso(n))
          call addfld(fldListTo(compice)%flds, 'Faxa_rain' //iso(n))
       else
          call addmap(fldListFr(compatm)%flds, 'Faxa_rainl'//iso(n), compocn, mapconsf, 'one', atm2ocn_fmap)
          call addmrg(fldListTo(compocn)%flds, 'Faxa_rain'//iso(n) , &
               mrg_from1=compatm, mrg_fld1='Faxa_rainc'//iso(n)//':'//'Faxa_rainl'//iso(n), &
               mrg_type1='sum_with_weights', mrg_fracname1='ofrac')
       end if
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rainl'//iso(n), rc=rc) .and.  &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rainc'//iso(n), rc=rc) .and. &
            fldchk(is_local%wrap%FBexp(compice)        , 'Faxa_rain' //iso(n), rc=rc)) then

          call addmap(fldListFr(compatm)%flds, 'Faxa_rainc'//iso(n), compice, mapconsf, 'one', atm2ice_fmap)
          call addmap(fldListFr(compatm)%flds, 'Faxa_rainl'//iso(n), compice, mapconsf, 'one', atm2ice_fmap)
          call addmrg(fldListTo(compice)%flds, 'Faxa_rain'//iso(n) , &
               mrg_from1=compatm, mrg_fld1='Faxa_rainc'//iso(n)//':'//'Faxa_rainl'//iso(n), &
               mrg_type1='sum')
       end if
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rainl'//iso(n), rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rainc'//iso(n), rc=rc) .and. &
            fldchk(is_local%wrap%FBexp(complnd)        , 'Faxa_rainl'//iso(n), rc=rc) .and. &
            fldchk(is_local%wrap%FBexp(complnd)        , 'Faxa_rainl'//iso(n), rc=rc)) then

          call addmap(fldListFr(compatm)%flds, 'Faxa_rainc'//iso(n), complnd, mapconsf, 'one', atm2lnd_fmap)
          call addmap(fldListFr(compatm)%flds, 'Faxa_rainl'//iso(n), complnd, mapconsf, 'one', atm2lnd_fmap)
          call addmrg(fldListTo(complnd)%flds, 'Faxa_rainc'//iso(n), &
               mrg_from1=compatm, mrg_fld1='Faxa_rainc'//iso(n), mrg_type1='copy')
          call addmrg(fldListTo(complnd)%flds, 'Faxa_rainl'//iso(n), &
               mrg_from1=compatm, mrg_fld1='Faxa_rainl'//iso(n), mrg_type1='copy')
       end if
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rain'//iso(n), rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compocn)        , 'Faxa_rain'//iso(n), rc=rc)) then

          call addmap(fldListFr(compatm)%flds, 'Faxa_rain'//iso(n), compocn, mapconsf, 'one', atm2ocn_fmap)
          call addmrg(fldListTo(compocn)%flds, 'Faxa_rain'//iso(n), mrg_from1=compatm, mrg_fld1='Faxa_rain'//iso(n), &
               mrg_type1='copy_with_weights', mrg_fracname1='ofrac')
       end if
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rain'//iso(n), rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compice)        , 'Faxa_rain'//iso(n), rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_rain'//iso(n), compice, mapconsf, 'one', atm2ice_fmap)
          call addmrg(fldListTo(compice)%flds, 'Faxa_rain'//iso(n), mrg_from1=compatm, mrg_fld1='Faxa_rain'//iso(n), &
               mrg_type1='copy')
       end if
    end do

    ! ---------------------------------------------------------------------
    !  From atm: Convective and large-scale (stable) snow rate
    ! ---------------------------------------------------------------------
    do n = 1,size(iso)
       if (phase == 'advertise') then
          call addfld(fldListFr(compatm)%flds, 'Faxa_snowc'//iso(n))
          call addfld(fldListFr(compatm)%flds, 'Faxa_snowl'//iso(n))
          call addfld(fldListFr(compatm)%flds, 'Faxa_snow' //iso(n))
          call addfld(fldListTo(complnd)%flds, 'Faxa_snowc'//iso(n))
          call addfld(fldListTo(complnd)%flds, 'Faxa_snowl'//iso(n))
          call addfld(fldListTo(compocn)%flds, 'Faxa_snow'//iso(n))
          call addfld(fldListTo(compice)%flds, 'Faxa_snow'//iso(n))
       else
          if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_snowl'//iso(n), rc=rc) .and.  &
               fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_snowc'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBexp(compocn)        , 'Faxa_snow' //iso(n), rc=rc)) then

             call addmap(fldListFr(compatm)%flds, 'Faxa_snowc'//iso(n), compocn, mapconsf, 'one', atm2ocn_fmap)
             call addmap(fldListFr(compatm)%flds, 'Faxa_snowl'//iso(n), compocn, mapconsf, 'one', atm2ocn_fmap)
             call addmrg(fldListTo(compocn)%flds, 'Faxa_snow'//iso(n) , &
                  mrg_from1=compatm, mrg_fld1='Faxa_snowc'//iso(n)//':'//'Faxa_snowl'//iso(n), &
                  mrg_type1='sum_with_weights', mrg_fracname1='ofrac')
          end if
          if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_snowl'//iso(n), rc=rc) .and.  &
               fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_snowc'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBexp(compice)        , 'Faxa_snow' //iso(n), rc=rc)) then

             call addmap(fldListFr(compatm)%flds, 'Faxa_snowc'//iso(n), compice, mapconsf, 'one', atm2ice_fmap)
             call addmap(fldListFr(compatm)%flds, 'Faxa_snowl'//iso(n), compice, mapconsf, 'one', atm2ice_fmap)
             call addmrg(fldListTo(compice)%flds, 'Faxa_snow'//iso(n) , &
                  mrg_from1=compatm, mrg_fld1='Faxa_snowc'//iso(n)//':'//'Faxa_snowl'//iso(n), mrg_type1='sum')
          end if
          if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_snowl'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_snowc'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBexp(complnd)        , 'Faxa_snowl'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBexp(complnd)        , 'Faxa_snowl'//iso(n), rc=rc)) then

             call addmap(fldListFr(compatm)%flds, 'Faxa_snowc'//iso(n), complnd, mapconsf, 'one', atm2lnd_fmap)
             call addmap(fldListFr(compatm)%flds, 'Faxa_snowl'//iso(n), complnd, mapconsf, 'one', atm2lnd_fmap)
             call addmrg(fldListTo(complnd)%flds, 'Faxa_snowc'//iso(n), &
                  mrg_from1=compatm, mrg_fld1='Faxa_snowc'//iso(n), mrg_type1='copy')
             call addmrg(fldListTo(complnd)%flds, 'Faxa_snowl'//iso(n), &
                  mrg_from1=compatm, mrg_fld1='Faxa_snowl'//iso(n), mrg_type1='copy')
          end if
          if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_snow'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBExp(compocn)        , 'Faxa_snow'//iso(n), rc=rc)) then

             call addmap(fldListFr(compatm)%flds, 'Faxa_snow'//iso(n), compocn, mapconsf, 'one', atm2ocn_fmap)
             call addmrg(fldListTo(compocn)%flds, 'Faxa_snow'//iso(n), mrg_from1=compatm, mrg_fld1='Faxa_snow'//iso(n), &
                  mrg_type1='copy_with_weights', mrg_fracname1='ofrac')
          end if
          if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_snow'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBExp(compice)        , 'Faxa_snow'//iso(n), rc=rc)) then
             call addmap(fldListFr(compatm)%flds, 'Faxa_snow'//iso(n), compice, mapconsf, 'one', atm2ice_fmap)
             call addmrg(fldListTo(compice)%flds, 'Faxa_snow'//iso(n), mrg_from1=compatm, mrg_fld1='Faxa_snow'//iso(n), &
                  mrg_type1='copy')
          end if
       end if
    end do

    ! ---------------------------------------------------------------------
    ! From atm: downward longwave heat flux
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_lwdn')
       call addfld(fldListTo(complnd)%flds, 'Faxa_lwdn')
       call addfld(fldListTo(compice)%flds, 'Faxa_lwdn')
       call addfld(fldListTo(compocn)%flds, 'Faxa_lwdn')
    else
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_lwdn', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_lwdn', complnd, mapconsf, 'one', atm2lnd_fmap)
          call addmap(fldListFr(compatm)%flds, 'Faxa_lwdn', compice, mapconsf, 'one', atm2ice_fmap)
          call addmap(fldListFr(compatm)%flds, 'Faxa_lwdn', compocn, mapconsf, 'one', atm2ocn_fmap)

          call addmrg(fldListTo(complnd)%flds, 'Faxa_lwdn', mrg_from1=compatm, mrg_fld1='Faxa_lwdn', mrg_type1='copy')
          call addmrg(fldListTo(compice)%flds, 'Faxa_lwdn', mrg_from1=compatm, mrg_fld1='Faxa_lwdn', mrg_type1='copy')
          call addmrg(fldListTo(compocn)%flds, 'Faxa_lwdn', mrg_from1=compatm, mrg_fld1='Faxa_lwdn', &
               mrg_type1='copy_with_weights', mrg_fracname1='ofrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! From atm: downward direct near-infrared incident solar radiation
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_swndr')
       call addfld(fldListTo(complnd)%flds, 'Faxa_swndr')
       call addfld(fldListTo(compice)%flds, 'Faxa_swndr')
       call addfld(fldListTo(compocn)%flds, 'Faxa_swndr')
    else
       call addmap(fldListFr(compatm)%flds, 'Faxa_swndr', complnd, mapconsf, 'one', atm2lnd_fmap)
       call addmap(fldListFr(compatm)%flds, 'Faxa_swndr', compice, mapconsf, 'one', atm2ice_fmap)
       call addmap(fldListFr(compatm)%flds, 'Faxa_swndr', compocn, mapconsf, 'one', atm2ocn_fmap)

       call addmrg(fldListTo(complnd)%flds, 'Faxa_swndr', mrg_from1=compatm, mrg_fld1='Faxa_swndr', mrg_type1='copy')
       call addmrg(fldListTo(compice)%flds, 'Faxa_swndr', mrg_from1=compatm, mrg_fld1='Faxa_swndr', mrg_type1='copy')
       call addmrg(fldListTo(compocn)%flds, 'Faxa_swndr', mrg_from1=compatm, mrg_fld1='Faxa_swndr', &
            mrg_type1='copy_with_weights', mrg_fracname1='ofrac')
    end if

    ! ---------------------------------------------------------------------
    ! From atm: downward direct visible incident solar radiation
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_swvdr')
       call addfld(fldListTo(complnd)%flds, 'Faxa_swvdr')
       call addfld(fldListTo(compice)%flds, 'Faxa_swvdr')
       call addfld(fldListTo(compocn)%flds, 'Faxa_swvdr')
    else
       call addmap(fldListFr(compatm)%flds, 'Faxa_swvdr', complnd, mapconsf, 'one', atm2lnd_fmap)
       call addmap(fldListFr(compatm)%flds, 'Faxa_swvdr', compice, mapconsf, 'one', atm2ice_fmap)
       call addmap(fldListFr(compatm)%flds, 'Faxa_swvdr', compocn, mapconsf, 'one', atm2ocn_fmap)

       call addmrg(fldListTo(complnd)%flds, 'Faxa_swvdr', mrg_from1=compatm, mrg_fld1='Faxa_swvdr', mrg_type1='copy')
       call addmrg(fldListTo(compice)%flds, 'Faxa_swvdr', mrg_from1=compatm, mrg_fld1='Faxa_swvdr', mrg_type1='copy')
       call addmrg(fldListTo(compocn)%flds, 'Faxa_swvdr', mrg_from1=compatm, mrg_fld1='Faxa_swvdr', &
            mrg_type1='copy_with_weights', mrg_fracname1='ofrac')
    end if

    ! ---------------------------------------------------------------------
    ! From atm: downward diffuse near-infrared incident solar radiation
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_swndf')
       call addfld(fldListTo(complnd)%flds, 'Faxa_swndf')
       call addfld(fldListTo(compice)%flds, 'Faxa_swndf')
       call addfld(fldListTo(compocn)%flds, 'Faxa_swndf')
    else
       call addmap(fldListFr(compatm)%flds, 'Faxa_swndf', complnd, mapconsf, 'one', atm2lnd_fmap)
       call addmap(fldListFr(compatm)%flds, 'Faxa_swndf', compice, mapconsf, 'one', atm2ice_fmap)
       call addmap(fldListFr(compatm)%flds, 'Faxa_swndf', compocn, mapconsf, 'one', atm2ocn_fmap)

       call addmrg(fldListTo(complnd)%flds, 'Faxa_swndf', mrg_from1=compatm, mrg_fld1='Faxa_swndf', mrg_type1='copy')
       call addmrg(fldListTo(compice)%flds, 'Faxa_swndf', mrg_from1=compatm, mrg_fld1='Faxa_swndf', mrg_type1='copy')
       call addmrg(fldListTo(compocn)%flds, 'Faxa_swndf', mrg_from1=compatm, mrg_fld1='Faxa_swndf', &
            mrg_type1='copy_with_weights', mrg_fracname1='ofrac')
    end if

    ! ---------------------------------------------------------------------
    ! From atm: downward Diffuse visible incident solar radiation
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_swvdf')
       call addfld(fldListTo(complnd)%flds, 'Faxa_swvdf')
       call addfld(fldListTo(compice)%flds, 'Faxa_swvdf')
       call addfld(fldListTo(compocn)%flds, 'Faxa_swvdf')
    else
       call addmap(fldListFr(compatm)%flds, 'Faxa_swvdf', complnd, mapconsf, 'one', atm2lnd_fmap)
       call addmap(fldListFr(compatm)%flds, 'Faxa_swvdf', compice, mapconsf, 'one', atm2ice_fmap)
       call addmap(fldListFr(compatm)%flds, 'Faxa_swvdf', compocn, mapconsf, 'one', atm2ocn_fmap)

       call addmrg(fldListTo(complnd)%flds, 'Faxa_swvdf', mrg_from1=compatm, mrg_fld1='Faxa_swvdf', mrg_type1='copy')
       call addmrg(fldListTo(compice)%flds, 'Faxa_swvdf', mrg_from1=compatm, mrg_fld1='Faxa_swvdf', mrg_type1='copy')
       call addmrg(fldListTo(compocn)%flds, 'Faxa_swvdf', mrg_from1=compatm, mrg_fld1='Faxa_swvdf', &
            mrg_type1='copy_with_weights', mrg_fracname1='ofrac')
    end if

    ! ---------------------------------------------------------------------
    ! From atm: hydrophylic black carbon dry deposition flux
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_bcphidry')
       call addfld(fldListTo(complnd)%flds, 'Faxa_bcphidry')
       call addfld(fldListTo(compice)%flds, 'Faxa_bcphidry')
    else
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_bcphidry', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_bcphidry', complnd, mapconsf, 'one', atm2lnd_fmap)
          call addmap(fldListFr(compatm)%flds, 'Faxa_bcphidry', compice, mapconsf, 'one', atm2ice_fmap)
          call addmap(fldListFr(compatm)%flds, 'Faxa_bcphidry', compocn, mapconsf, 'one', atm2ocn_fmap)

          call addmrg(fldListTo(complnd)%flds, 'Faxa_bcphidry', mrg_from1=compatm, mrg_fld1='Faxa_bcphidry', mrg_type1='copy')
          call addmrg(fldListTo(compice)%flds, 'Faxa_bcphidry', mrg_from1=compatm, mrg_fld1='Faxa_bcphidry', mrg_type1='copy')
          call addmrg(fldListTo(compocn)%flds, 'Faxa_bcphidry', mrg_from1=compatm, mrg_fld1='Faxa_bcphidry', &
               mrg_type1='copy_with_weights', mrg_fracname1='ofrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! From atm: hydrophobic black carbon dry deposition flux
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_bcphodry')
       call addfld(fldListTo(complnd)%flds, 'Faxa_bcphodry')
       call addfld(fldListTo(compice)%flds, 'Faxa_bcphodry')
    else
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_bcphodry', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_bcphodry', complnd, mapconsf, 'one', atm2lnd_fmap)
          call addmap(fldListFr(compatm)%flds, 'Faxa_bcphodry', compice, mapconsf, 'one', atm2ice_fmap)
          call addmap(fldListFr(compatm)%flds, 'Faxa_bcphodry', compocn, mapconsf, 'one', atm2ocn_fmap)

          call addmrg(fldListTo(complnd)%flds, 'Faxa_bcphodry', mrg_from1=compatm, mrg_fld1='Faxa_bcphodry', mrg_type1='copy')
          call addmrg(fldListTo(compice)%flds, 'Faxa_bcphodry', mrg_from1=compatm, mrg_fld1='Faxa_bcphodry', mrg_type1='copy')
          call addmrg(fldListTo(compocn)%flds, 'Faxa_bcphodry', mrg_from1=compatm, mrg_fld1='Faxa_bcphodry', &
               mrg_type1='copy_with_weights', mrg_fracname1='ofrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! From atm: hydrophylic black carbon wet deposition flux
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_bcphiwet')
       call addfld(fldListTo(complnd)%flds, 'Faxa_bcphiwet')
       call addfld(fldListTo(compice)%flds, 'Faxa_bcphiwet')
    else
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_bcphiwet', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_bcphiwet', complnd, mapconsf, 'one', atm2lnd_fmap)
          call addmap(fldListFr(compatm)%flds, 'Faxa_bcphiwet', compice, mapconsf, 'one', atm2ice_fmap)
          call addmap(fldListFr(compatm)%flds, 'Faxa_bcphiwet', compocn, mapconsf, 'one', atm2ocn_fmap)

          call addmrg(fldListTo(complnd)%flds, 'Faxa_bcphiwet', mrg_from1=compatm, mrg_fld1='Faxa_bcphiwet', mrg_type1='copy')
          call addmrg(fldListTo(compice)%flds, 'Faxa_bcphiwet', mrg_from1=compatm, mrg_fld1='Faxa_bcphiwet', mrg_type1='copy')
          call addmrg(fldListTo(compocn)%flds, 'Faxa_bcphiwet', mrg_from1=compatm, mrg_fld1='Faxa_bcphiwet', &
               mrg_type1='copy_with_weights', mrg_fracname1='ofrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! From atm: hydrophylic organic carbon dry deposition flux
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_ocphidry')
       call addfld(fldListTo(complnd)%flds, 'Faxa_ocphidry')
       call addfld(fldListTo(compice)%flds, 'Faxa_ocphidry')
    else
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_ocphidry', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_ocphidry', complnd, mapconsf, 'one', atm2lnd_fmap)
          call addmap(fldListFr(compatm)%flds, 'Faxa_ocphidry', compice, mapconsf, 'one', atm2ice_fmap)
          call addmap(fldListFr(compatm)%flds, 'Faxa_ocphidry', compocn, mapconsf, 'one', atm2ocn_fmap)

          call addmrg(fldListTo(complnd)%flds, 'Faxa_ocphidry', mrg_from1=compatm, mrg_fld1='Faxa_ocphidry', mrg_type1='copy')
          call addmrg(fldListTo(compice)%flds, 'Faxa_ocphidry', mrg_from1=compatm, mrg_fld1='Faxa_ocphidry', mrg_type1='copy')
          call addmrg(fldListTo(compocn)%flds, 'Faxa_ocphidry', mrg_from1=compatm, mrg_fld1='Faxa_ocphidry', &
               mrg_type1='copy_with_weights', mrg_fracname1='ofrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! From atm: hydrophobic organic carbon dry deposition flux
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_ocphodry')
       call addfld(fldListTo(complnd)%flds, 'Faxa_ocphodry')
       call addfld(fldListTo(compice)%flds, 'Faxa_ocphodry')
    else
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_ocphodry', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_ocphodry', complnd, mapconsf, 'one', atm2lnd_fmap)
          call addmap(fldListFr(compatm)%flds, 'Faxa_ocphodry', compice, mapconsf, 'one', atm2ice_fmap)
          call addmap(fldListFr(compatm)%flds, 'Faxa_ocphodry', compocn, mapconsf, 'one', atm2ocn_fmap)

          call addmrg(fldListTo(complnd)%flds, 'Faxa_ocphodry', mrg_from1=compatm, mrg_fld1='Faxa_ocphodry', mrg_type1='copy')
          call addmrg(fldListTo(compice)%flds, 'Faxa_ocphodry', mrg_from1=compatm, mrg_fld1='Faxa_ocphodry', mrg_type1='copy')
          call addmrg(fldListTo(compocn)%flds, 'Faxa_ocphodry', mrg_from1=compatm, mrg_fld1='Faxa_ocphodry', &
               mrg_type1='copy_with_weights', mrg_fracname1='ofrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! From atm: hydrophylic organic carbon wet deposition flux
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_ocphiwet')
       call addfld(fldListTo(complnd)%flds, 'Faxa_ocphiwet')
       call addfld(fldListTo(compice)%flds, 'Faxa_ocphiwet')
    else
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_ocphiwet', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_ocphiwet', complnd, mapconsf, 'one', atm2lnd_fmap)
          call addmap(fldListFr(compatm)%flds, 'Faxa_ocphiwet', compice, mapconsf, 'one', atm2ice_fmap)
          call addmap(fldListFr(compatm)%flds, 'Faxa_ocphiwet', compocn, mapconsf, 'one', atm2ocn_fmap)

          call addmrg(fldListTo(complnd)%flds, 'Faxa_ocphiwet', mrg_from1=compatm, mrg_fld1='Faxa_ocphiwet', mrg_type1='copy')
          call addmrg(fldListTo(compice)%flds, 'Faxa_ocphiwet', mrg_from1=compatm, mrg_fld1='Faxa_ocphiwet', mrg_type1='copy')
          call addmrg(fldListTo(compocn)%flds, 'Faxa_ocphiwet', mrg_from1=compatm, mrg_fld1='Faxa_ocphiwet', &
               mrg_type1='copy_with_weights', mrg_fracname1='ofrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to lnd and ice: dust wet deposition flux (size 1)
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_dstwet1')
       call addfld(fldListTo(complnd)%flds, 'Faxa_dstwet1')
       call addfld(fldListTo(compice)%flds, 'Faxa_dstwet1')
    else
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_dstwet1', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_dstwet1', complnd, mapconsf, 'one', atm2lnd_fmap)
          call addmap(fldListFr(compatm)%flds, 'Faxa_dstwet1', compice, mapconsf, 'one', atm2ice_fmap)
          call addmap(fldListFr(compatm)%flds, 'Faxa_dstwet1', compocn, mapconsf, 'one', atm2ocn_fmap)

          call addmrg(fldListTo(complnd)%flds, 'Faxa_dstwet1', mrg_from1=compatm, mrg_fld1='Faxa_dstwet1', mrg_type1='copy')
          call addmrg(fldListTo(compice)%flds, 'Faxa_dstwet1', mrg_from1=compatm, mrg_fld1='Faxa_dstwet1', mrg_type1='copy')
          call addmrg(fldListTo(compocn)%flds, 'Faxa_dstwet1', mrg_from1=compatm, mrg_fld1='Faxa_dstwet1', &
               mrg_type1='copy_with_weights', mrg_fracname1='ofrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to lnd and ice: dust wet deposition flux (size 2)
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_dstwet2')
       call addfld(fldListTo(complnd)%flds, 'Faxa_dstwet2')
       call addfld(fldListTo(compice)%flds, 'Faxa_dstwet2')
    else
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_dstwet2', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_dstwet2', complnd, mapconsf, 'one', atm2lnd_fmap)
          call addmap(fldListFr(compatm)%flds, 'Faxa_dstwet2', compice, mapconsf, 'one', atm2ice_fmap)
          call addmap(fldListFr(compatm)%flds, 'Faxa_dstwet2', compocn, mapconsf, 'one', atm2ocn_fmap)

          call addmrg(fldListTo(complnd)%flds, 'Faxa_dstwet2', mrg_from1=compatm, mrg_fld1='Faxa_dstwet2', mrg_type1='copy')
          call addmrg(fldListTo(compice)%flds, 'Faxa_dstwet2', mrg_from1=compatm, mrg_fld1='Faxa_dstwet2', mrg_type1='copy')
          call addmrg(fldListTo(compocn)%flds, 'Faxa_dstwet2', mrg_from1=compatm, mrg_fld1='Faxa_dstwet2', &
               mrg_type1='copy_with_weights', mrg_fracname1='ofrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to lnd and ice: dust wet deposition flux (size 3)
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_dstwet3')
       call addfld(fldListTo(complnd)%flds, 'Faxa_dstwet3')
       call addfld(fldListTo(compice)%flds, 'Faxa_dstwet3')
    else
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_dstwet3', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_dstwet3', complnd, mapconsf, 'one', atm2lnd_fmap)
          call addmap(fldListFr(compatm)%flds, 'Faxa_dstwet3', compice, mapconsf, 'one', atm2ice_fmap)
          call addmap(fldListFr(compatm)%flds, 'Faxa_dstwet3', compocn, mapconsf, 'one', atm2ocn_fmap)

          call addmrg(fldListTo(complnd)%flds, 'Faxa_dstwet3', mrg_from1=compatm, mrg_fld1='Faxa_dstwet3', mrg_type1='copy')
          call addmrg(fldListTo(compice)%flds, 'Faxa_dstwet3', mrg_from1=compatm, mrg_fld1='Faxa_dstwet3', mrg_type1='copy')
          call addmrg(fldListTo(compocn)%flds, 'Faxa_dstwet3', mrg_from1=compatm, mrg_fld1='Faxa_dstwet3', &
               mrg_type1='copy_with_weights', mrg_fracname1='ofrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to lnd and ice: dust wet deposition flux (size 4)
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_dstwet4')
       call addfld(fldListTo(complnd)%flds, 'Faxa_dstwet4')
       call addfld(fldListTo(compice)%flds, 'Faxa_dstwet4')
    else
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_dstwet4', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_dstwet4', complnd, mapconsf, 'one', atm2lnd_fmap)
          call addmap(fldListFr(compatm)%flds, 'Faxa_dstwet4', compice, mapconsf, 'one', atm2ice_fmap)
          call addmap(fldListFr(compatm)%flds, 'Faxa_dstwet4', compocn, mapconsf, 'one', atm2ocn_fmap)

          call addmrg(fldListTo(complnd)%flds, 'Faxa_dstwet4', mrg_from1=compatm, mrg_fld1='Faxa_dstwet4', mrg_type1='copy')
          call addmrg(fldListTo(compice)%flds, 'Faxa_dstwet4', mrg_from1=compatm, mrg_fld1='Faxa_dstwet4', mrg_type1='copy')
          call addmrg(fldListTo(compocn)%flds, 'Faxa_dstwet4', mrg_from1=compatm, mrg_fld1='Faxa_dstwet4', &
               mrg_type1='copy_with_weights', mrg_fracname1='ofrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to lnd and ice: dust dry deposition flux (size 1)
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_dstdry1')
       call addfld(fldListTo(complnd)%flds, 'Faxa_dstdry1')
       call addfld(fldListTo(compice)%flds, 'Faxa_dstdry1')
    else
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_dstdry1', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_dstdry1', complnd, mapconsf, 'one', atm2lnd_fmap)
          call addmap(fldListFr(compatm)%flds, 'Faxa_dstdry1', compice, mapconsf, 'one', atm2ice_fmap)
          call addmap(fldListFr(compatm)%flds, 'Faxa_dstdry1', compocn, mapconsf, 'one', atm2ocn_fmap)

          call addmrg(fldListTo(complnd)%flds, 'Faxa_dstdry1', mrg_from1=compatm, mrg_fld1='Faxa_dstdry1', mrg_type1='copy')
          call addmrg(fldListTo(compice)%flds, 'Faxa_dstdry1', mrg_from1=compatm, mrg_fld1='Faxa_dstdry1', mrg_type1='copy')
          call addmrg(fldListTo(compocn)%flds, 'Faxa_dstdry1', mrg_from1=compatm, mrg_fld1='Faxa_dstdry1', &
               mrg_type1='copy_with_weights', mrg_fracname1='ofrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to lnd and ice: dust dry deposition flux (size 2)
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_dstdry2')
       call addfld(fldListTo(complnd)%flds, 'Faxa_dstdry2')
       call addfld(fldListTo(compice)%flds, 'Faxa_dstdry2')
    else
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_dstdry2', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_dstdry2', complnd, mapconsf, 'one', atm2lnd_fmap)
          call addmap(fldListFr(compatm)%flds, 'Faxa_dstdry2', compice, mapconsf, 'one', atm2ice_fmap)
          call addmap(fldListFr(compatm)%flds, 'Faxa_dstdry2', compocn, mapconsf, 'one', atm2ocn_fmap)

          call addmrg(fldListTo(complnd)%flds, 'Faxa_dstdry2', mrg_from1=compatm, mrg_fld1='Faxa_dstdry2', mrg_type1='copy')
          call addmrg(fldListTo(compice)%flds, 'Faxa_dstdry2', mrg_from1=compatm, mrg_fld1='Faxa_dstdry2', mrg_type1='copy')
          call addmrg(fldListTo(compocn)%flds, 'Faxa_dstdry2', mrg_from1=compatm, mrg_fld1='Faxa_dstdry2', &
               mrg_type1='copy_with_weights', mrg_fracname1='ofrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to lnd and ice: dust dry deposition flux (size 3)
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_dstdry3')
       call addfld(fldListTo(complnd)%flds, 'Faxa_dstdry3')
       call addfld(fldListTo(compice)%flds, 'Faxa_dstdry3')
    else
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_dstdry3', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_dstdry3', complnd, mapconsf, 'one', atm2lnd_fmap)
          call addmap(fldListFr(compatm)%flds, 'Faxa_dstdry3', compice, mapconsf, 'one', atm2ice_fmap)
          call addmap(fldListFr(compatm)%flds, 'Faxa_dstdry3', compocn, mapconsf, 'one', atm2ocn_fmap)

          call addmrg(fldListTo(complnd)%flds, 'Faxa_dstdry3', mrg_from1=compatm, mrg_fld1='Faxa_dstdry3', mrg_type1='copy')
          call addmrg(fldListTo(compice)%flds, 'Faxa_dstdry3', mrg_from1=compatm, mrg_fld1='Faxa_dstdry3', mrg_type1='copy')
          call addmrg(fldListTo(compocn)%flds, 'Faxa_dstdry3', mrg_from1=compatm, mrg_fld1='Faxa_dstdry3', mrg_type1='copy_with_weights', mrg_fracname1='ofrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to lnd and ice: dust dry deposition flux (size 4)
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_dstdry4')
       call addfld(fldListTo(complnd)%flds, 'Faxa_dstdry4')
       call addfld(fldListTo(compice)%flds, 'Faxa_dstdry4')
    else
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_dstdry4', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_dstdry4', complnd, mapconsf, 'one', atm2lnd_fmap)
          call addmap(fldListFr(compatm)%flds, 'Faxa_dstdry4', compice, mapconsf, 'one', atm2ice_fmap)
          call addmap(fldListFr(compatm)%flds, 'Faxa_dstdry4', compocn, mapconsf, 'one', atm2ocn_fmap)

          call addmrg(fldListTo(complnd)%flds, 'Faxa_dstdry4', mrg_from1=compatm, mrg_fld1='Faxa_dstdry4', mrg_type1='copy')
          call addmrg(fldListTo(compice)%flds, 'Faxa_dstdry4', mrg_from1=compatm, mrg_fld1='Faxa_dstdry4', mrg_type1='copy')
          call addmrg(fldListTo(compocn)%flds, 'Faxa_dstdry4', mrg_from1=compatm, mrg_fld1='Faxa_dstdry4', &
               mrg_type1='copy_with_weights', mrg_fracname1='ofrac')
       end if
    end if

    !=====================================================================
    ! to Atmosphere
    !=====================================================================

    ! ---------------------------------------------------------------------
    ! to atm: Direct albedo (visible radiation)
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(complnd)%flds , 'Sl_avsdr')
       call addfld(fldListFr(compice)%flds , 'Si_avsdr')
       call addfld(fldListMed_ocnalb_o%flds, 'So_avsdr')
       call addfld(fldListTo(compatm)%flds , 'Sx_avsdr')
    else
       if ( fldchk(is_local%wrap%FBexp(compatm), 'Sx_avsdr', rc=rc)) then
          call addmap(fldListFr(complnd)%flds, 'Sl_avsdr', compatm, mapconsf, 'lfrin', lnd2atm_smap)
          call addmap(fldListFr(compice)%flds, 'Si_avsdr', compatm, mapconsf, 'ifrac', ice2atm_smap)
          call addmap(fldListMed_ocnalb_o%flds,'So_avsdr', compatm, mapconsf, 'ofrac', ocn2atm_smap)

          call addmrg(fldListTo(compatm)%flds, 'Sx_avsdr', &
               mrg_from1=complnd, mrg_fld1='Sl_avsdr', mrg_type1='merge', mrg_fracname1='lfrac', &
               mrg_from2=compice, mrg_fld2='Si_avsdr', mrg_type2='merge', mrg_fracname2='ifrac', &
               mrg_from3=compmed, mrg_fld3='So_avsdr', mrg_type3='merge', mrg_fracname3='ofrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to atm: Direct albedo (near-infrared radiation)
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(complnd)%flds , 'Sl_anidr')
       call addfld(fldListFr(compice)%flds , 'Si_anidr')
       call addfld(fldListMed_ocnalb_o%flds, 'So_anidr')
       call addfld(fldListTo(compatm)%flds , 'Sx_anidr')
    else
       if ( fldchk(is_local%wrap%FBexp(compatm), 'Sx_anidr', rc=rc)) then
          call addmap(fldListFr(complnd)%flds, 'Sl_anidr', compatm, mapconsf, 'lfrin', lnd2atm_smap)
          call addmap(fldListFr(compice)%flds, 'Si_anidr', compatm, mapconsf, 'ifrac', ice2atm_smap)
          call addmap(fldListMed_ocnalb_o%flds,'So_anidr', compatm, mapconsf, 'ofrac', ocn2atm_smap)

          call addmrg(fldListTo(compatm)%flds, 'Sx_anidr', &
               mrg_from1=complnd, mrg_fld1='Sl_anidr', mrg_type1='merge', mrg_fracname1='lfrac', &
               mrg_from2=compice, mrg_fld2='Si_anidr', mrg_type2='merge', mrg_fracname2='ifrac', &
               mrg_from3=compmed, mrg_fld3='So_anidr', mrg_type3='merge', mrg_fracname3='ofrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to atm: Diffuse albedo (visible radiation)
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(complnd)%flds , 'Sl_avsdf')
       call addfld(fldListFr(compice)%flds , 'Si_avsdf')
       call addfld(fldListMed_ocnalb_o%flds, 'So_avsdf')
       call addfld(fldListTo(compatm)%flds , 'Sx_avsdf')
    else
       if ( fldchk(is_local%wrap%FBexp(compatm), 'Sx_avsdf', rc=rc)) then
          call addmap(fldListFr(complnd)%flds, 'Sl_avsdf', compatm, mapconsf, 'lfrin', lnd2atm_smap)
          call addmap(fldListFr(compice)%flds, 'Si_avsdf', compatm, mapconsf, 'ifrac', ice2atm_smap)
          call addmap(fldListMed_ocnalb_o%flds,'So_avsdf', compatm, mapconsf, 'ofrac', ocn2atm_smap)

          call addmrg(fldListTo(compatm)%flds, 'Sx_avsdf', &
               mrg_from1=complnd, mrg_fld1='Sl_avsdf', mrg_type1='merge', mrg_fracname1='lfrac', &
               mrg_from2=compice, mrg_fld2='Si_avsdf', mrg_type2='merge', mrg_fracname2='ifrac', &
               mrg_from3=compmed, mrg_fld3='So_avsdf', mrg_type3='merge', mrg_fracname3='ofrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to atm: Diffuse albedo (near-infrared radiation)
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(complnd)%flds , 'Sl_anidf')
       call addfld(fldListFr(compice)%flds , 'Si_anidf')
       call addfld(fldListMed_ocnalb_o%flds, 'So_anidf')
       call addfld(fldListTo(compatm)%flds , 'Sx_anidf')
    else
       if ( fldchk(is_local%wrap%FBexp(compatm), 'Sx_anidf', rc=rc)) then
          call addmap(fldListFr(complnd)%flds, 'Sl_anidf', compatm, mapconsf, 'lfrin', lnd2atm_smap)
          call addmap(fldListFr(compice)%flds, 'Si_anidf', compatm, mapconsf, 'ifrac', ice2atm_smap)
          call addmap(fldListMed_ocnalb_o%flds,'So_anidf', compatm, mapconsf, 'ofrac', ocn2atm_smap)

          call addmrg(fldListTo(compatm)%flds, 'Sx_anidf', &
               mrg_from1=complnd, mrg_fld1='Sl_anidf', mrg_type1='merge', mrg_fracname1='lfrac', &
               mrg_from2=compice, mrg_fld2='Si_anidf', mrg_type2='merge', mrg_fracname2='ifrac', &
               mrg_from3=compmed, mrg_fld3='So_anidf', mrg_type3='merge', mrg_fracname3='ofrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to atm: Reference temperature at 2 meters
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(complnd)%flds , 'Sl_tref')
       call addfld(fldListFr(compice)%flds , 'Si_tref')
       call addfld(fldListMed_aoflux_a%flds, 'So_tref')
       call addfld(fldListMed_aoflux_o%flds, 'So_tref')
       call addfld(fldListTo(compatm)%flds , 'Sx_tref')
    else
       if ( fldchk(is_local%wrap%FBexp(compatm)         , 'Sx_tref', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(complnd,complnd ), 'Sl_tref', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compice,compice ), 'Si_tref', rc=rc) .and. &
            fldchk(is_local%wrap%FBMed_aoflux_o         , 'So_tref', rc=rc)) then
          call addmap(fldListFr(complnd)%flds , 'Sl_tref', compatm, mapconsf, 'lfrin', lnd2atm_fmap)
          call addmap(fldListFr(compice)%flds , 'Si_tref', compatm, mapconsf, 'ifrac', ice2atm_fmap)
          call addmap(fldListMed_aoflux_a%flds, 'So_tref', compocn, mapbilnr, 'one'  , atm2ocn_fmap) ! map atm->ocn
          call addmap(fldListMed_aoflux_o%flds, 'So_tref', compatm, mapconsf, 'ofrac', ocn2atm_fmap) ! map ocn->atm
          call addmrg(fldListTo(compatm)%flds , 'Sx_tref', &
               mrg_from1=complnd, mrg_fld1='Sl_tref', mrg_type1='merge', mrg_fracname1='lfrac', &
               mrg_from2=compice, mrg_fld2='Si_tref', mrg_type2='merge', mrg_fracname2='ifrac', &
               mrg_from3=compmed, mrg_fld3='So_tref', mrg_type3='merge', mrg_fracname3='ofrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to atm: reference specific humidity at 2 meters
    ! ---------------------------------------------------------------------
    do n = 1,size(iso)
       if (phase == 'advertise') then
          call addfld(fldListFr(complnd)%flds , 'Sl_qref'//iso(n))
          call addfld(fldListFr(compice)%flds , 'Si_qref'//iso(n))
          call addfld(fldListMed_aoflux_a%flds, 'So_qref'//iso(n))
          call addfld(fldListMed_aoflux_o%flds, 'So_qref'//iso(n))
          call addfld(fldListTo(compatm)%flds , 'Sx_qref'//iso(n))
       else
          if ( fldchk(is_local%wrap%FBexp(compatm)         , 'Sx_qref'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(complnd,complnd ), 'Sl_qref'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compice,compice ), 'Si_qref'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBMed_aoflux_o         , 'So_qref'//iso(n), rc=rc)) then
             call addmap(fldListFr(complnd)%flds , 'Sl_qref'//iso(n), compatm, mapconsf, 'lfrin', lnd2atm_fmap)
             call addmap(fldListFr(compice)%flds , 'Si_qref'//iso(n), compatm, mapconsf, 'ifrac', ice2atm_fmap)
             call addmap(fldListMed_aoflux_a%flds, 'So_qref'//iso(n), compocn, mapbilnr, 'one'  , atm2ocn_fmap) ! map atm->ocn
             call addmap(fldListMed_aoflux_o%flds, 'So_qref'//iso(n), compatm, mapconsf, 'ofrac', ocn2atm_fmap) ! map ocn->atm
             call addmrg(fldListTo(compatm)%flds , 'Sx_qref'//iso(n), &
                  mrg_from1=complnd, mrg_fld1='Sl_qref'//iso(n), mrg_type1='merge', mrg_fracname1='lfrac', &
                  mrg_from2=compice, mrg_fld2='Si_qref'//iso(n), mrg_type2='merge', mrg_fracname2='ifrac', &
                  mrg_from3=compmed, mrg_fld3='So_qref'//iso(n), mrg_type3='merge', mrg_fracname3='ofrac')
          end if
       end if
    end do

    ! ---------------------------------------------------------------------
    ! to atm: surface fraction velocity in land
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(complnd)%flds, 'Sl_fv')
       call addfld(fldListTo(compatm)%flds, 'Sl_fv')
    else
       if ( fldchk(is_local%wrap%FBexp(compatm)         , 'Sl_fv', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(complnd,complnd ), 'Sl_fv', rc=rc)) then
          call addmap(fldListFr(complnd)%flds, 'Sl_fv', compatm, mapconsf, 'lfrin', lnd2atm_fmap)
          call addmrg(fldListTo(compatm)%flds, 'Sl_fv', mrg_from1=complnd, mrg_fld1='Sl_fv', mrg_type1='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to atm: Aerodynamic resistance
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(complnd)%flds, 'Sl_ram1')
       call addfld(fldListTo(compatm)%flds, 'Sl_ram1')
    else
       if ( fldchk(is_local%wrap%FBexp(compatm)         , 'Sl_ram1', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(complnd,complnd ), 'Sl_ram1', rc=rc)) then
          call addmap(fldListFr(complnd)%flds, 'Sl_ram1', compatm, mapconsf, 'lfrin', lnd2atm_fmap)
          call addmrg(fldListTo(compatm)%flds, 'Sl_ram1', mrg_from1=complnd, mrg_fld1='Sl_ram1', mrg_type1='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to atm: surface snow water equivalent
    ! ---------------------------------------------------------------------
    do n = 1,size(iso)
       if (phase == 'advertise') then
          call addfld(fldListFr(complnd)%flds, trim(fldname))
          call addfld(fldListTo(compatm)%flds, trim(fldname))
       else
          if ( fldchk(is_local%wrap%FBexp(compatm)         , trim(fldname), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(complnd,complnd ), trim(fldname), rc=rc)) then
             call addmap(fldListFr(complnd)%flds, trim(fldname), compatm, mapconsf, 'lfrin', lnd2atm_fmap)
             call addmrg(fldListTo(compatm)%flds, trim(fldname), mrg_from1=complnd, mrg_fld1=trim(fldname), mrg_type1='copy')
          end if
       end if
    end do

    ! ---------------------------------------------------------------------
    ! to atm: Surface snow depth
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compice)%flds, 'Si_snowh')
       call addfld(fldListTo(compatm)%flds, 'Si_snowh')
    else
       if ( fldchk(is_local%wrap%FBexp(compatm)         , 'Si_snowh', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compice,compice ), 'Si_snowh', rc=rc)) then
          call addmap(fldListFr(compice)%flds, 'Si_snowh', compatm, mapconsf, 'ifrac', ice2atm_fmap)
          call addmrg(fldListTo(compatm)%flds, 'Si_snowh', mrg_from1=compice, mrg_fld1='Si_snowh', mrg_type1='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to atm: Mean ice volume per unit area
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compice)%flds, 'Si_vice')
       call addfld(fldListTo(compatm)%flds, 'Si_vice')
    else
       if ( fldchk(is_local%wrap%FBexp(compatm)         , 'Si_vice', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compice,compice ), 'Si_vice', rc=rc)) then
          call addmap(fldListFr(compice)%flds, 'Si_vice', compatm, mapconsf, 'ifrac', ice2atm_fmap)
          call addmrg(fldListTo(compatm)%flds, 'Si_vice', mrg_from1=compice, mrg_fld1='Si_vice', mrg_type1='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to atm: mean snow volume per unit area
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compice)%flds, 'Si_vsno')
       call addfld(fldListTo(compatm)%flds, 'Si_vsno')
    else
       if ( fldchk(is_local%wrap%FBexp(compatm)         , 'Si_vsno', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compice,compice ), 'Si_vsno', rc=rc)) then

          call addmap(fldListFr(compice)%flds, 'Si_vsno', compatm, mapconsf, 'ifrac', ice2atm_fmap)
          call addmrg(fldListTo(compatm)%flds, 'Si_vsno', mrg_from1=compice, mrg_fld1='Si_vsno', mrg_type1='copy')
       endif
    end if

    ! ---------------------------------------------------------------------
    ! to atm: surface saturation specific humidity in ocean
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListMed_aoflux_a%flds, 'So_ssq')
       call addfld(fldListMed_aoflux_o%flds, 'So_ssq')
       call addfld(fldListTo(compatm)%flds , 'So_ssq')
    else
       if (fldchk(is_local%wrap%FBexp(compatm), 'So_ssq', rc=rc) .and. &
           fldchk(is_local%wrap%FBMed_aoflux_o, 'So_ssq', rc=rc)) then
          call addmap(fldListMed_aoflux_a%flds, 'So_ssq', compocn, mapbilnr, 'one'  , atm2ocn_fmap) ! map atm->ocn
          call addmap(fldListMed_aoflux_o%flds, 'So_ssq', compatm, mapconsf, 'ofrac', ocn2atm_fmap) ! map ocn->atm
          call addmrg(fldListTo(compatm)%flds , 'So_ssq', mrg_from1=compmed, mrg_fld1='So_ssq', mrg_type1='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to atm: square of exch. coeff (tracers)
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListMed_aoflux_a%flds, 'So_re')
       call addfld(fldListMed_aoflux_o%flds, 'So_re')
       call addfld(fldListTo(compatm)%flds , 'So_re')
    else
       if (fldchk(is_local%wrap%FBexp(compatm), 'So_re', rc=rc) .and. &
           fldchk(is_local%wrap%FBMed_aoflux_o, 'So_re', rc=rc)) then
          call addmap(fldListMed_aoflux_a%flds, 'So_re', compocn, mapbilnr, 'one'  , atm2ocn_fmap) ! map atm->ocn
          call addmap(fldListMed_aoflux_o%flds, 'So_re', compatm, mapconsf, 'ofrac', ocn2atm_fmap) ! map ocn->atm
          call addmrg(fldListTo(compatm)%flds , 'So_re', mrg_from1=compmed, mrg_fld1='So_re', mrg_type1='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to atm: 10m wind
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(complnd)%flds , 'Sl_u10')
       call addfld(fldListFr(compice)%flds , 'Si_u10')
       call addfld(fldListMed_aoflux_a%flds, 'So_u10')
       call addfld(fldListMed_aoflux_o%flds, 'So_u10')
       call addfld(fldListTo(compatm)%flds , 'Sx_u10')
    else
       if (fldchk(is_local%wrap%FBexp(compatm), 'Sx_u10', rc=rc)) then
          call addmap(fldListMed_aoflux_a%flds, 'So_u10', compocn, mapbilnr, 'one'  , atm2ocn_fmap) ! map atm->ocn
          call addmap(fldListMed_aoflux_o%flds, 'So_u10', compatm, mapconsf, 'ofrac', ocn2atm_fmap) ! map ocn->atm
          call addmap(fldListFr(complnd)%flds , 'Sl_u10', compatm, mapconsf, 'lfrin', lnd2atm_fmap)
          call addmap(fldListFr(compice)%flds , 'Si_u10', compatm, mapconsf, 'ifrac', ice2atm_fmap)
          call addmrg(fldListTo(compatm)%flds , 'Sx_u10', &
               mrg_from1=complnd, mrg_fld1='Sl_u10', mrg_type1='merge', mrg_fracname1='lfrac', &
               mrg_from2=compice, mrg_fld2='Si_u10', mrg_type2='merge', mrg_fracname2='ifrac', &
               mrg_from3=compmed, mrg_fld3='So_u10', mrg_type3='merge', mrg_fracname3='ofrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to atm: surface fraction velocity from med
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListMed_aoflux_o%flds, 'So_ustar')
       call addfld(fldListTo(compatm)%flds , 'So_ustar')
    else
       if (fldchk(is_local%wrap%FBMed_aoflux_o, 'So_ustar', rc=rc) .and. &
           fldchk(is_local%wrap%FBExp(compocn), 'So_ustar', rc=rc)) then
          call addmap(fldListMed_aoflux_o%flds, 'So_ustar', compatm, mapconsf, 'ofrac', ocn2atm_fmap) ! map ocn->atm
          call addmrg(fldListTo(compocn)%flds , 'So_ustar', &
               mrg_from1=compmed, mrg_fld1='So_ustar', mrg_type1='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to atm, ice and wav: surface temperature
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(complnd)%flds, 'Sl_t')
       call addfld(fldListFr(compice)%flds, 'Si_t')
       call addfld(fldListFr(compocn)%flds, 'So_t')
       call addfld(fldListTo(compice)%flds, 'So_t')
       call addfld(fldListTo(compwav)%flds, 'So_t')
       call addfld(fldListTo(compatm)%flds, 'Si_t') ! nems-frac
       call addfld(fldListTo(compatm)%flds, 'So_t') ! cesm, nems-frac
       call addfld(fldListTo(compatm)%flds, 'Sx_t') ! cesm, nems-orig
    else
       if (fldchk(is_local%wrap%FBexp(compatm)        , 'Sx_t', rc=rc) .and. &
           fldchk(is_local%wrap%FBImp(complnd,complnd), 'Sl_t', rc=rc) .and. &
           fldchk(is_local%wrap%FBImp(compice,compice), 'Si_t', rc=rc) .and. &
           fldchk(is_local%wrap%FBImp(compocn,compocn), 'So_t', rc=rc)) then

          ! CESM
          call addmap(fldListFr(complnd)%flds, 'Sl_t', compatm, mapconsf , 'lfrin', lnd2atm_fmap)
          call addmap(fldListFr(compice)%flds, 'Si_t', compatm, mapconsf , 'ifrac', ice2atm_fmap)
          call addmap(fldListFr(compocn)%flds, 'So_t', compatm, mapconsf , 'ofrac', ocn2atm_fmap)
          call addmrg(fldListTo(compatm)%flds, 'Sx_t', &
               mrg_from1=complnd, mrg_fld1='Sl_t', mrg_type1='merge', mrg_fracname1='lfrac', &
               mrg_from2=compice, mrg_fld2='Si_t', mrg_type2='merge', mrg_fracname2='ifrac', &
               mrg_from3=compocn, mrg_fld3='So_t', mrg_type3='merge', mrg_fracname3='ofrac')

       else if (fldchk(is_local%wrap%FBexp(compatm)        , 'Sx_t', rc=rc) .and. &
                fldchk(is_local%wrap%FBImp(compocn,compocn), 'So_t', rc=rc) .and. &
                fldchk(is_local%wrap%FBImp(compice,compice), 'Si_t', rc=rc)) then

          ! NEMS-orig
          call addmap(fldListFr(compice)%flds, 'Si_t', compatm, mapconsf , 'ifrac', ice2atm_fmap)
          call addmap(fldListFr(compocn)%flds, 'So_t', compatm, mapconsf , 'ofrac', ocn2atm_fmap)
          call addmrg(fldListTo(compatm)%flds, 'Sx_t', &
               mrg_from1=compice, mrg_fld1='Si_t', mrg_type1='merge', mrg_fracname1='ifrac', &
               mrg_from2=compocn, mrg_fld2='So_t', mrg_type2='merge', mrg_fracname2='ofrac')
       end if

       if ( fldchk(is_local%wrap%FBexp(compatm)        , 'Si_t', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compice,compice), 'Si_t', rc=rc)) then
            ! NEMS-frac
          call addmap(fldListFr(compice)%flds, 'Si_t', compatm, mapconsf , 'ifrac', ice2atm_fmap)
          call addmrg(fldListTo(compatm)%flds, 'Si_t', mrg_from1=compice, mrg_fld1='Si_t', mrg_type1='copy')
       end if
       if ( fldchk(is_local%wrap%FBexp(compatm)        , 'So_t', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compocn,compocn), 'So_t', rc=rc)) then
          ! NEMS-frac and CESM
          call addmap(fldListFr(compocn)%flds, 'So_t', compatm, mapconsf , 'ofrac', ocn2atm_fmap)
          call addmrg(fldListTo(compatm)%flds, 'So_t', mrg_from1=compocn, mrg_fld1='So_t', mrg_type1='copy')
       end if
       if ( fldchk(is_local%wrap%FBexp(compice)        , 'So_t', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compocn,compocn), 'So_t', rc=rc)) then
          call addmap(fldListFr(compocn)%flds, 'So_t', compice, mapfcopy , 'unset', 'unset')
          call addmrg(fldListTo(compice)%flds, 'So_t', mrg_from1=compocn, mrg_fld1='So_t', mrg_type1='copy')
       end if
       if ( fldchk(is_local%wrap%FBexp(compwav), 'So_t', rc=rc)) then
          call addmap(fldListFr(compocn)%flds,' So_t', compwav, mapbilnr , 'one'  , ocn2wav_smap)
          call addmrg(fldListTo(compwav)%flds, 'So_t', mrg_from1=compocn, mrg_fld1='So_t', mrg_type1='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to atm: zonal surface stress
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(complnd)%flds , 'Fall_taux')
       call addfld(fldListFr(compice)%flds , 'Faii_taux')
       call addfld(fldListMed_aoflux_a%flds, 'Faox_taux')
       call addfld(fldListMed_aoflux_o%flds, 'Faox_taux')
       call addfld(fldListTo(compatm)%flds , 'Faii_taux') ! nems-frac
       call addfld(fldListTo(compatm)%flds , 'Faxx_taux') ! cesm, nems-orig
    else
       if (fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_taux', rc=rc) .and. &
           fldchk(is_local%wrap%FBMed_aoflux_a, 'Faox_taux', rc=rc)) then
          call addmap(fldListMed_aoflux_a%flds, 'Faox_taux', compocn, mapconsf, 'one'  , atm2ocn_fmap) ! atm->ocn
          call addmap(fldListMed_aoflux_o%flds, 'Faox_taux', compatm, mapconsf, 'ofrac', ocn2atm_fmap) ! ocn->atm
       end if
       if (fldchk(is_local%wrap%FBImp(complnd,complnd), 'Fall_taux', rc=rc)) then
          call addmap(fldListFr(complnd)%flds , 'Fall_taux', compatm, mapconsf, 'lfrin', lnd2atm_fmap)
       end if
       if ( fldchk(is_local%wrap%FBImp(complnd,complnd), 'Faii_taux', rc=rc)) then
          call addmap(fldListFr(compice)%flds , 'Faii_taux', compatm, mapconsf, 'ifrac', ice2atm_fmap)
       end if

       if ( fldchk(is_local%wrap%FBImp(complnd,complnd), 'Fall_taux', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compice,compice), 'Faii_taux', rc=rc) .and. &
            fldchk(is_local%wrap%FBMed_aoflux_o        , 'Faox_taux', rc=rc) .and. &
            fldchk(is_local%wrap%FBexp(compatm)        , 'Faxx_taux', rc=rc)) then

          ! CESM
          call addmrg(fldListTo(compatm)%flds , 'Faxx_taux', &
               mrg_from1=complnd, mrg_fld1='Fall_taux', mrg_type1='merge', mrg_fracname1='lfrac', &
               mrg_from2=compice, mrg_fld2='Faii_taux', mrg_type2='merge', mrg_fracname2='ifrac', &
               mrg_from3=compmed, mrg_fld3='Faox_taux', mrg_type3='merge', mrg_fracname3='ofrac')

       else if ( fldchk(is_local%wrap%FBImp(compice,compice), 'Faii_taux', rc=rc) .and. &
                 fldchk(is_local%wrap%FBMed_aoflux_o        , 'Faox_taux', rc=rc) .and. &
                 fldchk(is_local%wrap%FBexp(compatm)        , 'Faxx_taux', rc=rc)) then

          ! NEMS orig (here ofrac = 1.-ifrac)
          call addmrg(fldListTo(compatm)%flds , 'Faxx_taux', &
               mrg_from1=compice, mrg_fld1='Faii_taux', mrg_type1='merge', mrg_fracname1='ifrac', &
               mrg_from2=compmed, mrg_fld2='Faox_taux', mrg_type2='merge', mrg_fracname2='ofrac')

       else if ( fldchk(is_local%wrap%FBImp(compice,compice), 'Faii_taux', rc=rc) .and. &
                 fldchk(is_local%wrap%FBexp(compatm)        , 'Faii_taux', rc=rc)) then

          ! NEMS frac
          call addmrg(fldListTo(compatm)%flds, 'Faii_taux', &
               mrg_from1=compice, mrg_fld1='Faii_taux', mrg_type1='merge', mrg_fracname1='ifrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to atm: meridional surface stress'
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(complnd)%flds , 'Fall_tauy')
       call addfld(fldListFr(compice)%flds , 'Faii_tauy')
       call addfld(fldListMed_aoflux_a%flds, 'Faox_tauy')
       call addfld(fldListMed_aoflux_o%flds, 'Faox_tauy')
       call addfld(fldListTo(compatm)%flds , 'Faii_tauy') ! nems-frac
       call addfld(fldListTo(compatm)%flds , 'Faxx_tauy') ! cesm, nems-orig
    else
       if (fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_tauy', rc=rc) .and. &
           fldchk(is_local%wrap%FBMed_aoflux_a, 'Faox_tauy', rc=rc)) then
          call addmap(fldListMed_aoflux_a%flds, 'Faox_tauy', compocn, mapconsf, 'one'  , atm2ocn_fmap) ! map atm->ocn
          call addmap(fldListMed_aoflux_o%flds, 'Faox_tauy', compatm, mapconsf, 'ofrac', ocn2atm_fmap) ! map ocn->atm
       end if
       if ( fldchk(is_local%wrap%FBImp(complnd,complnd), 'Fall_tauy', rc=rc)) then
          call addmap(fldListFr(complnd)%flds , 'Fall_tauy', compatm, mapconsf, 'lfrin', lnd2atm_fmap)
       end if
       if ( fldchk(is_local%wrap%FBImp(complnd,complnd), 'Faii_tauy', rc=rc)) then
          call addmap(fldListFr(compice)%flds , 'Faii_tauy', compatm, mapconsf, 'ifrac', ice2atm_fmap)
       end if

       if (fldchk(is_local%wrap%FBImp(complnd,complnd), 'Fall_tauy', rc=rc) .and. &
           fldchk(is_local%wrap%FBImp(compice,compice), 'Faii_tauy', rc=rc) .and. &
           fldchk(is_local%wrap%FBMed_aoflux_o        , 'Faox_tauy', rc=rc) .and. &
           fldchk(is_local%wrap%FBexp(compatm)        , 'Faxx_tauy', rc=rc)) then

          ! CESM
          call addmrg(fldListTo(compatm)%flds , 'Faxx_tauy', &
               mrg_from1=complnd, mrg_fld1='Fall_tauy', mrg_type1='merge', mrg_fracname1='lfrac', &
               mrg_from2=compice, mrg_fld2='Faii_tauy', mrg_type2='merge', mrg_fracname2='ifrac', &
               mrg_from3=compmed, mrg_fld3='Faox_tauy', mrg_type3='merge', mrg_fracname3='ofrac')

       else if ( fldchk(is_local%wrap%FBImp(compice,compice), 'Faii_tauy', rc=rc) .and. &
                 fldchk(is_local%wrap%FBMed_aoflux_o        , 'Faox_tauy', rc=rc) .and. &
                 fldchk(is_local%wrap%FBexp(compatm)        , 'Faxx_tauy', rc=rc)) then

          ! NEMS orig (here ofrac = 1.-ifrac)
          call addmrg(fldListTo(compatm)%flds , 'Faxx_tauy', &
               mrg_from1=compice, mrg_fld1='Faii_tauy', mrg_type1='merge', mrg_fracname1='ifrac', &
               mrg_from2=compmed, mrg_fld2='Faox_tauy', mrg_type2='merge', mrg_fracname2='ofrac')

       else if ( fldchk(is_local%wrap%FBImp(compice,compice), 'Faii_tauy', rc=rc) .and. &
                 fldchk(is_local%wrap%FBexp(compatm)        , 'Faii_tauy', rc=rc)) then

          ! NEMS frac
          call addmrg(fldListTo(compatm)%flds, 'Faii_tauy', &
               mrg_from1=compice, mrg_fld1='Faii_tauy', mrg_type1='merge', mrg_fracname1='ifrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to atm: latent heat flux
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(complnd)%flds , 'Fall_lat')
       call addfld(fldListFr(compice)%flds , 'Faii_lat')
       call addfld(fldListMed_aoflux_a%flds, 'Faox_lat')
       call addfld(fldListMed_aoflux_o%flds, 'Faox_lat')
       call addfld(fldListTo(compatm)%flds , 'Faii_lat') ! nems-frac
       call addfld(fldListTo(compatm)%flds , 'Faxx_lat') ! cesm, nems-orig
    else
       if (fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_lat', rc=rc) .and. &
           fldchk(is_local%wrap%FBMed_aoflux_a, 'Faox_lat', rc=rc)) then
          call addmap(fldListMed_aoflux_a%flds,'Faox_lat', compocn, mapconsf, 'one'  , atm2ocn_fmap) ! map atm->ocn
          call addmap(fldListMed_aoflux_o%flds,'Faox_lat', compatm, mapconsf, 'ofrac', ocn2atm_fmap) ! map ocn->atm
       end if
       if ( fldchk(is_local%wrap%FBImp(complnd,complnd), 'Fall_lat', rc=rc)) then
          call addmap(fldListFr(complnd)%flds, 'Fall_lat', compatm, mapconsf, 'lfrin', lnd2atm_fmap)
       end if
       if ( fldchk(is_local%wrap%FBImp(complnd,complnd), 'Faii_lat', rc=rc)) then
          call addmap(fldListFr(compice)%flds , 'Faii_lat', compatm, mapconsf, 'ifrac', ice2atm_fmap)
       end if

       if ( fldchk(is_local%wrap%FBImp(complnd,complnd), 'Fall_lat', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compice,compice), 'Faii_lat', rc=rc) .and. &
            fldchk(is_local%wrap%FBMed_aoflux_o        , 'Faox_lat', rc=rc) .and. &
            fldchk(is_local%wrap%FBexp(compatm)        , 'Faxx_lat', rc=rc)) then

          ! CESM
          call addmrg(fldListTo(compatm)%flds , 'Faxx_lat', &
               mrg_from1=complnd, mrg_fld1='Fall_lat', mrg_type1='merge', mrg_fracname1='lfrac', &
               mrg_from2=compice, mrg_fld2='Faii_lat', mrg_type2='merge', mrg_fracname2='ifrac', &
               mrg_from3=compmed, mrg_fld3='Faox_lat', mrg_type3='merge', mrg_fracname3='ofrac')

       else if ( fldchk(is_local%wrap%FBImp(compice,compice), 'Faii_lat', rc=rc) .and. &
                 fldchk(is_local%wrap%FBMed_aoflux_o        , 'Faox_lat', rc=rc) .and. &
                 fldchk(is_local%wrap%FBexp(compatm)        , 'Faxx_lat', rc=rc)) then

          ! NEMS orig (here ofrac = 1.-ifrac)
          call addmrg(fldListTo(compatm)%flds , 'Faxx_lat', &
               mrg_from1=compice, mrg_fld1='Faii_lat', mrg_type1='merge', mrg_fracname1='ifrac', &
               mrg_from2=compmed, mrg_fld2='Faox_lat', mrg_type2='merge', mrg_fracname2='ofrac')

       else if ( fldchk(is_local%wrap%FBImp(compice,compice), 'Faii_lat', rc=rc) .and. &
                 fldchk(is_local%wrap%FBexp(compatm)        , 'Faii_lat', rc=rc)) then
          ! NEMS frac
          call addmap(fldListFr(compice)%flds, 'Faii_lat', compatm, mapconsf, 'one', ice2atm_fmap)
          call addmrg(fldListTo(compatm)%flds, 'Faii_lat', &
               mrg_from1=compice, mrg_fld1='Faii_lat', mrg_type1='merge', mrg_fracname1='ifrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to atm: surface sensible heat flux
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(complnd)%flds , 'Fall_sen')
       call addfld(fldListFr(compice)%flds , 'Faii_sen')
       call addfld(fldListMed_aoflux_a%flds, 'Faox_sen')
       call addfld(fldListMed_aoflux_o%flds, 'Faox_sen')
       call addfld(fldListTo(compatm)%flds , 'Faii_sen') ! nems-frac
       call addfld(fldListTo(compatm)%flds , 'Faxx_sen') ! cesm, nems-orig
    else

       if ( fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_sen', rc=rc) .and. &
            fldchk(is_local%wrap%FBMed_aoflux_a, 'Faox_sen', rc=rc)) then
          call addmap(fldListMed_aoflux_a%flds,'Faox_sen', compocn, mapconsf, 'one'  , atm2ocn_fmap) ! map atm->ocn
          call addmap(fldListMed_aoflux_o%flds,'Faox_sen', compatm, mapconsf, 'ofrac', ocn2atm_fmap) ! map ocn->atm
       end if
       if ( fldchk(is_local%wrap%FBImp(complnd,complnd), 'Fall_sen', rc=rc)) then
          call addmap(fldListFr(complnd)%flds, 'Fall_sen', compatm, mapconsf, 'lfrin', lnd2atm_fmap)
       end if
       if ( fldchk(is_local%wrap%FBImp(complnd,complnd), 'Faii_sen', rc=rc)) then
          call addmap(fldListFr(compice)%flds , 'Faii_sen', compatm, mapconsf, 'ifrac', ice2atm_fmap)
       end if

       if ( fldchk(is_local%wrap%FBImp(complnd,complnd), 'Fall_sen', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compice,compice), 'Faii_sen', rc=rc) .and. &
            fldchk(is_local%wrap%FBMed_aoflux_o        , 'Faox_sen', rc=rc) .and. &
            fldchk(is_local%wrap%FBexp(compatm)        , 'Faxx_sen', rc=rc)) then

          ! CESM
          call addmrg(fldListTo(compatm)%flds , 'Faxx_sen', &
               mrg_from1=complnd, mrg_fld1='Fall_sen', mrg_type1='merge', mrg_fracname1='lfrac', &
               mrg_from2=compice, mrg_fld2='Faii_sen', mrg_type2='merge', mrg_fracname2='ifrac', &
               mrg_from3=compmed, mrg_fld3='Faox_sen', mrg_type3='merge', mrg_fracname3='ofrac')

       else if ( fldchk(is_local%wrap%FBImp(compice,compice), 'Faii_sen', rc=rc) .and. &
                 fldchk(is_local%wrap%FBMed_aoflux_o        , 'Faox_sen', rc=rc) .and. &
                 fldchk(is_local%wrap%FBexp(compatm)        , 'Faxx_sen', rc=rc)) then

          ! NEMS orig (here ofrac = 1.-ifrac)
          call addmrg(fldListTo(compatm)%flds , 'Faxx_sen', &
               mrg_from1=compice, mrg_fld1='Faii_sen', mrg_type1='merge', mrg_fracname1='ifrac', &
               mrg_from2=compmed, mrg_fld2='Faox_sen', mrg_type2='merge', mrg_fracname2='ofrac')

       else if ( fldchk(is_local%wrap%FBImp(compice,compice), 'Faii_sen', rc=rc) .and. &
                 fldchk(is_local%wrap%FBexp(compatm)        , 'Faii_sen', rc=rc)) then
          ! NEMS frac
          call addmap(fldListFr(compice)%flds, 'Faii_sen', compatm, mapconsf, 'one', ice2atm_fmap)
          call addmrg(fldListTo(compatm)%flds, 'Faii_sen', &
               mrg_from1=compice, mrg_fld1='Faii_sen', mrg_type1='merge', mrg_fracname1='ifrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to atm: surface upward longwave heat flux
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(complnd)%flds , 'Fall_lwup')
       call addfld(fldListFr(compice)%flds , 'Faii_lwup')
       call addfld(fldListMed_aoflux_a%flds, 'Faox_lwup')
       call addfld(fldListMed_aoflux_o%flds, 'Faox_lwup')
       call addfld(fldListTo(compatm)%flds , 'Faii_lwup') ! nems-frac
       call addfld(fldListTo(compatm)%flds , 'Faxx_lwup') ! cesm, nems-orig
    else
       if ( fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_lwup', rc=rc) .and. &
            fldchk(is_local%wrap%FBMed_aoflux_a, 'Faox_lwup', rc=rc)) then
          call addmap(fldListMed_aoflux_a%flds,'Faox_lwup', compocn, mapconsf, 'one'  , atm2ocn_fmap) ! map atm->ocn
          call addmap(fldListMed_aoflux_o%flds,'Faox_lwup', compatm, mapconsf, 'ofrac', ocn2atm_fmap) ! map ocn->atm
       end if
       if ( fldchk(is_local%wrap%FBImp(complnd,complnd), 'Fall_lwup', rc=rc)) then
          call addmap(fldListFr(complnd)%flds, 'Fall_lwup', compatm, mapconsf, 'lfrin', lnd2atm_fmap)
       end if
       if ( fldchk(is_local%wrap%FBImp(complnd,complnd), 'Faii_lwup', rc=rc)) then
          call addmap(fldListFr(compice)%flds , 'Faii_lwup', compatm, mapconsf, 'ifrac', ice2atm_fmap)
       end if

       if ( fldchk(is_local%wrap%FBMed_aoflux_o        , 'Faox_lwup', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compice,compice), 'Faii_lwup', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(complnd,complnd), 'Fall_lwup', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compatm)        , 'Faxx_lwup', rc=rc)) then

          ! CESM
          call addmrg(fldListTo(compatm)%flds, 'Faxx_lwup', &
               mrg_from1=complnd, mrg_fld1='Fall_lwup', mrg_type1='merge', mrg_fracname1='lfrac', &
               mrg_from2=compice, mrg_fld2='Faii_lwup', mrg_type2='merge', mrg_fracname2='ifrac', &
               mrg_from3=compmed, mrg_fld3='Faox_lwup', mrg_type3='merge', mrg_fracname3='ofrac')

       else if ( fldchk(is_local%wrap%FBMed_aoflux_o        , 'Faox_lwup', rc=rc) .and. &
                 fldchk(is_local%wrap%FBImp(compice,compice), 'Faii_lwup', rc=rc) .and. &
                 fldchk(is_local%wrap%FBExp(compatm)        , 'Faxx_lwup', rc=rc)) then

          ! NEMS-orig (Note here ofrac = 1.-ifrac)
          call addmrg(fldListTo(compatm)%flds, 'Faxx_lwup', &
               mrg_from1=compice, mrg_fld1='Faii_lwup', mrg_type1='merge', mrg_fracname1='ifrac', &
               mrg_from2=compmed, mrg_fld2='Faox_lwup', mrg_type2='merge', mrg_fracname2='ofrac')

       else if ( fldchk(is_local%wrap%FBImp(compice,compice), 'Faii_lwup', rc=rc) .and. &
                 fldchk(is_local%wrap%FBexp(compatm)        , 'Faii_lwup', rc=rc)) then

          ! NEMS frac (merge done in fv3)
          call addmrg(fldListTo(compatm)%flds, 'Faii_lwup', &
               mrg_from1=compice, mrg_fld1='Faii_lwup', mrg_type1='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to atm: evaporation water flux
    ! ---------------------------------------------------------------------
    do n = 1,size(iso)
       if (phase == 'advertise') then
          call addfld(fldListFr(complnd)%flds , 'Fall_evap'//iso(n))
          call addfld(fldListFr(compice)%flds , 'Faii_evap'//iso(n))
          call addfld(fldListMed_aoflux_a%flds, 'Faox_evap'//iso(n))
          call addfld(fldListMed_aoflux_o%flds, 'Faox_evap'//iso(n))
          call addfld(fldListTo(compatm)%flds , 'Faii_evap'//iso(n)) ! nems-frac
          call addfld(fldListTo(compatm)%flds , 'Faxx_evap'//iso(n)) ! cesm, nems-orig
       else

          if ( fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_evap'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBMed_aoflux_a, 'Faox_evap'//iso(n), rc=rc)) then
             call addmap(fldListMed_aoflux_a%flds,'Faox_evap'//iso(n), compocn, mapconsf, 'one'  , atm2ocn_fmap) ! map atm->ocn
             call addmap(fldListMed_aoflux_o%flds,'Faox_evap'//iso(n), compatm, mapconsf, 'ofrac', ocn2atm_fmap) ! map ocn->atm
          end if
          if ( fldchk(is_local%wrap%FBImp(complnd,complnd), 'Fall_evap'//iso(n), rc=rc)) then
             call addmap(fldListFr(complnd)%flds, 'Fall_evap'//iso(n), compatm, mapconsf, 'lfrin', lnd2atm_fmap)
          end if
          if ( fldchk(is_local%wrap%FBImp(complnd,complnd), 'Faii_evap'//iso(n), rc=rc)) then
             call addmap(fldListFr(compice)%flds , 'Faii_evap'//iso(n), compatm, mapconsf, 'ifrac', ice2atm_fmap)
          end if

          if ( fldchk(is_local%wrap%FBImp(complnd,complnd), 'Fall_evap'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compice,compice), 'Faii_evap'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBMed_aoflux_o        , 'Faox_evap'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBexp(compatm)        , 'Faxx_evap'//iso(n), rc=rc)) then

             ! CESM
             call addmrg(fldListTo(compatm)%flds , 'Faxx_evap'//iso(n), &
                  mrg_from1=complnd, mrg_fld1='Fall_evap'//iso(n), mrg_type1='merge', mrg_fracname1='lfrac', &
                  mrg_from2=compice, mrg_fld2='Faii_evap'//iso(n), mrg_type2='merge', mrg_fracname2='ifrac', &
                  mrg_from3=compmed, mrg_fld3='Faox_evap'//iso(n), mrg_type3='merge', mrg_fracname3='ofrac')

          else if ( fldchk(is_local%wrap%FBImp(compice,compice), 'Faii_evap'//iso(n), rc=rc) .and. &
                    fldchk(is_local%wrap%FBMed_aoflux_o        , 'Faox_evap'//iso(n), rc=rc) .and. &
                    fldchk(is_local%wrap%FBexp(compatm)        , 'Faxx_evap'//iso(n), rc=rc)) then

             ! NEMS orig (here ofrac = 1.-ifrac)
             call addmrg(fldListTo(compatm)%flds , 'Faxx_evap'//iso(n), &
                  mrg_from1=compice, mrg_fld1='Faii_evap'//iso(n), mrg_type1='merge', mrg_fracname1='ifrac', &
                  mrg_from2=compmed, mrg_fld2='Faox_evap'//iso(n), mrg_type2='merge', mrg_fracname2='ofrac')

          else if ( fldchk(is_local%wrap%FBImp(compice,compice), 'Faii_evap'//iso(n), rc=rc) .and. &
                    fldchk(is_local%wrap%FBexp(compatm)        , 'Faii_evap'//iso(n), rc=rc)) then
             ! NEMS frac
             call addmap(fldListFr(compice)%flds, 'Faii_evap'//iso(n), compatm, mapconsf, 'one', ice2atm_fmap)
             call addmrg(fldListTo(compatm)%flds, 'Faii_evap'//iso(n), &
                  mrg_from1=compice, mrg_fld1='Faii_evap'//iso(n), mrg_type1='merge', mrg_fracname1='ifrac')
          end if
       end if
    end do

    ! ---------------------------------------------------------------------
    ! to atm: dust fluxes from land
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(complnd)%flds, 'Fall_flxdst1')
       call addfld(fldListFr(complnd)%flds, 'Fall_flxdst2')
       call addfld(fldListFr(complnd)%flds, 'Fall_flxdst3')
       call addfld(fldListFr(complnd)%flds, 'Fall_flxdst4')
       call addfld(fldListTo(compatm)%flds, 'Fall_flxdst1')
       call addfld(fldListTo(compatm)%flds, 'Fall_flxdst2')
       call addfld(fldListTo(compatm)%flds, 'Fall_flxdst3')
       call addfld(fldListTo(compatm)%flds, 'Fall_flxdst4')
    else
       if ( fldchk(is_local%wrap%FBImp(complnd, complnd), 'Fall_dstflx1', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compatm)         , 'Fall_dstflx1', rc=rc)) then

          ! Assume if the first bin is available from land then all bins are available
          call addmap(fldListFr(complnd)%flds, 'Fall_flxdst1', compatm, mapconsf, 'lfrin', lnd2atm_fmap)
          call addmap(fldListFr(complnd)%flds, 'Fall_flxdst2', compatm, mapconsf, 'lfrin', lnd2atm_fmap)
          call addmap(fldListFr(complnd)%flds, 'Fall_flxdst3', compatm, mapconsf, 'lfrin', lnd2atm_fmap)
          call addmap(fldListFr(complnd)%flds, 'Fall_flxdst4', compatm, mapconsf, 'lfrin', lnd2atm_fmap)

          call addmrg(fldListTo(compatm)%flds, 'Fall_flxdst1', &
               mrg_from1=complnd, mrg_fld1='Fall_flxdst1', mrg_type1='copy_with_weights', mrg_fracname1='lfrac')
          call addmrg(fldListTo(compatm)%flds, 'Fall_flxdst2', &
               mrg_from1=complnd, mrg_fld1='Fall_flxdst2', mrg_type1='copy_with_weights', mrg_fracname1='lfrac')
          call addmrg(fldListTo(compatm)%flds, 'Fall_flxdst3', &
               mrg_from1=complnd, mrg_fld1='Fall_flxdst3', mrg_type1='copy_with_weights', mrg_fracname1='lfrac')
          call addmrg(fldListTo(compatm)%flds, 'Fall_flxdst4', &
               mrg_from1=complnd, mrg_fld1='Fall_flxdst4', mrg_type1='copy_with_weights', mrg_fracname1='lfrac')
       end if
    end if

    !-----------------------------------------------------------------------------
    ! to atm: MEGAN emissions fluxes from land
    !-----------------------------------------------------------------------------
    if (phase == 'advertise') then
       do num = 1, max_megan
          write(cnum,'(i3.3)') num
          fldname = 'Fall_voc' // cnum
          call addfld(fldListFr(complnd)%flds, trim(fldname))
          call addfld(fldListTo(compatm)%flds, trim(fldname))
       end do
    else
       do num = 1, max_megan
          write(cnum,'(i3.3)') num
          fldname = 'Fall_voc' // cnum
          if ( fldchk(is_local%wrap%FBImp(complnd, complnd), trim(fldname), rc=rc) .and. &
               fldchk(is_local%wrap%FBExp(compatm)         , trim(fldname), rc=rc)) then

             call addmap(fldListFr(complnd)%flds, trim(fldname), compatm, mapconsf, 'one', atm2lnd_fmap)
             call addmrg(fldListTo(compatm)%flds, trim(fldname), &
                  mrg_from1=complnd, mrg_fld1=trim(fldname), mrg_type1='merge', mrg_fracname1='lfrac')
          end if
       end do
    end if

    !-----------------------------------------------------------------------------
    ! to atm: fire emissions fluxes from land
    !-----------------------------------------------------------------------------

    ! 'wild fire emission fluxes'
    if (phase == 'advertise') then
       do num = 1, max_fire
          write(cnum,'(i2.2)') num
          fldname  = 'Fall_fire' // cnum
          call addfld(fldListFr(complnd)%flds, trim(fldname))
          call addfld(fldListTo(compatm)%flds, trim(fldname))
       end do
    else
       do num = 1, max_fire
          write(cnum,'(i2.2)') num
          fldname  = 'Fall_fire' // cnum
          if ( fldchk(is_local%wrap%FBImp(complnd, complnd), trim(fldname), rc=rc) .and. &
               fldchk(is_local%wrap%FBExp(compatm)         , trim(fldname), rc=rc)) then
             call addmap(fldListFr(complnd)%flds, trim(fldname), compatm, mapconsf, 'one', lnd2atm_fmap)
             call addmrg(fldListTo(compatm)%flds, trim(fldname), &
                  mrg_from1=complnd, mrg_fld1=trim(fldname), mrg_type1='merge', mrg_fracname1='lfrac')
          end if
       end do
    end if

    ! 'wild fire plume height'
    if (phase == 'advertise') then
       call addfld(fldListFr(complnd)%flds, 'Sl_fztop')
       call addfld(fldListTo(compatm)%flds, 'Sl_fztop')
    else
       if ( fldchk(is_local%wrap%FBImp(complnd, complnd), 'Sl_fztop', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compatm)         , 'Sl_fztop', rc=rc)) then
          call addmap(fldListFr(complnd)%flds, 'Sl_fztop', compatm, mapconsf, 'one', lnd2atm_smap)
          call addmrg(fldListTo(compatm)%flds, 'Sl_fztop', &
               mrg_from1=complnd, mrg_fld1=trim(fldname), mrg_type1='copy')
       end if
    end if

    !-----------------------------------------------------------------------------
    ! to atm: dry deposition fields from land
    !-----------------------------------------------------------------------------
    if (phase == 'advertise') then
       do num = 1, max_ddep
          write(cnum,'(i2.2)') num
          fldname  = 'Sl_dd' // cnum
          call addfld(fldListFr(complnd)%flds, trim(fldname))
          call addfld(fldListTo(compatm)%flds, trim(fldname))
       end do
    else
       do num = 1, max_ddep
          write(cnum,'(i2.2)') num
          fldname  = 'Sl_dd' // cnum
          if ( fldchk(is_local%wrap%FBImp(complnd, complnd), trim(fldname), rc=rc) .and. &
               fldchk(is_local%wrap%FBExp(compatm)         , trim(fldname), rc=rc)) then
             call addmap(fldListFr(complnd)%flds, trim(fldname), compatm, mapconsf, 'one', lnd2atm_smap)
             call addmrg(fldListTo(compatm)%flds, trim(fldname), &
                  mrg_from1=complnd, mrg_fld1=trim(fldname), mrg_type1='copy')
          end if
       end do
    end if

    !=====================================================================
    ! to Ocean
    !=====================================================================

    ! ---------------------------------------------------------------------
    ! to ocn: zonal surface stress
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListMed_aoflux_o%flds, 'Faox_taux')
       call addfld(fldListFr(compice)%flds , 'Fioi_taux')
       call addfld(fldListFr(compatm)%flds , 'Faxa_taux')
       call addfld(fldListTo(compocn)%flds , 'Foxx_taux')
    else
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_taux', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compice,compice), 'Fioi_taux', rc=rc) .and. &
            fldchk(is_local%wrap%FBexp(compocn)        , 'Foxx_taux', rc=rc)) then

          ! NEMS orig
          call addmap(fldListFr(compatm)%flds, 'Faxa_taux', compocn, mapconsf, 'one'  , atm2ocn_fmap) ! map atm->ocn
          call addmap(fldListFr(compice)%flds, 'Fioi_taux', compocn, mapfcopy, 'unset', 'unset')
          ! custom merge calculation in med_phases_prep_ocn

       else if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_taux', rc=rc) .and. &
                 fldchk(is_local%wrap%FBexp(compocn)        , 'Foxx_taux', rc=rc)) then
          ! NEMS frac
          call addmap(fldListFr(compatm)%flds, 'Faxa_taux', compocn, mapconsf, 'one'  , atm2ocn_fmap) ! map atm->ocn
          call addmap(fldListFr(compice)%flds, 'Fioi_taux', compocn, mapfcopy, 'unset', 'unset')

          call addmrg(fldListTo(compocn)%flds, 'Foxx_taux', &
               mrg_from1=compatm, mrg_fld1='Faxa_taux', mrg_type1='merge', mrg_fracname1='ofrac', & ! ofrac=1-ifrac
               mrg_from2=compice, mrg_fld2='Fioi_taux', mrg_type2='merge', mrg_fracname2='ifrac')

       else if ( fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_taux', rc=rc) .and. &
                 fldchk(is_local%wrap%FBexp(compocn), 'Foxx_taux', rc=rc)) then
          ! CESM
          call addmap(fldListFr(compice)%flds, 'Fioi_taux', compocn, mapfcopy, 'unset', 'unset')

          call addmrg(fldListTo(compocn)%flds, 'Foxx_taux', &
               mrg_from1=compmed, mrg_fld1='Faox_taux', mrg_type1='merge', mrg_fracname1='ofrac', &
               mrg_from2=compice, mrg_fld2='Fioi_taux', mrg_type2='merge', mrg_fracname2='ifrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ocn: Meridional surface stress'
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListMed_aoflux_o%flds, 'Faox_tauy')
       call addfld(fldListFr(compice)%flds , 'Fioi_tauy')
       call addfld(fldListFr(compatm)%flds , 'Faxa_tauy')
       call addfld(fldListTo(compocn)%flds , 'Foxx_tauy')
    else
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_tauy', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compice,compice), 'Fioi_tauy', rc=rc) .and. &
            fldchk(is_local%wrap%FBexp(compocn)        , 'Foxx_tauy', rc=rc)) then
          ! NEMS orig
          call addmap(fldListFr(compatm)%flds, 'Faxa_tauy', compocn, mapconsf, 'one'  , atm2ocn_fmap) ! map atm->ocn
          call addmap(fldListFr(compice)%flds, 'Fioi_tauy', compocn, mapfcopy, 'unset', 'unset')
          ! TODO: implement a custom calculatino in med_phases_prep_ocn
          ! call addmrg(fldListTo(compocn)%flds, 'Foxx_tauy', &
          !      mrg_from1=compatm, mrg_fld1='Faxa_tauy', mrg_type1='merge', mrg_fracname1='atmwgt1', &
          !      mrg_from2=compice, mrg_fld2='Fioi_tauy', mrg_type2='merge', mrg_fracname2='icewgt1', &
          !      mrg_from3=compmed, mrg_fld3='Faox_tauy', mrg_type3='merge', mrg_fracname3='wgtm01'

       else if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_tauy', rc=rc) .and. &
                 fldchk(is_local%wrap%FBexp(compocn)        , 'Foxx_tauy', rc=rc)) then
          ! NEMS frac
          call addmap(fldListFr(compatm)%flds, 'Faxa_tauy', compocn, mapconsf, 'one'  , atm2ocn_fmap) ! map atm->ocn
          call addmap(fldListFr(compice)%flds, 'Fioi_tauy', compocn, mapfcopy, 'unset', 'unset')

          call addmrg(fldListTo(compocn)%flds, 'Foxx_tauy', &
               mrg_from1=compatm, mrg_fld1='Faxa_tauy', mrg_type1='merge', mrg_fracname1='ofrac', & ! ofrac=1-ifrac
               mrg_from2=compice, mrg_fld2='Fioi_tauy', mrg_type2='merge', mrg_fracname2='ifrac')

       else if ( fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_tauy', rc=rc) .and. &
                 fldchk(is_local%wrap%FBexp(compocn), 'Foxx_tauy', rc=rc)) then
          ! CESM
          call addmap(fldListFr(compice)%flds, 'Fioi_tauy', compocn, mapfcopy, 'unset', 'unset')

          call addmrg(fldListTo(compocn)%flds, 'Foxx_tauy', &
               mrg_from1=compmed, mrg_fld1='Faox_tauy', mrg_type1='merge', mrg_fracname1='ofrac', &
               mrg_from2=compice, mrg_fld2='Fioi_tauy', mrg_type2='merge', mrg_fracname2='ifrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ocn: surface latent heat flux
    ! ---------------------------------------------------------------------
    do n = 1,size(iso)
       if (phase == 'advertise') then
          call addfld(fldListMed_aoflux_o%flds, 'Faox_lat'//iso(n))
          call addfld(fldListTo(compocn)%flds , 'Foxx_lat'//iso(n))
       else
          if ( fldchk(is_local%wrap%FBexp(compocn), 'Foxx_lat'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_lat'//iso(n), rc=rc)) then
             ! CESM
             call addmrg(fldListTo(compocn)%flds, 'Foxx_lat'//iso(n), &
                  mrg_from1=compmed, mrg_fld1='Faox_lat'//iso(n), mrg_type1='merge', mrg_fracname1='ofrac')
          else
             ! Not used by NEMS orig and NEMS frac - latent derived from evap in med_phases_prep_ocn
          end if
       end if
    end do

    ! ---------------------------------------------------------------------
    ! to ocn: sensible heat flux
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds , 'Faxa_sen')
       call addfld(fldListMed_aoflux_o%flds, 'Faox_sen')
       call addfld(fldListFr(compice)%flds , 'Fioi_melth')
       call addfld(fldListTo(compocn)%flds , 'Foxx_sen')
    else
       if ( fldchk(is_local%wrap%FBexp(compocn)         , 'Foxx_sen'  , rc=rc) .and. &
            fldchk(is_local%wrap%FBMed_aoflux_o         , 'Faox_sen'  , rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compice, compice), 'Fioi_melth', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm, compatm), 'Faxa_sen'  , rc=rc)) then
          ! NEMS orig
          call addmap(fldListFr(compatm)%flds, 'Faxa_sen'  , compocn, mapconsf, 'one'  , atm2ocn_fmap) ! map atm->ocn
          call addmap(fldListFr(compice)%flds, 'Fioi_melth', compocn, mapfcopy, 'unset', 'unset')

       else if ( fldchk(is_local%wrap%FBexp(compocn)         , 'Foxx_sen', rc=rc) .and. &
                 fldchk(is_local%wrap%FBImp(compatm, compatm), 'Faxa_sen', rc=rc)) then
          ! NEMS frac
          call addmap(fldListFr(compatm)%flds, 'Faxa_sen', compocn, mapconsf, 'one'  , atm2ocn_fmap) ! map atm->ocn
          call addmrg(fldListTo(compocn)%flds, 'Foxx_sen', &
               mrg_from1=compatm, mrg_fld1='Faxa_sen'  , mrg_type1='merge', mrg_fracname1='ofrac', &
               mrg_from2=compice, mrg_fld2='Fioi_melth', mrg_type2='merge', mrg_fracname2='ifrac')

       else if ( fldchk(is_local%wrap%FBexp(compocn), 'Foxx_sen', rc=rc) .and. &
                 fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_sen', rc=rc)) then
          ! CESM
          call addmrg(fldListTo(compocn)%flds, 'Foxx_sen', &
               mrg_from1=compmed, mrg_fld1='Faox_sen', mrg_type1='merge', mrg_fracname1='ofrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ocn: heat flux from melting ice (CESM only for now)
    ! ---------------------------------------------------------------------
    ! TODO (mvertens, 2019-01-07): is fioi_melth being handled correctly here? Is fd.yaml correctly aliasing Fioi_melth?
    do n = 1,size(iso)
       if (phase == 'advertise') then
          call addfld(fldListFr(compice)%flds , 'Fioi_melth'//iso(n))
          call addfld(fldListTo(compocn)%flds , 'Fioi_melth'//iso(n))
       else
          if ( fldchk(is_local%wrap%FBexp(compocn)         , 'Fioi_melth'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compice, compice), 'Fioi_melth'//iso(n), rc=rc)) then
             call addmap(fldListFr(compice)%flds, 'Fioi_melth'//iso(n),    compocn,  mapfcopy, 'unset', 'unset')
             call addmrg(fldListTo(compocn)%flds, 'Fioi_melth'//iso(n), &
                  mrg_from1=compice, mrg_fld1='Fioi_melth'//iso(n), mrg_type1='copy_with_weights', mrg_fracname1='ifrac')
          end if
       end if
    end do

    ! ---------------------------------------------------------------------
    ! to ocn: net shortwave radiation
    ! ---------------------------------------------------------------------

    ! for budgets only
    if (phase == 'advertise') then
       call addfld(fldListFr(complnd)%flds, 'Fall_swnet')
       call addfld(fldListFr(compice)%flds, 'Faii_swnet')
       call addfld(fldListFr(compatm)%flds, 'Faxa_swnet')
    else
       if (fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_swnet', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_swnet', compice, mapconsf, 'one'  , atm2ice_fmap)
          call addmap(fldListFr(compatm)%flds, 'Faxa_swnet', compocn, mapconsf, 'one'  , atm2ocn_fmap)
       end if
       if (fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faii_swnet', rc=rc)) then
          call addmap(fldListFr(compice)%flds, 'Faii_swnet', compocn, mapfcopy, 'unset', 'unset')
       end if
    end if

    ! Shortwave radiation penetrating into ocean from ice (only mapping)
    if (phase == 'advertise') then
       call addfld(fldListFr(compice)%flds, 'Fioi_swpen')
       call addfld(fldListFr(compice)%flds, 'Fioi_swpen_vdr')
       call addfld(fldListFr(compice)%flds, 'Fioi_swpen_vdf')
       call addfld(fldListFr(compice)%flds, 'Fioi_swpen_idr')
       call addfld(fldListFr(compice)%flds, 'Fioi_swpen_idf')
    else
       if (fldchk(is_local%wrap%FBImp(compatm,compatm), 'Fioi_swpen', rc=rc)) then
          call addmap(fldListFr(compice)%flds, 'Fioi_swpen'    ,    compocn, mapfcopy, 'unset', 'unset')
       end if
       if (fldchk(is_local%wrap%FBImp(compatm,compatm), 'Fioi_swpen_vdr', rc=rc)) then
          call addmap(fldListFr(compice)%flds, 'Fioi_swpen_vdr', compocn, mapfcopy, 'unset', 'unset')
       end if
       if (fldchk(is_local%wrap%FBImp(compatm,compatm), 'Fioi_swpen_vdf', rc=rc)) then
          call addmap(fldListFr(compice)%flds, 'Fioi_swpen_vdf', compocn, mapfcopy, 'unset', 'unset')
       end if
       if (fldchk(is_local%wrap%FBImp(compatm,compatm), 'Fioi_swpen_idr', rc=rc)) then
          call addmap(fldListFr(compice)%flds, 'Fioi_swpen_idr', compocn, mapfcopy, 'unset', 'unset')
       end if
       if (fldchk(is_local%wrap%FBImp(compatm,compatm), 'Fioi_swpen_idf', rc=rc)) then
          call addmap(fldListFr(compice)%flds, 'Fioi_swpen_idf', compocn, mapfcopy, 'unset', 'unset')
       end if
    end if

    ! Net shortwave ocean (custom calculation in prep_phases_ocn_mod.F90)
    if (phase == 'advertise') then
       call addfld(fldListTo(compocn)%flds, 'Foxx_swnet')
       call addfld(fldListTo(compocn)%flds, 'Foxx_swnet_vdr')
       call addfld(fldListTo(compocn)%flds, 'Foxx_swnet_vdf')
       call addfld(fldListTo(compocn)%flds, 'Foxx_swnet_idr')
       call addfld(fldListTo(compocn)%flds, 'Foxx_swnet_idf')
    end if

    ! Per ice thickness fraction and sw penetrating into ocean from ice
    if (flds_i2o_per_cat) then
       if (phase == 'advertise') then
          ! 'fractional ice coverage wrt ocean for each thickness category '
          call addfld(fldListFr(compice)%flds, 'Si_ifrac_n')
          ! net shortwave radiation penetrating into ocean for each thickness category
          call addfld(fldListFr(compice)%flds, 'Fioi_swpen_ifrac_n')
          ! 'fractional atmosphere coverage wrt ocean'
          call addfld(fldListTo(compocn)%flds, 'Sf_afrac')
          ! 'net shortwave radiation times atmosphere fraction'
          call addfld(fldListTo(compocn)%flds, 'Foxx_swnet_afracr')
          ! 'fractional atmosphere coverage used in radiation computations wrt ocean'
          call addfld(fldListTo(compocn)%flds, 'Sf_afracr')
       else
          call addmap(fldListFr(compice)%flds, 'Si_ifrac_n', compocn, mapfcopy, 'unset', 'unset')
          call addmap(fldListFr(compice)%flds, 'Fioi_swpen_ifrac_n', compocn, mapfcopy, 'unset', 'unset')
          ! TODO (mvertens, 2018-12-21): add mapping and merging
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ocn: longwave heat fluxes
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListMed_aoflux_o%flds, 'Faox_lwup' )
       call addfld(fldListFr(compatm)%flds , 'Faxa_lwndn')
       call addfld(fldListFr(compatm)%flds , 'Foxx_lwnet')
       call addfld(fldListTo(compocn)%flds , 'Foxx_lwnet') ! mom6
       call addfld(fldListTo(compocn)%flds , 'Foxx_lwup' ) ! docn
       call addfld(fldListTo(compocn)%flds , 'Foxx_lwdn' ) ! docn
    else
       if ( fldchk(is_local%wrap%FBMed_aoflux_o        , 'Faox_lwup' , rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_lwdn' , rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compocn)        , 'Foxx_lwnet', rc=rc)) then
          ! CESM (mom6)
          call addmap(fldListFr(compatm)%flds, 'Faxa_lwdn', compocn, mapconsf, 'one'  , atm2ocn_fmap)
          call addmrg(fldListTo(compocn)%flds, 'Foxx_lwnet', &
               mrg_from1=compmed, mrg_fld1='Faox_lwup', mrg_type1='merge', mrg_fracname1='ofrac', &
               mrg_from2=compatm, mrg_fld2='Faxa_lwdn', mrg_type2='merge', mrg_fracname2='ofrac')

       else if ( fldchk(is_local%wrap%FBMed_aoflux_o        , 'Faox_lwup' , rc=rc) .and. &
                 fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_lwdn' , rc=rc)) then
          ! CESM (docn)
          call addmap(fldListFr(compatm)%flds, 'Faxa_lwdn', compocn, mapconsf, 'one'  , atm2ocn_fmap)
          call addmrg(fldListTo(compocn)%flds, 'Foxx_lwup', &
               mrg_from1=compmed, mrg_fld1='Faox_lwup', mrg_type1='merge', mrg_fracname1='ofrac')
          call addmrg(fldListTo(compocn)%flds, 'Faxa_lwdn', &
               mrg_from1=compmed, mrg_fld1='Faxa_lwdn', mrg_type1='merge', mrg_fracname1='ofrac')

       else if ( fldchk(is_local%wrap%FBMed_aoflux_o        , 'Faox_lwup' , rc=rc) .and. &
                 fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_lwdn' , rc=rc) .and. &
                 fldchk(is_local%wrap%FBImp(compatm,compatm), 'Foxx_lwnet', rc=rc) .and. &
                 fldchk(is_local%wrap%FBExp(compocn)        , 'Foxx_lwnet', rc=rc)) then
          ! NEMS-orig (mom6)
          call addmap(fldListFr(compatm)%flds, 'Faxa_lwdn' , compocn, mapconsf, 'one'  , atm2ocn_fmap)
          call addmap(fldListFr(compatm)%flds, 'Foxx_lwnet', compocn, mapconsf, 'one'  , atm2ocn_fmap)
          ! will have custom merge in med_phases_prep_ocn

       else if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Foxx_lwnet', rc=rc) .and. &
                 fldchk(is_local%wrap%FBExp(compocn)        , 'Foxx_lwnet', rc=rc)) then
          ! NEMS-frac (mom6)
          call addmap(fldListFr(compatm)%flds, 'Foxx_lwnet', compocn, mapconsf, 'one'  , atm2ocn_fmap)
          call addmrg(fldListTo(compocn)%flds, 'Foxx_lwnet', &
               mrg_from1=compatm, mrg_fld1='Foxx_lwnet', mrg_type1='merge', mrg_fracname1='ofrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ocn: evaporation water flux
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds , 'Faxx_evap') ! nems-orig, nems-frac
       call addfld(fldListTo(compocn)%flds , 'Foxx_evap') ! cesm
       call addfld(fldListMed_aoflux_o%flds, 'Faox_evap') ! cesm, nems-orig
    else
       if ( fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_evap', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compocn), 'Foxx_evap', rc=rc)) then
          ! CESM
          call addmrg(fldListTo(compocn)%flds, 'Foxx_evap', &
               mrg_from1=compmed, mrg_fld1='Faox_evap', mrg_type1='merge', mrg_fracname1='ofrac')
       else
          ! NEMS-frac and NEMS-orig have a custom calculation in med_phases_prep_ocn
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ocn: sea level pressure from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Sa_pslv')
       call addfld(fldListTo(compocn)%flds, 'Sa_pslv')
    else
       if ( fldchk(is_local%wrap%FBImp(compatm, compatm), 'Sa_pslv', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compocn)         , 'Sa_pslv', rc=rc)) then

          call addmap(fldListFr(compatm)%flds, 'Sa_pslv', compocn, mapbilnr, 'one', atm2ocn_smap)
          call addmap(fldListFr(compatm)%flds, 'Sa_pslv', compice, mapbilnr, 'one', atm2ocn_smap)

          call addmrg(fldListTo(compocn)%flds, 'Sa_pslv', &
               mrg_from1=compatm, mrg_fld1='Sa_pslv', mrg_type1='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ocn: wind speed squared at 10 meters from med
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListMed_aoflux_o%flds, 'So_duu10n')
       call addfld(fldListTo(compocn)%flds , 'So_duu10n')
    else
       if ( fldchk(is_local%wrap%FBMed_aoflux_o, 'So_duu10n', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compocn), 'So_duu10n', rc=rc)) then

          call addmap(fldListMed_aoflux_o%flds, 'So_duu10n', compatm, mapconsf, 'ofrac', ocn2atm_fmap) ! map ocn->atm
          call addmrg(fldListTo(compocn)%flds , 'So_duu10n', &
               mrg_from1=compmed, mrg_fld1='So_duu10n', mrg_type1='copy')
       end if
    end if


    ! ---------------------------------------------------------------------
    ! to ocn: water flux due to melting ice
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compice)%flds, 'Fioi_meltw')
       call addfld(fldListTo(compocn)%flds, 'Fioi_meltw')
    else
       if ( fldchk(is_local%wrap%FBImp(compice, compice), 'Fioi_meltw', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compocn)         , 'Fioi_meltw', rc=rc)) then

          call addmap(fldListFr(compice)%flds, 'Fioi_meltw', compocn,  mapfcopy, 'unset', 'unset')
          call addmrg(fldListTo(compocn)%flds, 'Fioi_meltw', &
               mrg_from1=compice, mrg_fld1='Fioi_meltw', mrg_type1='copy_with_weights', mrg_fracname1='ifrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ocn: salt flux from ice
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compice)%flds, 'Fioi_salt')
       call addfld(fldListTo(compocn)%flds, 'Fioi_salt')
    else
       if ( fldchk(is_local%wrap%FBImp(compice, compice), 'Fioi_salt', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compocn)         , 'Fioi_salt', rc=rc)) then

          call addmap(fldListFr(compice)%flds, 'Fioi_salt', compocn,  mapfcopy, 'unset', 'unset')
          call addmrg(fldListTo(compocn)%flds, 'Fioi_salt', &
               mrg_from1=compice, mrg_fld1='Fioi_salt', mrg_type1='copy_with_weights', mrg_fracname1='ifrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ocn: hydrophylic black carbon deposition flux from ice
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compice)%flds, 'Fioi_bcphi')
       call addfld(fldListTo(compocn)%flds, 'Fioi_bcphi')
    else
       if ( fldchk(is_local%wrap%FBImp(compice, compice), 'Fioi_bcphi', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compocn)         , 'Fioi_bcphi', rc=rc)) then

          call addmap(fldListFr(compice)%flds, 'Fioi_bcphi', compocn,  mapfcopy, 'unset', 'unset')
          call addmrg(fldListTo(compocn)%flds, 'Fioi_bcphi', &
               mrg_from1=compice, mrg_fld1='Fioi_bcphi', mrg_type1='copy_with_weights', mrg_fracname1='ifrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ocn: Hydrophobic black carbon deposition flux
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compice)%flds, 'Fioi_bcpho')
       call addfld(fldListTo(compocn)%flds, 'Fioi_bcpho')
    else
       if ( fldchk(is_local%wrap%FBImp(compice, compice), 'Fioi_bcpho', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compocn)         , 'Fioi_bcpho', rc=rc)) then

          call addmap(fldListFr(compice)%flds, 'Fioi_bcpho', compocn,  mapfcopy, 'unset', 'unset')
          call addmrg(fldListTo(compocn)%flds, 'Fioi_bcpho', &
               mrg_from1=compice, mrg_fld1='Fioi_bcpho', mrg_type1='copy_with_weights', mrg_fracname1='ifrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ocn: dust flux from ice
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compice)%flds, 'Fioi_flxdst')
       call addfld(fldListTo(compocn)%flds, 'Fioi_flxdst')
    else
       if ( fldchk(is_local%wrap%FBImp(compice, compice), 'Fioi_flxdst', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compocn)         , 'Fioi_flxdst', rc=rc)) then

          call addmap(fldListFr(compice)%flds, 'Fioi_flxdst', compocn,  mapfcopy, 'unset', 'unset')
          call addmrg(fldListTo(compocn)%flds, 'Fioi_flxdst', &
               mrg_from1=compice, mrg_fld1='Fioi_flxdst', mrg_type1='copy_with_weights', mrg_fracname1='ifrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ice: sea surface salinity
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compocn)%flds, 'So_s')
       call addfld(fldListTo(compice)%flds, 'So_s')
    else
       if ( fldchk(is_local%wrap%FBImp(compocn, compocn), 'So_s', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compocn)         , 'So_s', rc=rc)) then

          call addmap(fldListFr(compocn)%flds, 'So_s', compice,  mapfcopy, 'unset', 'unset')
          call addmrg(fldListTo(compice)%flds, 'So_s', &
               mrg_from1=compocn, mrg_fld1='So_s', mrg_type1='copy')
       end if
    end if

    !-----------------------------------------------------------------------------
    ! to ocn and land: nitrogen deposition fields
    !-----------------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_noy')
       call addfld(fldListFr(compatm)%flds, 'Faxa_nhx')
       call addfld(fldListTo(complnd)%flds, 'Faxa_noy')
       call addfld(fldListTo(complnd)%flds, 'Faxa_nhx')
       call addfld(fldListTo(compocn)%flds, 'Faxa_noy')
       call addfld(fldListTo(compocn)%flds, 'Faxa_nhx')
    else
       call addmap(fldListFr(compatm)%flds, 'Faxa_noy', complnd, mapbilnr, 'one', atm2lnd_smap)
       call addmap(fldListFr(compatm)%flds, 'Faxa_nhx', compocn, mapbilnr, 'one', atm2ocn_smap)

       call addmap(fldListFr(compatm)%flds, 'Faxa_noy', complnd, mapbilnr, 'one', atm2lnd_smap)
       call addmap(fldListFr(compatm)%flds, 'Faxa_nhx', compocn, mapbilnr, 'one', atm2ocn_smap)

       call addmrg(fldListTo(complnd)%flds, 'Faxa_noy', mrg_from1=compatm, mrg_fld1=trim(fldname), mrg_type1='copy_with_weights', mrg_fracname1='ofrac') !??? TODO
       call addmrg(fldListTo(complnd)%flds, 'Faxa_nhx', &
            mrg_from1=compatm, mrg_fld1=trim(fldname), mrg_type1='copy_with_weights', mrg_fracname1='ofrac') !??? TODO

       call addmrg(fldListTo(compocn)%flds, 'Faxa_noy', &
            mrg_from1=compatm, mrg_fld1=trim(fldname), mrg_type1='copy_with_weights', mrg_fracname1='ofrac')
       call addmrg(fldListTo(compocn)%flds, 'Faxa_nhx', &
            mrg_from1=compatm, mrg_fld1=trim(fldname), mrg_type1='copy_with_weights', mrg_fracname1='ofrac')
    end if

    !-----------------------------
    ! to ocn: liquid runoff from rof and glc components
    !-----------------------------
    do n = 1,size(iso)
       if (phase == 'advertise') then
          call addfld(fldListFr(compglc)%flds, 'Fogg_rofl'//iso(n))
          call addfld(fldListFr(comprof)%flds, 'Forr_rofl'//iso(n))
          call addmrg(fldListTo(compocn)%flds, 'Foxx_rofl'//iso(n))
       else
          if ( fldchk(is_local%wrap%FBImp(comprof, comprof), 'Forr_rofl'//iso(n), rc=rc)) then
             ! water flux into sea water due to runoff (liquid)
             call addmap(fldListFr(comprof)%flds, 'Forr_rofl'//iso(n), compocn, mapfiler, 'none', rof2ocn_liq_rmap)
          end if
          if ( fldchk(is_local%wrap%FBImp(comprof, comprof), 'Fogg_rofl'//iso(n), rc=rc)) then
             ! glc liquid runoff flux to ocean
             call addmap(fldListFr(compglc)%flds, 'Forr_rofl'//iso(n), compocn,  mapfiler, 'one', glc2ocn_liq_rmap)
          end if
          if ( fldchk(is_local%wrap%FBExp(comprof), 'Foxx_rofl'//iso(n), rc=rc)) then
             if ( fldchk(is_local%wrap%FBImp(comprof, comprof), 'Forr_rofl'//iso(n), rc=rc) .and. &
                  fldchk(is_local%wrap%FBImp(compglc, compglc), 'Fogg_rofl'//iso(n), rc=rc)) then
                call addmrg(fldListTo(compocn)%flds, 'Foxx_rofl'//iso(n), &
                     mrg_from1=comprof, mrg_fld1='Forr_rofl:Flrr_flood', mrg_type1='sum', &
                     mrg_from2=compglc, mrg_fld2='Fogg_rofl'//iso(n)           , mrg_type2='sum')

             else if ( fldchk(is_local%wrap%FBImp(comprof, comprof), 'Forr_rofl'//iso(n), rc=rc)) then
                call addmrg(fldListTo(compocn)%flds, 'Foxx_rofl'//iso(n), &
                     mrg_from1=comprof, mrg_fld1='Forr_rofl:Flrr_flood', mrg_type1='sum')

             else if (fldchk(is_local%wrap%FBImp(compglc, compglc), 'Fogg_rofl'//iso(n), rc=rc)) then
                call addmrg(fldListTo(compocn)%flds, 'Foxx_rofl'//iso(n), &
                     mrg_from1=compglc, mrg_fld1='Fogg_rofl'//iso(n), mrg_type1='sum')
             end if
          end if
       end if
    end do

    !-----------------------------
    ! to ocn: frozen runoff from rof and glc components
    !-----------------------------
    do n = 1,size(iso)
       if (phase == 'advertise') then
          call addfld(fldListFr(compglc)%flds, 'Fogg_rofi'//iso(n))
          call addfld(fldListFr(comprof)%flds, 'Forr_rofi'//iso(n)) ! rof ice runoff flux to ocean
          call addmrg(fldListTo(compocn)%flds, 'Foxx_rofi'//iso(n)) ! glc ice runoff flux to ocean
       else
          if ( fldchk(is_local%wrap%FBImp(comprof, comprof), 'Forr_rofi'//iso(n), rc=rc)) then
             call addmap(fldListFr(comprof)%flds, 'Forr_rofi'//iso(n), compocn, mapfiler, 'none', rof2ocn_liq_rmap)
          end if
          if ( fldchk(is_local%wrap%FBImp(comprof, comprof), 'Fogg_rofi'//iso(n), rc=rc)) then
             call addmap(fldListFr(compglc)%flds, 'Forr_rofi'//iso(n), compocn,  mapfiler, 'one', glc2ocn_liq_rmap)
          end if
          if ( fldchk(is_local%wrap%FBExp(comprof), 'Foxx_rofi'//iso(n), rc=rc)) then
             if ( fldchk(is_local%wrap%FBImp(comprof, comprof), 'Forr_rofi'//iso(n), rc=rc) .and. &
                  fldchk(is_local%wrap%FBImp(compglc, compglc), 'Fogg_rofi'//iso(n), rc=rc)) then
                call addmrg(fldListTo(compocn)%flds, 'Foxx_rofi'//iso(n), &
                     mrg_from1=comprof, mrg_fld1='Forr_rofi:Flrr_flood', mrg_type1='sum', &
                     mrg_from2=compglc, mrg_fld2='Fogg_rofi'//iso(n)           , mrg_type2='sum')

             else if ( fldchk(is_local%wrap%FBImp(comprof, comprof), 'Forr_rofi'//iso(n), rc=rc)) then
                call addmrg(fldListTo(compocn)%flds, 'Foxx_rofi'//iso(n), &
                     mrg_from1=comprof, mrg_fld1='Forr_rofi:Flrr_flood', mrg_type1='sum')

             else if (fldchk(is_local%wrap%FBImp(compglc, compglc), 'Fogg_rofi'//iso(n), rc=rc)) then
                call addmrg(fldListTo(compocn)%flds, 'Foxx_rofi'//iso(n), &
                     mrg_from1=compglc, mrg_fld1='Fogg_rofi'//iso(n), mrg_type1='sum')
             end if
          end if
       end if
    end do

    !-----------------------------
    ! to ocn: Langmuir multiplier from wave
    !-----------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compwav)%flds, 'Sw_lamult')
       call addfld(fldListTo(compocn)%flds, 'Sw_lamult')
    else
       if ( fldchk(is_local%wrap%FBExp(compocn)         , 'Sw_lamult', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compwav, compwav), 'Sw_lamult', rc=rc)) then

          call addmap(fldListFr(compwav)%flds, 'Sw_lamult', compocn,  mapbilnr, 'one', wav2ocn_smap)
          call addmrg(fldListTo(compocn)%flds, 'Sw_lamult', &
               mrg_from1=compwav, mrg_fld1='Sw_lamult', mrg_type1='copy')
       end if
    end if

    !-----------------------------
    ! to ocn: Stokes drift u component from wave
    !-----------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compwav)%flds, 'Sw_ustokes')
       call addfld(fldListTo(compocn)%flds, 'Sw_ustokes')
    else
       if ( fldchk(is_local%wrap%FBExp(compocn)         , 'Sw_ustokes', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compwav, compwav), 'Sw_ustokes', rc=rc)) then

          call addmap(fldListFr(compwav)%flds, 'Sw_ustokes', compocn,  mapbilnr, 'one', wav2ocn_smap)
          call addmrg(fldListTo(compocn)%flds, 'Sw_ustokes', &
               mrg_from1=compwav, mrg_fld1='Sw_ustokes', mrg_type1='copy')
       end if
    end if

    !-----------------------------
    ! to ocn: Stokes drift v component from wave
    !-----------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compwav)%flds, 'Sw_vstokes')
       call addfld(fldListTo(compocn)%flds, 'Sw_vstokes')
    else
       if ( fldchk(is_local%wrap%FBExp(compocn)         , 'Sw_vstokes', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compwav, compwav), 'Sw_vstokes', rc=rc)) then

          call addmap(fldListFr(compwav)%flds, 'Sw_vstokes', compocn,  mapbilnr, 'one', wav2ocn_smap)
          call addmrg(fldListTo(compocn)%flds, 'Sw_vstokes', &
               mrg_from1=compwav, mrg_fld1='Sw_vstokes', mrg_type1='copy')
       end if
    end if

    !-----------------------------
    ! to ocn: Stokes drift depth from wave
    !-----------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compwav)%flds, 'Sw_hstokes')
       call addfld(fldListTo(compocn)%flds, 'Sw_hstokes')
    else
       if ( fldchk(is_local%wrap%FBExp(compocn)         , 'Sw_hstokes', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compwav, compwav), 'Sw_hstokes', rc=rc)) then

          call addmap(fldListFr(compwav)%flds, 'Sw_hstokes', compocn,  mapbilnr, 'one', wav2ocn_smap)
          call addmrg(fldListTo(compocn)%flds, 'Sw_hstokes', &
               mrg_from1=compwav, mrg_fld1='Sw_hstokes', mrg_type1='copy')
       end if
    end if

    !=====================================================================
    ! to ice (and wave)
    !=====================================================================

    ! ---------------------------------------------------------------------
    ! to ice and wave: zonal sea water velocity from ocean
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compocn)%flds, 'So_u')
       call addfld(fldListTo(compice)%flds, 'So_u')
       call addfld(fldListTo(compwav)%flds, 'So_u')
    else
       if ( fldchk(is_local%wrap%FBImp(compocn, compocn), 'So_u', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compice)         , 'So_u', rc=rc)) then
          call addmap(fldListFr(compocn)%flds,' So_u', compice,  mapfcopy, 'unset', 'unset')
          call addmrg(fldListTo(compice)%flds, 'So_u', mrg_from1=compocn, mrg_fld1='So_u', mrg_type1='copy')
       end if
       if ( fldchk(is_local%wrap%FBImp(compocn, compocn), 'So_u', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compice)         , 'So_u', rc=rc)) then
          call addmap(fldListFr(compocn)%flds, 'So_u', compwav,  mapbilnr, 'one'  , 'ocn2wav_smap')
          call addmrg(fldListTo(compwav)%flds, 'So_u', mrg_from1=compocn, mrg_fld1='So_u', mrg_type1='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ice and wave: meridional sea water velocity from ocean
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compocn)%flds, 'So_v')
       call addfld(fldListTo(compice)%flds, 'So_v')
       call addfld(fldListTo(compwav)%flds, 'So_v')
    else
       if ( fldchk(is_local%wrap%FBImp(compocn, compocn), 'So_v', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compice)         , 'So_v', rc=rc)) then
          call addmap(fldListFr(compocn)%flds,' So_v', compice,  mapfcopy, 'unset', 'unset')
          call addmrg(fldListTo(compice)%flds, 'So_v', mrg_from1=compocn, mrg_fld1='So_v', mrg_type1='copy')
       end if
       if ( fldchk(is_local%wrap%FBImp(compocn, compocn), 'So_v', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compice)         , 'So_v', rc=rc)) then
          call addmap(fldListFr(compocn)%flds, 'So_v', compwav,  mapbilnr, 'one'  , 'ocn2wav_smap')
          call addmrg(fldListTo(compwav)%flds, 'So_v', mrg_from1=compocn, mrg_fld1='So_v', mrg_type1='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ice: zonal sea surface slope from ocean
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compocn)%flds, 'So_dhdx')
       call addfld(fldListTo(compice)%flds, 'So_dhdx')
    else
       if ( fldchk(is_local%wrap%FBImp(compocn, compocn), 'So_dhdx', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compice)         , 'So_dhdx', rc=rc)) then
          call addmap(fldListFr(compocn)%flds, 'So_dhdx', compice,  mapfcopy, 'unset', 'unset')
          call addmrg(fldListTo(compice)%flds, 'So_dhdx', mrg_from1=compocn, mrg_fld1='So_dhdx', mrg_type1='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ice: meridional sea surface slope from ocean
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compocn)%flds, 'So_dhdy')
       call addfld(fldListTo(compice)%flds, 'So_dhdy')
    else
       if ( fldchk(is_local%wrap%FBImp(compocn, compocn), 'So_dhdy', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compice)         , 'So_dhdy', rc=rc)) then
          call addmap(fldListFr(compocn)%flds, 'So_dhdy', compice,  mapfcopy, 'unset', 'unset')
          call addmrg(fldListTo(compice)%flds, 'So_dhdy', mrg_from1=compocn, mrg_fld1='So_dhdy', mrg_type1='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ice: ocean melt and freeze potential
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compocn)%flds, 'Fioo_q')
       call addfld(fldListTo(compice)%flds, 'Fioo_q')
    else
       if ( fldchk(is_local%wrap%FBImp(compocn, compocn), 'Fioo_q', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compice)         , 'Fioo_q', rc=rc)) then
          call addmap(fldListFr(compocn)%flds, 'Fioo_q', compice,  mapfcopy, 'unset', 'unset')
          call addmrg(fldListTo(compice)%flds, 'Fioo_q', mrg_from1=compocn, mrg_fld1='Fioo_q', mrg_type1='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to wave: ocean boundary layer depth
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compocn)%flds, 'So_bldepth')
       call addfld(fldListTo(compwav)%flds, 'So_bldepth')
    else
       if ( fldchk(is_local%wrap%FBImp(compocn, compocn), 'So_bldepth', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compwav)         , 'So_bldepth', rc=rc)) then
          call addmap(fldListFr(compocn)%flds, 'So_bldepth', compwav,  mapfcopy, 'unset', 'unset')
          call addmrg(fldListTo(compwav)%flds, 'So_bldepth', mrg_from1=compocn, mrg_fld1='So_bldepth', mrg_type1='copy')
       end if
    end if

    !=====================================================================
    ! to Runoff (river component)
    !=====================================================================

    ! ---------------------------------------------------------------------
    ! to rof: water flux from land (liquid surface)
    ! ---------------------------------------------------------------------
    do n = 1,size(iso)
       if (phase == 'advertise') then
          call addfld(fldListFr(complnd)%flds, 'Flrl_rofsur'//iso(n))
          call addfld(fldListTo(comprof)%flds, 'Flrl_rofsur'//iso(n))
       else
          if ( fldchk(is_local%wrap%FBImp(complnd, complnd), 'Flrl_rofsur'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBExp(comprof)         , 'Flrl_rofsur'//iso(n), rc=rc)) then
             call addmap(fldListFr(complnd)%flds, 'Flrl_rofsur'//iso(n), comprof, mapconsf, 'lfrin', lnd2rof_fmap)
             call addmrg(fldListTo(comprof)%flds, 'Flrl_rofsur'//iso(n), &
                  mrg_from1=complnd, mrg_fld1='Flrl_rofsur'//iso(n), mrg_type1='copy_with_weights', mrg_fracname1='lfrac')
          end if
       end if
    end do

    ! ---------------------------------------------------------------------
    ! to rof: water flux from land (liquid glacier, wetland, and lake)
    ! ---------------------------------------------------------------------
    do n = 1,size(iso)
       if (phase == 'advertise') then
          call addfld(fldListFr(complnd)%flds, 'Flrl_rofgwl'//iso(n))
          call addfld(fldListTo(comprof)%flds, 'Flrl_rofgwl'//iso(n))
       else
          if ( fldchk(is_local%wrap%FBImp(complnd, complnd), 'Flrl_rofgwl'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBExp(comprof)         , 'Flrl_rofgwl'//iso(n), rc=rc)) then
             call addmap(fldListFr(complnd)%flds, 'Flrl_rofgwl'//iso(n), comprof, mapconsf, 'lfrin', lnd2rof_fmap)
             call addmrg(fldListTo(comprof)%flds, 'Flrl_rofgwl'//iso(n), &
                  mrg_from1=complnd, mrg_fld1='Flrl_rofgwl'//iso(n), mrg_type1='copy_with_weights', mrg_fracname1='lfrac')
          end if
       end if
    end do

    ! ---------------------------------------------------------------------
    ! to rof: water flux from land (liquid subsurface)
    ! ---------------------------------------------------------------------
    do n = 1,size(iso)
       if (phase == 'advertise') then
          call addfld(fldListFr(complnd)%flds, 'Flrl_rofsub'//iso(n))
          call addfld(fldListTo(comprof)%flds, 'Flrl_rofsub'//iso(n))
       else
          if ( fldchk(is_local%wrap%FBImp(complnd, complnd), 'Flrl_rofsub'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBExp(comprof)         , 'Flrl_rofsub'//iso(n), rc=rc)) then
             call addmap(fldListFr(complnd)%flds, 'Flrl_rofsub'//iso(n), comprof, mapconsf, 'lfrin', lnd2rof_fmap)
             call addmrg(fldListTo(comprof)%flds, 'Flrl_rofsub'//iso(n), &
                  mrg_from1=complnd, mrg_fld1='Flrl_rofsub'//iso(n), mrg_type1='copy_with_weights', mrg_fracname1='lfrac')
          end if
       end if
    end do

    ! ---------------------------------------------------------------------
    ! to rof: water flux from land direct to ocean
    ! ---------------------------------------------------------------------
    do n = 1,size(iso)
       if (phase == 'advertise') then
          call addfld(fldListFr(complnd)%flds, 'Flrl_rofdto'//iso(n))
          call addfld(fldListTo(comprof)%flds, 'Flrl_rofdto'//iso(n))
       else
          if ( fldchk(is_local%wrap%FBImp(complnd, complnd), 'Flrl_rofdto'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBExp(comprof)         , 'Flrl_rofdto'//iso(n), rc=rc)) then
             call addmap(fldListFr(complnd)%flds, 'Flrl_rofdto'//iso(n), comprof, mapconsf, 'lfrin', lnd2rof_fmap)
             call addmrg(fldListTo(comprof)%flds, 'Flrl_rofdto'//iso(n), &
                  mrg_from1=complnd, mrg_fld1='Flrl_rofdto'//iso(n), mrg_type1='copy_with_weights', mrg_fracname1='lfrac')
          end if
       end if
    end do

    ! ---------------------------------------------------------------------
    ! to rof: water flux from land (frozen)
    ! ---------------------------------------------------------------------
    do n = 1,size(iso)
       if (phase == 'advertise') then
          call addfld(fldListFr(complnd)%flds, 'Flrl_rofi'//iso(n))
          call addfld(fldListTo(comprof)%flds, 'Flrl_rofi'//iso(n))
       else
          if ( fldchk(is_local%wrap%FBImp(complnd, complnd), 'Flrl_rofi'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBExp(comprof)         , 'Flrl_rofi'//iso(n), rc=rc)) then
             call addmap(fldListFr(complnd)%flds, 'Flrl_rofi'//iso(n), comprof, mapconsf, 'lfrin', lnd2rof_fmap)
             call addmrg(fldListTo(comprof)%flds, 'Flrl_rofi'//iso(n), &
                  mrg_from1=complnd, mrg_fld1='Flrl_rofi'//iso(n), mrg_type1='copy_with_weights', mrg_fracname1='lfrac')
          end if
       end if
    end do

    ! ---------------------------------------------------------------------
    ! to rof: irrigation flux (withdrawal from rivers)
    ! ---------------------------------------------------------------------
    do n = 1,size(iso)
       if (phase == 'advertise') then
          call addfld(fldListFr(complnd)%flds, 'Flrl_irrig'//iso(n))
          call addfld(fldListTo(comprof)%flds, 'Flrl_irrig'//iso(n))
       else
          if ( fldchk(is_local%wrap%FBImp(complnd, complnd), 'Flrl_irrig'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBExp(comprof)         , 'Flrl_irrig'//iso(n), rc=rc)) then
             call addmap(fldListFr(complnd)%flds, 'Flrl_irrig'//iso(n), comprof, mapconsf, 'lfrin', lnd2rof_fmap)
             call addmrg(fldListTo(comprof)%flds, 'Flrl_irrig'//iso(n), &
                  mrg_from1=complnd, mrg_fld1='Flrl_irrig'//iso(n), mrg_type1='copy_with_weights', mrg_fracname1='lfrac')
          end if
       end if
    end do

    ! ---------------------------------------------------------------------
    ! to lnd and to ocn: waterflux back to land due to flooding
    ! ---------------------------------------------------------------------
    do n = 1,size(iso)
       if (phase == 'advertise') then
          call addfld(fldListFr(comprof)%flds, 'Flrr_flood'//iso(n))
          call addfld(fldListTo(complnd)%flds, 'Flrr_flood'//iso(n))
          call addfld(fldListTo(compocn)%flds, 'Flrr_flood'//iso(n))
       else
          if ( fldchk(is_local%wrap%FBImp(comprof, comprof), 'Flrr_flood'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBExp(comprof)         , 'Flrr_flood'//iso(n), rc=rc)) then
             call addmap(fldListFr(comprof)%flds, 'Flrr_flood'//iso(n), complnd, mapconsf, 'one', rof2lnd_fmap)
             call addmrg(fldListTo(complnd)%flds, 'Flrr_flood'//iso(n), &
                  mrg_from1=comprof, mrg_fld1='Flrr_flood'//iso(n), mrg_type1='copy')
          end if
          if ( fldchk(is_local%wrap%FBImp(comprof, comprof), 'Flrr_flood'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBExp(compocn)         , 'Flrr_flood'//iso(n), rc=rc)) then
             call addmap(fldListFr(comprof)%flds, 'Flrr_flood'//iso(n), compocn, mapconsf, 'one', rof2ocn_fmap)
             ! TODO: custom merge for flood - how to handle consistency for flood back to ocean
          end if
       end if
    end do

    ! ---------------------------------------------------------------------
    ! to lnd: 'river channel total water volume
    ! ---------------------------------------------------------------------
    do n = 1,size(iso)
       if (phase == 'advertise') then
          call addfld(fldListFr(comprof)%flds, 'Flrr_volr'//iso(n))
          call addfld(fldListTo(complnd)%flds, 'Flrr_volr'//iso(n))
          call addfld(fldListTo(compocn)%flds, 'Flrr_volr'//iso(n))
       else
          if ( fldchk(is_local%wrap%FBImp(comprof, comprof), 'Flrr_volr'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBExp(comprof)         , 'Flrr_volr'//iso(n), rc=rc)) then
             call addmap(fldListFr(comprof)%flds, 'Flrr_volr'//iso(n), complnd, mapconsf, 'one', rof2lnd_fmap)
             call addmrg(fldListTo(complnd)%flds, 'Flrr_volr'//iso(n), &
                  mrg_from1=comprof, mrg_fld1='Flrr_volr'//iso(n), mrg_type1='copy')
          end if
       end if
    end do

    ! ---------------------------------------------------------------------
    ! to lnd: River channel main channel water volume
    ! ---------------------------------------------------------------------
    do n = 1,size(iso)
       if (phase == 'advertise') then
          call addfld(fldListFr(comprof)%flds, 'Flrr_volrmch'//iso(n))
          call addfld(fldListTo(complnd)%flds, 'Flrr_volrmch'//iso(n))
          call addfld(fldListTo(compocn)%flds, 'Flrr_volrmch'//iso(n))
       else
          if ( fldchk(is_local%wrap%FBImp(comprof, comprof), 'Flrr_volrmch'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBExp(comprof)         , 'Flrr_volrmch'//iso(n), rc=rc)) then
             call addmap(fldListFr(comprof)%flds, 'Flrr_volrmch'//iso(n), complnd, mapconsf, 'one', rof2lnd_fmap)
             call addmrg(fldListTo(complnd)%flds, 'Flrr_volrmch'//iso(n), &
                  mrg_from1=comprof, mrg_fld1='Flrr_volrmch'//iso(n), mrg_type1='copy')
          end if
       end if
    end do

    ! ---------------------------------------------------------------------
    ! to ice and glc: frozen runoff
    ! ---------------------------------------------------------------------
    do n = 1,size(iso)
       if (phase == 'advertise') then
          call addfld(fldListFr(comprof)%flds, 'Firr_rofi'//iso(n)) ! 'Water flux into sea ice due to runoff (frozen)'
          call addfld(fldListFr(compglc)%flds, 'Figg_rofi'//iso(n)) ! 'glc frozen runoff_iceberg flux to ice'
          call addfld(fldListTo(compice)%flds, 'Fixx_rofi'//iso(n)) ! 'Total frozen water flux into sea ice '
       else
          if ( fldchk(is_local%wrap%FBExp(compice)         , 'Fixx_rofi'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(comprof, comprof), 'Forr_rofi'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compglc, compglc), 'Figg_rofi'//iso(n), rc=rc)) then

             call addmap(fldListFr(comprof)%flds, 'Forr_rofi'//iso(n), compice, mapfiler, 'none', rof2ocn_ice_rmap)
             call addmap(fldListFr(compglc)%flds, 'Figg_rofi'//iso(n), compice, mapfiler, 'one' , glc2ice_rmap)
             call addmrg(fldListTo(compice)%flds, 'Fixx_rofi'//iso(n), &
                  mrg_from1=comprof, mrg_fld1='Firr_rofi'//iso(n), mrg_type1='sum', &
                  mrg_from2=compglc, mrg_fld2='Figg_rofi'//iso(n), mrg_type2='sum')

          else if ( fldchk(is_local%wrap%FBExp(compice)    , 'Fixx_rofi'//iso(n), rc=rc) .and. &
                    fldchk(is_local%wrap%FBImp(comprof, comprof), 'Forr_rofi'//iso(n), rc=rc)) then

             call addmap(fldListFr(comprof)%flds, 'Forr_rofi'//iso(n), compice, mapfiler, 'none', rof2ocn_ice_rmap)
             call addmrg(fldListTo(compice)%flds, 'Fixx_rofi'//iso(n), &
                  mrg_from1=comprof, mrg_fld1='Firr_rofi'//iso(n), mrg_type1='sum')
          end if
       end if
    end do

    !-----------------------------
    ! to lnd: glc -> lnd
    !-----------------------------

    ! ---------------------------------------------------------------------
    ! to lnd: ice sheet grid coverage on global grid
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compglc)%flds   , 'Sg_icemask')
       call addfld(fldListTo(complnd)%flds   , 'Sg_icemask')
       call addfld(fldListMed_x2l_fr_glc%flds, 'Sg_icemask') ! Needed for FB initialization
    else
       if ( fldchk(is_local%wrap%FBExp(complnd)         , 'Sg_icemask', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compglc, compglc), 'Sg_icemask', rc=rc)) then

          call addmap(fldListFr(compglc)%flds, 'Sg_icemask', complnd,  mapconsf, 'one', glc2lnd_smap)
          call addmrg(fldListTo(complnd)%flds, 'Sg_icemask', &
               mrg_from1=compglc, mrg_fld1='Sg_icemask', mrg_type1='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to lnd: ice sheet mask where we are potentially sending non-zero fluxes
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compglc)%flds    , 'Sg_icemask_coupled_fluxes')
       call addfld(fldListTo(complnd)%flds    , 'Sg_icemask_coupled_fluxes')
       call addfld(fldListMed_x2l_fr_glc%flds , 'Sg_icemask_coupled_fluxes')
    else
       if ( fldchk(is_local%wrap%FBExp(complnd)        , 'Sg_icemask_coupled_fluxes', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compglc,compglc), 'Sg_icemask_coupled_fluxes', rc=rc)) then
          call addmap(fldListFr(compglc)%flds, 'Sg_icemask_coupled_fluxes', complnd,  mapconsf, 'one', glc2lnd_smap)
          call addmrg(fldListTo(complnd)%flds, 'Sg_icemask_coupled_fluxes', &
               mrg_from1=compglc, mrg_fld1='Sg_icemask_coupled_fluxes', mrg_type1='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to lnd: fields with multiple elevation classes
    ! ---------------------------------------------------------------------

    ! for glc fields with multiple elevation classes in glc->lnd
    ! fields from glc->med do NOT have elevation classes
    ! fields from med->lnd are BROKEN into multiple elevation classes


    if (glc_nec > 0) then
       if (phase == 'advertise') then
          call addfld(fldListFr(compglc)%flds, 'Sg_ice_covered') ! fraction of glacier area
          call addfld(fldListFr(compglc)%flds, 'Sg_topo')        ! surface height of glacer
          call addfld(fldListFr(compglc)%flds, 'Flgg_hflx')      ! downward heat flux from glacier interior
          do num = 0, glc_nec
             cnum = glc_elevclass_as_string(num)
             call addfld(fldListMed_x2l_fr_glc%flds, 'Sg_ice_covered'//trim(cnum))
             call addfld(fldListMed_x2l_fr_glc%flds, 'Sg_topo'//trim(cnum))
             call addfld(fldListMed_x2l_fr_glc%flds, 'Flgg_hflx'//trim(cnum))
          end do
       else
          call addmap(FldListFr(compglc)%flds, 'Sg_ice_covered' , complnd, mapconsf, 'unset' , glc2lnd_fmap) ! TODO: normalization?
          call addmap(FldListFr(compglc)%flds, 'Sg_topo'        , compglc, mapconsf, 'custom', glc2lnd_fmap)
          call addmap(FldListFr(compglc)%flds, 'Flgg_hflx'      , compglc, mapconsf, 'custom', glc2lnd_fmap)
          do num = 0, glc_nec
             cnum = glc_elevclass_as_string(num)
             call addmrg(fldListTo(complnd)%flds, 'Sg_ice_covered'//trim(cnum), &
                  mrg_from1=compglc, mrg_fld1='Sg_ice_covered'//trim(cnum), mrg_type1='copy')
             call addmrg(fldListTo(complnd)%flds, 'Sg_topo' //trim(cnum),  &
                  mrg_from1=compglc, mrg_fld1='Sg_topo'//trim(cnum), mrg_type1='copy')
             call addmrg(fldListTo(complnd)%flds, 'Flgg_hflx'//trim(cnum), &
                  mrg_from1=compglc, mrg_fld1='Flgg_hflx'//trim(cnum), mrg_type1='copy')
          end do
       end if
    end if

    !-----------------------------
    ! to glc: from land fields with multiple elevation classes
    !-----------------------------

    ! glc fields with multiple elevation classes: lnd->glc
    ! - fields sent from lnd->med ARE in multiple elevation classes
    ! - fields sent from med->glc do NOT have elevation classes

    ! Sets a coupling field for all glc elevation classes (1:glc_nec) plus bare land (index 0).
    ! Note that, if glc_nec = 0, then we don't create any coupling fields (not even the bare land (0) fldindex)
    ! Note : Sl_topo is sent from lnd -> med, but is NOT sent to glc (only used for the remapping in the mediator)

    if (glc_nec > 0) then
       if (phase == 'advertise') then
          do num = 0, glc_nec
             cnum = glc_elevclass_as_string(num)
             call addfld(fldListFr(complnd)%flds, 'Flgl_qice'//trim(cnum)) ! glacier ice flux'
             call addfld(fldListFr(complnd)%flds, 'Sl_tsrf'  //trim(cnum)) ! surface temperature of glacier'
             call addfld(fldListFr(complnd)%flds, 'Sl_topo'  //trim(cnum)) ! surface height of glacier
          end do
          call addfld(fldListTo(compglc)%flds, 'Flgl_qice')
          call addfld(fldListTo(compglc)%flds, 'Sl_tsrf')
          call addfld(fldListTo(compglc)%flds, 'Sl_topo')
       else
          do num = 0, glc_nec
             cnum = glc_elevclass_as_string(num)
             call addmap(FldListFr(complnd)%flds, 'Flgl_qice', compglc, mapconsf, 'none', lnd2glc_fmap)
             call addmap(FldListFr(complnd)%flds, 'Sl_tsrf'  , compglc, mapbilnr, 'none', lnd2glc_smap)
             call addmap(FldListFr(complnd)%flds, 'Sl_topo'  , compglc, mapbilnr, 'none', lnd2glc_smap)
          end do
       end if
    end if

    ! ---------------------------------------------------------------------
    ! co2 exchange
    ! ---------------------------------------------------------------------

    if (flds_co2a) then
       ! ---------------------------------------------------------------------
       ! to lnd and ocn: prognostic CO2 at the lowest atm model level
       ! ---------------------------------------------------------------------
       if (phase == 'advertise') then
          call addfld(fldListFr(compatm)%flds, 'Sa_co2prog')
          call addfld(fldListTo(complnd)%flds, 'Sa_co2prog')
          call addfld(fldListTo(compocn)%flds, 'Sa_co2prog')
       else
          call addmap(fldListFr(compatm)%flds, 'Sa_co2prog', complnd, mapbilnr, 'one', atm2lnd_smap)
          call addmap(fldListFr(compatm)%flds, 'Sa_co2prog', compocn, mapbilnr, 'one', atm2ocn_smap)

          call addmrg(fldListTo(complnd)%flds, 'Sa_co2prog', &
               mrg_from1=compatm, mrg_fld1='Sa_co2prog', mrg_type1='copy')
          call addmrg(fldListTo(compocn)%flds, 'Sa_co2prog', &
               mrg_from1=compatm, mrg_fld1='Sa_co2prog', mrg_type1='copy')
       end if

       ! ---------------------------------------------------------------------
       ! to lnd and ocn: diagnostic CO2 at the lowest atm model level
       ! ---------------------------------------------------------------------
       if (phase == 'advertise') then
          call addfld(fldListFr(compatm)%flds, 'Sa_co2diag')
          call addfld(fldListTo(complnd)%flds, 'Sa_co2diag')
          call addfld(fldListTo(compocn)%flds, 'Sa_co2diag')
       else
          call addmap(fldListFr(compatm)%flds, 'Sa_co2diag', complnd, mapbilnr, 'one', atm2lnd_smap)
          call addmap(fldListFr(compatm)%flds, 'Sa_co2diag', compocn, mapbilnr, 'one', atm2ocn_smap)

          call addmrg(fldListTo(complnd)%flds, 'Sa_co2diag', &
               mrg_from1=compatm, mrg_fld1='Sa_co2diag', mrg_type1='copy')
          call addmrg(fldListTo(compocn)%flds, 'Sa_co2diag', &
               mrg_from1=compatm, mrg_fld1='Sa_co2diag', mrg_type1='copy')
       end if

    else if (flds_co2b) then

       ! ---------------------------------------------------------------------
       ! to lnd: prognostic CO2 at the lowest atm model level
       ! ---------------------------------------------------------------------
       if (phase == 'advertise') then
          call addfld(fldListFr(compatm)%flds, 'Sa_co2prog')
          call addfld(fldListTo(complnd)%flds, 'Sa_co2prog')
       else
          call addmap(fldListFr(compatm)%flds, 'Sa_co2prog', complnd, mapbilnr, 'one', atm2lnd_smap)
          call addmrg(fldListTo(complnd)%flds, 'Sa_co2prog', &
               mrg_from1=compatm, mrg_fld1='Sa_co2prog', mrg_type1='copy')
       end if

       ! ---------------------------------------------------------------------
       ! to lnd: diagnostic CO2 at the lowest atm model level
       ! ---------------------------------------------------------------------
       if (phase == 'advertise') then
          call addfld(fldListFr(compatm)%flds, 'Sa_co2diag')
          call addfld(fldListTo(complnd)%flds, 'Sa_co2diag')
       else
          call addmap(fldListFr(compatm)%flds, 'Sa_co2diag', complnd, mapbilnr, 'one', atm2lnd_smap)
          call addmrg(fldListTo(complnd)%flds, 'Sa_co2diag', &
               mrg_from1=compatm, mrg_fld1='Sa_co2diag', mrg_type1='copy')
       end if

       ! ---------------------------------------------------------------------
       ! to atm: surface flux of CO2 from land
       ! ---------------------------------------------------------------------
       if (phase == 'advertise') then
          call addfld(fldListFr(complnd)%flds, 'Fall_fco2_lnd')
          call addfld(fldListTo(compatm)%flds, 'Fall_fco2_lnd')
       else
          call addmap(fldListFr(complnd)%flds, 'Fall_fco2_lnd', compatm, mapconsf, 'one', atm2lnd_smap)
          call addmrg(fldListTo(compatm)%flds, 'Fall_fco2_lnd', &
               mrg_from1=complnd, mrg_fld1='Fall_fco2_lnd', mrg_type1='copy_with_weights', mrg_fracname1='lfrac')
       end if

    else if (flds_co2c) then

       ! ---------------------------------------------------------------------
       ! to lnd and ocn: prognostic CO2 at the lowest atm model level
       ! ---------------------------------------------------------------------
       if (phase == 'advertise') then
          call addfld(fldListFr(compatm)%flds, 'Sa_co2prog')
          call addfld(fldListTo(complnd)%flds, 'Sa_co2prog')
          call addfld(fldListTo(compocn)%flds, 'Sa_co2prog')
       else
          call addmap(fldListFr(compatm)%flds, 'Sa_co2prog', complnd, mapbilnr, 'one', atm2lnd_smap)
          call addmap(fldListFr(compatm)%flds, 'Sa_co2prog', compocn, mapbilnr, 'one', atm2ocn_smap)

          call addmrg(fldListTo(complnd)%flds, 'Sa_co2prog', &
               mrg_from1=compatm, mrg_fld1='Sa_co2prog', mrg_type1='copy')
          call addmrg(fldListTo(compocn)%flds, 'Sa_co2prog', &
               mrg_from1=compatm, mrg_fld1='Sa_co2prog', mrg_type1='copy')
       end if

       ! ---------------------------------------------------------------------
       ! to lnd and ocn: diagnostic CO2 at the lowest atm model level
       ! ---------------------------------------------------------------------
       if (phase == 'advertise') then
          call addfld(fldListFr(compatm)%flds, 'Sa_co2diag')
          call addfld(fldListTo(complnd)%flds, 'Sa_co2diag')
          call addfld(fldListTo(compocn)%flds, 'Sa_co2diag')
       else
          call addmap(fldListFr(compatm)%flds, 'Sa_co2diag', complnd, mapbilnr, 'one', atm2lnd_smap)
          call addmap(fldListFr(compatm)%flds, 'Sa_co2diag', compocn, mapbilnr, 'one', atm2ocn_smap)

          call addmrg(fldListTo(complnd)%flds, 'Sa_co2diag', &
               mrg_from1=compatm, mrg_fld1='Sa_co2diag', mrg_type1='copy')
          call addmrg(fldListTo(compocn)%flds, 'Sa_co2diag', &
               mrg_from1=compatm, mrg_fld1='Sa_co2diag', mrg_type1='copy')
       end if

       ! ---------------------------------------------------------------------
       ! to atm: surface flux of CO2 from land
       ! ---------------------------------------------------------------------
       if (phase == 'advertise') then
          call addfld(fldListFr(complnd)%flds, 'Fall_fco2_lnd')
          call addfld(fldListTo(compatm)%flds, 'Fall_fco2_lnd')
       else
          call addmap(fldListFr(complnd)%flds, 'Fall_fco2_lnd', compatm, mapconsf, 'one', atm2lnd_smap)
          call addmrg(fldListTo(compatm)%flds, 'Fall_fco2_lnd', &
               mrg_from1=complnd, mrg_fld1='Fall_fco2_lnd', mrg_type1='copy_with_weights', mrg_fracname1='lfrac')
       end if

       ! ---------------------------------------------------------------------
       ! to atm: surface flux of CO2 from ocn
       ! ---------------------------------------------------------------------
       if (phase == 'advertise') then
          call addfld(fldListFr(compocn)%flds, 'Faoo_fco2_ocn')
          call addfld(fldListTo(compatm)%flds, 'Faoo_fco2_ocn')
       else
          call addmap(fldListFr(complnd)%flds, 'Faoo_fco2_ocn', compatm, mapconsf, 'one', atm2lnd_smap)
          ! custom merge in med_phases_prep_atm
       end if
    endif

    !-----------------------------
    ! water isotope fields - TODO: add these to dictionary first
    !-----------------------------
    !   'Ratio of ocean surface level abund. H2_16O/H2O/Rstd'
    !    call fld_add(flds_o2x, "So_roce_16O")
    !    call fld_add(flds_x2i, "So_roce_16O")
    !    'Ratio of ocean surface level abund. HDO/H2O/Rstd'
    !    call fld_add(flds_o2x, "So_roce_HDO")
    !    call fld_add(flds_x2i, "So_roce_HDO")

    !-----------------------------------------------------------------------------
    ! CARMA fields (volumetric soil water)
    !-----------------------------------------------------------------------------
    ! TODO: add this
    ! if (carma_flds /= ' ') then
    !    do n = 1,shr_string_listGetNum(carma_flds)
    !       call addfld(fldListFr(complnd)%flds, trim(fldname))
    !       call addmap(fldListFr(complnd)%flds, trim(fldname), compatm, mapconsf, 'one',lnd2atm_smap)
    !       call addfld(fldListTo(compatm)%flds, trim(fldname), mrg_from1=complnd, mrg_fld1=trim(fldname), mrg_type1='copy')
    !    enddo
    ! endif


  end subroutine esmFldsExchange

end module esmFldsExchange_mod
