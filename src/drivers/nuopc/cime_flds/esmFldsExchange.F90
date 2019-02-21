module esmFldsExchange_mod

  !---------------------------------------------------------------------
  ! This is a mediator specific routine that determines ALL possible
  ! fields exchanged between components and their associated routing,
  ! mapping and merging
  !---------------------------------------------------------------------

  implicit none
  public

  public :: esmFldsExchange

  character(*), parameter :: u_FILE_u = &
       __FILE__

!================================================================================
contains
!================================================================================

  subroutine esmFldsExchange(gcomp, phase, rc)

    use ESMF
    use NUOPC
    use med_constants_mod     , only : CX, CS, CL
    use shr_nuopc_scalars_mod , only : flds_scalar_name
    use shr_nuopc_methods_mod , only : chkerr => shr_nuopc_methods_chkerr
    use shr_nuopc_methods_mod , only : fldchk => shr_nuopc_methods_FB_FldChk
    use med_internalstate_mod , only : InternalState
    use glc_elevclass_mod     , only : glc_elevclass_as_string
    use shr_sys_mod           , only : shr_sys_abort
    use esmFlds               , only : shr_nuopc_fldList_type
    use esmFlds               , only : addfld => shr_nuopc_fldList_AddFld
    use esmFlds               , only : addmap => shr_nuopc_fldList_AddMap
    use esmFlds               , only : addmrg => shr_nuopc_fldList_AddMrg
    use esmflds               , only : compmed, compatm, complnd, compocn
    use esmflds               , only : compice, comprof, compwav, compglc, ncomps
    use esmflds               , only : mapbilnr, mapconsf, mapconsd, mappatch
    use esmflds               , only : mapfcopy, mapfiler, mapnstod, mapnstod_consd, mapnstod_consf
    use esmflds               , only : fldListTo, fldListFr, fldListMed_aoflux, fldListMed_ocnalb
    use esmFlds               , only : coupling_mode

    ! input/output parameters:
    type(ESMF_GridComp)              :: gcomp
    character(len=*) , intent(in)    :: phase
    integer          , intent(inout) :: rc

    ! local variables:
    type(InternalState) :: is_local
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
    character(len=64), allocatable :: flds(:)
    character(len=64), allocatable :: suffix(:)
    character(len=*), parameter    :: subname='(esmFldsExchange)'
    !--------------------------------------

    rc = ESMF_SUCCESS

    !---------------------------------------
    ! Get the internal state
    !---------------------------------------

    if (phase /= 'advertise') then
       nullify(is_local%wrap)
       call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    !----------------------------------------------------------
    ! Determine supported coupling model
    !----------------------------------------------------------

    if (phase /= 'advertise') then

       ! CESM Default settings
       coupling_mode = 'cesm'

       if ( fldchk(is_local%wrap%FBexp(compatm)        , 'Faxx_taux', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compice,compice), 'Faii_taux', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_taux', rc=rc)) then

          ! NEMS orig
          ! atm receives merged atm/ocn fluxes computed in mediator and
          ! atm/ice and fluxes computed in ice. The atm/ocn fluxes are
          ! only used for gridcells that If no interpolated values can be
          ! obtained over ocn/ice gridcells on the atm grid (using
          ! bilinear or conservative methods), the interpolated values
          ! from the nearest neighbor method will be used.

          coupling_mode = 'nems_orig'

       else if ( fldchk(is_local%wrap%FBexp(compatm)        , 'Faii_taux', rc=rc) .and. &
                 fldchk(is_local%wrap%FBImp(compice,compice), 'Faii_taux', rc=rc) .and. &
                 fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_taux', rc=rc)) then

          ! NEMS frac
          ! atm receives atm/ice fluxes computed in the ice component and returns merged
          ! lnd/ice/ocn surface fluxes and states back to the mediator

          coupling_mode = 'nems_frac'

       end if
    end if

    !--------------------------------------
    ! Merging arguments:
    ! mrg_fromN = source component index that for the field to be merged
    ! mrg_fldN  = souce field name to be merged
    ! mrg_typeN = merge type ('copy', 'copy_with_weights', 'sum', 'sum_with_weights', 'merge')
    ! NOTE:
    ! mrg_from(compmed) can either be for mediator computed fields for atm/ocn fluxes or for ocn albedos
    !
    ! NOTE:
    ! FBMed_aoflux_o only refer to output fields to the atm/ocn that computed in the
    ! atm/ocn flux calculations. Input fields required from either the atm or the ocn for
    ! these computation will use the logical 'use_med_aoflux' below. This is used to determine
    ! mappings between the atm and ocn needed for these computations.
    !--------------------------------------

    !---------------------------
    ! For now hardwire these
    !---------------------------

    ! these must be less than or equal to the values in fd.yaml
    max_megan = 20
    max_ddep  = 80
    max_fire  = 10
    glc_nec   = 10
    flds_i2o_per_cat = .true.

    iso(1) = ''
    iso(2) = '_16O'
    iso(3) = '_18O'
    iso(4) = '_HDO'

    !----------------------------------------------------------
    ! Initialize mapping file names
    !----------------------------------------------------------

    ! to atm

    call NUOPC_CompAttributeGet(gcomp, name='ice2atm_fmapname', value=ice2atm_fmap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('ice2atm_fmapname = '// trim(ice2atm_fmap), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='ice2atm_smapname', value=ice2atm_smap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('ice2atm_smapname = '// trim(ice2atm_smap), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='lnd2atm_fmapname', value=lnd2atm_fmap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('lnd2atm_fmapname = '// trim(lnd2atm_fmap), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='ocn2atm_smapname', value=ocn2atm_smap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('ocn2atm_smapname = '// trim(ocn2atm_smap), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='ocn2atm_fmapname', value=ocn2atm_fmap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('ocn2atm_fmapname = '// trim(ocn2atm_fmap), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='lnd2atm_smapname', value=lnd2atm_smap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('lnd2atm_smapname = '// trim(lnd2atm_smap), ESMF_LOGMSG_INFO, rc=dbrc)

    ! to lnd

    call NUOPC_CompAttributeGet(gcomp, name='atm2lnd_fmapname', value=atm2lnd_fmap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('atm2lnd_fmapname = '// trim(atm2lnd_fmap), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='atm2lnd_smapname', value=atm2lnd_smap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('atm2lnd_smapname = '// trim(atm2lnd_smap), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='rof2lnd_fmapname', value=rof2lnd_fmap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('rof2lnd_fmapname = '// trim(rof2lnd_fmap), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='glc2lnd_fmapname', value=glc2lnd_fmap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('glc2lnd_smapname = '// trim(glc2lnd_fmap), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='glc2lnd_smapname', value=glc2lnd_smap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('glc2lnd_smapname = '// trim(glc2lnd_smap), ESMF_LOGMSG_INFO, rc=dbrc)

    ! to ice

    call NUOPC_CompAttributeGet(gcomp, name='atm2ice_fmapname', value=atm2ice_fmap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('atm2ice_fmapname = '// trim(atm2ice_fmap), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='atm2ice_smapname', value=atm2ice_smap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('atm2ice_smapname = '// trim(atm2ice_smap), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='atm2ice_vmapname', value=atm2ice_vmap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('atm2ice_vmapname = '// trim(atm2ice_vmap), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='glc2ice_rmapname', value=glc2ice_rmap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('glc2ice_rmapname = '// trim(glc2ice_rmap), ESMF_LOGMSG_INFO, rc=dbrc)

    ! to ocn

    call NUOPC_CompAttributeGet(gcomp, name='atm2ocn_fmapname', value=atm2ocn_fmap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('atm2ocn_fmapname = '// trim(atm2ocn_fmap), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='atm2ocn_smapname', value=atm2ocn_smap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('atm2ocn_smapname = '// trim(atm2ocn_smap), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='atm2ocn_vmapname', value=atm2ocn_vmap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('atm2ocn_vmapname = '// trim(atm2ocn_vmap), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='glc2ocn_liq_rmapname', value=glc2ocn_liq_rmap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('glc2ocn_liq_rmapname = '// trim(glc2ocn_liq_rmap), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='glc2ocn_ice_rmapname', value=glc2ocn_ice_rmap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('glc2ocn_ice_rmapname = '// trim(glc2ocn_ice_rmap), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='wav2ocn_smapname', value=wav2ocn_smap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('wav2ocn_smapname = '// trim(wav2ocn_smap), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='rof2ocn_fmapname', value=rof2ocn_fmap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('rof2ocn_fmapname = '// trim(rof2ocn_fmap), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='rof2ocn_liq_rmapname', value=rof2ocn_liq_rmap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('rof2ocn_liq_rmapname = '// trim(rof2ocn_liq_rmap), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='rof2ocn_ice_rmapname', value=rof2ocn_ice_rmap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('rof2ocn_ice_rmapname = '// trim(rof2ocn_ice_rmap), ESMF_LOGMSG_INFO, rc=dbrc)

    ! to rof

    call NUOPC_CompAttributeGet(gcomp, name='lnd2rof_fmapname', value=lnd2rof_fmap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('lnd2rof_fmapname = '// trim(lnd2rof_fmap), ESMF_LOGMSG_INFO, rc=dbrc)

    ! to glc

    call NUOPC_CompAttributeGet(gcomp, name='lnd2glc_fmapname', value=lnd2glc_fmap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('lnd2glc_fmapname = '// trim(lnd2glc_fmap), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='lnd2glc_smapname', value=lnd2glc_smap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('lnd2glc_smapname = '// trim(lnd2glc_smap), ESMF_LOGMSG_INFO, rc=dbrc)

    ! to wav

    call NUOPC_CompAttributeGet(gcomp, name='atm2wav_smapname', value=atm2wav_smap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('atm2wav_smapname = '// trim(atm2wav_smap), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='ice2wav_smapname', value=ice2wav_smap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('ice2wav_smapname = '// trim(ice2wav_smap), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='ocn2wav_smapname', value=ocn2wav_smap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('ocn2wav_smapname = '// trim(ocn2wav_smap), ESMF_LOGMSG_INFO, rc=dbrc)

    !=====================================================================
    ! scalar information
    !=====================================================================

    if (phase == 'advertise') then
       do n = 1,ncomps
          call addfld(fldListFr(n)%flds, trim(flds_scalar_name))
          call addfld(fldListTo(n)%flds, trim(flds_scalar_name))
       end do
    end if

    !=====================================================================
    ! FIELDS TO MEDIATOR component (for fractions and atm/ocn flux calculation)
    !=====================================================================

    !----------------------------------------------------------
    ! to med: masks from components
    !----------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(complnd)%flds, 'Sl_lfrin')
       call addfld(fldListFr(compocn)%flds, 'So_omask')
       call addfld(fldListFr(compice)%flds, 'Si_imask')
    else
       call addmap(fldListFr(compocn)%flds, 'So_omask', compice,  mapfcopy, 'unset', 'unset')
    end if

    ! ---------------------------------------------------------------------
    ! to med: atm and ocn fields required for atm/ocn flux calculation'
    ! ---------------------------------------------------------------------
    if (phase /= 'advertise') then
       if (trim(coupling_mode) == 'cesm' .or. trim(coupling_mode) == 'nems_orig') then

          call addfld(fldListFr(compatm)%flds, 'Sa_u')
          call addmap(fldListFr(compatm)%flds, 'Sa_u'   , compocn, mappatch, 'one', atm2ocn_vmap)

          call addfld(fldListFr(compatm)%flds, 'Sa_v')
          call addmap(fldListFr(compatm)%flds, 'Sa_v'   , compocn, mappatch, 'one', atm2ocn_vmap)

          call addfld(fldListFr(compatm)%flds, 'Sa_z')
          call addmap(fldListFr(compatm)%flds, 'Sa_z'   , compocn, mapconsf, 'one', atm2ocn_smap)

          call addfld(fldListFr(compatm)%flds, 'Sa_tbot')
          call addmap(fldListFr(compatm)%flds, 'Sa_tbot', compocn, mapconsf, 'one', atm2ocn_smap)

          call addfld(fldListFr(compatm)%flds, 'Sa_pbot')
          call addmap(fldListFr(compatm)%flds, 'Sa_pbot', compocn, mapconsf, 'one', atm2ocn_smap)

          do n = 1,size(iso)
             call addfld(fldListFr(compatm)%flds, 'Sa_shum'//iso(n))
             call addmap(fldListFr(compatm)%flds, 'Sa_shum'//iso(n), compocn, mapconsf, 'one', atm2ocn_smap)
          end do

          if (fldchk(is_local%wrap%FBImp(compatm,compatm), 'Sa_ptem', rc=rc)) then
             call addfld(fldListFr(compatm)%flds, 'Sa_ptem')
             call addmap(fldListFr(compatm)%flds, 'Sa_ptem', compocn, mapconsf, 'one', atm2ocn_smap)
          end if

          if (fldchk(is_local%wrap%FBImp(compatm,compatm), 'Sa_dens', rc=rc)) then
             call addfld(fldListFr(compatm)%flds, 'Sa_dens')
             call addmap(fldListFr(compatm)%flds, 'Sa_dens', compocn, mapconsf, 'one', atm2ocn_smap)
          end if
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to med: swnet fluxes used for budget calculation
    ! ---------------------------------------------------------------------
    ! TODO (mvertens, 2019-01-11): budget implemention needs to be done in CMEPS
    if (phase == 'advertise') then
       call addfld(fldListFr(complnd)%flds, 'Fall_swnet')
       call addfld(fldListFr(compice)%flds, 'Faii_swnet')
       call addfld(fldListFr(compatm)%flds, 'Faxa_swnet')
    else
       if (fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_swnet', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_swnet', compice, mapconsf, 'one'  , atm2ice_fmap)
          call addmap(fldListFr(compatm)%flds, 'Faxa_swnet', compocn, mapconsf, 'one'  , atm2ocn_fmap)
       end if
       if (fldchk(is_local%wrap%FBImp(compice,compice), 'Faii_swnet', rc=rc)) then
          call addmap(fldListFr(compice)%flds, 'Faii_swnet', compocn, mapfcopy, 'unset', 'unset')
       end if
    end if

    !=====================================================================
    ! FIELDS TO LAND
    !=====================================================================

    ! ---------------------------------------------------------------------
    ! from atm:
    ! to lnd: height at the lowest model level from atm
    ! to lnd: surface height from atm
    ! to lnd: zonal wind at the lowest model level from atm
    ! to lnd: meridional wind at the lowest model level from atm
    ! to lnd: Temperature at the lowest model level from atm
    ! to lnd: potential temperature at the lowest model level from atm
    ! to lnd: Pressure at the lowest model level from atm
    ! to lnd: specific humidity at the lowest model level from atm
    ! ---------------------------------------------------------------------

    allocate(flds(11))
    flds = (/'Sa_z', 'Sa_topo', 'Sa_u', 'Sa_v', 'Sa_tbot', 'Sa_ptem', &
             'Sa_pbot', 'Sa_shum', 'Sa_shum_16O', 'Sa_shum_18O', 'Sa_shum_HDO'/)

    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          call addfld(fldListFr(compatm)%flds, trim(fldname))
          call addfld(fldListTo(complnd)%flds, trim(fldname))
       else
          if ( fldchk(is_local%wrap%FBexp(complnd)         , trim(fldname), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm ), trim(fldname), rc=rc)) then
             call addmap(fldListFr(compatm)%flds, trim(fldname), complnd, mapbilnr, 'one', atm2lnd_smap)
             call addmrg(fldListTo(complnd)%flds, trim(fldname), &
                  mrg_from1=compatm, mrg_fld1=trim(fldname), mrg_type1='copy')
          end if
       end if
    end do
    deallocate(flds)

    ! ---------------------------------------------------------------------
    ! to lnd: convective and large scale precipitation rate water equivalent from atm
    ! to lnd: convective and large-scale (stable) snow rate from atm
    ! to lnd: downward longwave heat flux from atm
    ! to lnd: downward direct near-infrared incident solar radiation  from atm
    ! to lnd: downward direct visible incident solar radiation        from atm
    ! to lnd: downward diffuse near-infrared incident solar radiation from atm
    ! to lnd: downward Diffuse visible incident solar radiation       from atm
    ! to lnd: hydrophylic black carbon dry deposition flux   from atm
    ! to lnd: hydrophobic black carbon dry deposition flux   from atm
    ! to lnd: hydrophylic black carbon wet deposition flux   from atm
    ! to lnd: hydrophylic organic carbon dry deposition flux from atm
    ! to lnd: hydrophobic organic carbon dry deposition flux from atm
    ! to lnd: hydrophylic organic carbon wet deposition flux from atm
    ! to lnd: dust wet deposition flux (size 1) from atm
    ! to lnd: dust wet deposition flux (size 2) from atm
    ! to lnd: dust wet deposition flux (size 3) from atm
    ! to lnd: dust wet deposition flux (size 4) from atm
    ! to lnd: dust dry deposition flux (size 1) from atm
    ! to lnd: dust dry deposition flux (size 2) from atm
    ! to lnd: dust dry deposition flux (size 3) from atm
    ! to lnd: dust dry deposition flux (size 4) from atm
    ! to lnd: nitrogen deposition fields from atm
    ! ---------------------------------------------------------------------

    ! TODO (mvertens, 2019-12-13): the nitrogen deposition fluxes here
    ! are not treated the same was as in cesm2.0 release

    allocate(flds(25))
    flds = (/'Faxa_rainc'    , 'Faxa_rainl'   , 'Faxa_snowc'   , 'Faxa_snowl',               &
             'Faxa_lwdn'     , 'Faxa_swndr'   , 'Faxa_swvdr'   , 'Faxa_swndf', 'Faxa_swvdf', &
             'Faxa_bcphidry' , 'Faxa_bcphodry', 'Faxa_bcphiwet',                             &
             'Faxa_ocphidry' , 'Faxa_ocphodry', 'Faxa_ocphiwet',                             &
             'Faxa_dstwet1'  , 'Faxa_dstwet2' , 'Faxa_dstwet3' , 'Faxa_dstwet4',             &
             'Faxa_dstdry1'  , 'Faxa_dstdry2' , 'Faxa_dstdry3' , 'Faxa_dstdry4',             &
             'Faxa_noy'      , 'Faxa_nhx'/)

    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          call addfld(fldListFr(compatm)%flds, trim(fldname))
          call addfld(fldListTo(complnd)%flds, trim(fldname))
       else
          if ( fldchk(is_local%wrap%FBexp(complnd)         , trim(fldname), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm ), trim(fldname), rc=rc)) then
             call addmap(fldListFr(compatm)%flds, trim(fldname), complnd, mapconsf, 'one', atm2lnd_fmap)
             call addmrg(fldListTo(complnd)%flds, trim(fldname), &
                  mrg_from1=compatm, mrg_fld1=trim(fldname), mrg_type1='copy')
          end if
       end if
    end do
    deallocate(flds)

    ! ---------------------------------------------------------------------
    ! to lnd: river channel total water volume from rof
    ! to lnd: river channel main channel water volume from rof
    ! to lnd: river water flux back to land due to flooding
    ! ---------------------------------------------------------------------
    allocate(flds(12))
    flds = (/'Flrr_volr'   , 'Flrr_volr_16O'   , 'Flrr_volr_18O'   , 'Flrr_volr_HDO'   , &
             'Flrr_volrmch', 'Flrr_volrmch_16O', 'Flrr_volrmch_18O', 'Flrr_volrmch_HDO', &
             'Flrr_flood'  , 'Flrr_flood_16O'  , 'Flrr_flood_18O'  , 'Flrr_flood_HDO'  /)

    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          call addfld(fldListFr(comprof)%flds, trim(fldname))
          call addfld(fldListTo(complnd)%flds, trim(fldname))
       else
          if ( fldchk(is_local%wrap%FBExp(complnd)         , trim(fldname), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(comprof, comprof), trim(fldname), rc=rc)) then
             call addmap(fldListFr(comprof)%flds, trim(fldname), complnd, mapconsf, 'one', rof2lnd_fmap)
             call addmrg(fldListTo(complnd)%flds, trim(fldname), &
                  mrg_from1=comprof, mrg_fld1=trim(fldname), mrg_type1='copy')
          end if
       end if
    end do
    deallocate(flds)

    ! ---------------------------------------------------------------------
    ! to lnd: ice sheet grid coverage on global grid from glc
    ! to lnd: ice sheet mask where we are potentially sending non-zero fluxes from glc
    ! to lnd: fields with multiple elevation classes from glc
    ! ---------------------------------------------------------------------
    allocate(flds(2))
    flds = (/'Sg_icemask', 'Sg_icemask_coupled_fluxes'/)

    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          call addfld(fldListFr(compglc)%flds   , trim(fldname))
          call addfld(fldListTo(complnd)%flds   , trim(fldname))
       else
          if ( fldchk(is_local%wrap%FBExp(complnd)         , trim(fldname), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compglc, compglc), trim(fldname), rc=rc)) then
             call addmap(fldListFr(compglc)%flds, trim(fldname), complnd,  mapconsf, 'one', glc2lnd_smap)
             call addmrg(fldListTo(complnd)%flds, trim(fldname), &
                  mrg_from1=compglc, mrg_fld1=trim(fldname), mrg_type1='copy')
          end if
       end if
    end do
    deallocate(flds)

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
             call addfld(fldListTo(complnd)%flds, 'Sg_ice_covered'//trim(cnum))
             call addfld(fldListTo(complnd)%flds, 'Sg_topo'//trim(cnum))
             call addfld(fldListTo(complnd)%flds, 'Flgg_hflx'//trim(cnum))
          end do
       else
          do num = 0, glc_nec
             cnum = glc_elevclass_as_string(num)
             if ( fldchk(is_local%wrap%FBExp(complnd)         , 'Sg_ice_covered'//trim(cnum), rc=rc) .and. &
                  fldchk(is_local%wrap%FBExp(complnd)         , 'Sg_topo'//trim(cnum)       , rc=rc) .and. &
                  fldchk(is_local%wrap%FBExp(complnd)         , 'Flgg_hflx'//trim(cnum)     , rc=rc) .and. &
                  fldchk(is_local%wrap%FBImp(compglc,compglc) , 'Sg_ice_covered'            , rc=rc) .and. &
                  fldchk(is_local%wrap%FBImp(compglc,compglc) , 'Sg_topo'                   , rc=rc) .and. &
                  fldchk(is_local%wrap%FBImp(compglc,compglc) , 'Flgg_hflx'                 , rc=rc)) then
                if (num == 0) then
                   call addmap(FldListFr(compglc)%flds, 'Sg_ice_covered' , complnd, mapconsf, 'unset' , glc2lnd_fmap)
                   call addmap(FldListFr(compglc)%flds, 'Sg_topo'        , compglc, mapconsf, 'custom', glc2lnd_fmap)
                   call addmap(FldListFr(compglc)%flds, 'Flgg_hflx'      , compglc, mapconsf, 'custom', glc2lnd_fmap)
                end if
                call addmrg(fldListTo(complnd)%flds, 'Sg_ice_covered'//trim(cnum), &
                     mrg_from1=compglc, mrg_fld1='Sg_ice_covered'//trim(cnum), mrg_type1='copy')
                call addmrg(fldListTo(complnd)%flds, 'Sg_topo' //trim(cnum),  &
                     mrg_from1=compglc, mrg_fld1='Sg_topo'//trim(cnum), mrg_type1='copy')
                call addmrg(fldListTo(complnd)%flds, 'Flgg_hflx'//trim(cnum), &
                     mrg_from1=compglc, mrg_fld1='Flgg_hflx'//trim(cnum), mrg_type1='copy')
             end if
          end do
       end if
    end if

    !=====================================================================
    ! FIELDS TO ATMOSPHERE
    !=====================================================================

    !----------------------------------------------------------
    ! to atm: Fractions
    !----------------------------------------------------------
    if (phase == 'advertise') then
       ! the following are computed in med_phases_prep_atm
       call addfld(fldListTo(compatm)%flds, 'Sl_lfrac')
       call addfld(fldListTo(compatm)%flds, 'Si_ifrac')
       call addfld(fldListTo(compatm)%flds, 'So_ofrac')
    end if

    ! ---------------------------------------------------------------------
    ! to atm: merged direct  albedo (visible radiation)
    ! to atm: merged diffuse albedo (visible radiation)
    ! to atm: merged direct  albedo (near-infrared radiation)
    ! to atm: merged diffuse albedo (near-infrared radiation)
    ! ---------------------------------------------------------------------
    allocate(suffix(4))
    suffix = (/'avsdr', 'avsdf', 'anidr', 'anidf'/)

    do n = 1,size(suffix)
       if (phase == 'advertise') then
          call addfld(fldListFr(complnd)%flds, 'Sl_'//trim(suffix(n)))
          call addfld(fldListFr(compice)%flds, 'Si_'//trim(suffix(n)))
          call addfld(fldListMed_ocnalb%flds , 'So_'//trim(suffix(n)))
          call addfld(fldListTo(compatm)%flds, 'Sx_'//trim(suffix(n)))
       else
          ! CESM (cam, non-aqua-planet)
          if ( fldchk(is_local%wrap%FBexp(compatm)        , 'Sx_'//trim(suffix(n)), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(complnd,complnd), 'Sl_'//trim(suffix(n)), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compice,compice), 'Si_'//trim(suffix(n)), rc=rc) .and. &
               fldchk(is_local%wrap%FBMed_ocnalb_a        , 'So_'//trim(suffix(n)), rc=rc)) then
             call addmap(fldListFr(complnd)%flds, 'Sl_'//trim(suffix(n)), compatm, mapconsf, 'lfrin', lnd2atm_smap)
             call addmap(fldListFr(compice)%flds, 'Si_'//trim(suffix(n)), compatm, mapconsf, 'ifrac', ice2atm_smap)
             call addmap(fldListMed_ocnalb%flds , 'So_'//trim(suffix(n)), compatm, mapconsf, 'ofrac', ocn2atm_smap)
             call addmrg(fldListTo(compatm)%flds, 'Sx_'//trim(suffix(n)), &
                  mrg_from1=complnd, mrg_fld1='Sl_'//trim(suffix(n)), mrg_type1='merge', mrg_fracname1='lfrac', &
                  mrg_from2=compice, mrg_fld2='Si_'//trim(suffix(n)), mrg_type2='merge', mrg_fracname2='ifrac', &
                  mrg_from3=compmed, mrg_fld3='So_'//trim(suffix(n)), mrg_type3='merge', mrg_fracname3='ofrac')

          ! CESM (cam, aqua-planet)
          else if (fldchk(is_local%wrap%FBMed_ocnalb_a, 'So_'//trim(suffix(n)), rc=rc) .and. &
                   fldchk(is_local%wrap%FBexp(compatm), 'Sx_'//trim(suffix(n)), rc=rc)) then
             call addmap(fldListMed_ocnalb%flds , 'So_'//trim(suffix(n)), compatm, mapconsf, 'ofrac', ocn2atm_smap)
             call addmrg(fldListTo(compatm)%flds, 'Sx_'//trim(suffix(n)), &
                  mrg_from1=compmed, mrg_fld1='So_'//trim(suffix(n)), mrg_type1='merge', mrg_fracname1='ofrac')
          end if
       end if
    end do
    deallocate(suffix)

    ! ---------------------------------------------------------------------
    ! to atm: merged reference temperature at 2 meters
    ! to atm: merged 10m wind speed
    ! to atm: merged reference specific humidity at 2 meters
    ! to atm: merged reference specific water isoptope humidity at 2 meters
    ! ---------------------------------------------------------------------
    allocate(suffix(6))
    suffix = (/'tref', 'u10', 'qref', 'qref_16O', 'qref_18O', 'qref_HDO'/)

    do n = 1,size(suffix)
       if (phase == 'advertise') then
          call addfld(fldListFr(complnd)%flds , 'Sl_'//trim(suffix(n)))
          call addfld(fldListFr(compice)%flds , 'Si_'//trim(suffix(n)))
          call addfld(fldListMed_aoflux%flds  , 'So_'//trim(suffix(n)))
          call addfld(fldListTo(compatm)%flds , 'Sx_'//trim(suffix(n)))
       else
          ! CESM (cam, non-aqua-planet)
          if ( fldchk(is_local%wrap%FBexp(compatm)         , 'Sx_'//trim(suffix(n)), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(complnd,complnd ), 'Sl_'//trim(suffix(n)), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compice,compice ), 'Si_'//trim(suffix(n)), rc=rc) .and. &
               fldchk(is_local%wrap%FBMed_aoflux_o         , 'So_'//trim(suffix(n)), rc=rc)) then
             call addmap(fldListFr(complnd)%flds , 'Sl_'//trim(suffix(n)), compatm, mapconsf, 'lfrin', lnd2atm_fmap)
             call addmap(fldListFr(compice)%flds , 'Si_'//trim(suffix(n)), compatm, mapconsf, 'ifrac', ice2atm_fmap)
             call addmap(fldListMed_aoflux%flds  , 'So_'//trim(suffix(n)), compocn, mapbilnr, 'one'  , atm2ocn_fmap) ! map atm->ocn
             call addmap(fldListMed_aoflux%flds  , 'So_'//trim(suffix(n)), compatm, mapconsf, 'ofrac', ocn2atm_fmap) ! map ocn->atm
             call addmrg(fldListTo(compatm)%flds , 'Sx_'//trim(suffix(n)), &
                  mrg_from1=complnd, mrg_fld1='Sl_'//trim(suffix(n)), mrg_type1='merge', mrg_fracname1='lfrac', &
                  mrg_from2=compice, mrg_fld2='Si_'//trim(suffix(n)), mrg_type2='merge', mrg_fracname2='ifrac', &
                  mrg_from3=compmed, mrg_fld3='So_'//trim(suffix(n)), mrg_type3='merge', mrg_fracname3='ofrac')

          ! NEMS-orig - merged ocn temp
          else if (fldchk(is_local%wrap%FBexp(compatm)        , 'Sx_'//trim(suffix(n)), rc=rc) .and. &
                   fldchk(is_local%wrap%FBImp(compocn,compocn), 'So_'//trim(suffix(n)), rc=rc) .and. &
                   fldchk(is_local%wrap%FBImp(compice,compice), 'Si_'//trim(suffix(n)), rc=rc)) then
             call addmap(fldListFr(compice)%flds, 'Si_'//trim(suffix(n)), compatm, mapnstod_consf, 'ifrac', ice2atm_fmap)
             call addmap(fldListFr(compocn)%flds, 'So_'//trim(suffix(n)), compatm, mapnstod_consf, 'none' , ocn2atm_fmap)
             call addmrg(fldListTo(compatm)%flds, 'Sx_'//trim(suffix(n)), &
                  mrg_from1=compice, mrg_fld1='Si_'//trim(suffix(n)), mrg_type1='merge', mrg_fracname1='ifrac', &
                  mrg_from2=compocn, mrg_fld2='So_'//trim(suffix(n)), mrg_type2='merge', mrg_fracname2='ofrac')

          ! CESM (cam, aqua-planet)
          else if (fldchk(is_local%wrap%FBMed_aoflux_o, 'So_'//trim(suffix(n)), rc=rc)) then
             call addmap(fldListMed_aoflux%flds , 'So_'//trim(suffix(n)), compatm, mapconsf, 'ofrac', ocn2atm_fmap) ! map ocn->atm
             call addmrg(fldListTo(compatm)%flds, 'Sx_'//trim(suffix(n)), &
                  mrg_from1=compmed, mrg_fld1='So_'//trim(suffix(n)), mrg_type1='merge', mrg_fracname1='ofrac')
          end if
       end if
    end do
    deallocate(suffix)

    ! ---------------------------------------------------------------------
    ! to atm: merged zonal surface stress
    ! to atm: merged meridional surface stress
    ! to atm: merged surface latent heat flux
    ! to atm: merged surface sensible heat flux
    ! to atm: merged surface upward longwave heat flux
    ! to atm: evaporation water flux from water
    ! to atm: evaporation water flux from water isotopes
    ! ---------------------------------------------------------------------
    allocate(suffix(9))
    suffix = (/'taux', 'tauy', 'lat', 'sen', 'lwup', 'evap', 'evap_16O', 'evap_18O', 'evap_HDO' /)

    do n = 1,size(suffix)
       if (phase == 'advertise') then
          call addfld(fldListMed_aoflux%flds , 'Faox_'//trim(suffix(n)))
          call addfld(fldListFr(complnd)%flds, 'Fall_'//trim(suffix(n)))
          call addfld(fldListFr(compice)%flds, 'Faii_'//trim(suffix(n)))
          call addfld(fldListTo(compatm)%flds, 'Faii_'//trim(suffix(n))) ! nems-frac
          call addfld(fldListTo(compatm)%flds, 'Faxx_'//trim(suffix(n))) ! cesm, nems-orig
       else
          ! CESM (non aqua-planet)
          if ( fldchk(is_local%wrap%FBImp(complnd,complnd), 'Fall_'//trim(suffix(n)), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compice,compice), 'Faii_'//trim(suffix(n)), rc=rc) .and. &
               fldchk(is_local%wrap%FBMed_aoflux_o        , 'Faox_'//trim(suffix(n)), rc=rc) .and. &
               fldchk(is_local%wrap%FBexp(compatm)        , 'Faxx_'//trim(suffix(n)), rc=rc)) then
             call addmap(fldListMed_aoflux%flds  , 'Faox_'//trim(suffix(n)), compatm, mapconsf, 'ofrac', ocn2atm_fmap)
             call addmap(fldListFr(complnd)%flds , 'Fall_'//trim(suffix(n)), compatm, mapconsf, 'lfrin', lnd2atm_fmap)
             call addmap(fldListFr(compice)%flds , 'Faii_'//trim(suffix(n)), compatm, mapconsf, 'ifrac', ice2atm_fmap)
             call addmrg(fldListTo(compatm)%flds , 'Faxx_'//trim(suffix(n)), &
                  mrg_from1=complnd, mrg_fld1='Fall_'//trim(suffix(n)), mrg_type1='merge', mrg_fracname1='lfrac', &
                  mrg_from2=compice, mrg_fld2='Faii_'//trim(suffix(n)), mrg_type2='merge', mrg_fracname2='ifrac', &
                  mrg_from3=compmed, mrg_fld3='Faox_'//trim(suffix(n)), mrg_type3='merge', mrg_fracname3='ofrac')

          ! NEMS orig (here ofrac = 1.-ifrac)
          else if ( fldchk(is_local%wrap%FBImp(compice,compice), 'Faii_'//trim(suffix(n)), rc=rc) .and. &
                    fldchk(is_local%wrap%FBMed_aoflux_o        , 'Faox_'//trim(suffix(n)), rc=rc) .and. &
                    fldchk(is_local%wrap%FBexp(compatm)        , 'Faxx_'//trim(suffix(n)), rc=rc)) then
             call addmap(fldListMed_aoflux%flds  , 'Faox_'//trim(suffix(n)), compatm, mapnstod_consf, 'none' , ocn2atm_fmap)
             call addmap(fldListFr(compice)%flds , 'Faii_'//trim(suffix(n)), compatm, mapnstod_consf, 'ifrac', ice2atm_fmap)
             call addmrg(fldListTo(compatm)%flds , 'Faxx_'//trim(suffix(n)), &
                  mrg_from1=compice, mrg_fld1='Faii_'//trim(suffix(n)), mrg_type1='merge', mrg_fracname1='ifrac', &
                  mrg_from2=compmed, mrg_fld2='Faox_'//trim(suffix(n)), mrg_type2='merge', mrg_fracname2='ofrac')

          ! NEMS frac (merge done in fv3)
          else if ( fldchk(is_local%wrap%FBImp(compice,compice), 'Faii_'//trim(suffix(n)), rc=rc) .and. &
                    fldchk(is_local%wrap%FBexp(compatm)        , 'Faii_'//trim(suffix(n)), rc=rc)) then
             call addmap(fldListFr(compice)%flds, 'Faii_'//trim(suffix(n)), compatm, mapconsf, 'ifrac', ice2atm_fmap)
             call addmrg(fldListTo(compatm)%flds, 'Faii_'//trim(suffix(n)), &
                  mrg_from1=compice, mrg_fld1='Faii_'//trim(suffix(n)), mrg_type1='merge', mrg_fracname1='ifrac')

          ! CESM (cam, aqua-planet)
          else if (fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_'//trim(suffix(n)), rc=rc) .and. &
                   fldchk(is_local%wrap%FBexp(compatm), 'Faxx_'//trim(suffix(n)), rc=rc)) then
             call addmap(fldListMed_aoflux%flds , 'Faox_'//trim(suffix(n)), compatm, mapconsf, 'ofrac', ocn2atm_fmap)
             call addmrg(fldListTo(compatm)%flds, 'Faxx_'//trim(suffix(n)), &
                  mrg_from1=compmed, mrg_fld1='Faox_'//trim(suffix(n)), mrg_type1='merge', mrg_fracname1='ofrac')
          end if
       end if
    end do
    deallocate(suffix)

    ! ---------------------------------------------------------------------
    ! to atm: merged surface temperature and unmerged temperatures from ice and ocn
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(complnd)%flds, 'Sl_t')
       call addfld(fldListFr(compice)%flds, 'Si_t')
       call addfld(fldListFr(compocn)%flds, 'So_t')

       call addfld(fldListTo(compatm)%flds, 'So_t') ! cesm, nems-frac
       call addfld(fldListTo(compatm)%flds, 'Si_t') ! nems-frac
       call addfld(fldListTo(compatm)%flds, 'Sx_t') ! cesm, nems-orig
    else
       ! CESM - merged ocn/ice/lnd temp and unmerged ocn temp
       if (fldchk(is_local%wrap%FBexp(compatm)        , 'Sx_t', rc=rc) .and. &
           fldchk(is_local%wrap%FBImp(complnd,complnd), 'Sl_t', rc=rc) .and. &
           fldchk(is_local%wrap%FBImp(compice,compice), 'Si_t', rc=rc) .and. &
           fldchk(is_local%wrap%FBImp(compocn,compocn), 'So_t', rc=rc)) then
          call addmap(fldListFr(complnd)%flds, 'Sl_t', compatm, mapconsf , 'lfrin', lnd2atm_fmap)
          call addmap(fldListFr(compice)%flds, 'Si_t', compatm, mapconsf , 'ifrac', ice2atm_fmap)
          call addmap(fldListFr(compocn)%flds, 'So_t', compatm, mapconsf , 'ofrac', ocn2atm_fmap)
          call addmrg(fldListTo(compatm)%flds, 'Sx_t', &
               mrg_from1=complnd, mrg_fld1='Sl_t', mrg_type1='merge', mrg_fracname1='lfrac', &
               mrg_from2=compice, mrg_fld2='Si_t', mrg_type2='merge', mrg_fracname2='ifrac', &
               mrg_from3=compocn, mrg_fld3='So_t', mrg_type3='merge', mrg_fracname3='ofrac')
          call addmrg(fldListTo(compatm)%flds, 'So_t', &
               mrg_from1=compocn, mrg_fld1='So_t', mrg_type1='copy')

       ! NEMS-orig - merged ocn/ice temp and unmerged ocn temp
       else if (fldchk(is_local%wrap%FBexp(compatm)        , 'Sx_t', rc=rc) .and. &
                fldchk(is_local%wrap%FBImp(compocn,compocn), 'So_t', rc=rc) .and. &
                fldchk(is_local%wrap%FBImp(compice,compice), 'Si_t', rc=rc)) then
          call addmap(fldListFr(compice)%flds, 'Si_t', compatm, mapnstod_consf, 'ifrac', ice2atm_fmap)
          call addmap(fldListFr(compocn)%flds, 'So_t', compatm, mapnstod_consf, 'none' , ocn2atm_fmap)
          call addmrg(fldListTo(compatm)%flds, 'Sx_t', &
                mrg_from1=compice, mrg_fld1='Si_t', mrg_type1='merge', mrg_fracname1='ifrac', &
                mrg_from2=compocn, mrg_fld2='So_t', mrg_type2='merge', mrg_fracname2='ofrac')
          call addmrg(fldListTo(compatm)%flds, 'So_t', &
               mrg_from1=compocn, mrg_fld1='So_t', mrg_type1='copy')

       ! CESM aqua-planet - merged and unmerged ocn temp are the same
       else if ( fldchk(is_local%wrap%FBexp(compatm)        , 'Sx_t', rc=rc) .and. &
                 fldchk(is_local%wrap%FBImp(compocn,compocn), 'So_t', rc=rc)) then
          call addmap(fldListFr(compocn)%flds, 'So_t', compatm, mapconsf, 'ofrac', ocn2atm_fmap)
          call addmrg(fldListTo(compatm)%flds, 'Sx_t', &
               mrg_from1=compocn, mrg_fld1='So_t', mrg_type1='merge', mrg_fracname1='ofrac')
          call addmrg(fldListTo(compatm)%flds, 'So_t', &
               mrg_from1=compocn, mrg_fld1='So_t', mrg_type1='copy')
       end if

       ! NEMS-frac - unmerged ice temp
       if ( fldchk(is_local%wrap%FBexp(compatm)        , 'Si_t', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compice,compice), 'Si_t', rc=rc)) then
          call addmap(fldListFr(compice)%flds, 'Si_t', compatm, mapconsf , 'ifrac', ice2atm_fmap)
          call addmrg(fldListTo(compatm)%flds, 'Si_t', &
               mrg_from1=compice, mrg_fld1='Si_t', mrg_type1='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to atm: surface snow depth             from ice
    ! to atm: mean ice volume per unit area  from ice
    ! to atm: mean snow volume per unit area from ice
    ! ---------------------------------------------------------------------
    allocate(flds(3))
    flds = (/'Si_snowh', 'Si_vice', 'Si_vsno'/)

    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          call addfld(fldListFr(compice)%flds, trim(fldname))
          call addfld(fldListTo(compatm)%flds, trim(fldname))
       else
          if ( fldchk(is_local%wrap%FBexp(compatm)        , trim(fldname), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compice,compice), trim(fldname), rc=rc)) then
             if (trim(coupling_mode) == 'nems_orig') then
                call addmap(fldListFr(compice)%flds, trim(fldname), compatm, mapnstod_consf, 'ifrac', ice2atm_fmap)
             else
                call addmap(fldListFr(compice)%flds, trim(fldname), compatm, mapconsf      , 'ifrac', ice2atm_fmap)
             end if
             call addmrg(fldListTo(compatm)%flds, trim(fldname), &
                  mrg_from1=compice, mrg_fld1=trim(fldname), mrg_type1='copy')
          end if
       end if
    end do
    deallocate(flds)

    ! ---------------------------------------------------------------------
    ! to atm: surface saturation specific humidity in ocean from med aoflux
    ! to atm: square of exch. coeff (tracers)               from med aoflux
    ! to atm: surface fraction velocity                     from med aoflux
    ! ---------------------------------------------------------------------
    allocate(suffix(3))
    suffix = (/'ssq', 're', 'ustar'/)

    do n = 1,size(suffix)
       fldname = 'So_'//trim(suffix(n))
       if (phase == 'advertise') then
          call addfld(fldListMed_aoflux%flds  , trim(fldname))
          call addfld(fldListTo(compatm)%flds , trim(fldname))
       else
          if ( fldchk(is_local%wrap%FBexp(compatm) , trim(fldname), rc=rc) .and. &
               fldchk(is_local%wrap%FBMed_aoflux_o , trim(fldname), rc=rc)) then
             call addmap(fldListMed_aoflux%flds    , trim(fldname), compatm, mapconsf, 'ofrac', ocn2atm_fmap) ! map ocn->atm
             call addmrg(fldListTo(compatm)%flds   , trim(fldname), &
                  mrg_from1=compmed, mrg_fld1=trim(fldname), mrg_type1='copy')
          end if
       end if
    end do
    deallocate(suffix)

    ! ---------------------------------------------------------------------
    ! to atm: surface fraction velocity     from land
    ! to atm: aerodynamic resistance        from land
    ! to atm: surface snow water equivalent from land
    ! ---------------------------------------------------------------------
    allocate(suffix(3))
    suffix = (/'fv', 'ram1', 'snowh'/)

    do n = 1,size(suffix)
       fldname = 'Sl_'//trim(suffix(n))
       if (phase == 'advertise') then
          call addfld(fldListFr(complnd)%flds, trim(fldname))
          call addfld(fldListTo(compatm)%flds, trim(fldname))
       else
          if ( fldchk(is_local%wrap%FBexp(compatm)         , trim(fldname), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(complnd,complnd ), trim(fldname), rc=rc)) then
             call addmap(fldListFr(complnd)%flds, trim(fldname), compatm, mapconsf, 'lfrin', lnd2atm_fmap)
             call addmrg(fldListTo(compatm)%flds, trim(fldname), &
                  mrg_from1=complnd, mrg_fld1=trim(fldname), mrg_type1='copy')
          end if
       end if
    end do
    deallocate(suffix)

    ! ---------------------------------------------------------------------
    ! to atm: dust fluxes from land
    ! ---------------------------------------------------------------------
    allocate(suffix(4))
    suffix = (/'flxdst1', 'flxdst2', 'flxdst3', 'flxdst4'/)

    do n = 1,size(suffix)
       fldname = 'Fall_'//trim(suffix(n))
       if (phase == 'advertise') then
          call addfld(fldListFr(complnd)%flds, trim(fldname))
          call addfld(fldListTo(compatm)%flds, trim(fldname))
       else
          if ( fldchk(is_local%wrap%FBImp(complnd, complnd), trim(fldname), rc=rc) .and. &
               fldchk(is_local%wrap%FBExp(compatm)         , trim(fldname), rc=rc)) then
             call addmap(fldListFr(complnd)%flds, trim(fldname), compatm, mapconsf, 'lfrin', lnd2atm_fmap)
             call addmrg(fldListTo(compatm)%flds, trim(fldname), &
                  mrg_from1=complnd, mrg_fld1=trim(fldname), mrg_type1='copy_with_weights', mrg_fracname1='lfrac')
          end if
       end if
    end do
    deallocate(suffix)

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
    ! to atm: dry deposition from land
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
    ! FIELDS TO OCEAN (compocn)
    !=====================================================================

    !----------------------------------------------------------
    ! to ocn: fractional ice coverage wrt ocean from ice
    !----------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compice)%flds, 'Si_ifrac')
       call addfld(fldListTo(compocn)%flds, 'Si_ifrac')
    else
       call addmap(fldListFr(compice)%flds, 'Si_ifrac', compocn,  mapfcopy, 'unset', 'unset')
       call addmrg(fldListTo(compocn)%flds, 'Si_ifrac', mrg_from1=compice, mrg_fld1='Si_ifrac', mrg_type1='copy')
    end if

    ! ---------------------------------------------------------------------
    ! to ocn: downward longwave heat flux from atm
    ! to ocn: downward direct  near-infrared incident solar radiation from atm
    ! to ocn: downward diffuse near-infrared incident solar radiation from atm
    ! to ocn: downward dirrect visible incident solar radiation from atm
    ! to ocn: downward diffuse visible incident solar radiation from atm
    ! ---------------------------------------------------------------------
    allocate(flds(5))
    flds = (/'Faxa_lwdn', 'Faxa_swndr', 'Faxa_swndf', 'Faxa_swvdr', 'Faxa_swndf'/)

    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          call addfld(fldListFr(compatm)%flds, trim(fldname))
          call addfld(fldListTo(compocn)%flds, trim(fldname))
       else
          if ( fldchk(is_local%wrap%FBExp(compocn)        , trim(fldname), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), trim(fldname), rc=rc)) then
             call addmap(fldListFr(compatm)%flds, trim(fldname), compocn, mapconsf, 'one', atm2ocn_fmap)
             call addmrg(fldListTo(compocn)%flds, trim(fldname), mrg_from1=compatm, mrg_fld1=trim(fldname), &
                  mrg_type1='copy_with_weights', mrg_fracname1='ofrac')
          end if
       end if
    end do
    deallocate(flds)

    ! ---------------------------------------------------------------------
    ! to ocn: surface upward longwave heat flux from mediator
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListMed_aoflux%flds  , 'Faox_lwup')
       call addfld(fldListTo(compocn)%flds , 'Foxx_lwup') ! cesm, docn
    else
       if ( fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_lwup', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compocn), 'Foxx_lwup', rc=rc)) then
          call addmrg(fldListTo(compocn)%flds, 'Foxx_lwup', &
               mrg_from1=compmed, mrg_fld1='Faox_lwup', mrg_type1='merge', mrg_fracname1='ofrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ocn: merged longwave net heat flux
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds , 'Faxa_lwdn')
       call addfld(fldListFr(compatm)%flds , 'Faxa_lwnet')
       call addfld(fldListMed_aoflux%flds  , 'Faox_lwup' )
       call addfld(fldListTo(compocn)%flds , 'Foxx_lwnet')
    else
       ! CESM (mom6) (send longwave net to ocn via auto merge)
       if ( fldchk(is_local%wrap%FBExp(compocn)        , 'Foxx_lwnet', rc=rc) .and. &
            fldchk(is_local%wrap%FBMed_aoflux_o        , 'Faox_lwup' , rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_lwdn' , rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_lwdn', compocn, mapconsf, 'one'  , atm2ocn_fmap)
          call addmrg(fldListTo(compocn)%flds, 'Foxx_lwnet', &
               mrg_from1=compmed, mrg_fld1='Faox_lwup', mrg_type1='merge', mrg_fracname1='ofrac', &
               mrg_from2=compatm, mrg_fld2='Faxa_lwdn', mrg_type2='merge', mrg_fracname2='ofrac')

       ! NEMS-orig (mom6) (send longwave net to ocn via custom merge)
       else if ( fldchk(is_local%wrap%FBExp(compocn)        , 'Foxx_lwnet', rc=rc) .and. &
                 fldchk(is_local%wrap%FBMed_aoflux_o        , 'Faox_lwup' , rc=rc) .and. &
                 fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_lwdn' , rc=rc) .and. &
                 fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_lwnet', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_lwdn' , compocn, mapconsf, 'one'  , atm2ocn_fmap)
          call addmap(fldListFr(compatm)%flds, 'Faxa_lwnet', compocn, mapconsf, 'one'  , atm2ocn_fmap)

      ! NEMS-frac (mom6) (send longwave net to ocean via auto merge)
       else if ( fldchk(is_local%wrap%FBExp(compocn)        , 'Foxx_lwnet', rc=rc) .and. &
                 fldchk(is_local%wrap%FBImp(compatm,compatm), 'Foxx_lwnet', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_lwnet', compocn, mapconsf, 'one'  , atm2ocn_fmap)
          call addmrg(fldListTo(compocn)%flds, 'Foxx_lwnet', &
               mrg_from1=compatm, mrg_fld1='Foxx_lwnet', mrg_type1='merge', mrg_fracname1='ofrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ocn: net shortwave radiation from med
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_swvdr')
       call addfld(fldListFr(compatm)%flds, 'Faxa_swndr')
       call addfld(fldListFr(compatm)%flds, 'Faxa_swvdf')
       call addfld(fldListFr(compatm)%flds, 'Faxa_swndf')

       call addfld(fldListFr(compice)%flds, 'Fioi_swpen')
       call addfld(fldListFr(compice)%flds, 'Fioi_swpen_vdr')
       call addfld(fldListFr(compice)%flds, 'Fioi_swpen_vdf')
       call addfld(fldListFr(compice)%flds, 'Fioi_swpen_idr')
       call addfld(fldListFr(compice)%flds, 'Fioi_swpen_idf')

       call addfld(fldListTo(compocn)%flds, 'Foxx_swnet')
       call addfld(fldListTo(compocn)%flds, 'Foxx_swnet_vdr')
       call addfld(fldListTo(compocn)%flds, 'Foxx_swnet_vdf')
       call addfld(fldListTo(compocn)%flds, 'Foxx_swnet_idr')
       call addfld(fldListTo(compocn)%flds, 'Foxx_swnet_idf')
    else
       ! Net shortwave ocean (custom calculation in prep_phases_ocn_mod.F90)
       ! export swpent to ocn without bands
       if ( fldchk(is_local%wrap%FBExp(compocn)        , 'Foxx_swnet', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_swvdr', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_swvdf', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_swndr', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_swndr', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_swvdr', compocn, mapconsf, 'one', atm2ocn_fmap)
          call addmap(fldListFr(compatm)%flds, 'Faxa_swvdf', compocn, mapconsf, 'one', atm2ocn_fmap)
          call addmap(fldListFr(compatm)%flds, 'Faxa_swndr', compocn, mapconsf, 'one', atm2ocn_fmap)
          call addmap(fldListFr(compatm)%flds, 'Faxa_swndf', compocn, mapconsf, 'one', atm2ocn_fmap)

          ! import swpen from ice without bands
          if (fldchk(is_local%wrap%FBImp(compice,compice), 'Fioi_swpen', rc=rc)) then
             call addmap(fldListFr(compice)%flds, 'Fioi_swpen',  compocn, mapfcopy, 'unset', 'unset')
          end if
       end if

       ! export swnet to ocn by bands
       if ( fldchk(is_local%wrap%FBExp(compocn)        , 'Foxx_swnet_vdr', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compocn)        , 'Foxx_swnet_vdf', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compocn)        , 'Foxx_swnet_idr', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compocn)        , 'Foxx_swnet_idf', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_swvdr'    , rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_swvdf'    , rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_swndr'    , rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_swndf'    , rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_swvdr', compocn, mapconsf, 'one', atm2ocn_fmap)
          call addmap(fldListFr(compatm)%flds, 'Faxa_swvdf', compocn, mapconsf, 'one', atm2ocn_fmap)
          call addmap(fldListFr(compatm)%flds, 'Faxa_swndr', compocn, mapconsf, 'one', atm2ocn_fmap)
          call addmap(fldListFr(compatm)%flds, 'Faxa_swndf', compocn, mapconsf, 'one', atm2ocn_fmap)

          ! import swpen from ice by bands
          if ( fldchk(is_local%wrap%FBImp(compice,compice), 'Fioi_swpen_vdr', rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compice,compice), 'Fioi_swpen_vdf', rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compice,compice), 'Fioi_swpen_idr', rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compice,compice), 'Fioi_swpen_idf', rc=rc)) then
             call addmap(fldListFr(compice)%flds, 'Fioi_swpen_vdr', compocn, mapfcopy, 'unset', 'unset')
             call addmap(fldListFr(compice)%flds, 'Fioi_swpen_vdf', compocn, mapfcopy, 'unset', 'unset')
             call addmap(fldListFr(compice)%flds, 'Fioi_swpen_idr', compocn, mapfcopy, 'unset', 'unset')
             call addmap(fldListFr(compice)%flds, 'Fioi_swpen_idf', compocn, mapfcopy, 'unset', 'unset')
          end if
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ocn: per ice thickness fraction and sw penetrating into ocean from ice
    ! ---------------------------------------------------------------------
    if (flds_i2o_per_cat) then
       if (phase == 'advertise') then
          ! 'fractional ice coverage wrt ocean for each thickness category '
          call addfld(fldListFr(compice)%flds, 'Si_ifrac_n')
          call addfld(fldListTo(compocn)%flds, 'Si_ifrac_n')

          ! net shortwave radiation penetrating into ocean for each thickness category
          call addfld(fldListFr(compice)%flds, 'Fioi_swpen_ifrac_n')
          call addfld(fldListTo(compocn)%flds, 'Fioi_swpen_ifrac_n')

          ! 'fractional atmosphere coverage wrt ocean' (computed in med_phases_prep_ocn)
          call addfld(fldListTo(compocn)%flds, 'Sf_afrac')
          ! 'fractional atmosphere coverage used in radiation computations wrt ocean' (computed in med_phases_prep_ocn)
          call addfld(fldListTo(compocn)%flds, 'Sf_afracr')
          ! 'net shortwave radiation times atmosphere fraction' (computed in med_phases_prep_ocn)
          call addfld(fldListTo(compocn)%flds, 'Foxx_swnet_afracr')
       else
          call addmap(fldListFr(compice)%flds, 'Si_ifrac_n'        , compocn, mapfcopy, 'unset', 'unset')
          call addmap(fldListFr(compice)%flds, 'Fioi_swpen_ifrac_n', compocn, mapfcopy, 'unset', 'unset')
          ! Note that 'Sf_afrac, 'Sf_afracr' and 'Foxx_swnet_afracr' will have explicit merging in med_phases_prep_ocn
       end if
    end if

    ! ---------------------------------------------------------------------
    !  to ocn: precipitation rate water equivalent from atm
    !  to ocn: snow rate water equivalent from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       do n = 1,size(iso)
          call addfld(fldListFr(compatm)%flds, 'Faxa_rainc'//iso(n))
          call addfld(fldListFr(compatm)%flds, 'Faxa_rainl'//iso(n))
          call addfld(fldListFr(compatm)%flds, 'Faxa_rain' //iso(n))
          call addfld(fldListTo(compocn)%flds, 'Faxa_rain' //iso(n))

          call addfld(fldListFr(compatm)%flds, 'Faxa_snowc'//iso(n))
          call addfld(fldListFr(compatm)%flds, 'Faxa_snowl'//iso(n))
          call addfld(fldListFr(compatm)%flds, 'Faxa_snow' //iso(n))
          call addfld(fldListTo(compocn)%flds, 'Faxa_snow' //iso(n))
       end do
    else
       do n = 1,size(iso)
          ! Note that the mediator atm/ocn flux calculation needs Faxa_rainc for the gustiness parameterization
          ! which by default is not actually used
          if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rainl'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rainc'//iso(n), rc=rc) .and. &
              (fldchk(is_local%wrap%FBExp(compocn)        , 'Faxa_rain' //iso(n), rc=rc) &
               .or. trim(coupling_mode) == 'cesm' .or. trim(coupling_mode) == 'nems_orig')) then
             call addmap(fldListFr(compatm)%flds, 'Faxa_rainl'//iso(n), compocn, mapconsf, 'one', atm2ocn_fmap)
             call addmap(fldListFr(compatm)%flds, 'Faxa_rainc'//iso(n), compocn, mapconsf, 'one', atm2ocn_fmap)
             if (iso(n) == ' ') then
                call addmrg(fldListTo(compocn)%flds, 'Faxa_rain'//iso(n) , &
                     mrg_from1=compatm, mrg_fld1='Faxa_rainc:Faxa_rainl', &
                     mrg_type1='sum_with_weights', mrg_fracname1='ofrac')
             else
                call addmrg(fldListTo(compocn)%flds, 'Faxa_rain'//iso(n) , &
                     mrg_from1=compatm, mrg_fld1=trim('Faxa_rainc'//iso(n))//':'//trim('Faxa_rainl'//iso(n)), &
                     mrg_type1='sum_with_weights', mrg_fracname1='ofrac')
             end if
          else if ( fldchk(is_local%wrap%FBExp(compocn)        , 'Faxa_rain'//iso(n), rc=rc) .and. &
                    fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rain'//iso(n), rc=rc)) then
             call addmap(fldListFr(compatm)%flds, 'Faxa_rain'//iso(n), compocn, mapconsf, 'one', atm2ocn_fmap)
             call addmrg(fldListTo(compocn)%flds, 'Faxa_rain'//iso(n), mrg_from1=compatm, mrg_fld1='Faxa_rain'//iso(n), &
                  mrg_type1='copy')
          end if
          if ( fldchk(is_local%wrap%FBExp(compocn)        , 'Faxa_snow' //iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_snowl'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_snowc'//iso(n), rc=rc)) then
             call addmap(fldListFr(compatm)%flds, 'Faxa_snowl'//iso(n), compocn, mapconsf, 'one', atm2ocn_fmap)
             call addmap(fldListFr(compatm)%flds, 'Faxa_snowc'//iso(n), compocn, mapconsf, 'one', atm2ocn_fmap)
             if (iso(n) == ' ') then
                call addmrg(fldListTo(compocn)%flds, 'Faxa_snow' //iso(n) , &
                     mrg_from1=compatm, mrg_fld1='Faxa_snowc:Faxa_snowl', &
                     mrg_type1='sum_with_weights', mrg_fracname1='ofrac')
             else
                call addmrg(fldListTo(compocn)%flds, 'Faxa_snow' //iso(n) , &
                     mrg_from1=compatm, mrg_fld1=trim('Faxa_snowc'//iso(n))//':'//trim('Faxa_snowl'//iso(n)), &
                     mrg_type1='sum_with_weights', mrg_fracname1='ofrac')
             end if
          else if ( fldchk(is_local%wrap%FBExp(compocn)        , 'Faxa_snow'//iso(n), rc=rc) .and. &
                    fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_snow'//iso(n), rc=rc)) then
             call addmap(fldListFr(compatm)%flds, 'Faxa_snow'//iso(n), compocn, mapconsf, 'one', atm2ocn_fmap)
             call addmrg(fldListTo(compocn)%flds, 'Faxa_snow'//iso(n), mrg_from1=compatm, mrg_fld1='Faxa_snow'//iso(n), &
                  mrg_type1='copy')
          end if
       end do
    end if

    ! ---------------------------------------------------------------------
    ! to ocn: merged sensible heat flux
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds , 'Faxa_sen')
       call addfld(fldListMed_aoflux%flds  , 'Faox_sen')
       call addfld(fldListFr(compice)%flds , 'Fioi_melth')
       call addfld(fldListTo(compocn)%flds , 'Foxx_sen')
    else
       ! NEMS orig
       if ( fldchk(is_local%wrap%FBexp(compocn)         , 'Foxx_sen'  , rc=rc) .and. &
            fldchk(is_local%wrap%FBMed_aoflux_o         , 'Faox_sen'  , rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compice, compice), 'Fioi_melth', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm, compatm), 'Faxa_sen'  , rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_sen'  , compocn, mapconsf, 'one'  , atm2ocn_fmap) ! map atm->ocn
          call addmap(fldListFr(compice)%flds, 'Fioi_melth', compocn, mapfcopy, 'unset', 'unset')

       ! NEMS frac
       else if ( fldchk(is_local%wrap%FBexp(compocn)         , 'Foxx_sen', rc=rc) .and. &
                 fldchk(is_local%wrap%FBImp(compatm, compatm), 'Faxa_sen', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_sen', compocn, mapconsf, 'one'  , atm2ocn_fmap) ! map atm->ocn
          call addmrg(fldListTo(compocn)%flds, 'Foxx_sen', &
               mrg_from1=compatm, mrg_fld1='Faxa_sen'  , mrg_type1='merge', mrg_fracname1='ofrac', &
               mrg_from2=compice, mrg_fld2='Fioi_melth', mrg_type2='merge', mrg_fracname2='ifrac')

       ! CESM
       else if ( fldchk(is_local%wrap%FBexp(compocn), 'Foxx_sen', rc=rc) .and. &
                 fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_sen', rc=rc)) then
          call addmrg(fldListTo(compocn)%flds, 'Foxx_sen', &
               mrg_from1=compmed, mrg_fld1='Faox_sen', mrg_type1='merge', mrg_fracname1='ofrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ocn: surface latent heat flux and evaporation water flux
    ! ---------------------------------------------------------------------
    do n = 1,size(iso)
       if (phase == 'advertise') then
          call addfld(fldListMed_aoflux%flds , 'Faxa_lat' //iso(n))
          call addfld(fldListMed_aoflux%flds , 'Faox_lat' //iso(n))
          call addfld(fldListMed_aoflux%flds , 'Faox_evap'//iso(n))
          call addfld(fldListTo(compocn)%flds, 'Foxx_lat' //iso(n))
          call addfld(fldListTo(compocn)%flds, 'Foxx_evap'//iso(n))
       else
          ! CESM
          if ( fldchk(is_local%wrap%FBexp(compocn), 'Foxx_lat'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_lat'//iso(n), rc=rc)) then
             call addmrg(fldListTo(compocn)%flds, 'Foxx_lat'//iso(n), &
                  mrg_from1=compmed, mrg_fld1='Faox_lat'//iso(n), mrg_type1='merge', mrg_fracname1='ofrac')
          end if
          if ( fldchk(is_local%wrap%FBExp(compocn), 'Foxx_evap', rc=rc) .and. &
               fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_evap', rc=rc)) then
             call addmrg(fldListTo(compocn)%flds, 'Foxx_evap', &
                  mrg_from1=compmed, mrg_fld1='Faox_evap', mrg_type1='merge', mrg_fracname1='ofrac')
          end if
          ! NEMS orig
          if ( fldchk(is_local%wrap%FBexp(compocn)         , 'Foxx_lat'  , rc=rc) .and. &
               fldchk(is_local%wrap%FBMed_aoflux_o         , 'Foax_evap' , rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm, compatm), 'Faxa_lat'  , rc=rc)) then
             call addmap(fldListFr(compatm)%flds, 'Faxa_lat', compocn, mapconsf, 'one', atm2ocn_fmap)
          end if

          ! NEMS-frac and NEMS-orig
          ! Foxx_evap is passed to mom6 but but not the latent heat flux and  mom6 then computes
          ! the latent heat flux from the imported evaporative flux. However, the evap passed to mom6
          ! in med_phases_prep_ocn is in fact derived from  the latent heat flux obtained from the atm (fv3).
          ! TODO (mvertens, 2019-10-01): Can we unify this and have MOM6 use latent heat flux?
       end if
    end do

    ! ---------------------------------------------------------------------
    ! to ocn: wind speed squared at 10 meters from med
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListMed_aoflux%flds , 'So_duu10n')
       call addfld(fldListTo(compocn)%flds, 'So_duu10n')
    else
       if ( fldchk(is_local%wrap%FBMed_aoflux_o, 'So_duu10n', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compocn), 'So_duu10n', rc=rc)) then

          call addmap(fldListMed_aoflux%flds , 'So_duu10n', compatm, mapconsf, 'ofrac', ocn2atm_fmap) ! map ocn->atm
          call addmrg(fldListTo(compocn)%flds, 'So_duu10n', &
               mrg_from1=compmed, mrg_fld1='So_duu10n', mrg_type1='copy')
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
    ! to ocn: hydrophylic black carbon dry deposition flux from atm
    ! to ocn: hydrophobic black carbon dry deposition flux from atm
    ! to ocn: hydrophylic black carbon wet deposition flux from atm
    ! to ocn: hydrophylic organic carbon dry deposition flux from atm
    ! to ocn: hydrophobic organic carbon dry deposition flux from atm
    ! to ocn: hydrophylic organic carbon wet deposition flux to ice from atm
    ! to ocn: dust wet deposition flux (size 1) from atm
    ! to ocn: dust wet deposition flux (size 2) from atm
    ! to ocn: dust wet deposition flux (size 3) from atm
    ! to ocn: dust wet deposition flux (size 4) from atm
    ! to ocn: dust dry deposition flux (size 1) from atm
    ! to ocn: dust dry deposition flux (size 2) from atm
    ! to ocn: dust dry deposition flux (size 3) from atm
    ! to ocn: dust dry deposition flux (size 4) from atm
    ! ---------------------------------------------------------------------
    allocate(suffix(14))
    suffix = (/'bcphidry', 'bcphodry', 'bcphiwet',           &
               'ocphidry', 'ocphodry', 'ocphiwet',           &
               'dstwet1' , 'dstwet2' , 'dstwet3', 'dstwet4', &
               'dstdry1' , 'dstdry2' , 'dstdry3', 'dstdry4' /)

    do n = 1,size(suffix)
       fldname = 'Faxa_'//trim(suffix(n))
       if (phase == 'advertise') then
          call addfld(fldListFr(compatm)%flds, trim(fldname))
          call addfld(fldListTo(compocn)%flds, trim(fldname))
       else
          if ( fldchk(is_local%wrap%FBImp(compatm,compatm), trim(fldname), rc=rc) .and. &
               fldchk(is_local%wrap%FBExp(compocn)        , trim(fldname), rc=rc)) then
             call addmap(fldListFr(compatm)%flds, trim(fldname), compocn, mapconsf, 'one', atm2ocn_fmap)
             call addmrg(fldListTo(compocn)%flds, trim(fldname), &
                  mrg_from1=compatm, mrg_fld1=trim(fldname), mrg_type1='copy_with_weights', mrg_fracname1='ofrac')
          end if
       end if
    end do
    deallocate(suffix)

    !-----------------------------------------------------------------------------
    ! to ocn: nitrogen deposition fields from atm
    !-----------------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_noy')
       call addfld(fldListFr(compatm)%flds, 'Faxa_nhx')
       call addfld(fldListTo(compocn)%flds, 'Faxa_noy')
       call addfld(fldListTo(compocn)%flds, 'Faxa_nhx')
    else
       if ( fldchk(is_local%wrap%FBImp(compatm, compatm), 'Faxa_noy', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compocn)         , 'Faxa_noy', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_noy', compocn, mapbilnr, 'one', atm2ocn_smap)
          call addmrg(fldListTo(compocn)%flds, 'Faxa_noy', &
               mrg_from1=compatm, mrg_fld1=trim(fldname), mrg_type1='copy_with_weights', mrg_fracname1='ofrac')
       end if
       if ( fldchk(is_local%wrap%FBImp(compatm, compatm), 'Faxa_nhx', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compocn)         , 'Faxa_nhx', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_nhx', compocn, mapbilnr, 'one', atm2ocn_smap)
          call addmrg(fldListTo(compocn)%flds, 'Faxa_nhx', &
               mrg_from1=compatm, mrg_fld1=trim(fldname), mrg_type1='copy_with_weights', mrg_fracname1='ofrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ocn: merge zonal surface stress from ice and (atm or med)
    ! ---------------------------------------------------------------------
    allocate(suffix(2))
    suffix = (/'taux', 'tauy'/)

    do n = 1,size(suffix)
       if (phase == 'advertise') then
          call addfld(fldListMed_aoflux%flds  , 'Faox_'//trim(suffix(n)))
          call addfld(fldListFr(compice)%flds , 'Fioi_'//trim(suffix(n)))
          call addfld(fldListFr(compatm)%flds , 'Faxa_'//trim(suffix(n)))
          call addfld(fldListTo(compocn)%flds , 'Foxx_'//trim(suffix(n)))
       else
          ! NEMS orig and NEMS frac
          if ( fldchk(is_local%wrap%FBexp(compocn)        , 'Foxx_'//trim(suffix(n)), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_'//trim(suffix(n)), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compice,compice), 'Fioi_'//trim(suffix(n)), rc=rc)) then
             call addmap(fldListFr(compatm)%flds, 'Faxa_'//trim(suffix(n)), compocn, mapconsf, 'one'  , atm2ocn_fmap) ! map atm->ocn
             call addmap(fldListFr(compice)%flds, 'Fioi_'//trim(suffix(n)), compocn, mapfcopy, 'unset', 'unset')

             ! NEMS-frac
             call addmrg(fldListTo(compocn)%flds, 'Foxx_'//trim(suffix(n)), &
                  mrg_from1=compatm, mrg_fld1='Faxa_'//trim(suffix(n)), mrg_type1='merge', mrg_fracname1='ofrac', &
                  mrg_from2=compice, mrg_fld2='Fioi_'//trim(suffix(n)), mrg_type2='merge', mrg_fracname2='ifrac')
             ! NEMS-orig
             ! custom merge calculation in med_phases_prep_ocn will be done that will overwrite the auto-merge done above

          ! CESM
          else if (fldchk(is_local%wrap%FBexp(compocn), 'Foxx_'//trim(suffix(n)), rc=rc) .and. &
                   fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_'//trim(suffix(n)), rc=rc)) then
             call addmap(fldListFr(compice)%flds, 'Fioi_'//trim(suffix(n)), compocn, mapfcopy, 'unset', 'unset')
             call addmrg(fldListTo(compocn)%flds, 'Foxx_'//trim(suffix(n)), &
                  mrg_from1=compmed, mrg_fld1='Faox_'//trim(suffix(n)), mrg_type1='merge', mrg_fracname1='ofrac', &
                  mrg_from2=compice, mrg_fld2='Fioi_'//trim(suffix(n)), mrg_type2='merge', mrg_fracname2='ifrac')
          end if
       end if
    end do
    deallocate(suffix)

    ! ---------------------------------------------------------------------
    ! to ocn: water flux due to melting ice from ice
    ! ---------------------------------------------------------------------
    do n = 1,size(iso)
       if (phase == 'advertise') then
          call addfld(fldListFr(compice)%flds , 'Fioi_meltw'//iso(n))
          call addfld(fldListTo(compocn)%flds , 'Fioi_meltw'//iso(n))
       else
          if ( fldchk(is_local%wrap%FBexp(compocn)         , 'Fioi_meltw'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compice, compice), 'Fioi_meltw'//iso(n), rc=rc)) then
             call addmap(fldListFr(compice)%flds, 'Fioi_meltw'//iso(n),    compocn,  mapfcopy, 'unset', 'unset')
             call addmrg(fldListTo(compocn)%flds, 'Fioi_meltw'//iso(n), &
                  mrg_from1=compice, mrg_fld1='Fioi_meltw'//iso(n), mrg_type1='copy_with_weights', mrg_fracname1='ifrac')
          end if
       end if
    end do

    ! ---------------------------------------------------------------------
    ! to ocn: heat flux from melting ice from ice
    ! to ocn: salt flux from ice
    ! to ocn: hydrophylic black carbon deposition flux from ice
    ! to ocn: hydrophobic black carbon deposition flux from ice
    ! to ocn: dust flux from ice
    ! ---------------------------------------------------------------------
    ! TODO (mvertens, 2019-01-07): is fioi_melth being handled here?
    ! Is fd.yaml correctly aliasing Fioi_melth?

    allocate(flds(5))
    flds = (/'Fioi_melth', 'Fioi_salt', 'Fioi_bcphi', 'Fioi_bcpho', 'Fioi_flxdst'/)

    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          call addfld(fldListFr(compice)%flds, trim(fldname))
          call addfld(fldListTo(compocn)%flds, trim(fldname))
       else
          if ( fldchk(is_local%wrap%FBExp(compocn)         , trim(fldname), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compice, compice), trim(fldname), rc=rc)) then
             call addmap(fldListFr(compice)%flds, trim(fldname), compocn,  mapfcopy, 'unset', 'unset')
             call addmrg(fldListTo(compocn)%flds, trim(fldname), &
                  mrg_from1=compice, mrg_fld1=trim(fldname), mrg_type1='copy_with_weights', mrg_fracname1='ifrac')
          end if
       end if
    end do
    deallocate(flds)

    !-----------------------------
    ! to ocn: liquid runoff from rof and glc components
    ! to ocn: frozen runoff flux from rof and glc components
    ! to ocn: waterflux back to ocn due to flooding from rof
    !-----------------------------

    if (phase == 'advertise') then
       do n = 1,size(iso)
          call addfld(fldListFr(compglc)%flds, 'Fogg_rofl'//iso(n))
          call addfld(fldListFr(compglc)%flds, 'Fogg_rofi'//iso(n))
          call addfld(fldListFr(comprof)%flds, 'Forr_rofl'//iso(n))
          call addfld(fldListFr(comprof)%flds, 'Forr_rofi'//iso(n))
          call addfld(fldListTo(compocn)%flds, 'Foxx_rofl'//iso(n))
          call addfld(fldListTo(compocn)%flds, 'Foxx_rofi'//iso(n))
       end do
    else
       do n = 1,size(iso)
          ! from both rof and glc to con
          if ( fldchk(is_local%wrap%FBExp(compocn)         , 'Foxx_rofl'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(comprof, comprof), 'Forr_rofl'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compglc, compglc), 'Fogg_rofl'//iso(n), rc=rc)) then
             call addmap(fldListFr(comprof)%flds, 'Forr_rofl'//iso(n), compocn, mapconsf, 'none', rof2ocn_liq_rmap)
             call addmap(fldListFr(compglc)%flds, 'Fogg_rofl'//iso(n), compocn, mapconsf, 'one' , glc2ocn_liq_rmap)
             call addmrg(fldListTo(compocn)%flds, 'Foxx_rofl'//iso(n), &
                  mrg_from1=comprof, mrg_fld1='Forr_rofl:Flrr_flood', mrg_type1='sum', &
                  mrg_from2=compglc, mrg_fld2='Fogg_rofl'//iso(n)   , mrg_type2='sum')

          ! liquid runoff from rof and flood to ocn
          else if ( fldchk(is_local%wrap%FBExp(compocn)         , 'Foxx_rofl' //iso(n), rc=rc) .and. &
                    fldchk(is_local%wrap%FBImp(comprof, comprof), 'Forr_rofl' //iso(n), rc=rc) .and. &
                    fldchk(is_local%wrap%FBImp(comprof, comprof), 'Flrr_flood'//iso(n), rc=rc)) then
             call addmap(fldListFr(comprof)%flds, 'Flrr_flood'//iso(n), compocn, mapfiler, 'none', rof2ocn_fmap)
             call addmap(fldListFr(comprof)%flds, 'Forr_rofl' //iso(n), compocn, mapfiler, 'none', rof2ocn_liq_rmap)
             call addmrg(fldListTo(compocn)%flds, 'Foxx_rofl' //iso(n), &
                  mrg_from1=comprof, mrg_fld1='Forr_rofl:Flrr_flood', mrg_type1='sum')

          ! liquid from just rof to ocn
          else if ( fldchk(is_local%wrap%FBExp(compocn)         , 'Foxx_rofl'//iso(n), rc=rc) .and. &
                    fldchk(is_local%wrap%FBImp(comprof, comprof), 'Forr_rofl'//iso(n), rc=rc)) then
             call addmap(fldListFr(comprof)%flds, 'Forr_rofl'//iso(n), compocn, mapconsf, 'none', rof2ocn_liq_rmap)
             call addmrg(fldListTo(compocn)%flds, 'Foxx_rofl'//iso(n), &
                  mrg_from1=comprof, mrg_fld1='Forr_rofl:Flrr_flood', mrg_type1='sum')

          ! liquid runoff from just glc to ocn
          else if ( fldchk(is_local%wrap%FBExp(compocn)         , 'Foxx_rofl'//iso(n), rc=rc) .and. &
                    fldchk(is_local%wrap%FBImp(compglc, compglc), 'Fogg_rofl'//iso(n), rc=rc)) then
             call addmap(fldListFr(compglc)%flds, 'Fogg_rofl'//iso(n), compocn,  mapconsf, 'one', glc2ocn_liq_rmap)
             call addmrg(fldListTo(compocn)%flds, 'Foxx_rofl'//iso(n), &
                  mrg_from1=compglc, mrg_fld1='Fogg_rofl'//iso(n), mrg_type1='copy')
          end if

          ! ice runoff from both rof and glc to ocn
          if ( fldchk(is_local%wrap%FBExp(compocn)         , 'Foxx_rofi'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(comprof, comprof), 'Forr_rofi'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compglc, compglc), 'Fogg_rofi'//iso(n), rc=rc)) then
             call addmap(fldListFr(comprof)%flds, 'Forr_rofi'//iso(n), compocn, mapconsf, 'none', rof2ocn_ice_rmap)
             call addmap(fldListFr(compglc)%flds, 'Fogg_rofi'//iso(n), compocn, mapconsf, 'one' , glc2ocn_ice_rmap)
             call addmrg(fldListTo(compocn)%flds, 'Foxx_rofi'//iso(n), &
                  mrg_from1=comprof, mrg_fld1='Forr_rofi'//iso(n), mrg_type1='sum', &
                  mrg_from2=compglc, mrg_fld2='Fogg_rofi'//iso(n), mrg_type2='sum')

          ! ice runoff from just rof to ocn
          else if ( fldchk(is_local%wrap%FBExp(compocn)         , 'Foxx_rofi'//iso(n), rc=rc) .and. &
                    fldchk(is_local%wrap%FBImp(comprof, comprof), 'Forr_rofi'//iso(n), rc=rc)) then
             call addmap(fldListFr(comprof)%flds, 'Forr_rofi'//iso(n), compocn, mapconsf, 'none', rof2ocn_ice_rmap)
             call addmrg(fldListTo(compocn)%flds, 'Foxx_rofi'//iso(n), &
                  mrg_from1=comprof, mrg_fld1='Forr_rofi', mrg_type1='copy')

          ! ice runoff from just glc to ocn
          else if ( fldchk(is_local%wrap%FBExp(compocn)         , 'Foxx_rofi'//iso(n), rc=rc) .and. &
                    fldchk(is_local%wrap%FBImp(compglc, compglc), 'Fogg_rofi'//iso(n), rc=rc)) then
             call addmap(fldListFr(compglc)%flds, 'Fogg_rofi'//iso(n), compocn,  mapconsf, 'one', glc2ocn_ice_rmap)
             call addmrg(fldListTo(compocn)%flds, 'Foxx_rofi'//iso(n), &
                  mrg_from1=compglc, mrg_fld1='Fogg_rofi'//iso(n), mrg_type1='copy')
          end if
       end do
    end if

    !-----------------------------
    ! to ocn: Langmuir multiplier from wave
    ! to ocn: Stokes drift u component from wave
    ! to ocn: Stokes drift v component from wave
    ! to ocn: Stokes drift depth from wave
    !-----------------------------
    allocate(flds(4))
    flds = (/'Sw_lamult', 'Sw_ustokes', 'Sw_vstokes', 'Sw_hstokes'/)

    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          call addfld(fldListFr(compwav)%flds, trim(fldname))
          call addfld(fldListTo(compocn)%flds, trim(fldname))
       else
          if ( fldchk(is_local%wrap%FBExp(compocn)         , trim(fldname), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compwav, compwav), trim(fldname), rc=rc)) then
             call addmap(fldListFr(compwav)%flds, trim(fldname), compocn,  mapbilnr, 'one', wav2ocn_smap)
             call addmrg(fldListTo(compocn)%flds, trim(fldname), &
               mrg_from1=compwav, mrg_fld1=trim(fldname), mrg_type1='copy')
          end if
       end if
    end do
    deallocate(flds)


    !=====================================================================
    ! FIELDS TO ICE (compice)
    !=====================================================================

    ! ---------------------------------------------------------------------
    ! to ice: downward longwave heat flux from atm
    ! to ice: downward direct near-infrared incident solar radiation  from atm
    ! to ice: downward direct visible incident solar radiation        from atm
    ! to ice: downward diffuse near-infrared incident solar radiation from atm
    ! to ice: downward Diffuse visible incident solar radiation       from atm
    ! to ice: hydrophylic black carbon dry deposition flux   from atm
    ! to ice: hydrophobic black carbon dry deposition flux   from atm
    ! to ice: hydrophylic black carbon wet deposition flux   from atm
    ! to ice: hydrophylic organic carbon dry deposition flux from atm
    ! to ice: hydrophobic organic carbon dry deposition flux from atm
    ! to ice: hydrophylic organic carbon wet deposition flux from atm
    ! to ice: dust wet deposition flux (size 1) from atm
    ! to ice: dust wet deposition flux (size 2) from atm
    ! to ice: dust wet deposition flux (size 3) from atm
    ! to ice: dust wet deposition flux (size 4) from atm
    ! to ice: dust dry deposition flux (size 1) from atm
    ! to ice: dust dry deposition flux (size 2) from atm
    ! to ice: dust dry deposition flux (size 3) from atm
    ! to ice: dust dry deposition flux (size 4) from atm
    ! ---------------------------------------------------------------------
    allocate(flds(19))
    flds = (/'Faxa_lwdn'    , 'Faxa_swndr'   , 'Faxa_swvdr'   , 'Faxa_swndf' , 'Faxa_swvdf', &
             'Faxa_bcphidry', 'Faxa_bcphodry', 'Faxa_bcphiwet',                              &
             'Faxa_ocphidry', 'Faxa_ocphodry', 'Faxa_ocphiwet',                              &
             'Faxa_dstwet1' , 'Faxa_dstwet2' , 'Faxa_dstwet3' , 'Faxa_dstwet4',              &
             'Faxa_dstdry1' , 'Faxa_dstdry2' , 'Faxa_dstdry3' , 'Faxa_dstdry4'/)

    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          call addfld(fldListFr(compatm)%flds, trim(fldname))
          call addfld(fldListTo(compice)%flds, trim(fldname))
       else
          if ( fldchk(is_local%wrap%FBExp(compice)        , trim(fldname), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), trim(fldname), rc=rc)) then
             call addmap(fldListFr(compatm)%flds, trim(fldname), compice, mapconsf, 'one', atm2ice_fmap)
             call addmrg(fldListTo(compice)%flds, trim(fldname), &
                  mrg_from1=compatm, mrg_fld1=trim(fldname), mrg_type1='copy')
          end if
       end if
    end do
    deallocate(flds)

    ! ---------------------------------------------------------------------
    ! to ice: convective and large scale precipitation rate water equivalent from atm
    ! to ice: rain and snow rate from atm
    ! ---------------------------------------------------------------------
    do n = 1,size(iso)
       if (phase == 'advertise') then
          call addfld(fldListFr(compatm)%flds, 'Faxa_rainc'//iso(n))
          call addfld(fldListFr(compatm)%flds, 'Faxa_rainl'//iso(n))
          call addfld(fldListFr(compatm)%flds, 'Faxa_rain' //iso(n))
          call addfld(fldListTo(compice)%flds, 'Faxa_rain' //iso(n))

          call addfld(fldListFr(compatm)%flds, 'Faxa_snowc'//iso(n))
          call addfld(fldListFr(compatm)%flds, 'Faxa_snowl'//iso(n))
          call addfld(fldListFr(compatm)%flds, 'Faxa_snow' //iso(n))
          call addfld(fldListTo(compice)%flds, 'Faxa_snow' //iso(n))
       else
          if ( fldchk(is_local%wrap%FBexp(compice)        , 'Faxa_rain' //iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rainl'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rainc'//iso(n), rc=rc)) then
             call addmap(fldListFr(compatm)%flds, 'Faxa_rainc'//iso(n), compice, mapconsf, 'one', atm2ice_fmap)
             call addmap(fldListFr(compatm)%flds, 'Faxa_rainl'//iso(n), compice, mapconsf, 'one', atm2ice_fmap)
             if (iso(n) == ' ') then
                fldname = 'Faxa_rainc:Faxa_rainl'
             else
                fldname = trim('Faxa_rainc'//iso(n))//':'//trim('Faxa_rainl'//iso(n))
             end if
             call addmrg(fldListTo(compice)%flds, 'Faxa_rain' //iso(n) , &
                  mrg_from1=compatm, mrg_fld1=trim(fldname), mrg_type1='sum')
          else if ( fldchk(is_local%wrap%FBexp(compice)        , 'Faxa_rain'//iso(n), rc=rc) .and. &
                    fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rain'//iso(n), rc=rc)) then
             call addmap(fldListFr(compatm)%flds, 'Faxa_rain'//iso(n), compice, mapconsf, 'one', atm2ice_fmap)
             call addmrg(fldListTo(compice)%flds, 'Faxa_rain'//iso(n), &
                  mrg_from1=compatm, mrg_fld1='Faxa_rain'//iso(n), mrg_type1='copy')
          end if

          if ( fldchk(is_local%wrap%FBexp(compice)        , 'Faxa_snow' //iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_snowl'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_snowc'//iso(n), rc=rc)) then
             call addmap(fldListFr(compatm)%flds, 'Faxa_snowc'//iso(n), compice, mapconsf, 'one', atm2ice_fmap)
             call addmap(fldListFr(compatm)%flds, 'Faxa_snowl'//iso(n), compice, mapconsf, 'one', atm2ice_fmap)
             if (iso(n) == ' ') then
                fldname = 'Faxa_snowc:Faxa_snowl'
             else
                fldname = trim('Faxa_snowc'//iso(n))//':'//trim('Faxa_snowl'//iso(n))
             end if
             call addmrg(fldListTo(compice)%flds, 'Faxa_snow' //iso(n) , &
                  mrg_from1=compatm, mrg_fld1=trim(fldname), mrg_type1='sum')
          else if ( fldchk(is_local%wrap%FBexp(compice)        , 'Faxa_snow'//iso(n), rc=rc) .and. &
                    fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_snow'//iso(n), rc=rc)) then
             call addmap(fldListFr(compatm)%flds, 'Faxa_snow'//iso(n), compice, mapconsf, 'one', atm2ice_fmap)
             call addmrg(fldListTo(compice)%flds, 'Faxa_snow'//iso(n), &
                  mrg_from1=compatm, mrg_fld1='Faxa_snow'//iso(n), mrg_type1='copy')
          end if
       end if
    end do

    ! ---------------------------------------------------------------------
    ! to ice: height at the lowest model level from atm
    ! to ice: pressure at the lowest model level fromatm
    ! to ice: temperature at the lowest model level from atm
    ! to ice: potential temperature at the lowest model level from atm
    ! to ice: density at the lowest model level from atm
    ! to ice: zonal wind at the lowest model level from atm
    ! to ice: meridional wind at the lowest model level from atm
    ! to ice: specific humidity at the lowest model level from atm
    ! to ice: specific humidity for water isotopes at the lowest model level from atm
    ! ---------------------------------------------------------------------
    allocate(flds(11))
    flds = (/'Sa_z', 'Sa_pbot', 'Sa_tbot', 'Sa_ptem', 'Sa_dens', 'Sa_u', 'Sa_v', &
             'Sa_shum', 'Sa_shum_16O', 'Sa_shum_18O', 'Sa_shum_HDO'/)

    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          call addfld(fldListFr(compatm)%flds, trim(fldname))
          call addfld(fldListTo(compice)%flds, trim(fldname))
       else
          if ( fldchk(is_local%wrap%FBexp(compice)         , trim(fldname), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm ), trim(fldname), rc=rc)) then
             if (trim(fldname) == 'Sa_u' .or. trim(fldname) == 'Sa_v') then
                call addmap(fldListFr(compatm)%flds, trim(fldname), compice, mapconsf, 'one', atm2ice_vmap)
             else
                call addmap(fldListFr(compatm)%flds, trim(fldname), compice, mapconsf, 'one', atm2ice_smap)
             end if
             call addmrg(fldListTo(compice)%flds, trim(fldname), &
                  mrg_from1=compatm, mrg_fld1=trim(fldname), mrg_type1='copy')
          end if
       end if
    end do
    deallocate(flds)

    ! ---------------------------------------------------------------------
    ! to ice: sea surface temperature from ocn
    ! to ice: sea surface salinity from ocn
    ! to ice: zonal sea water velocity from ocn
    ! to ice: meridional sea water velocity from ocn
    ! to ice: zonal sea surface slope from ocean
    ! to ice: meridional sea surface slope from ocn
    ! ---------------------------------------------------------------------
    allocate(flds(6))
    flds = (/'So_t', 'So_s', 'So_u', 'So_v', 'So_dhdx', 'So_dhdy'/)

    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          call addfld(fldListFr(compocn)%flds, trim(fldname))
          call addfld(fldListTo(compice)%flds, trim(fldname))
       else
          if ( fldchk(is_local%wrap%FBexp(compice)        , trim(fldname), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compocn,compocn), trim(fldname), rc=rc)) then
             call addmap(fldListFr(compocn)%flds, trim(fldname), compice, mapfcopy , 'unset', 'unset')
             call addmrg(fldListTo(compice)%flds, trim(fldname), &
                  mrg_from1=compocn, mrg_fld1=trim(fldname), mrg_type1='copy')
          end if
       end if
    end do
    deallocate(flds)

    ! ---------------------------------------------------------------------
    ! to ice: ocean melt and freeze potential from ocn
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
    ! to ice: frozen runoff from rof and glc
    ! ---------------------------------------------------------------------
    do n = 1,size(iso)
       if (phase == 'advertise') then
          call addfld(fldListFr(comprof)%flds, 'Firr_rofi'//iso(n)) ! water flux into sea ice due to runoff (frozen)
          call addfld(fldListFr(compglc)%flds, 'Figg_rofi'//iso(n)) ! glc frozen runoff_iceberg flux to ice
          call addfld(fldListTo(compice)%flds, 'Fixx_rofi'//iso(n)) ! total frozen water flux into sea ice
       else
          if ( fldchk(is_local%wrap%FBExp(compice)         , 'Fixx_rofi'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(comprof, comprof), 'Forr_rofi'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compglc, compglc), 'Figg_rofi'//iso(n), rc=rc)) then

             call addmap(fldListFr(comprof)%flds, 'Forr_rofi'//iso(n), compice, mapconsf, 'none', rof2ocn_ice_rmap)
             call addmap(fldListFr(compglc)%flds, 'Figg_rofi'//iso(n), compice, mapconsf, 'one' , glc2ice_rmap)
             call addmrg(fldListTo(compice)%flds, 'Fixx_rofi'//iso(n), &
                  mrg_from1=comprof, mrg_fld1='Firr_rofi'//iso(n), mrg_type1='sum', &
                  mrg_from2=compglc, mrg_fld2='Figg_rofi'//iso(n), mrg_type2='sum')

          else if ( fldchk(is_local%wrap%FBExp(compice)    , 'Fixx_rofi'//iso(n), rc=rc) .and. &
                    fldchk(is_local%wrap%FBImp(comprof, comprof), 'Forr_rofi'//iso(n), rc=rc)) then

             call addmap(fldListFr(comprof)%flds, 'Forr_rofi'//iso(n), compice, mapconsf, 'none', rof2ocn_ice_rmap)
             call addmrg(fldListTo(compice)%flds, 'Fixx_rofi'//iso(n), &
                  mrg_from1=comprof, mrg_fld1='Firr_rofi'//iso(n), mrg_type1='sum')
          end if
       end if
    end do

    !=====================================================================
    ! FIELDS TO WAVE (compwav)
    !=====================================================================

    !----------------------------------------------------------
    ! to wav: fractional ice coverage wrt ocean from ice
    !----------------------------------------------------------
    if (phase == 'advertise') then
       ! the following is computed in med_phases_prep_wav
       call addfld(fldListFr(compice)%flds, 'Si_ifrac')
       call addfld(fldListTo(compwav)%flds, 'Si_ifrac')
    end if

    ! ---------------------------------------------------------------------
    ! to wav: ocean boundary layer depth from ocn
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compocn)%flds, 'So_bldepth')
       call addfld(fldListTo(compwav)%flds, 'So_bldepth')
    else
       if ( fldchk(is_local%wrap%FBImp(compocn, compocn), 'So_bldepth', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compwav)         , 'So_bldepth', rc=rc)) then
          call addmap(fldListFr(compocn)%flds, 'So_bldepth', compwav, mapbilnr, 'one', ocn2wav_smap)
          call addmrg(fldListTo(compwav)%flds, 'So_bldepth', mrg_from1=compocn, mrg_fld1='So_bldepth', mrg_type1='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to wav: zonal wind at the lowest model level from atm
    ! to wav: meridional wind at the lowest model level from atm
    ! ---------------------------------------------------------------------
    allocate(flds(2))
    flds = (/'Sa_u', 'Sa_v'/)

    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          call addfld(fldListFr(compatm)%flds, trim(fldname))
          call addfld(fldListTo(compwav)%flds, trim(fldname))
       else
          if ( fldchk(is_local%wrap%FBexp(compwav)         , trim(fldname), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm ), trim(fldname), rc=rc)) then
             call addmap(fldListFr(compatm)%flds, trim(fldname), compwav, mapbilnr, 'one', atm2wav_smap)
             call addmrg(fldListTo(compwav)%flds, trim(fldname), &
                  mrg_from1=compatm, mrg_fld1=trim(fldname), mrg_type1='copy')
          end if
       end if
    end do
    deallocate(flds)

    !=====================================================================
    ! FIELDS TO RIVER (comprof)
    !=====================================================================

    ! ---------------------------------------------------------------------
    ! to rof: water flux from land (liquid surface)
    ! to rof: water flux from land (liquid glacier, wetland, and lake)
    ! to rof: water flux from land (liquid subsurface)
    ! to rof: water flux from land direct to ocean
    ! to rof: irrigation flux from land (withdrawal from rivers)
    ! ---------------------------------------------------------------------
    ! TODO (mvertens, 2019-01-13): the following isotopes have not yet been defined in the NUOPC field dict
    ! allocate(flds(30))
    ! flds = (/'Flrl_rofsur', 'Flrl_rofsur_16O', 'Flrl_rofsur_18O', 'Flrl_rofsur_HDO', &
    !          'Flrl_rofgwl', 'Flrl_rofgwl_16O', 'Flrl_rofgwl_18O', 'Flrl_rofgwl_HDO', &
    !          'Flrl_rofsub', 'Flrl_rofsub_16O', 'Flrl_rofsub_18O', 'Flrl_rofsub_HDO', &
    !          'Flrl_rofdto', 'Flrl_rofdto_16O', 'Flrl_rofdto_18O', 'Flrl_rofdto_HDO', &
    !          'Flrl_rofi'  , 'Flrl_rofi_16O'  , 'Flrl_rofi_18O'  , 'Flrl_rofi_HDO'  , &
    !          'Flrl_irrig' , 'Flrl_irrig_16O' , 'Flrl_irrig_18O' , 'Flrl_irrig_HDO' /)

    allocate(flds(6))
    flds = (/'Flrl_rofsur', 'Flrl_rofgwl', 'Flrl_rofsub', 'Flrl_rofdto', 'Flrl_rofi', 'Flrl_irrig'/)

    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          call addfld(fldListFr(complnd)%flds, trim(fldname))
          call addfld(fldListTo(comprof)%flds, trim(fldname))
       else
          if ( fldchk(is_local%wrap%FBImp(complnd, complnd), trim(fldname), rc=rc) .and. &
               fldchk(is_local%wrap%FBExp(comprof)         , trim(fldname), rc=rc)) then
             call addmap(fldListFr(complnd)%flds, trim(fldname), comprof, mapconsd, 'lfrin', lnd2rof_fmap)
             call addmrg(fldListTo(comprof)%flds, trim(fldname), &
                  mrg_from1=complnd, mrg_fld1=trim(fldname), mrg_type1='copy_with_weights', mrg_fracname1='lfrac')
          end if
       end if
    end do
    deallocate(flds)

    !=====================================================================
    ! FIELDS TO LAND-ICE (compglc)
    !=====================================================================

    !-----------------------------
    ! to glc: from land
    !-----------------------------

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
          if ( fldchk(is_local%wrap%FBExp(complnd)         , 'Flgl_qice'//trim(cnum), rc=rc) .and. &
               fldchk(is_local%wrap%FBExp(complnd)         , 'Sl_tsrf'//trim(cnum)       , rc=rc) .and. &
               fldchk(is_local%wrap%FBExp(complnd)         , 'Sl_topo'//trim(cnum)    , rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compglc,compglc) , 'Sg_ice_covered'            , rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compglc,compglc) , 'Sg_topo'                   , rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compglc,compglc) , 'Flgg_hflx'                 , rc=rc)) then

             do num = 0, glc_nec
                cnum = glc_elevclass_as_string(num)
                call addmap(FldListFr(complnd)%flds, 'Flgl_qice'//trim(cnum), compglc, mapconsf, 'none', lnd2glc_fmap)
                call addmap(FldListFr(complnd)%flds, 'Sl_tsrf'//trim(cnum)  , compglc, mapbilnr, 'none', lnd2glc_smap)
                call addmap(FldListFr(complnd)%flds, 'Sl_topo'//trim(cnum)  , compglc, mapbilnr, 'none', lnd2glc_smap)
             end do
          end if
       end if
    end if

    !=====================================================================
    ! CO2 EXCHANGE
    !=====================================================================

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
