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
    use shr_sys_mod           , only : shr_sys_abort
    use esmflds

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
    character(len=4)    :: iso(1) ! TODO (mvertens, 2019-01-08): remove this
    character(len=4)    :: isos(3)
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
    logical             :: use_med_aoflux
    character(len=64), allocatable :: suffix(:)
    character(len=*), parameter    :: subname='(esmFlds_Init)'
    !--------------------------------------


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
    ! Determine if mediator computes aofluxes
    !---------------------------

    call NUOPC_CompAttributeGet(gcomp, name='compute_aofluxes_in_mediator', value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) use_med_aoflux
    call ESMF_LogWrite('compute_aofluxes_in_mediator = '// trim(cvalue), ESMF_LOGMSG_INFO, rc=dbrc)

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

    isos(1) = '_16O'
    isos(2) = '_18O'
    isos(3) = '_HDO'

    rc = ESMF_SUCCESS

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
    ! TO MEDIATOR component (for fractions and atm/ocn flux calculation)
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
       if ( ESMF_FieldBundleIsCreated(is_local%wrap%FBImp(compocn,compocn), rc=rc) .and. &
            ESMF_FieldBundleIsCreated(is_local%wrap%FBImp(compatm,compatm), rc=rc) .and. &
            use_med_aoflux) then 

          ! Following fields from the atm are mapped to the ocean grid in the mediator - BUT is not sent to the ocean
          if (.not. fldchk(is_local%wrap%FBImp(compatm,compatm), 'Sa_z', rc=rc)) then
             call shr_sys_abort(trim(subname)//'atm import Sa_z required for atm/ocn flux calculation') 
          end if
          if (.not. fldchk(is_local%wrap%FBImp(compatm,compatm), 'Sa_u', rc=rc)) then
             call shr_sys_abort(trim(subname)//'atm import Sa_u required for atm/ocn flux calculation') 
          end if
          if (.not. fldchk(is_local%wrap%FBImp(compatm,compatm), 'Sa_v', rc=rc)) then
             call shr_sys_abort(trim(subname)//'atm import Sa_v required for atm/ocn flux calculation') 
          end if
          if (.not. fldchk(is_local%wrap%FBImp(compatm,compatm), 'Sa_tbot', rc=rc)) then
             call shr_sys_abort(trim(subname)//'atm import Sa_tbot required for atm/ocn flux calculation') 
          end if
          if (.not. fldchk(is_local%wrap%FBImp(compatm,compatm), 'Sa_ptem', rc=rc)) then
             !call shr_sys_abort(trim(subname)//'atm import Sa_ptem required for atm/ocn flux calculation') 
          end if
          if (.not. fldchk(is_local%wrap%FBImp(compatm,compatm), 'Sa_dens', rc=rc)) then
             !call shr_sys_abort(trim(subname)//'atm import Sa_dens required for atm/ocn flux calculation') 
          end if
          do n = 1,size(iso)
             if (.not. fldchk(is_local%wrap%FBImp(compatm,compatm), trim('Sa_shum'//iso(n)), rc=rc)) then
                call shr_sys_abort(trim(subname)//trim('atm import Sa_shum'//iso(n))//' required for atm/ocn flux calculation') 
             end if
          end do
          call addmap(fldListFr(compatm)%flds, 'Sa_z'   , compocn, mapbilnr, 'one', atm2ocn_smap)
          call addmap(fldListFr(compatm)%flds, 'Sa_u'   , compocn, mappatch, 'one', atm2ocn_vmap)
          call addmap(fldListFr(compatm)%flds, 'Sa_v'   , compocn, mappatch, 'one', atm2ocn_vmap)
          call addmap(fldListFr(compatm)%flds, 'Sa_tbot', compocn, mapbilnr, 'one', atm2ocn_smap)
          call addmap(fldListFr(compatm)%flds, 'Sa_ptem', compocn, mapbilnr, 'one', atm2ocn_smap)
          call addmap(fldListFr(compatm)%flds, 'Sa_dens', compocn, mapbilnr, 'one', atm2ocn_smap)
          do n = 1,size(iso)
             call addmap(fldListFr(compatm)%flds, 'Sa_shum'//iso(n), compocn, mapbilnr, 'one', atm2ocn_smap)
          end do

          ! Following fields are not mapped since atm/ocn flux computation is assumed to be on ocean grid for now
          if (.not. fldchk(is_local%wrap%FBImp(compocn,compocn), 'So_t', rc=rc)) then
             call shr_sys_abort(trim(subname)//'ocn import So_t required for atm/ocn flux calculation') 
          end if
          if (.not. fldchk(is_local%wrap%FBImp(compocn,compocn), 'So_u', rc=rc)) then
             call shr_sys_abort(trim(subname)//'ocn import So_u required for atm/ocn flux calculation') 
          end if
          if (.not. fldchk(is_local%wrap%FBImp(compocn,compocn), 'So_v', rc=rc)) then
             call shr_sys_abort(trim(subname)//'ocn import So_v required for atm/ocn flux calculation') 
          end if
          if (.not. fldchk(is_local%wrap%FBImp(compocn,compocn), 'So_omask', rc=rc)) then
             call shr_sys_abort(trim(subname)//'ocn import So_omask required for atm/ocn flux calculation') 
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
    ! TO LAND
    !=====================================================================

    ! ---------------------------------------------------------------------
    ! to lnd: height at the lowest model level          from atm
    ! to lnd: Surface height                            from atm
    ! to lnd: zonal wind at the lowest model level      from atm
    ! to lnd: meridional wind at the lowest model level from atm
    ! ---------------------------------------------------------------------

    allocate(suffix(4))
    suffix = (/'z', 'topo', 'u', 'v'/) 

    do n = 1,size(suffix)
       fldname = 'Sa_'//trim(suffix(n))
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
    deallocate(suffix)

    ! ---------------------------------------------------------------------
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

    allocate(suffix(23))
    suffix = (/'lwdn', 'swndr', 'swvdr', 'swndf', 'swvdf', &
               'bcphidry', 'bcphodry', 'bcphiwet',         &
               'ocphidry', 'ocphodry', 'ocphiwet',         &
               'dstwet1', 'dstwet2', 'dstwet3', 'dstwet4', & 
               'dstdry1', 'dstdry2', 'dstdry3', 'dstdry4', &
               'noy', 'nhx'/) 

    do n = 1,size(suffix)
       fldname = 'Faxa_'//trim(suffix(n))
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
    deallocate(suffix)

    ! ---------------------------------------------------------------------
    ! to lnd: river channel total water volume from rof
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
    ! to lnd: river channel main channel water volume from rof
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
    ! to lnd: ice sheet grid coverage on global grid from glc
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compglc)%flds   , 'Sg_icemask')
       call addfld(fldListTo(complnd)%flds   , 'Sg_icemask')
    else
       if ( fldchk(is_local%wrap%FBExp(complnd)         , 'Sg_icemask', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compglc, compglc), 'Sg_icemask', rc=rc)) then
          call addmap(fldListFr(compglc)%flds, 'Sg_icemask', complnd,  mapconsf, 'one', glc2lnd_smap)
          call addmrg(fldListTo(complnd)%flds, 'Sg_icemask', &
               mrg_from1=compglc, mrg_fld1='Sg_icemask', mrg_type1='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to lnd: ice sheet mask where we are potentially sending non-zero fluxes from glc
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compglc)%flds, 'Sg_icemask_coupled_fluxes')
       call addfld(fldListTo(complnd)%flds, 'Sg_icemask_coupled_fluxes')
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
    ! FROM  ATMOSPHERE (this must be deleted)
    !=====================================================================

    ! ---------------------------------------------------------------------
    ! From atm: Temperature at the lowest model level
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Sa_tbot')
       call addfld(fldListTo(complnd)%flds, 'Sa_tbot')
       call addfld(fldListTo(compice)%flds, 'Sa_tbot')
    else
       if ( fldchk(is_local%wrap%FBexp(complnd)         , 'Sa_tbot', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm ), 'Sa_tbot', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Sa_tbot', complnd, mapbilnr, 'one', atm2lnd_smap)
          call addmrg(fldListTo(complnd)%flds, 'Sa_tbot', mrg_from1=compatm, mrg_fld1='Sa_tbot', mrg_type1='copy')
       end if
       if ( fldchk(is_local%wrap%FBexp(compice)         , 'Sa_tbot', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm ), 'Sa_tbot', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Sa_tbot', compice, mapbilnr, 'one', atm2ice_smap)
          call addmrg(fldListTo(compice)%flds, 'Sa_tbot', mrg_from1=compatm, mrg_fld1='Sa_tbot', mrg_type1='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to lnd and to ice: potential temperature at the lowest model level
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds , 'Sa_ptem')
       call addfld(fldListTo(complnd)%flds , 'Sa_ptem')
       call addfld(fldListTo(compice)%flds , 'Sa_ptem')
    else
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
    ! to lnd and ice: specific humidity at the lowest model level
    ! ---------------------------------------------------------------------
    do n = 1,size(iso)
       if (phase == 'advertise') then
          call addfld(fldListFr(compatm)%flds , 'Sa_shum'//iso(n))
          call addfld(fldListTo(complnd)%flds , 'Sa_shum'//iso(n))
          call addfld(fldListTo(compice)%flds , 'Sa_shum'//iso(n))
       else
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
          if (fldchk(is_local%wrap%FBImp(compatm,compatm), 'Sa_pbot', rc=rc) .and. &
              fldchk(is_local%wrap%FBExp(complnd)        , 'Sa_pbot', rc=rc)) then
             call addmap(fldListFr(compatm)%flds, 'Sa_pbot', complnd, mapbilnr, 'one', atm2lnd_smap)
             call addmrg(fldListTo(complnd)%flds, 'Sa_pbot', mrg_from1=compatm, mrg_fld1='Sa_pbot', mrg_type1='copy')
          end if
          if (fldchk(is_local%wrap%FBImp(compatm,compatm), 'Sa_pbot', rc=rc) .and. &
              fldchk(is_local%wrap%FBExp(compice)        , 'Sa_pbot', rc=rc)) then
             call addmap(fldListFr(compatm)%flds, 'Sa_pbot', compice, mapbilnr, 'one', atm2ice_smap)
             call addmrg(fldListTo(compice)%flds, 'Sa_pbot', mrg_from1=compatm, mrg_fld1='Sa_pbot', mrg_type1='copy')
          end if
    end if

    ! ---------------------------------------------------------------------
    ! From atm: Density at the lowest model level
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds , 'Sa_dens')
       call addfld(fldListTo(compice)%flds , 'Sa_dens')
    else
       if ( fldchk(is_local%wrap%FBexp(compice)         , 'Sa_dens', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm ), 'Sa_dens', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Sa_dens', compice, mapbilnr, 'one', atm2ice_smap)
          call addmrg(fldListTo(compice)%flds, 'Sa_dens', mrg_from1=compatm, mrg_fld1='Sa_dens', mrg_type1='copy')
       end if
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
          call addfld(fldListTo(compice)%flds, 'Faxa_rain' //iso(n))
       else
          if ( fldchk(is_local%wrap%FBexp(compice)        , 'Faxa_rain' //iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rainl'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rainc'//iso(n), rc=rc)) then
             call addmap(fldListFr(compatm)%flds, 'Faxa_rainc'//iso(n), compice, mapconsf, 'one', atm2ice_fmap)
             call addmap(fldListFr(compatm)%flds, 'Faxa_rainl'//iso(n), compice, mapconsf, 'one', atm2ice_fmap)
             call addmrg(fldListTo(compice)%flds, 'Faxa_rain' //iso(n) , &
                  mrg_from1=compatm, mrg_fld1=trim('Faxa_rainc'//iso(n))//':'//trim('Faxa_rainl'//iso(n)), &
                  mrg_type1='sum')
          end if
          if ( fldchk(is_local%wrap%FBexp(compice)        , 'Faxa_rain'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rain'//iso(n), rc=rc)) then
             call addmap(fldListFr(compatm)%flds, 'Faxa_rain'//iso(n), compice, mapconsf, 'one', atm2ice_fmap)
             call addmrg(fldListTo(compice)%flds, 'Faxa_rain'//iso(n), mrg_from1=compatm, mrg_fld1='Faxa_rain'//iso(n), &
                  mrg_type1='copy')
          end if
          if ( fldchk(is_local%wrap%FBexp(complnd)        , 'Faxa_rainl'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBexp(complnd)        , 'Faxa_rainl'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rainl'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rainc'//iso(n), rc=rc)) then
             call addmap(fldListFr(compatm)%flds, 'Faxa_rainc'//iso(n), complnd, mapconsf, 'one', atm2lnd_fmap)
             call addmap(fldListFr(compatm)%flds, 'Faxa_rainl'//iso(n), complnd, mapconsf, 'one', atm2lnd_fmap)
             call addmrg(fldListTo(complnd)%flds, 'Faxa_rainc'//iso(n), &
                  mrg_from1=compatm, mrg_fld1='Faxa_rainc'//iso(n), mrg_type1='copy')
             call addmrg(fldListTo(complnd)%flds, 'Faxa_rainl'//iso(n), &
                  mrg_from1=compatm, mrg_fld1='Faxa_rainl'//iso(n), mrg_type1='copy')
          end if
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
          call addfld(fldListTo(compice)%flds, 'Faxa_snow' //iso(n))
       else
          if ( fldchk(is_local%wrap%FBexp(compice)        , 'Faxa_snow' //iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_snowl'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_snowc'//iso(n), rc=rc)) then
             call addmap(fldListFr(compatm)%flds, 'Faxa_snowc'//iso(n), compice, mapconsf, 'one', atm2ice_fmap)
             call addmap(fldListFr(compatm)%flds, 'Faxa_snowl'//iso(n), compice, mapconsf, 'one', atm2ice_fmap)
             call addmrg(fldListTo(compice)%flds, 'Faxa_snow' //iso(n) , &
                  mrg_from1=compatm, mrg_fld1=trim('Faxa_snowc'//iso(n))//':'//trim('Faxa_snowl'//iso(n)), &
                  mrg_type1='sum')
          end if
          if ( fldchk(is_local%wrap%FBexp(compice)        , 'Faxa_snow' //iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_snow'//iso(n), rc=rc)) then
             call addmap(fldListFr(compatm)%flds, 'Faxa_snow'//iso(n), compice, mapconsf, 'one', atm2ice_fmap)
             call addmrg(fldListTo(compice)%flds, 'Faxa_snow'//iso(n), mrg_from1=compatm, mrg_fld1='Faxa_snow'//iso(n), &
                  mrg_type1='copy')
          end if
          if ( fldchk(is_local%wrap%FBexp(complnd)        , 'Faxa_snowl'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBexp(complnd)        , 'Faxa_snowl'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_snowl'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_snowc'//iso(n), rc=rc)) then
             call addmap(fldListFr(compatm)%flds, 'Faxa_snowc'//iso(n), complnd, mapconsf, 'one', atm2lnd_fmap)
             call addmap(fldListFr(compatm)%flds, 'Faxa_snowl'//iso(n), complnd, mapconsf, 'one', atm2lnd_fmap)
             call addmrg(fldListTo(complnd)%flds, 'Faxa_snowc'//iso(n), &
                  mrg_from1=compatm, mrg_fld1='Faxa_snowc'//iso(n), mrg_type1='copy')
             call addmrg(fldListTo(complnd)%flds, 'Faxa_snowl'//iso(n), &
                  mrg_from1=compatm, mrg_fld1='Faxa_snowl'//iso(n), mrg_type1='copy')
          end if
       end if
    end do
    ! ---------------------------------------------------------------------
    ! From atm: hydrophylic black carbon dry deposition flux
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_bcphidry')
       call addfld(fldListTo(compice)%flds, 'Faxa_bcphidry')
       call addfld(fldListTo(compocn)%flds, 'Faxa_bcphidry')
    else
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_bcphidry', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compice)        , 'Faxa_bcphidry', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_bcphidry', compice, mapconsf, 'one', atm2ice_fmap)
          call addmrg(fldListTo(compice)%flds, 'Faxa_bcphidry', mrg_from1=compatm, mrg_fld1='Faxa_bcphidry', &
               mrg_type1='copy')
       end if
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_bcphidry', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compice)        , 'Faxa_bcphidry', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_bcphidry', compocn, mapconsf, 'one', atm2ocn_fmap)
          call addmrg(fldListTo(compocn)%flds, 'Faxa_bcphidry', mrg_from1=compatm, mrg_fld1='Faxa_bcphidry', &
               mrg_type1='copy_with_weights', mrg_fracname1='ofrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! From atm: hydrophobic black carbon dry deposition flux
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_bcphodry')
       call addfld(fldListTo(compice)%flds, 'Faxa_bcphodry')
       call addfld(fldListTo(compocn)%flds, 'Faxa_bcphodry')
    else
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_bcphodry', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compice)        , 'Faxa_bcphodry', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_bcphodry', compice, mapconsf, 'one', atm2ice_fmap)
          call addmrg(fldListTo(compice)%flds, 'Faxa_bcphodry', mrg_from1=compatm, mrg_fld1='Faxa_bcphodry', &
               mrg_type1='copy')
       end if
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_bcphodry', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compice)        , 'Faxa_bcphodry', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_bcphodry', compocn, mapconsf, 'one', atm2ocn_fmap)
          call addmrg(fldListTo(compocn)%flds, 'Faxa_bcphodry', mrg_from1=compatm, mrg_fld1='Faxa_bcphodry', &
               mrg_type1='copy_with_weights', mrg_fracname1='ofrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! From atm: hydrophylic black carbon wet deposition flux
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_bcphiwet')
       call addfld(fldListTo(compice)%flds, 'Faxa_bcphiwet')
       call addfld(fldListTo(compocn)%flds, 'Faxa_bcphiwet')
    else
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_bcphiwet', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compice)        , 'Faxa_bcphiwet', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_bcphiwet', compice, mapconsf, 'one', atm2ice_fmap)
          call addmrg(fldListTo(compice)%flds, 'Faxa_bcphiwet', mrg_from1=compatm, mrg_fld1='Faxa_bcphiwet', &
               mrg_type1='copy')
       end if
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_bcphiwet', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compice)        , 'Faxa_bcphiwet', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_bcphiwet', compocn, mapconsf, 'one', atm2ocn_fmap)
          call addmrg(fldListTo(compocn)%flds, 'Faxa_bcphiwet', mrg_from1=compatm, mrg_fld1='Faxa_bcphiwet', &
               mrg_type1='copy_with_weights', mrg_fracname1='ofrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! From atm: hydrophylic organic carbon dry deposition flux
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_ocphidry')
       call addfld(fldListTo(compice)%flds, 'Faxa_ocphidry')
       call addfld(fldListTo(compocn)%flds, 'Faxa_ocphidry')
    else
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_ocphidry', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compice)        , 'Faxa_ocphidry', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_ocphidry', compice, mapconsf, 'one', atm2ice_fmap)
          call addmrg(fldListTo(compice)%flds, 'Faxa_ocphidry', mrg_from1=compatm, mrg_fld1='Faxa_ocphidry', &
               mrg_type1='copy')
       end if
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_ocphidry', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compice)        , 'Faxa_ocphidry', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_ocphidry', compocn, mapconsf, 'one', atm2ocn_fmap)
          call addmrg(fldListTo(compocn)%flds, 'Faxa_ocphidry', mrg_from1=compatm, mrg_fld1='Faxa_ocphidry', &
               mrg_type1='copy_with_weights', mrg_fracname1='ofrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! From atm: hydrophobic organic carbon dry deposition flux
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_ocphodry')
       call addfld(fldListTo(compice)%flds, 'Faxa_ocphodry')
       call addfld(fldListTo(compocn)%flds, 'Faxa_ocphodry')
    else
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_ocphodry', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compice)        , 'Faxa_ocphodry', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_ocphodry', compice, mapconsf, 'one', atm2ice_fmap)
          call addmrg(fldListTo(compice)%flds, 'Faxa_ocphodry', mrg_from1=compatm, mrg_fld1='Faxa_ocphodry', &
               mrg_type1='copy')
       end if
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_ocphodry', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compice)        , 'Faxa_ocphodry', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_ocphodry', compocn, mapconsf, 'one', atm2ocn_fmap)
          call addmrg(fldListTo(compocn)%flds, 'Faxa_ocphodry', mrg_from1=compatm, mrg_fld1='Faxa_ocphodry', &
               mrg_type1='copy_with_weights', mrg_fracname1='ofrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to lnd and ice and ocn: hydrophylic organic carbon wet deposition flux from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_ocphiwet')
       call addfld(fldListTo(compice)%flds, 'Faxa_ocphiwet')
       call addfld(fldListTo(compocn)%flds, 'Faxa_ocphiwet')
    else
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_ocphiwet', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compice)        , 'Faxa_ocphiwet', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_ocphiwet', compice, mapconsf, 'one', atm2ice_fmap)
          call addmrg(fldListTo(compice)%flds, 'Faxa_ocphiwet', mrg_from1=compatm, mrg_fld1='Faxa_ocphiwet', &
               mrg_type1='copy')
       end if
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_ocphiwet', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compice)        , 'Faxa_ocphiwet', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_ocphiwet', compocn, mapconsf, 'one', atm2ocn_fmap)
          call addmrg(fldListTo(compocn)%flds, 'Faxa_ocphiwet', mrg_from1=compatm, mrg_fld1='Faxa_ocphiwet', &
               mrg_type1='copy_with_weights', mrg_fracname1='ofrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to lnd and ice and ocn: dust wet deposition flux (size 1) from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_dstwet1')
       call addfld(fldListTo(compice)%flds, 'Faxa_dstwet1')
       call addfld(fldListTo(compocn)%flds, 'Faxa_dstwet1')
    else
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_dstwet1', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compice)        , 'Faxa_dstwet1', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_dstwet1', compice, mapconsf, 'one', atm2ice_fmap)
          call addmrg(fldListTo(compice)%flds, 'Faxa_dstwet1', mrg_from1=compatm, mrg_fld1='Faxa_dstwet1', &
               mrg_type1='copy')
       end if
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_dstwet1', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compice)        , 'Faxa_dstwet1', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_dstwet1', compocn, mapconsf, 'one', atm2ocn_fmap)
          call addmrg(fldListTo(compocn)%flds, 'Faxa_dstwet1', mrg_from1=compatm, mrg_fld1='Faxa_dstwet1', &
               mrg_type1='copy_with_weights', mrg_fracname1='ofrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to lnd and ice and ocn: dust wet deposition flux (size 2) from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_dstwet2')
       call addfld(fldListTo(compice)%flds, 'Faxa_dstwet2')
       call addfld(fldListTo(compocn)%flds, 'Faxa_dstwet2')
    else
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_dstwet2', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compice)        , 'Faxa_dstwet2', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_dstwet2', compice, mapconsf, 'one', atm2ice_fmap)
          call addmrg(fldListTo(compice)%flds, 'Faxa_dstwet2', mrg_from1=compatm, mrg_fld1='Faxa_dstwet2', &
               mrg_type1='copy')
       end if
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_dstwet2', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compice)        , 'Faxa_dstwet2', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_dstwet2', compocn, mapconsf, 'one', atm2ocn_fmap)
          call addmrg(fldListTo(compocn)%flds, 'Faxa_dstwet2', mrg_from1=compatm, mrg_fld1='Faxa_dstwet2', &
               mrg_type1='copy_with_weights', mrg_fracname1='ofrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to lnd and ice and ocn: dust wet deposition flux (size 3) from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_dstwet3')
       call addfld(fldListTo(compice)%flds, 'Faxa_dstwet3')
       call addfld(fldListTo(compocn)%flds, 'Faxa_dstwet3')
    else
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_dstwet3', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compice)        , 'Faxa_dstwet3', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_dstwet3', compice, mapconsf, 'one', atm2ice_fmap)
          call addmrg(fldListTo(compice)%flds, 'Faxa_dstwet3', mrg_from1=compatm, mrg_fld1='Faxa_dstwet3', &
               mrg_type1='copy')
       end if
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_dstwet3', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compice)        , 'Faxa_dstwet3', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_dstwet3', compocn, mapconsf, 'one', atm2ocn_fmap)
          call addmrg(fldListTo(compocn)%flds, 'Faxa_dstwet3', mrg_from1=compatm, mrg_fld1='Faxa_dstwet3', &
               mrg_type1='copy_with_weights', mrg_fracname1='ofrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to lnd and ice and ocn: dust wet deposition flux (size 4) from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_dstwet4')
       call addfld(fldListTo(compice)%flds, 'Faxa_dstwet4')
       call addfld(fldListTo(compocn)%flds, 'Faxa_dstwet4')
    else
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_dstwet4', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compice)        , 'Faxa_dstwet4', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_dstwet4', compice, mapconsf, 'one', atm2ice_fmap)
          call addmrg(fldListTo(compice)%flds, 'Faxa_dstwet4', mrg_from1=compatm, mrg_fld1='Faxa_dstwet4', &
               mrg_type1='copy')
       end if
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_dstwet4', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compice)        , 'Faxa_dstwet4', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_dstwet4', compocn, mapconsf, 'one', atm2ocn_fmap)
          call addmrg(fldListTo(compocn)%flds, 'Faxa_dstwet4', mrg_from1=compatm, mrg_fld1='Faxa_dstwet4', &
               mrg_type1='copy_with_weights', mrg_fracname1='ofrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to lnd and ice and ocn: dust dry deposition flux (size 1) from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_dstdry1')
       call addfld(fldListTo(compice)%flds, 'Faxa_dstdry1')
       call addfld(fldListTo(compocn)%flds, 'Faxa_dstdry1')
    else
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_dstdry1', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compice)        , 'Faxa_dstdry1', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_dstdry1', compice, mapconsf, 'one', atm2ice_fmap)
          call addmrg(fldListTo(compice)%flds, 'Faxa_dstdry1', mrg_from1=compatm, mrg_fld1='Faxa_dstdry1', &
               mrg_type1='copy')
       end if
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_dstdry1', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compice)        , 'Faxa_dstdry1', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_dstdry1', compocn, mapconsf, 'one', atm2ocn_fmap)
          call addmrg(fldListTo(compocn)%flds, 'Faxa_dstdry1', mrg_from1=compatm, mrg_fld1='Faxa_dstdry1', &
               mrg_type1='copy_with_weights', mrg_fracname1='ofrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to lnd and ice: dust dry deposition flux (size 2)
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_dstdry2')
       call addfld(fldListTo(compice)%flds, 'Faxa_dstdry2')
       call addfld(fldListTo(compocn)%flds, 'Faxa_dstdry2')
    else
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_dstdry2', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compice)        , 'Faxa_dstdry2', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_dstdry2', compice, mapconsf, 'one', atm2ice_fmap)
          call addmrg(fldListTo(compice)%flds, 'Faxa_dstdry2', mrg_from1=compatm, mrg_fld1='Faxa_dstdry2', &
               mrg_type1='copy')
       end if
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_dstdry2', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compice)        , 'Faxa_dstdry2', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_dstdry2', compocn, mapconsf, 'one', atm2ocn_fmap)
          call addmrg(fldListTo(compocn)%flds, 'Faxa_dstdry2', mrg_from1=compatm, mrg_fld1='Faxa_dstdry2', &
               mrg_type1='copy_with_weights', mrg_fracname1='ofrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to lnd and ice and ocn: dust dry deposition flux (size 3) from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_dstdry3')
       call addfld(fldListTo(compice)%flds, 'Faxa_dstdry3')
       call addfld(fldListTo(compocn)%flds, 'Faxa_dstdry3')
    else
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_dstdry3', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compice)        , 'Faxa_dstdry3', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_dstdry3', compice, mapconsf, 'one', atm2ice_fmap)
          call addmrg(fldListTo(compice)%flds, 'Faxa_dstdry3', mrg_from1=compatm, mrg_fld1='Faxa_dstdry3', &
               mrg_type1='copy')
       end if
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_dstdry3', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compice)        , 'Faxa_dstdry3', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_dstdry3', compocn, mapconsf, 'one', atm2ocn_fmap)
          call addmrg(fldListTo(compocn)%flds, 'Faxa_dstdry3', mrg_from1=compatm, mrg_fld1='Faxa_dstdry3', &
               mrg_type1='copy_with_weights', mrg_fracname1='ofrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to lnd and ice and ocn: dust dry deposition flux (size 4) from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_dstdry4')
       call addfld(fldListTo(compice)%flds, 'Faxa_dstdry4')
       call addfld(fldListTo(compocn)%flds, 'Faxa_dstdry4')
    else
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_dstdry4', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compice)        , 'Faxa_dstdry4', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_dstdry4', compice, mapconsf, 'one', atm2ice_fmap)
          call addmrg(fldListTo(compice)%flds, 'Faxa_dstdry4', mrg_from1=compatm, mrg_fld1='Faxa_dstdry4', &
               mrg_type1='copy')
       end if
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_dstdry4', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compice)        , 'Faxa_dstdry4', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_dstdry4', compocn, mapconsf, 'one', atm2ocn_fmap)
          call addmrg(fldListTo(compocn)%flds, 'Faxa_dstdry4', mrg_from1=compatm, mrg_fld1='Faxa_dstdry4', &
               mrg_type1='copy_with_weights', mrg_fracname1='ofrac')
       end if
    end if

    !=====================================================================
    ! TO ATMOSPHERE
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
    ! to atm: merged reference specific humidity at 2 meters
    ! to atm: merged 10m wind speed
    ! ---------------------------------------------------------------------
    allocate(suffix(3))
    suffix = (/'tref', 'qref', 'u10'/)

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
    ! to atm: merged reference specific humidity at 2 meters for isotopes from med aoflux
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       do n = 1,size(isos)
          call addfld(fldListFr(complnd)%flds , 'Sl_qref'//isos(n))
          call addfld(fldListFr(compice)%flds , 'Si_qref'//isos(n))
          call addfld(fldListMed_aoflux%flds  , 'So_qref'//isos(n))
          call addfld(fldListTo(compatm)%flds , 'Sx_qref'//isos(n))
       end do
    else
       do n = 1,size(isos)
          ! CESM (cam, non-aqua-planet)
          if ( fldchk(is_local%wrap%FBexp(compatm)         , 'Sx_qref'//isos(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(complnd,complnd ), 'Sl_qref'//isos(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compice,compice ), 'Si_qref'//isos(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBMed_aoflux_o         , 'So_qref'//isos(n), rc=rc)) then
             call addmap(fldListFr(complnd)%flds , 'Sl_qref'//isos(n), compatm, mapconsf, 'lfrin', lnd2atm_fmap)
             call addmap(fldListFr(compice)%flds , 'Si_qref'//isos(n), compatm, mapconsf, 'ifrac', ice2atm_fmap)
             call addmap(fldListMed_aoflux%flds  , 'So_qref'//isos(n), compocn, mapbilnr, 'one'  , atm2ocn_fmap) ! map atm->ocn
             call addmap(fldListMed_aoflux%flds  , 'So_qref'//isos(n), compatm, mapconsf, 'ofrac', ocn2atm_fmap) ! map ocn->atm
             call addmrg(fldListTo(compatm)%flds , 'Sx_qref'//isos(n), &
                  mrg_from1=complnd, mrg_fld1='Sl_qref'//isos(n), mrg_type1='merge', mrg_fracname1='lfrac', &
                  mrg_from2=compice, mrg_fld2='Si_qref'//isos(n), mrg_type2='merge', mrg_fracname2='ifrac', &
                  mrg_from3=compmed, mrg_fld3='So_qref'//isos(n), mrg_type3='merge', mrg_fracname3='ofrac')

          ! CESM (cam, aqua-planet)
          else if (fldchk(is_local%wrap%FBMed_aoflux_o, 'So_qref', rc=rc)) then
             call addmap(fldListMed_aoflux%flds , 'So_qref'//isos(n), compatm, mapconsf, 'ofrac', ocn2atm_fmap) ! map ocn->atm
             call addmrg(fldListTo(compatm)%flds, 'Sx_qref'//isos(n), &
                  mrg_from1=compmed, mrg_fld1='So_qref', mrg_type1='merge', mrg_fracname1='ofrac')
          end if
       end do
    end if

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
    ! to atm: surface snow depth             from ice
    ! to atm: mean ice volume per unit area  from ice
    ! to atm: mean snow volume per unit area from ice
    ! ---------------------------------------------------------------------
    allocate(suffix(3))
    suffix = (/'snowh', 'vice', 'vsno'/)

    do n = 1,size(suffix)
       fldname = 'Si_'//trim(suffix(n))
       if (phase == 'advertise') then
          call addfld(fldListFr(compice)%flds, trim(fldname))
          call addfld(fldListTo(compatm)%flds, trim(fldname))
       else
          if ( fldchk(is_local%wrap%FBexp(compatm)        , trim(fldname), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compice,compice), trim(fldname), rc=rc)) then
             call addmap(fldListFr(compice)%flds, trim(fldname), compatm, mapconsf, 'ifrac', ice2atm_fmap)
             call addmrg(fldListTo(compatm)%flds, trim(fldname), &
                  mrg_from1=compice, mrg_fld1=trim(fldname), mrg_type1='copy')
          end if
       end if
    end do
    deallocate(suffix)

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
          if (fldchk(is_local%wrap%FBexp(compatm) , trim(fldname), rc=rc) .and. &
               fldchk(is_local%wrap%FBMed_aoflux_o, trim(fldname), rc=rc)) then
             call addmap(fldListMed_aoflux%flds   , trim(fldname), compatm, mapconsf, 'ofrac', ocn2atm_fmap) ! map ocn->atm
             call addmrg(fldListTo(compatm)%flds  , trim(fldname), &
                  mrg_from1=compmed, mrg_fld1=trim(fldname), mrg_type1='copy')
          end if
       end if
    end do
    deallocate(suffix)

    ! ---------------------------------------------------------------------
    ! to atm: merged surface temperature and temperatures from ice and ocn
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(complnd)%flds, 'Sl_t')
       call addfld(fldListFr(compice)%flds, 'Si_t') 
       call addfld(fldListFr(compocn)%flds, 'So_t')

       call addfld(fldListTo(compatm)%flds, 'So_t') ! cesm, nems-frac
       call addfld(fldListTo(compatm)%flds, 'Si_t') ! nems-frac
       call addfld(fldListTo(compatm)%flds, 'Sx_t') ! cesm, nems-orig
    else
       ! CESM - merged ocn temp
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

       ! NEMS-orig - merged ocn temp
       else if (fldchk(is_local%wrap%FBexp(compatm)        , 'Sx_t', rc=rc) .and. &
                fldchk(is_local%wrap%FBImp(compocn,compocn), 'So_t', rc=rc) .and. &
                fldchk(is_local%wrap%FBImp(compice,compice), 'Si_t', rc=rc)) then
          call addmap(fldListFr(compice)%flds, 'Si_t', compatm, mapconsf , 'ifrac', ice2atm_fmap)
          call addmap(fldListFr(compocn)%flds, 'So_t', compatm, mapconsf , 'ofrac', ocn2atm_fmap)
          call addmrg(fldListTo(compatm)%flds, 'Sx_t', &
               mrg_from1=compice, mrg_fld1='Si_t', mrg_type1='merge', mrg_fracname1='ifrac', &
               mrg_from2=compocn, mrg_fld2='So_t', mrg_type2='merge', mrg_fracname2='ofrac')

       ! CESM aqua-planet - merged ocn temp
       else if ( fldchk(is_local%wrap%FBexp(compatm)        , 'Sx_t', rc=rc) .and. &
                 fldchk(is_local%wrap%FBImp(compocn,compocn), 'So_t', rc=rc)) then
          call addmap(fldListFr(compocn)%flds, 'So_t', compatm, mapconsf, 'ofrac', ocn2atm_fmap)
          call addmrg(fldListTo(compatm)%flds, 'Sx_t', &
               mrg_from1=compocn, mrg_fld1='So_t', mrg_type1='merge', mrg_fracname1='ofrac')
       end if

       ! NEMS-frac - unmerged ice temp
       if ( fldchk(is_local%wrap%FBexp(compatm)        , 'Si_t', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compice,compice), 'Si_t', rc=rc)) then
          call addmap(fldListFr(compice)%flds, 'Si_t', compatm, mapconsf , 'ifrac', ice2atm_fmap)
          call addmrg(fldListTo(compatm)%flds, 'Si_t', &
               mrg_from1=compice, mrg_fld1='Si_t', mrg_type1='copy')
       end if

       ! NEMS-frac and CESM - unmerged ocn tmep
       if ( fldchk(is_local%wrap%FBexp(compatm)        , 'So_t', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compocn,compocn), 'So_t', rc=rc)) then
          call addmap(fldListFr(compocn)%flds, 'So_t', compatm, mapconsf , 'ofrac', ocn2atm_fmap)
          call addmrg(fldListTo(compatm)%flds, 'So_t', mrg_from1=compocn, mrg_fld1='So_t', mrg_type1='copy')
       end if

    end if

    ! ---------------------------------------------------------------------
    ! to atm: merged zonal surface stress
    ! to atm: merged meridional surface stress
    ! to atm: merged surface latent heat flux
    ! to atm: merged surface sensible heat flux
    ! to atm: merged surface upward longwave heat flux
    ! ---------------------------------------------------------------------
    allocate(suffix(2))
    suffix = (/'taux', 'tauy', 'lat', 'sen', 'lwup', 'evap' /)

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
             call addmap(fldListMed_aoflux%flds  , 'Faox_'//trim(suffix(n)), compocn, mapconsf, 'one'  , atm2ocn_fmap) ! atm->ocn
             call addmap(fldListMed_aoflux%flds  , 'Faox_'//trim(suffix(n)), compatm, mapconsf, 'ofrac', ocn2atm_fmap) ! ocn->atm
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
             call addmap(fldListMed_aoflux%flds  , 'Faox_'//trim(suffix(n)), compocn, mapconsf, 'one'  , atm2ocn_fmap) ! atm->ocn
             call addmap(fldListMed_aoflux%flds  , 'Faox_'//trim(suffix(n)), compatm, mapconsf, 'ofrac', ocn2atm_fmap) ! ocn->atm
             call addmap(fldListFr(compice)%flds , 'Faii_'//trim(suffix(n)), compatm, mapconsf, 'ifrac', ice2atm_fmap)
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
             call addmap(fldListMed_aoflux%flds , 'Faox_'//trim(suffix(n)), compocn, mapconsf, 'one'  , atm2ocn_fmap) ! atm->ocn
             call addmap(fldListMed_aoflux%flds , 'Faox_'//trim(suffix(n)), compatm, mapconsf, 'ofrac', ocn2atm_fmap) ! ocn->atm
             call addmrg(fldListTo(compatm)%flds, 'Faxx_'//trim(suffix(n)), &
                  mrg_from1=compmed, mrg_fld1='Faox_'//trim(suffix(n)), mrg_type1='merge', mrg_fracname1='ofrac')
          end if
       end if
    end do
    deallocate(suffix)

    ! ---------------------------------------------------------------------
    ! to atm: evaporation water flux from water isotopes
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       do n = 1,size(isos)
          call addfld(fldListFr(complnd)%flds , 'Fall_evap'//isos(n))
          call addfld(fldListFr(compice)%flds , 'Faii_evap'//isos(n))
          call addfld(fldListMed_aoflux%flds  , 'Faox_evap'//isos(n))
          call addfld(fldListTo(compatm)%flds , 'Faii_evap'//isos(n)) ! nems-frac
          call addfld(fldListTo(compatm)%flds , 'Faxx_evap'//isos(n)) ! cesm, nems-orig
       end do
    else
       do n = 1,size(isos)
          ! CESM (cam, non-aqua-planet)
          if ( fldchk(is_local%wrap%FBImp(complnd,complnd), 'Fall_evap'//isos(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compice,compice), 'Faii_evap'//isos(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBMed_aoflux_o        , 'Faox_evap'//isos(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBexp(compatm)        , 'Faxx_evap'//isos(n), rc=rc)) then
             call addmap(fldListMed_aoflux%flds , 'Faox_evap'//isos(n), compatm, mapconsf, 'ofrac', ocn2atm_fmap) ! map ocn->atm
             call addmap(fldListMed_aoflux%flds , 'Faox_evap'//isos(n), compocn, mapconsf, 'one'  , atm2ocn_fmap) ! map atm->ocn
             call addmap(fldListFr(complnd)%flds, 'Fall_evap'//isos(n), compatm, mapconsf, 'lfrin', lnd2atm_fmap)
             call addmap(fldListFr(compice)%flds, 'Faii_evap'//isos(n), compatm, mapconsf, 'ifrac', ice2atm_fmap)
             call addmrg(fldListTo(compatm)%flds, 'Faxx_evap'//isos(n), &
                  mrg_from1=complnd, mrg_fld1='Fall_evap'//isos(n), mrg_type1='merge', mrg_fracname1='lfrac', &
                  mrg_from2=compice, mrg_fld2='Faii_evap'//isos(n), mrg_type2='merge', mrg_fracname2='ifrac', &
                  mrg_from3=compmed, mrg_fld3='Faox_evap'//isos(n), mrg_type3='merge', mrg_fracname3='ofrac')

          ! NEMS orig (here ofrac = 1.-ifrac)
          else if ( fldchk(is_local%wrap%FBImp(compice,compice), 'Faii_evap'//isos(n), rc=rc) .and. &
                    fldchk(is_local%wrap%FBMed_aoflux_o        , 'Faox_evap'//isos(n), rc=rc) .and. &
                    fldchk(is_local%wrap%FBexp(compatm)        , 'Faxx_evap'//isos(n), rc=rc)) then
             call addmap(fldListMed_aoflux%flds , 'Faox_evap'//isos(n), compatm, mapconsf, 'ofrac', ocn2atm_fmap) ! map ocn->atm
             call addmap(fldListMed_aoflux%flds , 'Faox_evap'//isos(n), compocn, mapconsf, 'one'  , atm2ocn_fmap) ! map atm->ocn
             call addmap(fldListFr(compice)%flds, 'Faii_evap'//isos(n), compatm, mapconsf, 'ifrac', ice2atm_fmap)
             call addmrg(fldListTo(compatm)%flds, 'Faxx_evap'//isos(n), &
                  mrg_from1=compice, mrg_fld1='Faii_evap'//isos(n), mrg_type1='merge', mrg_fracname1='ifrac', &
                  mrg_from2=compmed, mrg_fld2='Faox_evap'//isos(n), mrg_type2='merge', mrg_fracname2='ofrac')

          ! NEMS frac
          else if ( fldchk(is_local%wrap%FBImp(compice,compice), 'Faii_evap'//isos(n), rc=rc) .and. &
                    fldchk(is_local%wrap%FBexp(compatm)        , 'Faii_evap'//isos(n), rc=rc)) then
             call addmap(fldListFr(compice)%flds, 'Faii_evap'//isos(n), compatm, mapconsf, 'one', ice2atm_fmap)
             call addmrg(fldListTo(compatm)%flds, 'Faii_evap'//isos(n), &
                  mrg_from1=compice, mrg_fld1='Faii_evap'//isos(n), mrg_type1='merge', mrg_fracname1='ifrac')

          ! CESM (cam, aqua-planet)
          else if (fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_evap', rc=rc) .and. &
               fldchk(is_local%wrap%FBexp(compatm), 'Faxx_evap', rc=rc)) then
             call addmap(fldListMed_aoflux%flds , 'Faox_evap'//isos(n), compatm, mapconsf, 'ofrac', ocn2atm_fmap) ! map ocn->atm
             call addmap(fldListMed_aoflux%flds , 'Faox_evap'//isos(n), compocn, mapconsf, 'one'  , atm2ocn_fmap) ! map atm->ocn
             call addmrg(fldListTo(compatm)%flds, 'Faxx_evap', &
                  mrg_from1=compmed, mrg_fld1='Faox_evap', mrg_type1='merge', mrg_fracname1='ofrac')
          end if
       end do
    end if

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
    ! TO ocean component
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
    ! to ocn: downward longwave heat flux from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_lwdn')
       call addfld(fldListTo(compocn)%flds, 'Faxa_lwdn')
    else
       if ( fldchk(is_local%wrap%FBExp(compocn)        , 'Faxa_lwdn', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_lwdn', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_lwdn', compocn, mapconsf, 'one', atm2ocn_fmap)
          call addmrg(fldListTo(compocn)%flds, 'Faxa_lwdn', mrg_from1=compatm, mrg_fld1='Faxa_lwdn', &
               mrg_type1='copy_with_weights', mrg_fracname1='ofrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ocn: longwave net heat flux
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListMed_aoflux%flds  , 'Faox_lwup' )
       call addfld(fldListFr(compatm)%flds , 'Faxa_lwdn')
       call addfld(fldListFr(compatm)%flds , 'Foxx_lwnet')
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
                 fldchk(is_local%wrap%FBImp(compatm,compatm), 'Foxx_lwnet', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_lwdn' , compocn, mapconsf, 'one'  , atm2ocn_fmap)
          call addmap(fldListFr(compatm)%flds, 'Foxx_lwnet', compocn, mapconsf, 'one'  , atm2ocn_fmap)

      ! NEMS-frac (mom6) (send longwave net to ocean via auto merge)
       else if ( fldchk(is_local%wrap%FBExp(compocn)        , 'Foxx_lwnet', rc=rc) .and. &
                 fldchk(is_local%wrap%FBImp(compatm,compatm), 'Foxx_lwnet', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Foxx_lwnet', compocn, mapconsf, 'one'  , atm2ocn_fmap)
          call addmrg(fldListTo(compocn)%flds, 'Foxx_lwnet', &
               mrg_from1=compatm, mrg_fld1='Foxx_lwnet', mrg_type1='merge', mrg_fracname1='ofrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ocn: downward direct near-infrared incident solar radiation from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_swndr')
       call addfld(fldListTo(compocn)%flds, 'Faxa_swndr')
    else
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_swndr', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compocn)        , 'Faxa_swndr', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_swndr', compocn, mapconsf, 'one', atm2ocn_fmap)
          call addmrg(fldListTo(compocn)%flds, 'Faxa_swndr', mrg_from1=compatm, mrg_fld1='Faxa_swndr', &
               mrg_type1='copy_with_weights', mrg_fracname1='ofrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ocn: downward diffuse near-infrared incident solar radiation from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_swndf')
       call addfld(fldListTo(compocn)%flds, 'Faxa_swndf')
    else
       if ( fldchk(is_local%wrap%FBExp(compocn)        , 'Faxa_swndf', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_swndf', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_swndf', compocn, mapconsf, 'one', atm2ocn_fmap)
          call addmrg(fldListTo(compocn)%flds, 'Faxa_swndf', mrg_from1=compatm, mrg_fld1='Faxa_swndf', &
               mrg_type1='copy_with_weights', mrg_fracname1='ofrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ocn: downward Diffuse visible incident solar radiation from atm
    ! ---------------------------------------------------------------------

    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_swvdf')
       call addfld(fldListTo(compocn)%flds, 'Faxa_swvdf')
    else
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_swvdf', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compocn)        , 'Faxa_swvdf', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_swvdf', compocn, mapconsf, 'one', atm2ocn_fmap)
          call addmrg(fldListTo(compocn)%flds, 'Faxa_swvdf', mrg_from1=compatm, mrg_fld1='Faxa_swvdf', &
               mrg_type1='copy_with_weights', mrg_fracname1='ofrac')
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
    if (phase == 'advertise') then
       if (flds_i2o_per_cat) then
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
       end if
    else
       if (flds_i2o_per_cat) then
          call addmap(fldListFr(compice)%flds, 'Si_ifrac_n', compocn, mapfcopy, 'unset', 'unset')
          call addmap(fldListFr(compice)%flds, 'Fioi_swpen_ifrac_n', compocn, mapfcopy, 'unset', 'unset')
          ! TODO (mvertens, 2018-12-21): add mapping and merging
       end if
    end if

    ! ---------------------------------------------------------------------
    !  to ocn: precipitation rate water equivalent
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       do n = 1,size(iso)
          call addfld(fldListFr(compatm)%flds, 'Faxa_rainc'//iso(n))
          call addfld(fldListFr(compatm)%flds, 'Faxa_rainl'//iso(n))
          call addfld(fldListFr(compatm)%flds, 'Faxa_rain' //iso(n))
          call addfld(fldListTo(compocn)%flds, 'Faxa_rain' //iso(n))
       end do
    else
       do n = 1,size(iso)
          ! Note that the mediator atm/ocn flux calculation needs Faxa_rainc for the gustiness parameterization
          ! which by default is not actually used
          if ((fldchk(is_local%wrap%FBExp(compocn)        , 'Faxa_rain' //iso(n), rc=rc) .or. use_med_aoflux) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rainl'//iso(n), rc=rc) .and.  &
               fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rainc'//iso(n), rc=rc)) then

             call addmap(fldListFr(compatm)%flds, 'Faxa_rainl'//iso(n), compocn, mapconsf, 'one', atm2ocn_fmap)
             call addmap(fldListFr(compatm)%flds, 'Faxa_rainc'//iso(n), compocn, mapconsf, 'one', atm2ocn_fmap)
             call addmrg(fldListTo(compocn)%flds, 'Faxa_rain' //iso(n) , &
                  mrg_from1=compatm, mrg_fld1=trim('Faxa_rainc'//iso(n))//':'//trim('Faxa_rainl'//iso(n)), &
                  mrg_type1='sum_with_weights', mrg_fracname1='ofrac')
          end if
          if ( fldchk(is_local%wrap%FBExp(compocn)        , 'Faxa_rain'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rain'//iso(n), rc=rc)) then
             call addmap(fldListFr(compatm)%flds, 'Faxa_rain'//iso(n), compocn, mapconsf, 'one', atm2ocn_fmap)
             call addmrg(fldListTo(compocn)%flds, 'Faxa_rain'//iso(n), mrg_from1=compatm, mrg_fld1='Faxa_rain'//iso(n), &
                  mrg_type1='copy')
          end if
       end do
    end if

    ! ---------------------------------------------------------------------
    !  to ocn: snow rate water equivalent
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       do n = 1,size(iso)
          call addfld(fldListFr(compatm)%flds, 'Faxa_snowc'//iso(n))
          call addfld(fldListFr(compatm)%flds, 'Faxa_snowl'//iso(n))
          call addfld(fldListFr(compatm)%flds, 'Faxa_snow' //iso(n))
          call addfld(fldListTo(compocn)%flds, 'Faxa_snow' //iso(n))
       end do
    else
       do n = 1,size(iso)
          if ( fldchk(is_local%wrap%FBExp(compocn)        , 'Faxa_snow' //iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_snowl'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_snowc'//iso(n), rc=rc)) then
             call addmap(fldListFr(compatm)%flds, 'Faxa_snowl'//iso(n), compocn, mapconsf, 'one', atm2ocn_fmap)
             call addmap(fldListFr(compatm)%flds, 'Faxa_snowc'//iso(n), compocn, mapconsf, 'one', atm2ocn_fmap)
             call addmrg(fldListTo(compocn)%flds, 'Faxa_snow' //iso(n) , &
                  mrg_from1=compatm, mrg_fld1=trim('Faxa_snowc'//iso(n))//':'//trim('Faxa_snowl'//iso(n)), &
                  mrg_type1='sum_with_weights', mrg_fracname1='ofrac')
          end if
          if ( fldchk(is_local%wrap%FBExp(compocn)        , 'Faxa_snow' //iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_snow'//iso(n), rc=rc)) then
             call addmap(fldListFr(compatm)%flds, 'Faxa_snow'//iso(n), compocn, mapconsf, 'one', atm2ocn_fmap)
             call addmrg(fldListTo(compocn)%flds, 'Faxa_snow'//iso(n), mrg_from1=compatm, mrg_fld1='Faxa_snow'//iso(n), &
                  mrg_type1='copy')
          end if
       end do
    end if
    
    ! ---------------------------------------------------------------------
    ! to ocn: merge zonal surface stress from ice and (atm or med)
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListMed_aoflux%flds  , 'Faox_taux')
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
    ! to ocn: merge zonal meridional stress from ice and (atm or med)
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListMed_aoflux%flds  , 'Faox_tauy')
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
    ! to ocn: heat flux from melting ice form ice 
    ! ---------------------------------------------------------------------
    ! TODO (mvertens, 2019-01-07): is fioi_melth being handled correctly here? Is fd.yaml correctly aliasing Fioi_melth?
    do n = 1,size(iso)
       if (phase == 'advertise') then
          call addfld(fldListFr(compice)%flds , 'Fioi_melth'//iso(n))
          call addfld(fldListTo(compocn)%flds , 'Fioi_melth'//iso(n))
       else
          ! CESM
          if ( fldchk(is_local%wrap%FBexp(compocn)         , 'Fioi_melth'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compice, compice), 'Fioi_melth'//iso(n), rc=rc)) then
             call addmap(fldListFr(compice)%flds, 'Fioi_melth'//iso(n),    compocn,  mapfcopy, 'unset', 'unset')
             call addmrg(fldListTo(compocn)%flds, 'Fioi_melth'//iso(n), &
                  mrg_from1=compice, mrg_fld1='Fioi_melth'//iso(n), mrg_type1='copy_with_weights', mrg_fracname1='ifrac')
          end if
       end if
    end do

    ! ---------------------------------------------------------------------
    ! to ocn: sensible heat flux
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
    ! to ocn: surface latent heat flux from med
    ! ---------------------------------------------------------------------
    do n = 1,size(iso)
       if (phase == 'advertise') then
          call addfld(fldListMed_aoflux%flds , 'Faox_lat'//iso(n))
          call addfld(fldListTo(compocn)%flds, 'Foxx_lat'//iso(n))
       else
          ! CESM
          if ( fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_lat'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBexp(compocn), 'Foxx_lat'//iso(n), rc=rc)) then
             call addmrg(fldListTo(compocn)%flds, 'Foxx_lat'//iso(n), &
                  mrg_from1=compmed, mrg_fld1='Faox_lat'//iso(n), mrg_type1='merge', mrg_fracname1='ofrac')
          else
             ! Not used by NEMS orig and NEMS frac - latent derived from evap in med_phases_prep_ocn
          end if
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
    ! to ocn: evaporation water flux from atm and med
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds , 'Faxx_evap') ! nems-orig, nems-frac
       call addfld(fldListMed_aoflux%flds  , 'Faox_evap') ! cesm, nems-orig
       call addfld(fldListTo(compocn)%flds , 'Foxx_evap') ! cesm
    else
       if ( fldchk(is_local%wrap%FBExp(compocn), 'Foxx_evap', rc=rc) .and. &
            fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_evap', rc=rc)) then
          ! CESM
          call addmrg(fldListTo(compocn)%flds, 'Foxx_evap', &
               mrg_from1=compmed, mrg_fld1='Faox_evap', mrg_type1='merge', mrg_fracname1='ofrac')
       else
          ! NEMS-frac and NEMS-orig do not pass latent heat flux to mom6 - mom6 computes
          ! the latent heat flux from the imported evaporative flux
          ! TODO (mvertens, 2019-10-01): Can we unify this and have MOM6 use latent heat flux?
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
    allocate(suffix(16))
    suffix = (/'bcphidry', 'bcphodry', 'bcphiwet',         &
               'ocphidry', 'ocphodry', 'ocphiwet',         &
               'dstwet1', 'dstwet2', 'dstwet3', 'dstwet4', & 
               'dstdry1', 'dstdry2', 'dstdry3', 'dstdry4'/) 

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
    ! to ocn: water flux due to melting ice from ice
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
    ! to ocn: Hydrophobic black carbon deposition flux from ice
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

    !-----------------------------
    ! to ocn: liquid runoff from rof and glc components
    !-----------------------------
    if (phase == 'advertise') then
       do n = 1,size(iso)
          call addfld(fldListFr(compglc)%flds, 'Fogg_rofl'//iso(n))
          call addfld(fldListFr(comprof)%flds, 'Forr_rofl'//iso(n))
          call addfld(fldListTo(compocn)%flds, 'Foxx_rofl'//iso(n))
       end do
    else
       do n = 1,size(iso)
          if ( fldchk(is_local%wrap%FBExp(compocn)         , 'Foxx_rofl'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(comprof, comprof), 'Forr_rofl'//iso(n), rc=rc)) then
             ! water flux into sea water due to runoff (liquid)
             call addmap(fldListFr(comprof)%flds, 'Forr_rofl'//iso(n), compocn, mapfiler, 'none', rof2ocn_liq_rmap)
          end if
          if ( fldchk(is_local%wrap%FBExp(compocn)         , 'Fogg_rofl'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(comprof, comprof), 'Fogg_rofl'//iso(n), rc=rc)) then
             ! glc liquid runoff flux to ocean
             call addmap(fldListFr(compglc)%flds, 'Forr_rofl'//iso(n), compocn,  mapfiler, 'one', glc2ocn_liq_rmap)
          end if
          if ( fldchk(is_local%wrap%FBExp(compocn), 'Foxx_rofl'//iso(n), rc=rc)) then
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
       end do
    end if

    !-----------------------------
    ! to ocn: frozen runoff from rof and glc components
    !-----------------------------
    do n = 1,size(iso)
       if (phase == 'advertise') then
          call addfld(fldListFr(compglc)%flds, 'Fogg_rofi'//iso(n))
          call addfld(fldListFr(comprof)%flds, 'Forr_rofi'//iso(n))
          call addfld(fldListTo(compocn)%flds, 'Foxx_rofi'//iso(n))
       else
          if ( fldchk(is_local%wrap%FBExp(compocn)         , 'Foxx_rofi'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(comprof, comprof), 'Forr_rofi'//iso(n), rc=rc)) then
             ! water flux into sea water due to runoff (liquid)
             call addmap(fldListFr(comprof)%flds, 'Forr_rofi'//iso(n), compocn, mapfiler, 'none', rof2ocn_liq_rmap)
          end if
          if ( fldchk(is_local%wrap%FBExp(compocn)         , 'Fogg_rofi'//iso(n), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(comprof, comprof), 'Fogg_rofi'//iso(n), rc=rc)) then
             ! glc liquid runoff flux to ocean
             call addmap(fldListFr(compglc)%flds, 'Forr_rofi'//iso(n), compocn,  mapfiler, 'one', glc2ocn_liq_rmap)
          end if
          if ( fldchk(is_local%wrap%FBExp(compocn), 'Foxx_rofi'//iso(n), rc=rc)) then
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
    ! TO ice 
    !=====================================================================

    ! ---------------------------------------------------------------------
    ! to ice: downward longwave heat flux from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_lwdn')
       call addfld(fldListTo(compice)%flds, 'Faxa_lwdn')
    else
       if ( fldchk(is_local%wrap%FBExp(compice)        , 'Faxa_lwdn', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_lwdn', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_lwdn', compice, mapconsf, 'one', atm2ice_fmap)
          call addmrg(fldListTo(compice)%flds, 'Faxa_lwdn', mrg_from1=compatm, mrg_fld1='Faxa_lwdn', &
               mrg_type1='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ice: downward direct near-infrared incident solar radiation from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_swndr')
       call addfld(fldListTo(compice)%flds, 'Faxa_swndr')
    else
       if ( fldchk(is_local%wrap%FBExp(compice)        , 'Faxa_swndr', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_swndr', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_swndr', compice, mapconsf, 'one', atm2ice_fmap)
          call addmrg(fldListTo(compice)%flds, 'Faxa_swndr', mrg_from1=compatm, mrg_fld1='Faxa_swndr', &
               mrg_type1='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ice: downward direct visible incident solar radiation from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_swvdr')
       call addfld(fldListTo(compice)%flds, 'Faxa_swvdr')
    else
       if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_swvdr', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compice)        , 'Faxa_swvdr', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_swvdr', compice, mapconsf, 'one', atm2ocn_fmap)
          call addmrg(fldListTo(compice)%flds, 'Faxa_swvdr', mrg_from1=compatm, mrg_fld1='Faxa_swvdr', &
               mrg_type1='copy_with_weights', mrg_fracname1='ofrac')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ice: downward direct visible incident solar radiation from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_swvdr')
       call addfld(fldListTo(compice)%flds, 'Faxa_swvdr')
    else
       if ( fldchk(is_local%wrap%FBExp(compice)        , 'Faxa_swvdr', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_swvdr', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_swvdr', compice, mapconsf, 'one', atm2ice_fmap)
          call addmrg(fldListTo(compice)%flds, 'Faxa_swvdr', mrg_from1=compatm, mrg_fld1='Faxa_swvdr', &
               mrg_type1='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ice: downward diffuse near-infrared incident solar radiation from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_swndf')
       call addfld(fldListTo(compice)%flds, 'Faxa_swndf')
    else
       if ( fldchk(is_local%wrap%FBExp(compice)        , 'Faxa_swndf', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_swndf', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_swndf', compice, mapconsf, 'one', atm2ice_fmap)
          call addmrg(fldListTo(compice)%flds, 'Faxa_swndf', mrg_from1=compatm, mrg_fld1='Faxa_swndf', &
               mrg_type1='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ice: downward diffuse visible incident solar radiation from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Faxa_swvdf')
       call addfld(fldListTo(compice)%flds, 'Faxa_swvdf')
    else
       if ( fldchk(is_local%wrap%FBExp(compice)        , 'Faxa_swvdf', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_swvdf', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Faxa_swvdf', compice, mapconsf, 'one', atm2ice_fmap)
          call addmrg(fldListTo(compice)%flds, 'Faxa_swvdf', mrg_from1=compatm, mrg_fld1='Faxa_swvdf', &
               mrg_type1='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ice: height at the lowest model level from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Sa_z')
       call addfld(fldListTo(compice)%flds, 'Sa_z')
    else
       if ( fldchk(is_local%wrap%FBexp(compice)         , 'Sa_z', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm ), 'Sa_z', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Sa_z', compice, mapbilnr, 'one', atm2ice_smap)
          call addmrg(fldListTo(compice)%flds, 'Sa_z', mrg_from1=compatm, mrg_fld1='Sa_z', mrg_type1='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ice: zonal wind at the lowest model level from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Sa_u')
       call addfld(fldListTo(compice)%flds, 'Sa_u')
    else
       if ( fldchk(is_local%wrap%FBexp(compice)         , 'Sa_u', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm ), 'Sa_u', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Sa_u', compice, mappatch, 'one', atm2ice_vmap)
          call addmrg(fldListTo(compice)%flds, 'Sa_u', mrg_from1=compatm, mrg_fld1='Sa_u', mrg_type1='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ice: meridional wind at the lowest model level from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Sa_v')
       call addfld(fldListTo(compice)%flds, 'Sa_v')
    else
       if ( fldchk(is_local%wrap%FBexp(compice)         , 'Sa_v', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm ), 'Sa_v', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Sa_v', compice, mappatch, 'one', atm2ice_vmap)
          call addmrg(fldListTo(compice)%flds, 'Sa_v', mrg_from1=compatm, mrg_fld1='Sa_v', mrg_type1='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ice: sea surface temperature from ocn
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compocn)%flds, 'So_t')
       call addfld(fldListTo(compice)%flds, 'So_t')
    else
       if ( fldchk(is_local%wrap%FBexp(compice)        , 'So_t', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compocn,compocn), 'So_t', rc=rc)) then
          call addmap(fldListFr(compocn)%flds, 'So_t', compice, mapfcopy , 'unset', 'unset')
          call addmrg(fldListTo(compice)%flds, 'So_t', mrg_from1=compocn, mrg_fld1='So_t', mrg_type1='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ice: sea surface salinity from ocn
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compocn)%flds, 'So_s')
       call addfld(fldListTo(compice)%flds, 'So_s')
    else
       if ( fldchk(is_local%wrap%FBImp(compocn, compocn), 'So_s', rc=rc) .and. &
            fldchk(is_local%wrap%FBExp(compice)         , 'So_s', rc=rc)) then
          call addmap(fldListFr(compocn)%flds, 'So_s', compice,  mapfcopy, 'unset', 'unset')
          call addmrg(fldListTo(compice)%flds, 'So_s', mrg_from1=compocn, mrg_fld1='So_s', mrg_type1='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ice: zonal sea water velocity from ocn
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compocn)%flds, 'So_u')
       call addfld(fldListTo(compice)%flds, 'So_u')
    else
       if ( fldchk(is_local%wrap%FBExp(compice)         , 'So_u', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compocn, compocn), 'So_u', rc=rc)) then
          call addmap(fldListFr(compocn)%flds, 'So_u', compice,  mapfcopy, 'unset', 'unset')
          call addmrg(fldListTo(compice)%flds, 'So_u', mrg_from1=compocn, mrg_fld1='So_u', mrg_type1='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! to ice:  meridional sea water velocity from ocn
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compocn)%flds, 'So_v')
       call addfld(fldListTo(compice)%flds, 'So_v')
    else
       if ( fldchk(is_local%wrap%FBExp(compice)         , 'So_v', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compocn, compocn), 'So_v', rc=rc)) then
          call addmap(fldListFr(compocn)%flds, 'So_v', compice,  mapfcopy, 'unset', 'unset')
          call addmrg(fldListTo(compice)%flds, 'So_v', mrg_from1=compocn, mrg_fld1='So_v', mrg_type1='copy')
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
    ! to ice: meridional sea surface slope from ocn
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
          call addfld(fldListFr(comprof)%flds, 'Firr_rofi'//iso(n)) ! 'water flux into sea ice due to runoff (frozen)'
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

    !=====================================================================
    ! TO WAVE
    !=====================================================================

    !----------------------------------------------------------
    ! to wav: fractional ice coverage wrt ocean
    !----------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compice)%flds, 'Si_ifrac')
       call addfld(fldListTo(compwav)%flds, 'Si_ifrac')
    else
       if ( fldchk(is_local%wrap%FBexp(compwav)         , 'Si_ifrac', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compice,compice ), 'Si_ifrac', rc=rc)) then
          call addmrg(fldListTo(compwav)%flds, 'Si_ifrac', mrg_from1=compice, mrg_fld1='Si_ifrac', mrg_type1='copy')
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

    ! ---------------------------------------------------------------------
    ! To wav: zonal wind at the lowest model level from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Sa_u')
       call addfld(fldListTo(compwav)%flds, 'Sa_u')
    else
       if ( fldchk(is_local%wrap%FBexp(compwav)         , 'Sa_u', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm ), 'Sa_u', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Sa_u', compwav, mapbilnr, 'one', atm2wav_smap)
          call addmrg(fldListTo(compwav)%flds, 'Sa_u', mrg_from1=compatm, mrg_fld1='Sa_u', mrg_type1='copy')
       end if
    end if

    ! ---------------------------------------------------------------------
    ! To wav: meridional wind at the lowest model level from atm
    ! ---------------------------------------------------------------------
    if (phase == 'advertise') then
       call addfld(fldListFr(compatm)%flds, 'Sa_v')
       call addfld(fldListTo(compwav)%flds, 'Sa_v')
    else
       if ( fldchk(is_local%wrap%FBexp(compwav)         , 'Sa_v', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm,compatm ), 'Sa_v', rc=rc)) then
          call addmap(fldListFr(compatm)%flds, 'Sa_v', compwav, mapbilnr, 'one', atm2wav_smap)
          call addmrg(fldListTo(compwav)%flds, 'Sa_v', mrg_from1=compatm, mrg_fld1='Sa_v', mrg_type1='copy')
       end if
    end if

    !=====================================================================
    ! TO River
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

    !=====================================================================
    ! TO LAND ICE
    !=====================================================================

    !-----------------------------
    ! to glc: from lnd fields with multiple elevation classes
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
