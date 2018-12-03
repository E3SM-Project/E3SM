module esmFlds

  !---------------------------------------------------------------------
  ! This is a mediator specific routine that determines ALL possible
  ! fields exchanged between components and their associated routing,
  ! mapping and merging
  !---------------------------------------------------------------------

  use ESMF
  use NUOPC
  use med_constants_mod     , only : CX, CS, CL
  use shr_nuopc_fldList_mod , only : shr_nuopc_fldList_AddFld
  use shr_nuopc_fldList_mod , only : shr_nuopc_fldList_AddMap
  use shr_nuopc_fldList_mod , only : shr_nuopc_fldList_type
  use shr_nuopc_fldList_mod , only : mapbilnr, mapconsf, mapconsd, mappatch, mapfcopy, mapunset, mapfiler
  use shr_nuopc_scalars_mod , only : flds_scalar_name
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_chkerr
  use glc_elevclass_mod     , only : glc_elevclass_as_string

  implicit none
  public

  public :: esmFlds_Init

  !-----------------------------------------------
  ! Component and mapping array indices
  !-----------------------------------------------

  integer, parameter :: ncomps =8
  integer, parameter :: compmed=1
  integer, parameter :: compatm=2
  integer, parameter :: complnd=3
  integer, parameter :: compocn=4
  integer, parameter :: compice=5
  integer, parameter :: comprof=6
  integer, parameter :: compwav=7
  integer, parameter :: compglc=8
  character(len=*),parameter :: compname(ncomps) = (/'med','atm','lnd','ocn','ice','rof','wav','glc'/)

  type (shr_nuopc_fldList_type) :: fldListMed_aoflux_a
  type (shr_nuopc_fldList_type) :: fldListMed_aoflux_o
  type (shr_nuopc_fldList_type) :: fldListMed_ocnalb_o
  type (shr_nuopc_fldList_type) :: fldListMed_l2x_to_glc
  type (shr_nuopc_fldList_type) :: fldListMed_x2l_fr_glc
  type (shr_nuopc_fldList_type) :: fldListMed_g2x_to_lnd
  type (shr_nuopc_fldList_type) :: fldListTo(ncomps) ! advertise fields to components
  type (shr_nuopc_fldList_type) :: fldListFr(ncomps) ! advertise fields from components

  character(*), parameter :: u_FILE_u = &
       __FILE__

!================================================================================
contains
!================================================================================

  subroutine esmFlds_Init(gcomp, rc)

    ! input/output parameters:
    type(ESMF_GridComp)    :: gcomp
    integer, intent(inout) :: rc

    ! local variables:
    integer           :: ice_ncat                   ! number of sea ice thickness categories
    integer           :: glc_nec                    ! number of land-ice elevation classes
    integer           :: max_megan
    integer           :: max_ddep
    integer           :: max_fire
    logical           :: flds_i2o_per_cat
    integer           :: dbrc
    integer           :: num, i, n
    integer           :: n1, n2, n3, n4
    character(len=3)  :: cnum
    character(len=CL) :: cvalue
    character(len=CS) :: name, fldname
    character(len=CX) :: atm2ice_fmapname
    character(len=CX) :: atm2ice_smapname
    character(len=CX) :: atm2ice_vmapname
    character(len=CX) :: atm2lnd_fmapname
    character(len=CX) :: atm2lnd_smapname
    character(len=CX) :: atm2ocn_fmapname
    character(len=CX) :: atm2ocn_smapname
    character(len=CX) :: atm2ocn_vmapname
    character(len=CX) :: atm2wav_smapname
    character(len=CX) :: glc2lnd_fmapname
    character(len=CX) :: glc2lnd_smapname
    character(len=CX) :: glc2ice_rmapname
    character(len=CX) :: glc2ocn_liq_rmapname
    character(len=CX) :: glc2ocn_ice_rmapname
    character(len=CX) :: ice2atm_fmapname
    character(len=CX) :: ice2atm_smapname
    character(len=CX) :: ice2wav_smapname
    character(len=CX) :: lnd2atm_fmapname
    character(len=CX) :: lnd2atm_smapname
    character(len=CX) :: lnd2glc_fmapname
    character(len=CX) :: lnd2glc_smapname
    character(len=CX) :: lnd2rof_fmapname
    character(len=CX) :: ocn2atm_fmapname
    character(len=CX) :: ocn2atm_smapname
    character(len=CX) :: ocn2wav_smapname
    character(len=CX) :: rof2lnd_fmapname
    character(len=CX) :: rof2ocn_fmapname
    character(len=CX) :: rof2ocn_ice_rmapname
    character(len=CX) :: rof2ocn_liq_rmapname
    character(len=CX) :: wav2ocn_smapname
    logical           :: flds_co2a  ! use case
    logical           :: flds_co2b  ! use case
    logical           :: flds_co2c  ! use case
    character(len=*), parameter :: subname='(esmFlds_Init)'
    !--------------------------------------

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

    rc = ESMF_SUCCESS

    call NUOPC_CompAttributeGet(gcomp, name='flds_co2a', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) flds_co2a
    call ESMF_LogWrite('flds_co2a = '// trim(cvalue), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='flds_co2b', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) flds_co2b
    call ESMF_LogWrite('flds_co2b = '// trim(cvalue), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='flds_co2c', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) flds_co2c
    call ESMF_LogWrite('flds_co2c = '// trim(cvalue), ESMF_LOGMSG_INFO, rc=dbrc)

    !----------------------------------------------------------
    ! Initialize mapping file names
    !----------------------------------------------------------

    ! To atm

    call NUOPC_CompAttributeGet(gcomp, name='ice2atm_fmapname', value=ice2atm_fmapname, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('ice2atm_fmapname = '// trim(ice2atm_fmapname), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='ice2atm_smapname', value=ice2atm_smapname, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('ice2atm_smapname = '// trim(ice2atm_smapname), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='lnd2atm_fmapname', value=lnd2atm_fmapname, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('lnd2atm_fmapname = '// trim(lnd2atm_fmapname), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='ocn2atm_smapname', value=ocn2atm_smapname, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('ocn2atm_smapname = '// trim(ocn2atm_smapname), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='ocn2atm_fmapname', value=ocn2atm_fmapname, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('ocn2atm_fmapname = '// trim(ocn2atm_fmapname), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='lnd2atm_smapname', value=lnd2atm_smapname, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('lnd2atm_smapname = '// trim(lnd2atm_smapname), ESMF_LOGMSG_INFO, rc=dbrc)

    ! To lnd

    call NUOPC_CompAttributeGet(gcomp, name='atm2lnd_fmapname', value=atm2lnd_fmapname, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('atm2lnd_fmapname = '// trim(atm2lnd_fmapname), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='atm2lnd_smapname', value=atm2lnd_smapname, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('atm2lnd_smapname = '// trim(atm2lnd_smapname), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='rof2lnd_fmapname', value=rof2lnd_fmapname, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('rof2lnd_fmapname = '// trim(rof2lnd_fmapname), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='glc2lnd_fmapname', value=glc2lnd_fmapname, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('glc2lnd_smapname = '// trim(glc2lnd_fmapname), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='glc2lnd_smapname', value=glc2lnd_smapname, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('glc2lnd_smapname = '// trim(glc2lnd_smapname), ESMF_LOGMSG_INFO, rc=dbrc)

    ! To ice

    call NUOPC_CompAttributeGet(gcomp, name='atm2ice_fmapname', value=atm2ice_fmapname, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('atm2ice_fmapname = '// trim(atm2ice_fmapname), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='atm2ice_smapname', value=atm2ice_smapname, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('atm2ice_smapname = '// trim(atm2ice_smapname), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='atm2ice_vmapname', value=atm2ice_vmapname, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('atm2ice_vmapname = '// trim(atm2ice_vmapname), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='glc2ice_rmapname', value=glc2ice_rmapname, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('glc2ice_rmapname = '// trim(glc2ice_rmapname), ESMF_LOGMSG_INFO, rc=dbrc)

    ! To ocn

    call NUOPC_CompAttributeGet(gcomp, name='atm2ocn_fmapname', value=atm2ocn_fmapname, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('atm2ocn_fmapname = '// trim(atm2ocn_fmapname), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='atm2ocn_smapname', value=atm2ocn_smapname, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('atm2ocn_smapname = '// trim(atm2ocn_smapname), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='atm2ocn_vmapname', value=atm2ocn_vmapname, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('atm2ocn_vmapname = '// trim(atm2ocn_vmapname), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='glc2ocn_liq_rmapname', value=glc2ocn_liq_rmapname, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('glc2ocn_liq_rmapname = '// trim(glc2ocn_liq_rmapname), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='glc2ocn_ice_rmapname', value=glc2ocn_ice_rmapname, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('glc2ocn_ice_rmapname = '// trim(glc2ocn_ice_rmapname), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='wav2ocn_smapname', value=wav2ocn_smapname, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('wav2ocn_smapname = '// trim(wav2ocn_smapname), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='rof2ocn_fmapname', value=rof2ocn_fmapname, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('rof2ocn_fmapname = '// trim(rof2ocn_fmapname), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='rof2ocn_liq_rmapname', value=rof2ocn_liq_rmapname, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('rof2ocn_liq_rmapname = '// trim(rof2ocn_liq_rmapname), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='rof2ocn_ice_rmapname', value=rof2ocn_ice_rmapname, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('rof2ocn_ice_rmapname = '// trim(rof2ocn_ice_rmapname), ESMF_LOGMSG_INFO, rc=dbrc)

    ! To rof

    call NUOPC_CompAttributeGet(gcomp, name='lnd2rof_fmapname', value=lnd2rof_fmapname, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('lnd2rof_fmapname = '// trim(lnd2rof_fmapname), ESMF_LOGMSG_INFO, rc=dbrc)

    ! To glc

    call NUOPC_CompAttributeGet(gcomp, name='lnd2glc_fmapname', value=lnd2glc_fmapname, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('lnd2glc_fmapname = '// trim(lnd2glc_fmapname), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='lnd2glc_smapname', value=lnd2glc_smapname, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('lnd2glc_smapname = '// trim(lnd2glc_smapname), ESMF_LOGMSG_INFO, rc=dbrc)

    ! To wav

    call NUOPC_CompAttributeGet(gcomp, name='atm2wav_smapname', value=atm2wav_smapname, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('atm2wav_smapname = '// trim(atm2wav_smapname), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='ice2wav_smapname', value=ice2wav_smapname, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('ice2wav_smapname = '// trim(ice2wav_smapname), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='ocn2wav_smapname', value=ocn2wav_smapname, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('ocn2wav_smapname = '// trim(ocn2wav_smapname), ESMF_LOGMSG_INFO, rc=dbrc)

    !----------------------------------------------------------
    ! scalar information
    !----------------------------------------------------------

    do n = 1,ncomps
       call shr_nuopc_fldList_AddFld(fldListFr(n)%flds, trim(flds_scalar_name))
       call shr_nuopc_fldList_AddFld(fldListTo(n)%flds, trim(flds_scalar_name))
    end do

    !----------------------------------------------------------
    ! Masks from components
    !----------------------------------------------------------

    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds, 'Sl_lfrin')
    call shr_nuopc_fldList_AddFld(fldListFr(compocn)%flds, 'So_omask', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compocn)%flds(n1), compocn, compice,  mapfcopy, 'unset', 'unset')
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds, 'Si_imask')

    !----------------------------------------------------------
    ! Fractions sent to atm
    !----------------------------------------------------------

    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds, 'Sl_lfrac')
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds, 'Si_ifrac')
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds, 'So_ofrac')

    !----------------------------------------------------------
    ! Fractional ice coverage wrt ocean sent to ocn and wav
    !----------------------------------------------------------

    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds, 'Si_ifrac', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n1), compice, compocn,  mapfcopy, 'unset', 'unset')
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Si_ifrac', &
         merge_from1=compice, merge_field1='Si_ifrac', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compwav)%flds, 'Si_ifrac', &
         merge_from1=compice, merge_field1='Si_ifrac', merge_type1='copy')

    !----------------------------------------------------------
    ! Fields from atm
    !----------------------------------------------------------

    !  'Height at the lowest model level'
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Sa_z', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapbilnr, 'one', atm2lnd_smapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapbilnr, 'one', atm2ice_smapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapbilnr, 'one', atm2ocn_smapname)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Sa_z', &
         merge_from1=compatm, merge_field1='Sa_z', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Sa_z', &
         merge_from1=compatm, merge_field1='Sa_z', merge_type1='copy')

    ! 'Surface height'
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Sa_topo', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapbilnr, 'one', atm2lnd_smapname)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Sa_topo', &
         merge_from1=compatm, merge_field1='Sa_topo', merge_type1='copy')

    ! TODO: Sa_u and Sa_v are mapped to the ocean grid in the mediator - BUT are not sent to the ocean -
    ! They are only used in the atm/ocn flux calculation - so a special mapping will be done in the mediator
    ! for these fields

    !  'Zonal wind at the lowest model level'
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Sa_u', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapbilnr, 'one', atm2lnd_smapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mappatch, 'one', atm2ice_vmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mappatch, 'one', atm2ocn_vmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compwav, mapbilnr, 'one', atm2wav_smapname)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Sa_u', &
         merge_from1=compatm, merge_field1='Sa_u', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Sa_u', &
         merge_from1=compatm, merge_field1='Sa_u', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compwav)%flds, 'Sa_u', &
         merge_from1=compatm, merge_field1='Sa_u', merge_type1='copy')

    !  'Meridional wind at the lowest model level'
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Sa_v', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapbilnr, 'one', atm2lnd_smapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mappatch, 'one', atm2ice_vmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mappatch, 'one', atm2ocn_vmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compwav, mapbilnr, 'one', atm2wav_smapname)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Sa_v', &
         merge_from1=compatm, merge_field1='Sa_v', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Sa_v', &
         merge_from1=compatm, merge_field1='Sa_v', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compwav)%flds, 'Sa_v', &
         merge_from1=compatm, merge_field1='Sa_v', merge_type1='copy')

    !  'Temperature at the lowest model level'
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Sa_tbot', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapbilnr, 'one', atm2lnd_smapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapbilnr, 'one', atm2ice_smapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapbilnr, 'one', atm2ocn_smapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compwav, mapbilnr, 'one', atm2wav_smapname)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Sa_tbot', &
         merge_from1=compatm, merge_field1='Sa_tbot', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Sa_tbot', &
         merge_from1=compatm, merge_field1='Sa_tbot', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compwav)%flds, 'Sa_tbot', &
         merge_from1=compatm, merge_field1='Sa_tbot', merge_type1='copy')

    !  'Potential temperature at the lowest model level'
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Sa_ptem', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapbilnr, 'one', atm2lnd_smapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapbilnr, 'one', atm2ice_smapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapbilnr, 'one', atm2ocn_smapname)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Sa_ptem', &
         merge_from1=compatm, merge_field1='Sa_ptem', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Sa_ptem', &
         merge_from1=compatm, merge_field1='Sa_ptem', merge_type1='copy')

    !  'Specific humidity at the lowest model level'
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Sa_shum', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapbilnr, 'one', atm2lnd_smapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapbilnr, 'one', atm2ice_smapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapbilnr, 'one', atm2ocn_smapname)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Sa_shum', &
         merge_from1=compatm, merge_field1='Sa_shum', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Sa_shum', &
         merge_from1=compatm, merge_field1='Sa_shum', merge_type1='copy')

    !  'Pressure at the lowest model level'
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Sa_pbot', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapbilnr, 'one', atm2lnd_smapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapbilnr, 'one', atm2ice_smapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapbilnr, 'one', atm2ocn_smapname)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Sa_pbot', &
         merge_from1=compatm, merge_field1='Sa_pbot', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Sa_pbot', &
         merge_from1=compatm, merge_field1='Sa_pbot', merge_type1='copy')

    !  'Density at the lowest model level'
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Sa_dens', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapbilnr, 'one', atm2ice_smapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapbilnr, 'one', atm2ocn_smapname)
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Sa_dens', &
         merge_from1=compatm, merge_field1='Sa_dens', merge_type1='copy')

    !  'Convective and large scale precipitation rate'  water equivalent
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Faxa_rainc', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapconsf, 'one', atm2lnd_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapconsf, 'one', atm2ice_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapconsf, 'one', atm2ocn_fmapname)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Faxa_rainc', &
         merge_from1=compatm, merge_field1='Faxa_rainc', merge_type1='copy')

    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Faxa_rainl', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapconsf, 'one', atm2lnd_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapconsf, 'one', atm2ice_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapconsf, 'one', atm2ocn_fmapname)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Faxa_rainl', &
         merge_from1=compatm, merge_field1='Faxa_rainl', merge_type1='copy')

    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Faxa_rain', &
         merge_from1=compatm, merge_field1='Faxa_rainc:Faxa_rainl', merge_type1='accumulate', merge_fracname1='ofrac')
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Faxa_rain', &
         merge_from1=compatm, merge_field1='Faxa_rainc:Faxa_rainl', merge_type1='accumulate')

    !  'Convective and large-scale (stable) snow rate'
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Faxa_snowc', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapconsf, 'one', atm2lnd_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapconsf, 'one', atm2ice_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapconsf, 'one', atm2ocn_fmapname)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Faxa_snowc', &
         merge_from1=compatm, merge_field1='Faxa_snowc', merge_type1='copy')

    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Faxa_snowl', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapconsf, 'one', atm2lnd_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapconsf, 'one', atm2ice_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapconsf, 'one', atm2ocn_fmapname)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Faxa_snowl', &
         merge_from1=compatm, merge_field1='Faxa_snowl', merge_type1='copy')

    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Faxa_snow', &
         merge_from1=compatm, merge_field1='Faxa_snowc:Faxa_snowl', merge_type1='accumulate')
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Faxa_snow', &
         merge_from1=compatm, merge_field1='Faxa_snowc:Faxa_snowl', merge_type1='accumulate', merge_fracname1='ofrac')

    ! total precipitation to ocean (derived rain + snow, done AFTER mappings)
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Faxa_prec', &
         merge_from1=compatm, merge_field1='Faxa_rainc:Faxa_snowc:Faxa_rainl:Faxa_snowl', &
         merge_type1='accumulate', merge_fracname1='ofrac')

    !  'Downward longwave heat flux'
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Faxa_lwdn', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapconsf, 'one', atm2lnd_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapconsf, 'one', atm2ice_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapconsf, 'one', atm2ocn_fmapname)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Faxa_lwdn', &
         merge_from1=compatm, merge_field1='Faxa_lwdn', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Faxa_lwdn', &
         merge_from1=compatm, merge_field1='Faxa_lwdn', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Faxa_lwdn', &
         merge_from1=compatm, merge_field1='Faxa_lwdn', merge_type1='copy_with_weights', merge_fracname1='ofrac')

    !  'Direct near-infrared incident solar radiation'
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Faxa_swndr', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapconsf, 'one', atm2lnd_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapconsf, 'one', atm2ice_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapconsf, 'one', atm2ocn_fmapname)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Faxa_swndr', &
         merge_from1=compatm, merge_field1='Faxa_swndr', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Faxa_swndr', &
         merge_from1=compatm, merge_field1='Faxa_swndr', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Faxa_swndr', &
         merge_from1=compatm, merge_field1='Faxa_swndr', merge_type1='copy_with_weights', merge_fracname1='ofrac')

    !  'Direct visible incident solar radiation'
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Faxa_swvdr', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapconsf, 'one', atm2lnd_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapconsf, 'one', atm2ice_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapconsf, 'one', atm2ocn_fmapname)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Faxa_swvdr', &
         merge_from1=compatm, merge_field1='Faxa_swvdr', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Faxa_swvdr', &
         merge_from1=compatm, merge_field1='Faxa_swvdr', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Faxa_swvdr', &
         merge_from1=compatm, merge_field1='Faxa_swvdr', merge_type1='copy_with_weights', merge_fracname1='ofrac')

    !  'Diffuse near-infrared incident solar radiation'
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Faxa_swndf', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapconsf, 'one', atm2lnd_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapconsf, 'one', atm2ice_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapconsf, 'one', atm2ocn_fmapname)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Faxa_swndf', &
         merge_from1=compatm, merge_field1='Faxa_swndf', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Faxa_swndf', &
         merge_from1=compatm, merge_field1='Faxa_swndf', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Faxa_swndf', &
         merge_from1=compatm, merge_field1='Faxa_swndf', merge_type1='copy_with_weights', merge_fracname1='ofrac')

    ! 'Diffuse visible incident solar radiation'
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Faxa_swvdf', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapconsf, 'one', atm2lnd_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapconsf, 'one', atm2ice_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapconsf, 'one', atm2ocn_fmapname)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Faxa_swvdf', &
         merge_from1=compatm, merge_field1='Faxa_swvdf', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Faxa_swvdf', &
         merge_from1=compatm, merge_field1='Faxa_swvdf', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Faxa_swvdf', &
         merge_from1=compatm, merge_field1='Faxa_swvdf', merge_type1='copy_with_weights', merge_fracname1='ofrac')

    ! 'Net shortwave radiation'
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Faxa_swnet', fldindex=n1)  ! only diagnostic
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapconsf, 'one'  , atm2ice_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapconsf, 'one'  , atm2ocn_fmapname)
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds, 'Fall_swnet')
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds, 'Faii_swnet', fldindex=n1)  ! only diagnostic
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n1), compice, compocn, mapfcopy, 'unset', 'unset')
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Foxx_swnet') ! CUSTOM

    ! 'Net shortwave radiation penetrating into ice and ocean'
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds, 'Fioi_swpen', fldindex=n1) ! used for Foxx_swnet (see below)
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n1), compice, compocn, mapfcopy, 'unset', 'unset')

    ! 'Hydrophylic black carbon dry deposition flux'
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Faxa_bcphidry', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Faxa_bcphidry', &
         merge_from1=compatm, merge_field1='Faxa_bcphidry', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Faxa_bcphidry', &
         merge_from1=compatm, merge_field1='Faxa_bcphidry', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Faxa_bcphidry', &
         merge_from1=compatm, merge_field1='Faxa_bcphidry', merge_type1='copy_with_weights', merge_fracname1='ofrac')
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapconsf, 'one', atm2lnd_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapconsf, 'one', atm2ice_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapconsf, 'one', atm2ocn_fmapname)

    ! 'Hydrophobic black carbon dry deposition flux'
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Faxa_bcphodry', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Faxa_bcphodry', &
         merge_from1=compatm, merge_field1='Faxa_bcphodry', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Faxa_bcphodry', &
         merge_from1=compatm, merge_field1='Faxa_bcphodry', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Faxa_bcphodry', &
         merge_from1=compatm, merge_field1='Faxa_bcphodry', merge_type1='copy_with_weights', merge_fracname1='ofrac')
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapconsf, 'one', atm2lnd_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapconsf, 'one', atm2ice_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapconsf, 'one', atm2ocn_fmapname)

    ! 'Hydrophylic black carbon wet deposition flux'
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Faxa_bcphiwet', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Faxa_bcphiwet', &
         merge_from1=compatm, merge_field1='Faxa_bcphiwet', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Faxa_bcphiwet', &
         merge_from1=compatm, merge_field1='Faxa_bcphiwet', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Faxa_bcphiwet', &
         merge_from1=compatm, merge_field1='Faxa_bcphiwet', merge_type1='copy_with_weights', merge_fracname1='ofrac')
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapconsf, 'one', atm2lnd_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapconsf, 'one', atm2ice_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapconsf, 'one', atm2ocn_fmapname)

    ! 'Hydrophylic organic carbon dry deposition flux'
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Faxa_ocphidry', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Faxa_ocphidry', &
         merge_from1=compatm, merge_field1='Faxa_ocphidry', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Faxa_ocphidry', &
         merge_from1=compatm, merge_field1='Faxa_ocphidry', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Faxa_ocphidry', &
         merge_from1=compatm, merge_field1='Faxa_ocphidry', merge_type1='copy_with_weights', merge_fracname1='ofrac')
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapconsf, 'one', atm2lnd_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapconsf, 'one', atm2ice_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapconsf, 'one', atm2ocn_fmapname)

    ! 'Hydrophobic organic carbon dry deposition flux'
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Faxa_ocphodry', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Faxa_ocphodry', &
         merge_from1=compatm, merge_field1='Faxa_ocphodry', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Faxa_ocphodry', &
         merge_from1=compatm, merge_field1='Faxa_ocphodry', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Faxa_ocphodry', &
         merge_from1=compatm, merge_field1='Faxa_ocphodry', merge_type1='copy_with_weights', merge_fracname1='ofrac')
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapconsf, 'one', atm2lnd_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapconsf, 'one', atm2ice_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapconsf, 'one', atm2ocn_fmapname)

    ! 'Hydrophylic organic carbon wet deposition flux'
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Faxa_ocphiwet', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Faxa_ocphiwet', &
         merge_from1=compatm, merge_field1='Faxa_ocphiwet', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Faxa_ocphiwet', &
         merge_from1=compatm, merge_field1='Faxa_ocphiwet', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Faxa_ocphiwet', &
         merge_from1=compatm, merge_field1='Faxa_ocphiwet', merge_type1='copy_with_weights', merge_fracname1='ofrac')
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapconsf, 'one', atm2lnd_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapconsf, 'one', atm2ice_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapconsf, 'one', atm2ocn_fmapname)

    ! 'Dust wet deposition flux (size 1)'
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Faxa_dstwet1', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Faxa_dstwet1', &
         merge_from1=compatm, merge_field1='Faxa_dstwet1', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Faxa_dstwet1', &
         merge_from1=compatm, merge_field1='Faxa_dstwet1', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Faxa_dstwet1', &
         merge_from1=compatm, merge_field1='Faxa_dstwet1', merge_type1='copy_with_weights', merge_fracname1='ofrac')
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapconsf, 'one', atm2lnd_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapconsf, 'one', atm2ice_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapconsf, 'one', atm2ocn_fmapname)

    ! 'Dust wet deposition flux (size 2)'
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Faxa_dstwet2', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Faxa_dstwet2', &
         merge_from1=compatm, merge_field1='Faxa_dstwet2', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Faxa_dstwet2', &
         merge_from1=compatm, merge_field1='Faxa_dstwet2', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Faxa_dstwet2', &
         merge_from1=compatm, merge_field1='Faxa_dstwet2', merge_type1='copy_with_weights', merge_fracname1='ofrac')
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapconsf, 'one', atm2lnd_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapconsf, 'one', atm2ice_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapconsf, 'one', atm2ocn_fmapname)

    ! 'Dust wet deposition flux (size 3)'
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Faxa_dstwet3', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Faxa_dstwet3', &
         merge_from1=compatm, merge_field1='Faxa_dstwet3', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Faxa_dstwet3', &
         merge_from1=compatm, merge_field1='Faxa_dstwet3', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Faxa_dstwet3', &
         merge_from1=compatm, merge_field1='Faxa_dstwet3', merge_type1='copy_with_weights', merge_fracname1='ofrac')
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapconsf, 'one', atm2lnd_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapconsf, 'one', atm2ice_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapconsf, 'one', atm2ocn_fmapname)

    ! 'Dust wet deposition flux (size 4)'
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Faxa_dstwet4', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Faxa_dstwet4', &
         merge_from1=compatm, merge_field1='Faxa_dstwet4', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Faxa_dstwet4', &
         merge_from1=compatm, merge_field1='Faxa_dstwet4', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Faxa_dstwet4', &
         merge_from1=compatm, merge_field1='Faxa_dstwet4', merge_type1='copy_with_weights', merge_fracname1='ofrac')
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapconsf, 'one', atm2lnd_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapconsf, 'one', atm2ice_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapconsf, 'one', atm2ocn_fmapname)

    ! 'Dust dry deposition flux (size 1)'
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Faxa_dstdry1', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Faxa_dstdry1', &
         merge_from1=compatm, merge_field1='Faxa_dstdry1', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Faxa_dstdry1', &
         merge_from1=compatm, merge_field1='Faxa_dstdry1', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Faxa_dstdry1', &
         merge_from1=compatm, merge_field1='Faxa_dstdry1', merge_type1='copy_with_weights', merge_fracname1='ofrac')
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapconsf, 'one', atm2lnd_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapconsf, 'one', atm2ice_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapconsf, 'one', atm2ocn_fmapname)

    ! 'Dust dry deposition flux (size 2)'
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Faxa_dstdry2', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Faxa_dstdry2', &
         merge_from1=compatm, merge_field1='Faxa_dstdry2', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Faxa_dstdry2', &
         merge_from1=compatm, merge_field1='Faxa_dstdry2', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Faxa_dstdry2', &
         merge_from1=compatm, merge_field1='Faxa_dstdry2', merge_type1='copy_with_weights', merge_fracname1='ofrac')
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapconsf, 'one', atm2lnd_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapconsf, 'one', atm2ice_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapconsf, 'one', atm2ocn_fmapname)

    ! 'Dust dry deposition flux (size 3)'
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Faxa_dstdry3', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Faxa_dstdry3', &
         merge_from1=compatm, merge_field1='Faxa_dstdry3', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Faxa_dstdry3', &
         merge_from1=compatm, merge_field1='Faxa_dstdry3', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Faxa_dstdry3', &
         merge_from1=compatm, merge_field1='Faxa_dstdry3', merge_type1='copy_with_weights', merge_fracname1='ofrac')
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapconsf, 'one', atm2lnd_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapconsf, 'one', atm2ice_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapconsf, 'one', atm2ocn_fmapname)

    ! 'Dust dry deposition flux (size 4)'
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Faxa_dstdry4', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Faxa_dstdry4', &
         merge_from1=compatm, merge_field1='Faxa_dstdry4', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Faxa_dstdry4', &
         merge_from1=compatm, merge_field1='Faxa_dstdry4', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Faxa_dstdry4', &
         merge_from1=compatm, merge_field1='Faxa_dstdry4', merge_type1='copy_with_weights', merge_fracname1='ofrac')
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapconsf, 'one', atm2lnd_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapconsf, 'one', atm2ice_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapconsf, 'one', atm2ocn_fmapname)

    !----------------------------------------------------------
    ! states/fluxes to atm (and ocean)
    !----------------------------------------------------------

    ! 'Direct albedo (visible radiation)'
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds , 'Sl_avsdr', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1) , complnd, compatm, mapconsf, 'lfrin', lnd2atm_smapname)
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds , 'Si_avsdr', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n1) , compice, compatm, mapconsf, 'ifrac', ice2atm_smapname)
    call shr_nuopc_fldList_AddFld(fldListMed_ocnalb_o%flds, 'So_avsdr', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListMed_ocnalb_o%flds(n1), compocn, compatm, mapconsf, 'ofrac', ocn2atm_smapname)
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds , 'Sx_avsdr', &
         merge_from1=complnd, merge_field1='Sl_avsdr', merge_type1='merge', merge_fracname1='lfrac', &
         merge_from2=compice, merge_field2='Si_avsdr', merge_type2='merge', merge_fracname2='ifrac', &
         merge_from3=compmed, merge_field3='So_avsdr', merge_type3='merge', merge_fracname3='ofrac')

    ! 'Direct albedo (near-infrared radiation)'
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds , 'Sl_anidr', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1) , complnd, compatm, mapconsf, 'lfrin', lnd2atm_fmapname)
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds , 'Si_anidr', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n1) , compice, compatm, mapconsf, 'ifrac', ice2atm_fmapname)
    call shr_nuopc_fldList_AddFld(fldListMed_ocnalb_o%flds, 'So_anidr', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListMed_ocnalb_o%flds(n1), compocn, compatm, mapconsf, 'ofrac', ocn2atm_smapname)
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds , 'Sx_anidr', &
         merge_from1=complnd, merge_field1='Sl_anidr', merge_type1='merge', merge_fracname1='lfrac', &
         merge_from2=compice, merge_field2='Si_anidr', merge_type2='merge', merge_fracname2='ifrac', &
         merge_from3=compmed, merge_field3='So_anidr', merge_type3='merge', merge_fracname3='ofrac')

    ! 'Diffuse albedo (visible radiation)'
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds , 'Sl_avsdf', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1) , complnd, compatm, mapconsf, 'lfrin', lnd2atm_fmapname)
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds , 'Si_avsdf', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n1) , compice, compatm, mapconsf, 'ifrac', ice2atm_fmapname)
    call shr_nuopc_fldList_AddFld(fldListMed_ocnalb_o%flds, 'So_avsdf', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListMed_ocnalb_o%flds(n1), compocn, compatm, mapconsf, 'ofrac', ocn2atm_smapname)
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds , 'Sx_avsdf', &
         merge_from1=complnd, merge_field1='Sl_avsdf', merge_type1='merge', merge_fracname1='lfrac', &
         merge_from2=compice, merge_field2='Si_avsdf', merge_type2='merge', merge_fracname2='ifrac', &
         merge_from3=compmed, merge_field3='So_avsdf', merge_type3='merge', merge_fracname3='ofrac')

    ! 'Diffuse albedo (near-infrared radiation)'
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds , 'Sl_anidf', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1) , complnd, compatm, mapconsf, 'lfrin', lnd2atm_fmapname)
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds , 'Si_anidf', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n1) , compice, compatm, mapconsf, 'ifrac', ice2atm_fmapname)
    call shr_nuopc_fldList_AddFld(fldListMed_ocnalb_o%flds, 'So_anidf', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListMed_ocnalb_o%flds(n1), compocn, compatm, mapconsf, 'ofrac', ocn2atm_smapname)
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds , 'Sx_anidf', &
         merge_from1=complnd, merge_field1='Sl_anidf', merge_type1='merge', merge_fracname1='lfrac', &
         merge_from2=compice, merge_field2='Si_anidf', merge_type2='merge', merge_fracname2='ifrac', &
         merge_from3=compmed, merge_field3='So_anidf', merge_type3='merge', merge_fracname3='ofrac')

    ! 'Reference temperature at 2 meters'
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds , 'Sl_tref', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1) , complnd, compatm, mapconsf, 'lfrin', lnd2atm_fmapname)
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds , 'Si_tref', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n1) , compice, compatm, mapconsf, 'ifrac', ice2atm_fmapname)
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_a%flds, 'So_tref', fldindex=n1) ! Needed only for merging
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_a%flds(n1), compatm, compocn, mapbilnr, 'one'  , atm2ocn_fmapname) ! map atm->ocn
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_o%flds, 'So_tref', fldindex=n1) ! Needed only for merging
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_o%flds(n1), compocn, compatm, mapconsf, 'ofrac', ocn2atm_fmapname) ! map ocn->atm
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds , 'Sx_tref', &
         merge_from1=complnd, merge_field1='Sl_tref', merge_type1='merge', merge_fracname1='lfrac', &
         merge_from2=compice, merge_field2='Si_tref', merge_type2='merge', merge_fracname2='ifrac', &
         merge_from3=compmed, merge_field3='So_tref', merge_type3='merge', merge_fracname3='ofrac')

    ! 'Reference specific humidity at 2 meters'
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds , 'Sl_qref', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1) , complnd, compatm, mapconsf, 'lfrin', lnd2atm_fmapname)
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds , 'Si_qref', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n1) , compice, compatm, mapconsf, 'ifrac', ice2atm_fmapname)
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_a%flds, 'So_qref', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_a%flds(n1), compatm, compocn, mapbilnr, 'one'  , atm2ocn_fmapname) ! map atm->ocn
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_o%flds, 'So_qref', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_o%flds(n1), compocn, compatm, mapconsf, 'ofrac', ocn2atm_fmapname) ! map ocn->atm
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds , 'Sx_qref', &
         merge_from1=complnd, merge_field1='Sl_qref', merge_type1='merge', merge_fracname1='lfrac', &
         merge_from2=compice, merge_field2='Si_qref', merge_type2='merge', merge_fracname2='ifrac', &
         merge_from3=compmed, merge_field3='So_qref', merge_type3='merge', merge_fracname3='ofrac')

    ! 'Surface temperature'
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds, 'Sl_t', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1), complnd, compatm, mapconsf , 'lfrin', lnd2atm_fmapname)
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds, 'Si_t', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n1), compice, compatm, mapconsf , 'ifrac', ice2atm_fmapname)
    call shr_nuopc_fldList_AddFld(fldListFr(compocn)%flds, 'So_t', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compocn)%flds(n1), compocn, compatm, mapconsf , 'ofrac', ocn2atm_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compocn)%flds(n1), compocn, compwav, mapbilnr , 'one'  , ocn2wav_smapname) ! This will be a custom map - need to name it however
    call shr_nuopc_fldList_AddMap(fldListFr(compocn)%flds(n1), compocn, compice, mapfcopy , 'unset', 'unset')
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds, 'Sx_t', &
         merge_from1=complnd, merge_field1='Sl_t', merge_type1='merge', merge_fracname1='lfrac', &
         merge_from2=compice, merge_field2='Si_t', merge_type2='merge', merge_fracname2='ifrac', &
         merge_from3=compocn, merge_field3='So_t', merge_type3='merge', merge_fracname3='ofrac')
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds, 'So_t', &
         merge_from1=compocn, merge_field1='So_t', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'So_t', &
         merge_from1=compocn, merge_field1='So_t', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compwav)%flds, 'So_t', &
         merge_from1=compocn, merge_field1='So_t', merge_type1='copy')


    ! 'Surface fraction velocity in land'
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds, 'Sl_fv', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1), complnd, compatm, mapconsf, 'lfrin', lnd2atm_fmapname)
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds, 'Sl_fv', &
         merge_from1=complnd, merge_field1='Sl_fv', merge_type1='copy')

    ! 'Aerodynamic resistance'
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds, 'Sl_ram1', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1), complnd, compatm, mapconsf, 'lfrin', lnd2atm_fmapname)
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds, 'Sl_ram1',&
         merge_from1=complnd, merge_field1='Sl_ram1', merge_type1='copy')

    ! 'Surface snow water equivalent'
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds, 'Sl_snowh', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1), complnd, compatm, mapconsf, 'lfrin', lnd2atm_fmapname)
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds, 'Sl_snowh', &
         merge_from1=complnd, merge_field1='Sl_snowh', merge_type1='copy')

    ! 'Surface snow depth'
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds, 'Si_snowh', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n1), compice, compatm, mapconsf, 'ifrac', ice2atm_fmapname)
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds, 'Si_snowh', &
         merge_from1=compice, merge_field1='Si_snowh', merge_type1='copy')

    ! 'Surface saturation specific humidity in ocean'
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_a%flds, 'So_ssq', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_a%flds(n1), compatm, compocn, mapbilnr, 'one'  , atm2ocn_fmapname) ! map atm->ocn
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_o%flds, 'So_ssq', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_o%flds(n1), compocn, compatm, mapconsf, 'ofrac', ocn2atm_fmapname) ! map ocn->atm
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds , 'So_ssq', &
         merge_from1=compmed, merge_field1='So_ssq', merge_type1='copy')

    ! 'Square of exch. coeff (tracers)'
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_a%flds, 'So_re', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_a%flds(n1), compatm, compocn, mapbilnr, 'one'  , atm2ocn_fmapname) ! map atm->ocn
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_o%flds, 'So_re', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_o%flds(n1), compocn, compatm, mapconsf, 'ofrac', ocn2atm_fmapname) ! map ocn->atm
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds , 'So_re', &
         merge_from1=compmed, merge_field1='So_re', merge_type1='copy')

    ! '10m wind'
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds , 'Sl_u10', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1) , complnd, compatm, mapconsf, 'lfrin', lnd2atm_fmapname)
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds , 'Si_u10', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n1) , compice, compatm, mapconsf, 'ifrac', ice2atm_fmapname)
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_a%flds, 'So_u10', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_a%flds(n1), compatm, compocn, mapbilnr, 'one'  , atm2ocn_fmapname) ! map atm->ocn
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_o%flds, 'So_u10', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_o%flds(n1), compocn, compatm, mapconsf, 'ofrac', ocn2atm_fmapname) ! map ocn->atm
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds , 'Sx_u10', &
         merge_from1=complnd, merge_field1='Sl_u10', merge_type1='merge', merge_fracname1='lfrac', &
         merge_from2=compice, merge_field2='Si_u10', merge_type2='merge', merge_fracname2='ifrac', &
         merge_from3=compmed, merge_field3='So_u10', merge_type3='merge', merge_fracname3='ofrac')

    ! 'Zonal surface stress'
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds , 'Fall_taux', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1) , complnd, compatm, mapconsf, 'lfrin', lnd2atm_fmapname)
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds , 'Faii_taux', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n1) , compice, compatm, mapconsf, 'ifrac', ice2atm_fmapname)
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds , 'Fioi_taux', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n1) , compice, compocn, mapfcopy, 'unset', 'unset')
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_a%flds, 'Faox_taux', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_a%flds(n1), compatm, compocn, mapconsf, 'one'  , atm2ocn_fmapname) ! map atm->ocn
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_o%flds, 'Faox_taux', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_o%flds(n1), compocn, compatm, mapconsf, 'ofrac', ocn2atm_fmapname) ! map ocn->atm
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds , 'Faxx_taux', &
         merge_from1=complnd, merge_field1='Fall_taux', merge_type1='merge', merge_fracname1='lfrac', &
         merge_from2=compice, merge_field2='Faii_taux', merge_type2='merge', merge_fracname2='ifrac', &
         merge_from3=compmed, merge_field3='Faox_taux', merge_type3='merge', merge_fracname3='ofrac')
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds , 'Foxx_taux', &
         merge_from1=compmed, merge_field1='Faox_taux', merge_type1='merge', merge_fracname1='ofrac', &
         merge_from2=compice, merge_field2='Fioi_taux', merge_type2='merge', merge_fracname2='ifrac')

    ! 'Meridional surface stress'
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds , 'Fall_tauy', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1) , complnd, compatm, mapconsf, 'lfrin', lnd2atm_fmapname)
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds , 'Faii_tauy', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n1) , compice, compatm, mapconsf, 'ifrac', ice2atm_fmapname)
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds , 'Fioi_tauy', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n1) , compice, compocn, mapfcopy, 'unset', 'unset')
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_a%flds, 'Faox_tauy', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_a%flds(n1), compatm, compocn, mapconsf, 'one'  , atm2ocn_fmapname) ! map atm->ocn
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_o%flds, 'Faox_tauy', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_o%flds(n1), compocn, compatm, mapconsf, 'ofrac', ocn2atm_fmapname) ! map ocn->atm
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds , 'Faxx_tauy', &
         merge_from1=complnd, merge_field1='Fall_tauy', merge_type1='merge', merge_fracname1='lfrac', &
         merge_from2=compice, merge_field2='Faii_tauy', merge_type2='merge', merge_fracname2='ifrac', &
         merge_from3=compmed, merge_field3='Faox_tauy', merge_type3='merge', merge_fracname3='ofrac')
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Foxx_tauy', &
         merge_from1=compmed, merge_field1='Faox_tauy', merge_type1='merge', merge_fracname1='ofrac', &
         merge_from2=compice, merge_field2='Fioi_tauy', merge_type2='merge', merge_fracname2='ifrac')

    ! 'Surface latent heat flux'
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds , 'Fall_lat', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1) , complnd, compatm, mapconsf, 'lfrin', lnd2atm_fmapname)
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds , 'Faii_lat', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n1) , compice, compatm, mapconsf, 'ifrac', ice2atm_fmapname)
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_a%flds, 'Faox_lat', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_a%flds(n1), compatm, compocn, mapconsf, 'one'  , atm2ocn_fmapname) ! map atm->ocn
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_o%flds, 'Faox_lat', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_o%flds(n1), compocn, compatm, mapconsf, 'ofrac', ocn2atm_fmapname) ! map ocn->atm
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds , 'Faxx_lat', &
         merge_from1=complnd, merge_field1='Fall_lat', merge_type1='merge', merge_fracname1='lfrac', &
         merge_from2=compice, merge_field2='Faii_lat', merge_type2='merge', merge_fracname2='ifrac', &
         merge_from3=compmed, merge_field3='Faox_lat', merge_type3='merge', merge_fracname3='ofrac')
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds , 'Foxx_lat', &
         merge_from1=compmed, merge_field1='Faox_lat', merge_type1='merge', merge_fracname1='ofrac')

    ! 'Sensible heat flux'
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds , 'Fall_sen', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1) , complnd, compatm, mapconsf, 'lfrin', lnd2atm_fmapname)
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds , 'Faii_sen', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n1) , compice, compatm, mapconsf, 'ifrac', ice2atm_fmapname)
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_a%flds, 'Faox_sen', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_a%flds(n1), compatm, compocn, mapconsf, 'one'  , atm2ocn_fmapname) ! map atm->ocn
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_o%flds, 'Faox_sen', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_o%flds(n1), compocn, compatm, mapconsf, 'ofrac', ocn2atm_fmapname) ! map ocn->atm
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds , 'Faxx_sen', &
         merge_from1=complnd, merge_field1='Fall_sen', merge_type1='merge', merge_fracname1='lfrac', &
         merge_from2=compice, merge_field2='Faii_sen', merge_type2='merge', merge_fracname2='ifrac', &
         merge_from3=compmed, merge_field3='Faox_sen', merge_type3='merge', merge_fracname3='ofrac')
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds , 'Foxx_sen', &
         merge_from1=compmed, merge_field1='Faox_sen', merge_type1='merge', merge_fracname1='ofrac')

    ! 'Surface upward longwave heat flux'
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds , 'Fall_lwup', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1) , complnd, compatm, mapconsf, 'lfrin', lnd2atm_fmapname)
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds , 'Faii_lwup', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n1) , compice, compatm, mapconsf, 'ifrac', ice2atm_fmapname)
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_a%flds, 'Faox_lwup', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_a%flds(n1), compatm, compocn, mapconsf, 'one'  , atm2ocn_fmapname) ! map atm->ocn
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_o%flds, 'Faox_lwup', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_o%flds(n1), compocn, compatm, mapconsf, 'ofrac', ocn2atm_fmapname) ! map ocn->atm
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds , 'Faxx_lwup', &
         merge_from1=complnd, merge_field1='Fall_lwup', merge_type1='merge', merge_fracname1='lfrac', &
         merge_from2=compice, merge_field2='Faii_lwup', merge_type2='merge', merge_fracname2='ifrac', &
         merge_from3=compmed, merge_field3='Faox_lwup', merge_type3='merge', merge_fracname3='ofrac')
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds , 'Foxx_lwup', &
         merge_from1=compmed, merge_field1='Faox_lwup', merge_type1='merge', merge_fracname1='ofrac')

    ! 'Evaporation water flux'
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds , 'Fall_evap', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1) , complnd, compatm, mapconsf, 'lfrin', lnd2atm_fmapname)
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds , 'Faii_evap', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n1) , compice, compatm, mapconsf, 'ifrac', ice2atm_fmapname)
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_a%flds, 'Faox_evap', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_a%flds(n1), compatm, compocn, mapconsf, 'one'  , atm2ocn_fmapname) ! map atm->ocn
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_o%flds, 'Faox_evap', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_o%flds(n1), compocn, compatm, mapconsf, 'ofrac', ocn2atm_fmapname) ! map ocn->atm
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds , 'Faxx_evap', &
         merge_from1=complnd, merge_field1='Fall_evap', merge_type1='merge', merge_fracname1='lfrac', &
         merge_from2=compice, merge_field2='Faii_evap', merge_type2='merge', merge_fracname2='ifrac', &
         merge_from3=compmed, merge_field3='Faox_evap', merge_type3='merge', merge_fracname3='ofrac')
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds , 'Foxx_evap', &
         merge_from1=compmed, merge_field1='Faox_evap', merge_type1='merge', merge_fracname1='ofrac')

    ! 'Dust flux (particle bin number 1)'
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds, 'Fall_flxdst1', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1), complnd, compatm, mapconsf, 'lfrin', lnd2atm_fmapname)
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds, 'Fall_flxdst1', &
         merge_from1=complnd, merge_field1='Fall_flxdst1', merge_type1='copy_with_weights', merge_fracname1='lfrac')

    ! 'Dust flux (particle bin number 2)'
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds, 'Fall_flxdst2', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1), complnd, compatm, mapconsf, 'lfrin', lnd2atm_fmapname)
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds, 'Fall_flxdst2', &
         merge_from1=complnd, merge_field1='Fall_flxdst2', merge_type1='copy_with_weights', merge_fracname1='lfrac')

    ! 'Dust flux (particle bin number 3)'
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds, 'Fall_flxdst3', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1), complnd, compatm, mapconsf, 'lfrin', lnd2atm_fmapname)
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds, 'Fall_flxdst3', &
         merge_from1=complnd, merge_field1='Fall_flxdst3', merge_type1='copy_with_weights', merge_fracname1='lfrac')

    ! 'Dust flux (particle bin number 3)'
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds, 'Fall_flxdst4', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds, 'Fall_flxdst4', &
         merge_from1=complnd, merge_field1='Fall_flxdst4', merge_type1='copy_with_weights', merge_fracname1='lfrac')
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1), complnd, compatm, mapconsf, 'lfrin', lnd2atm_fmapname)

    !-----------------------------
    ! atm<->ocn only exchange
    !-----------------------------

    ! 'Sea level pressure'
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Sa_pslv', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapbilnr, 'one', atm2ocn_smapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapbilnr, 'one', atm2ocn_smapname)
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Sa_pslv', &
         merge_from1=compatm, merge_field1='Sa_pslv', merge_type1='copy')

    ! 'Wind speed squared at 10 meters'
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_a%flds, 'So_duu10n', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_a%flds(n1), compatm, compocn, mapbilnr, 'one'  , atm2ocn_fmapname) ! map atm->ocn
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_o%flds, 'So_duu10n')
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_o%flds(n1), compocn, compatm, mapconsf, 'ofrac', ocn2atm_fmapname) ! map ocn->atm
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds , 'So_duu10n', &
         merge_from1=compmed, merge_field1='So_duu10n', merge_type1='copy')

    ! 'Surface fraction velocity in ocean'
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_a%flds, 'So_ustar', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_a%flds(n1), compatm, compocn, mapbilnr, 'one'  , atm2ocn_fmapname) ! map atm->ocn
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_o%flds, 'So_ustar')
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_o%flds(n1), compocn, compatm, mapconsf, 'ofrac', ocn2atm_fmapname) ! map ocn->atm
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds , 'So_ustar', &
         merge_from1=compmed, merge_field1='So_ustar', merge_type1='copy')

    !-----------------------------
    ! ice->ocn exchange
    !-----------------------------

    ! 'Heat flux from melting'
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds, 'Fioi_melth', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Fioi_melth', &
         merge_from1=compice, merge_field1='Fioi_melth', merge_type1='copy_with_weights', merge_fracname1='ifrac')
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n1), compice, compocn,  mapfcopy, 'unset', 'unset')

    ! 'Water flux due to melting'
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds, 'Fioi_meltw', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Fioi_meltw', &
         merge_from1=compice, merge_field1='Fioi_meltw', merge_type1='copy_with_weights', merge_fracname1='ifrac')
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n1), compice, compocn,  mapfcopy, 'unset', 'unset')

    ! 'Salt flux'
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds, 'Fioi_salt', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Fioi_salt', &
         merge_from1=compice, merge_field1='Fioi_salt', merge_type1='copy_with_weights', merge_fracname1='ifrac')
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n1), compice, compocn,  mapfcopy, 'unset', 'unset')

    ! 'Hydrophylic black carbon deposition flux'
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds, 'Fioi_bcphi', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Fioi_bcphi', &
         merge_from1=compice, merge_field1='Fioi_bcphi', merge_type1='copy_with_weights', merge_fracname1='ifrac')
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n1), compice, compocn,  mapfcopy, 'unset', 'unset')

    ! 'Hydrophobic black carbon deposition flux'
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds, 'Fioi_bcpho', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Fioi_bcpho', &
         merge_from1=compice, merge_field1='Fioi_bcpho', merge_type1='copy_with_weights', merge_fracname1='ifrac')
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n1), compice, compocn,  mapfcopy, 'unset', 'unset')

    ! 'Dust flux'
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds, 'Fioi_flxdst', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Fioi_flxdst', &
         merge_from1=compice, merge_field1='Fioi_flxdst', merge_type1='copy_with_weights', merge_fracname1='ifrac')
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n1), compice, compocn,  mapfcopy, 'unset', 'unset')

    !-----------------------------
    ! ocn -> ice exchange (some of these fields are also used in the atm/ocn flux computation)
    !-----------------------------

    ! 'Sea surface salinity'
    call shr_nuopc_fldList_AddFld(fldListFr(compocn)%flds, 'So_s', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'So_s', &
         merge_from1=compocn, merge_field1='So_s', merge_type1='copy')
    call shr_nuopc_fldList_AddMap(fldListFr(compocn)%flds(n1), compocn, compice,  mapfcopy, 'unset', 'unset')

    ! 'Zonal sea water velocity'
    call shr_nuopc_fldList_AddFld(fldListFr(compocn)%flds, 'So_u', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'So_u', &
         merge_from1=compocn, merge_field1='So_u', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compwav)%flds, 'So_u', &
         merge_from1=compocn, merge_field1='So_u', merge_type1='copy')
    call shr_nuopc_fldList_AddMap(fldListFr(compocn)%flds(n1), compocn, compice,  mapfcopy, 'unset', 'unset')
    call shr_nuopc_fldList_AddMap(fldListFr(compocn)%flds(n1), compocn, compwav,  mapbilnr, 'one'  , 'ocn2wav_smapname')

    ! 'Fraction of sw penetrating surface layer for diurnal cycle'
    call shr_nuopc_fldList_AddFld(fldListFr(compocn)%flds , 'So_fswpen', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compocn)%flds(n1), compocn, compice,  mapfcopy, 'unset', 'unset')
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_a%flds, 'So_fswpen')
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_o%flds, 'So_fswpen')

    ! 'Meridional sea water velocity'
    call shr_nuopc_fldList_AddFld(fldListFr(compocn)%flds, 'So_v', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compocn)%flds(n1), compocn, compice,  mapfcopy, 'unset', 'unset')
    call shr_nuopc_fldList_AddMap(fldListFr(compocn)%flds(n1), compocn, compwav,  mapbilnr, 'one'  , 'ocn2wav_smapname')
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'So_v', &
         merge_from1=compocn, merge_field1='So_v', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compwav)%flds, 'So_v', &
         merge_from1=compocn, merge_field1='So_v', merge_type1='copy')

    ! 'Zonal sea surface slope'
    call shr_nuopc_fldList_AddFld(fldListFr(compocn)%flds, 'So_dhdx', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'So_dhdx', &
         merge_from1=compocn, merge_field1='So_dhdx', merge_type1='copy')
    call shr_nuopc_fldList_AddMap(fldListFr(compocn)%flds(n1), compocn, compice,  mapfcopy, 'unset', 'unset')

    ! 'Meridional sea surface slope'
    call shr_nuopc_fldList_AddFld(fldListFr(compocn)%flds, 'So_dhdy', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'So_dhdy', &
         merge_from1=compocn, merge_field1='So_dhdy', merge_type1='copy')
    call shr_nuopc_fldList_AddMap(fldListFr(compocn)%flds(n1), compocn, compice,  mapfcopy, 'unset', 'unset')

    ! 'Ocean Boundary Layer Depth'
    call shr_nuopc_fldList_AddFld(fldListFr(compocn)%flds, 'So_bldepth', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(compwav)%flds, 'So_bldepth', &
         merge_from1=compocn, merge_field1='So_bldepth', merge_type1='copy')
    call shr_nuopc_fldList_AddMap(fldListFr(compocn)%flds(n1), compocn, compwav,  mapbilnr, 'one', 'ocn2wav_smapname')

    ! 'Ocean melt and freeze potential'
    call shr_nuopc_fldList_AddFld(fldListFr(compocn)%flds, 'Fioo_q', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Fioo_q', &
         merge_from1=compocn, merge_field1='Fioo_q', merge_type1='copy')
    call shr_nuopc_fldList_AddMap(fldListFr(compocn)%flds(n1), compocn, compice,  mapfcopy, 'unset', 'unset')

    !-----------------------------
    ! lnd->rof exchange
    !-----------------------------

    ! 'Water flux from land (liquid surface)'
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds, 'Flrl_rofsur', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(comprof)%flds, 'Flrl_rofsur', &
         merge_from1=complnd, merge_field1='Flrl_rofsur', merge_type1='copy_with_weights', merge_fracname1='lfrac')
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1), complnd, comprof, mapconsf, 'lfrin', lnd2rof_fmapname)

    ! 'Water flux from land (liquid glacier, wetland, and lake)'
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds, 'Flrl_rofgwl', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(comprof)%flds, 'Flrl_rofgwl', &
         merge_from1=complnd, merge_field1='Flrl_rofgwl', merge_type1='copy_with_weights', merge_fracname1='lfrac')
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1), complnd, comprof, mapconsf, 'lfrin', lnd2rof_fmapname)

    ! 'Water flux from land (liquid subsurface)'
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds, 'Flrl_rofsub', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(comprof)%flds, 'Flrl_rofsub', &
         merge_from1=complnd, merge_field1='Flrl_rofsub', merge_type1='copy_with_weights', merge_fracname1='lfrac')
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1), complnd, comprof, mapconsf, 'lfrin', lnd2rof_fmapname)

    ! 'Water flux from land direct to ocean'
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds, 'Flrl_rofdto', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(comprof)%flds, 'Flrl_rofdto', &
         merge_from1=complnd, merge_field1='Flrl_rofdto', merge_type1='copy_with_weights', merge_fracname1='lfrac')
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1), complnd, comprof, mapconsf, 'lfrin', lnd2rof_fmapname)

    ! 'Water flux from land (frozen)'
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds, 'Flrl_rofi', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(comprof)%flds, 'Flrl_rofi', &
         merge_from1=complnd, merge_field1='Flrl_rofi', merge_type1='copy_with_weights', merge_fracname1='lfrac')
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1), complnd, comprof, mapconsf, 'lfrin', lnd2rof_fmapname)

    ! Irrigation flux (land/rof only)
    ! 'Irrigation flux (withdrawal from rivers)'
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds, 'Flrl_irrig', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(comprof)%flds, 'Flrl_irrig', &
         merge_from1=complnd, merge_field1='Flrl_irrig', merge_type1='copy_with_weights', merge_fracname1='lfrac')
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1), complnd, comprof, mapconsf, 'lfrin', lnd2rof_fmapname)

    !-----------------------------
    ! rof->lnd
    !-----------------------------

    ! 'Waterflux back to land due to flooding'
    call shr_nuopc_fldList_AddFld(fldListFr(comprof)%flds, 'Flrr_flood', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Flrr_flood', &
         merge_from1=comprof, merge_field1='Flrr_flood', merge_type1='copy')
    ! TODO: who should this be handled in terms of feeding back to the ocean
    !call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Flrr_flood') ! CUSTOM
    call shr_nuopc_fldList_AddMap(fldListFr(comprof)%flds(n1), comprof, complnd, mapconsf, 'one', rof2lnd_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(comprof)%flds(n1), comprof, compocn, mapconsf, 'one', rof2ocn_fmapname)

    ! 'River channel total water volume'
    call shr_nuopc_fldList_AddFld(fldListFr(comprof)%flds, 'Flrr_volr', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Flrr_volr', &
         merge_from1=comprof, merge_field1='Flrr_volr', merge_type1='copy')
    call shr_nuopc_fldList_AddMap(fldListFr(comprof)%flds(n1), comprof, complnd, mapconsf, 'one', rof2lnd_fmapname)

    ! 'River channel main channel water volume'
    call shr_nuopc_fldList_AddFld(fldListFr(comprof)%flds, 'Flrr_volrmch', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Flrr_volrmch', &
         merge_from1=comprof, merge_field1='Flrr_volrmch', merge_type1='copy')
    call shr_nuopc_fldList_AddMap(fldListFr(comprof)%flds(n1), comprof, complnd, mapconsf, 'one', rof2lnd_fmapname)

    !-----------------------------
    ! rof->ocn (liquid and frozen) and glc->ocn
    !-----------------------------

    ! 'glc liquid runoff flux to ocean'
    call shr_nuopc_fldList_AddFld(fldListFr(compglc)%flds, 'Fogg_rofl', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compglc)%flds(n1), compglc, compocn,  mapfiler, 'one', glc2ocn_liq_rmapname)

    ! 'Water flux into sea water due to runoff (liquid)'
    call shr_nuopc_fldList_AddFld(fldListFr(comprof)%flds, 'Forr_rofl', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(comprof)%flds(n1), comprof, compocn, mapfiler, 'none', rof2ocn_liq_rmapname)

    ! 'Total Water flux into sea water due to runoff (liquid)'
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Foxx_rofl', &
         merge_from1=comprof, merge_field1='Forr_rofl:Flrr_flood', merge_type1='accumulate', &
         merge_from2=compglc, merge_field2='Fogg_rofl'           , merge_type2='accumulate')

    ! 'glc frozen runoff flux to ocean'
    call shr_nuopc_fldList_AddFld(fldListFr(compglc)%flds, 'Fogg_rofi', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compglc)%flds(n1), compglc, compocn,  mapfiler, 'one', glc2ocn_ice_rmapname)

    ! Water flux into sea water due to runoff (frozen)'
    call shr_nuopc_fldList_AddFld(fldListFr(comprof)%flds, 'Forr_rofi', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(comprof)%flds(n1), comprof, compocn, mapfiler, 'none', rof2ocn_ice_rmapname)

    ! 'Total Water flux into sea water due to runoff (frozen)'
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Foxx_rofi', &
         merge_from1=comprof, merge_field1='Forr_rofi', merge_type1='accumulate', &
         merge_from2=compglc, merge_field2='Fogg_rofi', merge_type2='accumulate')

    !-----------------------------
    ! rof(frozen)->ice and glc->ice
    !-----------------------------

    ! 'Water flux into sea ice due to runoff (frozen)'
    call shr_nuopc_fldList_AddFld(fldListFr(comprof)%flds, 'Firr_rofi', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(comprof)%flds(n1), comprof, compice, mapfiler, 'none', rof2ocn_ice_rmapname)

    ! 'glc frozen runoff_iceberg flux to ice'
    call shr_nuopc_fldList_AddFld(fldListFr(compglc)%flds, 'Figg_rofi', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compglc)%flds(n1), compglc, compice,  mapfiler, 'one', glc2ice_rmapname)

    ! 'Total frozen water flux into sea ice '
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Fixx_rofi', &
         merge_from1=comprof, merge_field1='Firr_rofi', merge_type1='accumulate', &
         merge_from2=compglc, merge_field2='Figg_rofi', merge_type2='accumulate')

    !-----------------------------
    ! wav->ocn
    !-----------------------------

    !  'Langmuir multiplier'
    call shr_nuopc_fldList_AddFld(fldListFr(compwav)%flds, 'Sw_lamult', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Sw_lamult', &
         merge_from1=compwav, merge_field1='Sw_lamult', merge_type1='copy')
    call shr_nuopc_fldList_AddMap(fldListFr(compwav)%flds(n1), compwav, compocn,  mapbilnr, 'one', wav2ocn_smapname)

    !  'Stokes drift u component'
    call shr_nuopc_fldList_AddFld(fldListFr(compwav)%flds, 'Sw_ustokes', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Sw_ustokes', &
         merge_from1=compwav, merge_field1='Sw_ustokes', merge_type1='copy')
    call shr_nuopc_fldList_AddMap(fldListFr(compwav)%flds(n1), compwav, compocn,  mapbilnr, 'one', wav2ocn_smapname)

    !  'Stokes drift v component'
    call shr_nuopc_fldList_AddFld(fldListFr(compwav)%flds, 'Sw_vstokes', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Sw_vstokes', &
         merge_from1=compwav, merge_field1='Sw_vstokes', merge_type1='copy')
    call shr_nuopc_fldList_AddMap(fldListFr(compwav)%flds(n1), compwav, compocn, mapbilnr, 'one', wav2ocn_smapname)

    !  'Stokes drift depth'
    call shr_nuopc_fldList_AddFld(fldListFr(compwav)%flds, 'Sw_hstokes', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Sw_hstokes', &
         merge_from1=compwav, merge_field1='Sw_hstokes', merge_type1='copy')
    call shr_nuopc_fldList_AddMap(fldListFr(compwav)%flds(n1), compwav, compocn, mapbilnr, 'one', wav2ocn_smapname)

    !  'Downward solar radiation'
    call shr_nuopc_fldList_AddFld(FldListMed_ocnalb_o%flds, 'Faox_swdn', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListMed_ocnalb_o%flds(n1), compocn, compatm, mapconsf, 'ofrac', ocn2atm_smapname)

    !  'Upward solar radiation'
    call shr_nuopc_fldList_AddFld(FldListMed_ocnalb_o%flds, 'Faox_swup', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListMed_ocnalb_o%flds(n1), compocn, compatm, mapconsf, 'ofrac', ocn2atm_smapname)

    !-----------------------------
    ! glc -> ocn
    !-----------------------------

    !-----------------------------
    ! glc -> lnd
    !-----------------------------

    ! for glc fields with multiple elevation classes in glc->lnd
    ! fields from glc->med do NOT have elevation classes
    ! fields from med->lnd are BROKEN into multiple elevation classes

    !  'Ice sheet grid coverage on global grid'
    call shr_nuopc_fldList_AddFld(fldListFr(compglc)%flds   , 'Sg_icemask', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds   , 'Sg_icemask', &
         merge_from1=compglc, merge_field1='Sg_icemask', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListMed_x2l_fr_glc%flds, 'Sg_icemask') ! Needed for FB initialization
    call shr_nuopc_fldList_AddMap(fldListFr(compglc)%flds(n1), compglc, complnd,  mapconsf, 'one', glc2lnd_smapname)

    ! 'Ice sheet mask where we are potentially sending non-zero fluxes'
    call shr_nuopc_fldList_AddFld(fldListFr(compglc)%flds   , 'Sg_icemask_coupled_fluxes', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds   , 'Sg_icemask_coupled_fluxes', &
         merge_from1=compglc, merge_field1='Sg_icemask_coupled_fluxes', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListMed_x2l_fr_glc%flds, 'Sg_icemask_coupled_fluxes')
    call shr_nuopc_fldList_AddMap(fldListFr(compglc)%flds(n1), compglc, complnd,  mapconsf, 'one', glc2lnd_smapname)

    ! 'Fraction of glacier area'
    call shr_nuopc_fldList_AddFld(fldListFr(compglc)%flds, 'Sg_ice_covered', fldindex=n1)
    call shr_nuopc_fldList_AddMap(FldListFr(compglc)%flds(n1), compglc, complnd, mapconsf, 'unset', glc2lnd_fmapname) ! TODO: normalization?
    if (glc_nec > 0) then
       name = 'Sg_ice_covered'
       do num = 0, glc_nec
          cnum = glc_elevclass_as_string(num)
          call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Sg_ice_covered'//trim(cnum), &
               merge_from1=compglc, merge_field1='Sg_ice_covered'//trim(cnum), merge_type1='copy')
          call shr_nuopc_fldList_AddFld(fldListMed_x2l_fr_glc%flds, 'Sg_ice_covered'//trim(cnum))
       end do
    end if

    ! 'Surface height of glacier'
    call shr_nuopc_fldList_AddFld(fldListFr(compglc)%flds, 'Sg_topo', fldindex=n1)
    call shr_nuopc_fldList_AddMap(FldListFr(compglc)%flds(n1), compglc, compglc, mapconsf, 'custom', glc2lnd_fmapname)
    if (glc_nec > 0) then
       name = 'Sg_topo'
       do num = 0, glc_nec
          cnum = glc_elevclass_as_string(num)
          call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Sg_topo'//trim(cnum), &
               merge_from1=compglc, merge_field1='Sg_topo'//trim(cnum), merge_type1='copy')
          call shr_nuopc_fldList_AddFld(fldListMed_x2l_fr_glc%flds, 'Sg_topo'//trim(cnum))
       end do
    end if

    ! 'Downward heat flux from glacier interior'
    call shr_nuopc_fldList_AddFld(fldListFr(compglc)%flds, 'Flgg_hflx', fldindex=n1)
    call shr_nuopc_fldList_AddMap(FldListFr(compglc)%flds(n1), compglc, compglc, mapconsf, 'custom', glc2lnd_fmapname)
    if (glc_nec > 0) then
       name = 'Flgg_hflx'
       do num = 0, glc_nec
          cnum = glc_elevclass_as_string(num)
          call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Flgg_hflx'//trim(cnum), &
               merge_from1=compglc, merge_field1='Flgg_hflx'//trim(cnum), merge_type1='copy')
          call shr_nuopc_fldList_AddFld(fldListMed_x2l_fr_glc%flds, 'Flgg_hflx'//trim(cnum))
       end do
    end if

    !-----------------------------
    ! lnd -> glc
    !-----------------------------

    ! glc fields with multiple elevation classes: lnd->glc
    ! - fields sent from lnd->med are in multiple elevation classes
    ! - fields sent from med->glc do NOT have elevation classes
    ! - need to keep track of the l2x fields destined for glc in the
    !   additional variables, l2x_to_glc. This is needed so that can set up
    !   additional fields holding accumulated quantities of just these fields.

    ! Sets a coupling field for all glc elevation classes (1:glc_nec) plus bare land (index 0).
    ! Note that, if glc_nec = 0, then we don't create any coupling fields (not even the bare land (0) fldindex)

    ! 'New glacier ice flux'
    if (glc_nec > 0) then
       do num = 0, glc_nec
          name = 'Flgl_qice'
          cnum = glc_elevclass_as_string(num)
          call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds   , 'Flgl_qice'//trim(cnum), fldindex=n1)
          call shr_nuopc_fldList_AddFld(fldListMed_l2x_to_glc%flds, 'Flgl_qice'//trim(cnum))
       end do
    end if
    call shr_nuopc_fldList_AddFld(fldListTo(compglc)%flds, 'Flgl_qice')
    ! TODO: enter merging info
    call shr_nuopc_fldList_AddMap(FldListFr(complnd)%flds(n1), complnd, compglc, mapconsf, 'none', lnd2glc_fmapname)

    ! 'Surface temperature of glacier'
    if (glc_nec > 0) then
       name = 'Sl_tsrf'
       do num = 0, glc_nec
          cnum = glc_elevclass_as_string(num)
          call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds   , 'Sl_tsrf'//trim(cnum), fldindex=n1)
          call shr_nuopc_fldList_AddFld(fldListMed_l2x_to_glc%flds, 'Sl_tsrf'//trim(cnum))
       end do
    end if
    call shr_nuopc_fldList_AddFld(fldListTo(compglc)%flds, 'Sl_tsrf')
    ! TODO: enter merging info
    call shr_nuopc_fldList_AddMap(FldListFr(complnd)%flds(n1), complnd, compglc, mapconsf, 'none', lnd2glc_fmapname)

    ! Sl_topo is sent from lnd -> med, but is NOT sent to glc (it is only used for the remapping in the mediator)

    ! 'Surface height'
    if (glc_nec > 0) then
       name = 'Sl_topo'
       do num = 0, glc_nec
          cnum = glc_elevclass_as_string(num)
          call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds   , 'Sl_topo'//trim(cnum), fldindex=n1)
          call shr_nuopc_fldList_AddFld(fldListMed_l2x_to_glc%flds, 'Sl_topo'//trim(cnum))
       end do
    end if
    call shr_nuopc_fldList_AddFld(fldListTo(compglc)%flds, 'Sl_topo')
    call shr_nuopc_fldList_AddMap(FldListFr(complnd)%flds(n1), complnd, compglc, mapconsf, 'none', lnd2glc_fmapname)

    if (flds_co2a) then

       !  'Prognostic CO2 at the lowest model level'
       call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Sa_co2prog', fldindex=n1)
       call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Sa_co2prog', &
            merge_from1=compatm, merge_field1='Sa_co2prog', merge_type1='copy')
       call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Sa_co2prog', &
            merge_from1=compatm, merge_field1='Sa_co2prog', merge_type1='copy')
       call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapbilnr, 'one', atm2lnd_smapname)
       call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapbilnr, 'one', atm2ocn_smapname)

       !  'Diagnostic CO2 at the lowest model level'
       call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Sa_co2diag', fldindex=n1)
       call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Sa_co2diag', &
            merge_from1=compatm, merge_field1='Sa_co2diag', merge_type1='copy')
       call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Sa_co2diag', &
            merge_from1=compatm, merge_field1='Sa_co2diag', merge_type1='copy')
       call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapbilnr, 'one', atm2lnd_smapname)
       call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapbilnr, 'one', atm2ocn_smapname)

    else if (flds_co2b) then

       !  'Prognostic CO2 at the lowest model level'
       call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Sa_co2prog', fldindex=n1)
       call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Sa_co2prog', &
            merge_from1=compatm, merge_field1='Sa_co2prog', merge_type1='copy')
       call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapbilnr, 'one', atm2lnd_smapname)

       !  'Diagnostic CO2 at the lowest model level'
       call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Sa_co2diag', fldindex=n1)
       call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Sa_co2diag', &
            merge_from1=compatm, merge_field1='Sa_co2diag', merge_type1='copy')
       call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapbilnr, 'one', atm2lnd_smapname)

       !  'Surface flux of CO2 from land'
       call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds, 'Fall_fco2_lnd', fldindex=n1)
       call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds, 'Fall_fco2_lnd', &
            merge_from1=complnd, merge_field1='Fall_fco2_lnd', merge_type1='copy_with_weights', merge_fracname1='lfrac')
       call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1), complnd, compatm, mapconsf, 'one', atm2lnd_smapname)

    else if (flds_co2c) then

       !  'Prognostic CO2 at the lowest model level'
       call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Sa_co2prog', fldindex=n1)
       call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Sa_co2prog', &
            merge_from1=compatm, merge_field1='Sa_co2prog', merge_type1='copy')
       call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Sa_co2prog', &
            merge_from1=compatm, merge_field1='Sa_co2prog', merge_type1='copy')
       call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapbilnr, 'one', atm2lnd_smapname)
       call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapbilnr, 'one', atm2ocn_smapname)

       !  'Diagnostic CO2 at the lowest model level'
       call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Sa_co2diag', fldindex=n1)
       call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Sa_co2diag', &
            merge_from1=compatm, merge_field1='Sa_co2diag', merge_type1='copy')
       call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Sa_co2diag', &
            merge_from1=compatm, merge_field1='Sa_co2diag', merge_type1='copy')
       call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapbilnr, 'one', atm2lnd_smapname)
       call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapbilnr, 'one', atm2ocn_smapname)

       !  'Surface flux of CO2 from land'
       call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds, 'Fall_fco2_lnd', fldindex=n1)
       call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds, 'Fall_fco2_lnd', &
            merge_from1=complnd, merge_field1='Fall_fco2_lnd', merge_type1='copy_with_weights', merge_fracname1='lfrac')
       call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1), complnd, compatm, mapconsf, 'one', atm2lnd_smapname)

       ! 'Surface flux of CO2 from ocean'
       call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds, 'Faoo_fco2_ocn', fldindex=n1)
       call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds, 'Faoo_fco2_ocn') !CUSTOM
       call shr_nuopc_fldList_AddMap(fldListFr(compocn)%flds(n1), compocn, compatm, mapconsf, 'one', ocn2atm_smapname)
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

    !--------------------------------------------
    ! Atmospheric specific humidty at lowest level:
    !--------------------------------------------

    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Sa_shum_16O', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapbilnr, 'one', atm2lnd_smapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapbilnr, 'one', atm2ice_smapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapbilnr, 'one', atm2ocn_smapname)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Sa_shum_16O', &
         merge_from1=compatm, merge_field1='Sa_shum_16O', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Sa_shum_16O', &
         merge_from1=compatm, merge_field1='Sa_shum_16O', merge_type1='copy')

    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Sa_shum_18O', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapbilnr, 'one', atm2lnd_smapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapbilnr, 'one', atm2ice_smapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapbilnr, 'one', atm2ocn_smapname)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Sa_shum_18O', &
         merge_from1=compatm, merge_field1='Sa_shum_18O', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Sa_shum_18O', &
         merge_from1=compatm, merge_field1='Sa_shum_18O', merge_type1='copy')

    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Sa_shum_HDO', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapbilnr, 'one', atm2lnd_smapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapbilnr, 'one', atm2ice_smapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapbilnr, 'one', atm2ocn_smapname)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Sa_shum_HDO', &
         merge_from1=compatm, merge_field1='Sa_shum_HDO', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Sa_shum_HDO', &
         merge_from1=compatm, merge_field1='Sa_shum_HDO', merge_type1='copy')

    !--------------
    ! Isotopic surface snow water equivalent:
    !--------------

    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds, 'Sl_snowh_16O', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1), complnd, compatm, mapconsf, 'lfrin', lnd2atm_fmapname)
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds, 'Sl_snowh_16O', &
         merge_from1=complnd, merge_field1='Sl_snowh_16O', merge_type1='copy')

    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds, 'Sl_snowh_18O', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1), complnd, compatm, mapconsf, 'lfrin', lnd2atm_fmapname)
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds, 'Sl_snowh_18O', &
         merge_from1=complnd, merge_field1='Sl_snowh_18O', merge_type1='copy')

    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds, 'Sl_snowh_HDO', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1), complnd, compatm, mapconsf, 'lfrin', lnd2atm_fmapname)
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds, 'Sl_snowh_HDO', &
         merge_from1=complnd, merge_field1='Sl_snowh_HDO', merge_type1='copy')

    !--------------
    ! Isotopic Rain:
    !--------------

    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Faxa_rainc_16O', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapbilnr, 'one', atm2lnd_smapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapconsf, 'one', atm2ice_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapconsf, 'one', atm2ocn_fmapname)
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Faxa_rainl_16O', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapbilnr, 'one', atm2lnd_smapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapconsf, 'one', atm2ice_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapconsf, 'one', atm2ocn_fmapname)
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Faxa_rainl_HDO', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapbilnr, 'one', atm2lnd_smapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapconsf, 'one', atm2ice_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapconsf, 'one', atm2ocn_fmapname)

    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Faxa_rainc_18O', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapbilnr, 'one', atm2lnd_smapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapconsf, 'one', atm2ice_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapconsf, 'one', atm2ocn_fmapname)
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Faxa_rainl_18O', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapbilnr, 'one', atm2lnd_smapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapconsf, 'one', atm2ice_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapconsf, 'one', atm2ocn_fmapname)
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Faxa_rainc_HDO', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapbilnr, 'one', atm2lnd_smapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapconsf, 'one', atm2ice_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapconsf, 'one', atm2ocn_fmapname)

    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Faxa_rainl_16O', &
         merge_from1=compatm, merge_field1='Faxa_rainl_16O', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Faxa_rainl_18O', &
         merge_from1=compatm, merge_field1='Faxa_rainl_18O', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Faxa_rainl_HDO', &
         merge_from1=compatm, merge_field1='Faxa_rainl_HDO', merge_type1='copy')

    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Faxa_rainc_16O', &
         merge_from1=compatm, merge_field1='Faxa_rainc_16O', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Faxa_rainc_18O', &
         merge_from1=compatm, merge_field1='Faxa_rainc_18O', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Faxa_rainc_HDO', &
         merge_from1=compatm, merge_field1='Faxa_rainc_HDO', merge_type1='copy')

    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Faxa_rain_16O', &
         merge_from1=compatm, merge_field1='Faxa_rainc_16O:Faxa_rainl_16O',&
         merge_type1='accumulate', merge_fracname1='ofrac')
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Faxa_rain_18O', &
         merge_from1=compatm, merge_field1='Faxa_rainc_18O:Faxa_rainl_18O', &
         merge_type1='accumulate')
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Faxa_rain_HDO', &
         merge_from1=compatm, merge_field1='Faxa_rainc_HDO:Faxa_rainl_HDO', &
         merge_type1='accumulate')

    !-------------
    ! Isotopic snow:
    !-------------

    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Faxa_snowl_16O', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapbilnr, 'one', atm2lnd_smapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapconsf, 'one', atm2ice_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapconsf, 'one', atm2ocn_fmapname)
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Faxa_snowl_18O', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapbilnr, 'one', atm2lnd_smapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapconsf, 'one', atm2ice_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapconsf, 'one', atm2ocn_fmapname)
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Faxa_snowl_HDO', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapbilnr, 'one', atm2lnd_smapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapconsf, 'one', atm2ice_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapconsf, 'one', atm2ocn_fmapname)

    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Faxa_snowl_16O', &
         merge_from1=compatm, merge_field1='Faxa_snowl_16O', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Faxa_snowl_18O', &
         merge_from1=compatm, merge_field1='Faxa_snowl_18O', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Faxa_snowl_HDO', &
         merge_from1=compatm, merge_field1='Faxa_snowl_HDO', merge_type1='copy')

    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Faxa_snowc_16O', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapbilnr, 'one', atm2lnd_smapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapconsf, 'one', atm2ice_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapconsf, 'one', atm2ocn_fmapname)
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Faxa_snowc_18O', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapbilnr, 'one', atm2lnd_smapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapconsf, 'one', atm2ice_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapconsf, 'one', atm2ocn_fmapname)
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Faxa_snowc_HDO', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapbilnr, 'one', atm2lnd_smapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapconsf, 'one', atm2ice_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapconsf, 'one', atm2ocn_fmapname)

    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Faxa_snowc_16O', &
         merge_from1=compatm, merge_field1='Faxa_snowc_16O', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Faxa_snowc_18O', &
         merge_from1=compatm, merge_field1='Faxa_snowc_18O', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Faxa_snowc_HDO', &
         merge_from1=compatm, merge_field1='Faxa_snowc_HDO', merge_type1='copy')

    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Faxa_snow_16O', &
         merge_from1=compatm, merge_field1='Faxa_snowc_16O:Faxa_snowl_16O',&
         merge_type1='accumulate', merge_fracname1='ofrac')
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Faxa_snow_18O', &
         merge_from1=compatm, merge_field1='Faxa_snowc_18O:Faxa_snowl_18O', &
         merge_type1='accumulate')
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Faxa_snow_HDO', &
         merge_from1=compatm, merge_field1='Faxa_snowc_HDO:Faxa_snowl_HDO', &
         merge_type1='accumulate')

    !----------------------------------
    ! Isotopic precipitation (snow+snow):
    !----------------------------------

    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Faxa_prec_16O', &
         merge_from1=compatm, merge_field1='Faxa_rainc_16O:Faxa_snowc_16O:Faxa_rainl_16O:Faxa_snowl_16O', &
         merge_type1='accumulate', merge_fracname1='ofrac')

    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Faxa_prec_18O', &
         merge_from1=compatm, merge_field1='Faxa_rainc_18O:Faxa_snowc_18O:Faxa_rainl_18O:Faxa_snowl_18O', &
         merge_type1='accumulate', merge_fracname1='ofrac')

    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Faxa_prec_HDO', &
         merge_from1=compatm, merge_field1='Faxa_rainc_HDO:Faxa_snowc_HDO:Faxa_rainl_HDO:Faxa_snowl_HDO', &
         merge_type1='accumulate', merge_fracname1='ofrac')

    !-------------------------------------
    ! Isotopic two meter reference humidity:
    !-------------------------------------

    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds , 'Sl_qref_16O', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1) , complnd, compatm, mapconsf, 'lfrin', lnd2atm_fmapname)
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds , 'Si_qref_16O', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n1) , compice, compatm, mapconsf, 'ifrac', ice2atm_fmapname)
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_a%flds, 'So_qref_16O', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_a%flds(n1), compatm, compocn, mapbilnr, 'one'  , atm2ocn_fmapname) ! map atm->ocn
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_o%flds, 'So_qref_16O', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_o%flds(n1), compocn, compatm, mapconsf, 'ofrac', ocn2atm_fmapname) ! map ocn->atm
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds , 'Sx_qref_16O', &
         merge_from1=complnd, merge_field1='Sl_qref_16O', merge_type1='merge', merge_fracname1='lfrac', &
         merge_from2=compice, merge_field2='Si_qref_16O', merge_type2='merge', merge_fracname2='ifrac', &
         merge_from3=compmed, merge_field3='So_qref_16O', merge_type3='merge', merge_fracname3='ofrac')

    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds , 'Sl_qref_18O', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1) , complnd, compatm, mapconsf, 'lfrin', lnd2atm_fmapname)
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds , 'Si_qref_18O', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n1) , compice, compatm, mapconsf, 'ifrac', ice2atm_fmapname)
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_a%flds, 'So_qref_18O', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_a%flds(n1), compatm, compocn, mapbilnr, 'one'  , atm2ocn_fmapname) ! map atm->ocn
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_o%flds, 'So_qref_18O', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_o%flds(n1), compocn, compatm, mapconsf, 'ofrac', ocn2atm_fmapname) ! map ocn->atm
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds , 'Sx_qref_18O', &
         merge_from1=complnd, merge_field1='Sl_qref_18O', merge_type1='merge', merge_fracname1='lfrac', &
         merge_from2=compice, merge_field2='Si_qref_18O', merge_type2='merge', merge_fracname2='ifrac', &
         merge_from3=compmed, merge_field3='So_qref_18O', merge_type3='merge', merge_fracname3='ofrac')

    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds , 'Sl_qref_HDO', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1) , complnd, compatm, mapconsf, 'lfrin', lnd2atm_fmapname)
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds , 'Si_qref_HDO', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n1) , compice, compatm, mapconsf, 'ifrac', ice2atm_fmapname)
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_a%flds, 'So_qref_HDO', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_a%flds(n1), compatm, compocn, mapbilnr, 'one'  , atm2ocn_fmapname) ! map atm->ocn
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_o%flds, 'So_qref_HDO', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_o%flds(n1), compocn, compatm, mapconsf, 'ofrac', ocn2atm_fmapname) ! map ocn->atm
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds , 'Sx_qref_HDO', &
         merge_from1=complnd, merge_field1='Sl_qref_HDO', merge_type1='merge', merge_fracname1='lfrac', &
         merge_from2=compice, merge_field2='Si_qref_HDO', merge_type2='merge', merge_fracname2='ifrac', &
         merge_from3=compmed, merge_field3='So_qref_HDO', merge_type3='merge', merge_fracname3='ofrac')

    !-------------------------
    ! Isotopic Evaporation flux:
    !-------------------------

    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds , 'Fall_evap_16O', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1) , complnd, compatm, mapconsf, 'lfrin', lnd2atm_fmapname)
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds , 'Faii_evap_16O', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n1) , compice, compatm, mapconsf, 'ifrac', ice2atm_fmapname)
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_a%flds, 'Faox_evap_16O', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_a%flds(n1), compatm, compocn, mapconsf, 'one'  , atm2ocn_fmapname) ! map atm->ocn
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_o%flds, 'Faox_evap_16O', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_o%flds(n1), compocn, compatm, mapconsf, 'ofrac', ocn2atm_fmapname) ! map ocn->atm
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds , 'Faxx_evap_16O', &
         merge_from1=complnd, merge_field1='Fall_evap_16O', merge_type1='merge', merge_fracname1='lfrac', &
         merge_from2=compice, merge_field2='Faii_evap_16O', merge_type2='merge', merge_fracname2='ifrac', &
         merge_from3=compmed, merge_field3='Faox_evap_16O', merge_type3='merge', merge_fracname3='ofrac')
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds , 'Foxx_evap_16O', &
         merge_from1=compmed, merge_field1='Faox_evap_16O', merge_type1='merge', merge_fracname1='ofrac')

    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds , 'Fall_evap_18O', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1) , complnd, compatm, mapconsf, 'lfrin', lnd2atm_fmapname)
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds , 'Faii_evap_18O', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n1) , compice, compatm, mapconsf, 'ifrac', ice2atm_fmapname)
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_a%flds, 'Faox_evap_18O', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_a%flds(n1), compatm, compocn, mapconsf, 'one'  , atm2ocn_fmapname) ! map atm->ocn
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_o%flds, 'Faox_evap_18O', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_o%flds(n1), compocn, compatm, mapconsf, 'ofrac', ocn2atm_fmapname) ! map ocn->atm
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds , 'Faxx_evap_18O', &
         merge_from1=complnd, merge_field1='Fall_evap_18O', merge_type1='merge', merge_fracname1='lfrac', &
         merge_from2=compice, merge_field2='Faii_evap_18O', merge_type2='merge', merge_fracname2='ifrac', &
         merge_from3=compmed, merge_field3='Faox_evap_18O', merge_type3='merge', merge_fracname3='ofrac')
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds , 'Foxx_evap_18O', &
         merge_from1=compmed, merge_field1='Faox_evap_18O', merge_type1='merge', merge_fracname1='ofrac')

    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds , 'Fall_evap_HDO', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1) , complnd, compatm, mapconsf, 'lfrin', lnd2atm_fmapname)
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds , 'Faii_evap_HDO', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n1) , compice, compatm, mapconsf, 'ifrac', ice2atm_fmapname)
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_a%flds, 'Faox_evap_HDO', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_a%flds(n1), compatm, compocn, mapconsf, 'one'  , atm2ocn_fmapname) ! map atm->ocn
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_o%flds, 'Faox_evap_HDO', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_o%flds(n1), compocn, compatm, mapconsf, 'ofrac', ocn2atm_fmapname) ! map ocn->atm
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds , 'Faxx_evap_HDO', &
         merge_from1=complnd, merge_field1='Fall_evap_HDO', merge_type1='merge', merge_fracname1='lfrac', &
         merge_from2=compice, merge_field2='Faii_evap_HDO', merge_type2='merge', merge_fracname2='ifrac', &
         merge_from3=compmed, merge_field3='Faox_evap_HDO', merge_type3='merge', merge_fracname3='ofrac')
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds , 'Foxx_evap_HDO', &
         merge_from1=compmed, merge_field1='Faox_evap_HDO', merge_type1='merge', merge_fracname1='ofrac')

    !-----------------------------
    ! Isotopic sea ice melting flux:
    !-----------------------------

    ! 'Heat flux from melting'

    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds, 'Fioi_melth_16O', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Fioi_melth_16O', &
         merge_from1=compice, merge_field1='Fioi_melth_16O', merge_type1='copy_with_weights', merge_fracname1='ifrac')
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n1), compice, compocn,  mapfcopy, 'unset', 'unset')

    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds, 'Fioi_melth_18O', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Fioi_melth_18O', &
         merge_from1=compice, merge_field1='Fioi_melth_18O', merge_type1='copy_with_weights', merge_fracname1='ifrac')
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n1), compice, compocn,  mapfcopy, 'unset', 'unset')

    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds, 'Fioi_melth_HDO', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Fioi_melth_HDO', &
         merge_from1=compice, merge_field1='Fioi_melth_HDO', merge_type1='copy_with_weights', merge_fracname1='ifrac')
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n1), compice, compocn,  mapfcopy, 'unset', 'unset')

    ! 'Water flux due to melting'

    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds, 'Fioi_meltw_16O', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Fioi_meltw_16O', &
         merge_from1=compice, merge_field1='Fioi_meltw_16O', merge_type1='copy_with_weights', merge_fracname1='ifrac')
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n1), compice, compocn,  mapfcopy, 'unset', 'unset')

    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds, 'Fioi_meltw_18O', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Fioi_meltw_18O', &
         merge_from1=compice, merge_field1='Fioi_meltw_18O', merge_type1='copy_with_weights', merge_fracname1='ifrac')
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n1), compice, compocn,  mapfcopy, 'unset', 'unset')

    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds, 'Fioi_meltw_HDO', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Fioi_meltw_HDO', &
         merge_from1=compice, merge_field1='Fioi_meltw_HDO', merge_type1='copy_with_weights', merge_fracname1='ifrac')
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n1), compice, compocn,  mapfcopy, 'unset', 'unset')

    !-----------------------------
    !  Isotopic Runoff (r2o, l2x, x2r)
    !-----------------------------

    !    units    = 'kg m-2 s-1'
    !    longname = 'H2_16O Water flux from land (frozen)'
    !    stdname  = 'H2_16O_frozen_water_flux_into_runoff'
    !    call fld_add(flds_l2x, 'Flrl_rofi_16O')
    !    call fld_add(flds_x2r, flds_x2r_map,'Frxx_rofi_16O')

    !    longname = 'H2_18O Water flux from land (frozen)'
    !    stdname  = 'H2_18O_frozen_water_flux_into_runoff'
    !    call fld_add(flds_l2x, 'Flrl_rofi_18O')
    !    call fld_add(flds_x2r, flds_x2r_map,'Frxx_rofi_18O')

    !    longname = 'HDO Water flux from land (frozen)'
    !    stdname  = 'HDO_frozen_water_flux_into_runoff'
    !    call fld_add(flds_l2x, 'Flrl_rofi_HDO')
    !    call fld_add(flds_x2r, flds_x2r_map,'Frxx_rofi_HDO')

    !    longname = 'H2_16O Water flux from land (liquid)'
    !    stdname  = 'H2_16O_liquid_water_flux_into_runoff'
    !    call fld_add(flds_l2x, 'Flrl_rofl_16O')
    !    call fld_add(flds_x2r, flds_x2r_map,'Frxx_rofl_16O')

    !    longname = 'H2_18O Water flux from land (liquid)'
    !    stdname  = 'H2_18O_liquid_water_flux_into_runoff'
    !    call fld_add(flds_l2x, 'Flrl_rofl_18O')
    !    call fld_add(flds_x2r, flds_x2r_map,'Frxx_rofl_18O')

    !    longname = 'HDO Water flux from land (liquid)'
    !    stdname  = 'HDO_liquid_water_flux_into_runoff'
    !    call fld_add(flds_l2x, 'Flrl_rofl_HDO')
    !    call fld_add(flds_x2r, flds_x2r_map,'Frxx_rofl_HDO')

    !-----------------------------
    ! Isotopic r2x, x2o
    !-----------------------------

    !    units    = 'kg m-2 s-1'
    !    longname = 'H2_16O Water flux due to liq runoff '
    !    stdname  = 'H2_16O_water_flux_into_sea_water'
    !    call fld_add(flds_r2x, flds_r2x_map,'Forr_rofl_16O')
    !    call fld_add(flds_x2o, flds_x2o_map,'Foxx_rofl_16O')
    !    call fld_add(flds_r2o_liq, flds_r2o_liq_map,'Forr_rofl_16O')

    !    longname = 'H2_18O Water flux due to liq runoff '
    !    stdname  = 'H2_18O_water_flux_into_sea_water'
    !    call fld_add(flds_r2x, flds_r2x_map,'Forr_rofl_18O')
    !    call fld_add(flds_x2o, flds_x2o_map,'Foxx_rofl_18O')
    !    call fld_add(flds_r2o_liq, flds_r2o_liq_map,'Forr_rofl_18O')

    !    longname = 'HDO Water flux due to liq runoff '
    !    stdname  = 'HDO_water_flux_into_sea_water'
    !    call fld_add(flds_r2x, flds_r2x_map,'Forr_rofl_HDO')
    !    call fld_add(flds_x2o, flds_x2o_map,'Foxx_rofl_HDO')
    !    call fld_add(flds_r2o_liq, flds_r2o_liq_map,'Forr_rofl_HDO')

    !    longname = 'H2_16O Water flux due to ice runoff '
    !    stdname  = 'H2_16O_water_flux_into_sea_water'
    !    call fld_add(flds_r2x, flds_r2x_map,'Forr_rofi_16O')
    !    call fld_add(flds_x2o, flds_x2o_map,'Foxx_rofi_16O')
    !    call fld_add(flds_r2o_ice, flds_r2o_ice_map,'Forr_rofi_16O')

    !    longname = 'H2_18O Water flux due to ice runoff '
    !    stdname  = 'H2_18O_water_flux_into_sea_water'
    !    call fld_add(flds_r2x, flds_r2x_map,'Forr_rofi_18O')
    !    call fld_add(flds_x2o, flds_x2o_map,'Foxx_rofi_18O')
    !    call fld_add(flds_r2o_ice, flds_r2o_ice_map,'Forr_rofi_18O')

    !    longname = 'HDO Water flux due to ice runoff '
    !    stdname  = 'HDO_water_flux_into_sea_water'
    !    call fld_add(flds_r2x, flds_r2x_map,'Forr_rofi_HDO')
    !    call fld_add(flds_x2o, flds_x2o_map,'Foxx_rofi_HDO')
    !    call fld_add(flds_r2o_ice, flds_r2o_ice_map,'Forr_rofi_HDO')

    !    ! r2x, x2l

    !    units    = 'kg m-2 s-1'
    !    longname = 'H2_16O waterrflux due to flooding'
    !    stdname  = 'H2_16O_flodding_water_flux_back_to_land'
    !    call fld_add(flds_r2x, flds_r2x_map,'Flrr_flood_16O')
    !    call fld_add(flds_x2l, flds_x2l_map,'Flrr_flood_16O')

    !    longname = 'H2_18O waterrflux due to flooding'
    !    stdname  = 'H2_18O_flodding_water_flux_back_to_land'
    !    call fld_add(flds_r2x, flds_r2x_map,'Flrr_flood_18O')
    !    call fld_add(flds_x2l, flds_x2l_map,'Flrr_flood_18O')

    !    longname = 'HDO Waterrflux due to flooding'
    !    stdname  = 'HDO_flodding_water_flux_back_to_land'
    !    call fld_add(flds_r2x, flds_r2x_map,'Flrr_flood_HDO')
    !    call fld_add(flds_x2l, flds_x2l_map,'Flrr_flood_HDO')

    !    longname = 'H2_16O river channel water volume '
    !    stdname  = 'H2_16O_rtm_volr'
    !    call fld_add(flds_r2x, flds_r2x_map,'Flrr_volr_16O')
    !    call fld_add(flds_x2l, flds_x2l_map,'Flrr_volr_16O')

    !    longname = 'H2_18O river channel water volume '
    !    stdname  = 'H2_18O_rtm_volr'
    !    call fld_add(flds_r2x, flds_r2x_map,'Flrr_volr_18O')
    !    call fld_add(flds_x2l, flds_x2l_map,'Flrr_volr_18O')

    !    longname = 'HDO river channel water volume '
    !    stdname  = 'HDO_rtm_volr'
    !    call fld_add(flds_r2x, flds_r2x_map,'Flrr_volr_HDO')
    !    call fld_add(flds_x2l, flds_x2l_map,'Flrr_volr_HDO')

    !    ! longname = 'H2_18O Waterrflux due to flooding'
    !    ! stdname  = 'H2_18O_flodding_water_flux_back_to_land'
    !    ! call fld_add(flds_r2x, flds_r2x_map,'Flrr_flood_HDO')
    !    ! call fld_add(flds_x2l, flds_x2l_map,'Flrr_flood_HDO')

    !-----------------------------------------------------------------------------
    ! optional per thickness category fields
    !-----------------------------------------------------------------------------

    if (flds_i2o_per_cat) then
       do num = 1, ice_ncat
          write(cnum,'(i2.2)') num

          ! 'fractional ice coverage wrt ocean for thickness category ' // cnum
          name = 'Si_ifrac_' // cnum
          call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds, trim(name), fldindex=n1)
          call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, trim(name), &
               merge_from1=compice, merge_field1=trim(name), merge_type1='copy')
          call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n1), compice, compocn,  mapfcopy, 'unset', 'unset')

          ! Net shortwave radiation
          name = 'PFioi_swpen_ifrac_' // cnum
          call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds, trim(name), fldindex=n1)
          call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, trim(name), &
               merge_from1=compice, merge_field1=trim(name), merge_type1='copy')
          call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n1), compice, compocn,  mapfcopy, 'unset', 'unset')
       end do

       ! 'fractional atmosphere coverage wrt ocean'
       call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Sf_afrac')
       ! TODO: add mapping and merging

       ! 'fractional atmosphere coverage used in radiation computations wrt ocean'
       call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Sf_afracr')
       ! TODO: add mapping and merging

       ! 'net shortwave radiation times atmosphere fraction'
       call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Foxx_swnet_afracr')
       ! TODO: add mapping and merging
    end if

    !-----------------------------------------------------------------------------
    ! CARMA fields (volumetric soil water)
    !-----------------------------------------------------------------------------

    ! TODO: add this
    ! if (carma_fields /= ' ') then
    !    do n = 1,shr_string_listGetNum(carma_fields)
    !       call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds, trim(fldname), fldindex=n1)
    !       call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1), complnd, compatm, mapconsf, 'one', lnd2atm_smapname)
    !       call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds, trim(fldname), &
    !            merge_from1=complnd, merge_field1=trim(fldname), merge_type1='copy')
    !    enddo
    ! endif

    !-----------------------------------------------------------------------------
    ! MEGAN emissions fluxes fields
    !-----------------------------------------------------------------------------

    do num = 1, max_megan
       write(cnum,'(i3.3)') num
       fldname = 'Fall_voc' // cnum
       call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds, trim(fldname), fldindex=n1)
       call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1), complnd, compatm, mapconsf, 'one', atm2lnd_fmapname)
       call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds, trim(fldname), &
            merge_from1=complnd, merge_field1=trim(fldname), merge_type1='merge', merge_fracname1='lfrac')
    end do

    !-----------------------------------------------------------------------------
    ! Fire emissions fluxes fields
    !-----------------------------------------------------------------------------

    ! 'wild fire emission fluxes'
    do num = 1, max_fire
       write(cnum,'(i2.2)') num
       fldname  = 'Fall_fire' // cnum
       call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds, trim(fldname), fldindex=n1)
       call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1), complnd, compatm, mapconsf, 'one', lnd2atm_fmapname)
       call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds, trim(fldname), &
            merge_from1=complnd, merge_field1=trim(fldname), merge_type1='merge', merge_fracname1='lfrac')
    end do

    ! 'wild fire plume height'
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds, 'Sl_fztop', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1), complnd, compatm, mapconsf, 'one', lnd2atm_smapname)
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds, 'Sl_fztop', &
         merge_from1=complnd, merge_field1=trim(fldname), merge_type1='copy')

    !-----------------------------------------------------------------------------
    ! Dry Deposition fields
    !-----------------------------------------------------------------------------

    do num = 1, max_ddep
       write(cnum,'(i2.2)') num
       fldname  = 'Sl_dd' // cnum
       call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds, trim(fldname), fldindex=n1)
       call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1), complnd, compatm, mapconsf, 'one', lnd2atm_smapname)
       call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds, trim(fldname), &
            merge_from1=complnd, merge_field1=trim(fldname), merge_type1='copy')
    end do

    !-----------------------------------------------------------------------------
    ! Nitrogen Deposition fields
    !-----------------------------------------------------------------------------

    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Faxa_noy', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapbilnr, 'one', atm2lnd_smapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapbilnr, 'one', atm2ocn_smapname)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Faxa_noy', &
         merge_from1=compatm, merge_field1=trim(fldname), merge_type1='copy_with_weights', merge_fracname1='ofrac')
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Faxa_noy', &
         merge_from1=compatm, merge_field1=trim(fldname), merge_type1='copy_with_weights', merge_fracname1='ofrac')

    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Faxa_nhx', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapbilnr, 'one', atm2lnd_smapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapbilnr, 'one', atm2ocn_smapname)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Faxa_nhx', &
         merge_from1=compatm, merge_field1=trim(fldname), merge_type1='copy_with_weights', merge_fracname1='ofrac')
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Faxa_nhx', &
         merge_from1=compatm, merge_field1=trim(fldname), merge_type1='copy_with_weights', merge_fracname1='ofrac')

  end subroutine esmFlds_Init

end module esmFlds
