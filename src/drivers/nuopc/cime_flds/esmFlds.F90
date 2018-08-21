module esmFlds

  !---------------------------------------------------------------------
  ! This is a mediator specific routine that determines ALL possible
  ! fields exchanged between components and their associated routing,
  ! mapping and merging
  !---------------------------------------------------------------------

  use ESMF
  use NUOPC
  use med_constants_mod     , only : IN, R8, I8, CXX, CX, CS, CL
  use shr_nuopc_scalars_mod , only : flds_scalar_name
  use shr_nuopc_fldList_mod , only : mapbilnr, mapconsf, mapconsd, mappatch, mapfcopy, mapunset, mapfiler
  use shr_nuopc_fldList_mod , only : shr_nuopc_fldList_type
  use shr_nuopc_fldList_mod , only : shr_nuopc_fldList_AddMetaData
  use shr_nuopc_fldList_mod , only : shr_nuopc_fldList_AddFld
  use shr_nuopc_fldList_mod , only : shr_nuopc_fldList_AddMap
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_chkerr
  use seq_drydep_mod        , only : seq_drydep_init, seq_drydep_readnl, lnd_drydep
  use shr_megan_mod         , only : shr_megan_readnl, shr_megan_mechcomps_n
  use shr_fire_emis_mod     , only : shr_fire_emis_readnl, shr_fire_emis_mechcomps_n, shr_fire_emis_ztop_token
  use shr_carma_mod         , only : shr_carma_readnl
  use shr_ndep_mod          , only : shr_ndep_readnl
  use shr_flds_mod          , only : shr_flds_dom_coord, shr_flds_dom_other
  use shr_string_mod        , only : shr_string_listGetNum, shr_string_listGetName
  use glc_elevclass_mod     , only : glc_elevclass_as_string

  implicit none
  public

  ! IMPORTANT:
  ! TODO: the call to add metadata needs to be done on all pes
  ! all other calls are just on the mediator pes

  !----------------------------------------------------------------------------
  ! routines
  !----------------------------------------------------------------------------

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

  ! Initialize mediator field bundles in med.F90 routine DataInitialize

  type (shr_nuopc_fldList_type) :: fldListMed_aoflux_a
  type (shr_nuopc_fldList_type) :: fldListMed_aoflux_o
  type (shr_nuopc_fldList_type) :: fldListMed_ocnalb_o
  type (shr_nuopc_fldList_type) :: fldListMed_l2x_to_glc
  type (shr_nuopc_fldList_type) :: fldListMed_x2l_fr_glc
  type (shr_nuopc_fldList_type) :: fldListMed_g2x_to_lnd

  ! The following is used in the advertise AND realize fields in the components and mediator

  type (shr_nuopc_fldList_type) :: fldListTo(ncomps)
  type (shr_nuopc_fldList_type) :: fldListFr(ncomps)

  !----------------------------------------------------------------------------
  ! other field lists - column deliminated string
  !----------------------------------------------------------------------------

  character(len=CXX) :: drydep_fields       ! List of dry-deposition fields
  character(len=CXX) :: megan_voc_fields    ! List of MEGAN VOC emission fields
  character(len=CXX) :: fire_emis_fields    ! List of fire emission fields
  character(len=CX)  :: carma_fields        ! List of CARMA fields from lnd->atm
  character(len=CX)  :: ndep_fields         ! List of nitrogen deposition fields from atm->lnd/ocn
  integer            :: ice_ncat            ! number of sea ice thickness categories
  logical            :: add_ndep_fields     ! .true. => add ndep fields

  integer     , parameter :: CSS = 256           ! use longer short character
  character(*), parameter :: u_FILE_u = __FILE__

!================================================================================
contains
!================================================================================

  subroutine esmFlds_Init(gcomp, rc)

    ! input/output parameters:
    type(ESMF_GridComp)    :: gcomp
    integer, intent(inout) :: rc

    ! local variables:
    integer                :: dbrc
    type(ESMF_VM)          :: vm
    integer                :: localPet
    logical                :: mastertask
    character(len=CSS)     :: attname
    character(len=CSS)     :: units
    character(len=CSS)     :: longname
    character(len=CSS)     :: stdname
    integer                :: num, i, n
    integer                :: n1, n2, n3, n4
    character(len=  2)     :: cnum
    character(len=CSS)     :: name, fldname
    character(len=CL)      :: cvalue
    character(len=CX)      :: atm2ice_fmapname
    character(len=CX)      :: atm2ice_smapname
    character(len=CX)      :: atm2ice_vmapname
    character(len=CX)      :: atm2lnd_fmapname
    character(len=CX)      :: atm2lnd_smapname
    character(len=CX)      :: atm2ocn_fmapname
    character(len=CX)      :: atm2ocn_smapname
    character(len=CX)      :: atm2ocn_vmapname
    character(len=CX)      :: atm2wav_smapname
    character(len=CX)      :: glc2lnd_fmapname
    character(len=CX)      :: glc2lnd_smapname
    character(len=CX)      :: glc2ice_rmapname
    character(len=CX)      :: glc2ocn_liq_rmapname
    character(len=CX)      :: glc2ocn_ice_rmapname
    character(len=CX)      :: ice2atm_fmapname
    character(len=CX)      :: ice2atm_smapname
    character(len=CX)      :: ice2wav_smapname
    character(len=CX)      :: lnd2atm_fmapname
    character(len=CX)      :: lnd2atm_smapname
    character(len=CX)      :: lnd2glc_fmapname
    character(len=CX)      :: lnd2glc_smapname
    character(len=CX)      :: lnd2rof_fmapname
    character(len=CX)      :: ocn2atm_fmapname
    character(len=CX)      :: ocn2atm_smapname
    character(len=CX)      :: ocn2wav_smapname
    character(len=CX)      :: rof2lnd_fmapname
    character(len=CX)      :: rof2ocn_fmapname
    character(len=CX)      :: rof2ocn_ice_rmapname
    character(len=CX)      :: rof2ocn_liq_rmapname
    character(len=CX)      :: wav2ocn_smapname
    logical                :: flds_co2a  ! use case
    logical                :: flds_co2b  ! use case
    logical                :: flds_co2c  ! use case
    integer                :: glc_nec
    integer                :: mpicom
    character(len=*), parameter :: subname='(shr_nuopc_fldList_Init)'
    !--------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, localPet=localPet, mpiCommunicator=mpicom, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    mastertask = .false.
    if (localPet == 0) mastertask=.true.

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

    longname = trim(flds_scalar_name)
    stdname  = trim(flds_scalar_name)
    units    = 'unitless'
    call shr_nuopc_fldList_AddMetadata(fldname=trim(flds_scalar_name), longname=longname, stdname=stdname, units=units)
    do n = 1,ncomps
       call shr_nuopc_fldList_AddFld(fldListFr(n)%flds, trim(flds_scalar_name))
       call shr_nuopc_fldList_AddFld(fldListTo(n)%flds, trim(flds_scalar_name))
    end do

    !----------------------------------------------------------
    ! domain coordinates (appear in the share module shr_flds_mod)
    !----------------------------------------------------------

    shr_flds_dom_coord  = 'lat:lon:hgt'
    shr_flds_dom_other  = 'area:aream:mask:frac'

    longname = 'latitude'
    stdname  = 'latitude'
    units    = 'degrees north'
    call shr_nuopc_fldList_AddMetadata('lat', longname, stdname, units)

    longname = 'longitude'
    stdname  = 'longitude'
    units    = 'degrees east'
    call shr_nuopc_fldList_AddMetadata('lon', longname, stdname, units)

    longname = 'height'
    stdname  = 'height, depth, or levels'
    units    = 'unitless'
    call shr_nuopc_fldList_AddMetadata('hgt', longname, stdname, units)

    longname = 'cell_area_model'
    stdname  = 'cell area from model'
    units    = 'radian^2'
    call shr_nuopc_fldList_AddMetadata('area', longname, stdname, units)

    longname = 'cell_area_mapping'
    stdname  = 'cell area from mapping file'
    units    = 'radian^2'
    call shr_nuopc_fldList_AddMetadata('aream', longname, stdname, units)

    longname = 'mask'
    stdname  = 'mask'
    units    = '1'
    call shr_nuopc_fldList_AddMetadata('mask', longname, stdname, units)

    longname = 'area_fraction'
    stdname  = 'area fraction'
    units    = '1'
    call shr_nuopc_fldList_AddMetadata('frac', longname, stdname, units)

    !----------------------------------------------------------
    ! Masks from components
    !----------------------------------------------------------

    longname = 'Surface fraction in land'
    stdname  = 'land_fraction_from_land'
    units    = '1'
    call shr_nuopc_fldList_AddMetadata(fldname="Sl_lfrin", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds, 'Sl_lfrin')

    longname = 'Sea surface mask'
    stdname  = 'sea_surface_mask'
    units    = '1'
    call shr_nuopc_fldList_AddMetadata(fldname="So_omask", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(compocn)%flds, 'So_omask', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compocn)%flds(n1), compocn, compice,  mapfcopy, 'unset', 'unset')

    longname = 'Sea ice mask'
    stdname  = 'sea_ice_mask'
    units    = '1'
    call shr_nuopc_fldList_AddMetadata(fldname="Si_imask", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds, 'Si_imask')

    !----------------------------------------------------------
    ! Fractions sent to atm
    !----------------------------------------------------------

    longname = 'Surface land fraction'
    stdname  = 'land_area_fraction'
    units    = '1'
    call shr_nuopc_fldList_AddMetadata(fldname="Sl_lfrac", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds, 'Sl_lfrac')

    longname = 'Surface ice fraction'
    stdname  = 'sea_ice_area_fraction'
    call shr_nuopc_fldList_AddMetadata(fldname="Si_ifrac", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds, 'Si_ifrac')

    longname = 'Surface ocean fraction'
    stdname  = 'sea_area_fraction'
    units    = '1'
    call shr_nuopc_fldList_AddMetadata(fldname="So_ofrac", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds, 'So_ofrac')

    !----------------------------------------------------------
    ! Fractional ice coverage wrt ocean sent to ocn and wav
    !----------------------------------------------------------
    longname = 'Fractional ice coverage wrt ocean'
    stdname  = 'sea_ice_area_fraction'
    units    = '1'
    call shr_nuopc_fldList_AddMetadata(fldname="Si_ifrac", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds, 'Si_ifrac', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Si_ifrac', &
         merge_from1=compice, merge_field1='Si_ifrac', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compwav)%flds, 'Si_ifrac', &
         merge_from1=compice, merge_field1='Si_ifrac', merge_type1='copy')
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n1), compice, compocn,  mapfcopy, 'unset', 'unset')

    !----------------------------------------------------------
    ! Fields from atm
    !----------------------------------------------------------

    longname = 'Height at the lowest model level'
    stdname  = 'height'
    units    = 'm'
    call shr_nuopc_fldList_AddMetadata(fldname='Sa_z', longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Sa_z', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Sa_z', &
         merge_from1=compatm, merge_field1='Sa_z', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Sa_z', &
         merge_from1=compatm, merge_field1='Sa_z', merge_type1='copy')
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapbilnr, 'one', atm2lnd_smapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapbilnr, 'one', atm2ice_smapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapbilnr, 'one', atm2ocn_smapname)

    longname = 'Surface height'
    stdname  = 'height'
    units    = 'm'
    call shr_nuopc_fldList_AddMetadata(fldname='Sa_topo', longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Sa_topo', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Sa_topo', &
         merge_from1=compatm, merge_field1='Sa_topo', merge_type1='copy')
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapbilnr, 'one', atm2lnd_smapname)

    ! TODO: Sa_u and Sa_v are mapped to the ocean grid in the mediator - BUT are not sent to the ocean -
    ! They are only used in the atm/ocn flux calculation - so a special mapping will be done in the mediator
    ! for these fields

    longname = 'Zonal wind at the lowest model level'
    stdname  = 'eastward_wind'
    units    = 'm s-1'
    call shr_nuopc_fldList_AddMetadata(fldname='Sa_u', longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Sa_u', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Sa_u', &
         merge_from1=compatm, merge_field1='Sa_u', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Sa_u', &
         merge_from1=compatm, merge_field1='Sa_u', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compwav)%flds, 'Sa_u', &
         merge_from1=compatm, merge_field1='Sa_u', merge_type1='copy')
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapbilnr, 'one', atm2lnd_smapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mappatch, 'one', atm2ice_vmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mappatch, 'one', atm2ocn_vmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compwav, mapbilnr, 'one', atm2wav_smapname)

    longname = 'Meridional wind at the lowest model level'
    stdname  = 'northward_wind'
    units    = 'm s-1'
    call shr_nuopc_fldList_AddMetadata(fldname='Sa_v', longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Sa_v', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Sa_v', &
         merge_from1=compatm, merge_field1='Sa_v', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Sa_v', &
         merge_from1=compatm, merge_field1='Sa_v', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compwav)%flds, 'Sa_v', &
         merge_from1=compatm, merge_field1='Sa_v', merge_type1='copy')
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapbilnr, 'one', atm2lnd_smapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mappatch, 'one', atm2ice_vmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mappatch, 'one', atm2ocn_vmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compwav, mapbilnr, 'one', atm2wav_smapname)

    longname = 'Temperature at the lowest model level'
    stdname  = 'air_temperature'
    units    = 'K'
    call shr_nuopc_fldList_AddMetadata(fldname='Sa_tbot', longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Sa_tbot', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Sa_tbot', &
         merge_from1=compatm, merge_field1='Sa_tbot', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Sa_tbot', &
         merge_from1=compatm, merge_field1='Sa_tbot', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compwav)%flds, 'Sa_tbot', &
         merge_from1=compatm, merge_field1='Sa_tbot', merge_type1='copy')
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapbilnr, 'one', atm2lnd_smapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapbilnr, 'one', atm2ice_smapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapbilnr, 'one', atm2ocn_smapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compwav, mapbilnr, 'one', atm2wav_smapname)

    longname = 'Potential temperature at the lowest model level'
    stdname  = 'air_potential_temperature'
    units    = 'K'
    call shr_nuopc_fldList_AddMetadata(fldname='Sa_ptem', longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Sa_ptem', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Sa_ptem', &
         merge_from1=compatm, merge_field1='Sa_ptem', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Sa_ptem', &
         merge_from1=compatm, merge_field1='Sa_ptem', merge_type1='copy')
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapbilnr, 'one', atm2lnd_smapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapbilnr, 'one', atm2ice_smapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapbilnr, 'one', atm2ocn_smapname)

    longname = 'Specific humidity at the lowest model level'
    stdname  = 'specific_humidity'
    units    = 'kg kg-1'
    call shr_nuopc_fldList_AddMetadata(fldname='Sa_shum', longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Sa_shum', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Sa_shum', &
         merge_from1=compatm, merge_field1='Sa_shum', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Sa_shum', &
         merge_from1=compatm, merge_field1='Sa_shum', merge_type1='copy')
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapbilnr, 'one', atm2lnd_smapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapbilnr, 'one', atm2ice_smapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapbilnr, 'one', atm2ocn_smapname)

    longname = 'Pressure at the lowest model level'
    stdname  = 'air_pressure'
    units    = 'Pa'
    call shr_nuopc_fldList_AddMetadata(fldname='Sa_pbot', longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Sa_pbot', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Sa_pbot', &
         merge_from1=compatm, merge_field1='Sa_pbot', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Sa_pbot', &
         merge_from1=compatm, merge_field1='Sa_pbot', merge_type1='copy')
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapbilnr, 'one', atm2lnd_smapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapbilnr, 'one', atm2ice_smapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapbilnr, 'one', atm2ocn_smapname)

    longname = 'Density at the lowest model level'
    stdname  = 'air_density'
    units    = 'kg m-3'
    call shr_nuopc_fldList_AddMetadata(fldname='Sa_dens', longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Sa_dens', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Sa_dens', &
         merge_from1=compatm, merge_field1='Sa_dens', merge_type1='copy')
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice,  mapbilnr, 'one', atm2ice_smapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn,  mapbilnr, 'one', atm2ocn_smapname)

    units    = 'kg m-2 s-1'
    longname = 'Convective precipitation rate'
    stdname  = 'convective_precipitation_flux'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata(fldname='Faxa_rainc', longname=longname, stdname=stdname, units=units)
    longname = 'Large-scale (stable) precipitation rate' ! water equivalent
    stdname  = 'large_scale_precipitation_flux'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata(fldname='Faxa_rainl', longname=longname, stdname=stdname, units=units)
    longname = 'Water flux due to rain'
    stdname  = 'rainfall_flux'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata(fldname='Faxa_rain', longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Faxa_rainc', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Faxa_rainl', fldindex=n2)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Faxa_rainc', &
         merge_from1=compatm, merge_field1='Faxa_rainc', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Faxa_rainl', &
         merge_from1=compatm, merge_field1='Faxa_rainl', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Faxa_rain', &
         merge_from1=compatm, merge_field1='Faxa_rainc:Faxa_rainl', merge_type1='accumulate')
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Faxa_rain', &
         merge_from1=compatm, merge_field1='Faxa_rainc:Faxa_rainl', merge_type1='accumulate', merge_fracname1='ofrac')
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapconsf, 'one', atm2lnd_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n2), compatm, complnd, mapconsf, 'one', atm2lnd_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapconsf, 'one', atm2ice_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n2), compatm, compice, mapconsf, 'one', atm2ice_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapconsf, 'one', atm2ocn_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n2), compatm, compocn, mapconsf, 'one', atm2ocn_fmapname)

    longname = 'Convective snow rate'
    stdname  = 'convective_snowfall_flux'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata(fldname='Faxa_snowc', longname=longname, stdname=stdname, units=units)
    longname = 'Large-scale (stable) snow rate'
    stdname  = 'large_scale_snowfall_flux'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata(fldname='Faxa_snowl', longname=longname, stdname=stdname, units=units)
    longname = 'Water flux due to snow'
    stdname  = 'surface_snow_melt_flux'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata(fldname='Faxa_snow', longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Faxa_snowc', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Faxa_snowl', fldindex=n2)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Faxa_snowc', &
         merge_from1=compatm, merge_field1='Faxa_snowc', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Faxa_snowl', &
         merge_from1=compatm, merge_field1='Faxa_snowl', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Faxa_snow', &
         merge_from1=compatm, merge_field1='Faxa_snowc:Faxa_snowl', merge_type1='accumulate')
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Faxa_snow', &
         merge_from1=compatm, merge_field1='Faxa_snowc:Faxa_snowl', merge_type1='accumulate', merge_fracname1='ofrac')
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapconsf, 'one', atm2lnd_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n2), compatm, complnd, mapconsf, 'one', atm2lnd_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapconsf, 'one', atm2ice_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n2), compatm, compice, mapconsf, 'one', atm2ice_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapconsf, 'one', atm2ocn_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n2), compatm, compocn, mapconsf, 'one', atm2ocn_fmapname)

    ! total precipitation to ocean (derived rain + snow, done AFTER mappings)
    longname = 'Water flux (rain+snow)'
    stdname  = 'precipitation_flux'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata(fldname='Faxa_prec', longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Faxa_prec', &
         merge_from1=compatm, merge_field1='Faxa_rainc:Faxa_snowc:Faxa_rainl:Faxa_snowl', &
         merge_type1='accumulate', merge_fracname1='ofrac')

    longname = 'Downward longwave heat flux'
    stdname  = 'downwelling_longwave_flux'
    units    = 'W m-2'
    call shr_nuopc_fldList_AddMetadata(fldname='Faxa_lwdn', longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Faxa_lwdn', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Faxa_lwdn', &
         merge_from1=compatm, merge_field1='Faxa_lwdn', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Faxa_lwdn', &
         merge_from1=compatm, merge_field1='Faxa_lwdn', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Faxa_lwdn', &
         merge_from1=compatm, merge_field1='Faxa_lwdn', merge_type1='copy_with_weights', merge_fracname1='ofrac')
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapconsf, 'one', atm2lnd_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapconsf, 'one', atm2ice_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapconsf, 'one', atm2ocn_fmapname)

    longname = 'Direct near-infrared incident solar radiation'
    stdname  = 'surface_downward_direct_shortwave_flux_due_to_near_infrared_radiation'
    units    = 'W m-2'
    call shr_nuopc_fldList_AddMetadata(fldname="Faxa_swndr", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Faxa_swndr', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Faxa_swndr', &
         merge_from1=compatm, merge_field1='Faxa_swndr', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Faxa_swndr', &
         merge_from1=compatm, merge_field1='Faxa_swndr', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Faxa_swndr', &
         merge_from1=compatm, merge_field1='Faxa_swndr', merge_type1='copy_with_weights', merge_fracname1='ofrac')
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapconsf, 'one', atm2lnd_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapconsf, 'one', atm2ice_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapconsf, 'one', atm2ocn_fmapname)

    longname = 'Direct visible incident solar radiation'
    stdname  = 'surface_downward_direct_shortwave_flux_due_to_visible_radiation'
    units    = 'W m-2'
    call shr_nuopc_fldList_AddMetadata(fldname="Faxa_swvdr", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Faxa_swvdr', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Faxa_swvdr', &
         merge_from1=compatm, merge_field1='Faxa_swvdr', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Faxa_swvdr', &
         merge_from1=compatm, merge_field1='Faxa_swvdr', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Faxa_swvdr', &
         merge_from1=compatm, merge_field1='Faxa_swvdr', merge_type1='copy_with_weights', merge_fracname1='ofrac')
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapconsf, 'one', atm2lnd_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapconsf, 'one', atm2ice_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapconsf, 'one', atm2ocn_fmapname)

    longname = 'Diffuse near-infrared incident solar radiation'
    stdname  = 'surface_downward_diffuse_shortwave_flux_due_to_near_infrared_radiation'
    units    = 'W m-2'
    call shr_nuopc_fldList_AddMetadata(fldname="Faxa_swndf", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Faxa_swndf', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Faxa_swndf', &
         merge_from1=compatm, merge_field1='Faxa_swndf', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Faxa_swndf', &
         merge_from1=compatm, merge_field1='Faxa_swndf', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Faxa_swndf', &
         merge_from1=compatm, merge_field1='Faxa_swndf', merge_type1='copy_with_weights', merge_fracname1='ofrac')
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapconsf, 'one', atm2lnd_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapconsf, 'one', atm2ice_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapconsf, 'one', atm2ocn_fmapname)

    longname = 'Diffuse visible incident solar radiation'
    stdname  = 'surface_downward_diffuse_shortwave_flux_due_to_visible_radiation'
    units    = 'W m-2'
    call shr_nuopc_fldList_AddMetadata(fldname='Faxa_swvdf', longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Faxa_swvdf', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Faxa_swvdf', &
         merge_from1=compatm, merge_field1='Faxa_swvdf', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Faxa_swvdf', &
         merge_from1=compatm, merge_field1='Faxa_swvdf', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Faxa_swvdf', &
         merge_from1=compatm, merge_field1='Faxa_swvdf', merge_type1='copy_with_weights', merge_fracname1='ofrac')
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapconsf, 'one', atm2lnd_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapconsf, 'one', atm2ice_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapconsf, 'one', atm2ocn_fmapname)

    longname = 'Net shortwave radiation'
    stdname  = 'surface_net_shortwave_flux'
    units    = 'W m-2'
    call shr_nuopc_fldList_AddMetadata(fldname="Faxa_swnet", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddMetadata(fldname="Fall_swnet", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddMetadata(fldname="Faii_swnet", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddMetadata(fldname="Foxx_swnet", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Faxa_swnet', fldindex=n1)  ! only diagnostic
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds, 'Fall_swnet')
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds, 'Faii_swnet', fldindex=n2)  ! only diagnostic
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Foxx_swnet') ! CUSTOM
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapconsf, 'one'  , atm2ice_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapconsf, 'one'  , atm2ocn_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n2), compice, compocn, mapfcopy, 'unset', 'unset')

    longname = 'Net shortwave radiation penetrating into ice and ocean'
    stdname  = 'net_downward_shortwave_flux_in_sea_ice_due_to_penetration'
    units    = 'W m-2'
    call shr_nuopc_fldList_AddMetadata(fldname='Fioi_swpen', longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds, 'Fioi_swpen', fldindex=n1) ! used for Foxx_swnet (see below)
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n1), compice, compocn, mapfcopy, 'unset', 'unset')

    longname = 'Hydrophylic black carbon dry deposition flux'
    stdname  = 'dry_deposition_flux_of_hydrophylic_black_carbon'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata(fldname='Faxa_bcphidry', longname=longname, stdname=stdname, units=units)
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

    longname = 'Hydrophobic black carbon dry deposition flux'
    stdname  = 'dry_deposition_flux_of_hydrophobic_black_carbon'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata(fldname="Faxa_bcphodry", longname=longname, stdname=stdname, units=units)
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

    longname = 'Hydrophylic black carbon wet deposition flux'
    stdname  = 'wet_deposition_flux_of_hydrophylic_black_carbon'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata(fldname="Faxa_bcphiwet", longname=longname, stdname=stdname, units=units)
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

    longname = 'Hydrophylic organic carbon dry deposition flux'
    stdname  = 'dry_deposition_flux_of_hydrophylic_organic_carbon'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata(fldname="Faxa_ocphidry", longname=longname, stdname=stdname, units=units)
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

    longname = 'Hydrophobic organic carbon dry deposition flux'
    stdname  = 'dry_deposition_flux_of_hydrophobic_organic_carbon'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata(fldname="Faxa_ocphodry", longname=longname, stdname=stdname, units=units)
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

    longname = 'Hydrophylic organic carbon wet deposition flux'
    stdname  = 'wet_deposition_flux_of_hydrophylic_organic_carbon'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata(fldname="Faxa_ocphiwet", longname=longname, stdname=stdname, units=units)
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

    longname = 'Dust wet deposition flux (size 1)'
    stdname  = 'wet_deposition_flux_of_dust'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata(fldname="Faxa_dstwet1", longname=longname, stdname=stdname, units=units)
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

    longname = 'Dust wet deposition flux (size 2)'
    stdname  = 'wet_deposition_flux_of_dust'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata(fldname="Faxa_dstwet2", longname=longname, stdname=stdname, units=units)
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

    longname = 'Dust wet deposition flux (size 3)'
    stdname  = 'wet_deposition_flux_of_dust'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata(fldname="Faxa_dstwet3", longname=longname, stdname=stdname, units=units)
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

    longname = 'Dust wet deposition flux (size 4)'
    stdname  = 'wet_deposition_flux_of_dust'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata(fldname="Faxa_dstwet4", longname=longname, stdname=stdname, units=units)
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

    longname = 'Dust dry deposition flux (size 1)'
    stdname  = 'dry_deposition_flux_of_dust'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata(fldname="Faxa_dstdry1", longname=longname, stdname=stdname, units=units)
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

    longname = 'Dust dry deposition flux (size 2)'
    stdname  = 'dry_deposition_flux_of_dust'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata(fldname="Faxa_dstdry2", longname=longname, stdname=stdname, units=units)
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

    longname = 'Dust dry deposition flux (size 3)'
    stdname  = 'dry_deposition_flux_of_dust'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata(fldname="Faxa_dstdry3", longname=longname, stdname=stdname, units=units)
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

    longname = 'Dust dry deposition flux (size 4)'
    stdname  = 'dry_deposition_flux_of_dust'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata(fldname="Faxa_dstdry4", longname=longname, stdname=stdname, units=units)
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

    longname = 'Direct albedo (visible radiation)'
    stdname  = 'surface_direct_albedo_due_to_visible_radiation'
    units    = '1'
    call shr_nuopc_fldList_AddMetadata(fldname="Si_avsdr", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddMetadata(fldname="Sl_avsdr", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddMetadata(fldname="So_avsdr", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddMetadata(fldname="Sx_avsdr", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds , 'Sl_avsdr', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds , 'Si_avsdr', fldindex=n2)
    call shr_nuopc_fldList_AddFld(fldListMed_ocnalb_o%flds, 'So_avsdr', fldindex=n3)
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds , 'Sx_avsdr', &
         merge_from1=complnd, merge_field1='Sl_avsdr', merge_type1='merge', merge_fracname1='lfrac', &
         merge_from2=compice, merge_field2='Si_avsdr', merge_type2='merge', merge_fracname2='ifrac', &
         merge_from3=compmed, merge_field3='So_avsdr', merge_type3='merge', merge_fracname3='ofrac')
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1) , complnd, compatm, mapconsf, 'lfrin', lnd2atm_smapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n2) , compice, compatm, mapconsf, 'ifrac', ice2atm_smapname)
    call shr_nuopc_fldList_AddMap(fldListMed_ocnalb_o%flds(n3), compocn, compatm, mapconsf, 'ofrac', ocn2atm_smapname)

    longname = 'Direct albedo (near-infrared radiation)'
    stdname  = 'surface_direct_albedo_due_to_near_infrared_radiation'
    units    = '1'
    call shr_nuopc_fldList_AddMetadata(fldname="Si_anidr", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddMetadata(fldname="Sl_anidr", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddMetadata(fldname="So_anidr", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddMetadata(fldname="Sx_anidr", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds , 'Sl_anidr', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds , 'Si_anidr', fldindex=n2)
    call shr_nuopc_fldList_AddFld(fldListMed_ocnalb_o%flds, 'So_anidr', fldindex=n3)
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds , 'Sx_anidr', &
         merge_from1=complnd, merge_field1='Sl_anidr', merge_type1='merge', merge_fracname1='lfrac', &
         merge_from2=compice, merge_field2='Si_anidr', merge_type2='merge', merge_fracname2='ifrac', &
         merge_from3=compmed, merge_field3='So_anidr', merge_type3='merge', merge_fracname3='ofrac')
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1) , complnd, compatm, mapconsf, 'lfrin', lnd2atm_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n2) , compice, compatm, mapconsf, 'ifrac', ice2atm_fmapname)
    call shr_nuopc_fldList_AddMap(fldListMed_ocnalb_o%flds(n3), compocn, compatm, mapconsf, 'ofrac', ocn2atm_smapname)

    longname = 'Diffuse albedo (visible radiation)'
    stdname  = 'surface_diffuse_albedo_due_to_visible_radiation'
    units    = '1'
    call shr_nuopc_fldList_AddMetadata(fldname="Si_avsdf", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddMetadata(fldname="Sl_avsdf", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddMetadata(fldname="So_avsdf", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddMetadata(fldname="Sx_avsdf", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds , 'Sl_avsdf', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds , 'Si_avsdf', fldindex=n2)
    call shr_nuopc_fldList_AddFld(fldListMed_ocnalb_o%flds, 'So_avsdf', fldindex=n3)
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds , 'Sx_avsdf', &
         merge_from1=complnd, merge_field1='Sl_avsdf', merge_type1='merge', merge_fracname1='lfrac', &
         merge_from2=compice, merge_field2='Si_avsdf', merge_type2='merge', merge_fracname2='ifrac', &
         merge_from3=compmed, merge_field3='So_avsdf', merge_type3='merge', merge_fracname3='ofrac')
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1) , complnd, compatm, mapconsf, 'lfrin', lnd2atm_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n2) , compice, compatm, mapconsf, 'ifrac', ice2atm_fmapname)
    call shr_nuopc_fldList_AddMap(fldListMed_ocnalb_o%flds(n3), compocn, compatm, mapconsf, 'ofrac', ocn2atm_smapname)

    longname = 'Diffuse albedo (near-infrared radiation)'
    stdname  = 'surface_diffuse_albedo_due_to_near_infrared_radiation'
    units    = '1'
    call shr_nuopc_fldList_AddMetadata(fldname="Si_anidf", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddMetadata(fldname="Sl_anidf", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddMetadata(fldname="So_anidf", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddMetadata(fldname="Sx_anidf", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds , 'Sl_anidf', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds , 'Si_anidf', fldindex=n2)
    call shr_nuopc_fldList_AddFld(fldListMed_ocnalb_o%flds, 'So_anidf', fldindex=n3)
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds , 'Sx_anidf', &
         merge_from1=complnd, merge_field1='Sl_anidf', merge_type1='merge', merge_fracname1='lfrac', &
         merge_from2=compice, merge_field2='Si_anidf', merge_type2='merge', merge_fracname2='ifrac', &
         merge_from3=compmed, merge_field3='So_anidf', merge_type3='merge', merge_fracname3='ofrac')
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1) , complnd, compatm, mapconsf, 'lfrin', lnd2atm_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n2) , compice, compatm, mapconsf, 'ifrac', ice2atm_fmapname)
    call shr_nuopc_fldList_AddMap(fldListMed_ocnalb_o%flds(n3), compocn, compatm, mapconsf, 'ofrac', ocn2atm_smapname)

    longname = 'Reference temperature at 2 meters'
    stdname  = 'air_temperature'
    units    = 'K'
    call shr_nuopc_fldList_AddMetadata(fldname="Si_tref", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddMetadata(fldname="Sl_tref", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddMetadata(fldname="So_tref", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddMetadata(fldname="Sx_tref", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds , 'Sl_tref', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds , 'Si_tref', fldindex=n2)
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_a%flds, 'So_tref', fldindex=n3) ! Needed only for merging
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_o%flds, 'So_tref', fldindex=n3) ! Needed only for merging
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds , 'Sx_tref', &
         merge_from1=complnd, merge_field1='Sl_tref', merge_type1='merge', merge_fracname1='lfrac', &
         merge_from2=compice, merge_field2='Si_tref', merge_type2='merge', merge_fracname2='ifrac', &
         merge_from3=compmed, merge_field3='So_tref', merge_type3='merge', merge_fracname3='ofrac')
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1) , complnd, compatm, mapconsf, 'lfrin', lnd2atm_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n2) , compice, compatm, mapconsf, 'ifrac', ice2atm_fmapname)
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_o%flds(n3), compocn, compatm, mapconsf, 'ofrac', ocn2atm_fmapname) ! map ocn->atm
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_a%flds(n3), compatm, compocn, mapbilnr, 'one'  , atm2ocn_fmapname) ! map atm->ocn

    longname = 'Reference specific humidity at 2 meters'
    stdname  = 'specific_humidity'
    units    = 'kg kg-1'
    call shr_nuopc_fldList_AddMetadata(fldname="Si_qref", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddMetadata(fldname="Sl_qref", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddMetadata(fldname="So_qref", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddMetadata(fldname="Sx_qref", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds , 'Sl_qref', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds , 'Si_qref', fldindex=n2)
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_a%flds, 'So_qref', fldindex=n3)
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_o%flds, 'So_qref', fldindex=n3)
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds , 'Sx_qref', &
         merge_from1=complnd, merge_field1='Sl_qref', merge_type1='merge', merge_fracname1='lfrac', &
         merge_from2=compice, merge_field2='Si_qref', merge_type2='merge', merge_fracname2='ifrac', &
         merge_from3=compmed, merge_field3='So_qref', merge_type3='merge', merge_fracname3='ofrac')
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1) , complnd, compatm, mapconsf, 'lfrin', lnd2atm_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n2) , compice, compatm, mapconsf, 'ifrac', ice2atm_fmapname)
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_o%flds(n3), compocn, compatm, mapconsf, 'ofrac', ocn2atm_fmapname) ! map ocn->atm
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_a%flds(n3), compatm, compocn, mapbilnr, 'one'  , atm2ocn_fmapname) ! map atm->ocn

    longname = 'Surface temperature'
    stdname  = 'surface_temperature'
    units    = 'K'
    call shr_nuopc_fldList_AddMetadata(fldname="Si_t", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddMetadata(fldname="Sl_t", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddMetadata(fldname="So_t", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddMetadata(fldname="Sx_t", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds, 'Sl_t', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds, 'Si_t', fldindex=n2)
    call shr_nuopc_fldList_AddFld(fldListFr(compocn)%flds, 'So_t', fldindex=n3)
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
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1), complnd, compatm, mapconsf , 'lfrin', lnd2atm_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n2), compice, compatm, mapconsf , 'ifrac', ice2atm_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compocn)%flds(n3), compocn, compatm, mapconsf , 'ofrac', ocn2atm_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compocn)%flds(n3), compocn, compwav, mapbilnr , 'one'  , ocn2wav_smapname) ! This will be a custom map - need to name it however
    call shr_nuopc_fldList_AddMap(fldListFr(compocn)%flds(n3), compocn, compice, mapfcopy , 'unset', 'unset')

    longname = 'Surface fraction velocity in land'
    stdname  = 'fraction_velocity'
    units    = 'm s-1'
    call shr_nuopc_fldList_AddMetadata(fldname="Sl_fv", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds, 'Sl_fv', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds, 'Sl_fv', &
         merge_from1=complnd, merge_field1='Sl_fv', merge_type1='copy')
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1), complnd, compatm, mapconsf, 'lfrin', lnd2atm_fmapname)

    longname = 'Aerodynamic resistance'
    stdname  = 'aerodynamic_resistance'
    units    = 's/m'
    call shr_nuopc_fldList_AddMetadata(fldname="Sl_ram1", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds, 'Sl_ram1', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds, 'Sl_ram1',&
         merge_from1=complnd, merge_field1='Sl_ram1', merge_type1='copy')
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1), complnd, compatm, mapconsf, 'lfrin', lnd2atm_fmapname)

    longname = 'Surface snow water equivalent'
    stdname  = 'surface_snow_water_equivalent'
    units    = 'm'
    call shr_nuopc_fldList_AddMetadata(fldname="Sl_snowh", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds, 'Sl_snowh', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds, 'Sl_snowh', &
         merge_from1=complnd, merge_field1='Sl_snowh', merge_type1='copy')
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1), complnd, compatm, mapconsf, 'lfrin', lnd2atm_fmapname)

    longname = 'Surface snow depth'
    stdname  = 'surface_snow_thickness'
    units    = 'm'
    call shr_nuopc_fldList_AddMetadata(fldname="Si_snowh", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds, 'Si_snowh', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds, 'Si_snowh', &
         merge_from1=compice, merge_field1='Si_snowh', merge_type1='copy')
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n1), compice, compatm, mapconsf, 'ifrac', ice2atm_fmapname)

    longname = 'Surface saturation specific humidity in ocean'
    stdname  = 'specific_humidity_at_saturation'
    units    = 'kg kg-1'
    call shr_nuopc_fldList_AddMetadata(fldname="So_ssq", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_a%flds, 'So_ssq', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_o%flds, 'So_ssq')
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds , 'So_ssq', &
         merge_from1=compmed, merge_field1='So_ssq', merge_type1='copy')
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_o%flds(n1), compocn, compatm, mapconsf, 'ofrac', ocn2atm_fmapname) ! map ocn->atm
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_a%flds(n1), compatm, compocn, mapbilnr, 'one'  , atm2ocn_fmapname) ! map atm->ocn

    longname = 'Square of exch. coeff (tracers)'
    stdname  = 'square_of_exch_coeff'
    units    = '1'
    call shr_nuopc_fldList_AddMetadata(fldname="So_re", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_a%flds, 'So_re', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_o%flds, 'So_re')
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds , 'So_re', &
         merge_from1=compmed, merge_field1='So_re', merge_type1='copy')
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_o%flds(n1), compocn, compatm, mapconsf, 'ofrac', ocn2atm_fmapname) ! map ocn->atm
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_a%flds(n1), compatm, compocn, mapbilnr, 'one'  , atm2ocn_fmapname) ! map atm->ocn

    longname = '10m wind'
    stdname  = '10m_wind'
    units    = 'm'
    call shr_nuopc_fldList_AddMetadata(fldname="Sl_u10", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddMetadata(fldname="Si_u10", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddMetadata(fldname="So_u10", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddMetadata(fldname="Sx_u10", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds , 'Sl_u10', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds , 'Si_u10', fldindex=n2)
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_a%flds, 'So_u10', fldindex=n3)
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_o%flds, 'So_u10', fldindex=n3)
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds , 'Sx_u10', &
         merge_from1=complnd, merge_field1='Sl_u10', merge_type1='merge', merge_fracname1='lfrac', &
         merge_from2=compice, merge_field2='Si_u10', merge_type2='merge', merge_fracname2='ifrac', &
         merge_from3=compmed, merge_field3='So_u10', merge_type3='merge', merge_fracname3='ofrac')
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1) , complnd, compatm, mapconsf, 'lfrin', lnd2atm_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n2) , compice, compatm, mapconsf, 'ifrac', ice2atm_fmapname)
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_o%flds(n3), compocn, compatm, mapconsf, 'ofrac', ocn2atm_fmapname) ! map ocn->atm
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_a%flds(n3), compatm, compocn, mapbilnr, 'one'  , atm2ocn_fmapname) ! map atm->ocn

    longname = 'Zonal surface stress'
    stdname  = 'surface_downward_eastward_stress'
    units    = 'N m-2'
    call shr_nuopc_fldList_AddMetadata(fldname="Fall_taux", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddMetadata(fldname="Faox_taux", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddMetadata(fldname="Faii_taux", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddMetadata(fldname="Fioi_taux", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddMetadata(fldname="Faxx_taux", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddMetadata(fldname="Foxx_taux", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds , 'Fall_taux', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds , 'Faii_taux', fldindex=n2)
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds , 'Fioi_taux', fldindex=n3)
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_a%flds, 'Faox_taux', fldindex=n4)
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_o%flds, 'Faox_taux')
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds , 'Faxx_taux', &
         merge_from1=complnd, merge_field1='Fall_taux', merge_type1='merge', merge_fracname1='lfrac', &
         merge_from2=compice, merge_field2='Faii_taux', merge_type2='merge', merge_fracname2='ifrac', &
         merge_from3=compmed, merge_field3='Faox_taux', merge_type3='merge', merge_fracname3='ofrac')
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds , 'Foxx_taux', &
         merge_from1=compmed, merge_field1='Faox_taux', merge_type1='merge', merge_fracname1='ofrac', &
         merge_from2=compice, merge_field2='Fioi_taux', merge_type2='merge', merge_fracname2='ifrac')
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1) , complnd, compatm, mapconsf, 'lfrin', lnd2atm_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n2) , compice, compatm, mapconsf, 'ifrac', ice2atm_fmapname)
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_o%flds(n4), compocn, compatm, mapconsf, 'ofrac', ocn2atm_fmapname) ! map ocn->atm
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_a%flds(n4), compatm, compocn, mapconsf, 'one'  , atm2ocn_fmapname) ! map atm->ocn
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n3) , compice, compocn, mapfcopy, 'unset', 'unset')

    longname = 'Meridional surface stress'
    stdname  = 'surface_downward_northward_stress'
    units    = 'N m-2'
    call shr_nuopc_fldList_AddMetadata(fldname="Fall_tauy", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddMetadata(fldname="Faox_tauy", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddMetadata(fldname="Faii_tauy", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddMetadata(fldname="Fioi_tauy", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddMetadata(fldname="Faxx_tauy", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddMetadata(fldname="Foxx_tauy", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds , 'Fall_tauy', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds , 'Faii_tauy', fldindex=n2)
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds , 'Fioi_tauy', fldindex=n3)
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_a%flds, 'Faox_tauy', fldindex=n4)
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_o%flds, 'Faox_tauy')
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds , 'Faxx_tauy', &
         merge_from1=complnd, merge_field1='Fall_tauy', merge_type1='merge', merge_fracname1='lfrac', &
         merge_from2=compice, merge_field2='Faii_tauy', merge_type2='merge', merge_fracname2='ifrac', &
         merge_from3=compmed, merge_field3='Faox_tauy', merge_type3='merge', merge_fracname3='ofrac')
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Foxx_tauy', &
         merge_from1=compmed, merge_field1='Faox_tauy', merge_type1='merge', merge_fracname1='ofrac', &
         merge_from2=compice, merge_field2='Fioi_tauy', merge_type2='merge', merge_fracname2='ifrac')
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1) , complnd, compatm, mapconsf, 'lfrin', lnd2atm_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n2) , compice, compatm, mapconsf, 'ifrac', ice2atm_fmapname)
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_o%flds(n4), compocn, compatm, mapconsf, 'ofrac', ocn2atm_fmapname) ! map ocn->atm
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_a%flds(n4), compatm, compocn, mapconsf, 'one'  , atm2ocn_fmapname) ! map atm->ocn
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n3) , compice, compocn, mapfcopy, 'unset', 'unset')

    longname = 'Surface latent heat flux'
    stdname  = 'surface_upward_latent_heat_flux'
    units    = 'W m-2'
    call shr_nuopc_fldList_AddMetadata(fldname="Fall_lat", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddMetadata(fldname="Faox_lat", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddMetadata(fldname="Faii_lat", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddMetadata(fldname="Faxx_lat", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddMetadata(fldname="Foxx_lat", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds , 'Fall_lat', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds , 'Faii_lat', fldindex=n2)
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_a%flds, 'Faox_lat', fldindex=n3)
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_o%flds, 'Faox_lat')
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds , 'Faxx_lat', &
         merge_from1=complnd, merge_field1='Fall_lat', merge_type1='merge', merge_fracname1='lfrac', &
         merge_from2=compice, merge_field2='Faii_lat', merge_type2='merge', merge_fracname2='ifrac', &
         merge_from3=compmed, merge_field3='Faox_lat', merge_type3='merge', merge_fracname3='ofrac')
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds , 'Foxx_lat', &
         merge_from1=compmed, merge_field1='Faox_lat', merge_type1='merge', merge_fracname1='ofrac')
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1) , complnd, compatm, mapconsf, 'lfrin', lnd2atm_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n2) , compice, compatm, mapconsf, 'ifrac', ice2atm_fmapname)
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_o%flds(n3), compocn, compatm, mapconsf, 'ofrac', ocn2atm_fmapname) ! map ocn->atm
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_a%flds(n3), compatm, compocn, mapconsf, 'one'  , atm2ocn_fmapname) ! map atm->ocn

    longname = 'Sensible heat flux'
    stdname  = 'surface_upward_sensible_heat_flux'
    units    = 'W m-2'
    call shr_nuopc_fldList_AddMetadata(fldname="Fall_sen", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddMetadata(fldname="Faox_sen", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddMetadata(fldname="Faii_sen", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddMetadata(fldname="Faxx_sen", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddMetadata(fldname="Foxx_sen", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds , 'Fall_sen', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds , 'Faii_sen', fldindex=n2)
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_a%flds, 'Faox_sen', fldindex=n3)
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_o%flds, 'Faox_sen')
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds , 'Faxx_sen', &
         merge_from1=complnd, merge_field1='Fall_sen', merge_type1='merge', merge_fracname1='lfrac', &
         merge_from2=compice, merge_field2='Faii_sen', merge_type2='merge', merge_fracname2='ifrac', &
         merge_from3=compmed, merge_field3='Faox_sen', merge_type3='merge', merge_fracname3='ofrac')
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds , 'Foxx_sen', &
         merge_from1=compmed, merge_field1='Faox_sen', merge_type1='merge', merge_fracname1='ofrac')
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1) , complnd, compatm, mapconsf, 'lfrin', lnd2atm_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n2) , compice, compatm, mapconsf, 'ifrac', ice2atm_fmapname)
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_o%flds(n3), compocn, compatm, mapconsf, 'ofrac', ocn2atm_fmapname) ! map ocn->atm
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_a%flds(n3), compatm, compocn, mapconsf, 'one'  , atm2ocn_fmapname) ! map atm->ocn

    longname = 'Surface upward longwave heat flux'
    stdname  = 'surface_net_upward_longwave_flux'
    units    = 'W m-2'
    call shr_nuopc_fldList_AddMetadata(fldname="Fall_lwup", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddMetadata(fldname="Faox_lwup", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddMetadata(fldname="Faii_lwup", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddMetadata(fldname="Faxx_lwup", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddMetadata(fldname="Foxx_lwup", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds , 'Fall_lwup', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds , 'Faii_lwup', fldindex=n2)
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_a%flds, 'Faox_lwup', fldindex=n3)
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_o%flds, 'Faox_lwup')
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds , 'Faxx_lwup', &
         merge_from1=complnd, merge_field1='Fall_lwup', merge_type1='merge', merge_fracname1='lfrac', &
         merge_from2=compice, merge_field2='Faii_lwup', merge_type2='merge', merge_fracname2='ifrac', &
         merge_from3=compmed, merge_field3='Faox_lwup', merge_type3='merge', merge_fracname3='ofrac')
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds , 'Foxx_lwup', &
         merge_from1=compmed, merge_field1='Faox_lwup', merge_type1='merge', merge_fracname1='ofrac')
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1) , complnd, compatm, mapconsf, 'lfrin', lnd2atm_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n2) , compice, compatm, mapconsf, 'ifrac', ice2atm_fmapname)
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_o%flds(n3), compocn, compatm, mapconsf, 'ofrac', ocn2atm_fmapname) ! map ocn->atm
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_a%flds(n3), compatm, compocn, mapconsf, 'one'  , atm2ocn_fmapname) ! map atm->ocn

    longname = 'Evaporation water flux'
    stdname  = 'water_evaporation_flux'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata(fldname="Fall_evap", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddMetadata(fldname="Faox_evap", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddMetadata(fldname="Faii_evap", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddMetadata(fldname="Faxx_evap", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddMetadata(fldname="Foxx_evap", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds , 'Fall_evap', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds , 'Faii_evap', fldindex=n2)
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_a%flds, 'Faox_evap', fldindex=n3)
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_o%flds, 'Faox_evap')
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds , 'Faxx_evap', &
         merge_from1=complnd, merge_field1='Fall_evap', merge_type1='merge', merge_fracname1='lfrac', &
         merge_from2=compice, merge_field2='Faii_evap', merge_type2='merge', merge_fracname2='ifrac', &
         merge_from3=compmed, merge_field3='Faox_evap', merge_type3='merge', merge_fracname3='ofrac')
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds , 'Foxx_evap', &
         merge_from1=compmed, merge_field1='Faox_evap', merge_type1='merge', merge_fracname1='ofrac')
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1) , complnd, compatm, mapconsf, 'lfrin', lnd2atm_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n2) , compice, compatm, mapconsf, 'ifrac', ice2atm_fmapname)
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_o%flds(n3), compocn, compatm, mapconsf, 'ofrac', ocn2atm_fmapname) ! map ocn->atm
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_a%flds(n3), compatm, compocn, mapconsf, 'one'  , atm2ocn_fmapname) ! map atm->ocn

    longname = 'Dust flux (particle bin number 1)'
    stdname  = 'dust_flux'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata(fldname="Fall_flxdst1", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds, 'Fall_flxdst1', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds, 'Fall_flxdst1', &
         merge_from1=complnd, merge_field1='Fall_flxdst1', merge_type1='copy_with_weights', merge_fracname1='lfrac')
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1), complnd, compatm, mapconsf, 'lfrin', lnd2atm_fmapname)

    longname = 'Dust flux (particle bin number 2)'
    stdname  = 'dust_flux'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata(fldname="Fall_flxdst2", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds, 'Fall_flxdst2', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds, 'Fall_flxdst2', &
         merge_from1=complnd, merge_field1='Fall_flxdst2', merge_type1='copy_with_weights', merge_fracname1='lfrac')
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1), complnd, compatm, mapconsf, 'lfrin', lnd2atm_fmapname)

    longname = 'Dust flux (particle bin number 3)'
    stdname  = 'dust_flux'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata(fldname="Fall_flxdst3", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds, 'Fall_flxdst3', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds, 'Fall_flxdst3', &
         merge_from1=complnd, merge_field1='Fall_flxdst3', merge_type1='copy_with_weights', merge_fracname1='lfrac')
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1), complnd, compatm, mapconsf, 'lfrin', lnd2atm_fmapname)

    longname = 'Dust flux (particle bin number 4)'
    stdname  = 'dust_flux'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata(fldname="Fall_flxdst4", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds, 'Fall_flxdst4', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds, 'Fall_flxdst4', &
         merge_from1=complnd, merge_field1='Fall_flxdst4', merge_type1='copy_with_weights', merge_fracname1='lfrac')
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1), complnd, compatm, mapconsf, 'lfrin', lnd2atm_fmapname)

    !-----------------------------
    ! atm<->ocn only exchange
    !-----------------------------

    longname = 'Sea level pressure'
    stdname  = 'air_pressure_at_sea_level'
    units    = 'Pa'
    call shr_nuopc_fldList_AddMetadata(fldname="Sa_pslv", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Sa_pslv', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Sa_pslv', &
         merge_from1=compatm, merge_field1='Sa_pslv', merge_type1='copy')
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapbilnr, 'one', atm2ocn_smapname)
    call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compice, mapbilnr, 'one', atm2ocn_smapname)

    longname = 'Wind speed squared at 10 meters'
    stdname  = 'square_of_wind_speed'
    units    = 'm2 s-2'
    call shr_nuopc_fldList_AddMetadata(fldname="So_duu10n", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_a%flds, 'So_duu10n', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_o%flds, 'So_duu10n')
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds , 'So_duu10n', &
         merge_from1=compmed, merge_field1='So_duu10n', merge_type1='copy')
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_o%flds(n1), compocn, compatm, mapconsf, 'ofrac', ocn2atm_fmapname) ! map ocn->atm
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_a%flds(n1), compatm, compocn, mapbilnr, 'one'  , atm2ocn_fmapname) ! map atm->ocn

    longname = 'Surface fraction velocity in ocean'
    stdname  = 'fraction_velocity'
    units    = 'm s-1'
    call shr_nuopc_fldList_AddMetadata(fldname="So_ustar", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_a%flds, 'So_ustar', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_o%flds, 'So_ustar')
    call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds , 'So_ustar', &
         merge_from1=compmed, merge_field1='So_ustar', merge_type1='copy')
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_o%flds(n1), compocn, compatm, mapconsf, 'ofrac', ocn2atm_fmapname) ! map ocn->atm
    call shr_nuopc_fldList_AddMap(fldListMed_aoflux_a%flds(n1), compatm, compocn, mapbilnr, 'one'  , atm2ocn_fmapname) ! map atm->ocn

    !-----------------------------
    ! ice->ocn exchange
    !-----------------------------

    longname = 'Heat flux from melting'
    stdname  = 'surface_snow_melt_heat_flux'
    units    = 'W m-2'
    call shr_nuopc_fldList_AddMetadata(fldname="Fioi_melth", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds, 'Fioi_melth', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Fioi_melth', &
         merge_from1=compice, merge_field1='Fioi_melth', merge_type1='copy_with_weights', merge_fracname1='ifrac')
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n1), compice, compocn,  mapfcopy, 'unset', 'unset')

    longname = 'Water flux due to melting'
    stdname  = 'surface_snow_melt_flux'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata(fldname="Fioi_meltw", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds, 'Fioi_meltw', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Fioi_meltw', &
         merge_from1=compice, merge_field1='Fioi_meltw', merge_type1='copy_with_weights', merge_fracname1='ifrac')
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n1), compice, compocn,  mapfcopy, 'unset', 'unset')

    longname = 'Salt flux'
    stdname  = 'virtual_salt_flux_into_sea_water'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata(fldname="Fioi_salt", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds, 'Fioi_salt', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Fioi_salt', &
         merge_from1=compice, merge_field1='Fioi_salt', merge_type1='copy_with_weights', merge_fracname1='ifrac')
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n1), compice, compocn,  mapfcopy, 'unset', 'unset')

    longname = 'Hydrophylic black carbon deposition flux'
    stdname  = 'deposition_flux_of_hydrophylic_black_carbon'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata(fldname="Fioi_bcphi", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds, 'Fioi_bcphi', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Fioi_bcphi', &
         merge_from1=compice, merge_field1='Fioi_bcphi', merge_type1='copy_with_weights', merge_fracname1='ifrac')
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n1), compice, compocn,  mapfcopy, 'unset', 'unset')

    longname = 'Hydrophobic black carbon deposition flux'
    stdname  = 'deposition_flux_of_hydrophobic_black_carbon'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata(fldname="Fioi_bcpho", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds, 'Fioi_bcpho', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Fioi_bcpho', &
         merge_from1=compice, merge_field1='Fioi_bcpho', merge_type1='copy_with_weights', merge_fracname1='ifrac')
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n1), compice, compocn,  mapfcopy, 'unset', 'unset')

    longname = 'Dust flux'
    stdname  = 'dust_flux'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata(fldname="Fioi_flxdst", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds, 'Fioi_flxdst', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Fioi_flxdst', &
         merge_from1=compice, merge_field1='Fioi_flxdst', merge_type1='copy_with_weights', merge_fracname1='ifrac')
    call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n1), compice, compocn,  mapfcopy, 'unset', 'unset')

    !-----------------------------
    ! ocn -> ice exchange (some of these fields are also used in the atm/ocn flux computation)
    !-----------------------------

    longname = 'Sea surface salinity'
    stdname  = 'sea_surface_salinity'
    units    = 'g kg-1'
    call shr_nuopc_fldList_AddMetadata(fldname="So_s", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(compocn)%flds, 'So_s', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'So_s', &
         merge_from1=compocn, merge_field1='So_s', merge_type1='copy')
    call shr_nuopc_fldList_AddMap(fldListFr(compocn)%flds(n1), compocn, compice,  mapfcopy, 'unset', 'unset')

    longname = 'Zonal sea water velocity'
    stdname  = 'eastward_sea_water_velocity'
    units    = 'm s-1'
    call shr_nuopc_fldList_AddMetadata(fldname="So_u", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(compocn)%flds, 'So_u', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'So_u', &
         merge_from1=compocn, merge_field1='So_u', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compwav)%flds, 'So_u', &
         merge_from1=compocn, merge_field1='So_u', merge_type1='copy')
    call shr_nuopc_fldList_AddMap(fldListFr(compocn)%flds(n1), compocn, compice,  mapfcopy, 'unset', 'unset')
    call shr_nuopc_fldList_AddMap(fldListFr(compocn)%flds(n1), compocn, compwav,  mapbilnr, 'one'  , 'ocn2wav_smapname')

    longname = 'Fraction of sw penetrating surface layer for diurnal cycle'
    stdname  = 'Fraction_of_sw_penetrating_surface_layer'
    units    = '1'
    call shr_nuopc_fldList_AddMetadata(fldname="So_fswpen", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(compocn)%flds , 'So_fswpen', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_a%flds, 'So_fswpen')
    call shr_nuopc_fldList_AddFld(fldListMed_aoflux_o%flds, 'So_fswpen')
    call shr_nuopc_fldList_AddMap(fldListFr(compocn)%flds(n1), compocn, compice,  mapfcopy, 'unset', 'unset')

    longname = 'Meridional sea water velocity'
    stdname  = 'northward_sea_water_velocity'
    units    = 'm s-1'
    call shr_nuopc_fldList_AddMetadata(fldname="So_v", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(compocn)%flds, 'So_v', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'So_v', &
         merge_from1=compocn, merge_field1='So_v', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListTo(compwav)%flds, 'So_v', &
         merge_from1=compocn, merge_field1='So_v', merge_type1='copy')
    call shr_nuopc_fldList_AddMap(fldListFr(compocn)%flds(n1), compocn, compice,  mapfcopy, 'unset', 'unset')
    call shr_nuopc_fldList_AddMap(fldListFr(compocn)%flds(n1), compocn, compwav,  mapbilnr, 'one'  , 'ocn2wav_smapname')

    longname = 'Zonal sea surface slope'
    stdname  = 'sea_surface_eastward_slope'
    units    = 'm m-1'
    call shr_nuopc_fldList_AddMetadata(fldname="So_dhdx", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(compocn)%flds, 'So_dhdx', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'So_dhdx', &
         merge_from1=compocn, merge_field1='So_dhdx', merge_type1='copy')
    call shr_nuopc_fldList_AddMap(fldListFr(compocn)%flds(n1), compocn, compice,  mapfcopy, 'unset', 'unset')

    longname = 'Meridional sea surface slope'
    stdname  = 'sea_surface_northward_slope'
    units    = 'm m-1'
    call shr_nuopc_fldList_AddMetadata(fldname="So_dhdy", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(compocn)%flds, 'So_dhdy', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'So_dhdy', &
         merge_from1=compocn, merge_field1='So_dhdy', merge_type1='copy')
    call shr_nuopc_fldList_AddMap(fldListFr(compocn)%flds(n1), compocn, compice,  mapfcopy, 'unset', 'unset')

    longname = 'Ocean Boundary Layer Depth'
    stdname  = 'ocean_boundary_layer_depth'
    units    = 'm'
    call shr_nuopc_fldList_AddMetadata(fldname="So_bldepth", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(compocn)%flds, 'So_bldepth', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(compwav)%flds, 'So_bldepth', &
         merge_from1=compocn, merge_field1='So_bldepth', merge_type1='copy')
    call shr_nuopc_fldList_AddMap(fldListFr(compocn)%flds(n1), compocn, compwav,  mapbilnr, 'one', 'ocn2wav_smapname')

    longname = 'Ocean melt and freeze potential'
    stdname  = 'surface_snow_and_ice_melt_heat_flux'
    units    = 'W m-2'
    call shr_nuopc_fldList_AddMetadata(fldname="Fioo_q", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(compocn)%flds, 'Fioo_q', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Fioo_q', &
         merge_from1=compocn, merge_field1='Fioo_q', merge_type1='copy')
    call shr_nuopc_fldList_AddMap(fldListFr(compocn)%flds(n1), compocn, compice,  mapfcopy, 'unset', 'unset')

    !-----------------------------
    ! lnd->rof exchange
    !-----------------------------

    longname = 'Water flux from land (liquid surface)'
    stdname  = 'water_flux_into_runoff_surface'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata(fldname="Flrl_rofsur", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds, 'Flrl_rofsur', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(comprof)%flds, 'Flrl_rofsur', &
         merge_from1=complnd, merge_field1='Flrl_rofsur', merge_type1='copy_with_weights', merge_fracname1='lfrac')
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1), complnd, comprof, mapconsf, 'lfrin', lnd2rof_fmapname)

    longname = 'Water flux from land (liquid glacier, wetland, and lake)'
    stdname  = 'water_flux_into_runoff_from_gwl'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata(fldname="Flrl_rofgwl", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds, 'Flrl_rofgwl', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(comprof)%flds, 'Flrl_rofgwl', &
         merge_from1=complnd, merge_field1='Flrl_rofgwl', merge_type1='copy_with_weights', merge_fracname1='lfrac')
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1), complnd, comprof, mapconsf, 'lfrin', lnd2rof_fmapname)

    longname = 'Water flux from land (liquid subsurface)'
    stdname  = 'water_flux_into_runoff_subsurface'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata(fldname="Flrl_rofsub", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds, 'Flrl_rofsub', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(comprof)%flds, 'Flrl_rofsub', &
         merge_from1=complnd, merge_field1='Flrl_rofsub', merge_type1='copy_with_weights', merge_fracname1='lfrac')
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1), complnd, comprof, mapconsf, 'lfrin', lnd2rof_fmapname)

    longname = 'Water flux from land direct to ocean'
    stdname  = 'water_flux_direct_to_ocean'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata(fldname="Flrl_rofdto", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds, 'Flrl_rofdto', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(comprof)%flds, 'Flrl_rofdto', &
         merge_from1=complnd, merge_field1='Flrl_rofdto', merge_type1='copy_with_weights', merge_fracname1='lfrac')
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1), complnd, comprof, mapconsf, 'lfrin', lnd2rof_fmapname)

    longname = 'Water flux from land (frozen)'
    stdname  = 'frozen_water_flux_into_runoff'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata(fldname="Flrl_rofi", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds, 'Flrl_rofi', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(comprof)%flds, 'Flrl_rofi', &
         merge_from1=complnd, merge_field1='Flrl_rofi', merge_type1='copy_with_weights', merge_fracname1='lfrac')
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1), complnd, comprof, mapconsf, 'lfrin', lnd2rof_fmapname)

    ! Irrigation flux (land/rof only)
    longname = 'Irrigation flux (withdrawal from rivers)'
    stdname  = 'irrigation'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata(fldname="Flrl_irrig", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds, 'Flrl_irrig', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(comprof)%flds, 'Flrl_irrig', &
         merge_from1=complnd, merge_field1='Flrl_irrig', merge_type1='copy_with_weights', merge_fracname1='lfrac')
    call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1), complnd, comprof, mapconsf, 'lfrin', lnd2rof_fmapname)

    !-----------------------------
    ! rof->lnd
    !-----------------------------

    longname = 'Waterflux back to land due to flooding'
    stdname  = 'flooding_water_flux'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata(fldname="Flrr_flood", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(comprof)%flds, 'Flrr_flood', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Flrr_flood', &
         merge_from1=comprof, merge_field1='Flrr_flood', merge_type1='copy')
    ! TODO: who should this be handled in terms of feeding back to the ocean
    !call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Flrr_flood') ! CUSTOM
    call shr_nuopc_fldList_AddMap(fldListFr(comprof)%flds(n1), comprof, complnd, mapconsf, 'one', rof2lnd_fmapname)
    call shr_nuopc_fldList_AddMap(fldListFr(comprof)%flds(n1), comprof, compocn, mapconsf, 'one', rof2ocn_fmapname)

    longname = 'River channel total water volume'
    stdname  = 'rtm_volr'
    units    = 'm'
    call shr_nuopc_fldList_AddMetadata(fldname="Flrr_volr", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(comprof)%flds, 'Flrr_volr', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Flrr_volr', &
         merge_from1=comprof, merge_field1='Flrr_volr', merge_type1='copy')
    call shr_nuopc_fldList_AddMap(fldListFr(comprof)%flds(n1), comprof, complnd, mapconsf, 'one', rof2lnd_fmapname)

    longname = 'River channel main channel water volume'
    stdname  = 'rtm_volrmch'
    units    = 'm'
    call shr_nuopc_fldList_AddMetadata(fldname="Flrr_volrmch", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(comprof)%flds, 'Flrr_volrmch', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Flrr_volrmch', &
         merge_from1=comprof, merge_field1='Flrr_volrmch', merge_type1='copy')
    call shr_nuopc_fldList_AddMap(fldListFr(comprof)%flds(n1), comprof, complnd, mapconsf, 'one', rof2lnd_fmapname)

    !-----------------------------
    ! rof->ocn (liquid and frozen) and glc->ocn
    !-----------------------------

    longname = 'glc liquid runoff flux to ocean'
    stdname  = 'glacier_liquid_runoff_flux_to_ocean'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata(fldname='Fogg_rofl', longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(compglc)%flds, 'Fogg_rofl', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compglc)%flds(n1), compglc, compocn,  mapfiler, 'one', glc2ocn_liq_rmapname)

    longname = 'Water flux into sea water due to runoff (liquid)'
    stdname  = 'water_flux_into_sea_water'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata(fldname="Forr_rofl", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(comprof)%flds, 'Forr_rofl', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(comprof)%flds(n1), comprof, compocn, mapfiler, 'none', rof2ocn_liq_rmapname)

    longname = 'Total Water flux into sea water due to runoff (liquid)'
    stdname  = 'total_water_flux_into_sea_water'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata(fldname="Foxx_rofl", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Foxx_rofl', &
         merge_from1=comprof, merge_field1='Forr_rofl:Flrr_flood', merge_type1='accumulate', &
         merge_from2=compglc, merge_field2='Fogg_rofl'           , merge_type2='accumulate')

    longname = 'glc frozen runoff flux to ocean'
    stdname  = 'glacier_frozen_runoff_flux_to_ocean'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata(fldname='Fogg_rofi', longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(compglc)%flds, 'Fogg_rofi', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compglc)%flds(n1), compglc, compocn,  mapfiler, 'one', glc2ocn_ice_rmapname)

    longname = 'Water flux into sea water due to runoff (frozen)'
    stdname  = 'frozen_water_flux_into_sea_water'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata(fldname="Forr_rofi", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(comprof)%flds, 'Forr_rofi', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(comprof)%flds(n1), comprof, compocn, mapfiler, 'none', rof2ocn_ice_rmapname)

    longname = 'Total Water flux into sea water due to runoff (frozen)'
    stdname  = 'total_frozen_water_flux_into_sea_water'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata(fldname="Foxx_rofi", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Foxx_rofi', &
         merge_from1=comprof, merge_field1='Forr_rofi', merge_type1='accumulate', &
         merge_from2=compglc, merge_field2='Fogg_rofi', merge_type2='accumulate')

    !-----------------------------
    ! rof(frozen)->ice and glc->ice
    !-----------------------------

    longname = 'Water flux into sea ice due to runoff (frozen)'
    stdname  = 'frozen_water_flux_into_sea_ice'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata(fldname="Firr_rofi", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(comprof)%flds, 'Firr_rofi', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(comprof)%flds(n1), comprof, compice, mapfiler, 'none', rof2ocn_ice_rmapname)

    longname = 'glc frozen runoff_iceberg flux to ice'
    stdname  = 'glacier_frozen_runoff_flux_to_seaice'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata(fldname='Figg_rofi', longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(compglc)%flds, 'Figg_rofi', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListFr(compglc)%flds(n1), compglc, compice,  mapfiler, 'one', glc2ice_rmapname)

    longname = 'Total frozen water flux into sea ice '
    stdname  = 'total_frozen_water_flux_into_sea_ice'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata(fldname="Fixx_rofi", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListTo(compice)%flds, 'Fixx_rofi', &
         merge_from1=comprof, merge_field1='Firr_rofi', merge_type1='accumulate', &
         merge_from2=compglc, merge_field2='Figg_rofi', merge_type2='accumulate')

    !-----------------------------
    ! wav->ocn
    !-----------------------------

    longname = 'Langmuir multiplier'
    stdname  = 'wave_model_langmuir_multiplier'
    units    = '1'
    call shr_nuopc_fldList_AddMetadata(fldname='Sw_lamult', longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(compwav)%flds, 'Sw_lamult', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Sw_lamult', &
         merge_from1=compwav, merge_field1='Sw_lamult', merge_type1='copy')
    call shr_nuopc_fldList_AddMap(fldListFr(compwav)%flds(n1), compwav, compocn,  mapbilnr, 'one', wav2ocn_smapname)

    longname = 'Stokes drift u component'
    stdname  = 'wave_model_stokes_drift_eastward_velocity'
    units    = 'm/s'
    call shr_nuopc_fldList_AddMetadata(fldname='Sw_ustokes', longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(compwav)%flds, 'Sw_ustokes', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Sw_ustokes', &
         merge_from1=compwav, merge_field1='Sw_ustokes', merge_type1='copy')
    call shr_nuopc_fldList_AddMap(fldListFr(compwav)%flds(n1), compwav, compocn,  mapbilnr, 'one', wav2ocn_smapname)

    longname = 'Stokes drift v component'
    stdname  = 'wave_model_stokes_drift_northward_velocity'
    units    = 'm/s'
    call shr_nuopc_fldList_AddMetadata(fldname='Sw_vstokes', longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(compwav)%flds, 'Sw_vstokes', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Sw_vstokes', &
         merge_from1=compwav, merge_field1='Sw_vstokes', merge_type1='copy')
    call shr_nuopc_fldList_AddMap(fldListFr(compwav)%flds(n1), compwav, compocn, mapbilnr, 'one', wav2ocn_smapname)

    longname = 'Stokes drift depth'
    stdname  = 'wave_model_stokes_drift_depth'
    units    = 'm'
    call shr_nuopc_fldList_AddMetadata(fldname='Sw_hstokes', longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(compwav)%flds, 'Sw_hstokes', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Sw_hstokes', &
         merge_from1=compwav, merge_field1='Sw_hstokes', merge_type1='copy')
    call shr_nuopc_fldList_AddMap(fldListFr(compwav)%flds(n1), compwav, compocn, mapbilnr, 'one', wav2ocn_smapname)

    longname = 'Downward solar radiation'
    stdname  = 'surface_downward_shortwave_flux'
    units    = 'W m-2'
    call shr_nuopc_fldList_AddMetadata(fldname="Faox_swdn", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(FldListMed_ocnalb_o%flds, 'Faox_swdn', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListMed_ocnalb_o%flds(n1), compocn, compatm, mapconsf, 'ofrac', ocn2atm_smapname)

    longname = 'Upward solar radiation'
    stdname  = 'surface_upward_shortwave_flux'
    units    = 'W m-2'
    call shr_nuopc_fldList_AddMetadata(fldname="Faox_swup", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(FldListMed_ocnalb_o%flds, 'Faox_swup', fldindex=n1)
    call shr_nuopc_fldList_AddMap(fldListMed_ocnalb_o%flds(n1), compocn, compatm, mapconsf, 'ofrac', ocn2atm_smapname)

    !-----------------------------
    ! glc -> ocn
    !-----------------------------

    !-----------------------------
    ! glc -> lnd
    !-----------------------------

    ! initialize number of elevation classes
    call NUOPC_CompAttributeGet(gcomp, name='glc_nec', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) glc_nec
    call ESMF_LogWrite('glc_nec = '// trim(cvalue), ESMF_LOGMSG_INFO, rc=dbrc)

    ! for glc fields with multiple elevation classes in glc->lnd
    ! fields from glc->med do NOT have elevation classes
    ! fields from med->lnd are BROKEN into multiple elevation classes

    longname = 'Ice sheet grid coverage on global grid'
    stdname  = 'ice_sheet_grid_mask'
    units    = '1'
    call shr_nuopc_fldList_AddMetadata(fldname="Sg_icemask", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(compglc)%flds   , 'Sg_icemask', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds   , 'Sg_icemask', &
         merge_from1=compglc, merge_field1='Sg_icemask', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListMed_x2l_fr_glc%flds, 'Sg_icemask') ! Needed for FB initialization
    call shr_nuopc_fldList_AddMap(fldListFr(compglc)%flds(n1), compglc, complnd,  mapconsf, 'one', glc2lnd_smapname)

    longname = 'Ice sheet mask where we are potentially sending non-zero fluxes'
    stdname  = 'icemask_coupled'
    units    = '1'
    call shr_nuopc_fldList_AddMetadata(fldname="Sg_icemask_coupled_fluxes", longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(compglc)%flds   , 'Sg_icemask_coupled_fluxes', fldindex=n1)
    call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds   , 'Sg_icemask_coupled_fluxes', &
         merge_from1=compglc, merge_field1='Sg_icemask_coupled_fluxes', merge_type1='copy')
    call shr_nuopc_fldList_AddFld(fldListMed_x2l_fr_glc%flds, 'Sg_icemask_coupled_fluxes')
    call shr_nuopc_fldList_AddMap(fldListFr(compglc)%flds(n1), compglc, complnd,  mapconsf, 'one', glc2lnd_smapname)

    name = 'Sg_ice_covered'
    longname = 'Fraction of glacier area'
    stdname  = 'glacier_area_fraction'
    units    = '1'
    call shr_nuopc_fldList_AddMetadata(fldname='Sg_ice_covered', longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(compglc)%flds, 'Sg_ice_covered', fldindex=n1)
    call shr_nuopc_fldList_AddMap(FldListFr(compglc)%flds(n1), compglc, complnd, mapconsf, 'unset', glc2lnd_fmapname) ! TODO: normalization?
    if (glc_nec > 0) then
       do num = 0, glc_nec
          cnum = glc_elevclass_as_string(num)
          call shr_nuopc_fldList_AddMetadata(fldname= 'Sg_ice_covered'//trim(cnum), &
               longname = trim(longname)//' of elevation class '//trim(cnum), stdname =stdname, units=units)
          call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Sg_ice_covered'//trim(cnum), &
               merge_from1=compglc, merge_field1='Sg_ice_covered'//trim(cnum), merge_type1='copy')
          call shr_nuopc_fldList_AddFld(fldListMed_x2l_fr_glc%flds, 'Sg_ice_covered'//trim(cnum))
       end do
    end if

    name = 'Sg_topo'
    longname = 'Surface height of glacier'
    stdname  = 'height'
    units    = 'm'
    call shr_nuopc_fldList_AddMetadata(fldname='Sg_topo', longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(compglc)%flds, 'Sg_topo', fldindex=n1)
    call shr_nuopc_fldList_AddMap(FldListFr(compglc)%flds(n1), compglc, compglc, mapconsf, 'custom', glc2lnd_fmapname)
    if (glc_nec > 0) then
       do num = 0, glc_nec
          cnum = glc_elevclass_as_string(num)
          call shr_nuopc_fldList_AddMetadata(fldname= 'Sg_topo'//trim(cnum), &
               longname = trim(longname)//' of elevation class '//trim(cnum), stdname =stdname, units=units)
          call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Sg_topo'//trim(cnum), &
               merge_from1=compglc, merge_field1='Sg_topo'//trim(cnum), merge_type1='copy')
          call shr_nuopc_fldList_AddFld(fldListMed_x2l_fr_glc%flds, 'Sg_topo'//trim(cnum))
       end do
    end if

    name = 'Flgg_hflx'
    attname = name
    longname = 'Downward heat flux from glacier interior'
    stdname  = 'downward_heat_flux_in_glacier'
    units    = 'W m-2'
    call shr_nuopc_fldList_AddMetadata(fldname='Flgg_hflx', longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListFr(compglc)%flds, 'Flgg_hflx', fldindex=n1)
    call shr_nuopc_fldList_AddMap(FldListFr(compglc)%flds(n1), compglc, compglc, mapconsf, 'custom', glc2lnd_fmapname)
    if (glc_nec > 0) then
       do num = 0, glc_nec
          cnum = glc_elevclass_as_string(num)
          call shr_nuopc_fldList_AddMetadata(fldname= 'Flgg_hflx'//trim(cnum), &
               longname = trim(longname)//' of elevation class '//trim(cnum), stdname =stdname, units=units)
          call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds   , 'Flgg_hflx'//trim(cnum), &
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

    name = 'Flgl_qice'
    longname = 'New glacier ice flux'
    stdname  = 'ice_flux_out_of_glacier'
    units    = 'kg m-2 s-1'
    if (glc_nec > 0) then
       do num = 0, glc_nec
          cnum = glc_elevclass_as_string(num)
          call shr_nuopc_fldList_AddMetadata(fldname = 'Flgl_qice'//trim(cnum), &
               longname = trim(longname)//' of elevation class '//trim(cnum), stdname =stdname, units=units)
          call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds   , 'Flgl_qice'//trim(cnum), fldindex=n1)
          call shr_nuopc_fldList_AddFld(fldListMed_l2x_to_glc%flds, 'Flgl_qice'//trim(cnum))
       end do
    end if
    call shr_nuopc_fldList_AddMetadata(fldname= 'Flgl_qice', longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListTo(compglc)%flds, 'Flgl_qice')
    ! TODO: enter merging info
    call shr_nuopc_fldList_AddMap(FldListFr(complnd)%flds(n1), complnd, compglc, mapconsf, 'none', lnd2glc_fmapname)

    name = 'Sl_tsrf'
    longname = 'Surface temperature of glacier'
    stdname  = 'surface_temperature'
    units    = 'deg C'
    if (glc_nec > 0) then
       do num = 0, glc_nec
          cnum = glc_elevclass_as_string(num)
          call shr_nuopc_fldList_AddMetadata(fldname  = 'Sl_tsrf'//trim(cnum), &
               longname = trim(longname)//' of elevation class '//trim(cnum), stdname =stdname, units=units)
          call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds   , 'Sl_tsrf'//trim(cnum), fldindex=n1)
          call shr_nuopc_fldList_AddFld(fldListMed_l2x_to_glc%flds, 'Sl_tsrf'//trim(cnum))
       end do
    end if
    call shr_nuopc_fldList_AddMetadata(fldname= 'Sl_tsrf', longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListTo(compglc)%flds, 'Sl_tsrf')
    ! TODO: enter merging info
    call shr_nuopc_fldList_AddMap(FldListFr(complnd)%flds(n1), complnd, compglc, mapconsf, 'none', lnd2glc_fmapname)

    ! Sl_topo is sent from lnd -> med, but is NOT sent to glc (it is only used for the
    ! remapping in the mediator)
    name = 'Sl_topo'
    longname = 'Surface height'
    stdname  = 'height'
    units    = 'm'
    if (glc_nec > 0) then
       do num = 0, glc_nec
          cnum = glc_elevclass_as_string(num)
          call shr_nuopc_fldList_AddMetadata(fldname  = 'Sl_topo'//trim(cnum), &
               longname = trim(longname)//' of elevation class '//trim(cnum), stdname =stdname, units=units)
          call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds   , 'Sl_topo'//trim(cnum), fldindex=n1)
          call shr_nuopc_fldList_AddFld(fldListMed_l2x_to_glc%flds, 'Sl_topo'//trim(cnum))
       end do
    end if
    call shr_nuopc_fldList_AddMetadata(fldname= 'Sl_topo', longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListTo(compglc)%flds, 'Sl_topo')
    call shr_nuopc_fldList_AddMap(FldListFr(complnd)%flds(n1), complnd, compglc, mapconsf, 'none', lnd2glc_fmapname)

    if (flds_co2a) then

       longname = 'Prognostic CO2 at the lowest model level'
       stdname  = 'prognostic_CO2_lowest_level'
       units    = '1e-6 mol/mol'
       call shr_nuopc_fldList_AddMetadata(fldname='Sa_co2prog', longname=longname, stdname=stdname, units=units)
       call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Sa_co2prog', fldindex=n1)
       call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Sa_co2prog', &
            merge_from1=compatm, merge_field1='Sa_co2prog', merge_type1='copy')
       call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Sa_co2prog', &
            merge_from1=compatm, merge_field1='Sa_co2prog', merge_type1='copy')
       call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapbilnr, 'one', atm2lnd_smapname)
       call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapbilnr, 'one', atm2ocn_smapname)

       longname = 'Diagnostic CO2 at the lowest model level'
       stdname  = 'diagnostic_CO2_lowest_level'
       units    = '1e-6 mol/mol'
       call shr_nuopc_fldList_AddMetadata(fldname='Sa_co2diag', longname=longname, stdname=stdname, units=units)
       call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Sa_co2diag', fldindex=n1)
       call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Sa_co2diag', &
            merge_from1=compatm, merge_field1='Sa_co2diag', merge_type1='copy')
       call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Sa_co2diag', &
            merge_from1=compatm, merge_field1='Sa_co2diag', merge_type1='copy')
       call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapbilnr, 'one', atm2lnd_smapname)
       call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapbilnr, 'one', atm2ocn_smapname)

    else if (flds_co2b) then

       longname = 'Prognostic CO2 at the lowest model level'
       stdname  = 'prognostic_CO2_lowest_level'
       units    = '1e-6 mol/mol'
       call shr_nuopc_fldList_AddMetadata(fldname='Sa_co2prog', longname=longname, stdname=stdname, units=units)
       call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Sa_co2prog', fldindex=n1)
       call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Sa_co2prog', &
            merge_from1=compatm, merge_field1='Sa_co2prog', merge_type1='copy')
       call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapbilnr, 'one', atm2lnd_smapname)

       longname = 'Diagnostic CO2 at the lowest model level'
       stdname  = 'diagnostic_CO2_lowest_level'
       units    = '1e-6 mol/mol'
       call shr_nuopc_fldList_AddMetadata(fldname='Sa_co2diag', longname=longname, stdname=stdname, units=units)
       call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Sa_co2diag', fldindex=n1)
       call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Sa_co2diag', &
            merge_from1=compatm, merge_field1='Sa_co2diag', merge_type1='copy')
       call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapbilnr, 'one', atm2lnd_smapname)

       longname = 'Surface flux of CO2 from land'
       stdname  = 'surface_upward_flux_of_carbon_dioxide_where_land'
       units    = 'moles m-2 s-1'
       call shr_nuopc_fldList_AddMetadata(fldname='Fall_fco2_lnd', longname=longname, stdname=stdname, units=units)
       call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds, 'Fall_fco2_lnd', fldindex=n1)
       call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds, 'Fall_fco2_lnd', &
            merge_from1=complnd, merge_field1='Fall_fco2_lnd', merge_type1='copy_with_weights', merge_fracname1='lfrac')
       call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1), complnd, compatm, mapconsf, 'one', atm2lnd_smapname)

    else if (flds_co2c) then

       longname = 'Prognostic CO2 at the lowest model level'
       stdname  = 'prognostic_CO2_lowest_level'
       units    = '1e-6 mol/mol'
       call shr_nuopc_fldList_AddMetadata(fldname='Sa_co2prog', longname=longname, stdname=stdname, units=units)
       call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Sa_co2prog', fldindex=n1)
       call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Sa_co2prog', &
            merge_from1=compatm, merge_field1='Sa_co2prog', merge_type1='copy')
       call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Sa_co2prog', &
            merge_from1=compatm, merge_field1='Sa_co2prog', merge_type1='copy')
       call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapbilnr, 'one', atm2lnd_smapname)
       call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapbilnr, 'one', atm2ocn_smapname)

       longname = 'Diagnostic CO2 at the lowest model level'
       stdname  = 'diagnostic_CO2_lowest_level'
       units    = '1e-6 mol/mol'
       call shr_nuopc_fldList_AddMetadata(fldname='Sa_co2diag', longname=longname, stdname=stdname, units=units)
       call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, 'Sa_co2diag', fldindex=n1)
       call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, 'Sa_co2diag', &
            merge_from1=compatm, merge_field1='Sa_co2diag', merge_type1='copy')
       call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Sa_co2diag', &
            merge_from1=compatm, merge_field1='Sa_co2diag', merge_type1='copy')
       call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapbilnr, 'one', atm2lnd_smapname)
       call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapbilnr, 'one', atm2ocn_smapname)

       longname = 'Surface flux of CO2 from land'
       stdname  = 'surface_upward_flux_of_carbon_dioxide_where_land'
       units    = 'moles m-2 s-1'
       call shr_nuopc_fldList_AddMetadata(fldname='Fall_fco2_lnd', longname=longname, stdname=stdname, units=units)
       call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds, 'Fall_fco2_lnd', fldindex=n1)
       call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds, 'Fall_fco2_lnd', &
            merge_from1=complnd, merge_field1='Fall_fco2_lnd', merge_type1='copy_with_weights', merge_fracname1='lfrac')
       call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1), complnd, compatm, mapconsf, 'one', atm2lnd_smapname)

       longname = 'Surface flux of CO2 from ocean'
       stdname  = 'surface_upward_flux_of_carbon_dioxide_where_open_sea'
       units    = 'moles m-2 s-1'
       call shr_nuopc_fldList_AddMetadata(fldname='Faoo_fco2_ocn', longname=longname, stdname=stdname, units=units)
       call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds, 'Faoo_fco2_ocn', fldindex=n1)
       call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds, 'Faoo_fco2_ocn') !CUSTOM
       call shr_nuopc_fldList_AddMap(fldListFr(compocn)%flds(n1), compocn, compatm, mapconsf, 'one', ocn2atm_smapname)
    endif

    !-----------------------------
    ! water isotope fields
    !-----------------------------

    !    longname = 'Ratio of ocean surface level abund. H2_16O/H2O/Rstd'
    !    stdname  = 'ratio_ocean_surface_16O_abund'
    !    units    = '1'
    !    call fld_add(flds_o2x, flds_o2x_map, "So_roce_16O", longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_x2i, flds_x2i_map, "So_roce_16O")

    !    longname = 'Ratio of ocean surface level abund. HDO/H2O/Rstd'
    !    stdname  = 'ratio_ocean_surface_HDO_abund'
    !    call fld_add(flds_o2x, flds_o2x_map, "So_roce_HDO", longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_x2i, flds_x2i_map, "So_roce_HDO")

    !    !--------------------------------------------
    !    ! Atmospheric specific humidty at lowest level:
    !    !--------------------------------------------

    !    ! specific humidity of H216O at the lowest model level (kg/kg)
    !    longname = 'Specific humidty of H216O at the lowest model level'
    !    stdname  = 'H216OV'
    !    units    = 'kg kg-1'
    !    call fld_add(fldsfr_list(compatm) flds_a2x_map, "Sa_shum_16O", longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_x2l, flds_x2l_map, "Sa_shum_16O")
    !    call fld_add(flds_x2i, flds_x2i_map, "Sa_shum_16O")

    !    ! specific humidity of HD16O at the lowest model level (kg/kg)
    !    longname = 'Specific humidty of HD16O at the lowest model level'
    !    stdname  = 'HD16OV'
    !    call fld_add(fldsfr_list(compatm) flds_a2x_map, "Sa_shum_HDO", longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_x2l, flds_x2l_map, "Sa_shum_HDO")
    !    call fld_add(flds_x2i, flds_x2i_map, "Sa_shum_HDO")

    !    ! specific humidity of H218O at the lowest model level (kg/kg)
    !    longname = 'Specific humidty of H218O at the lowest model level'
    !    stdname  = 'H218OV'
    !    call fld_add(fldsfr_list(compatm) flds_a2x_map, "Sa_shum_18O", longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_x2l, flds_x2l_map, "Sa_shum_18O")
    !    call fld_add(flds_x2i, flds_x2i_map, "Sa_shum_18O")

    !    ! Surface snow water equivalent (land/atm only)
    !    longname = 'Isotopic surface snow water equivalent'
    !    stdname  = 'surface_snow_water_equivalent'
    !    units    = 'm'
    !    call fld_add(flds_l2x, flds_l2x_map, "Sl_snowh_16O", longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_l2x, flds_l2x_map, "Sl_snowh_18O", longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_l2x, flds_l2x_map, "Sl_snowh_HDO", longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_x2a, flds_x2a_map, "Sl_snowh_16O")
    !    call fld_add(flds_x2a, flds_x2a_map, "Sl_snowh_18O")
    !    call fld_add(flds_x2a, flds_x2a_map, "Sl_snowh_HDO")

    !    !--------------
    !    ! Isotopic Rain:
    !    !--------------

    !    !Isotopic Precipitation Fluxes:
    !    units    = 'kg m-2 s-1'
    !    longname = 'H216O Convective precipitation rate'
    !    stdname  = 'H2_16O_convective_precipitation_flux'
    !    call fld_add(fldsfr_list(compatm) flds_a2x_map, "Faxa_rainc_16O", longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_x2l, flds_x2l_map, "Faxa_rainc_16O")
    !    longname = 'H216O Large-scale (stable) precipitation rate'
    !    stdname  = 'H2_16O_large_scale_precipitation_flux'
    !    call fld_add(fldsfr_list(compatm) flds_a2x_map, "Faxa_rainl_16O", longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_x2l, flds_x2l_map, "Faxa_rainl_16O")
    !    longname = 'Water flux due to H216O rain' !equiv. to bulk
    !    stdname  = 'H2_16O_rainfall_flux'
    !    call fld_add(flds_x2o, flds_x2o_map, "Foxx_rain_16O", longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_x2i, flds_x2i_map, "Faxa_rain_16O")

    !    longname = 'H218O Convective precipitation rate'
    !    stdname  = 'H2_18O_convective_precipitation_flux'
    !    call fld_add(fldsfr_list(compatm) flds_a2x_map, "Faxa_rainc_18O", longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_x2l, flds_x2l_map, "Faxa_rainc_18O")
    !    longname = 'H218O Large-scale (stable) precipitation rate'
    !    stdname  = 'H2_18O_large_scale_precipitation_flux'
    !    call fld_add(flds_x2l, flds_x2l_map, "Faxa_rainl_18O", longname=longname, stdname=stdname, units=units)
    !    call fld_add(fldsfr_list(compatm) flds_a2x_map, "Faxa_rainl_18O")
    !    longname = 'Water flux due to H218O rain'
    !    stdname  = 'h2_18o_rainfall_flux'
    !    call fld_add(flds_x2i, flds_x2i_map, "Faxa_rain_18O", longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_x2o, flds_x2o_map, "Foxx_rain_18O", longname=longname, stdname=stdname, units=units)

    !    longname = 'HDO Convective precipitation rate'
    !    stdname  = 'HDO_convective_precipitation_flux'
    !    call fld_add(fldsfr_list(compatm) flds_a2x_map, "Faxa_rainc_HDO", longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_x2l, flds_x2l_map, "Faxa_rainc_HDO")
    !    longname = 'HDO Large-scale (stable) precipitation rate'
    !    stdname  = 'HDO_large_scale_precipitation_flux'
    !    call fld_add(flds_x2l, flds_x2l_map, "Faxa_rainl_HDO", longname=longname, stdname=stdname, units=units)
    !    call fld_add(fldsfr_list(compatm) flds_a2x_map, "Faxa_rainl_HDO")
    !    longname = 'Water flux due to HDO rain'
    !    stdname  = 'hdo_rainfall_flux'
    !    call fld_add(flds_x2i, flds_x2i_map, "Faxa_rain_HDO", longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_x2o, flds_x2o_map, "Foxx_rain_HDO", longname=longname, stdname=stdname, units=units)

    !    !-------------
    !    ! Isotopic snow:
    !    !-------------

    !    longname = 'H2_16O Convective snow rate (water equivalent)'
    !    stdname  = 'H2_16O_convective_snowfall_flux'
    !    call fld_add(fldsfr_list(compatm) flds_a2x_map, "Faxa_snowc_16O", longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_x2l, flds_x2l_map, "Faxa_snowc_16O")

    !    longname = 'H2_16O Large-scale (stable) snow rate (water equivalent)'
    !    stdname  = 'H2_16O_large_scale_snowfall_flux'
    !    call fld_add(fldsfr_list(compatm) flds_a2x_map, "Faxa_snowl_16O", longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_x2l, flds_x2l_map, "Faxa_snowl_16O")

    !    longname = 'Water equiv. H216O snow flux'
    !    stdname  = 'h2_16o_snowfall_flux'
    !    call fld_add(flds_x2o, flds_x2o_map, "Foxx_snow_16O", longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_x2i, flds_x2i_map, "Faxa_snow_16O", longname=longname, stdname=stdname, units=units)

    !    longname = 'H2_18O Convective snow rate (water equivalent)'
    !    stdname  = 'H2_18O_convective_snowfall_flux'
    !    call fld_add(fldsfr_list(compatm) flds_a2x_map, "Faxa_snowc_18O", longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_x2l, flds_x2l_map, "Faxa_snowc_18O")

    !    longname = 'H2_18O Large-scale (stable) snow rate (water equivalent)'
    !    stdname  = 'H2_18O_large_scale_snowfall_flux'
    !    call fld_add(fldsfr_list(compatm) flds_a2x_map, "Faxa_snowl_18O", longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_x2l, flds_x2l_map, "Faxa_snowl_18O")

    !    longname = 'Isotopic water equiv. snow flux of H218O'
    !    stdname  = 'h2_18o_snowfall_flux'
    !    call fld_add(flds_x2i, flds_x2i_map, "Faxa_snow_18O", longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_x2o, flds_x2o_map, "Foxx_snow_18O", longname=longname, stdname=stdname, units=units)

    !    longname = 'HDO Convective snow rate (water equivalent)'
    !    stdname  = 'HDO_convective_snowfall_flux'
    !    call fld_add(fldsfr_list(compatm) flds_a2x_map, "Faxa_snowc_HDO", longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_x2l, flds_x2l_map, "Faxa_snowc_HDO")

    !    longname = 'HDO Large-scale (stable) snow rate (water equivalent)'
    !    stdname  = 'HDO_large_scale_snowfall_flux'
    !    call fld_add(fldsfr_list(compatm) flds_a2x_map, "Faxa_snowl_HDO", longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_x2l, flds_x2l_map, "Faxa_snowl_HDO")

    !    longname = 'Isotopic water equiv. snow flux of HDO'
    !    stdname  = 'hdo_snowfall_flux'
    !    call fld_add(flds_x2i, flds_x2i_map, "Faxa_snow_HDO", longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_x2o, flds_x2o_map, "Foxx_snow_HDO", longname=longname, stdname=stdname, units=units)

    !    !----------------------------------
    !    ! Isotopic precipitation (rain+snow):
    !    !----------------------------------

    !    longname = 'Isotopic Water flux (rain+snow) for H2_16O'
    !    stdname  = 'h2_18o_precipitation_flux'
    !    call fld_add(flds_x2o, flds_x2o_map, "Foxx_prec_16O", longname=longname, stdname=stdname, units=units)  ! derived rain+snow

    !    longname = 'Isotopic Water flux (rain+snow) for H2_18O'
    !    stdname  = 'h2_18o_precipitation_flux'
    !    units    = 'kg m-2 s-1'
    !    call fld_add(flds_x2o, flds_x2o_map, "Foxx_prec_18O", longname=longname, stdname=stdname, units=units)  ! derived rain+snow

    !    longname = 'Isotopic Water flux (rain+snow) for HD_O'
    !    stdname  = 'hdo_precipitation_flux'
    !    units    = 'kg m-2 s-1'
    !    call fld_add(flds_x2o, flds_x2o_map, "Foxx_prec_HDO", longname=longname, stdname=stdname, units=units)  ! derived rain+snow

    !    !-------------------------------------
    !    ! Isotopic two meter reference humidity:
    !    !-------------------------------------

    !    ! H216O Reference specific humidity at 2 meters
    !    longname = 'Reference H216O specific humidity at 2 meters'
    !    stdname  = 'H216O_specific_humidity'
    !    units    = 'kg kg-1'
    !    call fld_add(flds_l2x, flds_l2x_map, "Sl_qref_16O", longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_i2x, flds_i2x_map, "Si_qref_16O", longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_xao, flds_xao_map, "So_qref_16O", longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_x2a, flds_x2a_map, "Sx_qref_16O", longname=longname, stdname=stdname, units=units)

    !    ! HD16O Reference specific humidity at 2 meters
    !    longname = 'Reference HD16O specific humidity at 2 meters'
    !    stdname  = 'HD16O_specific_humidity'
    !    units    = 'kg kg-1'
    !    call fld_add(flds_l2x, flds_l2x_map, "Sl_qref_HDO", longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_i2x, flds_i2x_map, "Si_qref_HDO", longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_xao, flds_xao_map, "So_qref_HDO", longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_x2a, flds_x2a_map, "Sx_qref_HDO", longname=longname, stdname=stdname, units=units)

    !    ! H218O Reference specific humidity at 2 meters
    !    longname = 'Reference H218O specific humidity at 2 meters'
    !    stdname  = 'H218O_specific_humidity'
    !    units    = 'kg kg-1'
    !    call fld_add(flds_l2x, flds_l2x_map, "Sl_qref_18O", longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_i2x, flds_i2x_map, "Si_qref_18O", longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_xao, flds_xao_map, "So_qref_18O", longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_x2a, flds_x2a_map, "Sx_qref_18O", longname=longname, stdname=stdname, units=units)

    !    !-------------------------
    !    ! Isotopic Evaporation flux:
    !    !-------------------------

    !    ! H216O Evaporation water flux
    !    longname = 'Evaporation H216O flux'
    !    stdname  = 'H216O_evaporation_flux'
    !    units    = 'kg m-2 s-1'
    !    call fld_add(flds_l2x "Fall_evap_16O", longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_i2x "Faii_evap_16O", longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_xao "Faox_evap_16O", longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_x2a "Faxx_evap_16O", longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_x2o "Foxx_evap_16O", longname=longname, stdname=stdname, units=units)

    !    ! HD16O Evaporation water flux
    !    longname = 'Evaporation HD16O flux'
    !    stdname  = 'HD16O_evaporation_flux'
    !    units    = 'kg m-2 s-1'
    !    call fld_add(flds_l2x, "Fall_evap_HDO", longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_i2x, "Faii_evap_HDO", longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_xao, "Faox_evap_HDO", longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_x2a, "Faxx_evap_HDO", longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_x2o, "Foxx_evap_HDO", longname=longname, stdname=stdname, units=units)

    !    ! H218O Evaporation water flux
    !    longname = 'Evaporation H218O flux'
    !    stdname  = 'H218O_evaporation_flux'
    !    units    = 'kg m-2 s-1'
    !    call fld_add(flds_l2x, "Fall_evap_18O", longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_i2x, "Faii_evap_18O", longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_xao, "Faox_evap_18O", longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_x2a, "Faxx_evap_18O", longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_x2o, "Foxx_evap_18O", longname=longname, stdname=stdname, units=units)

    !    !-----------------------------
    !    ! Isotopic sea ice melting flux:
    !    !-----------------------------

    !    ! H216O Water flux from melting
    !    units    = 'kg m-2 s-1'
    !    longname = 'H2_16O flux due to melting'
    !    stdname  = 'h2_16o_surface_snow_melt_flux'
    !    call fld_add(flds_i2x, flds_i2x_map, "Fioi_meltw_16O", longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_x2o, flds_x2o_map, "Foxx_meltw_16O", longname=longname, stdname=stdname, units=units)

    !    ! H218O Water flux from melting
    !    longname = 'H2_18O flux due to melting'
    !    stdname  = 'h2_18o_surface_snow_melt_flux'
    !    call fld_add(flds_i2x, flds_i2x_map, "Fioi_meltw_18O", longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_x2o, flds_x2o_map, "Foxx_meltw_18O", longname=longname, stdname=stdname, units=units)

    !    ! HDO Water flux from melting
    !    units    = 'kg m-2 s-1'
    !    longname = 'HDO flux due to melting'
    !    stdname  = 'hdo_surface_snow_melt_flux'
    !    call fld_add(flds_i2x, flds_i2x_map, "Fioi_meltw_HDO", longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_x2o, flds_x2o_map, "Foxx_meltw_HDO", longname=longname, stdname=stdname, units=units)

    !    !Iso-Runoff
    !    ! r2o, l2x, x2r

    !    units    = 'kg m-2 s-1'
    !    longname = 'H2_16O Water flux from land (frozen)'
    !    stdname  = 'H2_16O_frozen_water_flux_into_runoff'
    !    call fld_add(flds_l2x, flds_l2x_map,'Flrl_rofi_16O', longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_x2r, flds_x2r_map,'Frxx_rofi_16O', longname=longname, stdname=stdname, units=units)

    !    longname = 'H2_18O Water flux from land (frozen)'
    !    stdname  = 'H2_18O_frozen_water_flux_into_runoff'
    !    call fld_add(flds_l2x, flds_l2x_map,'Flrl_rofi_18O', longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_x2r, flds_x2r_map,'Frxx_rofi_18O', longname=longname, stdname=stdname, units=units)

    !    longname = 'HDO Water flux from land (frozen)'
    !    stdname  = 'HDO_frozen_water_flux_into_runoff'
    !    call fld_add(flds_l2x, flds_l2x_map,'Flrl_rofi_HDO', longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_x2r, flds_x2r_map,'Frxx_rofi_HDO', longname=longname, stdname=stdname, units=units)

    !    longname = 'H2_16O Water flux from land (liquid)'
    !    stdname  = 'H2_16O_liquid_water_flux_into_runoff'
    !    call fld_add(flds_l2x, flds_l2x_map,'Flrl_rofl_16O', longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_x2r, flds_x2r_map,'Frxx_rofl_16O', longname=longname, stdname=stdname, units=units)

    !    longname = 'H2_18O Water flux from land (liquid)'
    !    stdname  = 'H2_18O_liquid_water_flux_into_runoff'
    !    call fld_add(flds_l2x, flds_l2x_map,'Flrl_rofl_18O', longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_x2r, flds_x2r_map,'Frxx_rofl_18O', longname=longname, stdname=stdname, units=units)

    !    longname = 'HDO Water flux from land (liquid)'
    !    stdname  = 'HDO_liquid_water_flux_into_runoff'
    !    call fld_add(flds_l2x, flds_l2x_map,'Flrl_rofl_HDO', longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_x2r, flds_x2r_map,'Frxx_rofl_HDO', longname=longname, stdname=stdname, units=units)

    !    !-----------------------------
    !    ! Isotopic r2x, x2o
    !    !-----------------------------

    !    units    = 'kg m-2 s-1'
    !    longname = 'H2_16O Water flux due to liq runoff '
    !    stdname  = 'H2_16O_water_flux_into_sea_water'
    !    call fld_add(flds_r2x, flds_r2x_map,'Forr_rofl_16O', longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_x2o, flds_x2o_map,'Foxx_rofl_16O', longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_r2o_liq, flds_r2o_liq_map,'Forr_rofl_16O', longname=longname, stdname=stdname, units=units)

    !    longname = 'H2_18O Water flux due to liq runoff '
    !    stdname  = 'H2_18O_water_flux_into_sea_water'
    !    call fld_add(flds_r2x, flds_r2x_map,'Forr_rofl_18O', longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_x2o, flds_x2o_map,'Foxx_rofl_18O', longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_r2o_liq, flds_r2o_liq_map,'Forr_rofl_18O', longname=longname, stdname=stdname, units=units)

    !    longname = 'HDO Water flux due to liq runoff '
    !    stdname  = 'HDO_water_flux_into_sea_water'
    !    call fld_add(flds_r2x, flds_r2x_map,'Forr_rofl_HDO', longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_x2o, flds_x2o_map,'Foxx_rofl_HDO', longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_r2o_liq, flds_r2o_liq_map,'Forr_rofl_HDO', longname=longname, stdname=stdname, units=units)

    !    longname = 'H2_16O Water flux due to ice runoff '
    !    stdname  = 'H2_16O_water_flux_into_sea_water'
    !    call fld_add(flds_r2x, flds_r2x_map,'Forr_rofi_16O', longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_x2o, flds_x2o_map,'Foxx_rofi_16O', longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_r2o_ice, flds_r2o_ice_map,'Forr_rofi_16O', longname=longname, stdname=stdname, units=units)

    !    longname = 'H2_18O Water flux due to ice runoff '
    !    stdname  = 'H2_18O_water_flux_into_sea_water'
    !    call fld_add(flds_r2x, flds_r2x_map,'Forr_rofi_18O', longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_x2o, flds_x2o_map,'Foxx_rofi_18O', longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_r2o_ice, flds_r2o_ice_map,'Forr_rofi_18O', longname=longname, stdname=stdname, units=units)

    !    longname = 'HDO Water flux due to ice runoff '
    !    stdname  = 'HDO_water_flux_into_sea_water'
    !    call fld_add(flds_r2x, flds_r2x_map,'Forr_rofi_HDO', longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_x2o, flds_x2o_map,'Foxx_rofi_HDO', longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_r2o_ice, flds_r2o_ice_map,'Forr_rofi_HDO', longname=longname, stdname=stdname, units=units)

    !    ! r2x, x2l

    !    units    = 'kg m-2 s-1'
    !    longname = 'H2_16O waterrflux due to flooding'
    !    stdname  = 'H2_16O_flodding_water_flux_back_to_land'
    !    call fld_add(flds_r2x, flds_r2x_map,'Flrr_flood_16O', longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_x2l, flds_x2l_map,'Flrr_flood_16O')

    !    longname = 'H2_18O waterrflux due to flooding'
    !    stdname  = 'H2_18O_flodding_water_flux_back_to_land'
    !    call fld_add(flds_r2x, flds_r2x_map,'Flrr_flood_18O', longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_x2l, flds_x2l_map,'Flrr_flood_18O')

    !    longname = 'HDO Waterrflux due to flooding'
    !    stdname  = 'HDO_flodding_water_flux_back_to_land'
    !    call fld_add(flds_r2x, flds_r2x_map,'Flrr_flood_HDO', longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_x2l, flds_x2l_map,'Flrr_flood_HDO')

    !    longname = 'H2_16O river channel water volume '
    !    stdname  = 'H2_16O_rtm_volr'
    !    call fld_add(flds_r2x, flds_r2x_map,'Flrr_volr_16O', longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_x2l, flds_x2l_map,'Flrr_volr_16O')

    !    longname = 'H2_18O river channel water volume '
    !    stdname  = 'H2_18O_rtm_volr'
    !    call fld_add(flds_r2x, flds_r2x_map,'Flrr_volr_18O', longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_x2l, flds_x2l_map,'Flrr_volr_18O')

    !    longname = 'HDO river channel water volume '
    !    stdname  = 'HDO_rtm_volr'
    !    call fld_add(flds_r2x, flds_r2x_map,'Flrr_volr_HDO', longname=longname, stdname=stdname, units=units)
    !    call fld_add(flds_x2l, flds_x2l_map,'Flrr_volr_HDO')

    !    ! longname = 'H2_18O Waterrflux due to flooding'
    !    ! stdname  = 'H2_18O_flodding_water_flux_back_to_land'
    !    ! call fld_add(flds_r2x, flds_r2x_map,'Flrr_flood_HDO', longname=longname, stdname=stdname, units=units)
    !    ! call fld_add(flds_x2l, flds_x2l_map,'Flrr_flood_HDO')

    !-----------------------------------------------------------------------------
    ! optional per thickness category fields
    !-----------------------------------------------------------------------------

    call NUOPC_CompAttributeGet(gcomp, name='ice_ncat', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) ice_ncat
    call ESMF_LogWrite('ice_ncat = '// trim(cvalue), ESMF_LOGMSG_INFO, rc=dbrc)

    do num = 1, ice_ncat
       write(cnum,'(i2.2)') num

       ! Fractional ice coverage wrt ocean
       name = 'Si_ifrac_' // cnum
       longname = 'fractional ice coverage wrt ocean for thickness category ' // cnum
       stdname  = 'sea_ice_area_fraction'
       units    = '1'
       call shr_nuopc_fldList_AddMetadata(fldname=trim(name), longname=longname, stdname=stdname, units=units)
       call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds, trim(name), fldindex=n1)
       call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, trim(name), &
            merge_from1=compice, merge_field1=trim(name), merge_type1='copy')
       call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n1), compice, compocn,  mapfcopy, 'unset', 'unset')

       ! Net shortwave radiation
       name = 'PFioi_swpen_ifrac_' // cnum
       longname = 'net shortwave radiation penetrating into ice and ocean times ice fraction for thickness category ' // cnum
       stdname  = 'product_of_net_downward_shortwave_flux_at_sea_water_surface_and_sea_ice_area_fraction'
       units    = 'W m-2'
       call shr_nuopc_fldList_AddMetadata(fldname=trim(name), longname=longname, stdname=stdname, units=units)
       call shr_nuopc_fldList_AddFld(fldListFr(compice)%flds, trim(name), fldindex=n1)
       call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, trim(name), &
            merge_from1=compice, merge_field1=trim(name), merge_type1='copy')
       call shr_nuopc_fldList_AddMap(fldListFr(compice)%flds(n1), compice, compocn,  mapfcopy, 'unset', 'unset')
    end do

    longname = 'fractional atmosphere coverage wrt ocean'
    stdname  = 'atmosphere_area_fraction'
    units    = '1'
    call shr_nuopc_fldList_AddMetadata(fldname='Sf_afrac', longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Sf_afrac')
    ! TODO: add mapping and merging

    longname = 'fractional atmosphere coverage used in radiation computations wrt ocean'
    stdname  = 'atmosphere_area_fraction'
    units    = '1'
    call shr_nuopc_fldList_AddMetadata(fldname='Sf_afracr', longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Sf_afracr')
    ! TODO: add mapping and merging

    name = 'Foxx_swnet_afracr'
    longname = 'net shortwave radiation times atmosphere fraction'
    stdname = 'product_of_net_downward_shortwave_flux_at_sea_water_surface_and_atmosphere_area_fraction'
    units = 'W m-2'
    call shr_nuopc_fldList_AddMetadata(fldname='Foxx_swnet_afracr', longname=longname, stdname=stdname, units=units)
    call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, 'Foxx_swnet_afracr')
    ! TODO: add mapping and merging

    !-----------------------------------------------------------------------------
    ! CARMA fields
    ! if carma_flds are specified then setup fields for CLM to CAM communication
    !-----------------------------------------------------------------------------

    call shr_carma_readnl('drv_flds_in', mpicom, mastertask, carma_fields)
    if (carma_fields /= ' ') then
       longname = 'Volumetric soil water'
       stdname  = 'soil_water'
       units    = 'm3/m3'
       do n = 1,shr_string_listGetNum(carma_fields)
          call shr_string_listGetName(carma_fields, n, fldname)
          call shr_nuopc_fldList_AddMetadata(fldname=trim(fldname), longname=longname, stdname=stdname, units=units)
          call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds, trim(fldname), fldindex=n1)
          call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds, trim(fldname), &
               merge_from1=complnd, merge_field1=trim(fldname), merge_type1='copy')
          call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1), complnd, compatm, mapconsf, 'one', lnd2atm_smapname)
       enddo
    endif

    !-----------------------------------------------------------------------------
    ! MEGAN emissions fluxes fields
    ! if MEGAN emission are specified then setup fields for CLM to CAM communication
    !-----------------------------------------------------------------------------

    ! Note that shr_megan_readnl returns megan_voc_fields which is a
    ! colon deliminated string of the megan foc fields that will be
    ! exported by the land model

    call shr_megan_readnl('drv_flds_in', mpicom, mastertask, megan_voc_fields)
    if (shr_megan_mechcomps_n > 0) then
       longname = 'MEGAN emission fluxes'
       stdname  = 'megan'
       units    = 'molecules/m2/sec'
       do n = 1,shr_string_listGetNum(megan_voc_fields)
          call shr_string_listGetName(megan_voc_fields, n, fldname)
          call shr_nuopc_fldList_AddMetadata(fldname=trim(fldname), longname=longname, stdname=stdname, units=units)
          call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds, trim(fldname), fldindex=n1)
          call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds, trim(fldname), &
               merge_from1=complnd, merge_field1=trim(fldname), merge_type1='merge', merge_fracname1='lfrac')
          call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1), complnd, compatm, mapconsf, 'one', atm2lnd_smapname)
       enddo
    endif

    !-----------------------------------------------------------------------------
    ! Fire emissions fluxes fields
    ! if fire emission are specified then setup fields for CLM to CAM communication
    ! (emissions fluxes)
    !-----------------------------------------------------------------------------

    call shr_fire_emis_readnl('drv_flds_in', mpicom, mastertask, fire_emis_fields)
    if (shr_fire_emis_mechcomps_n>0) then
       longname = 'wild fire emission fluxes'
       stdname  = 'fire_emis'
       units    = 'kg/m2/sec'
       do n = 1,shr_string_listGetNum(fire_emis_fields)
          call shr_string_listGetName(fire_emis_fields, n, fldname)
          call shr_nuopc_fldList_AddMetadata(fldname=trim(fldname), longname=longname, stdname=stdname, units=units)
          call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds, trim(fldname), fldindex=n1)
          call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds, trim(fldname), &
               merge_from1=complnd, merge_field1=trim(fldname), merge_type1='merge', merge_fracname1='lfrac')
          call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1), complnd, compatm, mapconsf, 'one', lnd2atm_smapname)
       enddo

       longname = 'wild fire plume height'
       stdname  = 'fire_plume_top'
       units    = 'm'
       call shr_nuopc_fldList_AddMetadata(fldname=trim(shr_fire_emis_ztop_token), longname=longname, stdname=stdname, units=units)
       call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds, trim(shr_fire_emis_ztop_token), fldindex=n1)
       call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds, trim(shr_fire_emis_ztop_token), &
               merge_from1=complnd, merge_field1=trim(fldname), merge_type1='copy')
       call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1), complnd, compatm, mapconsf, 'one', lnd2atm_smapname)
    endif

    !-----------------------------------------------------------------------------
    ! Dry Deposition fields
    ! First read namelist and figure out the drydep field list to pass
    ! Then check if file exists and if not, n_drydep will be zero
    ! Then add dry deposition fields to land export and atmosphere import states
    ! Then initialize dry deposition fields
    ! Note: CAM and CLM will then call seq_drydep_setHCoeff
    !-----------------------------------------------------------------------------

    call seq_drydep_readnl("drv_flds_in", mpicom, mastertask, drydep_fields)
    if ( lnd_drydep ) then
       longname = 'dry deposition velocity'
       stdname  = 'drydep_vel'
       units    = 'cm/sec'
       do n = 1,shr_string_listGetNum(drydep_fields)
          call shr_string_listGetName(drydep_fields, n, fldname)
          call shr_nuopc_fldList_AddMetadata(fldname=fldname, longname=longname, stdname=stdname, units=units)
          call shr_nuopc_fldList_AddFld(fldListFr(complnd)%flds, trim(fldname), fldindex=n1)
          call shr_nuopc_fldList_AddFld(fldListTo(compatm)%flds, trim(fldname), &
               merge_from1=complnd, merge_field1=trim(fldname), merge_type1='copy')
          call shr_nuopc_fldList_AddMap(fldListFr(complnd)%flds(n1), complnd, compatm, mapconsf, 'one', lnd2atm_smapname)
       enddo
    endif
    call seq_drydep_init( )

    !-----------------------------------------------------------------------------
    ! Nitrogen Deposition fields
    ! First read namelist and figure out the ndepdep field list to pass
    ! Then check if file exists and if not, n_drydep will be zero
    ! Then add nitrogen deposition fields to atm export, lnd import and ocn import
    !-----------------------------------------------------------------------------

    call shr_ndep_readnl("drv_flds_in", mpicom, mastertask, ndep_fields, add_ndep_fields)
    if (add_ndep_fields) then
       longname = 'nitrogen deposition flux'
       stdname  = 'nitrogen_deposition'
       units    = 'kg(N)/m2/sec'
       do n = 1,shr_string_listGetNum(ndep_fields)
          call shr_string_listGetName(ndep_fields, n, fldname)
          call shr_nuopc_fldList_AddMetadata(fldname=trim(fldname), longname=longname, stdname=stdname, units=units)
          call shr_nuopc_fldList_AddFld(fldListFr(compatm)%flds, trim(fldname), fldindex=n1)
          call shr_nuopc_fldList_AddFld(fldListTo(complnd)%flds, trim(fldname), &
               merge_from1=compatm, merge_field1=trim(fldname), merge_type1='copy_with_weights', merge_fracname1='ofrac')
          call shr_nuopc_fldList_AddFld(fldListTo(compocn)%flds, trim(fldname), &
               merge_from1=compatm, merge_field1=trim(fldname), merge_type1='copy_with_weights', merge_fracname1='ofrac')
          call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, complnd, mapbilnr, 'one', atm2lnd_smapname)
          call shr_nuopc_fldList_AddMap(fldListFr(compatm)%flds(n1), compatm, compocn, mapbilnr, 'one', atm2ocn_smapname)
       enddo
    end if

  end subroutine esmFlds_Init

end module esmFlds
