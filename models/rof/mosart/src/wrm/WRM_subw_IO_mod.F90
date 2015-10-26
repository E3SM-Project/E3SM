!
MODULE WRM_subw_IO_mod
! Description: module to provide interface between WRM and other CLM components
! 
! Developed by Nathalie Voisin 2/1/2010
! REVISION HISTORY:
!-----------------------------------------------------------------------

! !USES:
  use RunoffMod     , only : Tctl, TUnit, rtmCTL
  use RtmSpmd       , only : masterproc
  use RtmVar        , only : iulog
  use RtmIO         , only : pio_subsystem, ncd_pio_openfile, ncd_pio_closefile
  use rof_cpl_indices, only : nt_rtm
  use shr_kind_mod  , only : r8 => shr_kind_r8, SHR_KIND_CL
  use shr_const_mod , only : SHR_CONST_REARTH, SHR_CONST_PI
  use shr_sys_mod   , only : shr_sys_flush, shr_sys_abort
  use WRM_type_mod  , only : ctlSubwWRM, WRMUnit, StorWater
  use WRM_start_op_year, only : WRM_init_StOp_FC
  use netcdf
  use pio
  
  implicit none
  private

  public WRM_init

!-----------------------------------------------------------------------
  contains
!-----------------------------------------------------------------------
  
  subroutine WRM_init
     ! !DESCRIPTION: initilization of WRM model
     implicit none

     integer :: nr, nd, iunit, j , mth, nsub, i, n      ! local loop indices
     integer :: cnt, nn, idam, ntotal
     integer :: begr, endr, maxnumdependentgrid
     integer :: ncdfid, did, dids(2), dsizes(2), ier, varid
     type(file_desc_t):: ncid       ! netcdf file
     type(var_desc_t) :: vardesc    ! netCDF variable description
     type(io_desc_t)  :: iodesc_int ! pio io desc
     type(io_desc_t)  :: iodesc_dbl ! pio io desc
     type(io_desc_t)  :: iodesc_int_dam ! pio io desc
     type(io_desc_t)  :: iodesc_dbl_dam ! pio io desc
     integer(kind=pio_offset) :: frame
     real(r8) :: peak, prorata , mn, mx                  ! peak value to define the start of operationalyr
     integer :: sgn,curr_sgn, nsc, ct, ct_mx, mth_op , idepend      ! number of sign change
     integer, pointer  :: compDOF(:)     ! pio decomp
     integer, pointer  :: idam_global(:) ! local to global dams mapping
     integer, pointer  :: idam_glo2dam(:) ! global to local dams mapping
     integer, pointer  :: UnitID_1D(:)   ! 1d to 2d dam mapping
!     integer, pointer  :: temp1D_int(:), temp2D_nc_int(:,:), temp2D_local_int(:,:), temp2D_int(:,:), temp3D_int(:,:,:) ! temporary
!     real(r8), pointer :: temp1D_dbl(:), temp2D_nc_dbl(:,:), temp2D_local_dbl(:,:), temp2D_dbl(:,:), temp3D_dbl(:,:,:) ! temporary
     integer, pointer  :: temp_gridID_from_Dam(:,:)
     character(len=*),parameter :: subname='WRM_init'
     character(len=*),parameter :: FORMI = '(2A,2i10)'
     character(len=*),parameter :: FORMR = '(2A,2g15.7)'

     begr = rtmCTL%begr
     endr = rtmCTL%endr

     if (masterproc) write(iulog,*) subname,' begr endr = ',begr,endr
     call shr_sys_flush(iulog)

     ctlSubwWRM%paraPath = '/pic/scratch/tcraig/wm_data/'
     ctlSubwWRM%paraFile = 'US_reservoir_8th_NLDAS2.nc'
     ctlSubwWRM%demandPath = '/pic/scratch/tcraig/wm_data/NLDAS2_GCAM_water_demand_'

     !-------------------
     ! dam count
     !-------------------

     ! start by opening and reading the dam IDs, can't use pio yet, need decomp first
     ier = nf90_open(trim(ctlSubwWRM%paraPath)//trim(ctlSubwWRM%paraFile), NF90_NOWRITE, ncdfid)
     ier = nf90_inq_dimid(ncdfid,'DependentGrids',did)
     ier = nf90_inquire_dimension(ncdfid,did,len=maxNumDependentGrid)
     ier = nf90_inq_dimid(ncdfid,'Dams',did)
     ier = nf90_inquire_dimension(ncdfid,did,len=CtlSubwWRM%NDam)

     allocate(UnitID_1D(ctlSubwWRM%NDam))
     ier = nf90_inq_varid(ncdfid,'unitID_1D',varid)
     ier = nf90_get_var(ncdfid,varid,UnitID_1d)
     allocate(temp_gridID_from_Dam(ctlSubwWRM%NDam, maxNumDependentGrid))
     ier = nf90_inq_varid(ncdfid,'gridID_from_Dam',varid)
     ier = nf90_get_var(ncdfid,varid,temp_gridID_from_Dam)
     ier = nf90_close(ncdfid)

     if (masterproc) write(iulog,*) subname,' dams = ',CtlSubwWRM%NDam
     if (masterproc) write(iulog,*) subname,' dependentgrids = ',maxNumDependentGrid
     if (masterproc) write(iulog,*) subname,' gridID_from_Dam = ',minval(temp_gridID_from_Dam),maxval(temp_gridID_from_Dam)

     ! now open same file with pio
     call ncd_pio_openfile(ncid, trim(ctlSubwWRM%paraPath)//trim(ctlSubwWRM%paraFile), 0)
     allocate(compdof(rtmCTL%lnumr))
     cnt = 0
     do n = begr,endr
        cnt = cnt + 1
!tcx, tcraig which of these is correct/better?
        compDOF(cnt) = rtmCTL%gindex(n)
!        compDOF(cnt) = TUnit%ID0(n)
     enddo
     ier = pio_inq_varid   (ncid, name='unit_ID', vardesc=vardesc)
     ier = pio_inq_vardimid(ncid, vardesc, dids)
     ier = pio_inq_dimlen  (ncid, dids(1),dsizes(1))
     ier = pio_inq_dimlen  (ncid, dids(2),dsizes(2))

     if (masterproc) write(iulog,*) subname,' dsizes = ',dsizes
     write(iulog,*) subname,' compDOF = ',minval(compDof),maxval(compDOF)
     call shr_sys_flush(iulog)

! decomp for lon/lat to local gridcells
     call pio_initdecomp(pio_subsystem, pio_double, dsizes, compDOF, iodesc_dbl)
     call pio_initdecomp(pio_subsystem, pio_int   , dsizes, compDOF, iodesc_int)
     deallocate(compdof)

! handled above in standard netcdf
!     ier = pio_inq_varid   (ncid, name='Dams', vardesc=vardesc)
!     ier = pio_inq_vardimid(ncid, vardesc, dids)
!     ier = pio_inq_dimlen  (ncid, dids(1), CtlSubwWRM%NDam)
!     ier = pio_inq_varid   (ncid, name='DependentGrids', vardesc=vardesc)
!     ier = pio_inq_vardimid(ncid, vardesc, dids)
!     ier = pio_inq_dimlen  (ncid, dids(1), maxNumDependentGrid)

     allocate (WRMUnit%isDam(begr:endr))
     WRMUnit%isDam = 0
     ier = pio_inq_varid (ncid, name='unit_ID', vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_int, WRMUnit%isDam, ier)
     if (masterproc) write(iulog,FORMI) trim(subname),' read unit_ID ',minval(WRMunit%isDam),maxval(WRMUnit%isDam)
     call shr_sys_flush(iulog)

     !--- count dams and set up dam indexing

     ctlSubwWRM%localNumDam = 0
     do nr=begr,endr
        if (WRMUnit%isDam(nr) > 0) then
           if (WRMUnit%isDam(nr) /= Tunit%ID0(nr)) then
              write(iulog,*) subname,' ERROR: Dam unit_ID does not match mosart ',nr,WRMUnit%isDam(nr),Tunit%ID0(nr)
              call shr_sys_abort(subname//' ERROR: Dam unit_ID does not match mosart')
           endif
           ctlSubwWRM%localNumDam = ctlSubwWRM%localNumDam + 1
        else
           WRMUnit%isDam(nr) = 0
        end if
     end do
     write(iulog,*) subname,' localNumDam = ',ctlSubwWRM%localNumDam

     !initialize INCVicell 
! tcx tcraig check this loop for correctness
     allocate(WRMUnit%INVicell(begr:endr))
     WRMUnit%INVicell =-99 
     allocate(WRMUnit%icell(ctlsubwWRM%localNumDam))
     WRMUnit%icell = 0
     allocate(idam_global(ctlSubwWRM%localNumDam))
     allocate(idam_glo2dam(ctlSubwWRM%localNumDam))
     idam_global = -99
     idam_glo2dam = -99
     nn = 0
     do nr=begr,endr
        if (WRMUnit%isDam(nr) > 0) then
           nn = nn + 1
           WRMUnit%INVicell(nr) = nn
           WRMUnit%icell(nn) = nr
           do nd=1,ctlSubwWRM%NDam
              if (WRMUnit%isDam(nr) .eq. UnitID_1D(nd)) then
!                 write(iulog,'(a,5i10)') subname//' idam_glo=',nr,nn,WRMUnit%isDam(nr),nd,UnitID_1D(nd)
                 idam_global(nn) = nd
                 idam_glo2dam(nn) = WRMUnit%isDam(nr)
                 exit
              end if
           end do
        end if
     end do

     write(iulog,*) subname,' icell    = ',minval(WRMUnit%icell),maxval(WRMUnit%icell)
     write(iulog,*) subname,' INVicell = ',minval(WRMUnit%INVicell),maxval(WRMUnit%INVicell)
     write(iulog,*) subname,' idam_global  = ',minval(idam_global),maxval(idam_global)
     write(iulog,*) subname,' idam_glo2dam = ',minval(idam_glo2dam),maxval(idam_glo2dam)

! decomp for ndam to local ndam
!     call pio_initdecomp(pio_subsystem, pio_double, ctlSubwWRM%NDam, idam_global, iodesc_dbl_dam)
!     call pio_initdecomp(pio_subsystem, pio_int   , ctlSubwWRM%NDam, idam_global, iodesc_int_dam)
! decomp for lon/lat to local ndam
     call pio_initdecomp(pio_subsystem, pio_double, dsizes, idam_glo2dam, iodesc_dbl_dam)
     call pio_initdecomp(pio_subsystem, pio_int   , dsizes, idam_glo2dam, iodesc_int_dam)

     !initialize dam dependencies
     allocate (WRMUnit%dam_Ndepend(ctlSubwWRM%localNumDam))
     WRMUnit%dam_Ndepend = 0        !
     allocate (WRMUnit%dam_depend(ctlSubwWRM%localNumDam,maxNumDependentGrid))
     WRMUnit%dam_depend = 0        !

     ier = pio_inq_varid (ncid, name='numGrid_from_Dam', vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_int_dam , WRMUnit%dam_Ndepend, ier)
     write(iulog,FORMI) trim(subname),' read numGrid_from_Dam',minval(WRMUnit%dam_Ndepend),maxval(WRMUnit%dam_Ndepend)
     call shr_sys_flush(iulog)

! initialize dam_depend
     do idam = 1,ctlSubwWRM%localNumDam
        nn = idam_global(idam)
        ntotal = 0
        do nd=1,maxNumDependentGrid
           if (temp_gridID_from_Dam(nn,nd) > 0) then
              ntotal = ntotal + 1
              WRMUnit%dam_depend(idam,ntotal) = temp_gridID_from_Dam(nn,nd)
           endif
        enddo
        write(iulog,*) subname,' compute dam_depend ntotal ',idam,ntotal,WRMUnit%dam_Ndepend(idam)
        if (ntotal > WRMUnit%dam_Ndepend(idam) .or. ntotal < WRMUnit%dam_Ndepend(idam)-1) then ! sometimes in the dependency data, the basin outlet is included for a dam as its dependent grid
           write(iulog,"(a,3i8)"), subname//" Attention check gridID_from_Dam ",idam, WRMUnit%dam_Ndepend(idam), ntotal
           !call endrun
        end if
        WRMUnit%dam_Ndepend(idam) = ntotal
     enddo
     deallocate(temp_gridID_from_Dam)

!tcx tcraig need to update this!
     !initialize INVsubw
! tcx deprecated by hongyi
!     allocate (WRMUnit%INVisubw(nsub))
!     WRMUnit%INVisubw = -99
!     WRMUnit%NUnitID = 0
!     do iunit=1,ctlSubWRM%NUnit
!        j = TUnit%icell(iunit, 1) ! cell number where dam is located, need indice
!        if ( j > WRMUnit%NUnitID ) then
!           WRMUnit%NUnitID = j
!        endif
!        WRMUnit%INVisubw(j)=iunit
!     end do
!     if ( WRMUnit%NUnitID > nsub ) then
!        write(iulog,*) subname,"ATTENTION, ID max number is larger than the hard coded value of %d, adjust in full knowledge, nsub "
!        call shr_sys_abort(subname//' ERROR: ID max number too large')
!     endif

     ! max subw id number, tcx
!     nsub = maxval(rglo2gdc)??

     allocate (WRMUnit%subw_Ndepend(begr:endr))
     WRMUnit%subw_Ndepend = 0        ! 
!tcx get rid of subw_depend?
!     allocate (WRMUnit%subw_depend(begr:endr,ctlSubwWRM%localNumDam))
!     WRMUnit%subw_depend = 0        !

     allocate (WRMUnit%DamName(ctlSubwWRM%localNumDam))
     allocate (WRMUnit%TotStorCapDepend(begr:endr))
     WRMUnit%TotStorCapDepend = 0._r8
     allocate (WRMUnit%TotInflowDepend(begr:endr))
     WRMUnit%TotInflowDepend = 0._r8

     allocate (WRMUnit%mask(ctlSubwWRM%localNumDam))
     WRMUnit%mask = 0
     allocate (WRMUnit%YEAR(ctlSubwWRM%localNumDam))
     WRMUnit%YEAR = 1900

     allocate (WRMUnit%SurfArea(ctlSubwWRM%localNumDam))
     WRMUnit%SurfArea = 0._r8
     allocate (WRMUnit%InstCap(ctlSubwWRM%localNumDam))
     WRMUnit%InstCap = 0._r8
     allocate (WRMUnit%StorCap(ctlSubwWRM%localNumDam))
     WRMUnit%StorCap = 0._r8
     allocate (WRMUnit%Height(ctlSubwWRM%localNumDam))
     WRMUnit%Height = 0._r8
     allocate (WRMUnit%Length(ctlSubwWRM%localNumDam))
     WRMUnit%Length = 0._r8
     allocate (WRMUnit%Depth(ctlSubwWRM%localNumDam))
     WRMUnit%Depth = 0._r8

     allocate (WRMUnit%MeanMthFlow(ctlSubwWRM%localNumDam,13))
     WRMUnit%MeanMthFlow = 0._r8
     allocate (WRMUnit%MeanMthDemand(ctlSubwWRM%localNumDam,13))
     WRMUnit%MeanMthDemand = 0._r8
     allocate (WRMUnit%StorTarget(ctlSubwWRM%localNumDam,13))
     WRMUnit%StorTarget = 0._r8
     allocate (WRMUnit%MinStorTarget(ctlSubwWRM%localNumDam))
     WRMUnit%MinStorTarget = 0._r8
     allocate (WRMUnit%MaxStorTarget(ctlSubwWRM%localNumDam))
     WRMUnit%MaxStorTarget = 0._r8
     allocate (WRMUnit%StorageCalibFlag(ctlSubwWRM%localNumDam))
     WRMUnit%StorageCalibFlag = 0

     allocate (WRMUnit%INVc(ctlSubwWRM%localNumDam))
     WRMUnit%INVc = 0._r8
     allocate (WRMUnit%use_Irrig(ctlSubwWRM%localNumDam))
     WRMUnit%use_Irrig = 0
     allocate (WRMUnit%use_Elec(ctlSubwWRM%localNumDam))
     WRMUnit%use_Elec = 0
     allocate (WRMUnit%use_Supp(ctlSubwWRM%localNumDam))
     WRMUnit%use_Supp = 0
     allocate (WRMUnit%use_FCon(ctlSubwWRM%localNumDam))
     WRMUnit%use_FCon = 0
     allocate (WRMUnit%use_Fish(ctlSubwWRM%localNumDam))
     WRMUnit%use_Fish = 0
     allocate (WRMUnit%use_Rec(ctlSubwWRM%localNumDam))
     WRMUnit%use_Rec = 0
     allocate (WRMUnit%use_Navi(ctlSubwWRM%localNumDam))
     WRMUnit%use_Navi = 0

     allocate (WRMUnit%Withdrawal(ctlSubwWRM%localNumDam))
     WRMUnit%Withdrawal = 0
     allocate (WRMUnit%Conveyance(ctlSubwWRM%localNumDam))
     WRMUnit%Conveyance = 0
     allocate (WRMUnit%MthStOp(ctlSubwWRM%localNumDam))
     WRMUnit%MthStOp = 0
     allocate (WRMUnit%StorMthStOp(ctlSubwWRM%localNumDam))
     WRMUnit%StorMthStOp = 0
     allocate (WRMUnit%MthStFC(ctlSubwWRM%localNumDam))
     WRMUnit%MthStFC = 0
     allocate (WRMUnit%MthNdFC(ctlSubwWRM%localNumDam))
     WRMUnit%MthNdFC = 0
     allocate (WRMUnit%MthFCtrack(ctlSubwWRM%localNumDam))
     WRMUnit%MthFCtrack = 0

     allocate (StorWater%demand(begr:endr))
     StorWater%demand=0._r8
     allocate (StorWater%supply(begr:endr))
     StorWater%supply=0._r8
     allocate (StorWater%deficit(begr:endr))
     StorWater%deficit=0._r8

     allocate (StorWater%WithDemIrrig(begr:endr))
     StorWater%WithDemIrrig=0._r8
     allocate (StorWater%WithDemNonIrrig(begr:endr))
     StorWater%WithDemNonIrrig=0._r8
     allocate (StorWater%ConDemIrrig(begr:endr))
     StorWater%ConDemIrrig=0._r8
     allocate (StorWater%ConDemNonIrrig(begr:endr))
     StorWater%ConDemNonIrrig=0._r8
     allocate (StorWater%SuppIrrig(begr:endr))
     StorWater%SuppIrrig=0._r8
     allocate (StorWater%SuppNonIrrig(begr:endr))
     StorWater%SuppNonIrrig=0._r8
     allocate (StorWater%ReturnIrrig(begr:endr))
     StorWater%ReturnIrrig=0._r8
     allocate (StorWater%ReturnNonIrrig(begr:endr))
     StorWater%ReturnNonIrrig=0._r8
     allocate (StorWater%GWShareIrrig(begr:endr))
     StorWater%GWShareIrrig=0._r8
     allocate (StorWater%GWShareNonIrrig(begr:endr))
     StorWater%GWShareNonIrrig=0._r8

     allocate (StorWater%pre_release(ctlSubwWRM%localNumDam, 13))
     StorWater%pre_release = 0._r8
     allocate (StorWater%storage(ctlSubwWRM%localNumDam))
     StorWater%storage = 0._r8
     allocate (StorWater%release(ctlSubwWRM%localNumDam))
     StorWater%release = 0._r8
     allocate (StorWater%FCrelease(ctlSubwWRM%localNumDam))
     StorWater%FCrelease = 0._r8
     allocate (StorWater%pot_evap(begr:endr))
     StorWater%pot_evap=0._r8

     !allocate (TmpStoRelease(ctlSubwWRM%localNumDam,:,:))

     ! read the wm para

     ier = pio_inq_varid (ncid, name='DamID_Spatial', vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_int_dam , WRMUnit%mask, ier)
     if (masterproc) write(iulog,FORMI) trim(subname),' read DamID_Spatial',minval(WRMUnit%mask),maxval(WRMUnit%mask)
     call shr_sys_flush(iulog)

     ier = pio_inq_varid (ncid, name='RUNOFF_CAP'   , vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_dbl_dam , WRMUnit%INVc, ier)
     if (masterproc) write(iulog,FORMR) trim(subname),' read RUNOFF_CAP',minval(WRMUnit%INVc),maxval(WRMUnit%INVc)
     call shr_sys_flush(iulog)

     ier = pio_inq_varid (ncid, name='Year'         , vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_int_dam , WRMUnit%YEAR, ier)
     if (masterproc) write(iulog,FORMI) trim(subname),' read Year',minval(WRMUnit%YEAR),maxval(WRMUnit%YEAR)
     call shr_sys_flush(iulog)

     ier = pio_inq_varid (ncid, name='dam_hgt'      , vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_dbl_dam , WRMUnit%Height, ier)
     if (masterproc) write(iulog,FORMR) trim(subname),' read dam_hgt',minval(WRMUnit%Height),maxval(WRMUnit%Height)
     call shr_sys_flush(iulog)

     ier = pio_inq_varid (ncid, name='dam_len'      , vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_dbl_dam , WRMUnit%Length, ier)
     if (masterproc) write(iulog,FORMR) trim(subname),' read dam_len',minval(WRMUnit%Length),maxval(WRMUnit%Length)
     call shr_sys_flush(iulog)

     ier = pio_inq_varid (ncid, name='area_skm'     , vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_dbl_dam , WRMUnit%SurfArea, ier)
     if (masterproc) write(iulog,FORMR) trim(subname),' read area_skm',minval(WRMUnit%SurfArea),maxval(WRMUnit%SurfArea)
     call shr_sys_flush(iulog)

     ier = pio_inq_varid (ncid, name='cap_mcm'      , vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_dbl_dam , WRMUnit%StorCap, ier)
     if (masterproc) write(iulog,FORMR) trim(subname),' read cap_mcm',minval(WRMUnit%StorCap),maxval(WRMUnit%StorCap)
     call shr_sys_flush(iulog)

     WRMUnit%SurfArea = WRMUnit%SurfArea * 1e6
     WRMUnit%StorCap=WRMUnit%StorCap*1e6

     !in MCM
     !WRMUnit%StorMthStOp = WRMUnit%StorCap*0.9

     ier = pio_inq_varid (ncid, name='depth_m'       , vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_dbl_dam  , WRMUnit%Depth, ier)
     if (masterproc) write(iulog,FORMR) trim(subname),' read depth_m',minval(WRMUnit%Depth),maxval(WRMUnit%Depth)
     call shr_sys_flush(iulog)

     ier = pio_inq_varid (ncid, name='use_irri'      , vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_int_dam  , WRMUnit%use_Irrig, ier)
     if (masterproc) write(iulog,FORMI) trim(subname),' read use_irri',minval(WRMUnit%use_Irrig),maxval(WRMUnit%use_Irrig)
     call shr_sys_flush(iulog)

     ier = pio_inq_varid (ncid, name='use_elec'      , vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_int_dam  , WRMUnit%use_Elec, ier)
     if (masterproc) write(iulog,FORMI) trim(subname),' read use_elec',minval(WRMUnit%use_Elec),maxval(WRMUnit%use_Elec)
     call shr_sys_flush(iulog)

     ier = pio_inq_varid (ncid, name='use_supp'      , vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_int_dam  , WRMUnit%use_Supp, ier)
     if (masterproc) write(iulog,FORMI) trim(subname),' read use_supp',minval(WRMUnit%use_Supp),maxval(WRMUnit%use_Supp)
     call shr_sys_flush(iulog)

     ier = pio_inq_varid (ncid, name='use_fcon'      , vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_int_dam  , WRMUnit%use_FCon, ier)
     if (masterproc) write(iulog,FORMI) trim(subname),' read use_fcon',minval(WRMUnit%use_FCon),maxval(WRMUnit%use_FCon)
     call shr_sys_flush(iulog)

     ier = pio_inq_varid (ncid, name='use_recr'      , vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_int_dam  , WRMUnit%use_Rec, ier)
     if (masterproc) write(iulog,FORMI) trim(subname),' read use_recr',minval(WRMUnit%use_Rec),maxval(WRMUnit%use_Rec)
     call shr_sys_flush(iulog)

     ier = pio_inq_varid (ncid, name='use_navi'      , vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_int_dam  , WRMUnit%use_Navi, ier)
     if (masterproc) write(iulog,FORMI) trim(subname),' read use_navi',minval(WRMUnit%use_Navi),maxval(WRMUnit%use_Navi)
     call shr_sys_flush(iulog)

     ier = pio_inq_varid (ncid, name='use_fish'      , vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_int_dam  , WRMUnit%use_Fish, ier)
     if (masterproc) write(iulog,FORMI) trim(subname),' read use_fish',minval(WRMUnit%use_Fish),maxval(WRMUnit%use_Fish)
     call shr_sys_flush(iulog)

     ier = pio_inq_varid (ncid, name='withdraw'      , vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_dbl_dam  , WRMUnit%Withdrawal, ier)
     if (masterproc) write(iulog,FORMR) trim(subname),' read withdraw',minval(WRMUnit%Withdrawal),maxval(WRMUnit%Withdrawal)
     call shr_sys_flush(iulog)

     ier = pio_inq_varid (ncid, name='conveyance'    , vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_dbl_dam  , WRMUnit%Conveyance, ier)
     if (masterproc) write(iulog,FORMR) trim(subname),' read conveyance',minval(WRMUnit%Conveyance),maxval(WRMUnit%Conveyance)
     call shr_sys_flush(iulog)

! tcx should this flag be here?
     if ( ctlSubwWRM%RegulationFlag > 0 ) then

! tcx, hongyi code, delete
!        ! reading mean monthly flow data
!        allocate(temp2D_local_dbl(ctlSubwWRM%NDam,12))
!        call mosart_wm_readnc_dbl_2D(ncid, 'Qmon', temp2D_local_dbl, ctlSubwWRM%NDam,12)
!        do nr=begr,endr
!           if(WRMUnit%isDam(nr)>0) then
!              idam = WRMUnit%INVicell(nr)
!              nn=idam_global(idam)
!              WRMUnit%MeanMthFlow(idam,1:12) = temp2D_local_dbl(nn,1:12)
!              WRMUnit%MeanMthFlow(idam,13) = sum(WRMUnit%MeanMthFlow(idam,1:12))/12.0_r8
!            end if
!        end do
!        deallocate(temp2D_local_dbl)
!
!        ! reading in mean monthly water demand
!        allocate(temp2D_local_dbl(ctlSubwWRM%NDam,12))
!        call mosart_wm_readnc_dbl_2D(ncid, 'demand', temp2D_local_dbl, ctlSubwWRM%NDam,12)
!        do nr=begr,endr
!           if(WRMUnit%isDam(nr)>0) then
!              idam = WRMUnit%INVicell(nr)
!              nn=idam_global(idam)
!              WRMUnit%MeanMthDemand(idam,1:12) = temp2D_local_dbl(nn,1:12)
!              WRMUnit%MeanMthDemand(idam,1:12) = WRMUnit%MeanMthDemand(idam,1:12) * WRMUnit%Withdrawal(idam)
!              WRMUnit%MeanMthDemand(idam,13) = sum(WRMUnit%MeanMthDemand(idam,1:12))/12.0_r8
!           end if
!        end do
!        deallocate(temp2D_local_dbl)

        ! reading mean monthly flow data
        do n = 1,12
           ier = pio_inq_varid (ncid, name='Qmon', vardesc=vardesc)
           frame = n
           call pio_setframe(vardesc,frame)
           call pio_read_darray(ncid, vardesc, iodesc_dbl_dam, WRMUnit%MeanMthFlow(:,n), ier)
           if (masterproc) write(iulog,FORMR) trim(subname),' read Qmon',minval(WRMUnit%MeanMthFlow(:,n)),maxval(WRMUnit%MeanMthFlow(:,n))
           call shr_sys_flush(iulog)
        enddo
        WRMUnit%MeanMthFlow(idam,13) = sum(WRMUnit%MeanMthFlow(idam,1:12))/12.0_r8
        if (masterproc) write(iulog,FORMR) trim(subname),' MeanMthFlow avg',minval(WRMUnit%MeanMthFlow(:,13)),maxval(WRMUnit%MeanMthFlow(:,13))

        ! reading in mean monthly water demand
        do n = 1,12
           ier = pio_inq_varid (ncid, name='demand', vardesc=vardesc)
           frame = n
           call pio_setframe(vardesc,frame)
           call pio_read_darray(ncid, vardesc, iodesc_dbl_dam, WRMUnit%MeanMthDemand(:,n), ier)
           if (masterproc) write(iulog,FORMR) trim(subname),' read demand',minval(WRMUnit%MeanMthDemand(:,n)),maxval(WRMUnit%MeanMthDemand(:,n))
           call shr_sys_flush(iulog)
        enddo
        WRMUnit%MeanMthDemand(idam,13) = sum(WRMUnit%MeanMthDemand(idam,1:12))/12.0_r8
        if (masterproc) write(iulog,FORMR) trim(subname),' MeanMthDemand avg',minval(WRMUnit%MeanMthDemand(:,13)),maxval(WRMUnit%MeanMthDemand(:,13))

        !initialize constant monthly pre-release based on longterm mean flow and demand (Biemans 2011) 
        do iunit=1,ctlSubwWRM%NDam
           do mth=1,12
              StorWater%pre_release(iunit,mth) = WRMUnit%MeanMthFlow(iunit,13)
           end do
           if ( WRMUnit%MeanMthDemand(iunit,13) >= (0.5_r8*WRMUnit%MeanMthFlow(iunit,13)) .and. WRMUnit%MeanMthFlow(iunit,13)>0._r8 ) then
              do mth=1,12
                 StorWater%pre_release(iunit,mth) = WRMUnit%MeanMthDemand(iunit,mth)/10._r8 + 9._r8/10._r8*WRMUnit%MeanMthFlow(iunit,13)*WRMUnit%MeanMthDemand(iunit,mth)/WRMUnit%MeanMthDemand(iunit, 13)
                 !TEST
                 !StorWater%pre_release(iunit,mth) = WRMUnit%MeanMthDemand(iunit,mth)/10._r8 + 9._r8/10._r8*WRMUnit%MeanMthFlow(iunit,13)*WRMUnit%MeanMthDemand(iunit,mth)/WRMUnit%MeanMthDemand(iunit, 13)*.5_r8
              end do
           else 
              do mth=1,12
                 if ( (WRMUnit%MeanMthFlow(iunit,13) + WRMUnit%MeanMthDemand(iunit,mth) - WRMUnit%MeanMthDemand(iunit,13))>0 ) then
                    StorWater%pre_release(iunit, mth) = WRMUnit%MeanMthFlow(iunit,13) + WRMUnit%MeanMthDemand(iunit,mth) - WRMUnit%MeanMthDemand(iunit,13)
                 endif 
                 ! test 2
                 !StorWater%pre_release(iunit, mth) = WRMUnit%MeanMthFlow(iunit,13)*0.5_r8 + WRMUnit%MeanMthDemand(iunit,mth) - WRMUnit%MeanMthDemand(iunit,13)
                 !TEST use pseudo regulated flow
                 !StorWater%pre_release(iunit, mth) = WRMUnit%MeanMthFlow(iunit,13)*.5_r8 + WRMUnit%MeanMthDemand(iunit,mth) - WRMUnit%MeanMthDemand(iunit,13)
              end do
           end if

           ! initialize storage in each reservoir - arbitrary 50%
           StorWater%storage(iunit) = 0.9_r8 * WRMUnit%StorCap(iunit)   
           if (WRMUnit%StorageCalibFlag(iunit).eq.1) then
              StorWater%storage(iunit)  = WRMUnit%StorTarget(iunit,13)
           endif
           if ( WRMUnit%StorCap(iunit) <= 0 ) then
              write(iulog,*) subname, "Error negative max cap for reservoir ", iunit, WRMUnit%StorCap(iunit)
              call shr_sys_abort(subname//' ERROR: negative max cap for reservoir')
           end if
        end do
        write(iulog,*) subname, "storage ",StorWater%pre_release(1,1), StorWater%storage(1)

        ! initialize start of the operationnal year based on long term simulation
        call WRM_init_StOp_FC

!tcx, move read of gridID_from_Dam and initialization to netcdf part, global data on all tasks
!        allocate(temp_gridID_from_Dam(ctlSubwWRM%NDam, maxNumDependentGrid))
!        call mosart_wm_readnc_int_2D(ncid, 'gridID_from_Dam', temp_gridID_from_Dam, ctlSubwWRM%NDam, maxNumDependentGrid)
!        do idam = 1,ctlSubwWRM%localNumDam
!           nn = idam_global(idam)
!           ntotal = 0
!           do nd=1,maxNumDependentGrid
!              if (temp_gridID_from_Dam(nn,nd) > 0) then
!                 ntotal = ntotal + 1
!                 WRMUnit%dam_depend(idam,ntotal) = temp_gridID_from_Dam(nn,nd)
!              endif
!           enddo
!           if (ntotal > WRMUnit%dam_Ndepend(idam) .or. ntotal < WRMUnit%dam_Ndepend(idam)-1) then ! sometimes in the dependency data, the basin outlet is included for a dam as its dependent grid
!              write(iulog,"(a,e12.6,i6,i6)"), "Attention reading gridID_from_Dam ", WRMUnit%StorCap(idam),WRMUnit%dam_Ndepend(idam), ntotal
!              !call endrun
!           end if
!           WRMUnit%dam_Ndepend(idam) = ntotal
!        enddo

! tcx old code looks wrong, seems to be grabbing local dams only
!        do nr=begr,endr
!           if(WRMUnit%isDam(nr)>0) then
!              idam = WRMUnit%INVicell(nr)
!              nn=idam_global(idam)
!              ntotal = 0
!              do nd=1,maxNumDependentGrid
!                 iunit = temp_gridID_from_Dam(nn,nd)
!                 if(iunit > 0) then
!                    n = rglo2gdc(iunit)
!                    if(n >= begr .and. n <= endr) then
!                       ntotal = ntotal + 1
!                       WRMUnit%dam_depend(idam,ntotal) = n
!                    end if
!                end if
!             end do
!             if(ntotal > WRMUnit%dam_Ndepend(idam) .or. ntotal < WRMUnit%dam_Ndepend(idam)-1) then ! sometimes in the dependency data, the basin outlet is included for a dam as its dependent grid
!                write(iulog,"(a,e12.6,i6,i6)"), "Attention reading gridID_from_Dam ", WRMUnit%StorCap(idam),WRMUnit%dam_Ndepend(idam), ntotal
!                !call endrun
!             end if
!             WRMUnit%dam_Ndepend(idam) = ntotal
!          end if
!        end do

! tcraig above is hongyi, below is what's left of Nathalie. reconcile this???
!! need to adjust for reservoir with zero inflow, do not  need to read the remaining
!           if ( WRMUnit%MeanMthFlow(iunit,13) <= 0._r8 ) then
!              WRMUnit%dam_Ndepend(iunit) = 0 ! this reservoir will not provide water to any subw, relieve database
!           end if
!
!           do j=1,WRMUnit%dam_Ndepend(iunit)
!              call split(stemp,' ',stmp1)
!              call str2num(stmp1, ctlSubwWRM%localNumDam, ierror) 
!!need additionnal check due to regionalization, need to remove non existing grid cell NV
!              if (  WRMUnit%INVisubw(ctlSubwWRM%localNumDam) .lt. 1 ) then
!                 WRMUnit%dam_Ndepend(iunit) = WRMUnit%dam_Ndepend(iunit) - 1
!              else
!                 WRMUnit%dam_depend(iunit,j) = nd 
!              endif
!           end do

        !initialize subw dependencies
        ier = pio_inq_varid (ncid, name='num_Dam2Grid', vardesc=vardesc)
        call pio_read_darray(ncid, vardesc, iodesc_int , WRMUnit%subw_Ndepend, ier)
        if (masterproc) write(iulog,FORMI) trim(subname),' read num_Dam2Grid',minval(WRMUnit%subw_Ndepend),maxval(WRMUnit%subw_Ndepend)
        call shr_sys_flush(iulog)

!        call MOSART_read_int(ncid, 'num_Dam2Grid', WRMUnit%subw_Ndepend)
!tcx this is highly suspicious
!tcx check this, it looks like it's operating only local decomp, should operate across all pes
!tcx this looks like an optimization to reduce cost later for gridpoints that have no MeanMthFlow
!tcx defer and review
        do idam=1,ctlSubwWRM%localNumDam
           ! need to adjust for reservoir with zero inflow, do not  need to read the remaining
           if (WRMUnit%MeanMthFlow(idam,13) <= 0._r8 .and. WRMUnit%dam_Ndepend(idam) > 1) then
              write(iulog,"(a,i8,i8)"), "Attention: Reservoir with zero inflow while non-zero dependent grids", WRMUnit%icell(idam), WRMUnit%dam_Ndepend(idam)
     !tcx         do nn=1,WRMUnit%dam_Ndepend(idam)
     !tcx            nr = WRMUnit%dam_depend(idam,nn)
     !tcx            if(nr >= begr .and. nr <= endr) then
     !tcx               WRMUnit%subw_Ndepend(nr) = WRMUnit%subw_Ndepend(nr) - 1
     !tcx            end if
     !tcx         end do
     !tcx         WRMUnit%dam_Ndepend(idam) = 0 ! this reservoir will not provide water to any subw, relieve database
           end if
        end do

!        allocate(temp1D_int(begr:endr))
!        temp1D_int = 0
!        do idam=1,ctlSubwWRM%localNumDam
!           do nn=1,WRMUnit%dam_Ndepend(idam)
!              nr = WRMUnit%dam_depend(idam,nn)
!              if(nr >= begr .and. nr <= endr) then
!                 temp1D_int(nr) = temp1D_int(nr) + 1
!! tcx not used anywhere as far as i can tell, so get rid of this block
!                 WRMUnit%subw_depend(nr,temp1D_int(nr)) = inverse_INVicell(idam)
!              else
!                 write(iulog,"(a15,i8,i6,i6)"), "Inconsistent dam_depend and subw_depend  ", begr, endr, nr
!              end if
!           end do
!        end do

!        !check the dependence database consistencies
!        do nr=begr,endr
!           if(.not.(temp1D_int(nr).eq.WRMUnit%subw_Ndepend(nr))) then
!              write(iulog,"(a35,i8,i6,i6)"), "Attention processing gridID_Dam2Grid ", TUnit%ID0(nr), WRMUnit%subw_Ndepend(nr), temp1D_int(nr)
!           end if
!        end do

        !check the dependence database consistencies
        do idam=1,ctlSubwWRM%localNumDam
           do j=1,WRMUnit%dam_Ndepend(idam)
              !if(WRMUnit%dam_depend(idam,j).eq.0) then
              !   WRMUnit%dam_depend(idam,j) = WRMUnit%dam_depend(idam,j) * 1
              !end if
              idepend = WRMUnit%dam_depend(idam,j)
              if ( idepend .lt. 0 ) then
                 write(iulog,*) subname,"Error checking dependency, zero idepend", idam, WRMUnit%dam_Ndepend(idam), j, idepend , WRMUnit%dam_depend(idam,j)
                 call shr_sys_abort(subname//' ERROR: zero idepend')
              endif
!tcx this should be accumlating global data, not just local data, idepend could be outside begr:endr
              WRMUnit%TotStorCapDepend(idepend) = WRMUnit%TotStorCapDepend(idepend) + WRMUnit%StorCap(idam)
              WRMUnit%TotInflowDepend(idepend) = WRMUnit%TotInflowDepend(idepend) + WRMUnit%MeanMthFlow(idam,13)
           end do
        end do

     end if !Regulation Flag

     call ncd_pio_closefile(ncid)

     ! check
     write(iulog,*) subname, "Done with WM init ..."
     !write(iulog,*) subname, WRMUnit%DamName(1), WRMUnit%Surfarea(1)
     !write(iulog,*) subname,WRMUnit%mask(1), WRMUnit%icell(1) 
     !write(iulog,*) subname, WRMUnit%dam_Ndepend(1), WRMUnit%dam_depend(1,2)
     !write(iulog,*) subname, "sub = 49",  TUnit%icell(49, 1),WRMUnit%subw_Ndepend(49),  WRMUnit%subw_depend(49,1) 
  end subroutine WRM_init

end MODULE WRM_subw_IO_mod
