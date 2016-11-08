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

     integer :: nr, nd, ng, mth        ! local loop indices
     integer :: cnt, idam, ntotal
     integer :: begr, endr, maxnumdependentgrid
     integer :: ncdfid, did, dids(2), dsize(1), dsizes(2), ier, varid
     type(file_desc_t):: ncid       ! netcdf file
     type(var_desc_t) :: vardesc    ! netCDF variable description
     type(io_desc_t)  :: iodesc_int_grd2grd ! pio io desc, global grid to local grid
     type(io_desc_t)  :: iodesc_dbl_grd2grd ! pio io desc, global grid to local grid
     type(io_desc_t)  :: iodesc_int_grd2dam  ! pio io desc, global grid to local dam
     type(io_desc_t)  :: iodesc_dbl_grd2dam  ! pio io desc, global grid to local dam
     type(io_desc_t)  :: iodesc_int_dam2dam   ! pio io desc, global dam to local dam
     type(io_desc_t)  :: iodesc_dbl_dam2dam   ! pio io desc, global dam to local dam
     integer(kind=PIO_OFFSET_KIND) :: frame
     real(r8) :: peak, prorata , mn, mx                  ! peak value to define the start of operationalyr
     integer :: sgn,curr_sgn, nsc, ct, ct_mx, mth_op , idepend      ! number of sign change
     integer, pointer  :: compDOF(:)           ! pio decomp
     integer, pointer  :: idam_ldam2gdam(:)    ! local dam index to global dam index
     integer, pointer  :: idam_ldam2ggrd(:)    ! local dam index to global grid index
     integer, pointer  :: UnitID_1D(:)         ! list of dam IDs, 1:ndams
     integer, pointer  :: temp_gridID_from_Dam(:,:)
     integer, pointer  :: temp_subw_Ndepend(:) 
     integer, pointer  :: dependcnt(:)
     character(len=*),parameter :: subname='WRM_init'
     character(len=*),parameter :: FORMI = '(2A,6i13)'
     character(len=*),parameter :: FORMR = '(2A,2g15.7)'

     begr = rtmCTL%begr
     endr = rtmCTL%endr

     if (masterproc) write(iulog,*) subname,' begr endr = ',begr,endr
     call shr_sys_flush(iulog)

     ctlSubwWRM%paraPath = '/pic/scratch/tcraig/wm_data/'
     ctlSubwWRM%paraFile = 'US_reservoir_8th_NLDAS2.nc'
     ctlSubwWRM%demandPath = '/pic/scratch/tcraig/wm_data/NLDAS2_GCAM_water_demand_'

     !-------------------
     ! read input dataset
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
     do nr = begr,endr
        cnt = cnt + 1
!tcx, read in data based on mosart decomposition
! THIS ASSUMES mosart and wrm use same grid indexing
        compDOF(cnt) = rtmCTL%gindex(nr)
     enddo
     ier = pio_inq_varid   (ncid, name='unit_ID', vardesc=vardesc)
     ier = pio_inq_vardimid(ncid, vardesc, dids)
     ier = pio_inq_dimlen  (ncid, dids(1),dsizes(1))
     ier = pio_inq_dimlen  (ncid, dids(2),dsizes(2))

     if (masterproc) write(iulog,*) subname,' dsizes = ',dsizes
     write(iulog,*) subname,' compDOF = ',minval(compDof),maxval(compDOF)
     call shr_sys_flush(iulog)

! decomp for lon/lat to local gridcells
     call pio_initdecomp(pio_subsystem, pio_double, dsizes, compDOF, iodesc_dbl_grd2grd)
     call pio_initdecomp(pio_subsystem, pio_int   , dsizes, compDOF, iodesc_int_grd2grd)
     deallocate(compdof)

     !--- read in DamIDs

     allocate (WRMUnit%isDam(begr:endr))
     WRMUnit%isDam = 0
     ier = pio_inq_varid (ncid, name='unit_ID', vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_int_grd2grd, WRMUnit%isDam, ier)
     if (masterproc) write(iulog,FORMI) trim(subname),' read unit_ID ',minval(WRMunit%isDam),maxval(WRMUnit%isDam)
     call shr_sys_flush(iulog)

     !--- count dams and set up dam indexing

     ctlSubwWRM%localNumDam = 0
     do nr = begr,endr
        if (WRMUnit%isDam(nr) > 0) then
           ! tcx tcraig, this is not the right check.  the damID should be independent of the grid at some future time
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

     !--- initialize icell, INVicell and dam indexing

     allocate(WRMUnit%INVicell(begr:endr))
     WRMUnit%INVicell =-99 
     allocate(WRMUnit%icell(ctlsubwWRM%localNumDam))
     WRMUnit%icell = 0
     allocate(idam_ldam2gdam(ctlSubwWRM%localNumDam))
     idam_ldam2gdam = -99
     allocate(idam_ldam2ggrd(ctlSubwWRM%localNumDam))
     idam_ldam2ggrd = -99
     cnt = 0
     do nr = begr,endr
        if (WRMUnit%isDam(nr) > 0) then
           cnt = cnt + 1
           WRMUnit%INVicell(nr) = cnt  ! local gridcell index to local dam index
           WRMUnit%icell(cnt) = nr     ! local dam index to local gridcell index
           do nd = 1,ctlSubwWRM%NDam
              if (WRMUnit%isDam(nr) .eq. UnitID_1D(nd)) then
!                 write(iulog,'(a,5i10)') subname//' idam_glo=',nr,cnt,WRMUnit%isDam(nr),nd,UnitID_1D(nd)
                 ! check if a dam was already found, if so abort
                 if (idam_ldam2ggrd(cnt) > 0) then
                    write(iulog,*) subname,' ERROR: dam ID appears twice in UnitID_1D'
                    call shr_sys_abort(subname//' ERROR dam ID appears twice in UnitID_1D')
                 endif
                 idam_ldam2gdam(cnt) = nd                 ! local dam index to global dam index
                 idam_ldam2ggrd(cnt) = WRMUnit%isDam(nr)  ! local dam index to global gridcell index
              end if
           end do
        end if
     end do

     write(iulog,*) subname,' icell    = ',minval(WRMUnit%icell),maxval(WRMUnit%icell)
     write(iulog,*) subname,' INVicell = ',minval(WRMUnit%INVicell),maxval(WRMUnit%INVicell)
     write(iulog,*) subname,' idam_ldam2gdam = ',minval(idam_ldam2gdam),maxval(idam_ldam2gdam)
     write(iulog,*) subname,' idam_ldam2ggrd = ',minval(idam_ldam2ggrd),maxval(idam_ldam2ggrd)

     ! --- initialize iodesc for grid2dam and dam2dam
     ! decomp for lon/lat to local dam
     call pio_initdecomp(pio_subsystem, pio_double, dsizes, idam_ldam2ggrd, iodesc_dbl_grd2dam)
     call pio_initdecomp(pio_subsystem, pio_int   , dsizes, idam_ldam2ggrd, iodesc_int_grd2dam)
     ! decomp for global dam to local dam
     dsize(1) = ctlSubwWRM%NDam
     call pio_initdecomp(pio_subsystem, pio_double, dsize , idam_ldam2gdam, iodesc_dbl_dam2dam)
     call pio_initdecomp(pio_subsystem, pio_int   , dsize , idam_ldam2gdam, iodesc_int_dam2dam)

     !--- initialize dam dependencies

     allocate (WRMUnit%dam_Ndepend(ctlSubwWRM%localNumDam))
     WRMUnit%dam_Ndepend = 0        !

     ier = pio_inq_varid (ncid, name='numGrid_from_Dam', vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_int_grd2dam , WRMUnit%dam_Ndepend, ier)
     write(iulog,FORMI) trim(subname),' read numGrid_from_Dam',minval(WRMUnit%dam_Ndepend),maxval(WRMUnit%dam_Ndepend)
     call shr_sys_flush(iulog)

     do idam = 1,ctlSubwWRM%localNumDam
        nd = idam_ldam2gdam(idam)
        ntotal = 0
        do ng = 1,maxNumDependentGrid
           if (temp_gridID_from_Dam(nd,ng) > 0) then
              ntotal = ntotal + 1
           endif
        enddo
!tcx        write(iulog,*) subname,' compute dam_depend ntotal ',idam,ntotal
        ! check that Ndepend and dam_depend are consistent
        ! sometimes in the dependency data, the basin outlet is included for a dam as its dependent grid
        if (ntotal > WRMUnit%dam_Ndepend(idam) .or. ntotal < WRMUnit%dam_Ndepend(idam)-1) then
           write(iulog,"(a,3i8)"), subname//"ERROR: Ndepend/dam_depend inconsistency1",idam, ntotal,WRMUnit%dam_Ndepend(idam)
           call shr_sys_abort(subname//' ERROR: numGrid_from_Dam not consistent with gridID_from_Dam list')
        end if
        WRMUnit%dam_Ndepend(idam) = ntotal
     enddo

     allocate (WRMUnit%dam_depend(ctlSubwWRM%localNumDam,maxval(WRMUnit%dam_Ndepend)))
     WRMUnit%dam_depend = 0        !

     do idam = 1,ctlSubwWRM%localNumDam
        nd = idam_ldam2gdam(idam)
        ntotal = 0
        do ng = 1,maxNumDependentGrid
           if (temp_gridID_from_Dam(nd,ng) > 0) then
              ntotal = ntotal + 1
              if (ntotal > WRMUnit%dam_Ndepend(idam)) then
                 write(iulog,"(a,3i8)"), subname//"ERROR: Ndepend/dam_depend inconsistency2",idam, WRMUnit%dam_Ndepend(idam), ntotal
                 call shr_sys_abort(subname//' ERROR: ntotal gt dam_Ndepend')
              end if
              WRMUnit%dam_depend(idam,ntotal) = temp_gridID_from_Dam(nd,ng)
           endif
        enddo
        if (ntotal /= WRMUnit%dam_Ndepend(idam)) then
           write(iulog,"(a,3i8)"), subname//"ERROR: Ndepend/dam_depend inconsistency3",idam, WRMUnit%dam_Ndepend(idam), ntotal
           call shr_sys_abort(subname//' ERROR: ntotal ne dam_Ndepend')
        end if
     enddo

     do idam = 1,ctlSubwWRM%localNumDam
        do ng = 1,WRMUnit%dam_Ndepend(idam)
           if (WRMUnit%dam_depend(idam,ng) <= 0) then
              write(iulog,'(2a,4i12)') subname,' ERROR dam_depend check1 ',idam,WRMUnit%dam_Ndepend(idam),ng,WRMUnit%dam_depend(idam,ng)
           endif
        enddo
     enddo

     !--- initialize inverse of dam_depend

     allocate(temp_subw_Ndepend(begr:endr))
     ier = pio_inq_varid (ncid, name='num_Dam2Grid', vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_int_grd2grd , temp_subw_Ndepend, ier)
     write(iulog,FORMI) trim(subname),' read num_Dam2Grid',minval(temp_subw_Ndepend),maxval(temp_subw_Ndepend)
     call shr_sys_flush(iulog)

     allocate (dependcnt(begr:endr))
     dependcnt = 0
     do nd = 1,ctlSubwWRM%NDam
        do ng = 1,maxNumDependentGrid
           if (temp_gridID_from_Dam(nd,ng) > 0) then
              do nr = begr,endr
                 if (temp_gridID_from_Dam(nd,ng) == rtmCTL%gindex(nr)) then
                    dependcnt(nr) = dependcnt(nr) + 1
                 endif
              enddo
           endif
        enddo
     enddo
     write(iulog,FORMI) trim(subname),' compute dependcnt',minval(dependcnt),maxval(dependcnt)

     allocate (WRMUnit%subw_Ndepend(begr:endr))
     WRMUnit%subw_Ndepend = 0
     allocate (WRMUnit%subw_depend(begr:endr,maxval(dependcnt)))
     WRMUnit%subw_depend = 0
     allocate (WRMUnit%subw_damlist(ctlSubwWRM%NDam))
     WRMUnit%subw_damlist = 0

     WRMUnit%subw_Ndepend = 0
     do nd = 1,ctlSubwWRM%NDam
        do ng = 1,maxNumDependentGrid
           if (temp_gridID_from_Dam(nd,ng) > 0) then
              do nr = begr,endr
                 if (temp_gridID_from_Dam(nd,ng) == rtmCTL%gindex(nr)) then
                    WRMUnit%subw_Ndepend(nr) = WRMUnit%subw_Ndepend(nr) + 1
                    if (WRMUnit%subw_Ndepend(nr) > dependcnt(nr)) then
                       write(iulog,*) subname,' ERROR: dependcnt ',nd,ng,nr,dependcnt(nr),WRMUnit%subw_Ndepend(nr)
                       call shr_sys_abort(subname//' ERROR: dependcnt')
                    endif
                    WRMUnit%subw_depend(nr,WRMUnit%subw_Ndepend(nr)) = nd
                    WRMUnit%subw_damlist(nd) = 1
                 endif
              enddo
           endif
        enddo
     enddo
     do nr = begr,endr
        if (WRMUnit%subw_Ndepend(nr) /= dependcnt(nr)) then
           write(iulog,*) subname,' ERROR: dependcnt,subw_Ndepend error ',nr,dependcnt(nr),WRMUnit%subw_Ndepend(nr)
           call shr_sys_abort(subname//' ERROR: dependcnt,subw_Ndepend')
        endif
     enddo
     deallocate(dependcnt)

     write(iulog,FORMI) trim(subname),' compute WRMUnit%subw_Ndepend',minval(WRMUnit%subw_Ndepend),maxval(WRMUnit%subw_Ndepend),sum(WRMUnit%subw_Ndepend)
     write(iulog,FORMI) trim(subname),' compute WRMUnit%subw_damlist',minval(WRMUnit%subw_damlist),maxval(WRMUnit%subw_damlist),sum(WRMUnit%subw_damlist)
     write(iulog,FORMI) trim(subname),' compute WRMUnit%subw_depend',minval(WRMUnit%subw_depend),maxval(WRMUnit%subw_depend),sum(WRMUnit%subw_depend)
     do nr = begr,endr
        if (WRMUnit%subw_Ndepend(nr) > 0) then
!tcx           write(iulog,*) trim(subname),' compute subw_depend',nr,WRMUnit%subw_Ndepend(nr),minval(WRMUnit%subw_depend(nr,1:WRMUnit%subw_Ndepend(nr))),maxval(WRMUnit%subw_depend(nr,1:WRMUnit%subw_Ndepend(nr)))
        else
!tcx           write(iulog,*) trim(subname),' compute subw_depend',nr,WRMUnit%subw_Ndepend(nr)
        endif
     enddo

     ! compare num_Dam2Grid to subw_Ndepend computed from temp_gridID_from_Dam
     ! tcx, tcraig this should probabaly be removed and num_Dam2Grid does not need to be on input file
     do nr = begr,endr
        if (temp_subw_Ndepend(nr) /= WRMUnit%subw_Ndepend(nr)) then
           write(iulog,*) subname,' ERROR: num_Dam2Grid not consistent with temp_gridID_from_Dam'
        endif
     enddo
     deallocate(temp_subw_Ndepend)
     deallocate(temp_gridID_from_Dam)

     !--- initial dam data

     allocate (WRMUnit%DamName(ctlSubwWRM%localNumDam))
     WRMUnit%DamName = 'undefined'
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
     call pio_read_darray(ncid, vardesc, iodesc_int_grd2dam , WRMUnit%mask, ier)
     if (masterproc) write(iulog,FORMI) trim(subname),' read DamID_Spatial',minval(WRMUnit%mask),maxval(WRMUnit%mask)
     call shr_sys_flush(iulog)

     ier = pio_inq_varid (ncid, name='RUNOFF_CAP'   , vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_dbl_grd2dam , WRMUnit%INVc, ier)
     if (masterproc) write(iulog,FORMR) trim(subname),' read RUNOFF_CAP',minval(WRMUnit%INVc),maxval(WRMUnit%INVc)
     call shr_sys_flush(iulog)

     ier = pio_inq_varid (ncid, name='Year'         , vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_int_grd2dam , WRMUnit%YEAR, ier)
     if (masterproc) write(iulog,FORMI) trim(subname),' read Year',minval(WRMUnit%YEAR),maxval(WRMUnit%YEAR)
     call shr_sys_flush(iulog)

     ier = pio_inq_varid (ncid, name='dam_hgt'      , vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_dbl_grd2dam , WRMUnit%Height, ier)
     if (masterproc) write(iulog,FORMR) trim(subname),' read dam_hgt',minval(WRMUnit%Height),maxval(WRMUnit%Height)
     call shr_sys_flush(iulog)

     ier = pio_inq_varid (ncid, name='dam_len'      , vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_dbl_grd2dam , WRMUnit%Length, ier)
     if (masterproc) write(iulog,FORMR) trim(subname),' read dam_len',minval(WRMUnit%Length),maxval(WRMUnit%Length)
     call shr_sys_flush(iulog)

     ier = pio_inq_varid (ncid, name='area_skm'     , vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_dbl_grd2dam , WRMUnit%SurfArea, ier)
     if (masterproc) write(iulog,FORMR) trim(subname),' read area_skm',minval(WRMUnit%SurfArea),maxval(WRMUnit%SurfArea)
     call shr_sys_flush(iulog)

     ier = pio_inq_varid (ncid, name='cap_mcm'      , vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_dbl_grd2dam , WRMUnit%StorCap, ier)
     if (masterproc) write(iulog,FORMR) trim(subname),' read cap_mcm',minval(WRMUnit%StorCap),maxval(WRMUnit%StorCap)
     call shr_sys_flush(iulog)

     WRMUnit%SurfArea = WRMUnit%SurfArea * 1e6
     WRMUnit%StorCap=WRMUnit%StorCap*1e6

     !in MCM
     !WRMUnit%StorMthStOp = WRMUnit%StorCap*0.9

     ier = pio_inq_varid (ncid, name='depth_m'       , vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_dbl_grd2dam  , WRMUnit%Depth, ier)
     if (masterproc) write(iulog,FORMR) trim(subname),' read depth_m',minval(WRMUnit%Depth),maxval(WRMUnit%Depth)
     call shr_sys_flush(iulog)

     ier = pio_inq_varid (ncid, name='use_irri'      , vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_int_grd2dam  , WRMUnit%use_Irrig, ier)
     if (masterproc) write(iulog,FORMI) trim(subname),' read use_irri',minval(WRMUnit%use_Irrig),maxval(WRMUnit%use_Irrig)
     call shr_sys_flush(iulog)

     ier = pio_inq_varid (ncid, name='use_elec'      , vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_int_grd2dam  , WRMUnit%use_Elec, ier)
     if (masterproc) write(iulog,FORMI) trim(subname),' read use_elec',minval(WRMUnit%use_Elec),maxval(WRMUnit%use_Elec)
     call shr_sys_flush(iulog)

     ier = pio_inq_varid (ncid, name='use_supp'      , vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_int_grd2dam  , WRMUnit%use_Supp, ier)
     if (masterproc) write(iulog,FORMI) trim(subname),' read use_supp',minval(WRMUnit%use_Supp),maxval(WRMUnit%use_Supp)
     call shr_sys_flush(iulog)

     ier = pio_inq_varid (ncid, name='use_fcon'      , vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_int_grd2dam  , WRMUnit%use_FCon, ier)
     if (masterproc) write(iulog,FORMI) trim(subname),' read use_fcon',minval(WRMUnit%use_FCon),maxval(WRMUnit%use_FCon)
     call shr_sys_flush(iulog)

     ier = pio_inq_varid (ncid, name='use_recr'      , vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_int_grd2dam  , WRMUnit%use_Rec, ier)
     if (masterproc) write(iulog,FORMI) trim(subname),' read use_recr',minval(WRMUnit%use_Rec),maxval(WRMUnit%use_Rec)
     call shr_sys_flush(iulog)

     ier = pio_inq_varid (ncid, name='use_navi'      , vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_int_grd2dam  , WRMUnit%use_Navi, ier)
     if (masterproc) write(iulog,FORMI) trim(subname),' read use_navi',minval(WRMUnit%use_Navi),maxval(WRMUnit%use_Navi)
     call shr_sys_flush(iulog)

     ier = pio_inq_varid (ncid, name='use_fish'      , vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_int_grd2dam  , WRMUnit%use_Fish, ier)
     if (masterproc) write(iulog,FORMI) trim(subname),' read use_fish',minval(WRMUnit%use_Fish),maxval(WRMUnit%use_Fish)
     call shr_sys_flush(iulog)

     ier = pio_inq_varid (ncid, name='withdraw'      , vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_dbl_grd2dam  , WRMUnit%Withdrawal, ier)
     if (masterproc) write(iulog,FORMR) trim(subname),' read withdraw',minval(WRMUnit%Withdrawal),maxval(WRMUnit%Withdrawal)
     call shr_sys_flush(iulog)

     ier = pio_inq_varid (ncid, name='conveyance'    , vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_dbl_grd2dam  , WRMUnit%Conveyance, ier)
     if (masterproc) write(iulog,FORMR) trim(subname),' read conveyance',minval(WRMUnit%Conveyance),maxval(WRMUnit%Conveyance)
     call shr_sys_flush(iulog)


! tcx should this flag be here?
!     if ( ctlSubwWRM%RegulationFlag > 0 ) then

        !--- read mean monthly flow data
        do mth = 1,12
           ier = pio_inq_varid (ncid, name='Qmon', vardesc=vardesc)
           frame = mth
           call pio_setframe(ncid,vardesc,frame)
           call pio_read_darray(ncid, vardesc, iodesc_dbl_dam2dam, WRMUnit%MeanMthFlow(:,mth), ier)
           if (masterproc) write(iulog,FORMR) trim(subname),' read Qmon',minval(WRMUnit%MeanMthFlow(:,mth)),maxval(WRMUnit%MeanMthFlow(:,mth))
           call shr_sys_flush(iulog)
        enddo
        do idam = 1,ctlSubwWRM%LocalNumDam
           WRMUnit%MeanMthFlow(idam,13) = sum(WRMUnit%MeanMthFlow(idam,1:12))/12.0_r8
        enddo
        if (masterproc) write(iulog,FORMR) trim(subname),' MeanMthFlow avg',minval(WRMUnit%MeanMthFlow(:,13)),maxval(WRMUnit%MeanMthFlow(:,13))

        !--- read in mean monthly water demand
        do mth = 1,12
           ier = pio_inq_varid (ncid, name='demand', vardesc=vardesc)
           frame = mth
           call pio_setframe(ncid,vardesc,frame)
           call pio_read_darray(ncid, vardesc, iodesc_dbl_dam2dam, WRMUnit%MeanMthDemand(:,mth), ier)
           if (masterproc) write(iulog,FORMR) trim(subname),' read demand',minval(WRMUnit%MeanMthDemand(:,mth)),maxval(WRMUnit%MeanMthDemand(:,mth))
           call shr_sys_flush(iulog)
        enddo
        do idam = 1,ctlSubwWRM%LocalNumDam
           WRMUnit%MeanMthDemand(idam,13) = sum(WRMUnit%MeanMthDemand(idam,1:12))/12.0_r8
        enddo
        if (masterproc) write(iulog,FORMR) trim(subname),' MeanMthDemand avg',minval(WRMUnit%MeanMthDemand(:,13)),maxval(WRMUnit%MeanMthDemand(:,13))

        !--- initialize constant monthly pre-release based on longterm mean flow and demand (Biemans 2011) 

        do idam=1,ctlSubwWRM%LocalNumDam
           do mth=1,12
              StorWater%pre_release(idam,mth) = WRMUnit%MeanMthFlow(idam,13)
           end do
           if ( WRMUnit%MeanMthDemand(idam,13) >= (0.5_r8*WRMUnit%MeanMthFlow(idam,13)) .and. WRMUnit%MeanMthFlow(idam,13) > 0._r8 ) then
              do mth=1,12
                 StorWater%pre_release(idam,mth) = WRMUnit%MeanMthDemand(idam,mth)/10._r8 + 9._r8/10._r8*WRMUnit%MeanMthFlow(idam,13)*WRMUnit%MeanMthDemand(idam,mth)/WRMUnit%MeanMthDemand(idam, 13)
                 !TEST
                 !StorWater%pre_release(idam,mth) = WRMUnit%MeanMthDemand(idam,mth)/10._r8 + 9._r8/10._r8*WRMUnit%MeanMthFlow(idam,13)*WRMUnit%MeanMthDemand(idam,mth)/WRMUnit%MeanMthDemand(idam, 13)*.5_r8
              end do
           else 
              do mth=1,12
                 if ( (WRMUnit%MeanMthFlow(idam,13) + WRMUnit%MeanMthDemand(idam,mth) - WRMUnit%MeanMthDemand(idam,13))>0 ) then
                    StorWater%pre_release(idam, mth) = WRMUnit%MeanMthFlow(idam,13) + WRMUnit%MeanMthDemand(idam,mth) - WRMUnit%MeanMthDemand(idam,13)
                 endif 
                 ! test 2
                 !StorWater%pre_release(idam, mth) = WRMUnit%MeanMthFlow(idam,13)*0.5_r8 + WRMUnit%MeanMthDemand(idam,mth) - WRMUnit%MeanMthDemand(idam,13)
                 !TEST use pseudo regulated flow
                 !StorWater%pre_release(idam, mth) = WRMUnit%MeanMthFlow(idam,13)*.5_r8 + WRMUnit%MeanMthDemand(idam,mth) - WRMUnit%MeanMthDemand(idam,13)
              end do
           end if

           !--- initialize storage in each reservoir - arbitrary 50%

           StorWater%storage(idam) = 0.9_r8 * WRMUnit%StorCap(idam)   
           if (WRMUnit%StorageCalibFlag(idam).eq.1) then
              StorWater%storage(idam)  = WRMUnit%StorTarget(idam,13)
           endif
           if ( WRMUnit%StorCap(idam) <= 0 ) then
              write(iulog,*) subname, "Error negative max cap for reservoir ", idam, WRMUnit%StorCap(idam)
              call shr_sys_abort(subname//' ERROR: negative max cap for reservoir')
           end if
           if (idam == 1) then
              write(iulog,*) subname, "storage ",StorWater%pre_release(idam,1), StorWater%storage(idam)
           endif
        end do

        !--- initialize start of the operational year based on long term simulation

        call WRM_init_StOp_FC

! tcraig below is from Nathalie. need to add this check back at some point for performance
!! need to adjust for reservoir with zero inflow, do not  need to read the remaining
!           if ( WRMUnit%MeanMthFlow(idam,13) <= 0._r8 ) then
!              WRMUnit%dam_Ndepend(idam) = 0 ! this reservoir will not provide water to any subw, relieve database
!           end if
!
!           do ng = 1,WRMUnit%dam_Ndepend(idam)
!              call split(stemp,' ',stmp1)
!              call str2num(stmp1, ctlSubwWRM%localNumDam, ierror) 
!!need additional check due to regionalization, need to remove non existing grid cell NV
!              if (  WRMUnit%INVisubw(ctlSubwWRM%localNumDam) .lt. 1 ) then
!                 WRMUnit%dam_Ndepend(idam) = WRMUnit%dam_Ndepend(idam) - 1
!              else
!                 WRMUnit%dam_depend(idam,ng) = nd 
!              endif
!           end do

        !check the dependence database consistencies
        do idam = 1,ctlSubwWRM%localNumDam
           do ng = 1,WRMUnit%dam_Ndepend(idam)
              !if(WRMUnit%dam_depend(idam,ng).eq.0) then
              !   WRMUnit%dam_depend(idam,ng) = WRMUnit%dam_depend(idam,ng) * 1
              !end if
              idepend = WRMUnit%dam_depend(idam,ng)
              if ( idepend <= 0 ) then
                 write(iulog,'(2a,4i12)') subname," Error: checking dependency, zero idepend", idam, WRMUnit%dam_Ndepend(idam), ng, idepend
                 call shr_sys_abort(subname//' ERROR: zero idepend')
              endif
!tcx this should be accumlating global data, not just local data, idepend could be outside begr:endr
!              WRMUnit%TotStorCapDepend(idepend) = WRMUnit%TotStorCapDepend(idepend) + WRMUnit%StorCap(idam)
!              WRMUnit%TotInflowDepend(idepend) = WRMUnit%TotInflowDepend(idepend) + WRMUnit%MeanMthFlow(idam,13)
           end do
        end do

!     end if !Regulation Flag

     call ncd_pio_closefile(ncid)

     ! check
     write(iulog,*) subname, "Done with WM init ..."
     !write(iulog,*) subname, WRMUnit%DamName(1), WRMUnit%Surfarea(1)
     !write(iulog,*) subname,WRMUnit%mask(1), WRMUnit%icell(1) 
     !write(iulog,*) subname, WRMUnit%dam_Ndepend(1), WRMUnit%dam_depend(1,2)
     !write(iulog,*) subname, "sub = 49",  TUnit%icell(49, 1),WRMUnit%subw_Ndepend(49),  WRMUnit%subw_depend(49,1) 
  end subroutine WRM_init

end MODULE WRM_subw_IO_mod
