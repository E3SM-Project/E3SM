!
MODULE WRM_subw_IO_mod
! Description: module to provide interface between WRM and other CLM components
! 
! Developed by Nathalie Voisin 2/1/2010
! REVISION HISTORY:
!-----------------------------------------------------------------------

! !USES:
  use RunoffMod     , only : Tctl, TUnit, rtmCTL
  use RtmSpmd       , only : masterproc, mpicom_rof, iam, ROFID, &
                             MPI_REAL8,MPI_INTEGER,MPI_CHARACTER,MPI_LOGICAL,MPI_MAX
  use RtmVar        , only : iulog, inst_suffix, smat_option
  use RtmFileUtils  , only : getfil, getavu, relavu
  use RtmIO         , only : pio_subsystem, ncd_pio_openfile, ncd_pio_closefile
  use RtmTimeManager, only : get_curr_date
  use rof_cpl_indices, only : nt_rtm
  use shr_kind_mod  , only : r8 => shr_kind_r8, SHR_KIND_CL
  use shr_const_mod , only : SHR_CONST_REARTH, SHR_CONST_PI
  use shr_sys_mod   , only : shr_sys_flush, shr_sys_abort
  use WRM_type_mod  , only : ctlSubwWRM, WRMUnit, StorWater, & 
                             gsMap_wg, gsMap_wd, sMatP_g2d, sMatP_d2g, &
                             aVect_wg, aVect_wd
  use WRM_modules   , only : RegulationRelease, WRM_storage_targets
  use WRM_start_op_year, only : WRM_init_StOp_FC
  use mct_mod
  use netcdf
  use pio
  use perf_mod      , only : t_startf, t_stopf

  implicit none
  private

  public WRM_init
  public WRM_ReadDemand
  public WRM_computeRelease

  type(io_desc_t)  :: iodesc_int_grd2grd ! pio io desc, global grid to local grid
  type(io_desc_t)  :: iodesc_dbl_grd2grd ! pio io desc, global grid to local grid
  type(io_desc_t)  :: iodesc_int_grd2dam ! pio io desc, global grid to local dam
  type(io_desc_t)  :: iodesc_dbl_grd2dam ! pio io desc, global grid to local dam
  type(io_desc_t)  :: iodesc_int_dam2dam ! pio io desc, global dam to local dam
  type(io_desc_t)  :: iodesc_dbl_dam2dam ! pio io desc, global dam to local dam
  character(len=*),parameter :: FORMI = '(2A,6i13)'
  character(len=*),parameter :: FORMR = '(2A,2g15.7)'
  character(len=*),parameter :: FORMR2= '(2A,i8,2g15.7)'

!-----------------------------------------------------------------------
  contains
!-----------------------------------------------------------------------
  
  subroutine WRM_init
     ! !DESCRIPTION: initilization of WRM model
     implicit none

     integer :: nr, nd, ng, mth        ! local loop indices
     integer :: cnt, idam, ntotal, cntw, cntg
     integer :: begr, endr, maxnumdependentgrid, lsize, gsize, lsized, gsized, ssize, dimsize
     integer :: iunit
     integer :: ncdfid, did, dids(2), ndims, dsize(1), dsizes(2), ier, varid
     integer :: start(2), count(2)
     integer(kind=pio_offset) :: frame
     real(r8) :: peak, prorata , mn, mx                  ! peak value to define the start of operationalyr
     integer :: sgn,curr_sgn, nsc, ct, ct_mx, mth_op , idepend      ! number of sign change
     integer :: igrow,igcol,iwgt
     integer :: unitn
     type(file_desc_t):: ncid       ! netcdf file
     type(var_desc_t) :: vardesc    ! netCDF variable description
     integer, pointer  :: compDOF(:)           ! pio decomp
     integer, pointer  :: temp_gridID_from_Dam(:,:)
     integer, allocatable :: gindex(:)         ! global index for gsmap
     integer, allocatable :: tmpmydam(:,:)     ! temporary copy of myDamm for resize
     logical           :: lexist               ! File exists
     type(mct_sMat)    :: sMat                 ! temporary sparse matrix, needed for sMatP
     character(len=256):: nlfilename_wrm

     character(len=350) :: paraFile, demandPath
     integer :: ExtractionFlag, ExtractionMainChannelFlag, RegulationFlag, &
        ReturnFlowFlag, TotalDemandFlag, GroundWaterFlag

     namelist /wrm_inparm/  &
        paraFile, demandPath, &
        ExtractionFlag, ExtractionMainChannelFlag, RegulationFlag, &
        ReturnFlowFlag, TotalDemandFlag, GroundWaterFlag

     character(len=*),parameter :: subname='(WRM_init)'

     begr = rtmCTL%begr
     endr = rtmCTL%endr

     if (masterproc) write(iulog,FORMI) subname,' begr endr = ',iam,begr,endr
     call shr_sys_flush(iulog)

     paraFile = '/pic/scratch/tcraig/wm_data/US_reservoir_8th_NLDAS2.nc'
     demandPath = '/pic/scratch/tcraig/wm_data/NLDAS2_GCAM_water_demand_'
     ExtractionFlag = 1
     ExtractionMainChannelFlag = 1
     RegulationFlag = 1
     ReturnFlowFlag = 0
     TotalDemandFlag = 1
     GroundWaterFlag = 0

     nlfilename_wrm = "mosart_in" // trim(inst_suffix)
     inquire (file = trim(nlfilename_wrm), exist = lexist)
     if ( .not. lexist ) then
        write(iulog,*) subname // ' ERROR: nlfilename_wrm does NOT exist:'&
             //trim(nlfilename_wrm)
        call shr_sys_abort(trim(subname)//' ERROR nlfilename_wrm does not exist')
     end if
     if (masterproc) then
        unitn = getavu()
        write(iulog,*) subname,' Read in wrm_inparm namelist from: ', trim(nlfilename_wrm)
        open( unitn, file=trim(nlfilename_wrm), status='old' )
        ier = 1
        do while ( ier /= 0 )
           read(unitn, wrm_inparm, iostat=ier)
           if (ier < 0) then
              call shr_sys_abort( subname//' encountered end-of-file on wrm_inparm read' )
           endif
        end do
        call relavu( unitn )
     endif

     call mpi_bcast(paraFile   ,len(paraFile)  , MPI_CHARACTER, 0, mpicom_rof, ier)
     call mpi_bcast(demandPath ,len(demandPath), MPI_CHARACTER, 0, mpicom_rof, ier)
     call mpi_bcast(ExtractionFlag,   1, MPI_INTEGER, 0, mpicom_rof, ier)
     call mpi_bcast(ExtractionMainChannelFlag, 1, MPI_INTEGER, 0, mpicom_rof, ier)
     call mpi_bcast(RegulationFlag,   1, MPI_INTEGER, 0, mpicom_rof, ier)
     call mpi_bcast(ReturnFlowFlag,   1, MPI_INTEGER, 0, mpicom_rof, ier)
     call mpi_bcast(TotalDemandFlag,  1, MPI_INTEGER, 0, mpicom_rof, ier)
     call mpi_bcast(GroundWaterFlag,  1, MPI_INTEGER, 0, mpicom_rof, ier)

     ctlSubwWRM%paraFile = paraFile
     ctlSubwWRM%demandPath = demandPath
     ctlSubwWRM%ExtractionFlag = ExtractionFlag
     ctlSubwWRM%ExtractionMainChannelFlag = ExtractionMainChannelFlag
     ctlSubwWRM%RegulationFlag = RegulationFlag
     ctlSubwWRM%ReturnFlowFlag = ReturnFlowFlag
     ctlSubwWRM%TotalDemandFlag = TotalDemandFlag
     ctlSubwWRM%GroundWaterFlag = GroundWaterFlag

     if (masterproc) then
        write(iulog,*) subname," paraFile        = ",trim(ctlSubwWRM%paraFile)
        write(iulog,*) subname," demandPath      = ",trim(ctlSubwWRM%demandPath)
        write(iulog,*) subname," ExtractionFlag  = ",ctlSubwWRM%ExtractionFlag
        write(iulog,*) subname," ExtractionMainChannelFlag = ",ctlSubwWRM%ExtractionMainChannelFlag
        write(iulog,*) subname," RegulationFlag  = ",ctlSubwWRM%RegulationFlag
        write(iulog,*) subname," ReturnFlowFlag  = ",ctlSubwWRM%ReturnFlowFlag
        write(iulog,*) subname," TotalDemandFlag = ",ctlSubwWRM%TotalDemandFlag
        write(iulog,*) subname," GroundWaterFlag = ",ctlSubwWRM%GroundWaterFlag
     endif

     !-------------------
     ! read input dataset
     !-------------------

     ! start by opening and reading the dam IDs, can't use pio yet, need decomp first
     ier = nf90_open(trim(ctlSubwWRM%paraFile), NF90_NOWRITE, ncdfid)
     ier = nf90_inq_dimid(ncdfid,'DependentGrids',did)
     ier = nf90_inquire_dimension(ncdfid,did,len=maxNumDependentGrid)
     ier = nf90_inq_dimid(ncdfid,'Dams',did)
     ier = nf90_inquire_dimension(ncdfid,did,len=ctlSubwWRM%NDam)

     if (masterproc) write(iulog,*) subname,' dams = ',ctlSubwWRM%NDam
     if (masterproc) write(iulog,*) subname,' dependentgrids = ',maxNumDependentGrid
     call shr_sys_flush(iulog)

     !-------------------
     ! now open same file with pio
     !-------------------

     call ncd_pio_openfile(ncid, trim(ctlSubwWRM%paraFile), 0)
     call shr_sys_flush(iulog)
     ier = pio_inq_varid   (ncid, name='DamID_Spatial', vardesc=vardesc)
     ier = pio_inq_vardimid(ncid, vardesc, dids)
     ier = pio_inq_dimlen  (ncid, dids(1),dsizes(1))
     ier = pio_inq_dimlen  (ncid, dids(2),dsizes(2))

     write(iulog,FORMI) subname,' lnumr = ',iam,rtmCTL%lnumr,begr,endr
     write(iulog,FORMI) subname,' gindex = ',iam,minval(rtmCTL%gindex),maxval(rtmCTL%gindex)
     if (masterproc) write(iulog,FORMI) subname,' dsizes = ',dsizes
     call shr_sys_flush(iulog)

     !-------------------
     ! decomp for lon/lat to local gridcells to dams
     !-------------------

     lsize = rtmCTL%lnumr
     allocate(gindex(lsize))
     cnt = 0
     do nr = begr,endr
        cnt = cnt + 1
!tcx, read in data based on mosart decomposition
! THIS ASSUMES mosart and wrm use same grid indexing
        gindex(cnt) = rtmCTL%gindex(nr)
     enddo
     write(iulog,FORMI) subname,' gindex = ',iam,minval(gindex),maxval(gindex)
     call shr_sys_flush(iulog)

     gsize = dsizes(1)*dsizes(2)
     call mct_gsMap_init(gsMap_wg, gindex, mpicom_rof, ROFID, lsize, gsize )
     call mct_aVect_init(aVect_wg, rList='fld1',lsize=lsize)
     write(iulog,FORMI) subname,' gsmap_wg lsize = ',iam,mct_gsMap_lsize(gsMap_wg,mpicom_rof)
     write(iulog,FORMI) subname,' avect_wg lsize = ',iam,mct_aVect_lsize(aVect_wg)
     call shr_sys_flush(iulog)

     call pio_initdecomp(pio_subsystem, pio_double, dsizes, gindex, iodesc_dbl_grd2grd)
     call pio_initdecomp(pio_subsystem, pio_int   , dsizes, gindex, iodesc_int_grd2grd)
     deallocate(gindex)

     !-------------------
     !--- read in DamIDs
     !-------------------

     allocate (WRMUnit%isDam(begr:endr))
     WRMUnit%isDam = 0
     ier = pio_inq_varid (ncid, name='DamID_Spatial', vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_int_grd2grd , WRMUnit%isDam, ier)
     if (masterproc) write(iulog,FORMI) trim(subname),' read DamID_Spatial',minval(WRMUnit%isDam),maxval(WRMUnit%isDam)
     call shr_sys_flush(iulog)

     !-------------------
     !--- count dams and set up dam indexing
     !-------------------

     ctlSubwWRM%localNumDam = 0
     do nr = begr,endr
        if (WRMUnit%isDam(nr) > 0) then
           ctlSubwWRM%localNumDam = ctlSubwWRM%localNumDam + 1
        else
           WRMUnit%isDam(nr) = 0
        end if
     end do
     write(iulog,FORMI) subname,' localNumDam = ',iam,ctlSubwWRM%localNumDam

     !-------------------
     !--- initialize icell, INVicell and dam indexing
     !-------------------

     allocate(WRMUnit%INVicell(begr:endr))
     WRMUnit%INVicell =-99 
     allocate(WRMUnit%icell(ctlsubwWRM%localNumDam))
     WRMUnit%icell = 0
     allocate(WRMUnit%damID(ctlSubwWRM%localNumDam))
     WRMUnit%damID = -99
     allocate(WRMUnit%grdID(ctlSubwWRM%localNumDam))
     WRMUnit%grdID = -99
     cnt = 0
     do nr = begr,endr
        if (WRMUnit%isDam(nr) > 0) then
           cnt = cnt + 1
           WRMUnit%INVicell(nr) = cnt  ! local gridcell index to local dam index
           WRMUnit%icell(cnt) = nr     ! local dam index to local gridcell index
           WRMUnit%damID(cnt) = WRMUnit%isDam(nr)  ! local dam index to global dam index
           WRMUnit%grdID(cnt) = rtmCTL%gindex(nr)  ! local dam index to global grid index
        endif
     enddo

     write(iulog,FORMI) subname,' icell    = ',iam,minval(WRMUnit%icell),maxval(WRMUnit%icell)
     write(iulog,FORMI) subname,' INVicell = ',iam,minval(WRMUnit%INVicell),maxval(WRMUnit%INVicell)
     write(iulog,FORMI) subname,' damID    = ',iam,minval(WRMUnit%damID),maxval(WRMUnit%damID)
     write(iulog,FORMI) subname,' grdID    = ',iam,minval(WRMUnit%grdID),maxval(WRMUnit%grdID)

     !-------------------
     !--- setup dam decomps
     !-------------------

     gsized = ctlSubwWRM%NDam
     lsized = ctlSubwWRM%localNumDam
     allocate(gindex(lsized))
     gindex(1:lsized) = WRMUnit%damID(1:lsized)
     call mct_gsMap_init(gsMap_wd, WRMUnit%damID, mpicom_rof, ROFID, lsized, gsized)
     call mct_aVect_init(aVect_wd, rList='fld1',lsize=lsized)
     write(iulog,FORMI) subname,' gsmap_wd lsize = ',iam,mct_gsMap_lsize(gsMap_wd,mpicom_rof)
     write(iulog,FORMI) subname,' avect_wd lsize = ',iam,mct_aVect_lsize(aVect_wd)
     deallocate(gindex)

     ! decomp for lon/lat to local dam
     call pio_initdecomp(pio_subsystem, pio_double, dsizes, WRMUnit%grdID, iodesc_dbl_grd2dam)
     call pio_initdecomp(pio_subsystem, pio_int   , dsizes, WRMUnit%grdID, iodesc_int_grd2dam)

     ! decomp for global dam to local dam
     dsize(1) = ctlSubwWRM%NDam
     call pio_initdecomp(pio_subsystem, pio_double, dsize , WRMUnit%damID, iodesc_dbl_dam2dam)
     call pio_initdecomp(pio_subsystem, pio_int   , dsize , WRMUnit%damID, iodesc_int_dam2dam)

     !-------------------
     !--- initialize dam/gridcell interactions
     !--- g2d will sum gridcell data to dams
     !--- d2g will sum dam data to gridcells
     !-------------------

     allocate (WRMUnit%dam_Ndepend(ctlSubwWRM%localNumDam))
     WRMUnit%dam_Ndepend = 0
     ier = pio_inq_varid (ncid, name='numGrid_from_Dam', vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_int_grd2dam , WRMUnit%dam_Ndepend, ier)
     write(iulog,FORMI) trim(subname),' read numGrid_from_Dam',iam,minval(WRMUnit%dam_Ndepend),maxval(WRMUnit%dam_Ndepend)
     call shr_sys_flush(iulog)

     !--- read in entire gridID_from_Dam field on all pes
     allocate(temp_gridID_from_Dam(ctlSubwWRM%NDam, maxNumDependentGrid))
     temp_gridID_from_Dam = -99
     ier = nf90_inq_varid(ncdfid,'gridID_from_Dam',varid)
     ier = nf90_get_var(ncdfid,varid,temp_gridID_from_Dam)
     ier = nf90_close(ncdfid)

     allocate (WRMUnit%myDamNum(begr:endr))
     WRMUnit%myDamNum = 0
     allocate (WRMUnit%myDam(10,begr:endr))
     WRMUnit%myDam = 0

     call t_startf('moswrm_init_mydamc')
     !--- total number of depend gridcells, cntg, myDam arrays
     !--- automatic resize of myDam because we don't want to do this search twice
     cntg = 0
     do nd = 1,ctlSubwWRM%NDam
        do ng = 1,maxNumDependentGrid
           if (temp_gridID_from_Dam(nd,ng) > 0) then
              cntg = cntg + 1
              do iunit = begr,endr
                 if (rtmCTL%gindex(iunit) == temp_gridID_from_Dam(nd,ng)) then
                    WRMUnit%myDamNum(iunit) = WRMUnit%myDamNum(iunit) + 1
                    dimsize = size(WRMUnit%myDam,dim=1)
                    if (WRMUnit%myDamNum(iunit) > dimsize) then
                       call t_startf('moswrm_init_mydamc_mem')
                       allocate(tmpmydam(dimsize,begr:endr))
                       tmpmydam(1:dimsize,begr:endr) = WRMUnit%myDam(1:dimsize,begr:endr)
                       deallocate(WRMUnit%myDam)
                       allocate(WRMUnit%myDam(WRMUnit%myDamNum(iunit)+10,begr:endr))
                       WRMUnit%myDam = 0
                       WRMUnit%myDam(1:dimsize,begr:endr) = tmpmydam(1:dimsize,begr:endr)
                       deallocate(tmpmydam)
                       call t_stopf('moswrm_init_mydamc_mem')
                    endif
                    WRMUnit%myDam(WRMUnit%myDamNum(iunit),iunit) = nd
                 endif
              enddo
           endif
        enddo
     enddo
     call t_stopf('moswrm_init_mydamc')

     write(iulog,FORMI) trim(subname),' myDamNumMax = ',iam,maxval(WRMUnit%myDamNum)

     !--- local number of depend gridcells
     cntw = 0
     do idam = 1,ctlSubwWRM%localNumDam
        nd = WRMUnit%damID(idam)
        ntotal = 0
        do ng = 1,maxNumDependentGrid
           if (temp_gridID_from_Dam(nd,ng) > 0) then
              ntotal = ntotal + 1
           endif
        enddo
        ! check that Ndepend and dam_depend are consistent
        ! sometimes in the dependency data, the basin outlet is included for a dam as its dependent grid
        if (ntotal > WRMUnit%dam_Ndepend(idam) .or. ntotal < WRMUnit%dam_Ndepend(idam)-1) then
           write(iulog,"(a,3i8)"), subname//"ERROR: Ndepend/gridID_from_Dam inconsistency1",idam,ntotal,WRMUnit%dam_Ndepend(idam)
           call shr_sys_abort(subname//' ERROR: numGrid_from_Dam not consistent with gridID_from_Dam list')
        end if
        WRMUnit%dam_Ndepend(idam) = ntotal
        cntw = cntw + ntotal
     enddo

     if (smat_option == 'opt') then

        call mct_sMat_init(sMat, gsized, gsize, cntw)  ! sMat, dst, src
        igcol = mct_sMat_indexIA(sMat,'gcol')  ! src
        igrow = mct_sMat_indexIA(sMat,'grow')  ! dst
        iwgt  = mct_sMat_indexRA(sMat,'weight')

        cnt = 0
        do idam = 1,ctlSubwWRM%localNumDam
           nd = WRMUnit%damID(idam)
           do ng = 1,WRMUnit%dam_Ndepend(idam)
              cnt = cnt + 1
              sMat%data%iAttr(igcol,cnt) = temp_gridID_from_Dam(nd,ng)
              sMat%data%iAttr(igrow,cnt) = nd
              sMat%data%rAttr(iwgt ,cnt) = 1.0_r8
           enddo
        enddo
        if (cnt /= cntw) then
           write(iulog,"(a,3i8)"), subname//'ERROR: sMat g2d cnt errora',cntw,cnt
           call shr_sys_abort(subname//' ERROR: sMat g2d cnt errora')
        endif

        call mct_sMatP_Init(sMatP_g2d, sMat, gsMap_wg, gsMap_wd, 0, mpicom_rof, ROFID)
        call mct_sMat_clean(sMat)

        call mct_sMat_init(sMat, gsize, gsized, cntw)
        igcol = mct_sMat_indexIA(sMat,'gcol')  ! src
        igrow = mct_sMat_indexIA(sMat,'grow')
        iwgt  = mct_sMat_indexRA(sMat,'weight')
        cnt = 0
        do idam = 1,ctlSubwWRM%localNumDam
           nd = WRMUnit%damID(idam)
           do ng = 1,WRMUnit%dam_Ndepend(idam)
              cnt = cnt + 1
              sMat%data%iAttr(igcol,cnt) = nd
              sMat%data%iAttr(igrow,cnt) = temp_gridID_from_Dam(nd,ng)
              sMat%data%rAttr(iwgt ,cnt) = 1.0_r8
           enddo
        enddo
        if (cnt /= cntw) then
           write(iulog,"(a,3i8)"), subname//'ERROR: sMat d2g cnt errora',cntw,cnt
           call shr_sys_abort(subname//' ERROR: sMat d2g cnt errora')
        endif

        call mct_sMatP_Init(sMatP_d2g, sMat, gsMap_wd, gsMap_wg, 0, mpicom_rof, ROFID)
        call mct_sMat_clean(sMat)

     elseif (smat_option == 'Xonly' .or. smat_option == 'Yonly') then

        if (masterproc) then
           call mct_sMat_init(sMat, gsized, gsize, cntg)
           igcol = mct_sMat_indexIA(sMat,'gcol')  ! src
           igrow = mct_sMat_indexIA(sMat,'grow')
           iwgt  = mct_sMat_indexRA(sMat,'weight')
           cnt = 0
           do nd = 1,ctlSubwWRM%NDam
              do ng = 1,maxNumDependentGrid
                 if (temp_gridID_from_Dam(nd,ng) > 0) then
                    cnt = cnt + 1
                    sMat%data%iAttr(igcol,cnt) = temp_gridID_from_Dam(nd,ng)
                    sMat%data%iAttr(igrow,cnt) = nd
                    sMat%data%rAttr(iwgt ,cnt) = 1.0_r8
                 endif
              enddo
           enddo
           if (cnt /= cntg) then
              write(iulog,"(a,3i8)"), subname//'ERROR: sMat g2d cnt errorb',cntg,cnt
              call shr_sys_abort(subname//' ERROR: sMat g2d cnt errorb')
           endif
        else
          call mct_sMat_init(sMat,1,1,1)
        endif
        call mct_sMatP_Init(sMatP_g2d, sMat, gsMap_wg, gsMap_wd, smat_option, 0, mpicom_rof, ROFID)
        call mct_sMat_clean(sMat)

        if (masterproc) then
           call mct_sMat_init(sMat, gsize, gsized, cntg)
           igcol = mct_sMat_indexIA(sMat,'gcol')  ! src
           igrow = mct_sMat_indexIA(sMat,'grow')
           iwgt  = mct_sMat_indexRA(sMat,'weight')
           cnt = 0
           do nd = 1,ctlSubwWRM%NDam
              do ng = 1,maxNumDependentGrid
                 if (temp_gridID_from_Dam(nd,ng) > 0) then
                    cnt = cnt + 1
                    sMat%data%iAttr(igcol,cnt) = nd
                    sMat%data%iAttr(igrow,cnt) = temp_gridID_from_Dam(nd,ng)
                    sMat%data%rAttr(iwgt ,cnt) = 1.0_r8
                 endif
              enddo
           enddo
           if (cnt /= cntg) then
              write(iulog,"(a,3i8)"), subname//'ERROR: sMat d2g cnt errorb',cntg,cnt
              call shr_sys_abort(subname//' ERROR: sMat d2g cnt errorb')
           endif
        else
          call mct_sMat_init(sMat,1,1,1)
        endif
        call mct_sMatP_Init(sMatP_d2g, sMat, gsMap_wd, gsMap_wg, smat_option, 0, mpicom_rof, ROFID)
        call mct_sMat_clean(sMat)

     else

        write(iulog,*) trim(subname),' WRM ERROR: invalid smat_option '//trim(smat_option)
        call shr_sys_abort(trim(subname)//' ERROR invald smat option')

     endif

     deallocate(temp_gridID_from_Dam)
     ssize = mct_smat_gNumEl(sMatP_d2g%Matrix,mpicom_rof)
     if (masterproc) write(iulog,*) subname," Done initializing SmatP_d2g, nElements = ",ssize
     ssize = mct_smat_gNumEl(sMatP_g2d%Matrix,mpicom_rof)
     if (masterproc) write(iulog,*) subname," Done initializing SmatP_g2d, nElements = ",ssize

     !-------------------
     !--- initial dam data
     !-------------------

     allocate (WRMUnit%DamName(ctlSubwWRM%localNumDam))
     WRMUnit%DamName = 'undefined'
     allocate (WRMUnit%TotStorCapDepend(begr:endr))
     WRMUnit%TotStorCapDepend = 0._r8
     allocate (WRMUnit%TotInflowDepend(begr:endr))
     WRMUnit%TotInflowDepend = 0._r8

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
     allocate (StorWater%demand0(begr:endr))
     StorWater%demand0=0._r8
     allocate (StorWater%supply(begr:endr))
     StorWater%supply=0._r8
     allocate (StorWater%deficit(begr:endr))
     StorWater%deficit=0._r8
     allocate (StorWater%storageG(begr:endr))
     StorWater%storageG=0._r8
     allocate (StorWater%releaseG(begr:endr))
     StorWater%releaseG=0._r8
     allocate (WRMUnit%StorMthStOpG(begr:endr))
     WRMUnit%StorMthStOpG = 0

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

     call WRM_readDemand()  ! initialize demand0

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
     ! NV uncommented out - needed for first year of simulation
     WRMUnit%StorMthStOp = WRMUnit%StorCap*0.85_r8

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
!NV : to be consistent with the flags, you can allow for extraction but not
!regulation, flag should stay although very  likely not used
     if ( ctlSubwWRM%RegulationFlag > 0 ) then

        !--- read mean monthly flow data
        do mth = 1,12
           ier = pio_inq_varid (ncid, name='Qmon', vardesc=vardesc)
           frame = mth
           call pio_setframe(vardesc,frame)
           call pio_read_darray(ncid, vardesc, iodesc_dbl_dam2dam, WRMUnit%MeanMthFlow(:,mth), ier)
           if (masterproc) write(iulog,FORMR) trim(subname),' read Qmon',minval(WRMUnit%MeanMthFlow(:,mth)),maxval(WRMUnit%MeanMthFlow(:,mth))
           call shr_sys_flush(iulog)
        enddo
!NV not okay at all but needed for checking as too much flow comng in with
!respect to set up parameters - need to be removed after test
!        WRMUnit%MeanMthFlow = WRMUnit%MeanMthFlow * 1.6_r8
! end not okay

        do idam = 1,ctlSubwWRM%LocalNumDam
           WRMUnit%MeanMthFlow(idam,13) = sum(WRMUnit%MeanMthFlow(idam,1:12))/12.0_r8
        enddo
        if (masterproc) write(iulog,FORMR) trim(subname),' MeanMthFlow avg',minval(WRMUnit%MeanMthFlow(:,13)),maxval(WRMUnit%MeanMthFlow(:,13))

        !--- read in mean monthly water demand
        do mth = 1,12
           ier = pio_inq_varid (ncid, name='demand', vardesc=vardesc)
           frame = mth
           call pio_setframe(vardesc,frame)
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

           !--- initialize storage in each reservoir - arbitrary 90%

           StorWater%storage(idam) = 0.9_r8 * WRMUnit%StorCap(idam)   
           if (WRMUnit%StorageCalibFlag(idam).eq.1) then
              StorWater%storage(idam)  = WRMUnit%StorTarget(idam,13)
           endif
           if ( WRMUnit%StorCap(idam) <= 0 ) then
              write(iulog,*) subname, "Error negative max cap for reservoir ", idam, WRMUnit%StorCap(idam)
              call shr_sys_abort(subname//' ERROR: negative max cap for reservoir')
           end if
           if (idam == 80) then
              write(iulog,*) subname, "storage ",StorWater%pre_release(idam,1), StorWater%storage(idam)
           endif
        end do

!NV
           if (masterproc) write(iulog,FORMR) trim(subname),'prerelease Jan',minval(StorWater%pre_release(:,1)),maxval(StorWater%pre_release(:,1)) 

           if (masterproc) write(iulog,FORMR) trim(subname),'prerelease Apr', minval(StorWater%pre_release(:,4)),maxval(StorWater%pre_release(:,4))
           if (masterproc) write(iulog,FORMR) trim(subname),'prerelease Jul',minval(StorWater%pre_release(:,7)),maxval(StorWater%pre_release(:,7))
           if (masterproc) write(iulog,FORMR) trim(subname),'prerelease Oct', minval(StorWater%pre_release(:,10)),maxval(StorWater%pre_release(:,10))
           if (masterproc) write(iulog,FORMR) trim(subname),'Coulee',StorWater%pre_release(80,1),StorWater%pre_release(80,3)
           if (masterproc) write(iulog,FORMR) trim(subname),'Coulee',StorWater%pre_release(80,5),StorWater%pre_release(80,9)
           if (masterproc) write(iulog,FORMR) trim(subname),'Coulee Flow',WRMUnit%MeanMthFlow(80,13), WRMUnit%MeanMthFlow(80,8)
           if (masterproc) write(iulog,FORMR) trim(subname),'Coulee Demand',WRMUnit%MeanMthFlow(80,1),WRMUnit%MeanMthFlow(80,3)
           if (masterproc) write(iulog,FORMR) trim(subname),'Coulee Demand',WRMUnit%MeanMthFlow(80,5),WRMUnit%MeanMthFlow(80,7)
           if (masterproc) write(iulog,FORMR) trim(subname),'Coulee Means Flow Demand',WRMUnit%MeanMthDemand(80,13),WRMUnit%MeanMthFlow(80,13)
           call shr_sys_flush(iulog)

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
!        do idam = 1,ctlSubwWRM%localNumDam
!           do ng = 1,WRMUnit%dam_Ndepend(idam)
!              !if (WRMUnit%dam_depend(idam,ng).eq.0) then
!              !   WRMUnit%dam_depend(idam,ng) = WRMUnit%dam_depend(idam,ng) * 1
!              !end if
!              idepend = WRMUnit%dam_depend(idam,ng)
!              if ( idepend <= 0 ) then
!                 write(iulog,'(2a,4i12)') subname," Error: checking dependency, zero idepend", idam, WRMUnit%dam_Ndepend(idam), ng, idepend
!                 call shr_sys_abort(subname//' ERROR: zero idepend')
!              endif
!!tcx this should be accumlating global data, not just local data, idepend could be outside begr:endr
!!              WRMUnit%TotStorCapDepend(idepend) = WRMUnit%TotStorCapDepend(idepend) + WRMUnit%StorCap(idam)
!!              WRMUnit%TotInflowDepend(idepend) = WRMUnit%TotInflowDepend(idepend) + WRMUnit%MeanMthFlow(idam,13)
!           end do
!        end do

!recommented out after discussion with tcx
! NV here I extracted piece of code without the databse check but needed for the
! remaining of the code
        !do idam = 1,ctlSubwWRM%localNumDam
        !   do ng = 1,WRMUnit%dam_Ndepend(idam)
        !      WRMUnit%TotStorCapDepend(idepend) =
!WRMUnit%TotStorCapDepend(idepend) + WRMUnit%StorCap(idam)
!              WRMUnit%TotInflowDepend(idepend) =
!WRMUnit%TotInflowDepend(idepend) + WRMUnit%MeanMthFlow(idam,13)
!           end do
!        end do

     end if !Regulation Flag

     call ncd_pio_closefile(ncid)

     ! check
     write(iulog,*) subname, "Done with WM init ..."
     !write(iulog,*) subname, WRMUnit%DamName(59), WRMUnit%Surfarea(59)
     !write(iulog,*) subname,WRMUnit%isDam(1), WRMUnit%icell(1) 
     !write(iulog,*) subname, WRMUnit%dam_Ndepend(1), WRMUnit%dam_depend(1,2)
     !write(iulog,*) subname, "sub = 49",  TUnit%icell(49, 1),WRMUnit%subw_Ndepend(49),  WRMUnit%subw_depend(49,1) 
  end subroutine WRM_init

!-----------------------------------------------------------------------

  subroutine WRM_readDemand
  ! !DESCRIPTION: read in the irrigation demand data for each time step
     implicit none
     character(len=250) :: demFileName  ! water demand file names
     integer :: ios, iunit, ilat, ilon    ! flag of IO status, local indices
     real(r8) :: ftemp1            ! tempory array

     integer :: yr, mon, day, tod
     character(len=4) :: strYear
     character(len=2) :: strMonth, strDay
     character(len=1000) :: fname
     integer  :: nr, ier
     type(file_desc_t):: ncid       ! netcdf file
     type(var_desc_t) :: vardesc    ! netCDF variable description
     character(len=*),parameter :: subname='(WRM_readDemand)'

     call get_curr_date(yr, mon, day, tod)
     write(iulog,'(2a,4i6)') subname,'at ',yr,mon,day,tod
    
     write(strYear,'(I4.4)') yr
     write(strMonth,'(I2.2)') mon
     fname = trim(ctlSubwWRM%demandPath)// strYear//'_'//strMonth//'.nc'

     write(iulog,*) subname, ' reading ',trim(fname)

     call ncd_pio_openfile(ncid, trim(fname), 0)
     ier = pio_inq_varid (ncid, name='totalDemand', vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_dbl_grd2grd , StorWater%demand0, ier)
     call ncd_pio_closefile(ncid)

     do nr=rtmCTL%begr,rtmCTL%endr
        if (StorWater%demand0(nr).lt.0._r8) then
           StorWater%demand0(nr) = 0._r8
        end if
     end do
     write(iulog,FORMR2) trim(subname),' read totalDemand',iam,minval(StorWater%demand0),maxval(StorWater%demand0)
     call shr_sys_flush(iulog)

  end subroutine WRM_readDemand

!-----------------------------------------------------------------------

  subroutine WRM_computeRelease()
     implicit none
     integer :: yr, mon, day, tod
     integer :: idam
     character(len=*),parameter :: subname = '(WRM_computeRelease)'

     call get_curr_date(yr, mon, day, tod)
     write(iulog,'(2a,4i6)') subname,'at ',yr,mon,day,tod
     do idam=1,ctlSubwWRM%localNumDam
        if ( mon .eq. WRMUnit%MthStOp(idam)) then
           WRMUnit%StorMthStOp(idam) = StorWater%storage(idam)
	   end if
     enddo
     call RegulationRelease()
     write(iulog,*) 'Start Coulee ',mon,day,tod,WRMUnit%MeanMthFlow(80,13)
     write(iulog,*) 'start Op mon, storage ', WRMUnit%MthStOp(80),WRMUnit%StorMthStOp(80)
     write(iulog,*)  'storage, release pre targets ',StorWater%storage(80), StorWater%release(80)
     call WRM_storage_targets()
     write(iulog,*) 'Coulee targets ',StorWater%release(80)

  end subroutine WRM_computeRelease

!-----------------------------------------------------------------------

end MODULE WRM_subw_IO_mod
