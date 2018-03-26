module decompInitMod

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module provides a descomposition into a clumped data structure which can
  ! be mapped back to atmosphere physics chunks.
  !
  ! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_sys_mod     , only : shr_sys_flush
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use spmdMod         , only : masterproc, iam, npes, mpicom, comp_id
  use abortutils      , only : endrun
  use clm_varctl      , only : iulog, use_fates
  use clm_varcon      , only : grlnd
  use GridcellType    , only : grc_pp
  use LandunitType    , only : lun_pp                
  use TopounitType    , only : top_pp                
  use ColumnType      , only : col_pp                
  use FatesInterfaceMod, only : fates_maxElementsPerSite
  use VegetationType  , only : veg_pp                
  use decompMod
  use mct_mod
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public decompInit_lnd          ! initializes lnd grid decomposition into clumps and processors
  public decompInit_clumps       ! initializes atm grid decomposition into clumps
  public decompInit_gtlcp         ! initializes g,l,c,p decomp info
  public decompInit_lnd_using_gp ! initialize lnd grid decomposition into clumps and processors using graph partitioning approach
  public decompInit_ghosts       ! initialize ghost/halo for land grid
  !
  ! !PRIVATE TYPES:
  private
  integer, pointer :: lcid(:)       ! temporary for setting ldecomp
  !------------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------------
  subroutine decompInit_lnd(lni,lnj,amask)
    !
    ! !DESCRIPTION:
    ! This subroutine initializes the land surface decomposition into a clump
    ! data structure.  This assumes each pe has the same number of clumps
    ! set by clump_pproc
    !
    ! !USES:
    use clm_varctl, only : nsegspc
    !
    ! !ARGUMENTS:
    implicit none
    integer , intent(in) :: amask(:)
    integer , intent(in) :: lni,lnj   ! domain global size
    !
    ! !LOCAL VARIABLES:
    integer :: lns                    ! global domain size
    integer :: ln,lj                  ! indices
    integer :: ag,an,ai,aj            ! indices
    integer :: numg                   ! number of land gridcells
    logical :: seglen1                ! is segment length one
    real(r8):: seglen                 ! average segment length
    real(r8):: rcid                   ! real value of cid
    integer :: cid,pid                ! indices
    integer :: n,m,ng                 ! indices
    integer :: ier                    ! error code
    integer :: beg,end,lsize,gsize    ! used for gsmap init
    integer, pointer :: gindex(:)     ! global index for gsmap init
    integer, pointer :: clumpcnt(:)   ! clump index counter
    integer, allocatable :: proc_ncell(:) ! number of cells assigned to a process
    integer, allocatable :: proc_begg(:)  ! beginning cell index assigned to a process
    !------------------------------------------------------------------------------

    lns = lni * lnj

    !--- set and verify nclumps ---
    if (clump_pproc > 0) then
       nclumps = clump_pproc * npes
       if (nclumps < npes) then
          write(iulog,*) 'decompInit_lnd(): Number of gridcell clumps= ',nclumps, &
               ' is less than the number of processes = ', npes
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if
    else
       write(iulog,*)'clump_pproc= ',clump_pproc,'  must be greater than 0'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    ! allocate and initialize procinfo (from decompMod.F90) and clumps 
    ! beg and end indices initialized for simple addition of cells later 

    allocate(procinfo%cid(clump_pproc), stat=ier)
    if (ier /= 0) then
       write(iulog,*) 'decompInit_lnd(): allocation error for procinfo%cid'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif
    procinfo%nclumps = clump_pproc
    procinfo%cid(:)  = -1
    procinfo%ncells  = 0
    procinfo%ntopounits  = 0
    procinfo%nlunits = 0
    procinfo%ncols   = 0
    procinfo%npfts   = 0
    procinfo%nCohorts = 0
    procinfo%begg    = 1
    procinfo%begt    = 1
    procinfo%begl    = 1
    procinfo%begc    = 1
    procinfo%begp    = 1
    procinfo%begCohort    = 1
    procinfo%endg    = 0
    procinfo%endt    = 0
    procinfo%endl    = 0
    procinfo%endc    = 0
    procinfo%endp    = 0
    procinfo%endCohort    = 0

    allocate(clumps(nclumps), stat=ier)
    if (ier /= 0) then
       write(iulog,*) 'decompInit_lnd(): allocation error for clumps'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if
    clumps(:)%owner   = -1
    clumps(:)%ncells  = 0
    clumps(:)%ntopounits = 0
    clumps(:)%nlunits = 0
    clumps(:)%ncols   = 0
    clumps(:)%npfts   = 0
    clumps(:)%nCohorts = 0
    clumps(:)%begg    = 1
    clumps(:)%begt    = 1
    clumps(:)%begl    = 1
    clumps(:)%begc    = 1
    clumps(:)%begp    = 1
    clumps(:)%begCohort    = 1
    clumps(:)%endg    = 0
    clumps(:)%endt    = 0
    clumps(:)%endl    = 0
    clumps(:)%endc    = 0
    clumps(:)%endp    = 0
    clumps(:)%endCohort    = 0

    ! assign clumps to proc round robin 
    cid = 0
    do n = 1,nclumps
       pid = mod(n-1,npes)
       if (pid < 0 .or. pid > npes-1) then
          write(iulog,*) 'decompInit_lnd(): round robin pid error ',n,pid,npes
          call endrun(msg=errMsg(__FILE__, __LINE__))
       endif
       clumps(n)%owner = pid
       if (iam == pid) then
          cid = cid + 1
          if (cid < 1 .or. cid > clump_pproc) then
             write(iulog,*) 'decompInit_lnd(): round robin pid error ',n,pid,npes
             call endrun(msg=errMsg(__FILE__, __LINE__))
          endif
          procinfo%cid(cid) = n
       endif
    enddo

    ! count total land gridcells
    numg = 0
    do ln = 1,lns
       if (amask(ln) == 1) then
          numg = numg + 1
       endif
    enddo

    if (npes > numg) then
       write(iulog,*) 'decompInit_lnd(): Number of processes exceeds number ', &
            'of land grid cells',npes,numg
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if
    if (nclumps > numg) then
       write(iulog,*) 'decompInit_lnd(): Number of clumps exceeds number ', &
            'of land grid cells',nclumps,numg
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    if (float(numg)/float(nclumps) < float(nsegspc)) then
       seglen1 = .true.
       seglen = 1.0_r8
    else
       seglen1 = .false.
       seglen = dble(numg)/(dble(nsegspc)*dble(nclumps))
    endif

    if (masterproc) then
       write(iulog,*) ' decomp precompute numg,nclumps,seglen1,avg_seglen,nsegspc=', &
            numg,nclumps,seglen1,&
            sngl(seglen),sngl(dble(numg)/(seglen*dble(nclumps)))
    end if

    ! Assign gridcells to clumps (and thus pes) ---

    allocate(lcid(lns), stat=ier)
    if (ier /= 0) then
       write(iulog,*) 'decompInit_lnd(): allocation error for lcid'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    lcid(:) = 0
    ng = 0
    do ln = 1,lns
       if (amask(ln) == 1) then
          ng = ng  + 1

          !--- give to clumps in order based on nsegspc
          if (seglen1) then
             cid = mod(ng-1,nclumps) + 1
          else
             rcid = (dble(ng-1)/dble(numg))*dble(nsegspc)*dble(nclumps)
             cid = mod(int(rcid),nclumps) + 1
          endif
          lcid(ln) = cid

          !--- give gridcell cell to pe that owns cid ---
          !--- this needs to be done to subsequently use function
          !--- get_proc_bounds(begg,endg) 
          if (iam == clumps(cid)%owner) then
             procinfo%ncells  = procinfo%ncells  + 1
          endif

          !--- give gridcell to cid ---
          clumps(cid)%ncells  = clumps(cid)%ncells  + 1

       end if
    enddo

    ! calculate number of cells per process
    allocate(proc_ncell(0:npes-1), stat=ier)
    if (ier /= 0) then
       write(iulog,*) 'decompInit_lnd(): allocation error for proc_ncell'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    proc_ncell(:) = 0
    do cid = 1,nclumps
       proc_ncell(clumps(cid)%owner) = proc_ncell(clumps(cid)%owner) + clumps(cid)%ncells
    enddo

    ! determine offset (begg) for all processes,
    ! and then procinfo%begg and procinfo%endg (for iam)
    allocate(proc_begg(0:npes-1), stat=ier)
    if (ier /= 0) then
       write(iulog,*) 'decompInit_lnd(): allocation error for proc_begg'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    proc_begg(0) = 1
    do pid = 1,npes-1
       proc_begg(pid) = proc_begg(pid-1) + proc_ncell(pid-1)
    enddo
    procinfo%begg = proc_begg(iam)
    procinfo%endg = (procinfo%begg-1) + procinfo%ncells

    ! determine offset for each clump assigned to each process
    ! (re-using proc_begg as work space)
    do cid = 1,nclumps
      clumps(cid)%begg = proc_begg(clumps(cid)%owner)
      proc_begg(clumps(cid)%owner) = proc_begg(clumps(cid)%owner) &
                                   + clumps(cid)%ncells
      clumps(cid)%endg = proc_begg(clumps(cid)%owner) - 1
    enddo

    ! free work space
    deallocate(proc_ncell, proc_begg)

    ! Set ldecomp

    allocate(ldecomp%gdc2glo(numg), stat=ier)
    if (ier /= 0) then
       write(iulog,*) 'decompInit_lnd(): allocation error1 for ldecomp, etc'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if
    allocate(clumpcnt(nclumps),stat=ier)
    if (ier /= 0) then
       write(iulog,*) 'decompInit_lnd(): allocation error1 for clumpcnt'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    ldecomp%gdc2glo(:) = 0
    ag = 0

    ! clumpcnt is the start gdc index of each clump

    clumpcnt = 0
    do cid = 1,nclumps
       clumpcnt(cid) = clumps(cid)%begg
    enddo

    ! now go through gridcells one at a time and increment clumpcnt
    ! in order to set gdc2glo

    do aj = 1,lnj
    do ai = 1,lni
       an = (aj-1)*lni + ai
       cid = lcid(an)
       if (cid > 0) then
          ag = clumpcnt(cid)
          ldecomp%gdc2glo(ag) = an
          clumpcnt(cid) = clumpcnt(cid) + 1
       end if
    end do
    end do

    deallocate(clumpcnt)

    ! Set gsMap_lnd_gdc2glo (the global index here includes mask=0 or ocean points)

    call get_proc_bounds(beg, end)
    allocate(gindex(beg:end))
    do n = beg,end
       gindex(n) = ldecomp%gdc2glo(n)
    enddo
    lsize = end-beg+1
    gsize = lni * lnj
    call mct_gsMap_init(gsMap_lnd_gdc2glo, gindex, mpicom, comp_id, lsize, gsize)
    deallocate(gindex)

    ! Diagnostic output

    if (masterproc) then
       write(iulog,*)' Surface Grid Characteristics'
       write(iulog,*)'   longitude points               = ',lni
       write(iulog,*)'   latitude points                = ',lnj
       write(iulog,*)'   total number of land gridcells = ',numg
       write(iulog,*)' Decomposition Characteristics'
       write(iulog,*)'   clumps per process             = ',clump_pproc
       write(iulog,*)' gsMap Characteristics'
       write(iulog,*) '  lnd gsmap glo num of segs      = ',mct_gsMap_ngseg(gsMap_lnd_gdc2glo)
       write(iulog,*)
    end if

    call shr_sys_flush(iulog)

  end subroutine decompInit_lnd

  !------------------------------------------------------------------------------
  subroutine decompInit_clumps(glcmask)
    !
    ! !DESCRIPTION:
    ! This subroutine initializes the land surface decomposition into a clump
    ! data structure.  This assumes each pe has the same number of clumps
    ! set by clump_pproc
    !
    ! !USES:
    use subgridMod, only : subgrid_get_gcellinfo
    use spmdMod
    !
    ! !ARGUMENTS:
    implicit none
    integer , pointer, optional   :: glcmask(:)  ! glc mask
    !
    ! !LOCAL VARIABLES:
    integer :: ln,an              ! indices
    integer :: i,g,l,k            ! indices
    integer :: cid,pid            ! indices
    integer :: n,m,np             ! indices
    integer :: anumg              ! lnd num gridcells
    integer :: icells             ! temporary
    integer :: begg, endg         ! temporary
    integer :: itopounits         ! temporary
    integer :: ilunits            ! temporary
    integer :: icols              ! temporary
    integer :: ipfts              ! temporary
    integer :: icohorts           ! temporary
    integer :: ier                ! error code
    integer :: glev, tlev, llev, clev, plev, hlev  ! order of subgrid levels in the allvec arrays
    integer :: nlev               ! number of subgrid levels
    integer, allocatable :: allvecg(:,:)  ! temporary vector "global"
    integer, allocatable :: allvecl(:,:)  ! temporary vector "local"
    integer, allocatable :: proc_nXXX(:) ! number of XXX assigned to a process
    integer, allocatable :: proc_begX(:) ! beginning XXX index assigned to a process
    integer :: ntest
    character(len=32), parameter :: subname = 'decompInit_clumps'
    !------------------------------------------------------------------------------
    
    !--- assign order of subgrid levels in allvecl and allvecg arrays ---
    nlev=6  ! number of subgrid levels
    glev=1  ! gridcell
    tlev=2  ! topounit
    llev=3  ! landunit
    clev=4  ! column
    plev=5  ! pft/patch
    hlev=6  ! cohort

    !--- assign gridcells to clumps (and thus pes) ---
    call get_proc_bounds(begg, endg)

    allocate(allvecl(nclumps,nlev))   ! local  clumps [gcells,topounits,lunits,cols,pfts,cohs]
    allocate(allvecg(nclumps,nlev))   ! global clumps [gcells,topounits,lunits,cols,pfts,cohs]

    ! Determine the number of gridcells, topounits, landunits, columns, pfts, and cohorts 
    ! on this processor 
    ! Determine number of topounits, landunits, columns and pfts for each global
    ! gridcell index (an) that is associated with the local gridcell index (ln)
    ! More detail: an is the row-major order 1d-index into the global ixj grid.

    itopounits=0
    ilunits=0
    icols=0
    ipfts=0
    icohorts=0 

    allvecg= 0
    allvecl= 0
    ! Loop through the gridcells on this proc
    do anumg = begg,endg
       ! an is the row-major order 1d-index into the global ixj grid.
       an  = ldecomp%gdc2glo(anumg)
       cid = lcid(an)
       ln  = anumg
       if (present(glcmask)) then
          call subgrid_get_gcellinfo (ln, ntopounits=itopounits, nlunits=ilunits, ncols=icols, npfts=ipfts, &
              ncohorts=icohorts, glcmask=glcmask(ln))
       else
          call subgrid_get_gcellinfo (ln, ntopounits=itopounits, nlunits=ilunits, ncols=icols, npfts=ipfts, &
               ncohorts=icohorts )
       endif
       allvecl(cid,glev) = allvecl(cid,glev) + 1           ! number of gridcells for local clump cid
       allvecl(cid,tlev) = allvecl(cid,tlev) + itopounits  ! number of topographic units for local clump cid
       allvecl(cid,llev) = allvecl(cid,llev) + ilunits     ! number of landunits for local clump cid
       allvecl(cid,clev) = allvecl(cid,clev) + icols       ! number of columns for local clump cid
       allvecl(cid,plev) = allvecl(cid,plev) + ipfts       ! number of pfts for local clump cid 
       allvecl(cid,hlev) = allvecl(cid,hlev) + icohorts    ! number of cohorts for local clump cid 
    enddo
    call mpi_allreduce(allvecl,allvecg,size(allvecg),MPI_INTEGER,MPI_SUM,mpicom,ier)

    ! Determine overall  total gridcells, landunits, columns and pfts and distribute
    ! gridcells over clumps

    numg = 0
    numt = 0
    numl = 0
    numc = 0
    nump = 0
    numCohort = 0

    do cid = 1,nclumps
       icells      = allvecg(cid,glev)  ! number of all clump cid gridcells (over all processors)
       itopounits  = allvecg(cid,tlev)  ! number of all clump cid topounits (over all processors)
       ilunits     = allvecg(cid,llev)  ! number of all clump cid landunits (over all processors)
       icols       = allvecg(cid,clev)  ! number of all clump cid columns (over all processors)
       ipfts       = allvecg(cid,plev)  ! number of all clump cid pfts (over all processors)
       icohorts    = allvecg(cid,hlev)  ! number of all clump cid cohorts (over all processors)

       !--- overall total ---
       numg = numg + icells         ! total number of gridcells
       numt = numt + itopounits     ! total number of landunits
       numl = numl + ilunits        ! total number of landunits
       numc = numc + icols          ! total number of columns
       nump = nump + ipfts          ! total number of pfts
       numCohort = numCohort + icohorts       ! total number of cohorts

       !--- give gridcell to cid ---
       clumps(cid)%ntopounits  = clumps(cid)%ntopounits  + itopounits  
       clumps(cid)%nlunits     = clumps(cid)%nlunits  + ilunits  
       clumps(cid)%ncols       = clumps(cid)%ncols    + icols
       clumps(cid)%npfts       = clumps(cid)%npfts    + ipfts
       clumps(cid)%nCohorts    = clumps(cid)%nCohorts + icohorts

       !--- give gridcell to the proc that owns the cid ---
       if (iam == clumps(cid)%owner) then
          procinfo%ntopounits  = procinfo%ntopounits  + itopounits
          procinfo%nlunits     = procinfo%nlunits  + ilunits
          procinfo%ncols       = procinfo%ncols    + icols
          procinfo%npfts       = procinfo%npfts    + ipfts
          procinfo%nCohorts    = procinfo%nCohorts + icohorts
       endif

    enddo

    ! determine offset for XXX (topounits/lunits/cols/pfts/cohorts) index
    ! for each process

    allocate(proc_nXXX(0:npes-1), stat=ier)
    if (ier /= 0) then
       write(iulog,*) 'decompInit_clumps(): allocation error for proc_nXXX'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    allocate(proc_begX(0:npes-1), stat=ier)
    if (ier /= 0) then
       write(iulog,*) 'decompInit_clumps(): allocation error for proc_begX'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    ! TOPOUNITS:
    ! calculate number of topographic units per process
    proc_nXXX(:) = 0
    do cid = 1,nclumps
       proc_nXXX(clumps(cid)%owner) = &
        proc_nXXX(clumps(cid)%owner) + clumps(cid)%ntopounits
    enddo

    ! determine offset (begt) for all processes,
    ! and then procinfo%begt and procinfo%endt (for iam)
    proc_begX(0) = 1
    do pid = 1,npes-1
       proc_begX(pid) = proc_begX(pid-1) + proc_nXXX(pid-1)
    enddo
    procinfo%begt = proc_begX(iam)
    procinfo%endt = (procinfo%begt-1) + procinfo%ntopounits

    ! determine topounit offset for each clump assigned to each process
    ! (re-using proc_begX as work space)
    do cid = 1,nclumps
      clumps(cid)%begt = proc_begX(clumps(cid)%owner)
      proc_begX(clumps(cid)%owner) = proc_begX(clumps(cid)%owner) &
                                   + clumps(cid)%ntopounits
      clumps(cid)%endt = proc_begX(clumps(cid)%owner) - 1
    enddo

    ! LUNITS:
    ! calculate number of lunits per process
    proc_nXXX(:) = 0
    do cid = 1,nclumps
       proc_nXXX(clumps(cid)%owner) = &
        proc_nXXX(clumps(cid)%owner) + clumps(cid)%nlunits
    enddo

    ! determine offset (begl) for all processes,
    ! and then procinfo%begl and procinfo%endl (for iam)
    proc_begX(0) = 1
    do pid = 1,npes-1
       proc_begX(pid) = proc_begX(pid-1) + proc_nXXX(pid-1)
    enddo
    procinfo%begl = proc_begX(iam)
    procinfo%endl = (procinfo%begl-1) + procinfo%nlunits

    ! determine lunit offset for each clump assigned to each process
    ! (re-using proc_begX as work space)
    do cid = 1,nclumps
      clumps(cid)%begl = proc_begX(clumps(cid)%owner)
      proc_begX(clumps(cid)%owner) = proc_begX(clumps(cid)%owner) &
                                   + clumps(cid)%nlunits
      clumps(cid)%endl = proc_begX(clumps(cid)%owner) - 1
    enddo

    ! COLS:
    ! calculate number of cols per process
    proc_nXXX(:) = 0
    do cid = 1,nclumps
       proc_nXXX(clumps(cid)%owner) = &
        proc_nXXX(clumps(cid)%owner) + clumps(cid)%ncols
    enddo

    ! determine offset (begc) for all processes,
    ! and then procinfo%begc and procinfo%endc (for iam)
    proc_begX(0) = 1
    do pid = 1,npes-1
       proc_begX(pid) = proc_begX(pid-1) + proc_nXXX(pid-1)
    enddo
    procinfo%begc = proc_begX(iam)
    procinfo%endc = (procinfo%begc-1) + procinfo%ncols

    ! determine col offset for each clump assigned to each process
    ! (re-using proc_begX as work space)
    do cid = 1,nclumps
      clumps(cid)%begc = proc_begX(clumps(cid)%owner)
      proc_begX(clumps(cid)%owner) = proc_begX(clumps(cid)%owner) &
                                   + clumps(cid)%ncols
      clumps(cid)%endc = proc_begX(clumps(cid)%owner) - 1
    enddo

    ! PFTS:
    ! calculate number of pfts per process
    proc_nXXX(:) = 0
    do cid = 1,nclumps
       proc_nXXX(clumps(cid)%owner) = &
        proc_nXXX(clumps(cid)%owner) + clumps(cid)%npfts
    enddo

    ! determine offset (begp) for all processes,
    ! and then procinfo%begp and procinfo%endp (for iam)
    proc_begX(0) = 1
    do pid = 1,npes-1
       proc_begX(pid) = proc_begX(pid-1) + proc_nXXX(pid-1)
    enddo
    procinfo%begp = proc_begX(iam)
    procinfo%endp = (procinfo%begp-1) + procinfo%npfts

    ! determine col offset for each clump assigned to each process
    ! (re-using proc_begX as work space)
    do cid = 1,nclumps
      clumps(cid)%begp = proc_begX(clumps(cid)%owner)
      proc_begX(clumps(cid)%owner) = proc_begX(clumps(cid)%owner) &
                                   + clumps(cid)%npfts
      clumps(cid)%endp = proc_begX(clumps(cid)%owner) - 1
    enddo

    ! COHORTS:
    ! calculate number of cohorts per process
    proc_nXXX(:) = 0
    do cid = 1,nclumps
       proc_nXXX(clumps(cid)%owner) = &
        proc_nXXX(clumps(cid)%owner) + clumps(cid)%nCohorts
    enddo

    ! determine offset (begCohort) for all processes,
    ! and then procinfo%begCohort and procinfo%endCohort (for iam)
    proc_begX(0) = 1
    do pid = 1,npes-1
       proc_begX(pid) = proc_begX(pid-1) + proc_nXXX(pid-1)
    enddo
    procinfo%begCohort = proc_begX(iam)
    procinfo%endCohort = (procinfo%begCohort-1) + procinfo%nCohorts

    ! determine col offset for each clump assigned to each process
    ! (re-using proc_begX as work space)
    do cid = 1,nclumps
      clumps(cid)%begCohort = proc_begX(clumps(cid)%owner)
      proc_begX(clumps(cid)%owner) = proc_begX(clumps(cid)%owner) &
                                   + clumps(cid)%nCohorts
      clumps(cid)%endCohort = proc_begX(clumps(cid)%owner) - 1
    enddo

    ! free work space
    deallocate(proc_nXXX, proc_begX)

    do n = 1,nclumps
       if (clumps(n)%ncells      /= allvecg(n,glev) .or. &
           clumps(n)%ntopounits  /= allvecg(n,tlev) .or. &
           clumps(n)%nlunits     /= allvecg(n,llev) .or. &
           clumps(n)%ncols       /= allvecg(n,clev) .or. &
           clumps(n)%npfts       /= allvecg(n,plev) .or. &
           clumps(n)%nCohorts    /= allvecg(n,hlev)) then

               write(iulog,*) 'decompInit_glcp(): allvecg error ncells ',iam,n,clumps(n)%ncells ,allvecg(n,glev)
               write(iulog,*) 'decompInit_glcp(): allvecg error topounits ',iam,n,clumps(n)%ntopounits,allvecg(n,tlev)
               write(iulog,*) 'decompInit_glcp(): allvecg error lunits ',iam,n,clumps(n)%nlunits,allvecg(n,llev)
               write(iulog,*) 'decompInit_glcp(): allvecg error ncols  ',iam,n,clumps(n)%ncols  ,allvecg(n,clev)
               write(iulog,*) 'decompInit_glcp(): allvecg error pfts   ',iam,n,clumps(n)%npfts  ,allvecg(n,plev)
               write(iulog,*) 'decompInit_glcp(): allvecg error cohorts ',iam,n,clumps(n)%nCohorts ,allvecg(n,hlev)

               call endrun(msg=errMsg(__FILE__, __LINE__))

       endif
    enddo

    deallocate(allvecg,allvecl)
    deallocate(lcid)

  end subroutine decompInit_clumps

  !------------------------------------------------------------------------------
  subroutine decompInit_gtlcp(lns,lni,lnj,glcmask)
    !
    ! !DESCRIPTION:
    ! Determine gsMaps for topounits, landunits, columns, pfts and cohorts
    !
    ! !USES:
    use spmdMod
    use spmdGathScatMod
    use subgridMod,       only : subgrid_get_gcellinfo
    use mct_mod
    !
    ! !ARGUMENTS:
    implicit none
    integer , intent(in) :: lns,lni,lnj ! land domain global size
    integer , pointer, optional   :: glcmask(:)  ! glc mask
    !
    ! !LOCAL VARIABLES:
    integer :: gi,ti,li,ci,pi,coi ! indices
    integer :: i,g,k,l,n,np       ! indices
    integer :: cid,pid            ! indices
    integer :: begg,endg          ! beg,end gridcells
    integer :: begt,endt          ! beg,end topographic units
    integer :: begl,endl          ! beg,end landunits
    integer :: begc,endc          ! beg,end columns
    integer :: begp,endp          ! beg,end pfts
    integer :: begCohort,endCohort    ! beg,end pfts
    integer :: numg               ! total number of gridcells across all processors
    integer :: numt               ! total number of topounits across all processors
    integer :: numl               ! total number of landunits across all processors
    integer :: numc               ! total number of columns across all processors
    integer :: nump               ! total number of pfts across all processors
    integer :: numCohort          ! ED cohorts
    integer :: icells             ! temporary
    integer :: itopounits         ! temporary
    integer :: ilunits            ! temporary
    integer :: icols              ! temporary
    integer :: ipfts              ! temporary
    integer :: icohorts           ! temporary
    integer :: ier                ! error code
    integer :: npmin,npmax,npint  ! do loop values for printing
    integer :: clmin,clmax        ! do loop values for printing
    integer :: locsize,globsize   ! used for gsMap init
    integer :: ng                 ! number of gridcells in gsMap_lnd_gdc2glo
    integer :: val1, val2         ! temporaries
    integer, pointer :: gindex(:) ! global index for gsMap init
    integer, pointer :: arrayglob(:) ! temporaroy
    integer, pointer :: gstart(:),  gcount(:)
    integer, pointer :: tstart(:),  tcount(:)
    integer, pointer :: lstart(:),  lcount(:)
    integer, pointer :: cstart(:),  ccount(:)
    integer, pointer :: pstart(:),  pcount(:)
    integer, pointer :: coStart(:), coCount(:)
    integer, pointer :: ioff(:)
    integer, parameter :: dbug=1      ! 0 = min, 1=normal, 2=much, 3=max
    character(len=32), parameter :: subname = 'decompInit_gtlcp'
    !------------------------------------------------------------------------------

    !init 

    call get_proc_bounds(begg, endg, begt, endt, begl, endl, begc, endc, begp, endp, &
         begCohort, endCohort)
    call get_proc_global(ng=numg, nt=numt, nl=numl, nc=numc, np=nump, nCohorts=numCohort)

    ! Determine global seg megs

    allocate(gstart(begg:endg))
    gstart(:) = 0
    allocate(gcount(begg:endg))
    gcount(:) = 0
    allocate(tstart(begg:endg))
    tstart(:) = 0
    allocate(tcount(begg:endg))
    tcount(:) = 0
    allocate(lstart(begg:endg))
    lstart(:) = 0
    allocate(lcount(begg:endg))
    lcount(:) = 0
    allocate(cstart(begg:endg))
    cstart(:) = 0
    allocate(ccount(begg:endg))
    ccount(:) = 0
    allocate(pstart(begg:endg))
    pstart(:) = 0
    allocate(pcount(begg:endg))
    pcount(:) = 0
    if ( use_fates ) then
       allocate(coStart(begg:endg))
       coStart(:) = 0
    endif
    allocate(coCount(begg:endg))
    coCount(:) = 0
    allocate(ioff(begg:endg)) 
    ioff(:) = 0

    ! Determine gcount, tcount, lcount, ccount and pcount

    do gi = begg,endg
       if (present(glcmask)) then
          call subgrid_get_gcellinfo (gi, ntopounits=itopounits, nlunits=ilunits, ncols=icols, npfts=ipfts, &
              ncohorts=icohorts, glcmask=glcmask(gi))
       else
          call subgrid_get_gcellinfo (gi, ntopounits=itopounits, nlunits=ilunits, ncols=icols, npfts=ipfts, &
               ncohorts=icohorts )
       endif
       gcount(gi)  = 1          ! number of gridcells for local gridcell index gi
       tcount(gi)  = itopounits ! number of topounits for local gridcell index gi
       lcount(gi)  = ilunits    ! number of landunits for local gridcell index gi
       ccount(gi)  = icols      ! number of columns for local gridcell index gi
       pcount(gi)  = ipfts      ! number of pfts for local gridcell index gi
       coCount(gi) = icohorts   ! number of ED cohorts for local gricell index gi
    enddo

    ! Determine gstart, lstart, cstart, pstart, coStart for the OUTPUT 1d data structures

    ! gather the gdc subgrid counts to masterproc in glo order
    ! compute glo ordered start indices from the counts
    ! scatter the subgrid start indices back out to the gdc gridcells
    ! set the local gindex array for the subgrid from the subgrid start and count arrays

    ng = mct_gsmap_gsize(gsmap_lnd_gdc2glo)
    allocate(arrayglob(ng))

    arrayglob(:) = 0
    call gather_data_to_master(gcount, arrayglob, grlnd)
    if (masterproc) then
       val1 = arrayglob(1)
       arrayglob(1) = 1
       do n = 2,ng
          val2 = arrayglob(n)
          arrayglob(n) = arrayglob(n-1) + val1
          val1 = val2
       enddo
    endif
    call scatter_data_from_master(gstart, arrayglob, grlnd)

    ! tstart for gridcell (n) is the total number of the topounits 
    ! over gridcells 1->n-1

    arrayglob(:) = 0
    call gather_data_to_master(tcount, arrayglob, grlnd)
    if (masterproc) then
       val1 = arrayglob(1)
       arrayglob(1) = 1
       do n = 2,ng
          val2 = arrayglob(n)
          arrayglob(n) = arrayglob(n-1) + val1
          val1 = val2
       enddo
    endif
    call scatter_data_from_master(tstart, arrayglob, grlnd)
    
    ! lstart for gridcell (n) is the total number of the landunits 
    ! over gridcells 1->n-1

    arrayglob(:) = 0
    call gather_data_to_master(lcount, arrayglob, grlnd)
    if (masterproc) then
       val1 = arrayglob(1)
       arrayglob(1) = 1
       do n = 2,ng
          val2 = arrayglob(n)
          arrayglob(n) = arrayglob(n-1) + val1
          val1 = val2
       enddo
    endif
    call scatter_data_from_master(lstart, arrayglob, grlnd)

    arrayglob(:) = 0
    call gather_data_to_master(ccount, arrayglob, grlnd)
    if (masterproc) then
       val1 = arrayglob(1)
       arrayglob(1) = 1
       do n = 2,ng
          val2 = arrayglob(n)
          arrayglob(n) = arrayglob(n-1) + val1
          val1 = val2
       enddo
    endif
    call scatter_data_from_master(cstart, arrayglob, grlnd)

    arrayglob(:) = 0
    call gather_data_to_master(pcount, arrayglob, grlnd)
    if (masterproc) then
       val1 = arrayglob(1)
       arrayglob(1) = 1
       do n = 2,ng
          val2 = arrayglob(n)
          arrayglob(n) = arrayglob(n-1) + val1
          val1 = val2
       enddo
    endif
    call scatter_data_from_master(pstart, arrayglob, grlnd)

    if ( use_fates ) then
       arrayglob(:) = 0
       call gather_data_to_master(coCount, arrayglob, grlnd)
       if (masterproc) then
          val1 = arrayglob(1)
          arrayglob(1) = 1
          do n = 2,ng
             val2 = arrayglob(n)
             arrayglob(n) = arrayglob(n-1) + val1
             val1 = val2
          enddo
       endif
       call scatter_data_from_master(coStart, arrayglob, grlnd)
    endif

    deallocate(arrayglob)

    ! Gridcell gsMap (compressed, no ocean points)

    allocate(gindex(begg:endg))
    i = begg-1
    do gi = begg,endg
       if (gcount(gi) <  1) then
          write(iulog,*) 'decompInit_gtlcp warning count g ',i,iam,gi,gcount(gi)
       endif
       do l = 1,gcount(gi)
          i = i + 1
          if (i < begg .or. i > endg) then
             write(iulog,*) 'decompInit_gtlcp error i ',i,begg,endg
             call endrun(msg=errMsg(__FILE__, __LINE__))
          endif
          gindex(i) = gstart(gi) + l - 1
       enddo
    enddo
    if (i /= endg) then
       write(iulog,*) 'decompInit_gtlcp error size ',i,begg,endg
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif
    locsize = endg-begg+1
    globsize = numg
    call mct_gsMap_init(gsMap_gce_gdc2glo, gindex, mpicom, comp_id, locsize, globsize)
    deallocate(gindex)

    ! Topounit gsmap

    allocate(gindex(begt:endt))
    ioff(:) = 0
    do ti = begt,endt
       gi = top_pp%gridcell(ti) !===this is determined internally from how landunits are spread out in memory
       gindex(ti) = tstart(gi) + ioff(gi) !=== the output gindex is ALWAYS the same regardless of how landuntis are spread out in memory
       ioff(gi)  = ioff(gi) + 1 
       ! check that this is less than [tstart(gi) + tcount(gi)]
    enddo
    locsize = endt-begt+1
    globsize = numt
    call mct_gsMap_init(gsMap_top_gdc2glo, gindex, mpicom, comp_id, locsize, globsize)
    deallocate(gindex)

    ! Landunit gsmap

    allocate(gindex(begl:endl))
    ioff(:) = 0
    do li = begl,endl
       gi = lun_pp%gridcell(li) !===this is determined internally from how landunits are spread out in memory
       gindex(li) = lstart(gi) + ioff(gi) !=== the output gindex is ALWAYS the same regardless of how landuntis are spread out in memory
       ioff(gi)  = ioff(gi) + 1 
       ! check that this is less than [lstart(gi) + lcount(gi)]
    enddo
    locsize = endl-begl+1
    globsize = numl
    call mct_gsMap_init(gsMap_lun_gdc2glo, gindex, mpicom, comp_id, locsize, globsize)
    deallocate(gindex)

    ! Column gsmap

    allocate(gindex(begc:endc))
    ioff(:) = 0
    do ci = begc,endc
       gi = col_pp%gridcell(ci)
       gindex(ci) = cstart(gi) + ioff(gi)
       ioff(gi) = ioff(gi) + 1 
       ! check that this is less than [cstart(gi) + ccount(gi)]
    enddo
    locsize = endc-begc+1
    globsize = numc
    call mct_gsMap_init(gsMap_col_gdc2glo, gindex, mpicom, comp_id, locsize, globsize)
    deallocate(gindex)

    ! PATCH gsmap

    allocate(gindex(begp:endp))
    ioff(:) = 0
    do pi = begp,endp
       gi = veg_pp%gridcell(pi)
       gindex(pi) = pstart(gi) + ioff(gi)
       ioff(gi) = ioff(gi) + 1 
       ! check that this is less than [pstart(gi) + pcount(gi)]
    enddo
    locsize = endp-begp+1
    globsize = nump
    call mct_gsMap_init(gsMap_patch_gdc2glo, gindex, mpicom, comp_id, locsize, globsize)
    deallocate(gindex)

    if ( use_fates ) then
       ! ED cohort gsMap
       allocate(gindex(begCohort:endCohort))
       ioff(:) = 0
       gi = begg
       do coi = begCohort,endCohort
          gindex(coi) = coStart(gi) + ioff(gi)
          ioff(gi) = ioff(gi) + 1
          if ( mod(coi, fates_maxElementsPerSite ) == 0 ) gi = gi + 1
       enddo
       locsize = endCohort-begCohort+1
       globsize = numCohort
       call mct_gsMap_init(gsMap_cohort_gdc2glo, gindex, mpicom, comp_id, locsize, globsize)
       deallocate(gindex)
    endif

    ! Deallocate start/count arrays
    deallocate(gstart, gcount)
    deallocate(tstart, tcount)
    deallocate(lstart, lcount)
    deallocate(cstart, ccount)
    deallocate(pstart, pcount)
    if ( use_fates ) then
       deallocate(coStart,coCount)
    endif
    deallocate(ioff)

    ! Diagnostic output

    if (masterproc) then
       write(iulog,*)' Surface Grid Characteristics'
       write(iulog,*)'   longitude points          = ',lni
       write(iulog,*)'   latitude points           = ',lnj
       write(iulog,*)'   total number of gridcells = ',numg
       write(iulog,*)'   total number of topounits = ',numt
       write(iulog,*)'   total number of landunits = ',numl
       write(iulog,*)'   total number of columns   = ',numc
       write(iulog,*)'   total number of pfts      = ',nump
       write(iulog,*)'   total number of cohorts   = ',numCohort
       write(iulog,*)' Decomposition Characteristics'
       write(iulog,*)'   clumps per process        = ',clump_pproc
       write(iulog,*)' gsMap Characteristics'
       write(iulog,*) '  lnd gsmap glo num of segs = ',mct_gsMap_ngseg(gsMap_lnd_gdc2glo)
       write(iulog,*) '  gce gsmap glo num of segs = ',mct_gsMap_ngseg(gsMap_gce_gdc2glo)
       write(iulog,*) '  top gsmap glo num of segs = ',mct_gsMap_ngseg(gsMap_top_gdc2glo)
       write(iulog,*) '  lun gsmap glo num of segs = ',mct_gsMap_ngseg(gsMap_lun_gdc2glo)
       write(iulog,*) '  col gsmap glo num of segs = ',mct_gsMap_ngseg(gsMap_col_gdc2glo)
       write(iulog,*) '  pft gsmap glo num of segs = ',mct_gsMap_ngseg(gsMap_patch_gdc2glo)
       write(iulog,*) '  coh gsmap glo num of segs = ',mct_gsMap_ngseg(gsMap_cohort_gdc2glo)
       write(iulog,*)
    end if

    ! Write out clump and proc info, one pe at a time, 
    ! barrier to control pes overwriting each other on stdout

    call shr_sys_flush(iulog)
    call mpi_barrier(mpicom,ier)
    npmin = 0
    npmax = npes-1
    npint = 1
    if (dbug == 0) then
       npmax = 0
    elseif (dbug == 1) then
       npmax = min(npes-1,4)
    elseif (dbug == 2) then
       npint = npes/8
    endif
    do np = npmin,npmax,npint
       pid = np
       if (dbug == 1) then
          if (np == 2) pid=npes/2-1
          if (np == 3) pid=npes-2
          if (np == 4) pid=npes-1
       endif
       pid = max(pid,0)
       pid = min(pid,npes-1)

       if (iam == pid) then
          write(iulog,*)
          write(iulog,*)'proc= ',pid,&
               ' beg gridcell= ',procinfo%begg, &
               ' end gridcell= ',procinfo%endg,                   &
               ' total gridcells per proc= ',procinfo%ncells
          write(iulog,*)'proc= ',pid,&
               ' beg topounit= ',procinfo%begt, &
               ' end topounit= ',procinfo%endt,                   &
               ' total topounits per proc= ',procinfo%ntopounits
          write(iulog,*)'proc= ',pid,&
               ' beg landunit= ',procinfo%begl, &
               ' end landunit= ',procinfo%endl,                   &
               ' total landunits per proc= ',procinfo%nlunits
          write(iulog,*)'proc= ',pid,&
               ' beg column  = ',procinfo%begc, &
               ' end column  = ',procinfo%endc,                   &
               ' total columns per proc  = ',procinfo%ncols
          write(iulog,*)'proc= ',pid,&
               ' beg pft     = ',procinfo%begp, &
               ' end pft     = ',procinfo%endp,                   &
               ' total pfts per proc     = ',procinfo%npfts
          write(iulog,*)'proc= ',pid,&
               ' beg coh     = ',procinfo%begCohort, &
               ' end coh     = ',procinfo%endCohort,                   &
               ' total coh per proc     = ',procinfo%nCohorts
          write(iulog,*)'proc= ',pid,&
               ' lnd ngseg   = ',mct_gsMap_ngseg(gsMap_lnd_gdc2glo), &
               ' lnd nlseg   = ',mct_gsMap_nlseg(gsMap_lnd_gdc2glo,iam)
          write(iulog,*)'proc= ',pid,&
               ' gce ngseg   = ',mct_gsMap_ngseg(gsMap_gce_gdc2glo), &
               ' gce nlseg   = ',mct_gsMap_nlseg(gsMap_gce_gdc2glo,iam)
          write(iulog,*)'proc= ',pid,&
               ' top ngseg   = ',mct_gsMap_ngseg(gsMap_top_gdc2glo), &
               ' top nlseg   = ',mct_gsMap_nlseg(gsMap_top_gdc2glo,iam)
          write(iulog,*)'proc= ',pid,&
               ' lun ngseg   = ',mct_gsMap_ngseg(gsMap_lun_gdc2glo), &
               ' lun nlseg   = ',mct_gsMap_nlseg(gsMap_lun_gdc2glo,iam)
          write(iulog,*)'proc= ',pid,&
               ' col ngseg   = ',mct_gsMap_ngseg(gsMap_col_gdc2glo), &
               ' col nlseg   = ',mct_gsMap_nlseg(gsMap_col_gdc2glo,iam)
          write(iulog,*)'proc= ',pid,&
               ' pft ngseg   = ',mct_gsMap_ngseg(gsMap_patch_gdc2glo), &
               ' pft nlseg   = ',mct_gsMap_nlseg(gsMap_patch_gdc2glo,iam)
          write(iulog,*)'proc= ',pid,&
               ' coh ngseg   = ',mct_gsMap_ngseg(gsMap_cohort_gdc2glo), &
               ' coh nlseg   = ',mct_gsMap_nlseg(gsMap_cohort_gdc2glo,iam)
          write(iulog,*)'proc= ',pid,' nclumps = ',procinfo%nclumps

          clmin = 1
          clmax = procinfo%nclumps
          if (dbug == 1) then
            clmax = 1
          elseif (dbug == 0) then
            clmax = -1
          endif
          do n = clmin,clmax
             cid = procinfo%cid(n)
             write(iulog,*)'proc= ',pid,' clump no = ',n, &
                  ' clump id= ',procinfo%cid(n),    &
                  ' beg gridcell= ',clumps(cid)%begg, &
                  ' end gridcell= ',clumps(cid)%endg, &
                  ' total gridcells per clump= ',clumps(cid)%ncells
             write(iulog,*)'proc= ',pid,' clump no = ',n, &
                  ' clump id= ',procinfo%cid(n),    &
                  ' beg topounit= ',clumps(cid)%begt, &
                  ' end topounit= ',clumps(cid)%endt, &
                  ' total topounits per clump = ',clumps(cid)%ntopounits
             write(iulog,*)'proc= ',pid,' clump no = ',n, &
                  ' clump id= ',procinfo%cid(n),    &
                  ' beg landunit= ',clumps(cid)%begl, &
                  ' end landunit= ',clumps(cid)%endl, &
                  ' total landunits per clump = ',clumps(cid)%nlunits
             write(iulog,*)'proc= ',pid,' clump no = ',n, &
                  ' clump id= ',procinfo%cid(n),    &
                  ' beg column  = ',clumps(cid)%begc, &
                  ' end column  = ',clumps(cid)%endc, &
                  ' total columns per clump  = ',clumps(cid)%ncols
             write(iulog,*)'proc= ',pid,' clump no = ',n, &
                  ' clump id= ',procinfo%cid(n),    &
                  ' beg pft     = ',clumps(cid)%begp, &
                  ' end pft     = ',clumps(cid)%endp, &
                  ' total pfts per clump     = ',clumps(cid)%npfts
             write(iulog,*)'proc= ',pid,' clump no = ',n, &
                  ' clump id= ',procinfo%cid(n),    &
                  ' beg cohort     = ',clumps(cid)%begCohort, &
                  ' end cohort     = ',clumps(cid)%endCohort, &
                  ' total cohorts per clump     = ',clumps(cid)%nCohorts
          end do
       end if
       call shr_sys_flush(iulog)
       call mpi_barrier(mpicom,ier)
    end do
    call shr_sys_flush(iulog)

  end subroutine decompInit_gtlcp

  !------------------------------------------------------------------------------
  subroutine decompInit_lnd_using_gp(lni, lnj, cellsOnCell, ncells_loc, maxEdges, &
       amask)
    !
    ! !DESCRIPTION:
    ! This subroutine initializes the land surface decomposition into a clump
    ! data structure using graph partitioning approach.  This assumes each pe
    ! has the same number of clumps set by clump_pproc.
    !
#ifdef USE_PETSC_LIB
#include <petsc/finclude/petsc.h>
#endif
    ! !USES:
    use clm_varctl, only : nsegspc
#ifdef USE_PETSC_LIB
    use petscsys
    use petscvec
    use petscmat
    use petscdm
#endif
    !
    ! !ARGUMENTS:
    implicit none
    !
    !
    integer , intent(in) :: amask(:)
    integer , intent(in) :: lni,lnj                     ! domain global size
    integer , intent(in) :: cellsOnCell(:,:)
    integer , intent(in) :: ncells_loc
    integer , intent(in) :: maxEdges
    !
    ! !LOCAL VARIABLES:
    integer            :: lns                           ! global domain size
    integer            :: ln,lj                         ! indices
    integer            :: numg                          ! number of land gridcells
    integer            :: cid,pid                       ! indices
    integer            :: n,i,m                         ! indices
    integer            :: ier                           ! error code
    integer            :: beg,end,lsize,gsize           ! used for gsmap init
    integer, pointer   :: gindex(:)                     ! global index for gsmap init
    integer            :: icell, iedge                  ! indices
    integer            :: offset                        ! temporary
    integer            :: cell_id_offset                ! temporary
    integer            :: count                         ! temporary
    integer            :: num_rows, num_cols            ! temporary
    integer            :: istart, iend                  ! temporary
    integer            :: ncells_tot                    ! total number of grid cells
    integer            :: ncells_owned                  ! number of grid cells owned by a processor after domain decomposition
    integer            :: ncells_per_clump              ! number of grid cells per clump
    integer            :: remainder                     ! temporary
    integer            :: cowner                        ! clump owner
    integer, pointer   :: i_index(:)                    ! temporary
    integer, pointer   :: j_index(:)                    ! temporary
    integer, pointer   :: local_conn_offsets(:)         ! temporary
    integer, pointer   :: local_conn(:)                 ! temporary
    integer, pointer   :: clump_ncells(:)               ! temporary
    integer, pointer   :: clump_begg(:)                 ! temporary
    integer, pointer   :: clump_endg(:)                 ! temporary
    integer, pointer   :: local_clump_info(:)           ! temporary
    integer, pointer   :: global_clump_info(:)          ! temporary
    integer, pointer   :: thread_count(:)               ! temporary
    integer, pointer   :: int_array(:)                  ! temporary
    integer, pointer   :: ncells_count(:)               ! temporary
#ifdef USE_PETSC_LIB
    PetscReal, pointer :: real_ptr(:)                   ! temporary
    PetscInt, pointer  :: int_ptr(:)                    ! temporary
    PetscBool          :: success                       ! temporary
    VecScatter         :: scatter                       ! temporary
    Vec                :: ids_old                       ! grid cell IDs before domain decomposition
    Vec                :: ids_new                       ! grid cell IDs after domain decomposition
    Vec                :: lcid_aft_decomp               ! local clump ID for grid cells owned after domain decomposition
    Vec                :: lcid_aft_decomp_for_all_procs ! local clump ID for all grid cells after domain decomposition
    Mat                :: Dual_mat                      ! dual matrix
    Mat                :: Dual_aij                      ! dual matrix in aij format
    Mat                :: Dual_aij_loc                  ! local part of the dual matrix in aij format
    MatPartitioning    :: part                          ! partitioning of dual matrix
    IS                 :: is_new_owner, is_new_id       ! index set that stores mpi rank and id of grid cell after domain decomposition
    IS                 :: is_from, is_to                ! temporary
    PetscErrorCode     :: ierr                          ! get error code from PETSc
#endif
    character(len=255) :: subname = 'decompInit_lnd_using_gp'
    !------------------------------------------------------------------------------

#ifndef USE_PETSC_LIB

    call endrun(msg='ERROR ' // trim(subname) //': Graph partitioning requires '//&
         'PETSc, but the code was compiled without -DUSE_PETSC_LIB')

#else

    lns = lni * lnj

    !--- set and verify nclumps ---
    if (clump_pproc > 0) then
       nclumps = clump_pproc * npes
       if (nclumps < npes) then
          write(iulog,*) trim(subname) // '(): Number of gridcell clumps= ',nclumps, &
               ' is less than the number of processes = ', npes
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if
    else
       write(iulog,*)trim(subname) // '(): clump_pproc= ',clump_pproc,'  must be greater than 0'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    ! allocate and initialize procinfo and clumps
    ! beg and end indices initialized for simple addition of cells later

    allocate(procinfo%cid(clump_pproc), stat=ier)
    if (ier /= 0) then
       write(iulog,*) trim(subname) // '(): allocation error for procinfo%cid'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    procinfo%nclumps   = clump_pproc
    procinfo%cid(:)    = -1
    procinfo%ncells    = 0
    procinfo%nlunits   = 0
    procinfo%ncols     = 0
    procinfo%npfts     = 0
    procinfo%nCohorts  = 0
    procinfo%begg      = 1
    procinfo%begl      = 1
    procinfo%begc      = 1
    procinfo%begp      = 1
    procinfo%begCohort = 1
    procinfo%endg      = 0
    procinfo%endl      = 0
    procinfo%endc      = 0
    procinfo%endp      = 0
    procinfo%endCohort = 0

    allocate(clumps(nclumps), stat=ier)
    if (ier /= 0) then
       write(iulog,*) trim(subname) // '(): allocation error for clumps'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if
    clumps(:)%owner     = -1
    clumps(:)%ncells    = 0
    clumps(:)%nlunits   = 0
    clumps(:)%ncols     = 0
    clumps(:)%npfts     = 0
    clumps(:)%nCohorts  = 0
    clumps(:)%begg      = 1
    clumps(:)%begl      = 1
    clumps(:)%begc      = 1
    clumps(:)%begp      = 1
    clumps(:)%begCohort = 1
    clumps(:)%endg      = 0
    clumps(:)%endl      = 0
    clumps(:)%endc      = 0
    clumps(:)%endp      = 0
    clumps(:)%endCohort = 0

    ! assign clumps to proc round robin
    cid = 0
    do n = 1,nclumps
       pid = mod(n-1,npes)
       if (pid < 0 .or. pid > npes-1) then
          write(iulog,*) trim(subname) // '(): round robin pid error ',n,pid,npes
          call endrun(msg=errMsg(__FILE__, __LINE__))
       endif
       clumps(n)%owner = pid
       if (iam == pid) then
          cid = cid + 1
          if (cid < 1 .or. cid > clump_pproc) then
             write(iulog,*) trim(subname) // '(): round robin pid error ',n,pid,npes
             call endrun(msg=errMsg(__FILE__, __LINE__))
          endif
          procinfo%cid(cid) = n
       endif
    enddo

    ! count total active land gridcells
    numg = 0
    do ln = 1,lns
       if (amask(ln) == 1) then
          numg = numg + 1
       endif
    enddo

    if (npes > numg) then
       write(iulog,*) trim(subname) // '(): Number of processes exceeds number ', &
            'of land grid cells',npes,numg
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    if (nclumps > numg) then
       write(iulog,*) trim(subname) // '(): Number of clumps exceeds number ', &
            'of land grid cells',nclumps,numg
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    if (numg /= lns) then
       write(iulog,*) trim(subname) // '(): Only implimented for numg == lns '
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    ! Determine the cell id offset on each processor
    cell_id_offset = 0
    call MPI_Scan(ncells_loc, cell_id_offset, 1, MPIU_INTEGER, MPI_SUM, mpicom, ierr)
    cell_id_offset = cell_id_offset - ncells_loc

    ! Determine the total number of grid cells
    call MPI_Allreduce(ncells_loc, ncells_tot, 1, MPIU_INTEGER, MPI_SUM, mpicom, ierr)

    ! Create a Dual matrix.
    call MatCreateAIJ(mpicom , & ! comm
         ncells_loc          , & ! m
         PETSC_DECIDE        , & ! n
         PETSC_DETERMINE     , & ! M
         ncells_tot          , & ! N
         maxEdges            , & ! d_nz
         PETSC_NULL_INTEGER  , & ! d_nnz
         maxEdges            , & ! o_nz
         PETSC_NULL_INTEGER  , & ! o_nnz
         Dual_aij            , & ! Mat
         ierr);CHKERRQ(ierr)

    !
    ! If cell_1 and cell_2 are neighbors, then set:
    !
    !  Dual_aij(cell_1, cell_2) = 1
    !  Dual_aij(cell_2, cell_1) = 1
    !
    do icell = 1, ncells_loc
       do iedge = 1, maxEdges

          if (cellsOnCell(iedge, icell) > 0) then

             call MatSetValue(Dual_aij        , &
                  icell+cell_id_offset-1      , &
                  cellsOnCell(iedge, icell)-1 , &
                  1.d0                        , &
                  INSERT_VALUES               , &
                  ierr);CHKERRQ(ierr)

             call MatSetValue(Dual_aij        , &
                  cellsOnCell(iedge, icell)-1 , &
                  icell+cell_id_offset-1      , &
                  1.d0                        , &
                  INSERT_VALUES               , &
                  ierr);CHKERRQ(ierr)
          end if
       end do
    end do

    call MatAssemblyBegin( Dual_aij,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(   Dual_aij,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

    !
    ! Create sparse matrix representing an adjacency list from Dual matrix
    !

    if ( npes > 1 ) then

       call MatMPIAIJGetLocalMat(Dual_aij , &
            MAT_INITIAL_MATRIX            , &
            Dual_aij_loc                  , &
            ierr);CHKERRQ(ierr)

       call MatGetRowIJF90(Dual_aij_loc   , &
            0                             , &
            PETSC_FALSE                   , &
            PETSC_FALSE                   , &
            num_rows                      , &
            i_index                       , &
            j_index                       , &
            success                       , &
            ierr);CHKERRQ(ierr)
    else

       call MatGetRowIJF90(Dual_aij       , &
            0                             , &
            PETSC_FALSE                   , &
            PETSC_FALSE                   , &
            num_rows                      , &
            i_index                       , &
            j_index                       , &
            success                       , &
            ierr);CHKERRQ(ierr)
    endif

    count = 0
    do icell = 1,num_rows
       istart   = i_index(icell)
       iend     = i_index(icell + 1) - 1
       num_cols = iend - istart + 1
       count    = count + num_cols
    enddo

    allocate (local_conn         (count      ))
    allocate (local_conn_offsets (num_rows+1 ))

    local_conn_offsets (1:num_rows+1) = i_index (1:num_rows+1)
    local_conn         (1:count)      = j_index (1:count     )

    call MatCreateMPIAdj(mpicom , &
         ncells_loc             , &
         ncells_tot             , &
         local_conn_offsets     , &
         local_conn             , &
         PETSC_NULL_INTEGER     , &
         Dual_mat               , &
         ierr); CHKERRQ(ierr)

    if ( npes > 1 ) then
       call MatRestoreRowIJF90(Dual_aij_loc , &
            0                               , &
            PETSC_FALSE                     , &
            PETSC_FALSE                     , &
            num_rows                        , &
            i_index                         , &
            j_index                         , &
            success                         , &
            ierr);CHKERRQ(ierr)
    else
       call MatGetRowIJF90(Dual_aij , &
            0                       , &
            PETSC_FALSE             , &
            PETSC_FALSE             , &
            num_rows                , &
            i_index                 , &
            j_index                 , &
            success                 , &
            ierr);CHKERRQ(ierr)
    endif
    call MatDestroy(Dual_aij,ierr);CHKERRQ(ierr)

    !
    ! Use graph partitioning to decompose the mesh
    !

    call MatPartitioningCreate(mpicom, part, ierr);CHKERRQ(ierr)

    call MatPartitioningSetAdjacency(part, Dual_mat, ierr);CHKERRQ(ierr)

    call MatPartitioningSetFromOptions(part, ierr);CHKERRQ(ierr)

    ! Now perform graph partioning.
    !   - After the call to MatPartitioningApply(),
    !     is_new_owner has entries that corresponds to the rank of
    !     processor which owns each grid cell
    call MatPartitioningApply(part, is_new_owner, ierr);CHKERRQ(ierr)

    ! Free up memory
    call MatDestroy(Dual_mat, ierr); CHKERRQ(ierr);
    call MatPartitioningDestroy(part, ierr);CHKERRQ(ierr)
    deallocate(local_conn        )
    deallocate(local_conn_offsets)

    ! Determine the number of grid cells owned after graph partitioning
    allocate(ncells_count(npes))
    call ISPartitioningCount(is_new_owner, npes, ncells_count, ierr);CHKERRQ(ierr)
    ncells_owned = ncells_count(iam + 1)
    deallocate(ncells_count)

    ! Determine the new ids of grid cells
    call ISPartitioningToNumbering(is_new_owner, is_new_id, ierr);CHKERRQ(ierr)

    !
    ! Compute information for each processor
    !
    procinfo%ncells = ncells_owned

    offset = 0
    call MPI_Scan(ncells_owned, offset, 1, MPIU_INTEGER, MPI_SUM, mpicom, ierr)
    procinfo%begg = offset + 1 - ncells_owned

    offset = 0
    call MPI_scan(ncells_owned, offset, 1, MPIU_INTEGER, MPI_SUM, mpicom, ierr)
    procinfo%endg = offset

    !
    ! Compute information about all clumps on all processors
    !

    ncells_per_clump = procinfo%ncells/clump_pproc
    remainder        = procinfo%ncells - ncells_per_clump*clump_pproc

    allocate (clump_ncells      (clump_pproc            ))
    allocate (clump_begg        (clump_pproc            ))
    allocate (clump_endg        (clump_pproc            ))

    allocate (local_clump_info  (0:3*clump_pproc-1      ))
    allocate (global_clump_info (0:3*clump_pproc*npes-1 ))

    allocate (thread_count      (0:npes-1               ))

    clump_ncells(:) = ncells_per_clump

    offset = procinfo%begg
    do m = 1,clump_pproc
       if (m-1 < remainder) clump_ncells(m) = clump_ncells(m) + 1

       clump_begg(m) = offset
       clump_endg(m) = offset + clump_ncells(m) - 1
       offset        = offset + clump_ncells(m)

       local_clump_info((m-1)*3 + 0) = clump_ncells(m)
       local_clump_info((m-1)*3 + 1) = clump_begg  (m)
       local_clump_info((m-1)*3 + 2) = clump_endg  (m)
    end do

    call MPI_Allgather(local_clump_info, 3*clump_pproc, MPI_INTEGER, &
         global_clump_info, 3*clump_pproc, MPI_INTEGER, mpicom, ier)

    count = 0
    thread_count(:) = 0
    do m = 1, nclumps
       cowner = clumps(m)%owner

       clumps(m)%ncells = global_clump_info(cowner*3 + thread_count(cowner)*3    )
       clumps(m)%begg   = global_clump_info(cowner*3 + thread_count(cowner)*3 + 1)
       clumps(m)%endg   = global_clump_info(cowner*3 + thread_count(cowner)*3 + 2)

       thread_count(cowner) = thread_count(cowner) + 1
    enddo

    deallocate (clump_ncells      )
    deallocate (clump_begg        )
    deallocate (clump_endg        )
    deallocate (local_clump_info  )
    deallocate (global_clump_info )
    deallocate (thread_count      )

    !
    ! Determine the natural ids of the grid cells that each processor owns after
    ! domain decomposition. This information will be used to set up ldecomp%gdc2glo(:).
    !
    call VecCreateMPI(mpicom, ncells_loc  , PETSC_DETERMINE, ids_old, ierr);CHKERRQ(ierr)
    call VecCreateMPI(mpicom, ncells_owned, PETSC_DETERMINE, ids_new, ierr);CHKERRQ(ierr)

    ! Create a vector containing old ids of grid cells (i.e. ids before
    ! domain decomposition)
    call VecGetArrayF90(ids_old, real_ptr, ierr)
    do i = 1, ncells_loc
       real_ptr(i) = i + cell_id_offset
    end do
    call VecRestoreArrayF90(ids_old, real_ptr, ierr)

    ! Create an IS to scatter data stored in ids_old
    allocate(int_array(ncells_loc))
    do m = 1, ncells_loc
       int_array(m) = m - 1 + cell_id_offset
    end do

    call ISCreateGeneral(mpicom, ncells_loc, int_array, &
         PETSC_COPY_VALUES, is_from, ierr);CHKERRQ(ierr);
    deallocate(int_array)

    ! Create a VectorScatter
    call VecScatterCreate(ids_old, is_from, ids_new, is_new_id, &
         scatter, ierr); CHKERRQ(ierr);

    call ISDestroy(is_from, ierr); CHKERRQ(ierr)
    call ISDestroy(is_new_id, ierr); CHKERRQ(ierr)

    ! Scatter data to get new ids of grid cells (i.e. ids after
    ! domain decomposition)
    call VecScatterBegin(scatter, ids_old, ids_new, &
         INSERT_VALUES, SCATTER_FORWARD, ierr); CHKERRQ(ierr)
    call VecScatterEnd(scatter, ids_old, ids_new, &
         INSERT_VALUES, SCATTER_FORWARD, ierr); CHKERRQ(ierr)

    ! Set ldecomp

    call get_proc_bounds(beg, end)
    allocate(ldecomp%gdc2glo(beg:end), stat=ier)
     if (ier /= 0) then
       write(iulog,*) trim(subname) // '(): allocation error1 for ldecomp, etc'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    ldecomp%gdc2glo(:) = 0

    call VecGetArrayF90(ids_new, real_ptr, ierr); CHKERRQ(ierr);
    do m = beg, end
       ldecomp%gdc2glo(m) = INT(real_ptr(m-beg+1))
    end do
    call VecRestoreArrayF90(ids_new, real_ptr, ierr); CHKERRQ(ierr);

    call VecDestroy(ids_old, ierr); CHKERRQ(ierr);
    call VecDestroy(ids_new, ierr); CHKERRQ(ierr);

    !
    ! For each grid cell, identify the processor that owns
    ! the grid cell after domain decomposition.
    !
    ! NOTE: This is not a scalable approach as each processor
    !       stores the information about the entire domain.
    !       But, decompInit_clumps() needs this information.
    !

    call VecCreateMPI(mpicom           , &
         ncells_loc                    , &
         PETSC_DETERMINE               , &
         lcid_aft_decomp               , &
         ierr); CHKERRQ(ierr)

    call VecCreateMPI(mpicom           , &
         numg                          , &
         PETSC_DETERMINE               , &
         lcid_aft_decomp_for_all_procs , &
         ierr); CHKERRQ(ierr)

    call ISGetIndicesF90(is_new_owner   , int_ptr , ierr); CHKERRQ(ierr)
    call VecGetArrayF90( lcid_aft_decomp, real_ptr, ierr); CHKERRQ(ierr)
    real_ptr(:) = real(int_ptr(:) + 1._r8)
    call ISRestoreIndicesF90(is_new_owner   , int_ptr , ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90( lcid_aft_decomp, real_ptr, ierr); CHKERRQ(ierr)

    allocate(int_array(numg))
    do m = 1, numg
       int_array(m) = m - 1
    end do
    call ISCreateGeneral(mpicom, numg, int_array, &
         PETSC_COPY_VALUES, is_from, ierr);CHKERRQ(ierr);
    deallocate(int_array)

    allocate(int_array(numg))
    do m = 1, numg
       int_array(m) = m - 1 + iam*numg
    end do
    call ISCreateGeneral(mpicom, numg, int_array, &
         PETSC_COPY_VALUES, is_to, ierr);CHKERRQ(ierr);
    deallocate(int_array)

    call VecScatterCreate(lcid_aft_decomp, is_from, &
         lcid_aft_decomp_for_all_procs, is_to, scatter, ierr);CHKERRQ(ierr)

    call ISDestroy(is_from, ierr); CHKERRQ(ierr)
    call ISDestroy(is_to  , ierr); CHKERRQ(ierr)

    call VecScatterBegin(scatter, lcid_aft_decomp, lcid_aft_decomp_for_all_procs, &
         INSERT_VALUES, SCATTER_FORWARD, ierr); CHKERRQ(ierr)
    call VecScatterEnd(scatter, lcid_aft_decomp, lcid_aft_decomp_for_all_procs, &
         INSERT_VALUES, SCATTER_FORWARD, ierr); CHKERRQ(ierr)
    call VecScatterDestroy(scatter, ierr); CHKERRQ(ierr)

    ! Set lcid

    allocate(lcid(lns))
    lcid(:) = 0

    call VecGetArrayF90(lcid_aft_decomp_for_all_procs, real_ptr, ierr); CHKERRQ(ierr)
    lcid(:) = INT(real_ptr(:))
    call VecRestoreArrayF90(lcid_aft_decomp_for_all_procs, real_ptr, ierr); CHKERRQ(ierr)

    call ISDestroy(is_new_owner, ierr);CHKERRQ(ierr)

    ! Set gsMap_lnd_gdc2glo (the global index here includes mask=0 or ocean points)

    call get_proc_bounds(beg, end)
    allocate(gindex(beg:end))
    do n = beg,end
       gindex(n) = ldecomp%gdc2glo(n)
    enddo
    lsize = end-beg+1
    gsize = lni * lnj

    call mct_gsMap_init(gsMap_lnd_gdc2glo, gindex, mpicom, comp_id, lsize, gsize)
    deallocate(gindex)

    ! Diagnostic output

    if (masterproc) then
       write(iulog,*)' Surface Grid Characteristics'
       write(iulog,*)'   longitude points               = ',lni
       write(iulog,*)'   latitude points                = ',lnj
       write(iulog,*)'   total number of land gridcells = ',numg
       write(iulog,*)' Decomposition Characteristics'
       write(iulog,*)'   clumps per process             = ',clump_pproc
       write(iulog,*)' gsMap Characteristics'
       write(iulog,*) '  lnd gsmap glo num of segs      = ',mct_gsMap_ngseg(gsMap_lnd_gdc2glo)
       write(iulog,*)
    end if

    call shr_sys_flush(iulog)

#endif

  end subroutine decompInit_lnd_using_gp

  !------------------------------------------------------------------------------
  subroutine decompInit_ghosts(glcmask)
    !
    ! !DESCRIPTION:
    ! On each proc, determine the number of following ghost/halo subgrid quantities:
    !  - topounits
    !  - landunits,
    !  - columns,
    !  - PFTs, and
    !  - cohorts.
    !
    ! !USES:
    use clm_varctl           , only : lateral_connectivity
    use subgridMod           , only : subgrid_get_gcellinfo
#ifdef USE_PETSC_LIB
    use domainLateralMod     , only : ldomain_lateral
    use UnstructuredGridType , only : ScatterDataG2L
#endif
    !
    implicit none
    !
    ! !ARGUMENTS:
    integer , pointer, optional  :: glcmask(:)             ! glc mask
    !
    ! !LOCAL VARIABLES:
    integer                      :: begg,endg              ! begin/end indices for grid
    integer                      :: anumg                  ! lnd num gridcells
    integer                      :: ln                     ! temporary
    integer                      :: nblocks                ! block size for PETSc vector
    integer                      :: itopounits             ! temporary
    integer                      :: ilunits                ! temporary
    integer                      :: icols                  ! temporary
    integer                      :: ipfts                  ! temporary
    integer                      :: icohorts               ! temporary
    integer                      :: ighost                 ! temporary
    integer                      :: ighost_beg, ighost_end ! temporary
    real(r8), pointer            :: data_send(:)
    real(r8), pointer            :: data_recv(:)
    integer                      :: ndata_send
    integer                      :: ndata_recv
    character(len=32), parameter :: subname = 'decompInit_ghosts'

    if (.not.lateral_connectivity) then

       ! No ghost cells
       procinfo%ncells_ghost    = 0
       procinfo%ntopounits_ghost   = 0
       procinfo%nlunits_ghost   = 0
       procinfo%ncols_ghost     = 0
       procinfo%npfts_ghost     = 0
       procinfo%nCohorts_ghost  = 0

       procinfo%begg_ghost      = 0
       procinfo%begt_ghost      = 0
       procinfo%begl_ghost      = 0
       procinfo%begc_ghost      = 0
       procinfo%begp_ghost      = 0
       procinfo%begCohort_ghost = 0
       procinfo%endg_ghost      = 0
       procinfo%endt_ghost      = 0
       procinfo%endl_ghost      = 0
       procinfo%endc_ghost      = 0
       procinfo%endp_ghost      = 0
       procinfo%endCohort_ghost = 0

       ! All = local (as no ghost cells)
       procinfo%ncells_all      = procinfo%ncells
       procinfo%ntopounits_all  = procinfo%ntopounits
       procinfo%nlunits_all     = procinfo%nlunits
       procinfo%ncols_all       = procinfo%ncols
       procinfo%npfts_all       = procinfo%npfts
       procinfo%nCohorts_all    = procinfo%nCohorts

       procinfo%begg_all        = procinfo%begg
       procinfo%begt_all        = procinfo%begt
       procinfo%begl_all        = procinfo%begl
       procinfo%begc_all        = procinfo%begc
       procinfo%begp_all        = procinfo%begp
       procinfo%begCohort_all   = procinfo%begCohort
       procinfo%endg_all        = procinfo%endg
       procinfo%endt_all        = procinfo%endt
       procinfo%endl_all        = procinfo%endl
       procinfo%endc_all        = procinfo%endc
       procinfo%endp_all        = procinfo%endp
       procinfo%endCohort_all   = procinfo%endCohort

    else

#ifndef USE_PETSC_LIB

    call endrun(msg='ERROR ' // trim(subname) //': decompInit_ghosts requires '//&
         'PETSc, but the code was compiled without -DUSE_PETSC_LIB')

#else
       call get_proc_bounds(begg, endg)

       ! Approach:
       ! 1) For a global PETSc vector, save the number of subgrid
       !    quantities for each grid cell.
       ! 2) Scatter the global PETSc vector to a local PETSc vector
       ! 3) Finally count the number of subgrid quantities for all
       !    ghost grid cells in the local PETSc vector

       nblocks = 5 ! topo + lun + col + pft + cohort

       ndata_send = nblocks*ldomain_lateral%ugrid%ngrid_local
       ndata_recv = nblocks*ldomain_lateral%ugrid%ngrid_ghosted

       allocate(data_send(ndata_send))
       allocate(data_recv(ndata_recv))

       data_send(:) = 0.d0

       ! Save information about number of subgrid categories for
       ! local grid cells

       do anumg = begg,endg
          ln  = anumg
          if (present(glcmask)) then
             call subgrid_get_gcellinfo (ln, ntopounits=itopounits, nlunits=ilunits, ncols=icols, npfts=ipfts, &
                  ncohorts=icohorts, glcmask=glcmask(ln))
          else
             call subgrid_get_gcellinfo (ln, ntopounits=itopounits, nlunits=ilunits, ncols=icols, npfts=ipfts, &
                  ncohorts=icohorts )
          endif

          data_send((anumg-begg)*nblocks + 1) = itopounits
          data_send((anumg-begg)*nblocks + 2) = ilunits
          data_send((anumg-begg)*nblocks + 3) = icols
          data_send((anumg-begg)*nblocks + 4) = ipfts
          data_send((anumg-begg)*nblocks + 5) = icohorts

       enddo

       ! Scatter: Global-to-Local
       call ScatterDataG2L(ldomain_lateral%ugrid, &
            nblocks, ndata_send, data_send, ndata_recv, data_recv)

       ! Get number of ghost quantites at all subgrid categories
       procinfo%ncells_ghost    = ldomain_lateral%ugrid%ngrid_ghost
       procinfo%ntopounits_ghost   = 0
       procinfo%nlunits_ghost   = 0
       procinfo%ncols_ghost     = 0
       procinfo%npfts_ghost     = 0
       procinfo%nCohorts_ghost  = 0

       ighost_beg = ldomain_lateral%ugrid%ngrid_local   + 1
       ighost_end = ldomain_lateral%ugrid%ngrid_ghosted

       do ighost = ighost_beg, ighost_end
          procinfo%ntopounits_ghost  = procinfo%ntopounits_ghost  + data_recv((ighost-1)*nblocks + 1)
          procinfo%nlunits_ghost  = procinfo%nlunits_ghost  + data_recv((ighost-1)*nblocks + 2)
          procinfo%ncols_ghost    = procinfo%ncols_ghost    + data_recv((ighost-1)*nblocks + 3)
          procinfo%npfts_ghost    = procinfo%npfts_ghost    + data_recv((ighost-1)*nblocks + 4)
          procinfo%ncohorts_ghost = procinfo%ncohorts_ghost + data_recv((ighost-1)*nblocks + 5)
       enddo

       ! Set 'begin' index for subgrid categories
       procinfo%begg_all        = procinfo%begg
       procinfo%begt_all        = procinfo%begt
       procinfo%begl_all        = procinfo%begl
       procinfo%begc_all        = procinfo%begc
       procinfo%begp_all        = procinfo%begp
       procinfo%begCohort_all   = procinfo%begCohort

       ! Set 'end' index for subgrid categories
       procinfo%endg_all        = procinfo%endg      + procinfo%ncells_ghost
       procinfo%endt_all        = procinfo%endt      + procinfo%ntopounits_ghost
       procinfo%endl_all        = procinfo%endl      + procinfo%nlunits_ghost
       procinfo%endc_all        = procinfo%endc      + procinfo%ncols_ghost
       procinfo%endp_all        = procinfo%endp      + procinfo%npfts_ghost
       procinfo%endCohort_all   = procinfo%endCohort + procinfo%nCohorts_ghost

#endif

    endif

  end subroutine decompInit_ghosts

end module decompInitMod
