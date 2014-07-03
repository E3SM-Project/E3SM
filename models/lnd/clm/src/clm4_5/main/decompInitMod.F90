module decompInitMod

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module provides a descomposition into a clumped data structure which can
  ! be mapped back to atmosphere physics chunks.
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use spmdMod     , only : masterproc, iam, npes, mpicom, comp_id
  use shr_sys_mod , only : shr_sys_flush
  use shr_log_mod , only : errMsg => shr_log_errMsg
  use abortutils  , only : endrun
  use clm_varctl  , only : iulog
  use mct_mod
  use decompMod
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public decompInit_lnd    ! initializes lnd grid decomposition into clumps and processors
  public decompInit_clumps ! initializes atm grid decomposition into clumps
  public decompInit_glcp   ! initializes g,l,c,p decomp info
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

    ! allocate and initialize procinfo and clumps 
    ! beg and end indices initialized for simple addition of cells later 

    allocate(procinfo%cid(clump_pproc), stat=ier)
    if (ier /= 0) then
       write(iulog,*) 'decompInit_lnd(): allocation error for procinfo%cid'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif
    procinfo%nclumps = clump_pproc
    procinfo%cid(:)  = -1
    procinfo%ncells  = 0
    procinfo%nlunits = 0
    procinfo%ncols   = 0
    procinfo%npfts   = 0
    procinfo%begg    = 1
    procinfo%begl    = 1
    procinfo%begc    = 1
    procinfo%begp    = 1
    procinfo%endg    = 0
    procinfo%endl    = 0
    procinfo%endc    = 0
    procinfo%endp    = 0

    allocate(clumps(nclumps), stat=ier)
    if (ier /= 0) then
       write(iulog,*) 'decompInit_lnd(): allocation error for clumps'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if
    clumps(:)%owner   = -1
    clumps(:)%ncells  = 0
    clumps(:)%nlunits = 0
    clumps(:)%ncols   = 0
    clumps(:)%npfts   = 0
    clumps(:)%begg    = 1
    clumps(:)%begl    = 1
    clumps(:)%begc    = 1
    clumps(:)%begp    = 1
    clumps(:)%endg    = 0
    clumps(:)%endl    = 0
    clumps(:)%endc    = 0
    clumps(:)%endp    = 0

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

    allocate(lcid(lns))
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
          if (iam >  clumps(cid)%owner) then
             procinfo%begg = procinfo%begg + 1
          endif
          if (iam >= clumps(cid)%owner) then
             procinfo%endg = procinfo%endg + 1
          endif

          !--- give gridcell to cid ---
          !--- increment the beg and end indices ---
          clumps(cid)%ncells  = clumps(cid)%ncells  + 1
          do m = 1,nclumps
             if ((clumps(m)%owner >  clumps(cid)%owner) .or. &
                 (clumps(m)%owner == clumps(cid)%owner .and. m > cid)) then
                clumps(m)%begg = clumps(m)%begg + 1
             endif
             
             if ((clumps(m)%owner >  clumps(cid)%owner) .or. &
                 (clumps(m)%owner == clumps(cid)%owner .and. m >= cid)) then
                clumps(m)%endg = clumps(m)%endg + 1
             endif
          enddo

       end if
    enddo

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
    ag = 1
    do pid = 0,npes-1
    do cid = 1,nclumps
       if (clumps(cid)%owner == pid) then
         clumpcnt(cid) = ag
         ag = ag + clumps(cid)%ncells
       endif
    enddo
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
  subroutine decompInit_clumps(lns,lni,lnj,glcmask)
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
    integer , intent(in) :: lns,lni,lnj ! land domain global size
    integer , pointer, optional   :: glcmask(:)  ! glc mask
    !
    ! !LOCAL VARIABLES:
    integer :: ln,an              ! indices
    integer :: i,g,l,k            ! indices
    integer :: cid,pid            ! indices
    integer :: n,m,np             ! indices
    integer :: anumg              ! lnd num gridcells
    integer :: begg,endg          ! beg,end gridcells
    integer :: begl,endl          ! beg,end landunits
    integer :: begc,endc          ! beg,end columns
    integer :: begp,endp          ! beg,end pfts
    integer :: icells             ! temporary
    integer :: ilunits            ! temporary
    integer :: icols              ! temporary
    integer :: ipfts              ! temporary
    integer :: ier                ! error code
    integer, allocatable :: allvecg(:,:)  ! temporary vector "global"
    integer, allocatable :: allvecl(:,:)  ! temporary vector "local"
    integer :: ntest
    character(len=32), parameter :: subname = 'decompInit_clumps'
    !------------------------------------------------------------------------------

    !--- assign gridcells to clumps (and thus pes) ---
    call get_proc_bounds(begg, endg)

    allocate(allvecl(nclumps,4))   ! local  clumps [gcells,lunit,cols,pfts]
    allocate(allvecg(nclumps,4))   ! global clumps [gcells,lunit,cols,pfts]

    ! Determine the number of gridcells, landunits, columns, and pfts 
    ! on this processor 
    ! Determine number of landunits, columns and pfts for each global
    ! gridcell index (an) that is associated with the local gridcell index (ln)

    allvecg= 0
    allvecl= 0
    do anumg = begg,endg
       an  = ldecomp%gdc2glo(anumg)
       cid = lcid(an)
       ln  = anumg
       if (present(glcmask)) then
          call subgrid_get_gcellinfo (ln, nlunits=ilunits, ncols=icols, npfts=ipfts, glcmask=glcmask(ln))
       else
          call subgrid_get_gcellinfo (ln, nlunits=ilunits, ncols=icols, npfts=ipfts)
       endif
       allvecl(cid,1) = allvecl(cid,1) + 1
       allvecl(cid,2) = allvecl(cid,2) + ilunits  ! number of landunits for local clump cid
       allvecl(cid,3) = allvecl(cid,3) + icols    ! number of columns for local clump cid
       allvecl(cid,4) = allvecl(cid,4) + ipfts    ! number of pfts for local clump cid 
    enddo
    call mpi_allreduce(allvecl,allvecg,size(allvecg),MPI_INTEGER,MPI_SUM,mpicom,ier)

    ! Determine overall  total gridcells, landunits, columns and pfts and distribute
    ! gridcells over clumps

    numg = 0
    numl = 0
    numc = 0
    nump = 0
    do cid = 1,nclumps
       icells  = allvecg(cid,1)  ! number of all clump cid gridcells (over all processors)
       ilunits = allvecg(cid,2)  ! number of all clump cid landunits (over all processors)
       icols   = allvecg(cid,3)  ! number of all clump cid columns (over all processors)
       ipfts   = allvecg(cid,4)  ! number of all clump cid pfts (over all processors)

       !--- overall total ---
       numg = numg + icells      ! total number of gridcells
       numl = numl + ilunits     ! total number of landunits
       numc = numc + icols       ! total number of columns
       nump = nump + ipfts       ! total number of pfts

       !--- give gridcell to cid ---
       !--- increment the beg and end indices ---
       clumps(cid)%nlunits = clumps(cid)%nlunits + ilunits  
       clumps(cid)%ncols   = clumps(cid)%ncols   + icols
       clumps(cid)%npfts   = clumps(cid)%npfts   + ipfts

       do m = 1,nclumps
          if ((clumps(m)%owner >  clumps(cid)%owner) .or. &
              (clumps(m)%owner == clumps(cid)%owner .and. m > cid)) then
             clumps(m)%begl = clumps(m)%begl + ilunits
             clumps(m)%begc = clumps(m)%begc + icols
             clumps(m)%begp = clumps(m)%begp + ipfts
          endif

          if ((clumps(m)%owner >  clumps(cid)%owner) .or. &
              (clumps(m)%owner == clumps(cid)%owner .and. m >= cid)) then
             clumps(m)%endl = clumps(m)%endl + ilunits
             clumps(m)%endc = clumps(m)%endc + icols
             clumps(m)%endp = clumps(m)%endp + ipfts
          endif
       enddo

       !--- give gridcell to the proc that owns the cid ---
       !--- increment the beg and end indices ---
       if (iam == clumps(cid)%owner) then
          procinfo%nlunits = procinfo%nlunits + ilunits
          procinfo%ncols   = procinfo%ncols   + icols
          procinfo%npfts   = procinfo%npfts   + ipfts
       endif

       if (iam >  clumps(cid)%owner) then
          procinfo%begl = procinfo%begl + ilunits
          procinfo%begc = procinfo%begc + icols
          procinfo%begp = procinfo%begp + ipfts
       endif

       if (iam >= clumps(cid)%owner) then
          procinfo%endl = procinfo%endl + ilunits
          procinfo%endc = procinfo%endc + icols
          procinfo%endp = procinfo%endp + ipfts
       endif
    enddo

    do n = 1,nclumps
       if (clumps(n)%ncells  /= allvecg(n,1) .or. &
           clumps(n)%nlunits /= allvecg(n,2) .or. &
           clumps(n)%ncols   /= allvecg(n,3) .or. &
           clumps(n)%npfts   /= allvecg(n,4)) then
          write(iulog,*) 'decompInit_glcp(): allvecg error ncells ',iam,n,clumps(n)%ncells ,allvecg(n,1)
          write(iulog,*) 'decompInit_glcp(): allvecg error lunits ',iam,n,clumps(n)%nlunits,allvecg(n,2)
          write(iulog,*) 'decompInit_glcp(): allvecg error ncols  ',iam,n,clumps(n)%ncols  ,allvecg(n,3)
          write(iulog,*) 'decompInit_glcp(): allvecg error pfts   ',iam,n,clumps(n)%npfts  ,allvecg(n,4)
          call endrun(msg=errMsg(__FILE__, __LINE__))
       endif
    enddo

    deallocate(allvecg,allvecl)
    deallocate(lcid)

  end subroutine decompInit_clumps

  !------------------------------------------------------------------------------
  subroutine decompInit_glcp(lns,lni,lnj,glcmask)
    !
    ! !DESCRIPTION:
    ! Determine gsmaps for landunits, columns and pfts
    !
    ! !USES:
    use clmtype   , only : lun, col, pft, grlnd
    use spmdMod
    use spmdGathScatMod
    use subgridMod, only : subgrid_get_gcellinfo
    use mct_mod
    !
    ! !ARGUMENTS:
    implicit none
    integer , intent(in) :: lns,lni,lnj ! land domain global size
    integer , pointer, optional   :: glcmask(:)  ! glc mask
    !
    ! !LOCAL VARIABLES:
    integer :: gi,li,ci,pi        ! indices
    integer :: i,g,k,l,n,np       ! indices
    integer :: cid,pid            ! indices
    integer :: begg,endg          ! beg,end gridcells
    integer :: begl,endl          ! beg,end landunits
    integer :: begc,endc          ! beg,end columns
    integer :: begp,endp          ! beg,end pfts
    integer :: numg               ! total number of gridcells across all processors
    integer :: numl               ! total number of landunits across all processors
    integer :: numc               ! total number of columns across all processors
    integer :: nump               ! total number of pfts across all processors
    integer :: icells             ! temporary
    integer :: ilunits            ! temporary
    integer :: icols              ! temporary
    integer :: ipfts              ! temporary
    integer :: ier                ! error code
    integer :: npmin,npmax,npint  ! do loop values for printing
    integer :: clmin,clmax        ! do loop values for printing
    integer :: locsize,globsize   ! used for gsmap init
    integer :: ng                 ! number of gridcells in gsmap_lnd_gdc2glo
    integer :: val1, val2         ! temporaries
    integer, pointer :: gindex(:) ! global index for gsmap init
    integer, pointer :: arrayglob(:) ! temporaroy
    integer, pointer :: gstart(:), gcount(:)
    integer, pointer :: lstart(:), lcount(:)
    integer, pointer :: cstart(:), ccount(:)
    integer, pointer :: pstart(:), pcount(:)
    integer, pointer :: ioff(:)
    integer, parameter :: dbug=1      ! 0 = min, 1=normal, 2=much, 3=max
    character(len=32), parameter :: subname = 'decompInit_glcp'
    !------------------------------------------------------------------------------

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
    call get_proc_global(numg, numl, numc, nump)

    ! Determine global seg megs

    allocate(gstart(begg:endg))
    gstart(:) = 0
    allocate(gcount(begg:endg))
    gcount(:) = 0
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
    allocate(ioff(begg:endg)) 
    ioff(:) = 0

    ! Determine gcount, lcount, ccount and pcount

    do gi = begg,endg
       if (present(glcmask)) then
          call subgrid_get_gcellinfo (gi, nlunits=ilunits, ncols=icols, npfts=ipfts, glcmask=glcmask(gi))
       else
          call subgrid_get_gcellinfo (gi, nlunits=ilunits, ncols=icols, npfts=ipfts)
       endif
       gcount(gi) = 1        ! number of gridcells for local gridcell index gi
       lcount(gi) = ilunits  ! number of landunits for local gridcell index gi
       ccount(gi) = icols    ! number of columns for local gridcell index gi
       pcount(gi) = ipfts    ! number of pfts for local gridcell index gi
    enddo

    ! Determine gstart, lstart, cstart and pstart for the OUTPUT 1d data structures

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

    deallocate(arrayglob)

    ! Gridcell gsmap (compressed, no ocean points)

    allocate(gindex(begg:endg))
    i = begg-1
    do gi = begg,endg
       if (gcount(gi) <  1) then
          write(iulog,*) 'decompInit_glcp warning count g ',k,iam,g,gcount(g)
       endif
       do l = 1,gcount(gi)
          i = i + 1
          if (i < begg .or. i > endg) then
             write(iulog,*) 'decompInit_glcp error i ',i,begg,endg
             call endrun(msg=errMsg(__FILE__, __LINE__))
          endif
          gindex(i) = gstart(gi) + l - 1
       enddo
    enddo
    if (i /= endg) then
       write(iulog,*) 'decompInit_glcp error size ',i,begg,endg
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif
    locsize = endg-begg+1
    globsize = numg
    call mct_gsMap_init(gsmap_gce_gdc2glo, gindex, mpicom, comp_id, locsize, globsize)
    deallocate(gindex)

    ! Landunit gsmap

    allocate(gindex(begl:endl))
    ioff(:) = 0
    do li = begl,endl
       gi = lun%gridcell(li) !===this is determined internally from how landunits are spread out in memory
       gindex(li) = lstart(gi) + ioff(gi) !=== the output gindex is ALWAYS the same regardless of how landuntis are spread out in memory
       ioff(gi)  = ioff(gi) + 1 
       ! check that this is less than [lstart(gi) + lcount(gi)]
    enddo
    locsize = endl-begl+1
    globsize = numl
    call mct_gsMap_init(gsmap_lun_gdc2glo, gindex, mpicom, comp_id, locsize, globsize)
    deallocate(gindex)

    ! Column gsmap

    allocate(gindex(begc:endc))
    ioff(:) = 0
    do ci = begc,endc
       gi = col%gridcell(ci)
       gindex(ci) = cstart(gi) + ioff(gi)
       ioff(gi) = ioff(gi) + 1 
       ! check that this is less than [cstart(gi) + ccount(gi)]
    enddo
    locsize = endc-begc+1
    globsize = numc
    call mct_gsMap_init(gsmap_col_gdc2glo, gindex, mpicom, comp_id, locsize, globsize)
    deallocate(gindex)

    ! PFT gsmap

    allocate(gindex(begp:endp))
    ioff(:) = 0
    do pi = begp,endp
       gi = pft%gridcell(pi)
       gindex(pi) = pstart(gi) + ioff(gi)
       ioff(gi) = ioff(gi) + 1 
       ! check that this is less than [pstart(gi) + pcount(gi)]
    enddo
    locsize = endp-begp+1
    globsize = nump
    call mct_gsMap_init(gsmap_pft_gdc2glo, gindex, mpicom, comp_id, locsize, globsize)
    deallocate(gindex)

    ! Deallocate start/count arrays
    deallocate(gstart,gcount)
    deallocate(lstart,lcount)
    deallocate(cstart,ccount)
    deallocate(pstart,pcount)
    deallocate(ioff)

    ! Diagnostic output

    if (masterproc) then
       write(iulog,*)' Surface Grid Characteristics'
       write(iulog,*)'   longitude points          = ',lni
       write(iulog,*)'   latitude points           = ',lnj
       write(iulog,*)'   total number of gridcells = ',numg
       write(iulog,*)'   total number of landunits = ',numl
       write(iulog,*)'   total number of columns   = ',numc
       write(iulog,*)'   total number of pfts      = ',nump
       write(iulog,*)' Decomposition Characteristics'
       write(iulog,*)'   clumps per process        = ',clump_pproc
       write(iulog,*)' gsMap Characteristics'
       write(iulog,*) '  lnd gsmap glo num of segs = ',mct_gsMap_ngseg(gsMap_lnd_gdc2glo)
       write(iulog,*) '  gce gsmap glo num of segs = ',mct_gsMap_ngseg(gsMap_gce_gdc2glo)
       write(iulog,*) '  lun gsmap glo num of segs = ',mct_gsMap_ngseg(gsMap_lun_gdc2glo)
       write(iulog,*) '  col gsmap glo num of segs = ',mct_gsMap_ngseg(gsMap_col_gdc2glo)
       write(iulog,*) '  pft gsmap glo num of segs = ',mct_gsMap_ngseg(gsMap_pft_gdc2glo)
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
               ' lnd ngseg   = ',mct_gsMap_ngseg(gsMap_lnd_gdc2glo), &
               ' lnd nlseg   = ',mct_gsMap_nlseg(gsMap_lnd_gdc2glo,iam)
          write(iulog,*)'proc= ',pid,&
               ' gce ngseg   = ',mct_gsMap_ngseg(gsMap_gce_gdc2glo), &
               ' gce nlseg   = ',mct_gsMap_nlseg(gsMap_gce_gdc2glo,iam)
          write(iulog,*)'proc= ',pid,&
               ' lun ngseg   = ',mct_gsMap_ngseg(gsMap_lun_gdc2glo), &
               ' lun nlseg   = ',mct_gsMap_nlseg(gsMap_lun_gdc2glo,iam)
          write(iulog,*)'proc= ',pid,&
               ' col ngseg   = ',mct_gsMap_ngseg(gsMap_col_gdc2glo), &
               ' col nlseg   = ',mct_gsMap_nlseg(gsMap_col_gdc2glo,iam)
          write(iulog,*)'proc= ',pid,&
               ' pft ngseg   = ',mct_gsMap_ngseg(gsMap_pft_gdc2glo), &
               ' pft nlseg   = ',mct_gsMap_nlseg(gsMap_pft_gdc2glo,iam)
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
          end do
       end if
       call shr_sys_flush(iulog)
       call mpi_barrier(mpicom,ier)
    end do
    call shr_sys_flush(iulog)

  end subroutine decompInit_glcp

end module decompInitMod
