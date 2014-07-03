module seq_mctext_mod

!---------------------------------------------------------------------
!
! Purpose:
!
! Shared routines for extension and computation of gsmaps, avs, and ggrids
!       
! Author: T Craig
!
!---------------------------------------------------------------------

  use shr_sys_mod
  use mct_mod
  use seq_flds_mod, only:seq_flds_dom_coord, seq_flds_dom_other 
  use seq_comm_mct
  use seq_diag_mct
  use shr_kind_mod, only: R8 => SHR_KIND_R8, CXX => SHR_KIND_CXX

  implicit none
  private  ! except
#include <mpif.h>
  save

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public  :: seq_mctext_gsmapInit
  public  :: seq_mctext_avInit
  public  :: seq_mctext_gGridInit
  public  :: seq_mctext_gsmapIdentical
  public  :: seq_mctext_gsmapExtend
  public  :: seq_mctext_avExtend
  private :: seq_mctext_gsmapCreate
  private :: seq_mctext_avCreate

!--------------------------------------------------------------------------
! Public data
!--------------------------------------------------------------------------

  integer,public :: seq_mctext_decomp

!--------------------------------------------------------------------------
! Private data
!--------------------------------------------------------------------------

  character(*),parameter :: subName = '(seq_mctext_mct)'
  real(r8),parameter :: c1 = 1.0_r8

!=======================================================================
   contains
!=======================================================================

  subroutine seq_mctext_gsmapInit( gsmap_old, ID_old, &
                                  gsmap_new, ID_new, ID_join)

    ! This routine initializes a gsmap based on another gsmap potentially
    ! on other pes.  It addresses non-overlap of pes.

    implicit none

    !-----------------------------------------------------
    ! 
    ! Arguments
    !
    type(mct_gsMap),intent(in)    :: gsmap_old
    integer        ,intent(in)    :: ID_old
    type(mct_gsMap),intent(inout) :: gsmap_new
    integer        ,intent(in)    :: ID_new
    integer        ,intent(in)    :: ID_join
    !
    ! Local Variables
    !
    character(len=*),parameter :: subname = "(seq_mctext_gsmapInit) "

    integer :: lsize
    integer :: mpicom_old
    integer :: mpicom_new
    integer :: mpicom_join
    integer :: ierr

    type(mct_gsMap)          :: gsmap_old_join     ! gsmap_old on joined id
    !-----------------------------------------------------

    call seq_comm_setptrs(ID_old ,mpicom=mpicom_old)
    call seq_comm_setptrs(ID_new ,mpicom=mpicom_new)
    call seq_comm_setptrs(ID_join,mpicom=mpicom_join)

    ! --- Set gsmaps
    ! ---   Extend the old one to now span all pes on ID_join
    ! ---   Create a new gsmap on pes associated with ID_new using info from the old one

    call seq_mctext_gsmapExtend(gsmap_old     , mpicom_old  , gsmap_old_join, mpicom_join, ID_join)
    call seq_mctext_gsmapCreate(gsmap_old_join, mpicom_join , gsmap_new     , mpicom_new , ID_new  )

    call mct_gsMap_clean(gsmap_old_join)

  end subroutine seq_mctext_gsmapInit

!=======================================================================

  subroutine seq_mctext_avInit( AV1_old, ID_old, AV1_new, ID_new, gsmap_new, ID_join)

    ! This routine initializes Avs that may need to be extended

    implicit none

    !-----------------------------------------------------
    ! 
    ! Arguments
    !
    type(mct_aVect),intent(inout) :: AV1_old
    integer        ,intent(in)    :: ID_old
    type(mct_aVect),intent(inout) :: AV1_new
    integer        ,intent(in)    :: ID_new
    type(mct_gsmap),intent(in)    :: gsmap_new
    integer        ,intent(in)    :: ID_join

    !
    ! Local Variables
    !
    character(len=*),parameter :: subname = "(seq_mctext_avInit) "

    integer :: lsize
    integer :: mpicom_new
    integer :: ierr

    !-----------------------------------------------------

    ! --- Setup data for use and make sure the old ID is ok

    call seq_comm_setptrs(ID_new ,mpicom=mpicom_new)

    ! --- Extend old avs and initialize new avs for use in the future

    lsize = 0
    if (seq_comm_iamin(ID_new)) then
       lsize = mct_gsMap_lsize(gsMap_new,mpicom_new)
    endif
    call seq_mctext_avExtend(AV1_old, ID_old, ID_join)
    call seq_mctext_avCreate(AV1_old, ID_old, AV1_new, ID_join, lsize)

  end subroutine seq_mctext_avInit

!=======================================================================

  subroutine seq_mctext_gGridInit( GG1_old, ID_old, GG1_new, ID_new, gsmap_new, ID_join)

    ! This routine initializes gGrids that may need to be extended

    implicit none

    !-----------------------------------------------------
    ! 
    ! Arguments
    !
    type(mct_gGrid),intent(inout) :: GG1_old
    integer        ,intent(in)    :: ID_old
    type(mct_gGrid),intent(inout) :: GG1_new
    integer        ,intent(in)    :: ID_new
    type(mct_gsmap),intent(in)    :: gsmap_new
    integer        ,intent(in)    :: ID_join

    !
    ! Local Variables
    !
    character(len=*),parameter :: subname = "(seq_mctext_gGridInit) "

    integer :: lsize
    integer :: mpicom_new
    integer :: ierr

    !-----------------------------------------------------

    ! --- Setup data for use and make sure the old ID is ok

    call seq_comm_setptrs(ID_new ,mpicom=mpicom_new)

    ! --- Extend old ggrids and initialize new ggrids for use in the future

    lsize = 0
    if (seq_comm_iamin(ID_new)) then
       lsize = mct_gsMap_lsize(gsMap_new,mpicom_new)
    endif
    call seq_mctext_avExtend(GG1_old%data, ID_old, ID_join)

    call mct_gGrid_init(GGrid=GG1_new, CoordChars=seq_flds_dom_coord, OtherChars=seq_flds_dom_other, lsize=lsize )
    call mct_avect_zero(GG1_new%data)

  end subroutine seq_mctext_gGridInit

!=======================================================================

  subroutine seq_mctext_gsmapExtend(gsmapi, mpicomi, gsmapo, mpicomo, compido)

    !----------------------------------------------------------------
    ! Extend/Convert a gsmap from one mpicom to another mpicom that contains
    ! at least all the pes that gsmap uses, but with different ranks
    !----------------------------------------------------------------
  
    implicit none
    type(mct_gsMap), intent(IN) :: gsmapi
    integer        , intent(IN) :: mpicomi
    type(mct_gsMap), intent(OUT):: gsmapo
    integer        , intent(IN) :: mpicomo
    integer        , intent(IN) :: compido

    character(len=*),parameter :: subname = "(seq_mctext_gsmapExtend) "
    integer :: n
    integer :: ngseg
    integer :: gsize
    integer :: msizei,msizeo
    integer :: mrank,mranko,mrankog   ! sets pe rank of root mpicomi pe in mpicomo
    integer :: mpigrpi,mpigrpo
    integer :: ierr
    integer, pointer :: pei(:),peo(:)
    integer, pointer :: start(:),length(:),peloc(:)

    mranko = -1

    ! --- create the new gsmap on the mpicomi root only

    if (mpicomi /= MPI_COMM_NULL) then
       call mpi_comm_rank(mpicomi,mrank,ierr)
       call shr_mpi_chkerr(ierr,subname//' gsm_cop mpi_comm_rank i')
       if (mrank == 0) then
          call mpi_comm_group(mpicomi,mpigrpi,ierr)
          call shr_mpi_chkerr(ierr,subname//' gsm_cop mpi_comm_group i')
          call mpi_comm_group(mpicomo,mpigrpo,ierr)
          call shr_mpi_chkerr(ierr,subname//' gsm_cop mpi_comm_group o')
          call mpi_comm_size(mpicomi,msizei,ierr)
          call shr_mpi_chkerr(ierr,subname//' gsm_cop mpi_comm_size i')
          call mpi_comm_size(mpicomo,msizeo,ierr)
          call shr_mpi_chkerr(ierr,subname//' gsm_cop mpi_comm_size o')

          ! --- setup the translation of pe numbers from the old gsmap(mpicom)
          ! --- to the new one, pei -> peo

          allocate(pei(0:msizei-1),peo(0:msizei-1))
          do n = 0,msizei-1
             pei(n) = n
          enddo

          peo = -1
          call mpi_group_translate_ranks(mpigrpi,msizei,pei,mpigrpo,peo,ierr)
          call shr_mpi_chkerr(ierr,subname//' gsm_cop mpi_group_translate_ranks')

          do n = 0,msizei-1
             if (peo(n) < 0 .or. peo(n) > msizeo-1) then
                write(logunit,*) subname,' peo out of bounds ',peo(n),msizeo
                call shr_sys_abort()
             endif
          enddo

          mranko = peo(0)

          ! --- compute the new gsmap which has the same start and length values
          ! --- but peloc is now the mapping of pei to peo

          ngseg = gsmapi%ngseg
          gsize = gsmapi%gsize
          allocate(start(ngseg),length(ngseg),peloc(ngseg))
          do n = 1,ngseg
             start(n)  = gsmapi%start(n)
             length(n) = gsmapi%length(n)
             peloc(n)  = peo(gsmapi%pe_loc(n))
          enddo

          ! --- initialize the gsmap on the root pe

          call mct_gsmap_init(gsmapo,compido,ngseg,gsize,start,length,peloc)

          deallocate(pei,peo,start,length,peloc)
       endif
    endif

    ! --- broadcast via allreduce the mpicomi root pe in mpicomo space
    ! --- mranko is -1 except on the root pe where is it peo of that pe

    call mpi_allreduce(mranko,mrankog,1,MPI_INTEGER,MPI_MAX,mpicomo,ierr)
    call shr_mpi_chkerr(ierr,subname//' gsm_cop mpi_allreduce max')

    ! --- broadcast the gsmap to all pes in mpicomo from mrankog

    call mct_gsmap_bcast(gsmapo, mrankog, mpicomo)

! tcx summarize decomp info
#if (1 == 0)
    write(logunit,*) trim(subname),'tcxa ',mpicomi,mpicomo
    call shr_sys_flush(logunit)
    call mpi_barrier(mpicomo,ierr)

    if (mpicomi /= MPI_COMM_NULL) then
       call mpi_comm_rank(mpicomi,mrank,ierr)
       write(logunit,*) 'tcxbi ',mrank
       if (mrank == 0) then
          write(logunit,*) 'tcxci ',gsmapi%ngseg,size(gsmapi%start),gsmapi%gsize,gsmapi%comp_id
          do n = 1,gsmapi%ngseg
             write(logunit,*) 'tcx gsmti ',n,gsmapi%start(n),gsmapi%length(n),gsmapi%pe_loc(n)
          enddo
          call shr_sys_flush(logunit)
      endif
    endif

    if (mpicomo /= MPI_COMM_NULL) then
       call mpi_comm_rank(mpicomo,mrank,ierr)
       write(logunit,*) 'tcxbo ',mrank
       if (mrank == 0) then
          write(logunit,*) 'tcxco ',gsmapo%ngseg,size(gsmapo%start),gsmapo%gsize,gsmapo%comp_id
          do n = 1,gsmapo%ngseg
             write(logunit,*) 'tcx gsmto ',n,gsmapo%start(n),gsmapo%length(n),gsmapo%pe_loc(n)
          enddo
          call shr_sys_flush(logunit)
       endif
    endif

    call shr_sys_flush(logunit)
    call mpi_barrier(mpicomo,ierr)
#endif


  end subroutine seq_mctext_gsmapExtend
!=======================================================================

  subroutine seq_mctext_gsmapCreate(gsmapi, mpicomi, gsmapo, mpicomo, compido)

    !---------------------------------------------------------------------
    ! creates a new gsmap on a subset of pes, requires setting a new decomp
    !---------------------------------------------------------------------
  
    implicit none
    type(mct_gsMap), intent(IN) :: gsmapi
    integer        , intent(IN) :: mpicomi
    type(mct_gsMap), intent(OUT):: gsmapo
    integer        , intent(IN) :: mpicomo
    integer        , intent(IN) :: compido

    character(len=*),parameter :: subname = "(seq_mctext_gsmapCreate) "
    integer :: n,m,k
    integer :: ktot            ! number of active cells in gsmap
    integer :: apesi, apeso    ! number of active pes in gsmap
    integer ::        lsizeo   ! local size for lindex
    integer :: ngsegi,ngsego   ! ngseg of mpicomi, mpicomo
    integer :: gsizei,gsizeo   ! gsize of mpicomi, mpicomo
    integer :: msizei,msizeo   ! size of mpicomi, mpicomo
    integer :: mranki,mranko   ! rank in mpicomi, mpicomo
    integer :: ierr
    integer :: decomp_type
    integer, pointer :: start(:),length(:),peloc(:),perm(:),gindex(:),lindex(:)
    real(r8):: rpeloc
    logical :: gsmap_bfbflag = .false. ! normally this should be set to false

    ! --- create a new gsmap on new pes based on the old gsmap
    ! --- gsmapi must be known on all mpicomo pes, compute the same 
    ! --- thing on all pes in parallel

    if (mpicomo /= MPI_COMM_NULL) then
       call mpi_comm_rank(mpicomi,mranki,ierr)
       call shr_mpi_chkerr(ierr,subname//' mpi_comm_rank i')
       call mpi_comm_size(mpicomi,msizei,ierr)
       call shr_mpi_chkerr(ierr,subname//' mpi_comm_size i')
       call mpi_comm_rank(mpicomo,mranko,ierr)
       call shr_mpi_chkerr(ierr,subname//' mpi_comm_rank o')
       call mpi_comm_size(mpicomo,msizeo,ierr)
       call shr_mpi_chkerr(ierr,subname//' mpi_comm_size o')

       ngsegi = gsmapi%ngseg
       gsizei = gsmapi%gsize
       gsizeo = gsizei
       call mct_gsMap_activepes(gsmapi,apesi)

       decomp_type = 0

       if (seq_mctext_decomp == 0) then
          if (msizeo == apesi) then      ! preserve segments and decomp
             ! For testing - set decomp_type to 1 - to have gsmapi and gsmapo identical
             if (gsmap_bfbflag) then
                decomp_type = 1     ! better in cpl to have all decomps "same-ish"
             else
                decomp_type = 2
             end if
          elseif (ngsegi >= msizeo) then ! preserve segments, new decomp
             decomp_type = 2
          else                           ! new segments
             decomp_type = 3
          endif
       else
          decomp_type = seq_mctext_decomp
       endif

!tcx       decomp_type = 3 ! over ride setting above for testing
!       if (mranko == 0) write(logunit,'(2A,4I)') trim(subname),' decomp_type =',decomp_type,ngsegi,msizeo,apesi

       select case (decomp_type)

       case(1)   ! --- preserve segments and decomp ---------------------

          ! -- copy the gsmap and translate the pes
          call mct_gsMap_copy(gsmapi,gsmapo)
          ngsego = ngsegi
          do n = 1,ngsego
             gsmapo%pe_loc(n) = mod(gsmapo%pe_loc(n),msizeo)    ! translate pes 1:1 from old to new
          enddo

       case(2)   ! --- preserve segments, new decomp --------------------

          ! --- preserve segments, sort the start and length, assign a new pe list
          ngsego = ngsegi
          allocate(start(ngsego),length(ngsego),peloc(ngsego),perm(ngsego))
          do n = 1,ngsego
             start(n)  = gsmapi%start(n)
             length(n) = gsmapi%length(n)
          enddo
          ! --- sort gsmap to minimize permute cost in mct
          call mct_indexset(perm)
          call mct_indexsort(ngsego,perm,start)
          call mct_permute(start,perm,ngsego)
          call mct_permute(length,perm,ngsego)
          ! --- give each pe "equal" number of segments, use reals to avoid integer overflow
          do n = 1,ngsego
             rpeloc = (((msizeo*c1)*((n-1)*c1))/(ngsego*c1))      ! give each pe "equal" number of segments, use reals to avoid integer overflow
             peloc(n) = int(rpeloc)
          enddo
          call mct_gsmap_init(gsmapo,ngsego,start,length,peloc,0,mpicomo,compido,gsizeo)
          deallocate(start,length,peloc,perm)

       case(3)   ! --- new segments, new decomp -------------------------

          ! --- new segments, compute gindex, then parse the gridcells out evenly

          k = 0
          do n = 1,ngsegi
          do m = 1,gsmapi%length(n)
             k = k + 1
             if (k > gsizei) then
                write(logunit,*) trim(subname),' ERROR in gindex ',k,gsizei
                call shr_sys_abort()
             endif
          enddo
          enddo
          ktot = k

          allocate(gindex(ktot),perm(ktot))  

          k = 0
          do n = 1,ngsegi
          do m = 1,gsmapi%length(n)
             k = k + 1
             gindex(k) = gsmapi%start(n) + m - 1
          enddo
          enddo
          call mct_indexset(perm)
          call mct_indexsort(ktot,perm,gindex)
          call mct_permute(gindex,perm,ktot)

          k = 0
          do m = 0,msizeo-1
             lsizeo = ktot/msizeo
             if (m < (ktot - lsizeo*msizeo)) lsizeo = lsizeo + 1
             if (mranko == m) then
                allocate(lindex(lsizeo))
                if (k+lsizeo > ktot) then
                   write(logunit,*) trim(subname),' ERROR: decomp out of bounds ',mranko,k,lsizeo,ktot
                   call shr_sys_abort()
                endif
                lindex(1:lsizeo) = gindex(k+1:k+lsizeo)
!                write(logunit,*) trim(subname),' decomp is ',mranko,lsizeo,k+1,k+lsizeo
             endif
             k = k + lsizeo
          enddo
          if (k /= ktot) then
             write(logunit,*) trim(subname),' ERROR: decomp incomplete ',k,ktot
             call shr_sys_abort()
          endif

          call mct_gsmap_init(gsmapo,lindex,mpicomo,compido,size(lindex),gsizeo)
          deallocate(gindex,perm,lindex)

       case default   ! --- unknown ---
          write(logunit,*) trim(subname),' ERROR decomp_type unknown ',decomp_type
          call shr_sys_abort(trim(subname)//' ERROR decomp_type unknown')

       end select

       if (mranko == 0) then
          write(logunit,102) trim(subname),' created new gsmap decomp_type =',decomp_type
          write(logunit,102) trim(subname),'   ngseg/gsize        = ', &
             mct_gsmap_ngseg(gsmapo),mct_gsmap_gsize(gsmapo)
          call mct_gsmap_activepes(gsmapo,apeso)
          write(logunit,102) trim(subname),'   mpisize/active_pes = ', &
             msizeo,apeso
          write(logunit,102) trim(subname),'   avg seg per pe/ape = ', &
             mct_gsmap_ngseg(gsmapo)/msizeo,mct_gsmap_ngseg(gsmapo)/apeso
          write(logunit,102) trim(subname),'   nlseg/maxnlsegs    = ', &
             mct_gsmap_nlseg(gsmapo,0),mct_gsmap_maxnlseg(gsmapo)
 102      format(2A,2I8)
       endif

!p       if (.not. mct_gsmap_increasing(gsmapo) ) then
!          write(logunit,*) trim(subname),' ERROR: gsmapo not increasing'
!          call shr_sys_abort()
!       endif

    endif

  end subroutine seq_mctext_gsmapCreate
    
!=======================================================================

subroutine seq_mctext_avExtend(AVin,IDin,ID)

  !-----------------------------------------------------------------------
  ! Extend an AV to a larger set of pes or
  ! Initialize an AV on another set of pes
  !-----------------------------------------------------------------------

  implicit none
  type(mct_aVect), intent(INOUT):: AVin
  integer         ,intent(IN)   :: IDin ! ID associated with AVin
  integer        , intent(IN)   :: ID   ! ID to initialize over

  ! Local variables

  character(len=*),parameter :: subname = "(seq_mctext_avExtend) "
  integer :: mpicom
  integer :: rank,rank2
  integer :: lsizei, lsizen
  integer :: srank,srankg
  integer :: ierr
  integer :: nints
  character(len=CXX) :: iList,rList


  call seq_comm_setptrs(ID,mpicom=mpicom,iam=rank)

  ! --- lsizen is the size of the newly initialized AV, zero is valid
  ! --- lsizei is -1 on any peszero on any pes where AV is not yet initialized

  lsizei = -1  
  if (seq_comm_iamin(IDin)) lsizei = mct_aVect_lsize(AVin)
  lsizen = 0

  ! --- find a pe that already has AVin allocated, use MPI_MAX to do so
  ! --- set the pe and broadcast it to all other pes using mpi_allreduce

  srank = -1
  srankg = -1
  if (lsizei > 0) srank = rank

  call mpi_allreduce(srank,srankg,1,MPI_INTEGER,MPI_MAX,mpicom,ierr)
  call shr_mpi_chkerr(ierr,subname//' mpi_allreduce max')

  if (srankg < 0) then
    write(logunit,*) subname,' WARNING AVin empty '
    return
  endif

  ! --- set the iList and rList from the broadcast pe (srankg) and 
  ! --- broadcast the lists

  iList = " "
  rList = " "
  if (rank == srankg) then
    if (mct_aVect_nIAttr(AVin) /= 0) iList = mct_aVect_ExportIList2c(AVin)
    if (mct_aVect_nRattr(AVin) /= 0) rList = mct_aVect_ExportRList2c(AVin)
  endif

  call mpi_bcast(iList,len(iList),MPI_CHARACTER,srankg,mpicom,ierr)
  call mpi_bcast(rList,len(rList),MPI_CHARACTER,srankg,mpicom,ierr)

  ! --- now allocate the AV on any pes where the orig size is zero.  those
  ! --- should be pes that either have no data and may have been allocated
  ! --- before (no harm in doing it again) or have never been allocated

  if (lsizei <= 0) then
    if(len_trim(iList) > 0 .and. len_trim(rList) > 0) then
      call mct_aVect_init(AVin,iList=iList,rList=rList,lsize=lsizen)
    elseif (len_trim(iList) > 0 .and. len_trim(rList) == 0) then
      call mct_aVect_init(AVin,iList=iList,lsize=lsizen)
    elseif (len_trim(iList) == 0 .and. len_trim(rList) > 0) then
      call mct_aVect_init(AVin,rList=rList,lsize=lsizen)
    endif
  endif

end subroutine seq_mctext_avExtend

!=======================================================================

subroutine seq_mctext_avCreate(AVin,IDin,AVout,ID,lsize)

  !-----------------------------------------------------------------------
  ! Extend an AV to a larger set of pes or
  ! Initialize an AV on another set of pes
  !-----------------------------------------------------------------------

  implicit none
  type(mct_aVect), intent(INOUT):: AVin
  integer         ,intent(IN)   :: IDin ! ID associated with AVin
  type(mct_aVect), intent(INOUT):: AVout
  integer        , intent(IN)   :: ID   ! ID to initialize over
  integer        , intent(IN)   :: lsize

  ! Local variables

  character(len=*),parameter :: subname = "(seq_mctext_avCreate) "
  integer :: mpicom
  integer :: rank,rank2
  integer :: lsizei, lsizen
  integer :: srank,srankg
  integer :: ierr
  integer :: nints
  character(len=CXX) :: iList,rList

  call seq_comm_setptrs(ID,mpicom=mpicom,iam=rank)

  ! --- lsizen is the size of the newly initialized AV, zero is valid

  lsizei = -1  
  if (seq_comm_iamin(IDin)) lsizei = mct_aVect_lsize(AVin)
  lsizen = lsize

  ! --- find a pe that already has AVin allocated, use MPI_MAX to do so
  ! --- set the pe and broadcast it to all other pes

  srank = -1
  srankg = -1
  if (lsizei > 0) srank = rank

  call mpi_allreduce(srank,srankg,1,MPI_INTEGER,MPI_MAX,mpicom,ierr)
  call shr_mpi_chkerr(ierr,subname//' mpi_allreduce max')

  if (srankg < 0) then
    write(logunit,*) subname,' ERROR AVin not initialized '
    call shr_sys_abort()
  endif

  ! --- set the iList and rList from the broadcast pe (srankg) and 
  ! --- broadcast the lists

  iList = " "
  rList = " "
  if (rank == srankg) then
    if (mct_aVect_nIAttr(AVin) /= 0) iList = mct_aVect_ExportIList2c(AVin)
    if (mct_aVect_nRattr(AVin) /= 0) rList = mct_aVect_ExportRList2c(AVin)
  endif

  call mpi_bcast(iList,len(iList),MPI_CHARACTER,srankg,mpicom,ierr)
  call mpi_bcast(rList,len(rList),MPI_CHARACTER,srankg,mpicom,ierr)

  ! --- now allocate the AV on all pes.  the AV should not exist before.
  ! --- If it does, mct should die.

  if(len_trim(iList) > 0 .and. len_trim(rList) > 0) then
    call mct_aVect_init(AVout,iList=iList,rList=rList,lsize=lsizen)
  elseif (len_trim(iList) > 0 .and. len_trim(rList) == 0) then
    call mct_aVect_init(AVout,iList=iList,lsize=lsizen)
  elseif (len_trim(iList) == 0 .and. len_trim(rList) > 0) then
    call mct_aVect_init(AVout,rList=rList,lsize=lsizen)
  endif

end subroutine seq_mctext_avCreate

!=======================================================================
logical function seq_mctext_gsmapIdentical(gsmap1,gsmap2)

  implicit none
  type(mct_gsMap), intent(IN):: gsmap1
  type(mct_gsMap), intent(IN):: gsmap2

  ! Local variables

  character(len=*),parameter :: subname = "(seq_mctext_gsmapIdentical) "
  integer :: n
  logical :: identical

  !-----------------------

  identical = .true.

  ! --- continue compare ---
  if (identical) then
     if (mct_gsMap_gsize(gsmap1) /= mct_gsMap_gsize(gsmap2)) identical = .false.
     if (mct_gsMap_ngseg(gsmap1) /= mct_gsMap_ngseg(gsmap2)) identical = .false.
  endif

  ! --- continue compare ---
  if (identical) then
     do n = 1,mct_gsMap_ngseg(gsmap1)
        if (gsmap1%start(n)  /= gsmap2%start(n) ) identical = .false.
        if (gsmap1%length(n) /= gsmap2%length(n)) identical = .false.
        if (gsmap1%pe_loc(n) /= gsmap2%pe_loc(n)) identical = .false.
     enddo
  endif

  seq_mctext_gsmapIdentical = identical

end function seq_mctext_gsmapIdentical
    
!=======================================================================

end module seq_mctext_mod
