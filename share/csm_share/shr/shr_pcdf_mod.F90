!===============================================================================
! SVN $Id: shr_pcdf_mod.F90 18683 2009-09-30 22:20:22Z kauff $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/drv/seq_mct/trunk_tags/drvseq3_0_36/driver/shr_pcdf_mod.F90 $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: shr_pcdf_mod -- generic pio file reader and writer
!
! !DESCRIPTION:
!
!    Reads & writes pio files
!
! !REMARKS:
!
!    supports aVect, 1d real and integer, and scalar real and integer fields
!    using a common decomp for all fields.  this is a heavily overloaded interface
!    that supports read and write of multiple fields/type to a file using a single call.
!
! !REVISION HISTORY:
!     2009-Oct-15 - T. Craig - initial implementation
!
! !INTERFACE: ------------------------------------------------------------------

module shr_pcdf_mod

  use shr_kind_mod,      only: R8 => SHR_KIND_R8, IN => SHR_KIND_IN
  use shr_kind_mod,      only: CL => SHR_KIND_CL, CS => SHR_KIND_CS
  use shr_sys_mod,       only: shr_sys_abort, shr_sys_flush
  use shr_const_mod,     only: shr_const_spval
  use shr_log_mod,       only: shr_log_unit, shr_log_level
  use mct_mod
  use pio
 
  implicit none

  private

 !PUBLIC TYPES:

  ! no public types

!!PUBLIC MEMBER FUNCTIONS

  public :: shr_pcdf_readwrite

!!PUBLIC DATA MEMBERS:

  ! no public data

!EOP

  character(len=*),parameter :: version   = 'shr_pcdf_v0_0_01'
  real(r8)        ,parameter :: fillvalue = SHR_CONST_SPVAL
  integer(in)     ,parameter :: ifillvalue = -999999
  
!===============================================================================
contains
!===============================================================================
subroutine shr_pcdf_readwrite(type,iosystem,pio_iotype,filename,mpicom,gsmap,dof,clobber,cdf64, &
                     id1,id1n,rs1,rs1n,is1,is1n,rf1,rf1n,if1,if1n,av1,av1n, &
                     id2,id2n,rs2,rs2n,is2,is2n,rf2,rf2n,if2,if2n,av2,av2n, &
                     id3,id3n,rs3,rs3n,is3,is3n,rf3,rf3n,if3,if3n,av3,av3n, &
                     id4,id4n,rs4,rs4n,is4,is4n,rf4,rf4n,if4,if4n,av4,av4n  )
  use pio, only : iosystem_desc_t
  implicit none

  character(len=*) , intent(in)    :: type      ! 'read' or 'write'
  type(iosystem_desc_t), intent(inout), target :: iosystem
  integer(IN), intent(in) :: pio_iotype
  character(len=*) , intent(in)    :: filename  ! filename
  integer(IN)      , intent(in)    :: mpicom    ! mpicom

  !--- one of these must be set ---
  type(mct_gsmap)  , optional, intent(in)    :: gsmap    ! decomp for all data
  integer(IN)      , optional, intent(in)    :: dof(:)   ! decomp for all data

  !--- optional settings ---
  logical          , optional, intent(in)    :: clobber
  logical          , optional, intent(in)    :: cdf64
  ! add root, stride, ntasks, netcdf/pnetcdf, etc

  !--- data to write ---

  !--- single scalar dimensions, assumed valid on the io root pe ---
  integer(IN)      , optional, intent(inout) :: id1      ! int field 1
  character(len=*) , optional, intent(in)    :: id1n     ! if1 name
  integer(IN)      , optional, intent(inout) :: id2      ! int field 2
  character(len=*) , optional, intent(in)    :: id2n     ! if2 name
  integer(IN)      , optional, intent(inout) :: id3      ! int field 3
  character(len=*) , optional, intent(in)    :: id3n     ! if3 name
  integer(IN)      , optional, intent(inout) :: id4      ! int field 4
  character(len=*) , optional, intent(in)    :: id4n     ! if4 name

  !--- single scalar variables, assumed valid on the io root pe ---
  real(R8)         , optional, intent(inout) :: rs1      ! real field 1
  character(len=*) , optional, intent(in)    :: rs1n     ! rf1 name
  real(R8)         , optional, intent(inout) :: rs2      ! real field 2
  character(len=*) , optional, intent(in)    :: rs2n     ! rf2 name
  real(R8)         , optional, intent(inout) :: rs3      ! real field 3
  character(len=*) , optional, intent(in)    :: rs3n     ! rf3 name
  real(R8)         , optional, intent(inout) :: rs4      ! real field 4
  character(len=*) , optional, intent(in)    :: rs4n     ! rf4 name
  integer(IN)      , optional, intent(inout) :: is1      ! int field 1
  character(len=*) , optional, intent(in)    :: is1n     ! if1 name
  integer(IN)      , optional, intent(inout) :: is2      ! int field 2
  character(len=*) , optional, intent(in)    :: is2n     ! if2 name
  integer(IN)      , optional, intent(inout) :: is3      ! int field 3
  character(len=*) , optional, intent(in)    :: is3n     ! if3 name
  integer(IN)      , optional, intent(inout) :: is4      ! int field 4
  character(len=*) , optional, intent(in)    :: is4n     ! if4 name

  !--- single field, decomposed f90 data in 1d arrays ---
  real(R8)         , optional, intent(inout) :: rf1(:)   ! real field 1
  character(len=*) , optional, intent(in)    :: rf1n     ! rf1 name
  real(R8)         , optional, intent(inout) :: rf2(:)   ! real field 2
  character(len=*) , optional, intent(in)    :: rf2n     ! rf2 name
  real(R8)         , optional, intent(inout) :: rf3(:)   ! real field 3
  character(len=*) , optional, intent(in)    :: rf3n     ! rf3 name
  real(R8)         , optional, intent(inout) :: rf4(:)   ! real field 4
  character(len=*) , optional, intent(in)    :: rf4n     ! rf4 name
  integer(IN)      , optional, intent(inout) :: if1(:)   ! int field 1
  character(len=*) , optional, intent(in)    :: if1n     ! if1 name
  integer(IN)      , optional, intent(inout) :: if2(:)   ! int field 2
  character(len=*) , optional, intent(in)    :: if2n     ! if2 name
  integer(IN)      , optional, intent(inout) :: if3(:)   ! int field 3
  character(len=*) , optional, intent(in)    :: if3n     ! if3 name
  integer(IN)      , optional, intent(inout) :: if4(:)   ! int field 4
  character(len=*) , optional, intent(in)    :: if4n     ! if4 name

  !--- attr vect, decomposed f90 data in av datatype ---
  type(mct_aVect)  , optional, intent(inout) :: av1      ! avect 1
  character(len=*) , optional, intent(in)    :: av1n     ! av1 name
  type(mct_aVect)  , optional, intent(inout) :: av2      ! avect 2
  character(len=*) , optional, intent(in)    :: av2n     ! av2 name
  type(mct_aVect)  , optional, intent(inout) :: av3      ! avect 3
  character(len=*) , optional, intent(in)    :: av3n     ! av3 name
  type(mct_aVect)  , optional, intent(inout) :: av4      ! avect 4
  character(len=*) , optional, intent(in)    :: av4n     ! av4 name

  !--- local ---
  integer(IN)   :: iam,ntasks
  integer(IN)   :: ier,rcode
  integer(IN)   :: loop,minloop,maxloop
  integer(IN)   :: n,nf
  logical       :: readtype
  integer(IN)   :: lsize,gsize
  logical       :: lclobber
  logical       :: lcdf64
  logical       :: exists
  integer       :: nmode
  character(CL) :: fname
  character(CL) :: vname
  type(mct_string) :: mstring     ! mct char type
  integer(IN)   :: dimid1(1)


  type(file_desc_t)     :: fid
  type(var_desc_t)      :: varid
  type(io_desc_t)       :: iodescd
  type(io_desc_t)       :: iodesci
  integer(IN), pointer  :: ldof(:)

  character(len=*),parameter :: subname = '(shr_pcdf_readwrite) '
  
  !-------------

  if (trim(type) == 'read') then
     readtype = .true.
  elseif (trim(type) == 'write') then
     readtype = .false.
  else
     call shr_sys_abort(subname//' ERROR: read write type invalid')
  endif

  lclobber = .false.
  if (present(clobber)) lclobber=clobber

  lcdf64 = .false.
  if (present(cdf64)) lcdf64=cdf64

  call mpi_comm_size(mpicom,ntasks,ier)
  call mpi_comm_rank(mpicom,iam,ier)

  if (iam == 0) then
     write(shr_log_unit,*) subname,' filename   = ',trim(filename)
     write(shr_log_unit,*) subname,' type       = ',trim(type)
     write(shr_log_unit,*) subname,' clobber    = ',lclobber
     write(shr_log_unit,*) subname,' cdf64      = ',lcdf64
     call shr_sys_flush(shr_log_unit)
  endif

  if (present(gsmap) .and. present(dof)) then
     call shr_sys_abort(trim(subname)//' ERROR: either gsmap OR dof must be an argument')
  endif
  if (present(gsmap)) then
     lsize = mct_gsmap_lsize(gsmap,mpicom)
     gsize = mct_gsmap_gsize(gsmap)
     call mct_gsmap_OrderedPoints(gsmap, iam, ldof)
     call pio_initdecomp(iosystem, pio_double, (/gsize/), ldof, iodescd)
     call pio_initdecomp(iosystem, pio_int   , (/gsize/), ldof, iodesci)
     deallocate(ldof)
  elseif (present(dof)) then
     lsize = size(dof)
     call shr_mpi_sum(lsize,gsize,mpicom,string=trim(subname),all=.true.)
     call pio_initdecomp(iosystem, pio_double, (/gsize/), ldof, iodescd)
     call pio_initdecomp(iosystem, pio_int   , (/gsize/), ldof, iodesci)
  else
     call shr_sys_abort(trim(subname)//' ERROR: either gsmap OR dof must be an argument')
  endif

  if (iam == 0) then
     if (len_trim(filename) == 0) then
        call shr_sys_abort(trim(subname)//' ERROR: filename is empty')
     endif
     inquire(file=trim(filename),exist=exists)
  endif
  call shr_mpi_bcast(exists,mpicom,trim(subname)//' exists')

  if (readtype) then
     if (.not.exists) then
        call shr_sys_abort(trim(subname)//' ERROR: '//trim(filename)//' doesnt exist')
     endif
     nmode = pio_nowrite
     rcode = pio_openfile(iosystem, fid, pio_iotype, trim(filename), nmode)
  else
     if (.not.lclobber .and. exists) then
        call shr_sys_abort(trim(subname)//' ERROR: '//trim(filename)//' exists, no clobber set')
     endif
     if (lclobber .or. .not.exists) then
        nmode = pio_clobber
        if (lcdf64) nmode = ior(nmode,PIO_64BIT_OFFSET)
        rcode = pio_createfile(iosystem, fid, pio_iotype, trim(filename), nmode)
     else
        nmode = pio_write
        if (lcdf64) nmode = ior(nmode,PIO_64BIT_OFFSET)
        rcode = pio_openfile(iosystem, fid, pio_iotype, trim(filename), nmode)
     endif
     rcode = pio_put_att(fid,pio_global,"file_version",version)
  endif
  call pio_seterrorhandling(fid,PIO_INTERNAL_ERROR)

  if (readtype) then
     minloop = 11
     maxloop = 11
  else
     minloop = 21
     maxloop = 22
  endif

  ! loop = 11 is read
  ! loop = 21 is define
  ! loop = 22 is write
  do loop = minloop,maxloop

     if (loop == 21) rcode = pio_def_dim(fid,'gsize',gsize,dimid1(1))

     if (present(id1)) then
        fname = 'id1'
        if (present(id1n)) fname = trim(id1n)
        if (loop == 11) call shr_pcdf_readdim(fid,trim(fname),id1)
        if (loop == 21) call shr_pcdf_writedim(fid,trim(fname),id1)
     endif

     if (present(id2)) then
        fname = 'id2'
        if (present(id2n)) fname = trim(id2n)
        if (loop == 11) call shr_pcdf_readdim(fid,trim(fname),id2)
        if (loop == 21) call shr_pcdf_writedim(fid,trim(fname),id2)
     endif

     if (present(id3)) then
        fname = 'id3'
        if (present(id3n)) fname = trim(id3n)
        if (loop == 11) call shr_pcdf_readdim(fid,trim(fname),id3)
        if (loop == 21) call shr_pcdf_writedim(fid,trim(fname),id3)
     endif

     if (present(id4)) then
        fname = 'id4'
        if (present(id4n)) fname = trim(id4n)
        if (loop == 11) call shr_pcdf_readdim(fid,trim(fname),id4)
        if (loop == 21) call shr_pcdf_writedim(fid,trim(fname),id4)
     endif

     if (present(rs1)) then
        fname = 'rs1'
        if (present(rs1n)) fname = trim(rs1n)
        if (loop == 11) call shr_pcdf_readr0d(fid,trim(fname),rs1)
        if (loop == 21) call shr_pcdf_defvar0d(fid,trim(fname),PIO_DOUBLE)
        if (loop == 22) call shr_pcdf_writer0d(fid,trim(fname),rs1)
     endif

     if (present(rs2)) then
        fname = 'rs2'
        if (present(rs2n)) fname = trim(rs2n)
        if (loop == 11) call shr_pcdf_readr0d(fid,trim(fname),rs2)
        if (loop == 21) call shr_pcdf_defvar0d(fid,trim(fname),PIO_DOUBLE)
        if (loop == 22) call shr_pcdf_writer0d(fid,trim(fname),rs2)
     endif

     if (present(rs3)) then
        fname = 'rs3'
        if (present(rs3n)) fname = trim(rs3n)
        if (loop == 11) call shr_pcdf_readr0d(fid,trim(fname),rs3)
        if (loop == 21) call shr_pcdf_defvar0d(fid,trim(fname),PIO_DOUBLE)
        if (loop == 22) call shr_pcdf_writer0d(fid,trim(fname),rs3)
     endif

     if (present(rs4)) then
        fname = 'rs4'
        if (present(rs4n)) fname = trim(rs4n)
        if (loop == 11) call shr_pcdf_readr0d(fid,trim(fname),rs4)
        if (loop == 21) call shr_pcdf_defvar0d(fid,trim(fname),PIO_DOUBLE)
        if (loop == 22) call shr_pcdf_writer0d(fid,trim(fname),rs4)
     endif

     if (present(is1)) then
        fname = 'is1'
        if (present(is1n)) fname = trim(is1n)
        if (loop == 11) call shr_pcdf_readi0d(fid,trim(fname),is1)
        if (loop == 21) call shr_pcdf_defvar0d(fid,trim(fname),PIO_INT)
        if (loop == 22) call shr_pcdf_writei0d(fid,trim(fname),is1)
     endif

     if (present(is2)) then
        fname = 'is2'
        if (present(is2n)) fname = trim(is2n)
        if (loop == 11) call shr_pcdf_readi0d(fid,trim(fname),is2)
        if (loop == 21) call shr_pcdf_defvar0d(fid,trim(fname),PIO_INT)
        if (loop == 22) call shr_pcdf_writei0d(fid,trim(fname),is2)
     endif
 
     if (present(is3)) then
        fname = 'is3'
        if (present(is3n)) fname = trim(is3n)
        if (loop == 11) call shr_pcdf_readi0d(fid,trim(fname),is3)
        if (loop == 21) call shr_pcdf_defvar0d(fid,trim(fname),PIO_INT)
        if (loop == 22) call shr_pcdf_writei0d(fid,trim(fname),is3)
     endif

     if (present(is4)) then
        fname = 'is4'
        if (present(is4n)) fname = trim(is4n)
        if (loop == 11) call shr_pcdf_readi0d(fid,trim(fname),is4)
        if (loop == 21) call shr_pcdf_defvar0d(fid,trim(fname),PIO_INT)
        if (loop == 22) call shr_pcdf_writei0d(fid,trim(fname),is4)
     endif

    if (present(rf1)) then
        fname = 'rf1'
        if (present(rf1n)) fname = trim(rf1n)
        if (loop == 11) call shr_pcdf_readr1d(fid,trim(fname),iodescd,rf1)
        if (loop == 21) call shr_pcdf_defvar1d(fid,trim(fname),PIO_DOUBLE,dimid1)
        if (loop == 22) call shr_pcdf_writer1d(fid,trim(fname),iodescd,rf1)
     endif

     if (present(rf2)) then
        fname = 'rf2'
        if (present(rf2n)) fname = trim(rf2n)
        if (loop == 11) call shr_pcdf_readr1d(fid,trim(fname),iodescd,rf2)
        if (loop == 21) call shr_pcdf_defvar1d(fid,trim(fname),PIO_DOUBLE,dimid1)
        if (loop == 22) call shr_pcdf_writer1d(fid,trim(fname),iodescd,rf2)
     endif

     if (present(rf3)) then
        fname = 'rf3'
        if (present(rf3n)) fname = trim(rf3n)
        if (loop == 11) call shr_pcdf_readr1d(fid,trim(fname),iodescd,rf3)
        if (loop == 21) call shr_pcdf_defvar1d(fid,trim(fname),PIO_DOUBLE,dimid1)
        if (loop == 22) call shr_pcdf_writer1d(fid,trim(fname),iodescd,rf3)
     endif

     if (present(rf4)) then
        fname = 'rf4'
        if (present(rf4n)) fname = trim(rf4n)
        if (loop == 11) call shr_pcdf_readr1d(fid,trim(fname),iodescd,rf4)
        if (loop == 21) call shr_pcdf_defvar1d(fid,trim(fname),PIO_DOUBLE,dimid1)
        if (loop == 22) call shr_pcdf_writer1d(fid,trim(fname),iodescd,rf4)
     endif

     if (present(if1)) then
        fname = 'if1'
        if (present(if1n)) fname = trim(if1n)
        if (loop == 11) call shr_pcdf_readi1d(fid,trim(fname),iodesci,if1)
        if (loop == 21) call shr_pcdf_defvar1d(fid,trim(fname),PIO_INT,dimid1)
        if (loop == 22) call shr_pcdf_writei1d(fid,trim(fname),iodesci,if1)
     endif

     if (present(if2)) then
        fname = 'if2'
        if (present(if2n)) fname = trim(if2n)
        if (loop == 11) call shr_pcdf_readi1d(fid,trim(fname),iodesci,if2)
        if (loop == 21) call shr_pcdf_defvar1d(fid,trim(fname),PIO_INT,dimid1)
        if (loop == 22) call shr_pcdf_writei1d(fid,trim(fname),iodesci,if2)
     endif

     if (present(if3)) then
        fname = 'if3'
        if (present(if3n)) fname = trim(if3n)
        if (loop == 11) call shr_pcdf_readi1d(fid,trim(fname),iodesci,if3)
        if (loop == 21) call shr_pcdf_defvar1d(fid,trim(fname),PIO_INT,dimid1)
        if (loop == 22) call shr_pcdf_writei1d(fid,trim(fname),iodesci,if3)
     endif

     if (present(if4)) then
        fname = 'if4'
        if (present(if4n)) fname = trim(if4n)
        if (loop == 11) call shr_pcdf_readi1d(fid,trim(fname),iodesci,if4)
        if (loop == 21) call shr_pcdf_defvar1d(fid,trim(fname),PIO_INT,dimid1)
        if (loop == 22) call shr_pcdf_writei1d(fid,trim(fname),iodesci,if4)
     endif

     if (present(av1)) then
        fname = 'av1_'
        if (present(av1n)) then
           if (trim(av1n) == '') then
              fname = trim(av1n)
           else
              fname = trim(av1n)//'_'
           endif
        endif
        nf = mct_aVect_nRattr(av1)
        do n = 1,nf
           call mct_aVect_getRList(mstring,n,av1)
           vname = trim(fname)//trim(mct_string_toChar(mstring))
           call mct_string_clean(mstring)
           if (loop == 11) call shr_pcdf_readr1d(fid,trim(vname),iodescd,av1%rAttr(n,:))
           if (loop == 21) call shr_pcdf_defvar1d(fid,trim(vname),PIO_DOUBLE,dimid1)
           if (loop == 22) call shr_pcdf_writer1d(fid,trim(vname),iodescd,av1%rAttr(n,:))
        enddo
        nf = mct_aVect_nIattr(av1)
        do n = 1,nf
           call mct_aVect_getIList(mstring,n,av1)
           vname = trim(fname)//trim(mct_string_toChar(mstring))
           call mct_string_clean(mstring)
           if (loop == 11) call shr_pcdf_readi1d(fid,trim(vname),iodesci,av1%iAttr(n,:))
           if (loop == 21) call shr_pcdf_defvar1d(fid,trim(vname),PIO_INT,dimid1)
           if (loop == 22) call shr_pcdf_writei1d(fid,trim(vname),iodesci,av1%iAttr(n,:))
        enddo
     endif

     if (present(av2)) then
        fname = 'av2_'
        if (present(av2n)) then
           if (trim(av2n) == '') then
              fname = trim(av2n)
           else
              fname = trim(av2n)//'_'
           endif
        endif
        nf = mct_aVect_nRattr(av2)
        do n = 1,nf
           call mct_aVect_getRList(mstring,n,av2)
           vname = trim(fname)//trim(mct_string_toChar(mstring))
           call mct_string_clean(mstring)
           if (loop == 11) call shr_pcdf_readr1d(fid,trim(vname),iodescd,av2%rAttr(n,:))
           if (loop == 21) call shr_pcdf_defvar1d(fid,trim(vname),PIO_DOUBLE,dimid1)
           if (loop == 22) call shr_pcdf_writer1d(fid,trim(vname),iodescd,av2%rAttr(n,:))
        enddo
        nf = mct_aVect_nIattr(av2)
        do n = 1,nf
           call mct_aVect_getIList(mstring,n,av2)
           vname = trim(fname)//trim(mct_string_toChar(mstring))
           call mct_string_clean(mstring)
           if (loop == 11) call shr_pcdf_readi1d(fid,trim(vname),iodesci,av2%iAttr(n,:))
           if (loop == 21) call shr_pcdf_defvar1d(fid,trim(vname),PIO_INT,dimid1)
           if (loop == 22) call shr_pcdf_writei1d(fid,trim(vname),iodesci,av2%iAttr(n,:))
        enddo
     endif

     if (present(av3)) then
        fname = 'av3_'
        if (present(av3n)) then
           if (trim(av3n) == '') then
              fname = trim(av3n)
           else
              fname = trim(av3n)//'_'
           endif
        endif
        nf = mct_aVect_nRattr(av3)
        do n = 1,nf
           call mct_aVect_getRList(mstring,n,av3)
           vname = trim(fname)//trim(mct_string_toChar(mstring))
           call mct_string_clean(mstring)
           if (loop == 11) call shr_pcdf_readr1d(fid,trim(vname),iodescd,av3%rAttr(n,:))
           if (loop == 21) call shr_pcdf_defvar1d(fid,trim(vname),PIO_DOUBLE,dimid1)
           if (loop == 22) call shr_pcdf_writer1d(fid,trim(vname),iodescd,av3%rAttr(n,:))
        enddo
        nf = mct_aVect_nIattr(av3)
        do n = 1,nf
           call mct_aVect_getIList(mstring,n,av3)
           vname = trim(fname)//trim(mct_string_toChar(mstring))
           call mct_string_clean(mstring)
           if (loop == 11) call shr_pcdf_readi1d(fid,trim(vname),iodesci,av3%iAttr(n,:))
           if (loop == 21) call shr_pcdf_defvar1d(fid,trim(vname),PIO_INT,dimid1)
           if (loop == 22) call shr_pcdf_writei1d(fid,trim(vname),iodesci,av3%iAttr(n,:))
        enddo
     endif

     if (present(av4)) then
        fname = 'av4_'
        if (present(av4n)) then
           if (trim(av4n) == '') then
              fname = trim(av4n)
           else
              fname = trim(av4n)//'_'
           endif
        endif
        nf = mct_aVect_nRattr(av4)
        do n = 1,nf
           call mct_aVect_getRList(mstring,n,av4)
           vname = trim(fname)//trim(mct_string_toChar(mstring))
           call mct_string_clean(mstring)
           if (loop == 11) call shr_pcdf_readr1d(fid,trim(vname),iodescd,av4%rAttr(n,:))
           if (loop == 21) call shr_pcdf_defvar1d(fid,trim(vname),PIO_DOUBLE,dimid1)
           if (loop == 22) call shr_pcdf_writer1d(fid,trim(vname),iodescd,av4%rAttr(n,:))
        enddo
        nf = mct_aVect_nIattr(av4)
        do n = 1,nf
           call mct_aVect_getIList(mstring,n,av4)
           vname = trim(fname)//trim(mct_string_toChar(mstring))
           call mct_string_clean(mstring)
           if (loop == 11) call shr_pcdf_readi1d(fid,trim(vname),iodesci,av4%iAttr(n,:))
           if (loop == 21) call shr_pcdf_defvar1d(fid,trim(vname),PIO_INT,dimid1)
           if (loop == 22) call shr_pcdf_writei1d(fid,trim(vname),iodesci,av4%iAttr(n,:))
        enddo
     endif

     if (loop == 21) rcode = pio_enddef(fid)
  enddo

  call pio_freedecomp(fid,iodesci)
  call pio_freedecomp(fid,iodescd)
  call pio_closefile(fid)

end subroutine shr_pcdf_readwrite

!===============================================================================
!===============================================================================
subroutine shr_pcdf_defvar0d(fid,fname,vtype)

  implicit none

  type(file_desc_t),intent(in) :: fid
  character(len=*) ,intent(in) :: fname
  integer(IN)      ,intent(in) :: vtype

  !--- local ---
  type(var_desc_t) :: varid
  integer(IN)      :: rcode
  character(len=*),parameter :: subname = '(shr_pcdf_defvar0d) '

  !-------------

  rcode = pio_def_var(fid,trim(fname),vtype,varid)

end subroutine shr_pcdf_defvar0d

!===============================================================================
subroutine shr_pcdf_defvar1d(fid,fname,vtype,dimid)

  implicit none

  type(file_desc_t),intent(in) :: fid
  character(len=*) ,intent(in) :: fname
  integer(IN)      ,intent(in) :: vtype
  integer(IN)      ,intent(in) :: dimid(:)

  !--- local ---
  type(var_desc_t) :: varid
  integer(IN)      :: rcode
  character(len=*),parameter :: subname = '(shr_pcdf_defvar1d) '

  !-------------

  rcode = pio_def_var(fid,trim(fname),vtype,dimid,varid)

end subroutine shr_pcdf_defvar1d

!===============================================================================
subroutine shr_pcdf_readr1d(fid,fname,iodesc,r1d)

  implicit none

  type(file_desc_t),intent(inout) :: fid
  character(len=*) ,intent(in)    :: fname
  type(io_desc_t)  ,intent(inout) :: iodesc
  real(R8)         ,intent(inout) :: r1d(:)

  !--- local ---
  type(var_desc_t) :: varid
  integer(IN)      :: dimid(4),ndims
  integer(IN)      :: vsize,fsize
  integer(IN)      :: rcode
  character(len=*),parameter :: subname = '(shr_pcdf_readr1d) '

  !-------------

  rcode = pio_inq_varid(fid,trim(fname),varid)

!--tcraig, here vsize is global, fsize is local, what check if any?
!  rcode = pio_inq_varndims(fid, varid, ndims)
!  rcode = pio_inq_vardimid(fid, varid, dimid(1:ndims))
!  rcode = pio_inq_dimlen(fid, dimid(1), vsize)
!  fsize = size(r1d)
!  if (vsize /= fsize) then
!     write(shr_log_unit,*) subname,' ERROR: vsize,fsize = ',vsize,fsize
!     call shr_sys_abort(trim(subname)//' ERROR: vsize,fsize')
!  endif

  call pio_read_darray(fid,varid,iodesc,r1d,rcode)

end subroutine shr_pcdf_readr1d

!===============================================================================
subroutine shr_pcdf_writer1d(fid,fname,iodesc,r1d)

  implicit none

  type(file_desc_t),intent(inout) :: fid
  character(len=*) ,intent(in)    :: fname
  type(io_desc_t)  ,intent(inout) :: iodesc
  real(R8)         ,intent(inout) :: r1d(:)

  !--- local ---
  type(var_desc_t) :: varid
  integer(IN)      :: dimid(4)
  integer(IN)      :: vsize,fsize
  real(R8)         :: lfillvalue
  integer(IN)      :: rcode
  character(len=*),parameter :: subname = '(shr_pcdf_writer1d) '

  !-------------

  lfillvalue = fillvalue

  rcode = pio_inq_varid(fid,trim(fname),varid)
  call pio_write_darray(fid, varid, iodesc, r1d, rcode, fillval=lfillvalue)

end subroutine shr_pcdf_writer1d
!===============================================================================
!===============================================================================
subroutine shr_pcdf_readi1d(fid,fname,iodesc,i1d)

  implicit none

  type(file_desc_t),intent(inout) :: fid
  character(len=*) ,intent(in)    :: fname
  type(io_desc_t)  ,intent(inout) :: iodesc
  integer(IN)      ,intent(inout) :: i1d(:)

  !--- local ---
  type(var_desc_t) :: varid
  integer(IN)      :: dimid(4),ndims
  integer(IN)      :: vsize,fsize
  integer(IN)      :: rcode
  character(len=*),parameter :: subname = '(shr_pcdf_readi1d) '

  !-------------

  rcode = pio_inq_varid(fid,trim(fname),varid)

!--tcraig, here vsize is global, fsize is local, what check if any?
!  rcode = pio_inq_varndims(fid, varid, ndims)
!  rcode = pio_inq_vardimid(fid, varid, dimid(1:ndims))
!  rcode = pio_inq_dimlen(fid, dimid(1), vsize)
!  fsize = size(i1d)
!  if (vsize /= fsize) then
!     write(shr_log_unit,*) subname,' ERROR: vsize,fsize = ',vsize,fsize
!     call shr_sys_abort(trim(subname)//' ERROR: vsize,fsize')
!  endif

  call pio_read_darray(fid,varid,iodesc,i1d,rcode)

end subroutine shr_pcdf_readi1d

!===============================================================================
subroutine shr_pcdf_writei1d(fid,fname,iodesc,i1d)

  implicit none

  type(file_desc_t),intent(inout) :: fid
  character(len=*) ,intent(in)    :: fname
  type(io_desc_t)  ,intent(inout) :: iodesc
  integer(IN)      ,intent(inout) :: i1d(:)

  !--- local ---
  type(var_desc_t) :: varid
  integer(IN)      :: dimid(4)
  integer(IN)      :: vsize,fsize
  integer(IN)      :: lfillvalue
  integer(IN)      :: rcode
  character(len=*),parameter :: subname = '(shr_pcdf_writei1d) '

  !-------------

  lfillvalue = ifillvalue

  rcode = pio_inq_varid(fid,trim(fname),varid)
  call pio_write_darray(fid, varid, iodesc, i1d, rcode, fillval=lfillvalue)

end subroutine shr_pcdf_writei1d
!===============================================================================
!===============================================================================
subroutine shr_pcdf_readr0d(fid,fname,r0d)

  implicit none

  type(file_desc_t),intent(inout) :: fid
  character(len=*) ,intent(in)    :: fname
  real(R8)         ,intent(inout) :: r0d

  !--- local ---
  type(var_desc_t) :: varid
  integer(IN)      :: rcode
  character(len=*),parameter :: subname = '(shr_pcdf_readr0d) '

  !-------------

  rcode = pio_inq_varid(fid,trim(fname),varid)
  rcode = pio_get_var(fid,varid,r0d)

end subroutine shr_pcdf_readr0d

!===============================================================================
subroutine shr_pcdf_writer0d(fid,fname,r0d)

  implicit none

  type(file_desc_t),intent(inout) :: fid
  character(len=*) ,intent(in)    :: fname
  real(R8)         ,intent(inout) :: r0d

  !--- local ---
  type(var_desc_t) :: varid
  integer(IN)      :: rcode
  character(len=*),parameter :: subname = '(shr_pcdf_writer0d) '

  !-------------

  rcode = pio_inq_varid(fid,trim(fname),varid)
  rcode = pio_put_var(fid, varid, r0d)

end subroutine shr_pcdf_writer0d
!===============================================================================
!===============================================================================
subroutine shr_pcdf_readi0d(fid,fname,i0d)

  implicit none

  type(file_desc_t),intent(inout) :: fid
  character(len=*) ,intent(in)    :: fname
  integer(IN)      ,intent(inout) :: i0d

  !--- local ---
  type(var_desc_t) :: varid
  integer(IN)      :: rcode
  character(len=*),parameter :: subname = '(shr_pcdf_readi0d) '

  !-------------

  rcode = pio_inq_varid(fid,trim(fname),varid)
  rcode = pio_get_var(fid,varid,i0d)

end subroutine shr_pcdf_readi0d

!===============================================================================
subroutine shr_pcdf_writei0d(fid,fname,i0d)

  implicit none

  type(file_desc_t),intent(inout) :: fid
  character(len=*) ,intent(in)    :: fname
  integer(IN)      ,intent(inout) :: i0d

  !--- local ---
  type(var_desc_t) :: varid
  integer(IN)      :: rcode
  character(len=*),parameter :: subname = '(shr_pcdf_writei0d) '

  !-------------

  rcode = pio_inq_varid(fid,trim(fname),varid)
  rcode = pio_put_var(fid, varid, i0d)

end subroutine shr_pcdf_writei0d
!===============================================================================
!===============================================================================
subroutine shr_pcdf_readdim(fid,fname,dim)

  implicit none

  type(file_desc_t),intent(inout) :: fid
  character(len=*) ,intent(in)    :: fname
  integer(IN)      ,intent(inout) :: dim

  !--- local ---
  integer(IN)      :: dimid
  integer(IN)      :: rcode
  character(len=*),parameter :: subname = '(shr_pcdf_readdim) '

  !-------------

  rcode = pio_inq_dimid(fid,trim(fname),dimid)
  rcode = pio_inq_dimlen(fid,dimid,dim)

end subroutine shr_pcdf_readdim

!===============================================================================
subroutine shr_pcdf_writedim(fid,fname,dim)

  implicit none

  type(file_desc_t),intent(inout) :: fid
  character(len=*) ,intent(in)    :: fname
  integer(IN)      ,intent(inout) :: dim

  !--- local ---
  integer(IN)      :: dimid
  integer(IN)      :: rcode
  character(len=*),parameter :: subname = '(shr_pcdf_writedim) '

  !-------------

  rcode = pio_def_dim(fid,trim(fname),dim,dimid)

end subroutine shr_pcdf_writedim
!===============================================================================
!===============================================================================
!===============================================================================

end module shr_pcdf_mod
