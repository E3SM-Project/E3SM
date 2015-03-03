!===============================================================================
! SVN $Id: seq_io_mod.F90 50621 2013-08-30 02:53:41Z mvertens $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/drv/seq_mct/branches/comptype/shr/seq_io_mod.F90 $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: seq_io_read_mod -- reads integer, real arrays and chacter of driver files
!
! !REMARKS:
!
! !REVISION HISTORY:
!    2007-Oct-26 - T. Craig first version
!    2007-Dec-06 - T. Craig update and improve
!    2008-Feb-16 - J. Edwards convert to PIO
!    2010-Nov    - J. Edwards move PIO init and namelists from components to driver
! Current Problems
!  - the original use of seq_io will now ONLY work with the cpl because
!    of hardwiring cpl_io_type and cpl_io_iosystem.  want the original
!    io capabilities to be usable by any component
!  - the init1 method depends on seq_comm for name consistency but seq_comm_init 
!    wants to be called after init1 so the global_comm can be modified for
!    async IO.  this needs to be reconciled.
!  - this routine stores information for all components but most methods are
!    hardwired to work only for the coupler.  should all the components info
!    be stored here or should this be more a general set of methods that are
!    reusable as it's original intent.
!
! !INTERFACE: ------------------------------------------------------------------

module seq_io_read_mod

  ! !USES:

  use shr_kind_mod, only: r8 => shr_kind_r8, in => shr_kind_in
  use shr_kind_mod, only: cl => shr_kind_cl, cs => shr_kind_cs
  use shr_pio_mod,  only: shr_pio_getiosys, shr_pio_getiotype
  use shr_sys_mod       ! system calls
  use seq_comm_mct
  use mct_mod           ! mct wrappers
  use pio

  implicit none
  private

  ! !PUBLIC TYPES:

  ! none

  ! !PUBLIC MEMBER FUNCTIONS:

  public seq_io_read

  ! !PUBLIC DATA MEMBERS

  ! none

  !EOP

  interface seq_io_read
     module procedure seq_io_read_int
     module procedure seq_io_read_int1d
     module procedure seq_io_read_r8
     module procedure seq_io_read_r81d
     module procedure seq_io_read_char
  end interface

!-------------------------------------------------------------------------------
! Local data
!-------------------------------------------------------------------------------

   character(*) , parameter :: prefix = "seq_io_"
   character(*) , parameter :: version ='cpl7v10'
   character(*) , parameter :: version0='cpl7v00'
   character(CL)            :: charvar   ! buffer for string read/write

!=================================================================================
contains
!=================================================================================

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_io_read_int - read scalar integer from netcdf file
  !
  ! !DESCRIPTION:
  !    Read scalar integer from netcdf file
  !
  ! !REVISION HISTORY:
  !    2007-Oct-26 - T. Craig - initial version
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_io_read_int(filename,pioid,idata,dname)

    ! !INPUT/OUTPUT PARAMETERS:
    implicit none
    character(len=*),intent(in)    :: filename ! file
    type(file_desc_t) :: pioid
    integer         ,intent(inout) :: idata    ! integer data
    character(len=*),intent(in)    :: dname    ! name of data

    !EOP

    integer :: i1d(1)
    character(*),parameter :: subName = '(seq_io_read_int) '

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    call seq_io_read_int1d(filename,pioid,i1d,dname)
    idata = i1d(1)

  end subroutine seq_io_read_int

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_io_read_int1d - read 1d integer from netcdf file
  !
  ! !DESCRIPTION:
  !    Read 1d integer array from netcdf file
  !
  ! !REVISION HISTORY:
  !    2007-Oct-26 - T. Craig - initial version
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_io_read_int1d(filename,pioid,idata,dname)

    ! !INPUT/OUTPUT PARAMETERS:
    implicit none
    character(len=*),intent(in)   :: filename ! file
    type(file_desc_t) :: pioid
    integer(in)     ,intent(inout):: idata(:)  ! integer data
    character(len=*),intent(in)   :: dname    ! name of data

    !EOP

    integer(in)                     :: rcode
    type(var_desc_t)                :: varid
    logical                         :: exists
    character(CL)                   :: name1
    character(*),parameter          :: subName = '(seq_io_read_int1d) '
    logical :: addprefix
    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------
    call seq_io_read_openfile(filename,pioid,addprefix) 

    if (addprefix) then
       name1 = trim(prefix)//trim(dname)
    else
       name1 = trim(dname)
    endif
    rcode = pio_inq_varid(pioid,trim(name1),varid)
    rcode = pio_get_var(pioid,varid,idata)

  end subroutine seq_io_read_int1d

  subroutine seq_io_read_openfile(filename,pioid,addprefix)
    character(len=*), intent(in) :: filename
    type(file_desc_t) :: pioid
    logical, intent(out) :: addprefix
    logical :: exists
    integer(in)                     :: iam,mpicom
    type(iosystem_desc_t) , pointer :: cpl_io_subsystem
    character(len=seq_comm_namelen) :: cpl_name
    integer(in)                     :: cpl_pio_iotype
    logical, save :: laddprefix
    integer :: rcode
    character(CL)                   :: lversion
    character(*),parameter          :: subName = '(seq_io_read_openfile) '

    if(.not. pio_file_is_open(pioid)) then
       cpl_name         =  seq_comm_name(CPLID)
       cpl_io_subsystem => shr_pio_getiosys(cpl_name)
       cpl_pio_iotype   =  shr_pio_getiotype(cpl_name)

       call seq_comm_setptrs(CPLID,iam=iam,mpicom=mpicom)
       if (iam==0) inquire(file=trim(filename),exist=exists)
       call shr_mpi_bcast(exists,mpicom,'seq_io_read_openfile')
       if (exists) then
          rcode = pio_openfile(cpl_io_subsystem, pioid, cpl_pio_iotype, trim(filename),pio_nowrite)
          call pio_seterrorhandling(pioid,PIO_BCAST_ERROR)
          rcode = pio_get_att(pioid,pio_global,"file_version",lversion)
          call pio_seterrorhandling(pioid,PIO_INTERNAL_ERROR)
          if (trim(lversion) == trim(version)) then
             laddprefix=.false.
          else
             laddprefix=.true.
          endif
       else
          if(iam==0) write(logunit,*) subname,' ERROR: file invalid ',trim(filename)
          call shr_sys_abort()
       endif
    endif
    addprefix = laddprefix

  end subroutine seq_io_read_openfile


  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_io_read_r8 - read scalar double from netcdf file
  !
  ! !DESCRIPTION:
  !    Read scalar double from netcdf file
  !
  ! !REVISION HISTORY:
  !    2007-Oct-26 - T. Craig - initial version
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_io_read_r8(filename,pioid,rdata,dname)

    ! !INPUT/OUTPUT PARAMETERS:
    implicit none
    character(len=*),intent(in)    :: filename ! file
    type(file_desc_t) :: pioid
    real(r8)        ,intent(inout) :: rdata    ! real data
    character(len=*),intent(in)    :: dname    ! name of data

    !EOP

    real(r8) :: r1d(1)
    character(*),parameter :: subName = '(seq_io_read_r8) '

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    call seq_io_read_r81d(filename,pioid,r1d,dname)
    rdata = r1d(1)

  end subroutine seq_io_read_r8

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_io_read_r81d - read 1d double array from netcdf file
  !
  ! !DESCRIPTION:
  !    Read 1d double array from netcdf file
  !
  ! !REVISION HISTORY:
  !    2007-Oct-26 - T. Craig - initial version
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_io_read_r81d(filename,pioid,rdata,dname)

    ! !INPUT/OUTPUT PARAMETERS:
    implicit none
    character(len=*),intent(in)    :: filename ! file
    type(file_desc_t) :: pioid
    real(r8)        ,intent(inout) :: rdata(:) ! real data
    character(len=*),intent(in)    :: dname    ! name of data

    !EOP

    type(var_desc_t)                :: varid
    character(CL)                   :: name1
    character(*),parameter          :: subName = '(seq_io_read_r81d) '
    logical :: addprefix
    integer :: rcode
    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------
    call seq_io_read_openfile(filename,pioid,addprefix) 

    if (addprefix) then
       name1 = trim(prefix)//trim(dname)
    else
       name1 = trim(dname)
    endif

    rcode = pio_inq_varid(pioid,trim(name1),varid)
    rcode = pio_get_var(pioid,varid,rdata)

  end subroutine seq_io_read_r81d

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_io_read_char - read char string from netcdf file
  !
  ! !DESCRIPTION:
  !    Read char string from netcdf file
  !
  ! !REVISION HISTORY:
  !    2010-July-06 - T. Craig - initial version
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_io_read_char(filename,pioid,rdata,dname)

    ! !INPUT/OUTPUT PARAMETERS:
    implicit none
    character(len=*),intent(in)    :: filename ! file
    type(file_desc_t) :: pioid
    character(len=*),intent(inout) :: rdata    ! character data
    character(len=*),intent(in)    :: dname    ! name of data

    !EOP

    type(var_desc_t)                :: varid
    character(CL)                   :: name1
    character(*),parameter          :: subName = '(seq_io_read_char) '
    logical :: addprefix
    integer :: rcode
    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------
    call seq_io_read_openfile(filename,pioid,addprefix) 

    if (addprefix) then
       name1 = trim(prefix)//trim(dname)
    else
       name1 = trim(dname)
    endif

    rcode = pio_inq_varid(pioid,trim(name1),varid)
    rcode = pio_get_var(pioid,varid,charvar)
    rdata = trim(charvar)

  end subroutine seq_io_read_char

  !===============================================================================
!===============================================================================
end module seq_io_read_mod
