!===============================================================================
! SVN $Id: shr_mct_mod.F90 18548 2009-09-26 23:55:51Z tcraig $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_091114/shr/shr_mct_mod.F90 $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: shr_mct_mod -- higher level mct type routines
!     needed to prevent some circular dependencies
!
! !REVISION HISTORY:
!     2009-Dec-16 - T. Craig - first prototype
!
! !INTERFACE: ------------------------------------------------------------------
module shr_mct_mod

! !USES:

   use shr_kind_mod, only : R8=>SHR_KIND_R8, IN=>SHR_KIND_IN, CL=>SHR_KIND_CL         ! shared kinds
   use shr_sys_mod          ! share system routines
   use shr_mpi_mod          ! mpi layer
   use shr_const_mod        ! constants
   use mct_mod

   use shr_log_mod          ,only: s_loglev               => shr_log_Level
   use shr_log_mod          ,only: s_logunit              => shr_log_Unit

   implicit none
   private

! PUBLIC: Public interfaces

   public :: shr_mct_sMatReadnc
   interface shr_mct_sMatPInitnc
      module procedure shr_mct_sMatPInitnc_mapfile
   end interface
   public :: shr_mct_sMatPInitnc
   public :: shr_mct_sMatReaddnc
   public :: shr_mct_sMatWritednc
   public :: shr_mct_queryConfigFile

!EOP

   !--- local use of kinds ---

   private :: R8, IN, CL

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_mct_sMatReadnc - read all mapping data from a NetCDF SCRIP file
!                              in to a full SparseMatrix
!
! !DESCRIPTION:
!   Read in mapping matrix data from a SCRIP netCDF data file so a sMat.
!
! !REMARKS:
!   Based on cpl_map_read
!
! !REVISION HISTORY:
!     2006 Nov 27: R. Jacob
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_mct_sMatReadnc(sMat,fileName)

  use netcdf

! !INPUT/OUTPUT PARAMETERS:

   type(mct_sMat),intent(inout)  :: sMat
   character(*),intent(in)  :: filename  ! netCDF file to read

!EOP

   !--- local ---
   integer(IN)           :: na      ! size of source domain
   integer(IN)           :: nb      ! size of destination domain
   integer(IN)           :: ns      ! number of non-zero elements in matrix
   integer(IN)           :: igrow   ! aVect index for matrix row
   integer(IN)           :: igcol   ! aVect index for matrix column
   integer(IN)           :: iwgt    ! aVect index for matrix element

   real(R8)   ,allocatable :: rtemp(:)  ! reals
   integer(IN),allocatable :: itemp(:)  ! ints

   integer(IN)           :: rcode   ! netCDF routine return code
   integer(IN)           :: fid     ! netCDF file      ID
   integer(IN)           :: vid     ! netCDF variable  ID
   integer(IN)           :: did     ! netCDF dimension ID

   character(*),parameter :: subName = '(shr_mct_sMatReadnc) '
   character(*),parameter :: F00 = "('(shr_mct_sMatReadnc) ',4a)"
   character(*),parameter :: F01 = '("(shr_mct_sMatReadnc) ",2(a,i9))'

   if (s_loglev > 0) write(s_logunit,F00) "reading mapping matrix data..."

   !----------------------------------------------------------------------------
   ! open & read the file
   !----------------------------------------------------------------------------
   if (s_loglev > 0) write(s_logunit,F00) "* file name                  : ",trim(fileName)
   rcode = nf90_open(filename,NF90_NOWRITE,fid)
   if (rcode /= NF90_NOERR) then
      write(s_logunit,F00) nf90_strerror(rcode)
      call mct_die(subName,"error opening Netcdf file")
   endif

   !--- allocate memory & get matrix data ----------
   rcode = nf90_inq_dimid (fid, 'n_s', did)  ! size of sparse matrix
   rcode = nf90_inquire_dimension(fid, did, len=ns)
   rcode = nf90_inq_dimid (fid, 'n_a', did)  ! size of  input vector
   rcode = nf90_inquire_dimension(fid, did, len=na)
   rcode = nf90_inq_dimid (fid, 'n_b', did)  ! size of output vector
   rcode = nf90_inquire_dimension(fid, did, len=nb)

   if (s_loglev > 0) write(s_logunit,F01) "* matrix dimensions src x dst: ",na,' x',nb
   if (s_loglev > 0) write(s_logunit,F01) "* number of non-zero elements: ",ns

   !----------------------------------------------------------------------------
   ! init the mct sMat data type
   !----------------------------------------------------------------------------
   ! mct_sMat_init must be given the number of rows and columns that
   ! would be in the full matrix.  Nrows= size of output vector=nb.
   ! Ncols = size of input vector = na.
   call mct_sMat_init(sMat, nb, na, ns)

   igrow = mct_sMat_indexIA(sMat,'grow')
   igcol = mct_sMat_indexIA(sMat,'gcol')
   iwgt  = mct_sMat_indexRA(sMat,'weight')

   !!!!!!!!!!!!!!!!!!!!!!!!!!
   ! read and load matrix weights
   allocate(rtemp(ns),stat=rcode)
   if (rcode /= 0) &
     call mct_die(subName,':: allocate weights',rcode)

   rcode = nf90_inq_varid(fid, 'S',vid)
   rcode = nf90_get_var(fid, vid, rtemp)
   if (rcode /= NF90_NOERR .and. s_loglev > 0) then
      write(s_logunit,F00) nf90_strerror(rcode)
   end if

   sMat%data%rAttr(iwgt ,:) =   rtemp(:)

   deallocate(rtemp, stat=rcode)
   if (rcode /= 0) call mct_perr_die(subName,':: deallocate weights',rcode)

   !!!!!!!!!!!!!!!!!!!!!!!!!!
   ! read and load rows
   allocate(itemp(ns),stat=rcode)
   if (rcode /= 0) call mct_perr_die(subName,':: allocate rows',rcode)

   rcode = nf90_inq_varid(fid, 'row', vid)
   rcode = nf90_get_var(fid, vid, itemp)
   if (rcode /= NF90_NOERR .and. s_loglev > 0) then
      write(s_logunit,F00) nf90_strerror(rcode)
   end if

   sMat%data%iAttr(igrow,:) = itemp(:)


   !!!!!!!!!!!!!!!!!!!!!!!!!!
   ! read and load columns
   itemp(:) = 0

   rcode = nf90_inq_varid(fid, 'col', vid)
   rcode = nf90_get_var(fid, vid, itemp)
   if (rcode /= NF90_NOERR .and. s_loglev > 0) then
      write(s_logunit,F00) nf90_strerror(rcode)
   end if

   sMat%data%iAttr(igcol,:) = itemp(:)

   deallocate(itemp, stat=rcode)
   if (rcode /= 0) call mct_perr_die(subName,':: deallocate cols',rcode)

   rcode = nf90_close(fid)

   if (s_loglev > 0) write(s_logunit,F00) "... done reading file"
   call shr_sys_flush(s_logunit)

end subroutine shr_mct_sMatReadnc

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_mct_queryConfigFile - get mct config file info
!
! !DESCRIPTION:
!   Query MCT config file variables
!
! !REMARKS:
!
! !REVISION HISTORY:
!     2013 Aug 17: T. Craig
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_mct_queryConfigFile(mpicom, ConfigFileName, &
           Label1,Value1,Label2,Value2,Label3,Value3)

! !INPUT/OUTPUT PARAMETERS:
   integer          ,intent(in)  :: mpicom
   character(len=*), intent(in)  :: ConfigFileName
   character(len=*), intent(in)  :: Label1
   character(len=*), intent(out) :: Value1
   character(len=*), intent(in) ,optional :: Label2
   character(len=*), intent(out),optional :: Value2
   character(len=*), intent(in) ,optional :: Label3
   character(len=*), intent(out),optional :: Value3

!EOP
   integer :: iret
   character(*),parameter :: subName = '(shr_mct_queryConfigFile) '

   call I90_allLoadF(ConfigFileName,0,mpicom,iret)
   if(iret /= 0) then
      write(s_logunit,*) trim(subname),"Cant find config file ",ConfigFileName
      call shr_sys_abort(trim(subname)//' File Not Found')
   endif

   call i90_label(trim(Label1),iret)
   if(iret /= 0) then
      write(s_logunit,*) trim(subname),"Cant find label ",Label1
      call shr_sys_abort(trim(subname)//' Label1 Not Found')
   endif

   call i90_gtoken(Value1,iret)
   if(iret /= 0) then
      write(s_logunit,*) trim(subname),"Error reading token ",Value1
      call shr_sys_abort(trim(subname)//' Error on read value1')
   endif

   if (present(Label2) .and. present(Value2)) then

      call i90_label(trim(Label2),iret)
      if(iret /= 0) then
         write(s_logunit,*) trim(subname),"Cant find label ",Label2
         call shr_sys_abort(trim(subname)//' Label2 Not Found')
      endif

      call i90_gtoken(Value2,iret)
      if(iret /= 0) then
         write(s_logunit,*)"Error reading token ",Value2
         call shr_sys_abort(trim(subname)//' Error on read value2')
      endif

   endif

   if (present(Label3) .and. present(Value3)) then

      call i90_label(trim(Label3),iret)
      if(iret /= 0) then
         write(s_logunit,*) trim(subname),"Cant find label ",Label3
         call shr_sys_abort(trim(subname)//' Label3 Not Found')
      endif

      call i90_gtoken(Value3,iret)
      if(iret /= 0) then
         write(s_logunit,*)"Error reading token ",Value3
         call shr_sys_abort(trim(subname)//' Error on read value3')
      endif

   endif

   call I90_Release(iret)

end subroutine shr_mct_queryConfigFile

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_mct_sMatPInitnc_mapfile - initialize a SparseMatrixPlus.
!
! !DESCRIPTION:
!   Read in mapping matrix data from a SCRIP netCDF data file in first an
!   Smat and then an SMatPlus
!
! !REMARKS:
!
! !REVISION HISTORY:
!     2012 Feb 27: M. Vertenstein
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_mct_sMatPInitnc_mapfile(sMatP, gsMapX, gsMapY, &
                                       filename, maptype, mpicom, &
                                       ni_i, nj_i, ni_o, nj_o, &
                                       areasrc, areadst)

! !INPUT/OUTPUT PARAMETERS:

   type(mct_sMatP),intent(inout)         :: sMatP
   type(mct_gsMap),intent(in)            :: gsMapX
   type(mct_gsMap),intent(in)            :: gsMapY
   character(*)   ,intent(in)            :: filename        ! scrip map file to read
   character(*)   ,intent(in)            :: maptype         ! map type
   integer        ,intent(in)            :: mpicom
   integer        ,intent(out), optional :: ni_i            ! number of longitudes on input grid
   integer        ,intent(out), optional :: nj_i            ! number of latitudes  on input grid
   integer        ,intent(out), optional :: ni_o            ! number of longitudes on output grid
   integer        ,intent(out), optional :: nj_o            ! number of latitudes  on output grid
   type(mct_Avect),intent(out), optional :: areasrc         ! area of src grid from mapping file
   type(mct_Avect),intent(out), optional :: areadst         ! area of src grid from mapping file

!EOP
   type(mct_sMat ) :: sMati    ! initial sMat from read (either root or decomp)
   type(mct_Avect) :: areasrc_map ! area of src grid from mapping file
   type(mct_Avect) :: areadst_map ! area of dst grid from mapping file

   integer          :: lsize
   integer          :: pe_loc
   logical          :: usevector
   character(len=3) :: Smaptype
   character(*),parameter :: areaAV_field = 'aream'
   character(*),parameter :: F00 = "('(shr_mct_sMatPInitnc) ',4a)"
   character(*),parameter :: F01 = "('(shr_mct_sMatPInitnc) ',a,i10)"

   call shr_mpi_commrank(mpicom,pe_loc)

   if (s_loglev > 0) write(s_logunit,*) " "
   if (s_loglev > 0) write(s_logunit,F00) "Initializing SparseMatrixPlus"
   if (s_loglev > 0) write(s_logunit,F00) "SmatP mapname ",trim(filename)
   if (s_loglev > 0) write(s_logunit,F00) "SmatP maptype ",trim(maptype)

   if (maptype == "X") then
      Smaptype = "src"
   else if(maptype == "Y") then
      Smaptype = "dst"
   end if

   call shr_mpi_commrank(mpicom, pe_loc)

   lsize = mct_gsMap_lsize(gsMapX, mpicom)
   call mct_aVect_init(areasrc_map, rList=areaAV_field, lsize=lsize)

   lsize = mct_gsMap_lsize(gsMapY, mpicom)
   call mct_aVect_init(areadst_map, rList=areaAV_field, lsize=lsize)

   if (present(ni_i) .and. present(nj_i) .and. present(ni_o) .and. present(nj_o)) then
      call shr_mct_sMatReaddnc(sMati, gsMapX, gsMapY, Smaptype, areasrc_map, areadst_map, &
           fileName, pe_loc, mpicom, ni_i, nj_i, ni_o, nj_o)
   else
      call shr_mct_sMatReaddnc(sMati, gsMapX, gsMapY, Smaptype, areasrc_map, areadst_map, &
           fileName, pe_loc, mpicom)
   end if
   call mct_sMatP_Init(sMatP, sMati, gsMapX, gsMapY, 0, mpicom, gsMapX%comp_id)

#ifdef CPP_VECTOR
   !--- initialize the vector parts of the sMat ---
   call mct_sMatP_Vecinit(sMatP)
#endif

   lsize = mct_smat_gNumEl(sMatP%Matrix,mpicom)
   if (s_loglev > 0) write(s_logunit,F01) "Done initializing SmatP, nElements = ",lsize

#ifdef CPP_VECTOR
   usevector = .true.
#else
   usevector = .false.
#endif
   if (present(areasrc)) then
      call mct_aVect_copy(aVin=areasrc_map, aVout=areasrc, vector=usevector)
   end if
   if (present(areadst)) then
      call mct_aVect_copy(aVin=areadst_map, aVout=areadst, vector=usevector)
   end if

   call mct_aVect_clean(areasrc_map)
   call mct_aVect_clean(areadst_map)

   call mct_sMat_Clean(sMati)

end subroutine shr_mct_sMatPInitnc_mapfile

!BOP ===========================================================================
!
! !IROUTINE:  shr_mct_sMatReaddnc - Do a distributed read of a NetCDF SCRIP file and
!                                return weights in a distributed SparseMatrix
!
! !DESCRIPTION:
!     Read in mapping matrix data from a SCRIP netCDF data file using
!     a low memory method and then scatter to all pes.
!
! !REMARKS:
!   This routine leverages gsmaps to determine scatter pattern
!   The scatter is implemented as a bcast of all weights then a local
!     computation on each pe to determine with weights to keep based
!     on gsmap information.
!   The algorithm to determine whether a weight belongs on a pe involves
!     creating a couple local arrays (lsstart and lscount) which are
!     the local values of start and length from the gsmap.  these are
!     sorted via a bubble sort and then searched via a binary search
!     to check whether a global index is on the local pe.
!   The local buffer sizes are estimated up front based on ngridcell/npes
!     plus 20% (see 1.2 below).  If the local buffer size fills up, then
!     the buffer is reallocated 50% large (see 1.5 below) and the fill
!     continues.  The idea is to trade off memory reallocation and copy
!     with memory usage.  1.2 and 1.5 are arbitary, other values may
!     result in better performance.
!   Once all the matrix weights have been read, the sMat is initialized,
!     the values from the buffers are copied in, and everything is deallocated.

! !SEE ALSO:
!    mct/m_SparseMatrix.F90 (MCT source code)
!
! !REVISION HISTORY:
!    2007-Jan-18 - T. Craig -- first version
!    2007-Mar-20 - R. Jacob -- rename to shr_mct_sMatReaddnc.  Remove use of cpl_
!                  variables and move to shr_mct_mod
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine shr_mct_sMatReaddnc(sMat,SgsMap,DgsMap,newdom,areasrc,areadst, &
                            fileName,mytask, mpicom, ni_i,nj_i,ni_o,nj_o )

! !USES:

  use netcdf

! !INPUT/OUTPUT PARAMETERS:

   type(mct_sMat)  ,intent(out)           :: sMat    ! mapping data
   type(mct_gsMap) ,intent(in) ,target    :: SgsMap  ! src gsmap
   type(mct_gSMap) ,intent(in) ,target    :: DgsMap  ! dst gsmap
   character(*)    ,intent(in)            :: newdom  ! type of sMat (src or dst)
   type(mct_Avect) ,intent(out), optional :: areasrc ! area of src grid from mapping file
   type(mct_Avect) ,intent(out), optional :: areadst ! area of dst grid from mapping file
   character(*)    ,intent(in)            :: filename! netCDF file to read
   integer(IN)     ,intent(in)            :: mytask   ! processor id
   integer(IN)     ,intent(in)            :: mpicom  ! communicator
   integer(IN)     ,intent(out), optional :: ni_i    ! number of lons on input grid
   integer(IN)     ,intent(out), optional :: nj_i    ! number of lats on input grid
   integer(IN)     ,intent(out), optional :: ni_o    ! number of lons on output grid
   integer(IN)     ,intent(out), optional :: nj_o    ! number of lats on output grid

! !EOP

   !--- local ---
   integer(IN)           :: n,m     ! generic loop indicies
   integer(IN)           :: na      ! size of source domain
   integer(IN)           :: nb      ! size of destination domain
   integer(IN)           :: ns      ! number of non-zero elements in matrix
   integer(IN)           :: igrow   ! aVect index for matrix row
   integer(IN)           :: igcol   ! aVect index for matrix column
   integer(IN)           :: iwgt    ! aVect index for matrix element
   integer(IN)           :: rsize   ! size of read buffer
   integer(IN)           :: cnt     ! local num of wgts
   integer(IN)           :: cntold  ! local num of wgts, previous read
   integer(IN)           :: start(1)! netcdf read
   integer(IN)           :: count(1)! netcdf read
   integer(IN)           :: bsize   ! buffer size
   integer(IN)           :: nread   ! number of reads
   logical               :: mywt    ! does this weight belong on my pe

   !--- buffers for i/o ---
   real(R8)   ,allocatable :: Sbuf(:)  ! real weights
   integer(IN),allocatable :: Rbuf(:)  ! ints rows
   integer(IN),allocatable :: Cbuf(:)  ! ints cols

   !--- variables associated with local computation of global indices
   integer(IN)             :: lsize     ! size of local seg map
   integer(IN)             :: commsize  ! size of local communicator
   integer(IN),allocatable :: lsstart(:) ! local seg map info
   integer(IN),allocatable :: lscount(:) ! local seg map info
   type(mct_gsMap),pointer :: mygsmap ! pointer to one of the gsmaps
   integer(IN)             :: l1,l2     ! generice indices for sort
   logical                 :: found     ! for sort

   !--- variable assocaited with local data buffers and reallocation
   real(R8)   ,allocatable :: Snew(:),Sold(:)  ! reals
   integer(IN),allocatable :: Rnew(:),Rold(:)  ! ints
   integer(IN),allocatable :: Cnew(:),Cold(:)  ! ints

   integer(IN)           :: rcode   ! netCDF routine return code
   integer(IN)           :: fid     ! netCDF file      ID
   integer(IN)           :: vid     ! netCDF variable  ID
   integer(IN)           :: did     ! netCDF dimension ID
   !--- arbitrary size of read buffer, this is the chunk size weights reading
   integer(IN),parameter :: rbuf_size = 100000

   !--- global source and destination areas ---
   type(mct_Avect) :: areasrc0   ! area of src grid from mapping file
   type(mct_Avect) :: areadst0   ! area of src grid from mapping file

   character(*),parameter :: areaAV_field = 'aream'

   !--- formats ---
   character(*),parameter :: subName = '(shr_mct_sMatReaddnc) '
   character(*),parameter :: F00 = '("(shr_mct_sMatReaddnc) ",4a)'
   character(*),parameter :: F01 = '("(shr_mct_sMatReaddnc) ",2(a,i10))'

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

 call shr_mpi_commsize(mpicom,commsize)
 if (mytask == 0) then
   if (s_loglev > 0) write(s_logunit,F00) "reading mapping matrix data decomposed..."

   !----------------------------------------------------------------------------
   ! open & read the file
   !----------------------------------------------------------------------------
   if (s_loglev > 0) write(s_logunit,F00) "* file name                  : ",trim(fileName)
   rcode = nf90_open(filename,NF90_NOWRITE,fid)
   if (rcode /= NF90_NOERR) then
      print *,'Failed to open file ',trim(filename)
      call shr_sys_abort(trim(subName)//nf90_strerror(rcode))
   end if


   !--- get matrix dimensions ----------
   rcode = nf90_inq_dimid(fid, 'n_s', did)  ! size of sparse matrix
   rcode = nf90_inquire_dimension(fid, did, len=ns)
   rcode = nf90_inq_dimid(fid, 'n_a', did)  ! size of  input vector
   rcode = nf90_inquire_dimension(fid, did, len=na)
   rcode = nf90_inq_dimid(fid, 'n_b', did)  ! size of output vector
   rcode = nf90_inquire_dimension(fid, did, len=nb)

   if (present(ni_i) .and. present(nj_i) .and. present(ni_o) .and. present(nj_o)) then
      rcode = nf90_inq_dimid(fid, 'ni_a', did)  ! number of lons in input grid
      rcode = nf90_inquire_dimension(fid, did, len=ni_i)
      rcode = nf90_inq_dimid(fid, 'nj_a', did)  ! number of lats in input grid
      rcode = nf90_inquire_dimension(fid, did, len=nj_i)
      rcode = nf90_inq_dimid(fid, 'ni_b', did)  ! number of lons in output grid
      rcode = nf90_inquire_dimension(fid, did, len=ni_o)
      rcode = nf90_inq_dimid(fid, 'nj_b', did)  ! number of lats in output grid
      rcode = nf90_inquire_dimension(fid, did, len=nj_o)
   end if

   if (s_loglev > 0) write(s_logunit,F01) "* matrix dims src x dst      : ",na,' x',nb
   if (s_loglev > 0) write(s_logunit,F01) "* number of non-zero elements: ",ns

 endif

   !--- read and load area_a ---
   if (present(areasrc)) then
   if (mytask == 0) then
      call mct_aVect_init(areasrc0,' ',areaAV_field,na)
      rcode = nf90_inq_varid(fid, 'area_a', vid)
      if (rcode /= NF90_NOERR) write(6,F00) nf90_strerror(rcode)
      rcode = nf90_get_var(fid, vid, areasrc0%rAttr)
      if (rcode /= NF90_NOERR) write(6,F00) nf90_strerror(rcode)
   endif
   call mct_aVect_scatter(areasrc0, areasrc, SgsMap, 0, mpicom, rcode)
   if (rcode /= 0) call mct_die("shr_mct_sMatReaddnc","Error on scatter of areasrc0")
   if (mytask == 0) then
!      if (present(dbug)) then
!         if (dbug > 2) then
!            write(6,*) subName,'Size of src ',mct_aVect_lSize(areasrc0)
!            write(6,*) subName,'min/max src ',minval(areasrc0%rAttr(1,:)),maxval(areasrc0%rAttr(1,:))
!         endif
!      end if
      call mct_aVect_clean(areasrc0)
   end if
   end if

   !--- read and load area_b ---
   if (present(areadst)) then
   if (mytask == 0) then
      call mct_aVect_init(areadst0,' ',areaAV_field,nb)
      rcode = nf90_inq_varid(fid, 'area_b', vid)
      if (rcode /= NF90_NOERR) write(6,F00) nf90_strerror(rcode)
      rcode = nf90_get_var(fid, vid, areadst0%rAttr)
      if (rcode /= NF90_NOERR) write(6,F00) nf90_strerror(rcode)
   endif
   call mct_aVect_scatter(areadst0, areadst, DgsMap, 0, mpicom, rcode)
   if (rcode /= 0) call mct_die("shr_mct_sMatReaddnc","Error on scatter of areadst0")
   if (mytask == 0) then
!      if (present(dbug)) then
!         if (dbug > 2) then
!            write(6,*) subName,'Size of dst ',mct_aVect_lSize(areadst0)
!            write(6,*) subName,'min/max dst ',minval(areadst0%rAttr(1,:)),maxval(areadst0%rAttr(1,:))
!         endif
!      end if
      call mct_aVect_clean(areadst0)
   endif
   endif

   if (present(ni_i) .and. present(nj_i) .and. present(ni_o) .and. present(nj_o)) then
      call shr_mpi_bcast(ni_i,mpicom,subName//" MPI in ni_i bcast")
      call shr_mpi_bcast(nj_i,mpicom,subName//" MPI in nj_i bcast")
      call shr_mpi_bcast(ni_o,mpicom,subName//" MPI in ni_o bcast")
      call shr_mpi_bcast(nj_o,mpicom,subName//" MPI in nj_o bcast")
   end if

   call shr_mpi_bcast(ns,mpicom,subName//" MPI in ns bcast")
   call shr_mpi_bcast(na,mpicom,subName//" MPI in na bcast")
   call shr_mpi_bcast(nb,mpicom,subName//" MPI in nb bcast")

   !--- setup local seg map, sorted
   if (newdom == 'src') then
      mygsmap => DgsMap
   elseif (newdom == 'dst') then
      mygsmap => SgsMap
   else
      write(s_logunit,F00) 'ERROR: invalid newdom value = ',newdom
      call shr_sys_abort(trim(subName)//" invalid newdom value")
   endif
   lsize = 0
   do n = 1,size(mygsmap%start)
      if (mygsmap%pe_loc(n) == mytask) then
         lsize=lsize+1
      endif
   enddo
   allocate(lsstart(lsize),lscount(lsize),stat=rcode)
   if (rcode /= 0) call mct_perr_die(subName,':: allocate Lsstart',rcode)

   lsize = 0
   do n = 1,size(mygsmap%start)
      if (mygsmap%pe_loc(n) == mytask) then  ! on my pe
         lsize=lsize+1
         found = .false.
         l1 = 1
         do while (.not.found .and. l1 < lsize)         ! bubble sort copy
            if (mygsmap%start(n) < lsstart(l1)) then
               do l2 = lsize, l1+1, -1
                  lsstart(l2) = lsstart(l2-1)
                  lscount(l2) = lscount(l2-1)
               enddo
               found = .true.
            else
               l1 = l1 + 1
            endif
         enddo
         lsstart(l1) = mygsmap%start(n)
         lscount(l1) = mygsmap%length(n)
      endif
   enddo
   do n = 1,lsize-1
      if (lsstart(n) > lsstart(n+1)) then
         write(s_logunit,F00) ' ERROR: lsstart not properly sorted'
         call shr_sys_abort()
      endif
   enddo

   rsize = min(rbuf_size,ns)                     ! size of i/o chunks
   bsize = ((ns/commsize) + 1 ) * 1.2   ! local temporary buffer size
   if (ns == 0) then
      nread = 0
   else
      nread = (ns-1)/rsize + 1                      ! num of reads to do
   endif

   allocate(Sbuf(rsize),Rbuf(rsize),Cbuf(rsize),stat=rcode)
   if (rcode /= 0) call mct_perr_die(subName,':: allocate Sbuf',rcode)
   allocate(Snew(bsize),Cnew(bsize),Rnew(bsize),stat=rcode)
   if (rcode /= 0) call mct_perr_die(subName,':: allocate Snew1',rcode)

   cnt = 0
   do n = 1,nread
      start(1) = (n-1)*rsize + 1
      count(1) = min(rsize,ns-start(1)+1)

      !--- read data on root pe
      if (mytask== 0) then
         rcode = nf90_inq_varid(fid, 'S', vid)
         rcode = nf90_get_var(fid, vid, Sbuf, start, count)
         if (rcode /= NF90_NOERR .and. s_loglev > 0) then
            write(s_logunit,F00) nf90_strerror(rcode)
         end if

         rcode = nf90_inq_varid(fid, 'row', vid)
         rcode = nf90_get_var(fid, vid, Rbuf, start, count)
         if (rcode /= NF90_NOERR .and. s_loglev > 0) then
            write(s_logunit,F00) nf90_strerror(rcode)
         end if

         rcode = nf90_inq_varid(fid, 'col', vid)
         rcode = nf90_get_var(fid, vid, Cbuf, start, count)
         if (rcode /= NF90_NOERR .and. s_loglev > 0) then
            write(s_logunit,F00) nf90_strerror(rcode)
         end if
      endif

      !--- send S, row, col to all pes
      call shr_mpi_bcast(Sbuf,mpicom,subName//" MPI in Sbuf bcast")
      call shr_mpi_bcast(Rbuf,mpicom,subName//" MPI in Rbuf bcast")
      call shr_mpi_bcast(Cbuf,mpicom,subName//" MPI in Cbuf bcast")

      !--- now each pe keeps what it should
      do m = 1,count(1)
         !--- should this weight be on my pe
         if (newdom == 'src') then
            mywt = mct_myindex(Rbuf(m),lsstart,lscount)
         elseif (newdom == 'dst') then
            mywt = mct_myindex(Cbuf(m),lsstart,lscount)
         endif

         if (mywt) then
            cntold = cnt
            cnt = cnt + 1

            !--- new arrays need to be bigger
            if (cnt > bsize) then
               !--- allocate old arrays and copy new into old
               allocate(Sold(cntold),Rold(cntold),Cold(cntold),stat=rcode)
               if (rcode /= 0) call mct_perr_die(subName,':: allocate old',rcode)
               Sold(1:cntold) = Snew(1:cntold)
               Rold(1:cntold) = Rnew(1:cntold)
               Cold(1:cntold) = Cnew(1:cntold)

               !--- reallocate new to bigger size, increase buffer by 50% (arbitrary)
               deallocate(Snew,Rnew,Cnew,stat=rcode)
               if (rcode /= 0) call mct_perr_die(subName,':: allocate new',rcode)
               bsize = 1.5 * bsize
               if (s_loglev > 1) write(s_logunit,F01) ' reallocate bsize to ',bsize
               allocate(Snew(bsize),Rnew(bsize),Cnew(bsize),stat=rcode)
               if (rcode /= 0) call mct_perr_die(subName,':: allocate old',rcode)

               !--- copy data back into new
               Snew(1:cntold) = Sold(1:cntold)
               Rnew(1:cntold) = Rold(1:cntold)
               Cnew(1:cntold) = Cold(1:cntold)
               deallocate(Sold,Rold,Cold,stat=rcode)
               if (rcode /= 0) call mct_perr_die(subName,':: deallocate old',rcode)
            endif

            Snew(cnt) = Sbuf(m)
            Rnew(cnt) = Rbuf(m)
            Cnew(cnt) = Cbuf(m)
         endif
      enddo  ! count
   enddo   ! nread

   deallocate(Sbuf,Rbuf,Cbuf, stat=rcode)
   if (rcode /= 0) call mct_perr_die(subName,':: deallocate Sbuf',rcode)

   !----------------------------------------------------------------------------
   ! init the mct sMat data type
   !----------------------------------------------------------------------------
   ! mct_sMat_init must be given the number of rows and columns that
   ! would be in the full matrix.  Nrows= size of output vector=nb.
   ! Ncols = size of input vector = na.
   call mct_sMat_init(sMat, nb, na, cnt)

   igrow = mct_sMat_indexIA(sMat,'grow')
   igcol = mct_sMat_indexIA(sMat,'gcol')
   iwgt  = mct_sMat_indexRA(sMat,'weight')

   if (cnt /= 0) then
      sMat%data%rAttr(iwgt ,1:cnt) = Snew(1:cnt)
      sMat%data%iAttr(igrow,1:cnt) = Rnew(1:cnt)
      sMat%data%iAttr(igcol,1:cnt) = Cnew(1:cnt)
   endif
   deallocate(Snew,Rnew,Cnew, stat=rcode)
   deallocate(lsstart,lscount,stat=rcode)
   if (rcode /= 0) call mct_perr_die(subName,':: deallocate new',rcode)

   if (mytask == 0) then
      rcode = nf90_close(fid)
      if (s_loglev > 0) write(s_logunit,F00) "... done reading file"
      call shr_sys_flush(s_logunit)
   endif

end subroutine shr_mct_sMatReaddnc

!BOP ===========================================================================
!
! !IROUTINE:  shr_mct_sMatWritednc - Do a distributed write of a NetCDF SCRIP file
!                                 based on a distributed SparseMatrix
!
! !DESCRIPTION:
!     Write out mapping matrix data from a SCRIP netCDF data file using
!     a low memory method.
!
! !SEE ALSO:
!    mct/m_SparseMatrix.F90 (MCT source code)
!
! !REVISION HISTORY:
!    2009-Dec-15 - T. Craig -- first version
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine shr_mct_sMatWritednc(sMat,iosystem, io_type, io_format, fileName,compid, mpicom)

! !USES:
  use pio, only : iosystem_desc_t
   use shr_pcdf_mod, only : shr_pcdf_readwrite
   implicit none
#include <mpif.h>

! !INPUT/OUTPUT PARAMETERS:

   type(mct_sMat)  ,intent(in)   :: sMat     ! mapping data
   type(iosystem_desc_t)         :: iosystem ! PIO subsystem description
   integer(IN)     ,intent(in)   :: io_type  ! type of io interface for this file
   integer(IN)     ,intent(in)   :: io_format ! type of io netcdf3 format for this file
   character(*)    ,intent(in)   :: filename ! netCDF file to read
   integer(IN)     ,intent(in)   :: compid   ! processor id
   integer(IN)     ,intent(in)   :: mpicom   ! communicator

 ! !local
   integer(IN) :: na,nb,ns,lsize,npes,ierr,my_task,n
   integer(IN), pointer :: start(:),count(:),ssize(:),pe_loc(:)
   integer(IN), pointer :: expvari(:)
   real(R8)   , pointer :: expvarr(:)
   type(mct_gsmap) :: gsmap
   type(mct_avect) :: AV
   character(*),parameter :: subName = '(shr_mct_sMatWritednc) '

!----------------------------------------

   call MPI_COMM_SIZE(mpicom,npes,ierr)
   call MPI_COMM_RANK(mpicom,my_task,ierr)
   allocate(start(npes),count(npes),ssize(npes),pe_loc(npes))

   na = mct_sMat_ncols(sMat)
   nb = mct_sMat_nrows(sMat)
   ns = mct_sMat_gNumEl(sMat,mpicom)
   lsize = mct_sMat_lsize(sMat)

   count(:) = -999
   pe_loc(:) = -999
   ssize(:) = 1
   call MPI_GATHER(lsize,1,MPI_INTEGER,count,ssize,MPI_INTEGER,0,mpicom,ierr)

   if (my_task == 0) then
      if (minval(count) < 0) then
         call shr_sys_abort(subname//' ERROR: count invalid')
      endif

      start(1) = 1
      pe_loc(1) = 0
      do n = 2,npes
         start(n) = start(n-1)+count(n-1)
         pe_loc(n) = n-1
      enddo

   endif

   call mct_gsmap_init(gsmap,npes,start,count,pe_loc,0,mpicom,compid,ns)
   deallocate(start,count,ssize,pe_loc)

   call mct_aVect_init(AV,iList='row:col',rList='S',lsize=lsize)
   allocate(expvari(lsize))
   call mct_sMat_ExpGRowI(sMat,expvari)
   AV%iAttr(1,:) = expvari(:)
   call mct_sMat_ExpGColI(sMat,expvari)
   AV%iAttr(2,:) = expvari(:)
   deallocate(expvari)
   allocate(expvarr(lsize))
   call mct_sMat_ExpMatrix(sMat,expvarr)
   AV%rAttr(1,:) = expvarr(:)
   deallocate(expvarr)

   call shr_pcdf_readwrite('write',iosystem,io_type, trim(filename),mpicom,gsmap,clobber=.false.,io_format=io_format, &
      id1=na,id1n='n_a',id2=nb,id2n='n_b',id3=ns,id3n='n_s',av1=AV,av1n='')

   call mct_gsmap_clean(gsmap)
   call mct_avect_clean(AV)

end subroutine shr_mct_sMatWritednc
!===============================================================================

end module shr_mct_mod
