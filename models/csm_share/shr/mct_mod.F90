!===============================================================================
! SVN $Id: mct_mod.F90 56641 2014-01-15 23:10:15Z tcraig $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_140509/shr/mct_mod.F90 $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: mct_mod -- provides a standard API naming convention for MCT code
!
! !DESCRIPTION:
!    This module should be used instead of accessing mct modules directly.  
!    This module:
!    \begin{itemize}
!    \item Uses Fortran {\sf use} renaming of MCT routines and data types so that they
!          all have an mct\_ prefix and related data types and routines have related names. 
!    \item Provides easy and uniform access to 
!          all MCT routines and data types that must be accessed.
!    \item Provides a convienient list of 
!          all MCT routines and data types that can be accessed.
!    \item Blocks access to MCT routines that are not used in cpl6.
!    \end{itemize}
!    This module also includes some MCT-only functions to augment
!    the MCT library.
!
! !REVISION HISTORY:
!     2001-Aug-14 - B. Kauffman - first prototype
!     2006-Apr-13 - M. Vertenstein - modified for sequential mode
!     2007-Mar-01 - R. Jacob - moved to shr
!
! !INTERFACE: ------------------------------------------------------------------
module mct_mod

! !USES:

   use shr_kind_mod         ! shared kinds
   use shr_sys_mod          ! share system routines
   use shr_mpi_mod          ! mpi layer
   use shr_const_mod        ! constants
   use shr_string_mod       ! string functions

   use shr_log_mod          ,only: s_loglev               => shr_log_Level
   use shr_log_mod          ,only: s_logunit              => shr_log_Unit

   use m_MCTWorld           ,only: mct_world_init         => init

   use m_AttrVect           ,only: mct_aVect              => AttrVect
   use m_AttrVect           ,only: mct_aVect_init         => init
   use m_AttrVect           ,only: mct_aVect_clean        => clean
   use m_AttrVect           ,only: mct_aVect_zero         => zero
   use m_AttrVect           ,only: mct_aVect_lsize        => lsize
   use m_AttrVect           ,only: mct_aVect_indexIA      => indexIA
   use m_AttrVect           ,only: mct_aVect_indexRA      => indexRA
   use m_AttrVect           ,only: mct_aVect_importIattr  => importIattr
   use m_AttrVect           ,only: mct_aVect_exportIattr  => exportIattr
   use m_AttrVect           ,only: mct_aVect_importRattr  => importRattr
   use m_AttrVect           ,only: mct_aVect_exportRattr  => exportRattr
   use m_AttrVect           ,only: mct_aVect_getIList     => getIList
   use m_AttrVect           ,only: mct_aVect_getRList     => getRList
   use m_AttrVect           ,only: mct_aVect_getIList2c   => getIListToChar
   use m_AttrVect           ,only: mct_aVect_getRList2c   => getRListToChar
   use m_AttrVect           ,only: mct_aVect_exportIList2c=> exportIListToChar
   use m_AttrVect           ,only: mct_aVect_exportRList2c=> exportRListToChar
   use m_AttrVect           ,only: mct_aVect_nIAttr       => nIAttr
   use m_AttrVect           ,only: mct_aVect_nRAttr       => nRAttr
   use m_AttrVect           ,only: mct_aVect_copy         => Copy
   use m_AttrVect           ,only: mct_aVect_permute      => Permute
   use m_AttrVect           ,only: mct_aVect_unpermute    => Unpermute
   use m_AttrVect           ,only: mct_aVect_SharedIndices=> AVSharedIndices
   use m_AttrVect           ,only: mct_aVect_setSharedIndices=> SharedIndices
   use m_AttrVectComms      ,only: mct_aVect_scatter      => scatter
   use m_AttrVectComms      ,only: mct_aVect_gather       => gather 
   use m_AttrVectComms      ,only: mct_aVect_bcast        => bcast  

   use m_GeneralGrid        ,only: mct_gGrid              => GeneralGrid
   use m_GeneralGrid        ,only: mct_gGrid_init         => init
   use m_GeneralGrid        ,only: mct_gGrid_clean        => clean
   use m_GeneralGrid        ,only: mct_gGrid_dims         => dims
   use m_GeneralGrid        ,only: mct_gGrid_lsize        => lsize
   use m_GeneralGrid        ,only: mct_ggrid_indexIA      => indexIA
   use m_GeneralGrid        ,only: mct_gGrid_indexRA      => indexRA
   use m_GeneralGrid        ,only: mct_gGrid_exportRattr  => exportRattr
   use m_GeneralGrid        ,only: mct_gGrid_importRattr  => importRattr
   use m_GeneralGrid        ,only: mct_gGrid_exportIattr  => exportIattr
   use m_GeneralGrid        ,only: mct_gGrid_importIattr  => importIattr
   use m_GeneralGrid        ,only: mct_gGrid_permute      => permute
   use m_GeneralGridComms   ,only: mct_gGrid_scatter      => scatter
   use m_GeneralGridComms   ,only: mct_gGrid_gather       => gather 
   use m_GeneralGridComms   ,only: mct_gGrid_bcast        => bcast  

   use m_Transfer           ,only: mct_send               => Send
   use m_Transfer           ,only: mct_recv               => Recv
      
   use m_GlobalSegMap       ,only: mct_gsMap              => GlobalSegMap
   use m_GlobalSegMap       ,only: mct_gsMap_init         => init
   use m_GlobalSegMap       ,only: mct_gsMap_clean        => clean
   use m_GlobalSegMap       ,only: mct_gsMap_lsize        => lsize
   use m_GlobalSegMap       ,only: mct_gsMap_gsize        => gsize
   use m_GlobalSegMap       ,only: mct_gsMap_gstorage     => GlobalStorage
   use m_GlobalSegMap       ,only: mct_gsMap_ngseg        => ngseg
   use m_GlobalSegMap       ,only: mct_gsMap_nlseg        => nlseg
   use m_GlobalSegMap       ,only: mct_gsMap_OP           => OrderedPoints
   use m_GlobalSegMap       ,only: mct_gsMap_maxnlseg     => max_nlseg
   use m_GlobalSegMap       ,only: mct_gsMap_activepes    => active_pes
   use m_GlobalSegMap       ,only: mct_gsMap_copy         => copy
   use m_GlobalSegMap       ,only: mct_gsMap_increasing   => increasing
   use m_GlobalSegMap       ,only: mct_gsMap_orderedPoints=> OrderedPoints
   use m_GlobalSegMapComms  ,only: mct_gsMap_bcast        => bcast 

   use m_Rearranger         ,only: mct_rearr              => Rearranger
   use m_Rearranger         ,only: mct_rearr_init         => init
   use m_Rearranger         ,only: mct_rearr_clean        => clean
   use m_Rearranger         ,only: mct_rearr_print        => print
   use m_Rearranger         ,only: mct_rearr_rearrange    => rearrange

   use m_Router             ,only: mct_router             => Router
   use m_Router             ,only: mct_router_init        => init

   use m_SparseMatrixToMaps ,only: mct_sMat_2XgsMap       => SparseMatrixToXGlobalSegMap
   use m_SparseMatrixToMaps ,only: mct_sMat_2YgsMap       => SparseMatrixToYGlobalSegMap
   use m_SparseMatrix       ,only: mct_sMat               => SparseMatrix
   use m_SparseMatrix       ,only: mct_sMat_Init          => init
   use m_SparseMatrix       ,only: mct_sMat_Vecinit       => vecinit
   use m_SparseMatrix       ,only: mct_sMat_Clean         => clean
   use m_SparseMatrix       ,only: mct_sMat_indexIA       => indexIA
   use m_SparseMatrix       ,only: mct_sMat_indexRA       => indexRA
   use m_SparseMatrix       ,only: mct_sMat_lsize         => lsize
   use m_SparseMatrix       ,only: mct_sMat_nrows         => nRows
   use m_SparseMatrix       ,only: mct_sMat_ncols         => nCols
   use m_SparseMatrix       ,only: mct_sMat_SortPermute   => SortPermute
   use m_SparseMatrix       ,only: mct_sMat_GNumEl        => GlobalNumElements
   use m_SparseMatrix       ,only: mct_sMat_ImpGRowI      => ImportGlobalRowIndices
   use m_SparseMatrix       ,only: mct_sMat_ImpGColI      => ImportGlobalColumnIndices
   use m_SparseMatrix       ,only: mct_sMat_ImpLRowI      => ImportLocalRowIndices
   use m_SparseMatrix       ,only: mct_sMat_ImpLColI      => ImportLocalColumnIndices
   use m_SparseMatrix       ,only: mct_sMat_ImpMatrix     => ImportMatrixElements
   use m_SparseMatrix       ,only: mct_sMat_ExpGRowI      => ExportGlobalRowIndices
   use m_SparseMatrix       ,only: mct_sMat_ExpGColI      => ExportGlobalColumnIndices
   use m_SparseMatrix       ,only: mct_sMat_ExpLRowI      => ExportLocalRowIndices
   use m_SparseMatrix       ,only: mct_sMat_ExpLColI      => ExportLocalColumnIndices
   use m_SparseMatrix       ,only: mct_sMat_ExpMatrix     => ExportMatrixElements
   use m_SparseMatrixComms  ,only: mct_sMat_ScatterByRow  => ScatterByRow
   use m_SparseMatrixComms  ,only: mct_sMat_ScatterByCol  => ScatterByColumn
   use m_SparseMatrixPlus   ,only: mct_sMatP              => SparseMatrixPlus
   use m_SparseMatrixPlus   ,only: mct_sMatP_Init         => init
   use m_SparseMatrixPlus   ,only: mct_sMatP_Vecinit      => vecinit
   use m_SparseMatrixPlus   ,only: mct_sMatP_clean        => clean
   use m_MatAttrVectMul     ,only: mct_sMat_avMult        => sMatAvMult
   use m_GlobalToLocal      ,only: mct_sMat_g2lMat        => GlobalToLocalMatrix

   use m_List               ,only: mct_list               => list     
   use m_List               ,only: mct_list_init          => init
   use m_List               ,only: mct_list_get           => get 
   use m_List               ,only: mct_list_nitem         => nitem 
   use m_List               ,only: mct_list_clean         => clean
   use m_string             ,only: mct_string             => string 
   use m_string             ,only: mct_string_clean       => clean
   use m_string             ,only: mct_string_toChar      => toChar 
   use m_die                ,only: mct_perr_die           => mp_perr_die
   use m_die                ,only: mct_die                => die
   use m_inpak90

   use m_Permuter           ,only: mct_permute            => Permute

   use m_MergeSorts         ,only: mct_indexset           => IndexSet
   use m_MergeSorts         ,only: mct_indexsort          => IndexSort

   implicit none

   public :: mct_aVect_info
   public :: mct_aVect_fldIndex
   public :: mct_aVect_sharedFields
   public :: mct_aVect_initSharedFields
   public :: mct_aVect_getRAttr
   public :: mct_aVect_putRAttr
   public :: mct_aVect_accum
   public :: mct_aVect_avg
   public :: mct_avect_mult
   public :: mct_avect_vecmult
   public :: mct_rearr_rearrange_fldlist
   public :: mct_gsmap_identical

   logical,public :: mct_usealltoall = .false.
   logical,public :: mct_usevector = .false.

!EOP

   !--- local kinds ---
   integer,parameter,private :: R8 = SHR_KIND_R8
   integer,parameter,private :: IN = SHR_KIND_IN
   integer,parameter,private :: CL = SHR_KIND_CL
   integer,parameter,private :: CX = SHR_KIND_CX
   integer,parameter,private :: CXX = SHR_KIND_CXX

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: mct_aVect_info - print out aVect info for debugging
!
! !DESCRIPTION:
!     Print out information about the input MCT {\it AttributeVector}
!     {\tt aVect} to stdout. {\tt flag} sets the level of information:
!     \begin{enumerate}
!     \item  print out names of attributes in {\tt aVect}.
!     \item  also print out local max and min of data in {\tt aVect}.
!     \item  also print out global max and min of data in {\tt aVect}.
!     \item  Same as 3 but include name of this routine.
!     \end{enumerate}
!     If {\tt flag} is 3 or higher, then optional argument {\tt comm}
!     must be provided.
!     If optional argument {\tt fld} is present, only information for
!     that field will be printed.
!     If optional argument {\tt istr} is present, it will be output
!     before any of the information.
!
!
! !REVISION HISTORY:
!     2003 Jul 01 - B. Kauffman, T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine mct_aVect_info(flag,aVect,comm,pe,fld,istr)

! !USES:
  
! !INPUT/OUTPUT PARAMETERS:

   integer(IN)    ,intent(in)           :: flag  ! info level flag
   type(mct_aVect),intent(in)           :: aVect ! Attribute vector
   integer(IN)    ,intent(in),optional  :: comm  ! MPI communicator
   integer(IN)    ,intent(in),optional  :: pe    ! processor number
   character(*)   ,intent(in),optional  :: fld   ! fld
   character(*)   ,intent(in),optional  :: istr  ! string for print

!EOP

   !--- local ---
   integer(IN)          :: i,j,k,n      ! generic indicies
   integer(IN)          :: ks,ke        ! start and stop k indices
   integer(IN)          :: nflds        ! number of flds in AV to diagnose
   integer(IN)          :: nsize        ! grid point size of AV
   type(mct_string)     :: item         ! mct string
   character(CL)        :: itemc        ! item converted to char
   integer(IN)          :: comm_loc     ! local variable for comm
   integer(IN)          :: pe_loc       ! local variable for pe
   logical              :: commOK       ! is comm available
   logical              :: peOK         ! is pe available
   real(R8),allocatable :: minl(:)      ! local  min
   real(R8),allocatable :: ming(:)      ! global min
   real(R8),allocatable :: maxl(:)      ! local  max
   real(R8),allocatable :: maxg(:)      ! global max

   !--- formats ---
   character(*),parameter :: subName = '(mct_aVect_info) '
   character(*),parameter :: F00 = "('(mct_aVect_info) ',8a)"
   character(*),parameter :: F01 = "('(mct_aVect_info) ',a,i9)"
   character(*),parameter :: F02 = "('(mct_aVect_info) ',240a)"
   character(*),parameter :: F03 = "('(mct_aVect_info) ',a,2es11.3,i4,2x,a)"

!-------------------------------------------------------------------------------
! NOTE: has hard-coded knowledge/assumptions about mct aVect data type internals
!-------------------------------------------------------------------------------

   commOK = .false.
   peOK   = .false.

   if (present(pe)) then
     peOK = .true.
     pe_loc = pe
   endif
   if (present(comm)) then
     commOK = .true.
     comm_loc = comm
     if (.not.PEOK) then
       call shr_mpi_commrank(comm,pe_loc,subName)
       peOK = .true.
     endif
   endif

   nsize = mct_aVect_lsize(aVect)

   if (present(fld)) then
     nflds = 1
     ks = mct_aVect_indexRA(aVect,fld,perrWith=subName)
     ke = ks
   else
     nflds = mct_aVect_nRAttr(aVect)
     ks = 1
     ke = nflds
   endif

   if (flag >= 1) then
     if (present(istr)) then
        if (s_loglev > 0) write(s_logunit,*) trim(istr)
     endif
     if (s_loglev > 0) write(s_logunit,F01) "local size =",nsize
     if (associated(aVect%iList%bf)) then
        if (s_loglev > 0) write(s_logunit,F02) "iList = ",aVect%iList%bf
     endif
     if (associated(aVect%rList%bf)) then
        if (s_loglev > 0) write(s_logunit,F02) "rList = ",aVect%rList%bf
     endif
   endif

   if (flag >= 2) then

     allocate(minl(nflds))
     allocate(maxl(nflds))

     do k=ks,ke
       minl(k) = minval(aVect%rAttr(k,:))
       maxl(k) = maxval(aVect%rAttr(k,:))
     enddo

     if (flag >= 4 .and. commOK) then
       allocate(ming(nflds))
       allocate(maxg(nflds))
       ming = 0._R8
       maxg = 0._R8
       call shr_mpi_min(minl,ming,comm,subName)
       call shr_mpi_max(maxl,maxg,comm,subName)
     endif

     do k=ks,ke
       call mct_aVect_getRList(item,k,aVect)
       itemc = mct_string_toChar(item)
       call mct_string_clean(item)
       if (s_loglev > 0) write(s_logunit,F03) 'l min/max ',minl(k),maxl(k),k,trim(itemc)
       if (flag >= 3 .and. commOK) then
         if ((peOK .and. pe_loc == 0) .or. .not.peOK) then
           if (s_loglev > 0) write(s_logunit,F03) 'g min/max ',ming(k),maxg(k),k,trim(itemc)
         endif
       endif
       if (flag >= 4 .and. commOK) then
         if ((peOK .and. pe_loc == 0) .or. .not.peOK) then
           if (s_loglev > 0) write(s_logunit,*) trim(subName),'g min/max ',ming(k),maxg(k),k,trim(itemc)
         endif
       endif
     enddo

      deallocate(minl)
      deallocate(maxl)
      if (flag >= 4 .and. commOK) then
         deallocate(ming)
         deallocate(maxg)
      endif

   endif

   call shr_sys_flush(s_logunit)

end subroutine mct_aVect_info

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: mct_aVect_fldIndex - get a real fld index from an AVect
!
! !DESCRIPTION:
!     Get the field index for a real field in an attribute vector.
!     This is like mct_aVect_indexRA but with a calling interface
!     that returns the index without any error messages.
!
! !REMARKS:
!   This is like the MCT routine indexRA
!
! !REVISION HISTORY:
!     2010 Oct 27 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

integer function mct_aVect_fldIndex(aVect,fld)

! !INPUT/OUTPUT PARAMETERS:

   type(mct_aVect),intent(in)  :: aVect    ! an Attribute vector
   character(*)   ,intent(in)  :: fld      ! field name string

!EOP

   !--- local ---

   !--- formats ---
   character(*),parameter :: subName = "(mct_aVect_fldIndex) "
   character(*),parameter :: F00 = "('(mct_aVect_fldIndex) ',8a)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   mct_aVect_fldIndex = mct_aVect_indexRA(aVect,trim(fld),perrWith='quiet')

end function mct_aVect_fldIndex

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: mct_aVect_sharedFields - get a shared real fld index from two AVects
!
! !DESCRIPTION:
!     Get the shared field index for a real field in two attribute vectors.
!
! !REMARKS:
!
! !REVISION HISTORY:
!     2013 Jul 17 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine mct_aVect_sharedFields(aVect1, aVect2, rlistout, ilistout)

! !INPUT/OUTPUT PARAMETERS:

   type(mct_aVect),intent(in)  :: aVect1   ! an Attribute vector
   type(mct_aVect),intent(in)  :: aVect2   ! an Attribute vector
   character(*)   ,intent(inout),optional  :: rlistout      ! field name string
   character(*)   ,intent(inout),optional  :: ilistout      ! field name string

!EOP

   !--- local ---
   integer(IN) :: nflds1,nflds2
   character(len=CXX) :: list1,list2

   !--- formats ---
   character(*),parameter :: subName = "(mct_aVect_sharedFields) "
   character(*),parameter :: F00 = "('(mct_aVect_sharedFields) ',8a)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   if (present(rlistout)) then
      nflds1 = mct_aVect_nRAttr(aVect1)
      nflds2 = mct_aVect_nRAttr(aVect2)
      rlistout = ''
      list1 = ''
      list2 = ''
      if (nflds1 > 0 .and. nflds2 > 0) then
         list1 = mct_aVect_exportRList2c(aVect1)
         list2 = mct_aVect_exportRlist2c(aVect2)
         call shr_string_listIntersect(list1,list2,rlistout)
      endif
   endif

   if (present(ilistout)) then
      nflds1 = mct_aVect_nIAttr(aVect1)
      nflds2 = mct_aVect_nIAttr(aVect2)
      ilistout = ''
      list1 = ''
      list2 = ''
      if (nflds1 > 0 .and. nflds2 > 0) then
         list1 = mct_aVect_exportIList2c(aVect1)
         list2 = mct_aVect_exportIlist2c(aVect2)
         call shr_string_listIntersect(list1,list2,ilistout)
      endif
   endif

end subroutine mct_aVect_sharedFields

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: mct_aVect_initSharedFields - init new AVect based on shared fields 
!     from two input aVects
!
! !DESCRIPTION:
!     Init new AVect based on shared fields of two input AVects
!
! !REMARKS:
!
! !REVISION HISTORY:
!     2013 Jul 17 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine mct_aVect_initSharedFields(aVect1, aVect2, aVect3, lsize)

! !INPUT/OUTPUT PARAMETERS:

   type(mct_aVect),intent(in)  :: aVect1   ! an Attribute vector
   type(mct_aVect),intent(in)  :: aVect2   ! an Attribute vector
   type(mct_aVect),intent(inout)  :: aVect3   ! new Attribute vector
   integer(IN)    ,intent(in)  :: lsize    ! aVect3 size

!EOP

   !--- local ---
   character(len=CXX) :: rlist,ilist

   !--- formats ---
   character(*),parameter :: subName = "(mct_aVect_initSharedFields) "
   character(*),parameter :: F00 = "('(mct_aVect_initSharedFields) ',8a)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   call mct_aVect_sharedFields(aVect1,aVect2,rlist,ilist)
   call mct_aVect_init(aVect3,ilist,rlist,lsize)

end subroutine mct_aVect_initSharedFields

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: mct_aVect_getRAttr - get real F90 array data out of an aVect
!
! !DESCRIPTION:
!     Get the data associated with attribute {\tt str} in 
!     {\it AttributeVector} {\tt aVect} and return in the
!     real F90 array data {\tt data}.
!     {\tt rcode} will be 0 if succesful, 1 if size of {\tt data}
!     does not match size  of {\tt aVect} and 2 if {\tt str} is
!     not found.
!
! !REMARKS:
!   This is like the MCT routine exportRAttr except the output argument
!   is not a pointer.
!
! !REVISION HISTORY:
!     2002 Apr xx - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine mct_aVect_getRAttr(aVect,str,data,rcode)

! !INPUT/OUTPUT PARAMETERS:

   type(mct_aVect)    ,intent(in)  :: aVect    ! an Attribute vector
   character(*)       ,intent(in)  :: str      ! field name string
   real(R8)           ,intent(out) :: data(:)  ! an F90 array
   integer(IN)        ,intent(out) :: rcode    ! return code

!EOP

   !--- local ---
   integer(IN) :: k,n,m
   integer(IN) :: aVsize

   !--- formats ---
   character(*),parameter :: subName = "(mct_aVect_getRAttr) "
   character(*),parameter :: F00 = "('(mct_aVect_getRAttr) ',8a)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   rcode = 0

   n = mct_aVect_lsize(aVect)
   m = size(data)
   if (n /= m) then
      if (s_loglev > 0) write(s_logunit,*) subName,"ERROR: size aV,data,attr = ",n,m,trim(str)
      data = SHR_CONST_SPVAL
      rcode = 1
      return
   end if
   
   k = mct_aVect_indexRA(aVect,trim(str) ,perrWith=subName)
   if ( k < 1) then
      if (s_loglev > 0) write(s_logunit,*) subName,"ERROR: attribute not found, var = ",trim(str),", k=",k
      data = SHR_CONST_SPVAL
      rcode = 2
      return
   end if

   data(:) = aVect%rAttr(k,:)

end subroutine mct_aVect_getRAttr

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: mct_aVect_putRAttr - put real F90 array data into an aVect
!
! !DESCRIPTION:
!     Put the data in array {\tt data} into the  {\it AttributeVector}
!     {\tt aVect} under the attribute {\tt str}.
!     {\tt rcode} will be 0 if succesful, 1 if size of {\tt data}
!     does not match size  of {\tt aVect} and 2 if {\tt str} is not
!     found.
!
! !REMARKS:
!   This is like the MCT routine importRAttr except the output argument
!   is not a pointer.

! !REVISION HISTORY:
!     2002 Apr xx - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine mct_aVect_putRAttr(aVect,str,data,rcode)

! !INPUT/OUTPUT PARAMETERS:

   type(mct_aVect),intent(inout) :: aVect ! Attribute vector
   character(*)   ,intent(in)  :: str
   real(R8)       ,intent(in)  :: data(:)
   integer(IN)    ,intent(out) :: rcode

!EOP

   !--- local ---
   integer(IN) :: k,n,m
   integer(IN) :: aVsize

   !--- formats ---
   character(*),parameter :: subName = "(mct_aVect_putRAttr) "
   character(*),parameter :: F00 = "('(mct_aVect_putRAttr) ',8a)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   rcode = 0

   n = mct_aVect_lsize(aVect)
   m = size(data)
   if (n /= m) then
      if (s_loglev > 0) write(s_logunit,*) subName,"ERROR: size aV,data,attr = ",n,m,trim(str)
      rcode = 1
      return
   end if
   
   k = mct_aVect_indexRA(aVect,trim(str) ,perrWith=subName)
   if ( k < 1) then
      if (s_loglev > 0) write(s_logunit,*) subName,"ERROR: attribute not found, var = ",trim(str),", k=",k
      rcode = 2
      return
   end if

   aVect%rAttr(k,:) = data(:) 

end subroutine mct_aVect_putRAttr

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: mct_aVect_accum - accumulate attributes from one aVect to another
!
! !DESCRIPTION:
! This routine accumulates from input argment {\tt aVin} into the output 
! {\it AttrVect} argument {\tt aVout} the real and integer attributes specified in 
! input {\tt CHARACTER} argument {\tt iList} and {\tt rList}. The attributes can
! be listed in any order.  If neither {\tt iList} nor {\tt rList} are provided, 
! all attributes shared between {\tt aVin} and {\tt aVout} will be copied.
!
! If any attributes in {\tt aVout} have different names but represent the
! the same quantity and should still be copied, you must provide a translation
! argument {\tt TrList} and/or {\tt TiList}.  The translation arguments should
! be identical to the {\tt rList} or {\tt iList} but with the correct {\tt aVout}
! name subsititued at the appropriate place.
!
! This routine leverages the mct copy routines directly
!
! {\bf N.B.:}  This routine will fail if the {\tt aVout} is not initialized or
! if any of the specified attributes are not present in either {\tt aVout} or {\tt aVin}.
!
! !REVISION HISTORY:
!    2002 Sep 15 - ? - initial version.
!     2013-Jul-20 - T. Craig -- updated
!
! !INTERFACE: ------------------------------------------------------------------

 subroutine mct_avect_accum(aVin, aVout, rList, TrList, iList, TiList, vector, sharedIndices,counter)

      implicit none

! !INPUT PARAMETERS: 

      type(mct_avect),            intent(in)    :: aVin
      character(len=*), optional, intent(in)    :: iList
      character(len=*), optional, intent(in)    :: rList
      character(len=*), optional, intent(in)    :: TiList
      character(len=*), optional, intent(in)    :: TrList
      logical, optional,          intent(in)    :: vector 
      type(mct_avect_SharedIndices), optional, intent(in) :: sharedIndices

! !OUTPUT PARAMETERS: 

      type(mct_avect),         intent(inout) :: aVout
      integer, optional,       intent(inout) :: counter


! !REVISION HISTORY:

!EOP ___________________________________________________________________

   !--- local ---
  logical :: usevector
  integer(IN) :: lsize,nflds,npts,i,j
  type(mct_avect) :: avotmp  ! temporary aVout copy
  character(*),parameter :: subName = '(mct_aVect_accum) '

!-----------------------------------------------------------------

  usevector = .false.
  if (present(vector)) then
     usevector = vector
  endif

  if (present(counter)) then
     counter = counter + 1
  endif

  ! --- allocate avotmp, a duplciate of aVout

  lsize = mct_aVect_lsize(aVout)
  call mct_avect_init(avotmp,aVout,lsize)
  call mct_avect_zero(avotmp)

  ! --- copy aVin fields into avotmp

  if (present(sharedIndices)) then

     if (present(rList) .and. present(iList)) then
        if (present(trList) .and. present(tilist)) then
           call mct_avect_copy(aVin, avotmp, rList, TrList, iList, tiList, vector = usevector, sharedIndices=sharedIndices)
        elseif (present(trList)) then
           call mct_avect_copy(aVin, avotmp, rList, TrList, iList, vector = usevector, sharedIndices=sharedIndices)
        elseif (present(tiList)) then
           call mct_avect_copy(aVin, avotmp, rList, iList=iList, tiList=tiList, vector = usevector, sharedIndices=sharedIndices)
        else
           call mct_avect_copy(aVin, avotmp, rList=rList, iList=iList, vector = usevector, sharedIndices=sharedIndices)
        endif
     else if (present(rList)) then
        if (present(trList)) then
           call mct_avect_copy(aVin, avotmp, rList, TrList, vector = usevector, sharedIndices=sharedIndices)
        else
           call mct_avect_copy(aVin, avotmp, rList, vector = usevector, sharedIndices=sharedIndices)
        endif

     else if (present(iList)) then
        if (present(tiList)) then
           call mct_avect_copy(aVin, avotmp, ilist=iList, tiList=tiList, vector = usevector, sharedIndices=sharedIndices)
        else
           call mct_avect_copy(aVin, avotmp, ilist=iList, vector = usevector, sharedIndices=sharedIndices)
        endif

     else
        call mct_avect_copy(aVin, avotmp, vector=usevector, sharedIndices=sharedIndices)

     endif

  else   ! sharedIndices

     if (present(rList) .and. present(iList)) then
        if (present(trList) .and. present(tilist)) then
           call mct_avect_copy(aVin, avotmp, rList, TrList, iList, tiList, vector = usevector)
        elseif (present(trList)) then
           call mct_avect_copy(aVin, avotmp, rList, TrList, iList, vector = usevector)
        elseif (present(tiList)) then
           call mct_avect_copy(aVin, avotmp, rList, iList=iList, tiList=tiList, vector = usevector)
        else
           call mct_avect_copy(aVin, avotmp, rList=rList, iList=iList, vector = usevector)
        endif
     else if (present(rList)) then
        if (present(trList)) then
           call mct_avect_copy(aVin, avotmp, rList, TrList, vector = usevector)
        else
           call mct_avect_copy(aVin, avotmp, rList, vector = usevector)
        endif

     else if (present(iList)) then
        if (present(tiList)) then
           call mct_avect_copy(aVin, avotmp, ilist=iList, tiList=tiList, vector = usevector)
        else
           call mct_avect_copy(aVin, avotmp, ilist=iList, vector = usevector)
        endif

     else
        call mct_avect_copy(aVin, avotmp, vector=usevector)

     endif

  endif ! shared indices

  ! --- accumulate avotmp into avout

  nflds = mct_aVect_nRAttr(aVout)
  npts  = mct_aVect_lsize (aVout)
!DIR$ CONCURRENT
!DIR$ PREFERVECTOR
  do i=1,npts 
  do j=1,nflds
     aVout%rattr(j,i) = aVout%rattr(j,i) + avotmp%rattr(j,i)
  enddo
  enddo

  ! --- clean avotmp

  call mct_avect_clean(avotmp)

 end subroutine mct_avect_accum

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: mct_aVect_avg - averages an accumulated attribute vector
!
! !DESCRIPTION:
!     Average the data in attribute vector {\tt aVect}.  Divides all fields in 
!     the attribute vector {\tt aVect} by the value of the input counter.
!
! !REVISION HISTORY:
!     2002-Sep-15 - T. Craig -- initial version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine mct_aVect_avg(aVect, counter)

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   type(mct_aVect),intent(inout) :: aVect   ! bundle to read
   integer        ,intent(in)    :: counter ! counter 

!EOP

   !--- local ---
   integer(IN) :: i,j    ! generic indicies
   integer(IN) :: npts   ! number of points (local) in an aVect field
   integer(IN) :: nflds  ! number of aVect fields (real)
   real(R8)    :: ravg   ! accumulation count

   !--- formats ---
   character(*),parameter :: subName = '(mct_aVect_avg) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   if (counter == 0 .or. counter == 1) return

   ravg = 1.0_R8/real(counter,R8)

   nflds = mct_aVect_nRAttr(aVect)
   npts  = mct_aVect_lsize (aVect)
!DIR$ CONCURRENT
!DIR$ PREFERVECTOR
   do i=1,npts 
   do j=1,nflds
      aVect%rattr(j,i) = aVect%rattr(j,i)*ravg
   enddo
   enddo

end subroutine mct_aVect_avg

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: mct_avect_mult - multiply an attribute vector by a field.
!
! !DESCRIPTION:
!     Replace each field in {\tt av} by the product of that field and the
!     field {\tt fld1} from input argument {\tt av1}.
!
!     If optional argument {\tt bunlist} is present, only those attributes 
!     in {\tt bun} will be replaced.
!
!     If optional argument {\tt initav} is present, then the data in {\tt av}
!     is replaced by the product of the data in {\tt initav} and {\tt fld1}
!     from {\tt av1}. NOTE:  this assume {\tt initav} has the exact same
!     attributes in the same order as {\tt av}.
!
!
! !REVISION HISTORY:
!     2007-Jun-11 - M. Vertenstein -- initial version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine mct_avect_mult(av,av1,fld1,avlist)

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   type(mct_aVect)      ,intent(inout) :: av       ! attribute vector output
   type(mct_aVect)      ,intent(in)    :: av1      ! attribute vector input
   character(*)         ,intent(in)    :: fld1     ! av1 field name
   character(*),optional,intent(in)    :: avlist   ! sublist of field in av

!EOP

   !--- local ---
   integer(IN) :: n,m            ! generic indicies
   integer(IN) :: npts           ! number of points (local) in an aVect field
   integer(IN) :: nfld           ! number of fields (local) in an aVect field
   integer(IN) :: nfldi          ! number of fields (local) in an aVect field
   integer(IN) :: nptsx          ! number of points (local) in an aVect field
   integer(IN) :: nptsi          ! number of points (local) in an aVect field
   integer(IN) :: kfld           ! field number of fld1 in av1
   integer(IN),dimension(:),allocatable :: kfldin   ! field numbers of avlist in av
   type(mct_list)   :: blist     ! avlist as a List
   type(mct_string) :: tattr     ! an attribute

   !--- formats ---
   character(*),parameter :: subName = '(mct_aVect_mult) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   nptsx = mct_aVect_lsize(av1)
   npts  = mct_aVect_lsize(av)
   if (nptsx /= npts .and. s_loglev > 0) write(s_logunit,*) subName,' ERROR: npts error1 ',npts,nptsx

   kfld  = mct_aVect_indexRA(av1,fld1,perrWith=subName)

   if (present(avlist)) then

     call mct_list_init(blist,avlist)

     nfld=mct_list_nitem(blist)

     allocate(kfldin(nfld))
     do m=1,nfld
       call mct_list_get(tattr,m,blist)
       kfldin(m) = mct_aVect_indexRA(av,mct_string_toChar(tattr))
       call mct_string_clean(tattr)
     enddo
     call mct_list_clean(blist)

#ifdef CPP_VECTOR
     do m=1,nfld
!CDIR SELECT(VECTOR)
!DIR$ CONCURRENT
     do n=1,npts
#else
     do n=1,npts
     do m=1,nfld
#endif
        av%rAttr(kfldin(m),n) = av%rAttr(kfldin(m),n)*av1%rAttr(kfld,n)
     enddo
     enddo

     deallocate(kfldin)

   else

     nfld  = mct_aVect_nRAttr(av)

#ifdef CPP_VECTOR
     do m=1,nfld
!CDIR SELECT(VECTOR)
!DIR$ CONCURRENT
     do n=1,npts
#else
     do n=1,npts
     do m=1,nfld
#endif
        av%rAttr(m,n) = av%rAttr(m,n)*av1%rAttr(kfld,n)
     enddo
     enddo

   endif

end subroutine mct_aVect_mult

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: mct_avect_vecmult - multiply an attribute vector by a field.
!
! !DESCRIPTION:
!     Replace each field in {\tt av} by the product of that field and the
!     field {\tt fld1} from input argument {\tt av1}.
!
!     If optional argument {\tt bunlist} is present, only those attributes 
!     in {\tt bun} will be replaced.
!
!     If optional argument {\tt initav} is present, then the data in {\tt av}
!     is replaced by the product of the data in {\tt initav} and {\tt fld1}
!     from {\tt av1}. NOTE:  this assume {\tt initav} has the exact same
!     attributes in the same order as {\tt av}.
!
!
! !REVISION HISTORY:
!     2007-Jun-11 - M. Vertenstein -- initial version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine mct_avect_vecmult(av,vec,avlist)

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   type(mct_aVect)      ,intent(inout) :: av       ! attribute vector output
   real(R8)             ,intent(in)    :: vec(:)
   character(*),optional,intent(in)    :: avlist   ! sublist of field in av

!EOP

   !--- local ---
   integer(IN) :: n,m            ! generic indicies
   integer(IN) :: npts           ! number of points (local) in an aVect field
   integer(IN) :: nfld           ! number of fields (local) in an aVect field
   integer(IN) :: nfldi          ! number of fields (local) in an aVect field
   integer(IN) :: nptsx          ! number of points (local) in an aVect field
   integer(IN) :: nptsi          ! number of points (local) in an aVect field
   integer(IN),dimension(:),allocatable :: kfldin   ! field numbers of avlist in av
   type(mct_list)   :: blist     ! avlist as a List
   type(mct_string) :: tattr     ! an attribute

   !--- formats ---
   character(*),parameter :: subName = '(mct_aVect_vecmult) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   nptsx = size(vec,1)
   npts  = mct_aVect_lsize(av)
   if (nptsx /= npts .and. s_loglev > 0) write(s_logunit,*) subName,' ERROR: npts error1 ',npts,nptsx


   if (present(avlist)) then

     call mct_list_init(blist,avlist)

     nfld=mct_list_nitem(blist)

     allocate(kfldin(nfld))
     do m=1,nfld
       call mct_list_get(tattr,m,blist)
       kfldin(m) = mct_aVect_indexRA(av,mct_string_toChar(tattr))
       call mct_string_clean(tattr)
     enddo
     call mct_list_clean(blist)

#ifdef CPP_VECTOR
     do m=1,nfld
!CDIR SELECT(VECTOR)
!DIR$ CONCURRENT
     do n=1,npts
#else
     do n=1,npts
     do m=1,nfld
#endif
        av%rAttr(kfldin(m),n) = av%rAttr(kfldin(m),n)*vec(n)
     enddo
     enddo

     deallocate(kfldin)

   else

     nfld  = mct_aVect_nRAttr(av)

#ifdef CPP_VECTOR
     do m=1,nfld
!CDIR SELECT(VECTOR)
!DIR$ CONCURRENT
     do n=1,npts
#else
     do n=1,npts
     do m=1,nfld
#endif
        av%rAttr(m,n) = av%rAttr(m,n)*vec(n)
     enddo
     enddo

   endif

end subroutine mct_aVect_vecmult

!===============================================================================
! !BOP ===========================================================================
!
! !IROUTINE:  subroutine mct_rearr_rearrange_fldlst - rearrange on a fieldlist
!
! !DESCRIPTION: 
!     Perform regarranger between two attribute vectors only on the fieldlist
!     that is provided
!
!
! !REVISION HISTORY: 
!    2007-Jun-22 - M. Vertenstein - first version
! 
! !INTERFACE:  -----------------------------------------------------------------

subroutine mct_rearr_rearrange_fldlist(avi, avo, Rearr, vector, alltoall, fldlist, tag)

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   type(mct_aVect) , intent(in)  :: avi
   type(mct_aVect) , intent(inout):: avo
   type(mct_rearr) , intent(in)  :: Rearr
   logical         , intent(in)  :: vector
   logical         , intent(in)  :: alltoall
   character(len=*), intent(in)  :: fldlist   
   integer(IN)     , intent(in),optional :: tag
! !EOP

   !---local ---
   type(mct_aVect) :: avi_fl
   type(mct_aVect) :: avo_fl
   integer(IN)     :: lsize
   integer(IN)     :: ltag

   !--- formats ---
   character(*),parameter :: subName = '(mct_rearr_rearrange_fldlist) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   if (present(tag)) then
      ltag = tag
   else
      ltag = 3000
   endif

   lsize = mct_aVect_lsize(avi)
   call mct_aVect_init (avi_fl, rlist=fldlist, lsize=lsize)
   call mct_aVect_zero (avi_fl)

   lsize = mct_aVect_lsize(avo)
   call mct_aVect_init (avo_fl, rlist=fldlist, lsize=lsize)
   call mct_aVect_zero (avo_fl)
   
   call mct_aVect_copy (aVin=avi, aVout=avi_fl)
   call mct_rearr_rearrange(avi_fl, avo_fl, Rearr, VECTOR=vector, ALLTOALL=alltoall, tag=ltag)
   call mct_aVect_copy (aVin=avo_fl, aVout=avo, vector=vector)

   call mct_aVect_clean(avi_fl)
   call mct_aVect_clean(avo_fl)

end subroutine mct_rearr_rearrange_fldlist

!=======================================================================
logical function mct_gsmap_Identical(gsmap1,gsmap2)

  implicit none
  type(mct_gsMap), intent(IN):: gsmap1
  type(mct_gsMap), intent(IN):: gsmap2

  ! Local variables

  character(len=*),parameter :: subname = "(mct_gsmap_Identical) "
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

  mct_gsmap_Identical = identical

end function mct_gsmap_Identical
    
!===============================================================================
! !BOP ===========================================================================
!
! !IROUTINE:  mct_myindex - binary search for index in list
!
! !DESCRIPTION: 
!     Do a binary search to see if a value is contained in a list of
!     values.  return true or false.  starti must be monotonically
!     increasing, function does NOT check this.
!
!
! !REVISION HISTORY: 
!    2007-Jan-17 - T. Craig -- first version
!    2007-Mar-20 - R. Jacob - move to mct_mod
! 
! !INTERFACE:  -----------------------------------------------------------------

logical function mct_myindex(index,starti,counti)

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   integer(IN) :: index       ! is this index in start/count list
   integer(IN) :: starti(:)   ! start list
   integer(IN) :: counti(:)   ! count list

! !EOP

   !--- local ---
   integer(IN)    :: nl,nc,nr,ncprev 
   integer(IN)    :: lsize
   logical        :: stopnow

   !--- formats ---
   character(*),parameter :: subName = '(mct_myindex) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   mct_myindex = .false.

   lsize = size(starti)
   if (lsize < 1) return

   nl = 0
   nr = lsize + 1
   nc = (nl+nr)/2
   stopnow = .false.
   do while (.not.stopnow)
      if (index < starti(nc)) then
         nr = nc
      elseif (index > (starti(nc) + counti(nc) - 1)) then
         nl = nc
      else
         mct_myindex = .true.
         return
      endif
      ncprev = nc
      nc = (nl + nr)/2
      if (nc == ncprev .or. nc < 1 .or. nc > lsize) stopnow = .true.
   enddo

   mct_myindex = .false.
   return

end function mct_myindex
!===============================================================================

end module mct_mod

