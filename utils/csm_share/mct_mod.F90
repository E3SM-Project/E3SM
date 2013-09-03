!===============================================================================
! SVN $Id: $
! SVN $URL: $
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


   use m_MCTWorld           ,only: mct_world_init         => init

   use m_AttrVect           ,only: mct_aVect              => AttrVect
   use m_AttrVect           ,only: mct_aVect_init         => init
   use m_AttrVect           ,only: mct_aVect_clean        => clean
   use m_AttrVect           ,only: mct_aVect_zero         => zero
   use m_AttrVect           ,only: mct_aVect_lsize        => lsize
   use m_AttrVect           ,only: mct_aVect_indexIA      => indexIA
   use m_AttrVect           ,only: mct_aVect_indexRA      => indexRA
   use m_AttrVect           ,only: mct_aVect_importRattr  => importRattr
   use m_AttrVect           ,only: mct_aVect_exportRattr  => exportRattr
   use m_AttrVect           ,only: mct_aVect_getIList     => getIList
   use m_AttrVect           ,only: mct_aVect_getRList     => getRList
   use m_AttrVect           ,only: mct_aVect_exportIList2c=> exportIListToChar
   use m_AttrVect           ,only: mct_aVect_exportRList2c=> exportRListToChar
   use m_AttrVect           ,only: mct_aVect_nIAttr       => nIAttr
   use m_AttrVect           ,only: mct_aVect_nRAttr       => nRAttr
   use m_AttrVect           ,only: mct_aVect_copy         => Copy
   use m_AttrVect           ,only: mct_aVect_permute      => Permute
   use m_AttrVect           ,only: mct_aVect_unpermute    => Unpermute
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
   use m_GlobalSegMap       ,only: mct_gsMap_orderedPoints=> OrderedPoints

   use m_Rearranger         ,only: mct_rearr              => Rearranger
   use m_Rearranger         ,only: mct_rearr_init         => init
   use m_Rearranger         ,only: mct_rearr_clean        => clean
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
   use m_SparseMatrixComms  ,only: mct_sMat_ScatterByRow  => ScatterByRow
   use m_SparseMatrixComms  ,only: mct_sMat_ScatterByCol  => ScatterByColumn
   use m_SparseMatrixPlus   ,only: mct_sMatP              => SparseMatrixPlus
   use m_SparseMatrixPlus   ,only: mct_sMatP_Init         => init
   use m_SparseMatrixPlus   ,only: mct_sMatP_Vecinit      => vecinit
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

!EOP

   !--- local kinds ---
   integer,parameter,private :: R8 = SHR_KIND_R8
   integer,parameter,private :: IN = SHR_KIND_IN
   integer,parameter,private :: CL = SHR_KIND_CL

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
     if (present(istr)) write(6,*) trim(istr)
     write(6,F01) "local size =",nsize
     if (associated(aVect%iList%bf)) write(6,F02) "iList = ",aVect%iList%bf
     if (associated(aVect%rList%bf)) write(6,F02) "rList = ",aVect%rList%bf
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
       write(6,F03) 'l min/max ',minl(k),maxl(k),k,trim(itemc)
       if (flag >= 3 .and. commOK) then
         if ((peOK .and. pe_loc == 0) .or. .not.peOK) then
           write(6,F03) 'g min/max ',ming(k),maxg(k),k,trim(itemc)
         endif
       endif
       if (flag >= 4 .and. commOK) then
         if ((peOK .and. pe_loc == 0) .or. .not.peOK) then
           write(6,*) trim(subName),'g min/max ',ming(k),maxg(k),k,trim(itemc)
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

   call shr_sys_flush(6)

end subroutine mct_aVect_info

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

   type(mct_aVect),intent(in)  :: aVect    ! an Attribute vector
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
      write(6,*) subName,"ERROR: size aV,data,attr = ",n,m,trim(str)
      data = SHR_CONST_SPVAL
      rcode = 1
      return
   end if
   
   k = mct_aVect_indexRA(aVect,trim(str) ,perrWith=subName)
   if ( k < 1) then
      write(6,*) subName,"ERROR: attribute not found, var = ",trim(str),", k=",k
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

   type(mct_aVect),intent(out) :: aVect ! Attribute vector
   character(*)       ,intent(in)  :: str
   real(R8)           ,intent(in)  :: data(:)
   integer(IN)        ,intent(out) :: rcode

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
      write(6,*) subName,"ERROR: size aV,data,attr = ",n,m,trim(str)
      rcode = 1
      return
   end if
   
   k = mct_aVect_indexRA(aVect,trim(str) ,perrWith=subName)
   if ( k < 1) then
      write(6,*) subName,"ERROR: attribute not found, var = ",trim(str),", k=",k
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
! {\bf N.B.:}  This routine will fail if the {\tt aVout} is not initialized or
! if any of the specified attributes are not present in either {\tt aVout} or {\tt aVin}.
!
! !REVISION HISTORY:
!    2002 Sep 15 - ? - initial version.
!
! !INTERFACE: ------------------------------------------------------------------

subroutine mct_aVect_accum(aVin, rList, TrList, iList, TiList, aVout)


! !USES:

   use m_die ,          only : die
   use m_stdio ,        only : stderr
   use m_String ,       only : String_toChar => toChar
   use m_String ,       only : String
   use m_String ,       only : String_init
   use m_String ,       only : String_clean => clean
   use m_List ,         only : List
   use m_List,          only : List_get => get
   use m_List,          only : List_nullify => nullify
   use m_List,          only : List_clean => clean
   use m_List,          only : init,nitem
   use m_AttrVect,      only : AttrVect
   use m_AttrVect,      only : lsize
   use m_AttrVect,      only : SharedAttrIndexList

   implicit none

!INPUT/OUTPUT PARAMETERS

   type(AttrVect)        ,intent(in)    :: aVin
   character(*), optional,intent(in)    :: iList
   character(*), optional,intent(in)    :: rList
   character(*), optional,intent(in)    :: TiList
   character(*), optional,intent(in)    :: TrList
   type(AttrVect)        ,intent(inout) :: aVout

!EOP

   !--- local ---
   type(List)   :: rcpList       !  The list of real attributes to accum
   type(List)   :: icpList       !  The list of integer attributes to accum
   type(List)   :: TrcpList      !  Names of output attributes corresponding to input
   type(List)   :: TicpList      !  Names of output attributes corresponding to input
   type(String) :: attr          !  an individual attribute
   type(String) :: attr2         !  an individual attribute
   integer(IN)  :: i,j           ! generic indicies
   integer(IN)  :: rcode         ! return code
   integer(IN)  :: inx,outx
   integer(IN)  :: num_indices   ! Overlapping attribute index number

   !--- Overlapping attribute index storage arrays: ---
   integer(IN), dimension(:), pointer :: aVinindices, aVoutindices

   character(7) :: data_flag      ! character variable used as data type flag

   !--- formats ---
   character(*),parameter :: myname_='mct_accum'

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   call List_nullify(rcpList)
   call List_nullify(icpList)
   call List_nullify(TrcpList)
   call List_nullify(TicpList)

   if (lsize(aVin) .ne. lsize(aVout)) then
      write(stderr,'(2a)')myname_, &
      'MCTERROR:  Input aV and output aV do not have the same size'
      write(stderr,*)myname_, &
      'MCTERROR: ',lsize(aVin),lsize(aVout)
      call die(myname_,'lsize check',rcode)
   endif

   !----------------------------------------------------------------------------
   ! Accum the listed real attributes
   !----------------------------------------------------------------------------
   if ( present(rList)) then
     if( len_trim(rList)>0 ) then

      call init(rcpList,rList)    ! init.List()

      !--- check translation list ---
      if ( present(TrList) ) then 
       if(len_trim(TrList)>0 ) then
         call init(TrcpList,TrList)
         if ( nitem(rcpList) .ne. nitem(TrcpList)) then
            write(stderr,'(2a)')myname_, &
            'MCTERROR:  Input rList and TrList do not have the same size'
            call die(myname_,'nitem TrList check',rcode)
         endif
       endif
      endif

      if (nitem(rcpList) .ge. 1) then
         do i=1,lsize(aVin)
         do j=1,nitem(rcpList)
            call List_get(attr,j,rcpList)
            if (present(TrList)) then
               call List_get(attr2,j,TrcpList)
            else
               call String_init(attr2,attr)
            endif
            inx=mct_aVect_indexRA(aVin,String_toChar(attr),dieWith=myname_//'real aVin')
            outx=mct_aVect_indexRA(aVout,String_toChar(attr2),dieWith=myname_//'real aVout')
            aVout%rAttr(outx,i)=aVout%rAttr(outx,i)+aVin%rAttr(inx,i)
            call String_clean(attr)
            call String_clean(attr2)
         enddo
         enddo
      endif

     call List_clean(rcpList)
     if (present(TrList)) call List_clean(TrcpList)

     endif
   endif

   !----------------------------------------------------------------------------
   ! Accum the listed integer attributes 
   !----------------------------------------------------------------------------
   if ( present(iList) ) then
     if (len_trim(iList)>0 ) then

      call init(icpList,iList)    ! init.List()

      !--- check translation list ---
      if ( present(TiList) ) then
        if (len_trim(TiList)>0 ) then
         call init(TicpList,TiList)
         if ( nitem(icpList) .ne. nitem(TicpList)) then
             write(stderr,'(2a)')myname_, &
            'MCTERROR:  Input iList and TiList do not have the same size'
            call die(myname_,'nitem TiList check',rcode)
         endif
        endif
      endif

     if (nitem(icpList) .ge. 1) then
        do i=1,lsize(aVin)
        do j=1,nitem(icpList)
           call List_get(attr,j,icpList)
           if (present(TiList)) then
              call List_get(attr2,j,TicpList)
           else
              call String_init(attr2,attr)
           endif
           inx =mct_aVect_indexIA(aVin ,String_toChar(attr) ,dieWith=myname_//'int aVin')
           outx=mct_aVect_indexIA(aVout,String_toChar(attr2),dieWith=myname_//'int aVout')
           aVout%iAttr(outx,i)=aVout%iAttr(outx,i)+aVin%iAttr(inx,i)
           call String_clean(attr)
           call String_clean(attr2)
        enddo
        enddo
     endif

     call List_clean(icpList)
     if (present(TrList)) call List_clean(TicpList)

     endif
   endif

   !----------------------------------------------------------------------------
   ! if neither rList nor iList is present, accum shared attibutes from in to out
   !----------------------------------------------------------------------------
   if ( .not.present(rList) .and. .not.present(iList)) then

      data_flag = 'REAL'
      call SharedAttrIndexList(aVin, aVout, data_flag, num_indices, &
                                aVinindices, aVoutindices)
      if (num_indices .gt. 0) then
#ifdef CPP_VECTOR
         do j=1,num_indices
!CDIR SELECT(VECTOR)
!DIR$ CONCURRENT
         do i=1,lsize(aVin)
#else
         do i=1,lsize(aVin)
         do j=1,num_indices
#endif
            aVout%rAttr(aVoutindices(j),i)= &
            & aVout%rAttr(aVoutindices(j),i)+aVin%rAttr(aVinindices(j),i)
         enddo
         enddo
      endif
      deallocate(aVinindices, aVoutindices,stat=rcode)
      if (rcode /= 0) call die(myname_,'deallocate real(Vinindices...',rcode)

      data_flag = 'INTEGER'
      call SharedAttrIndexList(aVin, aVout, data_flag, num_indices, &
                                aVinindices, aVoutindices)
      if (num_indices .gt. 0) then
#ifdef CPP_VECTOR
         do j=1,num_indices
!CDIR SELECT(VECTOR)
!DIR$ CONCURRENT
         do i=1,lsize(aVin)
#else
         do i=1,lsize(aVin)
         do j=1,num_indices
#endif
            aVout%iAttr(aVoutindices(j),i)= &
            & aVout%iAttr(aVoutindices(j),i)+aVin%iAttr(aVinindices(j),i)
         enddo
         enddo
      endif
      deallocate(aVinindices, aVoutindices,stat=rcode)
      if (rcode /= 0) call die(myname_,'deallocate int(Vinindices...',rcode)

   endif

end subroutine mct_aVect_accum

!===============================================================================

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: mct_sMat_readnc - read all mapping data from a NetCDF SCRIP file
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

subroutine mct_sMat_readnc(sMat,fileName)

#include <netcdf.inc>

! !INPUT/OUTPUT PARAMETERS:

   type(mct_sMat),intent(inout)  :: sMat
   character(*),intent(in)  :: filename  ! netCDF file to read

!EOP

   !--- local ---
   integer(IN)           :: n       ! generic loop indicies
   integer(IN)           :: na      ! size of source domain
   integer(IN)           :: nb      ! size of destination domain
   integer(IN)           :: ns      ! number of non-zero elements in matrix
   integer(IN)           :: ni,nj   ! number of row and col in the matrix
   integer(IN)           :: igrow   ! aVect index for matrix row
   integer(IN)           :: igcol   ! aVect index for matrix column
   integer(IN)           :: iwgt    ! aVect index for matrix element

   real(R8)   ,allocatable :: rtemp(:)  ! reals
   integer(IN),allocatable :: itemp(:)  ! ints

   integer(IN)           :: rcode   ! netCDF routine return code
   integer(IN)           :: fid     ! netCDF file      ID
   integer(IN)           :: vid     ! netCDF variable  ID
   integer(IN)           :: did     ! netCDF dimension ID

   character(*),parameter :: subName = '(mct_sMat_readnc) '
   character(*),parameter :: F00 = "('(mct_sMat_readnc) ',4a)"
   character(*),parameter :: F01 = '("(mct_sMat_readnc) ",2(a,i7))'

   write(6,F00) "reading mapping matrix data..."

   !----------------------------------------------------------------------------
   ! open & read the file
   !----------------------------------------------------------------------------
   write(6,F00) "* file name                  : ",trim(fileName)
   rcode = nf_open(filename,NF_NOWRITE,fid)
   if (rcode /= NF_NOERR) then
      write(6,F00) nf_strerror(rcode)
      call mct_die(subName,"error opening Netcdf file")
   endif

   !--- allocate memory & get matrix data ----------
   rcode = nf_inq_dimid (fid, 'n_s', did)  ! size of sparse matrix
   rcode = nf_inq_dimlen(fid, did  , ns)
   rcode = nf_inq_dimid (fid, 'n_a', did)  ! size of  input vector
   rcode = nf_inq_dimlen(fid, did  , na)
   rcode = nf_inq_dimid (fid, 'n_b', did)  ! size of output vector
   rcode = nf_inq_dimlen(fid, did  , nb)

   write(6,F01) "* matrix dimensions rows x cols :",na,' x',nb
   write(6,F01) "* number of non-zero elements: ",ns

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

   rcode = nf_inq_varid     (fid,'S'  ,vid)
   rcode = nf_get_var_double(fid,vid  ,rtemp  )
   if (rcode /= NF_NOERR) write(6,F00) nf_strerror(rcode)

   sMat%data%rAttr(iwgt ,:) =   rtemp(:)

   deallocate(rtemp, stat=rcode)
   if (rcode /= 0) call mct_perr_die(subName,':: deallocate weights',rcode)

   !!!!!!!!!!!!!!!!!!!!!!!!!!
   ! read and load rows
   allocate(itemp(ns),stat=rcode)
   if (rcode /= 0) call mct_perr_die(subName,':: allocate rows',rcode)

   rcode = nf_inq_varid     (fid,'row',vid)
   rcode = nf_get_var_int   (fid,vid  ,itemp)
   if (rcode /= NF_NOERR) write(6,F00) nf_strerror(rcode)

   sMat%data%iAttr(igrow,:) = itemp(:)


   !!!!!!!!!!!!!!!!!!!!!!!!!!
   ! read and load columns
   itemp(:) = 0

   rcode = nf_inq_varid     (fid,'col',vid)
   rcode = nf_get_var_int   (fid,vid  ,itemp)
   if (rcode /= NF_NOERR) write(6,F00) nf_strerror(rcode)

   sMat%data%iAttr(igcol,:) = itemp(:)

   deallocate(itemp, stat=rcode)
   if (rcode /= 0) call mct_perr_die(subName,':: deallocate cols',rcode)

   rcode = nf_close(fid)

   write(6,F00) "... done reading file"
   call shr_sys_flush(6)

end subroutine mct_sMat_readnc
!===============================================================================

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: mct_sMatP_initnc - initialize a SparseMatrixPlus.
!
! !DESCRIPTION:
!   Read in mapping matrix data from a SCRIP netCDF data file in first an
!   Smat and then an SMatPlus
!
! !REMARKS:
!
! !REVISION HISTORY:
!     2006 Nov 27: R. Jacob
!
! !INTERFACE: ------------------------------------------------------------------

subroutine mct_sMatP_initnc(sMatP,gsMapX,gsMapY,ConfigFileName, &
                           MapLabel,MapTypeLabel,mpicom)

! !INPUT/OUTPUT PARAMETERS:

   type(mct_sMatP),intent(inout)  :: sMatP
   type(mct_gsMap),intent(in)    :: gsMapX
   type(mct_gsMap),intent(in)    :: gsMapY
   character(*),intent(in)  :: ConfigFileName  ! config file to read
   character(*),intent(in)  :: MapLabel        ! map name
   character(*),intent(in)  :: MapTypeLabel    ! map type
   integer,intent(in)  :: mpicom

!EOP
   type(mct_sMat ) :: sMati    ! initial sMat from read (either root or decomp)
   type(mct_Avect) :: areasrc0 ! area of src grid from mapping file
   type(mct_Avect) :: areadst0 ! area of dst grid from mapping file

   integer :: iret
   integer :: pe_loc
   character(256)    :: fileName
   character(1) :: maptype
   character(*),parameter :: F00 = "('(mct_sMatP_initnc) ',4a)"

   call shr_mpi_commrank(mpicom,pe_loc)

   write(6,F00) "Initializing SparseMatrixPlus"
   call I90_allLoadF(ConfigFileName,0,mpicom,iret)
      if(iret /= 0) then
        write(6,*)"Cant find config file ",ConfigFileName
        call mct_die("mct_sMapP_initnc","File Not Found")
      endif

   call i90_label(trim(MapLabel),iret)
      if(iret /= 0) then
        write(6,*)"Cant find label ",MapLabel
        call mct_die("mct_sMapP_initnc","Label Not Found")
      endif

   call i90_gtoken(fileName,iret)
      if(iret /= 0) then
        write(6,*)"Error reading token ",MapLabel
        call mct_die("mct_sMapP_initnc","Error on read")
      endif
   write(6,F00) "map weight filename ",trim(fileName)

   call i90_label(trim(MapTypeLabel),iret)
   call i90_gtoken(maptype,iret)
   write(6,F00) "SmatP maptype ",maptype
   if(maptype == "X") then
      call mct_sMat_readdnc(sMati,gsMapX,gsMapY,'src', areasrc0, &
           areadst0, fileName, pe_loc, mpicom)
      call mct_sMatP_Init(sMatP,sMati,gsMapX,gsMapY, &
           0,mpicom,gsMapX%comp_id)
   else if(maptype == "Y") then
      call mct_sMat_readdnc(sMati,gsMapX,gsMapY,'dst', areasrc0, &
           areadst0, fileName, pe_loc,  mpicom)
      call mct_sMatP_Init(sMatP,sMati,gsMapX,gsMapY, &
          0,mpicom,gsMapY%comp_id)
   endif
#ifdef CPP_VECTOR
   !--- initialize the vector parts of the sMat ---
      call mct_sMatP_Vecinit(sMatP)
#endif


   write(6,F00) "Done initializing SmatP"



   call mct_sMat_Clean(sMati)
   call I90_Release(iret)

end subroutine mct_sMatP_initnc

!BOP ===========================================================================
!
! !IROUTINE:  mct_sMat_readdnc - Do a distributed read of a NetCDF SCRIP file and
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
!    2007-Mar-20 - R. Jacob -- rename to mct_sMat_readdnc.  Remove use of cpl_
!                  variables and move to mct_mod
! 
! !INTERFACE:  -----------------------------------------------------------------

subroutine mct_sMat_readdnc(sMat,SgsMap,DgsMap,newdom,areasrc,areadst, &
                            fileName,my_id, mycomm)

! !USES:

#include <netcdf.inc>

! !INPUT/OUTPUT PARAMETERS:

   type(mct_sMat)  ,intent(out)        :: sMat    ! mapping data
   type(mct_gsMap) ,intent(in) ,target :: SgsMap  ! src gsmap
   type(mct_gSMap) ,intent(in) ,target :: DgsMap  ! dst gsmap
   character(*)        ,intent(in)     :: newdom  ! type of sMat (src or dst)
   type(mct_Avect) ,intent(out)        :: areasrc ! area of src grid from mapping file
   type(mct_Avect) ,intent(out)        :: areadst ! area of dst grid from mapping file
   character(*)        ,intent(in)     :: filename! netCDF file to read
   integer(IN)        ,intent(in)      :: my_id   ! processor id
   integer(IN)        ,intent(in)      :: mycomm  ! communicator

! !EOP

   !--- local ---
   integer(IN)           :: n,m     ! generic loop indicies
   integer(IN)           :: na      ! size of source domain
   integer(IN)           :: nb      ! size of destination domain
   integer(IN)           :: ns      ! number of non-zero elements in matrix
   integer(IN)           :: ni,nj   ! number of row and col in the matrix
   integer(IN)           :: igrow   ! aVect index for matrix row
   integer(IN)           :: igcol   ! aVect index for matrix column
   integer(IN)           :: iwgt    ! aVect index for matrix element
 integer(IN)           :: iarea   ! aVect index for area
   integer(IN)           :: rsize   ! size of read buffer
   integer(IN)           :: cnt     ! local num of wgts
   integer(IN)           :: cntold  ! local num of wgts, previous read
   integer(IN)           :: start(1)! netcdf read
   integer(IN)           :: count(1)! netcdf read
   integer(IN)           :: bsize   ! buffer size
   integer(IN)           :: nread   ! number of reads 
   logical               :: mywt    ! does this weight belong on my pe

   !--- buffers for i/o ---
   real(R8)   ,allocatable :: rtemp(:) ! real temporary
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

   character,allocatable :: str(:)  ! variable length char string
   character(CL)         :: attstr  ! netCDF attribute name string
   integer(IN)           :: rcode   ! netCDF routine return code
   integer(IN)           :: fid     ! netCDF file      ID
   integer(IN)           :: vid     ! netCDF variable  ID
   integer(IN)           :: did     ! netCDF dimension ID
   !--- arbitrary size of read buffer, this is the chunk size weights reading
   integer(IN),parameter :: rbuf_size = 100000

   character(*),parameter :: areaAV_field = 'aream'

   !--- formats ---
   character(*),parameter :: subName = '(mct_sMat_readdnc) '
   character(*),parameter :: F00 = '("(mct_sMat_readdnc) ",4a)'
   character(*),parameter :: F01 = '("(mct_sMat_readdnc) ",2(a,i7))'

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

 call shr_mpi_commsize(mycomm,commsize)
 if (my_id == 0) then
   write(6,F00) "reading mapping matrix data decomposed..."

   !----------------------------------------------------------------------------
   ! open & read the file
   !----------------------------------------------------------------------------
   write(6,F00) "* file name                  : ",trim(fileName)
   rcode = nf_open(filename,NF_NOWRITE,fid)
   if (rcode /= NF_NOERR) write(6,F00) nf_strerror(rcode)

   !--- get matrix dimensions ----------
   rcode = nf_inq_dimid (fid, 'n_s', did)  ! size of sparse matrix
   rcode = nf_inq_dimlen(fid, did  , ns)
   rcode = nf_inq_dimid (fid, 'n_a', did)  ! size of  input vector
   rcode = nf_inq_dimlen(fid, did  , na)
   rcode = nf_inq_dimid (fid, 'n_b', did)  ! size of output vector
   rcode = nf_inq_dimlen(fid, did  , nb)

   write(6,F01) "* matrix dimensions rows x cols :",na,' x',nb
   write(6,F01) "* number of non-zero elements: ",ns

   !--- read and load area_a ---
   allocate(rtemp(na),stat=rcode)
   if (rcode /= 0) call mct_perr_die(subName,':: allocate area_a',rcode)

   rcode = nf_inq_varid     (fid,'area_a',vid)
   rcode = nf_get_var_double(fid,vid     ,rtemp)
   if (rcode /= NF_NOERR) write(6,F00) nf_strerror(rcode)

   call mct_aVect_init(areasrc,' ',areaAV_field,na)

   iarea = mct_aVect_indexRA(areasrc,areaAV_field)
   areasrc%rAttr(iarea,1:na) = rtemp(1:na)

   deallocate(rtemp, stat=rcode)
   if (rcode /= 0) call mct_perr_die(subName,':: deallocate area_a',rcode)

   !--- read and load area_b ---
   allocate(rtemp(nb),stat=rcode)
   if (rcode /= 0) call mct_perr_die(subName,':: allocate area_b',rcode)

   rcode = nf_inq_varid     (fid,'area_b',vid)
   rcode = nf_get_var_double(fid,vid     ,rtemp)
   if (rcode /= NF_NOERR) write(6,F00) nf_strerror(rcode)

   call mct_aVect_init(areadst,' ',areaAV_field,nb)

   iarea = mct_aVect_indexRA(areadst,areaAV_field)
   areadst%rAttr(iarea,1:nb) = rtemp(1:nb)

   deallocate(rtemp, stat=rcode)
   if (rcode /= 0) call mct_perr_die(subName,':: deallocate area_b',rcode)

  endif

   call shr_mpi_bcast(ns,mycomm,subName//" MPI in ns bcast")
   call shr_mpi_bcast(na,mycomm,subName//" MPI in na bcast")
   call shr_mpi_bcast(nb,mycomm,subName//" MPI in nb bcast")

   !--- setup local seg map, sorted
   if (newdom == 'src') then
      mygsmap => DgsMap
   elseif (newdom == 'dst') then
      mygsmap => SgsMap
   else
      write(6,F00) 'ERROR: invalid newdom value = ',newdom
      call shr_sys_abort(trim(subName)//" invalid newdom value")
   endif
   lsize = 0
   do n = 1,size(mygsmap%start)
      if (mygsmap%pe_loc(n) == my_id) then
         lsize=lsize+1
      endif
   enddo
   allocate(lsstart(lsize),lscount(lsize),stat=rcode)
   if (rcode /= 0) call mct_perr_die(subName,':: allocate Lsstart',rcode)

   lsize = 0
   do n = 1,size(mygsmap%start)
      if (mygsmap%pe_loc(n) == my_id) then  ! on my pe
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
        write(6,*) trim(subName),' lsstart not properly sorted'
         call shr_sys_abort()
      endif
   enddo

   rsize = min(rbuf_size,ns)                     ! size of i/o chunks
   bsize = ((ns/commsize) + 1 ) * 1.2   ! local temporary buffer size
   nread = (ns-1)/rsize + 1                      ! num of reads to do

   allocate(Sbuf(rsize),Rbuf(rsize),Cbuf(rsize),stat=rcode)
   if (rcode /= 0) call mct_perr_die(subName,':: allocate Sbuf',rcode)
   allocate(Snew(bsize),Cnew(bsize),Rnew(bsize),stat=rcode)
   if (rcode /= 0) call mct_perr_die(subName,':: allocate Snew1',rcode)

   cnt = 0
   do n = 1,nread
      start(1) = (n-1)*rsize + 1
      count(1) = min(rsize,ns-start(1)+1)

      !--- read data on root pe
      if (my_id== 0) then
         rcode = nf_inq_varid      (fid,'S'  ,vid)
         rcode = nf_get_vara_double(fid,vid,start,count,Sbuf)
         if (rcode /= NF_NOERR) write(6,F00) nf_strerror(rcode)

         rcode = nf_inq_varid      (fid,'row',vid)
         rcode = nf_get_vara_int   (fid,vid,start,count,Rbuf)
         if (rcode /= NF_NOERR) write(6,F00) nf_strerror(rcode)

         rcode = nf_inq_varid      (fid,'col',vid)
         rcode = nf_get_vara_int   (fid,vid,start,count,Cbuf)
         if (rcode /= NF_NOERR) write(6,F00) nf_strerror(rcode)
      endif

      !--- send S, row, col to all pes
      call shr_mpi_bcast(Sbuf,mycomm,subName//" MPI in Sbuf bcast")
      call shr_mpi_bcast(Rbuf,mycomm,subName//" MPI in Rbuf bcast")
      call shr_mpi_bcast(Cbuf,mycomm,subName//" MPI in Cbuf bcast")

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
               write(6,*) trim(subName),' reallocate bsize to ',bsize
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

   if (my_id == 0) then
      rcode = nf_close(fid)
      write(6,F00) "... done reading file"
      call shr_sys_flush(6)
   endif

end subroutine mct_sMat_readdnc


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

   if (counter == 0) return

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

end module mct_mod

