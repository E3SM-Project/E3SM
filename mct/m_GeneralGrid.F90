!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_GeneralGrid -- General Grid representation class and methods.
!
! !DESCRIPTION:
! The GeneralGrid data type is the representation of a generalized 
! coordinate grid.  The number of dimensions of the grid is limited 
! only by storage capacity.  The grid is representated by a literal 
! listing of the gridpoints, which makes this datatype suitable for 
! anything up to an unstructured grid.  Consider an example 
! {\tt GeneralGrid} variable {\tt GGrid}.  The coordinate names for the 
! grid are stored in the {\tt List} component {\tt GGrid\%coordinate\_list}.
! Multidimensional area/volume weight names are stored in the {\tt List} 
! component {\tt GGrid\%weight\_list}.  The names of global indices for 
! the gridpoints are stored in the {\tt List} ! component 
! {\tt GGrid\%index\_list}.  The {\tt INTEGER} index data and the {\tt REAL} 
! point location and area/volume weight data are stored in the {\tt AttrVect} 
! component {\tt GGrid\%data}.  The list of {\tt REAL} attributes for 
! {\tt GGrid\%data} (that is, {\tt GGrid\%data\%rList}) is formed by 
! concatenating {\tt GGrid\%coordinate\_list} and {\tt GGrid\%weight\_list} 
! (with the coordinates before the weights).  The list of  The list of 
! {\tt INTEGER} attributes for {\tt GGrid\%data} (that is, 
! {\tt GGrid\%data\%iList}) is defined by {\tt GGrid\%index\_list}, along 
! with the ever-present {\em global gridpoint number attribute} {\tt GlobGridNum}.
!
! !INTERFACE:

 module m_GeneralGrid

!
! !USES:
!
      use m_List, only : List   ! Support for List components.

      use m_AttrVect, only : AttrVect ! Support for AttrVect component.

      implicit none

      private   ! except

      public :: GeneralGrid  ! The class data structure

      public :: init         ! Create a GeneralGrid
      public :: clean        ! Destroy a GeneralGrid

                             ! Query functions-----------------
      public :: dims         ! Return dimensionality of the GeneralGrid
      public :: indexIA      ! Index integer attribute (indices)
      public :: indexRA      ! Index integer attribute (coords/weights)
      public :: lsize        ! Return local number of points

                             ! Manipulation--------------------
      public :: Sort         ! Sort point data by coordinates -> permutation
      public :: Permute      ! Rearrange point data using input permutation
      public :: SortPermute  ! Sort and Permute point data

    type GeneralGrid
      type(List)                     :: coordinate_list
      type(List)                     :: coordinate_sort_order
      logical, dimension(:), pointer :: descend
      type(List)                     :: weight_list
      type(List)                     :: other_list
      type(List)                     :: index_list
      type(AttrVect), pointer        :: data
    end type GeneralGrid

    interface init  ; module procedure &
	 init_, &
	 initgg_
    end interface
    interface clean ; module procedure clean_ ; end interface
    interface dims ; module procedure dims_ ; end interface
    interface indexIA ; module procedure indexIA_ ; end interface
    interface indexRA ; module procedure indexRA_ ; end interface
    interface lsize   ; module procedure lsize_   ; end interface
    interface Sort    ; module procedure Sort_    ; end interface
    interface Permute ; module procedure Permute_ ; end interface
    interface SortPermute ; module procedure SortPermute_ ; end interface

! !REVISION HISTORY:
!       25Sep00 - J.W. Larson <larson@mcs.anl.gov> - initial prototype
!       31Oct00 - J.W. Larson <larson@mcs.anl.gov> - modified the
!                 GeneralGrid type to allow inclusion of grid cell
!                 dimensions (lengths) and area/volume weights.
!       15Jan01 - J.W. Larson implemented new GeneralGrid type 
!                 definition and added numerous APIs.
!       17Jan01 - J.W. Larson fixed minor bug in module header use
!                 statement.
!       19Jan01 - J.W. Larson added other_list and coordinate_sort_order
!                 components to the GeneralGrid type.
!       21Mar01 - J.W. Larson - deleted the initv_ API (more study
!                 needed before implementation.
!       02May01 - J.W. Larson - added initgg_ API (replaces old initv_).
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_GeneralGrid'

 contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_ - initialize an empty GeneralGrid (no data stored).
!
! !DESCRIPTION:
! The routine {\tt init\_()} creates the storage space for grid point
! coordinates, area/volume weights, and other coordinate data ({\em e.g.}, 
! nearest-neighbor coordinates).  These data are referenced by {\tt List}
! components that are also created by this routine.  Finally, a {\tt List}
! defining coordinate sorting order is created, along with a set of flags
! that dictate whether each sorting is done in ascending or descending 
! order.
!
! !INTERFACE:

 subroutine init_(GGrid, CoordList, CoordSortOrder, descend, WeightList, &
                  OtherList, IndexList, lsize )
!
! !USES:
!
      use m_die

      use m_String,   only : String, char
      use m_List,     only : List
      use m_List,     only : List_init => init
      use m_List,     only : List_nitem => nitem

      use m_AttrVect, only : AttrVect
      use m_AttrVect, only : AttrVect_init => init

      implicit none

! !INPUT PARAMETERS:
!
      character(len=*),             intent(in) :: CoordList
      character(len=*), optional,   intent(in) :: CoordSortOrder
      character(len=*), optional,   intent(in) :: WeightList
      logical, dimension(:), optional, pointer :: descend
      character(len=*), optional,   intent(in) :: OtherList
      character(len=*), optional,   intent(in) :: IndexList
      integer,          optional,   intent(in) :: lsize

! !OUTPUT PARAMETERS:
!
      type(GeneralGrid), intent(out)   :: GGrid

! !REVISION HISTORY:
!       25Sep00 - Jay Larson <larson@mcs.anl.gov> - initial prototype
!       15Jan01 - Jay Larson <larson@mcs.anl.gov> - modified to fit
!                 new GeneralGrid definition.  
!       19Mar01 - Jay Larson <larson@mcs.anl.gov> - added OtherList
!       25Apr01 - Jay Larson <larson@mcs.anl.gov> - added GlobGridNum
!                 as a mandatory integer attribute.
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//'::init_'

  character*128 :: RAList, IAList
  integer :: n, list_len, ier, RA_len, IA_len
  integer :: ierr, ncoord, nsort

        ! if lsize is present, use it to set n; if not, set n=0

  n = 0
  if(present(lsize)) n=lsize

        ! Concatenate onto CoordList the values of WeightList and
        ! OtherList (if present)--store the result in RAList

  list_len = len(CoordList)

  if(present(WeightList)) then
     list_len = list_len + 1 + len(WeightList)
  endif

  if(present(OtherList)) then
     list_len = list_len + 1 + len(OtherList)
  endif

  RAList = trim(CoordList)
  RA_len = len(CoordList)

  if(present(WeightList)) then
     RAList = RAList(1:RA_len) // ':' // trim(WeightList)
     RA_len = RA_len + 1 + len(WeightList)
  endif

  if(present(OtherList)) then
     RAList = RAList(1:RA_len) // ':' // trim(OtherList)
     RA_len = RA_len + 1 + len(OtherList)
  endif

        ! Create the integer Attribute list IAList.  At a minimum,
        ! IAList contains the attribute GlobGridNum

  if(present(IndexList)) then
     IAList = trim(IndexList) // ':GlobGridNum'
  else
     IAList = 'GlobGridNum'
  endif

        ! Initialize the AttrVect storage component GGrid%data

  if(present(IndexList)) then
     call AttrVect_init(aV=GGrid%data, iList=trim(IAList), &
	                rList=trim(RAList), lsize=n)
  else
     call AttrVect_init(aV=GGrid%data, rList=trim(RAList), &
	                lsize=n)
  endif

        ! Initialize GGrid%coordinate_list

  call List_init(GGrid%coordinate_list, CoordList)

        ! If Initialize GGrid%weight_list

  if(present(WeightList)) then
     call List_init(GGrid%weight_list, WeightList)
  else
     call List_init(GGrid%weight_list, ' ')
  endif

        ! If Initialize GGrid%other_list

  if(present(OtherList)) then
     call List_init(GGrid%other_list, OtherList)
  else
     call List_init(GGrid%other_list, ' ')
  endif

        ! If Initialize GGrid%coordinate_sort_order.  Check 
        ! the string CoordSortOrder for validity.

  if(present(CoordSortOrder)) then
     call List_init(GGrid%coordinate_sort_order, CoordSortOrder)
  endif

        ! Check the number of coordinates versus the number of
        ! coordinate grid sort keys...they should be equal.

     ncoord = List_nitem(GGrid%coordinate_list)
     nsort = List_nitem(GGrid%coordinate_sort_order)

  if(ncoord /= nsort) then
     call perr_die(myname_,'ncoord-nsort /=0',ncoord-nsort)
  endif

        ! If the LOGICAL argument descend is present, check the
        ! number of entries to ensure they match the grid dimensionality.
        ! If descend is not present, assume all coordinate grid point
        ! sortings will be in ascending order.

  if(present(descend)) then
     if(size(descend) /= nsort) then
	call perr_die(myname_,'size(descend)-nsort /=0',&
	              size(descend)-nsort)
     endif
  endif

  allocate(GGrid%descend(nsort), stat=ierr)
  if(ierr /= 0) then
     call MP_perr_die(myname_,"allocate(GGrid%descend...",ierr)
  endif

  if(present(descend)) then
     do n=1,nsort
	GGrid%descend(n) = descend(n)
     end do
  else
     do n=1,nsort
	GGrid%descend(n) = .false.
     end do
  endif

 end subroutine init_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: initgg_ - initialize GeneralGrid from another GeneralGrid.
!
! !DESCRIPTION:
! The routine {\tt initgg\_()} creates the storage space for grid point
! coordinates, area/volume weights, and other coordinate data ({\em e.g.}, 
! nearest-neighbor coordinates).  These data are all copied from the 
! already initialized input {\tt GeneralGrid} argument {\tt iGGrid}.  This 
! routine initializes the output {\tt GeneralGrid} argument {\tt oGGrid} 
! with the same {\tt List} data as {\tt iGGrid}, but with storage space
! for {\tt lsize} gridpoints.
!
! {\bf N.B.}:  It is assumed that {\tt iGGrid} has been initialized.
!
! {\bf N.B.}:  The output {\tt GeneralGrid} {\tt oGGrid} is dynamically
! allocated memory.  When one no longer needs {\tt oGGrid}, one should
! release this space by invoking {\tt GeneralGrid\_clean()}.
!
! !INTERFACE:

 subroutine initgg_(oGGrid, iGGrid, lsize)
!
! !USES:
!
      use m_die,  only : MP_perr_die

      use m_List, only : List
      use m_List, only : assignment(=)

      use m_AttrVect, only:  AttrVect_init => init

      implicit none

! !INPUT PARAMETERS: 
!
      type(GeneralGrid), intent(in)  :: iGGrid
      integer, optional, intent(in)  :: lsize

! !OUTPUT PARAMETERS:
!
      type(GeneralGrid), intent(out) :: oGGrid

! !REVISION HISTORY:
!       02May01 - Jay Larson <larson@mcs.anl.gov> - Initial version.
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//'::initgg_'
! Number of grid points, number of grid dimensions
  integer :: n, ncoord
! Loop index and Error Flag
  integer :: i, ierr

        ! if lsize is present, use it to set n; if not, set n=0

  n = 0
  if(present(lsize)) n=lsize

        ! Dimensionality of the GeneralGrid

  ncoord = dims_(iGGrid)

        ! Argument check:

  if(associated(iGGrid%descend)) then
     if(size(iGGrid%descend) /= ncoord) then ! size mismatch
	call MP_perr_die(myname,"iGGrid dims/size of descend mismatch", &
	                 ncoord-size(iGGrid%descend))
     endif
  endif

        ! If iGGrid%descend has been allocated, copy its contents

  if(associated(iGGrid%descend)) then ! allocate and fill oGGrid%descend

     allocate(oGGrid%descend(ncoord), stat=ierr)
     if(ierr /= 0) then
	call MP_perr_die(myname,"allocate(oGGrid%descend...", ierr)
     endif

     do i=1,ncoord
	oGGrid%descend(i) = iGGrid%descend(i)
     end do

  endif

       ! Copy list data from iGGrid to oGGrid:

  oGGrid%coordinate_list = iGGrid%coordinate_list
  oGGrid%coordinate_sort_order = iGGrid%coordinate_sort_order
  oGGrid%weight_list = iGGrid%weight_list
  oGGrid%other_list = iGGrid%other_list
  oGGrid%index_list = iGGrid%index_list

       ! Now, initialize oGGrid%data from iGGrid%data, but 
       ! with length n.

  call AttrVect_init(oGGrid%data, iGGrid%data, n)

 end subroutine initgg_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - Destroy a GeneralGrid.
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean_(GGrid)
!
! !USES:
!
      use m_die

      use m_List,only : List_clean => clean
      use m_AttrVect,only : AttrVect_clean => clean

      implicit none

      type(GeneralGrid), intent(inout) :: GGrid

! !REVISION HISTORY:
!       25Sep00 - J.W. Larson <larson@mcs.anl.gov> - initial prototype
!       20Mar01 - J.W. Larson <larson@mcs.anl.gov> - complete version.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'
  integer :: ierr

  call AttrVect_clean(GGrid%data)

  call List_clean(GGrid%coordinate_list)
  call List_clean(GGrid%coordinate_sort_order)
  call List_clean(GGrid%weight_list)
  call List_clean(GGrid%other_list)
  call List_clean(GGrid%index_list)

  deallocate(GGrid%descend, stat=ierr)
  if(ierr /= 0) then
     call MP_perr_die(myname_,"deallocate(GGrid%descend...",ierr)
  endif
  
 end subroutine clean_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: dims_ -- Returns the dimensionality of the grid.
!
! !DESCRIPTION:
!
! !INTERFACE:

 integer function dims_(GGrid)
!
! !USES:
!
      use m_List,     only : List_nitem => nitem
      use m_die

      implicit none

      type(GeneralGrid), intent(in)  :: GGrid

! !REVISION HISTORY:
!       15Jan01 - Jay Larson <larson@mcs.anl.gov> - initial version
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//'::dims_'

 dims_ = List_nitem(GGrid%coordinate_list)

 end function dims_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: indexIA - return storage index of named member of index_list
!
! !DESCRIPTION:
!
! !INTERFACE:

 integer function indexIA_(GGrid, item, perrWith, dieWith)
!
! !USES:
!
      use m_AttrVect,     only : AttrVect_indexIA => indexIA
      use m_die
      use m_stdio,        only : stderr

      implicit none

      type(GeneralGrid), intent(in)  :: GGrid
      character(len=*),  intent(in)  :: item
      character(len=*), optional, intent(in) :: perrWith
      character(len=*), optional, intent(in) :: dieWith

! !REVISION HISTORY:
!       15Jan01 - Jay Larson <larson@mcs.anl.gov> - Initial version.
!EOP ___________________________________________________________________
!
 character(len=*),parameter :: myname_=myname//'::indexIA_'

 
 indexIA_ = AttrVect_indexIA(GGrid%data, item)

  if(indexIA_==0) then
     if(.not.present(dieWith)) then
        if(present(perrWith)) write(stderr,'(4a)') perrWith, &
          '" indexIA_() error, not found "',trim(item),'"'
     else
        write(stderr,'(4a)') dieWith,       &
         '" indexIA_() error, not found "',trim(item),'"'
        call die(dieWith)
     endif
  endif

 end function indexIA_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: indexRA - return storage index of named member of index_list
!
! !DESCRIPTION:
!
! !INTERFACE:

 integer function indexRA_(GGrid, item, perrWith, dieWith)
!
! !USES:
!
      use m_AttrVect,     only : AttrVect_indexRA => indexRA
      use m_die
      use m_stdio,        only : stderr

      implicit none

      type(GeneralGrid), intent(in)  :: GGrid
      character(len=*),  intent(in)  :: item
      character(len=*), optional, intent(in) :: perrWith
      character(len=*), optional, intent(in) :: dieWith

! !REVISION HISTORY:
!       15Jan01 - Jay Larson <larson@mcs.anl.gov> - Initial version.
!EOP ___________________________________________________________________
!
 character(len=*),parameter :: myname_=myname//'::indexRA_'

 
 indexRA_ = AttrVect_indexRA(GGrid%data, item)

  if(indexRA_==0) then
     if(.not.present(dieWith)) then
        if(present(perrWith)) write(stderr,'(4a)') perrWith, &
          '" indexRA_() error, not found "',trim(item),'"'
     else
        write(stderr,'(4a)') dieWith,       &
         '" indexRA_() error, not found "',trim(item),'"'
        call die(dieWith)
     endif
  endif

 end function indexRA_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: lsize - returns local length of coordinate data storage.
!
! !DESCRIPTION:
!
! !INTERFACE:

 integer function lsize_(GGrid)
!
! !USES:
!

      use m_AttrVect,     only : AttrVect_lsize => lsize

      implicit none

      type(GeneralGrid), intent(in)  :: GGrid

! !REVISION HISTORY:
!       15Jan01 - Jay Larson <larson@mcs.anl.gov> - Initial version.
!EOP ___________________________________________________________________
!
 character(len=*),parameter :: myname_=myname//'::lsize_'

 lsize_ = 0

 if(associated(GGrid%data%rAttr)) then
    lsize_ = AttrVect_lsize( GGrid%data )
 endif

 end function lsize_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: Sort_ - returns sort permutation defined by arbitrary keys.
!
! !DESCRIPTION:
! The subroutine {\tt Sort\_()} uses the list of keys present in the
! input  {\tt List} variable {\tt key\_List}.  This list of keys is 
! checked to ensure that {\em only} coordinate attributes are present 
! in the sorting keys, and that there are no redundant keys.  Once 
! checked, this list is used to find the appropriate real attributes
! referenced by the items in {\tt key\_list} ( that is, it identifies the 
! appropriate entries in {\tt GGrid\%data\%rList}), and then uses these 
! keys to generate a an output permutation {\tt perm} that will put
! the entries of the attribute vector {\tt GGrid\%data} in lexicographic 
! order as defined by {\tt key\_list} (the ordering in {\tt key\_list} 
! being from left to right.
!
! !INTERFACE:

 subroutine Sort_(GGrid, key_List, perm, descend)

!
! !USES:
!
      use m_die

      use m_AttrVect,     only : AttrVect_Sort => Sort
      use m_List,        only : List_nitem => nitem

      implicit none

      type(GeneralGrid), intent(in)  :: GGrid
      type(List),        intent(in)  :: key_list
      integer, dimension(:), pointer              :: perm
      logical, dimension(:), optional, intent(in) :: descend

! !REVISION HISTORY:
!       15Jan01 - Jay Larson <larson@mcs.anl.gov> - Initial version.
!       20Mar01 - Jay Larson <larson@mcs.anl.gov> - Final working version.
!EOP ___________________________________________________________________
!
 character(len=*),parameter :: myname_=myname//'::Sort_'
 logical, dimension(:), allocatable :: descending
 integer :: n, ierr

        ! Here is how we transmit the sort order keys stored
        ! in descending (if present):

  n = List_nitem(key_list)
  allocate(descending(n), stat=ierr)
  if(ierr /= 0) then
     call MP_perr_die(myname_,"allocate(descending...",ierr)
  endif

  if(present(descend)) then
     descending = descend
  else
     descending = .false.
  endif

        ! This is a straightforward call to AttrVect_Sort().
  
  call AttrVect_Sort(GGrid%data, key_list, perm, descending) 

        ! Clean up...

  deallocate(descending, stat=ierr)

  if(ierr /= 0) then
     call MP_perr_die(myname_,"deallocate(descending...",ierr)
  endif


 end subroutine Sort_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: Sortg_ - returns sort permutation based on GeneralGrid keys.
!
! !DESCRIPTION:
! The subroutine {\tt Sortg\_()} uses the list of sorting keys present in 
! the input  {\tt GeneralGrid} variable {\tt GGrid\%coordinate\_sort\_order}} 
! to create a sort permutation {\tt perm(:)}.  Sorting is either in ascending
! or descending order based on the entries of {\tt GGrid\%descend(:)}.
! The output index permutation is stored in the array {\tt perm(:)} that 
! will put the entries of the attribute vector {\tt GGrid\%data} in 
! lexicographic order as defined by {\tt GGrid\%coordinate\_sort\_order}. The 
! ordering in {\tt GGrid\%coordinate\_sort\_order} being from left to right.
!
! {\bf N.B.:}  This routine returnss an allocatable array perm(:).  This 
! allocated array must be deallocated when the user no longer needs it.  
! Failure to do so will cause a memory leak.
!
! !INTERFACE:

 subroutine Sortg_(GGrid, perm)

!
! !USES:
!
      implicit none

      type(GeneralGrid), intent(in)  :: GGrid
      integer, dimension(:), pointer              :: perm

! !REVISION HISTORY:
!       22Mar01 - Jay Larson <larson@mcs.anl.gov> - Initial version.
!EOP ___________________________________________________________________
!
 character(len=*),parameter :: myname_=myname//'::Sortg_'

        ! This is a straightforward call to the GGrid method Sort_():

     call Sort_(GGrid, GGrid%coordinate_sort_order, &
                             perm, GGrid%descend)

 end subroutine Sortg_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: Permute_ - apply input permutation to re-order grid points
!
! !DESCRIPTION:
! The subroutine {\tt Permute\_()} uses an input index permutation {\tt perm} 
! to re-order the coordinate data stored in the {\tt GeneralGrid} argument 
! {\tt GGrid}.  This permutation can be generated by the routine 
! {\tt Sort\_()} contained in this module.
!
! !INTERFACE:

 subroutine Permute_(GGrid, perm)

!
! !USES:
!

      use m_die ,          only : die
      use m_stdio ,        only : stderr

      use m_AttrVect,     only : AttrVect
      use m_AttrVect, only : AttrVect_Permute => Permute

      implicit none

      type(GeneralGrid),     intent(inout)   :: GGrid
      integer, dimension(:), intent(in)      :: perm

! !REVISION HISTORY:
!       15Jan01 - Jay Larson <larson@mcs.anl.gov> - API specification.
!       10Apr01 - Jay Larson <larson@mcs.anl.gov> - API modified, working
!                 code.
!EOP ___________________________________________________________________
!
 character(len=*),parameter :: myname_=myname//'::Permute_'

        ! This is a straightforward call to AttrVect_Permute:

  call AttrVect_Permute(GGrid%data, perm)

 end subroutine Permute_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: SortPermute_ - sort and re-order coordinate data
!
! !DESCRIPTION:
! The subroutine {\tt SortPermute\_()} uses the list of keys defined in 
! {\tt GGrid\%coordinate\_sort\_order} to create an index permutation 
! {\tt perm}, which is then applied to re-order the coordinate data stored 
! in the {\tt GeneralGrid} argument {\tt GGrid} (more specifically, the 
! gridpoint data stored in {\tt GGrid\%data}.  This permutation is generated  
! by the routine {\tt Sort\_()} contained in this module.  The permutation 
! is carried out by the routine {\tt Permute\_()} contained in this module.
!
! !INTERFACE:

 subroutine SortPermute_(GGrid)

!
! !USES:
!

      use m_die ,          only : die
      use m_stdio ,        only : stderr

      implicit none

      type(GeneralGrid),     intent(inout)   :: GGrid

! !REVISION HISTORY:
!       15Jan01 - Jay Larson <larson@mcs.anl.gov> - API specification.
!       10Apr01 - Jay Larson <larson@mcs.anl.gov> - API modified, working
!                 code.
!       13Apr01 - Jay Larson <larson@mcs.anl.gov> - Simplified API and
!                 code (Thanks to Tony Craig of NCAR for detecting the
!                 bug that inspired these changes).
!EOP ___________________________________________________________________
!
 character(len=*),parameter :: myname_=myname//'::SortPermute_'

 integer, dimension(:), pointer :: perm
 integer :: ierr

  call Sortg_(GGrid, perm)

  call Permute_(GGrid, perm)

! Clean up--deallocate temporary permutation array:

  deallocate(perm, stat=ierr)

 end subroutine SortPermute_

 end module m_GeneralGrid
















































