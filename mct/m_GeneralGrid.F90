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
! with the ever-present {\em global gridpoint number attribute} 
! {\tt GlobGridNum}.
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

! !PUBLIC TYPES:

    type GeneralGrid
      type(List)                     :: coordinate_list
      type(List)                     :: coordinate_sort_order
      logical, dimension(:), pointer :: descend
      type(List)                     :: weight_list
      type(List)                     :: other_list
      type(List)                     :: index_list
      type(AttrVect)                 :: data
    end type GeneralGrid

! !PUBLIC MEMBER FUNCTIONS:

      public :: GeneralGrid      ! The class data structure

      public :: init             ! Create a GeneralGrid
      public :: initCartesian    !
      public :: initUnstructured !
      public :: clean            ! Destroy a GeneralGrid

                             ! Query functions-----------------
      public :: dims         ! Return dimensionality of the GeneralGrid
      public :: indexIA      ! Index integer attribute (indices)
      public :: indexRA      ! Index integer attribute (coords/weights)
      public :: lsize        ! Return local number of points
      public :: exportIAttr  ! Return INTEGER attribute as a vector
      public :: exportRAttr  ! Return REAL attribute as a vector

                             ! Manipulation--------------------
      public :: importIAttr  ! Insert INTEGER vector as attribute
      public :: importRAttr  ! Insert REAL vector as attribute
      public :: Sort         ! Sort point data by coordinates -> permutation
      public :: Permute      ! Rearrange point data using input permutation
      public :: SortPermute  ! Sort and Permute point data

    interface init  ; module procedure &
	 init_, &
	 initl_, &
	 initgg_
    end interface
    interface initCartesian ; module procedure &
	 initCartesian_
    end interface
    interface initUnstructured ; module procedure &
	 initUnstructured_
    end interface
    interface clean ; module procedure clean_ ; end interface
    interface dims ; module procedure dims_ ; end interface
    interface indexIA ; module procedure indexIA_ ; end interface
    interface indexRA ; module procedure indexRA_ ; end interface
    interface lsize   ; module procedure lsize_   ; end interface
    interface exportIAttr ; module procedure exportIAttr_ ; end interface
    interface exportRAttr ; module procedure exportRAttr_ ; end interface
    interface importIAttr ; module procedure importIAttr_ ; end interface
    interface importRAttr ; module procedure importRAttr_ ; end interface
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
!       13Dec01 - J.W. Larson - added import and export methods.
!       27Mar02 - J.W. Larson <larson@mcs.anl.gov> - Corrected usage of
!                 m_die routines throughout this module.
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

 subroutine init_(GGrid, CoordChars, CoordSortOrder, descend, WeightChars, &
                  OtherChars, IndexChars, lsize )
!
! !USES:
!
      use m_stdio
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
      character(len=*),             intent(in) :: CoordChars
      character(len=*),             intent(in) :: CoordSortOrder
      character(len=*), optional,   intent(in) :: WeightChars
      logical, dimension(:), optional, pointer :: descend
      character(len=*), optional,   intent(in) :: OtherChars
      character(len=*), optional,   intent(in) :: IndexChars
      integer,          optional,   intent(in) :: lsize

! !OUTPUT PARAMETERS:
!
      type(GeneralGrid), intent(out)   :: GGrid

! !REVISION HISTORY:
!       25Sep00 - Jay Larson <larson@mcs.anl.gov> - initial prototype
!       15Jan01 - Jay Larson <larson@mcs.anl.gov> - modified to fit
!                 new GeneralGrid definition.  
!       19Mar01 - Jay Larson <larson@mcs.anl.gov> - added OtherChars
!       25Apr01 - Jay Larson <larson@mcs.anl.gov> - added GlobGridNum
!                 as a mandatory integer attribute.
!       13Jun01 - Jay Larson <larson@mcs.anl.gov> - No longer define 
!                 blank List attributes of the GeneralGrid.  Previous
!                 versions of this routine had this feature, and this
!                 caused problems with the GeneralGrid Send and Receive
!                 operations on the AIX platform.
!       13Jun01 - R. Jacob <jacob@mcs.anl.gov> - nullify any pointers
!                 for lists not declared.
!       15Feb02 - Jay Larson <larson@mcs.anl.gov> - made the input 
!                 argument CoordSortOrder mandatory (rather than
!                 optional).
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//'::init_'

  character*128 :: RAList, IAList
  integer :: n, list_len, ier, RA_len, IA_len
  integer :: ierr, ncoord, nsort

        ! if lsize is present, use it to set n; if not, set n=0

  n = 0
  if(present(lsize)) n=lsize

        ! Concatenate onto CoordList the values of WeightChars and
        ! OtherChars (if present)--store the result in RAList

  list_len = len(CoordChars)

  if(present(WeightChars)) then
     list_len = list_len + 1 + len(WeightChars)
  endif

  if(present(OtherChars)) then
     list_len = list_len + 1 + len(OtherChars)
  endif

  RAList = trim(CoordChars)
  RA_len = len(CoordChars)

  if(present(WeightChars)) then
     RAList = RAList(1:RA_len) // ':' // trim(WeightChars)
     RA_len = RA_len + 1 + len(WeightChars)
  endif

  if(present(OtherChars)) then
     RAList = RAList(1:RA_len) // ':' // trim(OtherChars)
     RA_len = RA_len + 1 + len(OtherChars)
  endif

        ! Create the integer Attribute list IAList.  At a minimum,
        ! IAList contains the attribute GlobGridNum

  if(present(IndexChars)) then
     IAList = trim(IndexChars) // ':GlobGridNum'
  else
     IAList = 'GlobGridNum'
  endif

  call List_init(GGrid%index_list, IAList)

        ! Initialize the AttrVect storage component GGrid%data

  call AttrVect_init(aV=GGrid%data, iList=trim(IAList), &
	                rList=trim(RAList), lsize=n)

        ! Initialize GGrid%coordinate_list

  call List_init(GGrid%coordinate_list, CoordChars)

        ! If Initialize GGrid%weight_list

  if(present(WeightChars)) then
     call List_init(GGrid%weight_list, WeightChars)
  else
     nullify(GGrid%weight_list%bf)
  endif

        ! If Initialize GGrid%other_list

  if(present(OtherChars)) then
     call List_init(GGrid%other_list, OtherChars)
  else
     nullify(GGrid%other_list%bf)
  endif

        ! Initialize GGrid%coordinate_sort_order.

  call List_init(GGrid%coordinate_sort_order, CoordSortOrder)

        ! Check the number of coordinates versus the number of
        ! coordinate grid sort keys...they should be equal.

     ncoord = List_nitem(GGrid%coordinate_list)
     nsort =  List_nitem(GGrid%coordinate_sort_order)

  if(ncoord /= nsort) then
     write(stderr,*) myname_,':: ERROR Arguments ncoord and nsort must be equal.',&
	  ' ncoord = ',ncoord,' nsort = ',nsort
     call die(myname_,'ncoord-nsort /=0',ncoord-nsort)
  endif

        ! If the LOGICAL argument descend is present, check the
        ! number of entries to ensure they match the grid dimensionality.
        ! If descend is not present, assume all coordinate grid point
        ! sortings will be in ascending order.

  if(present(descend)) then
     if(size(descend) /= nsort) then
     write(stderr,*) myname_, &
     ':: ERROR Number of elements in descend(:) must equal nsort.',&
	  ' size(descend) = ',size(descend),' nsort = ',nsort
	call die(myname_, 'size(descend)-nsort /=0', size(descend)-nsort)
     endif
  endif

  allocate(GGrid%descend(nsort), stat=ierr)
  if(ierr /= 0) then
     call die(myname_,"allocate(GGrid%descend) failed.",ierr)
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
! !IROUTINE: initl_ - initialize an empty GeneralGrid from Lists.
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

 subroutine initl_(GGrid, CoordList, CoordSortOrder, descend, WeightList, &
                   OtherList, IndexList, lsize )
!
! !USES:
!

      use m_stdio
      use m_die

      use m_List,     only : List
      use m_List,     only : List_init => init
      use m_List,     only : List_nitem => nitem
      use m_List,     only : List_concatenate => concatenate
      use m_List,     only : List_copy => copy
      use m_List,     only : List_clean => clean

      use m_AttrVect, only : AttrVect
      use m_AttrVect, only : AttrVect_init => init

      implicit none

! !INPUT PARAMETERS:
!
      Type(List),                      intent(in)  :: CoordList
      Type(List),                      intent(in)  :: CoordSortOrder
      Type(List),            optional, intent(in)  :: WeightList
      logical, dimension(:), optional, pointer     :: descend
      Type(List),            optional, intent(in)  :: OtherList
      Type(List),            optional, intent(in)  :: IndexList
      integer,               optional, intent(in)  :: lsize

! !OUTPUT PARAMETERS:
!
      type(GeneralGrid),               intent(out) :: GGrid

! !REVISION HISTORY:
!       10May01 - Jay Larson <larson@mcs.anl.gov> - initial version
!       08Aug01 - E.T. Ong <eong@mcs.anl.gov> - changed list assignment(=)
!                 to list copy to avoid compiler bugs with pgf90
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//'::initl_'

  type(List) :: TempList1, Templist2
  type(List) :: RAList, IAList
  integer :: i, n, ierr

       ! Copy array descend(:) (if present) to GGrid%descend.

  n = List_nitem(CoordList)

  allocate(GGrid%descend(n), stat=ierr)
  if(ierr /= 0) then
     call die(myname_,"allocate GGrid%descend...",ierr)
  endif

  if(present(descend)) then
     GGrid%descend(1:n) = descend(1:n)
  else
     GGrid%descend(1:n) = .FALSE.
  endif

       ! Process input lists and create the appropriate GeneralGrid
       ! List components

  call List_copy(GGrid%coordinate_list,CoordList)
  call List_copy(GGrid%coordinate_sort_order,CoordSortOrder)

  call List_copy(TempList1,CoordList)

       ! Concatenate present input Lists to create RAList, and
       ! at the same time assign the List components of GGrid

  if(present(WeightList)) then
     call List_copy(GGrid%weight_list,WeightList)
     call List_concatenate(TempList1, WeightList, TempList2)
     call List_clean(TempList1)
  else
     call List_init(GGrid%weight_list,'')
     call List_copy(TempList2,TempList1)
     call List_clean(TempList1)
  endif

  if(present(OtherList)) then
     call List_copy(GGrid%other_list,OtherList)
     call List_concatenate(TempList2, OtherList, RAList)
     call List_clean(TempList2)
  else
     call List_init(GGrid%other_list,'')
     call List_copy(RAList,TempList2)
     call List_clean(TempList2)
  endif

       ! Concatenate present input Lists to create IAList

  call List_init(TempList1,'GlobGridNum')

  if(present(IndexList)) then
     call List_concatenate(TempList1,IndexList,IAList)
     call List_clean(TempList1)
  else
     call List_copy(IAList,TempList1)
     call List_clean(TempList1)
  endif

  call List_copy(GGrid%index_list,IAList)


       ! Initialize GGrid%data using IAList, RAList, and lsize (if
       ! present).

  n = 0
  if(present(lsize)) n = lsize

  call AttrVect_init(GGrid%data, IAList, RAList, n)

       ! Finally, initialize GGrid%descend--first, count the coordinates:

  n = List_nitem(GGrid%coordinate_list)

       ! Allocate GGrid%descend

  allocate(GGrid%descend(n), stat=ierr)

       ! Initialize GGrid%descend from descend(:), if present.  If
       ! the argument descend(:) was not passed, set all the entries
       ! of GGrid%descend(:) to FALSE.

  if(present(descend)) then
     do i=1,n
	GGrid%descend(i) = descend(i)
     end do
  else
     GGrid%descend = .FALSE.
  endif

 end subroutine initl_

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
      use m_stdio
      use m_die

      use m_List, only : List
      use m_List, only : List_copy => copy

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
!       13Jun01 - Jay Larson <larson@mcs.anl.gov> - Now, undefined List
!                 components of the GeneralGrid iGGrid are no longer 
!                 copied to oGGrid.
!       08Aug01 - E.T. Ong <eong@mcs.anl.gov> - changed list assignment(=)
!                 to list copy to avoid compiler bugs with pgf90
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
	write(stderr,*) myname_,':: ERROR size(iGGrid%descend) must equal ncoord,',&
	     ' size(iGGrid%descend)=',size(iGGrid%descend),' ncoord=',ncoord
	call die(myname,"iGGrid dims/size of descend mismatch", &
	                 ncoord-size(iGGrid%descend))
     endif
  endif

        ! If iGGrid%descend has been allocated, copy its contents

  if(associated(iGGrid%descend)) then ! allocate and fill oGGrid%descend

     allocate(oGGrid%descend(ncoord), stat=ierr)
     if(ierr /= 0) then
	call die(myname,"allocate(oGGrid%descend...", ierr)
     endif

     do i=1,ncoord
	oGGrid%descend(i) = iGGrid%descend(i)
     end do

  endif

       ! Copy list data from iGGrid to oGGrid.  We know from
       ! interface and GeneralGrid type definition that 
       ! iGGrid%coordinate_list, iGGrid%coordinate_sort_order, 
       ! and iGGrid%index_list _should_ be defined.  Test all the
       ! list attributes before copying them.

  if(associated(oGGrid%coordinate_list%bf)) then
     call List_copy(oGGrid%coordinate_list,iGGrid%coordinate_list)
  endif

  if(associated(oGGrid%coordinate_sort_order%bf)) then
     call List_copy(oGGrid%coordinate_sort_order,iGGrid%coordinate_sort_order)
  endif

  if(associated(oGGrid%weight_list%bf)) then
     call List_copy(oGGrid%weight_list,iGGrid%weight_list)
  endif

  if(associated(oGGrid%other_list%bf)) then
     call List_copy(oGGrid%other_list,iGGrid%other_list)
  endif

  call List_copy(oGGrid%index_list,iGGrid%index_list)

       ! Now, initialize oGGrid%data from iGGrid%data, but 
       ! with length n.

  call AttrVect_init(oGGrid%data, iGGrid%data, n)

 end subroutine initgg_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: initCartesian_ - initialize a Cartesian GeneralGrid.
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
! Once the storage space in {\tt GGrid} is initialized, The gridpoint 
! coordinates are evaluated using the input arguments {\tt Dims} (the 
! number of points on each coordinate axis) and {\tt AxisData} (the 
! coordinate values on all of the points of all of the axes).
!
! !INTERFACE:

 subroutine initCartesian_(GGrid, CoordChars, CoordSortOrder, descend, &
                           WeightChars, OtherChars, IndexChars, Dims, AxisData)
!
! !USES:
!
      use m_stdio
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
      character(len=*),             intent(in) :: CoordChars
      character(len=*), optional,   intent(in) :: CoordSortOrder
      character(len=*), optional,   intent(in) :: WeightChars
      logical, dimension(:), optional, pointer :: descend
      character(len=*), optional,   intent(in) :: OtherChars
      character(len=*), optional,   intent(in) :: IndexChars
      integer, dimension(:),        pointer    :: Dims
      real, dimension(:),           pointer    :: AxisData

! !OUTPUT PARAMETERS:
!
      type(GeneralGrid), intent(out)   :: GGrid

! !REVISION HISTORY:
!       07Jun01 - Jay Larson <larson@mcs.anl.gov> - initial version.
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//'::initCartesian_'

 end subroutine initCartesian_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: initUnstructured_ - initialize an Unstructured GeneralGrid.
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
! Once the storage space in {\tt GGrid} is initialized, The gridpoint 
! coordinates are evaluated using the input arguments {\tt Dims} (the 
! number of points on each coordinate axis) and {\tt AxisData} (the 
! coordinate values on all of the points of all of the axes).
!
! !INTERFACE:

 subroutine initUnstructured_(GGrid, CoordChars, CoordSortOrder, descend, &
                           WeightChars, OtherChars, IndexChars, nPoints, &
                           PointData)
!
! !USES:
!
      use m_stdio
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
      character(len=*),             intent(in) :: CoordChars
      character(len=*), optional,   intent(in) :: CoordSortOrder
      character(len=*), optional,   intent(in) :: WeightChars
      logical, dimension(:), optional, pointer :: descend
      character(len=*), optional,   intent(in) :: OtherChars
      character(len=*), optional,   intent(in) :: IndexChars
      integer,                      intent(in) :: nPoints
      real, dimension(:),           pointer    :: PointData

! !OUTPUT PARAMETERS:
!
      type(GeneralGrid), intent(out)   :: GGrid

! !REVISION HISTORY:
!       07Jun01 - Jay Larson <larson@mcs.anl.gov> - initial version.
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//'::initUnstructured_'

 end subroutine initUnstructured_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - Destroy a GeneralGrid.
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean_(GGrid,stat)
!
! !USES:
!
      use m_stdio
      use m_die

      use m_List,only : List_clean => clean
      use m_AttrVect,only : AttrVect_clean => clean

      implicit none

! !INPUT/OUTPUT PARAMETERS: 
!
      type(GeneralGrid), intent(inout) :: GGrid
      integer, optional, intent(out)   :: stat

! !REVISION HISTORY:
!       25Sep00 - J.W. Larson <larson@mcs.anl.gov> - initial prototype
!       20Mar01 - J.W. Larson <larson@mcs.anl.gov> - complete version.
!       01Mar01 - E.T. Ong <eong@mcs.anl.gov> - removed dies to prevent
!                 crashes when cleaning uninitialized attrvects. Added
!                 optional stat argument.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'
  integer :: ierr

  if(present(stat)) then

     stat=0
     call AttrVect_clean(GGrid%data,ierr)
     if(ierr/=0) stat=ierr

     call List_clean(GGrid%coordinate_list,ierr)
     if(ierr/=0) stat=ierr
     call List_clean(GGrid%coordinate_sort_order,ierr)
     if(ierr/=0) stat=ierr
     call List_clean(GGrid%weight_list,ierr)
     if(ierr/=0) stat=ierr
     call List_clean(GGrid%other_list,ierr)
     if(ierr/=0) stat=ierr
     call List_clean(GGrid%index_list,ierr)
     if(ierr/=0) stat=ierr

  else

     call AttrVect_clean(GGrid%data)

     call List_clean(GGrid%coordinate_list)
     call List_clean(GGrid%coordinate_sort_order)
     call List_clean(GGrid%weight_list)
     call List_clean(GGrid%other_list)
     call List_clean(GGrid%index_list)

  endif

  deallocate(GGrid%descend, stat=ierr)

  if(ierr /= 0) then
     if(present(stat)) then
	stat=ierr
     else
	call warn(myname_,'deallocate(GGrid%descend...',ierr)
     endif
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
      use m_stdio
      use m_die

      use m_List,     only : List_nitem => nitem

      implicit none

! !INPUT PARAMETERS: 
!
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
      use m_die
      use m_stdio

      use m_AttrVect,     only : AttrVect_indexIA => indexIA

      implicit none

! !INPUT PARAMETERS: 
!
      type(GeneralGrid), intent(in)  :: GGrid
      character(len=*),  intent(in)  :: item
      character(len=*), optional, intent(in) :: perrWith
      character(len=*), optional, intent(in) :: dieWith

! !REVISION HISTORY:
!       15Jan01 - Jay Larson <larson@mcs.anl.gov> - Initial version.
!       27Mar02 - Jay Larson <larson@mcs.anl.gov> - Cleaned up error
!                 handling logic.
!EOP ___________________________________________________________________
!
 character(len=*),parameter :: myname_=myname//'::indexIA_'

 
 indexIA_ = AttrVect_indexIA(GGrid%data, item)

  if(indexIA_<=0) then

     if((.not. present(perrWith)) .and. (.not. present(dieWith))) then
	write(stderr,'(4a)') myname,       &
		'" :: indexIA_() error, not found "',trim(item),'"'
	call die(myname_)
     else
	if(present(perrWith)) then
	   write(stderr,'(4a)') perrWith, &
		'" indexIA_() error, not found "',trim(item),'"'
	endif
	if(present(dieWith)) then
	   write(stderr,'(4a)') dieWith,       &
		'" indexIA_() error, not found "',trim(item),'"'
	   call die(dieWith)
	endif
     endif ! if((.not. present(perrWith)) .and. ...

  endif ! if(indexIA_ <= 0)

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
      use m_stdio
      use m_die

      use m_AttrVect,     only : AttrVect_indexRA => indexRA

      implicit none

! !INPUT PARAMETERS: 
!
      type(GeneralGrid),          intent(in)  :: GGrid
      character(len=*),           intent(in)  :: item
      character(len=*), optional, intent(in) :: perrWith
      character(len=*), optional, intent(in) :: dieWith

! !REVISION HISTORY:
!       15Jan01 - Jay Larson <larson@mcs.anl.gov> - Initial version.
!       27Mar02 - Jay Larson <larson@mcs.anl.gov> - Cleaned up error
!                 handling logic.
!EOP ___________________________________________________________________
!
 character(len=*),parameter :: myname_=myname//'::indexRA_'

 
 indexRA_ = AttrVect_indexRA(GGrid%data, item)

  if(indexRA_<=0) then

     if((.not. present(perrWith)) .and. (.not. present(dieWith))) then
	write(stderr,'(4a)') myname,       &
		'" :: indexRA_() error, not found "',trim(item),'"'
	call die(myname_)
     else
	if(present(perrWith)) then
	   write(stderr,'(4a)') perrWith, &
		'" indexRA_() error, not found "',trim(item),'"'
	endif
	if(present(dieWith)) then
	   write(stderr,'(4a)') dieWith,       &
		'" indexRA_() error, not found "',trim(item),'"'
	   call die(dieWith)
	endif
     endif ! if((.not. present(perrWith)) .and. ...

  endif ! if(indexRA_ <= 0)

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

! !INPUT PARAMETERS: 
!
      type(GeneralGrid), intent(in)  :: GGrid

! !REVISION HISTORY:
!       15Jan01 - Jay Larson <larson@mcs.anl.gov> - Initial version.
!       27Mar02 - Jay Larson <larson@mcs.anl.gov> - slight logic change.
!EOP ___________________________________________________________________
!
 character(len=*),parameter :: myname_=myname//'::lsize_'

 integer :: lsizeI, lsizeR

 lsizeI = 0
 lsizeR = 0

 lsize_ = 0

 if(associated(GGrid%data%rAttr) .or. associated(GGrid%data%iAttr)) then

    if(associated(GGrid%data%rAttr)) then
       lsizeR = AttrVect_lsize( GGrid%data )
    endif
    if(associated(GGrid%data%rAttr)) then
       lsizeR = AttrVect_lsize( GGrid%data )
    endif

    lsize_ = max(lsizeI,lsizeR)

 endif

 end function lsize_

!BOP -------------------------------------------------------------------
!
! !IROUTINE: exportIAttr_ - return GeneralGrid INTEGER attribute as a vector
!
! !DESCRIPTION:
! This routine extracts from the input {\tt GeneralGrid} argument 
! {\tt GGrid} the integer attribute corresponding to the tag defined in 
! the input {\tt CHARACTER} argument {\tt AttrTag}, and returns it in 
! the {\tt INTEGER} output array {\tt outVect}, and its length in the 
! output {\tt INTEGER} argument {\tt lsize}.
!
! {\bf N.B.:}  This routine will fail if the {\tt AttrTag} is not in 
! the {\tt GeneralGrid} {\tt List} component {\tt GGrid\%data\%iList}.
!
! {\bf N.B.:}  The flexibility of this routine regarding the pointer 
! association status of the output argument {\tt outVect} means the
! user must invoke this routine with care.  If the user wishes this
! routine to fill a pre-allocated array, then obviously this array
! must be allocated prior to calling this routine.  If the user wishes
! that the routine {\em create} the output argument array {\tt outVect},
! then the user must ensure this pointer is not allocated (i.e. the user
! must nullify this pointer) at the time this routine is invoked.
!
! {\bf N.B.:}  If the user has relied on this routine to allocate memory
! associated with the pointer {\tt outVect}, then the user is responsible 
! for deallocating this array once it is no longer needed.  Failure to 
! do so will result in a memory leak.
!
! !INTERFACE:

 subroutine exportIAttr_(GGrid, AttrTag, outVect, lsize)
!
! !USES:
!
      use m_die 
      use m_stdio

      use m_AttrVect,      only : AttrVect_exportIAttr => exportIAttr

      implicit none

! !INPUT PARAMETERS: 

      type(GeneralGrid),      intent(in)  :: GGrid
      character(len=*),       intent(in)  :: AttrTag

! !OUTPUT PARAMETERS: 

      integer,  dimension(:), pointer     :: outVect
      integer,                intent(out) :: lsize

! !REVISION HISTORY:
!       13Dec01 - J.W. Larson <larson@mcs.anl.gov> - initial prototype.
!
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::exportIAttr_'

       ! Export the data (inheritance from AttrVect)

  call AttrVect_exportIAttr(GGrid%data, AttrTag, outVect, lsize)

 end subroutine exportIAttr_

!BOP -------------------------------------------------------------------
!
! !IROUTINE: exportRAttr_ - return GeneralGrid REAL attribute as a vector
!
! !DESCRIPTION:
! This routine extracts from the input {\tt GeneralGrid} argument 
! {\tt GGrid} the real attribute corresponding to the tag defined in 
! the input {\tt CHARACTER} argument {\tt AttrTag}, and returns it in 
! the {\tt REAL} output array {\tt outVect}, and its length in the 
! output {\tt INTEGER} argument {\tt lsize}.
!
! {\bf N.B.:}  This routine will fail if the {\tt AttrTag} is not in 
! the {\tt GeneralGrid} {\tt List} component {\tt GGrid\%data\%rList}.
!
! {\bf N.B.:}  The flexibility of this routine regarding the pointer 
! association status of the output argument {\tt outVect} means the
! user must invoke this routine with care.  If the user wishes this
! routine to fill a pre-allocated array, then obviously this array
! must be allocated prior to calling this routine.  If the user wishes
! that the routine {\em create} the output argument array {\tt outVect},
! then the user must ensure this pointer is not allocated (i.e. the user
! must nullify this pointer) at the time this routine is invoked.
!
! {\bf N.B.:}  If the user has relied on this routine to allocate memory
! associated with the pointer {\tt outVect}, then the user is responsible 
! for deallocating this array once it is no longer needed.  Failure to 
! do so will result in a memory leak.
!
! !INTERFACE:

 subroutine exportRAttr_(GGrid, AttrTag, outVect, lsize)
!
! !USES:
!
      use m_die
      use m_stdio

      use m_AttrVect,      only : AttrVect_exportRAttr => exportRAttr

      implicit none

! !INPUT PARAMETERS: 

      type(GeneralGrid),      intent(in)  :: GGrid
      character(len=*),       intent(in)  :: AttrTag

! !OUTPUT PARAMETERS: 

      real,  dimension(:),    pointer     :: outVect
      integer,                intent(out) :: lsize

! !REVISION HISTORY:
!       13Dec01 - J.W. Larson <larson@mcs.anl.gov> - initial prototype.
!
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::exportRAttr_'

       ! Export the data (inheritance from AttrVect)

  call AttrVect_exportRAttr(GGrid%data, AttrTag, outVect, lsize)

 end subroutine exportRAttr_

!BOP -------------------------------------------------------------------
!
! !IROUTINE: importIAttr_ - import GeneralGrid INTEGER attribute
!
! !DESCRIPTION:
! This routine imports data provided in the input {\tt INTEGER} vector 
! {\tt inVect} into the {\tt GeneralGrid} argument {\tt GGrid}, storing 
! it as the integer attribute corresponding to the tag defined in 
! the input {\tt CHARACTER} argument {\tt AttrTag}.  The input 
! {\tt INTEGER} argument {\tt lsize} is used to ensure there is 
! sufficient space in the {\tt GeneralGrid} to store the data.
!
! {\bf N.B.:}  This routine will fail if the {\tt AttrTag} is not in 
! the {\tt GeneralGrid} {\tt List} component {\tt GGrid\%data\%iList}.
!
! !INTERFACE:

 subroutine importIAttr_(GGrid, AttrTag, inVect, lsize)
!
! !USES:
!
      use m_die
      use m_stdio

      use m_AttrVect,      only : AttrVect_importIAttr => importIAttr

      implicit none

! !INPUT PARAMETERS: 

      character(len=*),       intent(in)    :: AttrTag
      integer,  dimension(:), pointer       :: inVect
      integer,                intent(in)    :: lsize

! !INPUT/OUTPUT PARAMETERS: 

      type(GeneralGrid),      intent(inout) :: GGrid

! !REVISION HISTORY:
!       13Dec01 - J.W. Larson <larson@mcs.anl.gov> - initial prototype.
!       27Mar02 - Jay Larson <larson@mcs.anl.gov> - improved error handling.
!
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::importIAttr_'

       ! Argument Check:

  if(lsize > lsize_(GGrid)) then
     write(stderr,*) myname_,':: ERROR, lsize > lsize_(GGrid).', &
	  'lsize = ',lsize,'lsize_(GGrid) = ',lsize_(GGrid)
     call die(myname_)
  endif

       ! Import the data (inheritance from AttrVect)

  call AttrVect_importIAttr(GGrid%data, AttrTag, inVect, lsize)

 end subroutine importIAttr_

!BOP -------------------------------------------------------------------
!
! !IROUTINE: importRAttr_ - import GeneralGrid REAL attribute
!
! !DESCRIPTION:
! This routine imports data provided in the input {\tt REAL} vector 
! {\tt inVect} into the {\tt GeneralGrid} argument {\tt GGrid}, storing 
! it as the real attribute corresponding to the tag defined in 
! the input {\tt CHARACTER} argument {\tt AttrTag}.  The input 
! {\tt INTEGER} argument {\tt lsize} is used to ensure there is 
! sufficient space in the {\tt GeneralGrid} to store the data.
!
! {\bf N.B.:}  This routine will fail if the {\tt AttrTag} is not in 
! the {\tt GeneralGrid} {\tt List} component {\tt GGrid\%data\%rList}.
!
! !INTERFACE:

 subroutine importRAttr_(GGrid, AttrTag, inVect, lsize)
!
! !USES:
!
      use m_die ,          only : die
      use m_die ,          only : MP_perr_die
      use m_stdio ,        only : stderr

      use m_AttrVect,      only : AttrVect_importRAttr => importRAttr

      implicit none

! !INPUT PARAMETERS: 

      character(len=*),       intent(in)    :: AttrTag
      real, dimension(:),     pointer       :: inVect
      integer,                intent(in)    :: lsize

! !INPUT/OUTPUT PARAMETERS: 

      type(GeneralGrid),      intent(inout) :: GGrid

! !REVISION HISTORY:
!       13Dec01 - J.W. Larson <larson@mcs.anl.gov> - initial prototype.
!       27Mar02 - Jay Larson <larson@mcs.anl.gov> - improved error handling.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::importRAttr_'

       ! Argument Check:

  if(lsize > lsize_(GGrid)) then
     write(stderr,*) myname_,':: ERROR, lsize > lsize_(GGrid).', &
	  'lsize = ',lsize,'lsize_(GGrid) = ',lsize_(GGrid)
     call die(myname_)
  endif

       ! Import the data (inheritance from AttrVect)

  call AttrVect_importRAttr(GGrid%data, AttrTag, inVect, lsize)

 end subroutine importRAttr_

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
      use m_stdio
      use m_die

      use m_AttrVect,     only : AttrVect_Sort => Sort
      use m_List,        only : List_nitem => nitem

      implicit none

! !INPUT PARAMETERS: 
!
      type(GeneralGrid),               intent(in) :: GGrid
      type(List),                      intent(in) :: key_list
      logical, dimension(:), optional, intent(in) :: descend

! !OUTPUT PARAMETERS: 
!
      integer, dimension(:), pointer              :: perm


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
     call die(myname_,"allocate(descending...",ierr)
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
     call die(myname_,"deallocate(descending...",ierr)
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
! the input  {\tt GeneralGrid} variable {\tt GGrid\%coordinate\_sort\_order}
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

! !INPUT PARAMETERS: 
!
      type(GeneralGrid),     intent(in) :: GGrid

! !OUTPUT PARAMETERS: 
!
      integer, dimension(:), pointer    :: perm

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

      use m_stdio
      use m_die

      use m_AttrVect,     only : AttrVect
      use m_AttrVect, only : AttrVect_Permute => Permute

      implicit none

! !INPUT PARAMETERS: 
!
      integer, dimension(:), intent(in)    :: perm

! !INPUT/OUTPUT PARAMETERS: 
!
      type(GeneralGrid),     intent(inout) :: GGrid


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
      use m_stdio
      use m_die

      implicit none

! !INPUT/OUTPUT PARAMETERS: 
!
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
  if(ierr /= 0) then
     call die(myname_,"deallocate(perm)",ierr)
  endif

 end subroutine SortPermute_

 end module m_GeneralGrid
















































