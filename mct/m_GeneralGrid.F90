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
! {\tt GGrid\%data\%iList}) is defined by {\tt GGrid\%index\_list}.
!
! !INTERFACE:

 module m_GeneralGrid

!
! !USES:
!
      use m_List, only : List   ! Support for List components.

      use m_AttrVect, only : AttrVect
      use m_AttrVect, only : AttrVect_init => init
      use m_AttrVect, only : AttrVect_clean => clean

      use m_AttrVect, only : AttrVect_Sort => Sort
      use m_AttrVect, only : AttrVect_Permute => Permute
      use m_AttrVect, only : AttrVect_SortPermute => SortPermute

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
      type(List) :: coordinate_list
      type(List) :: weight_list
      type(List) :: index_list
      type(AttrVect) :: data
    end type GeneralGrid

    interface init  ; module procedure &
	      init_, &
	      initv_
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
!                 definition and added numerous API's.
!       17Jan01 - J.W. Larson fixed minor bug in module header use
!                 statement.
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
!
! !INTERFACE:

 subroutine init_(GGrid, CoordList, WeightList, IndexList, lsize )
!
! !USES:
!
      use m_String,   only : String, char

      use m_AttrVect, only : AttrVect_init => init
      use m_die

      implicit none

      type(GeneralGrid), intent(out)   :: GGrid
      character(len=*),  intent(in)    :: CoordList
      character(len=*), optional, intent(in) :: WeightList
      character(len=*), optional, intent(in) :: IndexList
      integer,    optional,intent(in)  :: lsize

! !REVISION HISTORY:
!       25Sep00 - Jay Larson <larson@mcs.anl.gov> - initial prototype
!       15Jan01 - Jay Larson <larson@mcs.anl.gov> - modified to fit
!                 new GeneralGrid definition.  
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//'::init_'

  character, dimension(:), pointer :: RAList

  integer :: n, list_len, ier

        ! if lsize is present, use it to set n; if not, set n=0

  n = 0
  if(present(lsize)) n=lsize

        ! Initialize sMat using AttrVect_init

  if(present(WeightList)) then
     list_len = len(CoordList) + len(WeightList) + 1
     allocate(RAList(list_len), stat=ier)
     RAList = CoordList // ':' // WeightList
  endif

  if(present(IndexList)) then
     if(present(WeightList)) then
!        call AttrVect_init(GGrid%data, iList=IndexList, rList=RAList, &
!	  lsize=n)
     else
        call AttrVect_init(GGrid%data, iList=IndexList, rList=CoordList, &
	  lsize=n)
     endif
  else
     if(present(WeightList)) then
!        call AttrVect_init(GGrid%data, iList=' ', rList=RAList, &
!	  lsize=n)
     else
        call AttrVect_init(GGrid%data, iList=' ', rList=CoordList, &
	  lsize=n)
     endif
  endif

  if(present(WeightList)) deallocate(RAList, stat=ier)

 end subroutine init_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: initv_ - initialize the GeneralGrid using coordinate data
!
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine initv_(GGrid, CoordList, WeightList, IndexList, CoordData, &
                   WeightData, IndexData, lsize )
!
! !USES:
!
      use m_AttrVect, only : AttrVect_init => init
      use m_die

      implicit none

      type(GeneralGrid), intent(out)   :: GGrid
      character(len=*),  intent(in)    :: CoordList
      character(len=*), optional, intent(in) :: WeightList
      character(len=*), optional, intent(in) :: IndexList
      integer,    optional,intent(in)  :: lsize
      real, dimension(:,:),          intent(in) :: CoordData
      real, dimension(:,:), optional,intent(in) :: WeightData
      integer, dimension(:,:), optional,intent(in) :: IndexData

! !REVISION HISTORY:
!       15Jan01 - Jay Larson <larson@mcs.anl.gov> - API specification.
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//'::initv_'

 end subroutine initv_

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
      use m_AttrVect,only : AttrVect_clean => clean

      implicit none

      type(GeneralGrid), intent(inout) :: GGrid

! !REVISION HISTORY:
!       25Sep00 - J.W. Larson <larson@mcs.anl.gov> - initial prototype
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'

  call AttrVect_clean(GGrid%data)

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
! !IROUTINE: Sort_ - returns sort permutation defined by coordinate keys.
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

 subroutine Sort_(GGrid, key_List, perm, descend, perrWith, dieWith)

!
! !USES:
!

      use m_AttrVect,     only : AttrVect_Sort => Sort

      implicit none

      type(GeneralGrid), intent(in)  :: GGrid
      type(List),        intent(in)  :: key_list
      integer, dimension(:), pointer              :: perm
      logical, dimension(:), optional, intent(in) :: descend
      character(len=*), optional, intent(in)      :: perrWith
      character(len=*), optional, intent(in)      :: dieWith

! !REVISION HISTORY:
!       15Jan01 - Jay Larson <larson@mcs.anl.gov> - Initial version.
!EOP ___________________________________________________________________
!
 character(len=*),parameter :: myname_=myname//'::Sort_'

 end subroutine Sort_

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

 subroutine Permute_(GGrid, perm, perrWith, dieWith)

!
! !USES:
!

      use m_die ,          only : die
      use m_stdio ,        only : stderr
      use m_AttrVect, only : AttrVect_Permute => Permute

      implicit none

      type(GeneralGrid),     intent(inout)   :: GGrid
      integer, dimension(:), intent(in)      :: perm
      character(len=*), optional, intent(in) :: perrWith
      character(len=*), optional, intent(in) :: dieWith

! !REVISION HISTORY:
!       15Jan01 - Jay Larson <larson@mcs.anl.gov> - API specification.
!EOP ___________________________________________________________________
!
 character(len=*),parameter :: myname_=myname//'::Permute_'

 end subroutine Permute_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: SortPermute_ - sort and re-order coordinate data
!
! !DESCRIPTION:
! The subroutine {\tt SortPermute\_()} uses the list of keys defined in 
! {\tt key\_list} to creat an index permutation {\tt perm}, which is 
! then applied to re-order the coordinate data stored in the {\tt GeneralGrid} 
! argument {\tt GGrid}.  This permutation is generated  by the routine 
! {\tt Sort\_()} contained in this module.  The permutation is carried 
! out by the routine {\tt Permute\_()} contained in this module.
!
! !INTERFACE:

 subroutine SortPermute_(GGrid, key_list, perrWith, dieWith)

!
! !USES:
!

      use m_die ,          only : die
      use m_stdio ,        only : stderr
      use m_AttrVect, only : AttrVect_Permute => Permute

      implicit none

      type(GeneralGrid),     intent(inout)   :: GGrid
      type(List),            intent(in)      :: key_list
      character(len=*), optional, intent(in) :: perrWith
      character(len=*), optional, intent(in) :: dieWith

! !REVISION HISTORY:
!       15Jan01 - Jay Larson <larson@mcs.anl.gov> - API specification.
!EOP ___________________________________________________________________
!
 character(len=*),parameter :: myname_=myname//'::SortPermute_'

 end subroutine SortPermute_

 end module m_GeneralGrid

