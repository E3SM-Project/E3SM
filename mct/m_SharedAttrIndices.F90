!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_SharedAttrIndices - Mutual indexing tools
!
! !DESCRIPTION:
! This module contains tools for determining lists of shared attributes 
! between pairs of {\tt AttrVect} objects, pairs of {\tt Accumulator}
! objects, and combinations of the two types.
!
!
! !INTERFACE:

 module m_SharedAttrIndices

      implicit none

      private   ! except

      public :: SharedAttrIndexList  ! Returns the number of shared 
                                     ! attributes, and lists of the
                                     ! respective locations of these
                                     ! shared attributes
                                     
    interface SharedAttrIndexList ; module procedure   &
        aVaVSharedAttrIndexList_,  &   
        aCaCSharedAttrIndexList_,  &   
        aVaCSharedAttrIndexList_
    end interface

! !REVISION HISTORY:
!       07Feb01 - J.W. Larson <larson@mcs.anl.gov> - initial version
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='SharedAttrIndices'

 contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: aVaVSharedAttrIndexList_ - AttrVect shared attributes.
!
! !DESCRIPTION:  {\tt aVaVSharedAttrIndexList\_()} takes a pair of 
! user-supplied {\tt AttrVect} variables {\tt aV1} and {\tt aV2}, 
! and for choice of either {\tt REAL} or {\tt INTEGER} attributes (as
! specified literally in the input {\tt CHARACTER} argument {\tt attrib})
! returns the number of shared attributes {\tt NumShared}, and arrays of
! indices {\tt Indices1} and {\tt Indices2} to their storage locations
! in {\tt aV1} and {\tt aV2}, respectively.
!
! {\bf N.B.:}  This routine returns two allocated arrays---{\tt Indices1(:)} 
! and {\tt Indices2(:)}---which must be deallocated once the user no longer
! needs them.  Failure to do this will create a memory leak.
!
! !INTERFACE:

 subroutine aVaVSharedAttrIndexList_(aV1, aV2, attrib, NumShared, &
                                     Indices1, Indices2)

!
! !USES:
!
      use m_stdio
      use m_die,      only : MP_perr_die

      use m_AttrVect, only : AttrVect
      use m_AttrVect, only : AttrVect_lsize => lsize

      implicit none

! !INPUT PARAMETERS: 
!
      type(AttrVect),    intent(in)  :: aV1   
      type(AttrVect),    intent(in)  :: aV2
      character*7,       intent(in)  :: attrib

! !OUTPUT PARAMETERS:   
!
      integer,           intent(out) :: NumShared

      integer,dimension(:), pointer  :: Indices1
      integer,dimension(:), pointer  :: Indices2

! !REVISION HISTORY:
!       07Feb01 - J.W. Larson <larson@mcs.anl.gov> - initial version
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::aVaVSharedAttrIndexList_'

  integer :: ierr

       ! Based on the value of the argument attrib, pass the 
       ! appropriate pair of Lists for comparison...

  select case(trim(attrib))
  case('REAL','real')
     call GetSharedListIndices(aV1%rList, aV2%rList, NumShared, &
	                         Indices1, Indices2)
  case('INTEGER','integer')
     call GetSharedListIndices(aV1%iList, aV2%iList, NumShared, &
	                         Indices1, Indices2)
  case default
     write(stderr,'(4a)') myname_,":: value of argument attrib=",attrib, &
	  " not recognized.  Allowed values: REAL, real, INTEGER, integer"
     ierr = 1
     call MP_perr_die(myname_, 'invalid value for attrib', ierr)
  end select

 end subroutine aVaVSharedAttrIndexList_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: aCaCSharedAttrIndexList_ - AttrVect shared attributes.
!
! !DESCRIPTION:  {\tt aCaCSharedAttrIndexList\_()} takes a pair of 
! user-supplied {\tt Accumulator} variables {\tt aC1} and {\tt aC2}, 
! and for choice of either {\tt REAL} or {\tt INTEGER} attributes (as
! specified literally in the input {\tt CHARACTER} argument {\tt attrib})
! returns the number of shared attributes {\tt NumShared}, and arrays of
! indices {\tt Indices1} and {\tt Indices2} to their storage locations
! in {\tt aC1} and {\tt aC2}, respectively.
!
! {\bf N.B.:}  This routine returns two allocated arrays---{\tt Indices1(:)} 
! and {\tt Indices2(:)}---which must be deallocated once the user no longer
! needs them.  Failure to do this will create a memory leak.
!
! !INTERFACE:

 subroutine aCaCSharedAttrIndexList_(aC1, aC2, attrib, NumShared, &
                                     Indices1, Indices2)

!
! !USES:
!
      use m_stdio
      use m_die,         only : MP_perr_die

      use m_Accumulator, only : Accumulator
      use m_Accumulator, only : Accumulator_lsize => lsize

      implicit none

! !INPUT PARAMETERS: 
!
      type(Accumulator),    intent(in)  :: aC1   
      type(Accumulator),    intent(in)  :: aC2
      character*7,          intent(in)  :: attrib

! !OUTPUT PARAMETERS:   
!
      integer,           intent(out) :: NumShared

      integer,dimension(:), pointer  :: Indices1
      integer,dimension(:), pointer  :: Indices2

! !REVISION HISTORY:
!       07Feb01 - J.W. Larson <larson@mcs.anl.gov> - initial version
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::aCaCSharedAttrIndexList_'

  integer :: ierr

       ! Based on the value of the argument attrib, pass the 
       ! appropriate pair of Lists for comparison...

  select case(trim(attrib))
  case('REAL','real')
     call GetSharedListIndices(aC1%av%rList, aC2%av%rList, NumShared, &
	                         Indices1, Indices2)
  case('INTEGER','integer')
     call GetSharedListIndices(aC1%av%iList, aC2%av%iList, NumShared, &
	                         Indices1, Indices2)
  case default
     write(stderr,'(4a)') myname_,":: value of argument attrib=",attrib, &
	  " not recognized.  Allowed values: REAL, real, INTEGER, integer"
     ierr = 1
     call MP_perr_die(myname_, 'invalid value for attrib', ierr)
  end select

 end subroutine aCaCSharedAttrIndexList_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: aVaCSharedAttrIndexList_ - AttrVect/Accumulator shared attributes.
!
! !DESCRIPTION:  {\tt aVaCSharedAttrIndexList\_()} a user-supplied 
! {\tt AttrVect} variable {\tt aV} and an {\tt Accumulator} variable 
! {\tt aC}, and for choice of either {\tt REAL} or {\tt INTEGER} 
! attributes (as ! specified literally in the input {\tt CHARACTER} 
! argument {\tt attrib}) returns the number of shared attributes 
! {\tt NumShared}, and arrays of indices {\tt Indices1} and {\tt Indices2} 
! to their storage locations in {\tt aV} and {\tt aC}, respectively.
!
! {\bf N.B.:}  This routine returns two allocated arrays---{\tt Indices1(:)} 
! and {\tt Indices2(:)}---which must be deallocated once the user no longer
! needs them.  Failure to do this will create a memory leak.
!
! !INTERFACE:

 subroutine aVaCSharedAttrIndexList_(aV, aC, attrib, NumShared, &
                                     Indices1, Indices2)

!
! !USES:
!
      use m_stdio
      use m_die,         only : MP_perr_die

      use m_AttrVect,    only : AttrVect
      use m_AttrVect,    only : AttrVect_lsize => lsize
      use m_Accumulator, only : Accumulator
      use m_Accumulator, only : Accumulator_lsize => lsize
 
      implicit none

! !INPUT PARAMETERS: 
!
      type(AttrVect),    intent(in)  :: aV   
      type(Accumulator), intent(in)  :: aC
      character*7,       intent(in)  :: attrib

! !OUTPUT PARAMETERS:   
!
      integer,           intent(out) :: NumShared

      integer,dimension(:), pointer  :: Indices1
      integer,dimension(:), pointer  :: Indices2

! !REVISION HISTORY:
!       07Feb01 - J.W. Larson <larson@mcs.anl.gov> - initial version
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::aVaCSharedAttrIndexList_'

  integer :: ierr

       ! Based on the value of the argument attrib, pass the 
       ! appropriate pair of Lists for comparison...

  select case(trim(attrib))
  case('REAL','real')
     call GetSharedListIndices(aV%rList, aC%av%rList, NumShared, &
	                         Indices1, Indices2)
  case('INTEGER','integer')
     call GetSharedListIndices(aV%iList, aC%av%iList, NumShared, &
	                         Indices1, Indices2)
  case default
     write(stderr,'(4a)') myname_,":: value of argument attrib=",attrib, &
	  " not recognized.  Allowed values: REAL, real, INTEGER, integer"
     ierr = 1
     call MP_perr_die(myname_, 'invalid value for attrib', ierr)
  end select

 end subroutine aVaCSharedAttrIndexList_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: GetSharedListIndices_
!
! !DESCRIPTION:  {\tt GetSharedListIndices\_()} compares two user-
! supplied {\tt List} arguments {\tt List1} and {\tt Lis2} to determine:  
! the number of shared items {\tt NumShared}, and arrays of the locations 
! {\tt Indices1} and {\tt Indices2} in {\tt List1} and {\tt List2}, 
! respectively.
!
! {\bf N.B.:}  This routine returns two allocated arrays---{\tt Indices1(:)} 
! and {\tt Indices2(:)}---which must be deallocated once the user no longer
! needs them.  Failure to do this will create a memory leak.
!
! !INTERFACE:

 subroutine GetSharedListIndices_(List1, List2, NumShared, Indices1, &
                                   Indices2)

!
! !USES:
!
      use m_die,  only : MP_perr_die

      use m_String, only : String
      use m_String, only : String_clean => clean

      use m_List, only : List
      use m_List, only : List_nitem => nitem
      use m_List, only : List_get => get
      use m_List, only : List_index => index
 
      implicit none

! !INPUT PARAMETERS: 
!
      type(List),    intent(in)  :: List1
      type(List),    intent(in)  :: List2

! !OUTPUT PARAMETERS:   
!
      integer,           intent(out) :: NumShared

      integer,dimension(:), pointer  :: Indices1
      integer,dimension(:), pointer  :: Indices2

! !REVISION HISTORY:
!       07Feb01 - J.W. Larson <larson@mcs.anl.gov> - initial version
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::GetSharedListIndices_'

! Error flag
  integer :: ierr

! number of items in List1 and List2, respectively:
  integer :: nitem1, nitem2

! MAXIMUM number of matches possible:
  integer :: NumSharedMax

! Temporary storage for a string tag retrieved from a list:
  type(String) :: tag

! Loop counters / temporary indices:
  integer :: n1, n2

       ! Determine the number of items in each list:

  nitem1 = List_nitem(List1)
  nitem2 = List_nitem(List2)

       ! The maximum number of list item matches possible
       ! is the minimum(nitem1,nitem2):

  NumSharedMax = min(nitem1,nitem2)

       ! Allocate sufficient space for the matches we may find:

  allocate(Indices1(NumSharedMax), Indices2(NumSharedMax), stat=ierr)

       ! Initialize the counter for the number of matches found:

  NumShared = 0

       ! Scan through the two lists.  For the sake of speed, loop 
       ! over the shorter of the two lists...

  if(nitem1 <= nitem2) then ! List1 is shorter--scan it...

     do n1=1,NumSharedMax

       ! Retrieve string tag n1 from List1:
	call List_get(tag, n1, List1)

       ! Index this tag WRT List2--a nonzero value signifies a match
	n2 = List_index(List2, tag)

       ! Clear out tag for the next iteration...
	call String_clean(tag)

       ! If we have a hit, update NumShared, and load the indices
       ! n1 and n2 in Indices1 and Indices2, respectively...

	if((0 < n2) .and. (n2 <= nitem2)) then
	   NumShared = NumShared + 1
	   Indices1(NumShared) = n1
	   Indices2(NumShared) = n2
	endif

     end do ! do n1=1,NumSharedMax

  else ! List1 is shorter--scan it...

     do n2=1,NumSharedMax

       ! Retrieve string tag n2 from List2:
	call List_get(tag, n2, List2)

       ! Index this tag WRT List1--a nonzero value signifies a match
	n1 = List_index(List1, tag)

       ! Clear out tag for the next iteration...
	call String_clean(tag)

       ! If we have a hit, update NumShared, and load the indices
       ! n1 and n2 in Indices1 and Indices2, respectively...

	if((0 < n1) .and. (n1 <= nitem1)) then
	   NumShared = NumShared + 1
	   Indices1(NumShared) = n1
	   Indices2(NumShared) = n2
	endif

     end do ! do n2=1,NumSharedMax

  endif ! if(nitem1 <= nitem2)...

 end subroutine GetSharedListIndices_
 
 end module m_SharedAttrIndices
