!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_GlobalIntegral - Spatial Integrals and Averaging.
!
! !DESCRIPTION:
!
! !INTERFACE:

 module m_GlobalIntegral

      implicit none

      private	! except

! !PUBLIC MEMBER FUNCTIONS:

      public :: GlobalIntegral        ! Global Integral
      public :: GlobalAverage         ! Global Area Average

      public :: MaskedGlobalIntegral  ! Masked Global Integral
      public :: MaskedGlobalAverage   ! MaskedGlobal Area Average

      public :: PairedGlobalIntegrals ! A Pair of Global 
                                      ! Integrals 

      public :: PairedGlobalAverages  ! A Pair of Global 
                                      ! Area Averages

      interface GlobalIntegral ; module procedure &
	   GlobalIntegralRAttrGG_, &
	   GlobalIntegralRAttrV_
      end interface
      interface GlobalAverage ; module procedure &
	   GlobalAverageRAttrGG_, &
	   GlobalAverageRAttrV_
      end interface
      interface MaskedGlobalIntegral ; module procedure &
	   MaskedGlobalIntegralRAttrV_, &
	   MaskedGlobalIntegralRAttrGG_
      end interface
      interface MaskedGlobalAverage ; module procedure &
	   MaskedGlobalAverageRAttrV_, &
	   MaskedGlobalAverageRAttrGG_
      end interface
      interface PairedGlobalIntegrals ; module procedure &
	    PairedGlobalIntegralRAttrGG_, &
	    PairedGlobalIntegralRAttrV_
      end interface
      interface PairedGlobalAverages ; module procedure &
	   PairedGlobalAverageRAttrGG_, &
	   PairedGlobalAverageRAttrV_
      end interface

! !REVISION HISTORY:
! 	25Oct01 - J.W. Larson <larson@mcs.anl.gov> - Initial version
!        9May02 - J.W. Larson <larson@mcs.anl.gov> - Massive Refactoring.
!
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_GlobalIntegral'

 contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: GlobalIntegralRAttrGG_ - Compute global spatial integral.
!
! !DESCRIPTION:
! This routine computes spatial integrals of the {\tt REAL} attributes
! of the {\tt REAL} attributes of the input {\tt AttrVect} argument 
! {\tt inAv}.  {\tt GlobalIntegralRAttrGG\_()} takes the input 
! {\tt AttrVect} argument {\tt inAv} and computes the global spatial 
! integral using weights stored in the {\tt GeneralGrid} argument 
! {\tt GGrid} and identified by the {\tt CHARACTER} tag {\tt WeightTag}.  
! The integral of each {\tt REAL} attribute is returned in the output 
! {\tt AttrVect} argument {\tt outAv}. If {\tt GlobalIntegralRAttrGG\_()} 
! is invoked with the optional {\tt LOGICAL} input argument 
! {\tt SumWeights} set as {\tt .TRUE.}, then the weights are also summed 
! and stored in {\tt outAv} (and can be referenced with the attribute 
! tag defined by the argument{\tt WeightTag}.  If 
! {\tt GlobalIntegralRAttrGG\_()} is invoked with the optional {\tt INTEGER} 
! argument {\tt comm} (a Fortran MPI communicator handle), the summation
! operations for the integral are completed on the local process, then 
! reduced across the communicator, with all processes receiving the result.
!
! {\bf N.B.:  } The local lengths of the {\tt AttrVect} argument {\tt inAv} 
! and the {\tt GeneralGrid} {\tt GGrid} must be equal.  That is, there 
! must be a one-to-one correspondence between the field point values stored 
! in {\tt inAv} and the point weights stored in {\tt GGrid}.
!
! {\bf N.B.:  }  If {\tt GlobalIntegralRAttrGG\_()} is invoked with the 
! optional {\tt LOGICAL} input argument {\tt SumWeights} set as {\tt .TRUE.}, 
! then the value of {\tt WeightTag} must not conflict with any of the 
! {\tt REAL} attribute tags in {\tt inAv}.
!
! {\bf N.B.:  } The output {\tt AttrVect} argument {\tt outAv} is an 
! allocated data structure.  The user must deallocate it using the routine 
! {\tt AttrVect\_clean()} when it is no longer needed.  Failure to do so 
! will result in a memory leak.
!
! !INTERFACE:

 subroutine GlobalIntegralRAttrGG_(inAv, outAv, GGrid, WeightTag, &
                                   SumWeights, comm)
! ! USES:

      use m_stdio
      use m_die
      use m_mpif90

      use m_AttrVect, only : AttrVect
      use m_AttrVect, only : AttrVect_lsize => lsize

      use m_GeneralGrid, only : GeneralGrid
      use m_GeneralGrid, only : GeneralGrid_lsize => lsize
      use m_GeneralGrid, only : GeneralGrid_indexRA => indexRA
      use m_GeneralGrid, only : GeneralGrid_exportRAttr => exportRAttr

      implicit none

! !INPUT PARAMETERS:

      type(AttrVect),    intent(IN) :: inAv
      type(GeneralGrid), intent(IN) :: GGrid
      character(len=*),  intent(IN) :: WeightTag
      logical, optional, intent(IN) :: SumWeights
      integer, optional, intent(IN) :: comm

! !OUTPUT PARAMETERS:

      type(AttrVect),    intent(OUT) :: outAv

! !REVISION HISTORY:
! 	06Feb02 - J.W. Larson <larson@mcs.anl.gov> - initial version
! 	09May02 - J.W. Larson <larson@mcs.anl.gov> - Refactored and
!                 renamed GlobalIntegralRAttrGG_().
!       07Jun02 - J.W. Larson <larson@mcs.anl.gov> - Bug fix and further
!                 refactoring.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::GlobalIntegralRAttrGG_'

  integer :: ierr, length
  logical :: mySumWeights
  real, dimension(:), pointer :: gridWeights

       ! Argument Validity Checks

  if(AttrVect_lsize(inAv) /= GeneralGrid_lsize(GGrid)) then
     ierr = AttrVect_lsize(inAv) - GeneralGrid_lsize(GGrid)
     write(stderr,'(2a,i8,a,i8)') myname_, &
	  ':: inAv / GGrid length mismatch:  ', &
	  ' AttrVect_lsize(inAv) = ',AttrVect_lsize(inAv), &
	  ' GeneralGrid_lsize(GGrid) = ',GeneralGrid_lsize(GGrid)
     call die(myname_)
  endif

  if(present(SumWeights)) then
     mySumWeights = SumWeights
  else
     mySumWeights = .FALSE.
  endif

       ! ensure unambiguous pointer association status for gridWeights

  nullify(gridWeights) 

       ! Extract Grid Weights

  call GeneralGrid_exportRAttr(GGrid, WeightTag, gridWeights, length)

       ! 

  if(present(comm)) then ! do a distributed AllReduce-style integral:
     call GlobalIntegralRAttrV_(inAv, outAv, gridWeights, mySumWeights, &
	                        WeightTag, comm)
  else
     call GlobalIntegralRAttrV_(inAv, outAv, gridWeights, mySumWeights, &
                                WeightTag)
  endif

       ! Clean up temporary allocated space

  deallocate(gridWeights, stat=ierr)
  if(ierr /= 0) then
     write(stderr,'(2a,i8)') myname_, &
	              ':: deallocate(gridWeights...failed.  ierr=', ierr
     call die(myname_)
  endif

 end subroutine GlobalIntegralRAttrGG_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: GlobalIntegralRAttrV_ - Compute global spatial integral.
!
! !DESCRIPTION:
! This routine computes spatial integrals of the {\tt REAL} attributes
! of the {\tt REAL} attributes of the input {\tt AttrVect} argument 
! {\tt inAv}.  {\tt GlobalIntegralRAttrV\_()} takes the input 
! {\tt AttrVect} argument {\tt inAv} and computes the global spatial 
! integral using weights stored in the input {\tt REAL} array argument 
! {\tt Weights}.  The integral of each {\tt REAL} attribute is returned 
! in the output {\tt AttrVect} argument {\tt outAv}. If 
! {\tt GlobalIntegralRAttrV\_()} is invoked with the optional {\tt LOGICAL} 
! input argument {\tt SumWeights} set as {\tt .TRUE.}, then the weights 
! are also summed and stored in {\tt outAv} (and can be referenced with 
! the attribute name {\tt WeightTag}.  If {\tt GlobalIntegralRAttrV\_()} is 
! invoked with the optional {\tt INTEGER} argument {\tt comm} (a Fortran 
! MPI communicator handle), the summation operations for the integral are 
! completed on the local process, then reduced across the communicator, 
! with all processes receiving the result.
!
! {\bf N.B.:  } The local lengths of the {\tt AttrVect} argument {\tt inAv} 
! and the input array {\tt Weights} must be equal.  That is, there must be 
! a one-to-one correspondence between the field point values stored 
! in {\tt inAv} and the point weights stored in {\tt Weights}.
!
! {\bf N.B.:  }  If {\tt GlobalIntegralRAttrV\_()} is invoked with the 
! optional {\tt LOGICAL} input argument {\tt SumWeights} set as {\tt .TRUE.}. 
! In this case, the none of {\tt REAL} attribute tags in {\tt inAv} may be 
! named the same as the string contained in {\tt WeightTag}, which is an 
! attribute name reserved for the sum of the weights in the output {\tt AttrVect} 
! {\tt outAv}.
!
! {\bf N.B.:  } The output {\tt AttrVect} argument {\tt outAv} is an 
! allocated data structure.  The user must deallocate it using the routine 
! {\tt AttrVect\_clean()} when it is no longer needed.  Failure to do so 
! will result in a memory leak.
!
! !INTERFACE:

 subroutine GlobalIntegralRAttrV_(inAv, outAv, Weights, SumWeights, &
                                  WeightTag, comm)

! ! USES:

      use m_stdio
      use m_die
      use m_mpif90

      use m_AttrVect, only : AttrVect
      use m_AttrVect, only : AttrVect_lsize => lsize

      use m_AttrVectReduce, only : AttrVect_GlobalWeightedSumRAttr => &
	                                         GlobalWeightedSumRAttr
      use m_AttrVectReduce, only : AttrVect_LocalWeightedSumRAttr => &
	                                         LocalWeightedSumRAttr

      implicit none

! !INPUT PARAMETERS:

      type(AttrVect),               intent(IN) :: inAv
      real, dimension(:),           pointer    :: Weights
      logical,            optional, intent(IN) :: SumWeights
      character(len=*),   optional, intent(IN) :: WeightTag
      integer,            optional, intent(IN) :: comm

! !OUTPUT PARAMETERS:

      type(AttrVect),    intent(OUT) :: outAv

! !REVISION HISTORY:
! 	07Jun02 - J.W. Larson <larson@mcs.anl.gov> - initial version
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::GlobalIntegralRAttrV_'

  integer :: ierr, length
  logical :: mySumWeights

       ! Argument Validity Checks

  if(AttrVect_lsize(inAv) /= size(Weights)) then
     ierr = AttrVect_lsize(inAv) - size(Weights)
     write(stderr,'(2a,i8,a,i8)') myname_, &
	  ':: inAv / Weights array length mismatch:  ', &
	  ' AttrVect_lsize(inAv) = ',AttrVect_lsize(inAv), &
	  ' size(Weights) = ',size(Weights)
     call die(myname_)
  endif

  if(present(SumWeights)) then
     mySumWeights = SumWeights
     if(.not. present(WeightTag)) then
	write(stderr,*) myname_,':: if the input argument SumWeights=.TRUE.,', &
                        ' then the argument WeightTag must be provided.'
	call die(myname_)
     endif
  else
     mySumWeights = .FALSE.
  endif

       ! Compute the sum

  if(present(comm)) then ! compute distributed AllReduce-style sum:

     if(mySumWeights) then ! return the global sum of the weights in outAV
	call AttrVect_GlobalWeightedSumRAttr(inAV, outAV, Weights, &
	     comm, WeightTag)
     else
	call AttrVect_GlobalWeightedSumRAttr(inAV, outAV, Weights, comm)
     endif

  else ! compute local sum:

     if(mySumWeights) then ! return the global sum of the weights in outAV
	call AttrVect_LocalWeightedSumRAttr(inAV, outAV, Weights, &
                                            WeightTag)
     else
	call AttrVect_LocalWeightedSumRAttr(inAV, outAV, Weights)
     endif

  endif ! if(present(comm))...

 end subroutine GlobalIntegralRAttrV_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: GlobalAverageRAttrGG_ - Compute global spatial average.
!
! !DESCRIPTION:
! This routine computes spatial averages of the {\tt REAL} attributes
! of the input {\tt AttrVect} argument {\tt inAv}.  
! {\tt GlobalAverageRAttrGG\_()} takes the input {\tt AttrVect} argument 
! {\tt inAv} and computes the global spatial average using weights 
! stored in the {\tt GeneralGrid} argument {\tt GGrid} and identified by
! the {\tt CHARACTER} tag {\tt WeightTag}.  The average of each {\tt REAL}
! attribute is returned in the output {\tt AttrVect} argument {\tt outAv}.
! If {\tt GlobalAverageRAttrGG\_()} is invoked with the optional {\tt INTEGER} 
! argument {\tt comm} (a Fortran MPI communicator handle), the summation
! operations for the average are completed on the local process, then 
! reduced across the communicator, with all processes receiving the result.
!
! {\bf N.B.:  } The local lengths of the {\tt AttrVect} argument {\tt inAv} 
! and the {\tt GeneralGrid} {\tt GGrid} must be equal.  That is, there 
! must be a one-to-one correspondence between the field point values stored 
! in {\tt inAv} and the point weights stored in {\tt GGrid}.
!
! {\bf N.B.:  } The output {\tt AttrVect} argument {\tt outAv} is an 
! allocated data structure.  The user must deallocate it using the routine 
! {\tt AttrVect\_clean()} when it is no longer needed.  Failure to do so 
! will result in a memory leak.
!
! !INTERFACE:

 subroutine GlobalAverageRAttrGG_(inAv, outAv, GGrid, WeightTag, comm)

! ! USES:

      use m_stdio
      use m_die
      use m_mpif90

      use m_AttrVect, only : AttrVect
      use m_AttrVect, only : AttrVect_init => init
      use m_AttrVect, only : AttrVect_clean => clean
      use m_AttrVect, only : AttrVect_nRAttr => nRAttr
      use m_AttrVect, only : AttrVect_indexRA => indexRA
      use m_AttrVect, only : AttrVect_exportRListToChar => exportRListToChar

      use m_GeneralGrid, only : GeneralGrid

      implicit none

! !INPUT PARAMETERS:

      type(AttrVect),    intent(IN) :: inAv
      type(GeneralGrid), intent(IN) :: GGrid
      character(len=*),  intent(IN) :: WeightTag
      integer, optional, intent(IN) :: comm

! !OUTPUT PARAMETERS:

      type(AttrVect),    intent(OUT) :: outAv

! !REVISION HISTORY:
! 	08Feb02 - J.W. Larson <larson@mcs.anl.gov> - initial version
! 	08May02 - J.W. Larson <larson@mcs.anl.gov> - minor modifications:
!                 1) renamed the routine to GlobalAverageRAttrGG_
!                 2) changed calls to reflect new routine name
!                    GlobalIntegralRAttrGG_().
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::GlobalAverageRAtttrGG_'

  type(AttrVect) :: integratedAv
  integer :: i, ierr, iweight

       ! Compute the global integral:

  if(present(comm)) then
     call GlobalIntegralRAttrGG_(inAv, integratedAv, GGrid, WeightTag, &
	                         .TRUE., comm)
  else
     call GlobalIntegralRAttrGG_(inAv, integratedAv, GGrid, WeightTag, &
	                         .TRUE.)
  endif

       ! Check value of summed weights (to avoid division by zero):

  iweight = AttrVect_indexRA(integratedAv, WeightTag)
  if(integratedAv%rAttr(iweight, 1) == 0) then
     write(stderr,'(2a)') myname_, &
	  '::ERROR--Global sum of grid weights is zero.'
     call die(myname_)
  endif

       ! Initialize output AttrVect outAv:

  call AttrVect_init(outAv, rList=AttrVect_exportRListToChar(inAv), lsize=1)

       ! Divide by global weight sum to compute global averages from 
       ! global integrals.

  do i=1,AttrVect_nRAttr(outAv)
     outAv%rAttr(i,1) = integratedAv%rAttr(i,1) &
	                               / integratedAv%rAttr(iweight,1) 
  end do

       ! Clean up temporary AttrVect:

  call AttrVect_clean(integratedAv)

 end subroutine GlobalAverageRAttrGG_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: GlobalAverageRAttrV_ - Compute global spatial average.
!
! !DESCRIPTION:
! This routine computes spatial averages of the {\tt REAL} attributes
! of the input {\tt AttrVect} argument {\tt inAv}.  
! {\tt GlobalAverageRAttrV\_()} takes the input {\tt AttrVect} argument 
! {\tt inAv} and computes the global spatial average using weights 
! stored in the {\tt REAL} array {\tt Weights}.  The average of each 
! {\tt REAL} attribute is returned in the output {\tt AttrVect} argument 
! {\tt outAv}.  If {\tt GlobalAverageRAttrGG\_()} is invoked with the 
! optional {\tt INTEGER} argument {\tt comm} (a Fortran MPI communicator 
! handle), the summation operations for the average are completed on the 
! local process, then reduced across the communicator, with all processes 
! receiving the result.
!
! {\bf N.B.:  } The local lengths of the {\tt AttrVect} argument {\tt inAv} 
! and the input array {\tt Weights} must be equal.  That is, there must 
! be a one-to-one correspondence between the field point values stored 
! in {\tt inAv} and the point weights stored in {\tt Weights}.
!
! {\bf N.B.:  } The output {\tt AttrVect} argument {\tt outAv} is an 
! allocated data structure.  The user must deallocate it using the routine 
! {\tt AttrVect\_clean()} when it is no longer needed.  Failure to do so 
! will result in a memory leak.
!
! !INTERFACE:

 subroutine GlobalAverageRAttrV_(inAv, outAv, Weights, comm)

! ! USES:

      use m_stdio
      use m_die
      use m_mpif90

      use m_AttrVect, only : AttrVect
      use m_AttrVect, only : AttrVect_init => init
      use m_AttrVect, only : AttrVect_clean => clean
      use m_AttrVect, only : AttrVect_nRAttr => nRAttr
      use m_AttrVect, only : AttrVect_indexRA => indexRA
      use m_AttrVect, only : AttrVect_exportRListToChar => exportRListToChar

      use m_GeneralGrid, only : GeneralGrid

      implicit none

! !INPUT PARAMETERS:

      type(AttrVect),     intent(IN)  :: inAv
      real, dimension(:), pointer     :: Weights
      integer, optional,  intent(IN)  :: comm

! !OUTPUT PARAMETERS:

      type(AttrVect),     intent(OUT) :: outAv

! !REVISION HISTORY:
! 	10Jun02 - J.W. Larson <larson@mcs.anl.gov> - initial version
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::GlobalAverageRAtttrV_'

  type(AttrVect) :: integratedAv
  integer :: i, ierr, iweight

       ! Nullify all the pointers present in the AttrVect integratedAv:

  call AttrVect_clean(integratedAv)

       ! Compute the global integral:

  if(present(comm)) then
     call GlobalIntegralRAttrV_(inAv, integratedAv, Weights, &
	                         .TRUE., 'weights', comm)
  else
     call GlobalIntegralRAttrV_(inAv, integratedAv, Weights, .TRUE., 'weights')
  endif

       ! Check value of summed weights (to avoid division by zero):

  iweight = AttrVect_indexRA(integratedAv, 'weights')
  if(integratedAv%rAttr(iweight, 1) == 0.) then
     write(stderr,'(2a)') myname_, &
	  '::ERROR--Global sum of grid weights is zero.'
     call die(myname_)
  endif

       ! Initialize output AttrVect outAv:

  call AttrVect_init(outAv, rList=AttrVect_exportRListToChar(inAv), lsize=1)

       ! Divide by global weight sum to compute global averages from 
       ! global integrals.

  do i=1,AttrVect_nRAttr(outAv)
     outAv%rAttr(i,1) = integratedAv%rAttr(i,1) &
	                               / integratedAv%rAttr(iweight,1) 
  end do

       ! Clean up temporary AttrVect:

  call AttrVect_clean(integratedAv)

 end subroutine GlobalAverageRAttrV_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: MaskedGlobalIntegralRAttrGG_ - Masked global spatial integral.
!
! !DESCRIPTION: [NEEDS **LOTS** of work...]
! This routine computes masked spatial integrals of the {\tt REAL} 
! attributes of the input {\tt AttrVect} argument {\tt inAv}.  
!
!
! {\bf N.B.:  } The local lengths of the {\tt AttrVect} argument {\tt inAv} 
! and the input array {\tt Weights} must be equal.  That is, there must be 
! a one-to-one correspondence between the field point values stored 
! in {\tt inAv} and the point weights stored in {\tt Weights}.
!
! {\bf N.B.:  }  If {\tt GlobalIntegralRAttrV\_()} is invoked with the 
! optional {\tt LOGICAL} input argument {\tt SumWeights} set as {\tt .TRUE.}. 
! In this case, the none of {\tt REAL} attribute tags in {\tt inAv} may be 
! named the same as the string contained in {\tt WeightTag}, which is an 
! attribute name reserved for the sum of the weights in the output {\tt AttrVect} 
! {\tt outAv}.
!
! {\bf N.B.:  } The output {\tt AttrVect} argument {\tt outAv} is an 
! allocated data structure.  The user must deallocate it using the routine 
! {\tt AttrVect\_clean()} when it is no longer needed.  Failure to do so 
! will result in a memory leak.
!
! !INTERFACE:

 subroutine MaskedGlobalIntegralRAttrGG_(inAv, outAv, GGrid, SpatialWeightTag, &
                                         iMaskTags, rMaskTags, UseFastMethod,  &
					 SumWeights, WeightSumTag, comm)

! ! USES:

      use m_stdio
      use m_die
      use m_mpif90

      use m_String, only : String
      use m_String, only : String_toChar => toChar
      use m_String, only : String_clean => clean

      use m_List, only : List
      use m_List, only : List_init => init
      use m_List, only : List_clean => clean
      use m_List, only : List_nitem => nitem
      use m_List, only : List_get => get

      use m_AttrVect, only : AttrVect
      use m_AttrVect, only : AttrVect_lsize => lsize

      use m_GeneralGrid, only : GeneralGrid
      use m_GeneralGrid, only : GeneralGrid_lsize => lsize
      use m_GeneralGrid, only : GeneralGrid_indexRA => indexRA
      use m_GeneralGrid, only : GeneralGrid_exportIAttr => exportIAttr
      use m_GeneralGrid, only : GeneralGrid_exportRAttr => exportRAttr

      use m_AttrVectReduce, only : AttrVect_GlobalWeightedSumRAttr => &
	                                         GlobalWeightedSumRAttr
      use m_AttrVectReduce, only : AttrVect_LocalWeightedSumRAttr => &
	                                         LocalWeightedSumRAttr
      implicit none

! !INPUT PARAMETERS:

      type(AttrVect),                  intent(IN) :: inAv
      type(GeneralGrid),               intent(IN) :: GGrid
      character(len=*),                intent(IN) :: SpatialWeightTag
      character(len=*),      optional, intent(IN) :: iMaskTags
      character(len=*),      optional, intent(IN) :: rMaskTags
      logical,                         intent(IN) :: UseFastMethod
      logical,               optional, intent(IN) :: SumWeights
      character(len=*),      optional, intent(IN) :: WeightSumTag
      integer,               optional, intent(IN) :: comm

! !OUTPUT PARAMETERS:

      type(AttrVect),                  intent(OUT) :: outAv

! !REVISION HISTORY:
! 	11Jun02 - J.W. Larson <larson@mcs.anl.gov> - initial version
!
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::MaskedGlobalIntegralRAttrGG_'

  integer :: i, ierr, j, length
  logical :: mySumWeights

  type(List) :: iMaskList, rMaskList
  type(String) :: DummStr

  integer, dimension(:), pointer :: iMask, iMaskTemp
  real, dimension(:), pointer :: rMask, rMaskTemp
  integer :: TempMaskLength

  real, dimension(:), pointer :: SpatialWeights

  integer :: niM, nrM ! Number of iMasks and rMasks, respectively

       ! Argument Validity Checks

  if(AttrVect_lsize(inAv) /= GeneralGrid_lsize(GGrid)) then
     ierr = AttrVect_lsize(inAv) - GeneralGrid_lsize(GGrid)
     write(stderr,'(2a,i8,a,i8)') myname_, &
	  ':: inAv / GGrid length mismatch:  ', &
	  ' AttrVect_lsize(inAv) = ',AttrVect_lsize(inAv), &
	  ' GeneralGrid_lsize(GGrid) = ',GeneralGrid_lsize(GGrid)
     call die(myname_)
  endif

  if(present(SumWeights)) then
     mySumWeights = SumWeights
     if(.not. present(WeightSumTag)) then
	write(stderr,*) myname_,':: if the input argument SumWeights=.TRUE.,', &
                        ' then the argument WeightSumTag must be provided.'
	call die(myname_)
     endif
  else
     mySumWeights = .FALSE.
  endif

  if(present(iMaskTags)) then
     call List_init(iMaskList, iMaskTags)
     if(List_nitem(iMaskList) == 0) then
	write(stderr,'(3a)') myname_,':: ERROR--an INTEGER mask list with', &
	               'no valid items was provided.'
        call die(myname_)
     endif
  endif

  if(present(rMaskTags)) then
     call List_init(rMaskList, rMaskTags)
     if(List_nitem(iMaskList) == 0) then
	write(stderr,'(3a)') myname_,':: ERROR--an REAL mask list with', &
	               'no valid items was provided.'
        call die(myname_)
     endif
  endif

       ! Determine the on-processor vector length for use throughout
       ! this routine:

  length = AttrVect_lsize(inAv)

       !==========================================================
       ! Extract Spatial Weights from GGrid using SpatialWeightTag
       !==========================================================

  nullify(SpatialWeights)
  call GeneralGrid_exportRAttr(GGrid, SpatialWeightTag, SpatialWeights, &
                               TempMaskLength)
  if(TempMaskLength /= length) then
     write(stderr,'(3a,i8,a,i8)') myname_,&
	  ':: error on return from GeneralGrid_exportRAttr().' , &   
	  'Returned with SpatialWeights(:) length = ',TempMaskLength, &
	  ',which conflicts with AttrVect_lsize(inAv) = ',length
     call die(myname_)
  endif

       !==========================================================
       ! If the argument iMaskTags is present, create the combined
       ! iMask array:
       !==========================================================

  if(present(iMaskTags)) then ! assemble iMask(:) from all the integer
                              ! mask attributes stored in GGrid(:)

     allocate(iMask(length), iMaskTemp(length), stat=ierr)
     if(ierr /= 0) then
	write(stderr,'(3a,i8)') myname_,':: allocate(iMask(...) failed,', &
	     ' ierr=',ierr
	call die(myname_)
     endif

     niM = List_nitem(iMaskList)

     do i=1,niM

       ! Retrieve current iMask tag, and get this attribute from GGrid:
	call List_get(DummStr, i, iMaskList)
	call GeneralGrid_exportIAttr(GGrid, String_toChar(DummStr), &
	                             iMaskTemp, TempMaskLength)
	call String_clean(DummStr)
	if(TempMaskLength /= length) then
	   write(stderr,'(3a,i8,a,i8)') myname_,&
		 ':: error on return from GeneralGrid_exportIAttr().' , &   
		 'Returned with TempMaskLength = ',TempMaskLength, &
		 ',which conflicts with AttrVect_lsize(inAv) = ',length
	   call die(myname_)
	endif

	if(i == 1) then ! first pass--examine iMaskTemp(:) only

	   if(UseFastMethod) then ! straight copy of iMaskTemp(:)
	      do j=1,length
		 iMask(j) = iMaskTemp(j)
	      end do
	   else ! go through the entries of iMaskTemp(:) one-by-one
	      do j=1,length
		 select case(iMaskTemp(j))
		 case(0)
		    iMask(j) = 0
		 case(1)
		    iMask(j) = 1
		 case default
		    write(stderr,'(3a,i8,a,i8)') myname_, &
			 ':: FATAL--illegal INTEGER mask entry.  Integer mask ', &
			 'entries must be 0 or 1.  iMask(',j,') = ', iMask(j)
		    call die(myname_)
		 end select ! select case(iMaskTemp(j))...
	      end do ! do j=1,length
	   endif ! if(UseFastMethod)...

	else ! That is, i /= 1 ...

	   if(UseFastMethod) then ! straight product of iMask(:) 
                                  ! and iMaskTemp(:)
	      do j=1,length
		 iMask(j) = iMask(j) * iMaskTemp(j)
	      end do
	   else ! go through the entries of iMaskTemp(:) one-by-one
	      do j=1,length
		 select case(iMaskTemp(j))
		 case(0) ! zero out iMask(j)
		    iMask(j) = 0
		 case(1) ! do nothing
		 case default
		    write(stderr,'(3a,i8,a,i8)') myname_, &
			 ':: FATAL--illegal INTEGER mask entry.  Integer mask ', &
			 'entries must be 0 or 1.  iMask(',j,') = ', iMask(j)
		    call die(myname_)
		 end select ! select case(iMaskTemp(j))...
	      end do ! do j=1,length
	   endif ! if(UseFastMethod)...

	endif ! if(i == 1)...

     end do ! do i=1,niM...iMask retrievals

  endif ! if(present(iMaskTags))...

       !==========================================================
       ! If the argument rMaskTags is present, create the combined
       ! REAL mask rMask array:
       !==========================================================

  if(present(rMaskTags)) then ! assemble rMask(:) from all the integer
                              ! mask attributes stored in GGrid(:)

     allocate(rMask(length), rMaskTemp(length), stat=ierr)
     if(ierr /= 0) then
	write(stderr,'(3a,i8)') myname_,':: allocate(rMask(...) failed,', &
	     ' ierr=',ierr
	call die(myname_)
     endif

     nrM = List_nitem(rMaskList)

     do i=1,nrM

       ! Retrieve current rMask tag, and get this attribute from GGrid:
	call List_get(DummStr, i, rMaskList)
	call GeneralGrid_exportRAttr(GGrid, String_toChar(DummStr), &
	                             rMaskTemp, TempMaskLength)
	call String_clean(DummStr)
	if(TempMaskLength /= length) then
	   write(stderr,'(3a,i8,a,i8)') myname_,&
		 ':: error on return from GeneralGrid_exportRAttr().' , &   
		 'Returned with TempMaskLength = ',TempMaskLength, &
		 ',which conflicts with AttrVect_lsize(inAv) = ',length
	   call die(myname_)
	endif

	if(i == 1) then ! first pass--examine rMaskTemp(:) only

	   if(UseFastMethod) then ! straight copy of rMaskTemp(:)
	      do j=1,length
		 rMask(j) = rMaskTemp(j)
	      end do
	   else ! go through the entries of rMaskTemp(:) one-by-one
                ! to ensure they are in the range [0.,1.]
	      do j=1,length
		 if((rMaskTemp(j) >= 0.) .or. (rMaskTemp(j) <=1.)) then
		    rMask(j) = rMaskTemp(j)
		 else
		    write(stderr,'(3a,i8,a,i8)') myname_, &
			 ':: FATAL--illegal REAL mask entry.  Real mask ', &
			 'entries must be in [0.,1.]  rMask(',j,') = ', rMask(j)
		    call die(myname_)
		 endif ! if((rMaskTemp(j) >= 0.) .or. (rMaskTemp(j) <=1.))...
	      end do ! do j=1,length
	   endif ! if(UseFastMethod)...

	else ! That is, i /= 1 ...

	   if(UseFastMethod) then ! straight product of rMask(:) 
                                  ! and rMaskTemp(:)
	      do j=1,length
		 rMask(j) = rMask(j) * rMaskTemp(j)
	      end do
	   else ! go through the entries of rMaskTemp(:) one-by-one
                ! to ensure they are in the range [0.,1.]
	      do j=1,length
		 if((rMaskTemp(j) >= 0.) .or. (rMaskTemp(j) <=1.)) then
		    rMask(j) = rMask(j) * rMaskTemp(j)
		 else
		    write(stderr,'(3a,i8,a,i8)') myname_, &
			 ':: FATAL--illegal REAL mask entry.  Real mask ', &
			 'entries must be in [0.,1.]  rMask(',j,') = ', rMask(j)
		    call die(myname_)
		 endif ! if((rMaskTemp(j) >= 0.) .or. (rMaskTemp(j) <=1.))...
	      end do ! do j=1,length
	   endif ! if(UseFastMethod)...

	endif ! if(i == 1)...

     end do ! do i=1,niM...rMask retrievals

  endif ! if(present(rMaskTags))...

       !==========================================================
       ! Now that we have produced single INTEGER and REAL masks,
       ! compute the masked weighted sum.
       !==========================================================

  if(present(rMaskTags)) then ! We have a REAL Mask

     if(present(iMaskTags)) then ! and an INTEGER Mask

	if(present(comm)) then ! compute distributed AllReduce-style sum:
	   
	   if(mySumWeights) then ! return the global masked sum of the 
                                 ! weights in outAV
	      call MaskedGlobalIntegralRAttrV_(inAv, outAv, SpatialWeights, &
		                               iMask, rMask, UseFastMethod, &
					       SumWeights, WeightSumTag, comm)
	   else ! Do not return the masked sum of the weights
	      call MaskedGlobalIntegralRAttrV_(inAv, outAv, SpatialWeights, &
		                               iMask, rMask, UseFastMethod, &
					       comm=comm)
	   endif ! if(mySumWeights)...

	else ! compute local sum:

	   if(mySumWeights) then ! return the global masked sum of the 
                                 ! weights in outAV
	      call MaskedGlobalIntegralRAttrV_(inAv, outAv, SpatialWeights, &
		                               iMask, rMask, UseFastMethod, &
					       SumWeights, WeightSumTag)
	   else ! Do not return the masked sum of the weights
	      call MaskedGlobalIntegralRAttrV_(inAv, outAv, SpatialWeights, &
		                               iMask, rMask, UseFastMethod)
	   endif ! if(mySumWeights)...

	endif ! if(present(comm))...

     else ! REAL Mask Only Case...

	if(present(comm)) then ! compute distributed AllReduce-style sum:
	   
	   if(mySumWeights) then ! return the global masked sum of the 
                                 ! weights in outAV
	      call MaskedGlobalIntegralRAttrV_(inAv, outAv, SpatialWeights, &
		                               rMask=rMask, &
					       UseFastMethod=UseFastMethod, &
					       SumWeights=SumWeights, &
					       WeightSumTag=WeightSumTag, &
					       comm=comm)
	   else ! Do not return the masked sum of the weights
	      call MaskedGlobalIntegralRAttrV_(inAv, outAv, SpatialWeights, &
		                               rMask=rMask, &
					       UseFastMethod=UseFastMethod, &
					       comm=comm)
	   endif ! if(mySumWeights)...

	else ! compute local sum:

	   if(mySumWeights) then ! return the global masked sum of the 
                                 ! weights in outAV
	      call MaskedGlobalIntegralRAttrV_(inAv, outAv, SpatialWeights, &
		                               rMask=rMask, &
					       UseFastMethod=UseFastMethod, &
					       SumWeights=SumWeights, &
					       WeightSumTag=WeightSumTag)
	   else ! Do not return the masked sum of the weights
	      call MaskedGlobalIntegralRAttrV_(inAv, outAv, SpatialWeights, &
		                               rMask=rMask, &
					       UseFastMethod=UseFastMethod)
	   endif ! if(mySumWeights)...

	endif ! if(present(comm))...

     endif
  else ! no REAL Mask...

     if(present(iMaskTags)) then ! INTEGER Mask Only Case...

	if(present(comm)) then ! compute distributed AllReduce-style sum:
	   
	   if(mySumWeights) then ! return the global masked sum of the 
                                 ! weights in outAV
	      call MaskedGlobalIntegralRAttrV_(inAv, outAv, SpatialWeights, &
		                               iMask=iMask, &
					       UseFastMethod=UseFastMethod, &
					       SumWeights=SumWeights, &
					       WeightSumTag=WeightSumTag, &
					       comm=comm)
	   else ! Do not return the masked sum of the weights
	      call MaskedGlobalIntegralRAttrV_(inAv, outAv, SpatialWeights, &
		                               iMask=iMask, &
					       UseFastMethod=UseFastMethod, &
					       comm=comm)
	   endif ! if(mySumWeights)...

	else ! compute local sum:

	   if(mySumWeights) then ! return the global masked sum of the 
                                 ! weights in outAV
	      call MaskedGlobalIntegralRAttrV_(inAv, outAv, SpatialWeights, &
		                               iMask=iMask, &
					       UseFastMethod=UseFastMethod, &
					       SumWeights=SumWeights, &
					       WeightSumTag=WeightSumTag)
	   else ! Do not return the masked sum of the weights
	      call MaskedGlobalIntegralRAttrV_(inAv, outAv, SpatialWeights, &
		                               iMask=iMask, &
					       UseFastMethod=UseFastMethod)
	   endif ! if(mySumWeights)...

	endif ! if(present(comm))...

     else ! no INTEGER Mask / no REAL Mask Case...

	if(present(comm)) then ! compute distributed AllReduce-style sum:
	   
	   if(mySumWeights) then ! return the global masked sum of the 
                                 ! weights in outAV
	      call MaskedGlobalIntegralRAttrV_(inAv, outAv, SpatialWeights, &
					       UseFastMethod=UseFastMethod, &
					       SumWeights=SumWeights, &
					       WeightSumTag=WeightSumTag, &
					       comm=comm)
	   else ! Do not return the masked sum of the weights
	      call MaskedGlobalIntegralRAttrV_(inAv, outAv, SpatialWeights, &
					       UseFastMethod=UseFastMethod, &
					       comm=comm)
	   endif ! if(mySumWeights)...

	else ! compute local sum:

	   if(mySumWeights) then ! return the global masked sum of the 
                                 ! weights in outAV
	      call MaskedGlobalIntegralRAttrV_(inAv, outAv, SpatialWeights, &
					       UseFastMethod=UseFastMethod, &
					       SumWeights=SumWeights, &
					       WeightSumTag=WeightSumTag)
	   else ! Do not return the masked sum of the weights
	      call MaskedGlobalIntegralRAttrV_(inAv, outAv, SpatialWeights, &
					       UseFastMethod=UseFastMethod)
	   endif ! if(mySumWeights)...

	endif ! if(present(comm))...

     endif ! if(present(iMaskTags)...

  endif ! if(present(rMaskTags)...

       !==========================================================
       ! The masked spatial integral is now completed.
       ! Clean up the the various allocated mask structures.
       !==========================================================

  if(present(iMaskTags)) then ! clean up iMask and friends...
     call List_clean(iMaskList)
     deallocate(iMask, iMaskTemp, stat=ierr)
     if(ierr /= 0) then
	write(stderr,'(3a,i8)') myname_,':: deallocate(iMask(...) failed,', &
	     ' ierr=',ierr
	call die(myname_)
     endif
  endif

  if(present(rMaskTags)) then ! clean up rMask and co...
     call List_clean(rMaskList)
     deallocate(rMask, rMaskTemp, stat=ierr)
     if(ierr /= 0) then
	write(stderr,'(3a,i8)') myname_,':: deallocate(rMask(...) failed,', &
	     ' ierr=',ierr
	call die(myname_)
     endif
  endif

       ! Clean up SpatialWeights(:)

  deallocate(SpatialWeights, stat=ierr)
  if(ierr /= 0) then
     write(stderr,'(3a,i8)') myname_,':: deallocate(SpatialWeights(...) failed,', &
	                    ' ierr=',ierr
     call die(myname_)
  endif

 end subroutine MaskedGlobalIntegralRAttrGG_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: MaskedGlobalIntegralRAttrV_ - Masked global spatial integral.
!
! !DESCRIPTION: [NEEDS **LOTS** of work...]
! This routine computes masked spatial integrals of the {\tt REAL} 
! attributes of the input {\tt AttrVect} argument {\tt inAv}.  
!
! {\tt AttrVect} argument {\tt inAv} and computes the global spatial 
! integral using weights stored in the input {\tt REAL} array argument 
! {\tt Weights}.  The integral of each {\tt REAL} attribute is returned 
! in the output {\tt AttrVect} argument {\tt outAv}. If 
! {\tt GlobalIntegralRAttrV\_()} is invoked with the optional {\tt LOGICAL} 
! input argument {\tt SumWeights} set as {\tt .TRUE.}, then the weights 
! are also summed and stored in {\tt outAv} (and can be referenced with 
! the attribute name {\tt WeightTag}.  If {\tt GlobalIntegralRAttrV\_()} is 
! invoked with the optional {\tt INTEGER} argument {\tt comm} (a Fortran 
! MPI communicator handle), the summation operations for the integral are 
! completed on the local process, then reduced across the communicator, 
! with all processes receiving the result.
!
! {\bf N.B.:  } The local lengths of the {\tt AttrVect} argument {\tt inAv} 
! and the input array {\tt Weights} must be equal.  That is, there must be 
! a one-to-one correspondence between the field point values stored 
! in {\tt inAv} and the point weights stored in {\tt Weights}.
!
! {\bf N.B.:  }  If {\tt GlobalIntegralRAttrV\_()} is invoked with the 
! optional {\tt LOGICAL} input argument {\tt SumWeights} set as {\tt .TRUE.}. 
! In this case, the none of {\tt REAL} attribute tags in {\tt inAv} may be 
! named the same as the string contained in {\tt WeightTag}, which is an 
! attribute name reserved for the sum of the weights in the output {\tt AttrVect} 
! {\tt outAv}.
!
! {\bf N.B.:  } The output {\tt AttrVect} argument {\tt outAv} is an 
! allocated data structure.  The user must deallocate it using the routine 
! {\tt AttrVect\_clean()} when it is no longer needed.  Failure to do so 
! will result in a memory leak.
!
! !INTERFACE:

 subroutine MaskedGlobalIntegralRAttrV_(inAv, outAv, SpatialWeights, iMask, &
                                        rMask, UseFastMethod, SumWeights, &
                                        WeightSumTag, comm)

! ! USES:

      use m_stdio
      use m_die
      use m_mpif90

      use m_AttrVect, only : AttrVect
      use m_AttrVect, only : AttrVect_lsize => lsize

      use m_AttrVectReduce, only : AttrVect_GlobalWeightedSumRAttr => &
	                                         GlobalWeightedSumRAttr
      use m_AttrVectReduce, only : AttrVect_LocalWeightedSumRAttr => &
	                                         LocalWeightedSumRAttr
      implicit none

! !INPUT PARAMETERS:

      type(AttrVect),                  intent(IN) :: inAv
      real,    dimension(:),           pointer    :: SpatialWeights
      integer, dimension(:), optional, pointer    :: iMask
      real,    dimension(:), optional, pointer    :: rMask
      logical,                         intent(IN) :: UseFastMethod
      logical,               optional, intent(IN) :: SumWeights
      character(len=*),      optional, intent(IN) :: WeightSumTag
      integer,               optional, intent(IN) :: comm

! !OUTPUT PARAMETERS:

      type(AttrVect),                  intent(OUT) :: outAv

! !REVISION HISTORY:
! 	10Jun02 - J.W. Larson <larson@mcs.anl.gov> - initial version
!
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::MaskedGlobalIntegralRAttrV_'

  integer :: i, ierr, length
  logical :: mySumWeights
  real, dimension(:), pointer :: Weights

       ! Argument Validity Checks

  if(AttrVect_lsize(inAv) /= size(SpatialWeights)) then
     ierr = AttrVect_lsize(inAv) - size(SpatialWeights)
     write(stderr,'(2a,i8,a,i8)') myname_, &
	  ':: inAv / SpatialWeights array length mismatch:  ', &
	  ' AttrVect_lsize(inAv) = ',AttrVect_lsize(inAv), &
	  ' size(SpatialWeights) = ',size(SpatialWeights)
     call die(myname_)
  endif

  if(present(iMask)) then ! make sure it is the right length
     if(AttrVect_lsize(inAv) /= size(iMask)) then
	ierr = AttrVect_lsize(inAv) - size(iMask)
	write(stderr,'(2a,i8,a,i8)') myname_, &
	     ':: inAv / iMask array length mismatch:  ', &
	     ' AttrVect_lsize(inAv) = ',AttrVect_lsize(inAv), &
	     ' size(iMask) = ',size(iMask)
	call die(myname_)
     endif
  endif

  if(present(rMask)) then ! make sure it is the right length
     if(AttrVect_lsize(inAv) /= size(rMask)) then
	ierr = AttrVect_lsize(inAv) - size(rMask)
	write(stderr,'(2a,i8,a,i8)') myname_, &
	     ':: inAv / rMask array length mismatch:  ', &
	     ' AttrVect_lsize(inAv) = ',AttrVect_lsize(inAv), &
	     ' size(rMask) = ',size(rMask)
	call die(myname_)
     endif
  endif

  if(present(SumWeights)) then
     mySumWeights = SumWeights
     if(.not. present(WeightSumTag)) then
	write(stderr,*) myname_,':: if the input argument SumWeights=.TRUE.,', &
                        ' then the argument WeightSumTag must be provided.'
	call die(myname_)
     endif
  else
     mySumWeights = .FALSE.
  endif

       ! Create a common Weights(:) array...

  length = AttrVect_lsize(inAv)

  allocate(Weights(length), stat=ierr)
  if(ierr /= 0) then
     write(stderr,'(3a,i8)') myname_,':: allocate(Weights(...) failed,', &
	  ' ierr=',ierr
     call die(myname_)
  endif

       ! Combine weights and masks into a common Weights(:) array...

  if(UseFastMethod) then ! form the product of iMask, rMask, and SpatialWeights

     if(present(rMask)) then ! use it to form Weights(:)
	if(present(iMask)) then ! use it and rMask to form Weights(:)
	   do i=1,length
	      Weights(i) = rMask(i) * SpatialWeights(i) * iMask(i)
	   end do
	else
	   do i=1,length
	      Weights(i) = rMask(i) * SpatialWeights(i)
	   end do
	endif ! if(present(iMask))...
     else
	if(present(iMask)) then
	   do i=1,length
	      Weights(i) = SpatialWeights(i) * iMask(i)
	   end do
	else
	   do i=1,length
	      Weights(i) = SpatialWeights(i)
	   end do	   
	endif ! if(present(iMask))...
     endif ! if(present(rMask))...


  else ! Scan iMask and rMask carefully and set Weights(i) to zero
       ! when iMask(i) or rMask(i) is zero.  This avoids round-off
       ! effects from products and promotion of integers to reals.

     if(present(rMask)) then ! use it to form Weights(:)
	if(present(iMask)) then ! use it and rMask to form Weights(:)
	   do i=1,length
	      select case(iMask(i))
	      case(0)
		 Weights(i) = 0.
	      case(1)
		 if(rMask(i) == 1.) then
		    Weights(i) = SpatialWeights(i)
		 elseif(rMask(i) == 0.) then
		    Weights(i) = 0.
		 elseif((rMask(i) > 0.) .and. (rMask(i) < 1.)) then
		    Weights(i) = rMask(i) * SpatialWeights(i)
		 else ! rMask(i) < 0. or rMask(i) > 1.
		    write(stderr,'(3a,i8,a,f10.7)') myname_, &
			 ':: invalid value for real', &
			 'mask entry rMask(',i,') = ',rMask(i)
		    call die(myname_)
                 endif
	      case default
		 write(stderr,'(3a,i8,a,i8)') myname_, &
		      ':: invalid value for integer', &
		      'mask entry iMask(',i,') = ',iMask(i)
		 call die(myname_)
	      end select
	   end do
	else
	   do i=1,length
	      if(rMask(i) == 1.) then
		 Weights(i) = SpatialWeights(i)
	      elseif(rMask(i) == 0.) then
		 Weights(i) = 0.
	      elseif((rMask(i) > 0.) .and. (rMask(i) < 1.)) then
		 Weights(i) = rMask(i) * SpatialWeights(i)
	      else ! rMask(i) < 0. or rMask(i) > 1.
		 write(stderr,'(3a,i8,a,f10.7)') myname_, &
		      ':: invalid value for real', &
		      'mask entry rMask(',i,') = ',rMask(i)
		 call die(myname_)
	      endif
	   end do
	endif ! if(present(iMask))...
     else ! no rMask present...
	if(present(iMask)) then ! check iMask entries...
	   do i=1,length
	      select case(iMask(i))
	      case(0)
		 Weights(i) = 0.
	      case(1)
		 Weights(i) = SpatialWeights(i)
	      case default
		 write(stderr,'(3a,i8,a,i8)') myname_, &
		      ':: invalid value for integer', &
		      'mask entry iMask(',i,') = ',iMask(i)
		 call die(myname_)		 
	      end select
	   end do
	else ! straight assignment of SpatialWeights(:) 
	   do i=1,length
		 Weights(i) = SpatialWeights(i)
	   end do	   
	endif ! if(present(iMask))...
     endif ! if(present(rMask))...


  endif ! if(UseFastMethod)

       ! Now that the weights are combined into a common Weights(:),
       ! compute the masked weighted sum:

  if(present(comm)) then ! compute distributed AllReduce-style sum:
	
     if(mySumWeights) then ! return the global sum of the weights in outAV
	call AttrVect_GlobalWeightedSumRAttr(inAV, outAV, Weights, &
   	                                     comm, WeightSumTag)
     else
	call AttrVect_GlobalWeightedSumRAttr(inAV, outAV, Weights, comm)
     endif

  else ! compute local sum:

     if(mySumWeights) then ! return the global sum of the weights in outAV
	call AttrVect_LocalWeightedSumRAttr(inAV, outAV, Weights, &
	                                    WeightSumAttr=WeightSumTag)
     else
	call AttrVect_LocalWeightedSumRAttr(inAV, outAV, Weights)
     endif

  endif ! if(present(comm))...

       ! Clean up the allocated Weights(:) array

  deallocate(Weights, stat=ierr)
  if(ierr /= 0) then
     write(stderr,'(3a,i8)') myname_,':: deallocate(Weights(...) failed,', &
	  ' ierr=',ierr
     call die(myname_)
  endif

 end subroutine MaskedGlobalIntegralRAttrV_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: MaskedGlobalAverageRAttrGG_ - Masked global spatial average.
!
! !DESCRIPTION: [NEEDS **LOTS** of work...]
! This routine computes masked spatial averages of the {\tt REAL} 
! attributes of the input {\tt AttrVect} argument {\tt inAv}.  
!
!
! {\bf N.B.:  } The local lengths of the {\tt AttrVect} argument {\tt inAv} 
! and the input array {\tt Weights} must be equal.  That is, there must be 
! a one-to-one correspondence between the field point values stored 
! in {\tt inAv} and the point weights stored in {\tt Weights}.
!
! {\bf N.B.:  }  If {\tt GlobalAverageRAttrV\_()} is invoked with the 
! optional {\tt LOGICAL} input argument {\tt SumWeights} set as {\tt .TRUE.}. 
! In this case, the none of {\tt REAL} attribute tags in {\tt inAv} may be 
! named the same as the string contained in {\tt WeightTag}, which is an 
! attribute name reserved for the sum of the weights in the output {\tt AttrVect} 
! {\tt outAv}.
!
! {\bf N.B.:  } The output {\tt AttrVect} argument {\tt outAv} is an 
! allocated data structure.  The user must deallocate it using the routine 
! {\tt AttrVect\_clean()} when it is no longer needed.  Failure to do so 
! will result in a memory leak.
!
! !INTERFACE:

 subroutine MaskedGlobalAverageRAttrGG_(inAv, outAv, GGrid, SpatialWeightTag, &
                                        iMaskTags, rMaskTags, UseFastMethod,  &
  				        comm)

! ! USES:

      use m_stdio
      use m_die
      use m_mpif90

      use m_AttrVect, only : AttrVect
      use m_AttrVect, only : AttrVect_init => init
      use m_AttrVect, only : AttrVect_clean => clean
      use m_AttrVect, only : AttrVect_lsize => lsize
      use m_AttrVect, only : AttrVect_indexRA => indexRA
      use m_AttrVect, only : AttrVect_nRAttr => nRAttr
      use m_AttrVect, only : AttrVect_exportRListToChar => &
	                                                exportRListToChar 

      use m_GeneralGrid, only : GeneralGrid
      use m_GeneralGrid, only : GeneralGrid_lsize => lsize
      use m_GeneralGrid, only : GeneralGrid_indexRA => indexRA

      implicit none

! !INPUT PARAMETERS:

      type(AttrVect),                  intent(IN) :: inAv
      type(GeneralGrid),               intent(IN) :: GGrid
      character(len=*),                intent(IN) :: SpatialWeightTag
      character(len=*),      optional, intent(IN) :: iMaskTags
      character(len=*),      optional, intent(IN) :: rMaskTags
      logical,                         intent(IN) :: UseFastMethod
      integer,               optional, intent(IN) :: comm

! !OUTPUT PARAMETERS:

      type(AttrVect),                  intent(OUT) :: outAv

! !REVISION HISTORY:
! 	12Jun02 - J.W. Larson <larson@mcs.anl.gov> - initial version
!
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::MaskedGlobalAverageRAttrGG_'

  type(AttrVect) :: integratedAv
  character*9, parameter :: WeightSumTag = 'WeightSum'

  integer :: i, iweight

       !================================================================
       ! Do the integration using MaskedGlobalIntegralRAttrGG_(), which
       ! returns the intermediate integrals (including the masked weight
       ! sum) in the AttrVect integratedAv.
       !================================================================

  if(present(iMaskTags)) then

     if(present(rMaskTags)) then ! have both iMasks and rMasks

	if(present(comm)) then ! a distributed parallel sum
	   call MaskedGlobalIntegralRAttrGG_(inAv, integratedAv, GGrid, &
		                             SpatialWeightTag, iMaskTags, &
					     rMaskTags, UseFastMethod,  &
					     .TRUE., WeightSumTag, comm)
	else ! a purely local sum
	   call MaskedGlobalIntegralRAttrGG_(inAv, integratedAv, GGrid, &
		                             SpatialWeightTag, iMaskTags, &
					     rMaskTags, UseFastMethod,  &
					     .TRUE., WeightSumTag)
	endif ! if(present(comm))...

     else ! Only iMasks are in use

	if(present(comm)) then ! a distributed parallel sum
	   call MaskedGlobalIntegralRAttrGG_(inAv, integratedAv, GGrid, &
		                             SpatialWeightTag, iMaskTags, &
					     UseFastMethod=UseFastMethod,  &
					     SumWeights=.TRUE., &
					     WeightSumTag=WeightSumTag, &
					     comm=comm)

	else ! a purely local sum
	   call MaskedGlobalIntegralRAttrGG_(inAv, integratedAv, GGrid, &
		                             SpatialWeightTag, iMaskTags, &
					     UseFastMethod=UseFastMethod,  &
					     SumWeights=.TRUE., &
					     WeightSumTag=WeightSumTag)
	endif ! if(present(comm))...

     endif ! if(present(rMaskTags)...

  else ! no iMasks

     if(present(rMaskTags)) then ! Only rMasks are in use

	if(present(comm)) then ! a distributed parallel sum
	   call MaskedGlobalIntegralRAttrGG_(inAv, integratedAv, GGrid, &
		                             SpatialWeightTag, &
					     rMaskTags=rMaskTags, &
					     UseFastMethod=UseFastMethod,  &
					     SumWeights=.TRUE., &
					     WeightSumTag=WeightSumTag, &
					     comm=comm)
	else ! a purely local sum
	   call MaskedGlobalIntegralRAttrGG_(inAv, integratedAv, GGrid, &
		                             SpatialWeightTag, &
					     rMaskTags=rMaskTags, &
					     UseFastMethod=UseFastMethod,  &
					     SumWeights=.TRUE., &
					     WeightSumTag=WeightSumTag)
	endif

     else ! Neither iMasks nor rMasks are in use

	if(present(comm)) then ! a distributed parallel sum
	   call MaskedGlobalIntegralRAttrGG_(inAv, integratedAv, GGrid, &
		                             SpatialWeightTag, &
					     UseFastMethod=UseFastMethod,  &
					     SumWeights=.TRUE., &
					     WeightSumTag=WeightSumTag, &
					     comm=comm)
	else ! a purely local sum
	   call MaskedGlobalIntegralRAttrGG_(inAv, integratedAv, GGrid, &
		                             SpatialWeightTag, &
					     UseFastMethod=UseFastMethod,  &
					     SumWeights=.TRUE., &
					     WeightSumTag=WeightSumTag)
	endif ! if(present(comm))...

     endif ! if(present(rMaskTags))...

  endif ! if(present(iMaskTags))...

       !================================================================
       ! The masked integrals and masked weight sum now reside in 
       ! in the AttrVect integratedAv.  We now wish to compute the
       ! averages by dividing the integtrals by the masked weight sum.
       !================================================================

       ! Check value of summed weights (to avoid division by zero):

  iweight = AttrVect_indexRA(integratedAv, WeightSumTag)
  if(integratedAv%rAttr(iweight, 1) == 0.) then
     write(stderr,'(2a)') myname_, &
	  '::ERROR--Global sum of grid weights is zero.'
     call die(myname_)
  endif

       ! Initialize output AttrVect outAv:

  call AttrVect_init(outAv, rList=AttrVect_exportRListToChar(inAv), lsize=1)

       ! Divide by global weight sum to compute global averages from 
       ! global integrals.

  do i=1,AttrVect_nRAttr(outAv)
     outAv%rAttr(i,1) = integratedAv%rAttr(i,1) &
	                               / integratedAv%rAttr(iweight,1) 
  end do

       ! Clean up temporary AttrVect:

  call AttrVect_clean(integratedAv)

 end subroutine MaskedGlobalAverageRAttrGG_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: MaskedGlobalAverageRAttrV_ - Masked global spatial average.
!
! !DESCRIPTION: [NEEDS **LOTS** of work...]
! This routine computes spatial integrals of the {\tt REAL} attributes
! of the {\tt REAL} attributes of the input {\tt AttrVect} argument 
! {\tt inAv}.  {\tt GlobalIntegralRAttrV\_()} takes the input 
! {\tt AttrVect} argument {\tt inAv} and computes the global spatial 
! integral using weights stored in the input {\tt REAL} array argument 
! {\tt Weights}.  The integral of each {\tt REAL} attribute is returned 
! in the output {\tt AttrVect} argument {\tt outAv}. If 
! {\tt GlobalIntegralRAttrV\_()} is invoked with the optional {\tt LOGICAL} 
! input argument {\tt SumWeights} set as {\tt .TRUE.}, then the weights 
! are also summed and stored in {\tt outAv} (and can be referenced with 
! the attribute name {\tt WeightTag}.  If {\tt GlobalIntegralRAttrV\_()} is 
! invoked with the optional {\tt INTEGER} argument {\tt comm} (a Fortran 
! MPI communicator handle), the summation operations for the integral are 
! completed on the local process, then reduced across the communicator, 
! with all processes receiving the result.
!
! {\bf N.B.:  } The local lengths of the {\tt AttrVect} argument {\tt inAv} 
! and the input array {\tt Weights} must be equal.  That is, there must be 
! a one-to-one correspondence between the field point values stored 
! in {\tt inAv} and the point weights stored in {\tt Weights}.
!
! {\bf N.B.:  }  If {\tt GlobalIntegralRAttrV\_()} is invoked with the 
! optional {\tt LOGICAL} input argument {\tt SumWeights} set as {\tt .TRUE.}. 
! In this case, the none of {\tt REAL} attribute tags in {\tt inAv} may be 
! named the same as the string contained in {\tt WeightTag}, which is an 
! attribute name reserved for the sum of the weights in the output {\tt AttrVect} 
! {\tt outAv}.
!
! {\bf N.B.:  } The output {\tt AttrVect} argument {\tt outAv} is an 
! allocated data structure.  The user must deallocate it using the routine 
! {\tt AttrVect\_clean()} when it is no longer needed.  Failure to do so 
! will result in a memory leak.
!
! !INTERFACE:

 subroutine MaskedGlobalAverageRAttrV_(inAv, outAv, SpatialWeights, iMask, &
                                       rMask, UseFastMethod, comm)

! ! USES:

      use m_stdio
      use m_die
      use m_mpif90

      use m_AttrVect, only : AttrVect
      use m_AttrVect, only : AttrVect_init => init
      use m_AttrVect, only : AttrVect_clean => clean
      use m_AttrVect, only : AttrVect_lsize => lsize
      use m_AttrVect, only : AttrVect_nRAttr => nRAttr
      use m_AttrVect, only : AttrVect_indexRA => indexRA
      use m_AttrVect, only : AttrVect_exportRListToChar => exportRListToChar

      implicit none

! !INPUT PARAMETERS:

      type(AttrVect),                  intent(IN) :: inAv
      real,    dimension(:),           pointer    :: SpatialWeights
      integer, dimension(:), optional, pointer    :: iMask
      real,    dimension(:), optional, pointer    :: rMask
      logical,                         intent(IN) :: UseFastMethod
      integer,               optional, intent(IN) :: comm

! !OUTPUT PARAMETERS:

      type(AttrVect),                  intent(OUT) :: outAv

! !REVISION HISTORY:
! 	11Jun02 - J.W. Larson <larson@mcs.anl.gov> - initial version
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::MaskedGlobalAverageRAttrV_'

  type(AttrVect) :: integratedAv

  integer :: i, ierr, length, iweight
  logical :: mySumWeights

       ! Argument Validity Checks

  if(AttrVect_lsize(inAv) /= size(SpatialWeights)) then
     ierr = AttrVect_lsize(inAv) - size(SpatialWeights)
     write(stderr,'(2a,i8,a,i8)') myname_, &
	  ':: inAv / SpatialWeights array length mismatch:  ', &
	  ' AttrVect_lsize(inAv) = ',AttrVect_lsize(inAv), &
	  ' size(SpatialWeights) = ',size(SpatialWeights)
     call die(myname_)
  endif

  if(present(iMask)) then ! make sure it is the right length
     if(AttrVect_lsize(inAv) /= size(iMask)) then
	ierr = AttrVect_lsize(inAv) - size(iMask)
	write(stderr,'(2a,i8,a,i8)') myname_, &
	     ':: inAv / iMask array length mismatch:  ', &
	     ' AttrVect_lsize(inAv) = ',AttrVect_lsize(inAv), &
	     ' size(iMask) = ',size(iMask)
	call die(myname_)
     endif
  endif

  if(present(rMask)) then ! make sure it is the right length
     if(AttrVect_lsize(inAv) /= size(rMask)) then
	ierr = AttrVect_lsize(inAv) - size(rMask)
	write(stderr,'(2a,i8,a,i8)') myname_, &
	     ':: inAv / rMask array length mismatch:  ', &
	     ' AttrVect_lsize(inAv) = ',AttrVect_lsize(inAv), &
	     ' size(rMask) = ',size(rMask)
	call die(myname_)
     endif
  endif

       ! Nullify all the pointers present in the AttrVect integratedAv:

  call AttrVect_clean(integratedAv)

       ! Compute the masked weighted sum, including the sum of the 
       ! masked weights.

  if(present(comm)) then ! communicator handle present

     if(present(iMask)) then

	if(present(rMask)) then
	   call MaskedGlobalIntegralRAttrV_(inAv, integratedAv, SpatialWeights,  &
		                            iMask, rMask, UseFastMethod, .TRUE., &
					    'MaskedWeightsSum', comm)
	else ! no rMask
	   call MaskedGlobalIntegralRAttrV_(inAv, integratedAv, SpatialWeights,  &
		                            iMask=iMask, UseFastMethod=UseFastMethod, &
                                            SumWeights=.TRUE., &
					    WeightSumTag='MaskedWeightsSum', &
					    comm=comm)
	endif ! if(present(rMask))...

     else ! no iMask present...

	if(present(rMask)) then
	   call MaskedGlobalIntegralRAttrV_(inAv, integratedAv, SpatialWeights,  &
		                            rMask=rMask, UseFastMethod=UseFastMethod, &
					    SumWeights=.TRUE., &
					    WeightSumTag='MaskedWeightsSum', &
					    comm=comm)
	else ! neither rMask nor iMask present:
	   call MaskedGlobalIntegralRAttrV_(inAv, integratedAv, SpatialWeights,  &
		                            UseFastMethod=UseFastMethod, &
					    SumWeights=.TRUE., &
					    WeightSumTag='MaskedWeightsSum', &
					    comm=comm)
	endif ! if(present(rMask))...

     endif ! if(present(iMask))...

  else ! no communicator handle present

     if(present(iMask)) then

	if(present(rMask)) then
	   call MaskedGlobalIntegralRAttrV_(inAv, integratedAv, SpatialWeights,  &
		                            iMask, rMask, UseFastMethod, .TRUE., &
					    'MaskedWeightsSum')
	else ! no rMask
	   call MaskedGlobalIntegralRAttrV_(inAv, integratedAv, SpatialWeights,  &
		                            iMask=iMask, UseFastMethod=UseFastMethod, &
                                            SumWeights=.TRUE., &
					    WeightSumTag='MaskedWeightsSum')
	endif ! if(present(rMask))...

     else ! no iMask present...

	if(present(rMask)) then
	   call MaskedGlobalIntegralRAttrV_(inAv, integratedAv, SpatialWeights,  &
		                            rMask=rMask, UseFastMethod=UseFastMethod, &
					    SumWeights=.TRUE., &
					    WeightSumTag='MaskedWeightsSum')
	else ! neither rMask nor iMask present:
	   call MaskedGlobalIntegralRAttrV_(inAv, integratedAv, SpatialWeights,  &
		                            UseFastMethod=UseFastMethod, &
					    SumWeights=.TRUE., &
					    WeightSumTag='MaskedWeightsSum')
	endif ! if(present(rMask))...

     endif ! if(present(iMask))...

  endif ! if(present(comm))...

       ! At this point, integratedAv containes the masked global integrals
       ! of the REAL attributes of inAv, along with the sum of the weights.
       ! to compute the masked spatial average

       ! Check value of summed weights (to avoid division by zero):

  iweight = AttrVect_indexRA(integratedAv, 'MaskedWeightSum')
  if(integratedAv%rAttr(iweight, 1) == 0.) then
     write(stderr,'(2a)') myname_, &
	  '::ERROR--Global sum of grid weights is zero.'
     call die(myname_)
  endif

       ! Initialize output AttrVect outAv:

  call AttrVect_init(outAv, rList=AttrVect_exportRListToChar(inAv), lsize=1)

       ! Divide by global weight sum to compute global averages from 
       ! global integrals.

  do i=1,AttrVect_nRAttr(outAv)
     outAv%rAttr(i,1) = integratedAv%rAttr(i,1) &
	                               / integratedAv%rAttr(iweight,1) 
  end do

       ! Clean up temporary AttrVect:

  call AttrVect_clean(integratedAv)

 end subroutine MaskedGlobalAverageRAttrV_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: PairedGlobalIntegralRAttrGG_ - Two global spatial integrals.
!
! !DESCRIPTION:
! This routine computes spatial integrals of the {\tt REAL} attributes
! of the {\tt REAL} attributes of the input {\tt AttrVect} argument 
! {\tt inAv}.  {\tt GlobalIntegralRAttrGG\_()} takes the input 
! {\tt AttrVect} argument {\tt inAv} and computes the global spatial 
! integral using weights stored in the {\tt GeneralGrid} argument 
! {\tt GGrid} and identified by the {\tt CHARACTER} tag {\tt WeightTag}.  
! The integral of each {\tt REAL} attribute is returned in the output 
! {\tt AttrVect} argument {\tt outAv}. If {\tt GlobalIntegralRAttrGG\_()} 
! is invoked with the optional {\tt LOGICAL} input argument 
! {\tt SumWeights} set as {\tt .TRUE.}, then the weights are also summed 
! and stored in {\tt outAv} (and can be referenced with the attribute 
! tag defined by the argument {\tt WeightTag}.  If 
! {\tt GlobalIntegralRAttrGG\_()} is invoked with the optional {\tt INTEGER} 
! argument {\tt comm} (a Fortran MPI communicator handle), the summation
! operations for the integral are completed on the local process, then 
! reduced across the communicator, with all processes receiving the result.
!
! {\bf N.B.:  } The local lengths of the {\tt AttrVect} argument {\tt inAv1} 
! and the {\tt GeneralGrid} {\tt GGrid1} must be equal.  That is, there 
! must be a one-to-one correspondence between the field point values stored 
! in {\tt inAv1} and the point weights stored in {\tt GGrid1}.  The same
! relationship must apply between {\tt inAv2} and {\tt GGrid2}.
!
! {\bf N.B.:  }  If {\tt GlobalIntegralRAttrGG\_()} is invoked with the 
! optional {\tt LOGICAL} input argument {\tt SumWeights} set as {\tt .TRUE.}, 
! then the value of {\tt WeightTag1} must not conflict with any of the 
! {\tt REAL} attribute tags in {\tt inAv1} and the value of {\tt WeightTag2} 
! must not conflict with any of the {\tt REAL} attribute tags in {\tt inAv2}.
!
! {\bf N.B.:  } The output {\tt AttrVect} arguments {\tt outAv1} and 
! {\tt outAv2} are allocated data structures.  The user must deallocate them
!  using the routine {\tt AttrVect\_clean()} when they are no longer needed.  
! Failure to do so will result in a memory leak.
!
! !INTERFACE:

 subroutine PairedGlobalIntegralRAttrGG_(inAv1, outAv1, GGrid1, WeightTag1, &
                                         inAv2, outAv2, GGrid2, WeightTag2, &
					 SumWeights, comm)
! ! USES:

      use m_stdio
      use m_die
      use m_mpif90

      use m_AttrVect, only : AttrVect
      use m_AttrVect, only : AttrVect_init => init
      use m_AttrVect, only : AttrVect_lsize => lsize
      use m_AttrVect, only : AttrVect_nRAttr => nRAttr

      use m_GeneralGrid, only : GeneralGrid
      use m_GeneralGrid, only : GeneralGrid_lsize => lsize
      use m_GeneralGrid, only : GeneralGrid_indexRA => indexRA
      use m_GeneralGrid, only : GeneralGrid_exportRAttr => exportRAttr

      use m_AttrVectReduce, only : AttrVect_LocalWeightedSumRAttr => &
	                                         LocalWeightedSumRAttr

      implicit none

! !INPUT PARAMETERS:

      type(AttrVect),    intent(IN) :: inAv1
      type(GeneralGrid), intent(IN) :: GGrid1
      character(len=*),  intent(IN) :: WeightTag1
      type(AttrVect),    intent(IN) :: inAv2
      type(GeneralGrid), intent(IN) :: GGrid2
      character(len=*),  intent(IN) :: WeightTag2
      logical, optional, intent(IN) :: SumWeights
      integer,           intent(IN) :: comm

! !OUTPUT PARAMETERS:

      type(AttrVect),    intent(OUT) :: outAv1
      type(AttrVect),    intent(OUT) :: outAv2

! !REVISION HISTORY:
! 	09May02 - J.W. Larson <larson@mcs.anl.gov> - Initial version.
! 	10Jun02 - J.W. Larson <larson@mcs.anl.gov> - Refactored--now
!                 built on top of PairedIntegralRAttrV_().
!
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::PairedGlobalIntegralRAttrGG_'

       ! Argument Sanity Checks:

  integer :: ierr, length1, length2
  logical :: mySumWeights
  real, dimension(:), pointer :: gridWeights1, gridWeights2

       ! Argument Validity Checks

  if(AttrVect_lsize(inAv1) /= GeneralGrid_lsize(GGrid1)) then
     ierr = AttrVect_lsize(inAv1) - GeneralGrid_lsize(GGrid1)
     write(stderr,'(2a,i8,a,i8)') myname_, &
	  ':: inAv1 / GGrid1 length mismatch:  ', &
	  ' AttrVect_lsize(inAv1) = ',AttrVect_lsize(inAv1), &
	  ' GeneralGrid_lsize(GGrid1) = ',GeneralGrid_lsize(GGrid1)
     call die(myname_)
  endif

  if(AttrVect_lsize(inAv2) /= GeneralGrid_lsize(GGrid2)) then
     ierr = AttrVect_lsize(inAv2) - GeneralGrid_lsize(GGrid2)
     write(stderr,'(2a,i8,a,i8)') myname_, &
	  ':: inAv2 / GGrid2 length mismatch:  ', &
	  ' AttrVect_lsize(inAv2) = ',AttrVect_lsize(inAv2), &
	  ' GeneralGrid_lsize(GGrid2) = ',GeneralGrid_lsize(GGrid2)
     call die(myname_)
  endif

        ! Are we summing the integration weights for either input
        ! GeneralGrid?

  if(present(SumWeights)) then
     mySumWeights = SumWeights
  else
     mySumWeights = .FALSE.
  endif

       ! ensure unambiguous pointer association status for gridWeights1
       ! and gridWeights2

  nullify(gridWeights1) 
  nullify(gridWeights2) 

       ! Extract Grid Weights

  call GeneralGrid_exportRAttr(GGrid1, WeightTag1, gridWeights1, length1)
  call GeneralGrid_exportRAttr(GGrid2, WeightTag2, gridWeights2, length2)


  call PairedGlobalIntegralRAttrV_(inAv1, outAv1, gridweights1, WeightTag1, &
	                           inAv2, outAv2, gridweights2, WeightTag2, &
                                   mySumWeights, comm)

       ! Clean up allocated arrays:

  deallocate(gridWeights1, gridWeights2, stat=ierr)
  if(ierr /= 0) then
     write(stderr,'(2a,i8)') myname_, &
	  'ERROR--deallocate(gridWeights1,...) failed, ierr = ',ierr
     call die(myname_)
  endif

 end subroutine PairedGlobalIntegralRAttrGG_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: PairedGlobalIntegralRAttrV_ - Two global spatial integrals.
!
! !DESCRIPTION: [NEEDS MASSIVE REWRITE]
! This routine computes spatial integrals of the {\tt REAL} attributes
! of the {\tt REAL} attributes of the input {\tt AttrVect} argument 
! {\tt inAv}.  {\tt GlobalIntegralRAttrV\_()} takes the input 
! {\tt AttrVect} argument {\tt inAv} and computes the global spatial 
! integral using weights stored in the {\tt GeneralGrid} argument 
! {\tt GGrid} and identified by the {\tt CHARACTER} tag {\tt WeightName}.  
! The integral of each {\tt REAL} attribute is returned in the output 
! {\tt AttrVect} argument {\tt outAv}. If {\tt GlobalIntegralRAttrGG\_()} 
! is invoked with the optional {\tt LOGICAL} input argument 
! {\tt SumWeights} set as {\tt .TRUE.}, then the weights are also summed 
! and stored in {\tt outAv} (and can be referenced with the attribute 
! tag defined by the argument {\tt WeightName}.  If 
! {\tt GlobalIntegralRAttrGG\_()} is invoked with the optional {\tt INTEGER} 
! argument {\tt comm} (a Fortran MPI communicator handle), the summation
! operations for the integral are completed on the local process, then 
! reduced across the communicator, with all processes receiving the result.
!
! {\bf N.B.:  } The local lengths of the {\tt AttrVect} argument {\tt inAv1} 
! and the {\tt GeneralGrid} {\tt GGrid1} must be equal.  That is, there 
! must be a one-to-one correspondence between the field point values stored 
! in {\tt inAv1} and the point weights stored in {\tt GGrid1}.  The same
! relationship must apply between {\tt inAv2} and {\tt GGrid2}.
!
! {\bf N.B.:  }  If {\tt GlobalIntegralRAttrGG\_()} is invoked with the 
! optional {\tt LOGICAL} input argument {\tt SumWeights} set as {\tt .TRUE.}, 
! then the value of {\tt WeightName1} must not conflict with any of the 
! {\tt REAL} attribute tags in {\tt inAv1} and the value of {\tt WeightName2} 
! must not conflict with any of the {\tt REAL} attribute tags in {\tt inAv2}.
!
! {\bf N.B.:  } The output {\tt AttrVect} arguments {\tt outAv1} and 
! {\tt outAv2} are allocated data structures.  The user must deallocate them
!  using the routine {\tt AttrVect\_clean()} when they are no longer needed.  
! Failure to do so will result in a memory leak.
!
! !INTERFACE:

 subroutine PairedGlobalIntegralRAttrV_(inAv1, outAv1, Weights1, WeightName1, &
                                        inAv2, outAv2, Weights2, WeightName2, &
                                        SumWeights, comm)
! ! USES:

      use m_stdio
      use m_die
      use m_mpif90

      use m_AttrVect, only : AttrVect
      use m_AttrVect, only : AttrVect_lsize => lsize
      use m_AttrVect, only : AttrVect_nRAttr => nRAttr

      use m_GeneralGrid, only : GeneralGrid
      use m_GeneralGrid, only : GeneralGrid_lsize => lsize
      use m_GeneralGrid, only : GeneralGrid_indexRA => indexRA
      use m_GeneralGrid, only : GeneralGrid_exportRAttr => exportRAttr

      use m_AttrVectReduce, only : AttrVect_LocalWeightedSumRAttr => &
	                                         LocalWeightedSumRAttr

      implicit none

! !INPUT PARAMETERS:

      type(AttrVect),     intent(IN)  :: inAv1
      real, dimension(:), pointer     :: Weights1
      character(len=*),   intent(IN)  :: WeightName1
      type(AttrVect),     intent(IN)  :: inAv2
      real, dimension(:), pointer     :: Weights2
      character(len=*),   intent(IN)  :: WeightName2
      logical, optional,  intent(IN)  :: SumWeights
      integer,            intent(IN)  :: comm

! !OUTPUT PARAMETERS:

      type(AttrVect),     intent(OUT) :: outAv1
      type(AttrVect),     intent(OUT) :: outAv2

! !REVISION HISTORY:
! 	10Jun02 - J.W. Larson <larson@mcs.anl.gov> - Initial version.
!
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::PairedGlobalIntegralRAttrV_'

       ! Argument Sanity Checks:

  integer :: ierr, length1, length2, PairedBufferLength
  integer :: nRA1, nRA2
  logical :: mySumWeights
  real, dimension(:), pointer :: PairedBuffer, OutPairedBuffer

       ! Argument Validity Checks

  if(AttrVect_lsize(inAv1) /= size(Weights1)) then
     ierr = AttrVect_lsize(inAv1) - size(Weights1)
     write(stderr,'(2a,i8,a,i8)') myname_, &
	  ':: inAv1 / Weights1 length mismatch:  ', &
	  ' AttrVect_lsize(inAv1) = ',AttrVect_lsize(inAv1), &
	  ' size(Weights1) = ',size(Weights1)
     call die(myname_)
  endif

  if(AttrVect_lsize(inAv2) /= size(Weights2)) then
     ierr = AttrVect_lsize(inAv2) - size(Weights2)
     write(stderr,'(2a,i8,a,i8)') myname_, &
	  ':: inAv2 / Weights2 length mismatch:  ', &
	  ' AttrVect_lsize(inAv2) = ',AttrVect_lsize(inAv2), &
	  ' size(Weights2) = ',size(Weights2)
     call die(myname_)
  endif

        ! Are we summing the integration weights?

  if(present(SumWeights)) then
     mySumWeights = SumWeights
  else
     mySumWeights = .FALSE.
  endif

       ! Compute the local contributions to the two integrals:

  if(mySumWeights) then
     call AttrVect_LocalWeightedSumRAttr(inAv1, outAv1, Weights1, WeightName1)
     call AttrVect_LocalWeightedSumRAttr(inAv2, outAv2, Weights2, WeightName2)
  else
     call AttrVect_LocalWeightedSumRAttr(inAv1, outAv1, Weights1)
     call AttrVect_LocalWeightedSumRAttr(inAv2, outAv2, Weights2)
  endif

       ! Create the paired buffer for the Global Sum

  nRA1 = AttrVect_nRAttr(outAv1)
  nRA2 = AttrVect_nRAttr(outAv2)

  PairedBufferLength =  nRA1 + nRA2
  allocate(PairedBuffer(PairedBufferLength), OutPairedBuffer(PairedBufferLength), &
           stat=ierr)
  if(ierr /= 0) then
     write(stderr,'(2a,i8)') myname_, &
	  ':: Fatal error--allocate(PairedBuffer...failed, ierr = ',ierr
     call die(myname_)
  endif

       ! Load the paired buffer

  PairedBuffer(1:nRA1) = outAv1%rAttr(1:nRA1,1)
  PairedBuffer(nRA1+1:PairedBufferLength) = outAv2%rAttr(1:nRA2,1)

       ! Perform the global sum on the paired buffer

  call MPI_AllReduce(PairedBuffer, OutPairedBuffer, PairedBufferLength, &
                        MP_Type(PairedBuffer(1)), MP_SUM, comm, ierr)
  if(ierr /= 0) then
     write(stderr,'(2a,i8)') myname_, &
  	      ':: Fatal Error--MPI_ALLREDUCE() failed with ierror = ',ierr
     call MP_perr_die(myname_,'MPI_ALLREDUCE() failed',ierr)
  endif

       ! Unload OutPairedBuffer into outAv1 and outAv2:

  outAv1%rAttr(1:nRA1,1) = OutPairedBuffer(1:nRA1)
  outAv2%rAttr(1:nRA2,1) = OutPairedBuffer(nRA1+1:PairedBufferLength)

       ! Clean up allocated arrays:

  deallocate(PairedBuffer, OutPairedBuffer, stat=ierr)
  if(ierr /= 0) then
     write(stderr,'(2a,i8)') myname_, &
	  'ERROR--deallocate(PairedBuffer,...) failed, ierr = ',ierr
     call die(myname_)
  endif

 end subroutine PairedGlobalIntegralRAttrV_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: PairedGlobalAverageRAttrGG_ - Two global spatial averages.
!
! !DESCRIPTION:
! This routine computes spatial integrals of the {\tt REAL} attributes
! of the {\tt REAL} attributes of the input {\tt AttrVect} argument 
! {\tt inAv}.  {\tt GlobalIntegralRAttrGG\_()} takes the input 
! {\tt AttrVect} argument {\tt inAv} and computes the global spatial 
! integral using weights stored in the {\tt GeneralGrid} argument 
! {\tt GGrid} and identified by the {\tt CHARACTER} tag {\tt WeightTag}.  
! The integral of each {\tt REAL} attribute is returned in the output 
! {\tt AttrVect} argument {\tt outAv}. If {\tt GlobalIntegralRAttrGG\_()} 
! is invoked with the optional {\tt LOGICAL} input argument 
! {\tt SumWeights} set as {\tt .TRUE.}, then the weights are also summed 
! and stored in {\tt outAv} (and can be referenced with the attribute 
! tag defined by the argument {\tt WeightTag}.  If 
! {\tt GlobalIntegralRAttrGG\_()} is invoked with the optional {\tt INTEGER} 
! argument {\tt comm} (a Fortran MPI communicator handle), the summation
! operations for the integral are completed on the local process, then 
! reduced across the communicator, with all processes receiving the result.
!
! {\bf N.B.:  } The local lengths of the {\tt AttrVect} argument {\tt inAv1} 
! and the {\tt GeneralGrid} {\tt GGrid1} must be equal.  That is, there 
! must be a one-to-one correspondence between the field point values stored 
! in {\tt inAv1} and the point weights stored in {\tt GGrid1}.  The same
! relationship must apply between {\tt inAv2} and {\tt GGrid2}.
!
! {\bf N.B.:  }  If {\tt GlobalIntegralRAttrGG\_()} is invoked with the 
! optional {\tt LOGICAL} input argument {\tt SumWeights} set as {\tt .TRUE.}, 
! then the value of {\tt WeightTag1} must not conflict with any of the 
! {\tt REAL} attribute tags in {\tt inAv1} and the value of {\tt WeightTag2} 
! must not conflict with any of the {\tt REAL} attribute tags in {\tt inAv2}.
!
! {\bf N.B.:  } The output {\tt AttrVect} arguments {\tt outAv1} and 
! {\tt outAv2} are allocated data structures.  The user must deallocate them
!  using the routine {\tt AttrVect\_clean()} when they are no longer needed.  
! Failure to do so will result in a memory leak.
!
! !INTERFACE:

 subroutine PairedGlobalAverageRAttrGG_(inAv1, outAv1, GGrid1, WeightTag1, &
                                        inAv2, outAv2, GGrid2, WeightTag2, &
                                        comm)
! ! USES:

      use m_stdio
      use m_die
      use m_mpif90

      use m_AttrVect, only : AttrVect
      use m_AttrVect, only : AttrVect_init => init
      use m_AttrVect, only : AttrVect_clean => clean
      use m_AttrVect, only : AttrVect_lsize => lsize
      use m_AttrVect, only : AttrVect_nRAttr => nRAttr
      use m_AttrVect, only : AttrVect_indexRA => indexRA
      use m_AttrVect, only : AttrVect_exportRListToChar => exportRListToChar

      use m_GeneralGrid, only : GeneralGrid
      use m_GeneralGrid, only : GeneralGrid_lsize => lsize
      use m_GeneralGrid, only : GeneralGrid_indexRA => indexRA
      use m_GeneralGrid, only : GeneralGrid_exportRAttr => exportRAttr

      use m_AttrVectReduce, only : AttrVect_LocalWeightedSumRAttr => &
	                                         LocalWeightedSumRAttr

      implicit none

! !INPUT PARAMETERS:

      type(AttrVect),    intent(IN) :: inAv1
      type(GeneralGrid), intent(IN) :: GGrid1
      character(len=*),  intent(IN) :: WeightTag1
      type(AttrVect),    intent(IN) :: inAv2
      type(GeneralGrid), intent(IN) :: GGrid2
      character(len=*),  intent(IN) :: WeightTag2
      integer,           intent(IN) :: comm

! !OUTPUT PARAMETERS:

      type(AttrVect),    intent(OUT) :: outAv1
      type(AttrVect),    intent(OUT) :: outAv2

! !REVISION HISTORY:
! 	09May02 - J.W. Larson <larson@mcs.anl.gov> - Initial version.
! 	14Jun02 - J.W. Larson <larson@mcs.anl.gov> - Bug fix to reflect
!                 new interface to PairedGlobalIntegralRAttrGG_().
!
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::PairedGlobalAverageRAttrGG_'

  type(AttrVect) :: integratedAv1, integratedAv2
  integer :: i, ierr, iweight1, iweight2

       ! Compute the global integral:

     call PairedGlobalIntegralRAttrGG_(inAv1, integratedAv1, GGrid1, WeightTag1, &
	                               inAv2, integratedAv2, GGrid2,     &
				       WeightTag2, .TRUE., comm)


       ! Check value of summed weights (to avoid division by zero):

  iweight1 = AttrVect_indexRA(integratedAv1, WeightTag1)
  if(integratedAv1%rAttr(iweight1, 1) == 0) then   
     write(stderr,'(2a)') myname_, &
	  '::ERROR--Global sum of grid weights in first integral is zero.'
     call die(myname_)
  endif

  iweight2 = AttrVect_indexRA(integratedAv2, WeightTag2)
  if(integratedAv2%rAttr(iweight2, 1) == 0) then
     write(stderr,'(2a)') myname_, &
	  '::ERROR--Global sum of grid weights in second integral is zero.'
     call die(myname_)
  endif

       ! Initialize output AttrVects outAv1 and outAv2:

  call AttrVect_init(outAv1, rList=AttrVect_exportRListToChar(inAv1), lsize=1)
  call AttrVect_init(outAv2, rList=AttrVect_exportRListToChar(inAv2), lsize=1)

       ! Divide by global weight sum to compute global averages from 
       ! global integrals.

  do i=1,AttrVect_nRAttr(outAv1)
     outAv1%rAttr(i,1) = integratedAv1%rAttr(i,1) &
	                               / integratedAv1%rAttr(iweight1,1) 
  end do

  do i=1,AttrVect_nRAttr(outAv2)
     outAv2%rAttr(i,1) = integratedAv2%rAttr(i,1) &
	                               / integratedAv2%rAttr(iweight2,1) 
  end do

       ! Clean up temporary AttrVects:

  call AttrVect_clean(integratedAv1)
  call AttrVect_clean(integratedAv2)

 end subroutine PairedGlobalAverageRAttrGG_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: PairedGlobalAverageRAttrV_ - Two global spatial averages.
!
! !DESCRIPTION:
! This routine computes spatial integrals of the {\tt REAL} attributes
! of the {\tt REAL} attributes of the input {\tt AttrVect} argument 
! {\tt inAv}.  {\tt GlobalIntegralRAttrGG\_()} takes the input 
! {\tt AttrVect} argument {\tt inAv} and computes the global spatial 
! integral using weights stored in the {\tt REAL} array argument 
! {\tt Weights}.  The spatial average of each {\tt REAL} attribute is 
! returned in the output {\tt AttrVect} argument {\tt outAv}. The 
! input {\tt INTEGER} argument {\tt comm} is a Fortran MPI communicator 
! handle, and the summation operations for the spatial average are completed 
! on the local process, then reduced across the communicator, with all 
! processes receiving the result.
!
! {\bf N.B.:  } The local lengths of the {\tt AttrVect} argument {\tt inAv1} 
! and the array {\tt Weights} must be equal.  That is, there must be a 
! one-to-one correspondence between the field point values stored 
! in {\tt inAv1} and the spatial weights stored in {\tt Weights}
!
! {\bf N.B.:  } The output {\tt AttrVect} arguments {\tt outAv1} and 
! {\tt outAv2} are allocated data structures.  The user must deallocate them
!  using the routine {\tt AttrVect\_clean()} when they are no longer needed.  
! Failure to do so will result in a memory leak.
!
! !INTERFACE:

 subroutine PairedGlobalAverageRAttrV_(inAv1, outAv1, Weights1, inAv2, &
                                       outAv2, Weights2, comm)
! ! USES:

      use m_stdio
      use m_die
      use m_mpif90

      use m_AttrVect, only : AttrVect
      use m_AttrVect, only : AttrVect_init => init
      use m_AttrVect, only : AttrVect_clean => clean
      use m_AttrVect, only : AttrVect_lsize => lsize
      use m_AttrVect, only : AttrVect_nRAttr => nRAttr
      use m_AttrVect, only : AttrVect_indexRA => indexRA
      use m_AttrVect, only : AttrVect_exportRListToChar => exportRListToChar

      use m_AttrVectReduce, only : AttrVect_LocalWeightedSumRAttr => &
	                                         LocalWeightedSumRAttr

      implicit none

! !INPUT PARAMETERS:

      type(AttrVect),     intent(IN) :: inAv1
      real, dimension(:), pointer    :: Weights1
      type(AttrVect),     intent(IN) :: inAv2
      real, dimension(:), pointer    :: Weights2
      integer,            intent(IN) :: comm

! !OUTPUT PARAMETERS:

      type(AttrVect),    intent(OUT) :: outAv1
      type(AttrVect),    intent(OUT) :: outAv2

! !REVISION HISTORY:
! 	09May02 - J.W. Larson <larson@mcs.anl.gov> - Initial version.
!
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::PairedGlobalAverageRAttrV_'

  type(AttrVect) :: integratedAv1, integratedAv2
  integer :: i, ierr, iweight1, iweight2

       ! weight tags used to keep track of spatial weight sums
  character*8, parameter :: WeightName1='WeightSum1'
  character*8, parameter :: WeightName2='WeightSum2'

       ! Compute the paired global integral, including spatial weights:

  call PairedGlobalIntegralRAttrV_(inAv1, integratedAv1, Weights1, WeightName1, &
                                   inAv2, integratedAv2, Weights2, WeightName2, &
                                   .TRUE., comm)

       ! Check value of summed weights (to avoid division by zero):

  iweight1 = AttrVect_indexRA(integratedAv1, WeightName1)
  if(integratedAv1%rAttr(iweight1, 1) == 0) then   
     write(stderr,'(2a)') myname_, &
	  '::ERROR--Global sum of grid weights in first integral is zero.'
     call die(myname_)
  endif

  iweight2 = AttrVect_indexRA(integratedAv2, WeightName2)
  if(integratedAv2%rAttr(iweight2, 1) == 0) then
     write(stderr,'(2a)') myname_, &
	  '::ERROR--Global sum of grid weights in second integral is zero.'
     call die(myname_)
  endif

       ! Initialize output AttrVects outAv1 and outAv2:

  call AttrVect_init(outAv1, rList=AttrVect_exportRListToChar(inAv1), lsize=1)
  call AttrVect_init(outAv2, rList=AttrVect_exportRListToChar(inAv2), lsize=1)

       ! Divide by global weight sum to compute global averages from 
       ! global integrals.

  do i=1,AttrVect_nRAttr(outAv1)
     outAv1%rAttr(i,1) = integratedAv1%rAttr(i,1) &
	                               / integratedAv1%rAttr(iweight1,1) 
  end do

  do i=1,AttrVect_nRAttr(outAv2)
     outAv2%rAttr(i,1) = integratedAv2%rAttr(i,1) &
	                               / integratedAv2%rAttr(iweight2,1) 
  end do

       ! Clean up temporary AttrVects:

  call AttrVect_clean(integratedAv1)
  call AttrVect_clean(integratedAv2)

 end subroutine PairedGlobalAverageRAttrV_

 end module m_GlobalIntegral

