!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_GlobalIntegrals - Spatial Integrals and Averaging.
!
! !DESCRIPTION:
!
! !INTERFACE:

 module m_GlobalIntegrals

      implicit none

      private	! except

! !PUBLIC MEMBER FUNCTIONS:

      public :: GlobalIntegral        ! Global Integral
      public :: GlobalAverage         ! Global Area Average

!      public :: PairedGlobalIntegrals ! A Pair of Global 
                                      ! Integrals 

!      public :: PairedGlobalAverages  ! A Pair of Global 
                                      ! Area Averages

      interface GlobalIntegral ; module procedure &
	   GlobalIntegral_
      end interface
      interface GlobalAverage ; module procedure &
	   GlobalAverage_
      end interface
!      interface PariedGlobalIntegrals ; module procedure &
!	   PairedGlobalIntegrals_
!      end interface
!      interface PariedGlobalAveragess ; module procedure &
!	   PairedGlobalAverages_
!      end interface

! !REVISION HISTORY:
! 	25Oct01 - J.W. Larson <larson@mcs.anl.gov> - Initial version
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_GlobalIntegral'

 contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: GlobalIntegral_ - Compute global spatial integral.
!
! !DESCRIPTION:
! This routine computes spatial integrals of the {\tt REAL} attributes
! of the {\tt REAL} attributes of the input {\tt AttrVect} argument 
! {\tt inAv}.  {\tt GlobalIntegral\_()} takes the input {\tt AttrVect} 
! argument {\tt inAv} and computes the global spatial integral using weights 
! stored in the {\tt GeneralGrid} argument {\tt GGrid} and identified by
! the {\tt CHARACTER} tag {\tt WeightTag}.  The integral of each {\tt REAL}
! attribute is returned in the output {\tt AttrVect} argument {\tt outAv}.
! If {\tt GlobalIntegral\_()} is invoked with the optional {\tt LOGICAL} 
! input argument {\tt SumWeights} set as {\tt .TRUE.}, then the weights are 
! also summed and stored in {\tt outAv} (and can be referenced with the 
! attribute tag defined by the argument {\tt WeightTag}.  If 
! {\tt GlobalIntegral\_()} is invoked with the optional {\tt INTEGER} 
! argument {\tt comm} (a Fortran MPI communicator handle), the summation
! operations for the integral are completed on the local process, then 
! reduced across the communicator, with all processes receiving the result.
!
! {\bf N.B.:  } The local lengths of the {\tt AttrVect} argument {\tt inAv} 
! and the {\tt GeneralGrid} {\tt GGrid} must be equal.  That is, there 
! must be a one-to-one correspondence between the field point values stored 
! in {\tt inAv} and the point weights stored in {\tt GGrid}.
!
! {\bf N.B.:  }  If {\tt GlobalIntegral\_()} is invoked with the optional 
! {\tt LOGICAL} input argument {\tt SumWeights} set as {\tt .TRUE.}, then 
! the value of {\tt WeightTag} must not conflict with any of the {\tt REAL}
! attribute tags in {\tt inAv}.
!
! {\bf N.B.:  } The output {\tt AttrVect} argument {\tt outAv} is an 
! allocated data structure.  The user must deallocate it using the routine 
! {\tt AttrVect\_clean()} when it is no longer needed.  Failure to do so 
! will result in a memory leak.
!
! !INTERFACE:

 subroutine GlobalIntegral_(inAv, outAv, GGrid, WeightTag, SumWeights, &
                            comm)
! ! USES:

      use m_stdio
      use m_die
      use m_mpif90

      use m_List, only : List
      use m_List, only : List_init => init
      use m_List, only : List_clean => clean
      use m_List, only : List_nullify => nullify
      use m_List, only : List_concatenate => concatenate

      use m_AttrVect, only : AttrVect
      use m_AttrVect, only : AttrVect_init => init
      use m_AttrVect, only : AttrVect_lsize => lsize
      use m_AttrVect, only : AttrVect_nRAttr => nRAttr
      use m_AttrVect, only : AttrVect_indexRA => indexRA

      use m_GeneralGrid, only : GeneralGrid
      use m_GeneralGrid, only : GeneralGrid_lsize => lsize
      use m_GeneralGrid, only : GeneralGrid_indexRA => indexRA
      use m_GeneralGrid, only : GeneralGrid_exportRAttr => exportRAttr

      implicit none

! !INPUT PARAMETERS:

      type(AttrVect),    intent(in) :: inAv
      type(GeneralGrid), intent(in) :: GGrid
      character(len=*),  intent(in) :: WeightTag
      logical, optional, intent(in) :: SumWeights
      integer, optional, intent(in) :: comm

! !OUTPUT PARAMETERS:

      type(AttrVect),    intent(out) :: outAv

! !REVISION HISTORY:
! 	06Feb02 - J.W. Larson <larson@mcs.anl.gov> - initial version
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::GlobalIntegral_'

  integer :: ierr, length
  logical :: mySumWeights
  integer :: NumIntegrals, i, j
  real, dimension(:), allocatable :: LocalIntegrals, GlobalIntegrals
  real, dimension(:), pointer :: gridWeights
  type(List) :: WeightList, VarList, outAvIList, outAvRList

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

       ! Extract Grid Weights

  call GeneralGrid_exportRAttr(GGrid, WeightTag, gridWeights, length)

       ! Allocate space to compute integrals

  if(mySumWeights) then
     NumIntegrals = AttrVect_nRAttr(inAv) + 1
  else
     NumIntegrals = AttrVect_nRAttr(inAv)
  endif

  allocate(LocalIntegrals(NumIntegrals), GlobalIntegrals(NumIntegrals), & 
           stat=ierr)
  if(ierr /= 0) then
     call die(myname_, "allocate(LocalIntegrals... failed", ierr)
  endif

       ! Compute the desired integrals

  LocalIntegrals(1:NumIntegrals) = 0.
  GlobalIntegrals(1:NumIntegrals) = 0.

  do i=1, AttrVect_lsize(inAv)
     do j=1, AttrVect_nRAttr(inAv)
	LocalIntegrals(j) = LocalIntegrals(j) + inAv%rAttr(j,i)
     end do
  end do

       ! Sum GeneralGrid weights (if desired)

  if(mySumWeights) then
     do i=1,length
	LocalIntegrals(NumIntegrals) = LocalIntegrals(NumIntegrals) &
	                               + gridWeights(i)
     end do
  endif

       ! If this is a distributed operation, do an ALLREDUCE:

  if(present(comm)) then
     call MPI_ALLREDUCE(LocalIntegrals, GlobalIntegrals, NumIntegrals, &
	                MP_type(LocalIntegrals(1)), MP_SUM, comm, ierr)
     if(ierr /= 0) then
	call MP_perr_die(myname_, "MPI_ALLREDUCE(...", ierr)
     endif
  else
     GlobalIntegrals(1:NumIntegrals) = LocalIntegrals(1:NumIntegrals)
  endif

       ! Create output AttrVect outAv to store integrals

  if(mySumWeights) then ! include the Weights as an attribute

     call List_init(WeightList, WeightTag)
     call AttrVect_exportRList(inAv, VarList) ! creates VarList
     call List_concatenate(VarList, WeightList, outAvRList)
     call List_clean(WeightList)
     call List_clean(VarList)

  else ! we are returning only the integrals of the input variables

     call AttrVect_exportRList(inAv, outAvRList)

  endif

       ! Nullify the argument outAvIList (this is a work-around to
       ! cope with the fact the list-based AttrVect_init requires 
       ! both real and integer attribute lists.

  call List_nullify(outAvIList)

  call AttrVect_init(outAv, outAvIList, outAvRList, lsize=1)

  call List_clean(outAvRList)  ! outAvIList is already cleaned.

       ! Clean up allocated structures

  deallocate(LocalIntegrals, GlobalIntegrals, gridWeights, stat=ierr)
  if(ierr /= 0) then
     call die(myname_, "deallocate(LocalIntegrals...", ierr)
  endif

 end subroutine GlobalIntegral_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: GlobalAverage_ - Compute global spatial average.
!
! !DESCRIPTION:
! This routine computes spatial averages of the {\tt REAL} attributes
! of the input {\tt AttrVect} argument {\tt inAv}.  
! {\tt GlobalAverage\_()} takes the input {\tt AttrVect} argument 
! {\tt inAv} and computes the global spatial average using weights 
! stored in the {\tt GeneralGrid} argument {\tt GGrid} and identified by
! the {\tt CHARACTER} tag {\tt WeightTag}.  The average of each {\tt REAL}
! attribute is returned in the output {\tt AttrVect} argument {\tt outAv}.
! If {\tt GlobalAverage\_()} is invoked with the optional {\tt INTEGER} 
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

 subroutine GlobalAverage_(inAv, outAv, GGrid, WeightTag, comm)

! ! USES:

      use m_stdio
      use m_die
      use m_mpif90

      use m_AttrVect, only : AttrVect
      use m_AttrVect, only : AttrVect_init => init
      use m_AttrVect, only : AttrVect_nRAttr => nRAttr
      use m_AttrVect, only : AttrVect_indexRA => indexRA

      use m_GeneralGrid, only : GeneralGrid

      implicit none

! !INPUT PARAMETERS:

      type(AttrVect),    intent(in) :: inAv
      type(GeneralGrid), intent(in) :: GGrid
      character(len=*),  intent(in) :: WeightTag
      integer, optional, intent(in) :: comm

! !OUTPUT PARAMETERS:

      type(AttrVect),    intent(out) :: outAv

! !REVISION HISTORY:
! 	08Feb02 - J.W. Larson <larson@mcs.anl.gov> - initial version
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::GlobalAverage_'

  type(AttrVect) :: integratedAv
  integer :: i, ierr, iweight

       ! Compute the global integral:

  if(present(comm)) then
     call GlobalIntegral_(inAv, integratedAv, GGrid, WeightTag, &
	                  SumWeights=.TRUE., comm=comm)
  else
     call GlobalIntegral_(inAv, integratedAv, GGrid, WeightTag, &
	                  SumWeights=.TRUE.)
  endif

       ! Check value of summed weights (to avoid division by zero):

  iweight = AttrVect_indexRA(integratedAv, WeightTag)
  if(integratedAv%rAttr(iweight, 1) == 0) then
     write(stderr,'(2a)') myname_, &
	  '::ERROR--Global sum of grid weights is zero.'
     call die(myname_)
  endif

       ! Initialize output AttrVect outAv:

  call AttrVect_init(outAv, inAv, lsize=1)

       ! Divide by global weight sum to compute global averages from 
       ! global integrals.

  do i=1,AttrVect_nRAttr(outAv)
     outAv%rAttr(i,1) = integratedAv%rAttr(i,1) / integratedAv%rAttr(iweight,1) 
  end do

       ! Clean up temporary AttrVect:

  call AttrVect_clean(integratedAv)

 end subroutine GlobalAverage_

 end module m_GlobalIntegrals

