!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_SpatialIntegral - Spatial Integrals and Averages
!
! !DESCRIPTION:  This module provides spatial integration and averaging 
! services for the MCT.  For a field $\Phi$ sampled at a point ${\bf x}$ 
! in some multidimensional domain $\Omega$, the integral $I$ of 
! $\Phi({\bf x})$ is
! $$ I = \int_{\Omega} \Phi ({\bf x}) d\Omega .$$ 
! The spatial average $A$ of $\Phi({\bf x})$ over $\Omega$ is
! $$ A = {{ \int_{\Omega} \Phi ({\bf x}) d\Omega} \over 
! { \int_{\Omega} d\Omega} }. $$
! Since the {\tt AttrVect} represents a discretized field, the integrals 
! above are implemented as:
! $$ I = \sum_{i=1}^N \Phi_i \Delta \Omega_i $$
! and 
! $$ A = {{\sum_{i=1}^N \Phi_i \Delta \Omega_i } \over 
!{\sum_{i=1}^N \Delta \Omega_i } }, $$
! where $N$ is the number of physical locations, $\Phi_i$ is the value 
! of the field $\Phi$ at location $i$, and $\Delta \Omega_i$ is the spatial 
! weight (lenghth element, cross-sectional area element, volume element, 
! {\em et cetera}) at location $i$.
!
! MCT extends the concept of integrals and area/volume averages to include 
! {\em masked} integrals and averages.  MCT recognizes both {\em integer}
! and {\em real} masks.  An integer mask $M$ is a vector of integers (one
! corresponding to each physical location) with each element having value 
! either zero or one.  Integer masks are used to include/exclude data from
! averages or integrals.  For example, if one were to compute globally 
! averaged cloud amount over land (but not ocean nor sea-ice), one would 
! assign a $1$ to each location on the land and a $0$ to each non-land 
! location.  A {\em real} mask $F$ is a vector of real numbers (one corresponding 
! to each physical location) with each element having value within the 
! closed interval $[0,1]$.  .Real masks are used to represent fractional 
! area/volume coverage at a location by a given component model.  For 
! example, if one wishes to compute area averages over sea-ice, one must 
! include the ice fraction present at each point.  Masked Integrals and 
! averages are represented in the MCT by:
! $$ I = \sum_{i=1}^N {\prod_{j=1}^J M_i} {\prod_{k=1}^K F_i} 
! \Phi_i \Delta \Omega_i $$
! and 
! $$ A = {{\sum_{i=1}^N \bigg({\prod_{j=1}^J M_i}\bigg) \bigg( {\prod_{k=1}^K F_i}
! \bigg) \Phi_i 
! \Delta \Omega_i } \over 
!{\sum_{i=1}^N \bigg({\prod_{j=1}^J M_i}\bigg) \bigg( {\prod_{k=1}^K F_i} \bigg) 
!  \Delta \Omega_i } }, $$
! where $J$ is the number of integer masks and $K$ is the number of real masks.
!
! All of the routines in this module assume field data is stored in an 
! attribute vector ({\tt AttrVect}), and the integration/averaging is performed 
! only on the {\tt REAL} attributes.  Physical coordinate grid and mask 
! information is assumed to be stored as attributes in either a 
! {\tt GeneralGrid}, or pre-combined into a single integer mask and a single 
! real mask.  
!
! !INTERFACE:

 module m_SpatialIntegral

      implicit none

      private	! except

! !PUBLIC MEMBER FUNCTIONS:

      public :: SpatialIntegral        ! Spatial Integral
      public :: SpatialAverage         ! Spatial Area Average

      public :: MaskedSpatialIntegral  ! Masked Spatial Integral
      public :: MaskedSpatialAverage   ! MaskedSpatial Area Average

      public :: PairedSpatialIntegrals ! A Pair of Spatial 
                                      ! Integrals 

      public :: PairedSpatialAverages  ! A Pair of Spatial 
                                      ! Area Averages

      public :: PairedMaskedSpatialIntegrals ! A Pair of Masked
                                             ! Spatial Integrals 

      public :: PairedMaskedSpatialAverages ! A Pair of Masked
                                            ! Spatial Area Averages

      interface SpatialIntegral ; module procedure &
	   SpatialIntegralRAttrGG_, &
	   SpatialIntegralRAttrV_
      end interface
      interface SpatialAverage ; module procedure &
	   SpatialAverageRAttrGG_, &
	   SpatialAverageRAttrV_
      end interface
      interface MaskedSpatialIntegral ; module procedure &
	   MaskedSpatialIntegralRAttrV_, &
	   MaskedSpatialIntegralRAttrGG_
      end interface
      interface MaskedSpatialAverage ; module procedure &
	   MaskedSpatialAverageRAttrV_, &
	   MaskedSpatialAverageRAttrGG_
      end interface
      interface PairedSpatialIntegrals ; module procedure &
	    PairedSpatialIntegralRAttrGG_, &
	    PairedSpatialIntegralRAttrV_
      end interface
      interface PairedSpatialAverages ; module procedure &
	    PairedSpatialAverageRAttrGG_, &
	    PairedSpatialAverageRAttrV_
      end interface
      interface PairedMaskedSpatialIntegrals ; module procedure &
	    PairedMaskedIntegralRAttrGG_
      end interface
      interface PairedMaskedSpatialAverages ; module procedure &
	    PairedMaskedAverageRAttrGG_
      end interface

! !REVISION HISTORY:
! 	25Oct01 - J.W. Larson <larson@mcs.anl.gov> - Initial version
!        9May02 - J.W. Larson <larson@mcs.anl.gov> - Massive Refactoring.
!    10-14Jun02 - J.W. Larson <larson@mcs.anl.gov> - Added Masked methods.
!    17-18Jun02 - J.W. Larson <larson@mcs.anl.gov> - Added Paired/Masked 
!                 methods.
!       18Jun02 - J.W. Larson <larson@mcs.anl.gov> - Renamed module from 
!                 m_GlobalIntegral to m_SpatialIntegral.
!
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_SpatialIntegral'

 contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: SpatialIntegralRAttrGG_ - Compute spatial integral.
!
! !DESCRIPTION:
! This routine computes spatial integrals of the {\tt REAL} attributes
! of the {\tt REAL} attributes of the input {\tt AttrVect} argument 
! {\tt inAv}.  {\tt SpatialIntegralRAttrGG\_()} takes the input 
! {\tt AttrVect} argument {\tt inAv} and computes the spatial 
! integral using weights stored in the {\tt GeneralGrid} argument 
! {\tt GGrid} and identified by the {\tt CHARACTER} tag {\tt WeightTag}.  
! The integral of each {\tt REAL} attribute is returned in the output 
! {\tt AttrVect} argument {\tt outAv}. If {\tt SpatialIntegralRAttrGG\_()} 
! is invoked with the optional {\tt LOGICAL} input argument 
! {\tt SumWeights} set as {\tt .TRUE.}, then the weights are also summed 
! and stored in {\tt outAv} (and can be referenced with the attribute 
! tag defined by the argument{\tt WeightTag}.  If 
! {\tt SpatialIntegralRAttrGG\_()} is invoked with the optional {\tt INTEGER} 
! argument {\tt comm} (a Fortran MPI communicator handle), the summation
! operations for the integral are completed on the local process, then 
! reduced across the communicator, with all processes receiving the result.
!
! {\bf N.B.:  } The local lengths of the {\tt AttrVect} argument {\tt inAv} 
! and the {\tt GeneralGrid} {\tt GGrid} must be equal.  That is, there 
! must be a one-to-one correspondence between the field point values stored 
! in {\tt inAv} and the point weights stored in {\tt GGrid}.
!
! {\bf N.B.:  }  If {\tt SpatialIntegralRAttrGG\_()} is invoked with the 
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

 subroutine SpatialIntegralRAttrGG_(inAv, outAv, GGrid, WeightTag, &
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
!                 renamed SpatialIntegralRAttrGG_().
!       07Jun02 - J.W. Larson <larson@mcs.anl.gov> - Bug fix and further
!                 refactoring.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::SpatialIntegralRAttrGG_'

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
     call SpatialIntegralRAttrV_(inAv, outAv, gridWeights, mySumWeights, &
	                        WeightTag, comm)
  else
     call SpatialIntegralRAttrV_(inAv, outAv, gridWeights, mySumWeights, &
                                WeightTag)
  endif

       ! Clean up temporary allocated space

  deallocate(gridWeights, stat=ierr)
  if(ierr /= 0) then
     write(stderr,'(2a,i8)') myname_, &
	              ':: deallocate(gridWeights...failed.  ierr=', ierr
     call die(myname_)
  endif

 end subroutine SpatialIntegralRAttrGG_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: SpatialIntegralRAttrV_ - Compute spatial integral.
!
! !DESCRIPTION:
! This routine computes spatial integrals of the {\tt REAL} attributes
! of the {\tt REAL} attributes of the input {\tt AttrVect} argument 
! {\tt inAv}.  {\tt SpatialIntegralRAttrV\_()} takes the input 
! {\tt AttrVect} argument {\tt inAv} and computes the spatial 
! integral using weights stored in the input {\tt REAL} array argument 
! {\tt Weights}.  The integral of each {\tt REAL} attribute is returned 
! in the output {\tt AttrVect} argument {\tt outAv}. If 
! {\tt SpatialIntegralRAttrV\_()} is invoked with the optional {\tt LOGICAL} 
! input argument {\tt SumWeights} set as {\tt .TRUE.}, then the weights 
! are also summed and stored in {\tt outAv} (and can be referenced with 
! the attribute name {\tt WeightTag}.  If {\tt SpatialIntegralRAttrV\_()} is 
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
! {\bf N.B.:  }  If {\tt SpatialIntegralRAttrV\_()} is invoked with the 
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

 subroutine SpatialIntegralRAttrV_(inAv, outAv, Weights, SumWeights, &
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

  character(len=*),parameter :: myname_=myname//'::SpatialIntegralRAttrV_'

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

     if(mySumWeights) then ! return the spatial sum of the weights in outAV
	call AttrVect_GlobalWeightedSumRAttr(inAV, outAV, Weights, &
 	                                     comm, WeightTag)
     else
	call AttrVect_GlobalWeightedSumRAttr(inAV, outAV, Weights, comm)
     endif

  else ! compute local sum:

     if(mySumWeights) then ! return the spatial sum of the weights in outAV
	call AttrVect_LocalWeightedSumRAttr(inAV, outAV, Weights, &
                                            WeightTag)
     else
	call AttrVect_LocalWeightedSumRAttr(inAV, outAV, Weights)
     endif

  endif ! if(present(comm))...

 end subroutine SpatialIntegralRAttrV_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: SpatialAverageRAttrGG_ - Compute spatial average.
!
! !DESCRIPTION:
! This routine computes spatial averages of the {\tt REAL} attributes
! of the input {\tt AttrVect} argument {\tt inAv}.  
! {\tt SpatialAverageRAttrGG\_()} takes the input {\tt AttrVect} argument 
! {\tt inAv} and computes the spatial average using weights 
! stored in the {\tt GeneralGrid} argument {\tt GGrid} and identified by
! the {\tt CHARACTER} tag {\tt WeightTag}.  The average of each {\tt REAL}
! attribute is returned in the output {\tt AttrVect} argument {\tt outAv}.
! If {\tt SpatialAverageRAttrGG\_()} is invoked with the optional {\tt INTEGER} 
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

 subroutine SpatialAverageRAttrGG_(inAv, outAv, GGrid, WeightTag, comm)

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
!       18Jun02 - J.W. Larson <larson@mcs.anl.gov> - Renamed routine to
!                 SpatialAverageRAttrGG_().
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::SpatialAverageRAtttrGG_'

  type(AttrVect) :: integratedAv
  integer :: i, ierr, iweight

       ! Compute the spatial integral:

  if(present(comm)) then
     call SpatialIntegralRAttrGG_(inAv, integratedAv, GGrid, WeightTag, &
	                         .TRUE., comm)
  else
     call SpatialIntegralRAttrGG_(inAv, integratedAv, GGrid, WeightTag, &
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

       ! Divide by global weight sum to compute spatial averages from 
       ! spatial integrals.

  do i=1,AttrVect_nRAttr(outAv)
     outAv%rAttr(i,1) = integratedAv%rAttr(i,1) &
	                               / integratedAv%rAttr(iweight,1) 
  end do

       ! Clean up temporary AttrVect:

  call AttrVect_clean(integratedAv)

 end subroutine SpatialAverageRAttrGG_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: SpatialAverageRAttrV_ - Compute spatial average.
!
! !DESCRIPTION:
! This routine computes spatial averages of the {\tt REAL} attributes
! of the input {\tt AttrVect} argument {\tt inAv}.  
! {\tt SpatialAverageRAttrV\_()} takes the input {\tt AttrVect} argument 
! {\tt inAv} and computes the spatial average using weights 
! stored in the {\tt REAL} array {\tt Weights}.  The average of each 
! {\tt REAL} attribute is returned in the output {\tt AttrVect} argument 
! {\tt outAv}.  If {\tt SpatialAverageRAttrGG\_()} is invoked with the 
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

 subroutine SpatialAverageRAttrV_(inAv, outAv, Weights, comm)

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

  character(len=*),parameter :: myname_=myname//'::SpatialAverageRAtttrV_'

  type(AttrVect) :: integratedAv
  integer :: i, ierr, iweight

       ! Nullify all the pointers present in the AttrVect integratedAv:

  call AttrVect_clean(integratedAv)

       ! Compute the spatial integral:

  if(present(comm)) then
     call SpatialIntegralRAttrV_(inAv, integratedAv, Weights, &
	                         .TRUE., 'weights', comm)
  else
     call SpatialIntegralRAttrV_(inAv, integratedAv, Weights, .TRUE., 'weights')
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

       ! Divide by global weight sum to compute spatial averages from 
       ! spatial integrals.

  do i=1,AttrVect_nRAttr(outAv)
     outAv%rAttr(i,1) = integratedAv%rAttr(i,1) &
	                               / integratedAv%rAttr(iweight,1) 
  end do

       ! Clean up temporary AttrVect:

  call AttrVect_clean(integratedAv)

 end subroutine SpatialAverageRAttrV_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: MaskedSpatialIntegralRAttrGG_ - Masked spatial integral.
!
! !DESCRIPTION: 
! This routine computes masked spatial integrals of the {\tt REAL} 
! attributes of the input {\tt AttrVect} argument {\tt inAv}, returning 
! the masked integrals in the output {\tt AttrVect} {\tt outAv}.  All of 
! the masking data are assumed stored in the input {\tt GeneralGrid} 
! argument {\tt GGrid}.  If integer masks are to be used, their integer 
! attribute names in {\tt GGrid} are named as a colon-delimited list 
! in the optional {\tt CHARACTER} input argument {\tt iMaskTags}.  Real 
! masks (if desired) are referenced by their real attribute names in 
! {\tt GGrid} are named as a colon-delimited list in the optional 
! {\tt CHARACTER} input argument {\tt rMaskTags}.  The user specifies 
! a choice of mask combination method with the input {\tt LOGICAL} argument
! {\tt UseFastMethod}.  If ${\tt UseFastMethod} = {\tt .FALSE.}$ this 
! routine checks each mask entry to ensure that the integer masks contain
! only ones and zeroes, and that entries in the real masks are all in
! the closed interval $[0,1]$.  If ${\tt UseFastMethod} = {\tt .TRUE.}$,
! this routine performs direct products of the masks, assuming that the
! user has validated them in advance.  The optional {\tt LOGICAL} input 
! argument {\tt SumWeights} determines whether the masked sum of the spatial
! weights is computed and returned in {\tt outAv} with the real attribute 
! name supplied in the optional {\tt CHARACTER} input argument 
! {\tt WeightSumTag}.  This integral can either be a local (i.e. a global 
! memory space operation), or a global distributed integral.  The latter 
! is the case if the optional input {\tt INTEGER} argument {\tt comm} is
! supplied (which corresponds to a Fortran MPI communicatior handle).
!
! {\bf N.B.:  } The local lengths of the {\tt AttrVect} argument {\tt inAv} 
! and the input {\tt GeneralGrid} {\tt GGrid} must be equal.  That is, there 
! must be a one-to-one correspondence between the field point values stored 
! in {\tt inAv} and the point weights stored in {\tt GGrid}.
!
! {\bf N.B.:  }  If {\tt SpatialIntegralRAttrV\_()} is invoked with the 
! optional {\tt LOGICAL} input argument {\tt SumWeights} set as {\tt .TRUE.}. 
! In this case, the none of {\tt REAL} attribute tags in {\tt inAv} may be 
! named the same as the string contained in {\tt WeightSumTag}, which is an 
! attribute name reserved for the sum of the weights in the output {\tt AttrVect} 
! {\tt outAv}.
!
! {\bf N.B.:  } The output {\tt AttrVect} argument {\tt outAv} is an 
! allocated data structure.  The user must deallocate it using the routine 
! {\tt AttrVect\_clean()} when it is no longer needed.  Failure to do so 
! will result in a memory leak.
!
! !INTERFACE:

 subroutine MaskedSpatialIntegralRAttrGG_(inAv, outAv, GGrid, SpatialWeightTag, &
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

  character(len=*),parameter :: myname_=myname//'::MaskedSpatialIntegralRAttrGG_'

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
	      call MaskedSpatialIntegralRAttrV_(inAv, outAv, SpatialWeights, &
		                               iMask, rMask, UseFastMethod, &
					       SumWeights, WeightSumTag, comm)
	   else ! Do not return the masked sum of the weights
	      call MaskedSpatialIntegralRAttrV_(inAv, outAv, SpatialWeights, &
		                               iMask, rMask, UseFastMethod, &
					       comm=comm)
	   endif ! if(mySumWeights)...

	else ! compute local sum:

	   if(mySumWeights) then ! return the global masked sum of the 
                                 ! weights in outAV
	      call MaskedSpatialIntegralRAttrV_(inAv, outAv, SpatialWeights, &
		                               iMask, rMask, UseFastMethod, &
					       SumWeights, WeightSumTag)
	   else ! Do not return the masked sum of the weights
	      call MaskedSpatialIntegralRAttrV_(inAv, outAv, SpatialWeights, &
		                               iMask, rMask, UseFastMethod)
	   endif ! if(mySumWeights)...

	endif ! if(present(comm))...

     else ! REAL Mask Only Case...

	if(present(comm)) then ! compute distributed AllReduce-style sum:
	   
	   if(mySumWeights) then ! return the global masked sum of the 
                                 ! weights in outAV
	      call MaskedSpatialIntegralRAttrV_(inAv, outAv, SpatialWeights, &
		                               rMask=rMask, &
					       UseFastMethod=UseFastMethod, &
					       SumWeights=SumWeights, &
					       WeightSumTag=WeightSumTag, &
					       comm=comm)
	   else ! Do not return the masked sum of the weights
	      call MaskedSpatialIntegralRAttrV_(inAv, outAv, SpatialWeights, &
		                               rMask=rMask, &
					       UseFastMethod=UseFastMethod, &
					       comm=comm)
	   endif ! if(mySumWeights)...

	else ! compute local sum:

	   if(mySumWeights) then ! return the global masked sum of the 
                                 ! weights in outAV
	      call MaskedSpatialIntegralRAttrV_(inAv, outAv, SpatialWeights, &
		                               rMask=rMask, &
					       UseFastMethod=UseFastMethod, &
					       SumWeights=SumWeights, &
					       WeightSumTag=WeightSumTag)
	   else ! Do not return the masked sum of the weights
	      call MaskedSpatialIntegralRAttrV_(inAv, outAv, SpatialWeights, &
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
	      call MaskedSpatialIntegralRAttrV_(inAv, outAv, SpatialWeights, &
		                               iMask=iMask, &
					       UseFastMethod=UseFastMethod, &
					       SumWeights=SumWeights, &
					       WeightSumTag=WeightSumTag, &
					       comm=comm)
	   else ! Do not return the masked sum of the weights
	      call MaskedSpatialIntegralRAttrV_(inAv, outAv, SpatialWeights, &
		                               iMask=iMask, &
					       UseFastMethod=UseFastMethod, &
					       comm=comm)
	   endif ! if(mySumWeights)...

	else ! compute local sum:

	   if(mySumWeights) then ! return the global masked sum of the 
                                 ! weights in outAV
	      call MaskedSpatialIntegralRAttrV_(inAv, outAv, SpatialWeights, &
		                               iMask=iMask, &
					       UseFastMethod=UseFastMethod, &
					       SumWeights=SumWeights, &
					       WeightSumTag=WeightSumTag)
	   else ! Do not return the masked sum of the weights
	      call MaskedSpatialIntegralRAttrV_(inAv, outAv, SpatialWeights, &
		                               iMask=iMask, &
					       UseFastMethod=UseFastMethod)
	   endif ! if(mySumWeights)...

	endif ! if(present(comm))...

     else ! no INTEGER Mask / no REAL Mask Case...

	if(present(comm)) then ! compute distributed AllReduce-style sum:
	   
	   if(mySumWeights) then ! return the global masked sum of the 
                                 ! weights in outAV
	      call MaskedSpatialIntegralRAttrV_(inAv, outAv, SpatialWeights, &
					       UseFastMethod=UseFastMethod, &
					       SumWeights=SumWeights, &
					       WeightSumTag=WeightSumTag, &
					       comm=comm)
	   else ! Do not return the masked sum of the weights
	      call MaskedSpatialIntegralRAttrV_(inAv, outAv, SpatialWeights, &
					       UseFastMethod=UseFastMethod, &
					       comm=comm)
	   endif ! if(mySumWeights)...

	else ! compute local sum:

	   if(mySumWeights) then ! return the global masked sum of the 
                                 ! weights in outAV
	      call MaskedSpatialIntegralRAttrV_(inAv, outAv, SpatialWeights, &
					       UseFastMethod=UseFastMethod, &
					       SumWeights=SumWeights, &
					       WeightSumTag=WeightSumTag)
	   else ! Do not return the masked sum of the weights
	      call MaskedSpatialIntegralRAttrV_(inAv, outAv, SpatialWeights, &
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

 end subroutine MaskedSpatialIntegralRAttrGG_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: MaskedSpatialIntegralRAttrV_ - Masked spatial integral.
!
! !DESCRIPTION: 
! This routine computes masked spatial integrals of the {\tt REAL} 
! attributes of the input {\tt AttrVect} argument {\tt inAv}, returning 
! the masked integrals in the output {\tt AttrVect} argument {\tt outAv}.  
! The masked integral is computed using weights stored in the input 
! {\tt REAL} array argument {\tt SpatialWeights}.  Integer masking (if 
! desired) is provided in the optional input {\tt INTEGER} array {\tt iMask},
! and real masking (if desired) is provided in the optional input {\tt REAL} 
! array {\tt rMask}.  If {\tt SpatialIntegralRAttrV\_()} is invoked with the 
! optional {\tt LOGICAL} input argument {\tt SumWeights} set as {\tt .TRUE.}, 
! then the weights are also summed and stored in {\tt outAv} (and can be 
! referenced with the attribute name defined by the optional input 
! {\tt CHARACTER} argument {\tt WeightSumTag}.  If 
! {\tt SpatialIntegralRAttrV\_()} is invoked with the optional {\tt INTEGER} 
! argument {\tt comm} (a Fortran MPI communicator handle), the summation 
! operations for the integral are completed on the local process, then 
! reduced across the communicator, with all processes receiving the result.  
! Otherwise, the integral is assumed to be local (or equivalent to a global 
! address space).
!
! {\bf N.B.:  } The local lengths of the {\tt AttrVect} argument {\tt inAv} 
! and the input array {\tt Weights} must be equal.  That is, there must be 
! a one-to-one correspondence between the field point values stored 
! in {\tt inAv} and the point weights stored in {\tt SpatialWeights}.
!
! {\bf N.B.:  }  If {\tt SpatialIntegralRAttrV\_()} is invoked with the 
! optional {\tt LOGICAL} input argument {\tt SumWeights} set as {\tt .TRUE.}. 
! In this case, the none of {\tt REAL} attribute tags in {\tt inAv} may be 
! named the same as the string contained in {\tt WeightSumTag}, which is an 
! attribute name reserved for the sum of the weights in the output {\tt AttrVect} 
! {\tt outAv}.
!
! {\bf N.B.:  } The output {\tt AttrVect} argument {\tt outAv} is an 
! allocated data structure.  The user must deallocate it using the routine 
! {\tt AttrVect\_clean()} when it is no longer needed.  Failure to do so 
! will result in a memory leak.
!
! !INTERFACE:

 subroutine MaskedSpatialIntegralRAttrV_(inAv, outAv, SpatialWeights, iMask, &
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

  character(len=*),parameter :: myname_=myname//'::MaskedSpatialIntegralRAttrV_'

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

 end subroutine MaskedSpatialIntegralRAttrV_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: MaskedSpatialAverageRAttrGG_ - Masked spatial average.
!
! !DESCRIPTION: 
! This routine computes masked spatial averages of the {\tt REAL} 
! attributes of the input {\tt AttrVect} argument {\tt inAv}, returning 
! the masked averages in the output {\tt AttrVect} {\tt outAv}.  All of 
! the masking data are assumed stored in the input {\tt GeneralGrid} 
! argument {\tt GGrid}.  If integer masks are to be used, their integer 
! attribute names in {\tt GGrid} are named as a colon-delimited list 
! in the optional {\tt CHARACTER} input argument {\tt iMaskTags}.  Real 
! masks (if desired) are referenced by their real attribute names in 
! {\tt GGrid} are named as a colon-delimited list in the optional 
! {\tt CHARACTER} input argument {\tt rMaskTags}.  The user specifies 
! a choice of mask combination method with the input {\tt LOGICAL} argument
! {\tt UseFastMethod}.  If ${\tt UseFastMethod} = {\tt .FALSE.}$ this 
! routine checks each mask entry to ensure that the integer masks contain
! only ones and zeroes, and that entries in the real masks are all in
! the closed interval $[0,1]$.  If ${\tt UseFastMethod} = {\tt .TRUE.}$,
! this routine performs direct products of the masks, assuming that the
! user has validated them in advance.  This averaging can either be a 
! local (equivalent to a global memory space operation), or a global 
! distributed integral.  The latter is the case if the optional input 
! {\tt INTEGER} argument {\tt comm} is supplied (which corresponds to a 
! Fortran MPI communicatior handle).
!
! {\bf N.B.:  } The local lengths of the {\tt AttrVect} argument {\tt inAv} 
! and the input {\tt GeneralGrid} {\tt GGrid} must be equal.  That is, 
! there must be a one-to-one correspondence between the field point values 
! stored in {\tt inAv} and the point weights stored in {\tt GGrid}.
!
! {\bf N.B.:  } The output {\tt AttrVect} argument {\tt outAv} is an 
! allocated data structure.  The user must deallocate it using the routine 
! {\tt AttrVect\_clean()} when it is no longer needed.  Failure to do so 
! will result in a memory leak.
!
! !INTERFACE:

 subroutine MaskedSpatialAverageRAttrGG_(inAv, outAv, GGrid, SpatialWeightTag, &
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

  character(len=*),parameter :: myname_=myname//'::MaskedSpatialAverageRAttrGG_'

  type(AttrVect) :: integratedAv
  character*9, parameter :: WeightSumTag = 'WeightSum'

  integer :: i, iweight

       !================================================================
       ! Do the integration using MaskedSpatialIntegralRAttrGG_(), which
       ! returns the intermediate integrals (including the masked weight
       ! sum) in the AttrVect integratedAv.
       !================================================================

  if(present(iMaskTags)) then

     if(present(rMaskTags)) then ! have both iMasks and rMasks

	if(present(comm)) then ! a distributed parallel sum
	   call MaskedSpatialIntegralRAttrGG_(inAv, integratedAv, GGrid, &
		                             SpatialWeightTag, iMaskTags, &
					     rMaskTags, UseFastMethod,  &
					     .TRUE., WeightSumTag, comm)
	else ! a purely local sum
	   call MaskedSpatialIntegralRAttrGG_(inAv, integratedAv, GGrid, &
		                             SpatialWeightTag, iMaskTags, &
					     rMaskTags, UseFastMethod,  &
					     .TRUE., WeightSumTag)
	endif ! if(present(comm))...

     else ! Only iMasks are in use

	if(present(comm)) then ! a distributed parallel sum
	   call MaskedSpatialIntegralRAttrGG_(inAv, integratedAv, GGrid, &
		                             SpatialWeightTag, iMaskTags, &
					     UseFastMethod=UseFastMethod,  &
					     SumWeights=.TRUE., &
					     WeightSumTag=WeightSumTag, &
					     comm=comm)

	else ! a purely local sum
	   call MaskedSpatialIntegralRAttrGG_(inAv, integratedAv, GGrid, &
		                             SpatialWeightTag, iMaskTags, &
					     UseFastMethod=UseFastMethod,  &
					     SumWeights=.TRUE., &
					     WeightSumTag=WeightSumTag)
	endif ! if(present(comm))...

     endif ! if(present(rMaskTags)...

  else ! no iMasks

     if(present(rMaskTags)) then ! Only rMasks are in use

	if(present(comm)) then ! a distributed parallel sum
	   call MaskedSpatialIntegralRAttrGG_(inAv, integratedAv, GGrid, &
		                             SpatialWeightTag, &
					     rMaskTags=rMaskTags, &
					     UseFastMethod=UseFastMethod,  &
					     SumWeights=.TRUE., &
					     WeightSumTag=WeightSumTag, &
					     comm=comm)
	else ! a purely local sum
	   call MaskedSpatialIntegralRAttrGG_(inAv, integratedAv, GGrid, &
		                             SpatialWeightTag, &
					     rMaskTags=rMaskTags, &
					     UseFastMethod=UseFastMethod,  &
					     SumWeights=.TRUE., &
					     WeightSumTag=WeightSumTag)
	endif

     else ! Neither iMasks nor rMasks are in use

	if(present(comm)) then ! a distributed parallel sum
	   call MaskedSpatialIntegralRAttrGG_(inAv, integratedAv, GGrid, &
		                             SpatialWeightTag, &
					     UseFastMethod=UseFastMethod,  &
					     SumWeights=.TRUE., &
					     WeightSumTag=WeightSumTag, &
					     comm=comm)
	else ! a purely local sum
	   call MaskedSpatialIntegralRAttrGG_(inAv, integratedAv, GGrid, &
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

       ! Divide by global weight sum to compute spatial averages from 
       ! spatial integrals.

  do i=1,AttrVect_nRAttr(outAv)
     outAv%rAttr(i,1) = integratedAv%rAttr(i,1) &
	                               / integratedAv%rAttr(iweight,1) 
  end do

       ! Clean up temporary AttrVect:

  call AttrVect_clean(integratedAv)

 end subroutine MaskedSpatialAverageRAttrGG_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: MaskedSpatialAverageRAttrV_ - Masked spatial average.
!
! !DESCRIPTION: [NEEDS **LOTS** of work...]
! This routine computes spatial integrals of the {\tt REAL} attributes
! of the {\tt REAL} attributes of the input {\tt AttrVect} argument 
! {\tt inAv}.  {\tt SpatialIntegralRAttrV\_()} takes the input 
! {\tt AttrVect} argument {\tt inAv} and computes the spatial 
! integral using weights stored in the input {\tt REAL} array argument 
! {\tt Weights}.  The integral of each {\tt REAL} attribute is returned 
! in the output {\tt AttrVect} argument {\tt outAv}. If 
! {\tt SpatialIntegralRAttrV\_()} is invoked with the optional {\tt LOGICAL} 
! input argument {\tt SumWeights} set as {\tt .TRUE.}, then the weights 
! are also summed and stored in {\tt outAv} (and can be referenced with 
! the attribute name {\tt WeightTag}.  If {\tt SpatialIntegralRAttrV\_()} is 
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
! {\bf N.B.:  }  If {\tt SpatialIntegralRAttrV\_()} is invoked with the 
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

 subroutine MaskedSpatialAverageRAttrV_(inAv, outAv, SpatialWeights, iMask, &
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

  character(len=*),parameter :: myname_=myname//'::MaskedSpatialAverageRAttrV_'

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
	   call MaskedSpatialIntegralRAttrV_(inAv, integratedAv, SpatialWeights,  &
		                            iMask, rMask, UseFastMethod, .TRUE., &
					    'MaskedWeightsSum', comm)
	else ! no rMask
	   call MaskedSpatialIntegralRAttrV_(inAv, integratedAv, SpatialWeights,  &
		                            iMask=iMask, UseFastMethod=UseFastMethod, &
                                            SumWeights=.TRUE., &
					    WeightSumTag='MaskedWeightsSum', &
					    comm=comm)
	endif ! if(present(rMask))...

     else ! no iMask present...

	if(present(rMask)) then
	   call MaskedSpatialIntegralRAttrV_(inAv, integratedAv, SpatialWeights,  &
		                            rMask=rMask, UseFastMethod=UseFastMethod, &
					    SumWeights=.TRUE., &
					    WeightSumTag='MaskedWeightsSum', &
					    comm=comm)
	else ! neither rMask nor iMask present:
	   call MaskedSpatialIntegralRAttrV_(inAv, integratedAv, SpatialWeights,  &
		                            UseFastMethod=UseFastMethod, &
					    SumWeights=.TRUE., &
					    WeightSumTag='MaskedWeightsSum', &
					    comm=comm)
	endif ! if(present(rMask))...

     endif ! if(present(iMask))...

  else ! no communicator handle present

     if(present(iMask)) then

	if(present(rMask)) then
	   call MaskedSpatialIntegralRAttrV_(inAv, integratedAv, SpatialWeights,  &
		                            iMask, rMask, UseFastMethod, .TRUE., &
					    'MaskedWeightsSum')
	else ! no rMask
	   call MaskedSpatialIntegralRAttrV_(inAv, integratedAv, SpatialWeights,  &
		                            iMask=iMask, UseFastMethod=UseFastMethod, &
                                            SumWeights=.TRUE., &
					    WeightSumTag='MaskedWeightsSum')
	endif ! if(present(rMask))...

     else ! no iMask present...

	if(present(rMask)) then
	   call MaskedSpatialIntegralRAttrV_(inAv, integratedAv, SpatialWeights,  &
		                            rMask=rMask, UseFastMethod=UseFastMethod, &
					    SumWeights=.TRUE., &
					    WeightSumTag='MaskedWeightsSum')
	else ! neither rMask nor iMask present:
	   call MaskedSpatialIntegralRAttrV_(inAv, integratedAv, SpatialWeights,  &
		                            UseFastMethod=UseFastMethod, &
					    SumWeights=.TRUE., &
					    WeightSumTag='MaskedWeightsSum')
	endif ! if(present(rMask))...

     endif ! if(present(iMask))...

  endif ! if(present(comm))...

       ! At this point, integratedAv containes the masked spatial integrals
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

       ! Divide by global weight sum to compute spatial averages from 
       ! spatial integrals.

  do i=1,AttrVect_nRAttr(outAv)
     outAv%rAttr(i,1) = integratedAv%rAttr(i,1) &
	                               / integratedAv%rAttr(iweight,1) 
  end do

       ! Clean up temporary AttrVect:

  call AttrVect_clean(integratedAv)

 end subroutine MaskedSpatialAverageRAttrV_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: PairedSpatialIntegralRAttrGG_ - Two spatial integrals.
!
! !DESCRIPTION:
! This routine computes spatial integrals of the {\tt REAL} attributes
! of the {\tt REAL} attributes of the input {\tt AttrVect} arguments 
! {\tt inAv1} and {\tt inAv2}, returning the integrals in the output 
! {\tt AttrVect} arguments {\tt outAv1} and {\tt outAv2}, respectively .  
! The integrals of {\tt inAv1} and {\tt inAv2} are computed using 
! spatial weights stored in the input {\tt GeneralGrid} arguments  
! {\tt GGrid1} and {\tt GGrid2}, respectively.  The spatial weights in 
! in {\tt GGrid1} and {\tt GGrid2} are identified by the input {\tt CHARACTER} 
! arguments {\tt WeightTag1} and {\tt WeightTag2}, respectively.  
! If {\tt SpatialIntegralRAttrGG\_()} is invoked with the optional 
! {\tt LOGICAL} input argument 
! {\tt SumWeights} set as {\tt .TRUE.}, then the weights are also summed 
! and stored in {\tt outAv1} and {\tt outAv2}, and can be referenced with 
! the attribute tags defined by the arguments {\tt WeightTag1} and 
! {\tt WeightTag2}, respectively.  This paired integral is implicitly a 
! distributed operation (the whole motivation for pairing the integrals is 
! to reduce communication latency costs), and the Fortran MPI communicator
! handle is defined by the input {\tt INTEGER} argument {\tt comm}.  The 
! summation is an AllReduce operation, with all processes receiving the 
! global sum.
!
! {\bf N.B.:  } The local lengths of the {\tt AttrVect} argument {\tt inAv1} 
! and the {\tt GeneralGrid} {\tt GGrid1} must be equal.  That is, there 
! must be a one-to-one correspondence between the field point values stored 
! in {\tt inAv1} and the point weights stored in {\tt GGrid1}.  The same
! relationship must apply between {\tt inAv2} and {\tt GGrid2}.
!
! {\bf N.B.:  }  If {\tt SpatialIntegralRAttrGG\_()} is invoked with the 
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

 subroutine PairedSpatialIntegralRAttrGG_(inAv1, outAv1, GGrid1, WeightTag1, &
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

  character(len=*),parameter :: myname_=myname//'::PairedSpatialIntegralRAttrGG_'

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


  call PairedSpatialIntegralRAttrV_(inAv1, outAv1, gridweights1, WeightTag1, &
	                           inAv2, outAv2, gridweights2, WeightTag2, &
                                   mySumWeights, comm)

       ! Clean up allocated arrays:

  deallocate(gridWeights1, gridWeights2, stat=ierr)
  if(ierr /= 0) then
     write(stderr,'(2a,i8)') myname_, &
	  'ERROR--deallocate(gridWeights1,...) failed, ierr = ',ierr
     call die(myname_)
  endif

 end subroutine PairedSpatialIntegralRAttrGG_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: PairedSpatialIntegralRAttrV_ - Two spatial integrals.
!
! !DESCRIPTION: 
! This routine computes spatial integrals of the {\tt REAL} attributes
! of the {\tt REAL} attributes of the input {\tt AttrVect} arguments 
! {\tt inAv1} and {\tt inAv2}, returning the integrals in the output 
! {\tt AttrVect} arguments {\tt outAv1} and {\tt outAv2}, respectively .  
! The integrals of {\tt inAv1} and {\tt inAv2} are computed using 
! spatial weights stored in the input {\tt REAL} array arguments  
! {\tt Weights1} and {\tt Weights2}, respectively.  
! If {\tt SpatialIntegralRAttrGG\_()} is invoked with the optional 
! {\tt LOGICAL} input argument 
! {\tt SumWeights} set as {\tt .TRUE.}, then the weights are also summed 
! and stored in {\tt outAv1} and {\tt outAv2}, and can be referenced with 
! the attribute tags defined by the arguments {\tt WeightName1} and 
! {\tt WeightName2}, respectively.  This paired integral is implicitly a 
! distributed operation (the whole motivation for pairing the integrals is 
! to reduce communication latency costs), and the Fortran MPI communicator
! handle is defined by the input {\tt INTEGER} argument {\tt comm}.  The 
! summation is an AllReduce operation, with all processes receiving the 
! global sum.
!
! {\bf N.B.:  } The local lengths of the {\tt AttrVect} argument {\tt inAv1} 
! and the input {\tt REAL} array {\tt Weights1} must be equal.  That is, there 
! must be a one-to-one correspondence between the field point values stored 
! in {\tt inAv1} and the point weights stored in {\tt Weights}.  The same
! relationship must apply between {\tt inAv2} and {\tt Weights2}.
!
! {\bf N.B.:  }  If {\tt SpatialIntegralRAttrGG\_()} is invoked with the 
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

 subroutine PairedSpatialIntegralRAttrV_(inAv1, outAv1, Weights1, WeightName1, &
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

  character(len=*),parameter :: myname_=myname//'::PairedSpatialIntegralRAttrV_'

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

 end subroutine PairedSpatialIntegralRAttrV_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: PairedSpatialAverageRAttrGG_ - Paired spatial averages on 
! {\tt GeneralGrids}
!
! !DESCRIPTION:
! This routine computes spatial averages of the {\tt REAL} attributes
! of the {\tt REAL} attributes of the input {\tt AttrVect} arguments 
! {\tt inAv1} and {\tt inAv2}, returning the integrals in the output 
! {\tt AttrVect} arguments {\tt outAv1} and {\tt outAv2}, respectively .  
! The integrals of {\tt inAv1} and {\tt inAv2} are computed using 
! spatial weights stored in the input {\tt GeneralGrid} arguments  
! {\tt GGrid1} and {\tt GGrid2}, respectively.  The spatial weights in 
! in {\tt GGrid1} and {\tt GGrid2} are identified by the input {\tt CHARACTER} 
! arguments {\tt WeightTag1} and {\tt WeightTag2}, respectively.  
! This paired average is implicitly a 
! distributed operation (the whole motivation for pairing the averages is 
! to reduce communication latency costs), and the Fortran MPI communicator
! handle is defined by the input {\tt INTEGER} argument {\tt comm}.  The 
! summation is an AllReduce operation, with all processes receiving the 
! global sum.
!
! {\bf N.B.:  } The local lengths of the {\tt AttrVect} argument {\tt inAv1} 
! and the {\tt GeneralGrid} {\tt GGrid1} must be equal.  That is, there 
! must be a one-to-one correspondence between the field point values stored 
! in {\tt inAv1} and the point weights stored in {\tt GGrid1}.  The same
! relationship must apply between {\tt inAv2} and {\tt GGrid2}.
!
! {\bf N.B.:  } The output {\tt AttrVect} arguments {\tt outAv1} and 
! {\tt outAv2} are allocated data structures.  The user must deallocate them
!  using the routine {\tt AttrVect\_clean()} when they are no longer needed.  
! Failure to do so will result in a memory leak.
!
! !INTERFACE:

 subroutine PairedSpatialAverageRAttrGG_(inAv1, outAv1, GGrid1, WeightTag1, &
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
!                 new interface to PairedSpatialIntegralRAttrGG_().
!
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::PairedSpatialAverageRAttrGG_'

  type(AttrVect) :: integratedAv1, integratedAv2
  integer :: i, ierr, iweight1, iweight2

       ! Compute the spatial integral:

     call PairedSpatialIntegralRAttrGG_(inAv1, integratedAv1, GGrid1, WeightTag1, &
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

       ! Divide by global weight sum to compute spatial averages from 
       ! spatial integrals.

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

 end subroutine PairedSpatialAverageRAttrGG_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: PairedSpatialAverageRAttrV_ - Two spatial averages.
!
! !DESCRIPTION:
! This routine computes spatial averages of the {\tt REAL} attributes
! of the {\tt REAL} attributes of the input {\tt AttrVect} arguments 
! {\tt inAv1} and {\tt inAv2}, returning the integrals in the output 
! {\tt AttrVect} arguments {\tt outAv1} and {\tt outAv2}, respectively .  
! The averages of {\tt inAv1} and {\tt inAv2} are computed using 
! spatial weights stored in the input {\tt REAL} array arguments  
! {\tt Weights1} and {\tt Weights2}, respectively.  This paired average 
! is implicitly a 
! distributed operation (the whole motivation for pairing the integrals is 
! to reduce communication latency costs), and the Fortran MPI communicator
! handle is defined by the input {\tt INTEGER} argument {\tt comm}.  The 
! summation is an AllReduce operation, with all processes receiving the 
! global sum.
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

 subroutine PairedSpatialAverageRAttrV_(inAv1, outAv1, Weights1, inAv2, &
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

  character(len=*),parameter :: myname_=myname//'::PairedSpatialAverageRAttrV_'

  type(AttrVect) :: integratedAv1, integratedAv2
  integer :: i, ierr, iweight1, iweight2

       ! weight tags used to keep track of spatial weight sums
  character*8, parameter :: WeightName1='WeightSum1'
  character*8, parameter :: WeightName2='WeightSum2'

       ! Compute the paired spatial integral, including spatial weights:

  call PairedSpatialIntegralRAttrV_(inAv1, integratedAv1, Weights1, WeightName1, &
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

       ! Divide by global weight sum to compute spatial averages from 
       ! spatial integrals.

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

 end subroutine PairedSpatialAverageRAttrV_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: PairedMaskedIntegralRAttrGG_ - Paired masked integrals.
!
! !DESCRIPTION:
! This routine computes a pair of masked spatial integrals of the {\tt REAL} 
! attributes of the input {\tt AttrVect} arguments {\tt inAv} and 
! {\tt inAv2}, returning the masked integrals in the output {\tt AttrVect} 
! {\tt outAv1} and {\tt outAv2}, respectively.  All of the spatial weighting
! and masking data for each set of integrals are assumed stored in the input 
! {\tt GeneralGrid} arguments {\tt GGrid} and {\tt GGrid2}.  If integer 
! masks are to be used, their integer  attribute names in {\tt GGrid1} 
! and {\tt GGrid2} are named as a colon-delimited lists in the optional 
! {\tt CHARACTER} input arguments {\tt iMaskTags1} and {\tt iMaskTags2}, 
! respectively.  Real masks (if desired) are referenced by their real 
! attribute names in {\tt GGrid1} and {\tt GGrid2} are named as 
! colon-delimited lists in the optional {\tt CHARACTER} input arguments 
! {\tt rMaskTags1} and {\tt rMaskTags2}, respectively.  The user specifies 
! a choice of mask combination method with the input {\tt LOGICAL} argument
! {\tt UseFastMethod}.  If ${\tt UseFastMethod} = {\tt .FALSE.}$ this 
! routine checks each mask entry to ensure that the integer masks contain
! only ones and zeroes, and that entries in the real masks are all in
! the closed interval $[0,1]$.  If ${\tt UseFastMethod} = {\tt .TRUE.}$,
! this routine performs direct products of the masks, assuming that the
! user has validated them in advance.  The optional {\tt LOGICAL} input 
! argument {\tt SumWeights} determines whether the masked sum of the spatial
! weights is computed and returned in {\tt outAv1} and {\tt outAv2} with the 
! real attribute names supplied in the {\tt CHARACTER} input arguments
! {\tt SpatialWeightTag1}, and {\tt SpatialWeightTag2}, respectively.  
! This paired integral is implicitly a distributed operation (the whole 
! motivation for pairing the averages is to reduce communication latency 
! costs), and the Fortran MPI communicator handle is defined by the input 
! {\tt INTEGER} argument {\tt comm}.  The 
! summation is an AllReduce operation, with all processes receiving the 
! global sum.
!
! {\bf N.B.:  } The local lengths of the {\tt AttrVect} argument {\tt inAv1} 
! and the {\tt GeneralGrid} {\tt GGrid1} must be equal.  That is, there 
! must be a one-to-one correspondence between the field point values stored 
! in {\tt inAv1} and the point weights stored in {\tt GGrid1}.  The same
! relationship must apply between {\tt inAv2} and {\tt GGrid2}.
!
! {\bf N.B.:  }  If {\tt PairedMaskedIntegralRAttrGG\_()} is invoked with the 
! optional {\tt LOGICAL} input argument {\tt SumWeights} set as {\tt .TRUE.}, 
! then the value of {\tt SpatialWeightTag1} must not conflict with any of the 
! {\tt REAL} attribute tags in {\tt inAv1} and the value of 
! {\tt SpatialWeightTag2} must not conflict with any of the {\tt REAL} 
! attribute tags in {\tt inAv2}.
!
! {\bf N.B.:  } The output {\tt AttrVect} arguments {\tt outAv1} and 
! {\tt outAv2} are allocated data structures.  The user must deallocate them
!  using the routine {\tt AttrVect\_clean()} when they are no longer needed.  
! Failure to do so will result in a memory leak.
!
! !INTERFACE:

 subroutine PairedMaskedIntegralRAttrGG_(inAv1, outAv1, GGrid1, &
                                         SpatialWeightTag1, rMaskTags1, &
                                         iMaskTags1, inAv2, outAv2, GGrid2, &
                                         SpatialWeightTag2, rMaskTags2, &
                                         iMaskTags2, UseFastMethod, &
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

      type(AttrVect),              intent(IN) :: inAv1
      type(GeneralGrid),           intent(IN) :: GGrid1
      character(len=*),            intent(IN) :: SpatialWeightTag1
      character(len=*),  optional, intent(IN) :: iMaskTags1
      character(len=*),  optional, intent(IN) :: rMaskTags1
      type(AttrVect),              intent(IN) :: inAv2
      type(GeneralGrid),           intent(IN) :: GGrid2
      character(len=*),            intent(IN) :: SpatialWeightTag2
      character(len=*),  optional, intent(IN) :: iMaskTags2
      character(len=*),  optional, intent(IN) :: rMaskTags2
      logical,                     intent(IN) :: UseFastMethod
      logical,           optional, intent(IN) :: SumWeights
      integer,                     intent(IN) :: comm

! !OUTPUT PARAMETERS:

      type(AttrVect),    intent(OUT) :: outAv1
      type(AttrVect),    intent(OUT) :: outAv2

! !REVISION HISTORY:
! 	17Jun02 - J.W. Larson <larson@mcs.anl.gov> - Initial version.
! 	19Jun02 - J.W. Larson <larson@mcs.anl.gov> - Shortened the name
!                 for compatibility with the Portland Group f90 compiler
!
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_ = &
                            myname//'::PairedMaskedIntegralRAttrGG_'

  logical :: mySumWeights
  real, dimension(:), pointer :: PairedBuffer, OutPairedBuffer
  integer :: ierr, nRA1, nRA2, PairedBufferLength

        ! Basic Argument Validity Checks:

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

        ! Are we summing the integration weights for the input
        ! GeneralGrids?

  if(present(SumWeights)) then
     mySumWeights = SumWeights
  else
     mySumWeights = .FALSE.
  endif

        ! Begin by invoking MaskedSpatialIntegralRAttrGG_() for each
        ! AttrVect/GeneralGrid pair.  This is done LOCALLY to create
        ! integratedAv1 and integratedAv2, respectively.

        ! Local Masked Integral #1:

  if(present(iMaskTags1)) then

     if(present(rMaskTags1)) then ! both Integer and Real Masking
	call MaskedSpatialIntegralRAttrGG_(inAv1, outAv1, GGrid1, &
	                                  SpatialWeightTag1, iMaskTags1, &
					  rMaskTags1, UseFastMethod,  &
					  mySumWeights, SpatialWeightTag1)
     else ! Integer Masking Only
	call MaskedSpatialIntegralRAttrGG_(inAv1, outAv1, GGrid1, &
	                                  SpatialWeightTag1, &
					  iMaskTags=iMaskTags1, &
					  UseFastMethod=UseFastMethod,  &
					  SumWeights=mySumWeights, &
					  WeightSumTag=SpatialWeightTag1)
     endif ! if(present(rMaskTags1))...

  else ! No Integer Masking

     if(present(rMaskTags1)) then ! Real Masking Only
	call MaskedSpatialIntegralRAttrGG_(inAv1, outAv1, GGrid1, &
	                                  SpatialWeightTag=SpatialWeightTag1, &
					  rMaskTags=rMaskTags1, &
					  UseFastMethod=UseFastMethod,  &
					  SumWeights=mySumWeights, &
					  WeightSumTag=SpatialWeightTag1)
     else ! Neither Integer nor Real Masking
	call MaskedSpatialIntegralRAttrGG_(inAv1, outAv1, GGrid1, &
	                                  SpatialWeightTag=SpatialWeightTag1, &
					  UseFastMethod=UseFastMethod,  &
					  SumWeights=mySumWeights, &
					  WeightSumTag=SpatialWeightTag1)

     endif ! if(present(rMaskTags1))...

  endif ! if(present(iMaskTags1))...

        ! Local Masked Integral #2:

  if(present(iMaskTags2)) then

     if(present(rMaskTags2)) then ! both Integer and Real Masking
	call MaskedSpatialIntegralRAttrGG_(inAv2, outAv2, GGrid2, &
	                                  SpatialWeightTag2, iMaskTags2, &
					  rMaskTags2, UseFastMethod,  &
					  mySumWeights, SpatialWeightTag2)
     else ! Integer Masking Only
	call MaskedSpatialIntegralRAttrGG_(inAv2, outAv2, GGrid2, &
	                                  SpatialWeightTag2, &
					  iMaskTags=iMaskTags2, &
					  UseFastMethod=UseFastMethod,  &
					  SumWeights=mySumWeights, &
					  WeightSumTag=SpatialWeightTag2)
     endif ! if(present(rMaskTags2))...

  else ! No Integer Masking

     if(present(rMaskTags2)) then ! Real Masking Only
	call MaskedSpatialIntegralRAttrGG_(inAv2, outAv2, GGrid2, &
	                                  SpatialWeightTag=SpatialWeightTag2, &
					  rMaskTags=rMaskTags2, &
					  UseFastMethod=UseFastMethod,  &
					  SumWeights=mySumWeights, &
					  WeightSumTag=SpatialWeightTag2)
     else ! Neither Integer nor Real Masking
	call MaskedSpatialIntegralRAttrGG_(inAv2, outAv2, GGrid2, &
	                                  SpatialWeightTag=SpatialWeightTag2, &
					  UseFastMethod=UseFastMethod,  &
					  SumWeights=mySumWeights, &
					  WeightSumTag=SpatialWeightTag2)

     endif ! if(present(rMaskTags2))...

  endif ! if(present(iMaskTags2))...

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

  deallocate(PairedBuffer, OutPairedBuffer, stat=ierr)
  if(ierr /= 0) then
     write(stderr,'(2a,i8)') myname_, &
	  ':: Fatal error--deallocate(PairedBuffer...failed, ierr = ',ierr
     call die(myname_)
  endif

 end subroutine PairedMaskedIntegralRAttrGG_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: PairedMaskedAverageRAttrGG_ - Paired masked averages.
!
! !DESCRIPTION:
! This routine computes a pair of masked spatial averages of the {\tt REAL} 
! attributes of the input {\tt AttrVect} arguments {\tt inAv} and 
! {\tt inAv2}, returning the masked averagess in the output {\tt AttrVect} 
! {\tt outAv1} and {\tt outAv2}, respectively.  All of the spatial weighting
! and masking data for each set of averages are assumed stored in the input 
! {\tt GeneralGrid} arguments {\tt GGrid} and {\tt GGrid2}.  If integer 
! masks are to be used, their integer  attribute names in {\tt GGrid1} 
! and {\tt GGrid2} are named as a colon-delimited lists in the optional 
! {\tt CHARACTER} input arguments {\tt iMaskTags1} and {\tt iMaskTags2}, 
! respectively.  Real masks (if desired) are referenced by their real 
! attribute names in {\tt GGrid1} and {\tt GGrid2} are named as 
! colon-delimited lists in the optional {\tt CHARACTER} input arguments 
! {\tt rMaskTags1} and {\tt rMaskTags2}, respectively.  The user specifies 
! a choice of mask combination method with the input {\tt LOGICAL} argument
! {\tt UseFastMethod}.  If ${\tt UseFastMethod} = {\tt .FALSE.}$ this 
! routine checks each mask entry to ensure that the integer masks contain
! only ones and zeroes, and that entries in the real masks are all in
! the closed interval $[0,1]$.  If ${\tt UseFastMethod} = {\tt .TRUE.}$,
! this routine performs direct products of the masks, assuming that the
! user has validated them in advance.  This paired average is implicitly 
! a distributed operation (the whole motivation for pairing the averages 
! is to reduce communication latency costs), and the Fortran MPI communicator 
! handle is defined by the input {\tt INTEGER} argument {\tt comm}.  The 
! summation is an AllReduce operation, with all processes receiving the 
! global sum.
!
! {\bf N.B.:  } The local lengths of the {\tt AttrVect} argument {\tt inAv1} 
! and the {\tt GeneralGrid} {\tt GGrid1} must be equal.  That is, there 
! must be a one-to-one correspondence between the field point values stored 
! in {\tt inAv1} and the point weights stored in {\tt GGrid1}.  The same
! relationship must apply between {\tt inAv2} and {\tt GGrid2}.
!
! {\bf N.B.:  } The output {\tt AttrVect} arguments {\tt outAv1} and 
! {\tt outAv2} are allocated data structures.  The user must deallocate them
!  using the routine {\tt AttrVect\_clean()} when they are no longer needed.  
! Failure to do so will result in a memory leak.
!
! !INTERFACE:

 subroutine PairedMaskedAverageRAttrGG_(inAv1, outAv1, GGrid1, &
                                        SpatialWeightTag1, rMaskTags1, &
                                        iMaskTags1, inAv2, outAv2, GGrid2, &
                                        SpatialWeightTag2, rMaskTags2, &
                                        iMaskTags2, UseFastMethod, &
                                        comm)
! ! USES:

      use m_stdio
      use m_die
      use m_mpif90

      use m_AttrVect, only : AttrVect
      use m_AttrVect, only : AttrVect_init => init
      use m_AttrVect, only : AttrVect_clean => clean
      use m_AttrVect, only : AttrVect_exportRListToChar => exportRListToChar
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

      type(AttrVect),              intent(IN) :: inAv1
      type(GeneralGrid),           intent(IN) :: GGrid1
      character(len=*),            intent(IN) :: SpatialWeightTag1
      character(len=*),  optional, intent(IN) :: iMaskTags1
      character(len=*),  optional, intent(IN) :: rMaskTags1
      type(AttrVect),              intent(IN) :: inAv2
      type(GeneralGrid),           intent(IN) :: GGrid2
      character(len=*),            intent(IN) :: SpatialWeightTag2
      character(len=*),  optional, intent(IN) :: iMaskTags2
      character(len=*),  optional, intent(IN) :: rMaskTags2
      logical,                     intent(IN) :: UseFastMethod
      integer,                     intent(IN) :: comm

! !OUTPUT PARAMETERS:

      type(AttrVect),    intent(OUT) :: outAv1
      type(AttrVect),    intent(OUT) :: outAv2

! !REVISION HISTORY:
! 	17Jun02 - J.W. Larson <larson@mcs.anl.gov> - Initial version.
! 	19Jun02 - J.W. Larson <larson@mcs.anl.gov> - Shortened the name
!                 for compatibility with the Portland Group f90 compiler
! 	25Jul02 - J.W. Larson E.T. Ong - Bug fix.  This routine was 
!                 previously doing integrals rather than area averages.
!
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_ = &
                            myname//'::PairedMaskedAverageRAttrGG_'

  type(AttrVect) :: LocalIntegral1, LocalIntegral2
  real, dimension(:), pointer :: PairedBuffer, OutPairedBuffer
  integer :: i, ierr, nRA1, nRA2, PairedBufferLength
  real :: WeightSumInv

        ! Basic Argument Validity Checks:

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

        ! Begin by invoking MaskedSpatialIntegralRAttrGG_() for each
        ! AttrVect/GeneralGrid pair.  This is done LOCALLY to create
        ! LocalIntegral1 and LocalIntegral2, respectively.

        ! Local Masked Integral #1:

  if(present(iMaskTags1)) then

     if(present(rMaskTags1)) then ! both Integer and Real Masking
	call MaskedSpatialIntegralRAttrGG_(inAv1, LocalIntegral1, GGrid1, &
	                                  SpatialWeightTag1, iMaskTags1, &
					  rMaskTags1, UseFastMethod,  &
					  .TRUE., SpatialWeightTag1)
     else ! Integer Masking Only
	call MaskedSpatialIntegralRAttrGG_(inAv1, LocalIntegral1, GGrid1, &
	                                  SpatialWeightTag1, &
					  iMaskTags=iMaskTags1, &
					  UseFastMethod=UseFastMethod,  &
					  SumWeights=.TRUE., &
					  WeightSumTag=SpatialWeightTag1)
     endif ! if(present(rMaskTags1))...

  else ! No Integer Masking

     if(present(rMaskTags1)) then ! Real Masking Only
	call MaskedSpatialIntegralRAttrGG_(inAv1, LocalIntegral1, GGrid1, &
	                                  SpatialWeightTag=SpatialWeightTag1, &
					  rMaskTags=rMaskTags1, &
					  UseFastMethod=UseFastMethod,  &
					  SumWeights=.TRUE., &
					  WeightSumTag=SpatialWeightTag1)
     else ! Neither Integer nor Real Masking
	call MaskedSpatialIntegralRAttrGG_(inAv1, LocalIntegral1, GGrid1, &
	                                  SpatialWeightTag=SpatialWeightTag1, &
					  UseFastMethod=UseFastMethod,  &
					  SumWeights=.TRUE., &
					  WeightSumTag=SpatialWeightTag1)

     endif ! if(present(rMaskTags1))...

  endif ! if(present(iMaskTags1))...

        ! Local Masked Integral #2:

  if(present(iMaskTags2)) then

     if(present(rMaskTags2)) then ! both Integer and Real Masking
	call MaskedSpatialIntegralRAttrGG_(inAv2, LocalIntegral2, GGrid2, &
	                                  SpatialWeightTag2, iMaskTags2, &
					  rMaskTags2, UseFastMethod,  &
					  .TRUE., SpatialWeightTag2)
     else ! Integer Masking Only
	call MaskedSpatialIntegralRAttrGG_(inAv2, LocalIntegral2, GGrid2, &
	                                  SpatialWeightTag2, &
					  iMaskTags=iMaskTags2, &
					  UseFastMethod=UseFastMethod,  &
					  SumWeights=.TRUE., &
					  WeightSumTag=SpatialWeightTag2)
     endif ! if(present(rMaskTags2))...

  else ! No Integer Masking

     if(present(rMaskTags2)) then ! Real Masking Only
	call MaskedSpatialIntegralRAttrGG_(inAv2, LocalIntegral2, GGrid2, &
	                                  SpatialWeightTag=SpatialWeightTag2, &
					  rMaskTags=rMaskTags2, &
					  UseFastMethod=UseFastMethod,  &
					  SumWeights=.TRUE., &
					  WeightSumTag=SpatialWeightTag2)
     else ! Neither Integer nor Real Masking
	call MaskedSpatialIntegralRAttrGG_(inAv2, LocalIntegral2, GGrid2, &
	                                  SpatialWeightTag=SpatialWeightTag2, &
					  UseFastMethod=UseFastMethod,  &
					  SumWeights=.TRUE., &
					  WeightSumTag=SpatialWeightTag2)

     endif ! if(present(rMaskTags2))...

  endif ! if(present(iMaskTags2))...

       ! Create the paired buffer for the Global Sum

  nRA1 = AttrVect_nRAttr(LocalIntegral1)
  nRA2 = AttrVect_nRAttr(LocalIntegral2)

  PairedBufferLength =  nRA1 + nRA2
  allocate(PairedBuffer(PairedBufferLength), OutPairedBuffer(PairedBufferLength), &
           stat=ierr)
  if(ierr /= 0) then
     write(stderr,'(2a,i8)') myname_, &
	  ':: Fatal error--allocate(PairedBuffer...failed, ierr = ',ierr
     call die(myname_)
  endif

       ! Load the paired buffer

  PairedBuffer(1:nRA1) = LocalIntegral1%rAttr(1:nRA1,1)
  PairedBuffer(nRA1+1:PairedBufferLength) = LocalIntegral2%rAttr(1:nRA2,1)

       ! Perform the global sum on the paired buffer

  call MPI_AllReduce(PairedBuffer, OutPairedBuffer, PairedBufferLength, &
                        MP_Type(PairedBuffer(1)), MP_SUM, comm, ierr)
  if(ierr /= 0) then
     write(stderr,'(2a,i8)') myname_, &
  	      ':: Fatal Error--MPI_ALLREDUCE() failed with ierror = ',ierr
     call MP_perr_die(myname_,'MPI_ALLREDUCE() failed',ierr)
  endif

       ! Create outAv1 and outAv2 from inAv1 and inAv2, respectively:

  call AttrVect_init(outAv1, rList=AttrVect_exportRListToChar(inAv1), lsize=1)
  call AttrVect_init(outAv2, rList=AttrVect_exportRListToChar(inAv2), lsize=1)

       ! Unload/rescale OutPairedBuffer into outAv1 and outAv2:

  nRA1 = AttrVect_nRAttr(outAv1)
  nRA2 = AttrVect_nRAttr(outAv2)

       ! First outAv1:

  if(OutPairedBuffer(nRA1+1) /= 0.) then
     WeightSumInv = 1. / OutPairedBuffer(nRA1+1) ! Sum of weights on grid1
                                                 ! is the nRA1+1th element in
                                                 ! the paired buffer.
  else
     write(stderr,'(2a)') myname_, &
	  ':: FATAL ERROR--Sum of the Weights for integral #1 is zero!  Terminating...'
     call die(myname_)
  endif

       ! Rescale global integral to get global average:

  do i=1,nRA1
     outAv1%rAttr(i,1) = WeightSumInv * OutPairedBuffer(i)
  end do

       ! And then outAv2:

  if(OutPairedBuffer(PairedBufferLength) /= 0.) then
     WeightSumInv = 1. / OutPairedBuffer(PairedBufferLength) ! Sum of weights on grid2
                                                             ! is the last element in
                                                             ! the paired buffer.
  else
     write(stderr,'(2a)') myname_, &
	  ':: FATAL ERROR--Sum of the Weights for integral #2 is zero!  Terminating...'
     call die(myname_)
  endif

       ! Rescale global integral to get global average:

  do i=1,nRA2
     outAv2%rAttr(i,1) = WeightSumInv * OutPairedBuffer(i+nRA1+1)
  end do

       ! Clean up allocated structures

  call AttrVect_clean(LocalIntegral1)
  call AttrVect_clean(LocalIntegral2)

  deallocate(PairedBuffer, OutPairedBuffer, stat=ierr)
  if(ierr /= 0) then
     write(stderr,'(2a,i8)') myname_, &
	  ':: Fatal error--deallocate(PairedBuffer...failed, ierr = ',ierr
     call die(myname_)
  endif

 end subroutine PairedMaskedAverageRAttrGG_

 end module m_SpatialIntegral



