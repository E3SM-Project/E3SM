!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !ROUTINE: accumulate--Acumulate from an AttrVect to an Accumulator.
!
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine accumulate(aV, aC, action)

!
! !USES:
!
      use m_stdio, only : stdout,stderr
      use m_die,   only : MP_perr_die

      use m_AttrVect, only : AttrVect
      use m_AttrVect, only : AttrVect_lsize => lsize
      use m_AttrVect, only : AttrVect_nIAttr => nIAttr
      use m_AttrVect, only : AttrVect_nRAttr => nRAttr
      use m_AttrVect, only : AttrVect_indexRA => indexRA
      use m_AttrVect, only : AttrVect_indexIA => indexIA

      use m_Accumulator, only : Accumulator
      use m_Accumulator, only : Accumulator_lsize => lsize
      use m_Accumulator, only : Accumulator_nIAttr => nIAttr
      use m_Accumulator, only : Accumulator_nRAttr => nRAttr
      use m_Accumulator, only : Accumulator_indexRA => indexRA
      use m_Accumulator, only : Accumulator_indexIA => indexIA

      use m_SharedAttrIndices, only : SharedAttrIndexList

      implicit none

      type(AttrVect),     intent(in)    :: aV      ! Input AttrVect
      type(Accumulator),  intent(inout) :: aC      ! Output Accumulator
      character(len=*),   intent(in)    :: action  ! action to be taken.

! !REVISION HISTORY:
!       18Sep00 - J.W. Larson <larson@mcs.anl.gov> -- initial version.
!       07Feb01 - J.W. Larson <larson@mcs.anl.gov> -- General version.
!
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='accumulate'

! Overlapping attribute index number
  integer :: num_indices

! Overlapping attribute index storage arrays:
  integer, dimension(:), pointer :: aCindices, aVindices
  integer :: aCindex, aVindex

! Error flag and loop indices
  integer :: ierr, l, n

! Averaging time-weighting factor:
  real :: step_weight
  integer :: num_steps


        ! Sanity check of arguments:

  if(Accumulator_lsize(aC) /= AttrVect_lsize(aV)) then
     write(stderr,'(2a,2(1a,1i7))') myname, &
	  ":: mismatched Accumulator / AttrVect sizes.", &
	  "AttrVect_lsize(aV)=",AttrVect_lsize(aV), &
	  "Accumulator_lsize(aC)=",Accumulator_lsize(aC)
     ierr = 1
     call MP_perr_die(myname,'Accum. vs. AttrVect sizes.',ierr)
  endif

  if(aC%num_steps == 0) then
     write(stderr,'(2a)') myname, &
	  ":: Zero steps in Accumulation cycle."
     ierr = 2
     call MP_perr_die(myname,'Zero steps in accumulation cycle.',ierr)
  endif

        ! Accumulation of REAL attribute data:

  if(Accumulator_nRAttr(aC) /= 0) then ! Accumulate only if fields 
                                       ! are present

     call SharedAttrIndexList(aV, aC, 'REAL', num_indices, &
	                      aVindices, aCindices)

     if(num_indices /= 0) then

	do n=1,num_indices
	   aVindex = aVindices(n)
	   aCindex = aCindices(n)
	   do l=1,AttrVect_lsize(aV)
	      aC%av%rAttr(aCindex,l) = aC%av%rAttr(aCindex,l) + &
		      aV%rAttr(aVindex,l)
	   end do
	end do

	deallocate(aVindices, aCindices, stat=ierr)
	if(ierr /= 0) then
	   call MP_perr_die(myname,'first deallocate(aVindices...',ierr)
	endif

     endif ! if(num_indices /= 0)

  endif ! if(Accumulator_nRAttr(aC) /= 0)


        ! Accumulation of INTEGER attribute data:

  if(Accumulator_nIAttr(aC) /= 0) then ! Accumulate only if fields 
                                       ! are present

     call SharedAttrIndexList(aV, aC, 'INTEGER', num_indices, &
	                      aVindices, aCindices)

     if(num_indices /= 0) then

	do n=1,num_indices
	   aVindex = aVindices(n)
	   aCindex = aCindices(n)
	   do l=1,AttrVect_lsize(aV)
	      aC%av%iAttr(aCindex,l) = aC%av%iAttr(aCindex,l) + &
		      aV%iAttr(aVindex,l)
	   end do
	end do

	deallocate(aVindices, aCindices, stat=ierr)
	if(ierr /= 0) then
	   call MP_perr_die(myname,'second deallocate(aVindices...',ierr)
	endif

     endif ! if(num_indices /= 0)

  endif ! if(Accumulator_nIAttr(aC) /= 0)

        ! Increment aC%steps_done:

  aC%steps_done = aC%steps_done + 1

        ! If we are at the end of an averaging period, compute the
        ! average (if desired).

  if(aC%steps_done == aC%num_steps) then
     select case(action)
     case('average','AVERAGE')
	step_weight = 1 / float(aC%num_steps)
	do n=1,Accumulator_nRAttr(aC)
	   do l=1,Accumulator_lsize(aC)
	      aC%av%rAttr(n,l) = step_weight * aC%av%rAttr(n,l)
	   enddo
	enddo
	do n=1,Accumulator_nIAttr(aC)
	   do l=1,Accumulator_lsize(aC)
	      aC%av%iAttr(n,l) = aC%av%iAttr(n,l) / num_steps
	   enddo
	enddo
     case default
     end select
  endif

 end subroutine accumulate


