!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !ROUTINE: accumulate--Acumulate from an AttrVect to an Accumulator.
!
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine accumulate(aV, aC)

!
! !USES:
!
      use m_stdio, only : stdout,stderr
      use m_die,   only : die

      use m_AttrVect, only : AttrVect
      use m_AttrVect, only : AttrVect_lsize => lsize
      use m_AttrVect, only : AttrVect_nIAttr => nIAttr
      use m_AttrVect, only : AttrVect_nRAttr => nRAttr
      use m_AttrVect, only : AttrVect_indexRA => indexRA
      use m_AttrVect, only : AttrVect_indexIA => indexIA

      use m_Accumulator, only : Accumulator
      use m_Accumulator, only : MCT_SUM
      use m_Accumulator, only : MCT_AVG
      use m_Accumulator, only : Accumulator_lsize => lsize
      use m_Accumulator, only : Accumulator_nIAttr => nIAttr
      use m_Accumulator, only : Accumulator_nRAttr => nRAttr
      use m_Accumulator, only : Accumulator_indexRA => indexRA
      use m_Accumulator, only : Accumulator_indexIA => indexIA

      use m_SharedAttrIndices, only : SharedAttrIndexList

      implicit none

      type(AttrVect),     intent(in)    :: aV      ! Input AttrVect
      type(Accumulator),  intent(inout) :: aC      ! Output Accumulator

! !REVISION HISTORY:
!       18Sep00 - J.W. Larson <larson@mcs.anl.gov> -- initial version.
!       07Feb01 - J.W. Larson <larson@mcs.anl.gov> -- General version.
!       10Jun01 - E.T. Ong -- fixed divide-by-zero problem in integer
!                 attribute accumulation.
!       27Jul01 - E.T. Ong <eong@mcs.anl.gov> -- removed action argument.
!                 Make compatible with new Accumulator type.
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

! Character variable used as a data type flag:
  character*7 :: data_flag

        ! Sanity check of arguments:

  if(Accumulator_lsize(aC) /= AttrVect_lsize(aV)) then
     call die(myname, &
     "mismatched Accumulator  AttrVect sizes. AttrVect_lsize(aV) = ",&
     AttrVect_lsize(aV), "Accumulator_lsize(aC) = ", &
     Accumulator_lsize(aC))
  endif

  if(aC%num_steps == 0) then
     call die(myname, 'Zero steps in accumulation cycle.')
  endif

        ! Set num_steps from aC:

  num_steps = aC%num_steps

        ! Accumulation of REAL attribute data:

  if( aC%nrAction(MCT_SUM)+aC%nrAction(MCT_AVG) > 0 ) then 
        
        ! Accumulate only if fields are present 

     data_flag = 'REAL'
     call SharedAttrIndexList(aV, aC, data_flag, num_indices, &
	                      aVindices, aCindices)

     if(num_indices > 0) then
	do n=1,num_indices
	   aVindex = aVindices(n)
	   aCindex = aCindices(n)

	   ! Accumulate if the action is MCT_SUM or MCT_AVG
	   if( (aC%rAction(aCindex) == MCT_SUM).or. &
               (aC%rAction(aCindex) == MCT_AVG) ) then
              do l=1,AttrVect_lsize(aV)
		 aC%av%rAttr(aCindex,l) = aC%av%rAttr(aCindex,l) + &
		      aV%rAttr(aVindex,l)
	      end do
	   endif
	end do

	deallocate(aVindices, aCindices, stat=ierr)
	if(ierr /= 0) call die(myname,'first deallocate(aVindices...',ierr)

     endif ! if(num_indices > 0)

  endif ! if(aC%nrAction(MCT_SUM)+aC%nrAction(MCT_AVG) > 0)


        ! Accumulation of INTEGER attribute data:

  if( aC%niAction(MCT_SUM)+aC%niAction(MCT_AVG) > 0 ) then 

        ! Accumulate only if fields are present

     data_flag = 'INTEGER'
     call SharedAttrIndexList(aV, aC, data_flag, num_indices, &
	                      aVindices, aCindices)

     if(num_indices > 0) then

	do n=1,num_indices
	   aVindex = aVindices(n)
	   aCindex = aCindices(n)

	   ! Accumulate if the action is MCT_SUM or MCT_AVG
	   if( (aC%iAction(aCindex) == MCT_SUM) .or. &
               (aC%iAction(aCindex) == MCT_AVG) ) then
	      do l=1,AttrVect_lsize(aV)
		 aC%av%iAttr(aCindex,l) = aC%av%iAttr(aCindex,l) + &
		      aV%iAttr(aVindex,l)
	      end do
	   endif
	end do

	deallocate(aVindices, aCindices, stat=ierr)
	if(ierr /= 0) call die(myname,'second deallocate(aVindices...',ierr)

     endif ! if(num_indices > 0)

  endif ! if(aC%niAction(MCT_SUM)+aC%niAction(MCT_AVG) > 0 )

        ! Increment aC%steps_done:

  aC%steps_done = aC%steps_done + 1

        ! If we are at the end of an averaging period, compute the
        ! average (if desired).

  if(aC%steps_done == num_steps) then

     if( aC%nrAction(MCT_AVG) > 0 ) then
	step_weight = 1 / float(num_steps)
	do n=1,Accumulator_nRAttr(aC)
           if( aC%rAction(n) == MCT_AVG ) then
              do l=1,Accumulator_lsize(aC)
                 aC%av%rAttr(n,l) = step_weight * aC%av%rAttr(n,l)
              enddo
           endif
	enddo
     endif
     
     if( aC%niAction(MCT_AVG) > 0 ) then
	do n=1,Accumulator_nIAttr(aC)
           if( aC%iAction(n) == MCT_AVG ) then
              do l=1,Accumulator_lsize(aC)
                 aC%av%iAttr(n,l) = aC%av%iAttr(n,l) / num_steps
              enddo
           endif
	enddo
     endif

  endif

 end subroutine accumulate









