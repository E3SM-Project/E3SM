!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !ROUTINE: accumulate--Acumulate from an AttrVect to an Accumulator.
!
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine accumulate(aV, aC, action, ier)
!
! !USES:
!
      use m_stdio, only : stdout,stderr

      use m_AttrVect, only : AttrVect
      use m_AttrVect, only : AttrVect_lsize => lsize
      use m_AttrVect, only : AttrVect_indexRA => indexRA
      use m_AttrVect, only : AttrVect_indexIA => indexIA

      use m_Accumulator, only : Accumulator
      use m_Accumulator, only : Accumulator_lsize => lsize
      use m_Accumulator, only : Accumulator_indexRA => indexRA
      use m_Accumulator, only : Accumulator_indexIA => indexIA

      implicit none

      type(AttrVect), intent(in)     :: aV      ! Input AttrVect
      type(Accumulator), intent(out) :: aC      ! Output Accumulator
      character(len=*), intent(in)   :: action  ! Error Flag
      integer, intent(out)           :: ier     ! Error Flag

! !REVISION HISTORY:
!       18Sep00 - J.W. Larson <larson@mcs.anl.gov> -- initial version.
!
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='accumulate'    ! a stamp

        ! Sanity check of arguments:

  if(Accumulator_lsize(aC) /= AttrVect_lsize(aV)) then
     ier = -1
  endif

  if(aC%num_steps == 0) then
     ier = -1
  endif

  if(aC%av%rList /= ' ') then ! Accumulate only if fields are requested
     if(aC%av%rList == aV%rList) then
	do n=1,AttrVect_nRAttr(av)
	   do l=1,AttrVect_lsize(av)
	      aC%av%rList(n,l) = aC%av%rList(n,l) + aV%rList(n,l)
	   enddo
	enddo
     endif
  endif

  if(aC%av%rList /= ' ') then ! Accumulate only if fields are requested
     if(aC%av%iList == aV%iList) then
	do n=1,AttrVect_nIAttr(av)
	   do l=1,AttrVect_lsize(av)
	      aC%av%iList(n,l) = aC%av%iList(n,l) + aV%iList(n,l)
	   enddo
	enddo
     endif
  endif

        ! Increment aC%steps_done:

  aC%steps_done = aC%stepsp_done + 1

        ! If we're at the end of an averaging period, compute the
        ! average.

  if(aC%steps_done == aC%num_steps) then
     select case(action)
     case('average','AVERAGE')
	scale = 1 / float(aC%num_steps)
	do n=1,Accumulator_nRAttr(aC)
	   do l=1,Accumulator_lsize(aC)
	      aC%av%rList(n,l) = scale * aC%av%rList(n,l)
	   enddo
	enddo
	do n=1,AttrVect_nIAttr(aC)
	   do l=1,AttrVect_lsize(aC)
	      aC%av%rList(n,l) = aC%av%rList(n,l) / num_steps
	   enddo
	enddo
     case default
     end select
  endif

 end subroutine accumulate
