!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: MCT_Recv
!
! !DESCRIPTION:
! Recieve into the AttrVect the data coming from the component
! specified in the Router.  An Error will result if the
! attribute list of the incoming data doesn't match any of
! the attributes in the argument AttrVect.
! Requires a corresponding MCT_Send to be called on the other component.
!
! !INTERFACE:

 subroutine MCT_Recv(AtrVc, Rout)
!
! !USES:
!
      use m_Router,only  : Router
      use m_AttrVect,only : AttrVect
      use m_list,only:	List

      implicit none
      Type(AttrVect),intent(inout) :: 	AtrVc
      Type(Router),intent(in) ::	Rout

! !REVISION HISTORY:
!      07Feb01 - R. Jacob <jacob@mcs.anl.gov> - initial prototype
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='MCT_Recv'

end subroutine
