
      module MO_PRD_MAP

      CONTAINS

      subroutine PRD_MAP( template )
!-----------------------------------------------------------------------
!	... Form production indicies
!-----------------------------------------------------------------------
     
      use VAR_MOD, only : var_lim
      use RXT_MOD, only : prd_lim

      implicit none

!-----------------------------------------------------------------------
!	... Dummy args
!-----------------------------------------------------------------------
      integer, intent(inout) ::      template(:,:)
      
!-----------------------------------------------------------------------
!	... Local variables
!-----------------------------------------------------------------------
      integer  ::  i, clsno
      
      integer  ::  XLATE

      template(:,:2) = 0

      do i = 1,prd_lim
         if( template(i,3) < 0 ) then
            cycle
         else if( template(i,3) == 0 ) then
            exit
         else
            clsno = XLATE( template(i,3) )
            template(i,2) = clsno
            template(clsno,1) = template(clsno,1) + 1
         end if
      end do
      
      end subroutine PRD_MAP

      end module MO_PRD_MAP
