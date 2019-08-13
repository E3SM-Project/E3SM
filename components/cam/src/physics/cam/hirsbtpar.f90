
module hirsbtpar

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Initialize HIRS brightness temperature parameters for hirsrtm package
! 
! Author: 
! 
!-----------------------------------------------------------------------

  implicit none
  public
  save
 
! Parameter sizes for TOVS brightness temperature calculations

   integer, parameter :: pnb_hirs = 7    ! Number of TOVS/HIRS channels
   integer, parameter :: pnf_msu = 4     ! Number of TOVS/MSU channels
   integer, parameter :: msu_flag = 1    ! Flag to include MSU channels
                                         ! 1 = include; 0 = exclude 

   logical ::  dohirs                    ! Flag to call HIRSRTM
                                         ! .true. = call HIRSRTM
   integer :: ihirsfq                    ! Frequency to call the HIRSRTM routine
                                         ! Positive = time steps; Negative = hours
   character*8 hirsname(pnb_hirs)        ! Output field names for HIRS channels
   character*8 msuname(pnf_msu)          ! Output field names for MSU channels

end module hirsbtpar
