!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_column
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 12th January 2015
!> \details
!>  
!
!-----------------------------------------------------------------------

module cice_column

  use mpas_grid_types

  implicit none

  private
  save

  public :: &
       cice_column_test

contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_column_test
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 12th January 2015
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_column_test()

    use test_mod

    call test_sub()

  end subroutine cice_column_test

!-----------------------------------------------------------------------

end module cice_column


