! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copyright (c) 2015, Regents of the University of Colorado
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without modification, are 
! permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, this list of 
!    conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice, this list
!    of conditions and the following disclaimer in the documentation and/or other 
!    materials provided with the distribution.
!
! 3. Neither the name of the copyright holder nor the names of its contributors may be 
!    used to endorse or promote products derived from this software without specific prior
!    written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY 
! EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
! MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL 
! THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT 
! OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! History:
! 10/16/03  John Haynes   - Original version (haynes@atmos.colostate.edu)
! 01/31/06  John Haynes   - IDL to Fortran 90
! 01/01/15  Dustin Swales - Modified for COSPv2.0
! 
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
module array_lib
  USE COSP_KINDS, ONLY: wp
  implicit none
contains

  ! ############################################################################
  !                               function INFIND
  ! ############################################################################
  function infind(list,val)
    implicit none
    ! ##########################################################################
    ! Purpose:
    !   Finds the index of an array that is closest to a value, plus the
    !   difference between the value found and the value specified
    !
    ! Inputs:
    !   [list]   an array of sequential values
    !   [val]    a value to locate
    ! Optional input:
    !   [sort]   set to 1 if [list] is in unknown/non-sequential order
    !
    ! Returns:
    !   index of [list] that is closest to [val]
    !
    ! Optional output:
    !   [dist]   set to variable containing [list([result])] - [val]
    !
    ! Requires:
    !   mrgrnk library
    !
    ! ##########################################################################

    ! INPUTS
    real(wp), dimension(:), intent(in) :: &
         list   ! An array of sequential values
    real(wp), intent(in) :: &
         val    ! A value to locate
    ! OUTPUTS
    integer :: &
         infind ! Index of [list] that is closest to [val]

    ! Internal Variables
    real(wp), dimension(size(list)) :: lists
    integer :: nlist, result, tmp(1), sort_list
    integer, dimension(size(list)) :: mask
    
    sort_list = 0
    
    nlist = size(list)
    lists = list
    
    if (val >= lists(nlist)) then
       result = nlist
    else if (val <= lists(1)) then
       result = 1
    else
       mask(:) = 0
       where (lists < val) mask = 1
       tmp = minloc(mask,1)
       if (abs(lists(tmp(1)-1)-val) < abs(lists(tmp(1))-val)) then
          result = tmp(1) - 1
      else
         result = tmp(1)
      endif
   endif
   infind = result
 end function infind

end module array_lib
