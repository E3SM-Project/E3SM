subroutine diffuse_mom

!  Interface to the diffusion routines

use vars
implicit none
integer i,j,k

!call t_startf ('diffuse_mom')

if(RUN3D) then
!   call diffuse_mom3D()
   call diffuse_mom3D_xy()
   call diffuse_mom3D_z()
else
!   call diffuse_mom2D()
   call diffuse_mom2D_xy()
   call diffuse_mom2D_z()
endif

!call t_stopf ('diffuse_mom')

end subroutine diffuse_mom

