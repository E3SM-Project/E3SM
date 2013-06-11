! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=============================================================================================
 module mpas_atmphys_update
 use mpas_configure
 use mpas_grid_types

 use mpas_atmphys_driver_convection_deep
 use mpas_atmphys_vars

 implicit none
 private
 public:: physics_update, &
          update_radiation_diagnostics

 contains
 
!=============================================================================================
 subroutine physics_update(domain,dt)
!=============================================================================================

!input arguments:
 type(domain_type),intent(in):: domain
 real(kind=RKIND),intent(in) :: dt 

!local variables:
 type(block_type),pointer:: block
 integer:: i,j,k

!---------------------------------------------------------------------------------------------

 write(0,*)
 write(0,*) '--- begin physics_update:'
 block => domain % blocklist
 do while(associated(block))

    !parameterization of convection: update accumulated precipitation.
    !call update_convection_deep(dt,config_bucket_rainc,block%mesh,block%diag_physics)

    block => block % next
 end do
 write(0,*) '--- end physics_update:'

 end subroutine physics_update

!=============================================================================================
 subroutine update_radiation_diagnostics(bucket_radt,mesh,diag)
!=============================================================================================

!input arguments:
 real(kind=RKIND),intent(in):: bucket_radt
 type(mesh_type),intent(in):: mesh

!inout arguments:
 type(diag_physics_type),intent(inout):: diag

!local variables:
 integer:: iCell

!--------------------------------------------------------------------------------------------

 do iCell = 1, mesh%nCellsSolve
    !short-wave radiation:
    diag%acswdnb %array(iCell) = diag%acswdnb %array(iCell) + diag%swdnb %array(iCell)*dt_dyn
    diag%acswdnbc%array(iCell) = diag%acswdnbc%array(iCell) + diag%swdnbc%array(iCell)*dt_dyn
    diag%acswdnt %array(iCell) = diag%acswdnt %array(iCell) + diag%swdnt %array(iCell)*dt_dyn
    diag%acswdntc%array(iCell) = diag%acswdntc%array(iCell) + diag%swdntc%array(iCell)*dt_dyn
    diag%acswupb %array(iCell) = diag%acswupb %array(iCell) + diag%swupb %array(iCell)*dt_dyn
    diag%acswupbc%array(iCell) = diag%acswupbc%array(iCell) + diag%swupbc%array(iCell)*dt_dyn
    diag%acswupt %array(iCell) = diag%acswupt %array(iCell) + diag%swupt %array(iCell)*dt_dyn
    diag%acswuptc%array(iCell) = diag%acswuptc%array(iCell) + diag%swuptc%array(iCell)*dt_dyn
    !long-wave radiation:
    diag%aclwdnb %array(iCell) = diag%aclwdnb %array(iCell) + diag%lwdnb %array(iCell)*dt_dyn
    diag%aclwdnbc%array(iCell) = diag%aclwdnbc%array(iCell) + diag%lwdnbc%array(iCell)*dt_dyn
    diag%aclwdnt %array(iCell) = diag%aclwdnt %array(iCell) + diag%lwdnt %array(iCell)*dt_dyn
    diag%aclwdntc%array(iCell) = diag%aclwdntc%array(iCell) + diag%lwdntc%array(iCell)*dt_dyn
    diag%aclwupb %array(iCell) = diag%aclwupb %array(iCell) + diag%lwupb %array(iCell)*dt_dyn
    diag%aclwupbc%array(iCell) = diag%aclwupbc%array(iCell) + diag%lwupbc%array(iCell)*dt_dyn
    diag%aclwupt %array(iCell) = diag%aclwupt %array(iCell) + diag%lwupt %array(iCell)*dt_dyn
    diag%aclwuptc%array(iCell) = diag%aclwuptc%array(iCell) + diag%lwuptc%array(iCell)*dt_dyn
 enddo

 if(l_acradt .and. bucket_radt.gt.0._RKIND) then

    do iCell = 1, mesh%nCellsSolve
       !short-wave radiation:
       if(diag%acswdnb%array(iCell) .gt. bucket_radt) then
          diag%i_acswdnb%array(iCell) = diag%i_acswdnb%array(iCell) + 1
          diag%acswdnb%array(iCell) = diag%acswdnb%array(iCell) - bucket_radt
       endif   
       if(diag%acswdnbc%array(iCell) .gt. bucket_radt) then
          diag%i_acswdnbc%array(iCell) = diag%i_acswdnbc%array(iCell) + 1
          diag%acswdnbc%array(iCell) = diag%acswdnbc%array(iCell) - bucket_radt
       endif   
       if(diag%acswdnt%array(iCell) .gt. bucket_radt) then
          diag%i_acswdnt%array(iCell) = diag%i_acswdnt%array(iCell) + 1
          diag%acswdnt%array(iCell) = diag%acswdnt%array(iCell) - bucket_radt
       endif   
       if(diag%acswdntc%array(iCell) .gt. bucket_radt) then
          diag%i_acswdntc%array(iCell) = diag%i_acswdntc%array(iCell) + 1
          diag%acswdntc%array(iCell) = diag%acswdntc%array(iCell) - bucket_radt
       endif
       if(diag%acswupb%array(iCell) .gt. bucket_radt) then
          diag%i_acswupb%array(iCell) = diag%i_acswupb%array(iCell) + 1
          diag%acswupb%array(iCell) = diag%acswupb%array(iCell) - bucket_radt
       endif   
       if(diag%acswupbc%array(iCell) .gt. bucket_radt) then
          diag%i_acswupbc%array(iCell) = diag%i_acswupbc%array(iCell) + 1
          diag%acswupbc%array(iCell) = diag%acswupbc%array(iCell) - bucket_radt
       endif   
       if(diag%acswupt%array(iCell) .gt. bucket_radt) then
          diag%i_acswupt%array(iCell) = diag%i_acswupt%array(iCell) + 1
          diag%acswupt%array(iCell) = diag%acswupt%array(iCell) - bucket_radt
       endif   
       if(diag%acswuptc%array(iCell) .gt. bucket_radt) then
          diag%i_acswuptc%array(iCell) = diag%i_acswuptc%array(iCell) + 1
          diag%acswuptc%array(iCell) = diag%acswuptc%array(iCell) - bucket_radt
       endif
       !long-wave radiation:
       if(diag%aclwdnb%array(iCell) .gt. bucket_radt) then
          diag%i_aclwdnb%array(iCell) = diag%i_aclwdnb%array(iCell) + 1
          diag%aclwdnb%array(iCell) = diag%aclwdnb%array(iCell) - bucket_radt
       endif   
       if(diag%aclwdnbc%array(iCell) .gt. bucket_radt) then
          diag%i_aclwdnbc%array(iCell) = diag%i_aclwdnbc%array(iCell) + 1
          diag%aclwdnbc%array(iCell) = diag%aclwdnbc%array(iCell) - bucket_radt
       endif   
       if(diag%aclwdnt%array(iCell) .gt. bucket_radt) then
          diag%i_aclwdnt%array(iCell) = diag%i_aclwdnt%array(iCell) + 1
          diag%aclwdnt%array(iCell) = diag%aclwdnt%array(iCell) - bucket_radt
       endif   
       if(diag%aclwdntc%array(iCell) .gt. bucket_radt) then
          diag%i_aclwdntc%array(iCell) = diag%i_aclwdntc%array(iCell) + 1
          diag%aclwdntc%array(iCell) = diag%aclwdntc%array(iCell) - bucket_radt
       endif
       if(diag%aclwupb%array(iCell) .gt. bucket_radt) then
          diag%i_aclwupb%array(iCell) = diag%i_aclwupb%array(iCell) + 1
          diag%aclwupb%array(iCell) = diag%aclwupb%array(iCell) - bucket_radt
       endif   
       if(diag%aclwupbc%array(iCell) .gt. bucket_radt) then
          diag%i_aclwupbc%array(iCell) = diag%i_aclwupbc%array(iCell) + 1
          diag%aclwupbc%array(iCell) = diag%aclwupbc%array(iCell) - bucket_radt
       endif   
       if(diag%aclwupt%array(iCell) .gt. bucket_radt) then
          diag%i_aclwupt%array(iCell) = diag%i_aclwupt%array(iCell) + 1
          diag%aclwupt%array(iCell) = diag%aclwupt%array(iCell) - bucket_radt
       endif   
       if(diag%aclwuptc%array(iCell) .gt. bucket_radt) then
          diag%i_aclwuptc%array(iCell) = diag%i_aclwuptc%array(iCell) + 1
          diag%aclwuptc%array(iCell) = diag%aclwuptc%array(iCell) - bucket_radt
       endif
    enddo

 endif

 end subroutine update_radiation_diagnostics

!=============================================================================================
 end module mpas_atmphys_update
!=============================================================================================
