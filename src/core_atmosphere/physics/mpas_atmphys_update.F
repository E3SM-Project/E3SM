! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================
 module mpas_atmphys_update
 use mpas_configure
 use mpas_grid_types

 use mpas_atmphys_driver_convection_deep
 use mpas_atmphys_vars

 implicit none
 private
 public:: physics_update, &
          update_radiation_diagnostics


!>\brief update diagnostics.
!>\author Laura D. Fowler (send comments to laura@ucar.edu).
!>\date 2013-05-01.
!>
!>\details
!>
!> subroutines in mpas_atmphys_update:
!> -----------------------------------
!> physics_update              : not used.
!> update_radiation_diagnostics: update accumulated radiation diagnostics.
!>
!> add-ons and modifications to sourcecode:
!> ----------------------------------------
!>    * throughout the sourcecode, replaced all "var_struct" defined arrays by local pointers.
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2014-04-22.


 contains
 
!==================================================================================================
 subroutine physics_update(domain,dt)
!==================================================================================================

!input arguments:
 type(domain_type),intent(in):: domain
 real(kind=RKIND),intent(in) :: dt 

!local variables:
 type(block_type),pointer:: block
 integer:: i,j,k

!--------------------------------------------------------------------------------------------------

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

!==================================================================================================
 subroutine update_radiation_diagnostics(bucket_radt,mesh,diag_physics)
!==================================================================================================

!input arguments:
 real(kind=RKIND),intent(in):: bucket_radt
 type(mesh_type),intent(in):: mesh

!inout arguments:
 type(diag_physics_type),intent(inout):: diag_physics

!local variables and arrays:
 integer:: iCell,nCellsSolve
 integer,dimension(:),pointer  :: i_acswdnb,i_acswdnbc,i_acswdnt,i_acswdntc, &
                                  i_acswupb,i_acswupbc,i_acswupt,i_acswuptc, &
                                  i_aclwdnb,i_aclwdnbc,i_aclwdnt,i_aclwdntc, &
                                  i_aclwupb,i_aclwupbc,i_aclwupt,i_aclwuptc

 real(kind=RKIND),dimension(:),pointer:: swdnb,swdnbc,swdnt,swdntc,  &
                                         swupb,swupbc,swupt,swuptc,  &
                                         lwdnb,lwdnbc,lwdnt,lwdntc,  &
                                         lwupb,lwupbc,lwupt,lwuptc
 real(kind=RKIND),dimension(:),pointer:: acswdnb,acswdnbc,acswdnt,acswdntc,  &
                                         acswupb,acswupbc,acswupt,acswuptc,  &
                                         aclwdnb,aclwdnbc,aclwdnt,aclwdntc,  &
                                         aclwupb,aclwupbc,aclwupt,aclwuptc

!--------------------------------------------------------------------------------------------------

 nCellsSolve = mesh % nCellsSolve

 i_acswdnb    => diag_physics % i_acswdnb  % array
 i_acswdnbc   => diag_physics % i_acswdnbc % array
 i_acswdnt    => diag_physics % i_acswdnt  % array
 i_acswdntc   => diag_physics % i_acswdntc % array
 i_acswupb    => diag_physics % i_acswupb  % array
 i_acswupbc   => diag_physics % i_acswupbc % array
 i_acswupt    => diag_physics % i_acswupt  % array
 i_acswuptc   => diag_physics % i_acswuptc % array
 i_aclwdnb    => diag_physics % i_aclwdnb  % array
 i_aclwdnbc   => diag_physics % i_aclwdnbc % array
 i_aclwdnt    => diag_physics % i_aclwdnt  % array
 i_aclwdntc   => diag_physics % i_aclwdntc % array
 i_aclwupb    => diag_physics % i_aclwupb  % array
 i_aclwupbc   => diag_physics % i_aclwupbc % array
 i_aclwupt    => diag_physics % i_aclwupt  % array
 i_aclwuptc   => diag_physics % i_aclwuptc % array

 acswdnb      => diag_physics % acswdnb    % array
 acswdnbc     => diag_physics % acswdnbc   % array
 acswdnt      => diag_physics % acswdnt    % array
 acswdntc     => diag_physics % acswdntc   % array
 acswupb      => diag_physics % acswupb    % array
 acswupbc     => diag_physics % acswupbc   % array
 acswupt      => diag_physics % acswupt    % array
 acswuptc     => diag_physics % acswuptc   % array
 aclwdnb      => diag_physics % aclwdnb    % array
 aclwdnbc     => diag_physics % aclwdnbc   % array
 aclwdnt      => diag_physics % aclwdnt    % array
 aclwdntc     => diag_physics % aclwdntc   % array
 aclwupb      => diag_physics % aclwupb    % array
 aclwupbc     => diag_physics % aclwupbc   % array
 aclwupt      => diag_physics % aclwupt    % array
 aclwuptc     => diag_physics % aclwuptc   % array

 swdnb        => diag_physics % swdnb      % array
 swdnbc       => diag_physics % swdnbc     % array
 swdnt        => diag_physics % swdnt      % array
 swdntc       => diag_physics % swdntc     % array
 swupb        => diag_physics % swupb      % array
 swupbc       => diag_physics % swupbc     % array
 swupt        => diag_physics % swupt      % array
 swuptc       => diag_physics % swuptc     % array
 lwdnb        => diag_physics % lwdnb      % array
 lwdnbc       => diag_physics % lwdnbc     % array
 lwdnt        => diag_physics % lwdnt      % array
 lwdntc       => diag_physics % lwdntc     % array
 lwupb        => diag_physics % lwupb      % array
 lwupbc       => diag_physics % lwupbc     % array
 lwupt        => diag_physics % lwupt      % array
 lwuptc       => diag_physics % lwuptc     % array

 do iCell = 1, nCellsSolve
    !short-wave radiation:
    acswdnb(iCell)  = acswdnb (iCell) + swdnb (iCell)*dt_dyn
    acswdnbc(iCell) = acswdnbc(iCell) + swdnbc(iCell)*dt_dyn
    acswdnt(iCell)  = acswdnt (iCell) + swdnt (iCell)*dt_dyn
    acswdntc(iCell) = acswdntc(iCell) + swdntc(iCell)*dt_dyn
    acswupb(iCell)  = acswupb (iCell) + swupb (iCell)*dt_dyn
    acswupbc(iCell) = acswupbc(iCell) + swupbc(iCell)*dt_dyn
    acswupt(iCell)  = acswupt (iCell) + swupt (iCell)*dt_dyn
    acswuptc(iCell) = acswuptc(iCell) + swuptc(iCell)*dt_dyn
    !long-wave radiation:
    aclwdnb(iCell)  = aclwdnb (iCell) + lwdnb (iCell)*dt_dyn
    aclwdnbc(iCell) = aclwdnbc(iCell) + lwdnbc(iCell)*dt_dyn
    aclwdnt(iCell)  = aclwdnt (iCell) + lwdnt (iCell)*dt_dyn
    aclwdntc(iCell) = aclwdntc(iCell) + lwdntc(iCell)*dt_dyn
    aclwupb(iCell)  = aclwupb (iCell) + lwupb (iCell)*dt_dyn
    aclwupbc(iCell) = aclwupbc(iCell) + lwupbc(iCell)*dt_dyn
    aclwupt(iCell)  = aclwupt (iCell) + lwupt (iCell)*dt_dyn
    aclwuptc(iCell) = aclwuptc(iCell) + lwuptc(iCell)*dt_dyn
 enddo

 if(l_acradt .and. bucket_radt.gt.0._RKIND) then

    do iCell = 1, nCellsSolve
       !short-wave radiation:
       if(acswdnb(iCell) .gt. bucket_radt) then
          i_acswdnb(iCell) = i_acswdnb(iCell) + 1
          acswdnb(iCell) = acswdnb(iCell) - bucket_radt
       endif   
       if(acswdnbc(iCell) .gt. bucket_radt) then
          i_acswdnbc(iCell) = i_acswdnbc(iCell) + 1
          acswdnbc(iCell) = acswdnbc(iCell) - bucket_radt
       endif   
       if(acswdnt(iCell) .gt. bucket_radt) then
          i_acswdnt(iCell) = i_acswdnt(iCell) + 1
          acswdnt(iCell) = acswdnt(iCell) - bucket_radt
       endif   
       if(acswdntc(iCell) .gt. bucket_radt) then
          i_acswdntc(iCell) = i_acswdntc(iCell) + 1
          acswdntc(iCell) = acswdntc(iCell) - bucket_radt
       endif
       if(acswupb(iCell) .gt. bucket_radt) then
          i_acswupb(iCell) = i_acswupb(iCell) + 1
          acswupb(iCell) = acswupb(iCell) - bucket_radt
       endif   
       if(acswupbc(iCell) .gt. bucket_radt) then
          i_acswupbc(iCell) = i_acswupbc(iCell) + 1
          acswupbc(iCell) = acswupbc(iCell) - bucket_radt
       endif   
       if(acswupt(iCell) .gt. bucket_radt) then
          i_acswupt(iCell) = i_acswupt(iCell) + 1
          acswupt(iCell) = acswupt(iCell) - bucket_radt
       endif   
       if(acswuptc(iCell) .gt. bucket_radt) then
          i_acswuptc(iCell) = i_acswuptc(iCell) + 1
          acswuptc(iCell) = acswuptc(iCell) - bucket_radt
       endif
       !long-wave radiation:
       if(aclwdnb(iCell) .gt. bucket_radt) then
          i_aclwdnb(iCell) = i_aclwdnb(iCell) + 1
          aclwdnb(iCell) = aclwdnb(iCell) - bucket_radt
       endif   
       if(aclwdnbc(iCell) .gt. bucket_radt) then
          i_aclwdnbc(iCell) = i_aclwdnbc(iCell) + 1
          aclwdnbc(iCell) = aclwdnbc(iCell) - bucket_radt
       endif   
       if(aclwdnt(iCell) .gt. bucket_radt) then
          i_aclwdnt(iCell) = i_aclwdnt(iCell) + 1
          aclwdnt(iCell) = aclwdnt(iCell) - bucket_radt
       endif   
       if(aclwdntc(iCell) .gt. bucket_radt) then
          i_aclwdntc(iCell) = i_aclwdntc(iCell) + 1
          aclwdntc(iCell) = aclwdntc(iCell) - bucket_radt
       endif
       if(aclwupb(iCell) .gt. bucket_radt) then
          i_aclwupb(iCell) = i_aclwupb(iCell) + 1
          aclwupb(iCell) = aclwupb(iCell) - bucket_radt
       endif   
       if(aclwupbc(iCell) .gt. bucket_radt) then
          i_aclwupbc(iCell) = i_aclwupbc(iCell) + 1
          aclwupbc(iCell) = aclwupbc(iCell) - bucket_radt
       endif   
       if(aclwupt(iCell) .gt. bucket_radt) then
          i_aclwupt(iCell) = i_aclwupt(iCell) + 1
          aclwupt(iCell) = aclwupt(iCell) - bucket_radt
       endif   
       if(aclwuptc(iCell) .gt. bucket_radt) then
          i_aclwuptc(iCell) = i_aclwuptc(iCell) + 1
          aclwuptc(iCell) = aclwuptc(iCell) - bucket_radt
       endif
    enddo

 endif

 end subroutine update_radiation_diagnostics

!==================================================================================================
 end module mpas_atmphys_update
!==================================================================================================
