module  maml_module
!----------------------------------------------------------------------
! This module contains the MAML specific subroutines
!----------------------------------------------------------------------

use seq_comm_mct, only : num_inst_atm
use shr_kind_mod,  only: r8 => shr_kind_r8
!----------------------------------------------------------------------- 
! PRIVATE: Make default data and interfaces private
!----------------------------------------------------------------------- 
  private     ! By default all data is private to this module
!
! Public interfaces
!
  public cam_in_avg_mi_all ! Method to average surface interface fields(cam_in)  across instances, do the operation over all X_mi variables
  public cam_in_avg_mi     ! Method to average surface interface fields(cam_in)  across instances, do it for the specific X_mi variable
!  public cam_out_avg_mi    ! Method to average surface interface fields(cam_out) across instances
  public cam_out_rad_avg_mi    ! Method to average surface interface fields(cam_out) across instances
! public cam_in_copy_mi    
CONTAINS

subroutine cam_in_avg_mi_all(cam_in)

!----------------------------------------------------------------------- 
! Purpose: to make CRM domain average (across the multiple instances)
! for atmosphere model input. This subroutine needs to be called right before the subroutine
! surf_diag() inside the crm/physpkg.F90
!----------------------------------------------------------------------- 
   use camsrfexch,   only: cam_in_t
!-----------------------------------------------------------------------
!
! Input arguments
!
   type(cam_in_t),  intent(inout) :: cam_in

!
!---------------------------Local workspace-----------------------------
!
   real(r8) avgfac  
   integer  i, j 
!-----------------------------------------------------------------------
!
    ncol  = cam_in%ncol
    avgfac = 1._r8/real(num_inst_atm,r8)
    ! reinitialize *_avg, just in case
    cam_in%lhf = 0.
    cam_in%shf = 0.
    cam_in%wsx = 0.
    cam_in%wsy = 0.
    cam_in%snowhland = 0.
    cam_in%asdir = 0.
    cam_in%aldir = 0.
    cam_in%asdif = 0.
    cam_in%aldif = 0.
    cam_in%cflx(:,1) = 0.

    ! for non-MAML, num_inst_atm = avgfac = 1 
    do i = 1,ncol
      do j = 1,num_inst_atm
          cam_in%lhf(i)    = cam_in%lhf(i)+cam_in%lhf_mi(i,j)*avgfac
          cam_in%cflx(i,1) = cam_in%cflx(i,1)+cam_in%cflx1_mi(i,j)*avgfac
          cam_in%shf(i)    = cam_in%shf(i)+cam_in%shf_mi(i,j)*avgfac
          cam_in%wsx(i)    = cam_in%wsx(i)+cam_in%wsx_mi(i,j)*avgfac
          cam_in%wsy(i)    = cam_in%wsy(i)+cam_in%wsy_mi(i,j)*avgfac
          cam_in%snowhland(i) = cam_in%snowhland(i)+cam_in%snowhland_mi(i,j)*avgfac
          cam_in%asdir(i) = cam_in%asdir(i)+cam_in%asdir_mi(i,j)*avgfac
          cam_in%aldir(i) = cam_in%aldir(i)+cam_in%aldir_mi(i,j)*avgfac
          cam_in%asdif(i) = cam_in%asdif(i)+cam_in%asdif_mi(i,j)*avgfac
          cam_in%aldif(i) = cam_in%aldif(i)+cam_in%aldif_mi(i,j)*avgfac
      end do
   end do! i = 1, ncol
end subroutine cam_in_avg_mi_all


subroutine cam_in_avg_mi(cam_in_X, cam_in_X_mi, ncol)

!----------------------------------------------------------------------- 
! Purpose: to make CRM domain average (across the multiple instances)
! for atmosphere model input. This subroutine needs to be called right before the subroutine
! surf_diag() inside the crm/physpkg.F90
!----------------------------------------------------------------------- 
!
! Input arguments
!
   integer,   intent(in)     :: ncol
   real(r8),  intent(inout) :: cam_in_X(ncol)
   real(r8),  intent(inout) :: cam_in_X_mi(ncol, num_inst_atm)

!
!---------------------------Local workspace-----------------------------
!
   real(r8) avgfac  
   integer  i, j 
!-----------------------------------------------------------------------
!
    avgfac = 1._r8/real(num_inst_atm,r8)
    ! reinitialize *_avg, just in case
    cam_in_X = 0.

    ! for non-MAML, num_inst_atm = avgfac = 1 
    do i = 1,ncol
      do j = 1,num_inst_atm
          cam_in_X(i)    = cam_in_X(i)+cam_in_X_mi(i,j)*avgfac
      end do
   end do! i = 1, ncol
end subroutine cam_in_avg_mi

!subroutine cam_in_copy_mi(cam_in,ichnk,overwrite)

!----------------------------------------------------------------------- 
! Purpose: This is used for MMF configuration without MAML
! This subroutine copies cam_in%[X] to cam_in%[X]_mi after the data is transfered from surface through the coupler
! This is done because crm_input reads in cam_in%[X]_mi. For non-MAML, num_inst_atm = 1. 
!----------------------------------------------------------------------- 
!   use camsrfexch,   only: cam_in_t
!   use phys_grid ,   only: get_ncols_p
!-----------------------------------------------------------------------
!
! Input arguments
!
!   logical,         intent(in)    :: overwrite
!   integer,         intent(in)    :: ichnk
!   type(cam_in_t),  intent(inout) :: cam_in

!
!  Local Space
!
!   integer :: ncol
!   ncol = cam_in%ncol

   ! for non-MAML, num_inst_atm =  1 
!   if(overwrite) then
!     ! first timestep of the restart run, lhf, shf and cflx(:,1) are not overwritten 
!     cam_in(ichnk)%lhf_mi(1:ncol,1) = cam_in(ichnk)%lhf(1:ncol)
!     cam_in(ichnk)%shf_mi(1:ncol,1) = cam_in(ichnk)%shf(1:ncol)
!     cam_in(ichnk)%cflx_mi(1:ncol,1,1) = cam_in(ichnk)%cflx(1:ncol,1)
!   end if  
!     cam_in(ichnk)%wsx_mi(1:ncol,1) = cam_in(ichnk)%wsx(1:ncol)
!     cam_in(ichnk)%wsy_mi(1:ncol,1) = cam_in(ichnk)%wsy(1:ncol)
!     cam_in(ichnk)%snowhland_mi(1:ncol,1) = cam_in(ichnk)%snowhland(1:ncol)
!     cam_in(ichnk)%asdir_mi(1:ncol,1) = cam_in(ichnk)%asdir(1:ncol)
!     cam_in(ichnk)%aldir_mi(1:ncol,1) = cam_in(ichnk)%aldir(1:ncol)
!     cam_in(ichnk)%asdif_mi(1:ncol,1) = cam_in(ichnk)%asdif(1:ncol)
!     cam_in(ichnk)%aldif_mi(1:ncol,1) = cam_in(ichnk)%aldif(1:ncol)

!end subroutine cam_in_copy_mi

subroutine cam_out_avg_mi(cam_out)

!----------------------------------------------------------------------- 
! Purpose: to make CRM domain average (across the multiple instances)
! for cam output. This subroutine needs to be called at the end of cam_export()
!----------------------------------------------------------------------- 
   use camsrfexch,   only: cam_out_t
!-----------------------------------------------------------------------
!
! Input arguments
!
   type(cam_out_t),  intent(inout) :: cam_out

!
!---------------------------Local workspace-----------------------------
!
   real(r8) avgfac  
   integer  i, j 
!-----------------------------------------------------------------------
!
    ncol  = cam_out%ncol
    avgfac = 1._r8/real(num_inst_atm,r8)
    ! Initialize again.
    cam_out%tbot = 0.
    cam_out%precsc = 0.
    cam_out%precsl = 0.
    cam_out%precc = 0.
    cam_out%precl = 0.

    ! for non-MAML, num_inst_atm = avgfac = 1 
    do i = 1,ncol
      do j = 1,num_inst_atm
          cam_out%tbot(i)  = cam_out%tbot(i)+cam_out%tbot_mi(i,j)*avgfac
          cam_out%precsc(i) = cam_out%precsc(i)+cam_out%precsc_mi(i,j)*avgfac
          cam_out%precsl(i)= cam_out%precsl(i)+cam_out%precsl_mi(i,j)*avgfac
          cam_out%precc(i) = cam_out%precc(i)+cam_out%precc_mi(i,j)*avgfac
          cam_out%precl(i) = cam_out%precl(i)+cam_out%precl_mi(i,j)*avgfac
      end do
   end do! i = 1, ncol
end subroutine cam_out_avg_mi

subroutine cam_out_rad_avg_mi(cam_out)

!----------------------------------------------------------------------- 
! Purpose: to make CRM domain average (across the multiple instances)
! for cam output. This subroutine needs to be called at the end of cam_export()
!----------------------------------------------------------------------- 
   use camsrfexch,   only: cam_out_t
!-----------------------------------------------------------------------
!
! Input arguments
!
   type(cam_out_t),  intent(inout) :: cam_out

!
!---------------------------Local workspace-----------------------------
!
   real(r8) avgfac  
   integer  i, j 
!-----------------------------------------------------------------------
!
    ncol  = cam_out%ncol
    avgfac = 1._r8/real(num_inst_atm,r8)
    ! Initialize again.
    cam_out%sols = 0.
    cam_out%soll = 0.
    cam_out%solld = 0.
    cam_out%solsd = 0.
    cam_out%flwds = 0.

    ! for non-MAML, num_inst_atm = avgfac = 1 
    do i = 1,ncol
      do j = 1,num_inst_atm
          cam_out%sols(i)  = cam_out%sols(i)+cam_out%sols_mi(i,j)*avgfac
          cam_out%soll(i) = cam_out%soll(i)+cam_out%soll_mi(i,j)*avgfac
          cam_out%solld(i)= cam_out%solld(i)+cam_out%solld_mi(i,j)*avgfac
          cam_out%solsd(i) = cam_out%solsd(i)+cam_out%solsd_mi(i,j)*avgfac
          cam_out%flwds(i) = cam_out%flwds(i)+cam_out%flwds_mi(i,j)*avgfac
      end do
   end do! i = 1, ncol
end subroutine cam_out_rad_avg_mi
end module maml_module
