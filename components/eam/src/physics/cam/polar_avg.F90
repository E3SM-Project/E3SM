module polar_avg
!----------------------------------------------------------------------- 
! 
! Purpose: 
!  These routines are used by the fv dycore to set the collocated 
!  pole points at the limits of the latitude dimension to the same 
!  value. 
!
! Methods: 
!  The reprosum reproducible distributed sum is used for these 
!  calculations.
!
! Author: A. Mirin
! 
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!- use statements ------------------------------------------------------
!-----------------------------------------------------------------------
   use shr_kind_mod,  only: r8 => shr_kind_r8
   use dycore,        only: dycore_is
   use dyn_grid,      only: get_dyn_grid_parm
   use phys_grid,     only: get_ncols_p, get_lat_all_p
   use ppgrid,        only: begchunk, endchunk, pcols
   use shr_reprosum_mod, only: shr_reprosum_calc
#if ( defined SPMD )
   use mpishorthand,  only: mpicom
#endif

!-----------------------------------------------------------------------
!- module boilerplate --------------------------------------------------
!-----------------------------------------------------------------------
   implicit none
   private
   save

!-----------------------------------------------------------------------
! Public interfaces ----------------------------------------------------
!-----------------------------------------------------------------------
   public :: &
      polar_average           ! support for LR dycore polar averaging

   interface polar_average
      module procedure polar_average2d, polar_average3d
   end interface

   CONTAINS
!
!========================================================================
!
   subroutine polar_average2d(field)
!----------------------------------------------------------------------- 
! Purpose: Set the collocated pole points at the limits of the latitude 
!          dimension to the same value. 
! Author: J. Edwards
!-----------------------------------------------------------------------
!
! Arguments
!
     real(r8), intent(inout) :: field(pcols,begchunk:endchunk)
!
! Local workspace
!
     integer :: i, c, ln, ls, ncols
     integer :: plat, plon
     integer, allocatable :: lats(:)
#if (! defined SPMD)
     integer  :: mpicom = 0
#endif

     real(r8) :: sum(2)
     real(r8), allocatable :: n_pole(:), s_pole(:)
!
!-----------------------------------------------------------------------
!
     if(.not. dycore_is('LR')) return

     plon = get_dyn_grid_parm('plon')
     plat = get_dyn_grid_parm('plat')
     allocate(lats(pcols), n_pole(plon), s_pole(plon))
     ln=0
     ls=0
     do c=begchunk,endchunk
        call get_lat_all_p(c,pcols,lats) 
	ncols = get_ncols_p(c)

	do i=1,ncols
           if(lats(i).eq.1) then
              ln=ln+1
              n_pole(ln) = field(i,c)
           else if(lats(i).eq.plat) then
              ls=ls+1
              s_pole(ls) = field(i,c)
           end if
	enddo
        
     end do
     
     call shr_reprosum_calc(n_pole, sum(1:1), ln, plon, 1, &
                    gbl_count=plon, commid=mpicom)

     call shr_reprosum_calc(s_pole, sum(2:2), ls, plon, 1, &
                    gbl_count=plon, commid=mpicom)

     ln=0
     ls=0
     do c=begchunk,endchunk
        call get_lat_all_p(c,pcols,lats) 
	ncols = get_ncols_p(c)

	do i=1,ncols
           if(lats(i).eq.1) then
              ln=ln+1
              field(i,c) = sum(1)/plon
           else if(lats(i).eq.plat) then
              ls=ls+1
              field(i,c) = sum(2)/plon
           end if
	enddo
        
     end do

     deallocate(lats, n_pole, s_pole)
   
   end subroutine polar_average2d

!
!========================================================================
!

   subroutine polar_average3d(nlev, field)
!----------------------------------------------------------------------- 
! Purpose: Set the collocated pole points at the limits of the latitude 
!          dimension to the same value. 
! Author: J. Edwards
!-----------------------------------------------------------------------
!
! Arguments
!
     integer, intent(in) :: nlev
     real(r8), intent(inout) :: field(pcols,nlev,begchunk:endchunk)
!
! Local workspace
!
     integer :: i, c, ln, ls, ncols, k
     integer :: plat, plon
     integer, allocatable :: lats(:)
#if (! defined SPMD)
     integer  :: mpicom = 0
#endif

     real(r8) :: sum(nlev,2)
     real(r8), allocatable :: n_pole(:,:), s_pole(:,:)
!
!-----------------------------------------------------------------------
!
     if(.not. dycore_is('LR')) return

     plon = get_dyn_grid_parm('plon')
     plat = get_dyn_grid_parm('plat')
     allocate(lats(pcols), n_pole(plon,nlev), s_pole(plon,nlev))
     ln=0
     ls=0
     do c=begchunk,endchunk
        call get_lat_all_p(c,pcols,lats) 
	ncols = get_ncols_p(c)

        do i=1,ncols
           if(lats(i).eq.1) then
              ln=ln+1
              do k=1,nlev
                 n_pole(ln,k) = field(i,k,c)
              end do
           else if(lats(i).eq.plat) then
              ls=ls+1
              do k=1,nlev                 
                 s_pole(ls,k) = field(i,k,c)
              end do
           end if
        enddo
     end do
     
     call shr_reprosum_calc(n_pole, sum(:,1), ln, plon, nlev, &
                    gbl_count=plon, commid=mpicom)

     call shr_reprosum_calc(s_pole, sum(:,2), ls, plon, nlev, &
                    gbl_count=plon, commid=mpicom)

     ln=0
     ls=0
     do c=begchunk,endchunk
        call get_lat_all_p(c,pcols,lats) 
	ncols = get_ncols_p(c)

	do i=1,ncols
           if(lats(i).eq.1) then
              ln=ln+1
              do k=1,nlev
                 field(i,k,c) = sum(k,1)/plon
              end do
           else if(lats(i).eq.plat) then
              ls=ls+1
              do k=1,nlev
                 field(i,k,c) = sum(k,2)/plon
              end do
           end if
	enddo
        
     end do
   
     deallocate(lats, n_pole, s_pole)

   end subroutine polar_average3d

end module polar_avg
