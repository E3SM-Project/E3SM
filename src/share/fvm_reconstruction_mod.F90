#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

!MODULE FVM_RECONSTRUCTION_MOD--------------------------------------------CE-for FVM!
! AUTHOR: CHRISTOPH ERATH, 17.October 2011                                          !
! This module contains everything  to do (ONLY) a CUBIC (3rd order) reconstruction  ! 
!                                                                                   !
! IMPORTANT: the implementation is done for a ncfl > 1, which is not working        !
!            but it works for ncfl=1                                                !
!
! phl: ibaseref goes out of range when using boundscheck compile option
!      have not checked if if-statements below changes answers yet!
!
!
!
!-----------------------------------------------------------------------------------!
module fvm_reconstruction_mod

  use kinds, only                  : int_kind, real_kind
  use dimensions_mod, only         : nc,nhc,nhe
  use coordinate_systems_mod, only : cartesian2D_t,cartesian3D_t
  use control_mod, only : north, south, east, west, neast, nwest, seast, swest
  use fvm_control_volume_mod, only: fvm_struct
  use parallel_mod, only: abortmp
  implicit none
  private
  public :: reconstruction
contains
! ----------------------------------------------------------------------------------!
!SUBROUTINE RECONSTRUCTION------------------------------------------------CE-for FVM!
! AUTHOR: CHRISTOPH ERATH, 17.October 2011                                          !
! DESCRIPTION: controls the cubic (3rd order) reconstructions:                      !
!                                                                                   !
! CALLS: fillhalo_cubic, reconstruction_cubic                                       !
! INPUT: fcube    ...  tracer values incl. the halo zone                            !
!        fvm    ...  structure incl. tracer values aso                            !                                   ! 
! OUTPUT:recons   ...  has the reconstruction coefficients (5) for the 3rd order    !
!                      reconstruction: dx, dy, dx^2, dy^2, dxdy                     !
!-----------------------------------------------------------------------------------!
subroutine reconstruction(fcube,fvm,recons)
  use fvm_control_volume_mod, only: fvm_struct

  implicit none
  real (kind=real_kind), dimension(1-nhc:nc+nhc, 1-nhc:nc+nhc), intent(in) :: fcube
  type (fvm_struct), intent(in)                                     :: fvm
  real (kind=real_kind), dimension(5,1-nhe:nc+nhe,1-nhe:nc+nhe), &
                                                    intent(out)       :: recons
  
  real (kind=real_kind),dimension(-1:nc+2,2,2)                        :: fnewval
  real (kind=real_kind), dimension(-1:nc+2,2,2)                       :: fhalo
  real (kind=real_kind), dimension(0:nhe+1,2,2)                       :: fhaloex

  ! if element lies on a cube edge, recalculation the values in the halo zone
  if (fvm%cubeboundary > 0) then    
    call fillhalo_cubic(fcube,fvm,fnewval,fhalo,fhaloex)               
  end if

  call reconstruction_cubic(fcube,fvm,fnewval,fhalo,&
                            fhaloex, recons)
end subroutine reconstruction
!END SUBROUTINE RECONSTRUCTION--------------------------------------------CE-for FVM!


! ----------------------------------------------------------------------------------!
!SUBROUTINE FILLHALO_CUBIC------------------------------------------------CE-for FVM!
! AUTHOR: CHRISTOPH ERATH, 17.October 2011                                          !
! DESCRIPTION: uses the interpolation points in the fvm structure,                !
!        on different cubic faces, also in the halo region:                         !
!        because we also need the reconstruction coefficients in the halo zone,     !
!        which is basically calculated twice, on the original cell of an element    !
!        on face A and on a cell in the halo region of an element of face B         !
!        The crux is, that the interpolation has to be the same to ensure           !
!        conservation of the scheme                                                 !
!        SYMMETRY of the CUBE is used for calucaltion the interpolation_point       !
!                                                                                   !
! CALLS: interpolation_point, interpolation_value                                   !
! INPUT: fcube    ...  tracer values incl. the halo zone                            !
!        fvm    ...  fvm structure                                              !
! OUTPUT:fnewval ... contains the interpolationvalues for cell on the actual element!
!        fhalo ... contains interpolation values for the reconstruction in the halo !
!                  zone                                                             !
!        fhaloex...contains interpolation values for the reconstruction in the halo !
!                  zone, only for corner elements                                   !                  
!-----------------------------------------------------------------------------------!
subroutine fillhalo_cubic(fcube,fvm,fnewval,fhalo,fhaloex)

  use parallel_mod, only : haltmp
  use cube_mod, only     : cube_xstart, cube_xend, cube_ystart, cube_yend
  
  implicit none
  real (kind=real_kind), dimension(1-nhc:nc+nhc, 1-nhc:nc+nhc), intent(in) :: fcube
  type (fvm_struct), intent(in)                                          :: fvm
  
  real (kind=real_kind),  dimension(-1:nc+2,2,2), intent(out)        :: fnewval
  real (kind=real_kind),  dimension(-1:nc+2,2,2), intent(out)        :: fhalo
  real (kind=real_kind),  dimension(0:nhe+1,2,2), intent(out)        :: fhaloex
  
  
  integer                                       :: i, halo, ibaseref
  real (kind=real_kind), dimension(1:nhc)       :: tmpfval
  
  ! element is not on a corner, but shares a cube edge (call of subroutine)
!  write(*,*) "this should be fixed"
!  write(*,*) "arrays are out of bounds"
!  write(*,*) "Fortran runtime error: Index '-4' of dimension 2 of array 'fcube' outside of expected range (-3:8)"
!  write(*,*) "currently code does not fill halo"
!  return
  select case (fvm%cubeboundary)
    ! NOTE case(0) is excluded in the subroutine call
!-----------------------------------------------------------------------------------!
    !CASE WEST
    case(west) 
      ! west zone
      do halo=1,2
        do i=-1,nc+2
          ibaseref=fvm%ibase(i,halo,1)                                      
          fnewval(i,halo,1)=cubic_equispace_interp(fvm%dbeta, fvm%interp(i,halo,1),&
                                                  fcube(1-halo,ibaseref:ibaseref+3))                                      
        end do
      end do
      !westhalo: east zone (because of the reconstruction in the halo, conservation!)
      do halo=1,2
        do i=-1,nc+2
          ibaseref=fvm%ibasehalo(i,halo,1)                                      
          fhalo(i,halo,1)=cubic_equispace_interp(fvm%dbeta, fvm%interphalo(i,halo,1),&
                                                  fcube(halo,ibaseref:ibaseref+3))                                    
        end do
      end do
!-----------------------------------------------------------------------------------!
    !CASE EAST  
    case(east)
      ! east zone
      do halo=1,2
        do i=-1,nc+2
          ibaseref=fvm%ibase(i,halo,1)                                      
          fnewval(i,halo,1)=cubic_equispace_interp(fvm%dbeta, fvm%interp(i,halo,1),&
                                                  fcube(nc+halo,ibaseref:ibaseref+3))                                      
        end do 
      end do  
      !easthalo: west zone (because of the reconstruction in the halo, conservation!)
      do halo=1,2
        do i=-1,nc+2                                    
          ibaseref=fvm%ibasehalo(i,halo,1)                                      
          fhalo(i,halo,1)=cubic_equispace_interp(fvm%dbeta, fvm%interphalo(i,halo,1),&
                                                  fcube(nc+1-halo,ibaseref:ibaseref+3))                                    
        end do
      end do
!-----------------------------------------------------------------------------------!  
    !CASE NORTH  
    case(north)
      ! north zone
      do halo=1,2 
        do i=-1,nc+2
          ibaseref=fvm%ibase(i,halo,2)                                      
!          if (ibaseref<1-nhc) write(*,*) "ibaseref is out of range 1"
          fnewval(i,halo,2)=cubic_equispace_interp(fvm%dalpha, fvm%interp(i,halo,2),&
               fcube(ibaseref:ibaseref+3,nc+halo))  
        end do
      end do
      !northhalo: south zone (because of the reconstruction in the halo, conservation!)
      do halo=1,2
        do i=-1,nc+2
          ibaseref=fvm%ibasehalo(i,halo,2)                         
!          if (ibaseref<1-nhc) write(*,*) "ibaseref is out of range 2"             
          fhalo(i,halo,2)=cubic_equispace_interp(fvm%dalpha, fvm%interphalo(i,halo,2),&
                                                   fcube(ibaseref:ibaseref+3,nc+1-halo))
        end do
      end do   
!-----------------------------------------------------------------------------------!         
    !CASE SOUTH   
    case(south)
     !south zone
     do halo=1,2
        do i=-1,nc+2
          ibaseref=fvm%ibase(i,halo,2)                                      
!          if (ibaseref<1-nhc) write(*,*) "ibaseref is out of range 3"
          fnewval(i,halo,2)=cubic_equispace_interp(fvm%dalpha, fvm%interp(i,halo,2),&
                                                   fcube(ibaseref:ibaseref+3,1-halo))
        end do
      end do      

      !southhalo: north zone (because of the reconstruction in the halo, conservation!) 
      do halo=1,2
        do i=-1,nc+2
          ibaseref=fvm%ibasehalo(i,halo,2)                                      
!          if (ibaseref<1-nhc) write(*,*) "ibaseref is out of range 4"
          fhalo(i,halo,2)=cubic_equispace_interp(fvm%dalpha, fvm%interphalo(i,halo,2),&
                                                   fcube(ibaseref:ibaseref+3,halo))
        end do
      end do
!-----------------------------------------------------------------------------------!      
!CORNER TREATMENT
    !CASE SOUTH WEST
    case(swest) 
      ! west zone
      do halo=1,2
        do i=0,nc+2
          ibaseref=fvm%ibase(i,halo,1)                                      
          fnewval(i,halo,1)=cubic_equispace_interp(fvm%dbeta, fvm%interp(i,halo,1),&
                                                   fcube(1-halo,ibaseref:ibaseref+3))
        end do
      ! south zone
        do i=0,nc+2
          ibaseref=fvm%ibase(i,halo,2)                          
!          if (ibaseref<1-nhc) write(*,*) "ibaseref is out of range 5"            
          fnewval(i,halo,2)=cubic_equispace_interp(fvm%dalpha, fvm%interp(i,halo,2),&
                                                   fcube(ibaseref:ibaseref+3,1-halo))
        end do
      end do
      !westhalo: east zone & south zone
      !(because of the reconstruction in the halo, conservation!) 
      do halo=1,2
        do i=0,nc+2
          ibaseref=fvm%ibasehalo(i,halo,1)                                      
          fhalo(i,halo,1)=cubic_equispace_interp(fvm%dbeta, fvm%interphalo(i,halo,1),&
                                                   fcube(halo,ibaseref:ibaseref+3))
        end do
        do i=1,nhc
          tmpfval(i)=fcube(halo,-nhc+i)
        enddo
        do i=0,nhe+1
          ibaseref=fvm%ibasehaloex(i,halo,1)                                      
          fhaloex(nhe+1-i,halo,1)=cubic_equispace_interp(fvm%dalpha, fvm%interphaloex(i,halo,1),&
                                                         tmpfval)
        end do
      end do
      !southhalo: north zone & west zone
      !(because of the reconstruction in the halo, conservation!) 
      do halo=1,2
        do i=0,nc+2
          ibaseref=fvm%ibasehalo(i,halo,2)                                      
!          if (ibaseref<1-nhc) write(*,*) "ibaseref is out of range 6"
          fhalo(i,halo,2)=cubic_equispace_interp(fvm%dalpha, fvm%interphalo(i,halo,2),&
                                                   fcube(ibaseref:ibaseref+3,halo))
        end do
        do i=1,nhc
          tmpfval(i)=fcube(i-nhc,halo)
        enddo
        do i=0,nhe+1
          ibaseref=fvm%ibasehaloex(i,halo,2)                                      
          fhaloex(nhe+1-i,halo,2)=cubic_equispace_interp(fvm%dbeta, fvm%interphaloex(i,halo,2),&
                                                   tmpfval)
        end do
      end do
!-----------------------------------------------------------------------------------!      
    !CASE SOUTH EAST  
    case(seast)
      ! east zone
      do halo=1,2
        do i=0,nc+2
          ibaseref=fvm%ibase(i,halo,1)                                      
          fnewval(i,halo,1)=cubic_equispace_interp(fvm%dbeta, fvm%interp(i,halo,1),&
                                                   fcube(nc+halo,ibaseref:ibaseref+3))
        end do
      ! south zone
        do i=-1,nc+1
          ibaseref=fvm%ibase(i,halo,2)                             
          if (ibaseref<1-nhc) then 
!             write(*,*) "ibaseref is out of range 7"         
             ibaseref=1-nhc
          end if
          fnewval(i,halo,2)=cubic_equispace_interp(fvm%dalpha, fvm%interp(i,halo,2),&
                                                   fcube(ibaseref:ibaseref+3,1-halo))
        end do
      end do      
      !easthalo: west zone & south zone
      !(because of the reconstruction in the halo, conservation!) 
      do halo=1,2
        do i=0,nc+2
          ibaseref=fvm%ibasehalo(i,halo,1)                                      
          fhalo(i,halo,1)=cubic_equispace_interp(fvm%dbeta, fvm%interphalo(i,halo,1),&
                                                   fcube(nc+1-halo,ibaseref:ibaseref+3))
        end do
        do i=1,nhc
          tmpfval(i)=fcube(nc+1-halo,1-i)
        enddo
        do i=0,nhe+1
          ibaseref=fvm%ibasehaloex(i,halo,1)                                      
          fhaloex(i,halo,1)=cubic_equispace_interp(fvm%dalpha, fvm%interphaloex(i,halo,1),&
                                                   tmpfval)
        end do
      end do
      !(because of the reconstruction in the halo, conservation!) 
      !southhalo: north zone & west zone
      do halo=1,2
        do i=-1,nc+1
          ibaseref=fvm%ibasehalo(i,halo,2)                         
          if (ibaseref<1-nhc) then
!             write(*,*) "ibaseref is out of range 8"             
             ibaseref=1-nhc
          end if
          fhalo(i,halo,2)=cubic_equispace_interp(fvm%dalpha, fvm%interphalo(i,halo,2),&
                                                   fcube(ibaseref:ibaseref+3,halo))
        end do
        do i=1,nhc
          tmpfval(i)=fcube(nc+nhc+1-i,halo)
        enddo
        do i=0,nhe+1
          ibaseref=fvm%ibasehaloex(i,halo,2)                                      
          fhaloex(nhe+1-i,halo,2)=cubic_equispace_interp(fvm%dbeta, fvm%interphaloex(i,halo,2),&
                                                   tmpfval)
        end do
      end do  
!-----------------------------------------------------------------------------------!          
    !CASE NORTH EAST     
    case(neast)
      ! east zone
      do halo=1,2
        do i=-1,nc+1 
          ibaseref=fvm%ibase(i,halo,1)                                      
          if (ibaseref<1-nhc) then
!             write(*,*) "ibaseref is out of range 8b"
             ibaseref=1-nhc
          end if
          fnewval(i,halo,1)=cubic_equispace_interp(fvm%dbeta, fvm%interp(i,halo,1),&
                                                   fcube(nc+halo,ibaseref:ibaseref+3))
        end do
      ! north zone
        do i=-1,nc+1
          ibaseref=fvm%ibase(i,halo,2)                                      
          if (ibaseref<1-nhc) then
!             write(*,*) "ibaseref is out of range 9"
             ibaseref=1-nhc
          end if
          fnewval(i,halo,2)=cubic_equispace_interp(fvm%dalpha, fvm%interp(i,halo,2),&
                                                   fcube(ibaseref:ibaseref+3,nc+halo))
        end do
      end do
      !easthalo: west zone & north zone
      !(because of the reconstruction in the halo, conservation!) 
      do halo=1,2
        do i=-1,nc+1
          ibaseref=fvm%ibasehalo(i,halo,1)                                      
          if (ibaseref<1-nhc) then
!             write(*,*) "ibaseref is out of range 9b"
             ibaseref=1-nhc
          end if

          fhalo(i,halo,1)=cubic_equispace_interp(fvm%dbeta, fvm%interphalo(i,halo,1),&
                                                   fcube(nc+1-halo,ibaseref:ibaseref+3))
        end do
        do i=1,nhc
          tmpfval(i)=fcube(nc+1-halo,nc+i)
        enddo
        do i=0,nhe+1
          ibaseref=fvm%ibasehaloex(i,halo,1)                                      
          fhaloex(i,halo,1)=cubic_equispace_interp(fvm%dalpha, fvm%interphaloex(i,halo,1),&
                                                   tmpfval)
        end do
      end do
      !northhalo: south zone & east zone
      !(because of the reconstruction in the halo, conservation!)         
      do halo=1,2
        do i=-1,nc+1
          ibaseref=fvm%ibasehalo(i,halo,2)          
          if (ibaseref<1-nhc) then
!             write(*,*) "ibaseref is out of range 10"
             ibaseref=1-nhc
          end if
          fhalo(i,halo,2)=cubic_equispace_interp(fvm%dalpha, fvm%interphalo(i,halo,2),&
                                                   fcube(ibaseref:ibaseref+3,nc+1-halo))
        end do
        do i=1,nhc
          tmpfval(i)=fcube(nc+i,nc+1-halo)
        enddo
        do i=0,nhe+1
          ibaseref=fvm%ibasehaloex(i,halo,2)                                      
          fhaloex(i,halo,2)=cubic_equispace_interp(fvm%dbeta, fvm%interphaloex(i,halo,2),&
                                                   tmpfval)
        end do
      end do  
!-----------------------------------------------------------------------------------!                
    !CASE NORTH WEST   
    case(nwest)
      ! west zone
      do halo=1,2
        do i=-1,nc+1
          ibaseref=fvm%ibase(i,halo,1)                                      
          if (ibaseref<1-nhc) then
!             write(*,*) "ibaseref is out of range 10b"
             ibaseref=1-nhc
          end if

          fnewval(i,halo,1)=cubic_equispace_interp(fvm%dbeta, fvm%interp(i,halo,1),&
                                                   fcube(1-halo,ibaseref:ibaseref+3))
        end do
      ! north zone
        do i=0,nc+2
          ibaseref=fvm%ibase(i,halo,2)              
          if (ibaseref<1-nhc) then
!             write(*,*) "ibaseref is out of range 11"
             ibaseref=1-nhc
          end if                        
          fnewval(i,halo,2)=cubic_equispace_interp(fvm%dalpha, fvm%interp(i,halo,2),&
                                                   fcube(ibaseref:ibaseref+3,nc+halo))
        end do
      end do
      !westhalo: east zone & north zone
      !(because of the reconstruction in the halo, conservation!) 
      do halo=1,2
        do i=-1,nc+1
          ibaseref=fvm%ibasehalo(i,halo,1)                                      
          if (ibaseref<1-nhc) then
!             write(*,*) "ibaseref is out of range 11b"
             ibaseref=1-nhc
          end if
          fhalo(i,halo,1)=cubic_equispace_interp(fvm%dbeta, fvm%interphalo(i,halo,1),&
                                                   fcube(halo,ibaseref:ibaseref+3))
        end do
        do i=1,nhc
          tmpfval(i)=fcube(halo,nc+nhc-i+1)
        enddo
        do i=0,nhe+1
          ibaseref=fvm%ibasehaloex(i,halo,1)                                      
          fhaloex(nhe+1-i,halo,1)=cubic_equispace_interp(fvm%dalpha, fvm%interphaloex(i,halo,1),&
                                                   tmpfval)
        end do
      end do        
      !northhalo: south zone & west zone
      !(because of the reconstruction in the halo, conservation!) 
      do halo=1,2
        do i=0,nc+2
          ibaseref=fvm%ibasehalo(i,halo,2)                         
!          if (ibaseref<1-nhc) write(*,*) "ibaseref is out of range 12"             
          fhalo(i,halo,2)=cubic_equispace_interp(fvm%dalpha, fvm%interphalo(i,halo,2),&
                                                   fcube(ibaseref:ibaseref+3,nc+1-halo))
        end do
        do i=1,nhc
          tmpfval(i)=fcube(-i+1,nc+1-halo)
        enddo
        do i=0,nhe+1
          ibaseref=fvm%ibasehaloex(i,halo,2)                                      
          fhaloex(i,halo,2)=cubic_equispace_interp(fvm%dbeta, fvm%interphaloex(i,halo,2),&
                                                   tmpfval)
        end do
      end do
!-----------------------------------------------------------------------------------!      
      !THIS CASE SHOULD NOT HAPPEN!     
      case default
        print *, 'Fatal Error in second select statement:'
        print *, 'fvm_reconstruction_mod.F90 subroutine fillhalo_cubic!' 
        stop
    end select
end subroutine fillhalo_cubic
!END SUBROUTINE FILLHALO_CUBIC--------------------------------------------CE-for FVM!

! ----------------------------------------------------------------------------------!
!FUNCTION CUBIC_EQUISPACE_INTERP------------------------------------------CE-for FVM!
! AUTHOR: CHRISTOPH ERATH, 17.October 2011                                          !
! DESCRIPTION: apply cubic interpolation of Lagrange type on the specified array of !
!        values, where all points are equally spaced.                               !
!        alpha/beta coordinates                                                     !
!                                                                                   !
! CALLS: none                                                                       !
! INPUT: dx ... spacing of points, alpha/beta                                       !
!        x  ... X coordinate where interpolation is to be applied, alpha/beta       !
!        y  ... array of 4 values = f(x + k * dx) where k = 0,1,2,3                 !
! OUTPUT/RETURN:                                                                    !
!        cubic_equispace_interp ... interpolation value                             !                 
!-----------------------------------------------------------------------------------!
function cubic_equispace_interp(dx, x, y)

  implicit none
  real (kind=real_kind)                             :: cubic_equispace_interp
  real (kind=real_kind), intent(in)                 :: dx, x
  real (kind=real_kind), intent(in)                 :: y(:)

  cubic_equispace_interp =                                                   &
    (-y(1) / (6 * dx**3)) * (x - dx) * (x - 2 * dx) * (x - 3 * dx) +         &
    ( y(2) / (2 * dx**3)) * (x) * (x - 2 * dx) * (x - 3 * dx) +              &
    (-y(3) / (2 * dx**3)) * (x) * (x - dx) * (x - 3 * dx) +                  &
    ( y(4) / (6 * dx**3)) * (x) * (x - dx) * (x - 2 * dx)
end function cubic_equispace_interp
!END FUNCTION CUBIC_EQUISPACE_INTERP--------------------------------------CE-for FVM!

! ----------------------------------------------------------------------------------!
!SUBROUTINE RECONSTRUCTION_CUBIC------------------------------------------CE-for FVM!
! AUTHOR: CHRISTOPH ERATH, 17.October 2011                                          !
! DESCRIPTION: construct a piecewise cubic reconstruction in the alpha/beta         !
!              coordinate system                                                    !
!                                                                                   !
! CALLS: reconstruct_cubic_onface, reconstruct_cubic_x, reconstruct_cubic_y         !
! INPUT: fcube    ... tracer values incl. halo region                               !
!        fvm    ... fvm structure                                               !
!        fnewval  ... interpolated values for halo regions on cube edges            !
!        fhalo    ... interpolated values for reconstruction in the halo region     !
!        fhaloex  ... interpolated values for reconstruction in the halo region     !
!                     if it is a corner element                                     !
! OUTPUT:recons ... reconstruction coefficients df/dx, df/dy, d^2f/dx^2, d^2/dy^2   !                 
!-----------------------------------------------------------------------------------!
subroutine reconstruction_cubic(fcube,fvm,fnewval,fhalo,fhaloex,recons)
  implicit none
  real (kind=real_kind),   &
       dimension(1-nhc:nc+nhc, 1-nhc:nc+nhc), intent(in)          :: fcube
  type (fvm_struct), intent(in)                                 :: fvm
  real (kind=real_kind), dimension(-1:nc+2,2,2), intent(in)       :: fnewval
  real (kind=real_kind), dimension(-1:nc+2,2,2), intent(in)       :: fhalo
  real (kind=real_kind), dimension(0:nhe+1,2,2), intent(in)       :: fhaloex   
  
  real (kind=real_kind), &
       dimension(1:5,1-nhe:nc+nhe,1-nhe:nc+nhe), intent(out)      :: recons
        
  real (kind=real_kind),   &
        dimension(1-nhc:nc+nhc, 1-nhc:nc+nhc)                     :: fcubenew
  real (kind=real_kind), dimension(-2:nc+3,nhe+4)                 :: fotherface
  integer (kind=int_kind)                                         :: i, j

  fcubenew=fcube
  ! for an interior element (cubeboundary==0) the halo is correct
  !special cases, calculate the reconstruction coefficients in the halo, where the
  !halo region is on a different face as the element (elements, which have a edge
  !or a corner on the cube)  
  !the construction for the corner elements is tricky
  ! we construct the reconstruction coefficients df/dx, df/dy, d^2f/dx^2, d^2/dy^2
  ! and d^2/dxdy for the interior cells and for cells in the halo, which are on the
  ! same face as the original cell, note that below we calculate also the coefficients
  ! for cells in the halo on different faces, which we will overwrite below
  select case (fvm%cubeboundary) 
!-----------------------------------------------------------------------------------!    
    case (0)
      call reconstruct_cubic_onface(fvm, fcubenew, recons)
!-----------------------------------------------------------------------------------!
    case(west) 
      do i=-1,nc+2    
        fcubenew(0,i)=fnewval(i,1,1)
        fcubenew(-1,i)=fnewval(i,2,1)       
      end do
      call reconstruct_cubic_onface(fvm, fcubenew, recons)
      ! calculation for the halo zone
      do j=-2,nc+3
        do i=1,nhe+2
          fotherface(j,i)=fcube(-nhe-2+i,j)                  
        end do
      end do
      do j=-1,nc+2
        fotherface(j,nhe+3)=fhalo(j,1,1)
        fotherface(j,nhe+4)=fhalo(j,2,1)
      end do
      
      call reconstruct_cubic_haloy(fotherface, recons,fvm%dalpha, fvm%dbeta,&
            fvm%spherecentroid, fvm%invx1,fvm%invy1,fvm%jx_min1,fvm%jx_max1-1, &
             fvm%jy_min1,fvm%jy_max1-1,fvm%swap1)    
!-----------------------------------------------------------------------------------!        
    case(east)
      do i=-1,nc+2
        fcubenew(nc+1,i)=fnewval(i,1,1)
        fcubenew(nc+2,i)=fnewval(i,2,1)
      end do  
      call reconstruct_cubic_onface(fvm, fcubenew, recons)  
      !calculation for the halo zone
      do j=-2,nc+3
        do i=1,nhe+2
          fotherface(j,i+2)=fcube(nc+i,j)                  
        end do
      end do
      do j=-1,nc+2
        fotherface(j,1)=fhalo(j,2,1)
        fotherface(j,2)=fhalo(j,1,1)
      end do
      call reconstruct_cubic_haloy(fotherface, recons,fvm%dalpha, fvm%dbeta,&
            fvm%spherecentroid, fvm%invx1,fvm%invy1,fvm%jx_min1,fvm%jx_max1-1, &
             fvm%jy_min1,fvm%jy_max1-1,fvm%swap1)
!-----------------------------------------------------------------------------------!        
    case(north)
      do i=-1,nc+2
        fcubenew(i,nc+1)=fnewval(i,1,2)
        fcubenew(i,nc+2)=fnewval(i,2,2) 
      end do
      call reconstruct_cubic_onface(fvm, fcubenew, recons)
!       !calculation for the halo zone
      do i=-2,nc+3
        do j=1,nhe+2
          fotherface(i,j+2)=fcube(i,nc+j)                  
        end do
      end do
      do i=-1,nc+2
        fotherface(i,2)=fhalo(i,1,2)
        fotherface(i,1)=fhalo(i,2,2)
      end do
      call reconstruct_cubic_halox(fotherface, recons,fvm%dalpha, fvm%dbeta,&
            fvm%spherecentroid, fvm%invx1,fvm%invy1,fvm%jx_min1,fvm%jx_max1-1, &
             fvm%jy_min1,fvm%jy_max1-1,fvm%swap1)
!-----------------------------------------------------------------------------------! 
    case(south)
      do i=-1,nc+2    
        fcubenew(i,0)=fnewval(i,1,2)
        fcubenew(i,-1)=fnewval(i,2,2)        
      end do 
      call reconstruct_cubic_onface(fvm, fcubenew, recons)
      !calculation
      do i=-2,nc+3
        do j=1,nhe+2
          fotherface(i,j)=fcube(i,-nhe-2+j)                  
        end do
      end do
      do i=-1,nc+2
        fotherface(i,nhe+3)=fhalo(i,1,2)
        fotherface(i,nhe+4)=fhalo(i,2,2)
      end do
      call reconstruct_cubic_halox(fotherface, recons,fvm%dalpha, fvm%dbeta,&
            fvm%spherecentroid, fvm%invx1,fvm%invy1,fvm%jx_min1,fvm%jx_max1-1, &
             fvm%jy_min1,fvm%jy_max1-1,fvm%swap1)
!-----------------------------------------------------------------------------------!       
    case(swest) 
      do i=0,nc+1    !note: fcube(0,0) is not needed, we will overwrite it
        !west
        fcubenew(0,i)=fnewval(i,1,1)
        fcubenew(-1,i)=fnewval(i,2,1)
        !south
        fcubenew(i,0)=fnewval(i,1,2)
        fcubenew(i,-1)=fnewval(i,2,2)        
      end do
      fcubenew(0,nc+2)=fnewval(nc+2,1,1)
      fcubenew(nc+2,0)=fnewval(nc+2,1,2)
      ! the soutwest corner is an extrapolation, we take an average of the neighbours
      fcubenew(0,0)=(fcubenew(0,1)+fcubenew(1,0)+fcubenew(-1,0)+fcubenew(0,-1))/4.0D0 
      call reconstruct_cubic_onface(fvm, fcubenew, recons)
      !calculation for the halo zone
       do i=1,nc+3
         do j=1,nhe+2
           fotherface(i,j)=fcube(i,-nhe-2+j)                  
         end do
       end do
       do i=0,nc+2
         fotherface(i,nhe+3)=fhalo(i,1,2)
         fotherface(i,nhe+4)=fhalo(i,2,2)
       end do
       do j=0,nhe+1
         fotherface(0,nhe+3-j)=fhaloex(j,1,2)
         fotherface(-1,nhe+3-j)=fhaloex(j,2,2)
       end do
       fotherface(0,nhe+3)=(fotherface(-1,nhe+3)+fotherface(0,nhe+2)+ &
                               fotherface(1,nhe+3)+fotherface(0,nhe+4))/4.0D0
       call reconstruct_cubic_halox(fotherface, recons,fvm%dalpha, fvm%dbeta,&
             fvm%spherecentroid, fvm%invx1,fvm%invy1,fvm%jx_min1,fvm%jx_max1-1, &
              fvm%jy_min1,fvm%jy_max1-1,fvm%swap1)
      !west
      do j=1,nc+3
        do i=1,nhe+2
          fotherface(j,i)=fcube(-nhe-2+i,j)                  
        end do
      end do
      do j=0,nc+2
        fotherface(j,nhe+3)=fhalo(j,1,1)
        fotherface(j,nhe+4)=fhalo(j,2,1)
      end do
      do j=0,nhe+1
        fotherface(0,nhe+3-j)=fhaloex(j,1,1)
        fotherface(-1,nhe+3-j)=fhaloex(j,2,1)
      end do
      fotherface(0,nhe+3)=(fotherface(0,nhe+4)+fotherface(1,nhe+3)+ &
                              fotherface(-1,nhe+3)+fotherface(0,nhe+2))/4.0D0                  
      call reconstruct_cubic_haloy(fotherface, recons,fvm%dalpha, fvm%dbeta,&
            fvm%spherecentroid, fvm%invx2,fvm%invy2,fvm%jx_min2,&
             fvm%jx_max2-1, fvm%jy_min2,fvm%jy_max2-1,fvm%swap2)                               
!-----------------------------------------------------------------------------------!        
    case(seast)
      do i=0,nc+1
        !east
        fcubenew(nc+1,i)=fnewval(i,1,1)
        fcubenew(nc+2,i)=fnewval(i,2,1)
        !south
        fcubenew(i,0)=fnewval(i,1,2)
        fcubenew(i,-1)=fnewval(i,2,2)       
      end do
      fcubenew(nc+1,nc+2)=fnewval(nc+2,1,1)
      fcubenew(-1,0)=fnewval(-1,1,2)
      ! the souteast corner is an extrapolation, we take an average of the neighbours
      fcubenew(nc+1,0)=(fcubenew(nc+1,1)+fcubenew(nc,0)+&
                        fcubenew(nc+2,0)+fcubenew(nc+1,-1))/4.0D0
      call reconstruct_cubic_onface(fvm, fcubenew, recons)
      !calculation for the halo zone
      do i=-2,nc
        do j=1,nhe+2
          fotherface(i,j)=fcube(i,-nhe-2+j)                  
        end do
      end do
      do i=-1,nc+1
        fotherface(i,nhe+3)=fhalo(i,1,2)
        fotherface(i,nhe+4)=fhalo(i,2,2)
      end do
      do j=0,nhe+1
        fotherface(nc+1,nhe+3-j)=fhaloex(j,1,2)
        fotherface(nc+2,nhe+3-j)=fhaloex(j,2,2)
      end do
      fotherface(nc+1,nhe+3)=(fotherface(nc,nhe+3)+fotherface(nc+1,nhe+2)+ &
                             fotherface(nc+2,nhe+3)+fotherface(nc+1,nhe+4))/4.0D0

      call reconstruct_cubic_halox(fotherface, recons,fvm%dalpha, fvm%dbeta,&
             fvm%spherecentroid, fvm%invx1,fvm%invy1,fvm%jx_min1,fvm%jx_max1-1, &
             fvm%jy_min1,fvm%jy_max1-1,fvm%swap1)
      !east
      do j=1,nc+3
        do i=1,nhe+2
          fotherface(j,i+2)=fcube(nc+i,j)                  
        end do
      end do
      do j=0,nc+2
        fotherface(j,1)=fhalo(j,2,1)
        fotherface(j,2)=fhalo(j,1,1)
      end do
      do j=0,nhe+1
        fotherface(0,2+j)=fhaloex(j,1,1)
        fotherface(-1,2+j)=fhaloex(j,2,1)
      end do
      fotherface(0,2)=(fotherface(0,1)+fotherface(1,2)+ &
                              fotherface(0,3)+fotherface(-1,2))/4.0D0
      call reconstruct_cubic_haloy(fotherface, recons,fvm%dalpha, fvm%dbeta,&
           fvm%spherecentroid, fvm%invx2,fvm%invy2,fvm%jx_min2,&
           fvm%jx_max2-1, fvm%jy_min2,fvm%jy_max2-1,fvm%swap2)
!-----------------------------------------------------------------------------------!                          
    case(neast)
      do i=0,nc+1
        !east
        fcubenew(nc+1,i)=fnewval(i,1,1)
        fcubenew(nc+2,i)=fnewval(i,2,1)
        !north
        fcubenew(i,nc+1)=fnewval(i,1,2)
        fcubenew(i,nc+2)=fnewval(i,2,2)       
      end do
      fcubenew(nc+1,-1)=fnewval(-1,1,1)
      fcubenew(-1,nc+1)=fnewval(-1,1,2)
      ! the northeast corner is an extrapolation, we take an average of the neighbours
      fcubenew(nc+1,nc+1)=(fcubenew(nc,nc+1)+fcubenew(nc+1,nc)+&
                           fcubenew(nc+1,nc+2)+fcubenew(nc+2,nc+1))/4.0D0
      call reconstruct_cubic_onface(fvm, fcubenew, recons)
      !calculation for the halo zone
      do i=-2,nc
        do j=1,nhe+2
          fotherface(i,j+2)=fcube(i,nc+j)                  
        end do
      end do
      do i=-1,nc+1
        fotherface(i,2)=fhalo(i,1,2)
        fotherface(i,1)=fhalo(i,2,2)
      end do
      do j=0,nhe+1
        fotherface(nc+1,2+j)=fhaloex(j,1,2)
        fotherface(nc+2,2+j)=fhaloex(j,2,2)
      end do
      fotherface(nc+1,2)=(fotherface(nc,2)+fotherface(nc+1,1)+ &
                             fotherface(nc+2,2)+fotherface(nc+1,3))/4.0D0
      
      call reconstruct_cubic_halox(fotherface, recons,fvm%dalpha, fvm%dbeta,&
            fvm%spherecentroid, fvm%invx1,fvm%invy1,fvm%jx_min1,fvm%jx_max1-1, &
             fvm%jy_min1,fvm%jy_max1-1,fvm%swap1)
      
      !east
      do j=-2,nc
        do i=1,nhe+2
          fotherface(j,i+2)=fcube(nc+i,j)                  
        end do
      end do
      do j=-1,nc+1
        fotherface(j,1)=fhalo(j,2,1)
        fotherface(j,2)=fhalo(j,1,1)
      end do
      do j=0,nhe+1
        fotherface(nc+1,2+j)=fhaloex(j,1,1)
        fotherface(nc+2,2+j)=fhaloex(j,2,1)
      end do
      fotherface(nc+1,2)=(fotherface(nc+1,1)+fotherface(nc+2,2)+ &
                             fotherface(nc+1,3)+fotherface(nc,2))/4.0D0
      call reconstruct_cubic_haloy(fotherface, recons,fvm%dalpha, fvm%dbeta,&
            fvm%spherecentroid, fvm%invx2,fvm%invy2,fvm%jx_min2,&
             fvm%jx_max2-1, fvm%jy_min2,fvm%jy_max2-1,fvm%swap2)
!-----------------------------------------------------------------------------------!       
    case(nwest)    
      do i=0,nc+1
        !west
        fcubenew(0,i)=fnewval(i,1,1)
        fcubenew(-1,i)=fnewval(i,2,1)
        !north
        fcubenew(i,nc+1)=fnewval(i,1,2)
        fcubenew(i,nc+2)=fnewval(i,2,2) 
      end do     
      fcubenew(0,-1)=fnewval(-1,1,1)
      fcubenew(nc+2,nc+1)=fnewval(nc+2,1,2)  
      ! the northwest corner is an extrapolation, we take an average of the neighbours
      fcubenew(0,nc+1)=(fcubenew(0,nc)+fcubenew(1,nc+1)+&
                        fcubenew(-1,nc+1)+fcubenew(0,nc+2))/4.0D0  
      call reconstruct_cubic_onface(fvm, fcubenew, recons)
      !calculation for the halo zone
      do i=1,nc+3
        do j=1,nhe+2
          fotherface(i,j+2)=fcube(i,nc+j)                  
        end do
      end do
      do i=0,nc+2
        fotherface(i,2)=fhalo(i,1,2)
        fotherface(i,1)=fhalo(i,2,2)
      end do
      do j=0,nhe+1
        fotherface(0,2+j)=fhaloex(j,1,2)
        fotherface(-1,2+j)=fhaloex(j,2,2)
      end do
      fotherface(0,2)=(fotherface(-1,2)+fotherface(0,1)+ &
                              fotherface(1,2)+fotherface(0,3))/4.0D0
      
      call reconstruct_cubic_halox(fotherface, recons,fvm%dalpha, fvm%dbeta,&
            fvm%spherecentroid, fvm%invx1,fvm%invy1,fvm%jx_min1,fvm%jx_max1-1, &
             fvm%jy_min1,fvm%jy_max1-1,fvm%swap1)
      !WEST
      do j=-2,nc
        do i=1,nhe+2
          fotherface(j,i)=fcube(-nhe-2+i,j)                  
        end do
      end do
      do j=-1,nc+1
        fotherface(j,nhe+3)=fhalo(j,1,1)
        fotherface(j,nhe+4)=fhalo(j,2,1)
      end do
      do j=0,nhe+1
        fotherface(nc+1,nhe+3-j)=fhaloex(j,1,1)
        fotherface(nc+2,nhe+3-j)=fhaloex(j,2,1)
      end do
      fotherface(nc+1,nhe+3)=(fotherface(nc+1,nhe+4)+fotherface(nc+2,nhe+3)+ &
                              fotherface(nc,nhe+3)+fotherface(nc+1,nhe+2))/4.0D0 
      call reconstruct_cubic_haloy(fotherface, recons,fvm%dalpha, fvm%dbeta,&
            fvm%spherecentroid, fvm%invx2,fvm%invy2,fvm%jx_min2,&
             fvm%jx_max2-1, fvm%jy_min2,fvm%jy_max2-1,fvm%swap2)                          
!-----------------------------------------------------------------------------------!            
    !THIS CASE SHOULD NOT HAPPEN!     
    case default
      print *,'Fatal Error in first select statement:'
      call abortmp('fvm_reconstruction_mod.F90 subroutine reconstruction_cubic!' )
      
  end select
end subroutine reconstruction_cubic
!END SUBROUTINE RECONSTRUCTION_CUBIC--------------------------------------CE-for FVM!

! ----------------------------------------------------------------------------------!
!SUBROUTINE RECONSTRUCT_CUBIC_ONFACE--------------------------------------CE-for FVM!
! AUTHOR: CHRISTOPH ERATH, 14.November 2011                                         !
! DESCRIPTION: construct a piecewise cubic reconstruction in the alpha/beta         !
!              coordinate system for the interior element and halo zone on the same !
!              face as the element                                                  !
!                                                                                   !
! INPUT: fvm    ... fvm structure                                               !
!        fcube    ... tracer values in the patch                                    !                                    
! INPUT/OUTPUT:                                                                     !
!          recons ... reconstruction coefficients df/dx, df/dy, d^2f/dx^2, d^2/dy^2 !                 
!-----------------------------------------------------------------------------------!
subroutine reconstruct_cubic_onface(fvm, fcube, recons)
  implicit none

  type (fvm_struct), intent(in)                               :: fvm
  real (kind=real_kind),   &
      dimension(1-nhc:nc+nhc, 1-nhc:nc+nhc)                     :: fcube
  real (kind=real_kind), &
      dimension(1:5,1-nhe:nc+nhe,1-nhe:nc+nhe), intent(inout)   :: recons
  
  integer  :: i, j
        
  do j = fvm%jy_min, fvm%jy_max-1
    do i = fvm%jx_min, fvm%jx_max-1
      ! df/dx
      recons(1,i,j) = (- fcube(i+2,j) + 8 * fcube(i+1,j)      &
                       - 8 * fcube(i-1,j) + fcube(i-2,j))/(12 * fvm%dalpha)     
      ! df/dy
      recons(2,i,j) = (- fcube(i,j+2) + 8 * fcube(i,j+1)      &
                       - 8 * fcube(i,j-1) + fcube(i,j-2)) /( 12 * fvm%dbeta)
      ! Stretching
      recons(1,i,j) = recons(1,i,j) / (1 + fvm%spherecentroid(1,i,j)**2)
      
      recons(2,i,j) = recons(2,i,j) / (1 + fvm%spherecentroid(2,i,j)**2)

      ! d^2f/dx^2
      recons(3,i,j) = (- fcube(i+2,j) + 16 * fcube(i+1,j)     &
                       - 30 * fcube(i,j) + 16 * fcube(i-1,j)  &
                       - fcube(i-2,j)) / (12 * fvm%dalpha**2)

      ! d^2f/dy^2
      recons(4,i,j) = (- fcube(i,j+2) + 16 * fcube(i,j+1)     &
                       - 30 * fcube(i,j) + 16 * fcube(i,j-1)  &
                       - fcube(i,j-2)) / (12 * fvm%dbeta**2)
      ! d^2f/dxdy
      recons(5,i,j) = (+ fcube(i+1,j+1) - fcube(i-1,j+1)  &
                       - fcube(i+1,j-1) + fcube(i-1,j-1)) / (4 * fvm%dalpha * fvm%dbeta)
      ! stretching
      recons(3,i,j) = (- 2 * fvm%spherecentroid(1,i,j) * &
                      (1 + fvm%spherecentroid(1,i,j)**2) * recons(1,i,j) &
                       + recons(3,i,j)) / (1 + fvm%spherecentroid(1,i,j)**2)**2
      recons(4,i,j) = (- 2 * fvm%spherecentroid(2,i,j) * &
                      (1 + fvm%spherecentroid(2,i,j)**2) * recons(2,i,j) &
                       + recons(4,i,j)) / (1 + fvm%spherecentroid(2,i,j)**2)**2
      recons(5,i,j) = recons(5,i,j) / ((1 + fvm%spherecentroid(1,i,j)**2) * (1 + fvm%spherecentroid(2,i,j)**2))
      ! scaling
      recons(3,i,j) = 0.5 * recons(3,i,j)
      recons(4,i,j) = 0.5 * recons(4,i,j)
    enddo
  enddo
end subroutine reconstruct_cubic_onface
!END SUBROUTINE RECONSTRUCTION_CUBIC_ONFACE-------------------------------CE-for FVM!

! ----------------------------------------------------------------------------------!
!SUBROUTINE RECONSTRUCTION_CUBIC_HALOX------------------------------------CE-for FVM!
! AUTHOR: CHRISTOPH ERATH, 17.October 2011                                          !
! DESCRIPTION: construct a piecewise cubic reconstruction in the alpha/beta         !
!              coordinate system, the loop goes over cells in beta (y) direction    !
!                                                                                   !
! CALLS: reconstruct_cubic_x, reconstruct_cubic_y                                   !
! INPUT: f        ... tracer values in the patch                                    !                                    
!        delta1   ... cell step in alpha (coordinate) direction                     !
!        delta2   ... cell step in beta (coordinate) direction                      !
!        centroid ... area averaged moments for each cell                           !
!        cubeboundary..value 1-4 if an edge or 5-8 if a corner of the element       !
!                      is on a cube edge/corner, note this subroutine will never    !
!                      called with 0 for this value                                 !
!        invx ... if the values in f for the alpha coordinate are in the opposite   !
!                 direction, take -1, otherwise 1                                   !
!        invy ... if the values in f for the beta coordinate are in the opposite    !
!                 direction, take -1, otherwise 1                                   !
!        idxa ... start index of cell to reconstruct                                !
!        idxb ... end index of cell to reconstruct                                  !
!        swap ... if the entries as to be swaped in recons                          !
! INPUT/OUTPUT:                                                                     !
!        recons ... reconstruction coefficients df/dx, df/dy, d^2f/dx^2, d^2/dy^2   !
!-----------------------------------------------------------------------------------!
subroutine reconstruct_cubic_halox(f, recons, delta1, delta2,centroid, invx,invy,&
                                   idxa,idxb,idya,idyb,swap)
  implicit none
  real (kind=real_kind),dimension(-2:nc+3,1:nhe+4), intent(in)      :: f
  real (kind=real_kind), intent(in)                                 :: delta1, delta2
  real (kind=real_kind), dimension(5,1-nhe:nc+nhe,1-nhe:nc+nhe), intent(in) :: centroid
  real (kind=real_kind), dimension(1:5,1-nhe:nc+nhe,1-nhe:nc+nhe), intent(inout) :: recons
  integer, intent(in)                     :: invx, invy,idxa,idxb,idya,idyb
  logical, intent(in)                     :: swap
    
  integer (kind=int_kind)                 :: i, j, tmpj,sw
  
  if (swap) then
    sw=1
  else
    sw=0
  endif
  
  do j = idya, idyb    
    tmpj=3+j-idya
    do i = idxa, idxb
      ! df/dx
      recons(1+sw,i,j) = invx * (- f(i+2,tmpj) + 8 * f(i+1,tmpj) - 8 * f(i-1,tmpj) + &
                                   f(i-2,tmpj)) /  ( 12 * delta1)
      ! df/dy
      recons(2-sw,i,j) = invy * (- f(i,tmpj+2) + 8 * f(i,tmpj+1) - 8 * f(i,tmpj-1) + & 
                                   f(i,tmpj-2)) / ( 12 * delta2)
      ! Stretching
      recons(1,i,j) = recons(1,i,j) / (1 + centroid(1,i,j)**2)
      recons(2,i,j) = recons(2,i,j) / (1 + centroid(2,i,j)**2)
      ! d^2f/dx^2
      recons(3+sw,i,j) = (- f(i+2,tmpj) + 16 * f(i+1,tmpj) - 30 * f(i,tmpj) + &
                       16 * f(i-1,tmpj) - f(i-2,tmpj)) / (12 * delta1**2)
      ! d^2f/dy^2
      recons(4-sw,i,j) = (- f(i,tmpj+2) + 16 * f(i,tmpj+1) - 30 * f(i,tmpj) + &
                       16 * f(i,tmpj-1) - f(i,tmpj-2)) / (12 * delta2**2)

      ! d^2f/dxdy
      recons(5,i,j) = invx * invy * (+ f(i+1,tmpj+1) - f(i-1,tmpj+1) - &
                       f(i+1,tmpj-1) + f(i-1,tmpj-1))  / (4 * delta1 * delta2)

      ! Stretching
      recons(3,i,j) = (- 2 * centroid(1,i,j) * (1 + centroid(1,i,j)**2) * recons(1,i,j)    &
                       + recons(3,i,j)) / (1 + centroid(1,i,j)**2)**2

      recons(4,i,j) = (- 2 * centroid(2,i,j) * (1 + centroid(2,i,j)**2) * recons(2,i,j)   &
                       + recons(4,i,j)) / (1 + centroid(2,i,j)**2)**2

      recons(5,i,j) = recons(5,i,j) / ((1 + centroid(1,i,j)**2) * (1 + centroid(2,i,j)**2))
      ! Scaling
      recons(3,i,j) = 0.5 * recons(3,i,j)
      recons(4,i,j) = 0.5 * recons(4,i,j)
    enddo
  enddo 
end subroutine reconstruct_cubic_halox
!END SUBROUTINE RECONSTRUCTION_CUBIC_HALOX--------------------------------CE-for FVM!

! ----------------------------------------------------------------------------------!
!SUBROUTINE RECONSTRUCTION_CUBIC_HALOY------------------------------------CE-for FVM!
! AUTHOR: CHRISTOPH ERATH, 17.October 2011                                          !
! DESCRIPTION: construct a piecewise cubic reconstruction in the alpha/beta         !
!              coordinate system, the loop goes over cells in beta (y) direction    !
!                                                                                   !
! CALLS: reconstruct_cubic_x, reconstruct_cubic_y                                   !
! INPUT: f        ... tracer values in the patch                                    !                                    
!        delta1   ... cell step in alpha (coordinate) direction                     !
!        delta2   ... cell step in beta (coordinate) direction                      !
!        centroid ... area averaged moments for each cell                           !
!        cubeboundary..value 1-4 if an edge or 5-8 if a corner of the element       !
!                      is on a cube edge/corner, note this subroutine will never    !
!                      called with 0 for this value                                 !
!        invx ... if the values in f for the alpha coordinate are in the opposite   !
!                 direction, take -1, otherwise 1                                   !
!        invy ... if the values in f for the beta coordinate are in the opposite    !
!                 direction, take -1, otherwise 1                                   !
!        idxa ... start index of cell to reconstruct                                !
!        idxb ... end index of cell to reconstruct                                  !
!        swap ... if the entries as to be swaped in recons                          !
! INPUT/OUTPUT:                                                                     !
!        recons ... reconstruction coefficients df/dx, df/dy, d^2f/dx^2, d^2/dy^2   !                 
!-----------------------------------------------------------------------------------!
subroutine reconstruct_cubic_haloy(f, recons, delta1, delta2,centroid, invx,invy,&
                                  idxa,idxb,idya,idyb,swap)
  implicit none
  real (kind=real_kind),dimension(-2:nc+3,1:nhe+4), intent(in)      :: f
  real (kind=real_kind), intent(in)                                 :: delta1, delta2
  real (kind=real_kind), dimension(5,1-nhe:nc+nhe,1-nhe:nc+nhe), intent(in) :: centroid
  real (kind=real_kind), dimension(1:5,1-nhe:nc+nhe,1-nhe:nc+nhe), intent(inout) :: recons
  integer, intent(in)                     :: invx, invy,idxa,idxb,idya,idyb
  logical, intent(in)                     :: swap
     
  integer (kind=int_kind)                 :: i, j, tmpi,sw
  
  if (swap) then
    sw=1
  else
    sw=0
  endif
  do i = idxa, idxb
    tmpi=3+i-idxa
    do j = idya, idyb
      ! df/dx
      recons(1+sw,i,j) = invx * (- f(j,tmpi+2) + 8 * f(j,tmpi+1) - 8 * f(j,tmpi-1) + &
                                   f(j,tmpi-2)) /  ( 12 * delta1)
      ! df/dy
      recons(2-sw,i,j) = invy * (- f(j+2,tmpi) + 8 * f(j+1,tmpi) - 8 * f(j-1,tmpi) + & 
                                   f(j-2,tmpi)) / ( 12 * delta2)
      ! Stretching
      recons(1,i,j) = recons(1,i,j) / (1 + centroid(1,i,j)**2)
      recons(2,i,j) = recons(2,i,j) / (1 + centroid(2,i,j)**2)
      ! d^2f/dx^2
      recons(3+sw,i,j) = (- f(j,tmpi+2) + 16 * f(j,tmpi+1) - 30 * f(j,tmpi) + &
                       16 * f(j,tmpi-1) - f(j,tmpi-2)) / (12 * delta1**2)
      ! d^2f/dy^2
      recons(4-sw,i,j) = (- f(j+2,tmpi) + 16 * f(j+1,tmpi) - 30 * f(j,tmpi) + &
                       16 * f(j-1,tmpi) - f(j-2,tmpi)) / (12 * delta2**2)

      ! d^2f/dxdy
      recons(5,i,j) = invx * invy * (+ f(j+1,tmpi+1) - f(j+1,tmpi-1) - &
                       f(j-1,tmpi+1) + f(j-1,tmpi-1))  / (4 * delta1 * delta2)

      ! Stretching
      recons(3,i,j) = (- 2 * centroid(1,i,j) * (1 + centroid(1,i,j)**2) * recons(1,i,j)    &
                       + recons(3,i,j)) / (1 + centroid(1,i,j)**2)**2

      recons(4,i,j) = (- 2 * centroid(2,i,j) * (1 + centroid(2,i,j)**2) * recons(2,i,j)   &
                       + recons(4,i,j)) / (1 + centroid(2,i,j)**2)**2

      recons(5,i,j) = recons(5,i,j) / ((1 + centroid(1,i,j)**2) * (1 + centroid(2,i,j)**2))
      ! Scaling
      recons(3,i,j) = 0.5 * recons(3,i,j)
      recons(4,i,j) = 0.5 * recons(4,i,j)
    enddo
  enddo 
end subroutine reconstruct_cubic_haloy
!END SUBROUTINE RECONSTRUCTION_CUBIC_HALOY--------------------------------CE-for FVM!

end module fvm_reconstruction_mod
