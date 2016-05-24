#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

! =========================================
! rotation_init:
!
! Initialize cube rotation terms resulting
! from changing cube face coordinate systems
!
! =========================================

  subroutine rotation_init(elem, rot_type)
  use kinds, only : real_kind
  use control_mod, only : north, south, east, west, neast, nwest, seast, swest
  use dimensions_mod, only : nv
  use element_mod, only : element_t
  use cube_mod, only : contravariant_rot,covariant_rot, vmap, convert_gbl_index 
  implicit none
  type(element_t), intent(inout) :: elem(:)
    character(len=*) rot_type

    ! =======================================
    ! Local variables
    ! =======================================

    integer :: ie,je
    integer :: myface_no        ! current element face number
    integer :: nbrface_no       ! neighbor element face number
    integer :: inbr
    integer :: nrot,irot
    integer :: ii,i,j
    integer :: ir,jr
    integer :: nbrnum

    real (kind=real_kind) :: Dloc(2,2,nv)
    real (kind=real_kind) :: Drem(2,2,nv)
    real (kind=real_kind) :: x1,x2

    do ii=1,size(elem)

       call convert_gbl_index(elem(ii)%vertex%number,ie,je,myface_no)

       nrot   = 0
       do inbr=1,8
          if (elem(ii)%vertex%nbrs(inbr) /= -1 ) then
              call convert_gbl_index(elem(ii)%vertex%nbrs(inbr),ie,je,nbrface_no)
              if (myface_no /= nbrface_no) nrot=nrot+1
          end if
       end do
       if(associated(elem(ii)%desc%rot)) then
         deallocate(elem(ii)%desc%rot)
         NULLIFY(elem(ii)%desc%rot)
       end if

       ! =====================================================
       ! If there are neighbors on other cube faces, allocate 
       ! an array of rotation matrix structs.
       ! =====================================================

       if (nrot > 0) then
          allocate(elem(ii)%desc%rot(nrot))
          elem(ii)%desc%use_rotation=1
          irot=0          
          do inbr=1,8

             call convert_gbl_index(elem(ii)%vertex%nbrs(inbr),ie,je,nbrface_no)
             
             ! The cube edge (myface_no,nbrface_no) and inbr defines 
             ! a unique rotation given by (D^-1) on myface_no x (D) on nbrface_no

             if (myface_no /= nbrface_no .and. elem(ii)%vertex%nbrs(inbr) /= -1 ) then           

                irot=irot+1

                if (inbr <= 4) then      
                   allocate(elem(ii)%desc%rot(irot)%R(2,2,nv))  ! edge
                else                     
                   allocate(elem(ii)%desc%rot(irot)%R(2,2,1 ))   ! corner
                end if

                ! must compute Dloc on my face, Drem on neighbor face, 
                ! for each point on edge or corner.

                ! ==================================== 
                ! Equatorial belt east/west neighbors
                ! ==================================== 

                if (nbrface_no <= 4 .and. myface_no <= 4) then

                   if (inbr == west) then
                      do j=1,nv
                         x1 = elem(ii)%cartv(1,j)%x
                         x2 = elem(ii)%cartv(1,j)%y
                         call Vmap(Dloc(1,1,j), x1,x2,myface_no)
                         call Vmap(Drem(1,1,j),-x1,x2,nbrface_no)
                      end do
                   else if (inbr == east) then
                      do j=1,nv
                         x1 = elem(ii)%cartv(nv,j)%x
                         x2 = elem(ii)%cartv(nv,j)%y
                         call Vmap(Dloc(1,1,j), x1,x2,myface_no)
                         call Vmap(Drem(1,1,j),-x1,x2,nbrface_no)
                      end do
                   else if (inbr == swest ) then
                      x1 = elem(ii)%cartv(1,1)%x
                      x2 = elem(ii)%cartv(1,1)%y
                      call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                      call Vmap(Drem(1,1,1),-x1,x2,nbrface_no)
                   else if (inbr == nwest ) then
                      x1 = elem(ii)%cartv(1,nv)%x
                      x2 = elem(ii)%cartv(1,nv)%y
                      call Vmap(Dloc(1,1,1), x1,x2,myface_no)
                      call Vmap(Drem(1,1,1),-x1,x2,nbrface_no)
                   else if (inbr == seast ) then
                      x1 = elem(ii)%cartv(nv,1)%x
                      x2 = elem(ii)%cartv(nv,1)%y
                      call Vmap(Dloc(1,1,1), x1,x2,myface_no)
                      call Vmap(Drem(1,1,1),-x1,x2,nbrface_no)
                   else if (inbr == neast ) then
                      x1 = elem(ii)%cartv(nv,nv)%x
                      x2 = elem(ii)%cartv(nv,nv)%y
                      call Vmap(Dloc(1,1,1), x1,x2,myface_no)
                      call Vmap(Drem(1,1,1),-x1,x2,nbrface_no)
                   end if
                   
                end if

                ! Northern Neighbors of Equatorial Belt

                if ( myface_no <= 4 .and. nbrface_no == 6 ) then
                   if (inbr == north) then
                      do i=1,nv
                         ir=nv+1-i
                         x1 = elem(ii)%cartv(i,nv)%x
                         x2 = elem(ii)%cartv(i,nv)%y
                         if ( myface_no == 1) then
                            call Vmap(Dloc(1,1,i), x1,x2,myface_no)
                            call Vmap(Drem(1,1,i),x1,-x2,nbrface_no)
                         end if
                         if ( myface_no == 2) then
                            call Vmap(Dloc(1,1,i),x1,x2,myface_no)
                            call Vmap(Drem(1,1,i),x2,x1,nbrface_no)
                            nbrnum = elem(ii)%vertex%nbrs(inbr)
                            
                         end if
                         if ( myface_no == 3) then
                            call Vmap(Dloc(1,1,ir), x1,x2,myface_no)
                            call Vmap(Drem(1,1,ir),-x1,x2,nbrface_no)
                         end if
                         if ( myface_no == 4) then
                            call Vmap(Dloc(1,1,ir), x1,x2,myface_no)
                            call Vmap(Drem(1,1,ir),-x2,-x1,nbrface_no)
                         end if
                      end do
                   else if (inbr == nwest) then
                      x1 = elem(ii)%cartv(1,nv)%x
                      x2 = elem(ii)%cartv(1,nv)%y
                      call Vmap(Dloc(1,1,1), x1,x2,myface_no)
                      if ( myface_no == 1) call Vmap(Drem(1,1,1),x1,-x2,nbrface_no)
                      if ( myface_no == 2) call Vmap(Drem(1,1,1),x2, x1,nbrface_no)
                      if ( myface_no == 3) call Vmap(Drem(1,1,1),-x1,x2,nbrface_no)
                      if ( myface_no == 4) call Vmap(Drem(1,1,1),-x2,-x1,nbrface_no)
                   else if (inbr == neast) then
                      x1 = elem(ii)%cartv(nv,nv)%x
                      x2 = elem(ii)%cartv(nv,nv)%y
                      call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                      if ( myface_no == 1) call Vmap(Drem(1,1,1),x1,-x2,nbrface_no)
                      if ( myface_no == 2) call Vmap(Drem(1,1,1),x2, x1,nbrface_no)
                      if ( myface_no == 3) call Vmap(Drem(1,1,1),-x1,x2,nbrface_no)
                      if ( myface_no == 4) call Vmap(Drem(1,1,1),-x2,-x1,nbrface_no)
                   end if

                end if

                ! Southern Neighbors of Equatorial Belt

                if ( myface_no <= 4 .and. nbrface_no == 5 ) then
                   if (inbr == south) then
                      do i=1,nv
                         ir=nv+1-i
                         x1 = elem(ii)%cartv(i,1)%x
                         x2 = elem(ii)%cartv(i,1)%y
                         if ( myface_no == 1) then
                            call Vmap(Dloc(1,1,i), x1, x2,myface_no)
                            call Vmap(Drem(1,1,i), x1,-x2,nbrface_no)
                         end if
                         if ( myface_no == 2) then
                            call Vmap(Dloc(1,1,ir),x1,x2,myface_no)
                            call Vmap(Drem(1,1,ir),-x2,-x1,nbrface_no)
                         end if
                         if ( myface_no == 3) then
                            call Vmap(Dloc(1,1,ir), x1,x2,myface_no)
                            call Vmap(Drem(1,1,ir),-x1,x2,nbrface_no)
                         end if
                         if ( myface_no == 4) then
                            call Vmap(Dloc(1,1,i), x1,x2,myface_no)
                            call Vmap(Drem(1,1,i), x2,x1,nbrface_no)
                         end if
                      end do
                   else if (inbr == swest) then
                      x1 = elem(ii)%cartv(1,1)%x
                      x2 = elem(ii)%cartv(1,1)%y
                      call Vmap(Dloc(1,1,1),x1,x2,myface_no)


                      if ( myface_no == 1) call Vmap(Drem(1,1,1),x1,-x2,nbrface_no)
                      if ( myface_no == 2) call Vmap(Drem(1,1,1),-x2,-x1,nbrface_no)
                      if ( myface_no == 3) call Vmap(Drem(1,1,1),-x1,x2,nbrface_no)
                      if ( myface_no == 4) call Vmap(Drem(1,1,1),x2,x1,nbrface_no)

                   else if (inbr == seast) then
                      x1 = elem(ii)%cartv(nv,1)%x
                      x2 = elem(ii)%cartv(nv,1)%y
                      call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                      if ( myface_no == 1) call Vmap(Drem(1,1,1),x1,-x2,nbrface_no)
                      if ( myface_no == 2) call Vmap(Drem(1,1,1),-x2,-x1,nbrface_no)
                      if ( myface_no == 3) call Vmap(Drem(1,1,1),-x1,x2,nbrface_no)
                      if ( myface_no == 4) call Vmap(Drem(1,1,1),x2,x1,nbrface_no)
                   end if

                end if

                ! Neighbors of Northern Capping Face Number 6
                 
                if ( myface_no == 6 ) then
                   if (nbrface_no == 1) then
                      if (inbr == south) then
                         do i=1,nv
                            x1 = elem(ii)%cartv(i,1)%x
                            x2 = elem(ii)%cartv(i,1)%y
                            call Vmap(Dloc(1,1,i),x1,x2,myface_no)
                            call Vmap(Drem(1,1,i),x1,-x2,nbrface_no)
                         end do
                      else if (inbr == swest) then
                         x1 = elem(ii)%cartv(1,1)%x
                         x2 = elem(ii)%cartv(1,1)%y
                         call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                         call Vmap(Drem(1,1,1),x1,-x2,nbrface_no)
                      else if (inbr == seast) then
                         x1 = elem(ii)%cartv(nv,1)%x
                         x2 = elem(ii)%cartv(nv,1)%y
                         call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                         call Vmap(Drem(1,1,1),x1,-x2,nbrface_no)
                      end if
                   else if (nbrface_no == 2) then
                      if (inbr == east) then
                         do j=1,nv
                            x1 = elem(ii)%cartv(nv,j)%x
                            x2 = elem(ii)%cartv(nv,j)%y
                            call Vmap(Dloc(1,1,j),x1,x2,myface_no)
                            call Vmap(Drem(1,1,j),x2,x1,nbrface_no)
                         end do
                      else if (inbr == seast) then
                         x1 = elem(ii)%cartv(nv,1)%x
                         x2 = elem(ii)%cartv(nv,1)%y
                         call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                         call Vmap(Drem(1,1,1),x2,x1,nbrface_no)
                      else if (inbr == neast) then
                         x1 = elem(ii)%cartv(nv,nv)%x
                         x2 = elem(ii)%cartv(nv,nv)%y
                         call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                         call Vmap(Drem(1,1,1),x2,x1,nbrface_no)
                      end if
                   else if (nbrface_no == 3) then
                      if (inbr == north) then
                         do i=1,nv
                            ir =nv+1-i
                            x1 = elem(ii)%cartv(i,nv)%x
                            x2 = elem(ii)%cartv(i,nv)%y
                            call Vmap(Dloc(1,1,ir),x1,x2,myface_no)
                            call Vmap(Drem(1,1,ir),-x1,x2,nbrface_no)
                         end do
                      else if (inbr == nwest) then
                         x1 = elem(ii)%cartv(1,nv)%x
                         x2 = elem(ii)%cartv(1,nv)%y
                         call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                         call Vmap(Drem(1,1,1),-x1,x2,nbrface_no)
                      else if (inbr == neast) then
                         x1 = elem(ii)%cartv(nv,nv)%x
                         x2 = elem(ii)%cartv(nv,nv)%y
                         call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                         call Vmap(Drem(1,1,1),-x1,x2,nbrface_no)
                      end if
                   else if (nbrface_no == 4) then
                      if (inbr == west) then
                         do j=1,nv
                            jr=nv+1-j
                            x1 = elem(ii)%cartv(1,j)%x
                            x2 = elem(ii)%cartv(1,j)%y
                            call Vmap(Dloc(1,1,jr), x1, x2,myface_no )
                            call Vmap(Drem(1,1,jr),-x2,-x1,nbrface_no)
                         end do
                      else if (inbr == swest) then
                         x1 = elem(ii)%cartv(1,1)%x
                         x2 = elem(ii)%cartv(1,1)%y
                         call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                         call Vmap(Drem(1,1,1),-x2,-x1,nbrface_no)
                      else if (inbr == nwest) then
                         x1 = elem(ii)%cartv(1,nv)%x
                         x2 = elem(ii)%cartv(1,nv)%y
                         call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                         call Vmap(Drem(1,1,1),-x2,-x1,nbrface_no)
                      end if
                   end if
                end if

                ! Neighbors of South Capping Face Number 5

                if ( myface_no == 5 ) then
                   if (nbrface_no == 1) then
                      if (inbr == north) then
                         do i=1,nv
                            x1 = elem(ii)%cartv(i,nv)%x
                            x2 = elem(ii)%cartv(i,nv)%y
                            call Vmap(Dloc(1,1,i),x1,x2,myface_no)
                            call Vmap(Drem(1,1,i),x1,-x2,nbrface_no)
                         end do
                      else if (inbr == nwest) then
                         x1 = elem(ii)%cartv(1,nv)%x
                         x2 = elem(ii)%cartv(1,nv)%y
                         call Vmap(Dloc(:,:,1),x1,x2,myface_no)
                         call Vmap(Drem(:,:,1),x1,-x2,nbrface_no)
                      else if (inbr == neast) then
                         x1 = elem(ii)%cartv(nv,nv)%x
                         x2 = elem(ii)%cartv(nv,nv)%y
                         call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                         call Vmap(Drem(1,1,1),x1,-x2,nbrface_no)
                      end if
                   else if (nbrface_no == 2) then
                      if (inbr == east) then
                         do j=1,nv
                            jr=nv+1-j
                            x1 = elem(ii)%cartv(nv,j)%x
                            x2 = elem(ii)%cartv(nv,j)%y
                            call Vmap(Dloc(1,1,jr),x1,  x2,myface_no)
                            call Vmap(Drem(1,1,jr),-x2,-x1,nbrface_no)
                         end do
                      else if (inbr == seast) then
                         x1 = elem(ii)%cartv(nv,1)%x
                         x2 = elem(ii)%cartv(nv,1)%y
                         call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                         call Vmap(Drem(1,1,1),-x2,-x1,nbrface_no)
                      else if (inbr == neast) then
                         x1 = elem(ii)%cartv(nv,nv)%x
                         x2 = elem(ii)%cartv(nv,nv)%y
                         call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                         call Vmap(Drem(1,1,1),-x2,-x1,nbrface_no)
                      end if
                   else if (nbrface_no == 3) then
                      if (inbr == south) then
                         do i=1,nv
                            ir=nv+1-i
                            x1 = elem(ii)%cartv(i,1)%x
                            x2 = elem(ii)%cartv(i,1)%y
                            call Vmap(Dloc(1,1,ir),x1,x2,myface_no)
                            call Vmap(Drem(1,1,ir),-x1,x2,nbrface_no)
                         end do
                      else if (inbr == swest) then
                         x1 = elem(ii)%cartv(1,1)%x
                         x2 = elem(ii)%cartv(1,1)%y
                         call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                         call Vmap(Drem(1,1,1),-x1,x2,nbrface_no)
                      else if (inbr == seast) then
                         x1 = elem(ii)%cartv(nv,1)%x
                         x2 = elem(ii)%cartv(nv,1)%y
                         call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                         call Vmap(Drem(1,1,1),-x1,x2,nbrface_no)
                      end if
                   else if (nbrface_no == 4) then
                      if (inbr == west) then
                         do j=1,nv
                            x1 = elem(ii)%cartv(1,j)%x
                            x2 = elem(ii)%cartv(1,j)%y
                            call Vmap(Dloc(1,1,j),x1,x2,myface_no)
                            call Vmap(Drem(1,1,j),x2,x1,nbrface_no)
                         end do
                      else if (inbr == swest) then
                         x1 = elem(ii)%cartv(1,1)%x
                         x2 = elem(ii)%cartv(1,1)%y
                         call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                         call Vmap(Drem(1,1,1),x2,x1,nbrface_no)
                      else if (inbr == nwest) then
                         x1 = elem(ii)%cartv(1,nv)%x
                         x2 = elem(ii)%cartv(1,nv)%y
                         call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                         call Vmap(Drem(1,1,1),x2,x1,nbrface_no)
                      end if
                   end if
                end if

                elem(ii)%desc%rot(irot)%nbr = inbr
                if (rot_type == "covariant") then
                   do i=1,SIZE(elem(ii)%desc%rot(irot)%R(:,:,:),3)
                      elem(ii)%desc%rot(irot)%R(:,:,i)=covariant_rot(Dloc(:,:,i),Drem(:,:,i))
                   end do           
                else if (rot_type == "contravariant") then
                   do i=1,SIZE(elem(ii)%desc%rot(irot)%R(:,:,:),3)
                      elem(ii)%desc%rot(irot)%R(:,:,i)=contravariant_rot(Dloc(:,:,i),Drem(:,:,i))
                   end do           
                end if                   

             endif

          end do
       end if

    end do

  end subroutine rotation_init

