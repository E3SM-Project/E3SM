!***********************************************************************************!
! Author: Christoph Erath                                                           !
! Date: May 2011                                                                    !
! Module for CSLAM exchange buffer aso                                              !
!***********************************************************************************!
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module cslam_ghostcell_mod
  use kinds, only               : real_kind
  use edge_mod, only            : EdgeBuffer_t, EdgeDescriptor_t
  use dimensions_mod, only      : nc, nhc
  use control_mod, only         : north, south, east, west, neast, nwest, seast, swest
  
implicit none
private

public :: ghostcellpack
public :: ghostcellunpack
public :: ghostpointspack
public :: ghostpointsunpack

logical, private :: threadsafe=.true.

contains
!-----------------------------------------------------------------------------------!
! Christoph Erath, May 2011                                                         !
! Pack cells on the boundary of c into buf, use structure from shallow water
! This should be replaced by CSLAM buffer system
!-----------------------------------------------------------------------------------!
subroutine ghostcellpack(edge,c,vlyr,kptr,desc)

    type (EdgeBuffer_t)                      :: edge
    integer,              intent(in)         :: vlyr
    real (kind=real_kind),intent(in)         :: c(1-nhc:nc+nhc,1-nhc:nc+nhc,vlyr)
    integer,              intent(in)         :: kptr
    type (EdgeDescriptor_t),intent(in)       :: desc

    ! Local variables
    integer :: i,k,ir

    integer :: is,ie,in,iw
    
    real (kind=real_kind)                       :: NaN=-1.0
    NaN=sqrt(NaN)

    if(.not. threadsafe) then
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
       threadsafe=.true.
    end if
    
    is = desc%putmapP(south)
    ie = desc%putmapP(east)
    in = desc%putmapP(north)
    iw = desc%putmapP(west)

   do k=1,vlyr
      ! dirty hack, cslam uses one less, provide NaN there
      edge%buf(kptr+k,is+1)   = c(1  ,1 ,k)
      edge%buf(kptr+k,ie+1)   = c(nc ,1 ,k)
      edge%buf(kptr+k,in+1)   = c(1  ,nc,k)
      edge%buf(kptr+k,iw+1)   = c(1  ,1 ,k)
      ! set always the second entry to NaN
      edge%buf(kptr+k,is+2)   = NaN
      edge%buf(kptr+k,ie+2)   = NaN
      edge%buf(kptr+k,in+2)   = NaN
      edge%buf(kptr+k,iw+2)   = NaN
      do i=2,nc
         edge%buf(kptr+k,is+i+1)   = c(i  ,1 ,k)
         edge%buf(kptr+k,ie+i+1)   = c(nc ,i ,k)
         edge%buf(kptr+k,in+i+1)   = c(i  ,nc,k)
         edge%buf(kptr+k,iw+i+1)   = c(1  ,i ,k)
      enddo
   end do

    !  This is really kludgy way to setup the index reversals
    !  But since it is so a rare event not real need to spend time optimizing
    if(desc%reverse(south)) then
       !is = desc%putmapP(south)  is is already set
       do k=1,vlyr
         ! dirty hack, cslam uses one less, provide NaN there
         edge%buf(kptr+k,is+1)   = c(nc  ,1 ,k)  
         edge%buf(kptr+k,is+2)   = NaN
          do i=1,nc-1
             ir = nc+1-i+1
             edge%buf(kptr+k,is+ir)=c(i,1,k)
          enddo
       enddo
    endif

    if(desc%reverse(east)) then
       !ie = desc%putmapP(east)  ie is already set
       do k=1,vlyr
         edge%buf(kptr+k,ie+1)   = c(nc  ,nc ,k)  
         edge%buf(kptr+k,ie+2)   = NaN
          do i=1,nc-1
             ir = nc+1-i+1
             edge%buf(kptr+k,ie+ir)=c(nc,i,k)
          enddo
       enddo
    endif

    if(desc%reverse(north)) then
       !in = desc%putmapP(north) in is already set
       do k=1,vlyr
         edge%buf(kptr+k,in+1)   = c(nc  ,nc ,k)  
         edge%buf(kptr+k,in+2)   = NaN
          do i=1,nc-1
             ir = nc+1-i+1
             edge%buf(kptr+k,in+ir)=c(i,nc,k)
          enddo
       enddo
    endif

    if(desc%reverse(west)) then
       !iw = desc%putmapP(west) iw is already set
       do k=1,vlyr
         edge%buf(kptr+k,iw+1)   = c(1  ,nc ,k)  
         edge%buf(kptr+k,iw+2)   = NaN
          do i=1,nc-1
             ir = nc+1-i+1
             edge%buf(kptr+k,iw+ir)=c(1,i,k)
          enddo
       enddo
    endif
    
!  This is really kludgy way to setup the index reversals
!  But since it is so a rare event not real need to spend time optimizing
    if (desc%putmapP(swest) /= -1) then
         do k=1,vlyr
            edge%buf(kptr+k,desc%putmapP(swest)+1)=c(1  ,1 ,k)
         end do
    end if

    if (desc%putmapP(seast) /= -1) then
       do k=1,vlyr
         edge%buf(kptr+k,desc%putmapP(seast)+1)=c(nc ,1 ,k)
       end do
    end if

    if (desc%putmapP(neast) /= -1) then
       do k=1,vlyr
         edge%buf(kptr+k,desc%putmapP(neast)+1)=c(nc ,nc,k)
       end do
    end if

    if (desc%putmapP(nwest) /= -1) then
       do k=1,vlyr
         edge%buf(kptr+k,desc%putmapP(nwest)+1)=c(1  ,nc,k)
       end do
    end if
end subroutine ghostcellpack

!-----------------------------------------------------------------------------------!
! Christoph Erath, May 2011                                                         !
! UnPack cells on the boundary of v into buf, use structure from shallow water
! This should be replaced by CSLAM buffer system
!-----------------------------------------------------------------------------------!
subroutine ghostcellunpack(edge,c,vlyr,kptr,desc)
  type (EdgeBuffer_t),         intent(in)   :: edge

  integer,               intent(in)         :: vlyr
  real (kind=real_kind), intent(inout)      :: c(1-nhc:nc+nhc,1-nhc:nc+nhc,vlyr)
  integer,               intent(in)         :: kptr
  type (EdgeDescriptor_t)                   :: desc

  ! Local
  integer :: i,k
  integer :: is,ie,in,iw

  threadsafe=.false.

  is=desc%getmapP(south)
  ie=desc%getmapP(east)
  in=desc%getmapP(north)
  iw=desc%getmapP(west)
  
  do k=1,vlyr 
    c(1  ,1-nhc  ,k) = edge%buf(kptr+k,is+1  )
    c(nc+nhc ,1  ,k) = edge%buf(kptr+k,ie+1  )
    c(1  ,nc+nhc ,k) = edge%buf(kptr+k,in+1  )
    c(1-nhc  ,1  ,k) = edge%buf(kptr+k,iw+1  )    
! second entry in buf is NaN (cslam is not node based)
    do i=2,nc 
        c(i  ,1-nhc  ,k) = edge%buf(kptr+k,is+i+1  )
        c(nc+nhc ,i  ,k) = edge%buf(kptr+k,ie+i+1  )
        c(i  ,nc+nhc ,k) = edge%buf(kptr+k,in+i+1  )
        c(1-nhc  ,i  ,k) = edge%buf(kptr+k,iw+i+1  )
    end do
 end do

   
  if(desc%getmapP(swest) /= -1) then 
    do k=1,vlyr
      c(1-nhc  ,1-nhc ,k)=edge%buf(kptr+k,desc%getmapP(swest)+1)
    enddo
  endif
  
  if(desc%getmapP(seast) /= -1) then 
    do k=1,vlyr
      c(nc+nhc ,1-nhc ,k)=edge%buf(kptr+k,desc%getmapP(seast)+1)
    enddo
  endif
  
  if(desc%getmapP(neast) /= -1) then 
    do k=1,vlyr
      c(nc+nhc ,nc+nhc,k)=edge%buf(kptr+k,desc%getmapP(neast)+1)
    enddo
  endif
  
  if(desc%getmapP(nwest) /= -1) then 
    do k=1,vlyr
      c(1-nhc  ,nc+nhc,k)=edge%buf(kptr+k,desc%getmapP(nwest)+1)
    enddo
  endif

end subroutine ghostcellunpack


!-----------------------------------------------------------------------------------!
! Christoph Erath, May 2011                                                         !
! Pack points on the boundary of c into buf, use structure from shallow water
! This should be replaced by CSLAM buffer system
!-----------------------------------------------------------------------------------!
subroutine ghostpointspack(edge,points,kptr,desc)

    type (EdgeBuffer_t)                      :: edge
    real (kind=real_kind),intent(in)         :: points(1-nhc:nc+1+nhc,1-nhc:nc+1+nhc)
    integer,              intent(in)         :: kptr
    type (EdgeDescriptor_t),intent(in)       :: desc

    ! Local variables
    integer :: i,k,ir

    integer :: is,ie,in,iw

    if(.not. threadsafe) then
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
       threadsafe=.true.
    end if
    
    is = desc%putmapP(south)
    ie = desc%putmapP(east)
    in = desc%putmapP(north)
    iw = desc%putmapP(west)
    
    do i=1,nc+1
       edge%buf(kptr,is+i)   = points(i  ,2)
       edge%buf(kptr,ie+i)   = points(nc ,i)
       edge%buf(kptr,in+i)   = points(i  ,nc)
       edge%buf(kptr,iw+i)   = points(2  ,i)
    enddo

    !  This is really kludgy way to setup the index reversals
    !  But since it is so a rare event not real need to spend time optimizing
    if(desc%reverse(south)) then
       !is = desc%putmapP(south)  is is already set
      do i=1,nc+1
         ir = nc+2-i
         edge%buf(kptr,is+ir)=points(i,2)
      enddo
    endif

    if(desc%reverse(east)) then
       !ie = desc%putmapP(east)  ie is already set
      do i=1,nc+1
         ir = nc+2-i
         edge%buf(kptr,ie+ir)=points(nc,i)
      enddo
    endif

    if(desc%reverse(north)) then
       !in = desc%putmapP(north) in is already set
      do i=1,nc+1
         ir = nc+2-i
         edge%buf(kptr,in+ir)=points(i,nc)
      enddo
    endif

    if(desc%reverse(west)) then
       !iw = desc%putmapP(west) iw is already set
      do i=1,nc+1
         ir = nc+2-i
         edge%buf(kptr,iw+ir)=points(2,i)
      enddo
    endif
    
!  This is really kludgy way to setup the index reversals
!  But since it is so a rare event not real need to spend time optimizing
    if (desc%putmapP(swest) /= -1) then
        edge%buf(kptr,desc%putmapP(swest)+1)=points(2  ,2)
    end if

    if (desc%putmapP(seast) /= -1) then
        edge%buf(kptr,desc%putmapP(seast)+1)=points(nc ,2)
    end if

    if (desc%putmapP(neast) /= -1) then
        edge%buf(kptr,desc%putmapP(neast)+1)=points(nc ,nc)
    end if

    if (desc%putmapP(nwest) /= -1) then
        edge%buf(kptr,desc%putmapP(nwest)+1)=points(2  ,nc)
    end if
end subroutine ghostpointspack

!-----------------------------------------------------------------------------------!
! Christoph Erath, May 2011                                                         !
! UnPack points on the boundary of v into buf, use structure from shallow water
! This should be replaced by CSLAM buffer system
!-----------------------------------------------------------------------------------!
subroutine ghostpointsunpack(edge,points,kptr,desc)
  type (EdgeBuffer_t),         intent(in)   :: edge

  real (kind=real_kind), intent(inout)      :: points(1-nhc:nc+1+nhc,1-nhc:nc+1+nhc)
  integer,               intent(in)         :: kptr
  type (EdgeDescriptor_t)                   :: desc

  ! Local
  integer :: i,k
  integer :: is,ie,in,iw
  
  real (kind=real_kind)                       :: NaN=-1.0
  NaN=sqrt(NaN)

  threadsafe=.false.

  is=desc%getmapP(south)
  ie=desc%getmapP(east)
  in=desc%getmapP(north)
  iw=desc%getmapP(west)
      
  do i=1,nc+1 
    points(i  ,1-nhc) = edge%buf(kptr,is+i  )
    points(nc+1+nhc ,i) = edge%buf(kptr,ie+i  )
    points(i  ,nc+nhc+1) = edge%buf(kptr,in+i  )
    points(1-nhc  ,i) = edge%buf(kptr,iw+i  )
  end do

   
  if(desc%getmapP(swest) /= -1) then 
    points(1-nhc  ,1-nhc)=edge%buf(kptr,desc%getmapP(swest)+1)
  else
      points(1-nhc  ,1-nhc)=NaN    
  endif
  
  if(desc%getmapP(seast) /= -1) then 
      points(nc+1+nhc ,1-nhc)=edge%buf(kptr,desc%getmapP(seast)+1)
  else
      points(nc+1+nhc ,1-nhc)=NaN
  endif
  
  if(desc%getmapP(neast) /= -1) then 
    points(nc+1+nhc ,nc+1+nhc)=edge%buf(kptr,desc%getmapP(neast)+1)
  else  
    points(nc+1+nhc ,nc+1+nhc)=NaN
  endif
  
  if(desc%getmapP(nwest) /= -1) then 
    points(1-nhc  ,nc+1+nhc)=edge%buf(kptr,desc%getmapP(nwest)+1)
  else
    points(1-nhc  ,nc+1+nhc)=NaN
  endif

end subroutine ghostpointsunpack

end module cslam_ghostcell_mod
