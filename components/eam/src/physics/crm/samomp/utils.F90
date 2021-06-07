module utils
  use task_util_mod
  implicit none

contains

  integer function lenstr (string)

    ! returns string's length ignoring the rightmost blank and null characters

    implicit none
    character *(*) string
    integer k
    lenstr = 0
    do k = 1,len(string)
      if (string(k:k).ne.' '.and.string(k:k).ne.char(0)) then
        lenstr = lenstr+1
      end if
    end do
    111  return
  end



  subroutine averageXY(f,dimx1,dimx2,dimy1,dimy2,dimz,fm)

    use grid
    use params, only: crm_rknd
    implicit none
    integer dimx1, dimx2, dimy1, dimy2, dimz
    real(crm_rknd) f(dimx1:dimx2, dimy1:dimy2, dimz),fm(nzm)
    real(8) ff,factor
    integer i,j,k
    factor = 1./dble(nx*ny)
    do k =1,nzm
      ff = 0.
      do j =1,ny
        do i =1,nx
          ff = ff + f(i,j,k)
        end do
      end do
      ff = ff*factor
      fm(k) = real(ff,crm_rknd)
    end do
  end


  subroutine averageXY_MPI(f,dimx1,dimx2,dimy1,dimy2,dimz,fm)

    use grid
    use params, only: crm_rknd
    implicit none
    integer dimx1, dimx2, dimy1, dimy2, dimz
    real(crm_rknd) f(dimx1:dimx2, dimy1:dimy2, dimz),fm(nzm)
    real(8) fm1(nzm),fm2(nzm),factor
    integer i,j,k
    factor = 1./dble(nx*ny)
    do k =1,nzm
      fm1(k) = 0.
      do j =1,ny
        do i =1,nx
          fm1(k) = fm1(k) + f(i,j,k)
        end do
      end do
      fm1(k) = fm1(k) * factor
    end do
    if(dompi) then
      do k =1,nzm
        fm2(k) = fm1(k)
      end do
      call task_sum_real8(fm2,fm1,nzm)
      do k=1,nzm
        fm(k)=real( fm1(k)/dble(nsubdomains) ,crm_rknd)
      end do
    else
      do k=1,nzm
        fm(k)=real( fm1(k) ,crm_rknd)
      end do
    endif
  end




  subroutine fminmax_print(name,f,dimx1,dimx2,dimy1,dimy2,dimz)
    use grid
    use params, only: crm_rknd
    implicit none
    integer dimx1, dimx2, dimy1, dimy2, dimz
    real(crm_rknd) f(dimx1:dimx2, dimy1:dimy2, dimz),fmn(nz),fmx(nz)
    character *(*) name
    real(crm_rknd) fmin(1),fmax(1),fff(1)
    integer i,j,k

    do k=1,dimz
      if(dimx2.eq.1.and.dimy2.eq.1) then
        fmn(k) = f(1,1,k)
        fmx(k) = f(1,1,k)
      else
        fmn(k) = 1.e30
        fmx(k) =-1.e30
        do j=1,ny
          do i=1,nx
            fmn(k) = min(fmn(k),f(i,j,k))
            fmx(k) = max(fmx(k),f(i,j,k))
          end do
        enddo
      end if
    enddo
    fmin(1) = 1.e30
    fmax(1) =-1.e30
    do k=1,dimz
      fmin(1) = min(fmin(1),fmn(k))
      fmax(1) = max(fmax(1),fmx(k))
    end do

    if(dompi) then
      fff(1)=fmax(1)
      call task_max_real(fff(1),fmax(1),1)
      fff(1)=fmin(1)
      call task_min_real(fff(1),fmin(1),1)
    end if
    if(masterproc) print *,name,fmin,fmax
  end




  subroutine setvalue(f,n,f0)
    use params, only: crm_rknd
    implicit none
    integer n
    real(crm_rknd) f(n), f0
    integer k
    do k=1,n
      f(k)=f0
    end do
  end

  ! determine number of byte in a record in direct access files (can be anything, from 1 to 8):
  ! can't assume 1 as it is compiler and computer dependent
  integer function bytes_in_rec()
    implicit none
    character*8 str
    integer n, err
    open(1,status ='scratch',access ='direct',recl=1)
    do n = 1,8
      write(1,rec=1,iostat=err) str(1:n)
      if (err.ne.0) exit
      bytes_in_rec = n
    enddo
    close(1,status='delete')
  end

end module utils
