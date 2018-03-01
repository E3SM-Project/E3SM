module cmparray_mod

  use shr_kind_mod, only : r8 => shr_kind_r8
  
  implicit none
  private
  save
  
  public expdaynite, cmpdaynite

  interface CmpDayNite
    module procedure CmpDayNite_1d_R
    module procedure CmpDayNite_2d_R
    module procedure CmpDayNite_3d_R
    module procedure CmpDayNite_1d_R_Copy
    module procedure CmpDayNite_2d_R_Copy
    module procedure CmpDayNite_3d_R_Copy
    module procedure CmpDayNite_1d_I
    module procedure CmpDayNite_2d_I
    module procedure CmpDayNite_3d_I
  end interface ! CmpDayNite

  interface ExpDayNite
    module procedure ExpDayNite_1d_R
    module procedure ExpDayNite_2d_R
    module procedure ExpDayNite_3d_R
    module procedure ExpDayNite_1d_I
    module procedure ExpDayNite_2d_I
    module procedure ExpDayNite_3d_I
  end interface ! ExpDayNite

  interface cmparray
    module procedure cmparray_1d_R
    module procedure cmparray_2d_R
    module procedure cmparray_3d_R
  end interface ! cmparray

  interface chksum
    module procedure chksum_1d_R
    module procedure chksum_2d_R
    module procedure chksum_3d_R
    module procedure chksum_1d_I
    module procedure chksum_2d_I
    module procedure chksum_3d_I
  end interface ! chksum

  contains

  subroutine CmpDayNite_1d_R(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1)
    integer, intent(in) :: Nday, Nnite
    integer, intent(in) :: il1, iu1
    integer, intent(in), dimension(Nday) :: IdxDay
    integer, intent(in), dimension(Nnite) :: IdxNite
    real(r8), intent(inout), dimension(il1:iu1) :: Array

    call CmpDayNite_3d_R(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1, 1, 1, 1, 1)

    return
  end subroutine CmpDayNite_1d_R

  subroutine CmpDayNite_2d_R(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1, il2, iu2)
    integer, intent(in) :: Nday, Nnite
    integer, intent(in) :: il1, iu1
    integer, intent(in) :: il2, iu2
    integer, intent(in), dimension(Nday) :: IdxDay
    integer, intent(in), dimension(Nnite) :: IdxNite
    real(r8), intent(inout), dimension(il1:iu1,il2:iu2) :: Array

    call CmpDayNite_3d_R(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1, il2, iu2, 1, 1)

    return
  end subroutine CmpDayNite_2d_R

  subroutine CmpDayNite_3d_R(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1, il2,iu2, il3, iu3)
    integer, intent(in) :: Nday, Nnite
    integer, intent(in) :: il1, iu1
    integer, intent(in) :: il2, iu2
    integer, intent(in) :: il3, iu3
    integer, intent(in), dimension(Nday) :: IdxDay
    integer, intent(in), dimension(Nnite) :: IdxNite
    real(r8), intent(inout), dimension(il1:iu1,il2:iu2,il3:iu3) :: Array

    real(r8), dimension(il1:iu1) :: tmp
    integer :: i, j, k


    do k = il3, iu3
      do j = il2, iu2

        tmp(1:Nnite) = Array(IdxNite(1:Nnite),j,k)
        Array(il1:il1+Nday-1,j,k) = Array(IdxDay(1:Nday),j,k)
        Array(il1+Nday:il1+Nday+Nnite-1,j,k) = tmp(1:Nnite)

      end do
    end do

    return
  end subroutine CmpDayNite_3d_R

  subroutine CmpDayNite_1d_R_Copy(InArray, OutArray, Nday, IdxDay, Nnite, IdxNite, il1, iu1)
    integer, intent(in) :: Nday, Nnite
    integer, intent(in) :: il1, iu1
    integer, intent(in), dimension(Nday) :: IdxDay
    integer, intent(in), dimension(Nnite) :: IdxNite
    real(r8), intent(in), dimension(il1:iu1) :: InArray
    real(r8), intent(out), dimension(il1:iu1) :: OutArray

    call CmpDayNite_3d_R_Copy(InArray, OutArray, Nday, IdxDay, Nnite, IdxNite, il1, iu1, 1, 1, 1, 1)

    return
  end subroutine CmpDayNite_1d_R_Copy

  subroutine CmpDayNite_2d_R_Copy(InArray, OutArray, Nday, IdxDay, Nnite, IdxNite, il1, iu1, il2, iu2)
    integer, intent(in) :: Nday, Nnite
    integer, intent(in) :: il1, iu1
    integer, intent(in) :: il2, iu2
    integer, intent(in), dimension(Nday) :: IdxDay
    integer, intent(in), dimension(Nnite) :: IdxNite
    real(r8), intent(in), dimension(il1:iu1,il2:iu2) :: InArray
    real(r8), intent(out), dimension(il1:iu1,il2:iu2) :: OutArray

    call CmpDayNite_3d_R_Copy(InArray, OutArray, Nday, IdxDay, Nnite, IdxNite, il1, iu1, il2, iu2, 1, 1)

    return
  end subroutine CmpDayNite_2d_R_Copy

  subroutine CmpDayNite_3d_R_Copy(InArray, OutArray, Nday, IdxDay, Nnite, IdxNite, il1, iu1, il2,iu2, il3, iu3)
    integer, intent(in) :: Nday, Nnite
    integer, intent(in) :: il1, iu1
    integer, intent(in) :: il2, iu2
    integer, intent(in) :: il3, iu3
    integer, intent(in), dimension(Nday) :: IdxDay
    integer, intent(in), dimension(Nnite) :: IdxNite
    real(r8), intent(in), dimension(il1:iu1,il2:iu2,il3:iu3) :: InArray
    real(r8), intent(out), dimension(il1:iu1,il2:iu2,il3:iu3) :: OutArray

    integer :: i, j, k


    do k = il3, iu3
      do j = il2, iu2

         do i=il1,il1+Nday-1
            OutArray(i,j,k) = InArray(IdxDay(i-il1+1),j,k)
         enddo
         do i=il1+Nday,il1+Nday+Nnite-1
            OutArray(i,j,k) = InArray(IdxNite(i-(il1+Nday)+1),j,k)
         enddo
        

      end do
    end do

    return
  end subroutine CmpDayNite_3d_R_Copy

  subroutine CmpDayNite_1d_I(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1)
    integer, intent(in) :: Nday, Nnite
    integer, intent(in) :: il1, iu1
    integer, intent(in), dimension(Nday) :: IdxDay
    integer, intent(in), dimension(Nnite) :: IdxNite
    integer, intent(inout), dimension(il1:iu1) :: Array

    call CmpDayNite_3d_I(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1, 1, 1, 1, 1)

    return
  end subroutine CmpDayNite_1d_I

  subroutine CmpDayNite_2d_I(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1, il2, iu2)
    integer, intent(in) :: Nday, Nnite
    integer, intent(in) :: il1, iu1
    integer, intent(in) :: il2, iu2
    integer, intent(in), dimension(Nday) :: IdxDay
    integer, intent(in), dimension(Nnite) :: IdxNite
    integer, intent(inout), dimension(il1:iu1,il2:iu2) :: Array

    call CmpDayNite_3d_I(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1, il2, iu2, 1, 1)

    return
  end subroutine CmpDayNite_2d_I

  subroutine CmpDayNite_3d_I(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1, il2,iu2, il3, iu3)
    integer, intent(in) :: Nday, Nnite
    integer, intent(in) :: il1, iu1
    integer, intent(in) :: il2, iu2
    integer, intent(in) :: il3, iu3
    integer, intent(in), dimension(Nday) :: IdxDay
    integer, intent(in), dimension(Nnite) :: IdxNite
    integer, intent(inout), dimension(il1:iu1,il2:iu2,il3:iu3) :: Array

    integer, dimension(il1:iu1) :: tmp
    integer :: i, j, k


    do k = il3, iu3
      do j = il2, iu2

        tmp(1:Nnite) = Array(IdxNite(1:Nnite),j,k)
        Array(il1:il1+Nday-1,j,k) = Array(IdxDay(1:Nday),j,k)
        Array(il1+Nday:il1+Nday+Nnite-1,j,k) = tmp(1:Nnite)

      end do
    end do

    return
  end subroutine CmpDayNite_3d_I

  subroutine ExpDayNite_1d_R(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1)
    integer, intent(in) :: Nday, Nnite
    integer, intent(in) :: il1, iu1
    integer, intent(in), dimension(Nday) :: IdxDay
    integer, intent(in), dimension(Nnite) :: IdxNite
    real(r8), intent(inout), dimension(il1:iu1) :: Array

    call ExpDayNite_3d_R(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1, 1, 1, 1, 1)

    return
  end subroutine ExpDayNite_1d_R

  subroutine ExpDayNite_2d_R(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1, il2, iu2)
    integer, intent(in) :: Nday, Nnite
    integer, intent(in) :: il1, iu1
    integer, intent(in) :: il2, iu2
    integer, intent(in), dimension(Nday) :: IdxDay
    integer, intent(in), dimension(Nnite) :: IdxNite
    real(r8), intent(inout), dimension(il1:iu1,il2:iu2) :: Array

    call ExpDayNite_3d_R(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1, il2, iu2, 1, 1)

    return
  end subroutine ExpDayNite_2d_R

  subroutine ExpDayNite_3d_R(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1, il2,iu2, il3, iu3)
    integer, intent(in) :: Nday, Nnite
    integer, intent(in) :: il1, iu1
    integer, intent(in) :: il2, iu2
    integer, intent(in) :: il3, iu3
    integer, intent(in), dimension(Nday) :: IdxDay
    integer, intent(in), dimension(Nnite) :: IdxNite
    real(r8), intent(inout), dimension(il1:iu1,il2:iu2,il3:iu3) :: Array

    real(r8), dimension(il1:iu1) :: tmp
    integer :: i, j, k


    do k = il3, iu3
      do j = il2, iu2

        tmp(1:Nday) = Array(1:Nday,j,k)
        Array(IdxNite(1:Nnite),j,k) = Array(il1+Nday:il1+Nday+Nnite-1,j,k)
        Array(IdxDay(1:Nday),j,k) = tmp(1:Nday)

      end do
    end do

    return
  end subroutine ExpDayNite_3d_R

  subroutine ExpDayNite_1d_I(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1)
    integer, intent(in) :: Nday, Nnite
    integer, intent(in) :: il1, iu1
    integer, intent(in), dimension(Nday) :: IdxDay
    integer, intent(in), dimension(Nnite) :: IdxNite
    integer, intent(inout), dimension(il1:iu1) :: Array

    call ExpDayNite_3d_I(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1, 1, 1, 1, 1)

    return
  end subroutine ExpDayNite_1d_I

  subroutine ExpDayNite_2d_I(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1, il2, iu2)
    integer, intent(in) :: Nday, Nnite
    integer, intent(in) :: il1, iu1
    integer, intent(in) :: il2, iu2
    integer, intent(in), dimension(Nday) :: IdxDay
    integer, intent(in), dimension(Nnite) :: IdxNite
    integer, intent(inout), dimension(il1:iu1,il2:iu2) :: Array

    call ExpDayNite_3d_I(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1, il2, iu2, 1, 1)

    return
  end subroutine ExpDayNite_2d_I

  subroutine ExpDayNite_3d_I(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1, il2,iu2, il3, iu3)
    integer, intent(in) :: Nday, Nnite
    integer, intent(in) :: il1, iu1
    integer, intent(in) :: il2, iu2
    integer, intent(in) :: il3, iu3
    integer, intent(in), dimension(Nday) :: IdxDay
    integer, intent(in), dimension(Nnite) :: IdxNite
    integer, intent(inout), dimension(il1:iu1,il2:iu2,il3:iu3) :: Array

    integer, dimension(il1:iu1) :: tmp
    integer :: i, j, k


    do k = il3, iu3
      do j = il2, iu2

        tmp(1:Nday) = Array(1:Nday,j,k)
        Array(IdxNite(1:Nnite),j,k) = Array(il1+Nday:il1+Nday+Nnite-1,j,k)
        Array(IdxDay(1:Nday),j,k) = tmp(1:Nday)

      end do
    end do

    return
  end subroutine ExpDayNite_3d_I

!******************************************************************************!
!                                                                              !
!                                 DEBUG                                        !
!                                                                              !
!******************************************************************************!

  subroutine cmparray_1d_R(name, Ref, New, id1, is1, ie1)
    character(*), intent(in) :: name
    integer,  intent(in) :: id1, is1, ie1
    real(r8), intent(in), dimension(id1) :: Ref
    real(r8), intent(in), dimension(id1) :: New

    call cmparray_3d_R(name, Ref, New, id1, is1, ie1, 1, 1, 1, 1, 1, 1)
  end subroutine cmparray_1d_R

  subroutine cmparray_2d_R(name, Ref, New, id1, is1, ie1, id2, is2, ie2)
    character(*), intent(in) :: name
    integer,  intent(in) :: id1, is1, ie1
    integer,  intent(in) :: id2, is2, ie2
    real(r8), intent(in), dimension(id1, id2) :: Ref
    real(r8), intent(in), dimension(id1, id2) :: New

    call cmparray_3d_R(name, Ref, New, id1, is1, ie1, id2, is2, ie2, 1, 1, 1)
  end subroutine cmparray_2d_R

  subroutine cmparray_3d_R(name, Ref, New, id1, is1, ie1, id2, is2, ie2, id3, is3, ie3)
    character(*), intent(in) :: name
    integer,  intent(in) :: id1, is1, ie1
    integer,  intent(in) :: id2, is2, ie2
    integer,  intent(in) :: id3, is3, ie3
    real(r8), intent(in), dimension(id1, id2, id3) :: Ref
    real(r8), intent(in), dimension(id1, id2, id3) :: New

    integer :: i, j, k
    integer :: nerr
    logical :: found
    real(r8):: rdiff
    real(r8), parameter :: rtol = 1.0e-13_r8

    nerr = 0

    do k = is3, ie3
      do j = is2, ie2

        found = .false.
        do i = is1, ie1
          rdiff = abs(New(i,j,k)-Ref(i,j,k))
          rdiff = rdiff / merge(abs(Ref(i,j,k)), 1.0_r8, Ref(i,j,k) /= 0.0_r8)
          if ( rdiff > rtol ) then
            found = .true.
            exit
          end if
        end do

        if ( found ) then
          do i = is1, ie1
            rdiff = abs(New(i,j,k)-Ref(i,j,k))
            rdiff = rdiff / merge(abs(Ref(i,j,k)), 1.0_r8, Ref(i,j,k) /= 0.0_r8)
            if ( rdiff > rtol ) then
              print 666, name, i, j, k, Ref(i, j, k), New(i, j, k), rdiff
              nerr = nerr + 1
              if ( nerr > 10 ) stop
            end if
          end do
        end if

      end do
    end do

    return
666 format('cmp3d: ', a10, 3(1x, i4), 3(1x, e20.14))

  end subroutine cmparray_3d_R

  subroutine chksum_1d_R(name, Ref, id1, is1, ie1)
    character(*), intent(in) :: name
    integer,  intent(in) :: id1, is1, ie1
    real(r8), intent(in), dimension(id1) :: Ref

    call chksum_3d_R(name, Ref, id1, is1, ie1, 1, 1, 1, 1, 1, 1)
  end subroutine chksum_1d_R

  subroutine chksum_1d_I(name, Ref, id1, is1, ie1)
    character(*), intent(in) :: name
    integer,  intent(in) :: id1, is1, ie1
    integer,  intent(in), dimension(id1) :: Ref

    call chksum_3d_I(name, Ref, id1, is1, ie1, 1, 1, 1, 1, 1, 1)
  end subroutine chksum_1d_I

  subroutine chksum_2d_R(name, Ref, id1, is1, ie1, id2, is2, ie2)
    character(*), intent(in) :: name
    integer,  intent(in) :: id1, is1, ie1
    integer,  intent(in) :: id2, is2, ie2
    real(r8), intent(in), dimension(id1, id2) :: Ref

    call chksum_3d_R(name, Ref, id1, is1, ie1, id2, is2, ie2, 1, 1, 1)
  end subroutine chksum_2d_R

  subroutine chksum_2d_I(name, Ref, id1, is1, ie1, id2, is2, ie2)
    character(*), intent(in) :: name
    integer,  intent(in) :: id1, is1, ie1
    integer,  intent(in) :: id2, is2, ie2
    integer,  intent(in), dimension(id1, id2) :: Ref

    call chksum_3d_I(name, Ref, id1, is1, ie1, id2, is2, ie2, 1, 1, 1)
  end subroutine chksum_2d_I

  subroutine chksum_3d_R(name, Ref, id1, is1, ie1, id2, is2, ie2, id3, is3, ie3)
    character(*), intent(in) :: name
    integer,  intent(in) :: id1, is1, ie1
    integer,  intent(in) :: id2, is2, ie2
    integer,  intent(in) :: id3, is3, ie3
!orig    real(r8), intent(in), dimension(id1, id2, id3) :: Ref
    real(r8), intent(in), dimension(is1:ie1, is2:ie2, is3:ie3) :: Ref

    real(r8) :: chksum
    real(r8) :: rmin, rmax
    integer :: i, j, k
    integer :: imin, jmin, kmin
    integer :: imax, jmax, kmax

    imin = is1 ; jmin = is2 ; kmin = is3
    imax = is1 ; jmax = is2 ; kmax = is3
    rmin = Ref(is1, is2, is3) ; rmax = rmin

    chksum = 0.0_r8

    do k = is3, ie3
      do j = is2, ie2
        do i = is1, ie1
          chksum = chksum + abs(Ref(i,j,k))
          if ( Ref(i,j,k) < rmin ) then
            rmin = Ref(i,j,k)
            imin = i ; jmin = j ; kmin = k
          end if
          if ( Ref(i,j,k) > rmax ) then
            rmax = Ref(i,j,k)
            imax = i ; jmax = j ; kmax = k
          end if
        end do
      end do
    end do

    print 666, name, chksum, imin, jmin, kmin, imax, jmax, kmax
666 format('chksum: ', a8, 1x, e20.14, 6(1x, i4))

  end subroutine chksum_3d_R

  subroutine chksum_3d_I(name, Ref, id1, is1, ie1, id2, is2, ie2, id3, is3, ie3)
    character(*), intent(in) :: name
    integer,  intent(in) :: id1, is1, ie1
    integer,  intent(in) :: id2, is2, ie2
    integer,  intent(in) :: id3, is3, ie3
    integer,  intent(in), dimension(id1, id2, id3) :: Ref

    integer :: i, j, k
    integer :: chksum
    chksum = 0

    do k = is3, ie3
      do j = is2, ie2
        do i = is1, ie1
          chksum = chksum + abs(Ref(i,j,k))
        end do
      end do
    end do

    print 666, name, chksum
666 format('chksum: ', a8, 1x, i8)

  end subroutine chksum_3d_I

end module cmparray_mod
