module fill_positions

  contains

  subroutine fill_xpos (ncid, dimname, nx, ny, xpos, numx)
    use shr_kind_mod, only: r8 => shr_kind_r8

    implicit none

    include 'netcdf.inc'
!
! Input arguments
!
    integer ncid
    integer nx, ny
    integer numx(ny)
    character*(nf_max_name) dimname
    real(r8) :: xpos(nx,ny)
!
! Local workspace
!
    integer varid
    integer j

    if (dimname == 'lon' .and. nf_inq_varid (ncid, 'rlon', varid) == nf_noerr) then
      call wrap_get_var_double (ncid, varid, xpos)
      call wrap_inq_varid (ncid, 'nlon', varid)
      call wrap_get_var_int (ncid, varid, numx)
    else
      call wrap_inq_varid (ncid, dimname, varid)
      call wrap_get_var_double (ncid, varid, xpos(1,1))
      do j=2,ny
        xpos(:,j) = xpos(:,1)
      end do
      numx(:) = nx
    end if

    return
  end subroutine fill_xpos
!-------------------------------------------------------------------------------
  subroutine fill_yzpos (ncid, dimname, nyz, yzpos)
    use shr_kind_mod, only: r8 => shr_kind_r8

    implicit none

    include 'netcdf.inc'
!
! Input arguments
!
    integer ncid
    integer nyz
    character*(nf_max_name) dimname
    real(r8) :: yzpos(nyz)
!
! Local workspace
!
    integer varid

    call wrap_inq_varid (ncid, dimname, varid)
    call wrap_get_var_double (ncid, varid, yzpos)

    return
  end subroutine fill_yzpos
!-------------------------------------------------------------------------------
  subroutine fillvar (ncid, name, xtype, vshape, dimnames, &
                      ewdim, nsdim, zdim, var, varid, &
                      nvdims, vardids)

    use dimensions, only: info, get_dimlen, maxdims
    use varspecs_mod, only: varspecs

    implicit none

    include 'netcdf.inc'
!
! Input arguments
!
    integer ncid
    integer varid
    character*(*) name
    integer xtype
    character*(nf_max_name) dimnames(3)
    character*8 vshape
    type(info) :: ewdim(maxdims), nsdim(maxdims), zdim(maxdims)
    type(varspecs) :: var
    integer nvdims
    integer vardids(*)
!
! Local workspace
!
    integer nx, ny, nz
    integer dimlen
    integer n
    integer ret
    integer unlimdimid

    var%name   = name
    var%xtype  = xtype
    var%vshape = vshape
    var%varid  = varid

    var%totsiz   = 1    ! ini to to size which can be multiplied
    var%tpos     = 1    ! init to 1 in case unlimited dim not present
    var%count(:) = 1    ! init to 1 in case unlimited dim not present

    nx    = 1 
    ny    = 1 
    nz    = 1

    if (vshape == 'xy') then

      nx = get_dimlen (ewdim, dimnames(1))
      ny = get_dimlen (nsdim, dimnames(2))

      allocate(var%numx(ny))
      allocate(var%xpos(nx,ny))
      allocate(var%ypos(ny))

      call fill_xpos (ncid, dimnames(1), nx, ny, var%xpos, var%numx)
      call fill_yzpos (ncid, dimnames(2), ny, var%ypos)

    else if (vshape == 'xyz') then

      nx = get_dimlen (ewdim, dimnames(1))
      ny = get_dimlen (nsdim, dimnames(2))
      nz = get_dimlen (zdim, dimnames(3))

      allocate(var%numx(ny))
      allocate(var%xpos(nx,ny))
      allocate(var%ypos(ny))
      allocate(var%zpos(nz))

      call fill_xpos (ncid, dimnames(1), nx, ny, var%xpos, var%numx)
      call fill_yzpos (ncid, dimnames(2), ny, var%ypos)
      call fill_yzpos (ncid, dimnames(3), nz, var%zpos)

    else if (vshape == 'xzy') then

      nx = get_dimlen (ewdim, dimnames(1))
      ny = get_dimlen (nsdim, dimnames(3))
      nz = get_dimlen (zdim, dimnames(2))

      allocate(var%numx(ny))
      allocate(var%xpos(nx,ny))
      allocate(var%ypos(ny))
      allocate(var%zpos(nz))

      call fill_xpos (ncid, dimnames(1), nx, ny, var%xpos, var%numx)
      call fill_yzpos (ncid, dimnames(3), ny, var%ypos)
      call fill_yzpos (ncid, dimnames(2), nz, var%zpos)

    end if

    var%nx = nx
    var%ny = ny
    var%nz = nz

    ret = nf_inq_unlimdim (ncid, unlimdimid)

    do n=1,nvdims
      if (vardids(n) == unlimdimid) then
        var%tpos = n
        var%count(n) = 1
      else
        call wrap_inq_dimlen (ncid, vardids(n), dimlen)
        var%count(n) = dimlen
        var%totsiz = var%totsiz * dimlen
      end if
    end do

    return
  end subroutine fillvar

end module fill_positions

