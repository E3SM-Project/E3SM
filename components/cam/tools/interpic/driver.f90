subroutine driver (ncidi, ncido, ncidt, nvars, ntime)
!
! $Id$
!
  use shr_kind_mod, only: r8 => shr_kind_r8
  use varspecs_mod, only: varspecs
  use control
  use dimensions
  use fill_positions

  implicit none

  include 'netcdf.inc'
!
! Input arguments
!
  integer ncidi, ncido, ncidt     ! input, output, template netcdf file ids
  integer nvars                   ! number of variables
  integer ntime                   ! size of unlimited dimension (if present)
!
! Local workspace
!
  character*8 :: shapei, shapeo
  character*(nf_max_name) :: name, namei
  character*(nf_max_name) :: dimnames(3)
  character*(nf_max_name) :: dimnamesi(3)
  character*(nf_max_name) :: attname

  integer natts, nattsi               ! number of attributes for a given variable
  integer nvdims, nvdimsi              ! number of dimensions for this variable
  integer vardids(nf_max_var_dims) ! variable dimension id's
  integer vardidsi(nf_max_var_dims) ! variable dimension id's
  integer j, k       ! spatial indices
  integer n      ! index
  integer t      ! index over unlimited dimension
  integer v                    ! loop index over variable id
  integer vi, vo         ! returned variable id on output file
  integer nxi, nyi, nzi
  integer nxo, nyo, nzo
  integer xtype, xtypei                ! variable type (netcdf)
  integer tpos               ! position of unlimited dimension
  integer ncp
  integer nintp
  integer start(nf_max_var_dims)
  integer count(nf_max_var_dims)
  integer counti(nf_max_var_dims)
  integer totsiz
  integer :: indx_cp(nvars)
  integer :: indx_intp(nvars)

  logical copy

  real(r8), allocatable :: arrxyzi(:,:,:), arrxzyi(:,:,:)
  real(r8), allocatable :: arrxyzo(:,:,:), arrxzyo(:,:,:)

  type (varspecs) :: vari(nvars)
  type (varspecs) :: varo(nvars)

  logical is_special_case
  external is_special_case
!
! Initialize indices to invalid values
!
  indx_cp(:) = -1
  indx_intp(:) = -1

  ncp = 0
  nintp = 0


  do v=1,nvars
    copy = .true.
    
    call wrap_inq_var (ncidt, v, name, xtype, nvdims, vardids, natts)
!
! Skip any special case variables: they have already been dealt with.
! Also skip the variable if it is not on the input tape
!
    if (is_special_case (name, ncidt)) then
      if (verbose) write(6,*)'skipping special case var ',trim(name)
      cycle
    end if

    if (nf_inq_varid (ncidi, name, vi) == nf_noerr) then
      call wrap_inq_var (ncidi, vi, namei, xtypei, nvdimsi, vardidsi, nattsi)
    else
      if (verbose) write(6,*)trim(name),' not found on input: skipping'
      cycle
    end if
!
! Variable is on both input and template file, and is not a special case
!
    call wrap_def_var (ncido, name, xtype, nvdims, vardids, vo)
!
! Copy attributes from input to output
!
    do n=1,nattsi
      call wrap_inq_attname (ncidi, vi, n, attname)
      call wrap_copy_att (ncidi, vi, attname, ncido, vo)
    end do

    if (xtype == NF_FLOAT .or. xtype == NF_DOUBLE) then

      shapeo = get_shape (ncidt, vardids, nvdims, dimnames)
!
! Interpolated variables must be of a floating point type and have dimensions
! xy, xyz, or xzy.
!
      if (shapeo == 'xy' .or. shapeo == 'xyz' .or. shapeo == 'xzy') then
        if (xtypei /= NF_FLOAT .and. xtypei /= NF_DOUBLE) then
          write(6,*)'driver: type of ', trim(namei), &
                    ' does not match between input and template files'
          stop 999
        end if

        shapei = get_shape (ncidi, vardidsi, nvdimsi, dimnamesi)

        copy = .false.
        nintp = nintp + 1
        indx_intp(nintp) = v
      end if
    end if
!
! Variables to be copied will not be interpolated.  Determine which do and
! do not have an unlimited dimension
!
    if (copy) then
      ncp = ncp + 1
      indx_cp(ncp) = v
    end if
!
! Copy useful information for copying or interpolating into the struct
!
    call fillvar (ncidi, namei, xtypei, shapei, dimnamesi, &
                  ewdimi, nsdimi, zdimi, vari(v), vi, &
                  nvdimsi, vardidsi)

    call fillvar (ncidt, name, xtype, shapeo, dimnames, &
                  ewdim,  nsdim,  zdim,  varo(v), vo, &
                  nvdims, vardids)

    call compare_var (vari(v), varo(v))
  end do                           ! loop over input variables

  if (nf_enddef (ncido) /= NF_NOERR) stop 999
!
! Now loop over the unlimited dimension.  First do copies
!
  do t=1,ntime
    if (.not.silent) then
      write(6,*)'Starting time sample ',t
    end if

    do n=1,ncp
      v           = indx_cp(n)

      name        = vari(v)%name
      totsiz      = vari(v)%totsiz
      xtype       = vari(v)%xtype
      vi          = vari(v)%varid
      vo          = varo(v)%varid

      start(:)    = 1
      count(:)    = vari(v)%count(:)

      if (vari(v)%tpos > 0) then
        tpos      = vari(v)%tpos
        start(tpos) = t
      end if

      call cpvar (ncidi, ncido, vi, vo, name, &
                  totsiz, xtype, start, count)
    end do
!
! Now the data which need to be interpolated
!
    do n=1,nintp
      v           = indx_intp(n)

      name        = vari(v)%name
      vi          = vari(v)%varid
      shapei      = vari(v)%vshape
      nxi         = vari(v)%nx
      nyi         = vari(v)%ny
      nzi         = vari(v)%nz
      counti(:)   = vari(v)%count(:)

      vo          = varo(v)%varid
      shapeo      = varo(v)%vshape
      nxo         = varo(v)%nx
      nyo         = varo(v)%ny
      nzo         = varo(v)%nz
      count(:)    = varo(v)%count(:)

      start(:)    = 1

      if (vari(v)%tpos > 0) then
        tpos      = vari(v)%tpos
        start(tpos) = t
      end if

      allocate (arrxzyi(nxi,nzi,nyi))
      allocate (arrxzyo(nxo,nzo,nyo))

      if (shapei(1:2) == 'xy') then

        allocate (arrxyzi(nxi,nyi,nzi))

        call wrap_get_vara_double (ncidi, vi, start, counti, arrxyzi)

        do j=1,nyi
          do k=1,nzi
            arrxzyi(:,k,j) = arrxyzi(:,j,k)
          end do
        end do

        deallocate (arrxyzi)

      else

        call wrap_get_vara_double (ncidi, vi, start, counti, arrxzyi)

      end if

      if (verbose) then
        write(6,*)'Interpolating ',trim(name),' 1st elem=',arrxzyi(1,1,1)
      end if

      call interp_driver (arrxzyi, arrxzyo, vari(v), varo(v), nxi, &
                          nzi, nyi, nxo, nzo, nyo)

      if (shapeo(1:2) == 'xy') then

        allocate (arrxyzo(nxo,nyo,nzo))

        do j=1,nyo
          do k=1,nzo
            arrxyzo(:,j,k) = arrxzyo(:,k,j)
          end do
        end do

        call wrap_put_vara_double (ncido, vo, start, count, arrxyzo)
        deallocate (arrxyzo)

      else if (shapeo == 'xzy') then

        call wrap_put_vara_double (ncido, vo, start, count, arrxzyo)

      else

        write(6,*)'Unknown shape=',shapeo,' for variable ',name
        stop 999

      end if

      deallocate (arrxzyi)
      deallocate (arrxzyo)

    end do
  end do    

  return
end subroutine driver
