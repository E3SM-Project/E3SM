subroutine driver(ncidi, ncido, ncidt, ntime)

   use shr_kind_mod,   only: r8 => shr_kind_r8
   use control,        only: verbose, prec_override, prec_out, is_user_skip
   use dimensions,     only: ncoldimid, get_shape
   use fill_positions, only: varspec_t, fillvar
   use interp,         only: interp_driver

   implicit none

   include 'netcdf.inc'

   ! Input arguments

   integer ncidi, ncido, ncidt     ! input, output, template netcdf file ids
   integer ntime                   ! size of unlimited dimension (if present) in "input" file

   ! Local workspace

   character*8 :: shapei                   ! dimension order of variable on input file
   character*8 :: shapeo                   ! dimension order of variable on output file
   character*(nf_max_name) :: name         ! name of variable on input and output file
   character*(nf_max_name) :: dimnamesi(4) ! dimension names of a variable on the input file
   character*(nf_max_name) :: dimnames(4)  ! dimension names for output variable
   character*(nf_max_name) :: attname

   integer :: ndims, nvars, ngatts, unlimdimid

   integer :: natts                        ! number of attributes for a given variable
   integer :: nvdims, nvdimsi              ! number of dimensions for this variable
   integer :: vardids_i(nf_max_var_dims) ! variable dimension id's
   integer :: vardids_o(nf_max_var_dims)   ! variable dimension id's on output file
   integer :: j, k       ! spatial indices
   integer :: n      ! index
   integer :: t      ! index over unlimited dimension
   integer :: v, vi, vo            ! loop index over variable id
   integer :: vidi                 ! variable id on input file
   integer :: vido                 ! returned variable id on output file
   integer :: nxi, nyi, nzi
   integer :: nxo, nyo, nzo
   integer :: xtype, xtypeo            ! variable type (netcdf)
   integer :: tpos                     ! position of unlimited dimension
   integer :: ncp   = 0                ! number of variables to copy from input to output
   integer :: nintp = 0                ! number of variables that require interpolation
   integer :: starti(nf_max_var_dims)
   integer :: starto(nf_max_var_dims)
   integer :: counto(nf_max_var_dims)
   integer :: counti(nf_max_var_dims)
   integer :: totsiz
   integer, allocatable :: indx_cp(:)           ! input file var IDs for variables to be copied
   integer, allocatable :: indx_intp(:)         ! input file var IDs for variables to be interpolated
   integer :: i
   logical :: copy

   real(r8), target, allocatable :: arr3d_i(:,:,:)
   real(r8), pointer     :: arr3d_o(:,:,:)
   real(r8), pointer     :: arrxyzi(:,:,:)
   real(r8), target, allocatable :: arrxyzo(:,:,:)
   real,     allocatable :: arr3d_o_r4(:,:,:)

   type(varspec_t), allocatable :: varspec_i(:)
   type(varspec_t), allocatable :: varspec_o(:)

   logical :: is_special_case
   !-----------------------------------------------------------------------------

   ! loop over variables in "input" dataset
   call wrap_inq(ncidi, ndims, nvars, ngatts, unlimdimid)

   allocate(indx_cp(nvars))
   allocate(indx_intp(nvars))
   allocate(varspec_i(nvars))
   allocate(varspec_o(nvars))
   indx_cp(:)   = -1
   indx_intp(:) = -1

   do vidi = 1, nvars

      copy = .true.
    
      call wrap_inq_var(ncidi, vidi, name, xtype, nvdimsi, vardids_i, natts)

      ! Skip any special case variables: they have already been copied to the output file.
      if (is_special_case(name)) then
         if (verbose) write(6,*)'skipping special case var ',trim(name)
         cycle
      end if

      ! Skip any variables requested by user
      if (is_user_skip(name)) then
         if (verbose) write(6,*)'skipping user requested var ',trim(name)
         cycle
      end if

      ! return name of each dimension (dimnamesi), and the axis ordering (shapei)
      shapei = get_shape(ncidi, vardids_i, nvdimsi, dimnamesi)

      vardids_o = 0
      nvdims    = nvdimsi
      shapeo    = ' '

      ! Set arrays containing information needed to define the variable
      ! in the output file.

      if (ncoldimid > 0 .and. (shapei(1:2)=='xz' .or. shapei(1:2)=='xy')) then

         ! If the template file contains an "ncol" dimension then, then the output
         ! file uses an unstructured grid.  If the input file contains variables on
         ! a rectangular grid, then the 'xy' dimensions need to be compressed to a
         ! single column dimension 'n'.

         if (shapei .eq. 'xzy') then

            dimnames(1)  = 'ncol'
            vardids_o(1) = ncoldimid
            dimnames(2)  = dimnamesi(2)
            vardids_o(2) = get_dimid_o(dimnamesi(2))

            if (nvdimsi > 3) then
               dimnames(3)  = dimnamesi(4)
               vardids_o(3) = get_dimid_o(dimnamesi(4))
            end if

            nvdims = nvdimsi-1
            shapeo = 'nz'

         else if (shapei .eq. 'xyz') then

            dimnames(1)  = 'ncol'
            vardids_o(1) = ncoldimid
            dimnames(2)  = dimnamesi(3)
            vardids_o(2) = get_dimid_o(dimnamesi(3))

            if (nvdimsi > 3) then
               dimnames(3)  = dimnamesi(4)
               vardids_o(3) = get_dimid_o(dimnamesi(4))
            end if

            nvdims = nvdimsi-1
            shapeo = 'nz'

         else if(shapei .eq. 'xy') then

            dimnames(1)  = 'ncol'
            vardids_o(1) = ncoldimid

            if (nvdimsi > 2) then
               dimnames(2)  = dimnamesi(3)
               vardids_o(2) = get_dimid_o(dimnamesi(3))
            end if

            nvdims = nvdimsi - 1
            shapeo = 'n'

         end if

      else

         ! Assume that the dimension names on the output file match the
         ! dimension names on the input file.  Get the dimension ID from the 
         ! output file since that may be different from the ID in the input
         ! file.
         do i = 1, nvdimsi
            vardids_o(i) = get_dimid_o(dimnamesi(i))
         end do

         dimnames = dimnamesi
         nvdims   = nvdimsi
         shapeo   = shapei

      end if

      ! Allow override of real precision for output fields on 2D or 3D spatial grid
      xtypeo = xtype
      if (shapei(1:2) == 'xy' .or. shapei(1:2) == 'xz') then
         if (prec_override) then
            xtypeo = prec_out
         end if
      end if

      ! Define variable on output file with same type and dimension order as variable
      ! on input file
      if (nf_def_var(ncido, name, xtypeo, nvdims, vardids_o, vido) /= nf_noerr) then
         print*, 'driver: ERROR defining variable ', trim(name), ' type=', xtypeo
         print*, 'nvdims=', nvdims, ' dids=',vardids_o(1:nvdims)
         stop
         !      call wrap_def_var (ncido, name, xtypeo, nvdims, vardids_o, vido)
      else
         print*, 'create var ', trim(name), ' ndims=', nvdims, ' dids=', vardids_o(1:nvdims), ' vid=', vido
      end if

      ! Copy variable's attributes from input to output
      do n = 1, natts
         call wrap_inq_attname (ncidi, vidi, n, attname)
         call wrap_copy_att (ncidi, vidi, attname, ncido, vido)
      end do

      if (xtype == NF_FLOAT .or. xtype == NF_DOUBLE) then

         ! Interpolated variables must be of a floating point type and have dimensions
         ! xy, xyz, or xzy.
         ! Assumption: All variables using a horizontal grid need to be interpolated.
         if (shapeo == 'xy' .or. shapeo == 'xyz' .or. shapeo == 'xzy' .or. shapeo == 'n' .or. shapeo =='nz') then

            copy = .false.
            nintp = nintp + 1
            indx_intp(nintp) = vidi
         end if
      end if

      ! Variables to be copied will not be interpolated.
      if (copy) then
         ncp = ncp + 1
         indx_cp(ncp) = vidi
      end if

      ! Set the varspec_i array with coordinate information for each variable
      ! from the input file.
      call fillvar(ncidi, name, xtype, shapei, dimnamesi, &
                   vidi, nvdimsi, vardids_i, varspec_i(vidi))

      
      ! Set the varspec_o array with coordinate information for each variable
      ! from the output file.
      call fillvar(ncido, name, xtypeo, shapeo, dimnames, &
                   vido, nvdims, vardids_o, varspec_o(vidi))

      call compare_var(varspec_i(vidi), varspec_o(vidi))

      shapeo=''
      shapei=''
      dimnames=''

   end do         ! loop over input variables

   ! Take the output file out of define mode
   if (nf_enddef(ncido) /= NF_NOERR) then
      print *,'ERROR:  nf_enddef(ncido) = ',nf_enddef(ncido)
      stop
   endif

   ! Now loop over the time dimension.  First do copies
   if (verbose) write(6,*)'Copy and interpolate from input to output.  ntime=',ntime
   do t = 1, ntime

      if (verbose) write(6,*)'Starting time sample ',t

      do n = 1, ncp
         v           = indx_cp(n)
         name        = varspec_i(v)%name
         totsiz      = varspec_i(v)%totsiz
         xtype       = varspec_i(v)%xtype
         vi          = varspec_i(v)%varid
         vo          = varspec_o(v)%varid
         starti(:)   = 1
         counto(:)    = varspec_i(v)%count(:)

         if (varspec_i(v)%tpos > 0) then
            tpos      = varspec_i(v)%tpos
            starti(tpos) = t
         end if

         if (verbose) then
            write(6,*)'Copying variable ',trim(name), vi, vo
         end if

         call cpvar(ncidi, ncido, vi, vo, name, &
                    totsiz, xtype, starti, counto)
      end do

      ! Now the data which need to be interpolated

      do n = 1, nintp

         v           = indx_intp(n)

         name        = varspec_i(v)%name

         if (verbose) then
            write(6,*) 'Interpolating '//trim(name)
         end if

         vi          = varspec_i(v)%varid
         shapei      = varspec_i(v)%vshape
         nxi         = varspec_i(v)%nx
         nyi         = varspec_i(v)%ny
         nzi         = varspec_i(v)%nz
         counti(:)   = varspec_i(v)%count(:)

         if (verbose) then
            write(6,'(a16,a3,3(a5,i6))') '  input:  shape=', trim(shapei), '  nx=', nxi, '  ny=', nyi, '  nz=', nzi
         end if

         vo          = varspec_o(v)%varid
         shapeo      = varspec_o(v)%vshape
         nxo         = varspec_o(v)%nx
         nyo         = varspec_o(v)%ny
         nzo         = varspec_o(v)%nz
         counto(:)    = varspec_o(v)%count

         if (verbose) then
            write(6,'(a16,a3,3(a5,i6))') '  output: shape=', trim(shapeo), '  nx=', nxo, '  ny=', nyo, '  nz=', nzo
         end if


         starti(:)    = 1
         starto(:)    = 1

         if (varspec_i(v)%tpos > 0) then
            tpos      = varspec_i(v)%tpos
            if(shapeo(1:1).eq.'n') then
               starto(tpos-1) = t
            else
               starto(tpos) = t
            end if
            starti(tpos) = t
         end if

         allocate(arrxyzo(nxo,nyo,nzo))
         if (trim(shapei) == 'xzy') then
            allocate(arr3d_i(nxi,nzi,nyi))
            allocate(arrxyzi(nxi,nyi,nzi))
         else
            allocate(arr3d_i(nxi,nyi,nzi))
         end if

         call wrap_get_vara_double(ncidi, vi, starti, counti, arr3d_i)

         if (trim(shapei) == 'xzy') then
            do k = 1, nzi
               do j = 1, nyi
                  do i = 1, nxi
                     arrxyzi(i,j,k) = arr3d_i(i,k,j)
                  end do
               end do
            end do
         else
            arrxyzi => arr3d_i
         end if

         if (verbose) then
            write(6,*)'  input min, avg, max = ', minval(arrxyzi), &
                      sum(arrxyzi)/real(size(arrxyzi), r8), maxval(arrxyzi)
         end if

         call interp_driver(ncidi, ncido, arrxyzi, arrxyzo, varspec_i(v), varspec_o(v), nxi, &
                            nzi, nyi, nxo, nzo, nyo)

         if (verbose) then
            write(6,*)'  output min, avg, max= ', minval(arrxyzo), &
                      sum(arrxyzo)/real(size(arrxyzo), r8), maxval(arrxyzo)
         end if

         ! Output uses same dimension ordering as input
         if (trim(shapei) == 'xzy') then
            allocate(arr3d_o(nxo,nzo,nyo))
            do k = 1, nzo
               do j = 1, nyo
                  do i = 1, nxo
                     arr3d_o(i,k,j) = arrxyzo(i,j,k)
                  end do
               end do
            end do
         else
            arr3d_o => arrxyzo
         end if

         if (varspec_o(v)%xtype == NF_DOUBLE) then
            call wrap_put_vara_double(ncido, vo, starto, counto, arr3d_o)
         else

            ! Convert to 4-byte reals before output
            if (trim(shapei) == 'xzy') then
               allocate(arr3d_o_r4(nxo,nzo,nyo))
            else
               allocate(arr3d_o_r4(nxo,nyo,nzo))
            end if

            arr3d_o_r4 = arr3d_o
            if (nf_put_vara_real(ncido, vo, starto, counto, arr3d_o_r4) /= NF_NOERR) then
               write(6,*) 'driver: ERROR from nf_put_vara_real'
               stop
            end if
            deallocate(arr3d_o_r4)
         end if

         deallocate(arr3d_i)
         deallocate(arrxyzo)
         if (trim(shapei) == 'xzy') then
            deallocate(arrxyzi)
            deallocate(arr3d_o)
         end if

      end do

   end do

   deallocate(indx_cp)
   deallocate(indx_intp)
   deallocate(varspec_i)
   deallocate(varspec_o)

contains

integer function get_dimid_o(dimname)

! Lookup the dimension ID for the requested name in the output file

   character(len=*), intent(in) :: dimname

   if (nf_inq_dimid(ncido, trim(dimname), get_dimid_o) /= nf_noerr) then
      print *, 'driver: ERROR: input variable ', trim(name), ' has dimension ', trim(dimname), &
               'which is not found in the output file'
      stop
   end if
end function get_dimid_o

end subroutine driver
