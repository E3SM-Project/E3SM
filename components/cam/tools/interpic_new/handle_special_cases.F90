subroutine handle_special_cases(ncidt, ncido)

   use shr_kind_mod, only: r8 => shr_kind_r8
   use control,      only: verbose, is_user_include

   implicit none

   include 'netcdf.inc'

   ! arguments
   integer, intent(in) :: ncidt
   integer, intent(in) :: ncido

   ! Local workspace
   integer :: ndims
   integer :: nvars
   integer :: ngatts
   integer :: unlimdimid

   character*(nf_max_name) :: name
   integer :: xtype
   integer :: nvdims
   integer :: vardids(nf_max_var_dims) ! variable dimension ids
   integer :: natts

   integer :: n
   integer :: vt, vo
   character*(nf_max_name) :: dim_name
   integer :: dimlen
   integer :: totsiz
   integer :: vardids_o(nf_max_var_dims)
   integer :: start(nf_max_var_dims)
   integer :: count(nf_max_var_dims)
   integer :: ret

   character*(nf_max_name) :: attname

   logical :: is_special_case

   ! Copy all the special case variables from template to output
   ! Probably should have separate "define" and "put" loops to optimize redef/enddef
   !
   ! ***N.B.*** This routine assumes that the required dimensions have already
   !            been defined.

   call wrap_inq(ncidt, ndims, nvars, ngatts, unlimdimid)

   start(:) = 1

   do vt = 1, nvars

      call wrap_inq_var(ncidt, vt, name, xtype, nvdims, vardids, natts)

      if (is_special_case(name) .or. is_user_include(name)) then

         totsiz = 1

         ! The dimension IDs from the template file need to be translated
         ! to the appropriate ID for the output file.
         ! ***N.B.*** Assume that the dimension names are the same on
         !            the template and output files
         do n = 1, nvdims

            ! Get dimension name and length from the template file
            if (nf_inq_dim(ncidt, vardids(n), dim_name, dimlen) /= NF_NOERR) then
               write(6,*)'handle_special_cases: ERROR from nf_inq_dim for variable='//&
                  trim(name)//', dimID=',vardids(n)
            end if

            ! Get dimension ID from the output file
            if (nf_inq_dimid(ncido, dim_name, vardids_o(n)) /= NF_NOERR) then
               write(6,*)'handle_special_cases: ERROR from nf_inq_dimid for dimension='//&
                  trim(dim_name)
            end if

            count(n) = dimlen
            totsiz = totsiz * dimlen

         end do

         ! special variables that have previously been defined won't be overwritten
         ! because this define request will fail.
         if( nf_def_var(ncido, name, xtype, nvdims, vardids_o, vo) == NF_NOERR) then

            ! Copy attributes from input to output, then the variable itself
            do n = 1, natts
               call wrap_inq_attname(ncidt, vt, n, attname)
               call wrap_copy_att(ncidt, vt, attname, ncido, vo)
            end do
           
            ! leave define mode
            if (nf_enddef(ncido) /= NF_NOERR) stop 999

            if(verbose) print *, "Copying special var ", trim(name)
                          
            call cpvar(ncidt, ncido, vt, vo, name, &
                       totsiz, xtype, start, count )

            ! return to define mode
            ret = nf_redef(ncido)
         else
            if(verbose) print *, "Cannot define special var ", trim(name)
         end if

      end if ! is special case

   end do ! loop over variables on input file

end subroutine handle_special_cases
