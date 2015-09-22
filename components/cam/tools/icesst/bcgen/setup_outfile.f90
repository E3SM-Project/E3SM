subroutine setup_outfile (ncid, outfil, tag, dateid, datesecid, &
                          timeid, icesst, nlon, nlat, xlon,     &
                          xlat, iyr1out, mon1out, history)
   use prec
   use types

   implicit none

   include 'netcdf.inc'

   integer, intent(out) :: ncid                   ! output file netcdf id
   character(len=*), intent(in) :: outfil         ! output file name
   character(len=*), intent(in) :: tag            ! 'CLIM' or 'AMIP'
   integer, intent(out) :: dateid, datesecid      ! date variable ids
   integer, intent(out) :: timeid                 ! time variable id
   type(icesstparms), intent(inout) :: icesst(2)  ! struct containing ice or sst-specific info
   integer, intent(in) :: nlon                    ! number of longitudes
   integer, intent(in) :: nlat                    ! number of latitudes
   real(r8), intent(in) :: xlon(nlon)             ! longitude coordinate values
   real(r8), intent(in) :: xlat(nlat)             ! latitude coordinate values
   integer, intent(in) :: iyr1out                 ! start year written to output file
   integer, intent(in) :: mon1out                 ! start month written to output file
   character(len=*), intent(in) :: history        ! history attribute

   integer :: ret                                 ! return code
   integer :: dimids(3)                           ! dimension ids for multi-dimensioned arrays
   integer :: lonid, latid                        ! lon, lat variable ids

   character(len=32) :: name                      ! variable name
   character(len=32) :: str                       ! attribute string
!
! Open file and define dimensions
!
   call wrap_nf_create (outfil, NF_CLOBBER, ncid)
   call wrap_nf_def_dim (ncid, 'lon', nlon, dimids(1))
   call wrap_nf_def_dim (ncid, 'lat', nlat, dimids(2))
   call wrap_nf_def_dim (ncid, 'time', NF_UNLIMITED, dimids(3))
!
! Model-readable date info
!
   call wrap_nf_def_var (ncid, 'date', NF_INT, 1, dimids(3), dateid)
   call wrap_nf_put_att_text (ncid, dateid, 'long_name', 'current date (YYYYMMDD)')

   call wrap_nf_def_var (ncid, 'datesec', NF_INT, 1, dimids(3), datesecid)
   call wrap_nf_put_att_text (ncid, datesecid, 'long_name', 'current seconds of current date')
!
! Dimension variables
!
   call wrap_nf_def_var (ncid, 'lon', NF_DOUBLE, 1, dimids(1), lonid)
   call wrap_nf_put_att_text (ncid, lonid, 'long_name', 'longitude')
   call wrap_nf_put_att_text (ncid, lonid, 'units', 'degrees_east')

   call wrap_nf_def_var (ncid, 'lat', NF_DOUBLE, 1, dimids(2), latid)
   call wrap_nf_put_att_text (ncid, latid, 'long_name', 'latitude')
   call wrap_nf_put_att_text (ncid, latid, 'units', 'degrees_north')

   call wrap_nf_def_var (ncid, 'time', NF_DOUBLE, 1, dimids(3), timeid)
   write(str,100) iyr1out, mon1out
100 format('days since ',i4.4,'-',i2.2,'-01 00:00:00')
   call wrap_nf_put_att_text (ncid, timeid, 'units', trim(str))
   call wrap_nf_put_att_text (ncid, timeid, 'calendar', '365_day')

   ! history attribute
   call wrap_nf_put_att_text (ncid, NF_GLOBAL, 'history', trim(history))

!
! Define variables
!
   select case (tag)
   case ('CLIM')
!
! ice
!
      name = icesst(1)%climname
      call wrap_nf_def_var (ncid, name, NF_FLOAT, 3, dimids, icesst(1)%varidclim)
      call wrap_nf_put_att_text (ncid, icesst(1)%varidclim, 'long_name', icesst(1)%climlongname)
      call wrap_nf_put_att_text (ncid, icesst(1)%varidclim, 'units', icesst(1)%units)

      name = icesst(1)%climname_pre
      call wrap_nf_def_var (ncid, name, NF_FLOAT, 3, dimids, icesst(1)%varidclim_pre)
      call wrap_nf_put_att_text (ncid, icesst(1)%varidclim_pre, 'long_name', icesst(1)%climlongname_pre)
      call wrap_nf_put_att_text (ncid, icesst(1)%varidclim_pre, 'units', icesst(1)%units)
!
! sst
!
      name = icesst(2)%climname
      call wrap_nf_def_var (ncid, name, NF_FLOAT, 3, dimids, icesst(2)%varidclim)
      call wrap_nf_put_att_text (ncid, icesst(2)%varidclim, 'long_name', icesst(2)%climlongname)
      call wrap_nf_put_att_text (ncid, icesst(2)%varidclim, 'units', icesst(2)%units)

      name = icesst(2)%climname_pre
      call wrap_nf_def_var (ncid, name, NF_FLOAT, 3, dimids, icesst(2)%varidclim_pre)
      call wrap_nf_put_att_text (ncid, icesst(2)%varidclim_pre, 'long_name', icesst(2)%climlongname_pre)
      call wrap_nf_put_att_text (ncid, icesst(2)%varidclim_pre, 'units', icesst(2)%units)
   case ('AMIP')
!
! ice
!
      name = icesst(1)%amipname
      call wrap_nf_def_var (ncid, name, NF_FLOAT, 3, dimids, icesst(1)%varidamip)
      call wrap_nf_put_att_text (ncid, icesst(1)%varidamip, 'long_name', icesst(1)%amiplongname)
      call wrap_nf_put_att_text (ncid, icesst(1)%varidamip, 'units', icesst(1)%units)
      
      name = icesst(1)%amipname_pre
      call wrap_nf_def_var (ncid, name, NF_FLOAT, 3, dimids, icesst(1)%varidamip_pre)
      call wrap_nf_put_att_text (ncid, icesst(1)%varidamip_pre, 'long_name', icesst(1)%amiplongname_pre)
      call wrap_nf_put_att_text (ncid, icesst(1)%varidamip_pre, 'units', icesst(1)%units)
!
! sst
!
      name = icesst(2)%amipname
      call wrap_nf_def_var (ncid, name, NF_FLOAT, 3, dimids, icesst(2)%varidamip)
      call wrap_nf_put_att_text (ncid, icesst(2)%varidamip, 'long_name', icesst(2)%amiplongname)
      call wrap_nf_put_att_text (ncid, icesst(2)%varidamip, 'units', icesst(2)%units)

      name = icesst(2)%amipname_pre
      call wrap_nf_def_var (ncid, name, NF_FLOAT, 3, dimids, icesst(2)%varidamip_pre)
      call wrap_nf_put_att_text (ncid, icesst(2)%varidamip_pre, 'long_name', icesst(2)%amiplongname_pre)
      call wrap_nf_put_att_text (ncid, icesst(2)%varidamip_pre, 'units', icesst(2)%units)
   case default
      write(6,*) 'tag ', tag, ' is an invalid tag: valid values are CLIM and AMIP'
      stop 999
   end select

   ret = nf_enddef (ncid)
   if (ret < 0) then
      call handle_error (ret)
   end if

   call wrap_nf_put_var_double (ncid, lonid, xlon)
   call wrap_nf_put_var_double (ncid, latid, xlat)

   return
end subroutine setup_outfile
   
   
