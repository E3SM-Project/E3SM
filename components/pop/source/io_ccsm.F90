!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module io_ccsm

!BOP
! !MODULE: io_ccsm
!
! !DESCRIPTION:
!  This module provides a kludge interface for writing nonstandard
!  ccsm fields (eg, variables that are not on lat/lon grids and
!  thus are unable to be defined via the construct_io_field function)
!  to ccsm netCDF output files
!
! !REVISION HISTORY:
!  SVN:$Id: io_ccsm.F90 27495 2011-03-25 21:46:00Z njn01 $
!  
! 

! !USES:

   use kinds_mod
   use blocks
   use communicate
   use broadcast
   use exit_mod
   use domain
   use constants
   use io_netcdf
   use io_binary
   use io_types
   use io_tools

   implicit none
   public  ! to get io_types without having to explicitly use io_types
           ! module directly
   save

! !PUBLIC MEMBER FUNCTIONS:

   public ::            &
     data_set_nstd_ccsm

!EOP
!BOC
 

!EOC
!***********************************************************************

contains

!***********************************************************************
!BOP
! !IROUTINE: data_set_nstd_ccsm
! !INTERFACE:

 subroutine data_set_nstd_ccsm (data_file,operation,field_id, &
                                ndims,io_dims,nftype,         &
                                short_name,long_name,units,   &
                                time_dim,                     &
                                coordinates,                  &
                                fill_value,                   &
                                method_string,                &
                                data_1d_r8,                   &
                                data_2d_r8,                   &
                                data_2d_r4,                   &
                                data_3d_r4,                   &
                                data_4d_r4,                   &
                                data_1d_ch,                   &
                                data_2d_ch                    )

! !DESCRIPTION:
!  This routine is kludge interface to defining and writing the nonstandard
!  ccsm fields (eg, those not constructed via construct_io_field) 
!  to a netCDF file
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      ndims

   character (*), intent(in)  ::  &
      operation,                  &
      short_name,                 &
      long_name,                  &
      units,                      &
      coordinates,                &
      nftype,                     &
      method_string

   real (r4), dimension (:,:,:,:), intent(in) ::  &
      data_4d_r4
   real (r4), dimension (:,:,:),   intent(in) ::  &
      data_3d_r4
   real (r4), dimension (:,:),     intent(in) ::  &
      data_2d_r4

   real (r8), dimension (:,:), intent(in) ::  &
      data_2d_r8
   real (r8), dimension (:),   intent(in) ::  &
      data_1d_r8

   character (*), dimension (:,:), intent(in) ::  &
      data_2d_ch
   character (*), dimension (:),   intent(in) ::  &
      data_1d_ch


   real (r4), intent(in)  ::  &
      fill_value

   type (datafile),intent(inout) ::  &
      data_file

! !INPUT/OUTPUT PARAMETERS:
   integer (i4), intent(inout)   ::  &
      field_id  

   type (io_dim), intent(inout)  ::  &
      io_dims(:)

   type (io_dim), optional ::   &
      time_dim         ! dimension descriptor for (unlimited) time dim

   optional ::            & 
      short_name,         &
      long_name,          &
      units,              &
      coordinates,        &
      fill_value,         &
      method_string,      &
      data_1d_r8,         &
      data_2d_r8,         &
      data_2d_r4,         &
      data_3d_r4,         &
      data_4d_r4,         &
      data_1d_ch,         &
      data_2d_ch

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   logical (log_kind) :: supported
   logical (log_kind) :: lactive_time_dim

!-----------------------------------------------------------------------
!
!  Must be netCDF format -- binary version does not exist
!
!-----------------------------------------------------------------------

   if (data_file%data_format=='bin') then
         call exit_POP(sigAbort, &
             '(data_set_nstd_ccsm) ERROR: cannot call this routine with bin format')
   endif

!-----------------------------------------------------------------------
!
!  Deal with optional time dimension
!
!-----------------------------------------------------------------------

   if (present (time_dim)) then
     if (time_dim%active) then
      lactive_time_dim = .true.
     else
      lactive_time_dim = .false.
     endif
   else
      lactive_time_dim = .false.
   endif

!-----------------------------------------------------------------------
!
!  select operation to perform
!
!-----------------------------------------------------------------------

   select case (trim(operation))

!-----------------------------------------------------------------------
!
!  define an io field
!
!-----------------------------------------------------------------------

   case ('define')

      !*** note: should test for presence of optional variables
      call define_nstd_netcdf(data_file, ndims, io_dims,&
                              field_id,                             &
                              short_name, long_name, units,         &
                              coordinates=coordinates,              &
                              fill_value=fill_value,                &
                              method_string=method_string,          &
                              nftype=nftype                         )

!-----------------------------------------------------------------------
!
!  write an io field
!
!-----------------------------------------------------------------------

   case ('write')

      call write_nstd_netcdf(                     &
           data_file, field_id,                   &
           ndims, io_dims,nftype,                 &
           lactive_time_dim=lactive_time_dim,     &
           indata_1d_r8=data_1d_r8,               & ! pass all data arrays 
           indata_2d_r8=data_2d_r8,               & ! to write_nstd_netcdf
           indata_2d_r4=data_2d_r4,               &
           indata_3d_r4=data_3d_r4,               &
           indata_4d_r4=data_4d_r4,               &
           indata_1d_ch=data_1d_ch,               &
           indata_2d_ch=data_2d_ch                )

!-----------------------------------------------------------------------
!
!  unknown operation
!
!-----------------------------------------------------------------------

   case default

      if (my_task == master_task) &
         write(stdout,*) 'data_set_nstd_ccsm operation: ',trim(operation)
      call exit_POP(sigAbort,'data_set_nstd_ccsm: Unknown operation')

   end select

!-----------------------------------------------------------------------
!EOC

 end subroutine data_set_nstd_ccsm

!***********************************************************************

 end module io_ccsm

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
