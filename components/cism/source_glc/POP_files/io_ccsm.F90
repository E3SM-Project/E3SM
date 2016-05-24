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
!  SVN:$Id: io_ccsm.F90 808 2006-04-28 17:06:38Z njn01 $
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
                                coordinates,missing_value,    &
                                fill_value,                   &
                                implied_time_dim,             &
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
      nftype

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
      fill_value,             &
      missing_value

   type (datafile),intent(inout) ::  &
      data_file

! !INPUT/OUTPUT PARAMETERS:
   integer (i4), intent(inout)   ::  &
      field_id  

   type (io_dim), intent(inout)  ::  &
      io_dims(:)
   
   logical (log_kind), intent(inout) ::  &
      implied_time_dim

   optional ::            & 
      implied_time_dim,   &
      short_name,         &
      long_name,          &
      units,              &
      coordinates,        &
      missing_value,      &
      fill_value,         &
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

   integer (int_kind) :: num_writes  ! place-holder until more than one
                                     ! time level written to a single file

   logical (log_kind) :: supported

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

      call define_nstd_netcdf(data_file, ndims, io_dims, field_id,  &
                              short_name, long_name, units,         &
                              coordinates=coordinates,              &
                              missing_value=missing_value,          &
                              fill_value=fill_value,                &
                              nftype=nftype                         )

!-----------------------------------------------------------------------
!
!  write an io field
!
!-----------------------------------------------------------------------

   case ('write')

      num_writes =  1  ! for now, only support one time value per output file
      
      call write_nstd_netcdf(                     &
           data_file, field_id,num_writes,        &
           ndims, io_dims,nftype,                 &
           implied_time_dim=implied_time_dim,     &
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
