module cam_history_support
  use shr_kind_mod, only : r8=>shr_kind_r8, i8=>shr_kind_i8, shr_kind_cl
  use pio,          only : var_desc_t, file_desc_t
  use abortutils,   only : endrun
  use cam_logfile,  only: iulog
  use spmd_utils, only : masterproc

  implicit none
  private

  integer, parameter, public :: max_string_len = 256   ! Length of strings
  integer, parameter, public :: max_chars = shr_kind_cl         ! max chars for char variables
  integer, parameter, public :: fieldname_len = 24   ! max chars for field name
  integer, parameter, public :: fieldname_suffix_len =  3 ! length of field name suffix ("&IC")
  integer, parameter, public :: fieldname_lenp2      = fieldname_len + 2 ! allow for extra characters
  integer, parameter, public :: max_fieldname_len    = fieldname_len + fieldname_suffix_len ! max chars for field name (including suffix)

  integer, parameter, public :: pflds = 750               ! max number of fields for namelist entries fincl and fexcl
                                                   ! also used in write restart
  integer, parameter, public :: ptapes    = 8             ! max number of tapes


  type, public :: column_info
     character(len=max_chars) :: lat_name ! latitude name for this column or columns
     character(len=max_chars) :: lon_name ! latitude name for this column or columns
     integer :: num_lats            ! number of lats in a group of contiguous columns
     integer :: num_lons            ! number of lons in a group of contiguous columns
     integer :: columnlat(2)       ! beginning and ending latitude (range) dimensioned by groups
     integer :: columnlon(2)       ! beginning and ending longitude (range) dimensioned by groups
     integer, pointer :: hmask(:,:)
     integer, pointer :: hmask_dyn(:,:)
     integer :: column_cnt
     integer :: htape             ! need to keep track of which tape we are writing for iodesc reuse
  end type column_info

  type, public :: field_info

     ! (i.e., how often "outfld" is called):  every other; only during LW/SW
     ! radiation calcs; etc.
     logical :: flag_xyfill                    ! non-applicable xy points flagged with fillvalue
     integer, pointer :: mdims(:)              ! indicies into hist_mdims list

     real(r8) :: fillvalue                     ! fillvalue for this variable, set to default if not explicit in addfld

     integer :: numlev                         ! vertical dimension (.nc file and internal arr)

     integer :: begdim1                        ! on-node dim1 start index
     integer :: enddim1                        ! on-node dim1 end index

     integer :: begdim2                        ! on-node dim2 start index
     integer :: enddim2                        ! on-node dim2 end index

     integer :: begdim3                        ! on-node chunk or lat start index
     integer :: enddim3                        ! on-node chunk or lat end index

     integer :: decomp_type                    ! type of decomposition (physics or dynamics)

     integer :: vector_compliment              ! id for vector compliment for interpolation
     integer, pointer :: colperdim3(:)         ! number of valid elements per chunk or lat

     character(len=max_fieldname_len) :: name     ! field name
     character(len=max_chars) :: long_name        ! long name
     character(len=max_chars) :: units            ! units
     character(len=max_chars) :: sampling_seq     ! sampling sequence - if not every timestep, how often field is sampled
  end type field_info

  real(r8), parameter, public :: fillvalue = 1.e36_r8     ! fill value for netcdf fields


!
! hentry: elements of an entry in the list of active fields on a single history file
!
   type, public:: hentry
      type (field_info)     :: field            ! field information
      character(len=1)      :: avgflag          ! averaging flag
      character(len=max_chars) :: time_op          ! time operator (e.g. max, min, avg)
      character(len=max_chars),pointer :: field_column_name(:) ! names of column groups written to tape

      integer :: hwrt_prec                      ! history output precision
      real(r8), pointer :: hbuf(:,:,:)
      type(var_desc_t), pointer :: varid(:)      ! variable ids
      integer, pointer :: nacs(:,:)             ! accumulation counter
      type(var_desc_t) :: nacs_varid 
   end type hentry



! active_entry: vehicle for producing a ragged array
!
   type, public:: active_entry
	

      type(hentry), pointer :: hlist(:)

      type (column_info),pointer :: column   (:) ! array of history tape column entries
      type (column_info),pointer :: column_st(:) ! array of history tape column entries for staggered grid (FV)


!
! PIO ids
!

      type(file_desc_t) :: File            ! PIO file id

      type(var_desc_t) :: mdtid            ! var id for timestep
      type(var_desc_t) :: ndbaseid         ! var id for base day
      type(var_desc_t) :: nsbaseid         ! var id for base seconds of base day
      type(var_desc_t) :: nbdateid         ! var id for base date
      type(var_desc_t) :: nbsecid          ! var id for base seconds of base date
      type(var_desc_t) :: ndcurid          ! var id for current day
      type(var_desc_t) :: nscurid          ! var id for current seconds of current day
      type(var_desc_t) :: dateid           ! var id for current date
      type(var_desc_t) :: co2vmrid         ! var id for co2 volume mixing ratio
      type(var_desc_t) :: ch4vmrid         ! var id for ch4 volume mixing ratio
      type(var_desc_t) :: n2ovmrid         ! var id for n2o volume mixing ratio
      type(var_desc_t) :: f11vmrid         ! var id for f11 volume mixing ratio
      type(var_desc_t) :: f12vmrid         ! var id for f12 volume mixing ratio
      type(var_desc_t) :: sol_tsiid        ! var id for total solar irradiance (W/m2)
      type(var_desc_t) :: datesecid        ! var id for curent seconds of current date
#if ( defined BFB_CAM_SCAM_IOP )
      type(var_desc_t) :: bdateid         ! var id for base date
      type(var_desc_t) :: tsecid        ! var id for curent seconds of current date
#endif
      type(var_desc_t) :: nstephid         ! var id for current timestep
      type(var_desc_t) :: timeid           ! var id for time
      type(var_desc_t) :: tbndid           ! var id for time_bnds
      type(var_desc_t) :: date_writtenid   ! var id for date time sample written
      type(var_desc_t) :: time_writtenid   ! var id for time time sample written
      type(var_desc_t) :: nlonid           ! var id for number of longitudes
      type(var_desc_t) :: wnummaxid        ! var id for cutoff fourier wavenumber (reduced grid)
      type(var_desc_t) :: f107id           ! var id for f107
      type(var_desc_t) :: f107aid          ! var id for f107a
      type(var_desc_t) :: kpid             ! var id for kp
      type(var_desc_t) :: apid             ! var id for ap

   end type active_entry

   type, public :: hist_mdim_t
      character(len=16) :: name
      integer :: value
   end type hist_mdim_t

   integer, public :: registeredmdims=0, maxvarmdims=1
   integer, parameter :: maxmdims = 25    ! arbitrary limit
   type(hist_mdim_t), public :: hist_mdims(maxmdims)

   public :: register_hist_mdim, lookup_hist_mdim_indicies, sec2hms, date2yyyymmdd

   interface assignment(=)
      module procedure field_copy
   end interface



 contains


   subroutine field_copy(f_out, f_in)
     type(field_info), intent(in) :: f_in
     type(field_info), intent(out) :: f_out

     f_out%flag_xyfill= f_in%flag_xyfill                    
     f_out%fillvalue= f_in%fillvalue                    
     f_out%numlev =  f_in%numlev                         ! vertical dimension (.nc file and internal arr)
     f_out%begdim1 = f_in%begdim1                        ! on-node dim1 start index
     f_out%enddim1 = f_in%enddim1                        ! on-node dim1 end index
     f_out%begdim2 = f_in%begdim2                        ! on-node dim2 start index
     f_out%enddim2 = f_in%enddim2                        ! on-node dim2 end index
     f_out%begdim3 = f_in%begdim3                        ! on-node chunk or lat start index
     f_out%enddim3 = f_in%enddim3                        ! on-node chunk or lat end index
     f_out%decomp_type = f_in%decomp_type                    ! type of decomposition (physics or dynamics)

     f_out%vector_compliment = f_in%vector_compliment              ! id for vector compliment for interpolation

     f_out%name = f_in%name     ! field name
     f_out%long_name = f_in%long_name        ! long name
     f_out%units = f_in%units            ! units
     f_out%sampling_seq =  f_in%sampling_seq     ! sampling sequence - if not every timestep, how often field is sampled     

     if(associated(f_in%mdims)) then
        f_out%mdims=>f_in%mdims
     else
        nullify(f_out%mdims)
     end if
     if(associated(f_in%colperdim3)) then
        f_out%colperdim3=>f_in%colperdim3
     else
        nullify(f_out%colperdim3)
     end if




   end subroutine field_copy


   subroutine register_hist_mdim(name, val)
     character(len=*), intent(in) :: name
     integer, intent(in) :: val
     character(len=120) :: errormsg
     integer :: i
     logical :: found

     found=.false.
    !  if(masterproc) write(iulog,*) 'Registering mdim name ',name,' with value: ',val, registeredmdims
     do i=1,registeredmdims
        if(name .eq. hist_mdims(i)%name) then
           if(val .ne. hist_mdims(i)%value) then
              write(errormsg,*) 'Name ',name,' already registered with different value ',val,hist_mdims(i)%value
              call endrun(errormsg)
           end if
           found = .true.
           exit
        end if
     end do

     if(.not. found) then
        registeredmdims=registeredmdims+1
        hist_mdims(registeredmdims)%name=name
        hist_mdims(registeredmdims)%value=val
     end if

   end subroutine register_hist_mdim

   subroutine lookup_hist_mdim_indicies(mdimnames, mdimindicies)
     character(len=*), intent(in) :: mdimnames(:)
     integer, intent(out) :: mdimindicies(:)

     integer :: i, j
     integer :: cnt
     character(len=120) :: errormsg
     character(len=16) :: name


     cnt = size(mdimnames)
     mdimindicies = -1


     do j=1,cnt
        name = mdimnames(j)
        do i=1,registeredmdims
           if(name .eq. hist_mdims(i)%name) then
              mdimindicies(j)=i
           end if
        end do
     end do
     do j=1,cnt
        if(mdimindicies(j)<0) then
           do i=1,registeredmdims		
              print *,__FILE__,__LINE__,i,hist_mdims(i)%name
           end do
           write(errormsg,*) 'Name ',mdimnames(j),' not in registered mdimnames'
           call endrun(errormsg)
        end if
     end do


   end subroutine lookup_hist_mdim_indicies

!#######################################################################

character(len=8) function sec2hms (seconds)

! Input arguments

   integer, intent(in) :: seconds

! Local workspace

   integer :: hours     ! hours of hh:mm:ss
   integer :: minutes   ! minutes of hh:mm:ss
   integer :: secs      ! seconds of hh:mm:ss

   if (seconds < 0 .or. seconds > 86400) then
      write(iulog,*)'SEC2HRS: bad input seconds:', seconds
      call endrun ()
   end if

   hours   = seconds / 3600
   minutes = (seconds - hours*3600) / 60
   secs    = (seconds - hours*3600 - minutes*60)

   if (minutes < 0 .or. minutes > 60) then
      write(iulog,*)'SEC2HRS: bad minutes = ',minutes
      call endrun ()
   end if

   if (secs < 0 .or. secs > 60) then
      write(iulog,*)'SEC2HRS: bad secs = ',secs
      call endrun ()
   end if

   write(sec2hms,80) hours, minutes, secs
80 format(i2.2,':',i2.2,':',i2.2)
   return
end function sec2hms
character(len=10) function date2yyyymmdd (date)

! Input arguments

   integer, intent(in) :: date

! Local workspace

   integer :: year    ! year of yyyy-mm-dd
   integer :: month   ! month of yyyy-mm-dd
   integer :: day     ! day of yyyy-mm-dd

   if (date < 0) then
      call endrun ('DATE2YYYYMMDD: negative date not allowed')
   end if

   year  = date / 10000
   month = (date - year*10000) / 100
   day   = date - year*10000 - month*100

   write(date2yyyymmdd,80) year, month, day
80 format(i4.4,'-',i2.2,'-',i2.2)
   return
end function date2yyyymmdd

!#######################################################################


end module cam_history_support
