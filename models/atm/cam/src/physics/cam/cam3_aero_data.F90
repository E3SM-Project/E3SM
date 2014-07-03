module cam3_aero_data
!----------------------------------------------------------------------- 
! 
! Purposes: 
!       read, store, interpolate, and return fields
!         of aerosols to CAM.  The initialization
!         file (mass.nc) is assumed to be a monthly climatology
!         of aerosols from MATCH (on a sigma pressure
!         coordinate system).
!       also provide a "background" aerosol field to correct
!         for any deficiencies in the physical parameterizations
!         This fields is a "tuning" parameter.
!       Public methods:
!       (1) - initialization
!          read aerosol masses from external file
!             also pressure coordinates
!          convert from monthly average values to mid-month values
!       (2) - interpolation (time and vertical)
!          interpolate onto pressure levels of CAM
!          interpolate to time step of CAM
!          return mass of aerosols 
!
!-----------------------------------------------------------------------

  use shr_kind_mod,   only: r8 => shr_kind_r8
  use shr_scam_mod,   only: shr_scam_GetCloseLatLon
  use spmd_utils,     only: masterproc
  use ppgrid,         only: pcols, pver, pverp, begchunk, endchunk
  use phys_grid,      only: get_ncols_p, scatter_field_to_chunk
  use time_manager,   only: get_curr_calday
  use infnan,         only: nan, assignment(=)
  use abortutils,     only: endrun
  use scamMod,        only: scmlon,scmlat,single_column
  use error_messages, only: handle_ncerr
  use physics_types,  only: physics_state
  use boundarydata,   only: boundarydata_init, boundarydata_type
  use perf_mod,       only: t_startf, t_stopf
  use cam_logfile,    only: iulog
  use netcdf

  implicit none
  private
  save

  public :: &
     cam3_aero_data_readnl,       & ! read namelist
     cam3_aero_data_register,     & ! register these aerosols with pbuf2d
     cam3_aero_data_init,         & ! read from file, interpolate onto horiz grid
     cam3_aero_data_timestep_init   ! update data-aerosols to this timestep

  ! namelist variables
  logical, public :: cam3_aero_data_on = .false.
  character(len=256) :: bndtvaer = 'bndtvaer'   ! full pathname for time-variant aerosol mass climatology dataset

  ! naer is number of species in climatology
  integer, parameter :: naer = 11

  real(r8), parameter :: wgt_sscm = 6.0_r8 / 7.0_r8 ! Fraction of total seasalt mass in coarse mode

  ! indices to aerosol array (species portion)
  integer, parameter :: &
      idxSUL   =  1, &
      idxSSLTA =  2, & ! accumulation mode
      idxSSLTC =  3, & ! coarse mode
      idxOCPHO =  8, &
      idxBCPHO =  9, &
      idxOCPHI =  10, &
      idxBCPHI = 11

  ! indices to sections of array that represent 
  ! groups of aerosols
  integer, parameter :: &
      idxSSLTfirst    = 2, numSSLT  = 2, &
      idxDUSTfirst    = 4, &
      numDUST         = 4, &
      idxCARBONfirst = 8, &
      numCARBON      = 4

  ! names of aerosols are they are represented in
  ! the climatology file.
  ! Appended '_V' indicates field has been vertically summed.
  character(len=8), parameter :: aerosol_name(naer) =  &
     (/"MSUL_V  "&
      ,"MSSLTA_V"&
      ,"MSSLTC_V"&
      ,"MDUST1_V"&
      ,"MDUST2_V"&
      ,"MDUST3_V"&
      ,"MDUST4_V"&
      ,"MOCPHO_V"&
      ,"MBCPHO_V"&
      ,"MOCPHI_V"&
      ,"MBCPHI_V"/)

  ! number of different "groups" of aerosols
  integer, parameter :: num_aer_groups=4

  ! which group does each bin belong to?
  integer, dimension(naer), parameter ::  &
      group =(/1,2,2,3,3,3,3,4,4,4,4/)

  ! name of each group
  character(len=10), dimension(num_aer_groups), parameter :: &
      aerosol_names = (/'sul  ','sslt ','dust ','car  '/)

  ! this boundarydata_type is used for datasets in the ncols format only.
  type(boundarydata_type) :: aerosol_datan

  integer :: aernid = -1           ! netcdf id for aerosol file (init to invalid)
  integer :: species_id(naer) = -1 ! netcdf_id of each aerosol species (init to invalid)
  integer :: Mpsid                 ! netcdf id for MATCH PS
  integer :: nm = 1                ! index to prv month in array. init to 1 and toggle between 1 and 2
  integer :: np = 2                ! index to nxt month in array. init to 2 and toggle between 1 and 2
  integer :: mo_nxt = huge(1)      ! index to nxt month in file

  real(r8) :: cdaym                ! calendar day of prv month
  real(r8) :: cdayp                ! calendar day of next month

  ! aerosol mass 
  real(r8), allocatable :: aer_mass(:, :, :, :)

  ! Days into year for mid month date
  ! This variable is dumb, the dates are in the dataset to be read in but they are
  ! slightly different than this so getting rid of it causes a change which 
  ! exceeds roundoff.
  real(r8) :: Mid(12) = (/16.5_r8,  46.0_r8,  75.5_r8, 106.0_r8, 136.5_r8, 167.0_r8, &
                         197.5_r8, 228.5_r8, 259.0_r8, 289.5_r8, 320.0_r8, 350.5_r8 /)
  
  !  values read from file and temporary values used for interpolation
  !
  !  aerosolc is:
  !  Cumulative Mass at midpoint of each month
  !    on CAM's horizontal grid (col)
  !    on MATCH's levels (lev)
  !  aerosolc
  integer, parameter :: paerlev = 28       ! number of levels for aerosol fields (MUST = naerlev)
  integer :: naerlev                       ! size of level dimension in MATCH data
  integer :: naerlon
  integer :: naerlat
  real(r8), pointer :: M_hybi(:)           ! MATCH hybi
  real(r8), pointer :: M_ps(:,:)           ! surface pressure from MATCH file
  real(r8), pointer :: aerosolc(:,:,:,:,:) ! Aerosol cumulative mass from MATCH
  real(r8), pointer :: M_ps_cam_col(:,:,:) ! PS from MATCH on Cam Columns

  ! indices for fields in the physics buffer
  integer :: cam3_sul_idx, cam3_ssam_idx, cam3_sscm_idx, &
      cam3_dust1_idx, cam3_dust2_idx, cam3_dust3_idx, cam3_dust4_idx,&
      cam3_ocpho_idx, cam3_bcpho_idx, cam3_ocphi_idx, cam3_bcphi_idx

!================================================================================================
contains
!================================================================================================

subroutine cam3_aero_data_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'cam3_aero_data_readnl'

   namelist /cam3_aero_data_nl/ cam3_aero_data_on, bndtvaer
   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'cam3_aero_data_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, cam3_aero_data_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(cam3_aero_data_on, 1, mpilog, 0, mpicom)
   call mpibcast(bndtvaer, len(bndtvaer), mpichar, 0, mpicom)
#endif

   ! Prevent using these before they are set.
   cdaym = nan
   cdayp = nan

end subroutine cam3_aero_data_readnl

!================================================================================================

subroutine cam3_aero_data_register

   ! register old prescribed aerosols with physics buffer
   
   use physics_buffer, only: pbuf_add_field, dtype_r8

   call pbuf_add_field('cam3_sul',  'physpkg',dtype_r8,(/pcols,pver/),cam3_sul_idx)
   call pbuf_add_field('cam3_ssam', 'physpkg',dtype_r8,(/pcols,pver/),cam3_ssam_idx)
   call pbuf_add_field('cam3_sscm', 'physpkg',dtype_r8,(/pcols,pver/),cam3_sscm_idx)
   call pbuf_add_field('cam3_dust1','physpkg',dtype_r8,(/pcols,pver/),cam3_dust1_idx)
   call pbuf_add_field('cam3_dust2','physpkg',dtype_r8,(/pcols,pver/),cam3_dust2_idx)
   call pbuf_add_field('cam3_dust3','physpkg',dtype_r8,(/pcols,pver/),cam3_dust3_idx)
   call pbuf_add_field('cam3_dust4','physpkg',dtype_r8,(/pcols,pver/),cam3_dust4_idx)
   call pbuf_add_field('cam3_ocpho','physpkg',dtype_r8,(/pcols,pver/),cam3_ocpho_idx)
   call pbuf_add_field('cam3_bcpho','physpkg',dtype_r8,(/pcols,pver/),cam3_bcpho_idx)
   call pbuf_add_field('cam3_ocphi','physpkg',dtype_r8,(/pcols,pver/),cam3_ocphi_idx)
   call pbuf_add_field('cam3_bcphi','physpkg',dtype_r8,(/pcols,pver/),cam3_bcphi_idx)

end subroutine cam3_aero_data_register

!================================================================================================

subroutine cam3_aero_data_init(phys_state)
!------------------------------------------------------------------
!  Reads in:
!     file from which to read aerosol Masses on CAM grid. Currently
!        assumed to be MATCH ncep runs, averaged by month.
!     NOTE (Data have been externally interpolated onto CAM grid 
!        and backsolved to provide Mid-month values)
!     
!  Populates:
!     module variables:
!       aerosolc(pcols,paerlev+1,begchunk:endchunk,naer,2))
!       aerosolc(  column_index
!                , level_index (match levels)
!                , chunk_index 
!                , species_index
!                , month = 1:2 )
!       M_hybi(level_index = Lev_MATCH) = pressure at mid-level.
!       M_ps_cam_col(column,chunk,month) ! PS from MATCH on Cam Columns
!
!  Method:
!    read data from file
!    allocate memory for storage of aerosol data on CAM horizontal grid
!    distribute data to remote nodes
!    populates the module variables
!
!------------------------------------------------------------------
   use ioFileMod,    only: getfil

#if ( defined SPMD )
   use mpishorthand
#endif
   type(physics_state), intent(in) :: phys_state(begchunk:endchunk)

! local variables

   integer :: naerlev

   integer dateid                       ! netcdf id for date variable
   integer secid                        ! netcdf id for seconds variable
   integer londimid                     ! netcdf id for longitude dimension
   integer latdimid                     ! netcdf id for latitude dimension
   integer levdimid                     ! netcdf id for level dimension

   integer timesiz                      ! number of time samples (=12) in netcdf file
   integer latid                        ! netcdf id for latitude variable
   integer Mhybiid                      ! netcdf id for MATCH hybi
   integer timeid                       ! netcdf id for time variable
   integer dimids(nf90_max_var_dims)      ! variable shape
   integer :: start(4)                  ! start vector for netcdf calls
   integer :: kount(4)                  ! count vector for netcdf calls
   integer mo                           ! month index
   integer m                            ! constituent index
   integer :: n                         ! loop index
   integer :: i,j,k                     ! spatial indices
   integer :: date_aer(12)              ! Date on aerosol dataset (YYYYMMDD)
   integer :: attnum                    ! attribute number
   integer :: ierr                      ! netcdf return code
   real(r8) ::  coldata(paerlev)    ! aerosol field read in from dataset
   integer :: ret
   integer mo_prv                       ! index to previous month
   integer latidx,lonidx

   character(len=8) :: aname                   ! temporary aerosol name
   character(len=8) :: tmp_aero_name(naer) ! name for input to boundary data

   character(len=256) :: locfn          ! netcdf local filename to open
!
! aerosol_data will be read in from the aerosol boundary dataset, then scattered to chunks
! after filling in the bottom level with zeros
! 
   real(r8), allocatable :: aerosol_data(:,:,:)    ! aerosol field read in from dataset
   real(r8), allocatable :: aerosol_field(:,:,:)   ! (plon,paerlev+1,plat)  aerosol field to be scattered
   real(r8) :: caldayloc                           ! calendar day of current timestep
   real(r8) :: closelat,closelon

   character(len=*), parameter :: subname = 'cam3_aero_data_init'
   !------------------------------------------------------------------

   call t_startf(subname)

   allocate (aer_mass(pcols, pver, naer, begchunk:endchunk) )

   ! set new aerosol names because input file has 1 seasalt bin
   do m = 1, naer
      tmp_aero_name(m)=aerosol_name(m)
      if (aerosol_name(m)=='MSSLTA_V') tmp_aero_name(m) = 'MSSLT_V'
      if (aerosol_name(m)=='MSSLTC_V') tmp_aero_name(m) = 'MSSLT_V'
   end do

   allocate (aerosolc(pcols,paerlev+1,begchunk:endchunk,naer,2))
   aerosolc(:,:,:,:,:) = 0._r8

   caldayloc = get_curr_calday ()
   
   if (caldayloc < Mid(1)) then
      mo_prv = 12
      mo_nxt =  1
   else if (caldayloc >= Mid(12)) then
      mo_prv = 12
      mo_nxt =  1
   else
      do i = 2 , 12
         if (caldayloc < Mid(i)) then
            mo_prv = i-1
            mo_nxt = i
            exit
         end if
      end do
   end if

   ! Set initial calendar day values
   cdaym = Mid(mo_prv)
   cdayp = Mid(mo_nxt)

   if (masterproc) &
      write(iulog,*) subname//': CAM3 prescribed aerosol dataset is: ', trim(bndtvaer)

   call getfil (bndtvaer, locfn, 0)

   call handle_ncerr( nf90_open (locfn, 0, aernid),&
      subname, __LINE__)

   if (single_column) &
      call shr_scam_GetCloseLatLon(aernid,scmlat,scmlon,closelat,closelon,latidx,lonidx)

   ! Check to see if this dataset is in ncol format. 
   aerosol_datan%isncol=.false.
   ierr = nf90_inq_dimid( aernid,  'ncol', londimid )
   if ( ierr==NF90_NOERR ) then

      aerosol_datan%isncol=.true.
      call handle_ncerr(nf90_close(aernid),subname, __LINE__)

      call boundarydata_init(bndtvaer, phys_state, tmp_aero_name, naer, &
                             aerosol_datan, 3)

      aerosolc(:,1:paerlev,:,:,:)=aerosol_datan%fields

      M_ps_cam_col=>aerosol_datan%ps
      M_hybi=>aerosol_datan%hybi

   else 

      ! Allocate memory for dynamic arrays local to this module
      allocate (M_ps_cam_col(pcols,begchunk:endchunk,2))
      allocate (M_hybi(paerlev+1))
      ! TBH:  HACK to avoid use of uninitialized values when ncols < pcols
      M_ps_cam_col(:,:,:) = 0._r8

      if (masterproc) then

         ! First ensure dataset is CAM-ready

         call handle_ncerr(nf90_inquire_attribute (aernid, nf90_global, 'cam-ready', attnum=attnum),&
              subname//': interpaerosols needs to be run to create a cam-ready aerosol dataset')

         ! Get and check dimension info

         call handle_ncerr( nf90_inq_dimid( aernid,  'lon', londimid ),&
              subname, __LINE__)
         call handle_ncerr( nf90_inq_dimid( aernid,  'lev', levdimid ),&
              subname, __LINE__)
         call handle_ncerr( nf90_inq_dimid( aernid, 'time',   timeid ),&
              subname, __LINE__)
         call handle_ncerr( nf90_inq_dimid( aernid,  'lat', latdimid ),&
              subname, __LINE__)
         call handle_ncerr( nf90_inquire_dimension( aernid, londimid, len=naerlon ),&
              subname, __LINE__)
         call handle_ncerr( nf90_inquire_dimension( aernid, levdimid, len=naerlev ),&
              subname, __LINE__)
         call handle_ncerr( nf90_inquire_dimension( aernid, latdimid, len=naerlat ),&
              subname, __LINE__)
         call handle_ncerr( nf90_inquire_dimension( aernid,   timeid, len=timesiz ),&
              subname, __LINE__)

         call handle_ncerr( nf90_inq_varid( aernid, 'date',   dateid ),&
              subname, __LINE__)
         call handle_ncerr( nf90_inq_varid( aernid, 'datesec', secid ),&
              subname, __LINE__)

         do m = 1, naer
            aname=aerosol_name(m)
            ! rename because file has only one seasalt field
            if (aname=='MSSLTA_V') aname = 'MSSLT_V'
            if (aname=='MSSLTC_V') aname = 'MSSLT_V'
            call handle_ncerr( nf90_inq_varid( aernid, TRIM(aname), species_id(m)), &
               subname, __LINE__)
         end do

         call handle_ncerr( nf90_inq_varid( aernid, 'lat', latid   ),&
              subname, __LINE__)

         ! quick sanity check on one field
         call handle_ncerr( nf90_inquire_variable (aernid, species_id(1), dimids=dimids),&
              subname, __LINE__)

         if ( (dimids(4) /= timeid) .or. &
              (dimids(3) /= levdimid) .or. &
              (dimids(2) /= latdimid) .or. &
              (dimids(1) /= londimid) ) then
            write(iulog,*) subname//': Data must be ordered time, lev, lat, lon'
            write(iulog,*) 'data are       ordered as', dimids(4), dimids(3), dimids(2), dimids(1)
            write(iulog,*) 'data should be ordered as', timeid, levdimid, latdimid, londimid
            call endrun ()
         end if

         ! use hybi,PS from MATCH
         call handle_ncerr( nf90_inq_varid( aernid, 'hybi', Mhybiid   ),&
              subname, __LINE__)
         call handle_ncerr( nf90_inq_varid( aernid, 'PS', Mpsid   ),&
              subname, __LINE__)

         ! check dimension order for MATCH's surface pressure
         call handle_ncerr( nf90_inquire_variable (aernid, Mpsid, dimids=dimids),&
              subname, __LINE__)
         if ( (dimids(3) /= timeid) .or. &
              (dimids(2) /= latdimid) .or. &
              (dimids(1) /= londimid) ) then
            write(iulog,*) subname//': Pressure must be ordered time, lat, lon'
            write(iulog,*) 'data are       ordered as', dimids(3), dimids(2), dimids(1)
            write(iulog,*) 'data should be ordered as', timeid, levdimid, latdimid, londimid
            call endrun ()
         end if

         ! read in hybi from MATCH
         call handle_ncerr( nf90_get_var (aernid, Mhybiid, M_hybi),&
              subname, __LINE__)

         ! Retrieve date and sec variables.
         call handle_ncerr( nf90_get_var (aernid, dateid, date_aer),&
              subname, __LINE__)
         if (timesiz < 12) then
            write(iulog,*) subname//': When cycling aerosols, dataset must have 12 consecutive ', &
                 'months of data starting with Jan'
            write(iulog,*) 'Current dataset has only ',timesiz,' months'
            call endrun ()
         end if
         do mo = 1,12
            if (mod(date_aer(mo),10000)/100 /= mo) then
               write(iulog,*) subname//': When cycling aerosols, dataset must have 12 consecutive ', &
                    'months of data starting with Jan'
               write(iulog,*)'Month ',mo,' of dataset says date=',date_aer(mo)
               call endrun ()
            end if
         end do
         if (single_column) then
            naerlat=1
            naerlon=1
         endif
         kount(:) = (/naerlon,naerlat,paerlev,1/)
      end if          ! masterproc

      ! broadcast hybi to nodes

#if ( defined SPMD )
      call mpibcast (M_hybi, paerlev+1, mpir8, 0, mpicom)
      call mpibcast (kount, 3, mpiint, 0, mpicom)
      naerlon = kount(1)
      naerlat = kount(2)
#endif
      allocate(aerosol_field(kount(1),kount(3)+1,kount(2)))
      allocate(M_ps(kount(1),kount(2)))
      if (masterproc) allocate(aerosol_data(kount(1),kount(2),kount(3)))

      ! Retrieve Aerosol Masses (kg/m^2 in each layer), transpose to model order (lon,lev,lat),
      ! then scatter to slaves.
      if (nm /= 1 .or. np /= 2) call endrun (subname//': bad nm or np value')
      do n=nm,np
         if (n == 1) then
            mo = mo_prv
         else
            mo = mo_nxt
         end if
         
         do m=1,naer
            if (masterproc) then
               if (single_column) then
                  start(:) = (/lonidx,latidx,1,mo/)
               else
                  start(:) = (/1,1,1,mo/)
               endif
               kount(:) = (/naerlon,naerlat,paerlev,1/)

               call handle_ncerr( nf90_get_var (aernid, species_id(m),aerosol_data, start, kount),&
                    subname, __LINE__)
               do j=1,naerlat
                  do k=1,paerlev
                     aerosol_field(:,k,j) = aerosol_data(:,j,k)
                  end do
                  aerosol_field(:,paerlev+1,j) = 0._r8   ! value at bottom
               end do
               
            end if
            call scatter_field_to_chunk (1, paerlev+1, 1, naerlon, aerosol_field, &
                 aerosolc(:,:,:,m,n))
         end do

         ! Retrieve PS from Match

         if (masterproc) then
            if (single_column) then
               start(:) = (/lonidx,latidx,mo,-1/)
            else
               start(:) = (/1,1,mo,-1/)
            endif
            kount(:) = (/naerlon,naerlat,1,-1/)
            call handle_ncerr( nf90_get_var(aernid, Mpsid, M_ps,start,kount),&
                 subname, __LINE__)
         end if
         call scatter_field_to_chunk (1, 1, 1, naerlon, M_ps(:,:), M_ps_cam_col(:,:,n))
      end do     ! n=nm,np (=1,2)

      if(masterproc) deallocate(aerosol_data)
      deallocate(aerosol_field)

   end if   ! Check to see if this dataset is in ncol format. 

   call t_stopf(subname)

end subroutine cam3_aero_data_init

!================================================================================================

subroutine cam3_aero_data_timestep_init(pbuf2d,  phys_state)
!------------------------------------------------------------------
!
!  Input:
!     time at which aerosol masses are needed (get_curr_calday())
!     chunk index
!     CAM's vertical grid (pint)
!
!  Output:
!     values for Aerosol Mass at time specified by get_curr_calday
!     on vertical grid specified by pint (aer_mass) :: aerosol at time t
!
!  Method:
!     first determine which indexs of aerosols are the bounding data sets
!     interpolate both onto vertical grid aerm(),aerp().
!     from those two, interpolate in time.
!
!------------------------------------------------------------------

   use interpolate_data, only: get_timeinterp_factors
   
   use physics_buffer, only: physics_buffer_desc, dtype_r8, pbuf_set_field, pbuf_get_chunk
   use cam_logfile,     only: iulog
   use ppgrid,          only: begchunk,endchunk
   use physconst,       only: gravit

!
! aerosol fields interpolated to current time step
!   on pressure levels of this time step.
! these should be made read-only for other modules
! Is allocation done correctly here?
!
   
   type(physics_buffer_desc), pointer :: pbuf2d(:,:)
   type(physics_state), intent(in), dimension(begchunk:endchunk) :: phys_state

!
! Local workspace
!
   type(physics_buffer_desc), pointer :: phys_buffer_chunk(:)
   real(r8) :: pint(pcols,pverp)  ! interface pres.
   integer :: c                           ! chunk index
   real(r8) caldayloc                     ! calendar day of current timestep
   real(r8) fact1, fact2                  ! time interpolation factors

   integer i, k, j                        ! spatial indices
   integer m                              ! constituent index
   integer lats(pcols),lons(pcols)        ! latitude and longitudes of column
   integer ncol                           ! number of columns
   integer lchnk                          ! chunk index
   
   real(r8) speciesmin(naer)              ! minimal value for each species
!
! values before current time step "the minus month"
! aerosolm(pcols,pver) is value of preceeding month's aerosol masses
! aerosolp(pcols,pver) is value of next month's aerosol masses
!  (think minus and plus or values to left and right of point to be interpolated)
!
   real(r8) aerosolm(pcols,pver,naer,begchunk:endchunk) ! aerosol mass from MATCH in column,level at previous (minus) month
!
! values beyond (or at) current time step "the plus month"
!
   real(r8) aerosolp(pcols,pver,naer,begchunk:endchunk) ! aerosol mass from MATCH in column,level at next (plus) month 
   real(r8) :: mass_to_mmr(pcols,pver)

   character(len=*), parameter :: subname = 'cam3_aero_data_timestep_init'

   logical error_found
   !------------------------------------------------------------------

   call aerint(phys_state)

   caldayloc = get_curr_calday ()

   ! Determine time interpolation factors.  1st arg says we are cycling 1 year of data
   call get_timeinterp_factors (.true., mo_nxt, cdaym, cdayp, caldayloc, &
                    fact1, fact2, 'GET_AEROSOL:')

   ! interpolate (prv and nxt month) bounding datasets onto cam vertical grid.
   ! compute mass mixing ratios on CAMS's pressure coordinate
   !  for both the "minus" and "plus" months
   !
   !  This loop over chunk could probably be removed by working with the whole
   !  begchunk:endchunk group at once.  It would require a slight generalization 
   !  in vert_interpolate.
   do c = begchunk,endchunk  
                                
      lchnk = phys_state(c)%lchnk
      pint = phys_state(c)%pint
      ncol = get_ncols_p(c)

      call vert_interpolate (M_ps_cam_col(:,c,nm), pint, nm, aerosolm(:,:,:,c), ncol, c)
      call vert_interpolate (M_ps_cam_col(:,c,np), pint, np, aerosolp(:,:,:,c), ncol, c)

      ! Time interpolate.
      do m=1,naer
         do k=1,pver
            do i=1,ncol
               aer_mass(i,k,m,c) = aerosolm(i,k,m,c)*fact1 + aerosolp(i,k,m,c)*fact2
            end do
         end do
         ! Partition seasalt aerosol mass
         if (m .eq. idxSSLTA) then
            aer_mass(:ncol,:,m,c) = (1._r8-wgt_sscm)*aer_mass(:ncol,:,m,c) ! fraction of seasalt mass in accumulation mode
         elseif (m .eq. idxSSLTC) then
            aer_mass(:ncol,:,m,c) = wgt_sscm*aer_mass(:ncol,:,m,c)      ! fraction of seasalt mass in coarse mode
         endif
      end do

      ! exit if mass is negative (we have previously set
      !  cumulative mass to be a decreasing function.)
      speciesmin(:) = 0._r8 ! speciesmin(m) = 0 is minimum mass for each species
 
      error_found = .false.
      do m=1,naer
         do k=1,pver
            do i=1,ncol
               if (aer_mass(i, k, m,c) < speciesmin(m)) error_found = .true.
            end do
         end do
      end do
      if (error_found) then
         do m=1,naer
            do k=1,pver
               do i=1,ncol
                  if (aer_mass(i, k, m,c) < speciesmin(m)) then
                     write(iulog,*) subname//': negative mass mixing ratio, exiting'
                     write(iulog,*) 'm, column, pver',m, i, k ,aer_mass(i, k, m,c)
                     call endrun ()
                  end if
               end do
            end do
         end do
      end if
      do k = 1, pver
         mass_to_mmr(1:ncol,k) = gravit/(pint(1:ncol,k+1)-pint(1:ncol,k))
      enddo

      phys_buffer_chunk => pbuf_get_chunk(pbuf2d, lchnk)

      call pbuf_set_field(phys_buffer_chunk, cam3_sul_idx,   aer_mass(1:ncol,:,        idxSUL,c)*mass_to_mmr(:ncol,:), &
           start=(/1,1/), kount=(/ncol,pver/))
      call pbuf_set_field(phys_buffer_chunk, cam3_ssam_idx,  aer_mass(1:ncol,:,      idxSSLTA,c)*mass_to_mmr(:ncol,:), &
           start=(/1,1/), kount=(/ncol,pver/))
      call pbuf_set_field(phys_buffer_chunk, cam3_sscm_idx,  aer_mass(1:ncol,:,      idxSSLTC,c)*mass_to_mmr(:ncol,:), &
           start=(/1,1/), kount=(/ncol,pver/))
      call pbuf_set_field(phys_buffer_chunk, cam3_dust1_idx, aer_mass(1:ncol,:,  idxDUSTfirst,c)*mass_to_mmr(:ncol,:), &
           start=(/1,1/), kount=(/ncol,pver/))
      call pbuf_set_field(phys_buffer_chunk, cam3_dust2_idx, aer_mass(1:ncol,:,idxDUSTfirst+1,c)*mass_to_mmr(:ncol,:), &
           start=(/1,1/), kount=(/ncol,pver/))
      call pbuf_set_field(phys_buffer_chunk, cam3_dust3_idx, aer_mass(1:ncol,:,idxDUSTfirst+2,c)*mass_to_mmr(:ncol,:), &
           start=(/1,1/), kount=(/ncol,pver/))
      call pbuf_set_field(phys_buffer_chunk, cam3_dust4_idx, aer_mass(1:ncol,:,idxDUSTfirst+3,c)*mass_to_mmr(:ncol,:), &
           start=(/1,1/), kount=(/ncol,pver/))
      call pbuf_set_field(phys_buffer_chunk, cam3_ocpho_idx, aer_mass(1:ncol,:,      idxOCPHO,c)*mass_to_mmr(:ncol,:), &
           start=(/1,1/), kount=(/ncol,pver/))
      call pbuf_set_field(phys_buffer_chunk, cam3_bcpho_idx, aer_mass(1:ncol,:,      idxBCPHO,c)*mass_to_mmr(:ncol,:), &
           start=(/1,1/), kount=(/ncol,pver/))
      call pbuf_set_field(phys_buffer_chunk, cam3_ocphi_idx, aer_mass(1:ncol,:,      idxOCPHI,c)*mass_to_mmr(:ncol,:), &
           start=(/1,1/), kount=(/ncol,pver/))
      call pbuf_set_field(phys_buffer_chunk, cam3_bcphi_idx, aer_mass(1:ncol,:,      idxBCPHI,c)*mass_to_mmr(:ncol,:), &
           start=(/1,1/), kount=(/ncol,pver/))

   enddo ! c = begchunk:endchunk

end subroutine cam3_aero_data_timestep_init

!================================================================================================

subroutine vert_interpolate (Match_ps, pint, n, aerosol_mass, ncol, c)
!--------------------------------------------------------------------
! Input: match surface pressure, cam interface pressure, 
!        month index, number of columns, chunk index
! 
! Output: Aerosol mass mixing ratio (aerosol_mass)
!
! Method:
!         interpolate column mass (cumulative) from match onto
!           cam's vertical grid (pressure coordinate)
!         convert back to mass mixing ratio
!
!--------------------------------------------------------------------

   real(r8), intent(out) :: aerosol_mass(pcols,pver,naer)  ! aerosol mass from MATCH
   real(r8), intent(in) :: Match_ps(pcols)                ! surface pressure at a particular month
   real(r8), intent(in) :: pint(pcols,pverp)              ! interface pressure from CAM

   integer, intent(in) :: ncol,c                          ! chunk index and number of columns
   integer, intent(in) :: n                               ! prv or nxt month index
!
! Local workspace
!
   integer m                           ! index to aerosol species
   integer kupper(pcols)               ! last upper bound for interpolation
   integer i, k, kk, kkstart, kount    ! loop vars for interpolation
   integer isv, ksv, msv               ! loop indices to save

   logical bad                         ! indicates a bad point found
   logical lev_interp_comp             ! interpolation completed for a level 
   logical error_found

   real(r8) aerosol(pcols,pverp,naer)  ! cumulative mass of aerosol in column beneath upper 
                                       ! interface of level in column at particular month
   real(r8) dpl, dpu                   ! lower and upper intepolation factors
   real(r8) v_coord                    ! vertical coordinate
   real(r8) AER_diff                   ! temp var for difference between aerosol masses

   character(len=*), parameter :: subname = 'cam3_aero_data.vert_interpolate'
   !-----------------------------------------------------------------------

   call t_startf ('vert_interpolate')
!
! Initialize index array 
!
   do i=1,ncol
      kupper(i) = 1
   end do
!
! assign total mass to topmost level
!
   aerosol(:,1,:) = aerosolc(:,1,c,:,n)
!
! At every pressure level, interpolate onto that pressure level
!
   do k=2,pver
!
! Top level we need to start looking is the top level for the previous k
! for all longitude points
!
      kkstart = paerlev+1
      do i=1,ncol
         kkstart = min0(kkstart,kupper(i))
      end do
      kount = 0
!
! Store level indices for interpolation
!
! for the pressure interpolation should be comparing
! pint(column,lev) with M_hybi(lev)*M_ps_cam_col(month,column,chunk)
!
      lev_interp_comp = .false.
      do kk=kkstart,paerlev
         if(.not.lev_interp_comp) then
         do i=1,ncol
            v_coord = pint(i,k)
            if (M_hybi(kk)*Match_ps(i) .lt. v_coord .and. v_coord .le. M_hybi(kk+1)*Match_ps(i)) then
               kupper(i) = kk
               kount = kount + 1
            end if
         end do
!
! If all indices for this level have been found, do the interpolation and
! go to the next level
!
! Interpolate in pressure.
!
         if (kount.eq.ncol) then
            do m=1,naer
               do i=1,ncol
                  dpu = pint(i,k) - M_hybi(kupper(i))*Match_ps(i)
                  dpl = M_hybi(kupper(i)+1)*Match_ps(i) - pint(i,k)
                  aerosol(i,k,m) = &
                     (aerosolc(i,kupper(i)  ,c,m,n)*dpl + &
                     aerosolc(i,kupper(i)+1,c,m,n)*dpu)/(dpl + dpu)
               enddo !i
            end do
            lev_interp_comp = .true.
         end if
         end if
      end do
!
! If we've fallen through the kk=1,levsiz-1 loop, we cannot interpolate and
! must extrapolate from the bottom or top pressure level for at least some
! of the longitude points.
!

      if(.not.lev_interp_comp) then
         do m=1,naer
            do i=1,ncol
               if (pint(i,k) .lt. M_hybi(1)*Match_ps(i)) then
                  aerosol(i,k,m) =  aerosolc(i,1,c,m,n)
               else if (pint(i,k) .gt. M_hybi(paerlev+1)*Match_ps(i)) then
                  aerosol(i,k,m) = 0.0_r8
               else
                  dpu = pint(i,k) - M_hybi(kupper(i))*Match_ps(i)
                  dpl = M_hybi(kupper(i)+1)*Match_ps(i) - pint(i,k)
                  aerosol(i,k,m) = &
                     (aerosolc(i,kupper(i)  ,c,m,n)*dpl + &
                     aerosolc(i,kupper(i)+1,c,m,n)*dpu)/(dpl + dpu)
               end if
            end do
         end do

         if (kount.gt.ncol) then
            call endrun (subname//': Bad data: non-monotonicity suspected in dependent variable')
         end if
      end if
   end do

!   call t_startf ('vi_checks')
!
! aerosol mass beneath lowest interface (pverp) must be 0
!
   aerosol(1:ncol,pverp,:) = 0._r8
!
! Set mass in layer to zero whenever it is less than 
!   1.e-40 kg/m^2 in the layer
!
   do m = 1, naer
      do k = 1, pver
         do i = 1, ncol
            if (aerosol(i,k,m) < 1.e-40_r8) aerosol(i,k,m) = 0._r8
         end do
      end do
   end do
!
! Set mass in layer to zero whenever it is less than 
!   10^-15 relative to column total mass
!
   error_found = .false.
   do m = 1, naer
      do k = 1, pver
         do i = 1, ncol
            AER_diff = aerosol(i,k,m) - aerosol(i,k+1,m)
            if( abs(AER_diff) < 1e-15_r8*aerosol(i,1,m)) then
               AER_diff = 0._r8
            end if
            aerosol_mass(i,k,m)= AER_diff 
            if (aerosol_mass(i,k,m) < 0) error_found = .true.
         end do
      end do
   end do
   if (error_found) then
      do m = 1, naer
         do k = 1, pver
            do i = 1, ncol
               if (aerosol_mass(i,k,m) < 0) then
                  write(iulog,*) subname//': mass < 0, m, col, lev, mass',m, i, k, aerosol_mass(i,k,m)
                  write(iulog,*) subname//': aerosol(k),(k+1)',aerosol(i,k,m),aerosol(i,k+1,m)
                  write(iulog,*) subname//': pint(k+1),(k)',pint(i,k+1),pint(i,k)
                  write(iulog,*)'n,c',n,c
                  call endrun()
               end if
            end do
         end do
      end do
   end if

   call t_stopf ('vert_interpolate')

   return
end subroutine vert_interpolate

!================================================================================================

subroutine aerint (phys_state)

   type(physics_state), intent(in) :: phys_state(begchunk:endchunk)

   integer :: ntmp                                ! used in index swapping
   integer :: start(4)                            ! start vector for netcdf calls
   integer :: kount(4)                            ! count vector for netcdf calls
   integer :: i,j,k                               ! spatial indices
   integer :: m                                   ! constituent index
   integer :: cols, cole
   integer :: lchnk, ncol
   real(r8) :: caldayloc                          ! calendar day of current timestep
   real(r8) :: aerosol_data(naerlon,naerlat,paerlev)    ! aerosol field read in from dataset
   real(r8) :: aerosol_field(naerlon,paerlev+1,naerlat) ! aerosol field to be scattered
   integer latidx,lonidx
   real(r8) closelat,closelon

   character(len=*), parameter :: subname = 'cam3_aero_data.aerint'
   !-----------------------------------------------------------------------

   if (single_column) &
      call shr_scam_GetCloseLatLon(aernid,scmlat,scmlon,closelat,closelon,latidx,lonidx)
 
!
! determine if need to read in next month data
! also determine time interpolation factors
!
   caldayloc = get_curr_calday ()  
!
! If model time is past current forward timeslice, then
! masterproc reads in the next timeslice for time interpolation.  Messy logic is 
! for interpolation between December and January (mo_nxt == 1).  Just like
! ozone_data_timestep_init, sstint.
!
   if (caldayloc > cdayp .and. .not. (mo_nxt == 1 .and. caldayloc >= cdaym)) then
      mo_nxt = mod(mo_nxt,12) + 1
      cdaym = cdayp
      cdayp = Mid(mo_nxt)
!
! Check for valid date info
!
      if (.not. (mo_nxt == 1 .or. caldayloc <= cdayp)) then
         call endrun (subname//': Non-monotonicity suspected in input aerosol data')
      end if

      ntmp = nm
      nm = np
      np = ntmp

      if(aerosol_datan%isncol) then
         do lchnk=begchunk,endchunk
            ncol=phys_state(lchnk)%ncol
            cols=1
            cole=cols+aerosol_datan%count(cols,lchnk)-1
            do while(cole<=ncol)
               start=(/aerosol_datan%start(cols,lchnk),mo_nxt,1,-1/)
               kount=(/aerosol_datan%count(cols,lchnk),1,-1,-1/)
               call handle_ncerr( nf90_get_var(aerosol_datan%ncid, aerosol_datan%psid , &
                    aerosol_datan%ps(cols:cole,lchnk,np), start(1:2), &
                    kount(1:2)),&
                    subname, __LINE__)
               start(2)=1
               start(3)=mo_nxt
               kount(2)=paerlev
               kount(3)=1
               do m=1,naer
                  call handle_ncerr( nf90_get_var(aerosol_datan%ncid, aerosol_datan%dataid(m) , &
                       aerosol_datan%fields(cols:cole,:,lchnk,m,np),  &
                       start(1:3), kount(1:3)),&
                       subname, __LINE__)

               end do
               if(cols==ncol) exit
               cols=cols+aerosol_datan%count(cols,lchnk)
               cole=cols+aerosol_datan%count(cols,lchnk)-1
            end do
         end do
         aerosolc(:,1:paerlev,:,:,np)=aerosol_datan%fields(:,:,:,:,np)
      else
         do m=1,naer
            if (masterproc) then
               if (single_column) then
                  naerlon=1
                  naerlat=1
                  start(:) = (/lonidx,latidx,1,mo_nxt/)
               else
                  start(:) = (/1,1,1,mo_nxt/)
               endif
               kount(:) = (/naerlon,naerlat,paerlev,1/)
               call handle_ncerr( nf90_get_var (aernid, species_id(m), aerosol_data, start, kount),&
                    subname, __LINE__)

               do j=1,naerlat
                  do k=1,paerlev
                     aerosol_field(:,k,j) = aerosol_data(:,j,k)
                  end do
                  aerosol_field(:,paerlev+1,j) = 0._r8   ! value at bottom
               end do
            end if
            call scatter_field_to_chunk (1, paerlev+1, 1, naerlon, aerosol_field, &
                 aerosolc(:,:,:,m,np))
         end do
!
! Retrieve PS from Match
!
         if (masterproc) then
               if (single_column) then
                  naerlon=1
                  naerlat=1
                  start(:) = (/lonidx,latidx,mo_nxt,-1/)
               else
                  start(:) = (/1,1,mo_nxt,-1/)
               endif
               kount(:) = (/naerlon,naerlat,1,-1/)
               call handle_ncerr( nf90_get_var (aernid, Mpsid, M_ps, start, kount),&
                    subname, __LINE__)
               write(iulog,*) subname//': Read aerosols data for julian day', Mid(mo_nxt)
            end if
            call scatter_field_to_chunk (1, 1, 1, naerlon, M_ps(:,:), M_ps_cam_col(:,:,np))
         end if
      end if

end subroutine aerint

end module cam3_aero_data
