subroutine read_patch_background_sounding(SoundingFileName, ptop_model, &
     npatch_start, npatch_end, nzsnd, psnd, tsnd, qsnd, masterproc)
  use netcdf
  
  ! inputs
  character(LEN=250), intent(in) :: SoundingFileName
  real, intent(in) :: ptop_model ! pressure at top of model in mb
  integer, intent(in) :: nzsnd ! max length of psnd, qsnd, tsnd
  logical :: masterproc ! true if MPI rank==0

  ! outputs
  integer, intent(out) :: npatch_start, npatch_end ! indices of start/end of
                                           ! patched sounding in psnd, qsnd, tsnd
  real, intent(out) :: psnd(nzsnd) ! pressure sounding from SoundingFileName, mb
  real, intent(out) :: tsnd(nzsnd) ! temperature sounding from SoundingFileName, K
  real, intent(out) :: qsnd(nzsnd) ! water vapor sounding from SoundingFileName, kg/kg

  !local variables
  integer :: nlev
  integer :: status(11)
  integer :: k
  integer :: ncid, dimIDp, dimIDt, varID
  real, allocatable, dimension(:) :: psnd_in
  real, allocatable, dimension(:,:,:,:) :: tsnd_in, qsnd_in

  character(LEN=nf90_max_name) :: tmpName

  ! Read profiles from SCAM IOP data file
  status(:)   = nf90_NoErr
  status(1)   = nf90_open(TRIM(SoundingFileName),nf90_nowrite,ncid)
	
  ! get number of pressure levels
  status(2)   = nf90_inq_dimid(ncid,"lev",dimIDp)
  status(3)   = nf90_inquire_dimension(ncid, dimIDp, tmpName, nlev)

  ! get number of time levels
  status(4)   = nf90_inq_dimid(ncid,"time",dimIDt)
  if(status(4).ne.nf90_NoErr) then
    status(4)   = nf90_inq_dimid(ncid,"tsec",dimIDt)
  end if
  status(5) = nf90_inquire_dimension(ncid, dimIDt, tmpName, ntime)

  if(nlev.gt.nzsnd) then
    write(*,*) 'ERROR in read_patch_background_sounding.f90:'
    write(*,*) '***** number of levels in ', TRIM(SoundingFileName)
    write(*,*) '***** exceeds nzsnd.  Reset nzsnd in rad_full.f90'
    call rad_error()
  end if
    
  ALLOCATE(psnd_in(nlev), tsnd_in(1,1,nlev,ntime), qsnd_in(1,1,nlev,ntime), &
       STAT=ierr)
  if(ierr.ne.0) then
    write(*,*) 'ERROR in read_patch_background_sounding.f90:'
    write(*,*) '***** Could not allocate arrays to read in sounding'
    call rad_error()
  end if

  ! get pressure levels (in Pascal)
  status(6)   = nf90_inq_varid(ncid,"lev",varID)
  status(7)   = nf90_get_var(ncid, varID, psnd_in)

  ! get temperature and moisture (K and kg/kg, respectively)
  status(8)   = nf90_inq_varid(ncid,"T",varID)
  status(9)   = nf90_get_var(ncid, varID, tsnd_in)
  status(10)   = nf90_inq_varid(ncid,"q",varID)
  status(11)   = nf90_get_var(ncid, varID, qsnd_in)

  psnd(:) = 0.
  tsnd(:) = 0.
  qsnd(:) = 0.
  ! reverse order from top-down to bottom-up.
  psnd(1:nlev) = psnd_in(nlev:1:-1)/100. ! convert from Pa to hPa
  tsnd(1:nlev) = tsnd_in(1,1,nlev:1:-1,1)
  qsnd(1:nlev) = qsnd_in(1,1,nlev:1:-1,1)
  
  ! find whether we need this sounding to be patched on top of model.
  npatch_start = nlev+1
  npatch_end = nlev+1
  do k = 1,nlev
    if(psnd(k).lt.ptop_model - 10.) then
      ! start patch ~10hPa above model top.
      npatch_start = k
      EXIT
    end if
  end do

  if(npatch_start.le.nlev) then
    npatch_end = nlev
  end if
  
  if(masterproc) then
    write(*,*)
    write(*,*) 'Background sounding, p (mb), T (K), q (kg/kg)'
    do k = 1,nlev
      if(k.eq.npatch_start) write(*,*) '**** Start Patch Here *****'
      write(*,998) psnd(k), tsnd(k), qsnd(k)
998   format(f8.3,f8.3,e12.4)
      if(k.eq.npatch_end) write(*,*) '**** End Patch Here *****'
    end do
    if(npatch_start.gt.nlev) write(*,*) '**** No patching required -- model top is deeper than sounding ****'
    write(*,*)
  end if

  DEALLOCATE(psnd_in, tsnd_in, qsnd_in, STAT=ierr)
  if(ierr.ne.0) then
    write(*,*) 'ERROR in read_patch_background_sounding.f90:'
    write(*,*) '***** Could not allocate arrays to read in sounding'
    call rad_error()
  end if

end subroutine read_patch_background_sounding
