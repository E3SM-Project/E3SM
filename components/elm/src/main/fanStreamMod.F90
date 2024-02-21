module FanStreamMod

  !----------------------------------------------------------------------- 
  ! !DESCRIPTION: 
  ! Contains methods for reading in FAN nitrogen deposition (in the form of manure) data file
  ! Also includes functions for dynamic ndep2 file handling and interpolation.
  !
  ! !USES
  use shr_kind_mod, only: r8 => shr_kind_r8, CL => shr_kind_cl
  use shr_strdata_mod
  use shr_stream_mod
  use shr_string_mod
  use shr_sys_mod
  use shr_mct_mod
  use mct_mod
  use spmdMod     , only: mpicom, masterproc, comp_id, iam
  use elm_varctl  , only: iulog
  use abortutils  , only: endrun
  use fileutils   , only: getavu, relavu
  use decompMod   , only: bounds_type, ldecomp, gsmap_lnd_gdc2glo 
  use domainMod   , only: ldomain
  use ndepStreamMod, only: elm_domain_mct

  ! !PUBLIC TYPES:
  implicit none
  private
  save

  ! !PUBLIC MEMBER FUNCTIONS:
  public :: fanstream_init     ! position datasets for dynamic ndep2
  public :: fanstream_interp   ! interpolates between two years of ndep2 file data

  ! ! PRIVATE TYPES
  type(shr_strdata_type)  :: sdat_past, sdat_mix, sdat_urea, sdat_nitr, sdat_soilph         ! input data streams
  integer :: stream_year_first_fan      ! first year in stream to use
  integer :: stream_year_last_fan       ! last year in stream to use
  integer :: model_year_align_fan       ! align stream_year_firstfan with 

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !==============================================================================

contains

  !==============================================================================


  subroutine fanstream_init(bounds, NLFilename)
   !    
   ! Initialize data stream information.  
   !
   ! Uses:
   use elm_varctl       , only : inst_name
   use elm_time_manager , only : get_calendar
   use ncdio_pio        , only : pio_subsystem
   use shr_pio_mod      , only : shr_pio_getiotype
   use shr_nl_mod       , only : shr_nl_find_group_name
   use shr_log_mod      , only : errMsg => shr_log_errMsg
   !
   ! arguments
   implicit none
   type(bounds_type), intent(in) :: bounds  
   character(len=*),  intent(in) :: NLFilename   ! Namelist filename
   !
   ! local variables
   integer            :: nu_nml    ! unit for namelist file
   integer            :: nml_error ! namelist i/o error flag
   type(mct_ggrid)    :: dom_clm   ! domain information 
   character(len=CL)  :: stream_fldFileName_fan
   character(len=CL)  :: fanmapalgo = 'bilinear'
   character(*), parameter :: shr_strdata_unset = 'NOT_SET'
   character(*), parameter :: subName = "('fandyn_init')"
   character(*), parameter :: F00 = "('(fandyn_init) ',4a)"
   !-----------------------------------------------------------------------

   namelist /fan_nml/        &
        stream_year_first_fan,  &
        stream_year_last_fan,   &
        model_year_align_fan,   &
        fanmapalgo,             &
        stream_fldFileName_fan

   ! Default values for namelist
    stream_year_first_fan  = 1                ! first year in stream to use
    stream_year_last_fan   = 1                ! last  year in stream to use
    model_year_align_fan   = 1                ! align stream_year_first_fan with this model year
    stream_fldFileName_fan = ' '

   ! Read fandyn_nml namelist
   if (masterproc) then
      nu_nml = getavu()
      open( nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
      call shr_nl_find_group_name(nu_nml, 'fan_nml', status=nml_error)
      if (nml_error == 0) then
         read(nu_nml, nml=fan_nml,iostat=nml_error)
         if (nml_error /= 0) then
            call endrun(msg=' ERROR reading fan_nml namelist'//errMsg(sourcefile, __LINE__))
         end if
      else
         call endrun(msg=' ERROR finding fan_nml namelist'//errMsg(sourcefile, __LINE__))
      end if
      close(nu_nml)
      call relavu( nu_nml )
   endif

   call shr_mpi_bcast(stream_year_first_fan, mpicom)
   call shr_mpi_bcast(stream_year_last_fan, mpicom)
   call shr_mpi_bcast(model_year_align_fan, mpicom)
   call shr_mpi_bcast(stream_fldFileName_fan, mpicom)

   if (masterproc) then
      write(iulog,*) ' '
      write(iulog,*) 'fan stream settings:'
      write(iulog,*) '  stream_year_first_fan  = ',stream_year_first_fan
      write(iulog,*) '  stream_year_last_fan   = ',stream_year_last_fan   
      write(iulog,*) '  model_year_align_fan   = ',model_year_align_fan   
      write(iulog,*) '  stream_fldFileName_fan = ',stream_fldFileName_fan
      write(iulog,*) ' '
   endif

   call elm_domain_mct (bounds, dom_clm)

   call shr_strdata_create(sdat_past,name="clmfanpast",    &
        pio_subsystem=pio_subsystem,                & 
        pio_iotype=shr_pio_getiotype(inst_name),    &
        mpicom=mpicom, compid=comp_id,              &
        gsmap=gsmap_lnd_gdc2glo, ggrid=dom_clm,     &
        nxg=ldomain%ni, nyg=ldomain%nj,             &
        yearFirst=stream_year_first_fan,          &
        yearLast=stream_year_last_fan,            &
        yearAlign=model_year_align_fan,           &
        offset=0,                                   &
        domFilePath='',                             &
        domFileName=trim(stream_fldFileName_fan), &
        domTvarName='time',                         &
        domXvarName='x' ,                           &
        domYvarName='y' ,                           &  
        domAreaName='area',                         &
        domMaskName='mask',                         &
        filePath='',                                &
        filename=(/trim(stream_fldFileName_fan)/),&
        fldListFile='Nmanure_pastures',             &
        fldListModel='Nmanure_pastures',            &
        fillalgo='none',                            &
        mapalgo=fanmapalgo,                       &
        calendar=get_calendar(),                    &
	taxmode='extend'                            )

   if (masterproc) then
      call shr_strdata_print(sdat_past,'CLMFAN data')
   endif

   call shr_strdata_create(sdat_mix,name="clmfanmixed",    &
        pio_subsystem=pio_subsystem,                & 
        pio_iotype=shr_pio_getiotype(inst_name),    &
        mpicom=mpicom, compid=comp_id,              &
        gsmap=gsmap_lnd_gdc2glo, ggrid=dom_clm,     &
        nxg=ldomain%ni, nyg=ldomain%nj,             &
        yearFirst=stream_year_first_fan,          &
        yearLast=stream_year_last_fan,            &
        yearAlign=model_year_align_fan,           &
        offset=0,                                   &
        domFilePath='',                             &
        domFileName=trim(stream_fldFileName_fan), &
        domTvarName='time',                         &
        domXvarName='x' ,                           &
        domYvarName='y' ,                           &  
        domAreaName='area',                         &
        domMaskName='mask',                         &
        filePath='',                                &
        filename=(/trim(stream_fldFileName_fan)/),&
        fldListFile='Nmanure_mixed',                &
        fldListModel='Nmanure_mixed',               &
        fillalgo='none',                            &
        mapalgo=fanmapalgo,                       &
        calendar=get_calendar(),                    &
	taxmode='extend'                            )

   if (masterproc) then
      call shr_strdata_print(sdat_mix,'CLMFAN data')
   endif

   call shr_strdata_create(sdat_urea,name="clmfanurea",    &
        pio_subsystem=pio_subsystem,                & 
        pio_iotype=shr_pio_getiotype(inst_name),    &
        mpicom=mpicom, compid=comp_id,              &
        gsmap=gsmap_lnd_gdc2glo, ggrid=dom_clm,     &
        nxg=ldomain%ni, nyg=ldomain%nj,             &
        yearFirst=stream_year_first_fan,          &
        yearLast=stream_year_last_fan,            &
        yearAlign=model_year_align_fan,           &
        offset=0,                                   &
        domFilePath='',                             &
        domFileName=trim(stream_fldFileName_fan), &
        domTvarName='time',                         &
        domXvarName='x' ,                           &
        domYvarName='y' ,                           &  
        domAreaName='area',                         &
        domMaskName='mask',                         &
        filePath='',                                &
        filename=(/trim(stream_fldFileName_fan)/),&
        fldListFile='fract_urea',                &
        fldListModel='fract_urea',               &
        fillalgo='none',                            &
        mapalgo='nn',                       &
        calendar=get_calendar(),                    &
	taxmode='extend'                            )

   if (masterproc) then
      call shr_strdata_print(sdat_urea,'CLMFAN data')
   endif

   call shr_strdata_create(sdat_nitr,name="clmfannitr",    &
        pio_subsystem=pio_subsystem,                & 
        pio_iotype=shr_pio_getiotype(inst_name),    &
        mpicom=mpicom, compid=comp_id,              &
        gsmap=gsmap_lnd_gdc2glo, ggrid=dom_clm,     &
        nxg=ldomain%ni, nyg=ldomain%nj,             &
        yearFirst=stream_year_first_fan,          &
        yearLast=stream_year_last_fan,            &
        yearAlign=model_year_align_fan,           &
        offset=0,                                   &
        domFilePath='',                             &
        domFileName=trim(stream_fldFileName_fan), &
        domTvarName='time',                         &
        domXvarName='x' ,                           &
        domYvarName='y' ,                           &  
        domAreaName='area',                         &
        domMaskName='mask',                         &
        filePath='',                                &
        filename=(/trim(stream_fldFileName_fan)/),&
        fldListFile='fract_nitr',                &
        fldListModel='fract_nitr',               &
        fillalgo='none',                            &
        mapalgo='nn',                       &
        calendar=get_calendar(),                    &
	taxmode='extend'                            )

   call shr_strdata_create(sdat_soilph,name="clmfanph",    &
        pio_subsystem=pio_subsystem,                & 
        pio_iotype=shr_pio_getiotype(inst_name),    &
        mpicom=mpicom, compid=comp_id,              &
        gsmap=gsmap_lnd_gdc2glo, ggrid=dom_clm,     &
        nxg=ldomain%ni, nyg=ldomain%nj,             &
        yearFirst=stream_year_first_fan,          &
        yearLast=stream_year_last_fan,            &
        yearAlign=model_year_align_fan,           &
        offset=0,                                   &
        domFilePath='',                             &
        domFileName=trim(stream_fldFileName_fan), &
        domTvarName='time',                         &
        domXvarName='x' ,                           &
        domYvarName='y' ,                           &  
        domAreaName='area',                         &
        domMaskName='mask',                         &
        filePath='',                                &
        filename=(/trim(stream_fldFileName_fan)/),&
        fldListFile='soilph',                &
        fldListModel='fansoilph',               &
        fillalgo='none',                            &
        mapalgo='nn',                       &
        calendar=get_calendar(),                    &
	taxmode='extend'                            )

   if (masterproc) then
      call shr_strdata_print(sdat_mix,'CLMFAN data')
   endif

   
 end subroutine fanstream_init
  
 !================================================================
 subroutine fanstream_interp(bounds, atm2lnd_vars)

   !-----------------------------------------------------------------------
   use elm_time_manager, only : get_curr_date, get_days_per_year
   use elm_varcon      , only : secspday
   use atm2lndType     , only : atm2lnd_type
   !
   ! Arguments
   type(bounds_type) , intent(in)    :: bounds  
   type(atm2lnd_type), intent(inout) :: atm2lnd_vars
   !
   ! Local variables
   integer :: g, ig 
   integer :: year    ! year (0, ...) for nstep+1
   integer :: mon     ! month (1, ..., 12) for nstep+1
   integer :: day     ! day of month (1, ..., 31) for nstep+1
   integer :: sec     ! seconds into current date for nstep+1
   integer :: mcdate  ! Current model date (yyyymmdd)
   integer :: dayspyr ! days per year
   !-----------------------------------------------------------------------

   call get_curr_date(year, mon, day, sec)
   mcdate = year*10000 + mon*100 + day
   dayspyr = get_days_per_year( )

   call shr_strdata_advance(sdat_past, mcdate, sec, mpicom, 'clmfanpasture')

   ig = 0
   do g = bounds%begg,bounds%endg
      ig = ig+1
      atm2lnd_vars%forc_ndep_past_grc(g) = sdat_past%avs(1)%rAttr(1,ig) / (secspday * dayspyr)
   end do

   call shr_strdata_advance(sdat_mix, mcdate, sec, mpicom, 'clmfanmixed')

   ig = 0
   do g = bounds%begg,bounds%endg
      ig = ig+1
      atm2lnd_vars%forc_ndep_mgrz_grc(g) = sdat_mix%avs(1)%rAttr(1,ig) / (secspday * dayspyr)
   end do

   call shr_strdata_advance(sdat_urea, mcdate, sec, mpicom, 'clmfanurea')

   ig = 0
   do g = bounds%begg,bounds%endg
      ig = ig+1
      atm2lnd_vars%forc_ndep_urea_grc(g) = sdat_urea%avs(1)%rAttr(1,ig)
   end do

   call shr_strdata_advance(sdat_nitr, mcdate, sec, mpicom, 'clmfannitr')

   ig = 0
   do g = bounds%begg,bounds%endg
      ig = ig+1
      atm2lnd_vars%forc_ndep_nitr_grc(g) = sdat_nitr%avs(1)%rAttr(1,ig)
   end do

   call shr_strdata_advance(sdat_soilph, mcdate, sec, mpicom, 'clmfanph')

   ig = 0
   do g = bounds%begg,bounds%endg
      ig = ig+1
      atm2lnd_vars%forc_soilph_grc(g) = sdat_soilph%avs(1)%rAttr(1,ig)
   end do

   
 end subroutine fanstream_interp
    
end module FanStreamMod

