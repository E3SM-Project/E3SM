module pdepStreamMod

  !----------------------------------------------------------------------- 
  ! !DESCRIPTION: 
  ! Contains methods for reading in phosphorus deposition data file
  ! Also includes functions for dynamic pdep file handling and 
  ! interpolation.
  ! X.SHI
  ! !USES
  use shr_kind_mod, only: r8 => shr_kind_r8, CL => shr_kind_cl
  use shr_strdata_mod
  use shr_stream_mod
  use shr_string_mod
  use shr_sys_mod
  use shr_mct_mod
  use mct_mod
  use spmdMod     , only: mpicom, masterproc, comp_id, iam
  use clm_varctl  , only: iulog
  use controlMod  , only: NLFilename
  use abortutils  , only: endrun
  use fileutils   , only: getavu, relavu
  use decompMod   , only: bounds_type, ldecomp, gsmap_lnd_gdc2glo 
  use domainMod   , only: ldomain

  ! !PUBLIC TYPES:
  implicit none
  private
  save

  ! !PUBLIC MEMBER FUNCTIONS:
  public :: pdep_init      ! position datasets for dynamic pdep
  public :: pdep_interp    ! interpolates between two years of pdep file data
  public :: clm_domain_mct ! Sets up MCT domain for this resolution

  ! ! PRIVATE TYPES
  type(shr_strdata_type)  :: sdat         ! input data stream
  integer :: stream_year_first_pdep       ! first year in stream to use
  integer :: stream_year_last_pdep        ! last year in stream to use
  integer :: model_year_align_pdep        ! align stream_year_firstpdep with 
  !==============================================================================

contains

  !==============================================================================

  subroutine pdep_init(bounds)
   !    
   ! Initialize data stream information.  
   !
   ! Uses:
   use clm_varctl       , only : inst_name
   use clm_time_manager , only : get_calendar
   use ncdio_pio        , only : pio_subsystem
   use shr_pio_mod      , only : shr_pio_getiotype
   use shr_nl_mod       , only : shr_nl_find_group_name
   use shr_log_mod      , only : errMsg => shr_log_errMsg
   !
   ! arguments
   implicit none
   type(bounds_type), intent(in) :: bounds  
   !
   ! local variables
   integer            :: nu_nml    ! unit for namelist file
   integer            :: nml_error ! namelist i/o error flag
   type(mct_ggrid)    :: dom_clm   ! domain information 
   character(len=CL)  :: stream_fldFileName_pdep
   character(len=CL)  :: pdepmapalgo = 'bilinear'
   character(*), parameter :: shr_strdata_unset = 'NOT_SET'
   character(*), parameter :: subName = "('pdepdyn_init')"
   character(*), parameter :: F00 = "('(pdepdyn_init) ',4a)"
   !-----------------------------------------------------------------------

   namelist /pdepdyn_nml/        &
        stream_year_first_pdep,  &
	stream_year_last_pdep,   &
        model_year_align_pdep,   &
        pdepmapalgo,             &
        stream_fldFileName_pdep

   ! Default values for namelist
    stream_year_first_pdep  = 1                ! first year in stream to use
    stream_year_last_pdep   = 1                ! last  year in stream to use
    model_year_align_pdep   = 1                ! align stream_year_first_pdep with this model year
    stream_fldFileName_pdep = ' '

   ! Read pdepdyn_nml namelist
   if (masterproc) then
      nu_nml = getavu()
      open( nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
      call shr_nl_find_group_name(nu_nml, 'pdepdyn_nml', status=nml_error)
      if (nml_error == 0) then
         read(nu_nml, nml=pdepdyn_nml,iostat=nml_error)
         if (nml_error /= 0) then
            call endrun(msg=' ERROR reading pdepdyn_nml namelist'//errMsg(__FILE__, __LINE__))
         end if
      end if
      close(nu_nml)
      call relavu( nu_nml )
   endif

   call shr_mpi_bcast(stream_year_first_pdep, mpicom)
   call shr_mpi_bcast(stream_year_last_pdep, mpicom)
   call shr_mpi_bcast(model_year_align_pdep, mpicom)
   call shr_mpi_bcast(stream_fldFileName_pdep, mpicom)

   if (masterproc) then
      write(iulog,*) ' '
      write(iulog,*) 'pdepdyn stream settings:'
      write(iulog,*) '  stream_year_first_pdep  = ',stream_year_first_pdep  
      write(iulog,*) '  stream_year_last_pdep   = ',stream_year_last_pdep   
      write(iulog,*) '  model_year_align_pdep   = ',model_year_align_pdep   
      write(iulog,*) '  stream_fldFileName_pdep = ',stream_fldFileName_pdep
      write(iulog,*) ' '
   endif

   call clm_domain_mct (bounds, dom_clm)

   call shr_strdata_create(sdat,name="clmpdep",    &
        pio_subsystem=pio_subsystem,               & 
        pio_iotype=shr_pio_getiotype(inst_name),   &
        mpicom=mpicom, compid=comp_id,             &
        gsmap=gsmap_lnd_gdc2glo, ggrid=dom_clm,    &
        nxg=ldomain%ni, nyg=ldomain%nj,            &
        yearFirst=stream_year_first_pdep,          &
        yearLast=stream_year_last_pdep,            &
        yearAlign=model_year_align_pdep,           &
        offset=0,                                  &
        domFilePath='',                            &
        domFileName=trim(stream_fldFileName_pdep), &
        domTvarName='time',                        &
        domXvarName='lon' ,                        &
        domYvarName='lat' ,                        &  
        domAreaName='area',                        &
        domMaskName='mask',                        &
        filePath='',                               &
        filename=(/trim(stream_fldFileName_pdep)/),&
        fldListFile='PDEP_year',                   &
        fldListModel='PDEP_year',                  &
        fillalgo='none',                           &
        mapalgo=pdepmapalgo,                       &
        calendar=get_calendar(),                   &
	taxmode='extend'                           )

   if (masterproc) then
      call shr_strdata_print(sdat,'CLMPDEP data')
   endif

 end subroutine pdep_init
  
 !================================================================
 subroutine pdep_interp(bounds, atm2lnd_vars)

   !-----------------------------------------------------------------------
   use clm_time_manager, only : get_curr_date, get_days_per_year
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

   call shr_strdata_advance(sdat, mcdate, sec, mpicom, 'pdepdyn')

   ig = 0
   dayspyr = get_days_per_year( )
   do g = bounds%begg,bounds%endg
      ig = ig+1
      atm2lnd_vars%forc_pdep_grc(g) = sdat%avs(1)%rAttr(1,ig) / (secspday * dayspyr)
   end do
   
 end subroutine pdep_interp

 !==============================================================================
  subroutine clm_domain_mct(bounds, dom_clm)

    !-------------------------------------------------------------------
    ! Set domain data type for internal clm grid
    use elm_varcon  , only : re
    use domainMod   , only : ldomain
    use seq_flds_mod
    implicit none
    ! 
    ! arguments
    type(bounds_type), intent(in) :: bounds  
    type(mct_ggrid), intent(out)   :: dom_clm     ! Output domain information for land model
    !
    ! local variables
    integer :: g,i,j              ! index
    integer :: lsize              ! land model domain data size
    real(r8), pointer :: data(:)  ! temporary
    integer , pointer :: idata(:) ! temporary
    !-------------------------------------------------------------------
    !
    ! Initialize mct domain type
    ! lat/lon in degrees,  area in radians^2, mask is 1 (land), 0 (non-land)
    ! Note that in addition land carries around landfrac for the purposes of domain checking
    ! 
    lsize = mct_gsMap_lsize(gsmap_lnd_gdc2glo, mpicom)
    call mct_gGrid_init( GGrid=dom_clm, CoordChars=trim(seq_flds_dom_coord), &
                         OtherChars=trim(seq_flds_dom_other), lsize=lsize )
    !
    ! Allocate memory
    !
    allocate(data(lsize))
    !
    ! Determine global gridpoint number attribute, GlobGridNum, which is set automatically by MCT
    !
    call mct_gsMap_orderedPoints(gsmap_lnd_gdc2glo, iam, idata)
    call mct_gGrid_importIAttr(dom_clm,'GlobGridNum',idata,lsize)
    !
    ! Determine domain (numbering scheme is: West to East and South to North to South pole)
    ! Initialize attribute vector with special value
    !
    data(:) = -9999.0_R8 
    call mct_gGrid_importRAttr(dom_clm,"lat"  ,data,lsize) 
    call mct_gGrid_importRAttr(dom_clm,"lon"  ,data,lsize) 
    call mct_gGrid_importRAttr(dom_clm,"area" ,data,lsize) 
    call mct_gGrid_importRAttr(dom_clm,"aream",data,lsize) 
    data(:) = 0.0_R8     
    call mct_gGrid_importRAttr(dom_clm,"mask" ,data,lsize) 
    !
    ! Determine bounds
    !
    ! Fill in correct values for domain components
    ! Note aream will be filled in in the atm-lnd mapper
    !
    do g = bounds%begg,bounds%endg
       i = 1 + (g - bounds%begg)
       data(i) = ldomain%lonc(g)
    end do
    call mct_gGrid_importRattr(dom_clm,"lon",data,lsize) 

    do g = bounds%begg,bounds%endg
       i = 1 + (g - bounds%begg)
       data(i) = ldomain%latc(g)
    end do
    call mct_gGrid_importRattr(dom_clm,"lat",data,lsize) 

    do g = bounds%begg,bounds%endg
       i = 1 + (g - bounds%begg)
       data(i) = ldomain%area(g)/(re*re)
    end do
    call mct_gGrid_importRattr(dom_clm,"area",data,lsize) 

    do g = bounds%begg,bounds%endg
       i = 1 + (g - bounds%begg)
       data(i) = real(ldomain%mask(g), r8)
    end do
    call mct_gGrid_importRattr(dom_clm,"mask",data,lsize) 

    do g = bounds%begg,bounds%endg
       i = 1 + (g - bounds%begg)
       data(i) = real(ldomain%frac(g), r8)
    end do
    call mct_gGrid_importRattr(dom_clm,"frac",data,lsize) 

    deallocate(data)
    deallocate(idata)

  end subroutine clm_domain_mct
    
end module PdepStreamMod

