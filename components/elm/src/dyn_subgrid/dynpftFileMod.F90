module dynpftFileMod

#include "shr_assert.h"

  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Handle reading of the pftdyn dataset, which specifies transient areas of natural Patches
  !
  ! !USES:
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use shr_log_mod         , only : errMsg => shr_log_errMsg
  use decompMod           , only : bounds_type, BOUNDS_LEVEL_PROC
  use dynFileMod          , only : dyn_file_type
  use dynVarTimeInterpMod , only : dyn_var_time_interp_type
  use clm_varctl          , only : iulog
  use abortutils          , only : endrun
  use spmdMod             , only : masterproc, mpicom
  use elm_varcon          , only : grlnd, nameg
  use LandunitType        , only : lun_pp                
  !DW  not use at all     !   use ColumnType          , only : col                
  use VegetationType           , only : veg_pp                
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  private
  save
  public :: dynpft_init     ! initialize information read from pftdyn dataset
  public :: dynpft_interp   ! interpolate pftdyn information to current time step
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: dynpft_check_consistency   ! check consistency with surface dataset
  private :: dynpft_read_consistency_nl ! read namelist associated with consistency checks
  !
  ! ! PRIVATE TYPES
  type(dyn_file_type), target    :: dynpft_file   ! information for the pftdyn file
  type(dyn_var_time_interp_type) :: wtpatch       ! weight of each patch relative to the natural veg landunit

  character(len=*), parameter    :: varname = 'PCT_NAT_PFT'  ! name of variable on file
  !---------------------------------------------------------------------------

contains


  !-----------------------------------------------------------------------
  subroutine dynpft_init(bounds, dynpft_filename)
    !
    ! !DESCRIPTION:
    ! Initialize dynamic pft dataset (position it to the right time samples
    ! that bound the initial model date)
    !
    ! This also calls dynpft_interp for the initial time
    !
    ! !USES:
    use clm_varpar     , only : numpft, maxpatch_pft, natpft_size
    use dynTimeInfoMod , only : YEAR_POSITION_START_OF_TIMESTEP
    use dynTimeInfoMod , only : YEAR_POSITION_END_OF_TIMESTEP
    use ncdio_pio
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds          ! proc-level bounds
    character(len=*) , intent(in) :: dynpft_filename ! name of file containing transient pft information
    !
    ! !LOCAL VARIABLES:
    integer  :: wtpatch_shape(2)                  ! shape of the wtpatch data

    character(len= 32)     :: subname='dynpft_init'! subroutine name
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL(bounds%level == BOUNDS_LEVEL_PROC, subname // ': argument must be PROC-level bounds')

    ! Error check

    if ( maxpatch_pft /= numpft+1 )then
       call endrun(msg=' maxpatch_pft does NOT equal numpft+1 -- this is invalid for dynamic PFT case'//&
            errMsg(__FILE__, __LINE__) )
    end if

    if (masterproc) then
       write(iulog,*) 'Attempting to read pft dynamic landuse data .....'
    end if

    !dynpft_file = dyn_file_type(dynpft_filename, YEAR_POSITION_START_OF_TIMESTEP)
    dynpft_file = dyn_file_type(dynpft_filename, YEAR_POSITION_END_OF_TIMESTEP)

    ! Consistency checks
    call check_dim(dynpft_file, 'natpft', natpft_size)
    call dynpft_check_consistency(bounds)

    ! read data PCT_NAT_PFT corresponding to correct year
    !
    ! Note: if you want to change PCT_NAT_PFT so that it is NOT interpolated, but instead
    ! jumps to each year's value on Jan 1 of that year, simply change wtpatch to be of type
    ! dyn_var_time_uninterp_type (rather than dyn_var_time_interp_type), and change the
    ! following constructor to construct a variable of dyn_var_time_uninterp_type. That's
    ! all you need to do.

    wtpatch_shape = [(bounds%endg-bounds%begg+1), natpft_size]
    wtpatch = dyn_var_time_interp_type( &
         dyn_file=dynpft_file, varname=varname, &
         dim1name=grlnd, conversion_factor=100._r8, &
         do_check_sums_equal_1=.true., data_shape=wtpatch_shape)

    call dynpft_interp(bounds)

  end subroutine dynpft_init

  !-----------------------------------------------------------------------
  subroutine dynpft_check_consistency(bounds)
    !
    ! !DESCRIPTION:
    ! Check consistency between dynpft file and surface dataset.
    !
    ! This is done by assuming that PCT_NAT_PFT at time 1 in the pftdyn file agrees with
    ! PCT_NAT_PFT on the surface dataset.
    !
    ! !USES:
    use clm_varsur, only : wt_nat_patch
    use clm_varpar, only : natpft_size
    use ncdio_pio
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds  ! proc-level bounds
    !
    ! !LOCAL VARIABLES:
    logical             :: check_dynpft_consistency ! whether to do the consistency check in this routine
    integer             :: g                        ! index
    real(r8), pointer   :: wtpatch_time1(:,:)       ! weight of each pft in each grid cell at first time
    logical             :: readvar                  ! whether variable was read
    real(r8), parameter :: tol = 1.e-13_r8          ! tolerance for checking equality

    character(len=*), parameter :: subname = 'dynpft_check_consistency'
    !-----------------------------------------------------------------------
    
    call dynpft_read_consistency_nl(check_dynpft_consistency)

    if (check_dynpft_consistency) then

       ! Read first time slice of PCT_NAT_PATCH

       allocate(wtpatch_time1(bounds%begg:bounds%endg, natpft_size))
       call ncd_io(ncid=dynpft_file, varname=varname, flag='read', data=wtpatch_time1, &
            dim1name=grlnd, nt=1, readvar=readvar)
       if (.not. readvar) then
          call endrun(msg=' ERROR: ' // trim(varname) // ' NOT on landuse_timeseries file'//&
               errMsg(__FILE__, __LINE__))
       end if

       ! Convert from PCT to weight on grid cell
       wtpatch_time1(bounds%begg:bounds%endg,:) = wtpatch_time1(bounds%begg:bounds%endg,:) / 100._r8
    
       ! Compare with values read from surface dataset
       do g = bounds%begg, bounds%endg
          if (any(abs(wtpatch_time1(g,:) - wt_nat_patch(g,:)) > tol)) then
             write(iulog,*) subname//' mismatch between PCT_NAT_PATCH at initial time and that obtained from surface dataset'
             write(iulog,*) 'On landuse_timeseries file: ', wtpatch_time1(g,:)
             write(iulog,*) 'On surface dataset: ', wt_nat_patch(g,:)
             write(iulog,*) ' '
             write(iulog,*) 'Confirm that the year of your surface dataset'
             write(iulog,*) 'corresponds to the first year of your landuse_timeseries file'
             write(iulog,*) '(e.g., for a landuse_timeseries file starting at year 1850, which is typical,'
             write(iulog,*) 'you should be using an 1850 surface dataset),'
             write(iulog,*) 'and that your landuse_timeseries file is compatible with the surface dataset.'
             write(iulog,*) ' '
             write(iulog,*) 'If you are confident that you are using the correct landuse_timeseries file'
             write(iulog,*) 'and the correct surface dataset, then you can bypass this check by setting:'
             write(iulog,*) '  check_dynpft_consistency = .false.'
             write(iulog,*) 'in user_nl_clm'
             write(iulog,*) ' '
             call endrun(decomp_index=g, clmlevel=nameg, msg=errMsg(__FILE__, __LINE__))
          end if
       end do

       deallocate(wtpatch_time1)

    end if

  end subroutine dynpft_check_consistency

  !-----------------------------------------------------------------------
  subroutine dynpft_read_consistency_nl(check_dynpft_consistency)
    !
    ! !DESCRIPTION:
    ! Read namelist settings related to pftdyn consistency checks
    !
    ! !USES:
    use fileutils      , only : getavu, relavu
    use clm_nlUtilsMod , only : find_nlgroup_name
    use controlMod     , only : NLFilename
    use shr_mpi_mod    , only : shr_mpi_bcast
    !
    ! !ARGUMENTS:
    logical, intent(out) :: check_dynpft_consistency ! whether to do the consistency check
    !
    ! !LOCAL VARIABLES:
    integer :: nu_nml    ! unit for namelist file
    integer :: nml_error ! namelist i/o error flag

    character(len=*), parameter :: subname = 'dynpft_read_consistency_nl'
    !-----------------------------------------------------------------------

    namelist /dynpft_consistency_checks/ &
         check_dynpft_consistency

    ! Set default namelist values
    check_dynpft_consistency = .true.

    ! Read namelist
    if (masterproc) then
       nu_nml = getavu()
       open( nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
       call find_nlgroup_name(nu_nml, 'dynpft_consistency_checks', status=nml_error)
       if (nml_error == 0) then
          read(nu_nml, nml=dynpft_consistency_checks,iostat=nml_error)
          if (nml_error /= 0) then
             call endrun(msg='ERROR reading dynpft_consistency_checks namelist'//errMsg(__FILE__, __LINE__))
          end if
       end if
       close(nu_nml)
       call relavu( nu_nml )
    endif

    call shr_mpi_bcast (check_dynpft_consistency, mpicom)

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) 'dynpft_consistency_checks settings:'
       write(iulog,nml=dynpft_consistency_checks)
       write(iulog,*) ' '
    end if

  end subroutine dynpft_read_consistency_nl



  !-----------------------------------------------------------------------
  subroutine dynpft_interp(bounds)
    !
    ! !DESCRIPTION:
    ! Time interpolate dynamic pft data to get pft weights for model time
    !
    ! !USES:
    use landunit_varcon , only : istsoil
    use clm_varpar      , only : natpft_lb, natpft_ub
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds  ! proc-level bounds
    !
    ! !LOCAL VARIABLES:
    integer               :: m,p,l,g          ! indices
    real(r8), allocatable :: wtpatch_cur(:,:)   ! current pft weights
    character(len=32) :: subname='dynpft_interp' ! subroutine name
    !-----------------------------------------------------------------------

    ! Interpolate pctpft to current time step - output in pctpft_intp
    ! Map interpolated pctpft to subgrid weights
    ! assumes that maxpatch_pft = numpft + 1, that each landunit has only 1 column, 
    ! SCAM has not been defined, and create_croplandunit = .false.

    SHR_ASSERT_ALL(bounds%level == BOUNDS_LEVEL_PROC, subname // ': argument must be PROC-level bounds')

    ! Get pft weights for this time step

    ! As a workaround for an internal compiler error with ifort 13.1.2 on goldbach, call
    ! the specific name of this procedure rather than using its generic name
    call dynpft_file%time_info%set_current_year_get_year()

    allocate(wtpatch_cur(bounds%begg:bounds%endg, natpft_lb:natpft_ub))
    call wtpatch%get_current_data(wtpatch_cur)

    do p = bounds%begp,bounds%endp
       g = veg_pp%gridcell(p)
       l = veg_pp%landunit(p)

       ! Note that we only deal with the istsoil landunit here, NOT the istcrop landunit
       ! (if there is one)
       ! (However, currently [as of 5-9-13] the code won't let you run with transient
       ! Patches combined with create_crop_landunit anyway, so it's a moot point.)
       if (lun_pp%itype(l) == istsoil) then
          m = veg_pp%itype(p)

          ! Note that the following assignment assumes that all Patches share a single column
          veg_pp%wtcol(p) = wtpatch_cur(g,m)
       end if

    end do

    deallocate(wtpatch_cur)

  end subroutine dynpft_interp

end module dynpftFileMod
